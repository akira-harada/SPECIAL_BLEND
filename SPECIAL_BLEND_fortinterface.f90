program SPECIAL_BLEND
  implicit none
  real*8 :: params(15)
  integer :: datsize,linenum,nparam
  integer :: ana_mode ! 1: unbinned, 2: full-binned, 3: time-binned, 4: Gaussian approximation
  integer :: tbinnumber
  integer, parameter :: ebinnumber = 30
  real*8 :: gbeta,dist,Mdet
  real*8 :: mmin,mmax,rmin,rmax,emin,emax
  real*8 :: dt_ini,tmin,tmax
  real*8 :: t_read,e_read
  real*8,allocatable :: tdata(:),edata(:)
  real*8,allocatable :: t_bin(:),dt(:)
  real*8,allocatable :: hist(:,:)
  real*8,allocatable :: thist(:),e_ave(:)
  real*8 :: e_bin(ebinnumber),de(ebinnumber)
  real*8,allocatable :: mass(:),rad(:),et(:)
  real*8,allocatable :: mlogLH(:,:,:) ! minus log likelihood
  real*8,allocatable :: LH_MR(:,:),LH_RE(:,:),LH_EM(:,:)
  real*8,allocatable :: LH_M(:),LH_R(:),LH_E(:)
  real*8 :: MR95,MR68,peakMR(2),RE95,RE68,peakRE(2),EM95,EM68,peakEM(2)
  real*8 :: CIandBFM(5),CIandBFR(5),CIandBFE(5),lower95,lower68,bestfit,upper68,upper95
  integer :: i,j,status
  
  real*8 :: time_i,time_f ! for test
  
  call cpu_time(time_i)
  
  open(11,file='parameters.dat',status='old')
  do i = 1, 15
    read(11,*)params(i)
  enddo
  close(11)
  gbeta      = params(1)
  dist       = params(2)
  Mdet       = params(3)
  nparam     = int(params(4))
  mmin       = params(5)
  mmax       = params(6)
  rmin       = params(7)
  rmax       = params(8)
  emin       = params(9)
  emax       = params(10)
  ana_mode   = int(params(11))
  tbinnumber = int(params(12))
  dt_ini     = params(13)
  tmin       = params(14)
  tmax       = params(15)
  
  allocate(mass(nparam),rad(nparam),et(nparam))
  allocate(mlogLH(nparam,nparam,nparam))
  allocate(LH_MR(nparam,nparam),LH_RE(nparam,nparam),LH_EM(nparam,nparam))
  allocate(LH_M(nparam),LH_R(nparam),LH_E(nparam))
  allocate(t_bin(tbinnumber),dt(tbinnumber),hist(tbinnumber,ebinnumber),thist(tbinnumber),e_ave(tbinnumber))

  call grid(mass, nparam, mmax, mmin)
  call grid( rad, nparam, rmax, rmin)
  call grid(  et, nparam, emax, emin)
  
  i=0
  j=0
  open(12,file='time_energy.dat',status='old')
  do while(1.eq.1)
    read(12,*,end=100)t_read,e_read
    if(t_read .ge. tmin .and. t_read .le. tmax) i = i+1
    j = j + 1
  enddo
100 continue
    datsize = i
    linenum = j
  close(12)
  
  allocate(tdata(datsize),edata(datsize))
  
  open(12,file='time_energy.dat',status='old')
  i = 1
  do j = 1, linenum
    read(12,*)t_read,e_read
    if(t_read .ge. tmin .and. t_read .le. tmax) then
      tdata(i) = t_read
      edata(i) = e_read
      i = i + 1
    endif
  enddo
  close(12)
  
  if(ana_mode .eq. 1) then
    write(*,*)'unbinned analysis mode'
    write(*,*)'used event number',datsize,'/total event number',linenum
    write(*,*)linenum - datsize,'events are outside [tmin,tmax] and neglected'
    call eval_unbinned_likelihood(tdata,edata,mass,rad,et,gbeta,dist,Mdet,mlogLH,datsize,nparam)
  elseif(ana_mode .eq. 2) then
    write(*,*)'full-binned analysis mode'
    call binning(tdata,edata,dt_ini,tmin,tmax,t_bin,dt,e_bin,de,hist,status,datsize,tbinnumber)
    write(*,*)'total event number in histogram',int(sum(hist)),'/true total number',linenum
    write(*,*)linenum - int(sum(hist)),'events are outside [tmin,tmax] and neglected'
    if(status .ge. 1) then
      write(*,*)'binning error, try another bin number'
      stop
    endif
    call eval_binned_likelihood(t_bin,dt,e_bin,de,hist,mass,rad,et,gbeta,dist,Mdet,mlogLH,tbinnumber,nparam)
  elseif(ana_mode .eq. 3) then
    write(*,*)'time-binned analysis mode'
    call time_binning(tdata,dt_ini,tmin,tmax,t_bin,dt,thist,status,datsize,tbinnumber)
    write(*,*)'total event number in histogram',int(sum(thist)),'/true total number',linenum
    write(*,*)linenum - int(sum(thist)),'events are outside [tmin,tmax] and neglected'
    if(status .ge. 1) then
      write(*,*)'binning error, try another bin number'
      stop
    endif
    call eval_timebinned_likelihood(t_bin,dt,thist,edata,mass,rad,et,gbeta,dist,Mdet,mlogLH,tbinnumber,datsize,nparam)
  elseif(ana_mode .eq. 4) then
    write(*,*)'Gaussian-likelihood analysis mode'
    call time_binning(tdata,dt_ini,tmin,tmax,t_bin,dt,thist,status,datsize,tbinnumber)
    write(*,*)'total event number in histogram',int(sum(thist)),'/true total number',linenum
    write(*,*)linenum - int(sum(thist)),'events are outside [tmin,tmax] and neglected'
    if(status .ge. 1) then
      write(*,*)'binning error, try another bin number'
      stop
    endif
    call event_energy_averaging(thist,edata,e_ave,tbinnumber,datsize)
    call eval_Gaussian_likelihood(t_bin,dt,thist,e_ave,mass,rad,et,gbeta,dist,Mdet,mlogLH,tbinnumber,nparam)
  endif
  
  call mr_marginalize(mass,rad,et,mlogLH,LH_MR,MR95,MR68,peakMR,nparam)
  call re_marginalize(mass,rad,et,mlogLH,LH_RE,RE95,RE68,peakRE,nparam)
  call em_marginalize(mass,rad,et,mlogLH,LH_EM,EM95,EM68,peakEM,nparam)
  call  m_marginalize(mass,rad,et,mlogLH,LH_M, CIandBFM, nparam)
  call  r_marginalize(mass,rad,et,mlogLH,LH_R, CIandBFR, nparam)
  call  e_marginalize(mass,rad,et,mlogLH,LH_E, CIandBFE, nparam)

  open(20,file='LH_MR.dat',status='replace')
  write(20,*)'# CI95% level =',MR95,', CI68% level =',MR68
  write(20,*)'# peak at (',peakMR(1),',',peakMR(2),')'
  do i=1,nparam
    do j=1,nparam
      write(20,'(99E25.15e3)')mass(i),rad(j),LH_MR(i,j)
    enddo
    write(20,*)
  enddo
  close(20)
  
  open(21,file='LH_RE.dat',status='replace')
  write(21,*)'# CI95% level =',RE95,', CI68% level =',RE68
  write(21,*)'# peak at (',peakRE(1),',',peakRE(2),')'
  do i=1,nparam
    do j=1,nparam
      write(21,'(99E25.15e3)')rad(i),et(j),LH_RE(i,j)
    enddo
    write(21,*)
  enddo
  close(21)
  
  open(22,file='LH_EM.dat',status='replace')
  write(22,*)'# CI95% level =',EM95,', CI68% level =',EM68
  write(22,*)'# peak at (',peakEM(1),',',peakEM(2),')'
  do i=1,nparam
    do j=1,nparam
      write(22,'(99E25.15e3)')et(i),mass(j),LH_EM(i,j)
    enddo
    write(22,*)
  enddo
  close(22)
  
  open(30,file='LH_M.dat',status='replace')
  lower95 = CIandBFM(1)
  lower68 = CIandBFM(2)
  bestfit = CIandBFM(3)
  upper68 = CIandBFM(4)
  upper95 = CIandBFM(5)
  write(30,*)'#mass best fit:',bestfit,'+',upper68-bestfit,'/-',bestfit-lower68,'(68%) +',&
            &                              upper95-bestfit,'/-',bestfit-lower95,'(95%)'
  do i=1,nparam
    write(30,'(99E25.15e3)')mass(i),LH_M(i)
  enddo
  close(30)

  open(31,file='LH_R.dat',status='replace')
  lower95 = CIandBFR(1)
  lower68 = CIandBFR(2)
  bestfit = CIandBFR(3)
  upper68 = CIandBFR(4)
  upper95 = CIandBFR(5)
  write(31,*)'#radius best fit:',bestfit,'+',upper68-bestfit,'/-',bestfit-lower68,'(68%) +',&
            &                                upper95-bestfit,'/-',bestfit-lower95,'(95%)'
  do i=1,nparam
    write(31,'(99E25.15e3)')rad(i),LH_R(i)
  enddo
  close(31)
  
  open(32,file='LH_E.dat',status='replace')
  lower95 = CIandBFE(1)
  lower68 = CIandBFE(2)
  bestfit = CIandBFE(3)
  upper68 = CIandBFE(4)
  upper95 = CIandBFE(5)
  write(32,*)'#energy best fit:',bestfit,'+',upper68-bestfit,'/-',bestfit-lower68,'(68%) +',&
            &                                upper95-bestfit,'/-',bestfit-lower95,'(95%)'
  do i=1,nparam
    write(32,'(99E25.15e3)')et(i),LH_E(i)
  enddo
  close(32)
  
  deallocate(mass,rad,et)
  deallocate(mlogLH)
  deallocate(LH_MR,LH_RE,LH_EM)
  deallocate(LH_M,LH_R,LH_E)
  deallocate(tdata,edata)
  deallocate(t_bin,dt,hist,thist,e_ave)
  
  call cpu_time(time_f)
  write(*,*)'elapse time=',time_f-time_i,'sec'

  stop
end program SPECIAL_BLEND

!!!!! parameter search grid !!!!!
subroutine grid(vec,n,max,min)
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: max, min
  integer :: i
  real*8, dimension(n) :: vec

  do i=1,n
    vec(i) = min+dble(i-1)/dble(n-1)*(max-min) 
  enddo
      
end subroutine grid



