!!!!! unbinned analysis !!!!!
!!!!! subroutine to give mlogLH = -log(likelihood) !!!!!
subroutine eval_unbinned_likelihood(tdata,edata,mass,rad,et,gbeta,dist,Mdet,mlogLH,datsize,nparam)
  implicit none
  real*8, dimension(datsize), intent(in) :: tdata, edata
  real*8, dimension(nparam), intent(in) :: mass,rad,et
  real*8, intent(in) :: gbeta,dist,Mdet
  real*8, dimension(nparam,nparam,nparam), intent(out) :: mlogLH
  integer, intent(in) :: datsize,nparam
  integer :: i, j, k, p
  real*8 :: anarate_coef, rate_guess, e_ave_guess, t0, temp_guess, integral
  real*8, parameter :: FDint4=23.3d0,FDint5=118.d0

  do i=1,nparam
    do j=1,nparam
      do k=1,nparam
        t0            = 210.d0*(mass(i)/1.4d0)**1.2/(rad(j)/10.d0)**1.2*(gbeta/3.d0)**0.8d0/(et(k)/1.d52)**0.2d0
        anarate_coef  = 720.d0*(Mdet/32.5d0)/(dist/10.d0)**2&
                      & *(mass(i)/1.4d0)**7.5/(rad(j)/10.d0)**8*(gbeta/3.)**5
        integral      = anarate_coef * 100.d0 * (t0/100.d0)**(-6.5d0)/6.5d0
        mlogLH(i,j,k) = integral
        do p = 1, datsize
          rate_guess    = anarate_coef/((tdata(p)+t0)/100.)**7.5
          e_ave_guess   = 25.3d0*(mass(i)/1.4d0)**1.5/(rad(j)/10.d0)**2*(gbeta/3.d0)/((tdata(p)+t0)/100.d0)**1.5
          temp_guess    = e_ave_guess*FDint4/FDint5
          mlogLH(i,j,k) = mlogLH(i,j,k) - log(rate_guess) - 4.d0*log(edata(p)) + log(1.d0+exp(edata(p)/temp_guess)) &
                        &           + log(FDint4) + 5.d0*log(temp_guess)
        enddo
      enddo
    enddo
  enddo

end subroutine eval_unbinned_likelihood

!!!!! full-binned analysis !!!!!
!!!!! binning time data !!!!!
subroutine binning(tdata,edata,dt_ini,tmin,tmax,t_bin,dt,e_bin,de,hist,status,datsize,tbinnumber)
  implicit none
  integer, parameter :: ebinnumber = 30
  real*8, parameter :: emax = 100.d0
  real*8, dimension(datsize), intent(in) :: tdata,edata
  real*8, intent(in) :: dt_ini,tmin,tmax
  real*8, dimension(tbinnumber), intent(out) :: t_bin,dt
  real*8, dimension(ebinnumber), intent(out) :: e_bin,de
  real*8, dimension(tbinnumber,ebinnumber), intent(out) :: hist
  integer, intent(out) :: status
  integer, intent(in) :: datsize,tbinnumber
  integer :: i,j,p,ite,itemax
  real*8 :: ratio,ratio_min,ratio_max,einterval

  itemax = 100
  ratio_min = 1.d0
  ratio_max = 2.d0
  dt(1) = dt_ini
  t_bin(1) = tmin+dt(1)
  status = 0

  do ite = 1, itemax
    ratio = 0.5d0*(ratio_min + ratio_max)
    do i=2,tbinnumber
      dt(i)  = dt(i-1)*ratio
      t_bin(i) = t_bin(i-1) + dt(i)
    enddo
    if(t_bin(tbinnumber) .gt. tmax) ratio_max = ratio
    if(t_bin(tbinnumber) .lt. tmax) ratio_min = ratio
    if(abs(t_bin(tbinnumber)-tmax)/tmax .lt. 1.d-4) exit
  enddo
  if(ite .ge. itemax) status = 1
  t_bin = t_bin - 0.5d0*dt

  einterval = emax/dble(ebinnumber)
  de = einterval
  e_bin(1) = de(1)
  do j = 2, ebinnumber
    e_bin(j) = e_bin(j-1)+de(j)
  enddo
  e_bin = e_bin - 0.5d0*de

  do p = 1, datsize
    do i = 1, tbinnumber
      if(tdata(p) .gt. t_bin(i)-0.5d0*dt(i) .and. tdata(p) .lt. t_bin(i)+0.5d0*dt(i)) then
        j = int(edata(p)/einterval) + 1
        hist(i,j) = hist(i,j) + 1.d0
      endif
    enddo
  enddo
  
end subroutine binning

!!!!! subroutine to give mlogLH = -log(likelihood) !!!!!
subroutine eval_binned_likelihood(t_bin,dt,e_bin,de,hist,mass,rad,et,gbeta,dist,Mdet,mlogLH,tbinnumber,nparam)
  implicit none
  integer, parameter :: ebinnumber = 30
  real*8, dimension(tbinnumber), intent(in) :: t_bin,dt
  real*8, dimension(ebinnumber), intent(in) :: e_bin,de
  real*8, dimension(tbinnumber,ebinnumber), intent(in) :: hist
  real*8, dimension(nparam), intent(in) :: mass,rad,et
  real*8, intent(in) :: gbeta,dist,Mdet
  real*8, dimension(nparam,nparam,nparam), intent(out) :: mlogLH
  integer, intent(in) :: tbinnumber,nparam
  integer :: i, j, k, l, m
  real*8 :: anarate_coef, rate_guess, e_ave_guess, t0, temp_guess, spectrum
  real*8 :: log_n_in_bin_factorial
  real*8, parameter :: FDint4=23.3d0,FDint5=118.d0,tiny=1.d-50

  mlogLH = 0.d0

  do i=1,nparam
    do j=1,nparam
      do k=1,nparam
        t0            = 210.d0*(mass(i)/1.4d0)**1.2/(rad(j)/10.d0)**1.2*(gbeta/3.d0)**0.8d0/(et(k)/1.d52)**0.2d0
        anarate_coef  = 720.d0*(Mdet/32.5d0)/(dist/10.d0)**2&
                      & *(mass(i)/1.4d0)**7.5/(rad(j)/10.d0)**8*(gbeta/3.)**5
        do l = 1, tbinnumber
          rate_guess    = anarate_coef/((t_bin(l)+t0)/100.)**7.5
          e_ave_guess   = 25.3d0*(mass(i)/1.4d0)**1.5/(rad(j)/10.d0)**2*(gbeta/3.d0)/((t_bin(l)+t0)/100.d0)**1.5
          temp_guess    = e_ave_guess*FDint4/FDint5
          do m = 1, ebinnumber
            spectrum = max(e_bin(m)**4.d0/FDint4/temp_guess**5.d0/(1.d0+exp(e_bin(m)/temp_guess)),tiny)
            if(hist(l,m) .le. 0.5d0) then
              log_n_in_bin_factorial = 0.d0
            else
              log_n_in_bin_factorial = hist(l,m)*log(hist(l,m)) - hist(l,m) + 0.5d0*log(2.d0*3.141592d0*dble(hist(l,m)))   !log(hist(l,m)!) by Stirling's formula
            endif
            mlogLH(i,j,k) = mlogLH(i,j,k) - hist(l,m)*log(rate_guess*spectrum*dt(l)*de(m)) + log_n_in_bin_factorial &
                          &               + rate_guess*spectrum*dt(l)*de(m)
          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine eval_binned_likelihood

!!!!! time-binned analysis !!!!!
!!!!! binning time data !!!!!
subroutine time_binning(tdata,dt_ini,tmin,tmax,t_bin,dt,thist,status,datsize,tbinnumber)
  implicit none
  real*8, dimension(datsize), intent(in) :: tdata
  real*8, intent(in) :: dt_ini,tmin,tmax
  real*8, dimension(tbinnumber), intent(out) :: t_bin,dt,thist
  integer, intent(out) :: status
  integer, intent(in) :: datsize,tbinnumber
  integer :: i,j,ite,itemax
  real*8 :: ratio,ratio_min,ratio_max

  itemax = 100
  ratio_min = 1.d0
  ratio_max = 2.d0
  dt(1) = dt_ini
  t_bin(1) = tmin+dt(1)
  status = 0

  do ite = 1, itemax
    ratio = 0.5d0*(ratio_min + ratio_max)
    do i=2,tbinnumber
      dt(i)  = dt(i-1)*ratio
      t_bin(i) = t_bin(i-1) + dt(i)
    enddo
    if(t_bin(tbinnumber) .gt. tmax) ratio_max = ratio
    if(t_bin(tbinnumber) .lt. tmax) ratio_min = ratio
    if(abs(t_bin(tbinnumber)-tmax)/tmax .lt. 1.d-4) exit
  enddo
  if(ite .ge. itemax) status = 1
  t_bin = t_bin - 0.5d0*dt

  do i = 1, datsize
    do j = 1, tbinnumber
      if(tdata(i) .gt. t_bin(j)-0.5d0*dt(j) .and. tdata(i) .lt. t_bin(j)+0.5d0*dt(j)) thist(j) = thist(j) + 1.d0
    enddo
  enddo

end subroutine time_binning

!!!!! subroutine to give mlogLH = -log(likelihood) !!!!!
subroutine eval_timebinned_likelihood(t_bin,dt,thist,edata,mass,rad,et,gbeta,dist,Mdet,mlogLH,tbinnumber,datsize,nparam)
  implicit none
  real*8, dimension(tbinnumber), intent(in) :: t_bin,dt,thist
  real*8, dimension(datsize), intent(in) :: edata
  real*8, dimension(nparam), intent(in) :: mass,rad,et
  real*8, intent(in) :: gbeta,dist,Mdet
  real*8, dimension(nparam,nparam,nparam), intent(out) :: mlogLH
  integer, intent(in) :: tbinnumber,datsize,nparam
  integer :: i, j, k, l, m, p
  integer :: int_thist,counted
  real*8 :: anarate_coef, rate_guess, e_ave_guess, t0, temp_guess
  real*8 :: log_n_in_bin_factorial
  real*8, parameter :: FDint4=23.3d0,FDint5=118.d0

  mlogLH = 0.d0

  do i=1,nparam
    do j=1,nparam
      do k=1,nparam
        t0            = 210.d0*(mass(i)/1.4d0)**1.2/(rad(j)/10.d0)**1.2*(gbeta/3.d0)**0.8d0/(et(k)/1.d52)**0.2d0
        anarate_coef  = 720.d0*(Mdet/32.5d0)/(dist/10.d0)**2&
                      & *(mass(i)/1.4d0)**7.5/(rad(j)/10.d0)**8*(gbeta/3.)**5
        counted = 0
        do l = 1, tbinnumber
          rate_guess    = anarate_coef/((t_bin(l)+t0)/100.)**7.5
          e_ave_guess   = 25.3d0*(mass(i)/1.4d0)**1.5/(rad(j)/10.d0)**2*(gbeta/3.d0)/((t_bin(l)+t0)/100.d0)**1.5
          temp_guess    = e_ave_guess*FDint4/FDint5
          log_n_in_bin_factorial = thist(l)*log(thist(l)) - thist(l) + 0.5d0*log(2.d0*3.141592d0*dble(thist(l)))   !log(thist(l)!) by Stirling's formula
          mlogLH(i,j,k) = mlogLH(i,j,k) - thist(l)*log(rate_guess*dt(l)) + log_n_in_bin_factorial + rate_guess*dt(l)
          int_thist = int(thist(l))
          do m = 1, int_thist
            p = counted + m
            mlogLH(i,j,k) = mlogLH(i,j,k) - 4.d0*log(edata(p)) + log(1.d0+exp(edata(p)/temp_guess)) &
                          &               + log(FDint4) + 5.d0*log(temp_guess)
          enddo
          counted = counted + int_thist
        enddo
      enddo
    enddo
  enddo

end subroutine eval_timebinned_likelihood

!!!!! Gaussian-approximation analysis !!!!!
subroutine event_energy_averaging(thist,edata,e_ave,tbinnumber,datsize)
  implicit none
  real*8, dimension(tbinnumber), intent(in) :: thist
  real*8, dimension(datsize), intent(in) :: edata
  real*8, dimension(tbinnumber), intent(out) :: e_ave
  integer, intent(in) :: tbinnumber,datsize
  integer :: l,m,p
  integer :: counted,int_thist
  
  counted = 0
  e_ave = 0.d0
  do l = 1, tbinnumber
    int_thist = int(thist(l))
    if(int_thist .eq. 0) then
      write(*,*)l,'-th time bin has no events'
      write(*,*)'try another time bin number'
      stop
    endif
    do m = 1, int_thist
      p = counted + m
      e_ave(l) = e_ave(l) + edata(p)
    enddo
    counted = counted + int_thist
    e_ave(l) = e_ave(l)/thist(l)
  enddo
  
  return
end subroutine event_energy_averaging

!!!!! subroutine to give mlogLH = -log(likelihood) !!!!!
subroutine eval_Gaussian_likelihood(t_bin,dt,thist,e_ave,mass,rad,et,gbeta,dist,Mdet,mlogLH,tbinnumber,nparam)
  implicit none
  real*8, dimension(tbinnumber), intent(in) :: t_bin,dt,thist,e_ave
  real*8, dimension(nparam), intent(in) :: mass,rad,et
  real*8, intent(in) :: gbeta,dist,Mdet
  real*8, dimension(nparam,nparam,nparam), intent(out) :: mlogLH
  integer, intent(in) :: tbinnumber,nparam
  integer :: i, j, k, l, m, p
  integer :: int_thist,counted
  real*8 :: anarate_coef, rate_guess, e_ave_guess, t0, temp_guess
  real*8 :: log_n_in_bin_factorial

  mlogLH = 0.d0

  do i=1,nparam
    do j=1,nparam
      do k=1,nparam
        t0            = 210.d0*(mass(i)/1.4d0)**1.2/(rad(j)/10.d0)**1.2*(gbeta/3.d0)**0.8d0/(et(k)/1.d52)**0.2d0
        anarate_coef  = 720.d0*(Mdet/32.5d0)/(dist/10.d0)**2&
                      & *(mass(i)/1.4d0)**7.5/(rad(j)/10.d0)**8*(gbeta/3.)**5
        do l = 1, tbinnumber
          rate_guess    = anarate_coef/((t_bin(l)+t0)/100.)**7.5
          e_ave_guess   = 25.3d0*(mass(i)/1.4d0)**1.5/(rad(j)/10.d0)**2*(gbeta/3.d0)/((t_bin(l)+t0)/100.d0)**1.5
          mlogLH(i,j,k) = mlogLH(i,j,k) + 0.5d0*(thist(l)-rate_guess*dt(l))**2/rate_guess/dt(l) &
                        &               + 0.5d0*(e_ave(l)-e_ave_guess)**2/(e_ave_guess*0.05d0)**2
!          mlogLH(i,j,k) = mlogLH(i,j,k) + 0.5d0*(thist(l)/dt(l)-rate_guess)**2/(rate_guess**2/thist(l)) &
!                        &               + 0.5d0*(e_ave(l)-e_ave_guess)**2/(e_ave_guess*0.05d0)**2
        enddo
      enddo
    enddo
  enddo

end subroutine eval_Gaussian_likelihood


!!!!! subroutines to marginalize 3D likelihood to 2D likelihood !!!!!
subroutine mr_marginalize(mass,rad,et,mlogLH,LH_MR,CI95level,CI68level,peak,nparam)
  implicit none
  real*8, dimension(nparam), intent(in) :: mass,rad,et
  real*8, dimension(nparam,nparam,nparam), intent(in) :: mlogLH
  real*8, dimension(nparam,nparam), intent(out) :: LH_MR
  real*8, intent(out) :: CI95level,CI68level,peak(2)
  integer, intent(in) :: nparam

  integer :: i_min, j_min, k_min, minimums(3)
  integer :: i,j,k,l
  real*8 :: deltamlogLH, sumLH
  real*8 :: LHnormal,LHmax,ddparam
  real*8, dimension(100) :: levelintegral
  integer :: search_channel
  integer :: peakloc(2)

  minimums = minloc(mlogLH)
  i_min = minimums(1)
  j_min = minimums(2)
  k_min = minimums(3)

  do i=1,nparam
    do j=1,nparam
      deltamlogLH = 0.d0
      sumLH = 0.d0
      do k=1,nparam
        deltamlogLH     = mlogLH(i,j,k) - mlogLH(i_min,j_min,k_min)
        sumLH  = sumLH + exp(-deltamlogLH)
      enddo
      if(sumLH .lt. 1.d-98) sumLH = 0.d0
      LH_MR(i,j) = sumLH
    enddo
  enddo
  LHmax    = maxval(LH_MR)
  LHnormal = sum(LH_MR)
  peakloc  = maxloc(LH_MR)
  peak(1)  = mass(peakloc(1))
  peak(2)  = rad(peakloc(2))
  
  levelintegral = 0.d0
  do l=1,100
    do i=1,nparam
      do j=1,nparam
        if(LH_MR(i,j)/LHmax .ge. dble(l)*1.d-2) then
          levelintegral(l) = levelintegral(l) + LH_MR(i,j)/LHnormal
        endif
      enddo
    enddo
  enddo
  search_channel = 1 ! 1:CI95level, 2:CI68level
  do l=1,100
    if(levelintegral(l) .lt. 0.95d0 .and. search_channel .eq. 1) then
      CI95level = (dble(l)*(0.95d0-levelintegral(l-1))+dble(l-1)*(levelintegral(l)-0.95d0))&
                & /(levelintegral(l)-levelintegral(l-1))*1.d-2*LHmax/LHnormal
      search_channel = search_channel + 1
    endif
    if(levelintegral(l) .lt. 0.68d0 .and. search_channel .eq. 2) then
      CI68level = (dble(l)*(0.68d0-levelintegral(l-1))+dble(l-1)*(levelintegral(l)-0.68d0))&
                & /(levelintegral(l)-levelintegral(l-1))*1.d-2*LHmax/LHnormal
      exit
    endif
  enddo
  
  ddparam = (mass(nparam)-mass(1))*(rad(nparam)-rad(1))/(nparam-1)**2
  LH_MR = LH_MR/LHnormal/ddparam

end subroutine mr_marginalize

subroutine re_marginalize(mass,rad,et,mlogLH,LH_RE,CI95level,CI68level,peak,nparam)
  implicit none
  real*8, dimension(nparam), intent(in) :: mass,rad,et
  real*8, dimension(nparam,nparam,nparam), intent(in) :: mlogLH
  real*8, dimension(nparam,nparam), intent(out) :: LH_RE
  real*8, intent(out) :: CI95level,CI68level,peak(2)
  integer, intent(in) :: nparam

  integer :: i_min, j_min, k_min, minimums(3)
  integer :: i,j,k,l
  real*8 :: deltamlogLH, sumLH
  real*8 :: LHnormal,LHmax,ddparam
  real*8, dimension(100) :: levelintegral
  integer :: search_channel
  integer :: peakloc(2)

  minimums = minloc(mlogLH)
  i_min = minimums(1)
  j_min = minimums(2)
  k_min = minimums(3)

  do j=1,nparam
    do k=1,nparam
      deltamlogLH = 0.d0
      sumLH = 0.d0
      do i=1,nparam
        deltamlogLH     = mlogLH(i,j,k) - mlogLH(i_min,j_min,k_min)
        sumLH  = sumLH + exp(-deltamlogLH)
      enddo
      if(sumLH .lt. 1.d-98) sumLH = 0.d0
      LH_RE(j,k) = sumLH
    enddo
  enddo
  LHmax = maxval(LH_RE)
  LHnormal = sum(LH_RE)
  peakloc  = maxloc(LH_RE)
  peak(1)  = rad(peakloc(1))
  peak(2)  = et(peakloc(2))
  
  levelintegral = 0.d0
  do l=1,100
    do j=1,nparam
      do k=1,nparam
        if(LH_RE(j,k)/LHmax .ge. dble(l)*1.d-2) then
          levelintegral(l) = levelintegral(l) + LH_RE(j,k)/LHnormal
        endif
      enddo
    enddo
  enddo
  search_channel = 1 ! 1:CI95level, 2:CI68level
  do l=1,100
    if(levelintegral(l) .lt. 0.95d0 .and. search_channel .eq. 1) then
      CI95level = (dble(l)*(0.95d0-levelintegral(l-1))+dble(l-1)*(levelintegral(l)-0.95d0))&
                & /(levelintegral(l)-levelintegral(l-1))*1.d-2*LHmax/LHnormal
      search_channel = search_channel + 1
    endif
    if(levelintegral(l) .lt. 0.68d0 .and. search_channel .eq. 2) then
      CI68level = (dble(l)*(0.68d0-levelintegral(l-1))+dble(l-1)*(levelintegral(l)-0.68d0))&
                & /(levelintegral(l)-levelintegral(l-1))*1.d-2*LHmax/LHnormal
      exit
    endif
  enddo
  
  ddparam = (rad(nparam)-rad(1))*(et(nparam)-et(1))/(nparam-1)**2
  LH_RE = LH_RE/LHnormal/ddparam

end subroutine re_marginalize

subroutine em_marginalize(mass,rad,et,mlogLH,LH_EM,CI95level,CI68level,peak,nparam)
  implicit none
  real*8, dimension(nparam), intent(in) :: mass,rad,et
  real*8, dimension(nparam,nparam,nparam), intent(in) :: mlogLH
  real*8, dimension(nparam,nparam), intent(out) :: LH_EM
  real*8, intent(out) :: CI95level,CI68level,peak(2)
  integer, intent(in) :: nparam
  

  integer :: i_min, j_min, k_min, minimums(3)
  integer :: i,j,k,l
  real*8 :: deltamlogLH, sumLH
  real*8 :: LHnormal,LHmax,ddparam
  real*8, dimension(100) :: levelintegral
  integer :: search_channel
  integer :: peakloc(2)

  minimums = minloc(mlogLH)
  i_min = minimums(1)
  j_min = minimums(2)
  k_min = minimums(3)

  do k=1,nparam
    do i=1,nparam
      deltamlogLH = 0.d0
      sumLH = 0.d0
      do j=1,nparam
        deltamlogLH     = mlogLH(i,j,k) - mlogLH(i_min,j_min,k_min)
        sumLH  = sumLH + exp(-deltamlogLH)
      enddo
      if(sumLH .lt. 1.d-98) sumLH = 0.d0
      LH_EM(k,i) = sumLH
    enddo
  enddo
  LHmax = maxval(LH_EM)
  LHnormal = sum(LH_EM)
  peakloc  = maxloc(LH_EM)
  peak(1)  = et(peakloc(1))
  peak(2)  = mass(peakloc(2))
  
  levelintegral = 0.d0
  do l=1,100
    do k=1,nparam
      do i=1,nparam
        if(LH_EM(k,i)/LHmax .ge. dble(l)*1.d-2) then
          levelintegral(l) = levelintegral(l) + LH_EM(k,i)/LHnormal
        endif
      enddo
    enddo
  enddo
  search_channel = 1 ! 1:CI95level, 2:CI68level
  do l=1,100
    if(levelintegral(l) .lt. 0.95d0 .and. search_channel .eq. 1) then
      CI95level = (dble(l)*(0.95d0-levelintegral(l-1))+dble(l-1)*(levelintegral(l)-0.95d0))&
                & /(levelintegral(l)-levelintegral(l-1))*1.d-2*LHmax/LHnormal
      search_channel = search_channel + 1
    endif
    if(levelintegral(l) .lt. 0.68d0 .and. search_channel .eq. 2) then
      CI68level = (dble(l)*(0.68d0-levelintegral(l-1))+dble(l-1)*(levelintegral(l)-0.68d0))&
                & /(levelintegral(l)-levelintegral(l-1))*1.d-2*LHmax/LHnormal
      exit
    endif
  enddo
  
  ddparam = (et(nparam)-et(1))*(mass(nparam)-mass(1))/(nparam-1)**2
  LH_EM = LH_EM/LHnormal/ddparam

end subroutine em_marginalize

!!!!! subroutines to marginalize 3D likelihood to 1D likelihood !!!!!
subroutine m_marginalize(mass,rad,et,mlogLH,LH_M,CIandBF,nparam)
  implicit none
  real*8, dimension(nparam), intent(in) :: mass,rad,et
  real*8, dimension(nparam,nparam,nparam), intent(in) :: mlogLH
  real*8, dimension(nparam), intent(out) :: LH_M
  real*8, intent(out) :: CIandBF(5)
  integer, intent(in) :: nparam

  integer :: i_min, j_min, k_min, minimums(3)
  integer :: i,j,k,l
  real*8 :: deltamlogLH, sumLH
  real*8 :: LHnormal,LHmax,dparam
  real*8, dimension(100) :: levelintegral
  integer :: paramax(1)
  integer :: search_channel !1: lower95, 2:lower68, 3:upper68, 4:upper95
  real*8 :: CI95level,CI68level,maxpara,lower68,upper68,lower95,upper95,peak
  real*8 :: bestfit,x1,x2,x3,y1,y2,y3

  minimums = minloc(mlogLH)
  i_min = minimums(1)
  j_min = minimums(2)
  k_min = minimums(3)

  do i=1,nparam
    deltamlogLH = 0.d0
    sumLH = 0.d0
    do j=1,nparam
      do k=1,nparam
        deltamlogLH     = mlogLH(i,j,k) - mlogLH(i_min,j_min,k_min)
        sumLH  = sumLH + exp(-deltamlogLH)
      enddo
      if(sumLH .lt. 1.d-98) sumLH = 0.d0
      LH_M(i) = sumLH
    enddo
  enddo
  LHmax = maxval(LH_M)
  LHnormal = sum(LH_M)
  levelintegral = 0.d0
  do l=1,100
    do i=1,nparam
      if(LH_M(i)/LHmax .ge. dble(l)*1.d-2) then
        levelintegral(l) = levelintegral(l) + LH_M(i)/LHnormal
      endif
    enddo
  enddo
  
  search_channel = 1 ! 1:CI95level, 2:CI68level
  do l=1,100
    if(levelintegral(l) .lt. 0.95d0 .and. search_channel .eq. 1) then
      CI95level = (dble(l)*(0.95d0-levelintegral(l-1))+dble(l-1)*(levelintegral(l)-0.95d0))&
                &/(levelintegral(l)-levelintegral(l-1))*1.d-2*LHmax/LHnormal
      search_channel = search_channel + 1
    endif
    if(levelintegral(l) .lt. 0.68d0 .and. search_channel .eq. 2) then
      CI68level = (dble(l)*(0.68d0-levelintegral(l-1))+dble(l-1)*(levelintegral(l)-0.68d0))&
                &/(levelintegral(l)-levelintegral(l-1))*1.d-2*LHmax/LHnormal
      exit
    endif
  enddo
  
  LH_M  = LH_M/LHnormal

  paramax = maxloc(LH_M)
  search_channel = 1 !1:lower95, 2:lower68, 3:upper68, 4:upper95
  do i=1,nparam
    if(LH_M(i) .ge. CI95level .and. search_channel .eq. 1) then
      lower95 = (mass(i)*(CI95level-LH_M(i-1))+mass(i-1)*(LH_M(i)-CI95level))/(LH_M(i) - LH_M(i-1))
      search_channel = search_channel + 1
    endif
    if(LH_M(i) .ge. CI68level .and. search_channel .eq. 2) then
      lower68 = (mass(i)*(CI68level-LH_M(i-1))+mass(i-1)*(LH_M(i)-CI68level))/(LH_M(i) - LH_M(i-1))
      search_channel = search_channel + 1
    endif
    if(LH_M(i) .le. CI68level .and. search_channel .eq. 3) then
      upper68 = (mass(i)*(CI68level-LH_M(i-1))+mass(i-1)*(LH_M(i)-CI68level))/(LH_M(i) - LH_M(i-1))
      search_channel = search_channel + 1
    endif
    if(LH_M(i) .le. CI95level .and. search_channel .eq. 4) then
      upper95 = (mass(i)*(CI95level-LH_M(i-1))+mass(i-1)*(LH_M(i)-CI95level))/(LH_M(i) - LH_M(i-1))
      exit
    endif
  enddo
  x1 = mass(paramax(1)-1)
  y1 = LH_M(paramax(1)-1)
  x2 = mass(paramax(1))
  y2 = LH_M(paramax(1))
  x3 = mass(paramax(1)+1)
  y3 = LH_M(paramax(1)+1)
  bestfit = 0.5d0*(y1*(x2**2-x3**2)-y2*(x1**2-x3**2)+y3*(x1**2-x2**2))/(y1*(x2-x3)-y2*(x1-x3)+y3*(x1-x2))
  
  CIandBF(1) = lower95
  CIandBF(2) = lower68
  CIandBF(3) = bestfit
  CIandBF(4) = upper68
  CIandBF(5) = upper95
  
  dparam = (mass(nparam)-mass(1))/(nparam-1)
  LH_M = LH_M/dparam
  
end subroutine m_marginalize

subroutine r_marginalize(mass,rad,et,mlogLH,LH_R,CIandBF,nparam)
  implicit none
  real*8, dimension(nparam), intent(in) :: mass,rad,et
  real*8, dimension(nparam,nparam,nparam), intent(in) :: mlogLH
  real*8, dimension(nparam), intent(out) :: LH_R
  real*8, intent(out) :: CIandBF(5)
  integer, intent(in) :: nparam

  integer :: i_min, j_min, k_min, minimums(3)
  integer :: i,j,k,l
  real*8 :: deltamlogLH, sumLH
  real*8 :: LHnormal,LHmax,dparam
  real*8, dimension(100) :: levelintegral
  integer :: paramax(1)
  integer :: search_channel !1: lower95, 2:lower68, 3:upper68, 4:upper95
  real*8 :: CI95level,CI68level,maxpara,lower68,upper68,lower95,upper95,peak
  real*8 :: bestfit,x1,x2,x3,y1,y2,y3

  minimums = minloc(mlogLH)
  i_min = minimums(1)
  j_min = minimums(2)
  k_min = minimums(3)

  do j=1,nparam
    deltamlogLH = 0.d0
    sumLH = 0.d0
    do k=1,nparam
      do i=1,nparam
        deltamlogLH     = mlogLH(i,j,k) - mlogLH(i_min,j_min,k_min)
        sumLH  = sumLH + exp(-deltamlogLH)
      enddo
      if(sumLH .lt. 1.d-98) sumLH = 0.d0
      LH_R(j) = sumLH
    enddo
  enddo
  LHmax = maxval(LH_R)
  LHnormal = sum(LH_R)
  levelintegral = 0.d0
  do l=1,100
    do i=1,nparam
      if(LH_R(i)/LHmax .ge. dble(l)*1.d-2) then
        levelintegral(l) = levelintegral(l) + LH_R(i)/LHnormal
      endif
    enddo
  enddo
  
  search_channel = 1 ! 1:CI95level, 2:CI68level
  do l=1,100
    if(levelintegral(l) .lt. 0.95d0 .and. search_channel .eq. 1) then
      CI95level = (dble(l)*(0.95d0-levelintegral(l-1))+dble(l-1)*(levelintegral(l)-0.95d0))&
                &/(levelintegral(l)-levelintegral(l-1))*1.d-2*LHmax/LHnormal
      search_channel = search_channel + 1
    endif
    if(levelintegral(l) .lt. 0.68d0 .and. search_channel .eq. 2) then
      CI68level = (dble(l)*(0.68d0-levelintegral(l-1))+dble(l-1)*(levelintegral(l)-0.68d0))&
                &/(levelintegral(l)-levelintegral(l-1))*1.d-2*LHmax/LHnormal
      exit
    endif
  enddo
 
  LH_R  = LH_R/LHnormal
  
  paramax = maxloc(LH_R)
  search_channel = 1 !1:lower95, 2:lower68, 3:upper68, 4:upper95
  do i=1,nparam
    if(LH_R(i) .ge. CI95level .and. search_channel .eq. 1) then
      lower95 = (rad(i)*(CI95level-LH_R(i-1))+rad(i-1)*(LH_R(i)-CI95level))/(LH_R(i) - LH_R(i-1))
      search_channel = search_channel + 1
    endif
    if(LH_R(i) .ge. CI68level .and. search_channel .eq. 2) then
      lower68 = (rad(i)*(CI68level-LH_R(i-1))+rad(i-1)*(LH_R(i)-CI68level))/(LH_R(i) - LH_R(i-1))
      search_channel = search_channel + 1
    endif
    if(LH_R(i) .le. CI68level .and. search_channel .eq. 3) then
      upper68 = (rad(i)*(CI68level-LH_R(i-1))+rad(i-1)*(LH_R(i)-CI68level))/(LH_R(i) - LH_R(i-1))
      search_channel = search_channel + 1
    endif
    if(LH_R(i) .le. CI95level .and. search_channel .eq. 4) then
      upper95 = (rad(i)*(CI95level-LH_R(i-1))+rad(i-1)*(LH_R(i)-CI95level))/(LH_R(i) - LH_R(i-1))
      exit
    endif
  enddo
  x1 = rad(paramax(1)-1)
  y1 = LH_R(paramax(1)-1)
  x2 = rad(paramax(1))
  y2 = LH_R(paramax(1))
  x3 = rad(paramax(1)+1)
  y3 = LH_R(paramax(1)+1)
  bestfit = 0.5d0*(y1*(x2**2-x3**2)-y2*(x1**2-x3**2)+y3*(x1**2-x2**2))/(y1*(x2-x3)-y2*(x1-x3)+y3*(x1-x2))
  
  CIandBF(1) = lower95
  CIandBF(2) = lower68
  CIandBF(3) = bestfit
  CIandBF(4) = upper68
  CIandBF(5) = upper95
  
  dparam = (rad(nparam)-rad(1))/(nparam-1)
  LH_R = LH_R/dparam

end subroutine r_marginalize

subroutine e_marginalize(mass,rad,et,mlogLH,LH_E,CIandBF,nparam)
  implicit none
  real*8, dimension(nparam), intent(in) :: mass,rad,et
  real*8, dimension(nparam,nparam,nparam), intent(in) :: mlogLH
  real*8, dimension(nparam), intent(out) :: LH_E
  real*8, intent(out) :: CIandBF(5)
  integer, intent(in) :: nparam

  integer :: i_min, j_min, k_min, minimums(3)
  integer :: i,j,k,l
  real*8 :: deltamlogLH, sumLH
  real*8 :: LHnormal,LHmax,dparam
  real*8, dimension(100) :: levelintegral
  integer :: paramax(1)
  integer :: search_channel !1: lower95, 2:lower68, 3:upper68, 4:upper95
  real*8 :: CI95level,CI68level,maxpara,lower68,upper68,lower95,upper95,peak
  real*8 :: bestfit,x1,x2,x3,y1,y2,y3

  minimums = minloc(mlogLH)
  i_min = minimums(1)
  j_min = minimums(2)
  k_min = minimums(3)

  do k=1,nparam
    deltamlogLH = 0.d0
    sumLH = 0.d0
    do i=1,nparam
      do j=1,nparam
        deltamlogLH     = mlogLH(i,j,k) - mlogLH(i_min,j_min,k_min)
        sumLH  = sumLH + exp(-deltamlogLH)
      enddo
      if(sumLH .lt. 1.d-98) sumLH = 0.d0
      LH_E(k) = sumLH
    enddo
  enddo
  LHmax = maxval(LH_E)
  LHnormal = sum(LH_E)
  levelintegral = 0.d0
  do l=1,100
    do i=1,nparam
      if(LH_E(i)/LHmax .ge. dble(l)*1.d-2) then
        levelintegral(l) = levelintegral(l) + LH_E(i)/LHnormal
      endif
    enddo
  enddo
  
  search_channel = 1 ! 1:CI95level, 2:CI68level
  do l=1,100
    if(levelintegral(l) .lt. 0.95d0 .and. search_channel .eq. 1) then
      CI95level = (dble(l)*(0.95d0-levelintegral(l-1))+dble(l-1)*(levelintegral(l)-0.95d0))&
                &/(levelintegral(l)-levelintegral(l-1))*1.d-2*LHmax/LHnormal
      search_channel = search_channel + 1
    endif
    if(levelintegral(l) .lt. 0.68d0 .and. search_channel .eq. 2) then
      CI68level = (dble(l)*(0.68d0-levelintegral(l-1))+dble(l-1)*(levelintegral(l)-0.68d0))&
                &/(levelintegral(l)-levelintegral(l-1))*1.d-2*LHmax/LHnormal
      exit
    endif
  enddo
  
  LH_E  = LH_E/LHnormal
  
  paramax = maxloc(LH_E)
  search_channel = 1 !1:lower95, 2:lower68, 3:upper68, 4:upper95
  do i=1,nparam
    if(LH_E(i) .ge. CI95level .and. search_channel .eq. 1) then
      lower95 = (et(i)*(CI95level-LH_E(i-1))+et(i-1)*(LH_E(i)-CI95level))/(LH_E(i) - LH_E(i-1))
      search_channel = search_channel + 1
    endif
    if(LH_E(i) .ge. CI68level .and. search_channel .eq. 2) then
      lower68 = (et(i)*(CI68level-LH_E(i-1))+et(i-1)*(LH_E(i)-CI68level))/(LH_E(i) - LH_E(i-1))
      search_channel = search_channel + 1
    endif
    if(LH_E(i) .le. CI68level .and. search_channel .eq. 3) then
      upper68 = (et(i)*(CI68level-LH_E(i-1))+et(i-1)*(LH_E(i)-CI68level))/(LH_E(i) - LH_E(i-1))
      search_channel = search_channel + 1
    endif
    if(LH_E(i) .le. CI95level .and. search_channel .eq. 4) then
      upper95 = (et(i)*(CI95level-LH_E(i-1))+et(i-1)*(LH_E(i)-CI95level))/(LH_E(i) - LH_E(i-1))
      exit
    endif
  enddo
  x1 = et(paramax(1)-1)
  y1 = LH_E(paramax(1)-1)
  x2 = et(paramax(1))
  y2 = LH_E(paramax(1))
  x3 = et(paramax(1)+1)
  y3 = LH_E(paramax(1)+1)
  bestfit = 0.5d0*(y1*(x2**2-x3**2)-y2*(x1**2-x3**2)+y3*(x1**2-x2**2))/(y1*(x2-x3)-y2*(x1-x3)+y3*(x1-x2))
  
  CIandBF(1) = lower95
  CIandBF(2) = lower68
  CIandBF(3) = bestfit
  CIandBF(4) = upper68
  CIandBF(5) = upper95
  
  dparam = (et(nparam)-et(1))/(nparam-1)
  LH_E = LH_E/dparam

end subroutine e_marginalize

