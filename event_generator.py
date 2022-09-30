# event generation code based on Suwa et al. PTEP, 013E01, 2021
# below, \mathcal{R} in comments means equation (55) of the paper

import numpy as np

def event_generator():
  # constants
  FDint4 = 23.3308744907
  FDint5 = 118.266130956
  emin = 0.01
  emax = 100.0
  enum = 1000
  
  # load parameters
  pnsparams = np.loadtxt('./sourceparams.dat')
  obsparams = np.loadtxt('./parameters.dat')
  
  Mpns  = pnsparams[0]
  Rpns  = pnsparams[1]
  Etot  = pnsparams[2]
  gbeta = pnsparams[3]
  
  dist = obsparams[1]
  Mdet = obsparams[2]
  
  # determine important parameters for analytic light curve
  t0              = 210.0 * (Mpns/1.4)**1.2 / (Rpns/10.0)**1.2 * (gbeta/3.0)**0.8 / (Etot/1.0e52)**0.2
  anarate_coef    = 720.0 * (Mdet/32.5)/(dist/10.0)**2 * (Mpns/1.4)**7.5 / (Rpns/10.0)**8.0 * (gbeta/3.0)**5.0
  total_exp_count = anarate_coef * 100.0 * (t0/100.0)**(-6.5)/6.5 # = \int_0^\infty \mathcal{R}(t) dt
  
  s = np.random.poisson(total_exp_count) # determine total event number
  tevents = np.sort(t0 * ( (1.0-np.random.rand(s))**(-2.0/13.0) - 1.0)) # get 1D array of event time t from normalized accumulated count rate at detector = \int_0^t \mathcal{R}(t') dt'/\int_0^\infty \mathcal{R}(t') dt' = ((t + t0)/t0)^(-13/2)
  temp = 25.3 * (Mpns/1.4)**1.5 / (Rpns/10.0)**2.0 * (gbeta/3.0) / ((tevents+t0)/100.0)**1.5 * FDint4/FDint5 # get 1D array of PNS temperature at each time
  ene = np.linspace(emin,emax,enum) # prepare energy grid
  de = ene[1]-ene[0] # prepare energy interval for energy integral
  FDdist = np.divide.outer(np.power(ene,4),np.power(temp,5))/(1.0+np.exp(np.divide.outer(ene,temp)))*de/FDint4 # make 2D (time x energy) array of Fermi--Dirac distribution. de/FDint4 is added for the next step integration
  FDcum = np.add.accumulate(FDdist) # get normalized accumulated energy distribution by integrating w.r.t. energy: FDcum(t,e) = \int_0^e FDdist(t,e') de'/\int_0^\infty FDdist(t,e') de'
  
  erand = np.random.rand(s) # dice roll to determine event energy
  eevents=de*(0.5+np.abs(FDcum-erand).argmin(axis=0)) # abs(FDcum-erand) gives the 'distance' of the value of FDcum from the rolled dice, and hence argmin() returns the nearest (minimum distance) energy grid ID. other factors convert ID into energy
  
  events = np.stack([tevents,eevents])
  np.savetxt('time_energy.dat',events.T)

event_generator()
