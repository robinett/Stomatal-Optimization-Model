SperryModel <- function(Ca = 40,Ps = 0, Q = 1000, TA = 25, humidity=0.5,VC = c(2,3),Vmax25=100,Pcrit=-3.5,Tleaf = 30, VPD = 2, kmax = 1){
  
  #####
  #set model functioning data
  
  nsteps = 300
  
  
  ### Environmental Drivers
  #Ca = 40 #atmospheric CO2, in Pa, table 1 Sperry
  #D = 1 #atmospheric water pressure deficit, in kPa, Table 1 Sperry
  Oa = 21000 #atmospheric O2 conc, in Pa, Table 1 Sperry
  Patm = 101.3 #atmospheric air pressure, in kPa, Table 1 Sperry
  #Ps = 0 #soil water potential, in MPa, Table 1 Sperry
  #Q = 2000 #PAR photon flux density, umol*s-1*m-2, Table 1 Sperry
  #TA = 25 #air temperature, degree C, Table 1 Sperry
  u = 2 #wind speed, ms-1, Table 1 Sperry
  #humidity = 0.5 #humidity, %/100, Specified in this model
  
  ### Hydraulic cost and photosynthetic gain parameters
  c = 0.9 #curvature of light responce curve, Table 1 Sperry
  cprime = 0.98 #curvature factor for Je versus Jc limited photosynthesis, Table 1 Sperry
  lwidth = 0.05 #leaf width, m, specified
  lwidthm = lwidth*100 #leaf width, m, specified
  d = lwidth*0.72 #characteristic dimension of the leaf width, Campbell & Norman
  #Vmax25 = 100 #maximum carboxylation rate at 25C, umol*s-1*m-2, Sperry Table 1
  Jmax25 = 1.67*Vmax25 #maximium electron transport rate at 25C, umol*s-1*m-2, Sperry Table 1
  Kc = 41 #at 25C, michaelis-menton constant for carboxylation, Pa, Table 1 Sperry
  Ko = 28202 #at 25C, michaelis-menton constant for oxygenation, Pa, Table 1 Sperry
  Rabs = 740 #absorbed long and short wave radiation, W*m-2, Specified
  #VC = c(2,3) #two parameter [b,c] Weibull vulnerability curve, Table 1 Sperry
  fancyA = 0.3 #Quantum yield of electron transport, mol*mol-1, Table 1 Sperry
  gammastar = 4.36 #at 25C, CO2 compensation point, Pa, Table 1 Sperry
  epsilon = 0.97 #emissivity, Table 1 Sperry
  rho = 1 #air density, kg*m-3, specified constant
  sigma = 5.67E-8 #Stefan-Boltzman constant, W*m-2*K-4, specified constant
  Cp = 1004 #specific heat capacity of dry air, J*K-1*kg-1, specified constant
  Ta = 273.15 + TA #temperature in degrees kelvin, specified constant
  lambda = 2.5E6 #latent heat of vaporization in J/kg, specified constant
  
  
  
  #####
  #Calculate E Curve
  
  Evect = rep(0,nsteps)
  Pcvect = rep(0,nsteps)
  negPcvect = rep(0,nsteps)
  smallFcrownvect = rep(0,nsteps)
  

  for (i in 1:nsteps){
    Pc = (i/nsteps)*(Pcrit-Ps) + Ps
    Pcvect[i] = Pc
    negPcvect[i] = -Pc
    E = integrate(E_canopy_transp_rate,lower=Pc,upper=Ps,kmax=kmax,b=VC[1],c=VC[2])
    E = E$value
    Evect[i] = E #in mmol*m-2*s-1
  }
  
  #plot(negPcvect,Evect)
  
  
  
  #####
  #calculate A curve
  
  
  #initialize vectors
  Avect = rep(0,nsteps)
  Tlvect = rep(0,nsteps)
  Dlvect = rep(0,nsteps)
  Gwvect = rep(0,nsteps)
  Gcvect = rep(0,nsteps)
  Jvect = rep(0,nsteps)
  Civect = rep(0,nsteps)
  Avect = rep(0,nsteps)
  
  
  #calculate Leaf Temp for Pc=Ps through Pc=Pcrit (THE ENERGY)
  for (i in 1:nsteps){
    Tlvect[i] = Tleaf # in C
  }
  #plot(negPcvect,Tlvect)

  
  #calculate leaf to air vapor pres. deficit for Pc=Ps through Pc=Pcrit
  for (i in 1:nsteps){
    #substitute in once function is working
    Dlvect[i] = VPD #in kPa
    #because issues
  }
  #plot(negPcvect,Dlvect)
  
  
  
  #calculate diffusive water conductance of leaf for Pc=Ps through Pc=Pcrit
  for (i in 1:nsteps){
    Gwvect[i] = Gw_water_diff_cond_leaf(Evect[i],Dlvect[i])
  }
  #plot(negPcvect,Gwvect)
  
  
  #calculate diffusive CO2 conductance for leaf for Pc=Ps through Pc=Pcrit
  for (i in 1:nsteps){
    Gcvect[i] = Gc_CO2_diff_cond_leaf(Gwvect[i])
  }
  #plot(negPcvect,Gcvect)
  
  
  #calculate rate of electron transport for Pc=Ps through Pc=Pcrit
  for (i in 1:nsteps){
    Jvect[i] = J_elect_transp_rate(Q,Jmax25,c)
  }
  #plot(negPcvect,Jvect)
  
  
  #calculate internal leaf CO2 conc. for Pc=Ps through Pc=Pcrit
  for (i in 1:nsteps){
    interval = c(3,70)
    Jint = Jvect[i]
    Gcint = Gcvect[i]
    Ci = uniroot(solveCiForFunction,Jint=Jint,gammastar=gammastar,Vmax25=Vmax25,Kc=Kc,Oa=Oa,Ko=Ko,cprime=cprime,Gcint=Gcint,Ca=Ca,Patm=Patm,interval)
    Cisolved = Ci$root #in Pa
    Civect[i] = Cisolved
  }
  #plot(negPcvect,Civect)
  
  
  #calculate gross assimilation rate for Pc=Ps through Pc=Pcrit
  for (i in 1:nsteps){
    Avect[i] = A_gross_assim_rate(Civect[i],Gcvect[i],Ca,Patm)
  }
  #plot(negPcvect,Avect)
  
  
  
  
  #####
  #calculate beta, theta, profit curves
  
  #initialize vectors
  betavect = rep(0,(nsteps-1))
  thetavect = rep(0,(nsteps-1))
  profitvect = rep(0,(nsteps-1))
  
  
  #calcualte beta (carbon gain function)
  Amax = max(Avect, na.rm = TRUE)
  
  for (i in 1:(nsteps-1)){
    if(Amax == 0 || is.na(Avect[i] == TRUE)){
      betavect[i] = 0
    }else{
      betavect[i] = Avect[i]/Amax
    }
  }
  #plot(negPcvect[1:(nsteps-1)],betavect)
  
  #calculate theta (cost function)
  
  dEvect = rep(0,(nsteps-1))
  dPcvect = rep(0,(nsteps-1))
  kcvect = rep(0,(nsteps-1))
  thetavect = rep(0,(nsteps-1))
  
  
  O = 1
  for (i in 1:(nsteps-1)){
    dEvect[i] = Evect[i+1] - Evect[i]
    dPcvect[i] = negPcvect[i+1] - negPcvect[i]
    kcvect[i] = dEvect[i]/dPcvect[i]
    if(kcvect[i] == 0 || is.na(kcvect[i] == TRUE)){
      kcvect[i] = NA
    }
  }
  
  
  kcmax = max(kcvect, na.rm = TRUE)
  kcrit = kcvect[nsteps-1]
  
  for(i in 1:(nsteps-1)){
    thetavect[i] = ((kcmax - kcvect[i])/(kcmax-kcrit))
  }
  #thetavect[(O+1)] = NA
  for(i in 1:(nsteps-1)){
    if(is.na(thetavect[i]) == TRUE){
      thetavect[i] = 0
    }
  }
  #plot(negPcvect[1:(nsteps-1)],thetavect)
  
  #Find profit maximization point
  
  profitmaxvect = rep(0,nsteps)
  
  for(i in 1:(nsteps-1)){
    profitvect[i] = betavect[i]-thetavect[i]
  }
  #plot(negPcvect[1:(nsteps-1)],profitvect)
  
  profitmax = max(profitvect)
  profitmax
  
  P = 1
  
  while(profitvect[P] < profitmax || is.na(Avect[P]) == TRUE){
    P = P+1
  }
  
  returnvals = c(Pcvect[P],Evect[P],Gwvect[P],Gcvect[P],Civect[P],Avect[P])
  return(returnvals)
}