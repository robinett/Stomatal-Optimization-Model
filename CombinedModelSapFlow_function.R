CombinedModel <- function(params){
  
  ###DISCUSS IMPLEMENTATION OF ZWEIFEL PARAMETERS
  ###Zweifel Parameters crown storage
  PCrownMax = params$PCrownMax #maximum water content of crown, g, Zweiful Table 2
  PsiCrownMin = params$Pcrit #minimum water potential of storage pool crown, MPa, Zweiful Table 2
  RXylem = params$RXylem #Zweifel Table 2
  lwidth = params$lwidth
  leafarea = params$leafarea
  PsiStoragevect = params$PsiStoragevect
  
  
  simsteps = 100
  
  testtimes = rep(0,simsteps)
  for(i in 1:simsteps){
    testtimes[i] = i
  }
  
  
  
  ### Environmental Drivers
  Ca = params$Ca #atmospheric CO2, in Pa, table 1 Sperry
  VPD = params$VPD #atmospheric water pressure deficit, in kPa, Table 1 Sperry
  Oa = 21000 #atmospheric O2 conc, in Pa, Table 1 Sperry
  Patm = params$Patm #atmospheric air pressure, in kPa, Table 1 Sperry
  Ps = params$Ps #soil water potential, in MPa, Table 1 Sperry
  Q = params$Q #PAR photon flux density, umol*s-1*m-2, Table 1 Sperry
  TA = params$TA #air temperature, degree C, Table 1 Sperry
  u = params$u #wind speed, ms-1, Table 1 Sperry
  
  ### Hydraulic cost and photosynthetic gain parameters
  c = 0.9 #curvature of light responce curve, Table 1 Sperry
  cprime = 0.98 #curvature factor for Je versus Jc limited photosynthesis, Table 1 Sperry
  d = lwidth*0.72 #characteristic dimension of the leaf width, Campbell & Norman
  Vmax25 = params$Vmax25 #maximum carboxylation rate at 25C, umol*s-1*m-2, Sperry Table 1
  Jmax25 = 1.67*Vmax25 #maximium electron transport rate at 25C, umol*s-1*m-2, Sperry Table 1
  Kc = 41 #at 25C, michaelis-menton constant for carboxylation, Pa, Table 1 Sperry
  Ko = 28202 #at 25C, michaelis-menton constant for oxygenation, Pa, Table 1 Sperry
  Rabs = 740 #absorbed long and short wave radiation, W*m-2, Specified
  fancyA = 0.3 #Quantum yield of electron transport, mol*mol-1, Table 1 Sperry
  gammastar = 4.36 #at 25C, CO2 compensation point, Pa, Table 1 Sperry
  epsilon = 0.97 #emissivity, Table 1 Sperry
  rho = 1 #air density, kg*m-3, specified constant
  sigma = 5.67E-8 #Stefan-Boltzman constant, W*m-2*K-4, specified constant
  Cp = 1004 #specific heat capacity of dry air, J*K-1*kg-1, specified constant
  Ta = 273.15 + TA #temperature in degrees kelvin, specified constant
  lambda = 2.5E6 #latent heat of vaporization in J/kg, specified constant
  kmax = params$kmax



  Evect = rep(0,simsteps)
  Pcvect = rep(0,simsteps)
  negPcvect = rep(0,simsteps)
  smallFcrownvect = rep(0,simsteps)
  storageflowgramsvect = rep(0,simsteps)
  storageflowmmolvect = rep(0,simsteps)
  Enostoragevect = rep(0,simsteps)
  
  
  
  for(i in 1:simsteps){
    Pc = (i/simsteps)*Pcrit
    Pcvect[i] = Pc
    negPcvect[i] = -Pc
    E = integrate(E_canopy_transp_rate_Aspinwall,Pc,Ps,kmax,a,P50)
    storageFlowgrams = ((Pc-PsiStoragevect)/RXylem) #if it is positive, flow INTO storage
    #convert from g/h to mmol*s-1*m-2
    storageFlowmmol = storageFlowgrams/18.02 #now in mol/h
    storageFlowmmol = storageFlowmmol*1000 #now in mmol/h
    storageFlowmmol = storageFlowmmol/3600 #now in mmol/s
    storageFlowmmol = storageFlowmmol/leafarea #now in mmol*s-1*m-2
    
    StorageE = E$value - storageFlowmmol #units converted to mmol*s-1*m-2
    Enostoragevect[i] = E$value
    if(StorageE<0){
      StorageE = NA
    }
    Evect[i] = StorageE
    storageflowgramsvect[i] = storageFlowgrams
    storageflowmmolvect[i] = storageFlowmmol
  }
  #plot(negPcvect,Evect)
  
  
  #####
  #calculate A curve
  
  
  #initialize vectors
  Avect = rep(0,simsteps)
  Tlvect = rep(0,simsteps)
  Dlvect = rep(0,simsteps)
  Gwvect = rep(0,simsteps)
  Gcvect = rep(0,simsteps)
  Jvect = rep(0,simsteps)
  Civect = rep(0,simsteps)
  Avect = rep(0,simsteps)
  
  
  #calculate Leaf Temp for Pc=Ps through Pc=Pcrit (THE ENERGY)
  for (i in 1:simsteps){
    if(is.na(Evect[i]) == FALSE){
      #convert E from mmol*s-1*m-2 to kg*s-1*m-2
      E = Evect[i]/1000 #E here is in mol*s-1*m-2
      E = E*18.02 #E here is in g*s-1*m-2
      E = E/1000 #E is now in kg*s-1*m-2
      interval = c(250,350)
      
      Tl = uniroot.all(solveTlForFunction,interval,lower=250,upper=350,rho=rho,Cp=Cp,u=u,Ta=Ta,lambda=lambda,E=E,epsilon=epsilon,sigma=sigma,Rabs=Rabs)
      Tleaf = Tl - 273.15
      Tlvect[i] = Tleaf
    }
    if(is.na(Evect[i]) == TRUE || Tlvect[i] < TA){
      Tlvect[i] = NA
    }
    
  }
  #plot(negPcvect,Tlvect)
  
  #calculate leaf to air vapor pres. deficit for Pc=Ps through Pc=Pcrit
  for (i in 1:simsteps){
    if(is.na(Tlvect[i]) == FALSE){
      Dlvect[i] = VPD
    }
    
    if(is.na(Tlvect[i]) == TRUE){
      Dlvect[i] = NA
    }
    
    # if(Dlvect[i]<.1 || is.na(Dlvect[i]) == TRUE){
    #   Dlvect[i] = NA
    # }
  }
  #plot(negPcvect,Dlvect)
  
  
  
  #calculate diffusive water conductance of leaf for Pc=Ps through Pc=Pcrit
  for (i in 1:simsteps){
    if(is.na(Dlvect[i]) == FALSE){
      Gwvect[i] = Gw_water_diff_cond_leaf(Evect[i],Dlvect[i])
    }
    
    if(is.na(Dlvect[i]) == TRUE){
      Gwvect[i] = NA
    }
  }
  #plot(negPcvect,Gwvect)
  
  #calculate diffusive CO2 conductance for leaf for Pc=Ps through Pc=Pcrit
  for (i in 1:simsteps){
    if(is.na(Gwvect[i]) == FALSE){
      Gcvect[i] = Gc_CO2_diff_cond_leaf(Gwvect[i])
    }
    
    if(is.na(Gwvect[i]) == TRUE){
      Gcvect[i] = NA
    }
  }
  #plot(negPcvect,Gcvect)
  
  #calculate rate of electron transport for Pc=Ps through Pc=Pcrit
  for (i in 1:simsteps){
    if(is.na(Gcvect[i]) == FALSE){
      Jvect[i] = J_elect_transp_rate(Q,Jmax25,c)
    }
    
    if(is.na(Gcvect[i]) == TRUE){
      Jvect[i] = NA
    }
  }
  #plot(negPcvect,Jvect)
  
  #calculate internal leaf CO2 conc. for Pc=Ps through Pc=Pcrit
  for (i in 1:simsteps){
    if(is.na(Jvect[i]) == FALSE){
      interval = c(2,200)
      Jint = Jvect[i]
      Gcint = Gcvect[i]
      Ci = optimize(solveCiForOptimize,interval,Jint=Jint,gammastar=gammastar,Vmax25=Vmax25,Kc=Kc,Oa=Oa,Ko=Ko,cprime=cprime,Gcint=Gcint,Ca=Ca,Patm=Patm)
      Cisolved = Ci$minimum
      Civect[i] = Cisolved
    }
    
    if(is.na(Jvect[i]) == TRUE){
      Civect[i] = NA
    }
    if(Civect[i] > (Ca-0.01) && Civect[i] < (Ca+0.01) && is.na(Civect[i]) == FALSE){
      Civect[i] = Ca
    }
  }
  #plot(negPcvect,Civect)
  
  #calculate gross assimilation rate for Pc=Ps through Pc=Pcrit
  for (i in 1:simsteps){
    if(is.na(Civect[i]) == FALSE){
      Avect[i] = A_gross_assim_rate(Civect[i],Gcvect[i],Ca,Patm)
    }
    
    if(is.na(Jvect[i]) == TRUE){
      Avect[i] = NA
    }
  }
  #plot(negPcvect,Avect)
  
  
  
  #####
  #calculate beta, theta, profit curves
  
  #initialize vectors
  betavect = rep(0,(simsteps-1))
  thetavect = rep(0,(simsteps-1))
  profitvect = rep(0,(simsteps-1))
  
  
  #calcualte beta (carbon gain function)
  Amax = max(Avect, na.rm = TRUE)
  
  for (i in 1:(simsteps-1)){
    if(Amax == 0 || is.na(Avect[i] == TRUE)){
      betavect[i] = 0
    }else{
      betavect[i] = Avect[i]/Amax
    }
  }
  #plot(negPcvect[1:(simsteps-1)],betavect)
  
  #calculate theta (cost function)
  
  dEvect = rep(0,(simsteps-1))
  dPcvect = rep(0,(simsteps-1))
  kcvect = rep(0,(simsteps-1))
  thetavect = rep(0,(simsteps-1))
  
  
  ###NA STEPS MAKING THE CODE MUCH SLOWER???
  O = 1
  for (i in 1:(simsteps-1)){
    dEvect[i] = Evect[i+1] - Evect[i]
    dPcvect[i] = negPcvect[i+1] - negPcvect[i]
    kcvect[i] = dEvect[i]/dPcvect[i]
    if(kcvect[i] == 0 || is.na(kcvect[i] == TRUE)){
      kcvect[i] = NA
    }
  }
  
  
  kcmax = max(kcvect, na.rm = TRUE)
  kcrit = kcvect[simsteps-1]
  
  for(i in 1:(simsteps-1)){
    thetavect[i] = ((kcmax - kcvect[i])/(kcmax-kcrit))
  }
  #thetavect[(O+1)] = NA
  for(i in 1:(simsteps-1)){
    if(is.na(thetavect[i]) == TRUE){
      thetavect[i] = 0
    }
  }
  #plot(negPcvect[1:(simsteps-1)],thetavect)
  
  #Find profit maximization point
  
  profitmaxvect = rep(0,simsteps)
  
  for(i in 1:(simsteps-1)){
    profitvect[i] = betavect[i]-thetavect[i]
  }
  #plot(negPcvect[1:(simsteps-1)],profitvect)
  
  profitmax = max(profitvect)
  profitmax
  
  P = 1
  
  while(profitvect[P] < profitmax || is.na(Avect[P]) == TRUE){
    P = P+1
  }
  
  storageflowgramschangevect = storageflowgramsvect[P]
  storageflowmmolchangevect = storageflowmmolvect[P]
  Echangevect = Evect[P]
  Tlchangevect = Tlvect[P]
  Dlchangevect = Dlvect[P]
  Gwchangevect = Gwvect[P]
  Gcchangevect = Gcvect[P]
  Jchangevect = Jvect[P]
  Cichangevect = Civect[P]
  Achangevect = Avect[P]
  thetachangevect = thetavect[P]
  betachangevect = betavect[P]
  profitchangevect = profitmax
  negPcchangevect = negPcvect[P]
  Pcchangevect = Pcvect[P]
  PercEchangevect = -storageflowmmolvect[P]/Enostoragevect[P]
  kcchangevect = kcvect[P]
  tempchangevect = Tlvect[P]
  Qchangevect = Q
  Enostoragechangevect = Enostoragevect[P]
  
  returnvals = data.frame("storageflowgramschangevect" = storageflowgramschangevect, 
                          "storageflowmmolchangevect" = storageflowmmolchangevect, 
                          "Echangevect" = Echangevect, 
                          "Tlchangevect" = Tlchangevect,
                          "Dlchangevect" = Dlchangevect,
                          "Gwchangevect" = Gwchangevect,
                          "Gcchangevect" = Gcchangevect,
                          "Cichangevect" = Cichangevect,
                          "Achangevect" = Achangevect,
                          "negPcchangevect" = negPcchangevect,
                          "PercEchangevect" = PercEchangevect)
  return(returnvals)
  
  
}