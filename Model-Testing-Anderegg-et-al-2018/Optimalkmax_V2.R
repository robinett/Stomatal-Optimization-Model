OptimalkmaxV2 = function(VC,Vmax25,Pcmin,Ci=28,D=1,Ps=0,TA=25){

  
  nsteps = 1000
  
  ### Environmental Drivers
  Ca = Ci/0.7 #atmospheric CO2, in Pa, table 1 Sperry
  #D = 1 #atmospheric water pressure deficit, in kPa, Table 1 Sperry
  Oa = 21000 #atmospheric O2 conc, in Pa, Table 1 Sperry
  Patm = 101.3 #atmospheric air pressure, in kPa, Table 1 Sperry
  #Ps = 0 #soil water potential, in MPa, Table 1 Sperry
  Q = 2000 #PAR photon flux density, umol*s-1*m-2, Table 1 Sperry
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
  
  alpha = 0.3 #quantum yeild of electron transport
  J = (alpha*Q + Jmax25 - ((alpha*Q + Jmax25)^(2) - 4*c*alpha*Q*Jmax25)^(0.5))/(2*c)
  Jc = (Vmax25*(Ci - gammastar))/(Ci + Kc*(1+(Oa/Ko)))
  Je = (J/4)*((Ci - gammastar)/(Ci + 2*gammastar))
  A = (Je + Jc - ((Je + Jc)^(2) - 4*cprime*Je*Jc)^(0.5))/(2*cprime)
  
  Gc = (A*Patm*1000)/(Ca-Ci) #umol*s-1*m-2
  
  Gw = Gc*1.6/1000 #mmol*s-1*m-2
  
  E = Gw*D/101.3 #mmol*s-1*m-2
  
  kmaxguessvect = rep(0,nsteps)
  Eguessvect = rep(0,nsteps)
  
  for(i in 1:nsteps){
    kmaxguess = 10*i/nsteps
    kmaxguessvect[i] = kmaxguess
    b = VC[1]
    c = VC[2]
    Eguess = integrate(E_canopy_transp_rate,Pcmin,Ps,kmaxguess,b,c)
    Eguessvect[i] = Eguess$value
  }
  
  z=1
  while(Eguessvect[z] < E){
    z = z+1
  }
  
  kmaxtrue = kmaxguessvect[z]
  
  return(kmaxtrue)
}