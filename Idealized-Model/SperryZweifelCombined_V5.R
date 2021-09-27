rm(list = ls())
if(!is.null(dev.list())) dev.off()
setwd("/Users/Trent Robinett/Documents/Notre Dame/Medvigy Lab/Sperry_etal_2017_code")
library(FAdist)
library(ggplot2)
library("rootSolve")
library("pracma")
source("/Users/Trent Robinett/Documents/Notre Dame/Medvigy Lab/Sperry_etal_2017_code/Sperryfuncs.R")

#####
#NOTES
#if storageFlow is positive, flow INTO storage


#####
#set Parameters for model

#Interpolated from figures
kmax = 4.3 #derived from graph in Sperry via Vmax25
kcmax = 4.3
Pcrit = -3.2 #derived from Figure 1a

### Environmental Drivers
Ca = 40 #atmospheric CO2, in Pa, table 1 Sperry
D = 1 #atmospheric water pressure deficit, in kPa, Table 1 Sperry
Oa = 21000 #atmospheric O2 conc, in Pa, Table 1 Sperry
Patm = 101.3 #atmospheric air pressure, in kPa, Table 1 Sperry
Ps = 0 #soil water potential, in MPa, Table 1 Sperry
Qday = 2000 #PAR photon flux density, umol*s-1*m-2, Table 1 Sperry
Qnight = 0
TA = 25 #air temperature, degree C, Table 1 Sperry
u = 2 #wind speed, ms-1, Table 1 Sperry
humidity = 0.5 #humidity, %/100, Specified in this model

### Hydraulic cost and photosynthetic gain parameters
c = 0.9 #curvature of light responce curve, Table 1 Sperry
cprime = 0.98 #curvature factor for Je versus Jc limited photosynthesis, Table 1 Sperry
lwidth = 0.05 #leaf width, m, specified
lwidthm = lwidth*100 #leaf width, m, specified
leafarea = (2*lwidth)*(lwidth) #approximating leaf area
d = lwidth*0.72 #characteristic dimension of the leaf width, Campbell & Norman
Vmax25 = 100 #maximum carboxylation rate at 25C, umol*s-1*m-2, Sperry Table 1
Jmax25 = 1.67*Vmax25 #maximium electron transport rate at 25C, umol*s-1*m-2, Sperry Table 1
Kc = 41 #at 25C, michaelis-menton constant for carboxylation, Pa, Table 1 Sperry
Ko = 28202 #at 25C, michaelis-menton constant for oxygenation, Pa, Table 1 Sperry
Rabs = 740 #absorbed long and short wave radiation, W*m-2, Specified
VC = c(2,3) #two parameter [b,c] Weibull vulnerability curve, Table 1 Sperry
fancyA = 0.3 #Quantum yield of electron transport, mol*mol-1, Table 1 Sperry
gammastar = 4.36 #at 25C, CO2 compensation point, Pa, Table 1 Sperry
epsilon = 0.97 #emissivity, Table 1 Sperry
rho = 1 #air density, kg*m-3, specified constant
sigma = 5.67E-8 #Stefan-Boltzman constant, W*m-2*K-4, specified constant
Cp = 1004 #specific heat capacity of dry air, J*K-1*kg-1, specified constant
Ta = 273.15 + TA #temperature in degrees kelvin, specified constant
lambda = 2.5E6 #latent heat of vaporization in J/kg, specified constant
Enorm = 2.74 #storage with at optimal Pc with no storage modification to Sperry


###Zweifel Parameters stem storage
PBarkMax = 9.8 #maximum water content of bark, g, Zweiful Table 2
PsiBarkMin = -3 #minimum water potential of storage pool bark, MPa, Zweiful Table 2
k1Bark = 8.0 #amount of stored water at point of inflection of Psi(Pbark), unitless, Zweiful Table 2
k2Bark = 1.1 #index for slope of Psi(Pbark) at poitn of inflection, unitless, Zweiful Table 2

###Zweifel Parameters crown storage
PCrownMax = 87.3 #maximum water content of crown, g, Zweiful Table 2
PsiCrownMin = -4 #minimum water potential of storage pool crown, MPa, Zweiful Table 2
k1Crown = 69.3 #amount of stored water at point of inflection of Psi(Pcrown), unitless, Zweiful Table 2
k2Crown = 4.3 #index for slope of Psi(Pcrown) at point of inflection, unitless, Zweiful Table 2

###Zweifel other parameters
PsiBase = 0 #Zweifel calibration discussion
RXylem = 1 #Zweifel Table 2
maxtranspiration = 2 #Zweiful graph of T(input)
transpirationconst = 0 #for testing the model


#intialize time conditions
nhours = 100 #specify the number of hours to simluate
nstepsperhour = 6 #determine the number of steps per hour
nsteps = nstepsperhour*nhours
times = linspace(0, nhours, nsteps)


###Define inputs for ODE Function
PBarkInitial = PBarkMax
PCrownInitial = PCrownMax






storageflowgramschangevect = rep(0,nsteps)
storageflowmmolchangevect = rep(0,nsteps)
Echangevect = rep(0,nsteps)
Tlchangevect = rep(0,nsteps)
Dlchangevect = rep(0,nsteps)
Gwchangevect = rep(0,nsteps)
Gcchangevect = rep(0,nsteps)
Jchangevect = rep(0,nsteps)
Cichangevect = rep(0,nsteps)
Achangevect = rep(0,nsteps)
thetachangevect = rep(0,nsteps)
betachangevect = rep(0,nsteps)
profitchangevect = rep(0,nsteps)
zvect = rep(0,nsteps)
Pvect = rep(0,nsteps)
Pcchangevect = rep(0,nsteps)
PercEchangevect = rep(0,nsteps)
PcRangevect = rep(0,nsteps)
Pcchangevect = rep(0,nsteps)
kcchangevect = rep(0,nsteps)
tempchangevect = rep(0,nsteps)
PStoragevect = rep(0,nsteps)
PsiStoragevect = rep(0,nsteps)
negPcchangevect = rep(0,nsteps)
Qchangevect = rep(0,nsteps)
timeschangevect = rep(0,nsteps)
Enostoragechangevect = rep(0,nsteps)

###Initialize vectors
PStoragevect[1] = PCrownMax
Pcchangevect[1] = -0.65
PsiStoragevect[1] = -0.6
R = 0.5 






for(z in 1:nsteps){
  
  time = times[z]
  timeschangevect[z] = time
  print(time)
  zvect[z] = z
  Evect = rep(0,nsteps)
  Pcvect = rep(0,nsteps)
  negPcvect = rep(0,nsteps)
  smallFcrownvect = rep(0,nsteps)
  storageflowgramsvect = rep(0,nsteps)
  storageflowmmolvect = rep(0,nsteps)
  Enostoragevect = rep(0,nsteps)
  
  
  if(((time%%24))>=(0) && ((time%%24))<(7)){
    Q = Qnight
  }
  if(((time%%24))>=(7) && ((time%%24))<(13)){
    Q = ((time%%24)-7)/6*Qday
  }
  if(((time%%24))>=(13) && ((time%%24))<(19)){
    Q = (19-(time%%24))/6*Qday
  }
  if(((time%%24))>=(19) && ((time%%24))<(24)){
    Q = Qnight
  }
  
  for (i in 1:nsteps){
    Pc = (i/nsteps)*Pcrit
    Pcvect[i] = Pc
    negPcvect[i] = -Pc
    E = integrate(E_canopy_transp_rate,Pc,Ps,kmax,VC[1],VC[2])
    storageFlowgrams = ((Pc-PsiStoragevect[z])/R) #if it is positive, flow INTO storage
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
    if(is.na(Evect[i]) == FALSE){
      #convert E from mmol*s-1*m-2 to kg*s-1*m-2
      E = Evect[i]/1000 #E here is in mol*s-1*m-2
      E = E*18.02 #E here is in g*s-1*m-2
      E = E/1000 #E is now in kg*s-1*m-2
      
      Tl = uniroot.all(solveTl,lower= 250 ,upper=350)
      Tleaf = Tl - 273.15
      Tlvect[i] = Tleaf
    }
    if(is.na(Evect[i]) == TRUE || Tlvect[i] < TA){
      Tlvect[i] = NA
    }
    
  }
  #plot(negPcvect,Tlvect)
  
  #calculate leaf to air vapor pres. deficit for Pc=Ps through Pc=Pcrit
  for (i in 1:nsteps){
    if(is.na(Tlvect[i]) == FALSE){
      Dlvect[i] = Dl_leaf_air_pres_def(Tlvect[i],TA,humidity)
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
  for (i in 1:nsteps){
    if(is.na(Dlvect[i]) == FALSE){
      Gwvect[i] = Gw_water_diff_cond_leaf(Evect[i],Dlvect[i])
    }
    
    if(is.na(Dlvect[i]) == TRUE){
      Gwvect[i] = NA
    }
  }
  #plot(negPcvect,Gwvect)
  
  #calculate diffusive CO2 conductance for leaf for Pc=Ps through Pc=Pcrit
  for (i in 1:nsteps){
    if(is.na(Gwvect[i]) == FALSE){
      Gcvect[i] = Gc_CO2_diff_cond_leaf(Gwvect[i])
    }
    
    if(is.na(Gwvect[i]) == TRUE){
      Gcvect[i] = NA
    }
  }
  #plot(negPcvect,Gcvect)
  
  #calculate rate of electron transport for Pc=Ps through Pc=Pcrit
  for (i in 1:nsteps){
    if(is.na(Gcvect[i]) == FALSE){
      Jvect[i] = J_elect_transp_rate(Q,Jmax25,c)
    }
    
    if(is.na(Gcvect[i]) == TRUE){
      Jvect[i] = NA
    }
  }
  #plot(negPcvect,Jvect)
  
  #calculate internal leaf CO2 conc. for Pc=Ps through Pc=Pcrit
  for (i in 1:nsteps){
    if(is.na(Jvect[i]) == FALSE){
      Ci = uniroot(solveCi,lower=0.1,upper=60)
      Cisolved = Ci$root
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
  for (i in 1:nsteps){
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
  
  
  ###NA STEPS MAKING THE CODE MUCH SLOWER???
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
  
  if(negPcvect[P] < 0.1){
    PcRangevect[z] = "0.1<Pc"
  }
  if(negPcvect[P] < 0.2 && negPcvect[P] > 0.1){
    PcRangevect[z] = "0.1<Pc<0.2"
  }
  if(negPcvect[P] < 0.3 && negPcvect[P] > 0.2){
    PcRangevect[z] = "0.2<Pc<0.3"
  }
  if(negPcvect[P] < 0.4 && negPcvect[P] > 0.3){
    PcRangevect[z] = "0.3<Pc<0.4"
  }
  if(negPcvect[P] < 0.5 && negPcvect[P] > 0.4){
    PcRangevect[z] = "0.5<Pc<0.4"
  }
  if(negPcvect[P] < 0.6 && negPcvect[P] > 0.5){
    PcRangevect[z] = "0.5<Pc<0.6"
  }
  if(negPcvect[P] < 0.7 && negPcvect[P] > 0.6){
    PcRangevect[z] = "0.6<Pc<0.7"
  }
  if(negPcvect[P] < 0.8 && negPcvect[P] > 0.7){
    PcRangevect[z] = "0.7<Pc<0.8"
  }
  if(negPcvect[P] < 0.9 && negPcvect[P] > 0.8){
    PcRangevect[z] = "0.8<Pc<0.9"
  }
  if(negPcvect[P] < 1.0 && negPcvect[P] > 0.9){
    PcRangevect[z] = "0.9<Pc<1.0"
  }
  if(negPcvect[P] < 2.0 && negPcvect[P] > 1.0){
    PcRangevect[z] = "1.0<Pc<2.0"
  }
  if(negPcvect[P] < 3.0 && negPcvect[P] > 2.0){
    PcRangevect[z] = "2.0<Pc<3.0"
  }
  if(negPcvect[P] > 3.0){
    PcRangevect[z] = "Pc>3.0"
  }
  
  
  storageflowgramschangevect[z] = storageflowgramsvect[P]
  storageflowmmolchangevect[z] = storageflowmmolvect[P]
  Echangevect[z] = Evect[P]
  Tlchangevect[z] = Tlvect[P]
  Dlchangevect[z] = Dlvect[P]
  Gwchangevect[z] = Gwvect[P]
  Gcchangevect[z] = Gcvect[P]
  Jchangevect[z] = Jvect[P]
  Cichangevect[z] = Civect[P]
  Achangevect[z] = Avect[P]
  thetachangevect[z] = thetavect[P]
  betachangevect[z] = betavect[P]
  profitchangevect[z] = profitmax
  negPcchangevect[z] = negPcvect[P]
  Pcchangevect[z] = Pcvect[i]
  PercEchangevect[z] = -storageflowmmolvect[P]/Enostoragevect[P]
  kcchangevect[z] = kcvect[P]
  tempchangevect[z] = Tlvect[P]
  Qchangevect[z] = Q
  Enostoragechangevect[z] = Enostoragevect[P]
  
  if(z<nsteps){
    PStoragevect[z+1] = PStoragevect[z] + storageflowgramsvect[P]/nstepsperhour
    
    PsiStoragevect[z+1] = ((PsiCrownMin)/(exp((-k1Crown+PStoragevect[z+1])/(k2Crown))+1))
  }
  
}


plot(timeschangevect,Echangevect)
abline(v=7, col="green")
abline(v=31, col="green")
abline(v=55, col="green")
abline(v=79, col="green")
abline(v=19, col="red")
abline(v=43, col="red")
abline(v=67, col="red")
abline(v=91, col="red")
#plot(timeschangevect,PStoragevect)
abline(v=7, col="green")
abline(v=31, col="green")
abline(v=55, col="green")
abline(v=79, col="green")
abline(v=19, col="red")
abline(v=43, col="red")
abline(v=67, col="red")
abline(v=91, col="red")
plot(timeschangevect,PsiStoragevect)
abline(v=7, col="green")
abline(v=31, col="green")
abline(v=55, col="green")
abline(v=79, col="green")
abline(v=19, col="red")
abline(v=43, col="red")
abline(v=67, col="red")
abline(v=91, col="red")
plot(timeschangevect,storageflowgramschangevect)
abline(v=7, col="green")
abline(v=31, col="green")
abline(v=55, col="green")
abline(v=79, col="green")
abline(v=19, col="red")
abline(v=43, col="red")
abline(v=67, col="red")
abline(v=91, col="red")
abline(h=0)
plot(timeschangevect,Qchangevect)
abline(v=7, col="green")
abline(v=31, col="green")
abline(v=55, col="green")
abline(v=19, col="red")
abline(v=43, col="red")
abline(v=67, col="red")
plot(timeschangevect,PercEchangevect)
abline(v=7, col="green")
abline(v=31, col="green")
abline(v=55, col="green")
abline(v=79, col="green")
abline(v=19, col="red")
abline(v=43, col="red")
abline(v=67, col="red")
abline(v=91, col="red")
abline(h=0.7)
abline(h=0)
plot(timeschangevect,negPcchangevect)
abline(v=7, col="green")
abline(v=31, col="green")
abline(v=55, col="green")
abline(v=79, col="green")
abline(v=19, col="red")
abline(v=43, col="red")
abline(v=67, col="red")
abline(v=91, col="red")
plot(timeschangevect,Achangevect)
plot(timeschangevect,Enostoragechangevect)
abline(v=7, col="green")
abline(v=31, col="green")
abline(v=55, col="green")
abline(v=79, col="green")
abline(v=19, col="red")
abline(v=43, col="red")
abline(v=67, col="red")
abline(v=91, col="red")

datacomparison = data.frame("PercentSapflow"=PercEchangevect,
                            "HoursElapsed"=timeschangevect)

jpeg("PercentFromStorage.jpg", width=15)
ggplot(datacomparison, aes(x=HoursElapsed, y=PercentSapflow)) +
  geom_point(size=3) +
  geom_vline(xintercept = 7, color = "green", size = .7) +
  geom_vline(xintercept = 31, color = "green", size = .7) +
  geom_vline(xintercept = 55, color = "green", size = .7) +
  geom_vline(xintercept = 79, color = "green", size = .7) +
  geom_vline(xintercept = 19, color = "red", size = .7) +
  geom_vline(xintercept = 43, color = "red", size = .7) +
  geom_vline(xintercept = 67, color = "red", size = .7) +
  geom_vline(xintercept = 91, color = "red", size = .7) +
  xlab("Hours Elapsed") + ylab("% of Sapflow from Storage") +
  scale_y_continuous(breaks=seq(-1.1,1.5,0.4)) +
  geom_hline(yintercept = 0.7, color = "blue", size = .7) +
  theme(
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14)
  )
dev.off()


