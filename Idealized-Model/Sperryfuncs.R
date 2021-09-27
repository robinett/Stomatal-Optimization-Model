E_canopy_transp_rate <- function(P,kmax,b,c){
  library(FAdist)
  unIntE = kmax*exp(-((-P/b)^c)) 
  return(unIntE)
}


E_canopy_transp_rate_Aspinwall <- function(P,kmax,a,P50){
  library(FAdist)
  unIntE = kmax*(exp(a*(P-P50)) / (1+exp(a*(P-P50)))) 
  return(unIntE)
}


Tl_temp_leaf <- function(TA,E,Rstar,u,d,epsilon,rho){
  sigma = 5.67E-8 #Stefan-Boltzman constant, W*m-2*K-4
  Cp = 1004 #specific heat capacity of dry air, J*K-1*kg-1
  Ta = 273.15 + TA #temperature in degrees kelvin
  lambda =  2.5E6 #latent heat of vaporization in J/kg
  gb = 0.189*(u/d)^(-0.5) #in mol*m-2*s-1
  top = Rabs - ((epsilon*sigma)*(Ta^4)) - lambda*(E/1000)/2
  bottom = Cp*(gr + gha)
  secondhalf = top/bottom
  Tl = TA + secondhalf
  return(Tl)
}

#E needs to be in molwater*m^-2*s-1


Dl_leaf_air_pres_def <- function(Tl,TA,D){
  satvapleaf = (6.1078*exp((17.269*Tl)/(237.3+Tl)))/10 #these are calculated in KPa
  #typical assumption is that air inside leaf is saturated
  satvapair = (6.1078*exp((17.269*TA)/(237.3+TA)))/10
  actvapair = satvapair - D
  #multiply air by the relitive humidity [0,1], mark this as an input parameter
  #can significantly chand Dl
  deficit = satvapleaf - actvapair
  return(deficit)
}


Gw_water_diff_cond_leaf <- function(E,Dl){
  Gw = (E)/Dl*101.3 #mmol*s-1*m-2
}



Gc_CO2_diff_cond_leaf <- function(Gw){
  Gc = Gw*1000/1.6 #umol*s-1*m-2
}


J_elect_transp_rate <- function(Q,Jmax,c){
  alpha = 0.3 #quantum yeild of electron transport
  J = (alpha*Q + Jmax -((alpha*Q + Jmax)^(2) - 4*c*alpha*Q*Jmax)^0.5)/(2*c)
}

# 0 = ((((((Jvect[i]/4)*((Ci-gammastar)/(Ci+(2*gammastar))))) + (((Vmax*(Ci-gammastar))/(Ci + Kc*(1+(Oa/Ko))))) - (((((((Jvect[i]/4)*((Ci-gammastar)/(Ci+(2*gammastar))))) + (((Vmax*(Ci-gammastar))/(Ci + Kc*(1+(Oa/Ko))))))^(2)) - (4*cprime*(((Jvect[i]/4)*((Ci-gammastar)/(Ci+(2*gammastar)))))*(((Vmax*(Ci-gammastar))/(Ci + Kc*(1+(Oa/Ko)))))))^(0.5)))/(2*cprime)) - ((Gc)*(Ca-Ci))/(Patm))
# A = (((Je) + (Jc) - (((((Je) + (Jc))^(2)) - (4*cprime*(Je)*(Jc)))^(0.5)))/(2*cprime))
# J = Jvect[i]
# Je = ((J/4)*((Ci-gammastar)/(Ci+(2*gammastar))))
# Jc = ((Vmax*(Ci-gammastar))/(Ci + Kc*(1+(Oa/Ko))))
# old: ((((Jvect[i]/4)*((Ci-gammastar)/(Ci+2*gammastar)))+(Vmax25*(Ci-gammastar)/(Ci + Kc*(1+(Oa/Ko))))-((((Jvect[i]/4)*((Ci-gammastar)/(Ci+2*gammastar)))+(Vmax25*(Ci-gammastar)/(Ci + Kc*(1+(Oa/Ko)))))^2 - 4*cprime*((Jvect[i]/4)*((Ci-gammastar)/(Ci+2*gammastar)))*(Vmax25*(Ci-gammastar)/(Ci + Kc*(1+(Oa/Ko)))))^0.5)/2*cprime - ((Gcvect[i]*(Ca-Ci)/Patm)))

solveCi <- function(Ci)  ((((((Jvect[i]/4)*((Ci-gammastar)/(Ci+(2*gammastar))))) + (((Vmax25*(Ci-gammastar))/(Ci + Kc*(1+(Oa/Ko))))) - (((((((Jvect[i]/4)*((Ci-gammastar)/(Ci+(2*gammastar))))) + (((Vmax25*(Ci-gammastar))/(Ci + Kc*(1+(Oa/Ko))))))^(2)) - (4*cprime*(((Jvect[i]/4)*((Ci-gammastar)/(Ci+(2*gammastar)))))*(((Vmax25*(Ci-gammastar))/(Ci + Kc*(1+(Oa/Ko)))))))^(0.5)))/(2*cprime)) - ((Gcvect[i])*(Ca-Ci))/(Patm*1000))

solveCiForFunction <- function(Ci,Jint,gammastar,Vmax25,Kc,Oa,Ko,cprime,Gcint,Ca,Patm)  ((((((Jint/4)*((Ci-gammastar)/(Ci+(2*gammastar))))) + (((Vmax25*(Ci-gammastar))/(Ci + Kc*(1+(Oa/Ko))))) - (((((((Jint/4)*((Ci-gammastar)/(Ci+(2*gammastar))))) + (((Vmax25*(Ci-gammastar))/(Ci + Kc*(1+(Oa/Ko))))))^(2)) - (4*cprime*(((Jint/4)*((Ci-gammastar)/(Ci+(2*gammastar)))))*(((Vmax25*(Ci-gammastar))/(Ci + Kc*(1+(Oa/Ko)))))))^(0.5)))/(2*cprime)) - ((Gcint)*(Ca-Ci))/(Patm*1000))

solveCiForOptimize <- function(Ci,Jint,gammastar,Vmax25,Kc,Oa,Ko,cprime,Gcint,Ca,Patm)  abs(((((((Jint/4)*((Ci-gammastar)/(Ci+(2*gammastar))))) + (((Vmax25*(Ci-gammastar))/(Ci + Kc*(1+(Oa/Ko))))) - (((((((Jint/4)*((Ci-gammastar)/(Ci+(2*gammastar))))) + (((Vmax25*(Ci-gammastar))/(Ci + Kc*(1+(Oa/Ko))))))^(2)) - (4*cprime*(((Jint/4)*((Ci-gammastar)/(Ci+(2*gammastar)))))*(((Vmax25*(Ci-gammastar))/(Ci + Kc*(1+(Oa/Ko)))))))^(0.5)))/(2*cprime)) - ((Gcint)*(Ca-Ci))/(Patm*1000)))

solvedelta <- function(delta)  (((Cp*(0.189*(u/d)^(0.5) + 0.23)*((Ta*(1+delta))-Ta)) + (lambda*E) + (epsilon*sigma*((Ta^4)*(1+(4*delta))) - Rabs)))

solveTlForFunction <- function(Tl,rho,Cp,u,Ta,lambda,E,epsilon,sigma,Rabs)  ((rho*Cp*(1/(200*(sqrt(u*.05))))*(Tl-Ta)) + (lambda*E) + (epsilon*sigma*((Tl)^(4))) - Rabs)

solveTl <- function(Tl)  ((rho*Cp*(1/(200*(sqrt(u*.05))))*(Tl-Ta)) + (lambda*E) + (epsilon*sigma*((Tl)^(4))) - Rabs)

solvePStorageInit <- function(PStorageInit,Pcrit,k1Crown,k2Crown,PsiSoil)  abs((Pcrit)/(exp((-k1Crown+PStorageInit)/(k2Crown)) + 1) - PsiSoil)


A_gross_assim_rate <- function(Ci,Gc,Ca,Patm){
  A = ((Gc*((Ca-Ci)))/(Patm*1000)) 
  return(A)
}



