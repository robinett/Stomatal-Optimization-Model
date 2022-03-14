rm(list = ls())
if(!is.null(dev.list())) dev.off()
setwd("/Users/Trent Robinett/Documents/Notre Dame/Medvigy Lab/Sperry_etal_2017_code")
library(FAdist)
library(ggplot2)
library(rhdf5)
library("rootSolve")
library("pracma")
library("MLmetrics")
library("yardstick")
library("bigleaf")
library("photosynthesis")
source("/Users/Trent Robinett/Documents/Notre Dame/Medvigy Lab/Sperry_etal_2017_code/Sperryfuncs.R")
source("/Users/Trent Robinett/Documents/Notre Dame/Medvigy Lab/Sperry_etal_2017_code/CombinedModelSapFlow_function.R")
source("/Users/Trent Robinett/Documents/Notre Dame/Medvigy Lab/Sperry_etal_2017_code/ZenithAngle_function.R")
source("/Users/Trent Robinett/Documents/Notre Dame/Medvigy Lab/Sperry_etal_2017_code/Optimalkmax_V2_Aspinwall.R")


nsteps = 1448
nstepsperhour = 2
nlayers = 10


#####LOAD IN METEOROLOGICAL DATA
## Set file name
fname="DF_R1_MET_2010JUL.h5"

## Read in CO2. Units are ppm. There are 1488 values, one for 
## each half-hour of the month (48 per day, 31 days, 48*31=1488) 
## Note: the variables are set up a bit funny, hence the [,1,1] at the end.
co2=h5read(fname,"co2")[,1,1]
dlwrf=h5read(fname,"dlwrf")[,1,1]
nbdsf=h5read(fname,"nbdsf")[,1,1]
nddsf=h5read(fname,"nddsf")[,1,1]
vbdsf=h5read(fname,"vbdsf")[,1,1]
vddsf=h5read(fname,"vddsf")[,1,1]
prate=h5read(fname,"prate")[,1,1]
pres=h5read(fname,"pres")[,1,1]
sh=h5read(fname,"sh")[,1,1]
tmp=h5read(fname,"tmp")[,1,1]
ugrd=h5read(fname,"ugrd")[,1,1]
vgrd=h5read(fname,"vgrd")[,1,1]
totwind = sqrt(ugrd**2+vgrd**2)

#converted atmospheric drivers
vbdsfumol = vbdsf*4.6
vddsfumol = vddsf*4.6
preskPa = pres/1000
tmpC = tmp-273.15

#vulnerability curve parameters
P50 = -2.75
a = 2.3


#duke forest specific variables
lat = 36.0201 #latitude and long. are internet values. need citable values
lon = 78.9830

## See all the variables, for informational purposes:
h5ls(fname)

## co2 is atmospheric CO2, Units are ppm
## dlwrf is the downward longwave radiation, W/m2
## nbdsf is the NIR band shortwave beam radiation, W/m2
## nddsf is the NIR band shortwave diffuse radiation, W/m2
## vbdsf is the PAR (visible) band shortwave downward beam radiation, W/m2
## vddsf is the PAR (visible) band shortwave downward diffuse radiation, W/m2
## prate is precipitation rate; I don't think you'll need this
## pres is pressure in Pa
## sh is specific humidity, grams water vapor per gram air
## tmp is temperature in K
## ugrd is east-west wind and vgrd is north-south wind. Thus, total wind
##     is sqrt(ugrd**2+vgrd**2).(ignoring up-down wind, which is usually small)



#data storage matrices (row is timepoint, column is layer)
storageflowgramschangevect = matrix(0,nsteps,nlayers)
storageflowmmolchangevect = matrix(0,nsteps,nlayers)
Echangevect = matrix(0,nsteps,nlayers)
Tlchangevect = matrix(0,nsteps,nlayers)
Dlchangevect = matrix(0,nsteps,nlayers)
Gwchangevect = matrix(0,nsteps,nlayers)
Gcchangevect = matrix(0,nsteps,nlayers)
Jchangevect = matrix(0,nsteps,nlayers)
Cichangevect = matrix(0,nsteps,nlayers)
Achangevect = matrix(0,nsteps,nlayers)
thetachangevect = matrix(0,nsteps,nlayers)
betachangevect = matrix(0,nsteps,nlayers)
profitchangevect = matrix(0,nsteps,nlayers)
zvect = matrix(0,nsteps,nlayers)
Pvect = matrix(0,nsteps,nlayers)
Pcchangevect = matrix(0,nsteps,nlayers)
PercEchangevect = matrix(0,nsteps,nlayers)
PcRangevect = matrix(0,nsteps,nlayers)
Pcchangevect = matrix(0,nsteps,nlayers)
kcchangevect = matrix(0,nsteps,nlayers)
tempchangevect = matrix(0,nsteps,nlayers)
PStoragevect = matrix(0,nsteps+1,nlayers)
PsiStoragevect = matrix(0,nsteps+1,nlayers)
negPcchangevect = matrix(0,nsteps,nlayers)
Qchangevect = matrix(0,nsteps,nlayers)
Enostoragechangevect = matrix(0,nsteps,nlayers)
daychangevect = matrix(0,nsteps,nlayers)
hourchangevect = matrix(0,nsteps,nlayers)


#Other necessary inputs
Pcrit = -4
PCrownMax = 3000
k1Crown = PCrownMax/2
k2Crown = PCrownMax/10
RXylem = 0.01
lwidth = 0.0025
Ps = -1
Vmax25 = 26

###Initialize the Amount in Storage

interval = c(0,PCrownMax)
PStorageInit = optimize(solvePStorageInit,interval,Pcrit=Pcrit,k1Crown=k1Crown,k2Crown=k2Crown,PsiSoil=Ps)
PStorageInit


#Initialize storage vectors
PStoragevect[1,1:10] = PStorageInit$minimum
Pcchangevect[1,1:10] = -1
PsiStoragevect[1,1:10] = -1

#find optimal kmax for species
kmax = OptimalkmaxV2(P50,a,Vmax25,Pcrit)


for(z in 1:nsteps){
  print(z)
  #####Calculate Radiation at Each Layer
  doy = 182 + (z/48 - (z%%48)/48)
  timeofday = (z%%48)/48*24 - 4
  Ibvect = rep(0,nlayers)
  Idvect = rep(0,nlayers)
  Iscvect = rep(0,nlayers)
  sunlitr = rep(0,nlayers)
  shadedr = rep(0,nlayers)
  fsl = rep(0,nlayers)
  radvect = rep(0,nlayers)
  
  
  x = 1 #ratio of average projected leaf area of canopy elements on the horizontal and vertical surfaces (value of 1 chosen from Schafer et al. 2003)
  tansqz = ZenithAngle(doy,lat,lon)
  tansqznextday = ZenithAngle((doy+1),lat,lon)
  PI = 0.75 #clumping factor
  LAI = 4 #given to me by Dr. Medvigy for loblolly pine
  Lti = LAI/nlayers 
  footprint = 10
  leafarea = LAI*footprint
  leafareai = leafarea/nlayers
  cshoot = 0.609 #from Therezien etal 2007
  alphap = 0.83
  VPD = q.to.VPD(sh[z],tmpC[z],preskPa[z])
  ppm <- set_units(co2[z], "umol/mol")
  P <- set_units(preskPa[z], "kPa")
  co2Pa = ppm2pa(ppm, P)
  co2Pa = as.numeric(co2Pa)
  
  #following loop and if statement convert tansqz from UTC to EST
  timecount = 1
  if(timecount < timeofday){
    timecount = timecount + 1
  }
  
  tansqztime = tansqz[timecount]
  
  #light extinction coefficient
  Kbe = (((x^2) + (tansqztime))^(1/2))/(x + 1.744*(x+1.182)^(-0.733)) ####need to figure out the UTCnconversion dealio
  
  #transmission coefficient for direct light
  Tb = exp((-Kbe*Lti*PI))
  
  #direct beam radiation
  for(i in 1:nlayers){
    Ibvect[i] = Kbe*vbdsfumol[z]*cshoot^(i-1)
  }
  
  #transmission coefficient diffuse light
  Kd = 0.7 #from figure 15.4 in Campbell & Norman 1998, given LAI = 4 and x = 1
  Td = exp((-sqrt(alphap))*Kd*Lti)
  
  #diffuse radiation
  for(i in 1:nlayers){
    Idvect[i] = vddsfumol[z]*Td
  }
  
  #transmission coefficient for scattered radiation
  Ts = exp((-sqrt(alphap))*(Kbe)*(Lti)*(PI))
  
  #scattered radiation
  for(i in 1:nlayers){
    Iscvect[i] = (Ts - Tb)*vbdsfumol[z]
  }
  
  #fraction of foliage receiving sunlight
  for(i in 1:nlayers){
    fsl[i] = ((Kbe*Lti*PI)^(i-1))/(((factorial(i-1))*(1-cshoot))^(i-1))*(exp(((-Kbe)*Lti*PI)/(1-cshoot)))
  }
  
  for(i in 1:nlayers){
    sunlitr[i] = Ibvect[i] + Idvect[i] + Iscvect[i]
    shadedr[i] = Idvect[i] + Iscvect[i]
    
  }
  
  #Schafer et al 2013 uses weighted fraction to determine amount of radiance in each layer
  for(i in 1:nlayers){
    radvect[i] = fsl[i]*sunlitr[i] + fsl[i]*shadedr[i]
  }
  
  ###to do list for parameters
  #convert humidity to VPD
  #make sure that all units are the same
  #input changed vulnerability curve (i think i already kinda set this up somewhere?)
  
  totalstorageflowgrams = 0
  totalE = 0
  avgTl = 0
  avgDl = 0
  totalGw = 0
  avgCi = 0
  totalA = 0
  avgnegPc = 0
  avgpercE = 0
  totalPstorage = 0
  
  
  for(i in 1:nlayers){
    params = data.frame("PCrownMax"=PCrownMax, "k1Crown"=k1Crown, "k2Crown"=k2Crown,
                        "RXylem"=RXylem, "lwidth"=lwidth, "Ca"=co2Pa, "VPD"=VPD, 
                        "Patm"=preskPa[z], "Ps"=Ps, "Q"=radvect[i], "TA"=tmpC[z], 
                        "u"=totwind[z], "Vmax25"=Vmax25, "P50"=P50, "a"=a, "Pcrit"=Pcrit,
                        "PStoragevect"=PStoragevect[z,i], 
                        "PsiStoragevect"=PsiStoragevect[z,i], "kmax"=kmax,
                        "leafarea" = leafareai)
    
    Out = CombinedModel(params)
    
    PStoragevectnext = PStoragevect[z,i] + Out$storageflowgramschangevect/nstepsperhour
    PsiStoragevectnext = ((params$Pcrit)/(exp((-params$k1Crown+PStoragevectnext)/(params$k2Crown))+1))
    
    storageflowgramschangevect[z,i] = Out$storageflowgramschangevect
    storageflowmmolchangevect[z,i] = Out$storageflowmmolchangevect
    Echangevect[z,i] = Out$Echangevect
    Tlchangevect[z,i] = Out$Tlchangevect
    Dlchangevect[z,i] = Out$Dlchangevect
    Gwchangevect[z,i] = Out$Gwchangevect
    Gcchangevect[z,i] = Out$Gcchangevect
    Cichangevect[z,i] = Out$Cichangevect
    Achangevect[z,i] = Out$Achangevect
    negPcchangevect[z,i] = Out$negPcchangevect
    PercEchangevect[z,i] = Out$PercEchangevect
    Qchangevect[z,i] = radvect[i]
    PStoragevect[z+1,i] = PStoragevectnext
    PsiStoragevect[z+1,i] = PsiStoragevectnext
  }
  
}


finalstorageflowgram = rep(0,nsteps)
finalstorageflowmmol = rep(0,nsteps)
finalE = rep(0,nsteps)
finalTl = rep(0,nsteps)
finalDl = rep(0,nsteps)
finalGw = rep(0,nsteps)
finalGc = rep(0,nsteps)
finalCi = rep(0,nsteps)
finalnegPc = rep(0,nsteps)
finalPercE = rep(0,nsteps)
finalPStorage = rep(0,nsteps)
finalPsiStorage = rep(0,nsteps)
finalSapflow = rep(0,nsteps)
finalQ = rep(0,nsteps)

for(j in 1:nsteps){
  finalstorageflowgram[j] = sum(storageflowgramschangevect[j,])
  finalstorageflowmmol[j] = sum(storageflowmmolchangevect[j,])
  finalE[j] = sum(Echangevect[j,])
  finalTl[j] = mean(Tlchangevect[j,])
  finalDl[j] = mean(Dlchangevect[j,])
  finalGw[j] = sum(Gwchangevect[j,])
  finalGc[j] = sum(Gcchangevect[j,])
  finalCi[j] = mean(Cichangevect[j,])
  finalnegPc[j] = mean(negPcchangevect[j,])
  if(finalstorageflowmmol[j]<0){
    finalPercE[j] = finalstorageflowmmol[j]/finalE[j]
  }else{
    finalPercE[j] = 0
  }
  finalPStorage[j] = sum(PStoragevect[j,])
  finalPsiStorage[j] = mean(PsiStoragevect[j,])
  finalSapflow[j] = finalE[j] + finalstorageflowmmol[j]
  finalQ[j] = sum(Qchangevect[j,])
}

hello = linspace(0,31,1448)

datacomparison = data.frame("day"=hello[245:435], "E"=finalE[245:435], 
			    "PercE"=-finalPercE[245:435], "Sapflow"=finalSapflow[245:435])

print('')

# plot(hello,finalE)
# plot(hello,finalstorageflowmmol)
# plot(hello,finalPercE)
# plot(hello,finalPStorage)
# plot(hello,finalPsiStorage)
# plot(hello,finalSapflow)

pdf("DukeForestPercentFromStorage.pdf")
ggplot(datacomparison, aes(x=day, y=PercE)) +
  geom_point(size=3) +
  geom_vline(xintercept = 5.57, color = "green", size = .7) +
  geom_vline(xintercept = 6.59, color = "green", size = .7) +
  geom_vline(xintercept = 7.6, color = "green", size = .7) +
  geom_vline(xintercept = 8.65, color = "green", size = .7) +
  geom_vline(xintercept = 6.17, color = "red", size = .7) +
  geom_vline(xintercept = 7.19, color = "red", size = .7) +
  geom_vline(xintercept = 8.22, color = "red", size = .7) +
  geom_vline(xintercept = 9.25, color = "red", size = .7) +
  xlab("Hours Elapsed") + ylab("% of Transpiration from Storage") +
  #scale_y_continuous(breaks=seq(-1.1,1.5,0.4)) +
  geom_hline(yintercept = 0.4, color = "blue", size = .7) +
  geom_hline(yintercept = 0, color = "black", size = .7) +
  theme(
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14)
  )
dev.off()











# 
# pdf("totalESapflow.pdf", width=15)
# plot(hello[500:1000],finalE[500:1000])
# abline(v=hello[501], col="green")
# abline(v=hello[549], col="green")
# abline(v=hello[597], col="green")
# abline(v=hello[645], col="green")
# abline(v=hello[693], col="green")
# abline(v=hello[741], col="green")
# abline(v=hello[789], col="green")
# abline(v=hello[837], col="green")
# abline(v=hello[885], col="green")
# abline(v=hello[933], col="green")
# abline(v=hello[981], col="green")
# abline(v=hello[530], col="red")
# abline(v=hello[578], col="red")
# abline(v=hello[626], col="red")
# abline(v=hello[674], col="red")
# abline(v=hello[722], col="red")
# abline(v=hello[770], col="red")
# abline(v=hello[818], col="red")
# abline(v=hello[866], col="red")
# abline(v=hello[914], col="red")
# abline(v=hello[962], col="red")
# dev.off()
# 
# pdf("storageflowmmolSapflow.pdf", width=15)
# plot(hello[500:1000],finalstorageflowmmol[500:1000])
# abline(v=hello[501], col="green")
# abline(v=hello[549], col="green")
# abline(v=hello[597], col="green")
# abline(v=hello[645], col="green")
# abline(v=hello[693], col="green")
# abline(v=hello[741], col="green")
# abline(v=hello[789], col="green")
# abline(v=hello[837], col="green")
# abline(v=hello[885], col="green")
# abline(v=hello[933], col="green")
# abline(v=hello[981], col="green")
# abline(v=hello[530], col="red")
# abline(v=hello[578], col="red")
# abline(v=hello[626], col="red")
# abline(v=hello[674], col="red")
# abline(v=hello[722], col="red")
# abline(v=hello[770], col="red")
# abline(v=hello[818], col="red")
# abline(v=hello[866], col="red")
# abline(v=hello[914], col="red")
# abline(v=hello[962], col="red")
# dev.off()
# 
# pdf("PercESapflow.pdf", width=15)
# plot(hello[500:1000],finalPercE[500:1000])
# abline(v=hello[501], col="green")
# abline(v=hello[549], col="green")
# abline(v=hello[597], col="green")
# abline(v=hello[645], col="green")
# abline(v=hello[693], col="green")
# abline(v=hello[741], col="green")
# abline(v=hello[789], col="green")
# abline(v=hello[837], col="green")
# abline(v=hello[885], col="green")
# abline(v=hello[933], col="green")
# abline(v=hello[981], col="green")
# abline(v=hello[530], col="red")
# abline(v=hello[578], col="red")
# abline(v=hello[626], col="red")
# abline(v=hello[674], col="red")
# abline(v=hello[722], col="red")
# abline(v=hello[770], col="red")
# abline(v=hello[818], col="red")
# abline(v=hello[866], col="red")
# abline(v=hello[914], col="red")
# abline(v=hello[962], col="red")
# dev.off()
# 
# pdf("PStorageSapflow.pdf", width=15)
# plot(hello[500:1000],finalPStorage[500:1000])
# abline(v=hello[501], col="green")
# abline(v=hello[549], col="green")
# abline(v=hello[597], col="green")
# abline(v=hello[645], col="green")
# abline(v=hello[693], col="green")
# abline(v=hello[741], col="green")
# abline(v=hello[789], col="green")
# abline(v=hello[837], col="green")
# abline(v=hello[885], col="green")
# abline(v=hello[933], col="green")
# abline(v=hello[981], col="green")
# abline(v=hello[530], col="red")
# abline(v=hello[578], col="red")
# abline(v=hello[626], col="red")
# abline(v=hello[674], col="red")
# abline(v=hello[722], col="red")
# abline(v=hello[770], col="red")
# abline(v=hello[818], col="red")
# abline(v=hello[866], col="red")
# abline(v=hello[914], col="red")
# abline(v=hello[962], col="red")
# dev.off()
# 
# pdf("PsiStorageSapflow.pdf", width=15)
# plot(hello[500:1000],finalPsiStorage[500:1000])
# abline(v=hello[501], col="green")
# abline(v=hello[549], col="green")
# abline(v=hello[597], col="green")
# abline(v=hello[645], col="green")
# abline(v=hello[693], col="green")
# abline(v=hello[741], col="green")
# abline(v=hello[789], col="green")
# abline(v=hello[837], col="green")
# abline(v=hello[885], col="green")
# abline(v=hello[933], col="green")
# abline(v=hello[981], col="green")
# abline(v=hello[530], col="red")
# abline(v=hello[578], col="red")
# abline(v=hello[626], col="red")
# abline(v=hello[674], col="red")
# abline(v=hello[722], col="red")
# abline(v=hello[770], col="red")
# abline(v=hello[818], col="red")
# abline(v=hello[866], col="red")
# abline(v=hello[914], col="red")
# abline(v=hello[962], col="red")
# dev.off()
# 
# pdf("SapflowSapflow.pdf", width=15)
# plot(hello[500:1000],finalSapflow[500:1000])
# abline(v=hello[501], col="green")
# abline(v=hello[549], col="green")
# abline(v=hello[597], col="green")
# abline(v=hello[645], col="green")
# abline(v=hello[693], col="green")
# abline(v=hello[741], col="green")
# abline(v=hello[789], col="green")
# abline(v=hello[837], col="green")
# abline(v=hello[885], col="green")
# abline(v=hello[933], col="green")
# abline(v=hello[981], col="green")
# abline(v=hello[530], col="red")
# abline(v=hello[578], col="red")
# abline(v=hello[626], col="red")
# abline(v=hello[674], col="red")
# abline(v=hello[722], col="red")
# abline(v=hello[770], col="red")
# abline(v=hello[818], col="red")
# abline(v=hello[866], col="red")
# abline(v=hello[914], col="red")
# abline(v=hello[962], col="red")
# dev.off()
# 
# 
