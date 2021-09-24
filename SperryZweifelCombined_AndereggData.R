rm(list = ls())
if(!is.null(dev.list())) dev.off()
setwd("/Users/Trent Robinett/Documents/Notre Dame/Medvigy Lab/Sperry_etal_2017_code")
library(FAdist)
library(ggplot2)
library("rootSolve")
library("pracma")
library("MLmetrics")
library("yardstick")
source("/Users/Trent Robinett/Documents/Notre Dame/Medvigy Lab/Sperry_etal_2017_code/Sperryfuncs.R")
source("/Users/Trent Robinett/Documents/Notre Dame/Medvigy Lab/Sperry_etal_2017_code/SperryMain_Function.R")
source("/Users/Trent Robinett/Documents/Notre Dame/Medvigy Lab/Sperry_etal_2017_code/Optimalkmax_V2.R")
source("/Users/Trent Robinett/Documents/Notre Dame/Medvigy Lab/Sperry_etal_2017_code/TimeMatcher_Function.R")
source("/Users/Trent Robinett/Documents/Notre Dame/Medvigy Lab/Sperry_etal_2017_code/IncreaseDataPoints_function_V3.R")
source("/Users/Trent Robinett/Documents/Notre Dame/Medvigy Lab/Sperry_etal_2017_code/CombinedModel_function_AndereggData.R")
treedata = read.csv(file="Anderegg_etal_2018_data_timeinc_V2.csv",stringsAsFactors = FALSE)
speciesdata = read.csv(file="Anderegg_etal_2018_data_timeinc_V2_speciesonly_V2.csv")
speciesdata
#set intervale to run model over
minperstep = 15
nstepsperhour = 60/minperstep

#Calculate kmax for each species

kmaxvect = rep(0,nrow(speciesdata))

for(i in 1:(nrow(speciesdata))){
  VC = c(speciesdata$Weib_b[i],speciesdata$Weib_c[i])
  Vmax25 = speciesdata$Vmax25[i]
  Pcmin = speciesdata$PcMin[i]
  Ci = 28
  D = 1
  Ps = 0
  TA = 25
  kmax = OptimalkmaxV2(VC,Vmax25,Pcmin)
  kmaxvect[i] = kmax
}

cbind(speciesdata,kmaxvect)

z=1
for(i in 1:((nrow(treedata)-1))){
  if(treedata$Species[i] == treedata$Species[(i+1)]){
    treedata$kmax[i] = kmaxvect[z]
  }else{
    treedata$kmax[i] = kmaxvect[z]
    z = z+1
  }
}
treedata$kmax[i+1] = kmaxvect[z]


S = 1
for(i in 1:(nrow(treedata) - 1)){
  if(treedata$Species[i] == treedata$Species[i+1]){
    treedata$SpeciesNumber[i] = S
  }else if (treedata$Species[i] != treedata$Species[i+1]){
    treedata$SpeciesNumber[i] = S
    S = S + 1
  }
}
treedata$SpeciesNumber[i+1] = S


#assign a day # to each row in treedata
M = 1
DateNumber = rep(0,nrow(treedata))
for(i in 1:((nrow(treedata))-1)){
  if(treedata$SpeciesNumber[i] == treedata$SpeciesNumber[i+1]){
    if(treedata$Date[i] == treedata$Date[i+1]){
      DateNumber[i] = M
    }else{
      DateNumber[i] = M
      M = M + 1
    }
  }else{
    if(treedata$Date[i] == treedata$Date[i+1]){
      DateNumber[i] = M
    }else{
      DateNumber[i] = M
      M = 1
    }
  }
}
DateNumber[nrow(treedata)] = M

treedata = cbind(treedata,DateNumber)


#convert times to hours


nrowstreedata = nrow(treedata)


TimeMatcher = TimeMatcher()


W = 0
for(i in 1:(nrow(treedata))){
  for(j in 1:(nrow(TimeMatcher))){
    if(treedata$Time[i] == TimeMatcher$String[j]){
      W = j
    }
  }
  treedata$TimeInt[i] = TimeMatcher$TimeinMins[W]
}


#remove any duplicate times
badrow = 0
for(i in 1:(nrow(treedata) - 1)){
  if(treedata$TimeInt[i] == treedata$TimeInt[i+1]){
    badrow = c(badrow,i)
  }
}
badrow = badrow[2:length(badrow)]
for(i in 1:length(badrow)){
  j = length(badrow) - (i-1)
  treedata <- treedata[-c(badrow[j]),]
}


#remove any duplicate measurments within same time frame
badrow = 0
for(i in 1:(nrow(treedata) - 1)){
  if((treedata$TimeInt[i]+minperstep) > treedata$TimeInt[i+1] && treedata$DateNumber[i] == treedata$DateNumber[i+1]){
    badrow = c(badrow,i)
  }
}
if(max(badrow) > 0){
  badrow = badrow[2:length(badrow)]
  for(i in 1:length(badrow)){
    j = length(badrow) - (i-1)
    treedata <- treedata[-c(badrow[j]),]
  }
}



Cavect = 0
Qvect = 0
Psvect = 0
TAvect = 0
humidityvect = 0
VC1vect = 0
VC2vect = 0
Vmax25vect = 0
Pcritvect = 0
Tleafvect = 0
VPDvect = 0
kmaxvect = 0
timevect = 0
datevect = 0
Jmax25vect = 0
cvect = 0
Patmvect = 0
speciesvect = 0


#####for the first species
for(i in 1:(nrow(treedata))){
  if(treedata$SpeciesNumber[i] == 1){
    if(treedata$Patm[i] != -9999){
      Cavect = c(Cavect,treedata$CO2S[i]/1000000*treedata$Patm[i]*1000)
    }else{
      Cavect = c(Cavect,treedata$CO2S[i]/1000000*101.300*1000) #is this an OK assumption??
    }
    
    #Psvect = c(Psvect,treedata$LWPpredawn[i])
    Psvect = c(Psvect,-1)
    
    Qvect = c(Qvect,treedata$PARin[i])
    
    if(treedata$Tair[i] != -9999){
      TAvect = c(TAvect,treedata$Tair)
    }else{
      TAvect = c(TAvect,25)
    }
    
    if(treedata$RH[i] != -9999){
      humidityvect = c(humidityvect,treedata$RH[i]/100)
    }else{
      humidityvect = c(humidityvect,0.5)
    }
    
    if(treedata$Patm[i] != -9999){
      Patmvect = c(humidityvect,treedata$Patm[i])
    }else{
      Patmvect = c(Patmvect,101.3)
    }
    
    VC1vect = c(VC1vect,treedata$Weib_b[i])
    
    VC2vect = c(VC2vect,treedata$Weib_c[i])
    
    Vmax25vect = c(Vmax25vect,treedata$Vmax25[i])
    
    Pcritvect = c(Pcritvect,treedata$PcMin[i])
    
    Tleafvect = c(Tleafvect,treedata$Tleaf[i])
    
    VPDvect = c(VPDvect,treedata$VPD[i])
    
    kmaxvect = c(kmaxvect,treedata$kmax[i])
    
    timevect = c(timevect,treedata$TimeInt[i])
    
    datevect = c(datevect,treedata$DateNumber[i])
    
    Jmax25vect = c(Jmax25vect,(treedata$Vmax25[i])*1.67)
    
    cvect = c(cvect,0.9)
    
    speciesvect = c(speciesvect,treedata$Species[i])
    
    
    measuredA = treedata$Photo[i]
    measuredE = treedata$Trmmol[i]
    measuredPc = treedata$LWP[i]
    measuredGw = treedata$Cond[i]
    measuredCi = treedata$Ci[i]/1000000*101.3*1000
    
  }
}

Cavect = Cavect[2:length(Cavect)]
Psvect = Psvect[2:length(Psvect)]
Qvect = Qvect[2:length(Qvect)]
TAvect = TAvect[2:length(TAvect)]
humidityvect = humidityvect[2:length(humidityvect)]
VC1vect = VC1vect[2:length(VC1vect)]
VC2vect = VC2vect[2:length(VC2vect)]
Vmax25vect = Vmax25vect[2:length(Vmax25vect)]
Pcritvect = Pcritvect[2:length(Pcritvect)]
Tleafvect = Tleafvect[2:length(Tleafvect)]
VPDvect = VPDvect[2:length(VPDvect)]
kmaxvect = kmaxvect[2:length(kmaxvect)]
timevect = timevect[2:length(timevect)]
datevect = datevect[2:length(datevect)]
Jmax25vect = Jmax25vect[2:length(Jmax25vect)]
Patmvect = Patmvect[2:length(Patmvect)]
cvect = cvect[2:length(cvect)]
speciesvect = speciesvect[2:length(speciesvect)]



IncreaseDataReturned = IncreaseDataPoints(Cavect,Psvect,Qvect,TAvect,humidityvect,VC1vect,VC2vect,Vmax25vect,Pcritvect,Tleafvect,VPDvect,kmaxvect,timevect,datevect,minperstep,Jmax25vect,cvect,Patmvect)


times = rep(0,nrow(IncreaseDataReturned))
w = 1
for(i in 1:nrow(IncreaseDataReturned)){
  times[i] = w*minperstep
  w = w+1
}

IncreaseDataReturned = cbind(IncreaseDataReturned,times)



###Test what parameters give reasonable results
#test k1, k2, R, totalleafarea


PCrownMax = speciesdata$PCrownMax[1]
LAI = speciesdata$LAI[1]
footprint = 3
totalleafarea = 2
k1 = PCrownMax/2
k2 = PCrownMax/10
RXylem = 1


OutTrent = data.frame("Enostoragechangevect" = 1, "PercEchangevect" = 1,
                 "storageflowmmolchangevect" = 1, "negPcchangevect" = 1,
                 "Echangevect" = 1,"Gwchangevect" = 1,"Gcchangevect" = 1,
                 "Cichangevect" = 1,"Achangevect" = 1,"daychangevect" = 1,
                 "hourchangevect" = 1,"k1" = 1,"k2" = 1, "R" = 1, "totalleafarea" = 1)

OutSperry = data.frame("Enostoragechangevect" = 1, "PercEchangevect" = 1,
                      "storageflowmmolchangevect" = 1, "negPcchangevect" = 1,
                      "Echangevect" = 1,"Gwchangevect" = 1,"Gcchangevect" = 1,
                      "Cichangevect" = 1,"Achangevect" = 1,"daychangevect" = 1,
                      "hourchangevect" = 1,"k1" = 1,"k2" = 1, "R" = 1, "totalleafarea" = 1)


Returned = CombinedModel(IncreaseDataReturned,nstepsperhour,PCrownMax,k1,k2,RXylem,totalleafarea)
OutTrent = rbind(OutTrent,Returned)
OutTrent = OutTrent[-1,]







# RXylem = 10000000
# Returned = CombinedModel(IncreaseDataReturned,nstepsperhour,PCrownMax,k1,k2,RXylem,totalleafarea)
# OutSperry = rbind(OutSperry,Returned)
# OutSperry = OutSperry[-1,]


hello = linspace(1,nrow(OutTrent),nrow(OutTrent))
#plot(hello,OutTrent$Echangevect)
#plot(hello,OutTrent$PercEchangevect)
#plot(hello,OutTrent$storageflowmmolchangevect)
#plot(hello,OutTrent$Enostoragechangevect)





###START OF ANALYSIS

ticker = 0
i = 1
while(ticker == 0){
  if(treedata$Species[i] == treedata$Species[i+1]){
    ticker = 0
    i = i+1
  }else{
    ticker = i
  }
}

treedataspecies = treedata[1:ticker,]

nrowspecies = nrow(treedataspecies)

lowpoint = 1
highpoint = nrow(OutTrent)



####for my model
ApredTrent = rep(0,nrowspecies)
EpredTrent = rep(0,nrowspecies)
PcpredTrent = rep(0,nrowspecies)
Aact = rep(0,nrowspecies)
Eact = rep(0,nrowspecies)
Pcact = rep(0,nrowspecies)

count = 1
for(i in lowpoint:highpoint){
  OutTime = (OutTrent$daychangevect[i]*1440) + (OutTrent$hourchangevect[i]*60)
  treeTime = ((treedataspecies$DateNumber[count]-1)*1440) + (treedataspecies$TimeInt[count])
  treeTimelow = treeTime - (60/nstepsperhour)
  
  
  if((OutTime >= (treeTimelow - 0.0001)) && (treeTime > OutTime) && count <= nrowspecies){
    ApredTrent[count] = OutTrent$Achangevect[i]
    Aact[count] = treedataspecies$Photo[count]
    EpredTrent[count] = OutTrent$Echangevect[i]
    Eact[count] = treedataspecies$Trmmol[count]
    PcpredTrent[count] = -OutTrent$negPcchangevect[i]
    Pcact[count] = treedataspecies$LWP[count]
    
    count = count + 1
  }
  
  
}


AllpredTrent = c(ApredTrent,EpredTrent,PcpredTrent)
Allact = c(Aact,Eact,Pcact)


mapePcTrent = mape_vec(Pcact,PcpredTrent)
mapeETrent = mape_vec(Eact,EpredTrent)
mapeATrent = mape_vec(Aact,ApredTrent)
mapeoverallTrent = mape_vec(Allact,AllpredTrent)

mapePc2ndTrent = mape_vec(Pcact[150:305],PcpredTrent[150:305])
mapeE2ndTrent = mape_vec(Eact[150:305],EpredTrent[150:305])
mapeA2ndTrent = mape_vec(Aact[150:305],ApredTrent[150:305])


AavgerrorTrent = sum((ApredTrent-Aact))/length(ApredTrent)
EavgerrorTrent = sum((EpredTrent-Eact))/length(EpredTrent)
PcavgerrorTrent = sum((PcpredTrent-Pcact))/length(PcpredTrent)


####for sperry model

out = data.frame("Pccalc"=1,"measuredPc"=2,"Ecalc"=3,"measuredE"=4,"Gwcalc"=5,"measuredGw"=6,"Cicalc"=7,"measuredCi"=8,"Acalc"=9,"measuredA"=10)



#run over all elements of Anderegg Data
for(i in 1:(nrowspecies)){
  if(treedataspecies$Patm[i] != -9999){
    Ca = treedataspecies$CO2S[i]/1000000*treedataspecies$Patm[i]*1000
  }else{
    Ca = treedataspecies$CO2S[i]/1000000*101.300*1000 #is this an OK assumption??
  }
  Ps = treedataspecies$LWPpredawn[i]
  Q = treedataspecies$PARin[1]
  if(treedataspecies$Tair[1] != -9999){
    TA = treedataspecies$Tair 
  }else{
    TA = 25
  }
  if(treedataspecies$RH[i] != -9999){
    humidity = treedataspecies$RH[i]/100
  }else{
    humidity = 0.5
  }
  
  VC = rep(0,2)
  VC[1] = treedataspecies$Weib_b[i]
  VC[2] = treedataspecies$Weib_c[i]
  
  Vmax25 = treedataspecies$Vmax25[i]
  
  Pcrit = treedataspecies$PcMin[i]
  
  Tleaf = treedataspecies$Tleaf[i]
  
  VPD = treedataspecies$VPD[i]
  
  kmax = treedataspecies$kmax[i]
  
  measuredA = treedataspecies$Photo[i]
  measuredE = treedataspecies$Trmmol[i]
  measuredPc = treedataspecies$LWP[i]
  measuredGw = treedataspecies$Cond[i]
  measuredCi = treedataspecies$Ci[i]/1000000*101.3*1000 #is this an OK assumption??
  
  SperryModelReturned = SperryModel(Ca,Ps,Q,TA,humidity,VC,Vmax25,Pcrit,Tleaf,VPD,kmax)
  
  Pccalc = SperryModelReturned[1]
  Ecalc = SperryModelReturned[2]
  Gwcalc = SperryModelReturned[3]
  Gccalc = SperryModelReturned[4]
  Cicalc = SperryModelReturned[5]
  Acalc = SperryModelReturned[6]
  
  out = rbind(out,list(Pccalc,measuredPc,Ecalc,measuredE,Gwcalc,measuredGw,Cicalc,measuredCi,Acalc,measuredA))
  print(i)
}

OutSperry <- out[-c(1),]



ApredSperry = OutSperry$Acalc
EpredSperry = OutSperry$Ecalc
PcpredSperry = OutSperry$Pccalc
AllpredSperry = c(ApredSperry,EpredSperry,PcpredSperry)

mapePcSperry = mape_vec(Pcact,PcpredSperry)
mapeESperry = mape_vec(Eact,EpredSperry)
mapeASperry = mape_vec(Aact,ApredSperry)
mapeoverallSperry = mape_vec(Allact,AllpredSperry)

mapePc2ndSperry = mape_vec(Pcact[150:305],PcpredSperry[150:305])
mapeE2ndSperry = mape_vec(Eact[150:305],EpredSperry[150:305])
mapeA2ndSperry = mape_vec(Aact[150:305],ApredSperry[150:305])


AavgerrorSperry = sum((ApredSperry-Aact))/length(ApredSperry)
EavgerrorSperry = sum((EpredSperry-Eact))/length(EpredSperry)
PcavgerrorSperry = sum((PcpredSperry-Pcact))/length(PcpredSperry)



#Plot
Sperry = rep("Sperry",nrowspecies)
Data = rep("Data",nrowspecies)
Trent = rep("Trent",nrowspecies)

datacomparison = data.frame("E"=c(treedataspecies$Trmmol,EpredSperry,EpredTrent),
                            "model"=c(Data,Sperry,Trent),
                            "A"=c(treedataspecies$Photo,ApredSperry,ApredTrent),
                            "Pc"=c(treedataspecies$LWP,PcpredSperry,PcpredTrent))


statcomparison = data.frame("Model"=c("Sperry","Sperry","Sperry","Sperry",
                                      "Trent","Trent","Trent","Trent"),
                            "Stat"=c("A","Pc","E","Overall",
                                     "A","Pc","E","Overall"),
                            "Value"= c(mapeASperry,mapePcSperry,mapeESperry,mapeoverallSperry,
                                      mapeATrent,mapePcTrent,mapeETrent,mapeoverallTrent))
                            

rows = linspace(1,nrowspecies,nrowspecies)
datacomparison = cbind(datacomparison,rows)

pdf("TrentSperryEComparison.pdf")
ggplot(datacomparison, aes(x=rows, y=E, color=model)) +
  geom_point() +
  xlab("Observation Number") + ylab("Transpiration (mmol/s/m^2)") +
  theme(
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    legend.title = element_text(size=20),
    legend.text = element_text(size=14)
  )
dev.off()

#pdf("TrentSperryAComparison.pdf",size=15)
ggplot(datacomparison, aes(x=rows, y=A, color=model)) +
  geom_point()
#dev.off()

#pdf("TrentSperryPcComparison.pdf",size=15)
ggplot(datacomparison, aes(x=rows, y=Pc, color=model)) +
  geom_point()
#dev.off()

pdf("TrentSperryOverallComparison.pdf")
ggplot(data=statcomparison, aes(x=Stat, y=Value, fill=Model)) +
  geom_bar(stat="identity", position=position_dodge()) +
  xlab("Model Output") + ylab("Mean Absolute Percent Error") +
  scale_x_discrete(limits=c("E", "A", "Pc", "Overall"))+
  theme(
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    legend.title = element_text(size=20),
    legend.text = element_text(size=14)
  )
dev.off()








