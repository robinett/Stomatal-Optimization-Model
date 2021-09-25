IncreaseDataPoints <- function(Cavect,Psvect,Qvect,TAvect,humidityvect,VC1vect,VC2vect,Vmax25vect,Pcritvect,Tleafvect,VPDvect,kmaxvect,timevect,datevect,minperstep,Jmax25,c,Patm){
  
  nstepsperhour = 60/minperstep
  nstepsperday = nstepsperhour*24
  
  ndays = (max(datevect) - min(datevect) + 1)

  filledCavect = rep(0,(ndays*nstepsperday))
  filledPsvect = rep(0,(ndays*nstepsperday))
  filledQvect = rep(0,(ndays*nstepsperday))
  filledTAvect = rep(0,(ndays*nstepsperday))
  filledhumidityvect = rep(0,(ndays*nstepsperday))
  filledVC1vect = rep(0,(ndays*nstepsperday))
  filledVC2vect = rep(0,(ndays*nstepsperday))
  filledVmax25vect = rep(0,(ndays*nstepsperday))
  filledPcritvect = rep(0,(ndays*nstepsperday))
  filledTleafvect = rep(0,(ndays*nstepsperday))
  filledVPDvect = rep(0,(ndays*nstepsperday))
  filledkmaxvect = rep(0,(ndays*nstepsperday))
  filledJmax25vect = rep(0,(ndays*nstepsperday))
  filledcvect = rep(0,(ndays*nstepsperday))
  filledPatmvect = rep(0,(ndays*nstepsperday))
  filledspeciesname = rep(0,(ndays*nstepsperday))
  
  datevect = c(datevect,0)
  
  
  TimeCount = 1
  for(j in 1:(max(datevect) - min(datevect[1:(length(datevect) - 1)]) + 1)){
    print(j)
    for(i in 1:nstepsperday){
      if(j == 1){#for first day
        if(timevect[TimeCount] >= (i*minperstep + minperstep) && TimeCount == 1){ #before the first recorded point of the day
          filledCavect[(j-1)*nstepsperday + i] = Cavect[1]
          filledPsvect[(j-1)*nstepsperday + i] = Psvect[1]
          filledQvect[(j-1)*nstepsperday + i] = 0
          filledTAvect[(j-1)*nstepsperday + i] = TAvect[1]
          filledhumidityvect[(j-1)*nstepsperday + i] = humidityvect[1]
          filledVC1vect[(j-1)*nstepsperday + i] = VC1vect[1]
          filledVC2vect[(j-1)*nstepsperday + i] = VC2vect[1]
          filledVmax25vect[(j-1)*nstepsperday + i] = Vmax25vect[1]
          filledPcritvect[(j-1)*nstepsperday + i] = Pcritvect[1]
          filledTleafvect[(j-1)*nstepsperday + i] = Tleafvect[1]
          filledVPDvect[(j-1)*nstepsperday + i] = VPDvect[1]
          filledkmaxvect[(j-1)*nstepsperday + i] = kmaxvect[1]
          filledJmax25vect[(j-1)*nstepsperday + i] = Jmax25vect[1]
          filledcvect[(j-1)*nstepsperday + i] = cvect[1]
          filledPatmvect[(j-1)*nstepsperday + i] = Patmvect[1]
          filledspeciesname[(j-1)*nstepsperday + i] = speciesvect[1]
        }else if((timevect[TimeCount] >= (i*minperstep)) && (timevect[TimeCount] < (i*minperstep + minperstep))){ #if there is a recorded point within desired time range
          filledCavect[(j-1)*nstepsperday + i] = Cavect[TimeCount]
          filledPsvect[(j-1)*nstepsperday + i] = Psvect[TimeCount]
          filledQvect[(j-1)*nstepsperday + i] = Qvect[TimeCount]
          filledTAvect[(j-1)*nstepsperday + i] = TAvect[TimeCount]
          filledhumidityvect[(j-1)*nstepsperday + i] = humidityvect[TimeCount]
          filledVC1vect[(j-1)*nstepsperday + i] = VC1vect[TimeCount]
          filledVC2vect[(j-1)*nstepsperday + i] = VC2vect[TimeCount]
          filledVmax25vect[(j-1)*nstepsperday + i] = Vmax25vect[TimeCount]
          filledPcritvect[(j-1)*nstepsperday + i] = Pcritvect[TimeCount]
          filledTleafvect[(j-1)*nstepsperday + i] = Tleafvect[TimeCount]
          filledVPDvect[(j-1)*nstepsperday + i] = VPDvect[TimeCount]
          filledkmaxvect[(j-1)*nstepsperday + i] = kmaxvect[TimeCount]
          filledJmax25vect[(j-1)*nstepsperday + i] = Jmax25vect[TimeCount]
          filledcvect[(j-1)*nstepsperday + i] = cvect[TimeCount]
          filledPatmvect[(j-1)*nstepsperday + i] = Patmvect[TimeCount]
          filledspeciesname[(j-1)*nstepsperday + i] = speciesvect[TimeCount]
          
          TimeCount = TimeCount + 1
        }else if(timevect[TimeCount] >= (i*minperstep + minperstep) && TimeCount != 1){ #linear fit when no recorded point in desire time range
          prevtime = timevect[TimeCount-1]
          nexttime = timevect[TimeCount]
          diff = nexttime - prevtime
          idiff = i*minperstep - prevtime
          
          Caslope = (Cavect[TimeCount] - Cavect[TimeCount-1])/(nexttime-prevtime)
          filledCavect[(j-1)*nstepsperday + i] = Cavect[TimeCount-1] + idiff*Caslope
          
          Psslope = (Psvect[TimeCount] - Psvect[TimeCount-1])/(nexttime-prevtime)
          filledPsvect[(j-1)*nstepsperday + i] = Psvect[TimeCount-1] + idiff*Psslope
          
          Qslope = (Qvect[TimeCount] - Qvect[TimeCount-1])/(nexttime-prevtime)
          filledQvect[(j-1)*nstepsperday + i] = Qvect[TimeCount-1] + idiff*Qslope
          
          TAslope = (TAvect[TimeCount] - TAvect[TimeCount-1])/(nexttime-prevtime)
          filledTAvect[(j-1)*nstepsperday + i] = TAvect[TimeCount-1] + idiff*TAslope
          
          humidityslope = (humidityvect[TimeCount] - humidityvect[TimeCount-1])/(nexttime-prevtime)
          filledhumidityvect[(j-1)*nstepsperday + i] = humidityvect[TimeCount-1] + idiff*humidityslope
          
          VC1slope = (VC1vect[TimeCount] - VC1vect[TimeCount-1])/(nexttime-prevtime)
          filledVC1vect[(j-1)*nstepsperday + i] = VC1vect[TimeCount-1] + idiff*VC1slope
          
          VC2slope = (VC2vect[TimeCount] - VC2vect[TimeCount-1])/(nexttime-prevtime)
          filledVC2vect[(j-1)*nstepsperday + i] = VC2vect[TimeCount-1] + idiff*VC2slope
          
          Vmax25slope = (Vmax25vect[TimeCount] - Vmax25vect[TimeCount-1])/(nexttime-prevtime)
          filledVmax25vect[(j-1)*nstepsperday + i] = Vmax25vect[TimeCount-1] + idiff*Vmax25slope
          
          Pcritslope = (Pcritvect[TimeCount] - Pcritvect[TimeCount-1])/(nexttime-prevtime)
          filledPcritvect[(j-1)*nstepsperday + i] = Pcritvect[TimeCount-1] + idiff*Pcritslope
          
          Tleafslope = (Tleafvect[TimeCount] - Tleafvect[TimeCount-1])/(nexttime-prevtime)
          filledTleafvect[(j-1)*nstepsperday + i] = Tleafvect[TimeCount-1] + idiff*Tleafslope
          
          VPDslope = (VPDvect[TimeCount] - VPDvect[TimeCount-1])/(nexttime-prevtime)
          filledVPDvect[(j-1)*nstepsperday + i] = VPDvect[TimeCount-1] + idiff*VPDslope
          
          kmaxslope = (kmaxvect[TimeCount] - kmaxvect[TimeCount-1])/(nexttime-prevtime)
          filledkmaxvect[(j-1)*nstepsperday + i] = kmaxvect[TimeCount-1] + idiff*kmaxslope
          
          Jmax25slope = (Jmax25vect[TimeCount] - Jmax25vect[TimeCount-1])/(nexttime-prevtime)
          filledJmax25vect[(j-1)*nstepsperday + i] = Jmax25vect[TimeCount-1] + idiff*Jmax25slope
          
          cslope = (cvect[TimeCount] - cvect[TimeCount-1])/(nexttime-prevtime)
          filledcvect[(j-1)*nstepsperday + i] = cvect[TimeCount-1] + idiff*cslope
          
          Patmslope = (Patmvect[TimeCount] - Patmvect[TimeCount-1])/(nexttime-prevtime)
          filledPatmvect[(j-1)*nstepsperday + i] = Patmvect[TimeCount-1] + idiff*Patmslope
          
          filledspeciesname[(j-1)*nstepsperday + i] = speciesvect[TimeCount]
          
        }else if(TimeCount !=1 && datevect[TimeCount] > datevect[TimeCount-1]){ #after final recorded point in the day
          prevtime = timevect[TimeCount-1]
          nexttime = timevect[TimeCount] + 1440
          diff = nexttime - prevtime
          idiff = (i*minperstep) - prevtime
          
          Caslope = (Cavect[TimeCount] - Cavect[TimeCount-1])/(nexttime-prevtime)
          filledCavect[(j-1)*nstepsperday + i] = Cavect[TimeCount-1] + idiff*Caslope
          
          Psslope = (Psvect[TimeCount] - Psvect[TimeCount-1])/(nexttime-prevtime)
          filledPsvect[(j-1)*nstepsperday + i] = Psvect[TimeCount-1] + idiff*Psslope
          
          filledQvect[(j-1)*nstepsperday + i] = 0
          
          TAslope = (TAvect[TimeCount] - TAvect[TimeCount-1])/(nexttime-prevtime)
          filledTAvect[(j-1)*nstepsperday + i] = TAvect[TimeCount-1] + idiff*TAslope
          
          humidityslope = (humidityvect[TimeCount] - humidityvect[TimeCount-1])/(nexttime-prevtime)
          filledhumidityvect[(j-1)*nstepsperday + i] = humidityvect[TimeCount-1] + idiff*humidityslope
          
          VC1slope = (VC1vect[TimeCount] - VC1vect[TimeCount-1])/(nexttime-prevtime)
          filledVC1vect[(j-1)*nstepsperday + i] = VC1vect[TimeCount-1] + idiff*VC1slope
          
          VC2slope = (VC2vect[TimeCount] - VC2vect[TimeCount-1])/(nexttime-prevtime)
          filledVC2vect[(j-1)*nstepsperday + i] = VC2vect[TimeCount-1] + idiff*VC2slope
          
          Vmax25slope = (Vmax25vect[TimeCount] - Vmax25vect[TimeCount-1])/(nexttime-prevtime)
          filledVmax25vect[(j-1)*nstepsperday + i] = Vmax25vect[TimeCount-1] + idiff*Vmax25slope
          
          Pcritslope = (Pcritvect[TimeCount] - Pcritvect[TimeCount-1])/(nexttime-prevtime)
          filledPcritvect[(j-1)*nstepsperday + i] = Pcritvect[TimeCount-1] + idiff*Pcritslope
          
          Tleafslope = (Tleafvect[TimeCount] - Tleafvect[TimeCount-1])/(nexttime-prevtime)
          filledTleafvect[(j-1)*nstepsperday + i] = Tleafvect[TimeCount-1] + idiff*Tleafslope
          
          VPDslope = (VPDvect[TimeCount] - VPDvect[TimeCount-1])/(nexttime-prevtime)
          filledVPDvect[(j-1)*nstepsperday + i] = VPDvect[TimeCount-1] + idiff*VPDslope
          
          kmaxslope = (kmaxvect[TimeCount] - kmaxvect[TimeCount-1])/(nexttime-prevtime)
          filledkmaxvect[(j-1)*nstepsperday + i] = kmaxvect[TimeCount-1] + idiff*kmaxslope
          
          Jmax25slope = (Jmax25vect[TimeCount] - Jmax25vect[TimeCount-1])/(nexttime-prevtime)
          filledJmax25vect[(j-1)*nstepsperday + i] = Jmax25vect[TimeCount-1] + idiff*Jmax25slope
          
          cslope = (cvect[TimeCount] - cvect[TimeCount-1])/(nexttime-prevtime)
          filledcvect[(j-1)*nstepsperday + i] = cvect[TimeCount-1] + idiff*cslope
          
          Patmslope = (Patmvect[TimeCount] - Patmvect[TimeCount-1])/(nexttime-prevtime)
          filledPatmvect[(j-1)*nstepsperday + i] = Patmvect[TimeCount-1] + idiff*Patmslope
          
          filledspeciesname[(j-1)*nstepsperday + i] = speciesvect[TimeCount]
          
        }
      }else if(j != 1 && datevect[TimeCount] != 0){ #all days that are not the first day
        if(timevect[TimeCount] >= (i*minperstep + minperstep) && timevect[TimeCount] < timevect[TimeCount - 1]){ #before first point in the day = final point in the previous day (besideds light)
          prevtime = timevect[TimeCount-1]
          nexttime = timevect[TimeCount] + 1440
          diff = nexttime - prevtime
          idiff = (i*minperstep + 1440) - prevtime
          
          Caslope = (Cavect[TimeCount] - Cavect[TimeCount-1])/(nexttime-prevtime)
          filledCavect[(j-1)*nstepsperday + i] = Cavect[TimeCount-1] + idiff*Caslope
          
          Psslope = (Psvect[TimeCount] - Psvect[TimeCount-1])/(nexttime-prevtime)
          filledPsvect[(j-1)*nstepsperday + i] = Psvect[TimeCount-1] + idiff*Psslope
          
          filledQvect[(j-1)*nstepsperday + i] = 0
          
          TAslope = (TAvect[TimeCount] - TAvect[TimeCount-1])/(nexttime-prevtime)
          filledTAvect[(j-1)*nstepsperday + i] = TAvect[TimeCount-1] + idiff*TAslope
          
          humidityslope = (humidityvect[TimeCount] - humidityvect[TimeCount-1])/(nexttime-prevtime)
          filledhumidityvect[(j-1)*nstepsperday + i] = humidityvect[TimeCount-1] + idiff*humidityslope
          
          VC1slope = (VC1vect[TimeCount] - VC1vect[TimeCount-1])/(nexttime-prevtime)
          filledVC1vect[(j-1)*nstepsperday + i] = VC1vect[TimeCount-1] + idiff*VC1slope
          
          VC2slope = (VC2vect[TimeCount] - VC2vect[TimeCount-1])/(nexttime-prevtime)
          filledVC2vect[(j-1)*nstepsperday + i] = VC2vect[TimeCount-1] + idiff*VC2slope
          
          Vmax25slope = (Vmax25vect[TimeCount] - Vmax25vect[TimeCount-1])/(nexttime-prevtime)
          filledVmax25vect[(j-1)*nstepsperday + i] = Vmax25vect[TimeCount-1] + idiff*Vmax25slope
          
          Pcritslope = (Pcritvect[TimeCount] - Pcritvect[TimeCount-1])/(nexttime-prevtime)
          filledPcritvect[(j-1)*nstepsperday + i] = Pcritvect[TimeCount-1] + idiff*Pcritslope
          
          Tleafslope = (Tleafvect[TimeCount] - Tleafvect[TimeCount-1])/(nexttime-prevtime)
          filledTleafvect[(j-1)*nstepsperday + i] = Tleafvect[TimeCount-1] + idiff*Tleafslope
          
          VPDslope = (VPDvect[TimeCount] - VPDvect[TimeCount-1])/(nexttime-prevtime)
          filledVPDvect[(j-1)*nstepsperday + i] = VPDvect[TimeCount-1] + idiff*VPDslope
          
          kmaxslope = (kmaxvect[TimeCount] - kmaxvect[TimeCount-1])/(nexttime-prevtime)
          filledkmaxvect[(j-1)*nstepsperday + i] = kmaxvect[TimeCount-1] + idiff*kmaxslope
          
          Jmax25slope = (Jmax25vect[TimeCount] - Jmax25vect[TimeCount-1])/(nexttime-prevtime)
          filledJmax25vect[(j-1)*nstepsperday + i] = Jmax25vect[TimeCount-1] + idiff*Jmax25slope
          
          cslope = (cvect[TimeCount] - cvect[TimeCount-1])/(nexttime-prevtime)
          filledcvect[(j-1)*nstepsperday + i] = cvect[TimeCount-1] + idiff*cslope
          
          Patmslope = (Patmvect[TimeCount] - Patmvect[TimeCount-1])/(nexttime-prevtime)
          filledPatmvect[(j-1)*nstepsperday + i] = Patmvect[TimeCount-1] + idiff*Patmslope
          
          filledspeciesname[(j-1)*nstepsperday + i] = speciesvect[TimeCount]

        }else if((timevect[TimeCount] >= (i*minperstep)) && (timevect[TimeCount] < (i*minperstep + minperstep))){ #if there is a recorded point within desired time range
          filledCavect[(j-1)*nstepsperday + i] = Cavect[TimeCount]
          filledPsvect[(j-1)*nstepsperday + i] = Psvect[TimeCount]
          filledQvect[(j-1)*nstepsperday + i] = Qvect[TimeCount]
          filledTAvect[(j-1)*nstepsperday + i] = TAvect[TimeCount]
          filledhumidityvect[(j-1)*nstepsperday + i] = humidityvect[TimeCount]
          filledVC1vect[(j-1)*nstepsperday + i] = VC1vect[TimeCount]
          filledVC2vect[(j-1)*nstepsperday + i] = VC2vect[TimeCount]
          filledVmax25vect[(j-1)*nstepsperday + i] = Vmax25vect[TimeCount]
          filledPcritvect[(j-1)*nstepsperday + i] = Pcritvect[TimeCount]
          filledTleafvect[(j-1)*nstepsperday + i] = Tleafvect[TimeCount]
          filledVPDvect[(j-1)*nstepsperday + i] = VPDvect[TimeCount]
          filledkmaxvect[(j-1)*nstepsperday + i] = kmaxvect[TimeCount]
          filledJmax25vect[(j-1)*nstepsperday + i] = Jmax25vect[TimeCount]
          filledcvect[(j-1)*nstepsperday + i] = cvect[TimeCount]
          filledPatmvect[(j-1)*nstepsperday + i] = Patmvect[TimeCount]
          
          filledspeciesname[(j-1)*nstepsperday + i] = speciesvect[TimeCount]
          
          TimeCount = TimeCount + 1
        }else if((timevect[TimeCount] >= (i*minperstep + minperstep) && timevect[TimeCount] > timevect[TimeCount-1])){ #linear fit when no recorded point in desire time range
          prevtime = timevect[TimeCount-1]
          nexttime = timevect[TimeCount]
          diff = nexttime - prevtime
          idiff = i*minperstep - prevtime
          
          Caslope = (Cavect[TimeCount] - Cavect[TimeCount-1])/(nexttime-prevtime)
          filledCavect[(j-1)*nstepsperday + i] = Cavect[TimeCount-1] + idiff*Caslope
          
          Psslope = (Psvect[TimeCount] - Psvect[TimeCount-1])/(nexttime-prevtime)
          filledPsvect[(j-1)*nstepsperday + i] = Psvect[TimeCount-1] + idiff*Psslope
          
          Qslope = (Qvect[TimeCount] - Qvect[TimeCount-1])/(nexttime-prevtime)
          filledQvect[(j-1)*nstepsperday + i] = Qvect[TimeCount-1] + idiff*Qslope
          
          TAslope = (TAvect[TimeCount] - TAvect[TimeCount-1])/(nexttime-prevtime)
          filledTAvect[(j-1)*nstepsperday + i] = TAvect[TimeCount-1] + idiff*TAslope
          
          humidityslope = (humidityvect[TimeCount] - humidityvect[TimeCount-1])/(nexttime-prevtime)
          filledhumidityvect[(j-1)*nstepsperday + i] = humidityvect[TimeCount-1] + idiff*humidityslope
          
          VC1slope = (VC1vect[TimeCount] - VC1vect[TimeCount-1])/(nexttime-prevtime)
          filledVC1vect[(j-1)*nstepsperday + i] = VC1vect[TimeCount-1] + idiff*VC1slope
          
          VC2slope = (VC2vect[TimeCount] - VC2vect[TimeCount-1])/(nexttime-prevtime)
          filledVC2vect[(j-1)*nstepsperday + i] = VC2vect[TimeCount-1] + idiff*VC2slope
          
          Vmax25slope = (Vmax25vect[TimeCount] - Vmax25vect[TimeCount-1])/(nexttime-prevtime)
          filledVmax25vect[(j-1)*nstepsperday + i] = Vmax25vect[TimeCount-1] + idiff*Vmax25slope
          
          Pcritslope = (Pcritvect[TimeCount] - Pcritvect[TimeCount-1])/(nexttime-prevtime)
          filledPcritvect[(j-1)*nstepsperday + i] = Pcritvect[TimeCount-1] + idiff*Pcritslope
          
          Tleafslope = (Tleafvect[TimeCount] - Tleafvect[TimeCount-1])/(nexttime-prevtime)
          filledTleafvect[(j-1)*nstepsperday + i] = Tleafvect[TimeCount-1] + idiff*Tleafslope
          
          VPDslope = (VPDvect[TimeCount] - VPDvect[TimeCount-1])/(nexttime-prevtime)
          filledVPDvect[(j-1)*nstepsperday + i] = VPDvect[TimeCount-1] + idiff*VPDslope
          
          kmaxslope = (kmaxvect[TimeCount] - kmaxvect[TimeCount-1])/(nexttime-prevtime)
          filledkmaxvect[(j-1)*nstepsperday + i] = kmaxvect[TimeCount-1] + idiff*kmaxslope
          
          Jmax25slope = (Jmax25vect[TimeCount] - Jmax25vect[TimeCount-1])/(nexttime-prevtime)
          filledJmax25vect[(j-1)*nstepsperday + i] = Jmax25vect[TimeCount-1] + idiff*Jmax25slope
          
          cslope = (cvect[TimeCount] - cvect[TimeCount-1])/(nexttime-prevtime)
          filledcvect[(j-1)*nstepsperday + i] = cvect[TimeCount-1] + idiff*cslope
          
          Patmslope = (Patmvect[TimeCount] - Patmvect[TimeCount-1])/(nexttime-prevtime)
          filledPatmvect[(j-1)*nstepsperday + i] = Patmvect[TimeCount-1] + idiff*Patmslope
          
          filledspeciesname[(j-1)*nstepsperday + i] = speciesvect[TimeCount]
          
        }else if(datevect[TimeCount] > datevect[TimeCount-1]){ #linear fit between last daylight point of previous day and first daylight point of next day
          prevtime = timevect[TimeCount-1]
          nexttime = timevect[TimeCount] + 1440
          diff = nexttime - prevtime
          idiff = i*minperstep - prevtime
          
          Caslope = (Cavect[TimeCount] - Cavect[TimeCount-1])/(nexttime-prevtime)
          filledCavect[(j-1)*nstepsperday + i] = Cavect[TimeCount-1] + idiff*Caslope
          
          Psslope = (Psvect[TimeCount] - Psvect[TimeCount-1])/(nexttime-prevtime)
          filledPsvect[(j-1)*nstepsperday + i] = Psvect[TimeCount-1] + idiff*Psslope
          
          filledQvect[(j-1)*nstepsperday + i] = 0
          
          TAslope = (TAvect[TimeCount] - TAvect[TimeCount-1])/(nexttime-prevtime)
          filledTAvect[(j-1)*nstepsperday + i] = TAvect[TimeCount-1] + idiff*TAslope
          
          humidityslope = (humidityvect[TimeCount] - humidityvect[TimeCount-1])/(nexttime-prevtime)
          filledhumidityvect[(j-1)*nstepsperday + i] = humidityvect[TimeCount-1] + idiff*humidityslope
          
          VC1slope = (VC1vect[TimeCount] - VC1vect[TimeCount-1])/(nexttime-prevtime)
          filledVC1vect[(j-1)*nstepsperday + i] = VC1vect[TimeCount-1] + idiff*VC1slope
          
          VC2slope = (VC2vect[TimeCount] - VC2vect[TimeCount-1])/(nexttime-prevtime)
          filledVC2vect[(j-1)*nstepsperday + i] = VC2vect[TimeCount-1] + idiff*VC2slope
          
          Vmax25slope = (Vmax25vect[TimeCount] - Vmax25vect[TimeCount-1])/(nexttime-prevtime)
          filledVmax25vect[(j-1)*nstepsperday + i] = Vmax25vect[TimeCount-1] + idiff*Vmax25slope
          
          Pcritslope = (Pcritvect[TimeCount] - Pcritvect[TimeCount-1])/(nexttime-prevtime)
          filledPcritvect[(j-1)*nstepsperday + i] = Pcritvect[TimeCount-1] + idiff*Pcritslope
          
          Tleafslope = (Tleafvect[TimeCount] - Tleafvect[TimeCount-1])/(nexttime-prevtime)
          filledTleafvect[(j-1)*nstepsperday + i] = Tleafvect[TimeCount-1] + idiff*Tleafslope
          
          VPDslope = (VPDvect[TimeCount] - VPDvect[TimeCount-1])/(nexttime-prevtime)
          filledVPDvect[(j-1)*nstepsperday + i] = VPDvect[TimeCount-1] + idiff*VPDslope
          
          kmaxslope = (kmaxvect[TimeCount] - kmaxvect[TimeCount-1])/(nexttime-prevtime)
          filledkmaxvect[(j-1)*nstepsperday + i] = kmaxvect[TimeCount-1] + idiff*kmaxslope
          
          Jmax25slope = (Jmax25vect[TimeCount] - Jmax25vect[TimeCount-1])/(nexttime-prevtime)
          filledJmax25vect[(j-1)*nstepsperday + i] = Jmax25vect[TimeCount-1] + idiff*Jmax25slope
          
          cslope = (cvect[TimeCount] - cvect[TimeCount-1])/(nexttime-prevtime)
          filledcvect[(j-1)*nstepsperday + i] = cvect[TimeCount-1] + idiff*cslope
          
          Patmslope = (Patmvect[TimeCount] - Patmvect[TimeCount-1])/(nexttime-prevtime)
          filledPatmvect[(j-1)*nstepsperday + i] = Patmvect[TimeCount-1] + idiff*Patmslope
          
          filledspeciesname[(j-1)*nstepsperday + i] = speciesvect[TimeCount]
          
        }
       }else if(datevect[TimeCount] == 0){
         filledCavect[(j-1)*nstepsperday + i] = Cavect[TimeCount-1]
         filledPsvect[(j-1)*nstepsperday + i] = Psvect[TimeCount-1]
         filledQvect[(j-1)*nstepsperday + i] = 0
         filledTAvect[(j-1)*nstepsperday + i] = TAvect[TimeCount-1]
         filledhumidityvect[(j-1)*nstepsperday + i] = humidityvect[TimeCount-1]
         filledVC1vect[(j-1)*nstepsperday + i] = VC1vect[TimeCount-1]
         filledVC2vect[(j-1)*nstepsperday + i] = VC2vect[TimeCount-1]
         filledVmax25vect[(j-1)*nstepsperday + i] = Vmax25vect[TimeCount-1]
         filledPcritvect[(j-1)*nstepsperday + i] = Pcritvect[TimeCount-1]
         filledTleafvect[(j-1)*nstepsperday + i] = Tleafvect[TimeCount-1]
         filledVPDvect[(j-1)*nstepsperday + i] = VPDvect[TimeCount-1]
         filledkmaxvect[(j-1)*nstepsperday + i] = kmaxvect[TimeCount-1]
         filledJmax25vect[(j-1)*nstepsperday + i] = Jmax25vect[TimeCount-1]
         filledcvect[(j-1)*nstepsperday + i] = cvect[TimeCount-1]
         filledPatmvect[(j-1)*nstepsperday + i] = Patmvect[TimeCount-1]
         filledspeciesname[(j-1)*nstepsperday + i] = speciesvect[TimeCount]
         
      }
    }  
  }
  
  ReturnFrame = data.frame("filledCavect" = filledCavect,"filledPsvect" = filledPsvect, "filledQvect" = filledQvect, "filledTAvect" = filledTAvect, "filledhumidityvect" = filledhumidityvect, "filledVC1vect" = filledVC1vect, "filledVC2vect" = filledVC2vect, "filledVmax25vect" = filledVmax25vect, "filledPcritvect" = filledPcritvect, "filledTleafvect" = filledTleafvect, "filledVPDvect" = filledVPDvect, "filledkmaxvect" = filledkmaxvect, "filledJmax25vect" = filledJmax25vect, "filledcvect" = filledcvect, "filledPatmvect" = filledPatmvect, "filledspeciesname" = filledspeciesname)
  
  return(ReturnFrame)
}
