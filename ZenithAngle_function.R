ZenithAngle <- function(doy,lat,lon){


  ## This script determine the incoming solar radiation based on geometry alone.  It does
  ## not account for atmospheric composition, cloudiness, sunspot variability, etc.
  
  ## solar constant
  solar <- 1.3533e3
  
  ## doy is the day of the year.
  #doy <- 60
  
  ## latitude
  #lat <- 25.7
  
  ## longitude
  #lon <- -80.2
  
  ## some constants
  pi2 <- 2*pi
  pio180 <- pi / 180.
  
  
  ################################################
  ## rshort is the incoming solar (shortwave) radiation.  There is one value for
  ## each hour of the day.
  rshort <- rep(0,24)
  tansqz = rep(0,24)
  
  ## Loop over hours
  for(itime in 1:24){
  
    ## utc.sec is the number of seconds that have elapsed in this day.  Note that
    ## it is referenced against UTC time, not local time.
    utc.sec <- (itime-1)*3600.
  
    ## This is converts lat-lon (spherical) coordinates to x-y-z coordinates.
    wnx <- cos(pio180*lat)*cos(pio180*lon)
    wny <- cos(pio180*lat)*sin(pio180*lon)
    wnz <- sin(pio180*lat)
  
    ## Get solar declination.  This is an approximate function.  Don't worry about the details.
    t1 <- pi2 * doy / 366
    declin <- .322003-22.971*cos(t1)-.357898*cos(t1*2)-.14398*cos(t1*3)+3.94638*sin(t1)+.019334*sin(t1*2)+.05928*sin(t1*3)
  
    ## Calculate the time-dependence.  No need to worry about the details.
    t2 <- (279.134+.985647*doy)*pio180
    eqn.of.time <- 5.0323 - 100.976 * sin(t2) + 595.275*sin(t2*2)+3.6858 * sin(t2*3)-12.47*sin(t2*4)-430.847*cos(t2)+12.5024*cos(t2*2)+18.25*cos(t2*3)
  
    sun.longitude <- 180. - 360. * (utc.sec + eqn.of.time) / 86400.
  
    ## This is a unit vector defining the position of the sun.
    ## isn't zenith angle (90-sunz)??? But sunz doesnt change throughout the day?
    sunx <- cos(declin*pio180) * cos(sun.longitude * pio180)
    suny <- cos(declin*pio180) * sin(sun.longitude * pio180)
    sunz <- sin(declin*pio180)
  
    ## This is the dot product between the solar unit vector and your unit vector.
    ## (90-cosz) is zenith angle????
    cosz <- sunx*wnx + suny*wny + sunz*wnz
    secsqz = (1/(cosz))^2
    tansqz[itime] = secsqz - 1
    
    #this is the cosine of the solar zenith angle
  
    ## Multiply by the solar constant.  Take the "max" because there is no such thing
    ## as negative solar radiation.
    rshort[itime] <- max(0,solar * cosz)
  }

  return(tansqz)
    
}




## Make the graph.  Note that 0000 UTC corresponds to 8pm Eastern Daylight Time.
# outdir <- 'Documents/Courses/GEO430/Spring2015/Lectures/Unit1/'
# tiff(paste(outdir,'cosz.tif',sep=''))
# par(cex=1.5,lwd=2,mgp=c(2.5,1,0))
# plot(rshort,axes=F,ylab=expression(paste('Insolation [W ',m^-2,']')),xlab='Time [EST]',type='l')
# axis(1,at=c(1,4,7,10,13,16,19,22),labels=c('7PM','10PM','1AM','4AM','7AM','10AM','1PM','4PM'))
# axis(2)
# box()
# dev.off()

