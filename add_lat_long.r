library(sp)
setwd('C:/Users/ryan.morse/Desktop/1_habitat_analysis_2017/Obs_data/Obs Data')
setwd('/home/ryan/1_habitat_analysis_2017/Observer_data/Obs Data')
      
x=readRDS('fishLengths.RDS')
g=readRDS('gearTable.RDS')
gt=g$data

## remove records without lat long
x2=x[complete.cases(x$LATHBEG),]
x2=x2[complete.cases(x2$LONHBEG),]
## testing format dms to decimal degrees
z=x2$LATHBEG
dd=substr(z, 1, 2)
mm=substr(z, 3, 4)
ss=substr(z, 5, 6)
a=paste(dd, 'd', mm,'\'', ss, '" ', 'N', sep='')
Lat=as.numeric(sp::char2dms(a))
x2$LAT=Lat
z=x2$LONHBEG
dd=substr(z, 1, 2)
mm=substr(z, 3, 4)
ss=substr(z, 5, 6)
a=paste(dd, 'd', mm,'\'', ss, '" ', 'W', sep='')
Lon=as.numeric(sp::char2dms(a))
x2$LONG=Lon

# bin to quarter or half degree breaks
# br=0.25
# lnbr=seq(from=min(floor(x2$LONG)), by=br, to=max(ceiling(x2$LONG)))
# ltbr=seq(from=min(floor(x2$LAT)), by=br, to=max(ceiling(x2$LAT)))
x2$binlat=round(x2$LAT/0.25)*0.25
x2$binlon=round(x2$LONG/0.25)*0.25
# library(data.table)
# survdat[, Datetime := as.POSIXct( Datetime, format = "%Y-%m-%d %H:%M:%S") ]
# survdat[, `:=`(Datetime_max = Datetime + 345600,
#                Datetime_min = Datetime - 345600,
#                LAT_max = LAT + 0.25,
#                LAT_min = LAT - 0.25,
#                LON_max = LON + 0.25,
#                LON_min = LON - 0.25) ]
#     
setwd('/home/ryan/1_habitat_analysis_2017/Observer_data/Obs Data')
save(x2, file='OBS.Rdata')
