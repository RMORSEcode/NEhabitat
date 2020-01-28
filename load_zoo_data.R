### RM load habitat modeling database from Jim Morley
# There are two dataframes within this file.  The first is "dat" - list of species observations including haul ID, wtcpue, 
# logwtcpue and presfit (presences/absence).  The second is "hauls" which also includes a haulid, lat/lon and a number of environmental variables.
# You can use the haulid to join the tables for the desired species.
## see https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0196127

library(maps)
library(mapdata)
library(sp)
library(maptools)
library(marmap)
library(lubridate)
library(readxl)
library(mgcv)
library(dplyr)

load("C:/Users/ryan.morse/Desktop/1_habitat_analysis_2017/hauls_catch_Dec2017.RData")
table(hauls$region)
NEhauls=hauls[hauls$region %in% c('DFO_Newfoundland', 'DFO_ScotianShelf', 'DFO_SoGulf', 'NEFSC_NEUS', 'VIMS_NEAMAP'),]

length(unique(dat$haulid))

load("C:/Users/ryan.morse/Desktop/1_habitat_analysis_2017/hauls_catch_Dec2017.RData")
load("C:/Users/ryan.morse/Desktop/1_habitat_analysis_2017/survdat/Survdat.RData")
# add to survdat
survdat$latbin=round(survdat$LAT/0.25)*0.25
survdat$lonbin=round(survdat$LON/0.25)*0.25
survdat$month=month(survdat$EST_TOWDATE)
survdat$doy=yday(survdat$EST_TOWDATE)


load('/home/ryan/1_habitat_analysis_2017/habitat_ws_20191009.RData')

### add NE Hauls with length data from SDM group and length at maturity from stock assessment reports
GF=readRDS('C:/Users/ryan.morse/Documents/GitHub/SDM-convergence/data/jude_groundfish_training.rds')
GF=readRDS('/home/ryan/Git/SDM-convergence/data/jude_groundfish_training.rds')
# lmd=read_excel('C:/Users/ryan.morse/Documents/GitHub/NEhabitat/AdultMaturityList.xlsx')
# lmd$spnm=tolower(lmd$Taxon)
# test=strsplit(unique(GF$sppocean), split="_Atl")
# t=data.frame(matrix(unlist(test), nrow=13, byrow=T),stringsAsFactors=FALSE)
# colnames(t)='nm'
# t=left_join(t, lmd, by=c('nm'='spnm'), )
# tt=!duplicated(t$`BTS #`)
# Lm=t[tt,]



### Load EcoMon data - ichthyoplankton for NE Groundfish habitat work ###
# ZPD=read_excel('C:/Users/ryan.morse/Desktop/Iomega Drive Backup 20171012/1 RM/JHare data/EcoMon_Plankton_Data_v3_5.xlsx', sheet='Data' , col_names=T) # newest data through 2015; *** NEW FORMAT
# ZPD=read_excel('/home/ryan/1_habitat_analysis_2017/EcoMon_Plankton_Data_v3_6.xlsx', sheet='Data' , col_names=T) 
ZPD=read_excel('C:/Users/ryan.morse/Desktop/1_habitat_analysis_2017/EcoMon_Plankton_Data_v3_6.xlsx', sheet='Data' , col_names=T)
dt=as_date(ZPD$date)#, origin = "1899-12-30")
ichnms=read.csv('C:/Users/ryan.morse/Desktop/Iomega Drive Backup 20171012/1 RM/JHare data/ichthyonames.csv', header = F, stringsAsFactors = F)
nms=c('nofish_100m3',	'urospp_100m3',	'gadmor_100m3',	'melaeg_100m3',
      'polvir_100m3',	'merbil_100m3',	'sebspp_100m3',	'anaspp_100m3',	'parden_100m3',
      'pseame_100m3',	'glycyn_100m3',	'scoaqu_100m3',	'lopame_100m3')
gnms=ichnms[ichnms$V1 %in% nms,] # subset to NE Groundfish
gnms$V4=gnms$V3
gnms$V4[8]="Anarhichas lupus"
gnms$V4[7]="Sebastes fasciatus"
gnms$V4[2]="Urophycis tenuis" #white hake FIX issue with non classified hake larvae - use for both in model

gnms=data.frame(rbind(as.matrix(gnms), as.matrix(gnms[2,])), stringsAsFactors = F)
gnms$V4[14]="Urophycis chuss" #red hake
gnms$spp=lapply(gnms$V4, FUN=toupper)
gnms$spp=toupper(gnms$V4)
gnms$sp2=tolower(gnms$spp)
# tt=left_join(gnms, lmd, by=c("sp2"="spnm"))
# tt2=!duplicated(tt$`BTS #`)
# Lm2=tt[tt2,]
# Lmf=full_join(Lm, Lm2, by=c("nm"="sp2"))
# write.csv(Lmf, file='lmf.csv', col.names = T)
Lmf=read_excel('C:/Users/ryan.morse/Documents/GitHub/NEhabitat/Lm_included.xlsx') # read in final length at maturity
Lmf=read_excel('/home/ryan/Git/NEhabitat/Lm_included.xlsx')

svspplu=read.csv('C:/Users/ryan.morse/Desktop/1_habitat_analysis_2017/svspp_lookup.csv', stringsAsFactors = F)
svnms=svspplu[svspplu$SCINAME %in% gnms$spp,]


DOY=yday(dt) #day of year
month=as.numeric(format(dt, '%m'))
year=as.numeric(format(dt, '%Y'))
ZPD$year=year
ZPD$month=month
ZPD$dt=dt
ZPD$DOY=DOY
ZPD$day=as.numeric(format(dt, '%d'))
ZPD$lat2=ceiling(ZPD$lat) #use for binning into 1 degree bins for removal of undersampled bins
ZPD$lon2=floor(ZPD$lon) #use for binning into 1 degree bins for removal of undersampled bins

#keep just records with salinity and temp, remove NA from dominant ich and zoo
test=ZPD[complete.cases(ZPD$sfc_temp),]
test=test[complete.cases(test$sfc_salt),]
test=test[complete.cases(test$btm_salt),]
test=test[complete.cases(test$btm_temp),]
test=test[complete.cases(test$melaeg_100m3 ),]
test=test[complete.cases(test$calfin_100m3),]
test=test[complete.cases(test$ctyp_100m3),]
test$latbin=round(test$lat/0.25)*0.25
test$lonbin=round(test$lon/0.25)*0.25
barplot(table(test$year))


#subset just icthyoplankton of interst
ich=test[,which(colnames(test)%in% nms)]
ich$date=test$date
ich$lat=test$lat
ich$lon=test$lon
ich$latbin=round(ich$lat/0.25)*0.25
ich$lonbin=round(ich$lon/0.25)*0.25
vars=test[,c('date', 'lat', 'lon', 'latbin', 'lonbin', 'station', 'depth', 'sfc_temp', 'sfc_salt', 'btm_temp', 'btm_salt', 'month', 'year')]


# subset fraction of data for testing/ traiing
df %>% sample_frac(0.33)
## remove zero presence
# df=ich[which(ich$melaeg_100m3 >0),]
# dfv=vars[which(ich$melaeg_100m3>0),]
# dfz=test[which(ich$melaeg_100m3>0),]
df=ich
dfv=vars
dfz=test

### limit dfz to dominant taxa
X=20 # criteria to use as minimum percent in samples
ZPDa=dfz
ZPDa=ZPDa[!is.na(ZPDa$ich_gear),] # Remove NA in zooplankton rows
p.a=ZPDa[,15:197]
p.a[p.a > 0]=1 # presence/absence
count=colSums(p.a)
pct=(count/dim(ZPDa)[1])*100
crit=which(pct>X)
crit2=crit[31:60]
ZPDa=ZPDa[c(1:14,crit2+14)] # data limited to taxa occurring in > X percent of samples
p.a=p.a[,crit2]
dfz=ZPDa

## try fitting gams
y=log10(df$merbil_100m3+1) 
x=log10(dfz$calfin_100m3+1) #log Calfin
Sample_data <- data.frame(y,x)
g1=gam(y~ s(x), method="REML")
g1=gam(y~ s(x)+s(dfv$month, bs='cc', k=10) +s(dfv$btm_temp), method="REML")

g1=gam(log10(df$melaeg_100m3+1)~ s(log10(dfz$calfin_100m3+1))+s(dfv$month, bs='cc', k=12) +s(dfv$btm_temp) + s(log10(dfv$depth+1)), method="REML")

ggplot(Sample_data, aes(x, y)) + geom_point() + geom_smooth(method = "gam", formula = y ~s(x))

y=log10(df$gadmor_100m3+1) # log fish sp
z1=log10(dfz$calfin_100m3+1) #log Calfin
z2=log10(dfz$pseudo_100m3+1)
z3=log10(dfz$mlucens_100m3+1)
z4=log10(dfz$tlong_100m3+1)
z5=log10(dfz$cham_100m3+1)
z6=log10(dfz$ctyp_100m3+1)
z7=log10(dfz$calminor_100m3+1)
# Sample_data <- data.frame(y,x)
# g1=gam(df$melaeg_100m3 ~ s(dfz$calfin_100m3), method="REML")
g1=gam(y~ s(z1), method="REML")
g1=gam(y~ s(z1)+s(dfv$sfc_temp), method="REML")
g1=gam(y~ s(z1)+s(dfv$sfc_temp)+s(dfv$sfc_salt), method="REML")
g1=gam(y~ s(z7)+s(dfv$month, bs='cc', k=10) +s(dfv$btm_temp), method="REML")

g1=gam(y~ s(z1)+s(z2)+s(z3)+s(dfv$sfc_temp), method="REML")

g1=gam(y~ s(z4)+s(z5)+s(z6)+s(dfv$sfc_temp), method="REML")

g1=gam(y~ +s(z2, k=10)+s(dfv$lon, dfv$lat, k=48)+s(dfv$sfc_temp, k=20), family=nb, method="REML")
g1=gam(y~ +s(z2, k=10)+s(dfv$lon, dfv$lat, k=48)+s(dfv$sfc_temp, k=20), family="poisson", method="REML")

# ggplot(Sample_data, aes(x, y)) + geom_point() + geom_smooth(method = "gam", formula = y ~s(x))


plot(g1)
gam.check(g1) #p values are for the test of the null hypothesis that the basis dimension used is of sufficient size
summary.gam(g1) #test null hypothesis of a zero effect of the indicated spline
# par(mfrow = c(2,2))



## Merge data sets, add plus minus time (4days) and lat/long (0.5 degrees)
library(data.table)
survdat[, Datetime := as.POSIXct( Datetime, format = "%Y-%m-%d %H:%M:%S") ]
survdat[, `:=`(Datetime_max = Datetime + 345600,
               Datetime_min = Datetime - 345600,
               LAT_max = LAT + 0.25,
               LAT_min = LAT - 0.25,
               LON_max = LON + 0.25,
               LON_min = LON - 0.25) ]
dfz$LAT=dfz$lat
dfz$LON=dfz$lon
dfz=as.data.table(dfz)
dfz[, Datetime := as.POSIXct( date, format = "%Y-%m-%d %H:%M:%S") ]
dfz[, `:=`(Datetime_max = Datetime + 345600,
           Datetime_min = Datetime - 345600,
           LAT_max = lat + 0.25,
           LAT_min = lat - 0.25,
           LON_max = lon + 0.25,
           LON_min = lon - 0.25) ]

dt3=survdat[ dfz, on = .( Datetime <= Datetime_max, 
                          Datetime >= Datetime_min, 
                          LAT <= LAT_max, 
                          LAT >= LAT_min, 
                          LON <= LON_max, 
                          LON >= LON_min ),
             nomatch = 0L]

# subset on species...
i=3
t=dt3[which(dt3$svspp==svnms$SVSPP[i]),]



# ZPDb=ZPD[,c(seq(1,14,1), seq(290,297,1), seq(106,197,1))] # check to make sure these are correct against 'nms' if data source changes!!!
ZPDb=ZPD[,c(seq(1,14,1), seq(198,296,1))] # check to make sure these are correct against 'nms' if data source changes!!!
ZPDb=ZPDb[order(ZPDb$date),]
ZPDb=ZPDb[which(ZPDb$year > 1976),] # remove NA data in years prior to 1977


#### Select only taxa present in yearly data > x percent of samples
X=15 # criteria to use as minimum percent in samples
ZPDa=ZPDb
ZPDa=ZPDa[!is.na(ZPDa$ich_gear),] # Remove NA in zooplankton rows
# Reduce to taxa occurrance > x percent in samples
p.a=ZPDa[,15:106]
p.a[is.na(p.a$nofish_10m2)]=0
p.a[is.na(p.a$nofish_100m3)]=0
p.a[p.a > 0]=1 # presence/absence
count=colSums(p.a)
pct=(count/dim(ZPDa)[1])*100
crit=which(pct>X)
# ZPDa=ZPDa[c(1:14,crit+14)] # data limited to taxa occurring in > X percent of samples
nms=c('nofish_100m3',	'urospp_100m3',	'gadmor_100m3',	'melaeg_100m3',
      'polvir_100m3',	'merbil_100m3',	'sebspp_100m3',	'anaspp_100m3',	'parden_100m3',
      'pseame_100m3',	'glycyn_100m3',	'scoaqu_100m3',	'lopame_100m3')


nms=c('nofish_10m2',	'urospp_10m2',	'gadmor_10m2',	'melaeg_10m2',
      'polvir_10m2',	'merbil_10m2',	'sebspp_10m2',	'anaspp_10m2',	'parden_10m2',
      'pseame_10m2',	'glycyn_10m2',	'scoaqu_10m2',	'lopame_10m2')
tt=colnames(ZPDa)[c(1:14, 107:113)]
test=ZPDa[,c(tt,nms)]
p.a=test[,22:34]
# p.a[is.na(p.a$nofish_10m2)]=0
p.a[is.na(p.a$nofish_100m3)]=0
p.a[p.a > 0]=1 # presence/absence
count=colSums(p.a)
pct=(count/dim(ZPDa)[1])*100
barplot(pct)
ich=test
ichpa=ich
ichpa[,22:34]=p.a

pdf(file='EcoMon_PA_ichthyoplankton.pdf')
for (i in 2:length((nms))){
tt=ichpa[ichpa[nms[i]]>0,]
barplot(table(tt$month), main=nms[i])
barplot(table(tt$year), main=nms[i])
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(tt$lon,tt$lat, pch=19)
}
dev.off()

pdf(file='EcoMon_ichthyoplankton.pdf')
for (i in 9:length((nms))){
  tt=ich[complete.cases(ich[nms[i]]),]
  tt=ich[ich[nms[i]]>0,]
  # tt=tt[which((tt$month==8)|(tt$month==9)|(tt$month==10)),]
  barplot(table(tt$month), main=nms[i])
  barplot(table(tt$year), main=nms[i])
  barplot(table(tt[nms[i]]), main=paste(nms[i], ' count', sep=''))
  ttt=as.numeric(unlist((tt[nms[i]])))
  barplot(table(cut(ttt, 20)), main=paste(nms[i]))
  map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
  map.axes(las=1)
  points(tt$lon,tt$lat, pch=19)
}
dev.off()

## survdat
sn="SPRING"
sn="FALL"
i=155
test=survdat[which(survdat$SVSPP==i),]
test=test[which(test$SEASON==sn),]
test=test[complete.cases(test$LENGTH),]
barplot(table(test$LENGTH), main=paste(sn, svnms$COMNAME[which(svnms$SVSPP==i)], sep=" "))
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(test$LON, test$LAT, pch=19)

#subset for maturity
xx=22.3
test2=test[which(test$LENGTH<xx),]
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(test2$LON, test2$LAT, pch=19)
test2=test[which(test$LENGTH>xx),]
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(test2$LON, test2$LAT, pch=19)

# add observer length data, processed with add_lat_long.r
# obl= readRDS('C:/Users/ryan.morse/Desktop/1_habitat_analysis_2017/Obs_data/Obs Data/fishLengths.RDS')
obl=load("/home/ryan/1_habitat_analysis_2017/Observer_data/Obs Data/OBS.Rdata")
obl=load('C:/Users/ryan.morse/Desktop/1_habitat_analysis_2017/Obs_data/Obs Data/OBS.Rdata')
obl <- lapply(x2, function(x) as.numeric(as.character(x)))
obl=data.frame(obl)


### Area weight sample abundance using 1 degree binned latitudes for aggregation
setwd("G:/1 RM/JHare data")
setwd('C:/Users/ryan.morse/Desktop/Iomega Drive Backup 20171012/1 RM/JHare data')
wt=read.csv('NESareaweights.csv')
wt=wt[9:19,]
ZPDa.test=ZPDa
for (ll in wt$lat) {
  ZPDa.test[which(ZPDa.test$lat2 == ll),24:length(ZPDa)]=ZPDa.test[which(ZPDa.test$lat2 == ll),24:length(ZPDa)]/wt$wt[which(wt$lat == ll)]
}
ZPDa=ZPDa.test
rm(ZPDa.test)

#####
# SPLIT DATA INTO SEASONS (or not - 'Yearly') AND DROP YEARS
#####
#_Select season and subset, remove poorly sampled years
SEASON='Spr'
ZPDa=ZPDa[which((ZPDa$month==2)|(ZPDa$month==3)|(ZPDa$month==4)),]
### Manual exclusion:
# excludes=c(1989, 1990, 1991, 1992, 1993, 1994) #based on 1977-2013 NES analysis, later years not subject to rigor yet (July 2016)
# ZPDa=subset(ZPDb[which(!(ZPDb$year %in% excludes)),])

SEASON='Fall'
ZPDa=ZPDa[which((ZPDa$month== 9)|(ZPDa$month==10)|(ZPDa$month==11)),]
### Manual exclusion:
# excludes=c(1989, 1990, 1991, 1992, 1993) # based on 1977-2013 NES analysis, later years not subject to rigor yet (July 2016)
# ZPDa=subset(ZPDb[which(!(ZPDb$year %in% excludes)),])

SEASON='Yearly'
# ZPDa=ZPDa



## Count occurrence of lat lon combinations over time series
#__ all lat lon combinaions:
yy=unique(ZPDa$year)
test7=data.frame(ZPDa[,c('lat2', 'lon2')])
t=table(test7)
library(Matrix)
nnzero(t) # 47 for seasonal, 54 for yearly
tt=t
tt[tt==0]=NA
sum(!is.na(tt))
samples.all=data.frame(tt)

# Dataframe of number of 1-degree binned lat/lon combinations by year
samples.yrs=samples.all
for (yr in 1:length(yy)){
  test8=data.frame(ZPDa[ZPDa$year==yy[yr],c('lat2', 'lon2')])
  t2=table(test8)
  test9=data.frame(t2)
  samples.yrs=merge(samples.yrs, test9, by=c('lat2', 'lon2'), all=T)
}
colnames(samples.yrs)[3]='all'
colnames(samples.yrs)[4:(length(yy)+3)]=yy

###__ Drop from data those 1-degree bins sampled < NNN times out of entire time series (40 yrs as of 2016)
if ((SEASON == 'Spr') | (SEASON =='Fall')){
  NNN=30 # seasonal cut off: bins sampled less than this number over the time series will be removed
} else {
  NNN=35 # yearly cut off: bins sampled less than this number over the time series will be removed
}
samples.yrs.bins=data.frame(rowSums(samples.yrs[,3:(length(yy)+3)] >0, na.rm=T)) # count number of years lat/lon bin was sampled over time series
max(samples.yrs.bins) #40
samples.yrs.bins.dropped=samples.yrs[which(samples.yrs.bins < NNN),] # all lat/lon bins to be dropped
ZPDa$merge=paste(ZPDa$lat2, ZPDa$lon2, sep='_') #ugly solution but works - combine lat,lon into comparable vector for removal key
to.drop=paste(samples.yrs.bins.dropped$lat2, samples.yrs.bins.dropped$lon2, sep='_') # just coords of lat/lon bins to be dropped
ZPDa=subset(ZPDa[which(!ZPDa$merge %in% to.drop),])

###__ Drop from data years where < NNN bins were sampled in a given year -> max: 47 total bins for NES, drop years with < NNN bins
if ((SEASON == 'Spr') | (SEASON =='Fall')){
  NNN=30 # set criteria to drop years where < NNN bins sampled in given year
} else {
  NNN=35 # set criteria to drop years where < NNN bins sampled in given year
}

samples.yrs[samples.yrs==0]=NA # change zero to NA
na_count <-data.frame(sapply(samples.yrs[,3:42], function(y) sum(length(which(!is.na(y)))))) # count of sampled bins by year
table(na_count)
excludes=as.numeric(rownames(na_count)[na_count < NNN])
ZPDa=subset(ZPDa[which(!(ZPDa$year %in% excludes)),])

#Take median date from each cruise and assign cruise to bimonth
cruises=unique(ZPDa$cruise_name)
for (i in 1:length(cruises)){
  ZPDa$medmonth[ZPDa$cruise_name == cruises[i]]=median(ZPDa$DOY[ZPDa$cruise_name == cruises[i]])
}
# For bi-monthly means aggregation
# ZPD$bm=NA
ZPDa$bmm=NA
ZPDa$bmm[which(as.integer(ZPDa$medmonth) %in% seq(0,59))]=1
ZPDa$bmm[which(as.integer(ZPDa$medmonth) %in% seq(60,120))]=3
ZPDa$bmm[which(as.integer(ZPDa$medmonth) %in% seq(121,181))]=5
ZPDa$bmm[which(as.integer(ZPDa$medmonth) %in% seq(182,243))]=7
ZPDa$bmm[which(as.integer(ZPDa$medmonth) %in% seq(244,304))]=9
ZPDa$bmm[which(as.integer(ZPDa$medmonth) %in% seq(305,366))]=11



table(ZPDa$month)


