# Kevin's code to calculate center of biomass and movement for trawl survey data
## RM using this to clean survdat to use in habitat models and merge with zooplankotn/ichtyoplankton data
## must load zooplankton data prior to running this with {load_zoo_data.R}
## RM updated 2021 to include haddock data corrected for catchability for the Bigelow time series from Sean Lucey (fixed Feb 2021)


require(raster)
# require(ncdf)
library(sp)
library(maptools)
library(marmap)
library(geosphere)
library(lubridate)
library(dplyr)
library(tidyr)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

## Subset survdat to only groundfish species with Length at maturity data
# load("C:/Users/ryan.morse/Downloads/SurvdatBio (3).RData") #loads survdat.bio
load("C:/Users/ryan.morse/Downloads/survdat_lw.RData") # loads survdat.lw ; from Sean Lucey data into 2019
load("/home/ryan/Downloads/survdat_lw.RData")# loads survdat.lw ; from Sean Lucey data into 2019
newsurvdat=loadRData('/home/ryan/Downloads/NEFSC_BTS_2021.RData') # new data from Sean 2021, possible updates to haddock due to code q correction
test=newsurvdat$survdat$SVSPP %in% Lmf$SVSPP
nsdat=newsurvdat$survdat[test,]

load('/home/ryan/Downloads/NEFSC_BTS_2021.RData')
survdat <- survey$survdat %>% filter(YEAR!=2020)
test=survdat$SVSPP %in% Lmf$SVSPP
survdat=survdat[test,]
nsurvday=survdat

## 20210326 conversion factor issue from Sean Lucey - checking on differences between survdats
nsdat2=nsdat %>% filter(SVSPP==74)
sdat=survdat.lw %>% filter(SVSPP==74)
colnames(sdat)
colnames(nsdat)
# test1=nsdat2 %>% filter(YEAR==2010)
# test2=sdat %>% filter(YEAR==2010)
# test3=test1$BIOMASS==test2$BIOMASS
# test3=(test1$BIOMASS < 1.025*test2$BIOMASS) & (test1$BIOMASS > 0.975*test2$BIOMASS)
# sum(test3)/length(test3)
yyyy=unique(nsdat2$YEAR)
csdat=matrix(NA, ncol=5, nrow=length(unique(nsdat2$YEAR)))
colnames(csdat)=c('Year', 'Abundance', 'Biomass', 'Length', 'NumLen')
for (i in 1:length(unique(nsdat2$YEAR))){
  test1=nsdat2 %>% filter(YEAR==yyyy[i])
  test2=sdat %>% filter(YEAR==yyyy[i])
  csdat[i,1]=yyyy[i]
  test3=test1$ABUNDANCE==test2$ABUNDANCE
  csdat[i,2]=sum(test3)/length(test3)
  test3=test1$BIOMASS==test2$BIOMASS
  csdat[i,3]=sum(test3)/length(test3)
  test3=test1$LENGTH==test2$LENGTH
  csdat[i,4]=sum(test3)/length(test3)
  test3=test1$NUMLEN==test2$NUMLEN
  csdat[i,5]=sum(test3)/length(test3)
}
# test3=test1$ABUNDANCE==test2$ABUNDANCE
# test3=(test1$ABUNDANCE < 1.025*test2$ABUNDANCE) & (test1$ABUNDANCE > 0.975*test2$ABUNDANCE)
# sum(test3)/length(test3)
# test3=test1[test1$ABUNDANCE!=test2$ABUNDANCE,]
# test3$newABUNDANCE=test2$ABUNDANCE[test1$ABUNDANCE!=test2$ABUNDANCE]




#### 20210326 IMPORTANT NOTE:
#### YES THERE IS A LARGE DIFFERENCE DUE TO WRONG CONVERVERSION FACTORS APPLIED FOR HADDOCK (074), WINDOWPANE (108), AND BUTTERFISH (131)


test=survdat$SVSPP %in% Lmf$SVSPP
survdat=survdat[test,]
sort(unique(survdat$SVSPP)) #verify
length(unique(survdat$SVSPP))

## choose surdat and remove older non-length version to avoid confusion
rm(survdat)
survdat=survdat.lw
survdat=survdat.bio
survdat=nsdat

## newest survdat missing INDWT, WGTLEN, SIZECAT from survdat.lw, remove 2020 and add back
test=survdat %>% filter(YEAR!=2020)
test2=test %>% filter(SVSPP==074, YEAR==2010)
test3=survdat.lw %>% filter(SVSPP==074, YEAR==2010)
## they are the same
test$INDWT=survdat.lw$INDWT
test$WGTLEN=survdat.lw$WGTLEN
test$SIZECAT=survdat.lw$SIZECAT
survdat=test

### loop to split taxa into juvenile or adult using Lm data
# survdat$stg=NA
# for (j in 1:length(survdat$stg)){
#   survdat$stg[j]=ifelse(survdat$LENGTH[j]<Lmf$Lm[which(Lmf$SVSPP==survdat$SVSPP[j])], "juv", "adt")
# }
### Better to do this with left_join
test=data.frame(survdat[,c("SVSPP", "LENGTH")])
colnames(test)=c("SVSPP", "LENGTH")
test=left_join(test, Lmf[,c("SVSPP", "Lm")], by="SVSPP")
test$stg=ifelse(test$LENGTH<test$Lm, "Juv", "Adt")
survdat$stg=test$stg
rm(test)

### count juv and adult to see if numbers match ABUNDANCE 
survdat=survdat %>% group_by_at(vars(SVSPP,CRUISE6,STATION,STRATUM,ABUNDANCE,stg)) %>% mutate(stgsum=sum(NUMLEN), stgwtsum=sum(WGTLEN))
### now addd up weights for all fish measured, check against biomass, total abundance
survdat=survdat %>% group_by_at(vars(SVSPP,CRUISE6,STATION,STRATUM,ABUNDANCE)) %>% mutate(totlen=sum(NUMLEN), totwgt=sum(WGTLEN), abndiff=ABUNDANCE-totlen, biodiff=BIOMASS-totwgt)
### now add corrected biomass and abundance to account for any missing staged values (eg. for large hauls)
survdat=survdat %>% group_by_at(vars(SVSPP,CRUISE6,STATION,STRATUM,ABUNDANCE, stg)) %>%
  mutate(stgpctabn=stgsum/totlen, stgpctwgt=stgwtsum/totwgt, corABN=round(stgpctabn*ABUNDANCE,digits = 0), corBIO=round(stgpctwgt*BIOMASS, digits=2))

# 20210406 ONLY for newest survdat - missing wgtlen, only correct abundance
survdat=survdat %>% group_by_at(vars(SVSPP,CRUISE6,STATION,STRATUM,ABUNDANCE,stg)) %>% mutate(stgsum=sum(NUMLEN)) #, stgwtsum=sum(WGTLEN))
survdat=survdat %>% group_by_at(vars(SVSPP,CRUISE6,STATION,STRATUM,ABUNDANCE)) %>% mutate(totlen=sum(NUMLEN), abndiff=ABUNDANCE-totlen)
survdat=survdat %>% group_by_at(vars(SVSPP,CRUISE6,STATION,STRATUM,ABUNDANCE, stg)) %>%
  mutate(stgpctabn=stgsum/totlen, corABN=round(stgpctabn*ABUNDANCE,digits = 0))


# barplot(table(round(survdat$biodiff[(abs(survdat$biodiff>12))],1)))
y=survdat$SVSPP
x=abs(survdat$biodiff)
y=y[x>120]
x=x[x>120]
barplot(table(round(y))) #which species have the highest mismatches in biomass
barplot(table(round(x))) #how many

### now create dataframe with unique tows only, separated by stage
test=survdat[,c("CRUISE6","STATION","STRATUM","YEAR", "SVSPP", "stg")] # needs SVSPP...
svdtunq=survdat[!duplicated(test),]

svdtunq=survdat ## keep all stations for now, remove later

### calc percent of biomass by stage(adt, juv) per unique tow to be applied to BIOMASS for corrected values
# svdtunq$stgwgtpct=round(svdtunq$stgwtsum/svdtunq$totwgt, 2)
# svdtunq$corBIOMASS=svdtunq$stgwgtpct*svdtunq$BIOMASS

## not working yet
# svdunqadt=svdtunq %>% filter(stg=="Adt") %>% 
#   select(-(LENGTH:SIZECAT)) %>% 
#   select(-(CATCHSEX)) %>% 
#   spread(SVSPP, ABUNDANCE:corBIO)

# barplot(table(svdtunq$stgwgtpct))
# barplot(table(svdtunq$stgwgtpct[svdtunq$SEASON=='SPRING']))
# barplot(table(svdtunq$stgwgtpct[svdtunq$SEASON=='FALL']))
# barplot(table(svdtunq$stgwgtpct[svdtunq$stgwgtpct<0.99 & svdtunq$SEASON=='FALL']),ylim=c(0,1000), main='FALL')
# barplot(table(svdtunq$stgwgtpct[svdtunq$stgwgtpct<0.99 & svdtunq$SEASON=='SPRING']),ylim=c(0,1000), main='SPRING')
# 
# svspp=73
# svdtunq$stgwgtpct[is.na(svdtunq$stgwgtpct)]=0 #[svdtunq$stgwgtpct<0.99 & svdtunq$SEASON=='FALL' & svdtunq$SVSPP==svspp]),]
# 
# barplot(table(svdtunq$stgwgtpct[svdtunq$stgwgtpct<0.99 & svdtunq$SEASON=='FALL' & svdtunq$SVSPP==svspp]),ylim=c(0,100), main=paste('FALL ', svspp))
# barplot(table(svdtunq$stgwgtpct[svdtunq$stgwgtpct<0.99 & svdtunq$SEASON=='SPRING' & svdtunq$SVSPP==svspp]),ylim=c(0,100), main=paste('SPRING ', svspp))
# 
# 
# library(modes)
# # restricting to <0.99 was masking adult population. Can check on change over years...
# bimodality_amplitude(svdtunq$stgwgtpct[svdtunq$SEASON=='SPRING' & svdtunq$SVSPP==svspp], fig=T)
# bimodality_amplitude(svdtunq$stgwgtpct[svdtunq$SEASON=='FALL' & svdtunq$SVSPP==svspp], fig=T)
# sum(is.na(svdtunq$stgwgtpct[svdtunq$stgwgtpct<0.99 & svdtunq$SEASON=='FALL' & svdtunq$SVSPP==svspp]))

### change survdat long to wide, select biomass based (corBIO) or abundance based (corABN) ###
svdtunq=ungroup(svdtunq) #not needed
svdwide.bio=svdtunq %>% dplyr::select(SVSPP, stg, CRUISE6, STATION, STRATUM,SVVESSEL,YEAR,EST_TOWDATE, SEASON,LAT,LON,DEPTH,
                           SURFTEMP,BOTTEMP,SURFSALIN,BOTSALIN,corBIO) %>%
  pivot_wider(names_from=c(SVSPP, stg), values_from = corBIO, values_fill = list(corBIO=0))

svdwide.abn=svdtunq %>% dplyr::select(SVSPP, stg, CRUISE6, STATION, STRATUM,SVVESSEL,YEAR,EST_TOWDATE, SEASON,LAT,LON,DEPTH,
                               SURFTEMP,BOTTEMP,SURFSALIN,BOTSALIN,corABN) %>%
  pivot_wider(names_from=c(SVSPP, stg), values_from = corABN, values_fill = list(corABN=0), values_fn = {mean})


## 20210406 not using corrected biomass, take mean of multiple stations
svdwide.bio=svdtunq %>% dplyr::select(SVSPP, stg, CRUISE6, STATION, STRATUM,SVVESSEL,YEAR,EST_TOWDATE, SEASON,LAT,LON,DEPTH,
                                      SURFTEMP,BOTTEMP,SURFSALIN,BOTSALIN,BIOMASS) %>%
  pivot_wider(names_from=c(SVSPP, stg), values_from = BIOMASS, values_fill = list(BIOMASS=0), values_fn = {mean})

svdwide.abn=svdtunq %>% dplyr::select(SVSPP, stg, CRUISE6, STATION, STRATUM,SVVESSEL,YEAR,EST_TOWDATE, SEASON,LAT,LON,DEPTH,
                                      SURFTEMP,BOTTEMP,SURFSALIN,BOTSALIN,corABN) %>%
  pivot_wider(names_from=c(SVSPP, stg), values_from = corABN, values_fill = list(corABN=0), values_fn = {mean})

#remove column names with 'x_NA'
svdwide.bio=svdwide.bio%>% dplyr::select(-`72_NA`, -`73_NA`, -`74_NA`,-`75_NA`,-`76_NA`,-`77_NA`,-`101_NA`,-`102_NA`,-`103_NA`,-`105_NA`,
                          -`106_NA`,-`107_NA`,-`108_NA`,-`155_NA`,-`193_NA`,-`197_NA`)
svdwide.abn=svdwide.abn%>% dplyr::select(-`72_NA`, -`73_NA`, -`74_NA`,-`75_NA`,-`76_NA`,-`77_NA`,-`101_NA`,-`102_NA`,-`103_NA`,-`105_NA`,
                              -`106_NA`,-`107_NA`,-`108_NA`,-`155_NA`,-`193_NA`,-`197_NA`)

library(lubridate)
xdt=today()
xdt=gsub('-','',xdt)
save(svdwide.bio, file=paste('/home/ryan/Git/NEhabitat/', xdt, '_unique_wide_format_corBIO_by_stage.Rda', sep=''))
save(svdwide.abn, file=paste('/home/ryan/Git/NEhabitat/', xdt, '_unique_wide_format_corABN_by_stage.Rda', sep=''))

### create date and binned lat lon for mathching to other surveys
# svdate=data.frame(month(survdat$EST_TOWDATE))
# colnames(svdate)[1]='M'
# svdate$Y=year(survdat$EST_TOWDATE)
# svdate$D=day(survdat$EST_TOWDATE)
# # svdate$lat=survdat$LAT
# # svdate$lon=survdat$LON
# svdate$lonbin=round(survdat$LON/0.5)*0.5 ## half degree
# svdate$latbin=round(survdat$LAT/0.5)*0.5
# svdate$doy=as.numeric(strftime(survdat$EST_TOWDATE, format = "%j"))
# svdate$sdoy=svdate$doy-8
# svdate$ldoy=svdate$doy+8
# svdate$index=seq(from=1, to=length(svdate$M), by=1)

### use wide, unique vals to merge with zoo and ich
svdate=data.frame(month(svdwide.bio$EST_TOWDATE))
colnames(svdate)[1]='M'
svdate$Y=year(svdwide.bio$EST_TOWDATE)
svdate$D=day(svdwide.bio$EST_TOWDATE)
# svdate$lat=survdat$LAT
# svdate$lon=survdat$LON
svdate$lonbin=round(svdwide.bio$LON/0.5)*0.5 ## half degree
svdate$latbin=round(svdwide.bio$LAT/0.5)*0.5
svdate$doy=as.numeric(strftime(svdwide.bio$EST_TOWDATE, format = "%j"))
svdate$sdoy=svdate$doy-8
svdate$ldoy=svdate$doy+8
svdate$index=seq(from=1, to=length(svdate$M), by=1)


## get index of unique stations from survdat, apply to long format svdate with index
# retvars <- c("CRUISE6","STATION","STRATUM","YEAR", "SVSPP")
# sdq=survdat[,retvars]
# sdq$index=seq(from=1, to=length(survdat$YEAR), by=1)
# survdat_stations <- survdat[retvars] 
# sdq2=sdq[!duplicated(survdat_stations),] # just keep retvars
# svdate2=svdate[sdq2$index,]

#clean up svdate to match zooplankton and Chl
# svdate1=svdate[which(svdate$Y>1976 & svdate$Y<1982),]
# svdate2=svdate[which(svdate$Y>1981 & svdate$Y<1987),]
# svdate3=svdate[which(svdate$Y>1986 & svdate$Y<1997),]
# svdate4=svdate[which(svdate$Y>1996 & svdate$Y<2006),]
# svdate5=svdate[which(svdate$Y>2005 & svdate$Y<2012),]
# svdate6=svdate[which(svdate$Y>2011 & svdate$Y<2020),]

# library(geosphere)
# distm(c(lon1, lat1), c(lon2, lat2), fun = distHaversine)

### now for zoo and ich data
dfzdate=data.frame(year(dfz$date))
colnames(dfzdate)[1]='zY'
dfzdate$zM=month(dfz$date)
dfzdate$zD=day(dfz$date)
dfzdate$zdoy=as.numeric(strftime(dfz$date, format = "%j"))
dfzdate$zlonbin=round(dfz$lon/.5)*0.5
dfzdate$zlatbin=round(dfz$lat/.5)*0.5
dfzdate$zindex=seq(from=1, to=length(dfz$date), by=1) #add index of original order just to keep track before sorting on time
dfzdate=dfzdate[order(dfzdate$zY, dfzdate$zdoy),]
dfzdate$zsdoy=dfzdate$zdoy
dfzdate$zldoy=dfzdate$zdoy

# dfzdate1=dfzdate[which(dfzdate$zY>1976 & dfzdate$zY<1982),]
# dfzdate2=dfzdate[which(dfzdate$zY>1981 & dfzdate$zY<1987),]
# dfzdate3=dfzdate[which(dfzdate$zY>1986 & dfzdate$zY<1997),]
# dfzdate4=dfzdate[which(dfzdate$zY>1996 & dfzdate$zY<2006),]
# dfzdate5=dfzdate[which(dfzdate$zY>2005 & dfzdate$zY<2012),]
# dfzdate6=dfzdate[which(dfzdate$zY>2011 & dfzdate$zY<2020),]

# try merge and filter (did this in 2 parts, before and after 1998, took long time...)
# dfzdate2=dfzdate[which(dfzdate$zY<1997),]
# dfzdate3=dfzdate[which(dfzdate$zY>1996),]

# ttx=merge(svdate2, dfzdate2, all=T)
colnames(dfzdate2)
colnames(svdate2)
# ttx=left_join(svdate2, dfzdate2, by=c("lonbin"="zlonbin", "latbin"="zlatbin", "Y"="zY", "sdoy"="zdoy", "ldoy"="zdoy"))
### this works, just takes a long time
library(fuzzyjoin)
xdt=today()
xdt=gsub('-','',xdt)
setwd('/home/ryan/Git/NEhabitat/')
# tt=fuzzy_left_join(svdate1, dfzdate1, by=c("lonbin"="zlonbin", "latbin"="zlatbin", "Y"="zY", "sdoy"="zdoy", "ldoy"="zdoy"),
#                     match_fun=list(`==`,`==`,`==`,`<=`,`>=`))
# save(tt, file="survdat_zoo_merge_1977_1981.Rda")
# 
# tt1=fuzzy_left_join(svdate2, dfzdate2, by=c("lonbin"="zlonbin", "latbin"="zlatbin", "Y"="zY", "sdoy"="zdoy", "ldoy"="zdoy"),
#                    match_fun=list(`==`,`==`,`==`,`<=`,`>=`))
# save(tt1, file="survdat_zoo_merge_1982_1986.Rda")
# 
# tt2=fuzzy_left_join(svdate3, dfzdate3, by=c("lonbin"="zlonbin", "latbin"="zlatbin", "Y"="zY", "sdoy"="zdoy", "ldoy"="zdoy"), 
#                    match_fun=list(`==`,`==`,`==`,`<=`,`>=`))
# save(tt2, file="survdat_zoo_merge_1987_1996.Rda")
# 
# tt3=fuzzy_left_join(svdate4, dfzdate4, by=c("lonbin"="zlonbin", "latbin"="zlatbin", "Y"="zY", "sdoy"="zdoy", "ldoy"="zdoy"), 
#                    match_fun=list(`==`,`==`,`==`,`<=`,`>=`))
# save(tt3, file="survdat_zoo_merge_1997_2005.Rda")
# 
# tt4=fuzzy_left_join(svdate5, dfzdate5, by=c("lonbin"="zlonbin", "latbin"="zlatbin", "Y"="zY", "sdoy"="zdoy", "ldoy"="zdoy"), 
#                     match_fun=list(`==`,`==`,`==`,`<=`,`>=`))
# save(tt4, file="survdat_zoo_merge_2006_2011.Rda")
# 
# tt5=fuzzy_left_join(svdate6, dfzdate6, by=c("lonbin"="zlonbin", "latbin"="zlatbin", "Y"="zY", "sdoy"="zdoy", "ldoy"="zdoy"), 
#                     match_fun=list(`==`,`==`,`==`,`<=`,`>=`))
# save(tt5, file="survdat_zoo_merge_2012_2019.Rda")

for(i in 1:length(yrs)){
  selyear=yrs[i]
  svdatex=svdate %>% filter(Y==selyear)
  dfzdatex=dfzdate %>% filter(zY==selyear)
  tt=fuzzy_left_join(svdatex, dfzdatex, by=c("lonbin"="zlonbin", "latbin"="zlatbin", "Y"="zY", "sdoy"="zdoy", "ldoy"="zdoy"),
                     match_fun=list(`==`,`==`,`==`,`<=`,`>=`))
  save(tt, file=paste(xdt,"_survdat_zoo_merge_", selyear,"_.Rda", sep=""))
}

wd=getwd()
mergelist=list.files(wd, pattern = glob2rx("*20210414_survdat_zoo_merge_*")) # 
tt=loadRData(mergelist[[1]])
for(i in 2:length(yrs)){
  tt1=loadRData(mergelist[[i]])
  tt=rbind(tt, tt1)
}
save(tt, file=paste(xdt, '_merged_survdat_zoo_1977_2019.Rda', sep=''))
dmrg4=tt# [unique(tt$index),]

# tt2=tt[complete.cases(tt$zY),]
# tt2$index=seq(from=13952, to=(13951+length(tt2$M)), by=1)
# merge 2 dataframes together
# dmrg=rbind(tt, tt1)
# dmrg1=rbind(dmrg, tt2)
# dmrg2=rbind(dmrg1, tt3)
# dmrg3=rbind(dmrg2, tt4)
# dmrg4=rbind(dmrg3, tt5)

# save(dmrg4, file='updated_final_merged_survdat_zoo_1977_2018.Rda')

## now subset original dataframes before joining together
fish1.bio=svdwide.bio[dmrg4$index,]
fish1.abn=svdwide.abn[dmrg4$index,]
# zoo1=dfz2[dmrg4$zindex,]
zoo1=dfz[dmrg4$zindex,]
ich1=ich[dmrg4$zindex,]

# Create index to merge on
fish1.bio$mrgidx=seq(from=1, to=length(fish1.bio$YEAR), by=1)
fish1.abn$mrgidx=seq(from=1, to=length(fish1.abn$YEAR), by=1)
zoo1$mrgidx=fish1.bio$mrgidx
ich1$mrgidx=fish1.bio$mrgidx

## need to fix this, order is wrong from Lmf
colnames(ich1)[1:13]
colnames(ich1)[1:12]=c("77_ich","73_ich","74_ich","75_ich","72_ich","155_ich","192_ich","103_ich","106_ich","107_ich","108_ich","197_ich")
colnames(ich1)[1:13]
# colnames(ich1)[1:12]=svnms2$ichnm

### now extract static variables from rasters
# rug=loadRData('/home/ryan/Git/NEhabitat/rasters/scaledrugosity.RData') # will not load???
# rug=loadRData("/home/ryan/1_habitat_analysis_2017/static_vars/rast_rugosity.rdata")
### load bottom temp file to use as template for non conforming rasters (remove before load dynamic files)
bt=loadRData('/home/ryan/Git/NEhabitat/rasters/test/RAST_NESREG_1977.04.03.BT.TEMP.YEAR.000066596.RData') #BT
# bt=masked.raster
## depth raster
gz=loadRData('/home/ryan/Git/NEhabitat/rasters/test/depth.RData') #gdepth
# gz=resample(gdepth, bt, 'bilinear')
gz2=(gz*-1)
# gd2=crop(gz2, bt)
# gd2=mask(gd2, bt)
gd3=gz2
gd3[gd3>375]=NA # set values > 375 m to NA
gd4=resample(gd3, bt, 'bilinear')
gd4=crop(gd4, bt)
gd4=mask(gd4, bt)
## grain size raster
phi2mm=loadRData('/home/ryan/Git/NEhabitat/rasters/test/grainsizeMM.RData') #phi2mm
phi2=resample(phi2mm, bt)
phi2=crop(phi2, bt)
phi2=mask(phi2, bt)
## sand_pct raster
sandpct=loadRData('/home/ryan/1_habitat_analysis_2017/static_vars/rast_sand_fraction.rdata') #sand_pct
sand2=resample(sandpct, bt)
sand2=crop(sand2, bt)
sand2=mask(sand2, bt)
## mud_pct raster
mudpct=loadRData('/home/ryan/1_habitat_analysis_2017/static_vars/rast_mud_fraction.rdata') #sand_pct
mud2=resample(mudpct, bt)
mud2=crop(mud2, bt)
mud2=mask(mud2, bt)
## rugosity raster
# load('/home/ryan/Git/NEhabitat/rasters/test/scaledrugosity.RData') #rugscl
rug=loadRData('/home/ryan/Git/NEhabitat/rast_rugosity.rdata')
mmin=cellStats(rug, 'min')
mmax=cellStats(rug, 'max')
rugscl=calc(rug, fun=function(x){(x-mmin)/(mmax-mmin)}) # rescale from -4:2 -> 0:1
rug2=resample(rugscl, bt)
ex=extent(bt)
rug2=crop(rug2, ex)
rug2=mask(rug2, bt)
# load Chlorophyll climatology rasters, resample, rename
chllist=list.files('/home/ryan/1_habitat_analysis_2017/chl/NESREG/clim/', pattern='.RData')
for (i in (c(2,4,10))){
  load(paste('/home/ryan/1_habitat_analysis_2017/chl/NESREG/clim/', chllist[i], sep=''))
  masked.raster=crop(masked.raster, bt)
  masked.raster=resample(masked.raster, bt)
  masked.raster=mask(masked.raster, bt)
  newname=paste("chl.", i, sep='')
  assign(newname, masked.raster)
  rm(masked.raster)
}

# coords=sp::SpatialPoints(FData.abn[,c(10,9)])
coords=FData.abn[,c(10,9)]
extrug=raster::extract(rug2, coords)
extz=raster::extract(gd4, coords)
extphi=raster::extract(phi2, coords)
extsand=raster::extract(sand2, coords)
extmud=raster::extract(mud2, coords)
extchl2=raster::extract(chl.2, coords)
extchl4=raster::extract(chl.4, coords)
extchl10=raster::extract(chl.10, coords)

## choose 1 - Either biomass or abundance as unit of measure from survdat to use in GAMs ##
FData=merge(fish1.bio, ich1[,c(1:12,19)], by="mrgidx"); fishtype='biomass'
FData=merge(FData, zoo1[,c(15:27)], by="mrgidx")
# OR #
FData=merge(fish1.abn, ich1[,c(1:12,19)], by="mrgidx"); fishtype='abundance'
FData=merge(FData, zoo1[,c(15:27)], by="mrgidx")


### set up year list to match files with
yrlist=seq(from=1977, to=2019, by=1)
FDss=FData %>% dplyr::select(LON, LAT, YEAR, SEASON)
FData$BT=NA
FData$ST=NA
FData$cty=NA
FDss$BT=NA
FDss$ST=NA
FDss$cty=NA
SEASONss=unique(FData$SEASON)
for (jj in 1:2){
  season=SEASONss[jj]
  if (season=="SPRING"){
    SEASON='Spr'}
  if (season=="FALL"){
    SEASON='Fall'
  }
  ## list data files in each folder
  btlist=list.files(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/BT2', sep=''), pattern = 'RAST_NESREG_')
  stlist=list.files(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/ST2', sep=''), pattern = 'RAST_NESREG_')
  # zlist=list.files(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/pseudo', sep=''))
  zlist=list.files(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/ctyp', sep=''))
  ## parse year from filenames e.g #"RAST_NESREG_1977.04.03.BT.TEMP.YEAR.000066596.RData"
  tb=strsplit(btlist, split=('RAST_NESREG_'))
  ttb=sapply(tb, function(x) strsplit(x, "[.]")[[2]][1], USE.NAMES=FALSE)
  ttb2=as.numeric(ttb)
  #Surface temp
  ts=strsplit(stlist, split=('RAST_NESREG_'))
  tts=sapply(ts, function(x) strsplit(x, "[.]")[[2]][1], USE.NAMES=FALSE)
  tts2=as.numeric(tts)
  # Zooplankton (pseudocal)
  # tz=strsplit(zlist, split=('RAST_NESREG_')) #old one for PSE
  tz=strsplit(zlist, split=('RAST_ctypZZ_')) # for new 7-year series
  ttz=sapply(tz, function(x) strsplit(x, "[.]")[[2]][1], USE.NAMES=FALSE)
  ttz2=as.numeric(ttz)
  ### NOW loop over files, load yearly dynamic raster files and extract
  # FDss$BT=NA
  for (i in 1:length(yrlist)){
    bi=which(yrlist[i]==ttb2) # index of year
    bt=loadRData(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/BT2/', btlist[[bi]], sep=''))
    bi=which(yrlist[i]==tts2) # index of year
    st=loadRData(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/ST2/', stlist[[bi]], sep=''))
    bi=which(yrlist[i]==ttz2) # index of year
    # /home/ryan/1_habitat_analysis_2017/new 7 yr zoo/spring_rasters_7yr/RAST_ctypZZ_2008.04.01.06.TEMP.YEAR.000066596.RData
    # pse=loadRData(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/pseudo/', zlist[[bi]], sep=''))
    cty=loadRData(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/ctyp/', zlist[[bi]], sep=''))
    FDss$BT[which(FDss$YEAR==yrlist[i])]=extract(bt, FDss[which(FDss$YEAR==yrlist[i]),c(1:2)])
    FDss$ST[which(FDss$YEAR==yrlist[i])]=extract(st, FDss[which(FDss$YEAR==yrlist[i]),c(1:2)])
    FDss$cty[which(FDss$YEAR==yrlist[i])]=extract(cty, FDss[which(FDss$YEAR==yrlist[i]),c(1:2)])
  }
  FData$BT[which(FData$SEASON==season)]=FDss$BT[which(FDss$SEASON==season)]
  FData$ST[which(FData$SEASON==season)]=FDss$ST[which(FDss$SEASON==season)]
  FData$cty[which(FData$SEASON==season)]=FDss$cty[which(FDss$SEASON==season)]
}
#### Verify and replace
# plot(FData$BOTTEMP~FData$BT, type='p')
# plot(FData$BOTTEMP[which(FData$SEASON=="SPRING")]~FData$BT[which(FData$SEASON=="SPRING")], type='p')
# plot(log10(FData$ctyp_100m3[which(FData$SEASON=="SPRING")]+1)~FData$cty[which(FData$SEASON=="SPRING")], type='p')
FData$cty3=10^FData$cty # convert from log scale back to N per 100 m3
# plot(FData$BOTTEMP[which(FData$SEASON=="FALL")]~FData$BT[which(FData$SEASON=="FALL")], type='p')
FData$BOTTEMP[which(is.na(FData$BOTTEMP))]=FData$BT[which(is.na(FData$BOTTEMP))]
FData$SURFTEMP[which(is.na(FData$SURFTEMP))]=FData$ST[which(is.na(FData$SURFTEMP))]
FData$ctyp_100m3[which(is.na(FData$ctyp_100m3))]=FData$cty3[which(is.na(FData$ctyp_100m3))]
### Now add extracted static variables
FData$DEPTH2=extz
FData$rug=extrug
FData$grnszmm=extphi
FData$sand_pct=extsand
FData$mud_pct=extmud
FData$chl2=extchl2
FData$chl4=extchl4
FData$chl10=extchl10
### Finally, save output for either abundance or biomass
if (fishtype=='biomass'){
  FData.bio=FData
save(FData.abn, file=paste(xdt, "_Final_merged_fish_corABN_Zoo_Ich.Rda", sep=""))}
if (fishtype=='abundance'){
  FData.abn=FData
save(FData.bio, file=paste(xdt, "_Final_merged_fish_corBIO_Zoo_Ich.Rda", sep=""))
}



### troubleshooting....###
# FDspr=FData.abn[which(fish$SEASON==slctseason),] # subset to season
## check on overall NA counts
test=is.na(FData)
colSums(test)

## check on NA in zoo where data exists for ich
test2=FData.abn[complete.cases(FData.abn$`74_ich`),]
test3=is.na(test2)
colSums(test3)

unique(FData.bio$YEAR) ## 1977-2019, 1963-1976 WTF
# x=unique(dmrg$Y)
# y=unique(dmrg$zY)
z=seq(1977,2019,1) # years that SHOULD be included
z[!(z %in% unique(FData.bio$YEAR))]
# [1] 1981 1986 1996

# x[!(x %in% y)] # survdat not in zoo
# [1] 1983 1984 1986 1987 1988 1989 1995 1996 # survdat not in zoo
# z[!(z %in% x)] # missing from survdat
# [1] 1997 1998 2018 # missing from survdat
# z[!(z %in% y)] # missing from zoo
# [1] 1983 1984 1986 1987 1988 1989 1995 1996 1997 1998 2018 # missing from zoo merged with survdat

### check svdwide for order?
x=unique(svdwide.bio$YEAR) 
z=seq(1963,2019,1) # years that should be included
z[!(z %in% x)]
# numeric(0)

## OK.... so check svdate....
z=seq(1977,2019,1) # years that should be included
# y=unique(svdate$Y)# survdat
y=unique(svdate2$Y)# survdat
z[!(z %in% y)]
# numeric(0)

## OK..., check dfzdate for zoo
x=unique(dfzdate$zY) #zoo
z=seq(1977,2018,1)
z[!(z %in% x)] # missing from zoo
# [1] 1995 1996 1997 1998 ### problem 1

## check tt and tt2
z=seq(1977,2018,1)
x=unique(tt$Y)
x2=unique(tt2$Y)
x3=c(x,x2)
z[!(z %in% x3)] # missing from zoo

### found error -> taking complete.cases in 'load_zoo_data.R' dropped zoo from 1995-1998
# fixed and updated code, but now have to rerun fuzzy join

library(rfUtilities)
trymc=multi.collinear(FData,n=99, na.rm=T)

#### For biomass trends in along shelf distance, depth, distance to the coast ####
# set wd
setwd("K:/1 RM/2 Plankton Spatial Plots/fish_Kevin")
setwd("C:/Users/ryan.morse/Desktop/Iomega Drive Backup 20171012/1 RM/2 Plankton Spatial Plots/fish_Kevin")
setwd("/home/ryan/Git/NEhabitat")
setwd("C:/Users/ryan.morse/Documents/GitHub/NEhabitat")


# read in depth grid
gdepth=raster("nes_bath_data.nc", band=1)
nescoast=read.csv("test_nes_coast_2.csv", header=TRUE) #evenly spaced, more dense, outside LI, ChesBay
# read in coordinate for along shelf diagnal  diag.csv
diag=read.csv('diag.csv', header = TRUE)

# constants
radt=pi/180
R <- 6371 # Earth mean radius [km]




# select season
selseaon="SPRING"
selseaon="FALL"


# CODE TO READ IN STRATA COMPUTE AREAS, NOW JUST READ IN STRATAREAS dataframe
# readin in strata.shp and compute areas of strata
# TrawlStrata=readShapeSpatial("BTS_strata.shp")
TrawlStrata=rgdal::readOGR("BTS_strata.shp")
## TrawlStrata<-shapefile(file.choose())
plot(TrawlStrata)
AREA<-areaPolygon(TrawlStrata, r=6371000)/10^6

# for array of strata and area and make into dataframe
# stratareas=cbind(TrawlStrata@data$STRATA, AREA)
# colnames(stratareas) <- c("STRATA","AREA")
# stratareas=data.frame(stratareas)
# save(stratareas, file="stratareas.rdata")

# load stratareas
load("stratareas.rdata")
# load("C:/Users/ryan.morse/Desktop/Iomega Drive Backup 20171012/1 RM/2 Plankton Spatial Plots/fish_Kevin/stratareas.rdata")

# get Survdat.RData file
# load(file.choose())
# load('Survdat.rdata')
# load("K:/1 RM/0 R workspaces/Survdat (1).RData")

# trim the data.... need to choose season
# retvars <- c("CRUISE6","STATION","STRATUM","SVSPP","YEAR","SEASON","LAT","LON","DEPTH","ABUNDANCE","BIOMASS")
# survdat <- survdat[retvars]  		
# survdat <- survdat[(survdat$SEASON==selseaon),]

# stata to use
# offshore strata to use
CoreOffshoreStrata<-c(seq(1010,1300,10),1340, seq(1360,1400,10),seq(1610,1760,10))
# inshore strata to use, still sampled by Bigelow
CoreInshore73to12=c(3020, 3050, 3080 ,3110 ,3140 ,3170, 3200, 3230, 3260, 3290, 3320, 3350 ,3380, 3410 ,3440)
# combine
strata_used=c(CoreOffshoreStrata,CoreInshore73to12)

# find records to keep based on core strata
rectokeep=survdat$STRATUM %in% strata_used
rectokeep=svdtunq$STRATUM %in% strata_used

#table(rectokeep)
# add rec to keep to survdat
survdat2=cbind(survdat,rectokeep)
# survdat2=cbind(svdtunq,rectokeep)

# delete record form non-core strata
survdat2=survdat2[rectokeep,]
survdat2=svdtunq[rectokeep,]

# trim the data.... to prepare to find stations only
retvars <- c("CRUISE6","STATION","STRATUM","YEAR")
survdat_stations <- survdat2[retvars]  
# unique reduces to a record per tow
survdat_stations <- unique(survdat_stations)


# find unique tow records only, since length and weight removed
# unique deletes to one record per species
# survdat <- unique(survdat)

# add field with rounded BIOMASS scaler used to adjust distributions
survdat2$LOGBIO <- round(log10(survdat2$BIOMASS *10+10))
survdat2$LOGBIO <- round(log10(survdat2$corBIOMASS *10+10))

# take a look go from 1 to 5?
table(survdat2$LOGBIO)


# make table of strata by year
numtowsstratyr=table(survdat_stations$STRATUM,survdat_stations$YEAR)

# find records to keep based on core strata
rectokeep=stratareas$STRATA %in% strata_used

# add rec to keep to survdat2
stratareas=cbind(stratareas,rectokeep)

# delete record form non-core strata
stratareas_usedonly=stratareas[!stratareas$rectokeep=="FALSE",]

# creat areapertow
areapertow=numtowsstratyr

#compute area covered per tow per strata per year
for(i in 1:length(unique(survdat2$YEAR))){
  areapertow[,i]=stratareas_usedonly$AREA/numtowsstratyr[,i]
}

# change inf to NA and round and out in DF
areapertow[][is.infinite(areapertow[])]=NA
areapertow=round(areapertow)
areapertow=data.frame(areapertow)
colnames(areapertow) <- c("STRATA","YEAR","AREAWT")

# add col to survdat2 for strata weight
survdat2$AREAPERTOW=NA

#fill AREAPERTOW
dimsurvdat2=dim(survdat2)
for (i in 1:dimsurvdat2[1]){
  survdat2$AREAPERTOW[i]=areapertow$AREAWT[which(survdat2$STRATUM[i]==areapertow$STRATA & survdat2$YEAR[i]==areapertow$YEAR)]
}

table(ceiling(survdat2$AREAPERTOW/1000))
table(survdat2$LOGBIO)
table(ceiling(survdat2$AREAPERTOW/1000*survdat2$LOGBIO/9))

# add col to survdat2 for PLOTWT
survdat2$PLOTWT=NA
survdat2$PLOTWT= ceiling(survdat2$AREAPERTOW/1000*survdat2$LOGBIO/9)

table(survdat2$PLOTWT)

# Plot stations
plot(survdat2$LON[survdat2$YEAR==1974],survdat2$LAT[survdat2$YEAR==1974])



# put in shorter name
sdat=survdat2

## find unique tow records only, since length and weight removed unique deletes to one record per species
# survdat <- unique(survdat)
retvars <- c("CRUISE6","STATION","STRATUM","YEAR")
survdat_stations <- sdat[retvars]
allstn=sdat[!duplicated(survdat_stations),] # just keep retvars
allstn$LAT=sdat$LAT[!duplicated(survdat_stations)]
allstn$LON=sdat$LON[!duplicated(survdat_stations)]
allstn$SEASON=sdat$SEASON[!duplicated(survdat_stations)]
allstn$DEPTH=sdat$DEPTH[!duplicated(survdat_stations)]
allstn$SURFTEMP=sdat$SURFTEMP[!duplicated(survdat_stations)]
allstn$SURFSALIN=sdat$SURFSALIN[!duplicated(survdat_stations)]
allstn$BOTTEMP=sdat$BOTTEMP[!duplicated(survdat_stations)]
allstn$BOTSALIN=sdat$BOTSALIN[!duplicated(survdat_stations)]
allstn$latbin=sdat$latbin[!duplicated(survdat_stations)]
allstn$lonbin=sdat$lonbin[!duplicated(survdat_stations)]
allstn$month=sdat$month[!duplicated(survdat_stations)]
allstn$AREAPERTOW=sdat$AREAPERTOW[!duplicated(survdat_stations)]
allstn$PLOTWT=sdat$PLOTWT[!duplicated(survdat_stations)]
## Note LON was changed to one of these above column names in this whole file, hit undo, seems OK but check!!
## seems OK...


# allstn=survdat_stations[!duplicated(survdat_stations),]
# allstn=subset(allstn, select=-c(SVSPP, BIOMASS, wt, LOGBIO)) #drop species for list of all stations
# allstn$stg=NA


### MAKE SELECTIONS FOR FISH SVSPP SEASON AND STAGE ####
SELFISH=74
FISHNAME=Lmf$`Common Name`[which(Lmf$SVSPP==SELFISH)]
selseason="FALL" #"SPRING"|"FALL"
SELSTG="Juv" #'Adt'|'Juv'

FISH=sdat[which(sdat$SVSPP==SELFISH),]
# FISH$wt=0.0069*(FISH$LENGTH^3.08) # individual wt in grams based on length in cm for FISH (fishbase)
FISH$wt=Lmf$a[which(Lmf$SVSPP==SELFISH)]*(FISH$LENGTH^Lmf$b[which(Lmf$SVSPP==SELFISH)]) # individual wt in grams based on length in cm for FISH (fishbase)
test=FISH %>% group_by(CRUISE6,STATION,STRATUM,YEAR) %>% summarise(wtsum=sum(wt, na.rm = T))
retvars <- c("CRUISE6","STATION","STRATUM","YEAR")
survdat_stations <- FISH[retvars]
test2=FISH[!duplicated(survdat_stations),]
FISH=left_join(test2, test)
# FISH$wt=FISH$wt/1000
FISH=left_join(allstn, FISH) # now has all stations, with NA for stn data where svspp not caught
table(FISH$stg) # check
unique(FISH$SVSPP)

#now fix NA for LOGBIO, etc NA -> 0
# x[c("a", "b")][is.na(x[c("a", "b")])] <- 0

FISH.juv=FISH[which(FISH$stg=="juv" | is.na(FISH$stg)),]
FISH.adt=FISH[which(FISH$stg=="adt" | is.na(FISH$stg)),]

### NOW do calc for ASDIST, COB, DTC, DPTH ###
# put in shorter name

if(SELSTG=="Juv") {
  stage=FISH.juv
} else {
  stage=FISH.adt
}


#_______________________________
### below is from EcoMon processing of this type of data ###
####  Geosphere package to calc distance to coastline from pts (lon,lat), returns meters
dd = array(data = NA, dim = nrow(stage))
pts = data.frame(stage$LON, stage$LAT)
#line = t(rbind(nescoast$lon, nescoast$lat))
dd=dist2Line(pts[,], nescoast)
stage$dtc=dd[,1]/1000 # convert meters to KM

# Find distance to diagonal line (diag), use coordinates of nearest point to find distance to NC outerbanks (min(diag))
dd2 = array(data = NA, dim = nrow(stage))
dd2 = dist2Line(pts[,], diag, distfun=distHaversine)
#Distance of closest point to data along diag line to NC coast
p1 = diag[1,] #start of line
p2 = data.frame(dd2[,2], dd2[,3])
distNC = distCosine(p1, p2, r=6378137) /1000 # convert to KM (Great circle distance)
stage$distNC = distNC



#### MISSING DEPTHS (look for 9999) # fill only those records in misdepth with depth from grid
# create column for missing depth data intially with depth data
stage$misdepth=stage$DEPTH
missingdepth=which(is.na(stage$misdepth)) # == 9999)# find cases with missing depth data# =which(is.na(stage$depth))
for(k in missingdepth){
  stage$misdepth[k] = extract(gdepth,cbind(stage$LON[k],stage$LAT[k])) # * -1
}

## split out to season ##
stagesn=stage[which(stage$SEASON==selseason),]



stagesn$bio=stagesn$wtsum / stagesn$AREAPERTOW
stagesn$lgbio=floor(log10(stagesn$bio+1))
taxa=data.frame(stagesn$lgbio)
taxa[is.na(taxa)]=0
colnames(taxa)=FISHNAME

out_data=array(NA,c((max(stagesn$YEAR)-min(stagesn$YEAR)+1)*length(taxa),7))
row_c=0

# setwd('K:/1 RM/2 Plankton Spatial Plots/fish_Kevin')
# block to calculate summary final distance of populations
for (ii in 1:length(taxa)){
  print(c(colnames(taxa)[ii], SELSTG, selseason))
  for(j in min(stagesn$YEAR):max(stagesn$YEAR)){
    row_c=row_c+1 
    sumdistA=sum(stagesn$distNC[stagesn$YEAR==j] *taxa[stagesn$YEAR==j,ii]) #ASDIST
    lendist=sum(taxa[stagesn$YEAR==j,ii])
    mdist =sumdistA / lendist  
    
    sumdistB=sum(stagesn$dtc[stagesn$YEAR==j] * taxa[stagesn$YEAR==j,ii]) #DTOC
    sdtoc =sumdistB / lendist
    
    sumdistC=sum(stagesn$misdepth[stagesn$YEAR==j] *taxa[stagesn$YEAR==j, ii]) #Depth
    mdepth =sumdistC / lendist
    
    sumdistD=sum(stagesn$lat[stagesn$YEAR==j] * taxa[stagesn$YEAR==j,ii]) #Lat
    mlat=sumdistD/lendist
    
    sumdistE=sum(stagesn$lon[stagesn$YEAR==j] * taxa[stagesn$YEAR==j,ii]) #Lon
    mlon=sumdistE/lendist
    
    # outline=paste(j,colnames(taxa)[ii],SELSTG,selseason,mdist,sdtoc,mdepth,mlat,mlon, sep=",")
    # write.table(outline,file=paste(selseaon, "test_dis_depth.csv", sep='_'),row.name=F,col.names=F,append=TRUE)
    
    out_data[row_c,1]=j
    out_data[row_c,2]=paste(colnames(taxa)[ii], SELSTG, selseason, sep=" ")
    out_data[row_c,3]=mdist
    out_data[row_c,4]=sdtoc
    out_data[row_c,5]=mdepth
    out_data[row_c,6]=mlat
    out_data[row_c,7]=mlon
  }  
}

out_data=data.frame(out_data)



names(out_data)[names(out_data)=="X1"] <- "YR"
names(out_data)[names(out_data)=="X2"] <- "SP"
names(out_data)[names(out_data)=="X3"] <- "ASD"
names(out_data)[names(out_data)=="X4"] <- "DTC"
names(out_data)[names(out_data)=="X5"] <- "DEP"
names(out_data)[names(out_data)=="X6"] <- "LAT"
names(out_data)[names(out_data)=="X7"] <- "LON"


write.csv(out_data,file=outfile )

#### plotting ASD, DTC, Z, etc
setwd("C:/Users/ryan.morse/Documents/GitHub/NEhabitat")
nesbath=getNOAA.bathy(lon1=-77,lon2=-65,lat1=35,lat2=45, resolution=10, keep=F)
df=read.csv('SPRING_test_dis_depth.csv', check.names = F, stringsAsFactors = F, col.names = c('Year', 'Species', 'Stg', 'Season', 'ASD', 'DTC', 'Z', 'Lat', 'Lon'))
library(dplyr)
i=unique(df$Species)
j=unique(df$Season)
k=unique(df$Stg)

l=2 #1 for Cod, 2 Haddock, etc
m=2 #1 for Spring
n=2 #1 for Adt, 2 for Juv

pal <- colorRampPalette(c("blue", "yellow", "red"))
# colr=pal(48)
ss=df[which(df$Species==i[l] & df$Season==j[m] & df$Stg==k[n]),] # ss=filter(df, df$Species==i[l], df$Season==j[m], df$Stg==k[n])
colr=pal(length(ss$Year)) #make sure colors are correct length

## Plot ASD
plot(x=ss$Year,y=ss$ASD,pch=16,col=colr,main=paste('ASD ', i[l], j[m], k[n]))
lines(x=ss$Year,y=ss$ASD,col = "gray50")
points(x=ss$Year,y=ss$ASD,pch=16,col=colr)
## DTC
plot(x=ss$Year,y=ss$DTC,pch=16,col=colr, main=paste('DTC ', i[l], j[m], k[n]))
lines(x=ss$Year,y=ss$DTC,col = "gray50")
points(x=ss$Year,y=ss$DTC,pch=16,col=colr)
## Depth
plot(x=ss$Year,y=ss$Z,pch=16,col=colr, main=paste('Z ', i[l], j[m], k[n]))
lines(x=ss$Year,y=ss$Z,col = "gray50")
points(x=ss$Year,y=ss$Z,pch=16,col=colr)
## Lat Lon
#Lat Lon
plot(x=ss$Lon,y=ss$Lat,pch=16,col=colr, main=paste('Loc ', i[l], j[m], k[n]))
lines(x=ss$Lon,y=ss$Lat,col = "gray50")
points(x=ss$Lon,y=ss$Lat,pch=16,col=colr)

# map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map("worldHires", xlim=c(-72,-65),ylim=c(40,45), fill=T,border=0,col="gray70")
map.axes(las=1)
plot(nesbath,deep=-100, shallow=-100, step=1,add=T,lwd=1,col="gray80",lty=1)
lines(x=ss$Lon,y=ss$Lat,col = "gray50")
points(x=ss$Lon,y=ss$Lat,pch=16,col=colr)
