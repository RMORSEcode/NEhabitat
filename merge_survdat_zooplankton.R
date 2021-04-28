# require(raster)
library(lubridate)
library(readxl)
library(dplyr)
library(tidyr)
# library(sp)
# library(maptools)
# library(marmap)
# library(geosphere)


loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

Lmf=read_excel('/home/ryan/Git/NEhabitat/Lm_included.xlsx')
### 20210330 read newest data for 2021 SOE, from Harvey Walsh
ZPD=read.csv('/home/ryan/Desktop/Z/SOE_2021_PlanktonData/EcoMon_Plankton_Data_v3_7_DoNotDistribute.csv', stringsAsFactors = F)
# dt=as_date(ZPD$date)#, origin = "1899-12-30")
dt=as.Date(ZPD$date, format="%d-%b-%y")
ichnms=read.csv('/home/ryan/1_habitat_analysis_2017/ichthyonames.csv', header = F, stringsAsFactors = F)
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
svspplu=read.csv('/home/ryan/1_habitat_analysis_2017/svspp_lookup.csv', stringsAsFactors = F)


svnms=svspplu[svspplu$SCINAME %in% gnms$spp,]
svnms2=left_join(svnms, gnms, by=c("SCINAME"="spp"))
svnms$ich=svnms2$V1
rm(svnms2)
# svnms2=svnms[,c(1,5)]
svnms2=svnms[,c(1,4)]

ZPD$date=dt
ZPD$DOY=yday(ZPD$date) #day of year
# month=as.numeric(format(dt, '%m'))
# year=as.numeric(format(dt, '%Y'))
ZPD$year=year(ZPD$date)
ZPD$month=month(ZPD$date)
ZPD$dt=dt
# ZPD$DOY=DOY
# ZPD$day=as.numeric(format(dt, '%d'))
ZPD$day=day(ZPD$date)
ZPD$lat2=ceiling(ZPD$lat) #use for binning into 1 degree bins for removal of undersampled bins
ZPD$lon2=floor(ZPD$lon) #use for binning into 1 degree bins for removal of undersampled bins

test=ZPD[complete.cases(ZPD$melaeg_100m3 ),]
# test=test[complete.cases(test$calfin_100m3),]
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
## add SVSPP and stage to ich, transform to long format
ich$stg="ich"
ich=select(ich, -c(nofish_100m3)) #drop nofish (all zeros)
ichlong=ich %>% gather(ichsp, abundance, urospp_100m3:lopame_100m3)
ichlong=left_join(ichlong, svnms2, by=c("ichsp"="ich"))

df=ich
dfv=vars
dfz=test



## Subset survdat to only groundfish species with Length at maturity data
# load("C:/Users/ryan.morse/Downloads/SurvdatBio (3).RData") #loads survdat.bio
load("C:/Users/ryan.morse/Downloads/survdat_lw.RData") # loads survdat.lw ; from Sean Lucey data into 2019
load("/home/ryan/Downloads/survdat_lw.RData")# loads survdat.lw ; from Sean Lucey data into 2019
newsurvdat=loadRData('/home/ryan/Downloads/NEFSC_BTS_2021.RData') # new data from Sean 2021, possible updates to haddock due to code q correction
Lmf=read_excel('/home/ryan/Git/NEhabitat/Lm_included.xlsx')
test=newsurvdat$survdat$SVSPP %in% Lmf$SVSPP
nsdat=newsurvdat$survdat[test,]

## choose surdat and remove older non-length version to avoid confusion
rm(survdat)
survdat=survdat.lw
survdat=survdat.bio
survdat=nsdat

# NEED TO CHECK ON THIS PART, THE LENGTHS ARE THE SAME, STATIONS MAY NOT BE
## newest survdat missing INDWT, WGTLEN, SIZECAT from survdat.lw, remove 2020 and add back
test=survdat %>% filter(YEAR!=2020)
test2=test %>% filter(SVSPP==074, YEAR==2010)
test3=survdat.lw %>% filter(SVSPP==074, YEAR==2010)
## they are the same
test$INDWT=survdat.lw$INDWT
test$WGTLEN=survdat.lw$WGTLEN
test$SIZECAT=survdat.lw$SIZECAT
survdat=test

test=data.frame(survdat[,c("SVSPP", "LENGTH")])
colnames(test)=c("SVSPP", "LENGTH")
test=left_join(test, Lmf[,c("SVSPP", "Lm")], by="SVSPP")
test$stg=ifelse(test$LENGTH<test$Lm, "Juv", "Adt")
survdat$stg=test$stg
rm(test)


# 20210406 ONLY for newest survdat - missing wgtlen, only correct abundance
survdat=survdat %>% group_by_at(vars(SVSPP,CRUISE6,STATION,STRATUM,ABUNDANCE,stg)) %>% mutate(stgsum=sum(NUMLEN)) #, stgwtsum=sum(WGTLEN))
survdat=survdat %>% group_by_at(vars(SVSPP,CRUISE6,STATION,STRATUM,ABUNDANCE)) %>% mutate(totlen=sum(NUMLEN), abndiff=ABUNDANCE-totlen)
survdat=survdat %>% group_by_at(vars(SVSPP,CRUISE6,STATION,STRATUM,ABUNDANCE, stg)) %>%
  mutate(stgpctabn=stgsum/totlen, corABN=round(stgpctabn*ABUNDANCE,digits = 0))

### now create dataframe with unique tows only, separated by stage
test=survdat[,c("CRUISE6","STATION","STRATUM","YEAR", "SVSPP", "stg")] # needs SVSPP...
svdtunq=survdat[!duplicated(test),]

svdtunq=survdat ## keep all stations for now, remove latre


### SKIP THIS IF USING nsdat and above numbers correction
### count juv and adult to see if numbers match ABUNDANCE 
# survdat=survdat %>% group_by_at(vars(SVSPP,CRUISE6,STATION,STRATUM,ABUNDANCE,stg)) %>% mutate(stgsum=sum(NUMLEN), stgwtsum=sum(WGTLEN))
### now addd up weights for all fish measured, check against biomass, total abundance
# survdat=survdat %>% group_by_at(vars(SVSPP,CRUISE6,STATION,STRATUM,ABUNDANCE)) %>% mutate(totlen=sum(NUMLEN), totwgt=sum(WGTLEN), abndiff=ABUNDANCE-totlen, biodiff=BIOMASS-totwgt)
### now add corrected biomass and abundance to account for any missing staged values (eg. for large hauls)
# survdat=survdat %>% group_by_at(vars(SVSPP,CRUISE6,STATION,STRATUM,ABUNDANCE, stg)) %>%
  # mutate(stgpctabn=stgsum/totlen, stgpctwgt=stgwtsum/totwgt, corABN=round(stgpctabn*ABUNDANCE,digits = 0), corBIO=round(stgpctwgt*BIOMASS, digits=2))

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





