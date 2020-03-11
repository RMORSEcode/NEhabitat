# Kevin's code to calculate center of biomass and movement for trawl survey data


require(raster)
# require(ncdf)
library(sp)
library(maptools)
library(marmap)
library(geosphere)
library(lubridate)
library(dplyr)
library(tidyr)

## Subset survdat to only groundfish species with Length at maturity data
# load("C:/Users/ryan.morse/Downloads/SurvdatBio (3).RData") #loads survdat.bio
load("C:/Users/ryan.morse/Downloads/survdat_lw.RData") # loads survdat.lw ; from Sean Lucey data into 2019
load("/home/ryan/Downloads/survdat_lw.RData")# loads survdat.lw ; from Sean Lucey data into 2019
test=survdat.lw$SVSPP %in% Lmf$SVSPP
survdat.lw=survdat.lw[test,]
sort(unique(survdat.lw$SVSPP)) #verify
length(unique(survdat.lw$SVSPP))

test=survdat$SVSPP %in% Lmf$SVSPP
survdat=survdat[test,]
sort(unique(survdat$SVSPP)) #verify
length(unique(survdat$SVSPP))

## choose surdat and remove older non-length version to avoid confusion
rm(survdat)
survdat=survdat.lw
survdat=survdat.bio


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
test=data.frame(survdat[,c("SVSPP", "CRUISE6","YEAR", "STATION", "STRATUM", "ABUNDANCE", "BIOMASS","NUMLEN", "WGTLEN", "stg")])
# test=data.frame(survdat[1:25,c("SVSPP", "CRUISE6","YEAR", "STATION", "STRATUM", "ABUNDANCE", "NUMLEN", "stg")])
test2=test %>% group_by_at(vars(SVSPP,CRUISE6,STATION, STRATUM,ABUNDANCE,stg)) %>% mutate(stgsum=sum(NUMLEN), wtsum=sum(WGTLEN))
survdat$stgsum=test2$stgsum
survdat$stgwtsum=test2$wtsum
### now addd up weights for all fish measured, check against biomass, total abundance
test2=test %>% group_by_at(vars(SVSPP,CRUISE6,STATION, YEAR, STRATUM,ABUNDANCE)) %>% mutate(totlen=sum(NUMLEN), totwgt=sum(WGTLEN))
test2$abndiff=test2$ABUNDANCE-test2$totlen
test2$wgtdiff=test2$BIOMASS-test2$totwgt

survdat$totlen=test2$totlen
survdat$totwgt=test2$totwgt
survdat$abndiff=survdat$ABUNDANCE - survdat$totlen
survdat$biodiff=survdat$BIOMASS - survdat$totwgt


barplot(table(round(survdat$biodiff[(abs(survdat$biodiff>12))],1)))
y=survdat$SVSPP
x=abs(survdat$biodiff)
y=y[x>120]
x=x[x>120]
barplot(table(round(y))) #which species have the highest mismatches in biomass
barplot(table(round(x))) #how many

### now create dataframe with unique tows only, separated by stage
test=survdat[,c("CRUISE6","STATION","STRATUM","YEAR", "SVSPP", "stg")] # needs SVSPP...
svdtunq=survdat[!duplicated(test),]
### calc percent of biomass by stage(adt, juv) per unique tow to be applied to BIOMASS for corrected values
svdtunq$stgwgtpct=round(svdtunq$stgwtsum/svdtunq$totwgt, 2)
svdtunq$corBIOMASS=svdtunq$stgwgtpct*svdtunq$BIOMASS

svdunqadt=svdtunq %>% filter(stg=="Adt") %>% 
  select(-(LENGTH:SIZECAT)) %>% 
  select(-(CATCHSEX)) %>% 
  spread(SVSPP, ABUNDANCE:corBIOMASS)

barplot(table(svdtunq$stgwgtpct))
barplot(table(svdtunq$stgwgtpct[svdtunq$SEASON=='SPRING']))
barplot(table(svdtunq$stgwgtpct[svdtunq$SEASON=='FALL']))
barplot(table(svdtunq$stgwgtpct[svdtunq$stgwgtpct<0.99 & svdtunq$SEASON=='FALL']),ylim=c(0,1000), main='FALL')
barplot(table(svdtunq$stgwgtpct[svdtunq$stgwgtpct<0.99 & svdtunq$SEASON=='SPRING']),ylim=c(0,1000), main='SPRING')

svspp=73
svdtunq$stgwgtpct[is.na(svdtunq$stgwgtpct)]=0 #[svdtunq$stgwgtpct<0.99 & svdtunq$SEASON=='FALL' & svdtunq$SVSPP==svspp]),]

barplot(table(svdtunq$stgwgtpct[svdtunq$stgwgtpct<0.99 & svdtunq$SEASON=='FALL' & svdtunq$SVSPP==svspp]),ylim=c(0,100), main=paste('FALL ', svspp))
barplot(table(svdtunq$stgwgtpct[svdtunq$stgwgtpct<0.99 & svdtunq$SEASON=='SPRING' & svdtunq$SVSPP==svspp]),ylim=c(0,100), main=paste('SPRING ', svspp))


library(modes)
# restricting to <0.99 was masking adult population. Can check on change over years...
bimodality_amplitude(svdtunq$stgwgtpct[svdtunq$SEASON=='SPRING' & svdtunq$SVSPP==svspp], fig=T)
bimodality_amplitude(svdtunq$stgwgtpct[svdtunq$SEASON=='FALL' & svdtunq$SVSPP==svspp], fig=T)
sum(is.na(svdtunq$stgwgtpct[svdtunq$stgwgtpct<0.99 & svdtunq$SEASON=='FALL' & svdtunq$SVSPP==svspp]))

### create date and binned lat lon for mathching to other surveys
svdate=data.frame(month(survdat$EST_TOWDATE))
colnames(svdate)[1]='M'
svdate$Y=year(survdat$EST_TOWDATE)
svdate$D=day(survdat$EST_TOWDATE)
# svdate$lat=survdat$LAT
# svdate$lon=survdat$LON
svdate$lonbin=round(survdat$LON/0.5)*0.5 ## half degree
svdate$latbin=round(survdat$LAT/0.5)*0.5
svdate$doy=as.numeric(strftime(survdat$EST_TOWDATE, format = "%j"))
svdate$sdoy=svdate$doy-8
svdate$ldoy=svdate$doy+8
svdate$index=seq(from=1, to=length(svdate$M), by=1)


## get index of unique stations from survdat, apply to long format svdate with index
retvars <- c("CRUISE6","STATION","STRATUM","YEAR")
sdq=survdat[,retvars]
sdq$index=seq(from=1, to=length(survdat$YEAR), by=1)
survdat_stations <- survdat[retvars] 
sdq2=sdq[!duplicated(survdat_stations),] # just keep retvars
svdate2=svdate[sdq2$index,]

#clean up svdate to match zooplankton and Chl
svdate2=svdate2[which(svdate2$Y>2017),]


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

# try merge and filter (did this in 2 parts, before and after 1998, took long time...)
dfzdate2=dfzdate[which(dfzdate$zY>2017),]
# ttx=merge(svdate2, dfzdate2, all=T)
colnames(dfzdate2)
colnames(svdate2)
# ttx=left_join(svdate2, dfzdate2, by=c("lonbin"="zlonbin", "latbin"="zlatbin", "Y"="zY", "sdoy"="zdoy", "ldoy"="zdoy"))
### this works, just takes a long time
library(fuzzyjoin)
tt=fuzzy_left_join(svdate2, dfzdate2, by=c("lonbin"="zlonbin", "latbin"="zlatbin", "Y"="zY", "sdoy"="zdoy", "ldoy"="zdoy"), 
                   match_fun=list(`==`,`==`,`==`,`<=`,`>=`))
tt2=tt[complete.cases(tt$zY),]
# merge 2 dataframes together
dmrg=rbind(tt3, tt2)

## now subset original dataframes before joining together
fish1=survdat[dmrg$index,]
zoo1=dfz[dmrg$zindex,]
ich1=ich[dmrg$zindex,]


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
