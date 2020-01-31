# Kevin's code to calculate center of biomass and movement for trawl survey data


require(raster)
# require(ncdf)
library(sp)
library(maptools)
library(geosphere)

## Subset survdat to only groundfish species with Length at maturity data
load("C:/Users/ryan.morse/Downloads/SurvdatBio (3).RData")
test=survdat.bio$SVSPP %in% Lmf$SVSPP
survdat.bio2=survdat.bio[test,]
sort(unique(survdat.bio2$SVSPP)) #verify
length(unique(survdat.bio2$SVSPP))

test=survdat$SVSPP %in% Lmf$SVSPP
survdat=survdat[test,]
sort(unique(survdat$SVSPP)) #verify
length(unique(survdat$SVSPP))

survdat$stg=NA
for (j in 1:length(survdat$stg)){
  survdat$stg[j]=ifelse(survdat$LENGTH[j]<Lmf$Lm[which(Lmf$SVSPP==survdat$SVSPP[j])], "juv", "adt")
}


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

#table(rectokeep)
# add rec to keep to survdat
survdat2=cbind(survdat,rectokeep)

# delete record form non-core strata
survdat2=survdat2[rectokeep,]

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


### now reduce to single species, then add back NA for entire record, then split -> adult and juvenile
cod=sdat[which(sdat$SVSPP==73),]
cod$wt=0.0069*(cod$LENGTH^3.08) # individual wt in grams based on length in cm for cod (fishbase)
# test=aggregate(cod$wt, by=list(c(cod$CRUISE6, cod$STATION, cod$STRATUM, cod$YEAR, cod$SVSPP)), FUN=sum)
# test=cod %>% group_by((cod$CRUISE6, cod$STATION, cod$STRATUM, cod$YEAR, cod$SVSPP)) %>% summarise(wtsum=sum(cod$wt))
test=cod %>% group_by(CRUISE6,STATION,STRATUM,YEAR) %>% summarise(wtsum=sum(wt, na.rm = T))
# test2=cod[unique(cod[c("CRUISE6","STATION","STRATUM","YEAR")]),]
retvars <- c("CRUISE6","STATION","STRATUM","YEAR")
survdat_stations <- cod[retvars]
test2=cod[!duplicated(survdat_stations),]
# stn <- as.numeric(rownames(survdat_stations)) # get index of unique stations
# cod=cod[stn,] # this is not working...
cod=left_join(test2, test)
# cod$wt=cod$wt/1000
cod=left_join(allstn, cod) # now has all stations, with NA for stn data where svspp not caught
table(cod$stg) # check
unique(cod$SVSPP)

#now fix NA for LOGBIO, etc NA -> 0
# x[c("a", "b")][is.na(x[c("a", "b")])] <- 0

cod.juv=cod[which(cod$stg=="juv" | is.na(cod$stg)),]
cod.adt=cod[which(cod$stg=="adt" | is.na(cod$stg)),]

had=survdat2[which(survdat2$SVSPP==74),]
had$wt=0.0059*(had$LENGTH^3.13) # individual wt in grams based on length in cm for had (fishbase)
had.juv=had[which(had$stg=="juv" | is.na(had$stg)),]
had.adt=had[which(had$stg=="adt" | is.na(had$stg)),]

red=survdat2[which(survdat2$SVSPP==155),]
red$wt=0.018*(red$LENGTH^2.966) # K. Duclos 2015 thesis UNH
red.juv=red[which(red$stg=="juv" | is.na(red$stg)),]
red.adt=red[which(red$stg=="adt" | is.na(red$stg)),]



### NOW do calc for ASDIST, COB, DTC, DPTH ###
# put in shorter name
stage=cod.juv
stage=cod.adt


# read species list  sps.csv
# sps=read.csv('sps.csv', header = TRUE)
# # sps=read.csv('C:/Users/ryan.morse/Desktop/for ryan/sps.csv') #updated for NEUS
# numsps=nrow(sps)
# 
# numrecs=nrow(stage)
# 
# d = array(data = NA, dim = 150)
# 
# 
# for (j in 1:numrecs) {
#   
#   print(numrecs-j)
#   
#   
#   lat1=stage$LAT[j]* radt
#   long1=stage$LON[j]* radt
#   for (i in 1:150){
#     lat2=diag$Latitude[i]* radt
#     long2=diag$Longitude[i]* radt
#     d[i] <- acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2) * cos(long2-long1)) * R
#   }
#   dindex=which(d==min(d))
#   
#   lat1=34.60* radt
#   long1=-76.53* radt
#   
#   lat2=diag$Latitude[dindex]* radt
#   long2=diag$Longitude[dindex]* radt
#   stage$ASDIST[j] = acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2) * cos(long2-long1)) * R
#   
# }
# 
# 
# 
# 
# 
# # create column for missing depth data intially with depth data
# stage$MISDEPTH=stage$DEPTH
# 
# # find cases with missing depth data
# missingdepth=which(is.na(stage$DEPTH))
# 
# # fill only those records in misdepth with depth from grid
# for(k in missingdepth){
#   stage$MISDEPTH[k] = extract(gdepth,cbind(stage$LON[k],stage$LAT[k])) * -1
# }
# 
# 
# 
# 
# 
# # for (i in 1:numsps){
# #   print (i)
#   for(j in min(stage$YEAR):max(stage$YEAR)){
#     
#     sumdist=sum(stage$ASDIST[stage$YEAR==j & stage$SVSPP==sps$SVSPP[i]] *stage$PLOTWT[stage$YEAR==j & stage$SVSPP==sps$SVSPP[i]])
#     lendist=sum(stage$PLOTWT[stage$YEAR==j & stage$SVSPP==sps$SVSPP[i]])
#     mdist =sumdist / lendist
#     
#     sumdepth=sum(stage$MISDEPTH[stage$YEAR==j & stage$SVSPP==sps$SVSPP[i]] *stage$PLOTWT[stage$YEAR==j & stage$SVSPP==sps$SVSPP[i]])
#     lendepth=sum(stage$PLOTWT[stage$YEAR==j & stage$SVSPP==sps$SVSPP[i]])
#     mdepth =sumdepth / lendepth
#     
#     outline=paste(j,",",sps$SVSPP[i],",",mdist,",",mdepth)
#     write.table(outline,file="disdepthtest.csv",row.name=F,col.names=F,append=TRUE)
#   }  
# # }
# 
# # missing depths ???
# 
# stage$DEPTH[stage$SVSPP==73 & stage$YEAR==1973]

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

# copes=c(9,10,11,13,14,17,19,20,23,29,30,31)
# KEEP=c(seq(9,31),33,36,40,57) #26 taxa similar to NES regime shift taxa
# keep.many=seq(9,57) # lots, including original

#For ii species:
# taxa=ZPDa[,KEEP]
# taxa=floor(log10(taxa+1))
selseason="SPRING"
stage=stage[which(stage$SEASON==selseaon),]

stage$bio=stage$wtsum / stage$AREAPERTOW
stage$lgbio=floor(log10(stage$bio+1))
taxa=data.frame(stage$lgbio)
taxa[is.na(taxa)]=0
colnames(taxa)="Adt_Cod"



# setwd('K:/1 RM/2 Plankton Spatial Plots/fish_Kevin')
# block to calculate summary final distance of populations
for (ii in 1:length(taxa)){
  print(colnames(taxa)[ii])
  for(j in min(stage$YEAR):max(stage$YEAR)){
    sumdistA=sum(stage$distNC[stage$YEAR==j] *taxa[stage$YEAR==j,ii]) #ASDIST
    lendist=sum(taxa[stage$YEAR==j,ii])
    mdist =sumdistA / lendist  
    
    sumdistB=sum(stage$dtc[stage$YEAR==j] * taxa[stage$YEAR==j,ii]) #DTOC
    sdtoc =sumdistB / lendist
    
    sumdistC=sum(stage$misdepth[stage$YEAR==j] *taxa[stage$YEAR==j, ii]) #Depth
    mdepth =sumdistC / lendist
    
    sumdistD=sum(stage$lat[stage$YEAR==j] * taxa[stage$YEAR==j,ii]) #Lat
    mlat=sumdistD/lendist
    
    sumdistE=sum(stage$lon[stage$YEAR==j] * taxa[stage$YEAR==j,ii]) #Lon
    mlon=sumdistE/lendist
    
    # outline=paste(j,",",colnames(taxa)[ii],",",mdist,",",sdtoc, ',',mdepth, ',',mlat,',',mlon)
    outline=paste(j,colnames(taxa)[ii],mdist,sdtoc,mdepth,mlat,mlon, sep=",")
    
    write.table(outline,file=paste(selseaon, "test_dis_depth.csv", sep='_'),row.name=F,col.names=F,append=TRUE)
  }  
}
