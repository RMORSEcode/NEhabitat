library(mgcv)
library(gratia)
# library(car)
library(MASS)
library(stringr)
library(tidyr)
library(lubridate)
library(raster)
library(ROCR)
library(dplyr)

saveAUC <- function(truth, predicted, ...){
  pred <- prediction(as.vector(abs(predicted)), as.vector(truth))    
  roc <- performance(pred,"tpr","fpr")
  auc <- performance(pred,"auc")
  # plot(roc, ...)
  # abline(a=0, b= 1)
  return(auc@y.values)
}

plotROC <- function(truth, predicted, ...){
  pred <- prediction(as.vector(abs(predicted)), as.vector(truth))    
  roc <- performance(pred,"tpr","fpr")
  auc <- performance(pred,"auc")
  plot(roc, ...)
  text(0.6, 0.1, paste('AUC=', auc@y.values, sep=''))
  abline(a=0, b= 1)
  # print(auc@y.values)
}

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

### Load the data and rename ichthyoplankton
load('/home/ryan/Git/NEhabitat/Final_merged_fish_corBIO_Zoo_Ich.Rda') # Biomass
load('/home/ryan/Git/NEhabitat/Final_merged_fish_corABN_Zoo_Ich.Rda') # Abundance
## newest with updated survdat values
load('/home/ryan/Git/NEhabitat/20210326_unique_wide_format_corABN_by_stage.Rda')# Abundance
load('/home/ryan/Git/NEhabitat/20210326_unique_wide_format_corBIO_by_stage.Rda')# Biomass
test=FData.abn %>% left_join(svdwide.abn, by = c("CRUISE6", "STATION", "STRATUM", "SVVESSEL", 
                                                 "YEAR", "EST_TOWDATE", "SEASON", "LAT", "LON"))

og=FData.abn %>% filter(YEAR>1976)
og$MONTH=month(og$EST_TOWDATE)
tunq=og %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()==1) %>% mutate(num=n())
tdup=og %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()>1) %>% mutate(num=n())
tdupmn=tdup %>% mutate_if(is.numeric, mean) %>% mutate_if(is.character, funs(paste(unique(.), collapse = "_"))) %>% slice(1)
og2=bind_rows(tunq, tdupmn) %>% ungroup() #reassign name
og2[with(og2, order("EST_TOWDATE", "YEAR", "MONTH", "STATION", "STRATUM", "LAT", "LON")),]
og2$EST_TOWDATE=as.character(og2$EST_TOWDATE)

new=svdwide.abn %>% filter(YEAR>1976)
new$MONTH=month(new$EST_TOWDATE)
tunq=new %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()==1) %>% mutate(num=n())
tdup=new %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()>1) %>% mutate(num=n())
tdupmn=tdup %>% mutate_if(is.numeric, mean) %>% mutate_if(is.character, funs(paste(unique(.), collapse = "_"))) %>% slice(1)
new2=bind_rows(tunq, tdupmn) %>% ungroup() #reassign name
new2[with(new2, order("EST_TOWDATE", "YEAR", "MONTH", "STATION", "STRATUM", "LAT", "LON")),]
new2$EST_TOWDATE=as.character(new2$EST_TOWDATE)
# keep original, merge in new values
test=og2 %>% left_join(new2, by = c("CRUISE6", "STATION", "STRATUM", "SVVESSEL", 
                                                 "YEAR", "EST_TOWDATE", "LAT", "LON"))
## keep new values from survdat, merge ich and zoo
test=new2 %>% left_join(og3, by =c("CRUISE6", "STATION", "STRATUM", "SVVESSEL", "YEAR",  
                                   "LAT", "LON", "DEPTH")) #, "SURFTEMP", "BOTTEMP", "SURFSALIN", 
                                   # "BOTSALIN", "72_Juv", "72_Adt", "73_Juv", "73_Adt", "74_Juv", "74_Adt", 
                                   # "75_Adt", "75_Juv", "76_Adt", "76_Juv", "77_Juv", "77_Adt", "101_Juv", 
                                   # "101_Adt", "102_Juv", "102_Adt", "103_Adt", "103_Juv", "105_Juv", 
                                   # "105_Adt", "106_Adt", "106_Juv", "107_Juv", "107_Adt", "108_Adt", 
                                   # "108_Juv", "155_Juv", "155_Adt", "192_Juv", "192_Adt", "193_Adt", 
                                   # "193_Juv", "197_Juv", "197_Adt"))
test=new2 %>% left_join(og2, by =c("CRUISE6", "STATION", "STRATUM", "SVVESSEL", "YEAR",  
                                   "LAT", "LON", "DEPTH", "SURFTEMP", "BOTTEMP", "SURFSALIN", 
"BOTSALIN", "72_Juv", "72_Adt", "73_Juv", "73_Adt", "74_Juv", "74_Adt",
"75_Adt", "75_Juv", "76_Adt", "76_Juv", "77_Juv", "77_Adt", "101_Juv",
"101_Adt", "102_Juv", "102_Adt", "103_Adt", "103_Juv", "105_Juv",
"105_Adt", "106_Adt", "106_Juv", "107_Juv", "107_Adt", "108_Adt",
"108_Juv", "155_Juv", "155_Adt", "192_Juv", "192_Adt", "193_Adt",
"193_Juv", "197_Juv", "197_Adt"))# %>% filter(complete.cases(`74_ich`))

## something is not right, many NAs where there should be matchups
og2[1,c("CRUISE6", "STATION", "STRATUM", "SVVESSEL", "YEAR", "EST_TOWDATE", 
       "SEASON", "LAT", "LON", "DEPTH")]
new2[1,c("CRUISE6", "STATION", "STRATUM", "SVVESSEL", "YEAR", "EST_TOWDATE", 
         "SEASON", "LAT", "LON", "DEPTH")]

new3=new2[1:15,]
og4=og3[1:15,]
og4$STATION[which(!(og4$STATION %in% new3$STATION))]
# new3$EST_TOWDATE=as.character(new3$EST_TOWDATE)
# og3$EST_TOWDATE=as.character(og3$EST_TOWDATE)
# og3[1,c("CRUISE6", "STATION", "STRATUM", "SVVESSEL", "YEAR", "EST_TOWDATE", "SEASON", "LAT", "LON", "DEPTH")]==new3[1,c("CRUISE6", "STATION", "STRATUM", "SVVESSEL", "YEAR", "EST_TOWDATE", "SEASON", "LAT", "LON", "DEPTH")]

cog2=lapply(og2, class)
cnew2=lapply(new2, class)

## clean up
og3=og2 %>% select(-mrgidx, date)
og3$YEAR=as.integer(og3$YEAR)
og3$SEASON=as.character(og3$SEASON)
og3$SVVESSEL=as.character(og3$SVVESSEL)
og3$EST_TOWDATE=as.character(og3$EST_TOWDATE)

new3$EST_TOWDATE=as.character(new3$EST_TOWDATE)


### 20210406 merge original FData.abn with new survdat
test=FData.abn %>% left_join(svdwide.abn, by =c("CRUISE6", "STATION", 
                                                "STRATUM", "SVVESSEL", "YEAR","LAT", "LON", "DEPTH"))
# tunq=test %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()==1) %>% mutate(num=n())
# tdup=test %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()>1) %>% mutate(num=n())

tunq=test %>% group_by(CRUISE6, STATION, STRATUM, SVVESSEL, YEAR,LAT, LON, DEPTH) %>% filter(n()==1) %>% mutate(num=n())
tdup=test %>% group_by(CRUISE6, STATION, STRATUM, SVVESSEL, YEAR,LAT, LON, DEPTH) %>% filter(n()>1) %>% mutate(num=n())

tdupmn=tdup %>% mutate_if(is.numeric, mean) %>% mutate_if(is.character, funs(paste(unique(.), collapse = "_"))) %>% slice(1)
test2=bind_rows(tunq, tdupmn) %>% ungroup() #reassign name
test2[with(test2, order("EST_TOWDATE", "YEAR", "MONTH", "STATION", "STRATUM", "LAT", "LON")),]
test2$EST_TOWDATE=as.character(test2$EST_TOWDATE)

test4=test2[complete.cases(test2$`74_Adt.y`),] # 15619 long, but NAs for zoo taxa...
test5=test4[complete.cases(test4$ctyp_100m3),] # 8183 long...not good
test6=test2[complete.cases(test2$ctyp_100m3),] # 14219
test7=test6[1:500,]
test8=test7 %>% select(EST_TOWDATE.x, EST_TOWDATE.y, STATION, STRATUM, LAT, LON, `74_Juv.x`, `74_Juv.y`, `74_Adt.x`, `74_Adt.y`, `74_ich`)
# why would some be different? NA values for stations
CoreOffshoreStrata<-c(seq(1010,1300,10),1340, seq(1360,1400,10),seq(1610,1760,10))
# inshore strata to use, still sampled by Bigelow
CoreInshore73to12=c(3020, 3050, 3080 ,3110 ,3140 ,3170, 3200, 3230, 3260, 3290, 3320, 3350 ,3380, 3410 ,3440)
# combine
strata_used=c(CoreOffshoreStrata,CoreInshore73to12)
sum(unique(test8$STRATUM[which(is.na(test8$`74_Juv.y`))]) %in% strata_used)

### check on station numbers, lat lon, est_towdate
load("/home/ryan/Downloads/survdat_lw.RData")# loads survdat.lw ; from Sean Lucey data into 2019
nsdat2=nsdat %>% filter(SVSPP==74, YEAR>1976)
nsdat3=survdat %>% filter(SVSPP==74, YEAR>1976)
sdat=survdat.lw %>% filter(SVSPP==74, YEAR>1976)
colnames(nsdat3)
colnanes(sdat)
nsdat3=ungroup(nsdat3)
class(nsdat3)
class(sdat)
sdat=as_tibble(sdat)
test=sdat%>% left_join(nsdat3, by =c("CRUISE6", "STATION", "STRATUM", "SVVESSEL", "YEAR","LAT", "LON", "DEPTH"))
test2=test[1:100,]

## clean up

cnsdat3=lapply(nsdat3, class)
csdat=lapply(sdat, class)

sdat$YEAR=as.integer(sdat$YEAR)
sdat$SEASON=as.character(sdat$SEASON)
sdat$SVVESSEL=as.character(sdat$SVVESSEL)
sdat$EST_TOWDATE=as.character(sdat$EST_TOWDATE)

# nsdat3[,c('YEAR', 'SVSPP', 'TOW', 'STRATUM', 'STATION', 'CRUISE6', 'CATCHSEX', 'NUMLEN', 'DEPTH')]=as.integer(nsdat3[,c('YEAR', 'SVSPP', 'TOW', 'STRATUM', 'STATION', 'CRUISE6', 'CATCHSEX', 'NUMLEN', 'DEPTH')])
# nsdat4=lapply(nsdat3[,c('YEAR', 'SVSPP', 'TOW', 'STRATUM', 'STATION', 'CRUISE6', 'CATCHSEX', 'NUMLEN', 'DEPTH')], as.integer)

nsdat3$YEAR=as.integer(nsdat3$YEAR)
nsdat3$SVSPP=as.integer(nsdat3$SVSPP)
nsdat3$TOW=as.integer(nsdat3$TOW)
nsdat3$STRATUM=as.integer(nsdat3$STRATUM)
nsdat3$STATION=as.integer(nsdat3$STATION)
nsdat3$CRUISE6=as.integer(nsdat3$CRUISE6)
nsdat3$CATCHSEX=as.integer(nsdat3$CATCHSEX)
nsdat3$NUMLEN=as.integer(nsdat3$NUMLEN)
nsdat3$DEPTH=as.integer(nsdat3$DEPTH)
nsdat3$SEASON=as.character(nsdat3$SEASON)
nsdat3$SVVESSEL=as.character(nsdat3$SVVESSEL)
nsdat3$EST_TOWDATE=as.character(nsdat3$EST_TOWDATE)
# nsdat3$SIZECAT=as.character(nsdat3$SIZECAT)

sdat$YEAR=as.integer(sdat$YEAR)
sdat$SVSPP=as.integer(sdat$SVSPP)
sdat$TOW=as.integer(sdat$TOW)
sdat$STRATUM=as.integer(sdat$STRATUM)
sdat$STATION=as.integer(sdat$STATION)
sdat$CRUISE6=as.integer(sdat$CRUISE6)
sdat$CATCHSEX=as.integer(sdat$CATCHSEX)
sdat$NUMLEN=as.integer(sdat$NUMLEN)
sdat$DEPTH=as.integer(sdat$DEPTH)
sdat$SEASON=as.character(sdat$SEASON)
sdat$SVVESSEL=as.character(sdat$SVVESSEL)
sdat$EST_TOWDATE=as.character(sdat$EST_TOWDATE)
sdat$SIZECAT=as.character(sdat$SIZECAT)

sdat[with(sdat, order("EST_TOWDATE", "YEAR", "MONTH", "STATION", "STRATUM", "LAT", "LON")),]
nsdat3[with(nsdat3, order("EST_TOWDATE", "YEAR", "MONTH", "STATION", "STRATUM", "LAT", "LON")),]


### now create dataframe with unique tows only, separated by stage
test=nsdat3[,c("CRUISE6","STATION","STRATUM","YEAR", "SVSPP", "stg")] # needs SVSPP...
svdtunq=nsdat3[!duplicated(test),]

sdat=survdat %>% filter(SVSPP==74, YEAR>1976)
test=sdat[,c("CRUISE6","STATION","STRATUM","YEAR", "SVSPP", "stg")] # needs SVSPP...
ogsvdtunq=sdat[!duplicated(test),]

test=svdtunq[1:200,]
test2=ogsvdtunq[1:200,]
test3=left_join(test, test2, by =c("CRUISE6", "STATION", "STRATUM", "YEAR","LAT", "LON", "DEPTH"))

test3=left_join(ogsvdtunq, svdtunq, by =c("CRUISE6", "STATION", "STRATUM", "YEAR","LAT", "LON", "DEPTH"))
test3=inner_join(ogsvdtunq, svdtunq, by =c("CRUISE6", "STATION", "STRATUM", "YEAR","LAT", "LON", "DEPTH"))
barplot(table(test3$YEAR[which(test3$SEASON.x=="SPRING")]))

## other computer
load('/home/ryan/Documents/Git/NEhabitat/Final_merged_fish_corABN_Zoo_Ich.Rda') # Abundance

### Select season for GAMs
slctseason="FALL"; fishseas="Fall"
slctseason="SPRING"; fishseas="Spr"

### Choose fish to run, adjust by searching on and changing " `74_ " to new ID **** ###
### Subset the formatted data for a single species (see Lmf)
# Lmf=read_excel('/home/ryan/Git/NEhabitat/Lm_included.xlsx')
# fishname='Cod' #73
fishname='Haddock' #74
# fishname='SilverHake' #72
# fishname='Pollock' #75

fish=FData.abn %>% dplyr::select(YEAR, SEASON, LAT:BOTTEMP, `74_Adt`, `74_Juv`, `74_ich`, volume_100m3:chl12)
fish$MONTH=month(FData.abn$EST_TOWDATE)
fish=fish[complete.cases(fish),]
fish$`74_ich`=ceiling(fish$`74_ich`) # make integer from numbers per 100 m/3
fish2=fish[which(fish$SEASON==slctseason),] # subset to season

### REMOVE replicate values from zooplankton for the entire data set
tunq=fish2 %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()==1) %>% mutate(num=n())
tdup=fish2 %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()>1) %>% mutate(num=n())
tdupmn=tdup %>% mutate_if(is.numeric, mean) %>% mutate_if(is.character, funs(paste(unique(.), collapse = "_"))) %>% slice(1)
fish2=bind_rows(tunq, tdupmn) %>% ungroup() #reassign name

# ### now split data into training and testing set (75% train 25% test, randomly chosen)
set.seed(101) # Set Seed so that same sample can be reproduced in future also
sample <- sample.int(n = nrow(fish2), size = floor(.75*nrow(fish2)), replace = F)
trainPA <- fish2[sample, ]
testPA  <- fish2[-sample, ]
###
### OR split by year (pre and post X)
bp=table(fish2$YEAR)
bp=scale(bp, FALSE, sum(bp))*100 # scale to percent samples
barplot(bp[,1], names.arg=unique(fish2$YEAR))
spltyr=2010
trainPA=fish2 %>% filter(YEAR <= spltyr)
testPA=fish2 %>% filter(YEAR > spltyr)
dim(trainPA)[1]/dim(fish2)[1]
dim(testPA)[1]/dim(fish2)[1]
###
### OR specific years to leave out (large recruitment, Temp, plus 2 others to get to 20 percent)
sum(bp[c(20,27,29,30, 32,34),1])
todrop=c(2003,2010,2012,2013,2015,2017)
trainPA=fish2 %>% filter(!(YEAR %in% todrop))
testPA=fish2 %>% filter(YEAR %in% todrop)
unique(testPA$YEAR)
unique(trainPA$YEAR)
dim(trainPA)[1]/dim(fish2)[1]
dim(testPA)[1]/dim(fish2)[1]


### pivot longer so stage is repeated
trainPA=trainPA %>% pivot_longer(c(`74_Adt`, `74_Juv`, `74_ich`), names_to = c("SVSPP", "Stg"), names_sep ="_", values_to = "Number")
trainPA$Stg=factor(trainPA$Stg, ordered=F)
logd=trainPA[,8:19]
logd=log10(logd+1)
trainPA[,8:19]=logd
trainPA$pa=ifelse(trainPA$Number > 0, 1, 0)
## now do likewise for testing data
testPA=testPA %>% pivot_longer(c(`74_Adt`, `74_Juv`, `74_ich`), names_to = c("SVSPP", "Stg"), names_sep ="_", values_to = "Number")
testPA$Stg=factor(testPA$Stg, ordered=F)
logd=testPA[,8:19]
logd=log10(logd+1)
testPA[,8:19]=logd
testPA$pa=ifelse(testPA$Number > 0, 1, 0)

#### IF USING ABUNDANCE RUN THIS AND STOP -
# (logbio is used in the HGAM code as the dependent variable name)
trainBIO=trainPA %>% filter(Number >0)
trainBIO$Abundance=trainBIO$Number
trainBIO$logbio=log10(trainBIO$Number)
testBIO=testPA %>% filter(Number >0)
testBIO$Abundance=testBIO$Number
testBIO$logbio=log10(testBIO$Number)

### UPDATE DATA to remove ichtyhyoplankton for "FALL HADDOCK" AND COD (only 6 positive samples 1977-2017)
trainPA=trainPA %>% filter(Stg!='ich') %>% droplevels(exclude="ich")
# levels(trainPA$Stg)
testPA=testPA %>% filter(Stg!='ich') %>% droplevels(exclude="ich")
trainBIO=trainBIO %>% filter(Stg!='ich') %>% droplevels(exclude="ich")
testBIO=testBIO %>% filter(Stg!='ich') %>% droplevels(exclude="ich")

###______________________________________________________________________-
### IF USING BIOMASS AS SECOND PART SKIP ABOVE SECTION AND RUN BELOW:
# Do likewise for biomass data, to be combined with PA later for final model
fish.b=FData.bio %>% dplyr::select(YEAR, SEASON, LAT:BOTTEMP, `74_Adt`, `74_Juv`, `74_ich`, volume_100m3:chl12)
fish.b$MONTH=month(FData.bio$EST_TOWDATE)
fish.b=fish.b[complete.cases(fish.b),]
fish.b$`74_ich`=ceiling(fish.b$`74_ich`) # make integer from numbers per 100 m/3
fish2.b=fish.b[which(fish.b$SEASON==slctseason),] # subset to season
### Transform numbers per 100_m3 to KG using below conversions
# For positive logbiomass model need to update the ichthyoplankton from (numbers per 100 m^3) to biomass
# From: Penglase et al. (2013) DOI: 10.7717/peerj.20
# Cod avg of experiments w/ good and poor diet: mean dry weight of fish 5-30days = 0.25 mg DW
# mg->kg = x* 1e-6
# https://www.int-res.com/articles/meps/32/m032p229.pdf -> 1 mg wet weight avg
fish.b$`74_ich`=fish.b$`74_ich`*1e-6

### REMOVE replicate values from zooplankton for the entire data set
tunq=fish2.b %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()==1) %>% mutate(num=n())
tdup=fish2.b %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()>1) %>% mutate(num=n())
tdupmn=tdup %>% mutate_if(is.numeric, mean) %>% mutate_if(is.character, funs(paste(unique(.), collapse = "_"))) %>% slice(1)
fish2.b=bind_rows(tunq, tdupmn) %>% ungroup() #reassign name


### Split data into training and testing set (75% train 25% test, randomly chosen)
set.seed(101) # Set Seed so that same sample can be reproduced in future also
sample <- sample.int(n = nrow(fish2.b), size = floor(.75*nrow(fish2.b)), replace = F)
trainBIO <- fish2.b[sample, ]
testBIO  <- fish2.b[-sample, ]
### OR split by year (pre and post X)
bp=table(fish2.b$YEAR)
bp=scale(bp, FALSE, sum(bp))*100 # scale to percent samples
barplot(bp[,1], names.arg=unique(fish2.b$YEAR))
spltyr=2010
trainBIO=fish2.b %>% filter(YEAR <= spltyr)
testBIO=fish2.b %>% filter(YEAR > spltyr)
dim(trainBIO)[1]/dim(fish2.b)[1]
dim(testBIO)[1]/dim(fish2.b)[1]
### OR specific years to leave out (large recruitment, Temp, plus 2 others to get to 20 percent)
sum(bp[c(20,27,29,30, 32,34),1])
todrop=c(2003,2010,2012,2013,2015,2017)
trainBIO=fish2.b %>% filter(!(YEAR %in% todrop))
testBIO=fish2.b %>% filter(YEAR %in% todrop)
unique(testBIO$YEAR)
unique(trainBIO$YEAR)
dim(trainBIO)[1]/dim(fish2.b)[1]
dim(testBIO)[1]/dim(fish2.b)[1]


### pivot longer so stage is repeated
trainBIO=trainBIO %>% pivot_longer(c(`74_Adt`, `74_Juv`, `74_ich`), names_to = c("SVSPP", "Stg"), names_sep ="_", values_to = "Biomass")
# trainBIO=trainBIO %>% pivot_longer(c(`74_Adt`, `74_Juv`, `74_ich`), names_to = c("SVSPP", "Stg"), names_sep ="_", values_to = "Abundance")

trainBIO$Stg=factor(trainBIO$Stg, ordered=F)
logd=trainBIO[,8:19]
logd=log10(logd+1)
trainBIO[,8:19]=logd
# IF USING BIOMASS
trainBIO=trainBIO %>%  
  filter(., Biomass >0) %>%
  mutate(., "logbio"=log10(Biomass+1))
## IF USING ABUNDANCE
# trainBIO=trainBIO %>%  
  # filter(., Abundance >0) %>%
  # mutate(., "logbio"=log10(Abundance))
## now do likewise for testing data
testBIO=testBIO %>% pivot_longer(c(`74_Adt`, `74_Juv`, `74_ich`), names_to = c("SVSPP", "Stg"), names_sep ="_", values_to = "Biomass")
# testBIO=testBIO %>% pivot_longer(c(`74_Adt`, `74_Juv`, `74_ich`), names_to = c("SVSPP", "Stg"), names_sep ="_", values_to = "Abundance")

testBIO$Stg=factor(testBIO$Stg, ordered=F)
logd=testBIO[,8:19]
logd=log10(logd+1)
testBIO[,8:19]=logd
# IF USING BIOMASS
testBIO=testBIO %>%  
  filter(., Biomass >0) %>%
  mutate(., "logbio"=log10(Biomass+1))
### IF USING ABUNDANCE ONLY
# testBIO=testBIO %>%  
  # filter(., Abundance >0) %>%
  # mutate(., "logbio"=log10(Abundance))
### verify data
test=trainPA[,c('YEAR', 'Stg', 'pa')]
traindatsum=test %>%  group_by(YEAR, Stg) %>%  add_count() %>% mutate('sum'=sum(pa)) %>% 
  select(-pa) %>% distinct() %>% arrange(.by_group=T)
traindatsum$pct=traindatsum$sum/traindatsum$n

barplot(table(trainPA$YEAR[which(trainPA$pa>0 & trainPA$Stg=='ich')]), main=paste(slctseason, ' Ich ', fishname))
barplot(table(trainPA$YEAR[which(trainPA$pa>0 & trainPA$Stg=='Adt')]), main=paste(slctseason, ' Adt',fishname))
barplot(table(trainPA$YEAR[which(trainPA$pa>0 & trainPA$Stg=='Juv')]), main=paste(slctseason, ' Juv',fishname))

barplot(traindatsum$pct[which(traindatsum$Stg=='Adt')], names.arg=traindatsum$YEAR[which(traindatsum$Stg=='Adt')],main=paste(slctseason, ' Adt',fishname))
barplot(traindatsum$sum[which(traindatsum$Stg=='Adt')], names.arg=traindatsum$YEAR[which(traindatsum$Stg=='Adt')],main=paste(slctseason, ' Adt',fishname))
barplot(traindatsum$pct[which(traindatsum$Stg=='Juv')], names.arg=traindatsum$YEAR[which(traindatsum$Stg=='Juv')],main=paste(slctseason, ' Juv',fishname))
barplot(traindatsum$sum[which(traindatsum$Stg=='Juv')], names.arg=traindatsum$YEAR[which(traindatsum$Stg=='Juv')],main=paste(slctseason, ' Juv',fishname))
barplot(traindatsum$pct[which(traindatsum$Stg=='ich')], names.arg=traindatsum$YEAR[which(traindatsum$Stg=='ich')],main=paste(slctseason, ' ich',fishname))
barplot(traindatsum$sum[which(traindatsum$Stg=='ich')], names.arg=traindatsum$YEAR[which(traindatsum$Stg=='ich')],main=paste(slctseason, ' ich',fishname))



