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


## LOAD ORIGINAL DATA
# ogFData.abn=loadRData('/home/ryan/Git/NEhabitat/Final_merged_fish_corABN_Zoo_Ich.Rda')
## LOAD NEW DATA
FData.abn=loadRData('/home/ryan/Git/NEhabitat/20210415_Final_merged_fish_corABN_Zoo_Ich.Rda')

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

fish=FData.abn %>% dplyr::select(YEAR, SEASON, LAT:BOTTEMP, `74_Adt`, `74_Juv`, `74_ich`, volume_100m3:chl10)

fish$MONTH=month(FData.abn$EST_TOWDATE)
fish=fish[complete.cases(fish),]
fish$`74_ich`=ceiling(fish$`74_ich`) # make integer from numbers per 100 m/3
fish2=fish[which(fish$SEASON==slctseason),] # subset to season

### REMOVE replicate values from zooplankton for the entire data set
tunq=fish2 %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()==1) %>% mutate(num=n())
tdup=fish2 %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()>1) %>% mutate(num=n())
tdupmn=tdup %>% mutate_if(is.numeric, mean) %>% mutate_if(is.character, funs(paste(unique(.), collapse = "_"))) %>% slice(1)
fish2=bind_rows(tunq, tdupmn) %>% ungroup() #reassign name

# count non-zero values
colSums(fish2 !=0)

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



