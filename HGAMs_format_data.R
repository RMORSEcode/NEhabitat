library(mgcv)
library(gratia)
library(car)
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

### Select season for GAMs
# slctseason="FALL"; fishseas="Fall"
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

# ### now split data into training and testing set (75% train 25% test, randomly chosen)
set.seed(101) # Set Seed so that same sample can be reproduced in future also
sample <- sample.int(n = nrow(fish2), size = floor(.75*nrow(fish2)), replace = F)
trainPA <- fish2[sample, ]
testPA  <- fish2[-sample, ]

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

### Do likewise for biomass data, to be combined with PA later for final model
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

### Split data into training and testing set (75% train 25% test, randomly chosen)
set.seed(101) # Set Seed so that same sample can be reproduced in future also
sample <- sample.int(n = nrow(fish2.b), size = floor(.75*nrow(fish2.b)), replace = F)
trainBIO <- fish2.b[sample, ]
testBIO  <- fish2.b[-sample, ]
### pivot longer so stage is repeated
trainBIO=trainBIO %>% pivot_longer(c(`74_Adt`, `74_Juv`, `74_ich`), names_to = c("SVSPP", "Stg"), names_sep ="_", values_to = "Biomass")
trainBIO$Stg=factor(trainBIO$Stg, ordered=F)
logd=trainBIO[,8:19]
logd=log10(logd+1)
trainBIO[,8:19]=logd
trainBIO=trainBIO %>%  
  filter(., Biomass >0) %>%
  mutate(., "logbio"=log10(Biomass+1))
## now do likewise for testing data
testBIO=testBIO %>% pivot_longer(c(`74_Adt`, `74_Juv`, `74_ich`), names_to = c("SVSPP", "Stg"), names_sep ="_", values_to = "Biomass")
testBIO$Stg=factor(testBIO$Stg, ordered=F)
logd=testBIO[,8:19]
logd=log10(logd+1)
testBIO[,8:19]=logd
testBIO=testBIO %>%  
  filter(., Biomass >0) %>%
  mutate(., "logbio"=log10(Biomass+1))





