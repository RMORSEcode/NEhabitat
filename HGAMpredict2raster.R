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
library(maps)
library(mapdata)
library(RColorBrewer)
library(colorspace)
library(marmap)


library(lubridate)
xdt=today()
xdt=gsub('-','',xdt)

### set up year list to match files with
yrlist=seq(from=1977, to=2019, by=1)

### These are the files to loop on for years
## load example raster data for spring 1977, use in predict mode for gam
# load('/home/ryan/Git/NEhabitat/rasters/test/RAST_NESREG_1977.04.03.BT.TEMP.YEAR.000066596.RData') #BT
# bt=masked.raster
# load('/home/ryan/Git/NEhabitat/rasters/test/RAST_NESREG_1977.04.03.ST.TEMP.YEAR.000066596.RData') #ST
# st=masked.raster
# load('/home/ryan/Git/NEhabitat/rasters/test/RAST_NESREG_1977.04.01.07.TEMP.YEAR.000066596.RData') #pseudo
# pse=masked.raster

### List of finished models
# fish_modG4_fall_had, fish_modG4bll_fall_had, 
# fish_modG_spr_had, fish_modG4_spr_had, fish_modG4bll_spr_had 

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

### SELECT SPR of FALL seasons to process
SEASON='Spr' # Fall
SEASON='Fall'

### NAME OF FISH
# fishnm='Haddock' #74
# fishnm='SilverHake' #72
fishnm='Cod' #73

## get path and list of models
path1=paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/', fishnm,'/', sep='') # Spr/Haddock'
# modlist=list.files(path1)

## CHOOSE MODEL AND VERIFY
# modlist
# modchoice=9
# modlist[modchoice]
# usemodel=loadRData(paste(path1,modlist[modchoice], sep='')) #fish_modS #_spr_had

#### VERIFY MODEL SPECIES AND SEASON ####
trainPA$SVSPP[1] #verify species
trainPA$SEASON[1] #verify season
modlistpa=list.files(path1, pattern = glob2rx("*_pa_*Rdata")) # presence-absence models
modlistpb=list.files(path1, pattern = glob2rx("*_pb_*Rdata")) # positive biomass models

## subset to newest models - (added date to sort)
modlistpa2=modlistpa[grepl("*2021*", modlistpa)]
modlistpb2=modlistpb[grepl("*2021*", modlistpb)]

modlistpa=modlistpa2
modlistpb=modlistpb2

## draw GAM smooths, save to pdf
# pdf(paste(path1, 'PAmodels_smooths.pdf', sep=''), height=8, width=12)
for (i in (1:length(modlistpa))){
  modchoice=i
  usemodel=loadRData(paste(path1,modlistpa[modchoice], sep=''))
  modnm=strsplit(modlistpa[modchoice], split ='.Rdata')[[1]]
  draw(usemodel)
  ggplot2::ggsave(paste(path1, modnm,'_smooths.pdf', sep=''), height=8, width=12)
}
# dev.off()

### Plot ROC curves
pdf(paste(path1, 'PAmodels_AUC_',xdt,'.pdf', sep=''), height=4, width=6)
for (i in (1:length(modlistpa))){
  modchoice=i
  modlistpa[modchoice]
  usemodel=loadRData(paste(path1,modlistpa[modchoice], sep=''))
  pred.test=predict(usemodel,testPA,type='response')  # make prediction for test set
  preds.obs=data.frame(pred.test=pred.test,testPA$pa) # data frame of preds and obs
  colnames(preds.obs)=c("predicted","observed")
  preds.obs2=preds.obs[complete.cases(preds.obs$predicted),]
  plotROC(preds.obs2$observed,preds.obs2$predicted, colorize = TRUE, main=paste('All Stages ', modlistpa[i], sep=''))
  sset=which(testPA$Stg=='Adt')
  preds.obs2=preds.obs[sset,]
  preds.obs2=preds.obs2[complete.cases(preds.obs2$predicted),]
  plotROC(preds.obs2$observed,preds.obs2$predicted, colorize = TRUE, main=paste('Adt only ',modlistpa[i], sep=''))
  sset=which(testPA$Stg=='Juv')
  preds.obs2=preds.obs[sset,]
  preds.obs2=preds.obs2[complete.cases(preds.obs2$predicted),]
  plotROC(preds.obs2$observed,preds.obs2$predicted, colorize = TRUE, main=paste('Juv only ',modlistpa[i], sep=''))
  sset=which(testPA$Stg=='ich')
  preds.obs2=preds.obs[sset,]
  preds.obs2=preds.obs2[complete.cases(preds.obs2$predicted),]
  plotROC(preds.obs2$observed,preds.obs2$predicted, colorize = TRUE, main=paste('Ich only ',modlistpa[i], sep=''))
}
dev.off()

# ### verify model as seperate files
# pdf(paste(path1, 'PAmodels_AUC_all_Stg.pdf', sep=''), height=4, width=6)
# for (i in (1:length(modlistpa))){
#   modchoice=i
#   modlistpa[modchoice]
#   usemodel=loadRData(paste(path1,modlistpa[modchoice], sep=''))
#   pred.test=predict(usemodel,testPA,type='response')
#   preds.obs=data.frame(pred.test=pred.test,testPA$pa)
#   colnames(preds.obs)=c("predicted","observed")
#   preds.obs2=preds.obs2[complete.cases(preds.obs$predicted),]
#   plotROC(preds.obs2$observed,preds.obs2$predicted, colorize = TRUE, main=paste('All Stages ', modlistpa[i], sep=''))
# }
# dev.off()

modauc=data.frame(matrix(nrow=length(modlistpa), ncol=2, data=NA))
pdf(paste(path1, 'PA_models_by_Stg_Adt_AUC.pdf', sep=''), height=4, width=6)
for (i in (1:length(modlistpa))){
  modchoice=i
  modlistpa[modchoice]
  usemodel=loadRData(paste(path1,modlistpa[modchoice], sep=''))
  pred.test=predict(usemodel,testPA,type='response')
  preds.obs=data.frame(pred.test=pred.test,testPA$pa)
  colnames(preds.obs)=c("predicted","observed")
  sset=which(testPA$Stg=='Adt')
  preds.obs2=preds.obs[sset,]
  preds.obs2=preds.obs2[complete.cases(preds.obs2$predicted),]
  modauc[i,1]=modlistpa[modchoice]
  modauc[i,2]=saveAUC(preds.obs2$observed,preds.obs2$predicted)[[1]]
  # plotROC(preds.obs2$observed,preds.obs2$predicted, colorize = TRUE, main=paste('Adt only ',modlistpa[i], sep=''))
}
colnames(modauc)=c('model', 'Adt.AUC')
write.csv(modauc, file=paste(path1,'Adt_model_AUC.csv', sep=""), row.names = F)
dev.off()

modauc=data.frame(matrix(nrow=length(modlistpa), ncol=2, data=NA))
pdf(paste(path1, 'PA_models_by_Stg_Juv_AUC.pdf', sep=''), height=4, width=6)
for (i in (1:length(modlistpa))){
  modchoice=i
  modlistpa[modchoice]
  usemodel=loadRData(paste(path1,modlistpa[modchoice], sep=''))
  pred.test=predict(usemodel,testPA,type='response')
  preds.obs=data.frame(pred.test=pred.test,testPA$pa)
  colnames(preds.obs)=c("predicted","observed")
  sset=which(testPA$Stg=='Juv')
  preds.obs2=preds.obs[sset,]
  preds.obs2=preds.obs2[complete.cases(preds.obs2$predicted),]
  modauc[i,1]=modlistpa[modchoice]
  modauc[i,2]=saveAUC(preds.obs2$observed,preds.obs2$predicted)[[1]]
  # plotROC(preds.obs2$observed,preds.obs2$predicted, colorize = TRUE, main=paste('Juv only ',modlistpa[i], sep=''))
}
colnames(modauc)=c('model', 'Juv.AUC')
write.csv(modauc, file=paste(path1,'Juv_model_AUC.csv', sep=""), row.names = F)
dev.off()

modauc=data.frame(matrix(nrow=length(modlistpa), ncol=2, data=NA))
pdf(paste(path1, 'PA_models_by_Stg_Ich_AUC.pdf', sep=''), height=4, width=6)
for (i in (1:length(modlistpa))){
  modchoice=i
  modlistpa[modchoice]
  usemodel=loadRData(paste(path1,modlistpa[modchoice], sep=''))
  pred.test=predict(usemodel,testPA,type='response')
  preds.obs=data.frame(pred.test=pred.test,testPA$pa)
  colnames(preds.obs)=c("predicted","observed")
  sset=which(testPA$Stg=='ich')
  preds.obs2=preds.obs[sset,]
  preds.obs2=preds.obs2[complete.cases(preds.obs2$predicted),]
  modauc[i,1]=modlistpa[modchoice]
  modauc[i,2]=saveAUC(preds.obs2$observed,preds.obs2$predicted)[[1]]
  # plotROC(preds.obs2$observed,preds.obs2$predicted, colorize = TRUE, main=paste('Ich only ',modlistpa[i], sep=''))
}
colnames(modauc)=c('model', 'Ich.AUC')
write.csv(modauc, file=paste(path1,'Ich_model_AUC.csv', sep=""), row.names = F)
dev.off()

## list data files in each folder
btlist=list.files(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/BT2', sep=''))
stlist=list.files(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/ST2', sep=''))
# zlist=list.files(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/pseudo', sep=''))
zlist=list.files(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/ctyp', sep=''))
sslist=list.files(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/SS2', sep=''))
bslist=list.files(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/BS2', sep=''))
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
#Surface salinity
# tss=strsplit(sslist, split=('RAST_NESREG_'))
# ttss=sapply(tss, function(x) strsplit(x, "[.]")[[2]][1], USE.NAMES=FALSE)
# ttss2=as.numeric(ttss)
#Bottom salinity
# tbs=strsplit(bslist, split=('RAST_NESREG_'))
# ttbs=sapply(tbs, function(x) strsplit(x, "[.]")[[2]][1], USE.NAMES=FALSE)
# ttbs2=as.numeric(ttbs)

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
## Remove bottom temp raster
rm(bt)
fl=levels=c("Adt", "Juv", "ich")
fl=levels=c("Adt", "Juv")
# fishnm='SilverHake'
wd2=paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/', fishnm, '/', sep='')
modchoice=2
modchoicepb=2
usemodel=loadRData(paste(path1,modlistpa[modchoice], sep=''))
usemodelbio=loadRData(paste(path1,modlistpb[modchoicepb], sep=''))
### NOW loop over files, load yearly dynamic raster files and predict habitat from HGAM models
for (jj in 1:length(fl)){
  for (i in 1:length(yrlist)){
    bi=which(yrlist[i]==ttb2) # index of year
    bt=loadRData(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/BT2/', btlist[[bi]], sep=''))
    bi=which(yrlist[i]==tts2) # index of year
    st=loadRData(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/ST2/', stlist[[bi]], sep=''))
    bi=which(yrlist[i]==ttz2) # index of year
    # pse=loadRData(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/pseudo/', zlist[[bi]], sep=''))
    # cty=loadRData(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/ctyp/', zlist[[bi]], sep=''))
    ef <- data.frame(coordinates(bt), val=values(bt))
    colnames(ef)=c("LON", "LAT", "BOTTEMP")
    ef$SURFTEMP=values(st)
    ef$DEPTH=values(gd4)
    ef$Stg=factor(fl[jj])#, levels=c('Adt', 'Juv', 'ich'))
    # ef$pseudo_100m3=values(pse)
    ef$ctyp_100m3=values(cty)
    # ef$chl1=values(chl.1)
    ef$chl2=values(chl.2)
    # ef$chl3=values(chl.3)
    ef$chl4=values(chl.4)
    # ef$chl5=values(chl.5)
    # ef$chl6=values(chl.6)
    # ef$chl7=values(chl.7)
    # ef$chl8=values(chl.8)
    # ef$chl9=values(chl.9)
    ef$chl10=values(chl.10)
    # ef$chl11=values(chl.11)
    # ef$chl12=values(chl.12)
    ef$grnszmm=values(phi2)
    ef$sand_pct=values(sand2)
    # ef$rug=values(rug2)
    ef2=ef
    # ef2=ef[complete.cases(ef),]
    # tt=ef2$Stg
    # tt2=fl[[1]][[1]][2][[1]][tt] # subsets to stage 
    # ef$Stg=factor(ef$Stg, levels=c("Adt", "Juv", "ich"))
    test1 <- predict.gam(usemodel, ef2, type='response')
    test2=predict.gam(usemodelbio, ef2, type='response')
    ef2$predpa=test1
    ef2$predbio=10^(test2) # biomass used was logbio -> log10(x+1)
    ef2$combinedout=ef2$predpa * ef2$predbio
    wd4=paste(yrlist[i], '_', 'matrix', '_', SEASON, '_', fishnm, '_',fl[jj], '.RData', sep="")
    # wd3=paste(zooyrlist[i], '_', SEASON, '_', zoosp, '_',fl[jj], '.RData', sep="")
    save(ef2, file=paste(wd2, wd4, sep=""))
    spg1=ef2[,c('LON', 'LAT', 'combinedout')]
    spg2=ef2[,c('LON','LAT','predpa')]
    # tes1=rasterFromXYZ(spg1[complete.cases(spg1$Stg),])
    # save(tes1, file=paste(wd2,wd4, sep=''))
    coordinates(spg1)= ~ LON + LAT
    coordinates(spg2)= ~ LON + LAT
    gridded(spg1)=T
    gridded(spg2)=T
    rastDF=raster(spg1)
    extent(rastDF)=c(-75.95, -65.45, 35.65, 44.65)
    rastDF2=raster(spg2)
    extent(rastDF2)=c(-75.95, -65.45, 35.65, 44.65)
    if (i == 1){
      keepstack=rastDF
      keepstackpa=rastDF2
    }
    # plot(rastDF)
    # save(rastDF, file=paste(wd2,wd4, sep=''))
    if (i >1){
      keepstack=stack(keepstack, rastDF)
      keepstackpa=stack(keepstackpa, rastDF2)
    }
  }
  save(keepstack, file=paste(wd2,'stacked_', SEASON, '_', fishnm, '_',fl[jj], '.RData', sep=""))
  save(keepstackpa, file=paste(wd2,'PA_only_stacked_', SEASON, '_', fishnm, '_',fl[jj], '.RData', sep=""))
}


# ## DO for zooplankton species to model abundance
# usemodelpa=zoo_modG_pa #loadRData(paste(path1,modlist[modchoice], sep='')) #fish_modS #_spr_had
# usemodelbio=zoo_modG_pb
# zooyrlist=seq(from=1992, to=2019, by=1)
# fishnm='zoop'
# zoosp='pseudocal' #'calfin' #'pseudocal'
# rm(bt)
# wd2=paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/', fishnm, '/', zoosp, '/', sep='')
# ### NOW loop over files, load yearly dynamic raster files and predict habitat from HGAM models
# for (i in 1:length(zooyrlist)){
#   bi=which(zooyrlist[i]==ttb2) # index of year
#   bt=loadRData(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/BT2/', btlist[[bi]], sep=''))
#   bi=which(zooyrlist[i]==tts2) # index of year
#   st=loadRData(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/ST2/', stlist[[bi]], sep=''))
#   bi=which(zooyrlist[i]==ttss2) # index of year
#   ss=loadRData(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/SS2/', sslist[[bi]], sep=''))
#   bi=which(zooyrlist[i]==ttbs2) # index of year
#   bs=loadRData(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/BS2/', bslist[[bi]], sep=''))
#   ef <- data.frame(coordinates(bt), val=values(bt))
#   colnames(ef)=c("LON", "LAT", "BOTTEMP")
#   ef$SURFTEMP=values(st)
#   ef$DEPTH=values(gd4)
#   ef$SURFSALIN=values(ss)
#   ef$BOTSALIN=values(bs)
#   ef$grnszmm=values(phi2)
#   ef$rug=values(rug2)
#   ef2=ef[complete.cases(ef),]
#   # tt=ef2$Stg
#   # tt2=fl[[1]][[1]][2][[1]][tt] # subsets to stage 
#   # ef$Stg=factor(ef$Stg, levels=c("Adt", "Juv", "ich"))
#   test1 <- predict.gam(usemodelpa, ef2, type='response')
#   test2 <- predict.gam(usemodelbio, ef2, type='response')
#   ef2$predpa=test1
#   ef2$predbio=test2
#   ef2$combinedout=test1*test2
#   wd3=paste(zooyrlist[i], '_', SEASON, '_', zoosp, '_', '.RData', sep="")
#   save(ef2, file=paste(wd2, wd3, sep=""))
#   spg1=ef2[,c('LON', 'LAT', 'combinedout')]
#   wd4=paste(zooyrlist[i], '_', 'RASTER', '_', SEASON, '_', zoosp, '_', '.RData', sep="")
#   # tes1=rasterFromXYZ(spg1[complete.cases(spg1$Stg),])
#   # save(tes1, file=paste(wd2,wd4, sep=''))
#   coordinates(spg1)= ~ LON + LAT
#   gridded(spg1)=T
#   rastDF=raster(spg1)
#   if (i == 1){
#     keepstack=rastDF
#   }
#   # plot(rastDF)
#   # save(rastDF, file=paste(wd2,wd4, sep=''))
#   if (i >1){
#     keepstack=stack(keepstack, rastDF)
#   }
#   save(keepstack, file=paste(wd2,'stacked_', SEASON, '_', zoosp, '_', '.RData', sep=""))
# }
### save raster layers to new file:
# t1=subset(keepstack, 26)
# save(t1, file='/home/ryan/Git/NEhabitat/rasters/Fall/pseudo/RAST_NESREG_2017.04.01.07.TEMP.YEAR.000066596.RData')
# t1=subset(keepstack, 27)
# save(t1, file='/home/ryan/Git/NEhabitat/rasters/Fall/pseudo/RAST_NESREG_2018.04.01.07.TEMP.YEAR.000066596.RData')
# t1=subset(keepstack, 28)
# save(t1, file='/home/ryan/Git/NEhabitat/rasters/Fall/pseudo/RAST_NESREG_2019.04.01.07.TEMP.YEAR.000066596.RData')

modlistpa[modchoice]




### SAVE out model performance data to CSV file - deviance explained, AIC
modeval=data.frame(matrix(nrow=length(modlistpa), ncol=11, data=NA))
for (i in 1:length(modlistpa)){
  modchoice=i
  usemodel=loadRData(paste(path1,modlistpa[modchoice], sep=''))
  usemodelbio=loadRData(paste(path1,modlistpb[modchoice], sep=''))
  modeval[i,1]=modlistpa[modchoice]
  modeval[i,2]=summary(usemodel)$dev.expl
  modeval[i,3]=summary(usemodelbio)$dev.expl
  modeval[i,4]=usemodel$aic
  modeval[i,5]=usemodelbio$aic
  modeval[i,6]=sum(usemodel$edf)
  modeval[i,7]=sum(usemodelbio$edf)
  modeval[i,8]=df.residual(usemodel)
  modeval[i,9]=df.residual(usemodelbio)
  modeval[i,10]=((df.residual(usemodel)+sum(usemodel$edf))/3)-sum(usemodel$edf)
  modeval[i,11]=((df.residual(usemodelbio)+sum(usemodelbio$edf))/3)-sum(usemodelbio$edf)
}
colnames(modeval)=c('model', 'PA.dev.exp','BIO.dev.exp','PA.aic','BIO.aic','PA.edf','BIO.edf','PA.res.df','BIO.res.df', 'PA.corr.res.df','BIO.corr.res.df')
write.csv(format(modeval, digits=2), file=paste(wd2,'model_evaluation_', SEASON, '_', fishnm, '_', '.csv', sep=""), row.names = F)

#### Save model hindcast output trends (mean, trend, variance)
## Load rasters
p1=paste('/home/ryan/Git/NEhabitat/rasters/',SEASON,'/', fishnm, '/', sep='')
p2='fish_modGI_Fall_cod' #fall_haddock' #'fish_modGI_spr_Haddock' 'fish_modGI_spr_Haddock_select_years_mod' 'fish_modGI_spr_Haddock_abundance'
p3=paste('/PA_only_stacked_', SEASON, '_', fishnm, '_', sep='') #'PA_only_stacked_Spr_Haddock_'
p4=paste('/stacked_', SEASON, '_', fishnm, '_', sep='') #'stacked_Spr_Haddock_'
ichpa=loadRData(paste(p1,p2,p3,'ich.RData', sep=''))
juvpa=loadRData(paste(p1,p2,p3,'Juv.RData', sep=''))
adtpa=loadRData(paste(p1,p2,p3,'Adt.RData', sep=''))
ich=loadRData(paste(p1,p2,p4,'ich.RData', sep=''))
juv=loadRData(paste(p1,p2,p4,'Juv.RData', sep=''))
adt=loadRData(paste(p1,p2,p4,'Adt.RData', sep=''))
### Run script and save as PDF
pdf(paste(path1, 'PA_Hindcast_',p2,'_Ich.pdf', sep=''), height=4, width=6)
plotRasterTrends(ichpa)
dev.off()
pdf(paste(path1, 'PA_Hindcast_',p2,'_Juv.pdf', sep=''), height=4, width=6)
plotRasterTrends(juvpa)
dev.off()
pdf(paste(path1, 'PA_Hindcast_',p2,'_Adt.pdf', sep=''), height=4, width=6)
plotRasterTrends(adtpa)
dev.off()
pdf(paste(path1, 'Bio_Hindcast_',p2,'_Ich.pdf', sep=''), height=4, width=6)
plotRasterTrends(ich)
dev.off()
pdf(paste(path1, 'Bio_Hindcast_',p2,'_Juv.pdf', sep=''), height=4, width=6)
plotRasterTrends(juv)
dev.off()
pdf(paste(path1, 'Bio_Hindcast_',p2,'_Adt.pdf', sep=''), height=4, width=6)
plotRasterTrends(adt)
dev.off()




### plot out raster stacks to visualize changes
pdf(paste(p1,p2,'/', 'Juv_hindcast.pdf', sep=''), height=4, width=6)
for (i in (1:length(yrlist))){
  plot(juvpa[[i]], zlim=c(0,1), col=viridis::viridis(64), main=paste(yrlist[i], SEASON, fishnm, sep=" "))
  }
dev.off()

pdf(paste(p1,p2,'/', 'Adt_hindcast.pdf', sep=''), height=4, width=6)
for (i in (1:length(yrlist))){
  plot(adtpa[[i]], zlim=c(0,1), col=viridis::viridis(64), main=paste(yrlist[i], SEASON, fishnm, sep=" "))
}
dev.off()

pdf(paste(p1,p2,'/', 'Ich_hindcast.pdf', sep=''), height=4, width=6)
for (i in (1:length(yrlist))){
  plot(ichpa[[i]], zlim=c(0,1), col=viridis::viridis(64), main=paste(yrlist[i], SEASON, fishnm, sep=" "))
}
dev.off()

## Biomass (or abundance) models
pdf(paste(p1,p2,'/', 'Adt_abund_hindcast.pdf', sep=''), height=4, width=6)
for (i in (1:length(yrlist))){
  plot(adt[[i]], zlim=c(0,65), main=paste(yrlist[i], SEASON, fishnm, sep=" "))
}
dev.off()
## Biomass (or abundance) models
pdf(paste(p1,p2,'/', 'Juv_abund_hindcast.pdf', sep=''), height=4, width=6)
for (i in (1:length(yrlist))){
  plot(juv[[i]], zlim=c(0,65), main=paste(yrlist[i], SEASON, fishnm, sep=" "))
}
dev.off()
## Biomass (or abundance) models
pdf(paste(p1,p2,'/', 'Ich_abund_hindcast.pdf', sep=''), height=4, width=6)
for (i in (1:length(yrlist))){
  plot(ich[[i]], zlim=c(0,65), main=paste(yrlist[i], SEASON, fishnm, sep=" "))
}
dev.off()


plot(keepstack[[28]], zlim=c(0,5), col=viridis::viridis(64))


plotrasterNES=function(lograx, mn, mx, titlex){
  # fun=function(x) { if (is.na(x[1])){ NA } else {10^x}}
  # rax=calc(lograx, fun)
  rax=lograx
  if (missing(mn)){
    mn=min(rax@data@values, na.rm = T)
  } else {
    mn=mn
  }
  if (missing(mx)){
    mx=max(rax@data@values, na.rm = T)
  } else {
    mx=mx
  }
  if (missing(titlex)){
    titlex=''
  }
br <- seq(mn, mx, by = (mx-mn)/9) 
# cl <- colorspace::diverge_hcl(length(br) - 1, power = 1)
cl=RColorBrewer::brewer.pal(9, 'YlOrRd')
# rng=range(rax[],na.rm=T)
# if (rng[1] < 0){
#   rng[1]=0
# }
rng=c(mn, mx)
arg=list(at=rng, labels=round(rng,2))
plot(rax, col=cl, breaks=br,axis.args=arg, las=1, main=paste(titlex)) 
}

pse=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/zoop/pseudocal/zoo_modG/stacked_Spr_pseudocal_ich.RData')
pse2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/zoop/pseudocal/zoo_modG4/stacked_Spr_pseudocal_ich.RData')
pse3=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/zoop/pseudocal/zoo_modG4bll/stacked_Spr_pseudocal_ich.RData')

pse4=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/pseudo/RAST_NESREG_2016.04.01.07.TEMP.YEAR.000066596.RData')

plotrasterNES(pse[[25]], mn=0, mx=5, titlex='New modG 2016')
plotrasterNES(pse3[[25]], mn=0, mx=5, titlex='modG4bll 2016')

plotrasterNES(pse2[[25]], mn=0, mx=5, titlex='modG4 2016')
plotrasterNES(pse4, mn=0, mx=5, titlex='Kevin 2016')

### checking for recruitment bumps in Ich output - 2013 and 2003 were big recruit years, need area index
plotrasterNES(ichpa[[20]], mn=0, mx=1, titlex=paste(yrlist[20]))
plotrasterNES(ichpa[[21]], mn=0, mx=1, titlex=paste(yrlist[21]))
plotrasterNES(ichpa[[22]], mn=0, mx=1, titlex=paste(yrlist[22]))
plotrasterNES(ichpa[[23]], mn=0, mx=1, titlex=paste(yrlist[23]))
plotrasterNES(ichpa[[24]], mn=0, mx=1, titlex=paste(yrlist[24]))
plotrasterNES(ichpa[[25]], mn=0, mx=1, titlex=paste(yrlist[25]))
plotrasterNES(ichpa[[26]], mn=0, mx=1, titlex=paste(yrlist[26]))
plotrasterNES(ichpa[[27]], mn=0, mx=1, titlex=paste(yrlist[27]))
plotrasterNES(ichpa[[28]], mn=0, mx=1, titlex=paste(yrlist[28]))
plotrasterNES(ichpa[[29]], mn=0, mx=1, titlex=paste(yrlist[29]))


plotrasterNES(ichpa[[34]], mn=0, mx=1, titlex=paste(yrlist[34]))
plotrasterNES(ichpa[[35]], mn=0, mx=1, titlex=paste(yrlist[35]))
plotrasterNES(ichpa[[36]], mn=0, mx=1, titlex=paste(yrlist[36]))
plotrasterNES(ichpa[[37]], mn=0, mx=1, titlex=paste(yrlist[37]))
plotrasterNES(ichpa[[38]], mn=0, mx=1, titlex=paste(yrlist[38]))

ichhab_gbk=raster::extract(ichpa, gbk, fun=mean, na.rm=T)
plot(ichhab_gbk[1,]~yrlist, type='l', las=1)
ichhab_gom=raster::extract(ichpa, gom, fun=mean, na.rm=T)
plot(ichhab_gom[1,]~yrlist, type='l', las=1)

adthab_gbk=raster::extract(adtpa, gbk, fun=mean, na.rm=T)
plot(adthab_gbk[1,]~yrlist, type='l', las=1)
adthab_gom=raster::extract(adtpa, gom, fun=mean, na.rm=T)
plot(adthab_gom[1,]~yrlist, type='l', las=1)

juvhab_gbk=raster::extract(juvpa, gbk, fun=mean, na.rm=T)
plot(juvhab_gbk[1,]~yrlist, type='l', las=1)
juvhab_gom=raster::extract(juvpa, gom, fun=mean, na.rm=T)
plot(juvhab_gom[1,]~yrlist, type='l', las=1)

### using trawl survey strata for haddock
pdf(paste(path1, 'PA_Hindcast_',p2,'_extracted_mean_habitat.pdf', sep=''), height=4, width=6)

ichhab_tsGBK=raster::extract(ichpa, tsGBK, fun=mean, na.rm=T)
plot(ichhab_tsGBK[1,]~yrlist, type='l', las=1, ylab='GBK ich habitat', xlab='')
ichhab_tsGOM=raster::extract(ichpa, tsGOM, fun=mean, na.rm=T)
plot(ichhab_tsGOM[1,]~yrlist, type='l', las=1, ylab='GOM ich habitat', xlab='')

adthab_tsGBK=raster::extract(adtpa, tsGBK, fun=mean, na.rm=T)
plot(adthab_tsGBK[1,]~yrlist, type='l', las=1, ylab='GBK adt habitat', xlab='')
adthab_tsGOM=raster::extract(adtpa, tsGOM, fun=mean, na.rm=T)
plot(adthab_tsGOM[1,]~yrlist, type='l', las=1, ylab='GOM adt habitat', xlab='')

juvhab_tsGBK=raster::extract(juvpa, tsGBK, fun=mean, na.rm=T)
plot(juvhab_tsGBK[1,]~yrlist, type='l', las=1, ylab='GBK juv habitat', xlab='')
juvhab_tsGOM=raster::extract(juvpa, tsGOM, fun=mean, na.rm=T)
plot(juvhab_tsGOM[1,]~yrlist, type='l', las=1, ylab='GOM juv habitat', xlab='')
dev.off()

### plot model residuals against variables
modchoice=5 # possible different lengths for PA and BIO (e.g. haddock)
modlistpb[modchoice]
usemodelbio=loadRData(paste(path1,modlistpb[modchoice], sep=''))
x=resid(usemodelbio, type = 'deviance')
test1 <- predict.gam(usemodelbio, trainBIO, type='response')
plot(x~test1)
plot(x~trainBIO$SURFTEMP)
plot(x~trainBIO$DEPTH)
plot(x~trainBIO$ctyp_100m3)
plot(x~trainBIO$grnszmm)
plot(x~trainBIO$chl2)
plot(x~trainBIO$chl10)
plot(x~trainBIO$Stg)
### for PA model, predict and residual is (obs-pred), using fish2(->fishtest) formatted like trainPA (long on stg, modified data)
## using `format_data_for_testing_HGAMS.R`
modchoice=5 # possible different lengths for PA and BIO (e.g. haddock)
modlistpa[modchoice]
usemodel=loadRData(paste(path1,modlistpa[modchoice], sep=''))
# x=resid(usemodel, type = 'deviance')
# test1 <- predict.gam(usemodel, trainPA, type='response')
test1 <- predict.gam(usemodel, fishtest, type='response')
fishtest$pred=test1
fishtest$resid=fishtest$pa-fishtest$pred
plot(fishtest$resid~fishtest$SURFTEMP, ylab='PA residual', main=modlistpa[modchoice])
plot(fishtest$resid~fishtest$DEPTH, ylab='PA residual', main=modlistpa[modchoice])
plot(fishtest$resid~fishtest$ctyp_100m3, ylab='PA residual', main=modlistpa[modchoice])
plot(fishtest$resid~fishtest$grnszmm, ylab='PA residual', main=modlistpa[modchoice])
plot(fishtest$resid~fishtest$chl2, ylab='PA residual', main=modlistpa[modchoice])
plot(fishtest$resid~fishtest$chl10, ylab='PA residual', main=modlistpa[modchoice])
plot(fishtest$resid~fishtest$Stg, ylab='PA residual', main=modlistpa[modchoice])
### PDF of residuals for each stage (rather than all together as above)
modused=strsplit(modlistpa[modchoice], '.Rdata')[[1]][1]
stage='Adt' #choose: 'Adt' 'Juv' 'ich'
pdf(paste(p1,stage,'_', modused,'_residuals.pdf', sep=''), height=4, width=6)
fishtest.stg=fishtest[which(fishtest$Stg==stage),]
plot(fishtest.stg$resid~fishtest.stg$SURFTEMP, ylab=paste(stage, 'PA residual'), main=modlistpa[modchoice])
plot(fishtest.stg$resid~fishtest.stg$DEPTH, ylab=paste(stage, 'PA residual'), main=modlistpa[modchoice])
plot(fishtest.stg$resid~fishtest.stg$ctyp_100m3, ylab=paste(stage, 'PA residual'), main=modlistpa[modchoice])
plot(fishtest.stg$resid~fishtest.stg$grnszmm, ylab=paste(stage, 'PA residual'), main=modlistpa[modchoice])
plot(fishtest.stg$resid~fishtest.stg$chl2, ylab=paste(stage, 'PA residual'), main=modlistpa[modchoice])
plot(fishtest.stg$resid~fishtest.stg$chl10, ylab=paste(stage, 'PA residual'), main=modlistpa[modchoice])
# plot(fishtest.stg$resid~fishtest.stg$Stg, ylab=paste(stage, 'PA residual'), main=modlistpa[modchoice])
dev.off()

show_palette <- function(colors) {
  image(1:n, 1, as.matrix(1:n), col = colors, 
        xlab = "", ylab = "", xaxt = "n", 
        yaxt = "n", bty = "n")
}

### plot residuals spatially using fishtest.stg (from above)
pdf(paste(p1,p2,'/', stage,'_PA_model_residuals.pdf', sep=''), height=4, width=6)
pdf(paste(p1,'/', stage,'Trawl_survey_catch.pdf', sep=''), height=4, width=6)

x2=fishtest.stg=fishtest[which(fishtest$Stg==stage),]
x=x2$resid
for(i in (1:length(yrlist))){
  j=yrlist[i]
  # par(mar = c(0,0,0,0))
  # par(oma = c(0,0,0,0))
  map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
  text(-74,44, j)
  map.axes(las=1)
  # rng=c(floor(min(x)), ceiling(max(x)))
  # myrng=modelr::seq_range(rng, by=1)
  # mylen=length(unique(round(x,1)))
  # myord=order(unique(round(x,1)))
  # mycolor=colorRampPalette(brewer.pal(length(myrng+1),"RdBu"))(mylen)
  # points(trainPA$LON[which(trainPA$YEAR==j)], trainPA$LAT[which(trainPA$YEAR==j)], col=mycolor[round(x,1)], pch=19)
  mn=floor(min(x))
  mx=ceiling(max(x))
  high=max(abs(mn), mx)
  br <- seq(-high, high, by = high/15) 
  cl=colorRampPalette(brewer.pal(9,"RdBu"))(length(br))
  # rng=c(floor(min(x)), ceiling(max(x)))
  # arg=list(at=rng, labels=round(rng,2))
  # points(x2$LON[which(x2$YEAR==j)], x2$LAT[which(x2$YEAR==j)], col=cl, pch=19)#, breaks=br,axis.args=arg,las=1)
  points(x2$LON[which(x2$YEAR==j)], x2$LAT[which(x2$YEAR==j)], col=ifelse(x2$pa[which(x2$YEAR==j)]==1, 'red', 'black'), pch=19)#, breaks=br,axis.args=arg,las=1)
  # points(trainBIO$LON[which(trainBIO$YEAR==j)], trainBIO$LAT[which(trainBIO$YEAR==j)], col=cl, pch=19)#, breaks=br,axis.args=arg,las=1)
  plot(nesbath,deep=-200, shallow=-200, step=1,add=T,lwd=1,col=addTrans('black',150),lty=1)
}
# show_palette(cl)
# text(1,2,'br')
z=matrix(1:length(br),nrow=1)
xx=1
y=br
image(xx,y,z,col=cl,axes=FALSE,xlab="",ylab="")
axis(2)
text(1,2, 'Deviance Residual')
dev.off()



### check on correlation between covariates
library(corrplot)
t=fishtest.stg[,5:35]
t=fish2[,c(5:7,12:38)]
t=fish2[,c(5:7,12:20,22:23,25,26,28,30,36)] #drop mud, keep chl 2,4,9 use for spring biomod all vars run

M=cor(t)
corrplot(M, method='circle')
t=fishtest.stg[,c(5,6,9,20,23,25,33)] # spring haddock
t=fishtest.stg[,c(5,6,7,20,24,25,32)] # fall haddock
M=cor(t)
corrplot(M, method='number', type='lower')



### haddock WG data for recruits and SSB
haddocksr=read.csv('/home/ryan/Downloads/SR.csv', header=T, stringsAsFactors = F)

## RSSB from Kevin (Peretti)
ssbr=read.csv('/home/ryan/1_habitat_analysis_2017/SSB and recruits.csv',header=T, stringsAsFactors = F)
ssbr2=read.csv('/home/ryan/1_habitat_analysis_2017/SSB and recruits.csv',header=T, skip=3, stringsAsFactors = F)
plot(ssbr2[10:44,1], ssbr2[10:44,13], type='l', las=1) #GOM cod rssb
plot(ssbr2[10:44,1], ssbr2[10:44,19], type='l', las=1) #GBK haddock rssb

# from Liz Brooks survey based RSSB data for GB haddock
lbhad=readxl::read_xlsx('/home/ryan/Desktop/Haddock/SSB.yr_and_R.xlsx', skip=4)
plot(lbhad$ssb.yr.Spring~lbhad$years, type='l', main='SSB')
plot(lbhad$Age0.Fall ~lbhad$years, type='l', main='Age 0 Fall')
lines(lbhad$`Age1.Spring(lagged)` ~lbhad$years, type='l', col='red')
lbhad$r0ssb=lbhad$Age0.Fall/lbhad$ssb.yr.Spring
lbhad$r1ssb=lbhad$`Age1.Spring(lagged)`/lbhad$ssb.yr.Spring

lbhad2=lbhad[complete.cases(lbhad),]
lbhad2$anom0=log10(lbhad2$r0ssb)-mean(log10(lbhad2$r0ssb))/sd(log10(lbhad2$r0ssb))
lbhad2$anom1=log10(lbhad2$r1ssb)-mean(log10(lbhad2$r1ssb))/sd(log10(lbhad2$r1ssb))

# log10(hdts$recr)-mean(log10(hdts$recr))/sd(log10(hdts$recr))
plot(log10(lbhad2$r0ssb)~lbhad2$years, type='l')
plot(lbhad2$anom0~lbhad2$years, type='l', main='age0 rssb'); abline(h=0, lty=2)
plot(lbhad2$anom1~lbhad2$years, type='l', main='age1 rssb'); abline(h=0, lty=2)

# trying dynamic time warping GBK RSSB vs Juv habitat
library(dtw)
juvts=(juvhab_gbk[1,2:36]-mean(juvhab_gbk[1,2:36]))/sd(juvhab_gbk[1,2:36])
adtts=(adthab_gbk[1,4:38]-mean(adthab_gbk[1,2:36]))/sd(adthab_gbk[1,2:36])
ichts=(ichhab_gbk[1,4:38]-mean(ichhab_gbk[1,2:36]))/sd(ichhab_gbk[1,2:36])
rssbts=ssbr2[10:44,19]

# cross correlation coefficient plots
ccf(ichts, rssbts)
ccf(juvts, rssbts)
ccf(adtts, rssbts)
# DTW plots
alignment<-dtw(ichts,rssbts,keep=TRUE)
plot(alignment,type="threeway", main='Ich')
dtwPlotDensity(alignment, main='Ich')
alignment<-dtw(juvts,rssbts,keep=TRUE)
plot(alignment,type="threeway", main='Juv')
dtwPlotDensity(alignment, main='Juv')
alignment<-dtw(adtts,rssbts,keep=TRUE)
plot(alignment,type="threeway", main='Adt')
dtwPlotDensity(alignment, main='Adt')

plot(dtw(ichts,rssbts,keep=TRUE,
      step=rabinerJuangStepPattern(6,"c")),
  type="twoway",offset=-2, main='Ich vs RSSB');
plot(dtw(juvts,rssbts,keep=TRUE,
         step=rabinerJuangStepPattern(6,"c")),
     type="twoway",offset=-2, main='Juv vs RSSB');
plot(dtw(adtts,rssbts,keep=TRUE,
         step=rabinerJuangStepPattern(6,"c")),
     type="twoway",offset=-2, main='Adt vs RSSB');


hdts=haddocksr[46:88,] #1977-2019
hdts$ra=log10(hdts$recr)-mean(log10(hdts$recr))/sd(log10(hdts$recr))
hdts$ssba=hdts$ssb-mean(hdts$ssb)/sd(hdts$ssb)
hdts$rssb=hdts$ra/hdts$ssba
# hdts=haddocksr[61:88,] #1992-2019
# plot(haddocksr$year, log(haddocksr$recr), type='l')
plot(hdts$year, log10(hdts$recr), type='b', ylab='log10 recruits', xlab='', pch=20) #log recruits
plot(hdts$year, hdts$recr/hdts$ssb, type='b', pch=20, ylab='recr/SSB', xlab='') # recruits per ssb
plot(hdts$year, log10(hdts$recr)/log10(hdts$ssb), type='l', pch=20, ylab='log10(R)/ log10(SSB)', xlab='', las=1) # recruits per ssb

plot(hdts$year, log10(hdts$recr-mean(hdts$recr))/log10(hdts$ssb-mean(hdts$ssb)), type='b', pch=20, ylab='log10(R)/ log10(SSB)', xlab='') # recruits per ssb

plot(hdts$year, hdts$ssb, type='b', pch=20, ylab='SSB', xlab='')
plot(hdts$year, log10(hdts$ssb), type='b', pch=20, ylab='log10 SSB', xlab='')

test=data.frame(log10(hdts$recr), ichhab_gbk[1,])
cor(test)

test=data.frame(ssbr2[10:44,13], juvhab_gom[1,4:38])
test=test[complete.cases(test),]
cor(test)


### load Kevin's RF model hindcast output for haddock biomass and PA 1976-2019 for comparison (drop 1976)
fbms=list.files('/home/ryan/1_habitat_analysis_2017/KEVIN_HADDOCK_2020/for ryan/',pattern = glob2rx("*BM.spri.*"))
fbmf=list.files('/home/ryan/1_habitat_analysis_2017/KEVIN_HADDOCK_2020/for ryan/',pattern = glob2rx("*BM.fall.*"))
fpas=list.files('/home/ryan/1_habitat_analysis_2017/KEVIN_HADDOCK_2020/for ryan/',pattern = glob2rx("*PA.spri.*"))
fpaf=list.files('/home/ryan/1_habitat_analysis_2017/KEVIN_HADDOCK_2020/for ryan/',pattern = glob2rx("*PA.fall.*"))
# RAST_haddoc_1978.00.00.BM.spri.0155.000000000.RData
tb=strsplit(fbms, split=('RAST_haddoc_'))
ttb=sapply(tb, function(x) strsplit(x, "[.]")[[2]][1], USE.NAMES=FALSE) # check order
wd='/home/ryan/1_habitat_analysis_2017/KEVIN_HADDOCK_2020/for ryan/'
i=2
kfbms=loadRData(paste(wd, fbms[i], sep=''))
for (i in 3:length(fbms)){
  r2=loadRData(paste(wd, fbms[i], sep=''))
  kfbms=stack(kfbms, r2)
}
i=2
kfpas=loadRData(paste(wd, fpas[i], sep=''))
for (i in 3:length(fpas)){
  r2=loadRData(paste(wd, fpas[i], sep=''))
  kfpas=stack(kfpas, r2)
}
i=2
kfbmf=loadRData(paste(wd, fbmf[i], sep=''))
for (i in 3:length(fbmf)){
  r2=loadRData(paste(wd, fbmf[i], sep=''))
  kfbmf=stack(kfbmf, r2)
}
i=2
kfpaf=loadRData(paste(wd, fpaf[i], sep=''))
for (i in 3:length(fpaf)){
  r2=loadRData(paste(wd, fpaf[i], sep=''))
  kfpaf=stack(kfpaf, r2)
}

#load spring stacked hindcasts select years only model
adtpa=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGI_spr_Haddock_select_years_mod/PA_only_stacked_Spr_Haddock_Adt.RData')
juvpa=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGI_spr_Haddock_select_years_mod/PA_only_stacked_Spr_Haddock_Juv.RData')
#load spring stacked hindcasts
adtpa=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGI_spr_Haddock/PA_only_stacked_Spr_Haddock_Adt.RData')
juvpa=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGI_spr_Haddock/PA_only_stacked_Spr_Haddock_Juv.RData')
#load fall stacked hindcasts
juvpa=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_fall_haddock/PA_only_stacked_Sstacked_Fall_Haddock_Juv.RData')
adtpa=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_fall_haddock/PA_only_stacked_Sstacked_Fall_Haddock_Adt.RData')
# r=stack(adtpa, juvpa) # may need to do it year by year
r=stack(adtpa[[1]], juvpa[[1]])
combr=(max(r))
for (i in 2:dim(adtpa)[[3]]){
r2=stack(adtpa[[i]], juvpa[[i]])
a2=(max(r2))
combr=stack(combr, a2)
}
kdiff=kfpas-combr

plotRasterMeanNegBinom=function(rastck){
  time <- 1:nlayers(rastck) 
  newrast.m=calc(rastck, fun=mean, na.rm=T)
  mn=cellStats(newrast.m, min)
  mx=cellStats(newrast.m, max)
  high=max(abs(mn), mx)
  br <- seq(-high, high, by = high/15) 
  cl=colorRampPalette(brewer.pal(9,"RdBu"))(length(br))
  rng=range(newrast.m[],na.rm=T)
  arg=list(at=rng, labels=round(rng,3))
  plot(newrast.m, col=cl, breaks=br,axis.args=arg,las=1, main=paste('mean')) # Yearly slope
  maps::map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="black", add=T)
}

i=1
plot(adtpa[[i]], zlim=c(0,1), main=paste('HGAM adt ', yrlist[i]))
plot(juvpa[[i]], zlim=c(0,1), main=paste('HGAM juv', yrlist[i]))
plot(combr[[i]], zlim=c(0,1), main=paste('HGAM comb max', yrlist[i]))
plot(kfpas[[i]], zlim=c(0,1), main=paste('RF', yrlist[i]))
plot(kdiff[[i]], zlim=c(-1,1), col=colorRampPalette(brewer.pal(9,"RdBu"))(30), main=paste('RF-HGAM comb', yrlist[i]))
plotRasterMeanNegBinom(kdiff)

### load and stack Bottom temperaure
btlist=list.files(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/BT2', sep=''))
wd3=paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/BT2', sep='')
rastBT=loadRData(paste(wd3,'/',btlist[1], sep=''))
for (i in 2:length(btlist)){
  rastDF=loadRData(paste(wd3,'/',btlist[i], sep=''))
  rastBT=stack(rastBT, rastDF)
}
pdf(paste(wd3, '/BT_trends.pdf', sep=''), height=4, width=6)
plotRasterTrends(rastBT)
dev.off()

BT_gbk=raster::extract(rastBT, gbk, fun=mean, na.rm=T)
plot(BT_gbk[1,]~yrlist, type='l')
BT_gom=raster::extract(rastBT, gom, fun=mean, na.rm=T)
plot(BT_gom[1,]~yrlist, type='l')

### load and stack Bottom temperaure
bslist=list.files(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/BS2', sep=''))
wd3=paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/BS2', sep='')
rastBS=loadRData(paste(wd3,'/',bslist[1], sep=''))
for (i in 2:length(bslist)){
  rastDF=loadRData(paste(wd3,'/',bslist[i], sep=''))
  rastBS=stack(rastBS, rastDF)
}
BS_gbk=raster::extract(rastBS, gbk, fun=mean, na.rm=T)
plot(BS_gbk[1,]~yrlist[16:43], type='l')
BT_gom=raster::extract(rastBT, gom, fun=mean, na.rm=T)
plot(BT_gom[1,]~yrlist, type='l')

### load and stack PSeudocalanus
zlist=list.files(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/pseudo', sep=''))
wd3=paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/pseudo', sep='')
rastPSE=loadRData(paste(wd3,'/',zlist[1], sep=''))
for (i in 2:length(btlist)){
  rastDF=loadRData(paste(wd3,'/',zlist[i], sep=''))
  rastPSE=stack(rastPSE, rastDF)
}
PSE_gbk=raster::extract(rastPSE, gbk, fun=mean, na.rm=T)
plot(PSE_gbk[1,]~yrlist, type='l')
PSE_gom=raster::extract(rastPSE, gom, fun=mean, na.rm=T)
plot(PSE_gom[1,]~yrlist, type='l')

### load and stack Ctyp
zlist=list.files(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/ctyp', sep=''))
wd3=paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/ctyp', sep='')
rastCTY=loadRData(paste(wd3,'/',zlist[1], sep=''))
for (i in 2:length(zlist)){
  rastDF=loadRData(paste(wd3,'/',zlist[i], sep=''))
  rastCTY=stack(rastCTY, rastDF)
}
CTY_gbk=raster::extract(rastCTY, gbk, fun=mean, na.rm=T)
plot(CTY_gbk[1,]~yrlist, type='l')
CTY_gom=raster::extract(rastCTY, gom, fun=mean, na.rm=T)
plot(CTY_gom[1,]~yrlist, type='l')

## check GAM pseudocal
PSE=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/zoop/pseudocal/zoo_modG/stacked_Spr_pseudocal_ich.RData')
PSE_gbk=raster::extract(PSE, gbk, fun=mean, na.rm=T)
plot(PSE_gbk[1,]~yrlist[16:43], type='l')
PSE_gom=raster::extract(PSE, gom, fun=mean, na.rm=T)
plot(PSE_gom[1,]~yrlist[16:43], type='l')

### check against EcoMon anomalies aggregated to EPUs
tt1=loadRData("/home/ryan/Downloads/zoo_abun_anom.rdata")
plot(tt1$value[which(tt1$variable=='Pse' & tt1$Region=='GB' & tt1$Time > 1991)] ~ 
       tt1$Time[which(tt1$variable=='Pse' & tt1$Region=='GB' & tt1$Time > 1991)], 
     type='b', ylab='Pse GB anom', xlab='')
abline(h=0, lty=3)

plot(tt1$value[which(tt1$variable=='Cty' & tt1$Region=='GB' & tt1$Time > 1991)] ~ 
       tt1$Time[which(tt1$variable=='Cty' & tt1$Region=='GB' & tt1$Time > 1991)], 
     type='b', ylab='Cty GB anom', xlab='')
abline(h=0, lty=3)

gbk=rgdal::readOGR('/home/ryan/Desktop/shapefiles/epu_shapes/EPU_GBKPoly.shp')
gom=rgdal::readOGR('/home/ryan/Desktop/shapefiles/epu_shapes/EPU_GOMPoly.shp')


### Load Bottom trawl survey strata - filter to offhsore areas N of Hatteras and subset to strata used for stocks
polygons <- sf::st_read('/home/ryan/Desktop/shapefiles/TrawlSurvey/strata.shp') %>% 
  dplyr::filter(STRATA > 0) %>% filter(SET_ == 1) %>%
  dplyr::mutate(STRATUM=stringr::str_sub(as.character(STRATA),start=2L,end=3L))
tsgb=c(13:25,29,30) #GBK strata for haddock
tsgom=c(26:28,36:40) #GOM strata for haddock
test2=polygons$STR2 %in% tsgb
tsGBK=polygons[test2,]
plot(tsGBK)
test2=polygons$STR2 %in% tsgom
tsGOM=polygons[test2,]
plot(tsGOM)


# CoreOffshoreStrata<-c(seq(1010,1300,10),1340, seq(1360,1400,10),seq(1610,1760,10))
# inshore strata to use, still sampled by Bigelow
# CoreInshore73to12=c(3020, 3050, 3080 ,3110 ,3140 ,3170, 3200, 3230, 3260, 3290, 3320, 3350 ,3380, 3410 ,3440)
# combine
# strata_used=c(CoreOffshoreStrata,CoreInshore73to12)
# test=tstratas$STRATA %in% strata_used
# subt=tstratas[test,]
# plot(subt)



NES5=rgdal::readOGR('/home/ryan/Desktop/nes_gbk_gome_gomw_mabn_mabsPoly.shp')


### Function to extract data using a shapefile
extract_calc=function(x, shp){
  v2=list()
  for(i in 1:dim(x)[3]){
    v=extract(x[[i]], shp)
    v1=lapply(v, function(xx) mean(xx, na.rm=T))
    v2[i]=list(v1)
  }
  m=matrix(unlist(v2), ncol=dim(x)[3], nrow=length(shp@polygons)) # box 0-29 =rows, years =cols
  # colnames(m)=seq(1998, 2016, by=1)
  # rownamse(m)=
  return(m)
}
### Create raster stacks and save them
# wd2=paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/', fishnm, '/', sep='')
# wd4=paste(yrlist[i], '_', 'RASTERpred', '_', SEASON, '_', fishnm, '_',fl[jj], '_', '.RData', sep="")
wd4=('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modG4bll_fall_had/Rasters')
wd4=('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modG4bll_spr_had/Rasters')

rastlist=list.files(wd4)
rl=strsplit(rastlist, split=('_'))
rlst=sapply(rl, "[",5) #subset to stage
rlyr=sapply(rl, "[",1) #subset to year
ryr=as.numeric(rlyr)

rastlistA=rastlist[which(rlst=='Adt')]
rastlistJ=rastlist[which(rlst=='Juv')]
rastlistI=rastlist[which(rlst=='ich')]

# # wd3='/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modG4bll_fall_had/Rasters/'
# load(paste(wd4,'/',rastlistA[1], sep=''))
# AdtHad=rastDF
# for (i in 2:length(rastlistA)){
#   load(paste(wd4,'/',rastlistA[i], sep=''))
#   AdtHad=stack(AdtHad, rastDF)
# }
# load(paste(wd4,'/',rastlistJ[1], sep=''))
# JuvHad=rastDF
# for (i in 2:length(rastlistA)){
#   load(paste(wd4,'/',rastlistJ[i], sep=''))
#   JuvHad=stack(JuvHad, rastDF)
# }
# load(paste(wd4,'/',rastlistI[1], sep=''))
# IchHad=rastDF
# for (i in 2:length(rastlistI)){
#   load(paste(wd4,'/',rastlistI[i], sep=''))
#   IchHad=stack(IchHad, rastDF)
# }
# names(IchHad)=yrlist
# names(AdtHad)=yrlist
# names(JuvHad)=yrlist
# save(IchHad, file=paste(wd4, 'tacks/IchHad.Rdata', sep='')) 
# save(AdtHad, file=paste(wd4, 'tacks/AdtHad.Rdata', sep=''))
# save(JuvHad, file=paste(wd4, 'tacks/JuvHad.Rdata', sep=''))

## just load rater stacks
AdtHad=loadRData("/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/stacked_Spr_Haddock_Adt.RData")
adtcod=loadRData("/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modI_spr_cod/stacked_Spr_Cod_Adt.RData")
juvcod=loadRData("/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modI_spr_cod/stacked_Spr_Cod_Juv.RData")

psefall=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/zoop/pseudocal/zoo_modG_fall_pse/stacked_Fall_pseudocal_.RData')
psespr=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/zoop/pseudocal/zoo_modG/stacked_Spr_pseudocal_ich.RData')
plotRasterTrends(psefall)
plotRasterTrends(psespr)

## this is not working, needs to call name of saved raster 'rastDF'
# loadNstack=function(x, namex){
#   namex=load(x[1])
#   for (i in 2:length(x)){
#     nameb=load(x[i])
#     namex=stack(namex, nameb)
#   }
#   return(namex)
# }

lisz=list.files('/home/ryan/1_habitat_analysis_2017/zoo_data/maps/spring_raster/calfin', pattern='RAST_NESREG_')
wd2='/home/ryan/1_habitat_analysis_2017/zoo_data/maps/spring_raster/calfin/'
for (i in (1:40)){
x=loadRData(paste(wd2, lisz[i], sep=''))
  if (i == 1){
    stk=x
  } else {
    stk=stack(stk, x)
  }
}
plot(cellStats(keepstack, 'mean'), type='b')
plot(cellStats(stk, 'mean'), type='b')

### see plot_map.R in {dropbox} for code
par(mar = c(0,0,0,0))
par(oma = c(0,0,0,0))
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(fish2$LON[which(fish2$`73_Adt`==0)], fish2$LAT[which(fish2$`73_Adt`==0)], col=addTrans('purple', 20), pch=19)
points(fish2$LON[which(fish2$`73_Adt`>0)], fish2$LAT[which(fish2$`73_Adt`>0)], col=addTrans('red', 20), pch=19)
plot(nesbath,deep=-200, shallow=-200, step=1,add=T,lwd=1,col=addTrans('black',150),lty=1)

map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(fish2$LON[which(fish2$`73_Juv`==0)], fish2$LAT[which(fish2$`73_Juv`==0)], col=addTrans('purple', 20), pch=19)
points(fish2$LON[which(fish2$`73_Juv`>0)], fish2$LAT[which(fish2$`73_Juv`>0)], col=addTrans('red', 20), pch=19)
plot(nesbath,deep=-200, shallow=-200, step=1,add=T,lwd=1,col=addTrans('black',150),lty=1)

map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(fish2$LON[which(fish2$`73_ich`==0)], fish2$LAT[which(fish2$`73_ich`==0)], col=addTrans('purple', 20), pch=19)
points(fish2$LON[which(fish2$`73_ich`>0)], fish2$LAT[which(fish2$`73_ich`>0)], col=addTrans('red', 20), pch=19)
plot(nesbath,deep=-200, shallow=-200, step=1,add=T,lwd=1,col=addTrans('black',150),lty=1)

map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(trainPA$LON[which(trainPA$pa==1 & trainPA$Stg=='Juv')], trainPA$LAT[which(trainPA$pa==1& trainPA$Stg=='Juv')], col=addTrans('blue', 5), pch=19)
plot(nesbath,deep=-200, shallow=-200, step=1,add=T,lwd=1,col=addTrans('black',150),lty=1)
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(trainPA$LON[which(trainPA$pa==1 & trainPA$Stg=='ich')], trainPA$LAT[which(trainPA$pa==1& trainPA$Stg=='ich')], col=addTrans('green', 255), pch=19)
plot(nesbath,deep=-200, shallow=-200, step=1,add=T,lwd=1,col=addTrans('black',150),lty=1)


map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(testPA$LON[which(testPA$pa==1 & testPA$Stg=='ich')], testPA$LAT[which(testPA$pa==1& testPA$Stg=='ich')], col=addTrans('green', 255), pch=19)
plot(nesbath,deep=-200, shallow=-200, step=1,add=T,lwd=1,col=addTrans('black',150),lty=1)

sbfish=fish2 %>% dplyr::filter(YEAR==2010)
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70", mar=c(4,0,0,0))
map.axes(las=1)
points(sbfish$LON[which(sbfish[,8]>0)], sbfish$LAT[which(sbfish[,8]>0)], col=addTrans('red', 255), pch=19)
points(sbfish$LON[which(sbfish[,8]==0)], sbfish$LAT[which(sbfish[,8]==0)], col=addTrans('black', 125), pch=19)
plot(nesbath,deep=-200, shallow=-200, step=1, add=T,lwd=1,col=addTrans('black',150),lty=1)
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70", mar=c(4,0,0,0))
map.axes(las=1)
points(sbfish$LON[which(sbfish[,9]>0)], sbfish$LAT[which(sbfish[,9]>0)], col=addTrans('red', 255), pch=19)
points(sbfish$LON[which(sbfish[,9]==0)], sbfish$LAT[which(sbfish[,9]==0)], col=addTrans('black', 125), pch=19)
plot(nesbath,deep=-200, shallow=-200, step=1, add=T,lwd=1,col=addTrans('black',150),lty=1)

##### _______________________________________________________________________________________________________
# Correlation routine
#####
source("/home/ryan/Desktop/8 R Functions/Rfuncs.R")
source("/home/ryan/Desktop/8 R Functions/Gridcorts.R")
source("/home/ryan/Desktop/8 R Functions/KDE_funcs.R")

addTrans <- function(color,trans)
{
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}
# cg1=gridcorts(rasterstack=AdtHad, method='kendall',type='corel') #spearman
# cg2=gridcorts(rasterstack=AdtHad, method='kendall',type='pval')
# cg3=cg1  
# plot(cg3, col=diverge_hcl(120),main=paste(spp,'\nCorrelation to', season, datalab), las=1)
# plot(nesbath,deep=-50, shallow=-50, step=1,add=T,lwd=1,col='gray40',lty=1)
# plot(nesbath,deep=-100, shallow=-100, step=1,add=T,lwd=1,col='gray60',lty=1)
# plot(nesbath,deep=-200, shallow=-200, step=1,add=T,lwd=1,col='gray60',lty=1)
# map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="black", add=T)
# 
# #mark significant correlations
# plot(cg3, col=diverge_hcl(120),main=paste(spp,'\nCorrelation to', season, datalab), las=1)
# cg3[cg2>0.05]=NA
# test=rasterToPoints(cg3)
# points(test, pch=12, col=addTrans('black',50), cex=0.5)
# plot(nesbath,deep=-50, shallow=-50, step=1,add=T,lwd=1,col='gray40',lty=1)
# plot(nesbath,deep=-100, shallow=-100, step=1,add=T,lwd=1,col='gray60',lty=1)
# plot(nesbath,deep=-200, shallow=-200, step=1,add=T,lwd=1,col='gray60',lty=1)
# map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="black", add=T)

#####
# Linear model of environmental variable change, get slope and stats
# drop first layer for spring wind stress
#####
# shp.dat3=crop(shp.dat2, small)

plotRasterTrends=function(rastck){
  time <- 1:nlayers(rastck) 
  newrast.m=calc(rastck, fun=mean, na.rm=T)
  mn=cellStats(newrast.m, min)
  mx=cellStats(newrast.m, max)
  high=max(0, mx)
  br <- seq(0, high, by = high/15) 
  cl=colorRampPalette(brewer.pal(9,"Reds"))(length(br))
  rng=range(newrast.m[],na.rm=T)
  arg=list(at=rng, labels=round(rng,3))
  plot(newrast.m, col=cl, breaks=br,axis.args=arg,las=1, main=paste(length(time),'yrs','\nmean distribution')) # Yearly slope
  maps::map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="black", add=T)
  fun=function(x) { if (is.na(x[1])){ NA } else { m = lm(x ~ time); summary(m)$coefficients[2] }}
  newrast=calc(rastck, fun)
  mn=cellStats(newrast, min)
  mx=cellStats(newrast, max)
  high=max(abs(mn), mx)
  br <- seq(mn, mx, by = high/15) 
  cl <- colorspace::diverge_hcl(length(br) - 1, power = 1) 
  rng=range(newrast[],na.rm=T)
  arg=list(at=rng, labels=round(rng,3))
  plot(newrast, col=cl, breaks=br,axis.args=arg,las=1, main=paste(length(time),'yrs','\nYearly Slope')) # Yearly slope
  maps::map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="black", add=T)
  # plot slope times length of time
  newrast.t=newrast*length(time)
  mn=cellStats(newrast.t, min) #min(newrast.t@data@values, na.rm = T)
  mx=cellStats(newrast.t, max) #max(newrast.t@data@values, na.rm = T)
  high=max(abs(mn), mx)
  br <- seq(-high, high, by = high/15) 
  cl <- colorspace::diverge_hcl(length(br) - 1, power = 1) 
  rng=range(newrast.t[],na.rm=T)
  arg=list(at=rng, labels=round(rng,3))
  plot(newrast.t, col=cl, breaks=br,axis.args=arg,las=1, main=paste(length(time),'yr change')) # Time series change
  maps::map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="black", add=T)
  ## Calc Significance
  fun=function(x) { if (is.na(x[1])){ NA } else { m = lm(x ~ time); summary(m)$coefficients[8] }}
  p <- calc(rastck, fun=fun)
  m = c(0, 0.05, 1, 0.05, 1, 0)
  rclmat = matrix(m, ncol=3, byrow=TRUE)
  p.mask = reclassify(p, rclmat)
  fun=function(x) { x[x<1] <- NA; return(x)}
  p.mask.NA = calc(p.mask, fun)
  trend.sig = mask(newrast, p.mask.NA)
  mn=cellStats(newrast, min)
  mx=cellStats(newrast, max)
  high=max(abs(mn), mx)
  br <- seq(-high, high, by = high/15) 
  cl <- colorspace::diverge_hcl(length(br) - 1, power = 1) 
  rng=range(newrast[],na.rm=T)
  arg=list(at=rng, labels=round(rng,3))
  plot(newrast, main=paste(length(time),'yrs','\nYearly Slope'),col=cl, breaks=br,axis.args=arg,las=1) # Yearly slope
  test=rasterToPoints(trend.sig)
  points(test, pch='+', col=addTrans('black',50), cex=0.5)
  maps::map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="black", add=T)
  ## Plot Variance
  fun=function(x) { if (is.na(x[1])){ NA } else { m = lm(x ~ time); (summary(m)$sigma)}} ; TITL='Std deviation of slope'
  newrast.v=calc(rastck, fun)
  mn=cellStats(newrast.v, min)
  mx=cellStats(newrast.v, max)
  high=max(abs(mn), mx)
  br <- seq(0, high, by = high/15) 
  cl=colorRampPalette(brewer.pal(9,"Reds"))(length(br))
  rng=c(0, mx) #range(newrast.v[],na.rm=T)
  arg=list(at=rng, labels=round(rng,3))
  plot(newrast.v, main=paste(length(time),'yrs','\n', TITL),col=cl, breaks=br,axis.args=arg,las=1) # 
  maps::map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="black", add=T)
}



time <- 1:nlayers(ichcod) 
fun=function(x) { if (is.na(x[1])){ NA } else { m = lm(x ~ time); summary(m)$coefficients[2] }}
shp.dat.slope=calc(juvcod, fun)
shp.dat.time=shp.dat.slope*length(time)
mn=cellStats(shp.dat.slope, min) #min(shp.dat.slope@data@values, na.rm = T)
mx=cellStats(shp.dat.slope, max) #max(shp.dat.slope@data@values, na.rm = T)
high=max(abs(mn), mx)
br <- seq(-high, high, by = high/15) 
cl <- colorspace::diverge_hcl(length(br) - 1, power = 1) 
rng=range(shp.dat.slope[],na.rm=T)
arg=list(at=rng, labels=round(rng,3))
# plot(shp.dat.slope, main=paste(datalab, season, length(time),'yrs','\nYearly Slope'),col=diverge_hcl(120), las=1) # Yearly slope
plot(shp.dat.slope, col=cl, breaks=br,axis.args=arg,las=1, main=paste(length(time),'yrs','\nYearly Slope')) # Yearly slope
mn=cellStats(shp.dat.slope, min) #min(shp.dat.time@data@values, na.rm = T)
mx=cellStats(shp.dat.slope, max) #max(shp.dat.time@data@values, na.rm = T)
high=max(abs(mn), mx)
br <- seq(-high, high, by = high/15) 
cl <- colorspace::diverge_hcl(length(br) - 1, power = 1) 
rng=range(shp.dat.time[],na.rm=T)
arg=list(at=rng, labels=round(rng,3))
plot(shp.dat.time, col=cl, breaks=br,axis.args=arg,las=1, main=paste(length(time),'yr change')) # Time series change
# plot(nesbath,deep=-50, shallow=-50, step=1,add=T,lwd=1,col='gray40',lty=1)
# plot(nesbath,deep=-200, shallow=-200, step=1,add=T,lwd=1,col='gray60',lty=1)
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="black", add=T)

## Calc Significance
fun=function(x) { if (is.na(x[1])){ NA } else { m = lm(x ~ time); summary(m)$coefficients[8] }}
p <- calc(pse, fun=fun)
m = c(0, 0.05, 1, 0.05, 1, 0)
rclmat = matrix(m, ncol=3, byrow=TRUE)
p.mask = reclassify(p, rclmat)
fun=function(x) { x[x<1] <- NA; return(x)}
p.mask.NA = calc(p.mask, fun)
trend.sig = mask(shp.dat.slope, p.mask.NA)
mn=min(shp.dat.slope@data@values, na.rm = T)
mx=max(shp.dat.slope@data@values, na.rm = T)
high=max(abs(mn), mx)
br <- seq(-high, high, by = high/15) 
cl <- colorspace::diverge_hcl(length(br) - 1, power = 1) 
rng=range(shp.dat.slope[],na.rm=T)
arg=list(at=rng, labels=round(rng,3))
plot(shp.dat.slope, main=paste(length(time),'yrs','\nYearly Slope'),col=cl, breaks=br,axis.args=arg,las=1) # Yearly slope
# plot(shp.dat.slope, main=paste(datalab, season, '\nTrend'), col=diverge_hcl(120), las=1)
test=rasterToPoints(trend.sig)
points(test, pch='+', col=addTrans('black',50), cex=0.5)
plot(nesbath,deep=-50, shallow=-50, step=1,add=T,lwd=1,col='gray40',lty=1)
plot(nesbath,deep=-200, shallow=-200, step=1,add=T,lwd=1,col='gray60',lty=1)
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="black", add=T)

## Plot Variance
# fun=function(x) { if (is.na(x[1])){ NA } else { m = lm(x ~ time); (summary(m)$coefficients[2,2])}} ;TITL='Standard error of slope'# (coefficient)
fun=function(x) { if (is.na(x[1])){ NA } else { m = lm(x ~ time); (summary(m)$sigma)}} ; TITL='Std deviation of slope'# square root of estimated variance of random error
shp.dat.var=calc(shp.dat3, fun)
mn=min(shp.dat.var@data@values, na.rm = T)
mx=max(shp.dat.var@data@values, na.rm = T)
high=max(abs(mn), mx)
magn=floor(log10(mx))
br <- seq(0, high, by = high/15) 
# cl <- diverge_hcl(length(br) - 1, power = 1) 
cl=colorRampPalette(brewer.pal(9,"Reds"))(length(br))
rng=range(shp.dat.var[],na.rm=T)
arg=list(at=rng, labels=round(rng,abs(min(log10(rng)))))
plot(shp.dat.var, main=paste(datalab, season, length(time),'yrs','\n', TITL),col=cl, breaks=br,axis.args=arg,las=1) # 
plot(nesbath,deep=-50, shallow=-50, step=1,add=T,lwd=1,col='gray40',lty=1)
plot(nesbath,deep=-200, shallow=-200, step=1,add=T,lwd=1,col='gray60',lty=1)
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="black", add=T)
