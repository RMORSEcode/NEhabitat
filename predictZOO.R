# Choose
slctseason="FALL"; fishseas="Fall"; SEASON='Fall'
slctseason="SPRING"; fishseas="Spr"; SEASON='Spr'

# Testing GAMs on zooplankton abundance
zoo=FData.abn %>% dplyr::select(YEAR, SEASON, LAT:BOTSALIN, calfin_100m3:rug)
zoo$MONTH=month(FData.abn$EST_TOWDATE)
zoo=zoo[complete.cases(zoo),]
zoo2=zoo[which(zoo$SEASON==slctseason),] # subset to season
logd=zoo2[,10:19]
zoopa=ifelse(logd>0, 1, 0) # presence absence
logd10=log10(logd+1)
zoo2[,10:19]=logd10 # use in biomass (subset to only present)
zoo3=zoo2
zoo3[,10:19]=zoopa # use in PA only

set.seed(101) # Set Seed so that same sample can be reproduced in future also
# # Now Selecting 75% of data as sample from total 'n' rows of the data  
## presence absence only model
sample <- sample.int(n = nrow(zoo3), size = floor(.75*nrow(zoo3)), replace = F)
trainzoo <- zoo3[sample, ]
testzoo <- zoo3[-sample, ]

zoonm='calfin' #'pseudocal' #'calfin'

## positive biomass model - select zoo taxa to use for model and subset to positive biomass only
zooposbio=zoo2[which(zoo3$calfin_100m3>0),]
sample <- sample.int(n = nrow(zooposbio), size = floor(.75*nrow(zooposbio)), replace = F)
trainzoobio <- zooposbio[sample, ]
testzoobio <- zooposbio[-sample, ]

# + s(grnszmm, k=30, bs='tp') + s(rug, k=30, bs='tp')
zoo_modG_pa = gam(calfin_100m3 ~ s(BOTTEMP, k=35, bs="tp") +s(SURFTEMP, k=35, bs='tp') + s(DEPTH, k=30, bs="tp") + s(BOTSALIN, k=30, bs='tp') + s(SURFSALIN, k=30, bs="tp"), data=trainzoo, method = "REML", family="binomial", select=T)
zoo_modG_pb = gam(calfin_100m3 ~ s(BOTTEMP, k=35, bs="tp") +s(SURFTEMP, k=35, bs='tp') + s(DEPTH, k=30, bs="tp") + s(BOTSALIN, k=30, bs='tp') + s(SURFSALIN, k=30, bs="tp"), data=trainzoobio, method = "REML", family="gaussian", select=T)

save(zoo_modG_pa, file=paste('/home/ryan/Git/NEhabitat/rasters/',fishseas,'/zoop/zoo_modG_pa_',fishseas,'_',zoonm,'.Rdata',sep=''))
save(zoo_modG_pb, file=paste('/home/ryan/Git/NEhabitat/rasters/',fishseas,'/zoop/zoo_modG_pb_',fishseas,'_',zoonm,'.Rdata',sep='')) 
# ______
zoo_modG4_pa = gam(calfin_100m3 ~ s(LON, LAT, bs='tp') + s(BOTTEMP, k=35, bs="tp") +s(SURFTEMP, k=35, bs='tp') + s(DEPTH, k=30, bs="tp") + s(BOTSALIN, k=30, bs='tp') + s(SURFSALIN, k=30, bs="tp"), data=trainzoo, method = "REML", family="binomial", select = T)

zoo_modG4_pb = gam(calfin_100m3 ~ s(LON, LAT, bs='tp') + s(BOTTEMP, k=35, bs="tp") +s(SURFTEMP, k=35, bs='tp') + s(DEPTH, k=30, bs="tp") + s(BOTSALIN, k=30, bs='tp') + s(SURFSALIN, k=30, bs="tp"), data=trainzoobio, method = "REML", family="gaussian", select=T)

save(zoo_modG4_pa, file=paste('/home/ryan/Git/NEhabitat/rasters/',fishseas,'/zoop/zoo_modG4_pa_',fishseas,'_',zoonm,'.Rdata',sep=''))
save(zoo_modG4_pb, file=paste('/home/ryan/Git/NEhabitat/rasters/',fishseas,'/zoop/zoo_modG4_pb_',fishseas,'_',zoonm,'.Rdata',sep='')) 
# ______
# run soap film smoother below with Lat Lon first
trainzoopts=SpatialPoints(cbind(trainzoo$LON, trainzoo$LAT),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
trainzoo.in=sp::over(trainzoopts, sps)
trainzoo.sub=trainzoo
trainzoo.sub=trainzoo.sub[complete.cases(trainzoo.in),]

trainzoobiopts=SpatialPoints(cbind(trainzoobio$LON, trainzoobio$LAT),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
trainzoo.inbio=sp::over(trainzoobiopts, sps)
trainzoo.subbio=trainzoobio
trainzoo.subbio=trainzoo.subbio[complete.cases(trainzoo.inbio),]

zoo_modG4bll_pa = gam(calfin_100m3 ~ s(LON, LAT, bs = "so", xt = list(bnd = boundll)) + s(BOTTEMP, k=20, bs="tp") +s(SURFTEMP, k=20, bs='tp') + s(DEPTH, k=10, bs="tp") + s(BOTSALIN, k=30, bs='tp') + s(SURFSALIN, k=30, bs="tp"), data=trainzoo.sub, method = "REML", knots=knotsll, family="binomial", select=T)

zoo_modG4bll_pb = gam(calfin_100m3 ~ s(LON, LAT, bs = "so", xt = list(bnd = boundll)) + s(BOTTEMP, k=20, bs="tp") +s(SURFTEMP, k=20, bs='tp') + s(DEPTH, k=10, bs="tp") + s(BOTSALIN, k=30, bs='tp') + s(SURFSALIN, k=30, bs="tp"), data=trainzoo.subbio, method = "REML", knots=knotsll, family="gaussian", select=T)

save(zoo_modG4bll_pa, file=paste('/home/ryan/Git/NEhabitat/rasters/',fishseas,'/zoop/zoo_modG4bll_pa_',fishseas,'_',zoonm,'.Rdata',sep=''))
save(zoo_modG4bll_pb, file=paste('/home/ryan/Git/NEhabitat/rasters/',fishseas,'/zoop/zoo_modG4bll_pb_',fishseas,'_',zoonm,'.Rdata',sep=''))

# draw(zoo_modG)
# gam.check(zoo_modG)
# plot(fish_modG, shade=T, rug=T, scale=0)
# AIC(fish_modG)

path1=paste('/home/ryan/Git/NEhabitat/rasters/',fishseas,'/zoop/',zoonm,'/',sep='')
#### VERIFY MODEL SPECIES AND SEASON ####
fishseas
zoonm
trainzoo$SEASON[1] #verify season
modlistpa=list.files(path1, pattern = '_pa_') # presence-absence models
modlistpb=list.files(path1, pattern = '_pb_') # positive biomass models

### verify model as seperate files
pdf(paste(path1, 'PAmodels_AUC.pdf', sep=''), height=4, width=6)
for (i in (1:length(modlistpa))){
  modchoice=i
  modlistpa[modchoice]
  usemodel=loadRData(paste(path1,modlistpa[modchoice], sep=''))
  pred.test=predict(usemodel,testPA,type='response')
  preds.obs=data.frame(pred.test=pred.test,testPA$pa)
  colnames(preds.obs)=c("predicted","observed")
  preds.obs2=preds.obs2[complete.cases(preds.obs$predicted),]
  plotROC(preds.obs2$observed,preds.obs2$predicted, colorize = TRUE, main=paste('All Stages ', modlistpa[i], sep=''))
}
dev.off()
pred.test=predict(usemodel,testzoo,type='response')
preds.obs=data.frame(pred.test=pred.test,testzoo$calfin_100m3)
colnames(preds.obs)=c("predicted","observed")
plotROC(preds.obs2$observed,preds.obs2$predicted, colorize = TRUE, main=paste('All Stages ', modlistpa[i], sep=''))

## list data files in each folder
btlist=list.files(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/BT2', sep=''))
stlist=list.files(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/ST2', sep=''))
zlist=list.files(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/pseudo', sep=''))
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
tz=strsplit(zlist, split=('RAST_NESREG_'))
ttz=sapply(tz, function(x) strsplit(x, "[.]")[[2]][1], USE.NAMES=FALSE)
ttz2=as.numeric(ttz)
#Surface salinity
tss=strsplit(sslist, split=('RAST_NESREG_'))
ttss=sapply(tss, function(x) strsplit(x, "[.]")[[2]][1], USE.NAMES=FALSE)
ttss2=as.numeric(ttss)
#Bottom salinity
tbs=strsplit(bslist, split=('RAST_NESREG_'))
ttbs=sapply(tbs, function(x) strsplit(x, "[.]")[[2]][1], USE.NAMES=FALSE)
ttbs2=as.numeric(ttbs)

## DO for zooplankton species to model abundance
usemodelpa=zoo_modG_pa #loadRData(paste(path1,modlist[modchoice], sep='')) #fish_modS #_spr_had
usemodelbio=zoo_modG_pb
zooyrlist=seq(from=1992, to=2019, by=1)
fishnm='zoop'
zoosp='pseudocal' #'calfin' #'pseudocal'
rm(bt)
wd2=paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/', fishnm, '/', zoonm, '/', sep='')
### NOW loop over files, load yearly dynamic raster files and predict habitat from HGAM models
for (i in 1:length(zooyrlist)){
  bi=which(zooyrlist[i]==ttb2) # index of year
  bt=loadRData(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/BT2/', btlist[[bi]], sep=''))
  bi=which(zooyrlist[i]==tts2) # index of year
  st=loadRData(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/ST2/', stlist[[bi]], sep=''))
  bi=which(zooyrlist[i]==ttss2) # index of year
  ss=loadRData(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/SS2/', sslist[[bi]], sep=''))
  bi=which(zooyrlist[i]==ttbs2) # index of year
  bs=loadRData(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/BS2/', bslist[[bi]], sep=''))
  ef <- data.frame(coordinates(bt), val=values(bt))
  colnames(ef)=c("LON", "LAT", "BOTTEMP")
  ef$SURFTEMP=values(st)
  ef$DEPTH=values(gd4)
  ef$SURFSALIN=values(ss)
  ef$BOTSALIN=values(bs)
  ef$grnszmm=values(phi2)
  ef$rug=values(rug2)
  ef2=ef[complete.cases(ef),]
  # tt=ef2$Stg
  # tt2=fl[[1]][[1]][2][[1]][tt] # subsets to stage 
  # ef$Stg=factor(ef$Stg, levels=c("Adt", "Juv", "ich"))
  test1 <- predict.gam(usemodelpa, ef2, type='response')
  test2 <- predict.gam(usemodelbio, ef2, type='response')
  ef2$predpa=test1
  ef2$predbio=test2
  ef2$combinedout=test1*test2
  wd3=paste(zooyrlist[i], '_', SEASON, '_', zoosp, '_', '.RData', sep="")
  save(ef2, file=paste(wd2, wd3, sep=""))
  spg1=ef2[,c('LON', 'LAT', 'combinedout')]
  wd4=paste(zooyrlist[i], '_', 'RASTER', '_', SEASON, '_', zoosp, '_', '.RData', sep="")
  # tes1=rasterFromXYZ(spg1[complete.cases(spg1$Stg),])
  # save(tes1, file=paste(wd2,wd4, sep=''))
  coordinates(spg1)= ~ LON + LAT
  gridded(spg1)=T
  rastDF=raster(spg1)
  if (i == 1){
    keepstack=rastDF
  }
  # plot(rastDF)
  # save(rastDF, file=paste(wd2,wd4, sep=''))
  if (i >1){
    keepstack=stack(keepstack, rastDF)
  }
  save(keepstack, file=paste(wd2,'stacked_', SEASON, '_', zoosp, '_', '.RData', sep=""))
}
### save raster layers to new file:
t1=subset(keepstack, 26)
save(t1, file='/home/ryan/Git/NEhabitat/rasters/Fall/pseudo/RAST_NESREG_2017.04.01.07.TEMP.YEAR.000066596.RData')
t1=subset(keepstack, 27)
save(t1, file='/home/ryan/Git/NEhabitat/rasters/Fall/pseudo/RAST_NESREG_2018.04.01.07.TEMP.YEAR.000066596.RData')
t1=subset(keepstack, 28)
save(t1, file='/home/ryan/Git/NEhabitat/rasters/Fall/pseudo/RAST_NESREG_2019.04.01.07.TEMP.YEAR.000066596.RData')


### SAVE out model performance data to CSV file - deviance explained, AIC
modeval=data.frame(matrix(nrow=length(modlistpa), ncol=9, data=NA))
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
}
colnames(modeval)=c('model', 'PA.dev.exp','BIO.dev.exp','PA.aic','BIO.aic','PA.edf','BIO.edf','PA.res.df','BIO.res.df')
write.csv(modeval, file=paste(path1,'model_evaluation_', fishseas, '.csv', sep=""), row.names = F)


#### Save model hindcast output trends (mean, trend, variance)
## Load rasters
p1=paste('/home/ryan/Git/NEhabitat/rasters/',fishseas,'/', fishnm, '/', sep='')
p2='fish_modGSe_fall_haddock' #'fish_modG4_spr_Haddock/'
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

map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
points(trainPA$LON[which(trainPA$calfin_100m3>3 & trainPA$calfin_100m3<6)], trainPA$LAT[which(trainPA$calfin_100m3>3 & trainPA$calfin_100m3<6)], col=addTrans('red', 5), pch=19)
points(trainPA$LON[which(trainPA$calfin_100m3>1 & trainPA$calfin_100m3<3)], trainPA$LAT[which(trainPA$calfin_100m3>1 & trainPA$calfin_100m3<3)], col=addTrans('green', 5), pch=19)
points(trainPA$LON[which(trainPA$calfin_100m3<1)], trainPA$LAT[which(trainPA$calfin_100m3<1)], col=addTrans('blue', 5), pch=19)