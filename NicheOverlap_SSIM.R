# compute the ssmi statistic for habitat similarity for species SEASONAL pairs
#######################################################################
library (raster)
library(SPUTNIK)
#library(SDMTools)

# as of 10/2020 sputnik need edgeR and it was not in CRAN
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")

library("BiocManager")
install("edgeR")

# since less species in spring get the species list here
# spring rf occ data
setwd("C:/3_rf_v5/output/spring")
modfiles = list.files(pattern="bm_mod_rf_")
species  = substr(modfiles,11,16)
remove(modfiles)

dataoutfile = "C:/3_rf_v5/output/seasonal_ssim.csv"

for (i in 1:length(species)){
  print(i)
  # get spring files
  spr_dir = "C:/3_rf_v5/output/spring"
  files.spr = intersect(list.files(path=spr_dir,pattern = "RAST"),
                        list.files(path=spr_dir,pattern = species[i]))
  files.spr = files.spr[grep(".PA.", files.spr)]  
  
  # get fall files
  fall_dir = "C:/3_rf_v5/output/fall"
  files.fall = intersect(list.files(path=fall_dir,pattern = "RAST"),
                         list.files(path=fall_dir,pattern = species[i]))
  files.fall = files.fall[grep(".PA.", files.fall)]  
  

for(j in 1:length(files.spr)){
  load(as.character(paste0(spr_dir,"/",files.spr[j])))
  mr1 = masked.raster
  load(as.character(paste0(fall_dir,"/",files.fall[j])))
  mr2 = masked.raster
  
  
  r1m = rasterToPoints(mr1)
  r2m = rasterToPoints(mr2)
  
  ssim = SSIM(r1m[,3], r2m[,3], numBreaks = 256)
  
  year  = substr(files.spr[j],13,16)
  dataline <- matrix(c(species[i], year, ssim),1,3)
  write.table(dataline,file=dataoutfile,sep=",",row.name=F,col.names=F,append=TRUE)  
  
  } # end j
} # end i

## niche overlap index (similar to I-statistic from SDMTools)
library(dismo)
# use on rasters x and y 
# nicheOverlap(x, y, stat='I', mask=TRUE, checkNegatives=TRUE)

# ## Load Cod spring GSe raster stacks of hindcasts
# ich1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_spr_cod/PA_only_stacked_Spr_Cod_ich.RData')
# juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_spr_cod/PA_only_stacked_Spr_Cod_Juv.RData')
# adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_spr_cod/PA_only_stacked_Spr_Cod_Adt.RData')  
# ## Load Haddock spring GI raster stacks of hindcasts
# ich2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGI_spr_Haddock/PA_only_stacked_Spr_Haddock_ich.RData')
# juv2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGI_spr_Haddock/PA_only_stacked_Spr_Haddock_Juv.RData')
# adt2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGI_spr_Haddock/PA_only_stacked_Spr_Haddock_Adt.RData')

## 20210317 add new models
#spring GSe haddock
# ich1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_haddock/PA_only_stacked_Spr_Haddock_ich.RData')
# juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_haddock/PA_only_stacked_Spr_Haddock_Juv.RData')
# adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_haddock/PA_only_stacked_Spr_Haddock_Adt.RData')
# #spring GSe cod
# ich2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod/PA_only_stacked_Spr_Cod_ich.RData')
# juv2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod/PA_only_stacked_Spr_Cod_Juv.RData')
# adt2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod/PA_only_stacked_Spr_Cod_Adt.RData')
# #spring GI haddock
# ich1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGI_Spr_haddock/PA_only_stacked_Spr_Haddock_ich.RData')
# juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGI_Spr_haddock/PA_only_stacked_Spr_Haddock_Juv.RData')
# adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGI_Spr_haddock/PA_only_stacked_Spr_Haddock_Adt.RData')
# #spring GI cod
# ich2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGI_Spr_cod/PA_only_stacked_Spr_Cod_ich.RData')
# juv2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGI_Spr_cod/PA_only_stacked_Spr_Cod_Juv.RData')
# adt2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGI_Spr_cod/PA_only_stacked_Spr_Cod_Adt.RData')

#### Niche Overlap ####_________________________________________________________________________________
### 1) Seasonal Inter-species habitat overlap: Cod and Haddock by season, stage ###

# 20210421 add new models for haddock
#spring GSe haddock
ich1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/PA_only_stacked_Spr_Haddock_ich.RData')
juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/PA_only_stacked_Spr_Haddock_Juv.RData')
adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/PA_only_stacked_Spr_Haddock_Adt.RData')
#spring GSe cod
# ich2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod/PA_only_stacked_Spr_Cod_ich.RData')
# juv2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod/PA_only_stacked_Spr_Cod_Juv.RData')
# adt2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod/PA_only_stacked_Spr_Cod_Adt.RData')
#spring GSe cod updated from above (new data from Sean survdat)
ich2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod_20210518/PA_only_stacked_Spr_Cod_ich.RData')
juv2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod_20210518/PA_only_stacked_Spr_Cod_Juv.RData')
adt2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod_20210518/PA_only_stacked_Spr_Cod_Adt.RData')

adtno=data.frame(matrix(NA, ncol=1, nrow=dim(adt1)[3])); colnames(adtno)='I'
for( i in 1:dim(adt1)[3]){
adtno[i,1]=nicheOverlap(adt1[[i]], adt2[[i]], stat='I', mask=TRUE, checkNegatives=TRUE)
}
tadt=as.matrix(unlist(adtno))
plot(tadt~yrlist, type='b',ylim=c(.75,1), ylab='adult niche overlap', xlab='')

juvno=data.frame(matrix(NA, ncol=1, nrow=dim(juv1)[3])); colnames(juvno)='I'
for( i in 1:dim(juv1)[3]){
  juvno[i,1]=nicheOverlap(juv1[[i]], juv2[[i]], stat='I', mask=TRUE, checkNegatives=TRUE)
}
tjuv=as.matrix(unlist(juvno))
plot(tjuv~yrlist, type='b', ylim=c(.75,1), ylab='juvenile niche overlap', xlab='')

ichno=data.frame(matrix(NA, ncol=1, nrow=dim(ich1)[3])); colnames(ichno)='I'
for( i in 1:dim(ich1)[3]){
  ichno[i,1]=nicheOverlap(ich1[[i]], ich2[[i]], stat='I', mask=TRUE, checkNegatives=TRUE)
}
tich=as.matrix(unlist(ichno))
plot(tich~yrlist, type='b',ylim=c(.85,1), ylab='ichthyo niche overlap', xlab='')

plot(NULL, ylim=c(0.85,1), xlim=c(1977,2019), ylab="niche overlap", xlab="")
points(tadt~yrlist, type='b', col='red', pch=19)
points(tjuv~yrlist, type='b', col='blue', pch=19)
points(tich~yrlist, type='b', col='green', pch=19)

### save data out
Spr.Niche.Overlap=data.frame(matrix(NA, ncol=4, nrow=dim(ich1)[3])); colnames(Spr.Niche.Overlap)=c('year', 'adt', 'juv', 'ich')
Spr.Niche.Overlap[,1]=yrlist
Spr.Niche.Overlap[,2]=adtno
Spr.Niche.Overlap[,3]=juvno
Spr.Niche.Overlap[,4]=ichno
save(Spr.Niche.Overlap, file='/home/ryan/Git/NEhabitat/habitat index/Spring_Niche_Overlap_Cod_Haddock_GSe_20210518.rda')

### SSIM ### 
#1) Seasonal Inter-species habitat overlap: Cod and Haddock by season, stage ###
adtssim=data.frame(matrix(NA, ncol=1, nrow=dim(adt1)[3])); colnames(adtssim)='SSIM'
for( i in 1:dim(adt1)[3]){
  r1m = rasterToPoints(adt1[[i]])
  r2m = rasterToPoints(adt2[[i]])
  adtssim[i,1]=SSIM(r1m[,3], r2m[,3], numBreaks = 256)
}
plot(adtssim[,1]~yrlist, type='b', main='Adt SSIM')

juvssim=data.frame(matrix(NA, ncol=1, nrow=dim(adt1)[3])); colnames(juvssim)='SSIM'
for( i in 1:dim(adt1)[3]){
  r1m = rasterToPoints(juv1[[i]])
  r2m = rasterToPoints(juv2[[i]])
  juvssim[i,1]=SSIM(r1m[,3], r2m[,3], numBreaks = 256)
}
plot(juvssim[,1]~yrlist, type='b', main='Juv SSIM')

ichssim=data.frame(matrix(NA, ncol=1, nrow=dim(adt1)[3])); colnames(ichssim)='SSIM'
for( i in 1:dim(adt1)[3]){
  r1m = rasterToPoints(ich1[[i]])
  r2m = rasterToPoints(ich2[[i]])
  ichssim[i,1]=SSIM(r1m[,3], r2m[,3], numBreaks = 256)
}
plot(ichssim[,1]~yrlist, type='b', main='Ich SSIM')

plot(NULL, ylim=c(0.5,1), xlim=c(1977,2019), ylab="SSIM", xlab="")
# plot(NULL, ylim=c(0.2,1), xlim=c(1977,2019), ylab="SSIM", xlab="")
points(adtssim[,1]~yrlist, type='b', col='red', pch=19)
points(juvssim[,1]~yrlist, type='b', col='blue', pch=19)
points(ichssim[,1]~yrlist, type='b', col='green', pch=19)

### save data out
Spr.SSIM=data.frame(matrix(NA, ncol=4, nrow=dim(ich1)[3])); colnames(Spr.SSIM)=c('year', 'adt', 'juv', 'ich')
Spr.SSIM[,1]=yrlist
Spr.SSIM[,2]=adtssim
Spr.SSIM[,3]=juvssim
Spr.SSIM[,4]=ichssim
save(Spr.SSIM, file='/home/ryan/Git/NEhabitat/habitat index/Spring_SSIM_Cod_Haddock_GSe_20210518.rda')


## 20210317 add new models
#fall GSe haddock
# juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock/PA_only_stacked_Fall_Haddock_Juv.RData')
# adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock/PA_only_stacked_Fall_Haddock_Adt.RData')
## 20210421 add new models for haddock
#fall GSe haddock
juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock_20210421/PA_only_stacked_Fall_Haddock_Juv.RData')
adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock_20210421/PA_only_stacked_Fall_Haddock_Adt.RData')
#fall GSe cod
# juv2=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod/PA_only_stacked_Fall_Cod_Juv.RData')
# adt2=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod/PA_only_stacked_Fall_Cod_Adt.RData')
#fall GSe cod updated from above (new data from Sean survdat)
juv2=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod_20210518/PA_only_stacked_Fall_Cod_Juv.RData')
adt2=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod_20210518/PA_only_stacked_Fall_Cod_Adt.RData')
# #fall GI haddock
# juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGI_Fall_haddock/PA_only_stacked_Fall_Haddock_Juv.RData')
# adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGI_Fall_haddock/PA_only_stacked_Fall_Haddock_Adt.RData')
# #fall GI cod
# juv2=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGI_Fall_cod/PA_only_stacked_Fall_Cod_Juv.RData')
# adt2=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGI_Fall_cod/PA_only_stacked_Fall_Cod_Adt.RData')

adtno=data.frame(matrix(NA, ncol=1, nrow=dim(adt1)[3])); colnames(adtno)='I'
for( i in 1:dim(adt1)[3]){
  adtno[i,1]=nicheOverlap(adt1[[i]], adt2[[i]], stat='I', mask=TRUE, checkNegatives=TRUE)
}
t1=as.matrix(unlist(adtno))
plot(t1~yrlist, type='b',ylim=c(.85,1), main='Adt I-stat')

juvno=data.frame(matrix(NA, ncol=1, nrow=dim(juv1)[3])); colnames(juvno)='I'
for( i in 1:dim(juv1)[3]){
  juvno[i,1]=nicheOverlap(juv1[[i]], juv2[[i]], stat='I', mask=TRUE, checkNegatives=TRUE)
}
t2=as.matrix(unlist(juvno))
plot(t2~yrlist, type='b', ylim=c(.7,1),main='Juv I-stat')

plot(NULL, ylim=c(0.85,1), xlim=c(1977,2019), ylab="niche overlap", xlab="", las=1)
points(t1~yrlist, type='b', col='red', pch=19)
points(t2~yrlist, type='b', col='blue', pch=19)
legend('topleft', legend = c('adt', 'juv'), pch=19, col=c('red', 'blue'), bty='n', horiz = T)

### save data out
Fall.Niche.Overlap=data.frame(matrix(NA, ncol=3, nrow=dim(adt1)[3])); colnames(Fall.Niche.Overlap)=c('year', 'adt', 'juv')
Fall.Niche.Overlap[,1]=yrlist
Fall.Niche.Overlap[,2]=adtno
Fall.Niche.Overlap[,3]=juvno
# Fall.Niche.Overlap[,4]=ichno
save(Fall.Niche.Overlap, file='/home/ryan/Git/NEhabitat/habitat index/Fall_Niche_Overlap_Cod_Haddock_GSe_20210518.rda')

### SSIM ### 
adtssim=data.frame(matrix(NA, ncol=1, nrow=dim(adt1)[3])); colnames(adtssim)='SSIM'
for( i in 1:dim(adt1)[3]){
  r1m = rasterToPoints(adt1[[i]])
  r2m = rasterToPoints(adt2[[i]])
  adtssim[i,1]=SSIM(r1m[,3], r2m[,3], numBreaks = 256)
}
plot(adtssim[,1]~yrlist, type='b', main='Adt SSIM')

juvssim=data.frame(matrix(NA, ncol=1, nrow=dim(adt1)[3])); colnames(juvssim)='SSIM'
for( i in 1:dim(adt1)[3]){
  r1m = rasterToPoints(juv1[[i]])
  r2m = rasterToPoints(juv2[[i]])
  juvssim[i,1]=SSIM(r1m[,3], r2m[,3], numBreaks = 256)
}
plot(juvssim[,1]~yrlist, type='b', main='Juv SSIM')

# plot(NULL, ylim=c(0.5,1), xlim=c(1977,2019), ylab="SSIM", xlab="")
plot(NULL, ylim=c(0.2,1), xlim=c(1977,2019), ylab="SSIM", xlab="")
points(adtssim[,1]~yrlist, type='b', col='red', pch=19)
points(juvssim[,1]~yrlist, type='b', col='blue', pch=19)

### save data out
Fall.SSIM=data.frame(matrix(NA, ncol=3, nrow=dim(ich1)[3])); colnames(Fall.SSIM)=c('year', 'adt', 'juv')
Fall.SSIM[,1]=yrlist
Fall.SSIM[,2]=adtssim
Fall.SSIM[,3]=juvssim
save(Fall.SSIM, file='/home/ryan/Git/NEhabitat/habitat index/Fall_SSIM_Cod_Haddock_GSe_20210518.rda')


### 2) Intra-species habitat overlap: between seasons by stage and species ###

## 20210317 add new models
# # fall GSe haddock
# juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock/PA_only_stacked_Fall_Haddock_Juv.RData')
# adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock/PA_only_stacked_Fall_Haddock_Adt.RData')
# #spring GSe haddock
# juv2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_haddock/PA_only_stacked_Spr_Haddock_Juv.RData')
# adt2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_haddock/PA_only_stacked_Spr_Haddock_Adt.RData')
#fall GSe haddock
juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock_20210421/PA_only_stacked_Fall_Haddock_Juv.RData')
adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock_20210421/PA_only_stacked_Fall_Haddock_Adt.RData')
#spring GSe haddock
juv2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/PA_only_stacked_Spr_Haddock_Juv.RData')
adt2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/PA_only_stacked_Spr_Haddock_Adt.RData')

adtno=data.frame(matrix(NA, ncol=1, nrow=dim(adt1)[3])); colnames(adtno)='I'
for( i in 1:dim(adt1)[3]){
  adtno[i,1]=nicheOverlap(adt1[[i]], adt2[[i]], stat='I', mask=TRUE, checkNegatives=TRUE)
}
t1=as.matrix(unlist(adtno))
plot(t1~yrlist, type='b',ylim=c(.85,1), main='Adt I-stat')

juvno=data.frame(matrix(NA, ncol=1, nrow=dim(juv1)[3])); colnames(juvno)='I'
for( i in 1:dim(juv1)[3]){
  juvno[i,1]=nicheOverlap(juv1[[i]], juv2[[i]], stat='I', mask=TRUE, checkNegatives=TRUE)
}
t2=as.matrix(unlist(juvno))
plot(t2~yrlist, type='b', ylim=c(.85,1),main='Juv I-stat')
plot(NULL, ylim=c(0.85,1), xlim=c(1977,2019), ylab="niche overlap", xlab="", las=1)
points(t1~yrlist, type='b', col='red', pch=19)
points(t2~yrlist, type='b', col='blue', pch=19)
legend('topleft', legend = c('adt', 'juv'), pch=19, col=c('red', 'blue'), bty='n', horiz = T)

### save data out
Haddock.Niche.Overlap=data.frame(matrix(NA, ncol=3, nrow=dim(adt1)[3])); colnames(Haddock.Niche.Overlap)=c('year', 'adt', 'juv')
Haddock.Niche.Overlap[,1]=yrlist
Haddock.Niche.Overlap[,2]=adtno
Haddock.Niche.Overlap[,3]=juvno
save(Haddock.Niche.Overlap, file='/home/ryan/Git/NEhabitat/habitat index/Haddock_Spr_Fall_Niche_Overlap_GSe_20210518.rda')

## SSIM
adtssim=data.frame(matrix(NA, ncol=1, nrow=dim(adt1)[3])); colnames(adtssim)='SSIM'
for( i in 1:dim(adt1)[3]){
  r1m = rasterToPoints(adt1[[i]])
  r2m = rasterToPoints(adt2[[i]])
  adtssim[i,1]=SSIM(r1m[,3], r2m[,3], numBreaks = 256)
}
plot(adtssim[,1]~yrlist, type='b', main='Adt SSIM')

juvssim=data.frame(matrix(NA, ncol=1, nrow=dim(adt1)[3])); colnames(juvssim)='SSIM'
for( i in 1:dim(adt1)[3]){
  r1m = rasterToPoints(juv1[[i]])
  r2m = rasterToPoints(juv2[[i]])
  juvssim[i,1]=SSIM(r1m[,3], r2m[,3], numBreaks = 256)
}
plot(juvssim[,1]~yrlist, type='b', main='Juv SSIM')

plot(NULL, ylim=c(0.5,1), xlim=c(1977,2019), ylab="SSIM", xlab="")
# plot(NULL, ylim=c(0.2,1), xlim=c(1977,2019), ylab="SSIM", xlab="")
points(adtssim[,1]~yrlist, type='b', col='red', pch=19)
points(juvssim[,1]~yrlist, type='b', col='blue', pch=19)
### save data out
Haddock.SSIM=data.frame(matrix(NA, ncol=3, nrow=dim(adt1)[3])); colnames(Haddock.SSIM)=c('year', 'adt', 'juv')
Haddock.SSIM[,1]=yrlist
Haddock.SSIM[,2]=adtssim
Haddock.SSIM[,3]=juvssim
save(Haddock.SSIM, file='/home/ryan/Git/NEhabitat/habitat index/Haddock_Spr_Fall_SSIM_GSe_20210518.rda')

#fall GSe cod updated from above (new data from Sean survdat)
juv2=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod_20210518/PA_only_stacked_Fall_Cod_Juv.RData')
adt2=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod_20210518/PA_only_stacked_Fall_Cod_Adt.RData')
#spring GSe cod
juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod_20210518/PA_only_stacked_Spr_Cod_Juv.RData')
adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod_20210518/PA_only_stacked_Spr_Cod_Adt.RData')

adtno=data.frame(matrix(NA, ncol=1, nrow=dim(adt1)[3])); colnames(adtno)='I'
for( i in 1:dim(adt1)[3]){
  adtno[i,1]=nicheOverlap(adt1[[i]], adt2[[i]], stat='I', mask=TRUE, checkNegatives=TRUE)
}
t1=as.matrix(unlist(adtno))
plot(t1~yrlist, type='b',ylim=c(.85,1), main='Adt I-stat')

juvno=data.frame(matrix(NA, ncol=1, nrow=dim(juv1)[3])); colnames(juvno)='I'
for( i in 1:dim(juv1)[3]){
  juvno[i,1]=nicheOverlap(juv1[[i]], juv2[[i]], stat='I', mask=TRUE, checkNegatives=TRUE)
}
t2=as.matrix(unlist(juvno))
plot(t2~yrlist, type='b', ylim=c(.85,1),main='Juv I-stat')
plot(NULL, ylim=c(0.85,1), xlim=c(1977,2019), ylab="niche overlap", xlab="", las=1)
points(t1~yrlist, type='b', col='red', pch=19)
points(t2~yrlist, type='b', col='blue', pch=19)
legend('topleft', legend = c('adt', 'juv'), pch=19, col=c('red', 'blue'), bty='n', horiz = T)

### save data out
Cod.Niche.Overlap=data.frame(matrix(NA, ncol=3, nrow=dim(adt1)[3])); colnames(Cod.Niche.Overlap)=c('year', 'adt', 'juv')
Cod.Niche.Overlap[,1]=yrlist
Cod.Niche.Overlap[,2]=adtno
Cod.Niche.Overlap[,3]=juvno
save(Cod.Niche.Overlap, file='/home/ryan/Git/NEhabitat/habitat index/Cod_Spr_Fall_Niche_Overlap_GSe_20210518.rda')

## SSIM
adtssim=data.frame(matrix(NA, ncol=1, nrow=dim(adt1)[3])); colnames(adtssim)='SSIM'
for( i in 1:dim(adt1)[3]){
  r1m = rasterToPoints(adt1[[i]])
  r2m = rasterToPoints(adt2[[i]])
  adtssim[i,1]=SSIM(r1m[,3], r2m[,3], numBreaks = 256)
}
plot(adtssim[,1]~yrlist, type='b', main='Adt SSIM')

juvssim=data.frame(matrix(NA, ncol=1, nrow=dim(adt1)[3])); colnames(juvssim)='SSIM'
for( i in 1:dim(adt1)[3]){
  r1m = rasterToPoints(juv1[[i]])
  r2m = rasterToPoints(juv2[[i]])
  juvssim[i,1]=SSIM(r1m[,3], r2m[,3], numBreaks = 256)
}
plot(juvssim[,1]~yrlist, type='b', main='Juv SSIM')

plot(NULL, ylim=c(0.5,1), xlim=c(1977,2019), ylab="SSIM", xlab="")
# plot(NULL, ylim=c(0.2,1), xlim=c(1977,2019), ylab="SSIM", xlab="")
points(adtssim[,1]~yrlist, type='b', col='red', pch=19)
points(juvssim[,1]~yrlist, type='b', col='blue', pch=19)
### save data out
Cod.SSIM=data.frame(matrix(NA, ncol=3, nrow=dim(adt1)[3])); colnames(Cod.SSIM)=c('year', 'adt', 'juv')
Cod.SSIM[,1]=yrlist
Cod.SSIM[,2]=adtssim
Cod.SSIM[,3]=juvssim
save(Cod.SSIM, file='/home/ryan/Git/NEhabitat/habitat index/Cod_Spr_Fall_SSIM_GSe_20210518.rda')


### 3) Life stage habitat overlaps: Cod and Haddock by season: adt-juv, adt-ich, juv-ich

## Choose season and fish, load files and run (as appropriate with ich)fSEAS='Spr'
fSEAS='Fall' #Spr Fall
fFISH='Haddock' # Cod Haddock
#spring GSe haddock
ich1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/PA_only_stacked_Spr_Haddock_ich.RData')
juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/PA_only_stacked_Spr_Haddock_Juv.RData')
adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/PA_only_stacked_Spr_Haddock_Adt.RData')
#spring GSe cod
ich1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod_20210518/PA_only_stacked_Spr_Cod_ich.RData')
juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod_20210518/PA_only_stacked_Spr_Cod_Juv.RData')
adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod_20210518/PA_only_stacked_Spr_Cod_Adt.RData')
#fall GSe haddock
juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock_20210421/PA_only_stacked_Fall_Haddock_Juv.RData')
adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock_20210421/PA_only_stacked_Fall_Haddock_Adt.RData')
#fall GSe cod updated from above (new data from Sean survdat)
juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod_20210518/PA_only_stacked_Fall_Cod_Juv.RData')
adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod_20210518/PA_only_stacked_Fall_Cod_Adt.RData')
## Choose above season and fish and run (as appropriate with ich)
adtjuvno=data.frame(matrix(NA, ncol=1, nrow=dim(adt1)[3])); colnames(adtjuvno)='I'
for( i in 1:dim(adt1)[3]){
  adtjuvno[i,1]=nicheOverlap(adt1[[i]], juv1[[i]], stat='I', mask=TRUE, checkNegatives=TRUE)
}
t1=as.matrix(unlist(adtjuvno))
# plot(t1~yrlist, type='b',ylim=c(.85,1), main='AdtJuv I-stat')
if (fSEAS=='Fall'){
  break
} else {
adtichno=data.frame(matrix(NA, ncol=1, nrow=dim(adt1)[3])); colnames(adtichno)='I'
for( i in 1:dim(adt1)[3]){
  adtichno[i,1]=nicheOverlap(adt1[[i]], ich1[[i]], stat='I', mask=TRUE, checkNegatives=TRUE)
}
t2=as.matrix(unlist(adtichno))
# plot(t2~yrlist, type='b',ylim=c(.85,1), main='AdtIch I-stat')
juvichno=data.frame(matrix(NA, ncol=1, nrow=dim(adt1)[3])); colnames(juvichno)='I'
for( i in 1:dim(adt1)[3]){
  juvichno[i,1]=nicheOverlap(juv1[[i]], ich1[[i]], stat='I', mask=TRUE, checkNegatives=TRUE)
}
t3=as.matrix(unlist(juvichno))
# plot(t3~yrlist, type='b',ylim=c(.85,1), main='JuvIch I-stat')
}
if (fSEAS=='Spr'){
plot(NULL, ylim=c(0.85,1), xlim=c(1977,2019), ylab="niche overlap", xlab="", las=1, main=paste(fSEAS, fFISH, sep=' '))
points(t1~yrlist, type='b', col='red', pch=19)
points(t2~yrlist, type='b', col='blue', pch=19)
points(t3~yrlist, type='b', col='green', pch=19)
legend('topleft', legend = c('adtjuv', 'adtich', 'juvich'), pch=19, col=c('red', 'blue', 'green'), bty='n', horiz = T)
} else if (fSEAS=='Fall'){
## FALL plots
plot(NULL, ylim=c(0.85,1), xlim=c(1977,2019), ylab="niche overlap", xlab="", las=1 ,main=paste(fSEAS, fFISH, sep=' '))
points(t1~yrlist, type='b', col='red', pch=19)
legend('topleft', legend = c('adtjuv'), pch=19, col=c('red'), bty='n', horiz = T)
}

### save data out 
# Spr haddock
Haddock.Stage.Niche.Overlap=data.frame(matrix(NA, ncol=4, nrow=dim(adt1)[3])); colnames(Haddock.Stage.Niche.Overlap)=c('year', 'adtjuv', 'adtich', 'juvich')
Haddock.Stage.Niche.Overlap[,1]=yrlist
Haddock.Stage.Niche.Overlap[,2]=adtjuvno
Haddock.Stage.Niche.Overlap[,3]=adtichno
Haddock.Stage.Niche.Overlap[,4]=juvichno
save(Haddock.Stage.Niche.Overlap, file='/home/ryan/Git/NEhabitat/habitat index/Haddock_Spr_stage_Niche_Overlap_GSe_20210518.rda')
# Fall Haddock
Haddock.Stage.Niche.Overlap=data.frame(matrix(NA, ncol=2, nrow=dim(adt1)[3])); colnames(Haddock.Stage.Niche.Overlap)=c('year', 'adtjuv')
Haddock.Stage.Niche.Overlap[,1]=yrlist
Haddock.Stage.Niche.Overlap[,2]=adtjuvno
save(Haddock.Stage.Niche.Overlap, file='/home/ryan/Git/NEhabitat/habitat index/Haddock_Fall_stage_Niche_Overlap_GSe_20210518.rda')
# Spring Cod
Cod.Stage.Niche.Overlap=data.frame(matrix(NA, ncol=4, nrow=dim(adt1)[3])); colnames(Cod.Stage.Niche.Overlap)=c('year', 'adtjuv', 'adtich', 'juvich')
Cod.Stage.Niche.Overlap[,1]=yrlist
Cod.Stage.Niche.Overlap[,2]=adtjuvno
Cod.Stage.Niche.Overlap[,3]=adtichno
Cod.Stage.Niche.Overlap[,4]=juvichno
save(Cod.Stage.Niche.Overlap, file='/home/ryan/Git/NEhabitat/habitat index/Cod_Spr_stage_Niche_Overlap_GSe_20210518.rda')
# Fall Cod
Cod.Stage.Niche.Overlap=data.frame(matrix(NA, ncol=2, nrow=dim(adt1)[3])); colnames(Cod.Stage.Niche.Overlap)=c('year', 'adtjuv')
Cod.Stage.Niche.Overlap[,1]=yrlist
Cod.Stage.Niche.Overlap[,2]=adtjuvno
save(Cod.Stage.Niche.Overlap, file='/home/ryan/Git/NEhabitat/habitat index/Cod_Fall_stage_Niche_Overlap_GSe_20210518.rda')

## SSIM
adtjuvssim=data.frame(matrix(NA, ncol=1, nrow=dim(adt1)[3])); colnames(adtjuvssim)='SSIM'
for( i in 1:dim(adt1)[3]){
  r1m = rasterToPoints(adt1[[i]])
  r2m = rasterToPoints(juv1[[i]])
  adtjuvssim[i,1]=SSIM(r1m[,3], r2m[,3], numBreaks = 256)
}
# plot(tadtssim~yrlist, type='b', main='Adt SSIM')
if (fSEAS=='Fall'){
  break
} else {
  adtichssim=data.frame(matrix(NA, ncol=1, nrow=dim(adt1)[3])); colnames(adtichssim)='SSIM'
  for( i in 1:dim(adt1)[3]){
    r1m = rasterToPoints(adt1[[i]])
    r2m = rasterToPoints(ich1[[i]])
    adtichssim[i,1]=SSIM(r1m[,3], r2m[,3], numBreaks = 256)
  }
  juvichssim=data.frame(matrix(NA, ncol=1, nrow=dim(adt1)[3])); colnames(juvichssim)='SSIM'
  for( i in 1:dim(adt1)[3]){
    r1m = rasterToPoints(juv1[[i]])
    r2m = rasterToPoints(ich1[[i]])
    juvichssim[i,1]=SSIM(r1m[,3], r2m[,3], numBreaks = 256)
  }
}
if (fSEAS=='Spr'){
  plot(NULL, ylim=c(0.5,1), xlim=c(1977,2019), ylab="SSIM", xlab="", las=1, main=paste(fSEAS, fFISH, sep=' '))
  points(adtjuvssim[,1]~yrlist, type='b', col='red', pch=19)
  points(adtichssim[,1]~yrlist, type='b', col='blue', pch=19)
  points(juvichssim[,1]~yrlist, type='b', col='green', pch=19)
  legend('topleft', legend = c('adtjuv', 'adtich', 'juvich'), pch=19, col=c('red', 'blue', 'green'), bty='n', horiz = T)
} else if (fSEAS=='Fall'){
  ## FALL plots
  plot(NULL, ylim=c(0.5,1), xlim=c(1977,2019), ylab="niche overlap", xlab="", las=1 ,main=paste(fSEAS, fFISH, sep=' '))
  points(adtjuvssim[,1]~yrlist, type='b', col='red', pch=19)
  legend('topleft', legend = c('adtjuv'), pch=19, col=c('red'), bty='n', horiz = T)
}
### save data out 
# Spr haddock
Haddock.Stage.SSIM=data.frame(matrix(NA, ncol=4, nrow=dim(adt1)[3])); colnames(Haddock.Stage.SSIM)=c('year', 'adtjuv', 'adtich', 'juvich')
Haddock.Stage.SSIM[,1]=yrlist
Haddock.Stage.SSIM[,2]=adtjuvssim
Haddock.Stage.SSIM[,3]=adtichssim
Haddock.Stage.SSIM[,4]=juvichssim
save(Haddock.Stage.SSIM, file='/home/ryan/Git/NEhabitat/habitat index/Haddock_Spr_stage_SSIM_GSe_20210518.rda')
# Fall Haddock
Haddock.Stage.SSIM=data.frame(matrix(NA, ncol=2, nrow=dim(adt1)[3])); colnames(Haddock.Stage.SSIM)=c('year', 'adtjuv')
Haddock.Stage.SSIM[,1]=yrlist
Haddock.Stage.SSIM[,2]=adtjuvssim
save(Haddock.Stage.SSIM, file='/home/ryan/Git/NEhabitat/habitat index/Haddock_Fall_stage_SSIM_GSe_20210518.rda')
# Spring Cod
Cod.Stage.SSIM=data.frame(matrix(NA, ncol=4, nrow=dim(adt1)[3])); colnames(Cod.Stage.SSIM)=c('year', 'adtjuv', 'adtich', 'juvich')
Cod.Stage.SSIM[,1]=yrlist
Cod.Stage.SSIM[,2]=adtjuvssim
Cod.Stage.SSIM[,3]=adtichssim
Cod.Stage.SSIM[,4]=juvichssim
save(Cod.Stage.SSIM, file='/home/ryan/Git/NEhabitat/habitat index/Cod_Spr_stage_SSIM_GSe_20210518.rda')
# Fall Cod
Cod.Stage.SSIM=data.frame(matrix(NA, ncol=2, nrow=dim(adt1)[3])); colnames(Cod.Stage.SSIM)=c('year', 'adtjuv')
Cod.Stage.SSIM[,1]=yrlist
Cod.Stage.SSIM[,2]=adtjuvssim
save(Cod.Stage.SSIM, file='/home/ryan/Git/NEhabitat/habitat index/Cod_Fall_stage_SSIM_GSe_20210518.rda')







# PLOTTING SAVED DATA

### Niche Overlap between stages (species, season)
# Spring Cod
load('/home/ryan/Git/NEhabitat/habitat index/Cod_Spr_stage_Niche_Overlap_GSe_20210518.rda')
plot(NULL, ylim=c(0.85,1), xlim=c(1977,2019), ylab="niche overlap", xlab="", las=1 ,main='Spr Cod')
points(Cod.Stage.Niche.Overlap$adtjuv~yrlist, type='b', col='red', pch=19)
points(Cod.Stage.Niche.Overlap$adtich~yrlist, type='b', col='blue', pch=19)
points(Cod.Stage.Niche.Overlap$juvich~yrlist, type='b', col='green', pch=19)
legend('topleft', legend = c('adt-juv','adt-ich', 'juv-ich'), pch=19, col=c('red', 'blue', 'green'), bty='n', horiz = T)
# Spr Haddock
load('/home/ryan/Git/NEhabitat/habitat index/Haddock_Spr_stage_Niche_Overlap_GSe_20210518.rda')
plot(NULL, ylim=c(0.85,1), xlim=c(1977,2019), ylab="niche overlap", xlab="", las=1 ,main='Spr Haddock')
points(Haddock.Stage.Niche.Overlap$adtjuv~yrlist, type='b', col='red', pch=19)
points(Haddock.Stage.Niche.Overlap$adtich~yrlist, type='b', col='blue', pch=19)
points(Haddock.Stage.Niche.Overlap$juvich~yrlist, type='b', col='green', pch=19)
legend('topleft', legend = c('adt-juv','adt-ich', 'juv-ich'), pch=19, col=c('red', 'blue', 'green'), bty='n', horiz = T)
# Fall (both)
load('/home/ryan/Git/NEhabitat/habitat index/Haddock_Fall_stage_Niche_Overlap_GSe_20210518.rda')
load('/home/ryan/Git/NEhabitat/habitat index/Cod_Fall_stage_Niche_Overlap_GSe_20210518.rda')
plot(NULL, ylim=c(0.85,1), xlim=c(1977,2019), ylab="niche overlap", xlab="", las=1 ,main='Fall Adt-Juv')
points(Cod.Stage.Niche.Overlap$adtjuv~yrlist, type='b', col='red', pch=19)
points(Haddock.Stage.Niche.Overlap$adtjuv~yrlist, type='b', col='blue', pch=19)
legend('topleft', legend = c('Cod','Haddock'), pch=19, col=c('red', 'blue'), bty='n', horiz = T)

### Niche Overlap between season (species, stage)
# Haddock
load('/home/ryan/Git/NEhabitat/habitat index/Haddock_Spr_Fall_Niche_Overlap_GSe_20210518.rda')
plot(NULL, ylim=c(0.85,1), xlim=c(1977,2019), ylab="niche overlap", xlab="", las=1 ,main='Haddock')
points(Haddock.Niche.Overlap$adt~yrlist, type='b', col='red', pch=19)
points(Haddock.Niche.Overlap$juv~yrlist, type='b', col='blue', pch=19)
legend('topleft', legend = c('Adt','Juv'), pch=19, col=c('red', 'blue'), bty='n', horiz = T)
# Cod
load('/home/ryan/Git/NEhabitat/habitat index/Cod_Spr_Fall_Niche_Overlap_GSe_20210518.rda')
plot(NULL, ylim=c(0.85,1), xlim=c(1977,2019), ylab="niche overlap", xlab="", las=1 ,main='Cod')
points(Cod.Niche.Overlap$adt~yrlist, type='b', col='red', pch=19)
points(Cod.Niche.Overlap$juv~yrlist, type='b', col='blue', pch=19)
legend('topleft', legend = c('Adt','Juv'), pch=19, col=c('red', 'blue'), bty='n', horiz = T)

### Niche Overlap between species (season, stage)
# Spr
load('/home/ryan/Git/NEhabitat/habitat index/Spring_Niche_Overlap_Cod_Haddock_GSe_20210518.rda')
plot(NULL, ylim=c(0.8,1), xlim=c(1977,2019), ylab="niche overlap", xlab="", las=1 ,main='Spr Cod-Haddock')
points(Spr.Niche.Overlap$adt~yrlist, type='b', col='red', pch=19)
points(Spr.Niche.Overlap$juv~yrlist, type='b', col='blue', pch=19)
points(Spr.Niche.Overlap$ich~yrlist, type='b', col='green', pch=19)
legend('topleft', legend = c('Adt', 'Juv','Ich'), pch=19, col=c('red', 'blue', 'green'), bty='n', horiz = T)
# Fall
load('/home/ryan/Git/NEhabitat/habitat index/Fall_Niche_Overlap_Cod_Haddock_GSe_20210518.rda')
plot(NULL, ylim=c(0.8,1), xlim=c(1977,2019), ylab="niche overlap", xlab="", las=1 ,main='Fall Cod-Haddock')
points(Fall.Niche.Overlap$adt~yrlist, type='b', col='red', pch=19)
points(Fall.Niche.Overlap$juv~yrlist, type='b', col='blue', pch=19)
legend('topleft', legend = c('Adt', 'Juv'), pch=19, col=c('red', 'blue'), bty='n', horiz = T)

#### Structural Similarity Index (SSIM) ####_________________________________________________________________________________
#### INTRASPECIES DIFFERENCES

#spring GSe haddock
ich1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/PA_only_stacked_Spr_Haddock_ich.RData')
juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/PA_only_stacked_Spr_Haddock_Juv.RData')
adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/PA_only_stacked_Spr_Haddock_Adt.RData')
#spring GSe cod updated from above (new data from Sean survdat)
ich2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod_20210518/PA_only_stacked_Spr_Cod_ich.RData')
juv2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod_20210518/PA_only_stacked_Spr_Cod_Juv.RData')
adt2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod_20210518/PA_only_stacked_Spr_Cod_Adt.RData')
#spring GSe haddock
# ich1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_haddock/PA_only_stacked_Spr_Haddock_ich.RData')
# juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_haddock/PA_only_stacked_Spr_Haddock_Juv.RData')
# adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_haddock/PA_only_stacked_Spr_Haddock_Adt.RData')
# #spring GSe cod
# ich2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod/PA_only_stacked_Spr_Cod_ich.RData')
# juv2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod/PA_only_stacked_Spr_Cod_Juv.RData')
# adt2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod/PA_only_stacked_Spr_Cod_Adt.RData')
#fall GSe cod
# juv2=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod/PA_only_stacked_Fall_Cod_Juv.RData')
# adt2=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod/PA_only_stacked_Fall_Cod_Adt.RData')
# #fall GSe haddock
# juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock/PA_only_stacked_Fall_Haddock_Juv.RData')
# adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock/PA_only_stacked_Fall_Haddock_Adt.RData')
#fall GSe haddock
juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock_20210421/PA_only_stacked_Fall_Haddock_Juv.RData')
adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock_20210421/PA_only_stacked_Fall_Haddock_Adt.RData')
#fall GSe cod
# juv2=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod/PA_only_stacked_Fall_Cod_Juv.RData')
# adt2=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod/PA_only_stacked_Fall_Cod_Adt.RData')
#fall GSe cod updated from above (new data from Sean survdat)
juv2=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod_20210518/PA_only_stacked_Fall_Cod_Juv.RData')
adt2=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod_20210518/PA_only_stacked_Fall_Cod_Adt.RData')


### Structural similarity index

##### Seasonal difference  within species spring vs fall GSe 
# ## HADDOCK
# #spring GSe haddock
# juv2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_haddock/PA_only_stacked_Spr_Haddock_Juv.RData')
# adt2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_haddock/PA_only_stacked_Spr_Haddock_Adt.RData')
# #fall GSe haddock
# juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock/PA_only_stacked_Fall_Haddock_Juv.RData')
# adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock/PA_only_stacked_Fall_Haddock_Adt.RData')
# 
# ## COD
# #spring GSe cod
# juv2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod/PA_only_stacked_Spr_Cod_Juv.RData')
# adt2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod/PA_only_stacked_Spr_Cod_Adt.RData')
# #fall GSe cod
# juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod/PA_only_stacked_Fall_Cod_Juv.RData')
# adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod/PA_only_stacked_Fall_Cod_Adt.RData')


adtssim=data.frame(matrix(NA, ncol=1, nrow=dim(adt1)[3])); colnames(adtssim)='SSIM'
for( i in 1:dim(adt1)[3]){
  r1m = rasterToPoints(adt1[[i]])
  r2m = rasterToPoints(adt2[[i]])
  adtssim[i,1]=SSIM(r1m[,3], r2m[,3], numBreaks = 256)
}
tadtssim=as.matrix(unlist(adtssim))
plot(tadtssim~yrlist, type='b', main='Adt SSIM')

juvssim=data.frame(matrix(NA, ncol=1, nrow=dim(adt1)[3])); colnames(juvssim)='SSIM'
for( i in 1:dim(adt1)[3]){
  r1m = rasterToPoints(juv1[[i]])
  r2m = rasterToPoints(juv2[[i]])
  juvssim[i,1]=SSIM(r1m[,3], r2m[,3], numBreaks = 256)
}
tjuvssim=as.matrix(unlist(juvssim))
plot(tjuvssim~yrlist, type='b', main='Juv SSIM')

ichssim=data.frame(matrix(NA, ncol=1, nrow=dim(adt1)[3])); colnames(ichssim)='SSIM'
for( i in 1:dim(adt1)[3]){
  r1m = rasterToPoints(ich1[[i]])
  r2m = rasterToPoints(ich2[[i]])
  ichssim[i,1]=SSIM(r1m[,3], r2m[,3], numBreaks = 256)
}
tichssim=as.matrix(unlist(ichssim))
plot(tichssim~yrlist, type='b', main='Ich SSIM')

plot(NULL, ylim=c(0.5,1), xlim=c(1977,2019), ylab="SSIM", xlab="")
plot(NULL, ylim=c(0.2,1), xlim=c(1977,2019), ylab="SSIM", xlab="")

points(tadtssim~yrlist, type='b', col='red', pch=19)
points(tjuvssim~yrlist, type='b', col='blue', pch=19)
points(tichssim~yrlist, type='b', col='green', pch=19)

# adtno=data.frame(matrix(NA, ncol=1, nrow=dim(adt1)[3])); colnames(adtno)='SSIM'
# for( i in 1:dim(adt1)[3]){
#   r1m = rasterToPoints(kfpaf[[i]])
#   r2m = rasterToPoints(kfpas[[i]])
#   adtno[i,1]=SSIM(r1m[,3], r2m[,3], numBreaks = 256)
#   # ssim = SSIM(r1m[,3], r2m[,3], numBreaks = 256)
# }
# t=as.matrix(unlist(adtno))
# plot(t~yrlist, type='b',ylim=c(.85,1), main='Adt SSIM')







