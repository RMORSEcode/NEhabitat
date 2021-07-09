
# Using numbers per tow (biomass) models instead of PA

#### Niche Overlap ####_________________________________________________________________________________
### 1) Seasonal Inter-species habitat overlap: Cod and Haddock by season, stage ###

# 20210421 add new models for haddock
#spring GSe haddock
ich1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/stacked_Spr_Haddock_ich.RData')
juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/stacked_Spr_Haddock_Juv.RData')
adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/stacked_Spr_Haddock_Adt.RData')
#spring GSe cod
# ich2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod/stacked_Spr_Cod_ich.RData')
# juv2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod/stacked_Spr_Cod_Juv.RData')
# adt2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod/stacked_Spr_Cod_Adt.RData')
#spring GSe cod updated from above (new data from Sean survdat)
ich2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod_20210518/stacked_Spr_Cod_ich.RData')
juv2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod_20210518/stacked_Spr_Cod_Juv.RData')
adt2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod_20210518/stacked_Spr_Cod_Adt.RData')

adtno=data.frame(matrix(NA, ncol=1, nrow=dim(adt1)[3])); colnames(adtno)='I'
for( i in 1:dim(adt1)[3]){
  adtno[i,1]=nicheOverlap(adt1[[i]], adt2[[i]], stat='I', mask=TRUE, checkNegatives=TRUE)
}
tadt=as.matrix(unlist(adtno))
plot(tadt~yrlist, type='b',ylim=c(.85,1), ylab='adult niche overlap', xlab='')

juvno=data.frame(matrix(NA, ncol=1, nrow=dim(juv1)[3])); colnames(juvno)='I'
for( i in 1:dim(juv1)[3]){
  juvno[i,1]=nicheOverlap(juv1[[i]], juv2[[i]], stat='I', mask=TRUE, checkNegatives=TRUE)
}
tjuv=as.matrix(unlist(juvno))
plot(tjuv~yrlist, type='b', ylim=c(.85,1), ylab='juvenile niche overlap', xlab='')

ichno=data.frame(matrix(NA, ncol=1, nrow=dim(ich1)[3])); colnames(ichno)='I'
for( i in 1:dim(ich1)[3]){
  ichno[i,1]=nicheOverlap(ich1[[i]], ich2[[i]], stat='I', mask=TRUE, checkNegatives=TRUE)
}
tich=as.matrix(unlist(ichno))
plot(tich~yrlist, type='b',ylim=c(.85,1), ylab='ichthyo niche overlap', xlab='')

plot(NULL, ylim=c(0.65,1), xlim=c(1977,2019), ylab="niche overlap", xlab="")
points(tadt~yrlist, type='b', col='red', pch=19)
points(tjuv~yrlist, type='b', col='blue', pch=19)
points(tich~yrlist, type='b', col='green', pch=19)
legend('topleft', legend = c('adt', 'juv', 'ich'), pch=19, col=c('red', 'blue', 'green'), bty='n', horiz = T)

### save data out
Spr.Niche.Overlap=data.frame(matrix(NA, ncol=4, nrow=dim(ich1)[3])); colnames(Spr.Niche.Overlap)=c('year', 'adt', 'juv', 'ich')
Spr.Niche.Overlap[,1]=yrlist
Spr.Niche.Overlap[,2]=adtno
Spr.Niche.Overlap[,3]=juvno
Spr.Niche.Overlap[,4]=ichno
save(Spr.Niche.Overlap, file='/home/ryan/Git/NEhabitat/habitat index/numbers/Spring_Niche_Overlap_Cod_Haddock_GSe_20210518.rda')

## 20210317 add new models
#fall GSe haddock
# juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock/stacked_Fall_Haddock_Juv.RData')
# adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock/stacked_Fall_Haddock_Adt.RData')
## 20210421 add new models for haddock
#fall GSe haddock
juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock_20210421/stacked_Fall_Haddock_Juv.RData')
adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock_20210421/stacked_Fall_Haddock_Adt.RData')
#fall GSe cod
# juv2=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod/stacked_Fall_Cod_Juv.RData')
# adt2=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod/stacked_Fall_Cod_Adt.RData')
#fall GSe cod updated from above (new data from Sean survdat)
juv2=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod_20210518/stacked_Fall_Cod_Juv.RData')
adt2=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod_20210518/stacked_Fall_Cod_Adt.RData')
# #fall GI haddock
# juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGI_Fall_haddock/stacked_Fall_Haddock_Juv.RData')
# adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGI_Fall_haddock/stacked_Fall_Haddock_Adt.RData')
# #fall GI cod
# juv2=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGI_Fall_cod/stacked_Fall_Cod_Juv.RData')
# adt2=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGI_Fall_cod/stacked_Fall_Cod_Adt.RData')

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

plot(NULL, ylim=c(0.65,1), xlim=c(1977,2019), ylab="niche overlap", xlab="", las=1)
points(t1~yrlist, type='b', col='red', pch=19)
points(t2~yrlist, type='b', col='blue', pch=19)
legend('topleft', legend = c('adt', 'juv'), pch=19, col=c('red', 'blue'), bty='n', horiz = T)

### save data out
Fall.Niche.Overlap=data.frame(matrix(NA, ncol=3, nrow=dim(adt1)[3])); colnames(Fall.Niche.Overlap)=c('year', 'adt', 'juv')
Fall.Niche.Overlap[,1]=yrlist
Fall.Niche.Overlap[,2]=adtno
Fall.Niche.Overlap[,3]=juvno
# Fall.Niche.Overlap[,4]=ichno
save(Fall.Niche.Overlap, file='/home/ryan/Git/NEhabitat/habitat index/numbers/Fall_Niche_Overlap_Cod_Haddock_GSe_20210518.rda')



### 2) Intra-species habitat overlap: between seasons by stage and species ###

## 20210317 add new models
# # fall GSe haddock
# juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock/stacked_Fall_Haddock_Juv.RData')
# adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock/stacked_Fall_Haddock_Adt.RData')
# #spring GSe haddock
# juv2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_haddock/stacked_Spr_Haddock_Juv.RData')
# adt2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_haddock/stacked_Spr_Haddock_Adt.RData')
#fall GSe haddock
juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock_20210421/stacked_Fall_Haddock_Juv.RData')
adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock_20210421/stacked_Fall_Haddock_Adt.RData')
#spring GSe haddock
juv2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/stacked_Spr_Haddock_Juv.RData')
adt2=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/stacked_Spr_Haddock_Adt.RData')

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
save(Haddock.Niche.Overlap, file='/home/ryan/Git/NEhabitat/habitat index/numbers/Haddock_Spr_Fall_Niche_Overlap_GSe_20210518.rda')


#fall GSe cod updated from above (new data from Sean survdat)
juv2=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod_20210518/stacked_Fall_Cod_Juv.RData')
adt2=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod_20210518/stacked_Fall_Cod_Adt.RData')
#spring GSe cod
juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod_20210518/stacked_Spr_Cod_Juv.RData')
adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod_20210518/stacked_Spr_Cod_Adt.RData')
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
save(Cod.Niche.Overlap, file='/home/ryan/Git/NEhabitat/habitat index/numbers/Cod_Spr_Fall_Niche_Overlap_GSe_20210518.rda')

### 3) Life stage habitat overlaps: Cod and Haddock by season: adt-juv, adt-ich, juv-ich

## Choose season and fish, load files and run (as appropriate with ich)fSEAS='Spr'
fSEAS='Spr' #Spr Fall
fFISH='Cod' # Cod Haddock
#spring GSe haddock
ich1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/stacked_Spr_Haddock_ich.RData')
juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/stacked_Spr_Haddock_Juv.RData')
adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/stacked_Spr_Haddock_Adt.RData')
#spring GSe cod
ich1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod_20210518/stacked_Spr_Cod_ich.RData')
juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod_20210518/stacked_Spr_Cod_Juv.RData')
adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod_20210518/stacked_Spr_Cod_Adt.RData')
#fall GSe haddock
juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock_20210421/stacked_Fall_Haddock_Juv.RData')
adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock_20210421/stacked_Fall_Haddock_Adt.RData')
#fall GSe cod updated from above (new data from Sean survdat)
juv1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod_20210518/stacked_Fall_Cod_Juv.RData')
adt1=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod_20210518/stacked_Fall_Cod_Adt.RData')
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
  plot(NULL, ylim=c(0.65,1), xlim=c(1977,2019), ylab="niche overlap", xlab="", las=1, main=paste(fSEAS, fFISH, sep=' '))
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
save(Haddock.Stage.Niche.Overlap, file='/home/ryan/Git/NEhabitat/habitat index/numbers/Haddock_Spr_stage_Niche_Overlap_GSe_20210518.rda')
# Fall Haddock
Haddock.Stage.Niche.Overlap=data.frame(matrix(NA, ncol=2, nrow=dim(adt1)[3])); colnames(Haddock.Stage.Niche.Overlap)=c('year', 'adtjuv')
Haddock.Stage.Niche.Overlap[,1]=yrlist
Haddock.Stage.Niche.Overlap[,2]=adtjuvno
save(Haddock.Stage.Niche.Overlap, file='/home/ryan/Git/NEhabitat/habitat index/numbers/Haddock_Fall_stage_Niche_Overlap_GSe_20210518.rda')
# Spring Cod
Cod.Stage.Niche.Overlap=data.frame(matrix(NA, ncol=4, nrow=dim(adt1)[3])); colnames(Cod.Stage.Niche.Overlap)=c('year', 'adtjuv', 'adtich', 'juvich')
Cod.Stage.Niche.Overlap[,1]=yrlist
Cod.Stage.Niche.Overlap[,2]=adtjuvno
Cod.Stage.Niche.Overlap[,3]=adtichno
Cod.Stage.Niche.Overlap[,4]=juvichno
save(Cod.Stage.Niche.Overlap, file='/home/ryan/Git/NEhabitat/habitat index/numbers/Cod_Spr_stage_Niche_Overlap_GSe_20210518.rda')
# Fall Cod
Cod.Stage.Niche.Overlap=data.frame(matrix(NA, ncol=2, nrow=dim(adt1)[3])); colnames(Cod.Stage.Niche.Overlap)=c('year', 'adtjuv')
Cod.Stage.Niche.Overlap[,1]=yrlist
Cod.Stage.Niche.Overlap[,2]=adtjuvno
save(Cod.Stage.Niche.Overlap, file='/home/ryan/Git/NEhabitat/habitat index/numbers/Cod_Fall_stage_Niche_Overlap_GSe_20210518.rda')


# PLOTTING SAVED DATA

### Niche Overlap between stages (species, season) USING NUMBERS model
# Spring Cod
load('/home/ryan/Git/NEhabitat/habitat index/numbers/Cod_Spr_stage_Niche_Overlap_GSe_20210518.rda')
plot(NULL, ylim=c(0.65,1), xlim=c(1977,2019), ylab="niche overlap", xlab="", las=1 ,main='Spr Cod')
points(Cod.Stage.Niche.Overlap$adtjuv~yrlist, type='b', col='red', pch=19)
points(Cod.Stage.Niche.Overlap$adtich~yrlist, type='b', col='blue', pch=19)
points(Cod.Stage.Niche.Overlap$juvich~yrlist, type='b', col='green', pch=19)
legend('topleft', legend = c('adt-juv','adt-ich', 'juv-ich'), pch=19, col=c('red', 'blue', 'green'), bty='n', horiz = T)
# Spr Haddock
load('/home/ryan/Git/NEhabitat/habitat index/numbers/Haddock_Spr_stage_Niche_Overlap_GSe_20210518.rda')
plot(NULL, ylim=c(0.65,1), xlim=c(1977,2019), ylab="niche overlap", xlab="", las=1 ,main='Spr Haddock')
points(Haddock.Stage.Niche.Overlap$adtjuv~yrlist, type='b', col='red', pch=19)
points(Haddock.Stage.Niche.Overlap$adtich~yrlist, type='b', col='blue', pch=19)
points(Haddock.Stage.Niche.Overlap$juvich~yrlist, type='b', col='green', pch=19)
legend('topleft', legend = c('adt-juv','adt-ich', 'juv-ich'), pch=19, col=c('red', 'blue', 'green'), bty='n', horiz = T)
# Fall (both)
load('/home/ryan/Git/NEhabitat/habitat index/numbers/Haddock_Fall_stage_Niche_Overlap_GSe_20210518.rda')
load('/home/ryan/Git/NEhabitat/habitat index/numbers/Cod_Fall_stage_Niche_Overlap_GSe_20210518.rda')
plot(NULL, ylim=c(0.85,1), xlim=c(1977,2019), ylab="niche overlap", xlab="", las=1 ,main='Fall Adt-Juv')
points(Cod.Stage.Niche.Overlap$adtjuv~yrlist, type='b', col='red', pch=19)
points(Haddock.Stage.Niche.Overlap$adtjuv~yrlist, type='b', col='blue', pch=19)
legend('topleft', legend = c('Cod','Haddock'), pch=19, col=c('red', 'blue'), bty='n', horiz = T)

### Niche Overlap between season (species, stage)
# Haddock
load('/home/ryan/Git/NEhabitat/habitat index/numbers/Haddock_Spr_Fall_Niche_Overlap_GSe_20210518.rda')
plot(NULL, ylim=c(0.75,1), xlim=c(1977,2019), ylab="niche overlap", xlab="", las=1 ,main='Haddock')
points(Haddock.Niche.Overlap$adt~yrlist, type='b', col='red', pch=19)
points(Haddock.Niche.Overlap$juv~yrlist, type='b', col='blue', pch=19)
legend('topleft', legend = c('Adt','Juv'), pch=19, col=c('red', 'blue'), bty='n', horiz = T)
# Cod
load('/home/ryan/Git/NEhabitat/habitat index/numbers/Cod_Spr_Fall_Niche_Overlap_GSe_20210518.rda')
plot(NULL, ylim=c(0.75,1), xlim=c(1977,2019), ylab="niche overlap", xlab="", las=1 ,main='Cod')
points(Cod.Niche.Overlap$adt~yrlist, type='b', col='red', pch=19)
points(Cod.Niche.Overlap$juv~yrlist, type='b', col='blue', pch=19)
legend('topleft', legend = c('Adt','Juv'), pch=19, col=c('red', 'blue'), bty='n', horiz = T)

### Niche Overlap between species (season, stage)
# Spr
load('/home/ryan/Git/NEhabitat/habitat index/numbers/Spring_Niche_Overlap_Cod_Haddock_GSe_20210518.rda')
plot(NULL, ylim=c(0.65,1), xlim=c(1977,2019), ylab="niche overlap", xlab="", las=1 ,main='Spr Cod-Haddock')
points(Spr.Niche.Overlap$adt~yrlist, type='b', col='red', pch=19)
points(Spr.Niche.Overlap$juv~yrlist, type='b', col='blue', pch=19)
points(Spr.Niche.Overlap$ich~yrlist, type='b', col='green', pch=19)
legend('topleft', legend = c('Adt', 'Juv','Ich'), pch=19, col=c('red', 'blue', 'green'), bty='n', horiz = T)
# Fall
load('/home/ryan/Git/NEhabitat/habitat index/numbers/Fall_Niche_Overlap_Cod_Haddock_GSe_20210518.rda')
plot(NULL, ylim=c(0.65,1), xlim=c(1977,2019), ylab="niche overlap", xlab="", las=1 ,main='Fall Cod-Haddock')
points(Fall.Niche.Overlap$adt~yrlist, type='b', col='red', pch=19)
points(Fall.Niche.Overlap$juv~yrlist, type='b', col='blue', pch=19)
legend('topleft', legend = c('Adt', 'Juv'), pch=19, col=c('red', 'blue'), bty='n', horiz = T)
#_________________________________________________________________
#### AS ABOVE BUT USING PA MODEL
### Niche Overlap between stages (species, season) USING NUMBERS model
# Spring Cod

par(oma=c(3,2,2,0))
par(mar=c(2,4,2,1) + 0.1)
# pdf(file=paste('/home/ryan/Git/NEhabitat/Niche_Overlap_plots_with_signif_slopes.pdf', sep=''))

load('/home/ryan/Git/NEhabitat/habitat index/Cod_Spr_stage_Niche_Overlap_GSe_20210518.rda')
plot(NULL, ylim=c(0.85,1), xlim=c(1977,2019), ylab="Stage niche overlap", xlab="", las=1)
lines(Cod.Stage.Niche.Overlap$adtjuv~yrlist, lty=1, col='black')
m = lm(Cod.Stage.Niche.Overlap$adtjuv~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=1)
  text(1980,0.95, paste('AdJv ', round(summary(m)$coefficients[8],4)))
}
lines(Cod.Stage.Niche.Overlap$adtich~yrlist, lty=2, col='black')
m = lm(Cod.Stage.Niche.Overlap$adtich~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=2)
  text(1980,0.91, paste('AdIc ', round(summary(m)$coefficients[8],4)))
}
lines(Cod.Stage.Niche.Overlap$juvich~yrlist, lty=3, col='black')
m = lm(Cod.Stage.Niche.Overlap$juvich~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=3)
  text(1980,0.87, paste('JvIc ', round(summary(m)$coefficients[8],4)))
}
legend('topleft', legend = c('Ad-Jv', 'Ad-Ic', 'Jv-Ic'), lty=c(1,2,3), col='black', bty='n', horiz = T) #col=c('red', 'blue', 'green')
text(2016,0.85, 'Spr cod')
# Spr Haddock
load('/home/ryan/Git/NEhabitat/habitat index/Haddock_Spr_stage_Niche_Overlap_GSe_20210518.rda')
plot(NULL, ylim=c(0.85,1), xlim=c(1977,2019), ylab="Stage niche overlap", xlab="", las=1)
lines(Haddock.Stage.Niche.Overlap$adtjuv~yrlist, lty=1, col='black')
m = lm(Haddock.Stage.Niche.Overlap$adtjuv~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=1)
  text(1980,0.97, paste('AdJv ', round(summary(m)$coefficients[8],4)))
}
lines(Haddock.Stage.Niche.Overlap$adtich~yrlist, lty=2, col='black')
m = lm(Haddock.Stage.Niche.Overlap$adtich~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=2)
  text(1980,0.85, paste('AdIc ', round(summary(m)$coefficients[8],4)))
}
lines(Haddock.Stage.Niche.Overlap$juvich~yrlist, lty=3, col='black')
m = lm(Haddock.Stage.Niche.Overlap$juvich~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=3)
  text(1980,0.92, paste('JvIc ', round(summary(m)$coefficients[8],4)))
}
legend('topleft', legend = c('Ad-Jv', 'Ad-Ic', 'Jv-Ic'), lty=c(1,2,3), col='black', bty='n', horiz = T)
text(2016,0.85, 'Spr had')
# Fall (both)
load('/home/ryan/Git/NEhabitat/habitat index/Haddock_Fall_stage_Niche_Overlap_GSe_20210518.rda')
load('/home/ryan/Git/NEhabitat/habitat index/Cod_Fall_stage_Niche_Overlap_GSe_20210518.rda')
plot(NULL, ylim=c(0.88,0.95), xlim=c(1977,2019), ylab="Stage niche overlap", xlab="", las=1)
lines(Cod.Stage.Niche.Overlap$adtjuv~yrlist, lty=1, col='black')
m = lm(Cod.Stage.Niche.Overlap$adtjuv~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=1)
  text(1980,0.93, paste('cod ', round(summary(m)$coefficients[8],4)))
}
lines(Haddock.Stage.Niche.Overlap$adtjuv~yrlist, lty=2, col='black')
m = lm(Haddock.Stage.Niche.Overlap$adtjuv~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=2)
  text(1980,0.9, paste('had ', round(summary(m)$coefficients[8],4)))
}
legend('topleft', legend = c('cod','had'), lty=c(1,2), col='black', bty='n', horiz = T)
text(2016,0.88, 'Fall Ad-Jv')

### Niche Overlap between season (species, stage)
# Haddock
load('/home/ryan/Git/NEhabitat/habitat index/Haddock_Spr_Fall_Niche_Overlap_GSe_20210518.rda')
plot(NULL, ylim=c(0.85,1), xlim=c(1977,2019), ylab="Seasonal niche overlap", xlab="", las=1)
lines(Haddock.Niche.Overlap$adt~yrlist, lty=1, col='black')
m = lm(Haddock.Niche.Overlap$adt~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=1)
  text(1980,0.87, paste('Ad ', round(summary(m)$coefficients[8],4)))
}
lines(Haddock.Niche.Overlap$juv~yrlist, lty=2, col='black')
m = lm(Haddock.Niche.Overlap$juv~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=2)
  text(1980,0.96, paste('Jv ', round(summary(m)$coefficients[8],4)))
}
legend('topleft', legend = c('Ad','Jv'), lty=c(1,2), col='black', bty='n', horiz = T)
text(2016,0.85, 'haddock')
# Cod
load('/home/ryan/Git/NEhabitat/habitat index/Cod_Spr_Fall_Niche_Overlap_GSe_20210518.rda')
plot(NULL, ylim=c(0.85,1), xlim=c(1977,2019), ylab="Seasonal niche overlap", xlab="", las=1)
lines(Cod.Niche.Overlap$adt~yrlist, lty=1, col='black')
m = lm(Cod.Niche.Overlap$adt~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=1)
  text(1980,0.96, paste('Ad ', round(summary(m)$coefficients[8],4)))
}
lines(Cod.Niche.Overlap$juv~yrlist, lty=2, col='black')
m = lm(Cod.Niche.Overlap$juv~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=2)
  text(1980,0.9, paste('Jv ', round(summary(m)$coefficients[8],5)))
}
legend('topleft', legend = c('Ad','Jv'), lty=c(1,2), col='black', bty='n', horiz = T)
text(2016,0.85, 'cod')

### Niche Overlap between species (season, stage)
# Spr
load('/home/ryan/Git/NEhabitat/habitat index/Spring_Niche_Overlap_Cod_Haddock_GSe_20210518.rda')
plot(NULL, ylim=c(0.8,1), xlim=c(1977,2019), ylab="Species niche overlap", xlab="", las=1)
lines(Spr.Niche.Overlap$adt~yrlist, lty=1, col='black')
m = lm(Spr.Niche.Overlap$adt~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=1)
  text(1980,0.9, paste('Ad ', round(summary(m)$coefficients[8],4)))
  }
lines(Spr.Niche.Overlap$juv~yrlist, lty=2, col='black')
m = lm(Spr.Niche.Overlap$juv~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=2)
  text(1980,0.82, paste('Jv ', round(summary(m)$coefficients[8],4)))
}
lines(Spr.Niche.Overlap$ich~yrlist, lty=3, col='black')
m = lm(Spr.Niche.Overlap$ich~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=3)
  text(1980,0.96, paste('Ic ', round(summary(m)$coefficients[8],4)))
  }
legend('topleft', legend = c('Ad', 'Jv','Ic'), lty=c(1,2,3), col='black', bty='n', horiz = T)
text(2016,0.8, 'Spr')
# Fall
load('/home/ryan/Git/NEhabitat/habitat index/Fall_Niche_Overlap_Cod_Haddock_GSe_20210518.rda')
plot(NULL, ylim=c(0.78,0.93), xlim=c(1977,2019), ylab="Species niche overlap", xlab="", las=1)
lines(Fall.Niche.Overlap$adt~yrlist, lty=1, col='black')
m = lm(Fall.Niche.Overlap$adt~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=1)
  text(1980,0.9, paste('Ad ', round(summary(m)$coefficients[8],5)))
}
lines(Fall.Niche.Overlap$juv~yrlist, lty=2, col='black')
m = lm(Fall.Niche.Overlap$juv~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=2)
  text(1980,0.82, paste('Jv ', round(summary(m)$coefficients[8],5)))
}
legend('topleft', legend = c('Ad', 'Jv'), lty=c(1,2), col='black', bty='n', horiz = T)
text(2016,0.78, 'Fall')
# dev.off()
