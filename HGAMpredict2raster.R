### set up year list to match files with
yrlist=seq(from=1977, to=2016, by=1)

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

### SELECT SPR of FALL seasons to process
SEASON='Spr' # Fall
SEASON='Fall'
usemodel=fish_modG4bll_fall_had

## list data files in each folder
btlist=list.files(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/BT', sep=''))
stlist=list.files(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/ST', sep=''))
zlist=list.files(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/pseudo', sep=''))
## parse year from filenames e.g #"RAST_NESREG_1977.04.03.BT.TEMP.YEAR.000066596.RData"
tb=strsplit(btlist, split=('RAST_NESREG_'))
ttb=sapply(tb, function(x) strsplit(x, "[.]")[[2]][1], USE.NAMES=FALSE)
ttb2=as.numeric(ttb)
#Surface temp
ts=strsplit(btlist, split=('RAST_NESREG_'))
tts=sapply(ts, function(x) strsplit(x, "[.]")[[2]][1], USE.NAMES=FALSE)
tts2=as.numeric(tts)
# Zooplankton (pseudocal)
tz=strsplit(btlist, split=('RAST_NESREG_'))
ttz=sapply(tz, function(x) strsplit(x, "[.]")[[2]][1], USE.NAMES=FALSE)
ttz2=as.numeric(ttz)

### load bottom temp file to use as template for non conforming rasters (remove before load dynamic files)
load('/home/ryan/Git/NEhabitat/rasters/test/RAST_NESREG_1977.04.03.BT.TEMP.YEAR.000066596.RData') #BT
bt=masked.raster
## depth raster
load('/home/ryan/Git/NEhabitat/rasters/test/depth.RData') #gdepth
gz=resample(gdepth, bt, 'bilinear')
gz2=(gz*-1)
gd2=crop(gz2, bt)
gd2=mask(gd2, bt)
## grain size raster
load('/home/ryan/Git/NEhabitat/rasters/test/grainsizeMM.RData') #phi2mm
phi2=resample(phi2mm, bt)
phi2=crop(phi2, bt)
phi2=mask(phi2, bt)
## rugosity raster
load('/home/ryan/Git/NEhabitat/rasters/test/scaledrugosity.RData') #rugscl
# load('/home/ryan/Git/NEhabitat/rast_rugosity.rdata')
# rug=masked.raster
# mmin=cellStats(rug, 'min')
# mmax=cellStats(rug, 'max')
# rugscl=calc(rug, fun=function(x){(x-mmin)/(mmax-mmin)}) # rescale from -4:2 -> 0:1
rug2=resample(rugscl, bt)
ex=extent(bt)
rug2=crop(rug2, ex)
rug2=mask(rug2, bt)

## Remove bottom temp raster
rm(bt)
fl=levels=c("Adt", "Juv", "ich")
fishnm='Haddock'
wd2=paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/', fishnm, '/', sep='')
### NOW loop over files, load yearly dynamic raster files and predict habitat from HGAM models
for (jj in 1:3){
  for (i in 1:length(yrlist)){
    bi=which(yrlist[i]==ttb2) # index of year
    load(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/BT/', btlist[[bi]], sep=''))
    bt=masked.raster # rename
    bi=which(yrlist[i]==tts2) # index of year
    load(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/ST/', stlist[[bi]], sep=''))
    st=masked.raster
    bi=which(yrlist[i]==ttz2) # index of year
    load(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/pseudo/', zlist[[bi]], sep=''))
    pse=masked.raster
    ef <- data.frame(coordinates(bt), val=values(bt))
    colnames(ef)=c("LON", "LAT", "BOTTEMP")
    ef$SURFTEMP=values(st)
    ef$DEPTH=values(gd2)
    ef$Stg=factor(fl[jj], levels=c('Adt', 'Juv', 'ich'))
    ef$pseudo_100m3=values(pse)
    ef$grnszmm=values(phi2)
    ef$rug=values(rug2)
    ef2=ef[complete.cases(ef),]
    # tt=ef2$Stg
    # tt2=fl[[1]][[1]][2][[1]][tt] # subsets to stage 
    # ef$Stg=factor(ef$Stg, levels=c("Adt", "Juv", "ich"))
    test1 <- predict.gam(usemodel, ef2, type='response')
    ef2$pred=test1
    wd3=paste(yrlist[i], '_', SEASON, '_', fishnm, '_',fl[jj], '.RData', sep="")
    save(ef2, file=paste(wd2, wd3, sep=""))
    spg1=ef2[,c('LON', 'LAT', 'pred')]
    wd4=paste(yrlist[i], '_', 'RASTERpred', '_', SEASON, '_', fishnm, '_',fl[jj], '_', '.RData', sep="")
    # tes1=rasterFromXYZ(spg1[complete.cases(spg1$Stg),])
    # save(tes1, file=paste(wd2,wd4, sep=''))
    coordinates(spg1)= ~ LON + LAT
    gridded(spg1)=T
    rastDF=raster(spg1)
    # plot(rastDF)
    save(rastDF, file=paste(wd2,wd4, sep=''))
  }
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

# wd3='/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modG4bll_fall_had/Rasters/'
load(paste(wd4,'/',rastlistA[1], sep=''))
AdtHad=rastDF
for (i in 2:length(rastlistA)){
  load(paste(wd4,'/',rastlistA[i], sep=''))
  AdtHad=stack(AdtHad, rastDF)
}
load(paste(wd4,'/',rastlistJ[1], sep=''))
JuvHad=rastDF
for (i in 2:length(rastlistA)){
  load(paste(wd4,'/',rastlistJ[i], sep=''))
  JuvHad=stack(JuvHad, rastDF)
}
load(paste(wd4,'/',rastlistI[1], sep=''))
IchHad=rastDF
for (i in 2:length(rastlistI)){
  load(paste(wd4,'/',rastlistI[i], sep=''))
  IchHad=stack(IchHad, rastDF)
}
names(IchHad)=yrlist
names(AdtHad)=yrlist
names(JuvHad)=yrlist
save(IchHad, file=paste(wd4, 'tacks/IchHad.Rdata', sep='')) 
save(AdtHad, file=paste(wd4, 'tacks/AdtHad.Rdata', sep=''))
save(JuvHad, file=paste(wd4, 'tacks/JuvHad.Rdata', sep=''))


## this is not working, needs to call name of saved raster 'rastDF'
# loadNstack=function(x, namex){
#   namex=load(x[1])
#   for (i in 2:length(x)){
#     nameb=load(x[i])
#     namex=stack(namex, nameb)
#   }
#   return(namex)
# }


##### _______________________________________________________________________________________________________
# Correlation routine
#####
source("/home/ryan/Desktop/8 R Functions/Rfuncs.R")
source("/home/ryan/Desktop/8 R Functions/Gridcorts.R")
source("/home/ryan/Desktop/8 R Functions/KDE_funcs.R")
library(RColorBrewer)
library(colorspace)
cg1=gridcorts(rasterstack=AdtHad, method='kendall',type='corel') #spearman
cg2=gridcorts(rasterstack=AdtHad, method='kendall',type='pval')
cg3=cg1  
plot(cg3, col=diverge_hcl(120),main=paste(spp,'\nCorrelation to', season, datalab), las=1)
plot(nesbath,deep=-50, shallow=-50, step=1,add=T,lwd=1,col='gray40',lty=1)
plot(nesbath,deep=-100, shallow=-100, step=1,add=T,lwd=1,col='gray60',lty=1)
plot(nesbath,deep=-200, shallow=-200, step=1,add=T,lwd=1,col='gray60',lty=1)
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="black", add=T)

#mark significant correlations
plot(cg3, col=diverge_hcl(120),main=paste(spp,'\nCorrelation to', season, datalab), las=1)
cg3[cg2>0.05]=NA
test=rasterToPoints(cg3)
points(test, pch=12, col=addTrans('black',50), cex=0.5)
plot(nesbath,deep=-50, shallow=-50, step=1,add=T,lwd=1,col='gray40',lty=1)
plot(nesbath,deep=-100, shallow=-100, step=1,add=T,lwd=1,col='gray60',lty=1)
plot(nesbath,deep=-200, shallow=-200, step=1,add=T,lwd=1,col='gray60',lty=1)
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="black", add=T)

#####
# Linear model of environmental variable change, get slope and stats
# drop first layer for spring wind stress
#####
shp.dat3=crop(shp.dat2, small)
time <- 1:nlayers(AdtHad) 
fun=function(x) { if (is.na(x[1])){ NA } else { m = lm(x ~ time); summary(m)$coefficients[2] }}
shp.dat.slope=calc(IchHad, fun)
mn=min(shp.dat.slope@data@values, na.rm = T)
mx=max(shp.dat.slope@data@values, na.rm = T)
high=max(abs(mn), mx)
br <- seq(-high, high, by = high/15) 
cl <- diverge_hcl(length(br) - 1, power = 1) 
rng=range(shp.dat.slope[],na.rm=T)
arg=list(at=rng, labels=round(rng,3))
# plot(shp.dat.slope, main=paste(datalab, season, length(time),'yrs','\nYearly Slope'),col=diverge_hcl(120), las=1) # Yearly slope
plot(shp.dat.slope, col=cl, breaks=br,axis.args=arg,las=1) # Yearly slope
plot(nesbath,deep=-50, shallow=-50, step=1,add=T,lwd=1,col='gray40',lty=1)
plot(nesbath,deep=-200, shallow=-200, step=1,add=T,lwd=1,col='gray60',lty=1)
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="black", add=T)

## Calc Significance
fun=function(x) { if (is.na(x[1])){ NA } else { m = lm(x ~ time); summary(m)$coefficients[8] }}
p <- calc(IchHad, fun=fun)
m = c(0, 0.05, 1, 0.05, 1, 0)
rclmat = matrix(m, ncol=3, byrow=TRUE)
p.mask = reclassify(p, rclmat)
fun=function(x) { x[x<1] <- NA; return(x)}
p.mask.NA = calc(p.mask, fun)
trend.sig = mask(shp.dat.slope, p.mask.NA)
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
