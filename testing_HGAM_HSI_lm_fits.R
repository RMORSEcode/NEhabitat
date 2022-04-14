### Testing whether linear models are appropriate for determining change in HSI hindcasts per HGAM paper review 20220327

## this does not work as expected.... try next step
# resfun = function(x) { if (is.na(x[1])){ NA } else { m <- lm(x ~ time, na.action=na.exclude); r <- residuals.lm(m);return (r)}}
# newrast=calc(rastck, resfun)
# mn=cellStats(newrast, min)
# mx=cellStats(newrast, max)
# high=max(abs(mn), mx)
# br <- seq(-high, high, by = high/15) 
# cl <- colorspace::diverge_hcl(length(br) - 1, power = 1) 
# rng=range(newrast[],na.rm=T)
# arg=list(at=rng, labels=round(rng,3))
# plot(newrast, col=cl, breaks=br,axis.args=arg,las=1, main=paste(length(time),'yrs','\nYearly Slope'), legend=F) # Yearly slope
# plot(newrast, legend.only=TRUE, col=cl, breaks=br, axis.args=arg, las=1, legend.width=1, legend.shrink=0.75, smallplot=c(0.6, 0.61, 0.2, 0.5)); par(mar = par("mar"))
# maps::map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="black", add=T)


## nor does this
# https://stackoverflow.com/questions/46786053/extracting-residuals-from-pixel-by-pixel-regression
# residualfun <- function(x) {
#   nly=dim(x)[3]
#   r <- rep(NA, nly)
#   obs <- x[1:nly]
#   cov <- x[nly + 1:nly]
#   # Remove NA values before model
#   x.nona <- which(!is.na(obs) & !is.na(cov))
#   # If more than 2 points proceed to lm
#   if (length(x.nona) > 2) {
#     m <- NA
#     try(m <- lm(obs[x.nona] ~ cov[x.nona]))
#     # If model worked, calculate residuals
#     if (is(m)[1] == "lm") {
#       r[x.nona] <- residuals.lm(m)
#     } else {
#       # alternate value to find where model did not work
#       r[x.nona] <- -1e32
#     }
#   }
#   return(r)
# }
# 
# res <- calc(adtpa, residualfun)
# res <- calc(adtpa, resfun)


### raster to matrix of vectors, do lm stuff, back to raster and stack -> visualize spatial residuals from lm
## choose data to use from 
# SPRING HADDOCK
obs=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/PA_only_stacked_Spr_Haddock_Adt.RData');SEASON='Spr'; fishnm='Haddock'
obs=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/PA_only_stacked_Spr_Haddock_Juv.RData');SEASON='Spr'; fishnm='Haddock'
obs=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/PA_only_stacked_Spr_Haddock_ich.RData');SEASON='Spr'; fishnm='Haddock'
# adt=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/stacked_Spr_Haddock_Adt.RData')
# juv=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/stacked_Spr_Haddock_Juv.RData')
# ich=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/stacked_Spr_Haddock_ich.RData')
# ;SEASON='Spr'; fishnm='Haddock'

# SPRING COD
obs=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod_20210518/PA_only_stacked_Spr_Cod_Adt.RData');SEASON='Spr';fishnm='Cod'
obs=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod_20210518/PA_only_stacked_Spr_Cod_Juv.RData');SEASON='Spr';fishnm='Cod'
obs=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod_20210518/PA_only_stacked_Spr_Cod_ich.RData');SEASON='Spr';fishnm='Cod'
# adt=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod_20210518/stacked_Spr_Cod_Adt.RData')
# juv=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod_20210518/stacked_Spr_Cod_Juv.RData')
# ich=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod_20210518/stacked_Spr_Cod_ich.RData')
# ;SEASON='Spr';fishnm='Cod'

# FALL HADDOCK
obs=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock_20210421/PA_only_stacked_Fall_Haddock_Adt.RData');SEASON='Fall';fishnm='Haddock'
obs=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock_20210421/PA_only_stacked_Fall_Haddock_Juv.RData');SEASON='Fall';fishnm='Haddock'
# adt=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock_20210421/stacked_Fall_Haddock_Adt.RData')
# juv=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock_20210421/stacked_Fall_Haddock_Juv.RData')
# ;SEASON='Fall';fishnm='Haddock' #74

# FALL COD
obs=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod_20210518/PA_only_stacked_Fall_Cod_Adt.RData');SEASON='Fall';fishnm='Cod'
obs=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod_20210518/PA_only_stacked_Fall_Cod_Juv.RData');SEASON='Fall';fishnm='Cod'
# adt=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod_20210518/stacked_Fall_Cod_Adt.RData')
# juv=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod_20210518/stacked_Fall_Cod_Juv.RData')
# ;SEASON='Fall';fishnm='Cod'


# obs=adtpa
crs(obs)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
x=getValues(obs)
dim(x)
x.res=matrix(data = NA, nrow = nrow(x), ncol = ncol(x), byrow = FALSE,
             dimnames = NULL)
x.fit=matrix(data = NA, nrow = nrow(x), ncol = ncol(x), byrow = FALSE,
             dimnames = NULL)
## create residuals in form of vector
for(i in 1:nrow(x)){
# x1=x[i,]
  if (is.na(x[i,1])){
    next
  }
lmx=lm(as.numeric(x[i,])~time, na.action = na.exclude)
x.res[i,]=resid(lmx)
x.fit[i,]=fitted(lmx)
}
xco=coordinates(obs[[1]])
#  convert back to raster using rasterFromXYZ
x1rm=matrix(data=c(xco[,1],xco[,2],x.res[,1]), ncol=3, nrow=9450)
x.resid=rasterFromXYZ(x1rm, res=c(0.1,0.1), crs=xcr)
for(i in 2:ncol(x)){
  x1rm=matrix(data=c(xco[,1],xco[,2],x.res[,i]), ncol=3, nrow=9450)
  x1r2=rasterFromXYZ(x1rm, res=c(0.1,0.1), crs=xcr)
  x.resid=stack(x.resid, x1r2)
}

#plot all years in one pixel
sum(!is.na(x.fit[,1]))
j=seq(1,nrow(x.fit),1)[!is.na(x.fit[,1])]
# i=j[1000:1010]
i=j

## plot all pixels in all years
plot(x.res[i,]~x.fit[i,], ylim=c(-0.2,0.2)); abline(0,0)
# plot all pixels in a year (1-43)
i=15
plot(x.res[,i]~x.fit[,i], ylim=c(-0.2,0.2)); abline(0,0)

barplot(table(round(x.fit[j,],2)))
quantile(x.fit, na.rm=T)

### plot histogram of residuals
hist(x.res[j,])
### plot histogram of model fit values, and quantiles with 50% split for high qualtity HSI
hist(x.fit[j,]); abline(v=0.5, col='red'); abline(v=quantile(x.fit[j,], na.rm=T), col='blue')


## plot spatial residuals
# plot(x.resid[[1]])
wd="/home/ryan/Git/NEhabitat"
pdf(file=paste(wd,SEASON,fishnm,'Yearly_Spatial_Residuals.pdf', sep='_'))
for(i in 1:dim(x.resid)[[3]]){
newrast=x.resid[[i]]
mn=cellStats(newrast, min)
mx=cellStats(newrast, max)
high=max(abs(mn), mx)
br <- seq(mn, mx, by = high/15) 
cl <- colorspace::diverge_hcl(length(br) - 1, power = 1) 
rng=range(newrast[],na.rm=T)
arg=list(at=rng, labels=round(rng,3))
plot(newrast, col=cl, breaks=br,axis.args=arg,las=1, main=paste(yrlist[i], 'spatial residual'))
maps::map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="black", add=T)
}
dev.off()
