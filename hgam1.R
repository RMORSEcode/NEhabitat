#### Hierarchical GAMs fitting life history stages for NE groundfish
# RM 2020
# 1) Train model w/ data through 1999, predict 2000-2018
# 2) Train 1977-2011, predict 2012-2018
# 3) Train 1977-2012, predict 2013-2018
# 4) Train 1977-2014, predict 2015-2018

library(mgcv)
library(gratia)
library(car)
library(MASS)
library(stringr)
# library(gamm4)
library(tidyr)
# library(ggplot2)
# library(ggthemes)
# library(viridis)
# library(cowplot)
# library(kableExtra)
# library(docxtools)
# library(knitr)
# library(tibble)
library(dplyr)
# library(gratia)
# library(latex2exp)
library(lubridate)
library(raster)


rqresiduals = function (gam.obj) {
  if(!"gam" %in% attr(gam.obj,"class")){
    stop('"gam.obj has to be of class "gam"')
  }
  if (!grepl("^Tweedie|^Negative Binomial|^poisson|^binomial|^gaussian|^Gamma|^inverse.gaussian",
             gam.obj$family$family)){
    stop(paste("family " , gam.obj$family$family, 
               " is not currently supported by the statmod library, 
               and any randomized quantile residuals would be inaccurate.",
               sep=""))
  }
  if (grepl("^Tweedie", gam.obj$family$family)) {
    if (is.null(environment(gam.obj$family$variance)$p)) {
      p.val <- gam.obj$family$getTheta(TRUE)
      environment(gam.obj$family$variance)$p <- p.val
    }
    qres <- qres.tweedie(gam.obj)
  }
  else if (grepl("^Negative Binomial", gam.obj$family$family)) {
    if ("extended.family" %in% class(gam.obj$family)) {
      gam.obj$theta <- gam.obj$family$getTheta(TRUE)
    }
    else {
      gam.obj$theta <- gam.obj$family$getTheta()
    }
    qres <- qres.nbinom(gam.obj)
  }
  else {
    qres <- qresid(gam.obj)
  }
  return(qres)
}

load('/home/ryan/Git/NEhabitat/Final_merged_fish_corABN_Zoo_Ich.Rda') # Biomass
load('/home/ryan/Git/NEhabitat/Final_merged_fish_corABN_Zoo_Ich.Rda') # Abundance
### update colnames for ichythoplankton ### Note: urospp assigned to 76 (WHK) but is not distinguishable from 77 (RHK)
# colnames(FData.bio[50:61])
# [1] "urospp_100m3" "gadmor_100m3" "melaeg_100m3" "polvir_100m3" "merbil_100m3" "sebspp_100m3"
# [7] "anaspp_100m3" "parden_100m3" "pseame_100m3" "glycyn_100m3" "scoaqu_100m3" "lopame_100m3"
colnames(FData.bio)[50:61]=c("76_ich", "73_ich", "74_ich", "75_ich", "72_ich", "155_ich", "192_ich", "103_ich","106_ich","107_ich", "108_ich", "197_ich")
colnames(FData.abn)[50:61]=c("76_ich", "73_ich", "74_ich", "75_ich", "72_ich", "155_ich", "192_ich", "103_ich","106_ich","107_ich", "108_ich", "197_ich")


### LOAD RASTER5 and extract static vars

# /home/ryan/Git/NEhabitat/rast_rugosity.rdata

### extract from rasters at fish points
locs=FData.abn[,c('LON', 'LAT')]
coordinates(locs) <- ~LON+LAT
mypoints = SpatialPoints(locs,proj4string = CRS("+init=epsg:4326"))
# myproj = CRS(masked.raster)
# points.proj = spTransform(mypoints, myproj)
load('/home/ryan/Git/NEhabitat/rast_phi_fraction.rdata') 
## convert phi (with negative vals) size to mm
phi=masked.raster
phi2mm=calc(phi, fun=function(x){2^-x})
plot(phi2mm)
FData.abn$grnszmm = extract(phi2mm, mypoints)
load('/home/ryan/Git/NEhabitat/rast_mud_fraction.rdata')
FData.abn$mud_pct = extract(masked.raster, mypoints)
load('/home/ryan/Git/NEhabitat/rast_sand_fraction.rdata')
FData.abn$sand_pct = extract(masked.raster, mypoints)
# load('/home/ryan/Git/NEhabitat/rast_gdepth.rdata')
# FData.abn$gdepth = extract(masked.raster, mypoints) # odd scale neg vals start at 1000m
load('/home/ryan/Git/NEhabitat/rast_rugosity.rdata')
rug=masked.raster
mmin=cellStats(rug, 'min')
mmax=cellStats(rug, 'max')
rugscl=calc(rug, fun=function(x){(x-mmin)/(mmax-mmin)}) # rescale from -4:2 -> 0:1
FData.abn$rug = extract(masked.raster, mypoints)
plot(masked.raster, xlim=c(-80,-60), ylim=c(36,48))
plot(rugscl)
FData.abn$rug = extract(rugscl, mypoints)

### select data
# fish=FData.abn %>% dplyr::select(YEAR, SEASON, LAT:BOTSALIN, `73_Adt`, `73_Juv`, `73_ich`, calfin_100m3:chaeto_100m3) # need to fill Salinity
fish=FData.abn %>% dplyr::select(YEAR, SEASON, LAT:BOTTEMP, `73_Adt`, `73_Juv`, `73_ich`, calfin_100m3:grnszmm)
fish$MONTH=month(FData.abn$EST_TOWDATE)
fish=fish[complete.cases(fish),]
fish$`73_ich`=ceiling(fish$`73_ich`) # make integer from numbers per 100 m/3
# fish[`73_ich`]=round(fish[`73_ich`], digts=0)
table(fish$YEAR)

## fitting distribution (not working...)
qqp(fish$`73_Adt`, "norm")
qqp(log(fish$`73_Adt`+1), "norm")

qqp(fish$`73_Adt`, "lnorm")
nbinom <- fitdistr(fish$`73_Adt`, "Negative Binomial")
qqp(FData.bio$`73_Adt`, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])
qqnorm(log(fish$`73_Adt`+1)); qqline(log(fish$`73_Adt`+1))

### Need to make data frame wide -> long with just species of interest...
# e.g.: df %>% pivot_longer(c(x, y, z), names_to = "key", values_to = "value")
# pivot_longer(data, cols, names_to = "name", names_prefix = NULL,
#              names_sep = NULL, names_pattern = NULL, names_ptypes = list(),
#              names_repair = "check_unique", values_to = "value",
#              values_drop_na = FALSE, values_ptypes = list())
fish2=fish %>% pivot_longer(c(`73_Adt`, `73_Juv`, `73_ich`), names_to = c("SVSPP", "Stg"), names_sep ="_", values_to = "Number")
fish2$Stg=factor(fish2$Stg, ordered=F)
table(fish2$MONTH)

### try negative binomial abd tweedie for counts

# Global smoother
fishG=gam(log10(fish2$Number+1)~te(fish2$MONTH, fish2$BOTTEMP, log10(fish2$calfin_100m3+1), bs=c("cc", "tp", "tp"), 
                                  k=c(8, 10, 5), m=2), family=Gamma(link="log"), method="REML")

fishG=gam(Number~te(MONTH, BOTTEMP, calfin_100m3, bs=c("cc", "tp", "tp"), k=c(8, 10, 5), m=2), data=fish2, family=nb(), method="REML")



# Error in eval(family$initialize) : 
  # non-positive values not allowed for the 'gamma' family

fishG=gam(fish2$Number~te(fish2$MONTH, fish2$BOTTEMP, log10(fish2$calfin_100m3+1), bs=c("cc", "tp", "tp"), 
                                   k=c(8, 10, 5), m=2), family=gaussian(link="identity"), method="REML")
# stage was not included explicitly in above, does this make it GI model?
fishGI=gam(log10(fish2$Number+1)~te(fish2$MONTH, fish2$BOTTEMP, fish2$Stg, bs=c("cc", "tp", "re"), 
                          k=c(8, 10, 5), m=2) + s(log10(fish2$calfin_100m3+1), bs="tp"), family=gaussian(link="identity"), method="REML")

fishGI=gam(fish2$Number~te(fish2$MONTH, fish2$BOTTEMP, bs=c("cc", "tp"), k=c(8, 10)) + s(fish2$Stg, bs="re", k=3) +
             s(log10(fish2$calfin_100m3+1), bs="tp", k=5), family=gaussian(link="identity"), method="REML")


plot(fishGI)
draw(fishGI)
gam.check(fishGI) #p values are for the test of the null hypothesis that the basis dimension used is of sufficient size
summary.gam(fishGI) #test null hypothesis of a zero effect of the indicated spline
vis.gam(fishGI)
plot.gam(fishGI)

fishGS=gam(log10(fish2$Number+1)~te(fish2$MONTH, fish2$BOTTEMP, fish2$Stg, bs=c("cc", "tp", "tp"), 
                                   k=c(8, 10, 3), m=2), method="REML")
fishGS2=gam(log10(fish2$Number+1)~te(fish2$MONTH, fish2$BOTTEMP, bs=c("cc", "tp"), k=c(8, 10)) + 
             s(fish2$Stg, bs="tp", k=3), method="REML")

### updated models below -> working from models.r in pederson folder on desktop
fish_modG = gam(Number ~s(BOTTEMP, k=100, bs="tp") + s(DEPTH, k=60, bs="tp") + s(Stg, k=4, bs='re') 
                + s(pseudo_100m3, k=6, bs="cr"), data=fish4, method = "REML", family="nb")

fish_modG2 = gam(Number ~s(BOTTEMP, k=100, bs="tp") + s(DEPTH, k=60, bs="tp") + s(Stg, k=3, bs='re') 
                 + s(pseudo_100m3, k=6, bs="cr"), data=fish4, method = "REML", family="tw")


### biomod2 attempt ####
library(biomod2)
resp.var=fish2$Number
resp.var[resp.var > 0]=1 # presence absesnce
latlon=as.matrix((data.frame(fish2$LON, fish2$LAT)))
expl.var=fish2 %>% select(DEPTH, SURFTEMP, BOTTEMP, calfin_100m3:chaeto_100m3)
expl.var=data.frame(expl.var)
myBiomodData=BIOMOD_FormatingData(resp.var,
                     expl.var,
                     resp.xy = latlon,
                     resp.name = "Cod",
                     eval.resp.var = NULL,
                     eval.expl.var = NULL,
                     eval.resp.xy = NULL,
                     PA.nb.rep = 0,
                     PA.nb.absences = 1000,
                     PA.strategy = 'random',
                     PA.dist.min = 0,
                     PA.dist.max = NULL,
                     PA.sre.quant = 0.025,
                     PA.table = NULL,
                     na.rm = TRUE)

myBiomodData
plot(myBiomodData)
myBiomodOption <- BIOMOD_ModelingOptions()

myBiomodModelOut <- BIOMOD_Modeling( 
  myBiomodData, 
  models = c('SRE','CTA','RF','MARS','FDA'), 
  models.options = myBiomodOption, 
  NbRunEval=3, 
  DataSplit=80, 
  Prevalence=0.5, 
  VarImport=3,
  models.eval.meth = c('TSS','ROC'),
  SaveObj = TRUE,
  rescal.all.models = TRUE,
  do.full.models = FALSE,
  modeling.id = paste('Cod',"FirstModeling",sep=""))

myBiomodModelOut 
myBiomodModelEval <- get_evaluations(myBiomodModelOut)
dimnames(myBiomodModelEval)
myBiomodModelEval["TSS","Testing.data","RF",,]
myBiomodModelEval["TSS","Testing.data",,,]

myBiomodModelEval["ROC","Testing.data",,,]
get_variables_importance(myBiomodModelOut)

### Plot model performance:
### by models
gg1 <- models_scores_graph( myBiomodModelOut,
                            by = 'models',
                            metrics = c('ROC','TSS') )

### by cross validation run
gg2 <- models_scores_graph( myBiomodModelOut,
                            by = 'cv_run',
                            metrics = c('ROC','TSS') )





## try fitting gams
y=log10(df$merbil_100m3+1) 
x=log10(dfz$calfin_100m3+1) #log Calfin
Sample_data <- data.frame(y,x)
g1=gam(y~ s(x), method="REML")
g1=gam(y~ s(x)+s(dfv$month, bs='cc', k=10) +s(dfv$btm_temp), method="REML")

g1=gam(log10(df$melaeg_100m3+1)~ s(log10(dfz$calfin_100m3+1))+s(dfv$month, bs='cc', k=12) +s(dfv$btm_temp) + s(log10(dfv$depth+1)), method="REML")

ggplot(Sample_data, aes(x, y)) + geom_point() + geom_smooth(method = "gam", formula = y ~s(x))

y=log10(df$gadmor_100m3+1) # log fish sp
z1=log10(dfz$calfin_100m3+1) #log Calfin
z2=log10(dfz$pseudo_100m3+1)
z3=log10(dfz$mlucens_100m3+1)
z4=log10(dfz$tlong_100m3+1)
z5=log10(dfz$cham_100m3+1)
z6=log10(dfz$ctyp_100m3+1)
z7=log10(dfz$calminor_100m3+1)
# Sample_data <- data.frame(y,x)
# g1=gam(df$melaeg_100m3 ~ s(dfz$calfin_100m3), method="REML")
g1=gam(y~ s(z1), method="REML")
g1=gam(y~ s(z1)+s(dfv$sfc_temp), method="REML")
g1=gam(y~ s(z1)+s(dfv$sfc_temp)+s(dfv$sfc_salt), method="REML")
g1=gam(y~ s(z7)+s(dfv$month, bs='cc', k=10) +s(dfv$btm_temp), method="REML")

g1=gam(y~ s(z1)+s(z2)+s(z3)+s(dfv$sfc_temp), method="REML")

g1=gam(y~ s(z4)+s(z5)+s(z6)+s(dfv$sfc_temp), method="REML")

g1=gam(y~ +s(z2, k=10)+s(dfv$lon, dfv$lat, k=48)+s(dfv$sfc_temp, k=20), family=nb, method="REML")
g1=gam(y~ +s(z2, k=10)+s(dfv$lon, dfv$lat, k=48)+s(dfv$sfc_temp, k=20), family="poisson", method="REML")

# ggplot(Sample_data, aes(x, y)) + geom_point() + geom_smooth(method = "gam", formula = y ~s(x))


plot(g1)
gam.check(g1) #p values are for the test of the null hypothesis that the basis dimension used is of sufficient size
summary.gam(g1) #test null hypothesis of a zero effect of the indicated spline
# par(mfrow = c(2,2))



