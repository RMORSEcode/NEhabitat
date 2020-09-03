library(biomod2)
library(tidyr)
library(RColorBrewer)
library(dismo)
library(dplyr)
library(lubridate)
library(raster)
load('/home/ryan/Documents/Git/NEhabitat/Final_merged_fish_corBIO_Zoo_Ich.Rda') # Biomass
load('/home/ryan/Documents/Git/NEhabitat/Final_merged_fish_corABN_Zoo_Ich.Rda') # Abundance

colnames(FData.bio)[50:61]=c("76_ich", "73_ich", "74_ich", "75_ich", "72_ich", "155_ich", "192_ich", "103_ich","106_ich","107_ich", "108_ich", "197_ich")
colnames(FData.abn)[50:61]=c("76_ich", "73_ich", "74_ich", "75_ich", "72_ich", "155_ich", "192_ich", "103_ich","106_ich","107_ich", "108_ich", "197_ich")

setwd('/home/ryan/Biomod models')

# Add variables from rasters to both biomass and numbers database
### extract from rasters at fish points
locs=FData.abn[,c('LON', 'LAT')]
coordinates(locs) <- ~LON+LAT
mypoints = SpatialPoints(locs,proj4string = CRS("+init=epsg:4326"))
# myproj = CRS(masked.raster)
# points.proj = spTransform(mypoints, myproj)
load('/home/ryan/Documents/Git/NEhabitat/rast_phi_fraction.rdata') 
## convert phi (with negative vals) size to mm
phi=masked.raster
phi2mm=calc(phi, fun=function(x){2^-x})
plot(phi2mm)
FData.abn$grnszmm = extract(phi2mm, mypoints)
FData.bio$grnszmm = extract(phi2mm, mypoints)
load('/home/ryan/Documents/Git/NEhabitat/rast_mud_fraction.rdata')
FData.abn$mud_pct = extract(masked.raster, mypoints)
FData.bio$mud_pct = extract(masked.raster, mypoints)
load('/home/ryan/Documents/Git/NEhabitat/rast_sand_fraction.rdata')
FData.abn$sand_pct = extract(masked.raster, mypoints)
FData.bio$sand_pct = extract(masked.raster, mypoints)
# load('/home/ryan/Git/NEhabitat/rast_gdepth.rdata')
# FData.abn$gdepth = extract(masked.raster, mypoints) # odd scale neg vals start at 1000m
load('/home/ryan/Documents/Git/NEhabitat/rast_rugosity.rdata')
rug=masked.raster
mmin=cellStats(rug, 'min')
mmax=cellStats(rug, 'max')
rugscl=calc(rug, fun=function(x){(x-mmin)/(mmax-mmin)}) # rescale from -4:2 -> 0:1
# FData.abn$rug = extract(masked.raster, mypoints)
# FData.bio$rug = extract(masked.raster, mypoints)
plot(masked.raster, xlim=c(-80,-60), ylim=c(36,48))
plot(rugscl)
FData.abn$rug = extract(rugscl, mypoints)
FData.bio$rug = extract(rugscl, mypoints)

# Select season for GAMs
#slctseason="SPRING" #"FALL"
slctseason="FALL"


# Subset the formatted data for a single species to use in models (need to update manually as these are column names)
### select data

fname='Haddock'
fish=FData.abn %>% dplyr::select(YEAR, SEASON, LAT:BOTTEMP, `74_Adt`, `74_Juv`, `74_ich`, calfin_100m3:rug)
fish$MONTH=month(FData.abn$EST_TOWDATE)
fish=fish[complete.cases(fish),]
fish$`74_ich`=ceiling(fish$`74_ich`) # make integer from numbers per 100 m/3

### Need to make data frame wide -> long with just species of interest...
# e.g.: df %>% pivot_longer(c(x, y, z), names_to = "key", values_to = "value")
# pivot_longer(data, cols, names_to = "name", names_prefix = NULL,
#              names_sep = NULL, names_pattern = NULL, names_ptypes = list(),
#              names_repair = "check_unique", values_to = "value",
#              values_drop_na = FALSE, values_ptypes = list())
fish2=fish %>% pivot_longer(c(`74_Adt`, `74_Juv`, `74_ich`), names_to = c("SVSPP", "Stg"), names_sep ="_", values_to = "Number")
fish2$Stg=factor(fish2$Stg, ordered=F)
table(fish2$MONTH)

## log transform plankton
fish3=fish2
logd=fish3[,8:17]
logd=log10(logd+1)
fish3[,8:17]=logd

# Now set up traing and testing data; split data sets
# 1) Train model w/ data through 1999, predict 2000-2018
# 2) Train 1977-2011, predict 2012-2018
# 3) Train 1977-2012, predict 2013-2018
# 4) Train 1977-2014, predict 2015-2018

# Select Season and Fish species
# Training/testing year split
dat.split.year<- 2018

## select season and split training/testing; Change to Presence-Absence for fish Numbers
fish5=fish3[which(fish3$SEASON==slctseason),]
fish4<- fish5 %>%
  filter(., YEAR <= dat.split.year) %>%
  mutate(., "pa" = ifelse(Number > 0, 1, 0))
# dat.test<- fish5 %>%
  # filter(., YEAR > dat.split.year) %>%
  # mutate(., "pa" = ifelse(Number > 0, 1, 0),
         # "pred.id" = paste("ID.", seq(from = 1, to = nrow(.), by = 1), sep = ""))

### Subset to stage for individual models ----> choose one only
fish4.ss=fish4 %>% filter(., Stg=='Adt'); stage='Adt'
fish4.ss=fish4 %>% filter(., Stg=='Juv'); stage='Juv'
fish4.ss=fish4 %>% filter(., Stg=='ich'); stage='Ich'

modname=paste(slctseason, ' ', stage, ' ',fname)
modname2=paste(stage, '.',fname,'.', slctseason, sep='')

### biomod2 ####
resp.var=fish4.ss$pa
# resp.var[resp.var > 0]=1 # presence absence
latlon=as.matrix((data.frame(fish4.ss$LON, fish4.ss$LAT)))
expl.var=fish4.ss %>% select(DEPTH, SURFTEMP, BOTTEMP, calfin_100m3:rug)
expl.var=data.frame(expl.var)
myBiomodData=BIOMOD_FormatingData(resp.var,
                                  expl.var,
                                  resp.xy = latlon,
                                  resp.name = modname2,
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

# 3-fold cross validation with 80/20 random split
myBiomodModelOut <- BIOMOD_Modeling( 
  myBiomodData, 
  models = c('CTA', 'RF', 'FDA', 'ANN', 'GAM', 'GBM'), 
  models.options = myBiomodOption, 
  NbRunEval=3, 
  DataSplit=80, 
  Prevalence=0.5, 
  VarImport=3,
  models.eval.meth = c('TSS','ROC'),
  SaveObj = TRUE,
  rescal.all.models = FALSE,
  do.full.models = FALSE,
  modeling.id = paste(modname,"FirstModeling",sep=""))

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
