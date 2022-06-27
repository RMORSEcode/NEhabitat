library(biomod2)
library(tidyr)
library(RColorBrewer)
library(dismo)
library(dplyr)
library(lubridate)
library(raster)
library(ggplot2)
load('/home/ryan/Documents/Git/NEhabitat/Final_merged_fish_corBIO_Zoo_Ich.Rda') # Biomass
load('/home/ryan/Documents/Git/NEhabitat/Final_merged_fish_corABN_Zoo_Ich.Rda') # Abundance

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

setwd('/home/ryan/Biomod models')
# get date to save in file names
xdt=today()
xdt=gsub('-','',xdt)

# Add variables from rasters to both biomass and numbers database -NOT NEEDED ANYMORE
### extract from rasters at fish points
# locs=FData.abn[,c('LON', 'LAT')]
# coordinates(locs) <- ~LON+LAT
# mypoints = SpatialPoints(locs,proj4string = CRS("+init=epsg:4326"))
# myproj = CRS(masked.raster)
# points.proj = spTransform(mypoints, myproj)
# load('/home/ryan/Documents/Git/NEhabitat/rast_phi_fraction.rdata') 
# ## convert phi (with negative vals) size to mm
# phi=masked.raster
# phi2mm=calc(phi, fun=function(x){2^-x})
# # plot(phi2mm)
# FData.abn$grnszmm = extract(phi2mm, mypoints)
# FData.bio$grnszmm = extract(phi2mm, mypoints)
# load('/home/ryan/Documents/Git/NEhabitat/rast_mud_fraction.rdata')
# FData.abn$mud_pct = extract(masked.raster, mypoints)
# FData.bio$mud_pct = extract(masked.raster, mypoints)
# load('/home/ryan/Documents/Git/NEhabitat/rast_sand_fraction.rdata')
# FData.abn$sand_pct = extract(masked.raster, mypoints)
# FData.bio$sand_pct = extract(masked.raster, mypoints)
# # load('/home/ryan/Git/NEhabitat/rast_gdepth.rdata')
# # FData.abn$gdepth = extract(masked.raster, mypoints) # odd scale neg vals start at 1000m
# load('/home/ryan/Documents/Git/NEhabitat/rast_rugosity.rdata')
# rug=masked.raster
# mmin=cellStats(rug, 'min')
# mmax=cellStats(rug, 'max')
# rugscl=calc(rug, fun=function(x){(x-mmin)/(mmax-mmin)}) # rescale from -4:2 -0:1
# # FData.abn$rug = extract(masked.raster, mypoints)
# # FData.bio$rug = extract(masked.raster, mypoints)
# # plot(masked.raster, xlim=c(-80,-60), ylim=c(36,48))
# # plot(rugscl)
# FData.abn$rug = extract(rugscl, mypoints)
# FData.bio$rug = extract(rugscl, mypoints)

# Select season for GAMs
slctseason="SPRING" #"FALL"
slctseason="FALL"


### Choose fish to run, adjust by searching on and changing " `73_ " to new ID ****
fname='Cod' #73
# fname='Haddock' #74
# fname='SilverHake' #72
# fname='Pollock' #75

###____________________________________________________________________
## use for zooplankton only
# fish=FData.abn %>% dplyr::select(YEAR, SEASON:`197_ich`, volume_100m3:chl12)
# fish$MONTH=month(FData.abn$EST_TOWDATE)
# fish=fish[complete.cases(fish),]
# ### change this ###
# fish$pa=ifelse(fish$ctyp_100m3>0, 1, 0)
# # fish2=fish %>% dplyr::select(-`calfin_100m3`) # change to PA and drop column
# # fish2=fish %>% dplyr::select(-`pseudo_100m3`) # change to PA and drop column
# fish2=fish %>% dplyr::select(-`ctyp_100m3`) # change to PA and drop column
# 
# fish2=fish2[which(fish2$SEASON==slctseason),]
# logd=fish2[,56:66]
# logd=log10(logd+1)
# fish2[,56:66]=logd
# trainPA=fish2[complete.cases(fish2),]
# trainPA.ss=trainPA[,c(3:9, 57:66, 71:82)]
# stage='zoo'
# fname='ctyp' #'pseudo' #'calfin'
###_______________________________________________________________________


### change fish to whatever species you are modeling!!!
fish=FData.abn %>% dplyr::select(YEAR, SEASON, LAT:BOTTEMP, `73_Adt`, `73_Juv`, `73_ich`, volume_100m3:chl12)
fish$MONTH=month(FData.abn$EST_TOWDATE)
fish=fish[complete.cases(fish),]
fish$`73_ich`=ceiling(fish$`73_ich`) # make integer from numbers per 100 m/3
fish2=fish[which(fish$SEASON==slctseason),] # subset to season

### REMOVE replicate values from zooplankton for the entire data set
tunq=fish2 %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()==1) %>% mutate(num=n())
tdup=fish2 %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()>1) %>% mutate(num=n())
tdupmn=tdup %>% mutate_if(is.numeric, mean) %>% mutate_if(is.character, funs(paste(unique(.), collapse = "_"))) %>% slice(1)
fish2=bind_rows(tunq, tdupmn) %>% ungroup()#reassign name

## subset to models used for HGAMs (after running all initially)
# fish2=fish2 %>% select(YEAR:BOTTEMP, DEPTH, chl10, chl2, grnszmm, ctyp_100m3, `73_Adt`, `73_Juv`, `73_ich`)
# for spring haddock limited models 20210225
# fish2=fish2 %>% select(YEAR:`73_ich`, ctyp_100m3,grnszmm,chl2, chl4,chl9, num)




# # ### now split data into training and testing set (75% train 25% test, randomly chosen)
# set.seed(101) # Set Seed so that same sample can be reproduced in future also
# # # Now Selecting 75% of data as sample from total 'n' rows of the data  
# sample <- sample.int(n = nrow(fish2), size = floor(.75*nrow(fish2)), replace = F)
# trainPA <- fish2[sample, ]
# testPA  <- fish2[-sample, ]

### pivot longer so stage is repeated
trainPA=fish2 %>% pivot_longer(c(`73_Adt`, `73_Juv`, `73_ich`), names_to = c("SVSPP", "Stg"), names_sep ="_", values_to = "Number")
trainPA$Stg=factor(trainPA$Stg, ordered=F)
# table(trainPA$MONTH)
## log transform plankton for training set
logd=trainPA[,8:19]
logd=log10(logd+1)
trainPA[,8:19]=logd
trainPA$pa=ifelse(trainPA$Number>0, 1, 0)
trainPA=trainPA[complete.cases(trainPA),]

### only used if running select vars for ensemble approach
# trainPA$ctyp_100m3=log10(trainPA$ctyp_100m3+1)

# Now set up traing and testing data; split data sets
# 1) Train model w/ data through 1999, predict 2000-2018
# 2) Train 1977-2011, predict 2012-2018
# 3) Train 1977-2012, predict 2013-2018
# 4) Train 1977-2014, predict 2015-2018
## Select Season and Fish species
## Training/testing year split
# dat.split.year<- 2018
#  

### Subset to stage for individual models ----choose one only
trainPA.ss=trainPA %>% filter(., Stg=='Adt'); stage='Adt'
trainPA.ss=trainPA %>% filter(., Stg=='Juv'); stage='Juv'
trainPA.ss=trainPA %>% filter(., Stg=='ich'); stage='Ich'
# trainPA.ss$DEPTH=as.numeric(trainPA.ss$DEPTH)

# trainPA.ss=trainPA.ss %>% select(., -SEASON)
# ### now split data into training and testing set (75% train 25% test, randomly chosen)
set.seed(101) # Set Seed so that same sample can be reproduced in future also
sample <- sample.int(n = nrow(trainPA.ss), size = floor(.75*nrow(trainPA.ss)), replace = F)
testPA.ss  <- trainPA.ss[-sample, ]
trainPA.ss <- trainPA.ss[sample, ]
explvardata='full'# 'limited' #'full' ## whether all or limited data used for explanatory vars
modname=paste(slctseason, ' ', stage, ' ',explvardata,' ', fname)
modname2=paste(stage, '_',fname,'_', explvardata,'_', slctseason, sep='')

### biomod2 #### 
resp.var=trainPA.ss$pa
# resp.var[resp.var 0]=1 # presence absence
latlon=as.matrix((data.frame(trainPA.ss$LON, trainPA.ss$LAT)))
# expl.var=trainPA.ss %>% select(DEPTH, SURFTEMP, BOTTEMP, volume_100m3:chl12)
# expl.var=trainPA.ss %>% select(-LAT, -LON,-YEAR, -SEASON)
# eval.expl.var.ss= data.frame(testPA.ss %>% select(-LAT, -LON, -YEAR, -SEASON))

### subset to only those found important from full suite (run that first)
# expl.var=expl.var %>% select(DEPTH:BOTSALIN, chl1:chl12)
# xdt=paste(xdt,'selectedvarsonly', sep='')
# eval.expl.var.ss= data.frame(testPA.ss %>% select(DEPTH:BOTSALIN, chl1:chl12)) #data.frame(testPA.ss %>% select(DEPTH, SURFTEMP, BOTTEMP, volume_100m3:chl12))

### limited models -- for spring haddock models
# expl.var=trainPA.ss %>% select(DEPTH, SURFTEMP, ctyp_100m3, grnszmm, chl10, chl2)
# expl.var=trainPA.ss %>% select(DEPTH:chl9) # 20210225 spring haddock adt limited
# expl.var$DEPTH=as.numeric(expl.var$DEPTH)
# eval.expl.var.ss=testPA.ss %>% select(DEPTH, SURFTEMP, ctyp_100m3, grnszmm, chl10, chl2)
# eval.expl.var.ss=testPA.ss %>% select(DEPTH:chl9)
# eval.expl.var.ss$DEPTH=as.numeric(eval.expl.var.ss$DEPTH)
# xdt=paste(xdt,'selectedvarsonly', sep='')

# USE FOR INITIAL RUN full suite of vars for initial models
expl.var=trainPA.ss %>% select(DEPTH:grnszmm,sand_pct,rug,chl2,chl10,chl4) #drop co-linear vars
eval.expl.var.ss=testPA.ss %>% select(DEPTH:grnszmm,sand_pct,rug,chl2,chl10,chl4)

# expl.var=data.frame(expl.var)
eval.resp.var.ss=testPA.ss$pa
# eval.expl.var.ss= eval.expl.var.ss  #data.frame(testPA.ss %>% select(DEPTH, SURFTEMP, BOTTEMP, volume_100m3:chl12))
latlon2=as.matrix((data.frame(testPA.ss$LON, testPA.ss$LAT)))

class(expl.var)="data.frame"
class(eval.expl.var.ss)="data.frame"
myBiomodData=BIOMOD_FormatingData(resp.var,
                                  expl.var,
                                  resp.xy = latlon,
                                  resp.name = modname2,
                                  eval.resp.var = eval.resp.var.ss,
                                  eval.expl.var = eval.expl.var.ss,
                                  eval.resp.xy = latlon2,
                                  PA.nb.rep = 0,
                                  PA.nb.absences = 1000,
                                  PA.strategy = 'random',
                                  PA.dist.min = 0,
                                  PA.dist.max = NULL,
                                  PA.sre.quant = 0.025,
                                  PA.table = NULL,
                                  na.rm = TRUE)

myBiomodData
# plot(myBiomodData)
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
  modeling.id = paste(modname2,"_",xdt,sep=""))

myBiomodModelOut 
myBiomodModelEval <- get_evaluations(myBiomodModelOut)
dimnames(myBiomodModelEval)
myBiomodModelEval["TSS","Testing.data","RF",,]
myBiomodModelEval["TSS","Testing.data",,,]

myBiomodModelEval["ROC","Testing.data",,,]
get_variables_importance(myBiomodModelOut)
test=get_variables_importance(myBiomodModelOut)
test2=apply(test, c(1,2), 'mean')
write.csv(format(test2, digits=3), file=paste(modname2, '_', xdt,'_variable_importance.csv', sep=''))
test3=matrix(data=NA, nrow=dim(test2)[1], ncol=dim(test2)[2])
for (i in 1:dim(test2)[2]){
  test3[,i]=rownames(test2)[rev(order(test2[,i]))]
}
colnames(test3)=colnames(test2)
# table(test3[1:7,])
rev(sort(table(test3[1:7,])))
rev(sort(table(test3[1:5,])))


### open saved list of variable importance and select important factors:
wd2='/home/ryan/Biomod models' 
wd='/home/ryan/Git/NEhabitat'
varlist=list.files(wd, pattern="variable_importance.csv")
varlist
## select certain runs more easily
# library(stringr)
# varlist2=str_detect(varlist, 'Cod')
# varlist3=varlist[varlist2]
i=13 # 2 10 14 | 1 13
varlist[i]
# test2=read.csv('/home/ryan/Biomod models/Ich_Haddock_SPRINGvariable_importance.csv', stringsAsFactors = F, row.names = 1)
test2=read.csv(paste(wd,'/',varlist[i],sep=''), stringsAsFactors = F, row.names = 1)
# standardize values
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
test2_norm <- as.data.frame(lapply(test2[1:6], min_max_norm))
test2sn=test2_norm
test2sn[test2_norm<.5]=NA
rowSums(test2sn,na.rm=T)
m=length(rownames(test2)[order(rowSums(test2sn,na.rm=T), decreasing=T)][ sort(rowSums(test2sn,na.rm=T), decreasing=T)>0])
m2=matrix(ncol=2, nrow=m)
m2=data.frame(m2)
m2[,1]=rownames(test2)[order(rowSums(test2sn,na.rm=T), decreasing=T)][ sort(rowSums(test2sn,na.rm=T), decreasing=T)>0]
m2[,2]=sort(rowSums(test2sn,na.rm=T), decreasing=T)[ sort(rowSums(test2sn,na.rm=T), decreasing=T)>0]
colnames(m2)=c('Var', 'St_Mean_6')
write.csv(m2, file=paste(wd2,'/','FINAL_SUMMARY_', varlist[i], sep=''))

test3=matrix(data=NA, nrow=dim(test2)[1], ncol=dim(test2)[2])
for (i in 1:dim(test2)[2]){
  test3[,i]=rownames(test2)[rev(order(test2[,i]))]
}
colnames(test3)=colnames(test2)
# table(test3[1:7,])
rev(sort(table(test3[1:7,])))
rev(sort(table(test3[1:5,])))




### Plot model performance:
### by models
gg1 <- models_scores_graph( myBiomodModelOut,
                            by = 'models',
                            metrics = c('ROC','TSS') )
gg1 + ggtitle(paste(modname))

### by cross validation run
gg2 <- models_scores_graph( myBiomodModelOut,
                            by = 'cv_run',
                            metrics = c('ROC','TSS') )
gg2 + ggtitle(paste(modname))

### plot response curves
myRF <- BIOMOD_LoadModels(myBiomodModelOut, models = 'RF')
myRespPlot2D=response.plot2(
  myRF,
  Data = get_formal_data(myBiomodModelOut, 'expl.var'),
  show.variables = get_formal_data(myBiomodModelOut,'expl.var.names'),
  do.bivariate = FALSE,
  fixed.var.metric = 'median',
  col = c("blue", "red"),
  legend = TRUE,
  data_species = get_formal_data(myBiomodModelOut, 'resp.var'),
  plot = TRUE)

myEnsembleOut=BIOMOD_EnsembleModeling( modeling.output=myBiomodModelOut,
                         chosen.models = 'all',
                         em.by = 'PA_dataset+repet',
                         eval.metric = 'all',
                         eval.metric.quality.threshold = NULL,
                         models.eval.meth = c('KAPPA','TSS','ROC'),
                         prob.mean = TRUE,
                         prob.cv = FALSE,
                         prob.ci = FALSE,
                         prob.ci.alpha = 0.05,
                         prob.median = FALSE,
                         committee.averaging = FALSE,
                         prob.mean.weight = FALSE,
                         prob.mean.weight.decay = 'proportional',
                         VarImport = 0)

get_evaluations(myEnsembleOut)

### load saved model example:
# rm(myBiomodModelOut)
# rm(myBiomodModelEval)
# load("~/Biomod models/Adt.Haddock.SPRING/Adt.Haddock.SPRING.Adt_Haddock_SPRING20200928.models.out")
# myBiomodModelOut=Adt.Haddock.SPRING.Adt_Haddock_SPRING20200928.models.out
# myBiomodModelOut=loadRData("~/Biomod models/SilverHake/adult.haddock.spr/adult.haddock.Adult HaddockFirstModeling.models.out")
myBiomodModelEval <- get_evaluations(myBiomodModelOut)

myBiomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExpl,
  proj.name = 'current',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')
