### For fall haddock haddock
### based on biomod runs and variable importance for each of 3 stages run seperately
## add date to model names

library(lubridate)
xdt=today()
xdt=gsub('-','',xdt)
wdhome='/home/ryan/Git/NEhabitat/rasters/'
wdhome='C:/Users/ryan.morse/Documents/GitHub/NEhabitat/rasters/'

### open saved list of variable importance and select important factors:
wd='/home/ryan/Git/NEhabitat'
wd='/home/ryan/Documents/Git/NEhabitat' #second machine
impvarlist=list.files(wd, pattern="variable_importance.csv")
impvarlist
xxi=2
impvarlist[xxi]
# test2=read.csv('/home/ryan/Biomod models/Ich_Haddock_SPRINGvariable_importance.csv', stringsAsFactors = F, row.names = 1)
test2=read.csv(paste(wd,'/',impvarlist[xxi],sep=''), stringsAsFactors = F, row.names = 1)
test3=matrix(data=NA, nrow=dim(test2)[1], ncol=dim(test2)[2])
for (i in 1:dim(test2)[2]){
  test3[,i]=rownames(test2)[rev(order(test2[,i]))]
}
colnames(test3)=colnames(test2)
# table(test3[1:7,])
rev(sort(table(test3[1:7,])))
rev(sort(table(test3[1:5,])))
barplot(rev(sort(table(test3[1:5,]))))
### variables of importance from biomod2 runs:
## Juv [1:5]
# chl1    DEPTH     chl2  grnszmm  BOTTEMP SURFTEMP     chl9    chl12    chl10     chl3 
## Adt [1:7]
# SURFTEMP     chl9     chl1 sand_pct  grnszmm    DEPTH     chl7     chl2  BOTTEMP     chl8     chl5     chl4 
# 5        5        5        4        4        4        3        3        3        1        1        1 
# chl3    chl12    chl10 
# 1        1        1 
# Adt [1:5]
# SURFTEMP sand_pct     chl9     chl1    DEPTH     chl7  grnszmm     chl2     chl8     chl4  BOTTEMP 
# 5        4        4        4        3        3        2        2        1        1        1 
## Selected vars:
# chl1    DEPTH     chl2  grnszmm  BOTTEMP SURFTEMP     chl9
## to change:
# ctyp_100m3 -> chl9
# chl10 > chl1
# add BOTTEMP
## done

impvarlist2=impvarlist[grepl("*Haddock_full_FALL*", impvarlist)] # select all stages from a season
tt=matrix(NA, ncol=12, nrow=2)
for(ii in 1:2){
  xxi=ii
test2=read.csv(paste(wd,'/',impvarlist[xxi],sep=''), stringsAsFactors = F, row.names = 1)
test3=matrix(data=NA, nrow=dim(test2)[1], ncol=dim(test2)[2])
for (i in 1:dim(test2)[2]){
  test3[,i]=rownames(test2)[rev(order(test2[,i]))]
}
colnames(test3)=colnames(test2)
# table(test3[1:7,])
im1=rev(sort(table(test3[1:7,])))
im2=rev(sort(table(test3[1:7,])))
im3=rev(sort(table(test3[1:7,])))

test1=data.frame(im1) %>% filter(Freq>=3)
test2=data.frame(im2) %>% filter(Freq>=3)
test3=data.frame(im3) %>% filter(Freq>=3)

# Reduce(intersect, list(test1[,1], test3[,1], test2[,1])) # too limiting...
## this isnt working right, need to sort first I think
# t4=test1 %>% full_join(test2, keep=T, by = c("Var1", "Freq"))
# t5=t4 %>% full_join(test3, keep=T, by = c("Var1.x"="Var1", "Freq.x"="Freq"))

x1=sort(as.character(test1$Var1))
x2=sort(as.character(test2$Var1))
x3=sort(as.character(test3$Var1))
xf=paste(c(x1,x2,x3))
# unique(xf)
tt[ii,1:length(unique(xf))]=unique(xf)
fnm=strsplit(impvarlist2, split="_")[[ii]][1]
rownames(tt)[ii]=fnm
}
write.csv(tt, file=paste(wd,'/haddock_fall_variable_importance_final_selection.csv', sep=''))

#  NEW MODEL 20210312:
# [1] "BOTTEMP"    "chl2"       "chl4"       "DEPTH"      "grnszmm"    "sand_pct"   "SURFTEMP"   "chl10" "ctyp_100m3" ""   
# current model: SURFTEMP, BOTTEMP,DEPTH,Stg,chl9,grnszmm,chl2,chl1
#  To add:"chl4" "sand_pct" "chl10" "ctyp_100m3"
#  to remove: chl9, chl1 done



### UPDATE DATA to remove ichtyhyoplankton for fall haddock (only 6 positive samples 1977-2017)
trainPA=trainPA %>% filter(Stg!='ich') %>% droplevels(exclude="ich")
# levels(trainPA$Stg)
testPA=testPA %>% filter(Stg!='ich') %>% droplevels(exclude="ich")
trainBIO=trainBIO %>% filter(Stg!='ich') %>% droplevels(exclude="ich")
testBIO=testBIO %>% filter(Stg!='ich') %>% droplevels(exclude="ich")

### set directory of location to save runs
wd='/home/ryan/Git/NEhabitat/rasters/'
wd='/home/ryan/Documents/Git/NEhabitat/rasters/'


# + s(ctyp_100m3, k=20, bs="cs") + s(sand_pct, k=15, bs='ts') 
# 1) Model G: single global smoother for all stages
# y ~ s(x, bs='ts') + s(fac, bs='re')
# Use binomial distribution for presence absence model, gaussian for biomass model, combine later with predict in 'HGAMpredict2raster.R'
fish_modG_pa = gam(pa ~ s(SURFTEMP, k=30, bs='ts') + s(ctyp_100m3, k=20, bs="cs") + s(sand_pct, k=15, bs='ts') + s(BOTTEMP, k=30, bs='ts') + s(DEPTH, k=20, bs="ts") + s(Stg, k=3, bs='re') + s(chl4, k=30, bs="cs") + s(grnszmm, k=15, bs='ts') + s(chl2, k=30, bs='cs') + s(chl10, k=30, bs='cs'),  data=trainPA, method = "REML", family="binomial", select=T)
fish_modG_pb = gam(logbio ~ s(SURFTEMP, k=30, bs='ts') + s(ctyp_100m3, k=20, bs="cs") + s(sand_pct, k=15, bs='ts') + s(BOTTEMP, k=30, bs='ts') + s(DEPTH, k=20, bs="ts") + s(Stg, k=3, bs='re') + s(chl4, k=30, bs="cs") + s(grnszmm, k=15, bs='ts') + s(chl2, k=30, bs='cs') + s(chl10, k=30, bs='cs'),  data=trainBIO, method = "REML", family="gaussian", select=T)
save(fish_modG_pa, file=paste(wdhome,fishseas,'/',fishname,'/fish_modG_pa_',fishseas,'_',fishname,'_',xdt,'.Rdata',sep=''))
save(fish_modG_pb, file=paste(wdhome,fishseas,'/',fishname,'/fish_modG_pb_',fishseas,'_',fishname,'_',xdt,'.Rdata',sep=''))   


# Add latitude and longitude as interaction terms to 'Fish_modG'
fish_modG4_pa = gam(pa ~ s(LON, LAT) + s(SURFTEMP, k=30, bs='ts') + s(ctyp_100m3, k=20, bs="cs") + s(sand_pct, k=15, bs='ts') + s(ctyp_100m3, k=20, bs="cs") + s(sand_pct, k=15, bs='ts') + s(BOTTEMP, k=30, bs='ts') + s(DEPTH, k=20, bs="ts") + s(Stg, k=3, bs='re') + s(chl4, k=30, bs="cs") + s(grnszmm, k=15, bs='ts') + s(chl2, k=30, bs='ts') + s(chl10, k=30, bs='ts'), data=trainPA, method = "REML", family="binomial", select=T)
fish_modG4_pb = gam(logbio ~ s(LON, LAT) +s(SURFTEMP, k=30, bs='ts') + s(ctyp_100m3, k=20, bs="cs") + s(sand_pct, k=15, bs='ts') + s(ctyp_100m3, k=20, bs="cs") + s(sand_pct, k=15, bs='ts') + s(BOTTEMP, k=30, bs='ts') + s(DEPTH, k=20, bs="ts") + s(Stg, k=3, bs='re') + s(chl4, k=30, bs="cs") + s(grnszmm, k=15, bs='ts') + s(chl2, k=30, bs='ts') + s(chl10, k=30, bs='ts'), data=trainBIO, method = "REML", family="gaussian", select=T)
save(fish_modG4_pa, file=paste(wdhome,fishseas,'/',fishname,'/fish_modG4_pa_',fishseas,'_',fishname,'_',xdt,'.Rdata',sep=''))
save(fish_modG4_pb, file=paste(wdhome,fishseas,'/',fishname,'/fish_modG4_pb_',fishseas,'_',fishname,'_',xdt,'.Rdata',sep=''))  


# # TRY above model using soap film smoother on Lat an Lon
# # m2 <- gam(-depth ~ s(os_x, os_y, bs = "so", xt = list(bnd = bound)),data = depth, method = "REML", knots = knots)
# # Now use Lat and Lon for soapfilm smoother
# # Create soap film boundary for fitting GAMs
# library("rgdal")
# # gdepth=raster('/home/ryan/Desktop/nes_bath_data.nc', band=1)
# # gdepth[gdepth>0]=NA
# # gdepth[gdepth< -500]=NA
# outline <- read.csv("/home/ryan/Desktop/Final_NES_shelf_shape.csv")
# outl.crd=sp::coordinates(cbind(outline$longitude, outline$latitude))
# # change bound to use LAT LON, rerun
# boundll <- list(list(x = outl.crd[,1], y = outl.crd[,2], f = rep(0, nrow(outl.crd))))
# N <- 25
# gx <- seq(min(outl.crd[,1]), max(outl.crd[,1]), len = N)
# gy <- seq(min(outl.crd[,2]), max(outl.crd[,2]), len = N)
# gp <- expand.grid(gx, gy)
# names(gp) <- c("x","y")
# knotsll <- gp[with(gp, inSide(boundll, x, y)), ]
# names(knotsll) <- c("LON", "LAT")
# names(boundll[[1]]) <- c("LON", "LAT", "f")
# xym <- cbind(outline$longitude, outline$latitude)
# p = sp::Polygon(xym)
# ps = sp::Polygons(list(p),1)
# sps = sp::SpatialPolygons(list(ps))
# proj4string(sps) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
# trainPApts=SpatialPoints(cbind(trainPA$LON, trainPA$LAT),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# trainPA.in=sp::over(trainPApts, sps)
# trainPA.sub=trainPA
# trainPA.sub=trainPA.sub[complete.cases(trainPA.in),]
# # plot(boundll[[1]]$LON, boundll[[1]]$LAT, type='l', asp=1); points(knotsll, add=T, col='red')
# trainBIOpts=SpatialPoints(cbind(trainBIO$LON, trainBIO$LAT),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# trainBIO.in=sp::over(trainBIOpts, sps)
# trainBIO.sub=trainBIO
# trainBIO.sub=trainBIO.sub[complete.cases(trainBIO.in),]
# ##### Chunk to remove select knots on boundary (specific to N specified in above chunk)
# ## remove knots on the boundary
# ## for N=25
# knotsll=knotsll[-19,]
# knotsll=knotsll[-22,]
# knotsll=knotsll[-37,]
# knotsll=knotsll[-77,]
# knotsll=knotsll[-129,]
# ### Now run soap film smoother with LAT and LON values
# fish_modG4bll_pa = gam(pa ~ s(LON, LAT, bs = "so", xt = list(bnd = boundll)) + s(SURFTEMP, k=30, bs='ts') + s(BOTTEMP, k=30, bs='ts') + s(DEPTH, k=20, bs="ts") + s(Stg, k=3, bs='fs') + s(chl4, k=30, bs="cs") + s(grnszmm, k=15, bs='ts') + s(chl2, k=30, bs='ts') + s(chl10, k=30, bs='ts'), data=trainPA.sub, method = "REML", knots=knotsll, family="binomial", select=T)
# fish_modG4bll_pb = gam(logbio ~ s(LON, LAT, bs = "so", xt = list(bnd = boundll)) + s(SURFTEMP, k=30, bs='ts') + s(BOTTEMP, k=30, bs='ts') + s(DEPTH, k=20, bs="ts") + s(Stg, k=3, bs='fs') + s(chl4, k=30, bs="cs") + s(grnszmm, k=15, bs='ts') + s(chl2, k=30, bs='ts') + s(chl10, k=30, bs='ts'), data=trainBIO.sub, method = "REML", knots=knotsll, family="gaussian", select=T)
# save(fish_modG4bll_pa, file=paste('/home/ryan/Git/NEhabitat/rasters/',fishseas,'/',fishname,'/fish_modG4bll_pa_',fishseas,'_',fishname,'_',xdt,'.Rdata',sep=''))
# save(fish_modG4bll_pb, file=paste('/home/ryan/Git/NEhabitat/rasters/',fishseas,'/',fishname,'/fish_modG4bll_pb_',fishseas,'_',fishname,'_',xdt,'.Rdata',sep=''))  

# + s(ctyp_100m3, k=20, bs="cs") +s(ctyp_100m3, Stg, k=20, bs="fs") + s(sand_pct, k=15, bs='ts') + s(sand_pct, Stg, k=15, bs='fs') 
# GS model: Y~ s(conc, bs='ts') +s(conc, plant, bs='fs')
# remove lat and lon

fish_modGS_pa =  gam(pa ~ s(SURFTEMP, k=30, bs='ts') + s(SURFTEMP, Stg, k=30, bs='fs')  + s(ctyp_100m3, k=20, bs="cs") + 
                        s(ctyp_100m3, Stg, k=20, bs="fs") + s(sand_pct, k=15, bs='ts') + s(sand_pct, Stg, k=15, bs='fs') + 
                        s(BOTTEMP, k=30, bs='ts') + s(BOTTEMP, Stg, k=30, bs='fs') + s(DEPTH, k=20, bs='ts') + s(DEPTH, Stg, k=20, bs="fs") +
                        s(chl4, k=30, bs='ts') + s(chl4, Stg, k=30, bs="fs") + s(grnszmm, k=15, bs='ts') + s(grnszmm, Stg, k=15, bs='fs') + 
                        s(chl2, k=30, bs='ts') + s(chl2, Stg, k=30, bs='fs') +s(chl10, k=30, bs='ts') + s(chl10, Stg, k=30, bs='fs'), 
                      data=trainPA, method = "REML", family="binomial", select=T)
fish_modGS_pb =  gam(logbio ~ s(SURFTEMP, k=30, bs='ts') + s(SURFTEMP, Stg, k=30, bs='fs') + s(ctyp_100m3, k=20, bs="cs") +s(ctyp_100m3, Stg, k=20, bs="fs") + 
                        s(sand_pct, k=15, bs='ts') + s(sand_pct, Stg, k=15, bs='fs')+ s(BOTTEMP, k=30, bs='ts') + s(BOTTEMP, Stg, k=30, bs='fs') + 
                        s(DEPTH, k=20, bs='ts') + s(DEPTH, Stg, k=20, bs="fs") + s(chl4, k=30, bs='ts') + s(chl4, Stg, k=30, bs="fs") + s(grnszmm, k=15, bs='ts') +
                        s(grnszmm, Stg, k=15, bs='fs') + s(chl2, k=30, bs='ts') + s(chl2, Stg, k=30, bs='fs') +s(chl10, k=30, bs='ts') + s(chl10, Stg, k=30, bs='fs'), 
                      data=trainBIO, method = "REML", family="gaussian", select=T)
save(fish_modGS_pa, file=paste(wdhome,fishseas,'/',fishname,'/fish_modGSe_pa_',fishseas,'_',fishname,'_',xdt,'.Rdata',sep=''))
save(fish_modGS_pb, file=paste(wdhome,fishseas,'/',fishname,'/fish_modGSe_pb_',fishseas,'_',fishname,'_',xdt,'.Rdata',sep=''))

# + s(ctyp_100m3, k=20, bs="cs") +s(ctyp_100m3, by=Stg, k=20, m=1, bs="ts") + s(sand_pct, k=15, bs='ts') + s(sand_pct, by=Stg, k=15, m=1, bs='ts')
# GI model
# y ~ s(x, bs='ts') + s(x, by=fac, m=1, bs='ts') + s(fac, bs='re')
fish_modGI_pa =  gam(pa ~ s(SURFTEMP, k=30, m=2, bs='ts') + s(SURFTEMP, by=Stg, k=30, m=1, bs='ts') + s(ctyp_100m3, k=20, bs="cs") + 
                       s(ctyp_100m3, by=Stg, k=20, m=1, bs="ts") + s(sand_pct, k=15, bs='ts') + s(sand_pct, by=Stg, k=15, m=1, bs='ts')+ 
                       s(BOTTEMP, k=30, m=2, bs='ts') + s(BOTTEMP, by=Stg, k=30, m=1, bs='ts') + s(DEPTH, k=20, m=2, bs="ts") + 
                       s(DEPTH, by=Stg, k=20, m=1, bs="ts") + s(chl4, k=30, m=2, bs="ts") + s(chl4, m=1, k=30, by=Stg, bs="ts")+ 
                       s(grnszmm, m=2, k=15, bs='ts') + s(grnszmm, by=Stg, k=15, m=1, bs='ts') + s(chl2, m=2, k=30, bs='ts') + 
                       s(chl2, by=Stg, k=30, m=1, bs='ts') + s(chl10, m=2, k=30, bs='ts') +s(chl10, by=Stg, k=30, m=1, bs='ts') + 
                       s(Stg, k=30, bs='re'), data=trainPA, method = "REML", family="binomial", select=T)
fish_modGI_pb =  gam(logbio ~ s(SURFTEMP, k=30, m=2, bs='ts') + s(SURFTEMP, by=Stg, k=30, m=1, bs='ts') + s(ctyp_100m3, k=20, bs="cs") + 
                       s(ctyp_100m3, by=Stg, k=20, m=1, bs="ts") + s(sand_pct, k=15, bs='ts') + s(sand_pct, by=Stg, k=15, m=1, bs='ts')+ 
                       s(BOTTEMP, m=2, k=30, bs='ts') + s(BOTTEMP, by=Stg, m=1, k=30, bs='ts') + s(DEPTH, k=20, m=2, bs="ts") + 
                       s(DEPTH, by=Stg, k=20, m=1, bs="ts") + s(chl4, k=30, m=2, bs="ts") + s(chl4, m=1, k=30, by=Stg, bs="ts")+ 
                       s(grnszmm, m=2, k=15, bs='ts') + s(grnszmm, by=Stg, k=15, m=1, bs='ts') + s(chl2, m=2, k=30, bs='ts') + 
                       s(chl2, by=Stg, k=30, m=1, bs='ts') + s(chl10, m=2, k=30, bs='ts') +s(chl10, by=Stg, k=30, m=1, bs='ts') + 
                       s(Stg, k=30, bs='re'), data=trainBIO, method = "REML", family="gaussian", select=T)
save(fish_modGI_pa, file=paste(wdhome,fishseas,'/',fishname,'/fish_modGI_pa_',fishseas,'_',fishname,'_',xdt,'.Rdata',sep=''))
save(fish_modGI_pb, file=paste(wdhome,fishseas,'/',fishname,'/fish_modGI_pb_',fishseas,'_',fishname,'_',xdt,'.Rdata',sep=''))

# + s(ctyp_100m3, Stg, k=20, bs="fs") + s(sand_pct, Stg, k=15, bs='fs')
# S model
# y∼s(x, fac, bs="fs") or y∼t2(x1, x2, fac, bs=c("ts", "ts", "re")
fish_modS_pa =  gam(pa ~ s(SURFTEMP, Stg, k=30, bs='fs') + s(ctyp_100m3, Stg, k=20, bs="fs") + s(sand_pct, Stg, k=15, bs='fs') + 
                      s(BOTTEMP, Stg, k=30, bs='fs') + s(DEPTH, Stg, k=20, bs="fs") + s(chl4, Stg, k=30, bs="fs") + 
                      s(grnszmm, Stg, k=15, bs='fs') + s(chl2, Stg, k=30, bs='fs') +  s(chl10, Stg, k=30, bs='fs'), 
                    data=trainPA, method = "REML", family="binomial", select=T)
fish_modS_pb =  gam(logbio ~ s(SURFTEMP, Stg, k=30, bs='fs') + s(ctyp_100m3, Stg, k=20, bs="fs") + s(sand_pct, Stg, k=15, bs='fs') + 
                      s(BOTTEMP, Stg, k=30, bs='fs') + s(DEPTH, Stg, k=20, bs="fs") + s(chl4, Stg, k=30, bs="fs") + 
                      s(grnszmm, Stg, k=15, bs='fs') + s(chl2, Stg, k=30, bs='fs') +  s(chl10, Stg, k=30, bs='fs'), 
                    data=trainBIO, method = "REML", family="gaussian", select=T)
save(fish_modS_pa, file=paste(wdhome,fishseas,'/',fishname,'/fish_modS_pa_',fishseas,'_',fishname,'_',xdt,'.Rdata',sep=''))
save(fish_modS_pb, file=paste(wdhome,fishseas,'/',fishname,'/fish_modS_pb_',fishseas,'_',fishname,'_',xdt,'.Rdata',sep='')) 

# + s(ctyp_100m3, by=Stg, k=20, bs="ts") + s(sand_pct, by=Stg, k=15, bs='ts')
# I model
# y∼fac+s(x, by=fac) or y∼fac+te(x1,x2, by=fac)
fish_modI_pa =  gam(pa ~ Stg + s(SURFTEMP, by=Stg, k=30, bs='ts') + s(ctyp_100m3, by=Stg, k=20, bs="ts") + s(sand_pct, by=Stg, k=15, bs='ts') +
                      s(BOTTEMP, by=Stg, k=30, bs='ts') + s(DEPTH, by=Stg, k=20, bs="ts") + s(chl4, by=Stg, k=20, bs="ts") + 
                      s(grnszmm, by=Stg, k=15, bs='ts') + s(chl2, by=Stg, k=20, bs='ts') + s(chl10, by=Stg, k=20, bs='ts'), 
                    data=trainPA, method = "REML", family="binomial", select=T)
fish_modI_pb =  gam(logbio ~ Stg + s(SURFTEMP, by=Stg, k=30, bs='ts') + s(ctyp_100m3, by=Stg, k=20, bs="ts") + s(sand_pct, by=Stg, k=15, bs='ts') + 
                      s(BOTTEMP, by=Stg, k=30, bs='ts') + s(DEPTH, by=Stg, k=20, bs="ts") + s(chl4, by=Stg, k=20, bs="ts") + 
                      s(grnszmm, by=Stg, k=15, bs='ts') + s(chl2, by=Stg, k=20, bs='ts') + s(chl10, by=Stg, k=20, bs='ts'), 
                    data=trainBIO, method = "REML", family="gaussian", select=T)
save(fish_modI_pa, file=paste(wdhome,fishseas,'/',fishname,'/fish_modI_pa_',fishseas,'_',fishname,'_',xdt,'.Rdata',sep=''))
save(fish_modI_pb, file=paste(wdhome,fishseas,'/',fishname,'/fish_modI_pb_',fishseas,'_',fishname,'_',xdt,'.Rdata',sep=''))

