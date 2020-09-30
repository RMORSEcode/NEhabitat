### For spring haddock



# Now use Pederson et al (2019) as a model to format models:
#   
# 1. A single common smoother for all observations; We will refer to this as model *G*, as it only has a Global smoother.
# 2. A global smoother plus group-level smoothers that have the same wiggliness. We will refer to this as model *GS* (for Global smoother with individual effects that have a Shared penalty)
# 3. A global smoother plus group-level smoothers with differing wiggliness. We will refer to this as model *GI* (for Global smoother with individual effects that have Individual penalties)
# 4. Group-specific smoothers without a global smoother, but with all smoothers having the same wiggliness. We will refer to this as model *S*.
# 5. Group-specific smoothers with different wiggliness. We will refer to this as model *I*. 

# 1) Model G: single global smoother for all stages
# y ~ s(x, bs='tp') + s(fac, bs='re')
# Use binomial distribution for presence absence model, gaussian for biomass model, combine later with predict in 'HGAMpredict2raster.R'
fish_modG_pa = gam(pa ~ s(SURFTEMP, k=20, bs='tp') + s(DEPTH, k=10, bs="tp") + s(Stg, k=3, bs='re') + s(ctyp_100m3, k=20, bs="cr") + s(grnszmm, k=10, bs='tp') + s(chl2, k=20, bs='tp') + s(chl10, k=20, bs='tp'),  data=trainPA, method = "REML", family="binomial")
fish_modG_pb = gam(logbio ~ (SURFTEMP, k=20, bs='tp') + s(DEPTH, k=10, bs="tp") + s(Stg, k=3, bs='re') + s(ctyp_100m3, k=20, bs="cr") + s(grnszmm, k=10, bs='tp') + s(chl2, k=20, bs='tp') + s(chl10, k=20, bs='tp'), data=trainBIO, method = "REML", family="gaussian")
save(fish_modG_pa, file=paste('/home/ryan/Git/NEhabitat/rasters/',fishseas,'/',fishname,'/fish_modG_pa_',fishseas,'_',fishname,'.Rdata',sep=''))
save(fish_modG_pb, file=paste('/home/ryan/Git/NEhabitat/rasters/',fishseas,'/',fishname,'/fish_modG_pb_',fishseas,'_',fishname,'.Rdata',sep=''))   


# Add latitude and longitude as interaction terms to 'Fish_modG'
fish_modG4_pa = gam(pa ~ s(LON, LAT) + s(SURFTEMP, k=20, bs='tp') + s(DEPTH, k=10, bs="tp") + s(Stg, k=4, bs='re') + s(ctyp_100m3, k=20, bs="cr") + s(grnszmm, k=10, bs='tp') + , bs='tp'), data=trainPA, method = "REML", family="binomial")
fish_modG4_pb = gam(logbio ~ s(LON, LAT) +s(SURFTEMP, k=20, bs='tp') + s(DEPTH, k=10, bs="tp") + s(Stg, k=4, bs='re') + s(ctyp_100m3, k=20, bs="cr") + s(grnszmm, k=10, bs='tp') + s(chl2, k=20, bs='tp') + s(chl10, k=20, bs='tp'), data=trainBIO, method = "REML", family="gaussian")
save(fish_modG4_pa, file=paste('/home/ryan/Git/NEhabitat/rasters/',fishseas,'/',fishname,'/fish_modG4_pa_',fishseas,'_',fishname,'.Rdata',sep=''))
save(fish_modG4_pb, file=paste('/home/ryan/Git/NEhabitat/rasters/',fishseas,'/',fishname,'/fish_modG4_pb_',fishseas,'_',fishname,'.Rdata',sep=''))  


# TRY above model using soap film smoother on Lat an Lon
# m2 <- gam(-depth ~ s(os_x, os_y, bs = "so", xt = list(bnd = bound)),data = depth, method = "REML", knots = knots)
# Now use Lat and Lon for soapfilm smoother
# Create soap film boundary for fitting GAMs
library("rgdal")
# gdepth=raster('/home/ryan/Desktop/nes_bath_data.nc', band=1)
# gdepth[gdepth>0]=NA
# gdepth[gdepth< -500]=NA
outline <- read.csv("/home/ryan/Desktop/Final_NES_shelf_shape.csv")
outl.crd=sp::coordinates(cbind(outline$longitude, outline$latitude))
# change bound to use LAT LON, rerun
boundll <- list(list(x = outl.crd[,1], y = outl.crd[,2], f = rep(0, nrow(outl.crd))))
N <- 25
gx <- seq(min(outl.crd[,1]), max(outl.crd[,1]), len = N)
gy <- seq(min(outl.crd[,2]), max(outl.crd[,2]), len = N)
gp <- expand.grid(gx, gy)
names(gp) <- c("x","y")
knotsll <- gp[with(gp, inSide(boundll, x, y)), ]
names(knotsll) <- c("LON", "LAT")
names(boundll[[1]]) <- c("LON", "LAT", "f")
xym <- cbind(outline$longitude, outline$latitude)
p = sp::Polygon(xym)
ps = sp::Polygons(list(p),1)
sps = sp::SpatialPolygons(list(ps))
proj4string(sps) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
trainPApts=SpatialPoints(cbind(trainPA$LON, trainPA$LAT),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
trainPA.in=sp::over(trainPApts, sps)
trainPA.sub=trainPA
trainPA.sub=trainPA.sub[complete.cases(trainPA.in),]
# plot(boundll[[1]]$LON, boundll[[1]]$LAT, type='l', asp=1); points(knotsll, add=T, col='red')
trainBIOpts=SpatialPoints(cbind(trainBIO$LON, trainBIO$LAT),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
trainBIO.in=sp::over(trainBIOpts, sps)
trainBIO.sub=trainBIO
trainBIO.sub=trainBIO.sub[complete.cases(trainBIO.in),]


##### Chunk to remove select knots on boundary (specific to N specified in above chunk)
## remove knots on the boundary
## for N=25
knotsll=knotsll[-19,]
knotsll=knotsll[-22,]
knotsll=knotsll[-37,]
knotsll=knotsll[-77,]
knotsll=knotsll[-129,]

### Now run soap film smoother with LAT and LON values
fish_modG4bll_pa = gam(pa ~ s(LON, LAT, bs = "so", xt = list(bnd = boundll)) + s(SURFTEMP, k=20, bs='tp') + s(DEPTH, k=10, bs="tp") + s(Stg, k=4, bs='fs') + s(ctyp_100m3, k=20, bs="cr") + s(grnszmm, k=10, bs='tp') + s(chl2, k=20, bs='tp') + s(chl10, k=20, bs='tp'), data=trainPA.sub, method = "REML", knots=knotsll, family="binomial")
fish_modG4bll_pb = gam(logbio ~ s(LON, LAT, bs = "so", xt = list(bnd = boundll)) + s(SURFTEMP, k=20, bs='tp') + s(DEPTH, k=10, bs="tp") + s(Stg, k=4, bs='fs') + s(ctyp_100m3, k=20, bs="cr") + s(grnszmm, k=10, bs='tp') + s(chl2, k=20, bs='tp') + s(chl10, k=20, bs='tp'), data=trainBIO.sub, method = "REML", knots=knotsll, family="gaussian")
save(fish_modG4bll_pa, file=paste('/home/ryan/Git/NEhabitat/rasters/',fishseas,'/',fishname,'/fish_modG4bll_pa_',fishseas,'_',fishname,'.Rdata',sep=''))
save(fish_modG4bll_pb, file=paste('/home/ryan/Git/NEhabitat/rasters/',fishseas,'/',fishname,'/fish_modG4bll_pb_',fishseas,'_',fishname,'.Rdata',sep=''))  




# GS model: Y~ s(conc, bs='tp') +s(conc, plant, bs='fs')
# remove lat and lon
fish_modGSe_pa =  gam(pa ~ s(SURFTEMP, bs='tp') + s(SURFTEMP, Stg, k=10, bs='fs') + s(DEPTH, bs='tp') + s(DEPTH, Stg, k=10, bs="fs") + s(ctyp_100m3, bs='tp') + s(ctyp_100m3, Stg, k=10, bs="fs") + s(grnszmm, bs='tp') + s(grnszmm, Stg, k=10, bs='fs') + s(chl2, k=20, bs='tp') + s(chl2, Stg, k=20, bs='fs') +s(chl10, k=20, bs='tp') + s(chl10, Stg, k=20, bs='fs'), data=trainPA, method = "REML", family="binomial")
fish_modGSe_pb =  gam(logbio ~ s(SURFTEMP, bs='tp') + s(SURFTEMP, Stg, k=10, bs='fs') + s(DEPTH, bs='tp') + s(DEPTH, Stg, k=10, bs="fs") + s(ctyp_100m3, bs='tp') + s(ctyp_100m3, Stg, k=10, bs="fs") + s(grnszmm, bs='tp') + s(grnszmm, Stg, k=10, bs='fs') + s(chl2, k=20, bs='tp') + s(chl2, Stg, k=20, bs='fs') +s(chl10, k=20, bs='tp') + s(chl10, Stg, k=20, bs='fs'), data=trainBIO, method = "REML", family="gaussian")
save(fish_modGSe_pa, file=paste('/home/ryan/Git/NEhabitat/rasters/',fishseas,'/',fishname,'/fish_modGSe3_pa_',fishseas,'_',fishname,'.Rdata',sep=''))
save(fish_modGSe_pb, file=paste('/home/ryan/Git/NEhabitat/rasters/',fishseas,'/',fishname,'/fish_modGSe3_pb_',fishseas,'_',fishname,'.Rdata',sep=''))


# GI model
# y ~ s(x, bs='tp') + s(x, by=fac, m=1, bs='tp') + s(fac, bs='re')
fish_modGI_pa =  gam(pa ~ s(SURFTEMP, k=10, m=2, bs='tp') + s(SURFTEMP, by=Stg, k=10, m=1, bs='tp') + s(DEPTH, k=10, m=2, bs="tp") + s(DEPTH, by=Stg, k=10, m=1, bs="tp") + s(ctyp_100m3, k=10, m=2, bs="tp") + s(ctyp_100m3, m=1, k=10, by=Stg, bs="tp")+ s(grnszmm, m=2, k=10, bs='tp') + s(grnszmm, by=Stg, k=10, m=1, bs='tp') + s(chl2, m=2, k=20, bs='tp') +s(chl2, by=Stg, k=20, m=1, bs='tp') + s(chl10, m=2, k=20, bs='tp') +s(chl10, by=Stg, k=20, m=1, bs='tp') +s(Stg, k=10, bs='re'), data=trainPA, method = "REML", family="binomial")
fish_modGI_pb =  gam(logbio ~ s(SURFTEMP, k=10, m=2, bs='tp') + s(SURFTEMP, by=Stg, k=10, m=1, bs='tp') + s(DEPTH, k=10, m=2, bs="tp") + s(DEPTH, by=Stg, k=10, m=1, bs="tp") + s(ctyp_100m3, k=10, m=2, bs="tp") + s(ctyp_100m3, m=1, k=10, by=Stg, bs="tp")+ s(grnszmm, m=2, k=10, bs='tp') + s(grnszmm, by=Stg, k=10, m=1, bs='tp') + s(chl2, m=2, k=20, bs='tp') +s(chl2, by=Stg, k=20, m=1, bs='tp') + s(chl10, m=2, k=20, bs='tp') +s(chl10, by=Stg, k=20, m=1, bs='tp') +s(Stg, k=10, bs='re'), data=trainBIO, method = "REML", family="gaussian")
save(fish_modGI_pa, file=paste('/home/ryan/Git/NEhabitat/rasters/',fishseas,'/',fishname,'/fish_modGI_pa_',fishseas,'_',fishname,'.Rdata',sep=''))
save(fish_modGI_pb, file=paste('/home/ryan/Git/NEhabitat/rasters/',fishseas,'/',fishname,'/fish_modGI_pb_',fishseas,'_',fishname,'.Rdata',sep=''))


# S model
# y∼s(x, fac, bs="fs") or y∼t2(x1, x2, fac, bs=c("tp", "tp", "re")
fish_modS_pa =  gam(pa ~ s(SURFTEMP, Stg, k=10, bs='fs') + s(DEPTH, Stg, k=10, bs="fs") + s(ctyp_100m3, Stg, k=10, bs="fs") + s(grnszmm, Stg, k=10, bs='fs') + s(chl2, Stg, k=20, bs='fs') +  s(chl10, Stg, k=20, bs='fs'), data=trainPA, method = "REML", family="binomial")
fish_modS_pb =  gam(logbio ~ s(SURFTEMP, Stg, k=10, bs='fs') + s(DEPTH, Stg, k=10, bs="fs") + s(ctyp_100m3, Stg, k=10, bs="fs") + s(grnszmm, Stg, k=10, bs='fs') + s(chl2, Stg, k=20, bs='fs') +  s(chl10, Stg, k=20, bs='fs'), data=trainBIO, method = "REML", family="gaussian")
save(fish_modS_pa, file=paste('/home/ryan/Git/NEhabitat/rasters/',fishseas,'/',fishname,'/fish_modS_pa_',fishseas,'_',fishname,'.Rdata',sep=''))
save(fish_modS_pb, file=paste('/home/ryan/Git/NEhabitat/rasters/',fishseas,'/',fishname,'/fish_modS_pb_',fishseas,'_',fishname,'.Rdata',sep='')) 

                 
# I model
# y∼fac+s(x, by=fac) or y∼fac+te(x1,x2, by=fac)
fish_modI_pa =  gam(pa ~ Stg + s(SURFTEMP, by=Stg, k=10, bs='tp') + s(DEPTH, by=Stg, k=10, bs="tp") + s(ctyp_100m3, by=Stg, k=10, bs="tp") + s(grnszmm, by=Stg, k=10, bs='tp') + s(chl2, by=Stg, k=20, bs='tp') + s(chl10, by=Stg, k=20, bs='tp'), data=trainPA, method = "REML", family="binomial")
fish_modI_pb =  gam(logbio ~ Stg + s(SURFTEMP, by=Stg, k=10, bs='tp') + s(DEPTH, by=Stg, k=10, bs="tp") + s(ctyp_100m3, by=Stg, k=10, bs="tp") + s(grnszmm, by=Stg, k=10, bs='tp') + s(chl2, by=Stg, k=20, bs='tp') + s(chl10, by=Stg, k=20, bs='tp'), data=trainBIO, method = "REML", family="gaussian")
save(fish_modI_pa, file=paste('/home/ryan/Git/NEhabitat/rasters/',fishseas,'/',fishname,'/fish_modI_pa_',fishseas,'_',fishname,'.Rdata',sep=''))
save(fish_modI_pb, file=paste('/home/ryan/Git/NEhabitat/rasters/',fishseas,'/',fishname,'/fish_modI_pb_',fishseas,'_',fishname,'.Rdata',sep=''))

                             