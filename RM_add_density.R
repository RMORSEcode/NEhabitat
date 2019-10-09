### RM 20190829 addding water mass determination to habitat modeling efforts


library(oce)
library(raster)

### load bottom salinity
load("G:/1_habitat_analysis_2017/botsal/spring_spdf/rasters/RAST_NESREG_2012.04.03.BS.TEMP.YEAR.000066596.RData")
BS=masked.raster

### load bottom temperature
load("G:/1_habitat_analysis_2017/bottemp/spring_spdf/rasters/RAST_NESREG_2012.04.03.BT.TEMP.YEAR.000066596.RData")
BT=masked.raster

### load depth raster (note different crs and extent than NESREG data rasters)
load("G:/1_habitat_analysis_2017/static_vars/rast_gdepth.rdata")
z=masked.raster
z2=z
# z[z>0]=
z[z>=0]=0.001 # get rid of all positive altitudes - NA does not work with oce
z2[z2>=0]=NA
# p=swPressure(z2, latitude = 42)

### calc pressure from depth
zvec=as.data.frame.vector(z)
p=swPressure(z, latitude = 42)

mp=matrix(p, nrow=721, ncol=961, byrow=T)
rp=raster(mp)
tt=crs(z)
crs(rp)=tt

tt=extent(z)
extent(rp)=tt

s=as.vector(BS)
t=as.vector(BT)


tt=crs(rp)
crs(BS)=tt

### now resample fine grid to NESREG
rp2=resample(rp, BS)
mskrp2=mask(rp2, BS)
p2=as.vector(mskrp2)

### compute density using salinity, temp, and pressure
rho=swRho(s, t, p2)
rhomat=matrix(rho, nrow=90, ncol=105, byrow=T)
Rho=raster(rhomat)
plot(Rho)