load(file='/home/ryan/Git/NEhabitat/Final_merged_fish_corABN_Zoo_Ich.Rda')
library(dplyr)
library(lubridate)
library(ncdf4)
library(raster)
library(angstroms)
library(ncdump)       ## devtools::install_github("r-gris/ncdump")
library(tabularaster) ## devtools::install_github("r-gris/tabularaster")
library(rworldmap)

NAfill=data.frame(FData.abn$EST_TOWDATE)
NAfill$Y=FData.abn$YEAR
NAfill$M=month(FData.abn$EST_TOWDATE)
NAfill$D=day(FData.abn$EST_TOWDATE)
NAfill$dstr=format(NAfill$FData.abn.EST_TOWDATE, format='%Y-%m-%d')
NAfill$doy=yday(FData.abn$EST_TOWDATE)
NAfill$sday=dmrg4$sdoy
NAfill$lday=dmrg4$ldoy
NAfill$lat=FData.abn$LAT
NAfill$lon=FData.abn$LON
NAfill$lonbin=FData.abn$lonbin
NAfill$latbin=FData.abn$latbin
NAfill$zoona=is.na(FData.abn$calfin_100m3)
NAfill$ichna=is.na(FData.abn$gadmor_100m3)
save(NAfill, file='merged_coords_file.Rda')

## get list of roms daily files to open # RM_NWA-SZ.HCob05T_avg_1980-02-02.nc
b=NAfill$dstr[NAfill$Y>1979 & NAfill$Y<2015]
a=unique(b)
aa=paste('RM_NWA-SZ.HCob05T_avg_', a,'.nc', sep='')
romsfile='/home/ryan/roms/test2/1/RM_NWA-SZ.HCob05T_avg_1980-02-01.nc'
ncd <- NetCDF(romsfile)
vars=ncd$variable$name
vname <- "nlg"
dd <- romsdata(romsfile, varname = vname, slice = c(40, 1), transpose = T)
plot(dd)  ## this is pure 0-dimX, 0-dimY index space
longlat <- romscoords(romsfile, transpose = TRUE)
contour(longlat[[1]], add = TRUE, lty = 2)
contour(longlat[[2]], add = TRUE, lty = 2)
bound <- romsboundary(longlat)
# projection(bound) <- "+init=epsg:4326"
extent(bound)
# class      : Extent 
# xmin       : 274.8872 
# xmax       : 308.4482 
# ymin       : 32.22423 
# ymax       : 53.76019 
plot_cgrid(romsfile, include = c("u", "v", "rho","psi"), cell = TRUE)
coord <- romscoords(romsfile, spatial=c("lon_rho", "lat_rho")) #transpose=T)



for (i in 1:length(b)){
  romsfile=aa[i]
  flats=NAfill$lat[which(b %in% a[i])]
  flons=NAfill$lon[which(b %in% a[i])]
  fcoords=(data.frame(flons, flats))
  nc=nc_open(romsfile)
  # lon=ncvar_get(nc, 'lon')
  roms_nlg=ncvar_get(nc, 'nlg')
  m=raster(nc$var$nlg)
  extract(nc$var$nlg, SpatialPoints(fcoords), sp = T)
}



