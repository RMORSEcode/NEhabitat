## plot variable importance summary
varimp=readxl::read_xlsx('/home/ryan/Git/NEhabitat/Biomod_variable_importance_summary.xlsx')
varlist=readxl::read_xlsx('/home/ryan/Git/NEhabitat/Biomod_variable_importance_list.xlsx', col_names =T)
# table(unique(varlist))
sort(table(varlist))


modchoice=i
modlistpa[modchoice]
usemodel=loadRData(paste(path1,modlistpa[modchoice], sep=''))

par(oma=c(3,2,2,0))
par(mar=c(2,4,2,1) + 0.1)

### plotting and fitting for haddock
test=fish2
test=test %>% pivot_longer(c(`74_Adt`, `74_Juv`, `74_ich`), names_to = c("SVSPP", "Stg"), 
                           names_sep ="_", values_to = "Number")
test$Stg=factor(test$Stg, ordered=F)
logd=test[,8:19]
logd=log10(logd+1)
test[,8:19]=logd
test$pa=ifelse(test$Number > 0, 1, 0)

test=test %>% filter(Stg!='ich') %>% droplevels(exclude="ich")
test1=predict.gam(usemodel, test, type='response')
test$pred=test1

plotROC(test$pa,test$pred, colorize = TRUE, main=paste('All Stages '))
# AUC by year
sbtest=test %>% dplyr::filter(YEAR==2013)
plotROC(sbtest$pa,sbtest$pred, colorize = TRUE, main=paste('All Stages '))


saveAUC <- function(truth, predicted, ...){
  pred <- prediction(as.vector(abs(predicted)), as.vector(truth))    
  roc <- performance(pred,"tpr","fpr")
  auc <- performance(pred,"auc")
  # plot(roc, ...)
  # abline(a=0, b= 1)
  return(auc@y.values)
}

senSpec=function(truth, predicted, ...){
  pred <- prediction(as.vector(abs(predicted)), as.vector(truth))    
  perf1 <- performance(pred, "sens", "spec")
  plot(perf1)
  return(perf1)
}

### SAVE out model performance data to CSV file - deviance explained, AIC
modeval=data.frame(matrix(nrow=length(yrlist), ncol=3, data=NA))
for (i in 1:length(yrlist)){
  j=yrlist[i]
  sbtest=test %>% dplyr::filter(YEAR==j, Stg=='Adt')
  sbtest2=test %>% dplyr::filter(YEAR==j, Stg=='Juv')
  sbtest3=test %>% dplyr::filter(YEAR==j, Stg=='ich')
  # modchoice=i
  # usemodel=loadRData(paste(path1,modlistpa[modchoice], sep=''))
  # usemodelbio=loadRData(paste(path1,modlistpb[modchoice], sep=''))
  # modeval[i,1]=modlistpa[modchoice]
  if (sum(sbtest$pa)<3){
    next
  }
  modeval[i,1]=saveAUC(sbtest$pa,sbtest$pred)[[1]]
  modeval[i,2]=saveAUC(sbtest2$pa,sbtest2$pred)[[1]]
  modeval[i,3]=saveAUC(sbtest2$pa,sbtest3$pred)[[1]]
  
  # modeval[i,2]=summary(usemodel)$dev.expl
  # modeval[i,3]=summary(usemodelbio)$dev.expl
  # modeval[i,4]=usemodel$aic
  # modeval[i,5]=usemodelbio$aic
  # modeval[i,6]=sum(usemodel$edf)
  # modeval[i,7]=sum(usemodelbio$edf)
  # modeval[i,8]=df.residual(usemodel)
  # modeval[i,9]=df.residual(usemodelbio)
  # modeval[i,10]=((df.residual(usemodel)+sum(usemodel$edf))/3)-sum(usemodel$edf)
  # modeval[i,11]=((df.residual(usemodelbio)+sum(usemodelbio$edf))/3)-sum(usemodelbio$edf)
}
colnames(modeval)=c('adtAUC', 'juvAUC', 'ichAUC')
modeval$year=yrlist
barplot(modeval$adtAUC, names.arg=yrlist, las=1, main='Adt')
barplot(modeval$juvAUC, names.arg=yrlist, las=1, main="Juv")
plot(modeval$adtAUC~ yrlist, las=1, main="Adt", type='b', pch=19, ylim=c(0.8, 1))
plot(modeval$juvAUC~ yrlist, las=1, main="Juv", type='b', pch=19, ylim=c(0.7, 1))


senSpec(sbtest$pa,sbtest$pred)


# colnames(modeval)=c('model', 'PA.dev.exp','BIO.dev.exp','PA.aic','BIO.aic','PA.edf','BIO.edf','PA.res.df','BIO.res.df', 'PA.corr.res.df','BIO.corr.res.df')
write.csv(format(modeval, digits=2), file=paste(wd2,'Model_AUC_by_year_', SEASON, '_', fishnm, '_', '.csv', sep=""), row.names = F)


SEASON="Fall"
fishnm="Haddock"

p1=paste('/home/ryan/Git/NEhabitat/rasters/',SEASON,'/', fishnm, '/', sep='')
p2='fish_modGSe_Fall_haddock' #'fish_modG4_spr_Haddock/'
p3=paste('/PA_only_stacked_', SEASON, '_', fishnm, '_', sep='') #'PA_only_stacked_Spr_Haddock_'
p4=paste('/stacked_', SEASON, '_', fishnm, '_', sep='') #'stacked_Spr_Haddock_'
ichpa=loadRData(paste(p1,p2,p3,'ich.RData', sep=''))
juvpa=loadRData(paste(p1,p2,p3,'Juv.RData', sep=''))
adtpa=loadRData(paste(p1,p2,p3,'Adt.RData', sep=''))
ich=loadRData(paste(p1,p2,p4,'ich.RData', sep=''))
juv=loadRData(paste(p1,p2,p4,'Juv.RData', sep=''))
adt=loadRData(paste(p1,p2,p4,'Adt.RData', sep=''))

ichhab_gbk=raster::extract(ichpa, gbk, fun=mean, na.rm=T)
plot(ichhab_gbk[1,]~yrlist, type='l', las=1)
ichhab_gom=raster::extract(ichpa, gom, fun=mean, na.rm=T)
plot(ichhab_gom[1,]~yrlist, type='l', las=1)

adthab_gbk=raster::extract(adtpa, gbk, fun=mean, na.rm=T)
plot(adthab_gbk[1,]~yrlist, type='l', las=1)
adthab_gom=raster::extract(adtpa, gom, fun=mean, na.rm=T)
plot(adthab_gom[1,]~yrlist, type='l', las=1)

juvhab_gbk=raster::extract(juvpa, gbk, fun=mean, na.rm=T)
plot(juvhab_gbk[1,]~yrlist, type='l', las=1)
juvhab_gom=raster::extract(juvpa, gom, fun=mean, na.rm=T)
plot(juvhab_gom[1,]~yrlist, type='l', las=1)
## using trawl suvey stata
ichhab_tsGBK=raster::extract(ichpa, tsGBK, fun=mean, na.rm=T)
plot(ichhab_tsGBK[1,]~yrlist, type='l', las=1)
ichhab_tsGOM=raster::extract(ichpa, tsGOM, fun=mean, na.rm=T)
plot(ichhab_tsGOM[1,]~yrlist, type='l', las=1)

adthab_tsGBK=raster::extract(adtpa, tsGBK, fun=mean, na.rm=T)
plot(adthab_tsGBK[1,]~yrlist, type='l', las=1)
adthab_tsGOM=raster::extract(adtpa, tsGOM, fun=mean, na.rm=T)
plot(adthab_tsGOM[1,]~yrlist, type='l', las=1)

juvhab_tsGBK=raster::extract(juvpa, tsGBK, fun=mean, na.rm=T)
plot(juvhab_tsGBK[1,]~yrlist, type='l', las=1)
juvhab_tsGOM=raster::extract(juvpa, tsGOM, fun=mean, na.rm=T)
plot(juvhab_tsGOM[1,]~yrlist, type='l', las=1)

adthab_tsGBK=raster::extract(adt, tsGBK, fun=sum, na.rm=T)
plot(adthab_tsGBK[1,]~yrlist, type='l', las=1)
adthab_tsGOM=raster::extract(adt, tsGOM, fun=mean, na.rm=T)
plot(adthab_tsGOM[1,]~yrlist, type='l', las=1)

juvhab_tsGBK=raster::extract(juv, tsGBK, fun=mean, na.rm=T)
plot(juvhab_tsGBK[1,]~yrlist, type='l', las=1)
juvhab_tsGOM=raster::extract(juv, tsGOM, fun=mean, na.rm=T)
plot(juvhab_tsGOM[1,]~yrlist, type='l', las=1)

### using kevins RF results
adthab_tsGBK=raster::extract(kfpas, tsGBK, fun=mean, na.rm=T)
plot(adthab_tsGBK[1,]~yrlist, type='l', las=1)
adthab_tsGBK=raster::extract(kfpaf, tsGBK, fun=mean, na.rm=T)
plot(adthab_tsGBK[1,]~yrlist, type='l', las=1)
adthab_tsGOM=raster::extract(kfpas, tsGOM, fun=mean, na.rm=T)
plot(adthab_tsGOM[1,]~yrlist, type='l', las=1)
adthab_tsGOM=raster::extract(kfpaf, tsGOM, fun=mean, na.rm=T)
plot(adthab_tsGOM[1,]~yrlist, type='l', las=1)

adthab_tsGBK=raster::extract(kfbms, tsGBK, fun=mean, na.rm=T)
plot(adthab_tsGBK[1,]~yrlist, type='l', las=1)


t=fish2 %>% filter(SEASON=='Spring') %>% group_by(YEAR)


# from Liz Brooks survey based RSSB data for GB haddock
lbhad=readxl::read_xlsx('/home/ryan/Desktop/Haddock/SSB.yr_and_R.xlsx', skip=4)
plot(lbhad$ssb.yr.Spring~lbhad$years, type='l', main='SSB')
plot(lbhad$Age0.Fall ~lbhad$years, type='l', main='Age 0 Fall')
lines(lbhad$`Age1.Spring(lagged)` ~lbhad$years, type='l', col='red')
lbhad$r0ssb=lbhad$Age0.Fall/lbhad$ssb.yr.Spring
lbhad$r1ssb=lbhad$`Age1.Spring(lagged)`/lbhad$ssb.yr.Spring

lbhad2=lbhad[complete.cases(lbhad),]
lbhad2$anom0=(log10(lbhad2$r0ssb)-mean(log10(lbhad2$r0ssb)))/sd(log10(lbhad2$r0ssb))
lbhad2$anom1=(log10(lbhad2$r1ssb)-mean(log10(lbhad2$r1ssb)))/sd(log10(lbhad2$r1ssb))
lbhad3=left_join(lbhad, lbhad2, by=c('years', 'ssb.yr.Spring', 'Age0.Fall', 'Age1.Spring(lagged)', 'r0ssb', 'r1ssb'))

# log10(hdts$recr)-mean(log10(hdts$recr))/sd(log10(hdts$recr))
plot(log10(lbhad2$r0ssb)~lbhad2$years, type='l')
plot(lbhad3$anom0~lbhad3$years, type='l', main='age0 rssb'); abline(h=0, lty=2)
plot(lbhad3$anom1~lbhad3$years, type='l', main='age1 rssb'); abline(h=0, lty=2)

plot(lbhad3$anom0~lbhad3$years, type='l', col='black', main='rssb'); abline(h=0, lty=2)
lines(lbhad3$anom1~lbhad3$years, col='red')

### read in modeled recruits spawners from GOM cod stock assessment update report 2019
# M=0.2 for fixed mortality, or ramp from 0.2 to 0.4 for mramp
# recruits are age1 numbers on Jan 1 (not lagged yet)
codrssb=readxl::read_xlsx('/home/ryan/Desktop/GOM_cod_rssb_2019.xlsx', skip=2)
# codrssb$rssb=codrssb$`Age1 (000s)`/codrssb$`SSB_mfix (mt)`
# codrssb$anom=(log10(codrssb$rssb)-mean(log10(codrssb$rssb)))/sd(log10(codrssb$rssb))
codrssb$anom2=(log10(codrssb$rssb_mramp_lag2)-mean(log10(codrssb$rssb_mramp_lag2), na.rm=T))/sd(log10(codrssb$rssb_mramp_lag2),na.rm=T)
codrssb$anom=(log10(codrssb$rssb_mramp_lag1)-mean(log10(codrssb$rssb_mramp_lag1), na.rm=T))/sd(log10(codrssb$rssb_mramp_lag1),na.rm=T)
plot(codrssb$anom ~codrssb$Yr, type='l', main='age1 rssb'); abline(h=0, lty=2)


# trying dynamic time warping GOM RSSB vs Juv habitat
library(dtw)
juvts=(juvhab_gbk[1,2:36]-mean(juvhab_gbk[1,2:36]))/sd(juvhab_gbk[1,2:36])
adtts=(adthab_gbk[1,4:38]-mean(adthab_gbk[1,2:36]))/sd(adthab_gbk[1,2:36])
ichts=(ichhab_gbk[1,4:38]-mean(ichhab_gbk[1,2:36]))/sd(ichhab_gbk[1,2:36])
rssbts=ssbr2[10:44,19]

# cross correlation coefficient plots
ccf(ichts, rssbts)
ccf(juvts, rssbts)
ccf(adtts, rssbts)
# DTW plots
alignment<-dtw(ichts,rssbts,keep=TRUE)
plot(alignment,type="threeway", main='Ich')
dtwPlotDensity(alignment, main='Ich')
alignment<-dtw(juvts,rssbts,keep=TRUE)
plot(alignment,type="threeway", main='Juv')
dtwPlotDensity(alignment, main='Juv')
alignment<-dtw(adtts,rssbts,keep=TRUE)
plot(alignment,type="threeway", main='Adt')
dtwPlotDensity(alignment, main='Adt')

plot(dtw(ichts,rssbts,keep=TRUE,
         step=rabinerJuangStepPattern(6,"c")),
     type="twoway",offset=-2, main='Ich vs RSSB');
plot(dtw(juvts,rssbts,keep=TRUE,
         step=rabinerJuangStepPattern(6,"c")),
     type="twoway",offset=-2, main='Juv vs RSSB');
plot(dtw(adtts,rssbts,keep=TRUE,
         step=rabinerJuangStepPattern(6,"c")),
     type="twoway",offset=-2, main='Adt vs RSSB');


plotRasterMeanNegBinom=function(rastck){
  time <- 1:nlayers(rastck) 
  newrast.m=calc(rastck, fun=mean, na.rm=T)
  mn=cellStats(newrast.m, min)
  mx=cellStats(newrast.m, max)
  high=max(abs(mn), mx)
  br <- seq(-high, high, by = high/15) 
  cl=colorRampPalette(brewer.pal(9,"RdYlBu"))(length(br))
  rng=range(newrast.m[],na.rm=T)
  arg=list(at=rng, labels=round(rng,3))
  plot(newrast.m, col=cl, breaks=br,axis.args=arg,las=1, main=paste('mean')) # Yearly slope
  maps::map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="black", add=T)
}
  
  fun=function(x) { if (is.na(x[1])){ NA } else { m = lm(x ~ time); summary(m)$coefficients[2] }}
  newrast=calc(rastck, fun)
  mn=cellStats(newrast, min)
  mx=cellStats(newrast, max)
  high=max(abs(mn), mx)
  br <- seq(mn, mx, by = high/15) 
  cl <- colorspace::diverge_hcl(length(br) - 1, power = 1) 
  rng=range(newrast[],na.rm=T)
  arg=list(at=rng, labels=round(rng,3))
  plot(newrast, col=cl, breaks=br,axis.args=arg,las=1, main=paste(length(time),'yrs','\nYearly Slope')) # Yearly slope
  maps::map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="black", add=T)
  # plot slope times length of time
  newrast.t=newrast*length(time)
  mn=cellStats(newrast.t, min) #min(newrast.t@data@values, na.rm = T)
  mx=cellStats(newrast.t, max) #max(newrast.t@data@values, na.rm = T)
  high=max(abs(mn), mx)
  br <- seq(-high, high, by = high/15) 
  cl <- colorspace::diverge_hcl(length(br) - 1, power = 1) 
  rng=range(newrast.t[],na.rm=T)
  arg=list(at=rng, labels=round(rng,3))
  plot(newrast.t, col=cl, breaks=br,axis.args=arg,las=1, main=paste(length(time),'yr change')) # Time series change
  maps::map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="black", add=T)
  ## Calc Significance
  fun=function(x) { if (is.na(x[1])){ NA } else { m = lm(x ~ time); summary(m)$coefficients[8] }}
  p <- calc(rastck, fun=fun)
  m = c(0, 0.05, 1, 0.05, 1, 0)
  rclmat = matrix(m, ncol=3, byrow=TRUE)
  p.mask = reclassify(p, rclmat)
  fun=function(x) { x[x<1] <- NA; return(x)}
  p.mask.NA = calc(p.mask, fun)
  trend.sig = mask(newrast, p.mask.NA)
  mn=cellStats(newrast, min)
  mx=cellStats(newrast, max)
  high=max(abs(mn), mx)
  br <- seq(-high, high, by = high/15) 
  cl <- colorspace::diverge_hcl(length(br) - 1, power = 1) 
  rng=range(newrast[],na.rm=T)
  arg=list(at=rng, labels=round(rng,3))
  plot(newrast, main=paste(length(time),'yrs','\nYearly Slope'),col=cl, breaks=br,axis.args=arg,las=1) # Yearly slope
  test=rasterToPoints(trend.sig)
  points(test, pch='+', col=addTrans('black',50), cex=0.5)
  maps::map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="black", add=T)
  ## Plot Variance
  fun=function(x) { if (is.na(x[1])){ NA } else { m = lm(x ~ time); (summary(m)$sigma)}} ; TITL='Std deviation of slope'
  newrast.v=calc(rastck, fun)
  mn=cellStats(newrast.v, min)
  mx=cellStats(newrast.v, max)
  high=max(abs(mn), mx)
  br <- seq(0, high, by = high/15) 
  cl=colorRampPalette(brewer.pal(9,"Reds"))(length(br))
  rng=c(0, mx) #range(newrast.v[],na.rm=T)
  arg=list(at=rng, labels=round(rng,3))
  plot(newrast.v, main=paste(length(time),'yrs','\n', TITL),col=cl, breaks=br,axis.args=arg,las=1) # 
  maps::map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="black", add=T)
}


t=fish2 %>% select(YEAR, MONTH, LAT, LON, `74_ich`)
t2=
  plot(t$Value~t$year, type='b', main=paste(selepu), las=1, xlab='', ylab='anom')

t2=aggregate(t, by=list(t$YEAR), FUN=mean, na.rm=T)


xx=35
ich1=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/stacked_Spr_Haddock_ich.RData')
plot(ich1[[xx]])
yrlist[xx]

map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70")
map.axes(las=1)
test=fish2 %>% filter(YEAR==yrlist[xx])
points(test$LON[which(test$`73_ich`==0)], test$LAT[which(test$`73_ich`==0)], col='grey70', pch=19)
points(test$LON[which(test$`73_ich`>0)], test$LAT[which(test$`73_ich`>0)], col='black', pch=19)
plot(nesbath,deep=-200, shallow=-200, step=1,add=T,lwd=1,col=addTrans('black',150),lty=1)

#__________________ extract and plot
# haddock=1
# cod=2
# s=spring
# f=fall
#spring GSe haddock
ich1s=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/PA_only_stacked_Spr_Haddock_ich.RData')
juv1s=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/PA_only_stacked_Spr_Haddock_Juv.RData')
adt1s=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGSe_Spr_20210421_haddock/PA_only_stacked_Spr_Haddock_Adt.RData')
#spring GSe cod
ich2s=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod_20210518/PA_only_stacked_Spr_Cod_ich.RData')
juv2s=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod_20210518/PA_only_stacked_Spr_Cod_Juv.RData')
adt2s=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGSe_Spr_cod_20210518/PA_only_stacked_Spr_Cod_Adt.RData')
#fall GSe haddock
juv1f=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock_20210421/PA_only_stacked_Fall_Haddock_Juv.RData')
adt1f=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGSe_Fall_haddock_20210421/PA_only_stacked_Fall_Haddock_Adt.RData')
#fall GSe cod updated from above (new data from Sean survdat)
juv2f=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod_20210518/PA_only_stacked_Fall_Cod_Juv.RData')
adt2f=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGSe_Fall_cod_20210518/PA_only_stacked_Fall_Cod_Adt.RData')

# # Testing GI models GS models used for everything before Sep 2021! NOTE ALSO that ONLY GS models rerun after data issue in April/May 2021
# this means that all model comparison output is not correct (models use different data) : checked habitat area plots for cod in spring and fall
# and they were very similar to the updated model GSe, so no cause for concern as they are very similar
# #spring GI haddock
# ich1s=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGI_Spr_haddock/PA_only_stacked_Spr_Haddock_ich.RData')
# juv1s=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGI_Spr_haddock/PA_only_stacked_Spr_Haddock_Juv.RData')
# adt1s=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Haddock/fish_modGI_Spr_haddock/PA_only_stacked_Spr_Haddock_Adt.RData')
# #spring GI cod
# ich2s=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGI_Spr_cod/PA_only_stacked_Spr_Cod_ich.RData')
# juv2s=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGI_Spr_cod/PA_only_stacked_Spr_Cod_Juv.RData')
# adt2s=loadRData('/home/ryan/Git/NEhabitat/rasters/Spr/Cod/fish_modGI_Spr_cod/PA_only_stacked_Spr_Cod_Adt.RData')
# #fall GI haddock
# juv1f=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGI_Fall_haddock/PA_only_stacked_Fall_Haddock_Juv.RData')
# adt1f=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Haddock/fish_modGI_Fall_haddock/PA_only_stacked_Fall_Haddock_Adt.RData')
# #fall GI cod updated from above (new data from Sean survdat)
# juv2f=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGI_Fall_cod/PA_only_stacked_Fall_Cod_Juv.RData')
# adt2f=loadRData('/home/ryan/Git/NEhabitat/rasters/Fall/Cod/fish_modGI_Fall_cod/PA_only_stacked_Fall_Cod_Adt.RData')


selstg=ich1s; sellab='Ich'; selfis='haddock'
selstg=juv1s; sellab='Juv'; selfis='haddock'
selstg=adt1s; sellab='Adt'; selfis='haddock'
selstg=ich2s; sellab='Ich'; selfis='cod'
selstg=juv2s; sellab='Juv'; selfis='cod'
selstg=adt2s; sellab='Adt'; selfis='cod'
xhab_tsGBK=raster::extract(selstg, tsGBK, fun=mean, na.rm=T)
plot(colMeans(xhab_tsGBK)~yrlist, type='l', las=1, ylab=paste(sellab,' GOM habitat'), xlab='', main=selfis)
abline(v=c(2003,2010, 2013, 2016), lty=2)
xhab_gb=raster::extract(selstg, gbk, fun=mean, na.rm=T)
plot(xhab_gb[1,]~yrlist, type='l', las=1, ylab=paste(sellab,' GB habitat', xlab=''))
abline(v=c(2003,2010, 2013, 2016), lty=2)

# Haddock
xhabich1s_tsGBK=raster::extract(ich1s, tsGBK, fun=mean, na.rm=T)
xhabjuv1s_tsGBK=raster::extract(juv1s, tsGBK, fun=mean, na.rm=T)
xhabadt1s_tsGBK=raster::extract(adt1s, tsGBK, fun=mean, na.rm=T)
xhabjuv1f_tsGBK=raster::extract(juv1f, tsGBK, fun=mean, na.rm=T)
xhabadt1f_tsGBK=raster::extract(adt1f, tsGBK, fun=mean, na.rm=T)
xhabich1s_gbk=raster::extract(ich1s, gbk, fun=mean, na.rm=T)
xhabjuv1s_gbk=raster::extract(juv1s, gbk, fun=mean, na.rm=T)
xhabadt1s_gbk=raster::extract(adt1s, gbk, fun=mean, na.rm=T)
xhabjuv1f_gbk=raster::extract(juv1f, gbk, fun=mean, na.rm=T)
xhabadt1f_gbk=raster::extract(adt1f, gbk, fun=mean, na.rm=T)
# Cod
xhabich2s_tsGBK=raster::extract(ich2s, tsGBK, fun=mean, na.rm=T)
xhabjuv2s_tsGBK=raster::extract(juv2s, tsGBK, fun=mean, na.rm=T)
xhabadt2s_tsGBK=raster::extract(adt2s, tsGBK, fun=mean, na.rm=T)
xhabjuv2f_tsGBK=raster::extract(juv2f, tsGBK, fun=mean, na.rm=T)
xhabadt2f_tsGBK=raster::extract(adt2f, tsGBK, fun=mean, na.rm=T)
xhabich2s_gbk=raster::extract(ich2s, gbk, fun=mean, na.rm=T)
xhabjuv2s_gbk=raster::extract(juv2s, gbk, fun=mean, na.rm=T)
xhabadt2s_gbk=raster::extract(adt2s, gbk, fun=mean, na.rm=T)
xhabjuv2f_gbk=raster::extract(juv2f, gbk, fun=mean, na.rm=T)
xhabadt2f_gbk=raster::extract(adt2f, gbk, fun=mean, na.rm=T)

### GOM TS (GOM) and EPU (gom)
# Haddock
xhabich1s_tsGOM=raster::extract(ich1s, tsGOM, fun=mean, na.rm=T)
xhabjuv1s_tsGOM=raster::extract(juv1s, tsGOM, fun=mean, na.rm=T)
xhabadt1s_tsGOM=raster::extract(adt1s, tsGOM, fun=mean, na.rm=T)
xhabjuv1f_tsGOM=raster::extract(juv1f, tsGOM, fun=mean, na.rm=T)
xhabadt1f_tsGOM=raster::extract(adt1f, tsGOM, fun=mean, na.rm=T)
xhabich1s_gom=raster::extract(ich1s, gom, fun=mean, na.rm=T)
xhabjuv1s_gom=raster::extract(juv1s, gom, fun=mean, na.rm=T)
xhabadt1s_gom=raster::extract(adt1s, gom, fun=mean, na.rm=T)
xhabjuv1f_gom=raster::extract(juv1f, gom, fun=mean, na.rm=T)
xhabadt1f_gom=raster::extract(adt1f, gom, fun=mean, na.rm=T)
# Cod
xhabich2s_tsGOM=raster::extract(ich2s, tsGOM, fun=mean, na.rm=T)
xhabjuv2s_tsGOM=raster::extract(juv2s, tsGOM, fun=mean, na.rm=T)
xhabadt2s_tsGOM=raster::extract(adt2s, tsGOM, fun=mean, na.rm=T)
xhabjuv2f_tsGOM=raster::extract(juv2f, tsGOM, fun=mean, na.rm=T)
xhabadt2f_tsGOM=raster::extract(adt2f, tsGOM, fun=mean, na.rm=T)
xhabich2s_gom=raster::extract(ich2s, gom, fun=mean, na.rm=T)
xhabjuv2s_gom=raster::extract(juv2s, gom, fun=mean, na.rm=T)
xhabadt2s_gom=raster::extract(adt2s, gom, fun=mean, na.rm=T)
xhabjuv2f_gom=raster::extract(juv2f, gom, fun=mean, na.rm=T)
xhabadt2f_gom=raster::extract(adt2f, gom, fun=mean, na.rm=T)


# all stages, cod and haddock spring TS GBK strata
plot(NULL, ylim=c(0.15,0.6), xlim=c(1977,2019), ylab="GBK Habitat", xlab="", las=1)
lines(colMeans(xhabich1s_tsGBK)~yrlist, lty=1, col='green')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(colMeans(xhabich2s_tsGBK)~yrlist, lty=3, col='green')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(colMeans(xhabjuv1s_tsGBK)~yrlist, lty=1, col='blue')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(colMeans(xhabjuv2s_tsGBK)~yrlist, lty=3, col='blue')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(colMeans(xhabadt1s_tsGBK)~yrlist, lty=1, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(colMeans(xhabadt2s_tsGBK)~yrlist, lty=3, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
abline(v=c(2003,2010, 2013, 2016), lty=2)
# all stages, cod and haddock fall
plot(NULL, ylim=c(0,0.5), xlim=c(1977,2019), ylab="GBK Habitat", xlab="", las=1)
lines(colMeans(xhabjuv1f_tsGBK)~yrlist, lty=1, col='blue')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(colMeans(xhabjuv2f_tsGBK)~yrlist, lty=3, col='blue')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(colMeans(xhabadt1f_tsGBK)~yrlist, lty=1, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(colMeans(xhabadt2f_tsGBK)~yrlist, lty=3, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')

# all stages GBK TS, haddock spring TS GBK strata
plot(NULL, ylim=c(0.15,0.5), xlim=c(1977,2019), ylab="GBK spr haddock habitat", xlab="", las=1)
lines(colMeans(xhabich1s_tsGBK)~yrlist, lty=1, col='green')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(colMeans(xhabjuv1s_tsGBK)~yrlist, lty=1, col='blue')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(colMeans(xhabadt1s_tsGBK)~yrlist, lty=1, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
abline(v=c(2003,2010, 2013, 2016), lty=3)
legend('topleft', legend = c('adt', 'juv', 'ich'),lty=1, col=c('red', 'blue', 'green'), bty='n', horiz = T)
# all stages, haddock fall
plot(NULL, ylim=c(0.3,0.5), xlim=c(1977,2019), ylab="GBK fall haddock habitat", xlab="", las=1)
lines(colMeans(xhabjuv1f_tsGBK)~yrlist, lty=1, col='blue')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(colMeans(xhabadt1f_tsGBK)~yrlist, lty=1, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
legend('topleft', legend = c('adt', 'juv'),lty=1, col=c('red', 'blue'), bty='n', horiz = T)
abline(v=c(2003,2010, 2013, 2016), lty=3)
# all stages, cod spring TS GBK strata
plot(NULL, ylim=c(0.15,0.6), xlim=c(1977,2019), ylab="GBK spr cod habitat", xlab="", las=1)
lines(colMeans(xhabich2s_tsGBK)~yrlist, lty=1, col='green')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(colMeans(xhabjuv2s_tsGBK)~yrlist, lty=1, col='blue')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(colMeans(xhabadt2s_tsGBK)~yrlist, lty=1, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
# abline(v=c(2003,2010, 2013, 2016), lty=2)
legend('topleft', legend = c('adt', 'juv', 'ich'),lty=1, col=c('red', 'blue', 'green'), bty='n', horiz = T)
# all stages, cod fall
plot(NULL, ylim=c(0,0.5), xlim=c(1977,2019), ylab="GBK fall cod habitat", xlab="", las=1)
lines(colMeans(xhabjuv2f_tsGBK)~yrlist, lty=1, col='blue')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(colMeans(xhabadt2f_tsGBK)~yrlist, lty=1, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
legend('topleft', legend = c('adt', 'juv'),lty=1, col=c('red', 'blue'), bty='n', horiz = T)



# all stages, cod and haddock spring GB EPU
plot(NULL, ylim=c(0.2,0.7), xlim=c(1977,2019), ylab="GB Habitat", xlab="", las=1)
lines(xhabich1s_gbk[1,]~yrlist, lty=1, col='green')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(xhabich2s_gbk[1,]~yrlist, lty=3, col='green')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(xhabjuv1s_gbk[1,]~yrlist, lty=1, col='blue')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(xhabjuv2s_gbk[1,]~yrlist, lty=3, col='blue')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(xhabadt1s_gbk[1,]~yrlist, lty=1, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(xhabadt2s_gbk[1,]~yrlist, lty=3, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
abline(v=c(2003,2010, 2013, 2016), lty=2)
# all stages, cod and haddock fall
plot(NULL, ylim=c(0,0.4), xlim=c(1977,2019), ylab="GB Habitat", xlab="", las=1)
lines(xhabjuv1f_gbk[1,]~yrlist, lty=1, col='blue')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(xhabjuv2f_gbk[1,]~yrlist, lty=3, col='blue')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(xhabadt1f_gbk[1,]~yrlist, lty=1, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(xhabadt2f_gbk[1,]~yrlist, lty=3, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')

# all stages haddock spring GB EPU
plot(NULL, ylim=c(0.2,0.5), xlim=c(1977,2019), ylab="GB spr haddock habitat", xlab="", las=1)
lines(xhabich1s_gbk[1,]~yrlist, lty=1, col='green')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(xhabjuv1s_gbk[1,]~yrlist, lty=1, col='blue')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(xhabadt1s_gbk[1,]~yrlist, lty=1, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
abline(v=c(2003,2010, 2013, 2016), lty=3)
legend('topleft', legend = c('adt', 'juv', 'ich'),lty=1, col=c('red', 'blue', 'green'), bty='n', horiz = T)
# all stages haddock fall
plot(NULL, ylim=c(0,0.4), xlim=c(1977,2019), ylab="GB fall haddock habitat", xlab="", las=1)
lines(xhabjuv1f_gbk[1,]~yrlist, lty=1, col='blue')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(xhabadt1f_gbk[1,]~yrlist, lty=1, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
legend('topleft', legend = c('adt', 'juv'),lty=1, col=c('red', 'blue'), bty='n', horiz = T)
abline(v=c(2003,2010, 2013, 2016), lty=3)
# all stages, cod spring GB EPU
plot(NULL, ylim=c(0.2,0.7), xlim=c(1977,2019), ylab="GB spr cod habitat", xlab="", las=1)
lines(xhabich2s_gbk[1,]~yrlist, lty=1, col='green')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(xhabjuv2s_gbk[1,]~yrlist, lty=1, col='blue')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(xhabadt2s_gbk[1,]~yrlist, lty=1, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
legend('topleft', legend = c('adt', 'juv', 'ich'),lty=1, col=c('red', 'blue', 'green'), bty='n', horiz = T)
# abline(v=c(2003,2010, 2013, 2016), lty=2)
# all stages, cod and haddock fall
plot(NULL, ylim=c(0,0.3), xlim=c(1977,2019), ylab="GB fall cod habitat", xlab="", las=1)
lines(xhabjuv2f_gbk[1,]~yrlist, lty=1, col='blue')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(xhabadt2f_gbk[1,]~yrlist, lty=1, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
legend('topleft', legend = c('adt', 'juv'),lty=1, col=c('red', 'blue'), bty='n', horiz = T)

### GOM ###
# all stages TS, haddock spring TS GOM strata
plot(NULL, ylim=c(0,0.5), xlim=c(1977,2019), ylab="GOM spr haddock habitat", xlab="", las=1)
lines(colMeans(xhabich1s_tsGOM)~yrlist, lty=1, col='green')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(colMeans(xhabjuv1s_tsGOM)~yrlist, lty=1, col='blue')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(colMeans(xhabadt1s_tsGOM)~yrlist, lty=1, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
abline(v=c(2003,2010, 2013, 2016), lty=3)
legend('topleft', legend = c('adt', 'juv', 'ich'),lty=1, col=c('red', 'blue', 'green'), bty='n', horiz = T)
# all stages, haddock fall
plot(NULL, ylim=c(0.2,0.5), xlim=c(1977,2019), ylab="GOM fall haddock habitat", xlab="", las=1)
lines(colMeans(xhabjuv1f_tsGOM)~yrlist, lty=1, col='blue')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(colMeans(xhabadt1f_tsGOM)~yrlist, lty=1, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
legend('topleft', legend = c('adt', 'juv'),lty=1, col=c('red', 'blue'), bty='n', horiz = T)
abline(v=c(2003,2010, 2013, 2016), lty=3)
# all stages, cod spring TS GOM strata
plot(NULL, ylim=c(0.15,0.45), xlim=c(1977,2019), ylab="GOM spr cod habitat", xlab="", las=1)
lines(colMeans(xhabich2s_tsGOM)~yrlist, lty=1, col='green')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(colMeans(xhabjuv2s_tsGOM)~yrlist, lty=1, col='blue')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(colMeans(xhabadt2s_tsGOM)~yrlist, lty=1, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
# abline(v=c(2003,2010, 2013, 2016), lty=2)
legend('topleft', legend = c('adt', 'juv', 'ich'),lty=1, col=c('red', 'blue', 'green'), bty='n', horiz = T)
# all stages, cod fall
plot(NULL, ylim=c(0,0.5), xlim=c(1977,2019), ylab="GOM fall cod habitat", xlab="", las=1)
lines(colMeans(xhabjuv2f_tsGOM)~yrlist, lty=1, col='blue')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(colMeans(xhabadt2f_tsGOM)~yrlist, lty=1, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
legend('topleft', legend = c('adt', 'juv'),lty=1, col=c('red', 'blue'), bty='n', horiz = T)


# all stages haddock spring gom EPU
plot(NULL, ylim=c(0,0.3), xlim=c(1977,2019), ylab="GOM spr haddock habitat", xlab="", las=1)
lines(xhabich1s_gom[1,]~yrlist, lty=1, col='green')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(xhabjuv1s_gom[1,]~yrlist, lty=1, col='blue')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(xhabadt1s_gom[1,]~yrlist, lty=1, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
abline(v=c(2003,2010, 2013, 2016), lty=3)
legend('topleft', legend = c('adt', 'juv', 'ich'),lty=1, col=c('red', 'blue', 'green'), bty='n', horiz = T)
# all stages haddock fall
plot(NULL, ylim=c(0.2,0.5), xlim=c(1977,2019), ylab="GOM fall haddock habitat", xlab="", las=1)
lines(xhabjuv1f_gom[1,]~yrlist, lty=1, col='blue')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(xhabadt1f_gom[1,]~yrlist, lty=1, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
legend('topleft', legend = c('adt', 'juv'),lty=1, col=c('red', 'blue'), bty='n', horiz = T)
abline(v=c(2003,2010, 2013, 2016), lty=3)
# all stages, cod spring GOM EPU
plot(NULL, ylim=c(0.1,0.4), xlim=c(1977,2019), ylab="GOM spr cod habitat", xlab="", las=1)
lines(xhabich2s_gom[1,]~yrlist, lty=1, col='green')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(xhabjuv2s_gom[1,]~yrlist, lty=1, col='blue')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(xhabadt2s_gom[1,]~yrlist, lty=1, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
legend('topleft', legend = c('adt', 'juv', 'ich'),lty=1, col=c('red', 'blue', 'green'), bty='n', horiz = T)
# abline(v=c(2003,2010, 2013, 2016), lty=2)
# all stages, cod and haddock fall
plot(NULL, ylim=c(0,0.5), xlim=c(1977,2019), ylab="GOM fall cod habitat", xlab="", las=1)
lines(xhabjuv2f_gom[1,]~yrlist, lty=1, col='blue')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(xhabadt2f_gom[1,]~yrlist, lty=1, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
legend('topleft', legend = c('adt', 'juv'),lty=1, col=c('red', 'blue'), bty='n', horiz = T)


## Get areas of intervals for time series (non projected, 0.1 degree resolution)
# sum(ras[] >= 0.1 & ras[] <= 0.2) * res(ras)[1]^2
test=ich1s[[35]]
crs(test)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
crs(test)
tt=raster::area(test)
testproj=projectRaster(test, crs="+init=epsg:2163")
equalarea = "+proj=laea +lat_0=41 +lon_0=-70"
testproj=projectRaster(test, crs=equalarea, method = "bilinear")
plot(testproj)
res(testproj)
crs(testproj)
(res(testproj)[1]/1000) ^2
res(testproj)[1]*res(testproj)[2]/1e6


stest=sum(testproj[] >= 0.5, na.rm=T)* (res(testproj)[1]/1000) ^2
sum(testproj[] >= 0, na.rm=T)

result <- data.frame(raster::extract(tt, gbk, cellnumbers = T))


geosphere::distHaversine(c(-65,45), c(-65, 45.1))
geosphere::distHaversine(c(-65,35), c(-65, 35.1))
geosphere::distHaversine(c(-65,45), c(-65.1, 45))
geosphere::distHaversine(c(-65,35), c(-65.1, 35))
coordinates(test, 300)

result <- data.frame(raster::extract(test, gbk, cellnumbers = T))
raster::coordinates(test)[result[,1],]
t=cbind(result,coordinates(test)[result[,1],])


### solution: reclassify and extract area (meters^2)
RasterClassAreaSum=function(r, minval,lobrk,hibrk,maxval){
  # set crs to lat long non projected, non spec'd crs string
  if (is.na(crs(r))){
    crs(r)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  }
  if (missing(minval)){
    mnv=cellStats(r, min, na.rm=T)
  } else {
    mnv=minval
  }
  if (missing(lobrk)){
    lov=0.2
  } else {
    lov=lobrk
  }
  if (missing(hibrk)){
    mdv=0.5
  } else{
    mdv=hibrk
  }
  if (missing(maxval)){
    mxv=cellStats(r, max, na.rm=T)
  } else{
    mxv=maxval
  }
  reclass_df=c(mnv,lov,1,lov,mdv,2,mdv,mxv,3)
  # reclass_df=c(0,0.2,1,0.2,0.5,2,0.5,1,3)
  reclass_m <- matrix(reclass_df,
                    ncol = 3,
                    byrow = TRUE)
  # equalarea = "+proj=laea +lat_0=41 +lon_0=-70"
  # r=projectRaster(r, crs=equalarea, method = "bilinear")
  testclass=raster::reclassify(r, reclass_m)
  t=tapply(raster::area(r), testclass[], sum, na.rm=T)
  t2=as.matrix(t, nrow=1, ncol=3)
  return(t2)
}
## haddock area extractions
ichhads.area=matrix(NA, ncol=3, nrow=length(yrlist))
for(i in 1:length(yrlist)){
ichhads.area[i,1:3]=RasterClassAreaSum(ich1s[[i]], 0, 0.2, 0.5, 1)
}
wd1='/home/ryan/Git/NEhabitat/habitat index/'
save(ichhads.area, file=paste(wd1,'ichhads_areaclass.rda', sep=''))
plot(ichhads.area[,3]~yrlist, type='b', las=1,xlab='',ylab='', main='ich had spr')

juvhads.area=matrix(NA, ncol=3, nrow=length(yrlist))
for(i in 1:length(yrlist)){
  juvhads.area[i,1:3]=RasterClassAreaSum(juv1s[[i]], 0, 0.2, 0.5, 1)
}
plot(juvhads.area[,3]~yrlist, type='b', las=1,xlab='',ylab='',main='juv had spr')

adthads.area=matrix(NA, ncol=3, nrow=length(yrlist))
for(i in 1:length(yrlist)){
  adthads.area[i,1:3]=RasterClassAreaSum(adt1s[[i]], 0, 0.2, 0.5, 1)
}
plot(adthads.area[,3]~yrlist, type='b', las=1,xlab='',ylab='',main='adt had spr')
# par(mar=c(1,1,1,1))
par(oma=c(3,2,2,0))
par(mar=c(2,4,2,1) + 0.1)
ylab.txt=expression('Spring habitat area (km)'^2)
plot(NULL, ylim=c(10000,45000), xlim=c(1977,2019), ylab='', xlab="", las=1)
lines(juvhads.area[,3]~yrlist, lty=2, col='black')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
m = lm(juvhads.area[,3]~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=2, col='black')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
lines(adthads.area[,3]~yrlist, lty=1, col='black')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
m = lm(adthads.area[,3]~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=1, col='black')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
lines(ichhads.area[,3]~yrlist, lty=3, col='black')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
m = lm(ichhads.area[,3]~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=3, col='black')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
mtext(ylab.txt,side=2, line =3.5)
legend('topleft', legend = c('Ad had', 'Jv had', 'Ic had'),lty=c(1,2,3), col='black', bty='n', horiz = T)
abline(v=c(2003,2010, 2013, 2016), lty=3)

## cod area extractions
ichcods.area=matrix(NA, ncol=3, nrow=length(yrlist))
for(i in 1:length(yrlist)){
  ichcods.area[i,1:3]=RasterClassAreaSum(ich2s[[i]], 0, 0.2, 0.5, 1)
}
plot(ichcods.area[,3]~yrlist, type='b', las=1,xlab='',ylab='',main='ich cod spr')

juvcods.area=matrix(NA, ncol=3, nrow=length(yrlist))
for(i in 1:length(yrlist)){
  juvcods.area[i,1:3]=RasterClassAreaSum(juv2s[[i]], 0, 0.2, 0.5, 1)
}
plot(juvcods.area[,3]~yrlist, type='b', las=1,xlab='',ylab='',main='juv cod spr')

adtcods.area=matrix(NA, ncol=3, nrow=length(yrlist))
for(i in 1:length(yrlist)){
  adtcods.area[i,1:3]=RasterClassAreaSum(adt2s[[i]], 0, 0.2, 0.5, 1)
}
plot(adtcods.area[,3]~yrlist, type='b', las=1,xlab='',ylab='',main='adt cod spr')

par(oma=c(3,2,2,0))
par(mar=c(2,4,2,1) + 0.1)
ylab.txt=expression('Spring habitat area (km)'^2)
plot(NULL, ylim=c(10000,65000), xlim=c(1977,2019), ylab='', xlab="", las=1)
# plot(NULL, ylim=c(10000,65000), xlim=c(1977,2019), ylab='',main="Spr cod habitat area", xlab="", las=1)
lines(juvcods.area[,3]~yrlist, lty=2, col='black')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
m = lm(juvcods.area[,3]~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=2, col='black')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
lines(adtcods.area[,3]~yrlist, lty=1, col='black')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
m = lm(adtcods.area[,3]~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=1, col='black')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
lines(ichcods.area[,3]~yrlist, lty=3, col='black')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
m = lm(ichcods.area[,3]~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=3, col='black')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
mtext(ylab.txt,side=2, line =3.5)
legend('bottomleft', legend = c('Ad cod', 'Jv cod', 'Ic cod'),lty=c(1,2,3), col='black', bty='n', horiz = T)

# abline(v=c(2003,2010, 2013, 2016), lty=3)
plot(NULL, ylim=c(10000,110000), xlim=c(1977,2019), ylab='',main="Spr cod habitat area", xlab="", las=1)
lines(juvcods.area[,3]~yrlist, lty=1, col='blue')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(adtcods.area[,3]~yrlist, lty=1, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(ichcods.area[,3]~yrlist, lty=1, col='green')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(juvcods.area[,2]~yrlist, lty=2, col='blue')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(adtcods.area[,2]~yrlist, lty=2, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(ichcods.area[,2]~yrlist, lty=2, col='green')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
legend('bottomleft', legend = c('adt', 'juv', 'ich'),lty=1, col=c('red', 'blue', 'green'), bty='n', horiz = T)


### looks bad... do 2 figs instead
# par(oma=c(2,2,2,0))
# par(mar=c(1,4,1,1) + 0.1)
# ylab.txt=expression('Spring habitat area (km)'^2)
# plot(NULL, ylim=c(10000,65000), xlim=c(1977,2019), ylab='', xlab="", las=1)
# # plot(NULL, ylim=c(10000,65000), xlim=c(1977,2019), ylab='',main="Spr cod habitat area", xlab="", las=1)
# lines(juvcods.area[,3]~yrlist, lty=2, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
# lines(adtcods.area[,3]~yrlist, lty=1, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
# lines(ichcods.area[,3]~yrlist, lty=3, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
# mtext(ylab.txt,side=2, line =3.5)
# lines(juvhads.area[,3]~yrlist, lty=2, col='black')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
# lines(adthads.area[,3]~yrlist, lty=1, col='black')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
# lines(ichhads.area[,3]~yrlist, lty=3, col='black')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
# legend('bottomleft', legend = c('Ad cod', 'Jv cod', 'Ic cod'),lty=c(1,2,3), col='black', bty='n', horiz = T)


### Fall haddock
juvhadf.area=matrix(NA, ncol=3, nrow=length(yrlist))
for(i in 1:length(yrlist)){
  juvhadf.area[i,1:3]=RasterClassAreaSum(juv1f[[i]], 0, 0.2, 0.5, 1)
}
plot(juvhadf.area[,3]~yrlist, type='b', las=1,xlab='',ylab='',main='juv had fall')

adthadf.area=matrix(NA, ncol=3, nrow=length(yrlist))
for(i in 1:length(yrlist)){
  adthadf.area[i,1:3]=RasterClassAreaSum(adt1f[[i]], 0, 0.2, 0.5, 1)
}
plot(adthadf.area[,3]~yrlist, type='b', las=1,xlab='',ylab='',main='adt had fall')

par(oma=c(3,2,2,0))
par(mar=c(2,4,2,1) + 0.1)
ylab.txt=expression('Fall habitat area (km)'^2)
plot(NULL, ylim=c(30000,70000), xlim=c(1977,2019), ylab='', xlab="", las=1)
lines(juvhadf.area[,3]~yrlist, lty=2, col='black')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
m = lm(juvhadf.area[,3]~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=2, col='black')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
lines(adthadf.area[,3]~yrlist, lty=1, col='black')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
m = lm(adthadf.area[,3]~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=1, col='black')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
mtext(ylab.txt,side=2, line =3.5)
legend('topleft', legend = c('Ad had', 'Jv had'),lty=c(1,2), col='black', bty='n', horiz = T)
# abline(v=c(2003,2010, 2013, 2016), lty=3)

## Cod fall
juvcodf.area=matrix(NA, ncol=3, nrow=length(yrlist))
for(i in 1:length(yrlist)){
  juvcodf.area[i,1:3]=RasterClassAreaSum(juv2f[[i]], 0, 0.2, 0.5, 1)
}
plot(juvcodf.area[,3]~yrlist, type='b', las=1,xlab='',ylab='',main='juv cod fall')

adtcodf.area=matrix(NA, ncol=3, nrow=length(yrlist))
for(i in 1:length(yrlist)){
  adtcodf.area[i,1:3]=RasterClassAreaSum(adt2f[[i]], 0, 0.2, 0.5, 1)
}
plot(adtcodf.area[,3]~yrlist, type='b', las=1,xlab='',ylab='',main='adt cod fall')

ylab.txt=expression('Fall habitat area (km)'^2)
par(oma=c(3,2,2,0))
par(mar=c(2,4,2,1) + 0.1)
plot(NULL, ylim=c(0,55000), xlim=c(1977,2019), ylab='', xlab="", las=1)
# plot(NULL, ylim=c(0,55000), xlim=c(1977,2019), ylab='',main="Fall cod habitat area", xlab="", las=1)
lines(juvcodf.area[,3]~yrlist, lty=2, col='black')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
m = lm(juvcodf.area[,3]~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=2, col='black')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
## plot regimes
# qlr <- Fstats(juvcodf.area[,3]~yrlist) # 2007 sup.F = 14.528, p-value = 0.01485
# x=breakpoints(qlr)
# fm1 <- lm(juvcodf.area[,3] ~ breakfactor(x, breaks = 2))
# lines(ts(fitted(fm1), start = 1977), col = 'black', lty=2, lwd=2)
lines(adtcodf.area[,3]~yrlist, lty=1, col='black')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
m = lm(adtcodf.area[,3]~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=1, col='black')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
# qlr <- Fstats(adtcodf.area[,3]~yrlist) # 2007 sup.F = 14.528, p-value = 0.01485
# x=breakpoints(qlr)
# fm1 <- lm(adtcodf.area[,3] ~ breakfactor(x, breaks = 1))
# lines(ts(fitted(fm1), start = 1977), col = 'black', lty=2, lwd=2)
mtext(ylab.txt,side=2, line =3.5)
legend('bottomleft', legend = c('Ad cod', 'Jv cod'),lty=c(1,2), col=c('black', 'black'), bty='n', horiz = T)

# abline(v=c(2003,2010, 2013, 2016), lty=3)
plot(NULL, ylim=c(0,55000), xlim=c(1977,2019), ylab='',main="Fall cod habitat area", xlab="", las=1)
lines(juvcodf.area[,2]~yrlist, lty=2, col='black')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(adtcodf.area[,2]~yrlist, lty=2, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(juvcodf.area[,3]~yrlist, lty=1, col='black')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
lines(adtcodf.area[,3]~yrlist, lty=1, col='red')# las=1, ylab=paste(sellab,' GB habitat', xlab='')
legend('bottomleft', legend = c('adt', 'juv'),lty=1, col=c('red', 'black'), bty='n', horiz = T)
abline(v=2004, lty=3)

nes=rgdal::readOGR('/home/ryan/Desktop/shapefiles/epu_shapes/EPU_NESPoly.shp')
gbk=rgdal::readOGR("/home/ryan/Desktop/shapefiles/epu_shapes/EPU_GBKPoly.shp")
gom=rgdal::readOGR("/home/ryan/Desktop/shapefiles/epu_shapes/EPU_GOMPoly.shp")
mab=rgdal::readOGR("/home/ryan/Desktop/shapefiles/epu_shapes/EPU_MABPoly.shp")
scs=rgdal::readOGR("/home/ryan/Desktop/shapefiles/epu_shapes/EPU_SCSPoly.shp")
#combine shapefiles GOM and GBK
gom.scs.shp=rgeos::gUnion(gom, scs, byid=F, id=NULL)
mab.gbk.shp=rgeos::gUnion(mab, gbk, byid=F, id=NULL)
NES.shp=rgeos::gUnion(mab.gbk.shp, gom.scs.shp, byid=F, id=NULL)

#get area of NES
library(sf); library(dplyr)
p <- st_read("/home/ryan/Desktop/shapefiles/epu_shapes/EPU_NESPoly.shp")
p %>% 
  mutate(crown_area_m2 = st_area(p))

# totalarea=raster::extract(ich1s[[1]], NES.shp) 2631 pixels inside NES shape
reclass_df=c(mnv,lov,1,lov,mdv,2,mdv,mxv,3)
# reclass_df=c(0,0.2,1,0.2,0.5,2,0.5,1,3)
reclass_m <- matrix(reclass_df,
                    ncol = 3,
                    byrow = TRUE)
testclass=raster::reclassify(r, reclass_m)
t=tapply(raster::area(r), testclass[], sum, na.rm=T)

zlist=list.files(paste('/home/ryan/Git/NEhabitat/rasters/Spr/ctyp/', sep=''), pattern = 'RAST_ctypZZ_')
wd3=paste('/home/ryan/Git/NEhabitat/rasters/Spr/ctyp', sep='')
zlist=list.files(paste('/home/ryan/Git/NEhabitat/rasters/Fall/ctyp/', sep=''), pattern = 'RAST_ctypZZ_')
wd3=paste('/home/ryan/Git/NEhabitat/rasters/Fall/ctyp', sep='')
rastZT=loadRData(paste(wd3,'/',zlist[1], sep=''))
for (i in 2:length(zlist)){
  rastDF=loadRData(paste(wd3,'/',zlist[i], sep=''))
  rastZT=stack(rastZT, rastDF)
}
## save output
# Zoospring=rastZT
# Zoofall=rastZT

## try calfin to use on extractions of ich
zlist=list.files(paste('/home/ryan/Git/NEhabitat/rasters/Spr/cfin/', sep=''), pattern = 'RAST_calfin_')
wd3=paste('/home/ryan/Git/NEhabitat/rasters/Spr/cfin', sep='')
zlist=list.files(paste('/home/ryan/Git/NEhabitat/rasters/Fall/cfin/', sep=''), pattern = 'RAST_calfin_')
wd3=paste('/home/ryan/Git/NEhabitat/rasters/Fall/cfin', sep='')
rastZT=loadRData(paste(wd3,'/',zlist[1], sep=''))
for (i in 2:length(zlist)){
  rastDF=loadRData(paste(wd3,'/',zlist[i], sep=''))
  rastZT=stack(rastZT, rastDF)
}
## save output
Zoospring2=rastZT
Zoofall2=rastZT

### load and stack Surface temperaure
stlist=list.files(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/ST2', sep=''), pattern = 'RAST_NESREG_')
wd3=paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/ST2', sep='')
rastST=loadRData(paste(wd3,'/',stlist[1], sep=''))
for (i in 2:length(stlist)){
  rastDF=loadRData(paste(wd3,'/',stlist[i], sep=''))
  rastST=stack(rastST, rastDF)
}
## save output
STspring=rastST
STfall=rastST

### load and stack Bottom temperaure
# btlist=list.files(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/BT2', sep=''))
btlist=list.files(paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/BT2', sep=''), pattern = 'RAST_NESREG_')
wd3=paste('/home/ryan/Git/NEhabitat/rasters/', SEASON,'/BT2', sep='')
rastBT=loadRData(paste(wd3,'/',btlist[1], sep=''))
for (i in 2:length(btlist)){
  rastDF=loadRData(paste(wd3,'/',btlist[i], sep=''))
  rastBT=stack(rastBT, rastDF)
}
## save output
# BTspring=rastBT
# BTfall=rastBT
pdf(paste(wd3, '/BT_trends.pdf', sep=''), height=4, width=6)
plotRasterTrends(rastBT)
dev.off()

BT_gbk=raster::extract(rastBT, gbk, fun=mean, na.rm=T)
plot(BT_gbk[1,]~yrlist, type='l')
BT_gom=raster::extract(rastBT, gom, fun=mean, na.rm=T)
plot(BT_gom[1,]~yrlist, type='l')
BT_mab=raster::extract(rastBT, mab, fun=mean, na.rm=T)
plot(BT_mab[1,]~yrlist, type='l')

BT_gbk=raster::extract(rastBT, tsGBK, fun=mean, na.rm=T)
plot(BT_gbk[1,]~yrlist, type='l')
BT_gom=raster::extract(rastBT, tsGOM, fun=mean, na.rm=T)
plot(BT_gom[1,]~yrlist, type='l')
# BT_mab=raster::extract(rastBT, mab, fun=mean, na.rm=T)
# plot(BT_mab[1,]~yrlist, type='l')

BT_nes.fall=raster::extract(rastBT, nes, fun=mean, na.rm=T)
plot(BT_nes.fall[1,]~yrlist, type='l')

## subset BT based on habitat -> now using function
# sBTcodhab=matrix(NA, ncol=1, nrow=length(yrlist))
# for(i in 1:length(yrlist)){
# test=adt2s[[i]]
# test2=Which(test>0.5, cells=T) #cell numbers
# # test3=rasterFromCells(rastBT[[1]], test2, values=T)
# testXY=xyFromCell(test, test2)
# sBTcodhab[i,1]=mean(extract(rastBT[[i]],testXY, fun=mean, na.rm=T))
# }
# plot(sBTcodhab~yrlist, type='b')
extractrasterAbyrasterBvals=function(stckA, stckB, Bval){
matdata=matrix(NA, ncol=1, nrow=length(yrlist))
for(i in 1:dim(stckA)[3]){
  test=stckB[[i]] # the raster stack to use for subsetting
  test2=Which(test>Bval, cells=T) #cell numbers of rasterstack B > Bval
  testXY=xyFromCell(test, test2) # coords
  matdata[i,1]=mean(extract(stckA[[i]],testXY, fun=mean, na.rm=T)) # extract from raster
}
return(matdata)
}

sBTadtcodhab=extractrasterAbyrasterBvals(BTspring, adt2s, Bval=0.5)
sBTjuvcodhab=extractrasterAbyrasterBvals(BTspring, juv2s, Bval=0.5)
sBTichcodhab=extractrasterAbyrasterBvals(BTspring, ich2s, Bval=0.5)
sBTadthadhab=extractrasterAbyrasterBvals(BTspring, adt1s, Bval=0.5)
sBTjuvhadhab=extractrasterAbyrasterBvals(BTspring, juv1s, Bval=0.5)
sBTichhadhab=extractrasterAbyrasterBvals(BTspring, ich1s, Bval=0.5)
fBTadtcodhab=extractrasterAbyrasterBvals(BTfall, adt2f, Bval=0.5)
fBTjuvcodhab=extractrasterAbyrasterBvals(BTfall, juv2f, Bval=0.5)
fBTadthadhab=extractrasterAbyrasterBvals(BTfall, adt1f, Bval=0.5)
fBTjuvhadhab=extractrasterAbyrasterBvals(BTfall, juv1f, Bval=0.5)

## plot extracted T for high quality habitat
plot(sBTadtcodhab~yrlist, type='l', ylim=c(3,7), main='Spr Cod', ylab='BT')
lines(sBTjuvcodhab~yrlist, lty=3)
legend('topleft', legend = c('adt', 'juv'),lty=c(1,3), bty='n', horiz = T)

plot(fBTadtcodhab~yrlist, type='l', ylim=c(8.5,11.5), main="Fall Cod", ylab='BT', las=1)
lines(fBTjuvcodhab~yrlist, lty=3)
legend('topleft', legend = c('adt', 'juv'),lty=c(1,3), bty='n', horiz = T)

plot(sBTadthadhab~yrlist, type='l', ylim=c(4,8), main='Spr Haddock', ylab='BT')
lines(sBTjuvhadhab~yrlist, lty=3)
legend('topleft', legend = c('adt', 'juv'),lty=c(1,3), bty='n', horiz = T)

plot(fBTadthadhab~yrlist, type='l', ylim=c(7.5,11.5), main="Fall Haddock", ylab='BT', las=1)
lines(fBTjuvhadhab~yrlist, lty=3)
legend('topleft', legend = c('adt', 'juv'),lty=c(1,3), bty='n', horiz = T)

## both spr cod and haddock on same plot
plot(sBTadtcodhab~yrlist, type='l', ylim=c(3,8), main='Spr', ylab='BT', xlab='', las=1)
m = lm(sBTadtcodhab~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=1, col='black')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
lines(sBTjuvcodhab~yrlist, lty=2, col='black')
m = lm(sBTjuvcodhab~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=2, col='black')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
lines(sBTichcodhab~yrlist, lty=3, col='black')
m = lm(sBTichcodhab~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=3, col='black')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
lines(sBTadthadhab~yrlist, lty=1, col='red')
m = lm(sBTadthadhab~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=1, col='red')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
lines(sBTjuvhadhab~yrlist, lty=2, col='red')
m = lm(sBTjuvhadhab~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=2, col='red')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
lines(sBTichhadhab~yrlist, lty=3, col='red')
m = lm(sBTichhadhab~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=3, col='red')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
# legend('topleft', legend = c('Ad cod', 'Jv cod', 'Ad had', 'Jv had'),lty=c(1,3,1,3), col=c('black', 'black', 'red','red'), bty='n', horiz = T)
# legend('topleft', legend = c('Ad cod', 'Jv cod', 'Ic cod', 'Ad had', 'Jv had', 'Ic had'),lty=c(1,2,3,1,2,3), col=c('black', 'black', 'black', 'red','red','red'), cex=0.7, bty='n', horiz = T)
# legend('topleft', legend = c('Ad cod', 'Jv cod', 'Ic cod', 'Ad had', 'Jv had', 'Ic had'),lty=c(1,2,3,1,2,3), col=c('black', 'black', 'black', 'red','red','red'), cex=0.7, bty='n', ncol=3)
legend('topleft', legend = c('Ad cod', 'Ad had', 'Jv cod', 'Jv had', 'Ic cod', 'Ic had'),lty=c(1,1,2,2,3,3), col=c('black', 'red', 'black', 'red','black', 'red'), bty='n', ncol=3)

## both fall cod and haddock on same plot
plot(fBTadtcodhab~yrlist, type='l', ylim=c(7.5,11.5), main="Fall", ylab='BT', xlab='', las=1)
m = lm(fBTadtcodhab~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=1, col='black')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
  }
lines(fBTjuvcodhab~yrlist, lty=2)
m = lm(fBTjuvcodhab~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=2, col='black')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
lines(fBTadthadhab~yrlist, lty=1, col='red')
m = lm(fBTadthadhab~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=1, col='red')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
lines(fBTjuvhadhab~yrlist, lty=2, col='red')
m = lm(fBTjuvhadhab~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=2, col='red')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
# legend('topleft', legend = c('Ad cod', 'Jv cod', 'Ad had', 'Jv had'),lty=c(1,2,1,2), col=c('black', 'black', 'red','red'), cex=0.7, bty='n', horiz = T)
legend('topleft', legend = c('Ad cod', 'Ad had', 'Jv cod', 'Jv had'),lty=c(1,1,2,2), col=c('black', 'red', 'black', 'red'), bty='n', ncol=2)


sSTadtcodhab=extractrasterAbyrasterBvals(STspring, adt2s, Bval=0.5)
sSTjuvcodhab=extractrasterAbyrasterBvals(STspring, juv2s, Bval=0.5)
sSTichcodhab=extractrasterAbyrasterBvals(STspring, ich2s, Bval=0.5)
sSTadthadhab=extractrasterAbyrasterBvals(STspring, adt1s, Bval=0.5)
sSTjuvhadhab=extractrasterAbyrasterBvals(STspring, juv1s, Bval=0.5)
sSTichhadhab=extractrasterAbyrasterBvals(STspring, ich1s, Bval=0.5)
fSTadtcodhab=extractrasterAbyrasterBvals(STfall, adt2f, Bval=0.5)
fSTjuvcodhab=extractrasterAbyrasterBvals(STfall, juv2f, Bval=0.5)
fSTadthadhab=extractrasterAbyrasterBvals(STfall, adt1f, Bval=0.5)
fSTjuvhadhab=extractrasterAbyrasterBvals(STfall, juv1f, Bval=0.5)

## both spr cod and haddock on same plot
plot(sSTadtcodhab~yrlist, type='l', ylim=c(2,7), main='Spr', ylab='ST', xlab='', las=1)
m = lm(sSTadtcodhab~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=1, col='black')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
lines(sSTjuvcodhab~yrlist, lty=2, col='black')
m = lm(sSTjuvcodhab~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=2, col='black')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
lines(sSTichcodhab~yrlist, lty=3, col='black')
m = lm(sSTichcodhab~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=3, col='black')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
lines(sSTadthadhab~yrlist, lty=1, col='red')
m = lm(sSTadthadhab~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=1, col='red')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
lines(sSTjuvhadhab~yrlist, lty=2, col='red')
m = lm(sSTjuvhadhab~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=2, col='red')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
lines(sSTichhadhab~yrlist, lty=3, col='red')
m = lm(sSTichhadhab~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=3, col='red')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
# legend('topleft', legend = c('Ad cod', 'Jv cod', 'Ic cod', 'Ad had', 'Jv had', 'Ic had'),lty=c(1,2,3,1,2,3), col=c('black', 'black', 'black', 'red', 'red','red'), cex=0.7, bty='n', horiz = T)
legend('topleft', legend = c('Ad cod', 'Ad had', 'Jv cod', 'Jv had', 'Ic cod', 'Ic had'),lty=c(1,1,2,2,3,3), col=c('black', 'red', 'black', 'red','black', 'red'), bty='n', ncol=3)



## both fall cod and haddock on same plot
plot(fSTadtcodhab~yrlist, type='l', ylim=c(11,16), main="Fall", ylab='ST', xlab='', las=1)
m = lm(fSTadtcodhab~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=1, col='black')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
lines(fSTjuvcodhab~yrlist, lty=2)
m = lm(fSTjuvcodhab~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=2, col='black')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
lines(fSTadthadhab~yrlist, lty=1, col='red')
m = lm(fSTadthadhab~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=1, col='red')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
lines(fSTjuvhadhab~yrlist, lty=2, col='red')
m = lm(fSTjuvhadhab~yrlist)
if(summary(m)$coefficients[8] < 0.05){
  abline(m, lty=2, col='red')
  print(paste('slope= ', round(summary(m)$coefficients[2],3),sep=''))# slope
  print(paste('stde= ', round(summary(m)$coefficients[4],3),sep='')) # std error
  print(paste('p= ', round(summary(m)$coefficients[8],6), sep='')) # pvalue
}
# legend('topleft', legend = c('Ad cod', 'Jv cod', 'Ad had', 'Jv had'),lty=c(1,2,1,2), col=c('black', 'black', 'red','red'), cex=0.7, bty='n', horiz = T)
legend('topleft', legend = c('Ad cod', 'Ad had', 'Jv cod', 'Jv had'),lty=c(1,1,2,2), col=c('black', 'red', 'black', 'red'), bty='n', ncol=2)


sZooadtcodhab=extractrasterAbyrasterBvals(Zoospring, adt2s, Bval=0.5)
sZoojuvcodhab=extractrasterAbyrasterBvals(Zoospring, juv2s, Bval=0.5)
sZooichcodhab=extractrasterAbyrasterBvals(Zoospring, ich2s, Bval=0.5)
fZooadtcodhab=extractrasterAbyrasterBvals(Zoofall, adt2f, Bval=0.5)
fZoojuvcodhab=extractrasterAbyrasterBvals(Zoofall, juv2f, Bval=0.5)
sZooadthadhab=extractrasterAbyrasterBvals(Zoospring, adt1s, Bval=0.5)
sZoojuvhadhab=extractrasterAbyrasterBvals(Zoospring, juv1s, Bval=0.5)
sZooichhadhab=extractrasterAbyrasterBvals(Zoospring, ich1s, Bval=0.5)
fZooadthadhab=extractrasterAbyrasterBvals(Zoofall, adt1f, Bval=0.5)
fZoojuvhadhab=extractrasterAbyrasterBvals(Zoofall, juv1f, Bval=0.5)

### testing extractions on cfin out of curiosity...
# sZooadtcodhab=extractrasterAbyrasterBvals(Zoospring2, adt2s, Bval=0.5)
# sZoojuvcodhab=extractrasterAbyrasterBvals(Zoospring2, juv2s, Bval=0.5)
# sZooichcodhab=extractrasterAbyrasterBvals(Zoospring2, ich2s, Bval=0.5)
# fZooadtcodhab=extractrasterAbyrasterBvals(Zoofall2, adt2f, Bval=0.5)
# fZoojuvcodhab=extractrasterAbyrasterBvals(Zoofall2, juv2f, Bval=0.5)
# sZooadthadhab=extractrasterAbyrasterBvals(Zoospring2, adt1s, Bval=0.5)
# sZoojuvhadhab=extractrasterAbyrasterBvals(Zoospring2, juv1s, Bval=0.5)
# sZooichhadhab=extractrasterAbyrasterBvals(Zoospring2, ich1s, Bval=0.5)
# fZooadthadhab=extractrasterAbyrasterBvals(Zoofall2, adt1f, Bval=0.5)
# fZoojuvhadhab=extractrasterAbyrasterBvals(Zoofall2, juv1f, Bval=0.5)
# plot(sZooadtcodhab~yrlist, type='l', ylim=c(3.7,4.5), main="Spr", ylab='Cfin', xlab='', las=1)
# lines(sZoojuvcodhab~yrlist, lty=2)
# lines(sZooichcodhab~yrlist, lty=3)
# lines(sZooadthadhab~yrlist, lty=1, col='red')
# lines(sZoojuvhadhab~yrlist, lty=2, col='red')
# lines(sZooichhadhab~yrlist, lty=3, col='red')
# legend('topleft', legend = c('Ad cod', 'Jv cod', 'Ic cod', 'Ad had', 'Jv had', 'Ic had'),lty=c(1,2,3,1,2,3), col=c('black', 'black', 'black', 'red', 'red','red'), cex=0.7, bty='n', horiz = T)
# plot(fZooadtcodhab~yrlist, type='l', ylim=c(2.5,4.5), main="Fall", ylab='Cfin', xlab='', las=1)
# lines(fZoojuvcodhab~yrlist, lty=2)
# lines(fZooadthadhab~yrlist, lty=1, col='red')
# lines(fZoojuvhadhab~yrlist, lty=2, col='red')
# legend('topleft', legend = c('Ad cod', 'Jv cod', 'Ad had', 'Jv had'),lty=c(1,2,1,2), col=c('black', 'black', 'red','red'), cex=0.7, bty='n', horiz = T)


plot(sZooadtcodhab~yrlist, type='b')
plot(sZoojuvcodhab~yrlist, type='b')
plot(fZooadtcodhab~yrlist, type='b')
plot(fZoojuvcodhab~yrlist, type='b')
abline(v=2008, lty=3)

### plot extracted Ctypicus abundance for habitat areas
plot(sZooadtcodhab~yrlist, type='l', ylim=c(1,3.5), main="Spr", ylab='Ctyp', xlab='', las=1)
lines(sZoojuvcodhab~yrlist, lty=2)
lines(sZooichcodhab~yrlist, lty=3)
lines(sZooadthadhab~yrlist, lty=1, col='red')
lines(sZoojuvhadhab~yrlist, lty=2, col='red')
lines(sZooichhadhab~yrlist, lty=3, col='red')
# legend('topleft', legend = c('Ad cod', 'Jv cod', 'Ic cod', 'Ad had', 'Jv had', 'Ic had'),lty=c(1,2,3,1,2,3), col=c('black', 'black', 'black', 'red', 'red','red'), cex=0.7, bty='n', horiz = T)
legend('topleft', legend = c('Ad cod', 'Ad had', 'Jv cod', 'Jv had', 'Ic cod',  'Ic had'),lty=c(1,1,2,2,3,3), col=c('black', 'red', 'black', 'red', 'black', 'red'), bty='n', ncol=3)

plot(fZooadtcodhab~yrlist, type='l', ylim=c(3,5), main="Fall", ylab='Ctyp', xlab='', las=1)
lines(fZoojuvcodhab~yrlist, lty=2)
lines(fZooadthadhab~yrlist, lty=1, col='red')
lines(fZoojuvhadhab~yrlist, lty=2, col='red')
# legend('topleft', legend = c('Ad cod', 'Jv cod', 'Ad had', 'Jv had'),lty=c(1,2,1,2), col=c('black', 'black', 'red','red'), cex=0.7, bty='n', horiz = T)
legend('topleft', legend = c('Ad cod', 'Ad had', 'Jv cod', 'Jv had'),lty=c(1,1,2,2), col=c('black','red', 'black', 'red'), bty='n', ncol=2)

### bottom temperature thermal habitat (C) area selections --
## from values and plots above, use for range for (Spr/Fall):
# Adt had 4.5-7 / 8-10
# Juv had 4.5-7 / 9.5-11
# Adt cod 4.5-6.5 / 8-10
# Juv cod 3.0-7 / 9.5-11
## -> spring 3-7
## -> fall 8-11

isEmpty <- function(x) {
  return(identical(x, numeric(0)))
}
thermhabs.area=matrix(NA, ncol=3, nrow=length(yrlist))
for(i in 1:length(yrlist)){
  thermhabs.area2=RasterClassAreaSum(BTspring[[i]], 0, 3, 7, 30)
  thermhabs.area[i,1]=ifelse(isEmpty(thermhabs.area2[which(as.numeric(rownames(thermhabs.area2))==1)]), NA, thermhabs.area2[which(as.numeric(rownames(thermhabs.area2))==1)])
  thermhabs.area[i,2]=ifelse(isEmpty(thermhabs.area2[which(as.numeric(rownames(thermhabs.area2))==2)]), NA, thermhabs.area2[which(as.numeric(rownames(thermhabs.area2))==2)])
  thermhabs.area[i,3]=ifelse(isEmpty(thermhabs.area2[which(as.numeric(rownames(thermhabs.area2))==3)]), NA, thermhabs.area2[which(as.numeric(rownames(thermhabs.area2))==3)])
}
plot(thermhabs.area[,2]~yrlist, type='b', las=1,xlab='',ylab='',main='spr 3-7 area')

thermhabs.area=matrix(NA, ncol=3, nrow=length(yrlist))
for(i in 1:length(yrlist)){
  bbox=extent(-74, -65.45, 40, 44.65)
  test=crop(BTspring[[i]], bbox)
  thermhabs.area2=RasterClassAreaSum(test, 0, 3, 7, 30)
  thermhabs.area[i,1]=ifelse(isEmpty(thermhabs.area2[which(as.numeric(rownames(thermhabs.area2))==1)]), NA, thermhabs.area2[which(as.numeric(rownames(thermhabs.area2))==1)])
  thermhabs.area[i,2]=ifelse(isEmpty(thermhabs.area2[which(as.numeric(rownames(thermhabs.area2))==2)]), NA, thermhabs.area2[which(as.numeric(rownames(thermhabs.area2))==2)])
  thermhabs.area[i,3]=ifelse(isEmpty(thermhabs.area2[which(as.numeric(rownames(thermhabs.area2))==3)]), NA, thermhabs.area2[which(as.numeric(rownames(thermhabs.area2))==3)])
}
plot(thermhabs.area[,2]~yrlist, type='b', las=1,xlab='',ylab='',main='Spr habitat area')

thermhabf.area=matrix(NA, ncol=3, nrow=length(yrlist))
for(i in 1:length(yrlist)){
  thermhabf.area2=RasterClassAreaSum(BTspring[[i]], 0, 8, 11, 25)
  thermhabf.area[i,1]=ifelse(isEmpty(thermhabf.area2[which(as.numeric(rownames(thermhabf.area2))==1)]), NA, thermhabf.area2[which(as.numeric(rownames(thermhabf.area2))==1)])
  thermhabf.area[i,2]=ifelse(isEmpty(thermhabf.area2[which(as.numeric(rownames(thermhabf.area2))==2)]), NA, thermhabf.area2[which(as.numeric(rownames(thermhabf.area2))==2)])
  thermhabf.area[i,3]=ifelse(isEmpty(thermhabf.area2[which(as.numeric(rownames(thermhabf.area2))==3)]), NA, thermhabf.area2[which(as.numeric(rownames(thermhabf.area2))==3)])
}
plot(thermhabf.area[,2]~yrlist, type='b', las=1,xlab='',ylab='',main='fall 8-11 area')

thermhabf.area=matrix(NA, ncol=3, nrow=length(yrlist))
for(i in 1:length(yrlist)){
  bbox=extent(-74, -65.45, 40, 44.65)
  test=crop(BTfall[[i]], bbox)
  thermhabf.area2=RasterClassAreaSum(test, 0, 8, 11, 25)
  thermhabf.area[i,1]=ifelse(isEmpty(thermhabf.area2[which(as.numeric(rownames(thermhabf.area2))==1)]), NA, thermhabf.area2[which(as.numeric(rownames(thermhabf.area2))==1)])
  thermhabf.area[i,2]=ifelse(isEmpty(thermhabf.area2[which(as.numeric(rownames(thermhabf.area2))==2)]), NA, thermhabf.area2[which(as.numeric(rownames(thermhabf.area2))==2)])
  thermhabf.area[i,3]=ifelse(isEmpty(thermhabf.area2[which(as.numeric(rownames(thermhabf.area2))==3)]), NA, thermhabf.area2[which(as.numeric(rownames(thermhabf.area2))==3)])
}
plot(thermhabf.area[,2]~yrlist, type='b', las=1,xlab='',ylab='',main='fall 8-11 area')

par(oma=c(3,2,2,0))
par(mar=c(2,4,2,1) + 0.1)
ylab.txt=expression('Thermal habitat area (km)'^2)
plot(thermhabs.area[,2]~yrlist, type='l', ylim=c(5e4, 1.8e5),las=1,xlab='',ylab='')#,main='Thermal habitat area')
lines(thermhabf.area[,2]~yrlist, lty=2, las=1)
mtext(ylab.txt,side=2, line =3.5)
legend('bottomleft', legend = c('Spr: 3-7 C', 'Fall: 8-11 C'),lty=c(1,2), col=c('black', 'black'), bty='n', horiz = T)

### testing breakpoints for extracted BT
library(strucchange)
qlr <- Fstats(fBTadthadhab ~ yrlist) # 2007 sup.F = 14.528, p-value = 0.01485
qlr <- Fstats(fBTjuvhadhab ~ yrlist) # 2007 sup.F = 17.777, p-value = 0.003497
qlr <- Fstats(fBTadtcodhab ~ yrlist) # 2007 sup.F = 11.36, p-value = 0.05712
qlr <- Fstats(fBTjuvcodhab ~ yrlist) # 2007 sup.F = 8.1339, p-value = 0.203

x=breakpoints(qlr)
yrlist[x[[1]]]
sctest(qlr, type = "supF")
plot(qlr)

qlr <- Fstats(sBTadthadhab ~ yrlist) # 2002 sup.F = 11.864, p-value = 0.04636
qlr <- Fstats(sBTjuvhadhab ~ yrlist) # 2002 sup.F = 9.5983, p-value = 0.1161
qlr <- Fstats(sBTadtcodhab ~ yrlist) # 2002 sup.F = 8.92, p-value = 0.151
qlr <- Fstats(sBTjuvcodhab ~ yrlist) # 2002 sup.F = 7.7629, p-value = 0.2325


qlr <- Fstats(sSTadthadhab ~ yrlist) # 2011 sup.F = 5.0774, p-value = 0.5612
qlr <- Fstats(sSTjuvhadhab ~ yrlist) # 1991 sup.F = 5.3306, p-value = 0.5215
qlr <- Fstats(sSTadtcodhab ~ yrlist) # 1991 sup.F = 4.638, p-value = 0.6331
qlr <- Fstats(sSTjuvcodhab ~ yrlist) # 2011 sup.F = 4.7803, p-value = 0.6095

qlr <- Fstats(fSTadthadhab ~ yrlist) # 1987 sup.F = 21.124, p-value = 0.000753
qlr <- Fstats(fSTjuvhadhab ~ yrlist) # 1987 sup.F = 22.131, p-value = 0.0004713
qlr <- Fstats(fSTadtcodhab ~ yrlist) # 1987 sup.F = 25.765, p-value = 8.498e-05
qlr <- Fstats(fSTjuvcodhab ~ yrlist) # 1987 sup.F = 27.3, p-value = 4.089e-05


## plot regimes
qlr <- Fstats(fBTadthadhab ~ yrlist) # 2007 sup.F = 14.528, p-value = 0.01485
x=breakpoints(qlr)
# yrlist[x[[1]]]
# sctest(qlr, type = "supF")

## plot regimes Fall bottom temperature
plot(fBTjuvhadhab ~ yrlist, type='l', lty=3, las=1, ylab='', xlab='', ylim=c(7.5,11.5), col='red')
fm1 <- lm(fBTjuvhadhab ~ breakfactor(x, breaks = 1))
lines(ts(fitted(fm1), start = 1977), col = 'red', lty=3, lwd=2)
lines(fBTadthadhab ~ yrlist, type='l', lty=1, col='red')
fm1 <- lm(fBTadthadhab ~ breakfactor(x, breaks = 1))
lines(ts(fitted(fm1), start = 1977), lty=1, lwd=2, col = 'red')
# legend('topleft', legend = c('Ad had', 'Jv had'),lty=c(1,3), col=c('black', 'black'), bty='n', horiz = T)
lines(fBTadtcodhab ~ yrlist, type='l', lty=1, col='black')
fm1 <- lm(fBTadtcodhab ~ breakfactor(x, breaks = 1))
lines(ts(fitted(fm1), start = 1977), lwd=2, col = 'black')
lines(fBTjuvcodhab ~ yrlist, type='l', lty=3, col='black')
fm1 <- lm(fBTjuvcodhab ~ breakfactor(x, breaks = 1))
# lines(ts(fitted(fm1), start = 1977), lwd=2, lty=3, col = 'black')
legend('topleft', legend = c('Ad cod', 'Jv cod', 'Ad had', 'Jv had'),lty=c(1,3,1,3), col=c('black', 'black', 'red','red'), bty='n', horiz = T, cex=0.8)
abline(v=2008, lty=2)

qlr <- Fstats(fSTadthadhab ~ yrlist) # 1987 sup.F = 21.124, p-value = 0.000753
t=matrix(data=fSTadthadhab)
qlr <- Fstats(ts(t ~ yrlist)) # 1987 sup.F = 21.124, p-value = 0.000753

plot(qlr)
sctest(qlr)
plot(fSTadthadhab~yrlist)
lines(breakpoints(qlr))

qlr <- Fstats(fSTjuvhadhab ~ yrlist)
plot(qlr)
x=breakpoints(fSTjuvhadhab~1)
yrlist[x[[1]]]
plot(x)
plot(fSTjuvhadhab ~ yrlist, type='l', lty=3, las=1, ylab='', xlab='', ylim=c(11,16), col='red')
fm1 <- lm(fSTjuvhadhab ~ breakfactor(x, breaks = 2))
lines(ts(fitted(fm1), start = 1977), col = 'red', lty=3, lwd=2)
lines(confint(fm1, breaks =3))


plot(fSTjuvhadhab ~ yrlist, type='l', lty=3, las=1, ylab='', xlab='', ylim=c(11,16), col='red')
fm1 <- lm(fSTjuvhadhab ~ breakfactor(x, breaks = 1))
lines(ts(fitted(fm1), start = 1977), col = 'red', lty=3, lwd=2)
lines(fBTadthadhab ~ yrlist, type='l', lty=1, col='red')
fm1 <- lm(fBTadthadhab ~ breakfactor(x, breaks = 1))
lines(ts(fitted(fm1), start = 1977), lty=1, lwd=2, col = 'red')
# legend('topleft', legend = c('Ad had', 'Jv had'),lty=c(1,3), col=c('black', 'black'), bty='n', horiz = T)
lines(fBTadtcodhab ~ yrlist, type='l', lty=1, col='black')
fm1 <- lm(fBTadtcodhab ~ breakfactor(x, breaks = 1))
lines(ts(fitted(fm1), start = 1977), lwd=2, col = 'black')
lines(fBTjuvcodhab ~ yrlist, type='l', lty=3, col='black')
fm1 <- lm(fBTjuvcodhab ~ breakfactor(x, breaks = 1))
# lines(ts(fitted(fm1), start = 1977), lwd=2, lty=3, col = 'black')
legend('topleft', legend = c('Ad cod', 'Jv cod', 'Ad had', 'Jv had'),lty=c(1,3,1,3), col=c('black', 'black', 'red','red'), bty='n', horiz = T, cex=0.8)
abline(v=2008, lty=2)