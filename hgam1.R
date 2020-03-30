#### Hierarchical GAMs fitting life history stages for NE groundfish
# RM 2020
# 1) Train model w/ data through 1999, predict 2000-2018
# 2) Train 1977-2011, predict 2012-2018
# 3) Train 1977-2012, predict 2013-2018
# 4) Train 1977-2014, predict 2015-2018

library(mgcv)
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

load('/home/ryan/Git/NEhabitat/Final_merged_fish_corABN_Zoo_Ich.Rda') # Biomass
load('/home/ryan/Git/NEhabitat/Final_merged_fish_corABN_Zoo_Ich.Rda') # Abundance
### update colnames for ichythoplankton ### Note: urospp assigned to 76 (WHK) but is not distinguishable from 77 (RHK)
# colnames(FData.bio[50:61])
# [1] "urospp_100m3" "gadmor_100m3" "melaeg_100m3" "polvir_100m3" "merbil_100m3" "sebspp_100m3"
# [7] "anaspp_100m3" "parden_100m3" "pseame_100m3" "glycyn_100m3" "scoaqu_100m3" "lopame_100m3"
colnames(FData.bio)[50:61]=c("76_ich", "73_ich", "74_ich", "75_ich", "72_ich", "155_ich", "192_ich", "103_ich","106_ich","107_ich", "108_ich", "197_ich")
colnames(FData.abn)[50:61]=c("76_ich", "73_ich", "74_ich", "75_ich", "72_ich", "155_ich", "192_ich", "103_ich","106_ich","107_ich", "108_ich", "197_ich")

### select data
# fish=FData.abn %>% select(YEAR, SEASON, LAT:BOTSALIN, `73_Adt`, `73_Juv`, `73_ich`, calfin_100m3:chaeto_100m3) # need to fill Salinity
fish=FData.abn %>% select(YEAR, SEASON, LAT:BOTTEMP, `73_Adt`, `73_Juv`, `73_ich`, calfin_100m3:chaeto_100m3)
fish$MONTH=month(FData.abn$EST_TOWDATE)
fish=fish[complete.cases(fish),]
fish$`73_ich`=ceiling(fish$`73_ich`) # make integer from numbers per 100 m/3
# fish[`73_ich`]=round(fish[`73_ich`], digts=0)
table(fish$YEAR)

## fitting distribution (not working...)
qqp(fish$`73_Adt`, "norm")
qqp(fish$`73_Adt`, "lnorm")
nbinom <- fitdistr(fish$`73_Adt`, "Negative Binomial")
qqp(FData.bio$`73_Adt`, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])


### Need to make data frame wide -> long with just species of interest...
# e.g.: df %>% pivot_longer(c(x, y, z), names_to = "key", values_to = "value")
# pivot_longer(data, cols, names_to = "name", names_prefix = NULL,
#              names_sep = NULL, names_pattern = NULL, names_ptypes = list(),
#              names_repair = "check_unique", values_to = "value",
#              values_drop_na = FALSE, values_ptypes = list())
fish2=fish %>% pivot_longer(c(`73_Adt`, `73_Juv`, `73_ich`), names_to = c("SVSPP", "Stg"), names_sep ="_", values_to = "Number")
fish2$Stg=as.factor(fish2$Stg)
table(fish2$MONTH)

# Global smoother
fishG=gam(log10(fish2$Number+1)~te(fish2$MONTH, fish2$BOTTEMP, log10(fish2$calfin_100m3+1), bs=c("cc", "tp", "tp"), 
                                  k=c(8, 10, 5), m=2), family=Gamma(link="log"), method="REML")
fishG=gam(fish2$Number~te(fish2$MONTH, fish2$BOTTEMP, log10(fish2$calfin_100m3+1), bs=c("cc", "tp", "tp"), 
                                   k=c(8, 10, 5), m=2), family=gaussian(link="identity"), method="REML")
# stage was not included explicitly in above, does this make it GI model?
fishGI=gam(fish2$Number~te(fish2$MONTH, fish2$BOTTEMP, fish2$Stg, bs=c("cc", "tp", "re"), 
                          k=c(8, 10, 5), m=2) + s(log10(fish2$calfin_100m3+1), bs="tp"), family=gaussian(link="identity"), method="REML")

fishGI=gam(fish2$Number~te(fish2$MONTH, fish2$BOTTEMP, bs=c("cc", "tp"), k=c(8, 10)) + s(fish2$Stg, bs="re", k=3) + s(log10(fish2$calfin_100m3+1), bs="tp", k=5), family=gaussian(link="identity"), method="REML")


plot(fishGI)
gam.check(fishGI) #p values are for the test of the null hypothesis that the basis dimension used is of sufficient size
summary.gam(fishGI) #test null hypothesis of a zero effect of the indicated spline

fishGS=gam(log10(fish2$Number+1)~te(fish2$MONTH, fish2$BOTTEMP, fish2$Stg, bs=c("cc", "tp", "tp"), 
                                   k=c(8, 10, 5), m=2), method="REML")



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



