require(raster)
# require(ncdf)
library(sp)
library(maptools)
library(marmap)
library(geosphere)
library(lubridate)
library(dplyr)
library(tidyr)

## Subset survdat to only groundfish species with Length at maturity data
# load("C:/Users/ryan.morse/Downloads/SurvdatBio (3).RData") #loads survdat.bio
load("C:/Users/ryan.morse/Downloads/survdat_lw.RData") # loads survdat.lw ; from Sean Lucey data into 2019
load("/home/ryan/Downloads/survdat_lw.RData")# loads survdat.lw ; from Sean Lucey data into 2019
test=survdat.lw$SVSPP %in% Lmf$SVSPP
survdat.lw=survdat.lw[test,]
sort(unique(survdat.lw$SVSPP)) #verify
length(unique(survdat.lw$SVSPP))

# test=survdat$SVSPP %in% Lmf$SVSPP
# survdat=survdat[test,]
# sort(unique(survdat$SVSPP)) #verify
# length(unique(survdat$SVSPP))

## choose surdat and remove older non-length version to avoid confusion
rm(survdat)
survdat=survdat.lw
survdat=survdat.bio


### loop to split taxa into juvenile or adult using Lm data
# survdat$stg=NA
# for (j in 1:length(survdat$stg)){
#   survdat$stg[j]=ifelse(survdat$LENGTH[j]<Lmf$Lm[which(Lmf$SVSPP==survdat$SVSPP[j])], "juv", "adt")
# }
### Better to do this with left_join
test=data.frame(survdat[,c("SVSPP", "LENGTH")])
colnames(test)=c("SVSPP", "LENGTH")
test=left_join(test, Lmf[,c("SVSPP", "Lm")], by="SVSPP")
test$stg=ifelse(test$LENGTH<test$Lm, "Juv", "Adt")
survdat$stg=test$stg
rm(test)

### count juv and adult to see if numbers match ABUNDANCE 
survdat=survdat %>% group_by_at(vars(SVSPP,CRUISE6,STATION,STRATUM,ABUNDANCE,stg)) %>% mutate(stgsum=sum(NUMLEN), stgwtsum=sum(WGTLEN))
### now addd up weights for all fish measured, check against biomass, total abundance
survdat=survdat %>% group_by_at(vars(SVSPP,CRUISE6,STATION,STRATUM,ABUNDANCE)) %>% mutate(totlen=sum(NUMLEN), totwgt=sum(WGTLEN), abndiff=ABUNDANCE-totlen, biodiff=BIOMASS-totwgt)
### now add corrected biomass and abundance to account for any missing staged values (eg. for large hauls)
survdat=survdat %>% group_by_at(vars(SVSPP,CRUISE6,STATION,STRATUM,ABUNDANCE, stg)) %>%
  mutate(stgpctabn=stgsum/totlen, stgpctwgt=stgwtsum/totwgt, corABN=round(stgpctabn*ABUNDANCE,digits = 0), corBIO=round(stgpctwgt*BIOMASS, digits=2))


# barplot(table(round(survdat$biodiff[(abs(survdat$biodiff>12))],1)))
# y=survdat$SVSPP
# x=abs(survdat$biodiff)
# y=y[x>120]
# x=x[x>120]
# barplot(table(round(y))) #which species have the highest mismatches in biomass
# barplot(table(round(x))) #how many

### now create dataframe with unique tows only, separated by stage
test=survdat[,c("CRUISE6","STATION","STRATUM","YEAR", "SVSPP", "stg")] # needs SVSPP...
svdtunq=survdat[!duplicated(test),]
### calc percent of biomass by stage(adt, juv) per unique tow to be applied to BIOMASS for corrected values
# svdtunq$stgwgtpct=round(svdtunq$stgwtsum/svdtunq$totwgt, 2)
# svdtunq$corBIOMASS=svdtunq$stgwgtpct*svdtunq$BIOMASS

## not working yet
# svdunqadt=svdtunq %>% filter(stg=="Adt") %>% 
#   select(-(LENGTH:SIZECAT)) %>% 
#   select(-(CATCHSEX)) %>% 
#   spread(SVSPP, ABUNDANCE:corBIO)

# barplot(table(svdtunq$stgwgtpct))
# barplot(table(svdtunq$stgwgtpct[svdtunq$SEASON=='SPRING']))
# barplot(table(svdtunq$stgwgtpct[svdtunq$SEASON=='FALL']))
# barplot(table(svdtunq$stgwgtpct[svdtunq$stgwgtpct<0.99 & svdtunq$SEASON=='FALL']),ylim=c(0,1000), main='FALL')
# barplot(table(svdtunq$stgwgtpct[svdtunq$stgwgtpct<0.99 & svdtunq$SEASON=='SPRING']),ylim=c(0,1000), main='SPRING')
# 
# svspp=73
# svdtunq$stgwgtpct[is.na(svdtunq$stgwgtpct)]=0 #[svdtunq$stgwgtpct<0.99 & svdtunq$SEASON=='FALL' & svdtunq$SVSPP==svspp]),]
# 
# barplot(table(svdtunq$stgwgtpct[svdtunq$stgwgtpct<0.99 & svdtunq$SEASON=='FALL' & svdtunq$SVSPP==svspp]),ylim=c(0,100), main=paste('FALL ', svspp))
# barplot(table(svdtunq$stgwgtpct[svdtunq$stgwgtpct<0.99 & svdtunq$SEASON=='SPRING' & svdtunq$SVSPP==svspp]),ylim=c(0,100), main=paste('SPRING ', svspp))
# 
# 
# library(modes)
# # restricting to <0.99 was masking adult population. Can check on change over years...
# bimodality_amplitude(svdtunq$stgwgtpct[svdtunq$SEASON=='SPRING' & svdtunq$SVSPP==svspp], fig=T)
# bimodality_amplitude(svdtunq$stgwgtpct[svdtunq$SEASON=='FALL' & svdtunq$SVSPP==svspp], fig=T)
# sum(is.na(svdtunq$stgwgtpct[svdtunq$stgwgtpct<0.99 & svdtunq$SEASON=='FALL' & svdtunq$SVSPP==svspp]))

### change survdat long to wide, select biomass based (corBIO) or abundance based (corABN) ###
svdtunq=ungroup(svdtunq) #not needed
svdwide.bio=svdtunq %>% select(SVSPP, stg, CRUISE6, STATION, STRATUM,SVVESSEL,YEAR,EST_TOWDATE, SEASON,LAT,LON,DEPTH,
                               SURFTEMP,BOTTEMP,SURFSALIN,BOTSALIN,corBIO) %>%
  pivot_wider(names_from=c(SVSPP, stg), values_from = corBIO, values_fill = list(corBIO=0))

svdwide.abn=svdtunq %>% select(SVSPP, stg, CRUISE6, STATION, STRATUM,SVVESSEL,YEAR,EST_TOWDATE, SEASON,LAT,LON,DEPTH,
                               SURFTEMP,BOTTEMP,SURFSALIN,BOTSALIN,corABN) %>%
  pivot_wider(names_from=c(SVSPP, stg), values_from = corABN, values_fill = list(corABN=0))

#remove column names with 'x_NA'
svdwide.bio=svdwide.bio%>% select(-`72_NA`, -`73_NA`, -`74_NA`,-`75_NA`,-`76_NA`,-`77_NA`,-`101_NA`,-`102_NA`,-`103_NA`,-`105_NA`,
                                  -`106_NA`,-`107_NA`,-`108_NA`,-`155_NA`,-`193_NA`,-`197_NA`)
svdwide.abn=svdwide.abn%>% select(-`72_NA`, -`73_NA`, -`74_NA`,-`75_NA`,-`76_NA`,-`77_NA`,-`101_NA`,-`102_NA`,-`103_NA`,-`105_NA`,
                                  -`106_NA`,-`107_NA`,-`108_NA`,-`155_NA`,-`193_NA`,-`197_NA`)
save(svdwide.bio, file='unique_wide_format_corBIO_by_stage.Rda')
save(svdwide.abn, file='unique_wide_format_corABN_by_stage.Rda')

### now select just haddock and get full list by season
svdsmall=svdwide.abn %>% select(`74_Juv`, `74_Adt`, SEASON, YEAR, EST_TOWDATE, LAT, LON) %>% 
  filter(SEASON=='SPRING', YEAR > 1976)

### Load the data and rename ichthyoplankton
load('/home/ryan/Git/NEhabitat/Final_merged_fish_corBIO_Zoo_Ich.Rda') # Biomass
load('/home/ryan/Git/NEhabitat/Final_merged_fish_corABN_Zoo_Ich.Rda') # Abundance

### Select season for GAMs
slctseason="FALL"; fishseas="Fall"
slctseason="SPRING"; fishseas="Spr"

### Choose fish to run, adjust by searching on and changing " `74_ " to new ID **** ###
### Subset the formatted data for a single species (see Lmf)
# Lmf=read_excel('/home/ryan/Git/NEhabitat/Lm_included.xlsx')
# fishname='Cod' #73
fishname='Haddock' #74
# fishname='SilverHake' #72
# fishname='Pollock' #75

fish=FData.abn %>% dplyr::select(YEAR, SEASON, date, LAT:BOTTEMP, `74_Adt`, `74_Juv`, `74_ich`, volume_100m3:chl12)
fish$MONTH=month(FData.abn$EST_TOWDATE)
fish$DAY=day(FData.abn$EST_TOWDATE)
fish=fish[complete.cases(fish),]
fish$`74_ich`=ceiling(fish$`74_ich`) # make integer from numbers per 100 m/3
fish2=fish[which(fish$SEASON==slctseason),] # subset to season


## to use in HGAM testing residuals for PA
fishtest=fish2
### pivot longer so stage is repeated
fishtest=fishtest %>% pivot_longer(c(`74_Adt`, `74_Juv`, `74_ich`), names_to = c("SVSPP", "Stg"), names_sep ="_", values_to = "Number")
fishtest$Stg=factor(fishtest$Stg, ordered=F)
logd=fishtest[,9:20] #fishtest[,8:19]
logd=log10(logd+1)
# fishtest[,8:19]=logd
fishtest[,9:20]=logd
fishtest$pa=ifelse(fishtest$Number > 0, 1, 0)
# testBIO=testBIO %>%  
#   filter(., Biomass >0) %>%
#   mutate(., "logbio"=log10(Biomass+1))
### drop ich for periods when there are none avail (fall for haddock)
fishtest=fishtest %>% filter(Stg!='ich') %>% droplevels(exclude="ich")

test=left_join(fish2, svdsmall, by=c('LAT', 'LON', 'YEAR'))
## testing length of single year, number of tows expected...
x=survdat.lw
test=data.frame(x[,c("SVSPP", "LENGTH")])
colnames(test)=c("SVSPP", "LENGTH")
test=left_join(test, Lmf[,c("SVSPP", "Lm")], by="SVSPP")
test$stg=ifelse(test$LENGTH<test$Lm, "Juv", "Adt")
x$stg=test$stg
rm(test)

### count juv and adult to see if numbers match ABUNDANCE 
x=x %>% group_by_at(vars(SVSPP,CRUISE6,STATION,STRATUM,ABUNDANCE,stg)) %>% mutate(stgsum=sum(NUMLEN), stgwtsum=sum(WGTLEN))
### now addd up weights for all fish measured, check against biomass, total abundance
x=x %>% group_by_at(vars(SVSPP,CRUISE6,STATION,STRATUM,ABUNDANCE)) %>% mutate(totlen=sum(NUMLEN), totwgt=sum(WGTLEN), abndiff=ABUNDANCE-totlen, biodiff=BIOMASS-totwgt)
### now add corrected biomass and abundance to account for any missing staged values (eg. for large hauls)
x=x %>% group_by_at(vars(SVSPP,CRUISE6,STATION,STRATUM,ABUNDANCE, stg)) %>%
  mutate(stgpctabn=stgsum/totlen, stgpctwgt=stgwtsum/totwgt, corABN=round(stgpctabn*ABUNDANCE,digits = 0), corBIO=round(stgpctwgt*BIOMASS, digits=2))



tt=x %>% filter(YEAR==2005, SVSPP==74, SEASON=='SPRING')
test=tt[,c("CRUISE6","STATION","LAT", "LON", "STRATUM","YEAR", "SVSPP", "TOW", "stg")] # needs SVSPP...

test=tt[,c("LAT", "LON","YEAR", "TOW", "SVSPP", "stg")] # needs SVSPP...

tt2=tt[!duplicated(test),]

length(unique(c(tt2$LON, tt2$LAT)))
### plot samples
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70", mar=c(4,0,0,0))
map.axes(las=1)
points(tt2$LON, tt2$LAT)
# br=c(-Inf, 0,10,50,100,1000)
br=c(-1,0,1,10,100,1000)
# labs=c('0', '1', '>1', '>10', '>100')
labs2=c(-Inf,0,1,2,3) # log10 breaks
labs3=c(0.5, 1, 1.25, 1.5, 2) # plotting size breaks
labs3=c(0.25, 0.5, 1, 1.5, 2) # plotting size breaks

t=fish2 %>% filter(YEAR==2004)
t$bin=cut(t$`74_ich`, breaks=br, labels=labs3, right=T)
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70", mar=c(4,0,0,0))
points(t$LON, t$LAT, pch=1, cex=as.numeric(as.character(t$bin)), col=ifelse(as.numeric(as.character(t$bin))<0.5, 'black', 'red'))
legend(-70, 38, pch=1, pt.cex=labs3, legend=c('0', '1', '>1', '>10', '>100'), col=ifelse(labs3<0.5, 'black', 'red'))
# points(t$LON, t$LAT, pch=19, cex=as.numeric(t$bin)/2, col=ifelse(as.numeric(t$bin)/2<1, 'black', 'red'))
# legend(-70, 38, pch=19, cex=unique(as.numeric(t2))/2, legend=c('0', '< 10', '< 50'))

# checking on frequency years with higher positive tows in the MAB
t=fish2 %>% filter(`74_ich`>0 & LAT<40.5)
table(t$YEAR)
barplot(table(t$YEAR))
t=fish2 %>% filter(`74_ich`>0 & LAT<40)

m=t[rev(order(t$`74_ich`)),]
table(m$YEAR[1:20])

### checking on duplicated data points...
t5=fish2 %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n() ==1) # drops from 9470 to ~3600 samples!
table(t5$YEAR)

#CHECK on single year duplicate stations
t=fish2 %>% filter(YEAR==2004)
t$bin=cut(t$`74_ich`, breaks=br, labels=labs3, right=T)

tdup=t %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()>1) %>% mutate(num=n())
tdupmn=tdup %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()>1) %>% summarise_at(vars(`74_Adt`:chaeto_100m3), mean, na.rm=T) # take mean of all duplicate rows
tdupmd=t %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()>1) %>% summarise_all(funs(med)) # take mean of all duplicate rows

### this works to average duplicates, triplicates, (etc) and select only first (after taking mean of values)
tdupmn=t %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()>1) %>% mutate(num=n()) %>% mutate_if(is.numeric, mean) %>%
  mutate_if(is.character, funs(paste(unique(.), collapse = "_"))) %>% slice(1) #distinct()

tunq=t %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()==1) %>% mutate(num=n())
tdup=t %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()>1) %>% mutate(num=n())
tdupmn=tdup %>% mutate_if(is.numeric, median) %>% mutate_if(is.character, funs(paste(unique(.), collapse = "_"))) %>% slice(1)
tfin=bind_rows(tunq, tdupmn)
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70", mar=c(4,0,0,0))
points(tfin$LON, tfin$LAT, pch=1, cex=as.numeric(as.character(tfin$bin)), col=ifelse(as.numeric(as.character(tfin$bin))<0.5, 'black', 'red'))
legend(-70, 38, pch=1, pt.cex=labs3, legend=c('0', '1', '>1', '>10', '>100'),
       col=ifelse(labs3<0.5, 'black', 'red'))

### now try it on the entire data set
tunq=fish2 %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()==1) %>% mutate(num=n())
tdup=fish2 %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()>1) %>% mutate(num=n())
tdupmn=tdup %>% mutate_if(is.numeric, median) %>% mutate_if(is.character, funs(paste(unique(.), collapse = "_"))) %>% slice(1)
tfin=bind_rows(tunq, tdupmn)

tunq=t %>% group_by(LAT, LON, MONTH, YEAR) %>% filter(n()>1) %>% mutate(num=n()) %>% mutate_if(is.numeric, mean) %>%
  mutate_if(is.character, funs(paste(unique(.), collapse = "_"))) %>% slice(1)


map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70", mar=c(4,0,0,0))
points(tdup$LON, tdup$LAT, pch=1, cex=as.numeric(as.character(tdup$bin)), col=ifelse(as.numeric(as.character(tdup$bin))<0.5, 'black', 'red'))
legend(-70, 38, pch=1, pt.cex=labs3, legend=c('0', '1', '>1', '>10', '>100'),
       col=ifelse(labs3<0.5, 'black', 'red'))

tunq=t %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()==1) %>% mutate(num=n())
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70", mar=c(4,0,0,0))
points(tunq$LON, tunq$LAT, pch=1, cex=as.numeric(as.character(tunq$bin)), col=ifelse(as.numeric(as.character(tunq$bin))<0.5, 'black', 'red'))
legend(-70, 38, pch=1, pt.cex=labs3, legend=c('0', '1', '>1', '>10', '>100'),
       col=ifelse(labs3<0.5, 'black', 'red'))





