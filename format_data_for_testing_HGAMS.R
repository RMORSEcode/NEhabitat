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

### REMOVE replicate zooplankton data - added to formatting for HGAMS script 2/23/21
tunq=fish2 %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()==1) %>% mutate(num=n())
tdup=fish2 %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()>1) %>% mutate(num=n())
tdupmn=tdup %>% mutate_if(is.numeric, mean) %>% mutate_if(is.character, funs(paste(unique(.), collapse = "_"))) %>% slice(1)
fish2=bind_rows(tunq, tdupmn) %>% ungroup()


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



### plot map of sample magnitude
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70", mar=c(4,0,0,0))
map.axes(las=1)
points(tt2$LON, tt2$LAT)
# br=c(-Inf, 0,10,50,100,1000)
br=c(-1,0,1,10,100,1000) # set breaks for ichtyo
# labs=c('0', '1', '>1', '>10', '>100')
labs2=c(-Inf,0,1,2,3) # log10 breaks
# labs3=c(0.5, 1, 1.25, 1.5, 2) # plotting size breaks
labs3=c(0.25, 0.5, 1, 1.5, 2) # plotting size breaks for open circle routine below:

yy=2004 # select year to plot
t=fish2 %>% filter(YEAR==yy)
t$bin=cut(t$`74_ich`, breaks=br, labels=labs3, right=T)
map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70", mar=c(4,0,0,0))
points(t$LON, t$LAT, pch=1, cex=as.numeric(as.character(t$bin)), col=ifelse(as.numeric(as.character(t$bin))<0.5, 'black', 'red'))
legend(-70, 38, pch=1, pt.cex=labs3, legend=c('0', '1', '>1', '>10', '>100'), col=ifelse(labs3<0.5, 'black', 'red'))
legend(-76,45, legend=yy, bty='n')
# points(t$LON, t$LAT, pch=19, cex=as.numeric(t$bin)/2, col=ifelse(as.numeric(t$bin)/2<1, 'black', 'red'))
# legend(-70, 38, pch=19, cex=unique(as.numeric(t2))/2, legend=c('0', '< 10', '< 50'))

# checking on frequency years with higher positive tows in the MAB
t=fish2 %>% filter(`74_ich`>0 & LAT<40.5)
table(t$YEAR)
barplot(table(t$YEAR))
t=fish2 %>% filter(`74_ich`>0 & LAT<40)

m=t[rev(order(t$`74_ich`)),]
table(m$YEAR[1:20])

### checking on duplicated data points in fish2... (solved, keeping for posterity) ###
# t5=fish2 %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n() ==1) # drops from 9470 to ~3600 samples!
# table(t5$YEAR)
# 
# #CHECK on single year duplicate stations
# yy=2004
# t=fish2 %>% filter(YEAR==yy)
# t$bin=cut(t$`74_ich`, breaks=br, labels=labs3, right=T)
# 
# tdup=t %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()>1) %>% mutate(num=n())
# tdupmn=tdup %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()>1) %>% summarise_at(vars(`74_Adt`:chaeto_100m3), mean, na.rm=T) # take mean of all duplicate rows
# tdupmd=t %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()>1) %>% summarise_all(funs(med)) # take mean of all duplicate rows
# 
# ### this works to average duplicates, triplicates, (etc) and select only first (after taking mean of values)
# tdupmn=t %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()>1) %>% mutate(num=n()) %>% mutate_if(is.numeric, mean) %>%
#   mutate_if(is.character, funs(paste(unique(.), collapse = "_"))) %>% slice(1) #distinct()
# 
# tunq=t %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()==1) %>% mutate(num=n())
# tdup=t %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()>1) %>% mutate(num=n())
# tdupmn=tdup %>% mutate_if(is.numeric, median) %>% mutate_if(is.character, funs(paste(unique(.), collapse = "_"))) %>% slice(1)
# tfin=bind_rows(tunq, tdupmn)
# map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70", mar=c(4,0,0,0))
# points(tfin$LON, tfin$LAT, pch=1, cex=as.numeric(as.character(tfin$bin)), col=ifelse(as.numeric(as.character(tfin$bin))<0.5, 'black', 'red'))
# legend(-70, 38, pch=1, pt.cex=labs3, legend=c('0', '1', '>1', '>10', '>100'),
#        col=ifelse(labs3<0.5, 'black', 'red'))
# legend(-77,45, legend=yy)
# 
# ### now try it on the entire data set - works, added to formatting for HGAMS script 2/23/21
# tunq=fish2 %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()==1) %>% mutate(num=n())
# tdup=fish2 %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()>1) %>% mutate(num=n())
# tdupmn=tdup %>% mutate_if(is.numeric, median) %>% mutate_if(is.character, funs(paste(unique(.), collapse = "_"))) %>% slice(1)
# tfin=bind_rows(tunq, tdupmn)
# 
# tunq=t %>% group_by(LAT, LON, MONTH, YEAR) %>% filter(n()>1) %>% mutate(num=n()) %>% mutate_if(is.numeric, mean) %>%
#   mutate_if(is.character, funs(paste(unique(.), collapse = "_"))) %>% slice(1)
# 
## map unique and duplicated values from yearly 't' set above
# tdup=t %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()>1) %>% mutate(num=n())
# map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70", mar=c(4,0,0,0))
# points(tdup$LON, tdup$LAT, pch=1, cex=as.numeric(as.character(tdup$bin)), col=ifelse(as.numeric(as.character(tdup$bin))<0.5, 'black', 'red'))
# legend(-70, 38, pch=1, pt.cex=labs3, legend=c('0', '1', '>1', '>10', '>100'),
#        col=ifelse(labs3<0.5, 'black', 'red'))
# tunq=t %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()==1) %>% mutate(num=n())
# map("worldHires", xlim=c(-77,-65),ylim=c(35,45), fill=T,border=0,col="gray70", mar=c(4,0,0,0))
# points(tunq$LON, tunq$LAT, pch=1, cex=as.numeric(as.character(tunq$bin)), col=ifelse(as.numeric(as.character(tunq$bin))<0.5, 'black', 'red'))
# legend(-70, 38, pch=1, pt.cex=labs3, legend=c('0', '1', '>1', '>10', '>100'),
#        col=ifelse(labs3<0.5, 'black', 'red'))


#### copied from HGAMs_format_data.R 20210421 - used to try and merge zoo, ich, and survdat, ended up not using 
# and just following the fuzzy join method by individual year (faster, no crashes). Code saved here just in case

### Load the data and rename ichthyoplankton
load('/home/ryan/Git/NEhabitat/Final_merged_fish_corBIO_Zoo_Ich.Rda') # Biomass
load('/home/ryan/Git/NEhabitat/Final_merged_fish_corABN_Zoo_Ich.Rda') # Abundance
## newest with updated survdat values
load('/home/ryan/Git/NEhabitat/20210326_unique_wide_format_corABN_by_stage.Rda')# Abundance
load('/home/ryan/Git/NEhabitat/20210326_unique_wide_format_corBIO_by_stage.Rda')# Biomass
test=FData.abn %>% left_join(svdwide.abn, by = c("CRUISE6", "STATION", "STRATUM", "SVVESSEL", 
                                                 "YEAR", "EST_TOWDATE", "SEASON", "LAT", "LON"))

og=FData.abn %>% filter(YEAR>1976)
og$MONTH=month(og$EST_TOWDATE)
tunq=og %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()==1) %>% mutate(num=n())
tdup=og %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()>1) %>% mutate(num=n())
tdupmn=tdup %>% mutate_if(is.numeric, mean) %>% mutate_if(is.character, funs(paste(unique(.), collapse = "_"))) %>% slice(1)
og2=bind_rows(tunq, tdupmn) %>% ungroup() #reassign name
og2[with(og2, order("EST_TOWDATE", "YEAR", "MONTH", "STATION", "STRATUM", "LAT", "LON")),]
og2$EST_TOWDATE=as.character(og2$EST_TOWDATE)

new=svdwide.abn %>% filter(YEAR>1976)
new$MONTH=month(new$EST_TOWDATE)
tunq=new %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()==1) %>% mutate(num=n())
tdup=new %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()>1) %>% mutate(num=n())
tdupmn=tdup %>% mutate_if(is.numeric, mean) %>% mutate_if(is.character, funs(paste(unique(.), collapse = "_"))) %>% slice(1)
new2=bind_rows(tunq, tdupmn) %>% ungroup() #reassign name
new2[with(new2, order("EST_TOWDATE", "YEAR", "MONTH", "STATION", "STRATUM", "LAT", "LON")),]
new2$EST_TOWDATE=as.character(new2$EST_TOWDATE)
# keep original, merge in new values
test=og2 %>% left_join(new2, by = c("CRUISE6", "STATION", "STRATUM", "SVVESSEL", 
                                    "YEAR", "EST_TOWDATE", "LAT", "LON"))
## keep new values from survdat, merge ich and zoo
test=new2 %>% left_join(og3, by =c("CRUISE6", "STATION", "STRATUM", "SVVESSEL", "YEAR",  
                                   "LAT", "LON", "DEPTH")) #, "SURFTEMP", "BOTTEMP", "SURFSALIN", 
# "BOTSALIN", "72_Juv", "72_Adt", "73_Juv", "73_Adt", "74_Juv", "74_Adt", 
# "75_Adt", "75_Juv", "76_Adt", "76_Juv", "77_Juv", "77_Adt", "101_Juv", 
# "101_Adt", "102_Juv", "102_Adt", "103_Adt", "103_Juv", "105_Juv", 
# "105_Adt", "106_Adt", "106_Juv", "107_Juv", "107_Adt", "108_Adt", 
# "108_Juv", "155_Juv", "155_Adt", "192_Juv", "192_Adt", "193_Adt", 
# "193_Juv", "197_Juv", "197_Adt"))
test=new2 %>% left_join(og2, by =c("CRUISE6", "STATION", "STRATUM", "SVVESSEL", "YEAR",  
                                   "LAT", "LON", "DEPTH", "SURFTEMP", "BOTTEMP", "SURFSALIN", 
                                   "BOTSALIN", "72_Juv", "72_Adt", "73_Juv", "73_Adt", "74_Juv", "74_Adt",
                                   "75_Adt", "75_Juv", "76_Adt", "76_Juv", "77_Juv", "77_Adt", "101_Juv",
                                   "101_Adt", "102_Juv", "102_Adt", "103_Adt", "103_Juv", "105_Juv",
                                   "105_Adt", "106_Adt", "106_Juv", "107_Juv", "107_Adt", "108_Adt",
                                   "108_Juv", "155_Juv", "155_Adt", "192_Juv", "192_Adt", "193_Adt",
                                   "193_Juv", "197_Juv", "197_Adt"))# %>% filter(complete.cases(`74_ich`))

## something is not right, many NAs where there should be matchups
og2[1,c("CRUISE6", "STATION", "STRATUM", "SVVESSEL", "YEAR", "EST_TOWDATE", 
        "SEASON", "LAT", "LON", "DEPTH")]
new2[1,c("CRUISE6", "STATION", "STRATUM", "SVVESSEL", "YEAR", "EST_TOWDATE", 
         "SEASON", "LAT", "LON", "DEPTH")]

new3=new2[1:15,]
og4=og3[1:15,]
og4$STATION[which(!(og4$STATION %in% new3$STATION))]
# new3$EST_TOWDATE=as.character(new3$EST_TOWDATE)
# og3$EST_TOWDATE=as.character(og3$EST_TOWDATE)
# og3[1,c("CRUISE6", "STATION", "STRATUM", "SVVESSEL", "YEAR", "EST_TOWDATE", "SEASON", "LAT", "LON", "DEPTH")]==new3[1,c("CRUISE6", "STATION", "STRATUM", "SVVESSEL", "YEAR", "EST_TOWDATE", "SEASON", "LAT", "LON", "DEPTH")]

cog2=lapply(og2, class)
cnew2=lapply(new2, class)

## clean up
og3=og2 %>% select(-mrgidx, date)
og3$YEAR=as.integer(og3$YEAR)
og3$SEASON=as.character(og3$SEASON)
og3$SVVESSEL=as.character(og3$SVVESSEL)
og3$EST_TOWDATE=as.character(og3$EST_TOWDATE)

new3$EST_TOWDATE=as.character(new3$EST_TOWDATE)


### 20210406 merge original FData.abn with new survdat
test=FData.abn %>% left_join(svdwide.abn, by =c("CRUISE6", "STATION", 
                                                "STRATUM", "SVVESSEL", "YEAR","LAT", "LON", "DEPTH"))
# tunq=test %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()==1) %>% mutate(num=n())
# tdup=test %>% group_by(LAT, LON, MONTH, YEAR, SURFTEMP) %>% filter(n()>1) %>% mutate(num=n())

tunq=test %>% group_by(CRUISE6, STATION, STRATUM, SVVESSEL, YEAR,LAT, LON, DEPTH) %>% filter(n()==1) %>% mutate(num=n())
tdup=test %>% group_by(CRUISE6, STATION, STRATUM, SVVESSEL, YEAR,LAT, LON, DEPTH) %>% filter(n()>1) %>% mutate(num=n())

tdupmn=tdup %>% mutate_if(is.numeric, mean) %>% mutate_if(is.character, funs(paste(unique(.), collapse = "_"))) %>% slice(1)
test2=bind_rows(tunq, tdupmn) %>% ungroup() #reassign name
test2[with(test2, order("EST_TOWDATE", "YEAR", "MONTH", "STATION", "STRATUM", "LAT", "LON")),]
test2$EST_TOWDATE=as.character(test2$EST_TOWDATE)

test4=test2[complete.cases(test2$`74_Adt.y`),] # 15619 long, but NAs for zoo taxa...
test5=test4[complete.cases(test4$ctyp_100m3),] # 8183 long...not good
test6=test2[complete.cases(test2$ctyp_100m3),] # 14219
test7=test6[1:500,]
test8=test7 %>% select(EST_TOWDATE.x, EST_TOWDATE.y, STATION, STRATUM, LAT, LON, `74_Juv.x`, `74_Juv.y`, `74_Adt.x`, `74_Adt.y`, `74_ich`)
# why would some be different? NA values for stations
CoreOffshoreStrata<-c(seq(1010,1300,10),1340, seq(1360,1400,10),seq(1610,1760,10))
# inshore strata to use, still sampled by Bigelow
CoreInshore73to12=c(3020, 3050, 3080 ,3110 ,3140 ,3170, 3200, 3230, 3260, 3290, 3320, 3350 ,3380, 3410 ,3440)
# combine
strata_used=c(CoreOffshoreStrata,CoreInshore73to12)
sum(unique(test8$STRATUM[which(is.na(test8$`74_Juv.y`))]) %in% strata_used)

### check on station numbers, lat lon, est_towdate
load("/home/ryan/Downloads/survdat_lw.RData")# loads survdat.lw ; from Sean Lucey data into 2019
nsdat2=nsdat %>% filter(SVSPP==74, YEAR>1976)
nsdat3=survdat %>% filter(SVSPP==74, YEAR>1976)
sdat=survdat.lw %>% filter(SVSPP==74, YEAR>1976)
colnames(nsdat3)
colnanes(sdat)
nsdat3=ungroup(nsdat3)
class(nsdat3)
class(sdat)
sdat=as_tibble(sdat)
test=sdat%>% left_join(nsdat3, by =c("CRUISE6", "STATION", "STRATUM", "SVVESSEL", "YEAR","LAT", "LON", "DEPTH"))
test2=test[1:100,]

## clean up

cnsdat3=lapply(nsdat3, class)
csdat=lapply(sdat, class)

sdat$YEAR=as.integer(sdat$YEAR)
sdat$SEASON=as.character(sdat$SEASON)
sdat$SVVESSEL=as.character(sdat$SVVESSEL)
sdat$EST_TOWDATE=as.character(sdat$EST_TOWDATE)

# nsdat3[,c('YEAR', 'SVSPP', 'TOW', 'STRATUM', 'STATION', 'CRUISE6', 'CATCHSEX', 'NUMLEN', 'DEPTH')]=as.integer(nsdat3[,c('YEAR', 'SVSPP', 'TOW', 'STRATUM', 'STATION', 'CRUISE6', 'CATCHSEX', 'NUMLEN', 'DEPTH')])
# nsdat4=lapply(nsdat3[,c('YEAR', 'SVSPP', 'TOW', 'STRATUM', 'STATION', 'CRUISE6', 'CATCHSEX', 'NUMLEN', 'DEPTH')], as.integer)

nsdat3$YEAR=as.integer(nsdat3$YEAR)
nsdat3$SVSPP=as.integer(nsdat3$SVSPP)
nsdat3$TOW=as.integer(nsdat3$TOW)
nsdat3$STRATUM=as.integer(nsdat3$STRATUM)
nsdat3$STATION=as.integer(nsdat3$STATION)
nsdat3$CRUISE6=as.integer(nsdat3$CRUISE6)
nsdat3$CATCHSEX=as.integer(nsdat3$CATCHSEX)
nsdat3$NUMLEN=as.integer(nsdat3$NUMLEN)
nsdat3$DEPTH=as.integer(nsdat3$DEPTH)
nsdat3$SEASON=as.character(nsdat3$SEASON)
nsdat3$SVVESSEL=as.character(nsdat3$SVVESSEL)
nsdat3$EST_TOWDATE=as.character(nsdat3$EST_TOWDATE)
# nsdat3$SIZECAT=as.character(nsdat3$SIZECAT)

sdat$YEAR=as.integer(sdat$YEAR)
sdat$SVSPP=as.integer(sdat$SVSPP)
sdat$TOW=as.integer(sdat$TOW)
sdat$STRATUM=as.integer(sdat$STRATUM)
sdat$STATION=as.integer(sdat$STATION)
sdat$CRUISE6=as.integer(sdat$CRUISE6)
sdat$CATCHSEX=as.integer(sdat$CATCHSEX)
sdat$NUMLEN=as.integer(sdat$NUMLEN)
sdat$DEPTH=as.integer(sdat$DEPTH)
sdat$SEASON=as.character(sdat$SEASON)
sdat$SVVESSEL=as.character(sdat$SVVESSEL)
sdat$EST_TOWDATE=as.character(sdat$EST_TOWDATE)
sdat$SIZECAT=as.character(sdat$SIZECAT)

sdat[with(sdat, order("EST_TOWDATE", "YEAR", "MONTH", "STATION", "STRATUM", "LAT", "LON")),]
nsdat3[with(nsdat3, order("EST_TOWDATE", "YEAR", "MONTH", "STATION", "STRATUM", "LAT", "LON")),]


### now create dataframe with unique tows only, separated by stage
test=nsdat3[,c("CRUISE6","STATION","STRATUM","YEAR", "SVSPP", "stg")] # needs SVSPP...
svdtunq=nsdat3[!duplicated(test),]

sdat=survdat %>% filter(SVSPP==74, YEAR>1976)
test=sdat[,c("CRUISE6","STATION","STRATUM","YEAR", "SVSPP", "stg")] # needs SVSPP...
ogsvdtunq=sdat[!duplicated(test),]

test=svdtunq[1:200,]
test2=ogsvdtunq[1:200,]
test3=left_join(test, test2, by =c("CRUISE6", "STATION", "STRATUM", "YEAR","LAT", "LON", "DEPTH"))

test3=left_join(ogsvdtunq, svdtunq, by =c("CRUISE6", "STATION", "STRATUM", "YEAR","LAT", "LON", "DEPTH"))
test3=inner_join(ogsvdtunq, svdtunq, by =c("CRUISE6", "STATION", "STRATUM", "YEAR","LAT", "LON", "DEPTH"))
barplot(table(test3$YEAR[which(test3$SEASON.x=="SPRING")]))




