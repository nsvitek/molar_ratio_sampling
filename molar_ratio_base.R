# !diagnostics off

# Use molar rows of Recent Peromyscus gossypinus & others as comparison to fossil record
# dependencies -----
library(dplyr) #organize
library(ggplot2) #plot
# library(gridExtra) #plot
library(ggforce) #allows drawing circles and ellipses
# library(scales) #allows use of standard colors. May not be completely necessary
# library(ggthemes) #to get Paul Tol colors
library(readxl) #read raw data
library(reshape2) #for reformatting & comparing CV, measurement error data
library(lmodel2) #for RMA regression
library(MCMCglmm) #for Bayesian modelling
library(stringr) #using str_detect for sample size resampling in non-mice

# locateScripts<-"C:/cygwin/home/N.S/scripts/molar_ratio_sampling"
locateScripts<-"C:/scripts/molar_ratio_sampling"
# locateData<-"C:/Users/N.S/Dropbox/Documents/research/vitek-etal_inhibitory-cascade-isolated"
locateData<-"D:/Dropbox/Documents/research/vitek-etal_inhibitory-cascade-isolated"

setwd(locateScripts)
# source("../scripts/function_bootstrap.R") #may not need this anymore. 
#script-specific functions:
source("molar_ratio_functions.R")

# settings ------
replicates<-999 #bootstrap number, low for now, increase for final calcs
# n<-10 #resample number

# load data -----------
setwd(locateData)
#read in raw data
mouse.raw<-read_excel("input/peromyscus_gossypinus_icm.xlsx")
#read in comparative data. ICM comparative variance is mostly from various studies involving Asahara
compICM.raw<-read_excel("input/comparison_ICM_asahara.xlsx")
#and also Roseman and Delezene's dataset for comparison
compICM2.raw<-read_excel("input/comparison_ICM_roseman.xlsx", sheet = "Molar Size Data", na="NA")
#and also Labonne et al's dataset for comparison
compICM3.raw<-read_excel("input/comparison_ICM_labonne.xlsx")

#comparative MMC data from Monson et al. 2019
compMMC.raw<-read_excel("input/comparison_MMC_monson.xlsx")
#read in the data of macrocranion
fossil.raw<-read_excel("input/measurements_macrocranion.xlsx")
#read in data pertaining to CV of raw lengths and areas
# comp.areaCV.raw<-read_excel("inhibitory_results/input/comparison_CA_CV.xlsx")
# comp.lengthCV.raw

# format cotton mouse --------
mouse.raw$specimen<-paste(mouse.raw$institution,mouse.raw$catalog_number,sep="-")
mouse<-dplyr::select(mouse.raw,-c(institution,catalog_number,genus,species,filestart,level_fromfloor,
                                  voxels,locality))

mouse$Order<-"Rodentia"
mouse$Species<-"Peromyscus gossypinus"
#take means of raw measurement replicates
mouse$m1.length<-apply(cbind(mouse$m1.l1,mouse$m1.l2,mouse$m1.l3),1,mean)
mouse$m1.width<-apply(cbind(mouse$m1.w1,mouse$m1.w2,mouse$m1.w3),1,mean)
mouse$m2.length<-apply(cbind(mouse$m2.l1,mouse$m2.l2,mouse$m2.l3),1,mean)
mouse$m2.width<-apply(cbind(mouse$m2.w1,mouse$m2.w2,mouse$m2.w3),1,mean)
mouse$m3.length<-apply(cbind(mouse$m3.l1,mouse$m3.l2,mouse$m3.l3),1,mean)
mouse$m3.width<-apply(cbind(mouse$m3.w1,mouse$m3.w2,mouse$m3.w3),1,mean)

mouse$total.length<-mouse$m1.length+mouse$m2.length+mouse$m3.length
#clean data
mouse<-dplyr::select(mouse,specimen,state,county,sex,m1.length,m1.width,m2.length,m2.width,
                     m3.length,m3.width,total.length)
#calculate crown area:
mouse$m1.area<-mouse$m1.length*mouse$m1.width
mouse$m2.area<-mouse$m2.length*mouse$m2.width
mouse$m3.area<-mouse$m3.length*mouse$m3.width
mouse$total.area<-mouse$m1.area+mouse$m2.area+mouse$m3.area

#calculate proportional size of molars relative to total molar size
mouse$prop.m1A<-mouse$m1.area/mouse$total.area
mouse$prop.m2A<-mouse$m2.area/mouse$total.area
mouse$prop.m3A<-mouse$m3.area/mouse$total.area
mouse$prop.m1L<-mouse$m1.length/mouse$total.length
mouse$prop.m2L<-mouse$m2.length/mouse$total.length
mouse$prop.m3L<-mouse$m3.length/mouse$total.length

#calculate mouse "a.i" ratio (if it follows ICM):
mouse$a.iA<-mouse$m2.area/mouse$m1.area
mouse$a.iL<-mouse$m2.length/mouse$m1.length

#calculate ICM ratios:
mouse$m2.m1A<-mouse$m2.area/mouse$m1.area
mouse$m3.m1A<-mouse$m3.area/mouse$m1.area

#calculate MMC ratios:
mouse$m2.m1L<-mouse$m2.length/mouse$m1.length
mouse$m3.m1L<-mouse$m3.length/mouse$m1.length

#test assumption of normality
#In Kavanagh et al. 2007, the area metrics are not transformed
#Does that protocol seems okay in terms of normality, according to tests?
shapiro.test(mouse$m3.m1A)
shapiro.test(mouse$m2.m1A) #yes

#Do MMC ratios also meet assumptions of normality?
shapiro.test(mouse$m3.m1L)
shapiro.test(mouse$m2.m1L) #yes

# format ICM comparative data ------
#get ratios, other info for Roseman & Delezene data
compICM2.raw$m3.m1A<-compICM2.raw$M3.Area/compICM2.raw$M1.Area
compICM2.raw$m2.m1A<-compICM2.raw$M2.Area/compICM2.raw$M1.Area
compICM2<-compICM2.raw[!is.na(compICM2.raw$m3.m1A),]
compICM2<-compICM2[!is.na(compICM2$m2.m1A),]
compICM2$Order<-"Primates"
compICM2$Reference<-"Roseman and Delezene 2019"

#format Labonne data
colnames(compICM3.raw)<-c("family","subfamily","Species","p4","M1.Area","M2.Area","M3.Area",
                          "composite","m2.m1A","m3.m1A","fossilextant","ref")
compICM3<-compICM3.raw %>% filter(.,composite=="no",fossilextant=="E") %>% 
  filter(.,grepl(" sp\\.",Species,perl=TRUE)==FALSE) %>% 
  select(.,Species,M1.Area,M2.Area,M3.Area,m3.m1A,m2.m1A)
compICM3.spp.counts<-compICM3 %>% group_by(Species) %>% summarize (N=n()) 
exclude.me<-compICM3.spp.counts$Species[which(compICM3.spp.counts$N==1)]
compICM3<-filter(compICM3,Species %in% exclude.me == FALSE)
compICM3$Order<-"Rodentia"
compICM3$Sex<-"unknown"
compICM3$Reference<-"Labonne et al. 2012"

compICM.pub<-rbind(compICM2,compICM3)
#how many species remaining?
compICM.pub %>% group_by(Species) %>% summarize (N=n()) %>% nrow
nrow(compICM.pub)

#add in peromyscus gossypinus data
compICM2add<-dplyr::select(mouse,m2.m1A,m3.m1A)
compICM2add$Species<-"Peromyscus gossypinus area"
compICM2add$Order<-"Rodentia"
compICM2add$Sex<-"unknown"
compICM2add$M1.Area<-mouse$m1.area
compICM2add$M2.Area<-mouse$m2.area
compICM2add$M3.Area<-mouse$m3.area
compICM2add$Reference<-"this study"
compICM.indv<-rbind(compICM.pub,compICM2add)

#summarize individual-level data so it can be collated with population-level data
ICM.CV<- compICM.indv %>% group_by(Order, Species, Sex) %>%  
  summarize(N=n(), m3m1.mean=mean(m3.m1A), m3m1.SD=sd(m3.m1A),m3m1.CV=sd(m3.m1A)/mean(m3.m1A)*100,
            m2m1.mean=mean(m2.m1A), m2m1.SD=sd(m2.m1A),m2m1.CV=sd(m2.m1A)/mean(m2.m1A)*100,Reference=unique(Reference))
colnames(ICM.CV)[which(colnames(ICM.CV)=="Sex")]<-"condition"

#organize collated population-level ICM ratio data into standardized set of variables
compICM.raw$N<-as.integer(compICM.raw$N)
compICM.raw$m2m1.CV<-(compICM.raw$m2m1.SD/compICM.raw$m2m1.mean) *100
compICM.raw$m3m1.CV<-(compICM.raw$m3m1.SD/compICM.raw$m3m1.mean) *100
compICM.pop<-compICM.raw %>% bind_rows(.,ICM.CV)

compICM.pop %>% group_by(Species) %>% summarize (N=n()) %>% nrow -1 #-1 to not double-count new Peromyscus measures
write.csv(compICM.pop,"output/SI_TableX_AreasByPop.csv")
# format MMC comparative data -----
#clean out singleton species
MMC.spp.counts<-compMMC.raw %>% group_by(Species) %>% summarize (N=n()) 
exclude.me<-MMC.spp.counts$Species[which(MMC.spp.counts$N==1)]
compMMC.raw<-dplyr::filter(compMMC.raw,Species %in% exclude.me == FALSE)
#also calc m2:m1 ratio
compMMC.raw$MMC2<-compMMC.raw$`DM2L (mm)`/compMMC.raw$`DM1L (mm)`
#how many species remaining?
compMMC.raw %>% group_by(Species) %>% summarize (N=n()) %>% nrow

#MMC: add Peromyscus gossypinus to comparative dataset
mouse.MMC.2add<-dplyr::select(mouse,m3.m1L)
mouse.MMC.2add$Order<-"Rodentia"
mouse.MMC.2add$Family<-"Cricetidae"
mouse.MMC.2add$Species<-"Peromyscus gossypinus length"
mouse.MMC.2add$`Museum ID`<-mouse.MMC.2add$`Museum Repository`<-NA
mouse.MMC.2add$`DM3L (mm)`<-mouse$m3.length
mouse.MMC.2add$`DM2L (mm)`<-mouse$m2.length
mouse.MMC.2add$`DM1L (mm)`<-mouse$m1.length
mouse.MMC.2add$`DP4L (mm)`<-NA
mouse.MMC.2add$MMC<-mouse$m3.m1L
mouse.MMC.2add$MMC2<-mouse$m2.m1L
mouse.MMC.2add$PMM<-mouse.MMC.2add$`Measured by`<-NA

mouse.MMC.2add<-dplyr::select(mouse.MMC.2add,-m3.m1L)

#no need to write raw measures because these are only Monson SI + Peromyscus, 
#which are/will be published already in separate tabels
compMMC<-rbind(compMMC.raw,mouse.MMC.2add)

# format fossils -----
fossil.raw$length<-apply(cbind(fossil.raw$length1,fossil.raw$length2,fossil.raw$length3),1,
                         function(x) mean(as.numeric(x)))
fossil.raw$width<-apply(cbind(fossil.raw$talonid_width1,fossil.raw$talonid_width2,
                              fossil.raw$talonid_width3),1,function(x) mean(as.numeric(x)))
fossil.raw$area<-fossil.raw$length * fossil.raw$width
# calculate mouse measurement error ------
source(paste(locateScripts,"molar_ratio_measurement_error.R",sep="/"))
#output that you care about is avg.repeat

# ICM expectations: length or area? ----
source(paste(locateScripts,"molar_ratio_ICM_expectations.R",sep="/"))

# standing variation in species -----
#MMC: calculate CV by species from Monson, gossypinus data
MMC.CV<- compMMC %>% group_by(Order, Species) %>% 
  summarize(N=n(), MMC.mean=mean(MMC), MMC.SD=sd(MMC), MMC.CV=sd(MMC)/mean(MMC)*100,
            MMC2.mean=mean(MMC2), MMC2.SD=sd(MMC2), MMC2.CV=sd(MMC2)/mean(MMC2)*100)
MMC.CV<-MMC.CV[complete.cases(MMC.CV),]
write.csv(MMC.CV,"output/SI_TableX_LengthsByPop.csv")

#create object for plotting
MMC.CV$condition<-"unknown"
CV.survey<-bind_rows(select(compICM.pop,Order,Species,condition,m2m1.CV,m3m1.CV) %>% melt(id=c("Order","Species","condition")),
                     select(MMC.CV,Order,Species,condition,MMC.CV,MMC2.CV) %>% melt(id=c("Order","Species","condition")))
#remember to reset to numeric
CV.survey$value<-as.numeric(CV.survey$value)

#organize by which has a tooth and what kind of measurement
CV.survey$tooth<-"M3"
CV.survey$tooth[grepl("2",CV.survey$variable,perl=TRUE)]<-"M2"
CV.survey$dimension<-"area"
CV.survey$dimension[grepl("MMC",CV.survey$variable,perl=TRUE)]<-"length"

CV.survey2<-CV.survey %>% select(-variable) %>% dcast(., Order + Species + condition + dimension ~ tooth)

#summarize, report
MMC<-ungroup(MMC.CV) %>% summarize(mean=mean(MMC.CV),SD=sd(MMC.CV),median=median(MMC.CV),min=min(MMC.CV),max=max(MMC.CV))
MMC2<- ungroup(MMC.CV) %>% summarize(mean=mean(MMC2.CV),SD=sd(MMC2.CV),median=median(MMC2.CV),min=min(MMC2.CV),max=max(MMC2.CV))
M2.M1<-compICM.pop %>% summarize(mean=mean(m2m1.CV),SD=sd(m2m1.CV),median=median(m2m1.CV),min=min(m2m1.CV),max=max(m2m1.CV))
M2.M1.noRaccoonDog<-compICM.pop %>% dplyr::filter(.,Species!="Nyctereutes procyonoides") %>% 
  summarize(mean=mean(m2m1.CV),SD=sd(m2m1.CV),median=median(m2m1.CV),min=min(m2m1.CV),max=max(m2m1.CV))
M3.M1<-compICM.pop %>% summarize(mean=mean(m3m1.CV),SD=sd(m3m1.CV),median=median(m3m1.CV),min=min(m3m1.CV),max=max(m3m1.CV))
M3.M1.noRaccoonDog<-compICM.pop %>% dplyr::filter(.,Species!="Nyctereutes procyonoides") %>% 
  summarize(mean=mean(m3m1.CV),SD=sd(m3m1.CV),median=median(m3m1.CV),min=min(m3m1.CV),max=max(m3m1.CV))

CV.summary<-rbind(MMC,MMC2,M2.M1, M3.M1,M2.M1.noRaccoonDog,M3.M1.noRaccoonDog) %>% as.data.frame
CV.summary$CI.upper<-CV.summary$mean+(CV.summary$SD*1.96)
CV.summary$CI.lower<-CV.summary$mean-(CV.summary$SD*1.96)
row.names(CV.summary)<-c("M3.M1L","M2.M1L","M2.M1A", "M3.M1A","M2.M1A.noRaccoonDog","M3.M1A.noRaccoonDog")
write.csv(CV.summary,"output/CV_summary_stats.csv")

# geographic variation: u-tests -----
#gossypinus: compare the 4 state-level populations using U test and bootstrap.
#define testing variables
group.var<-mouse$state  #which variable are you using as your group identity
pop.sizes<-mouse %>% group_by(state) %>% summarize(n()) #calculate sample size for each group

#pairwise U-test with Bonferroni correction approach following Asahara 2014
U.MMC<-Pairwise.U(mouse$m3.m1L,group.var)
U.ICM31<-Pairwise.U(mouse$m3.m1A,group.var)
U.ICM21<-Pairwise.U(mouse$m2.m1A,group.var)

cbind(U.MMC, U.ICM21[,5:8], U.ICM31[,5:8]) %>% 
  write.csv(.,"output/gossypinus_population_Utest_M23.csv")

#Alternate approach, but provides similar results.
# #bootstrap ANOVA approach
# resam.distribution<-BootstrapANOVA(mouse$m3.m1L,group.var,pop.sizes,replicates)
# bootstrap.p<-which(resam.distribution>=resam.distribution[1]) %>% length()/replicates
# 
# resam.distribution.ICM1<-BootstrapANOVA(mouse$m3.m1A,group.var,pop.sizes,replicates)
# bootstrap.p.ICM1<-which(resam.distribution.ICM1>=resam.distribution.ICM1[1]) %>% length()/replicates
# 
# resam.distribution.ICM2<-BootstrapANOVA(mouse$m2.m1A,group.var,pop.sizes,replicates)
# bootstrap.p.ICM2<-which(resam.distribution.ICM2>=resam.distribution.ICM2[1]) %>% length()/replicates
# cbind(bootstrap.p,bootstrap.p.ICM1,bootstrap.p.ICM2) %>% 
#   write.csv(.,"output/gossypinus_population_bootstrapANOVA.csv")

# #visualize
# hist(resam.distribution.ICM2) #no significant difference from random
# abline(v=resam.distribution.ICM2[1],col="red")
# aov(lm(measurement.ICM1~group.var)) %>% summary

# master table of mouse stats ------
#make a table matching names to ratio
ratio<-matrix(c("m3.m1L","MMC.mean","MMC.SD",
                "m2.m1L","MMC2.mean","MMC2.SD",
                "m3.m1A","m3m1.mean","m3m1.SD",
                "m2.m1A","m2m1.mean","m2m1.SD"),byrow = TRUE,ncol=3)

#make table of summary stats for total mouse population, mouse.stats
mouse.stats<-data.frame(Mu=rep(NA,nrow(ratio)),row.names=ratio[,1])
mouse.stats$Mu<-mouse[,ratio[,1]] %>% apply(.,2,mean)
mouse.stats$Sigma<-mouse[,ratio[,1]] %>% apply(.,2,sd)
mouse.stats$CV<-mouse.stats$Sigma / mouse.stats$Mu * 100
#confidence intervals, assuming normality, from Sokal and Rohlf
mouse.stats$M.CI.U<-mouse.stats$Mu + mouse.stats$Sigma * 1.96
mouse.stats$M.CI.L<-mouse.stats$Mu - mouse.stats$Sigma * 1.96
mouse.stats$S.CI.L<-(((nrow(mouse)-1) * (mouse.stats$Sigma)^2) /
                            qchisq(c(0.975),df=nrow(mouse)-1)) %>% sqrt
mouse.stats$S.CI.U<-(((nrow(mouse)-1) * (mouse.stats$Sigma)^2) /
                            qchisq(c(0.025),df=nrow(mouse)-1)) %>% sqrt
# geographic variation: signs test ----------- 
#make a table of CVs by state
ratio.pop.CV<- mouse %>% group_by(state) %>% 
  summarize(N=n(),CV.MMC=sd(m3.m1L)/mean(m3.m1L)*100,mean.MMC=mean(m3.m1L),
            CV.m3m1=sd(m3.m1A)/mean(m3.m1A)*100,mean.m3m1=mean(m3.m1A),
            CV.m2m1=sd(m2.m1A)/mean(m2.m1A)*100,mean.m2m1=mean(m2.m1A)) 
#add total
ratio.pop.CV.2<-rbind(select(ratio.pop.CV,state,CV.MMC,CV.m3m1,CV.m2m1),
                      c("Total",mouse.stats$CV[1],mouse.stats$CV[2],mouse.stats$CV[3])) %>%
  melt(id="state")
#take back out of character
ratio.pop.CV.2$value<-as.numeric(ratio.pop.CV.2$value)

#does pooling localities result in a higher CV?
mice.high.CV<-(which(ratio.pop.CV$CV.m3m1 < mouse.stats$CV[which(ratio[,1]=="m3.m1A")]) %>% length) + 
  (which(ratio.pop.CV$CV.m2m1 < mouse.stats$CV[which(ratio[,1]=="m2.m1A")]) %>% length) +
  (which(ratio.pop.CV$CV.MMC < mouse.stats$CV[which(ratio[,1]=="m2.m1L")]) %>% length)
#here, "successes" are where pooled CV is higher than locality cV
binom.test(mice.high.CV,12,alternative = "greater") #no

# #check to see if size CV increases with pooled sample (which it could be based on plotted data)
# mouse.total.size.CV<-c(sd(mouse$m3.area)/mean(mouse$m3.area)*100,
#                        sd(mouse$m3.length)/mean(mouse$m3.length)*100)
# mouse.pop.size.CV<-mouse %>% group_by(state) %>% 
#   summarize(CV.area=sd(m3.area)/mean(m3.area)*100,CV.length=sd(m3.length)/mean(m3.length)*100)
# which(mouse.pop.size.CV$CV.area < mouse.total.size.CV[1]) %>% length
# which(mouse.pop.size.CV$CV.length < mouse.total.size.CV[2]) %>% length
# binom.test(6,8,alternative = "greater")
# sexual dimorphism  -------
ICM.dimorph<-compICM.indv %>% filter(.,compICM.indv$Sex !="unknown")
ICM.di.spp.CV<-ICM.dimorph %>% group_by(Species) %>% 
  summarize(N=n(),
            CV.m3m1=sd(m3.m1A)/mean(m3.m1A)*100,mean.m3m1=mean(m3.m1A),
            CV.m2m1=sd(m2.m1A)/mean(m2.m1A)*100,mean.m2m1=mean(m2.m1A),
            CV.m3.area=sd(M3.Area)/mean(M3.Area)*100,mean.m3.area=mean(M3.Area))

ICM.di.sex.CV<-ICM.dimorph %>% group_by(Species,Sex) %>% 
  summarize(N=n(),
            CV.m3m1=sd(m3.m1A)/mean(m3.m1A)*100,mean.m3m1=mean(m3.m1A),
            CV.m2m1=sd(m2.m1A)/mean(m2.m1A)*100,mean.m2m1=mean(m2.m1A),
            CV.m3.area=sd(M3.Area)/mean(M3.Area)*100,mean.m3.area=mean(M3.Area))

# sexual dimorphism: signs & U tests --------
#does pooling sexes result in a higher CV?
n.cv<-n.size<-m3m1.U<-m2m1.U<-m3.size.U<-NULL
n.spp<-nrow(ICM.di.spp.CV) #number of primate spp with sexual dimorphism info
for (i in 1:n.spp){
  spp.choice<-ICM.di.sex.CV %>% filter(Species %in% ICM.di.spp.CV$Species[i])
  n.cv<-c(n.cv,
          (which(spp.choice$CV.m3m1 < ICM.di.spp.CV$CV.m3m1[i]) %>% length) + 
            (which(spp.choice$CV.m2m1 < ICM.di.spp.CV$CV.m2m1[i]) %>% length))
  n.size<-c(n.size, 
            which(spp.choice$CV.m3.area < ICM.di.spp.CV$CV.m3.area[i]) %>% length)
  #also a U-test?
  sample.choice<-ICM.dimorph %>% filter(Species %in% ICM.di.spp.CV$Species[i])
  m3.size.U<-rbind(m3.size.U,
                   Pairwise.U(sample.choice$M3.Area,sample.choice$Sex))
  m3m1.U<-rbind(m3m1.U,
                Pairwise.U(sample.choice$m3.m1A,sample.choice$Sex))
  m2m1.U<-rbind(m2m1.U,
                Pairwise.U(sample.choice$m2.m1A,sample.choice$Sex))
}
row.names(m3m1.U)<-row.names(m2m1.U)<-ICM.di.spp.CV$Species
binom.test(sum(n.cv),n.spp*2*2,alternative = "greater") #no
# binom.test(sum(n.size),n.spp*2*2,alternative = "greater")

#make object for plotting
ICM.di.spp.CV$Sex<-"Total"
ratio.sex.CV<-bind_rows(ICM.di.sex.CV,ICM.di.spp.CV) %>% 
  select(Species, Sex, CV.m3m1, CV.m2m1)
ratio.sex.CV.2<-melt(ratio.sex.CV, id=c("Species","Sex"))

# sensitivity of mean and sd of ratios to sample size -----------

#which parameters to test
resampled.names<-c("N","MMC.mean","MMC.SD",#"MMC.max","MMC.min", #maxima and minima are here just in case you're curious
                   "MMC2.mean","MMC2.SD",#"MMC2.max","MMC2.min",
                   "m3m1.mean","m3m1.SD",#"m3m1.max","m3m1.min",
                   "m2m1.mean","m2m1.SD")#,"m2m1.max","m2m1.min")
species.name<-"Peromyscus gossypinus"
source(paste(locateScripts,"molar_ratio_sample_size.R",sep="/"))
# test for phylogenetic signal in CV --------------
# source(paste(locateScripts,"molar_ratio_physig.R",sep="/"))
# create simulated composite molar rows --------
#create resampled pseudoreplicates of molars sizes and ratios from composite cheek teeth
#both, simulate all ratios at once:
sim.isolate.names<-c("N","m1.area","m2.area","m3.area",
                     "m1.length","m2.length","m3.length",
                     "m2m1.SD","m3m1.SD","MMC.SD","MMC2.SD",
                     "m2m1.mean","m3m1.mean","MMC.mean","MMC2.mean")
sim.isolate.base<-matrix(NA,nrow=replicates,ncol=length(sim.isolate.names),
                         dimnames=list(c(NULL),sim.isolate.names)) %>% as.data.frame #empty holder frame
for (sample.n in 1:50){
  print(paste("starting sample size",sample.n))
  sim.isolate<-sim.isolate.base
  sim.isolate[,"N"]<-sample.n #sample size
  for (i in 1:nrow(sim.isolate.base)){ #resample each tooth position and take mean
    m1A.sample<-mouse$m1.area %>% sample(.,size=sample.n,replace=TRUE)
    m2A.sample<-mouse$m2.area %>% sample(.,size=sample.n,replace=TRUE)
    m3A.sample<-mouse$m3.area %>% sample(.,size=sample.n,replace=TRUE)
    m1L.sample<-mouse$m1.length %>% sample(.,size=sample.n,replace=TRUE)
    m2L.sample<-mouse$m2.length %>% sample(.,size=sample.n,replace=TRUE)
    m3L.sample<-mouse$m3.length %>% sample(.,size=sample.n,replace=TRUE)
    sim.isolate$m1.area[i]<-m1A.sample %>% mean
    sim.isolate$m2.area[i]<-m2A.sample %>% mean
    sim.isolate$m3.area[i]<-m3A.sample %>% mean
    sim.isolate$m1.length[i]<-m1L.sample %>% mean
    sim.isolate$m2.length[i]<-m2L.sample %>% mean
    sim.isolate$m3.length[i]<-m3L.sample %>% mean
    #make a pseudosample from the isolated molars by resampling 1 specimen each time
    pseudopop.names<-c("m2m1A","m3m1A","m3m1L","m2m1L")
    pseudopop<-matrix(NA,nrow=replicates,ncol=length(pseudopop.names),
                      dimnames=list(c(NULL),pseudopop.names)) %>% as.data.frame #empty holder frame
    for(j in 1:replicates){ #ceiling(sample.n/2)
      pseudopop$m2m1.mean[j]<-sample(m2A.sample,size=1)/sample(m1A.sample,size=1)
      pseudopop$m3m1.mean[j]<-sample(m3A.sample,size=1)/sample(m1A.sample,size=1)
      pseudopop$MMC.mean[j]<-sample(m3L.sample,size=1)/sample(m1L.sample,size=1)
      pseudopop$MMC2.mean[j]<-sample(m2L.sample,size=1)/sample(m1L.sample,size=1)
    }
    sim.isolate$m2m1.SD[i]<-pseudopop$m2m1A %>% sd
    sim.isolate$m3m1.SD[i]<-pseudopop$m3m1A %>% sd
    sim.isolate$MMC.SD[i]<-pseudopop$m3m1L %>% sd
    sim.isolate$MMC2.SD[i]<-pseudopop$m2m1L %>% sd
  }
  
  sim.isolate$m2.m1A<-sim.isolate$m2.area/sim.isolate$m1.area
  sim.isolate$m3.m1A<-sim.isolate$m3.area/sim.isolate$m1.area
  sim.isolate$m3.m1L<-sim.isolate$m3.length/sim.isolate$m1.length
  sim.isolate$m2.m1L<-sim.isolate$m2.length/sim.isolate$m1.length
  if(sample.n==1){simulation.stats<-as.data.frame(sim.isolate)}
  if(sample.n>=2){simulation.stats<-rbind(simulation.stats,sim.isolate)}
}
save(simulation.stats,file = "output/composite_simulation_stats.RData")

# compare stats between simulated and real samples ------
#objects are simulation.stats and mouse.stats and ratio
sim.metric.names<-c("too.big","too.small","too.little","too.much")

for(ratio.choice in 1:nrow(ratio)){
  sim.metrical<-matrix(NA,nrow=50,ncol=length(sim.metric.names),
                     dimnames=list(c(NULL),sim.metric.names)) %>% as.data.frame #empty holder frame
  sim.metrical$N<-rep(1:50)
  for (i in 1:50){
    #with respect to sample size
    simulation.subsample<-simulation.stats[which(simulation.stats$N==i),]
    #How often is your simulated "composite" molar row outside of the 95% CI of the true population mean?
    sim.metrical$too.big[i]<-which(simulation.subsample[,ratio[ratio.choice,2]] > 
                                   CI[ratio[ratio.choice,2],"upper"]) %>% length()/replicates
    sim.metrical$too.small[i]<-which(simulation.subsample[,ratio[ratio.choice,2]] < 
                                     CI[ratio[ratio.choice],"lower"]) %>% length()/replicates
    #Can you model a CI by resampling isolated teeth? Make the "pseudosample" of isolated teeth
    #and use it to model some kind of confidence interval, is that interval greater than, equal to, or
    #smaller than the true confidence interval?
    sim.metrical$too.much[i]<-which(simulation.subsample[,ratio[ratio.choice,3]] > 
                                    CI[ratio[ratio.choice,3],"upper"]) %>% length()/replicates
    sim.metrical$too.little[i]<-which(simulation.subsample[,ratio[ratio.choice,3]] < 
                                      CI[ratio[ratio.choice,3],"lower"]) %>% length()/replicates
  }
  sim.metrical$ratio.choice<-ratio[ratio.choice,1]
  if(ratio.choice==1){sim.metrics<-sim.metrical}
  if(ratio.choice>=2){sim.metrics<-rbind(sim.metrics,sim.metrical)}
}
#tally totals
sim.metrics$outside.mean<-(sim.metrics$too.big + sim.metrics$too.small) %>%  round(.,3)
sim.metrics$outside.var<-(sim.metrics$too.much + sim.metrics$too.little) %>% round(.,3)
# reformat results for reporting purposes --------
#recast to write to a table
sim.metrics.to.write<-sim.metrics %>%
  select(N,ratio.choice,outside.mean,outside.var) %>%
  melt(.,id=c("N","ratio.choice")) %>%
  dcast(.,N ~ ratio.choice + variable)
write.csv(sim.metrics.to.write,
          paste("output/",species.name,"_composite_size_metrics.csv",sep=""),
          row.names = FALSE)

#melt to to plot
sim.melted<-simulation.stats %>% select(-m1.area,-m2.area,-m3.area,-m1.length,-m2.length,-m3.length) %>% 
  melt(.,id=c("N"))
sim.melted$tooth<-"M3"
sim.melted$tooth[grepl("2",sim.melted$variable,perl=TRUE)]<-"M2"
sim.melted$dimension<-"Area"
sim.melted$dimension[grepl("MMC",sim.melted$variable,perl=TRUE)]<-"Length"
sim.melted$statistic<-"Mean"
sim.melted$statistic[grepl("SD",sim.melted$variable,perl=TRUE)]<-"Standard Deviation"
sim.melted$CI.upper<-apply(sim.melted,1, function(x) CI[x[2],2])
sim.melted$CI.lower<-apply(sim.melted,1, function(x) CI[x[2],1])
sim.melted$outside<-(sim.melted$CI.lower > sim.melted$value) |
  (sim.melted$CI.upper < sim.melted$value)

composite.size.plot<-ggplot(data=sim.melted, aes(x=N,y=value))+
  facet_grid(facets = statistic + tooth ~ dimension, scales="free_y", switch="y") +
  geom_point(alpha=0.25,aes(color=outside)) +
  geom_smooth(color="black",method="auto") + #make smoother, prettier?
  scale_color_manual(values=c("gray","red")) +
  ggtitle(species.name) + theme_minimal() + 
  theme(legend.position="none",axis.title.y=element_blank(),strip.placement = "outside",
                            plot.title=element_text(face="italic"))
ggsave(composite.size.plot, filename = paste("output/","composite_sample_size.pdf",sep=""), dpi = 300)

## #Tables ------
#show distribution of CV, SD, Mean, Range for sample sizes for each metric?
## #plots --------
# source(paste(locateScripts,"molar_ratio_extant_plots.R"))
