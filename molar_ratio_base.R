# !diagnostics off

# Use molar rows of Recent Peromyscus gossypinus & others as comparison to fossil record
# dependencies -----
library(dplyr) #organize
library(ggplot2) #plot
# library(gridExtra) #plot
library(ggforce) #allows drawing circles and ellipses
library(scales) #allows use of standard colors. May not be completely necessary
library(readxl) #read raw data
library(reshape2) #for reformatting & comparing CV, measurement error data
# library(ggthemes) #to get Paul Tol colors
library(lmodel2) #for RMA regression
library(moments) #for kurtosis

locateScripts<-"C:/cygwin/home/N.S/scripts/molar_ratio_sampling"
# locateScripts<-"C:/scripts/scripts"
locateData<-"C:/Users/N.S/Dropbox/Documents/research/vitek-etal_inhibitory-cascade-isolated"
# locateData<-"D:/Dropbox/Documents/research/vitek-etal_inhibitory-cascade-isolated"

setwd(locateScripts)
source("../scripts/function_bootstrap.R")
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
mouse$Species<-"Peromyscus"
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

#add in peromyscus gossypinus data
compICM2add<-dplyr::select(mouse,Order,m2.m1A,m3.m1A)
compICM2add$Species<-"Peromyscus gossypinus"
compICM2add$Sex<-"both"
compICM2add$M1.Area<-mouse$m1.area
compICM2add$M2.Area<-mouse$m2.area
compICM2add$M3.Area<-mouse$m3.area
compICM.indv<-rbind(compICM2,compICM2add)

#summarize individual-level data so it can be collated with population-level data
ICM.CV<- compICM.indv %>% group_by(Order, Species, Sex) %>%  
  summarize(N=n(), m3m1.mean=mean(m3.m1A), m3m1.SD=sd(m3.m1A),m3m1.CV=sd(m3.m1A)/mean(m3.m1A)*100,
            m2m1.mean=mean(m2.m1A), m2m1.SD=sd(m2.m1A),m2m1.CV=sd(m2.m1A)/mean(m2.m1A)*100)
colnames(ICM.CV)[which(colnames(ICM.CV)=="Sex")]<-"condition"

#organize collated population-level ICM ratio data into standardized set of variables
compICM.raw$N<-as.integer(compICM.raw$N)
compICM.raw$m2m1.CV<-(compICM.raw$m2m1.SD/compICM.raw$m2m1.mean) *100
compICM.raw$m3m1.CV<-(compICM.raw$m3m1.SD/compICM.raw$m3m1.mean) *100
compICM.pop<-dplyr::select(compICM.raw,-citation) %>% bind_rows(.,ICM.CV)

#confidence intervals?
compICM.pop$m2m1.CI<-compICM.pop$m2m1.SD*1.96
compICM.pop$m3m1.CI<-compICM.pop$m3m1.SD*1.96
# format MMC comparative data -----
compMMC.raw %>% head
#MMC: add Peromyscus gossypinus to comparative dataset

mouse.MMC.2add<-dplyr::select(mouse,m3.m1L)
mouse.MMC.2add$Order<-"Rodentia"
mouse.MMC.2add$Family<-"Cricetidae"
mouse.MMC.2add$Species<-"Peromyscus gossypinus"
mouse.MMC.2add$`Museum ID`<-mouse.MMC.2add$`Museum Repository`<-NA
mouse.MMC.2add$`DM3L (mm)`<-mouse$m3.length
mouse.MMC.2add$`DM2L (mm)`<-mouse$m2.length
mouse.MMC.2add$`DM1L (mm)`<-mouse$m1.length
mouse.MMC.2add$`DP4L (mm)`<-NA
mouse.MMC.2add$MMC<-mouse$m3.m1L
mouse.MMC.2add$PMM<-mouse.MMC.2add$`Measured by`<-NA

mouse.MMC.2add<-dplyr::select(mouse.MMC.2add,-m3.m1L)

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

# standing variation between populations -----
#gossypinus: compare the 4 state-level populations using U test and bootstrap.
# #visualize intraspecific gossypinus data
# # cairo_pdf("output/MMC_population_gossy.pdf")
# ggplot(mouse,aes(x=state, y=m3.m1L, fill=county))+
#   geom_hline(yintercept=mean(mouse$m3.m1L),linetype="dashed",color="black",size=1.5)+
#   geom_dotplot(binaxis="y")+
#   theme_minimal()
# # dev.off()
# 
# # cairo_pdf("output/ICM_population_gossy.pdf")
# ggplot(mouse,aes(x=m2.m1A, y=m3.m1A, colour=state))+
#   stat_chull(fill=NA)+
#   geom_point(size=1.5)+
#   theme_minimal()
# # dev.off()
# #plot to look at distribution. Structuring by site?
# ggplot() + #theme(legend.position=c(.7,.1)) +
#   geom_point(data=mouse,aes(x=m2.m1A, y=m3.m1A,fill=state),size=4,pch=21) 

#define testing variables
measurement<-mouse$m3.m1L #which variable are you using to measure size
measurement.ICM1<-mouse$m3.m1A
measurement.ICM2<-mouse$m2.m1A
group.var<-mouse$state  #which variable are you using as your group identity
pop.sizes<-mouse %>% group_by(state) %>% summarize(n()) #calculate sample size for each group

#pairwise U-test with Bonferroni correction approach following Asahara 2014
U.MMC<-Pairwise.U(measurement,group.var)
U.ICM1<-Pairwise.U(measurement.ICM1,group.var)
U.ICM2<-Pairwise.U(measurement.ICM1,group.var)

U.MMC$U.p.vals.ICM31<-U.ICM1$U.p.vals
U.MMC$U.p.vals.ICM21<-U.ICM2$U.p.vals
write.csv(U.MMC,"output/gossypinus_population_Utest.csv")

#bootstrap ANOVA approach
resam.distribution<-BootstrapANOVA(measurement,group.var,pop.sizes,replicates)
bootstrap.p<-which(resam.distribution>=resam.distribution[1]) %>% length()/replicates

resam.distribution.ICM1<-BootstrapANOVA(measurement.ICM1,group.var,pop.sizes,replicates)
bootstrap.p.ICM1<-which(resam.distribution.ICM1>=resam.distribution.ICM1[1]) %>% length()/replicates

resam.distribution.ICM2<-BootstrapANOVA(measurement.ICM2,group.var,pop.sizes,replicates)
bootstrap.p.ICM2<-which(resam.distribution.ICM2>=resam.distribution.ICM2[1]) %>% length()/replicates
cbind(bootstrap.p,bootstrap.p.ICM1,bootstrap.p.ICM2) %>% 
  write.csv(.,"out/gossypinus_population_bootstrapANOVA.csv")

# #visualize
# hist(resam.distribution.ICM2) #no significant difference from random
# abline(v=resam.distribution.ICM2[1],col="red")
# aov(lm(measurement.ICM1~group.var)) %>% summary

#calculate mean and coefficient of variation by population
ratio.pop.CV<- mouse %>% group_by(state) %>% 
  summarize(N=n(),CV.MMC=sd(m3.m1L)/mean(m3.m1L)*100,mean.MMC=mean(m3.m1L),
            CV.m3m1=sd(m3.m1A)/mean(m3.m1A)*100,mean.m3m1=mean(m3.m1A),
            CV.m2m1=sd(m2.m1A)/mean(m2.m1A)*100,mean.m3m1=mean(m2.m1A))


# standing variation in species -----
#MMC: calculate CV by species from Monson, gossypinus data
MMC.CV<- compMMC %>% group_by(Order, Species) %>% 
  summarize(N=n(), MMC.mean=mean(MMC), MMC.SD=sd(MMC), MMC.CV=sd(MMC)/mean(MMC)*100)
MMC.CV<-MMC.CV[complete.cases(MMC.CV),]

#summarize, report.
MMC<-ungroup(MMC.CV) %>% summarize(mean=mean(MMC.CV),SD=sd(MMC.CV),median=median(MMC.CV))
M2.M1<-compICM.pop %>% summarize(mean=mean(m2m1.CV),SD=sd(m2m1.CV),median=median(m2m1.CV))
M2.M1.noRaccoonDog<-compICM.pop %>% dplyr::filter(.,Species!="Nyctereutes procyonoides") %>% 
  summarize(mean=mean(m2m1.CV),SD=sd(m2m1.CV),median=median(m2m1.CV))
M3.M1<-compICM.pop %>% summarize(mean=mean(m3m1.CV),SD=sd(m3m1.CV),median=median(m3m1.CV))
M3.M1.noRaccoonDog<-compICM.pop %>% dplyr::filter(.,Species!="Nyctereutes procyonoides") %>% 
  summarize(mean=mean(m3m1.CV),SD=sd(m3m1.CV),median=median(m3m1.CV))

CV.summary<-rbind(MMC,M2.M1, M3.M1,M2.M1.noRaccoonDog,M3.M1.noRaccoonDog) %>% as.data.frame
CV.summary$CI.upper<-CV.summary$mean+CV.summary$SD*1.96
row.names(CV.summary)<-c("MMC","M2.M1", "M3.M1","M2.M1.noRaccoonDog","M3.M1.noRaccoonDog")
write.csv(CV.summary,"output/CV_summary_stats.csv")
# sensitivity of mean and sd of ratios to sample size -----------

#creates resampled.stats
#Both: resample gossypinus dataset for all 3 ratios, check for sensitivity to sample size
resampled.names<-c("N","MMC.CV","MMC.mean","MMC.SD",
                   "m3m1.CV","m3m1.mean","m3m1.SD","m2m1.CV","m2m1.mean","m2m1.SD",
                   "MMC.max","MMC.min","m3m1.max","m3m1.min","m2m1.max","m2m1.min")
resampled.base<-matrix(NA,nrow=replicates,ncol=length(resampled.names),
                       dimnames=list(c(NULL),resampled.names))
for (sample.n in 2:nrow(mouse)){
  resampled<-resampled.base
  resampled[,"N"]<-sample.n #sample size
  for(i in 1:replicates){
    mouse.resampled<-dplyr::sample_n(mouse,size=sample.n,replace=FALSE)
    resampled[i,"MMC.CV"]<-sd(mouse.resampled$m3.m1L)/mean(mouse.resampled$m3.m1L)*100
    resampled[i,"MMC.mean"]<-mean(mouse.resampled$m3.m1L)
    resampled[i,"MMC.SD"]<-sd(mouse.resampled$m3.m1L)
    resampled[i,"m3m1.CV"]<-sd(mouse.resampled$m3.m1A)/mean(mouse.resampled$m3.m1A)*100
    resampled[i,"m2m1.CV"]<-sd(mouse.resampled$m2.m1A)/mean(mouse.resampled$m2.m1A)*100
    resampled[i,"m3m1.mean"]<-mean(mouse.resampled$m3.m1A)
    resampled[i,"m2m1.mean"]<-mean(mouse.resampled$m2.m1A)
    resampled[i,"m3m1.SD"]<-sd(mouse.resampled$m3.m1A)
    resampled[i,"m2m1.SD"]<-sd(mouse.resampled$m2.m1A)
    resampled[i,"m3m1.max"]<-max(mouse.resampled$m3.m1A)
    resampled[i,"m2m1.max"]<-max(mouse.resampled$m2.m1A)
    resampled[i,"m3m1.min"]<-min(mouse.resampled$m3.m1A)
    resampled[i,"m2m1.min"]<-min(mouse.resampled$m2.m1A)
    resampled[i,"MMC.max"]<-max(mouse.resampled$m3.m1L)
    resampled[i,"MMC.min"]<-min(mouse.resampled$m3.m1L)
    
  }
  if(sample.n==2){resampled.stats<-as.data.frame(resampled)}
  if(sample.n>=3){resampled.stats<-rbind(resampled.stats,resampled)}
}

#make a table matching names to ratio
ratio<-matrix(c("m3.m1L","MMC.mean","MMC.SD","resampled.range.MMC","m3.m1L.SD","V.m3.m1L",
                "m3.m1A","m3m1.mean","m3m1.SD","resampled.range.m3m1","m3.m1A.SD","V.m3.m1A",
                "m2.m1A","m2m1.mean","m2m1.SD","resampled.range.m2m1","m2.m1A.SD","V.m2.m1A"),byrow = TRUE,ncol=6)

#calculate summary stats for original mouse dataset to use as standar
#do for all three ratios. Write to a table. 
for(ratio.choice in 1:nrow(ratio)){
  true.numbers<-mouse[,ratio[ratio.choice,1]] %>% unlist
  mouse.standard<-data.frame(Mu=round(mean(true.numbers),3),Sigma=round(sd(true.numbers),3),
                             CV=round(sd(true.numbers)/mean(true.numbers)*100,3),
                             CI.U=round(mean(true.numbers)+sd(true.numbers)*1.96,3),
                             CI.L=round(mean(true.numbers)-sd(true.numbers)*1.96,3),
                             Range=round(max(true.numbers)-min(true.numbers),3),
                             Kurtosis=round(kurtosis(true.numbers),3))
  mouse.standard$CI.range<-mouse.standard$CI.U-mouse.standard$CI.L
  #CI for mouse variance, assuming normality, from Sokal and Rohlf
  mouse.standard$var.CI.L<-(((nrow(mouse)-1) * (mouse.standard$Sigma)^2) /
              qchisq(c(0.975),df=nrow(mouse)-1)) %>% sqrt
  mouse.standard$var.CI.U<-(((nrow(mouse)-1) * (mouse.standard$Sigma)^2) /
              qchisq(c(0.025),df=nrow(mouse)-1)) %>% sqrt
  #create a vector of mean error
  vs.true<-data.frame(prop.mean.outside.95ci=rep(NA,nrow(mouse)-1),
                      prop.sd.above.95ci=rep(NA,nrow(mouse)-1),
                      prop.sd.below.95ci=rep(NA,nrow(mouse)-1),
                      upper.mean=rep(NA,nrow(mouse)-1),
                      lower.mean=rep(NA,nrow(mouse)-1),
                      N=rep(2:nrow(mouse)))
  for (i in 2:nrow(mouse)){
    resampled.means<-resampled.stats[which(resampled.stats$N==i),ratio[ratio.choice,2]] %>% unlist
    vs.true$prop.mean.outside.95ci[i-1]<-which(abs(resampled.means-mouse.standard$Mu)>=(mouse.standard$Sigma*1.96)) %>% 
      length(.)/replicates
    #get range of means
    vs.true$upper.mean[i-1]<-max(resampled.means)
    vs.true$lower.mean[i-1]<-min(resampled.means)
    resampled.standarddevs<-resampled.stats[which(resampled.stats$N==i),ratio[ratio.choice,3]] %>% unlist
    vs.true$prop.sd.above.95ci[i-1]<-which((resampled.standarddevs-mouse.standard$var.CI.U)>0) %>% 
      length(.)/replicates
    vs.true$prop.sd.below.95ci[i-1]<-which((resampled.standarddevs-mouse.standard$var.CI.L)<0) %>% 
      length(.)/replicates
  }
  vs.true$ratio.name<-ratio[ratio.choice,1]
  if(ratio.choice==1){    
    mouse.stats<-mouse.standard
    sample.size.indicators<-vs.true
    }
  if(ratio.choice>1){
    mouse.stats<-rbind(mouse.stats,mouse.standard)
    sample.size.indicators<-rbind(sample.size.indicators,vs.true)
    }
}
#melt to to plot
ssi.long<-sample.size.indicators %>% dplyr::select(.,-upper.mean, -lower.mean) %>% 
  melt(,id=c("N","ratio.name"))

#recast to write to a table
sample.size.indicators.to.write<-melt(sample.size.indicators,id=c("N","ratio.name")) %>% 
  dcast(.,N ~ ratio.name + variable)
write.csv(sample.size.indicators.to.write,"output/sample_size_metrics.csv")
# #plot
# source(paste(locateScripts,"molar_ratio_extant_plots.R"))

# # test for phylogenetic signal in CV --------------
# library(phytools)
# library(treeman)
# 
# data(mammals)
# str(mammals)
# MMC.CV$Species
# MMC.CV$MatchSpp<-gsub("([A-Za-z]*) ([A-Za-z]*)","\\1_\\2",MMC.CV$Species,perl=TRUE)
# 
# #prune the tree
# mammalsMMC<-rmTips(mammals, tids = mammals@tips[-which(mammals@tips %in% MMC.CV$MatchSpp)]) %>%
#   as(.,"phylo")
# plot(mammalsMMC) #chec: a reasonable tree?
# rm(rm="mammals") #is an 8mb dataset, remove when finished with it
# #create variable for testing phylogenetic signal
# phylo.MMC<-MMC.CV$MMC.CV[which(MMC.CV$MatchSpp %in% mammalsMMC$tip.label)]
# names(phylo.MMC)<-MMC.CV$MatchSpp[which(MMC.CV$MatchSpp %in% mammalsMMC$tip.label)]
# #test
# CV.K<-phylosig(mammalsMMC, phylo.MMC, method="K", test=TRUE, nsim=1000)
# CV.L<-phylosig(mammalsMMC, phylo.MMC, method="lambda", test=TRUE, nsim=1000)
# 
# #likelihood ratio test
# fitBrownian<-brownie.lite(paintSubTree(mammalsMMC,mammalsMMC$edge[1,1],state="1"),phylo.MMC)
# lam1<-phylosig(mammalsMMC,x1,method="lambda")
# LR<-2*(CV.L$logL-fitBrownian$logL1)
# pchisq(LR,df=1,lower.tail=FALSE) #significantly different from "phylogenetic signal", or pure Brownian motion
# create simulated composite molar rows --------
#create resampled pseudoreplicates of molars sizes and ratios from composite cheek teeth
#both, simulate all ratios at once:
sim.isolate.names<-c("N","m1.area","m2.area","m3.area","m1.length","m2.length","m3.length",
                     # "m1m2.area.C","m2m3.area.C","m1m3.area.C",
                     # "m1m2.length.C","m2m3.length.C","m1m3.length.C",
                     "m2.m1A.SD","m3.m1A.SD","m3.m1L.SD",
                     "m2.m1A","m3.m1A","m3.m1L")
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
    # sim.isolate$m1m2.area.C[i]<-cov(m1A.sample,m2A.sample)
    # sim.isolate$m2m3.area.C[i]<-cov(m2A.sample,m3A.sample)
    # sim.isolate$m1m3.area.C[i]<-cov(m1A.sample,m3A.sample)
    # sim.isolate$m1m2.length.C[i]<-cov(m1L.sample,m2L.sample)
    # sim.isolate$m2m3.length.C[i]<-cov(m2L.sample,m3L.sample)
    # sim.isolate$m1m3.length.C[i]<-cov(m1L.sample,m3L.sample)
    # print(paste("pseudopopulating rep", i))
    #make a pseudosample from the isolated molars by resampling 50% of specimen each time
    pseudopop.names<-c("m2m1A","m3m1A","m3m1L")
    pseudopop<-matrix(NA,nrow=replicates,ncol=length(pseudopop.names),
                             dimnames=list(c(NULL),pseudopop.names)) %>% as.data.frame #empty holder frame
    for(j in 1:replicates){ #ceiling(sample.n/2)
      pseudopop$m2m1A[j]<-sample(m2A.sample,size=1)/
        sample(m1A.sample,size=1)
      pseudopop$m3m1A[j]<-sample(m3A.sample,size=1)/
        sample(m1A.sample,size=1)
      pseudopop$m3m1L[j]<-sample(m3L.sample,size=1)/
        sample(m1L.sample,size=1)
    }
    sim.isolate$m2.m1A.SD[i]<-pseudopop$m2m1A %>% sd
    sim.isolate$m3.m1A.SD[i]<-pseudopop$m3m1A %>% sd
    sim.isolate$m3.m1L.SD[i]<-pseudopop$m3m1L %>% sd
  }
  
  sim.isolate$m2.m1A<-sim.isolate$m2.area/sim.isolate$m1.area
  sim.isolate$m3.m1A<-sim.isolate$m3.area/sim.isolate$m1.area
  sim.isolate$m3.m1L<-sim.isolate$m3.length/sim.isolate$m1.length
  if(sample.n==1){simulation.stats<-as.data.frame(sim.isolate)}
  if(sample.n>=2){simulation.stats<-rbind(simulation.stats,sim.isolate)}
}

## compare stats between simulated and real samples ------
#objects are simulation.stats and mouse.stats and ratio
sim.metric.names<-c("too.big","too.small","too.little","too.much")
#mouse.ratio.estimated.table comes from molar_ratio_ICM_expectations.R
SD.CI<-HPDinterval(as.mcmc(mouse.ratio.estimated.table)) %>% sqrt %>% round(.,5) %>% as.data.frame 

for(ratio.choice in 1:nrow(ratio)){
  sim.metric<-matrix(NA,nrow=50,ncol=length(sim.metric.names),
                     dimnames=list(c(NULL),sim.metric.names)) %>% as.data.frame #empty holder frame
  for (i in 1:50){
    #with respect to sample size
    simulation.subsample<-simulation.stats[which(simulation.stats$N==i),]
    #How often is your simulated "composite" molar row outside of the 95% CI of the true population mean?
    sim.metric$too.big[i]<-which(simulation.subsample[,ratio[ratio.choice,1]] > mouse.stats$CI.U[ratio.choice]) %>% 
      length()
    sim.metric$too.small[i]<-which(simulation.subsample[,ratio[ratio.choice,1]] < mouse.stats$CI.L[ratio.choice]) %>% 
      length()
    #Can you model a CI by resampling isolated teeth? Make the "pseudosample" of isolated teeth
    #and use it to model some kind of confidence interval, is that interval greater than, equal to, or
    #smaller than the true confidence interval?
    sim.metric$too.much[i]<-which(simulation.subsample[,ratio[ratio.choice,5]] > 
                                    SD.CI[ratio[ratio.choice,6],"upper"]) %>% length
    sim.metric$too.little[i]<-which(simulation.subsample[,ratio[ratio.choice,5]] < 
                                      SD.CI[ratio[ratio.choice,6],"lower"]) %>% 
      length
  }
  sim.metric$ratio.choice<-ratio[ratio.choice,1]
  if(ratio.choice==1){sim.metrics<-sim.metric}
  if(ratio.choice>=2){sim.metrics<-rbind(sim.metrics,sim.metric)}
}

#recast to write to a table
# sample.size.indicators.to.write<-melt(sample.size.indicators,id=c("N","ratio.name")) %>% 
#   dcast(.,N ~ ratio.name + variable)
write.csv(sim.metrics,"output/resample_metrics.csv")

ggplot(data=simulation.stats, aes_string(x="N",y="m3.m1A.SD"))+
  geom_point(alpha=0.25,color="gray") +
  # geom_smooth(color="black",method="loess") + #make smoother, prettier? or quantile?
  geom_hline(yintercept=SD.CI["V.m3.m1A","lower"],linetype="dashed",color="blue",size=1)+
  geom_hline(yintercept=SD.CI["V.m3.m1A","upper"],linetype="dashed",color="blue",size=1) 

ggplot() +
  geom_point(data=simulation.subsample,aes(x=m2.m1A, y=m3.m1A))+
  geom_point(data=pseudopop, aes(x=m2m1A,y=m3m1A),color="blue")+
  geom_point(data=mouse,aes(x=m2.m1A,y=m3.m1A),size=3,color="red")

#If yes, Does this creation of a "pseudosample" also mimic the covariance structure of the true sample?
#Use MCMCGlmm results of posterior distribution. Take mode/HPD for covariance for "true" sample
#Take the speudosample to calculate covariances.
#subtract that covariance from the "true" modal covariance and see if it's in the HPD
# HPDinterval(as.mcmc(mouseA.estimated.table)) %>% round(.,5)
# cbind(posterior.mode(as.mcmc(estimated.table)),
#       HPDinterval(as.mcmc(estimated.table)))

## #compare means -----
t.test(mouse$m2.m1A, simulation.stats$m2.m1A[which(simulation.stats$N==1)],
       alternative="two.sided",paired=FALSE,var.equal=FALSE) #%>% 
  # capture.output(.,file="inhibitory_results/gossypinus_ttest_m2m1.txt")
t.test(mouse$m3.m1, sim.isolate.area$m3.m1,alternative="two.sided",paired=FALSE,var.equal=FALSE) %>% 
  capture.output(.,file="inhibitory_results/gossypinus_ttest_m3m1.txt")

## #evaluate quality of fossil dataset -------
## #plot data ------

#make the two polygons of "forbidden space" in Kavanagh model
poly.right<-data.frame(id=rep("one",4),x=c(1,1,2.2,2.2),y=c(0,1,2,0))
poly.left<-data.frame(id=rep("two",4),x=c(0,1,1,0),y=c(0,1,2.2,2.2))


#plot complete molar rows vs. simulated composite molar rows
cairo_pdf("inhibitory_results/gossypinus_ICM.pdf",width=4,height=4)
ggplot()+ 
  # geom_polygon(data=poly.right,aes(x=x,y=y)) +
  # geom_polygon(data=poly.left,aes(x=x,y=y)) +
  geom_point(data=simulation.stats[which(simulation.stats$N==10),],
             aes(x=m2.m1A,y=m3.m1A),alpha=0.5,fill="gray20",color="gray20",pch=21) +
  geom_point(data=mouse,aes(x=m2.m1A, y=m3.m1A, pch=state),size=4,fill=hue_pal()(7)[1]) +
  # scale_shape_manual(values=c(21,22,23))+
  # coord_cartesian(xlim = c(0.7,0.95),ylim = c (0.5,0.75))+ 
  xlab("M/2:M/1")+ylab("M/3:M/1") +
  theme_minimal() +theme(legend.position = "none")
dev.off()
embedFonts("inhibitory_results/gossypinus_ICM.pdf")

#repeat, but with lengths instead of areas
ggplot()+ 
  # geom_polygon(data=poly.right,aes(x=x,y=y)) +
  # geom_polygon(data=poly.left,aes(x=x,y=y)) +
  geom_point(data=sim.isolate.length,aes(x=m2.m1,y=m3.m1),alpha=0.5,fill="gray20",color="gray20",pch=21) +
  geom_point(data=mouse,aes(x=m2.m1L, y=m3.m1L, pch=Institution),size=4,fill=hue_pal()(7)[1]) +
  scale_shape_manual(values=c(21,22,23))+
  # coord_cartesian(xlim = c(0.5,0.95),ylim = c (0.5,0.95))+
  xlab("M/2:M/1")+ylab("M/3:M/1") +
  theme_minimal() +theme(legend.position = "none")

# #compare confidence intervals to distribution
ggplot()+
  geom_polygon(data=poly.right,aes(x=x,y=y)) +
  geom_polygon(data=poly.left,aes(x=x,y=y)) +
  geom_circle(data=raw.ci,aes(x0=m2.m1,y0=m3.m1,r=ci95),fill="blue",alpha=0.4) +
  geom_circle(data=sim.a.ci,aes(x0=m2.m1,y0=m3.m1,r=ci95),fill="yellow",alpha=0.4) +
  geom_point(data=sim.isolate.area,aes(x=m2.m1,y=m3.m1),alpha=0.5,fill="gray20",color="gray20",pch=21) +
  geom_point(data=mouse,aes(x=m2.m1, y=m3.m1),size=4,pch=21,fill=hue_pal()(7)[1]) +
  coord_cartesian(xlim = c(0.6,1),ylim = c (0.5,0.8))+
  theme_minimal()

## #NOTES ------
#Test from Roseman and Delezene 2019: if ICM working, M3 = 2*M2 - M1., or M3 = 2 * M2 + (-1) * M1
#The other mathematical predictions involve covariances, which cannot do with isolated molars,
#could possibly use equation 3 if you have a bunch of two-tooth jaws.
#Schroer and Wood: this is the equation that ends up with a slope of 2 and an intercept of -1
#stated another way by Bernal: [M3/M1 = 2 * (M2/M1) - 1]
#original statement in Kavanagh is 1+[(a-i)/i](x-1), where a=activator, i=inhibitor, x=tooth position
#so M1 = 1
#and M2 = 1+(a-i)/i, or 1+(a/i)-1 or a/i
#and M3 = 1+(a/i)*2-(i/i)*2 or 2*(a/i) - 1, which seems to most nicely translate to R&D's substitution
#or Bernal's ratio substitutions: These equations are all how you get the line. 
#so how you get Polly's model-consistent space? In PDP's space, for example, you can have
#m2m1 = 1 and m3m1 = 0, or (0.5, 0.5), or (1, 0.5)
#that is, that "2" in the IC model could be any number...between...
#and the "1" can be between 0 and 1?
#it seem to be that m2m1 = 1 and m3m1 = m2m1, if you're between those in a certain way, then you're okay
#m3m1 < m2m1 < 1 or
#1 < m3m1 < m2m1 describes the space, which is to say a space consistent with 
#cumulative change (words of H&G) or directional change resulting from a developmental cascade (PDP)

## #Tables ------
#show distribution of CV, SD, Mena, Range for sample sizes for each metric?