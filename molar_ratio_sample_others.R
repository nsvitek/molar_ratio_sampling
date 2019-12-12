#re-do sample resampling on a cluster for previously published data of probably adequate sample size

# basic set-up ------------
# library(dplyr)
# library(ggplot2) #plot
# library(readxl) #read raw data
# library(reshape2) #for reformatting & comparing CV, measurement error data
# library(stringr) #using str_detect for sample size resampling in non-mice (and therefore in mice too)
# library(MCMCglmm) #for Bayesian modelling

#use format like:
library(dplyr,lib.loc="/gpfs/home/nvitek/R_packages")
library(ggplot2,lib.loc="/gpfs/home/nvitek/R_packages")
library(readxl,lib.loc="/gpfs/home/nvitek/R_packages")
library(reshape2,lib.loc="/gpfs/home/nvitek/R_packages")
library(stringr,lib.loc="/gpfs/home/nvitek/R_packages")
library(MCMCglmm,lib.loc="/gpfs/home/nvitek/R_packages")

MinSize<-41 #the result from Peromyscus dataset, how many you need to confidently reconstruct variance

replicates<-1000
ratio<-matrix(c("m3.m1L","MMC.mean","MMC.SD",
                "m2.m1L","MMC2.mean","MMC2.SD"),byrow = TRUE,ncol=3)
#reformat resampled.names to fit dataset
resampled.names<-c("N","MMC.mean","MMC.SD",#"MMC.max","MMC.min", #maxima and minima are here just in case you're curious
                   "MMC2.mean","MMC2.SD")#,"m2m1.max","m2m1.min")

# make new "mice" -----------
#comparative MMC data from Monson et al. 2019
compMMC.raw<-read_excel("input/comparison_MMC_monson.xlsx")
#calculate second ratio
compMMC.raw$m2.m1L<-compMMC.raw$`DM2L (mm)`/compMMC.raw$`DM1L (mm)`
compMMC.raw$m3.m1L<-compMMC.raw$MMC

#get sample size
MMC.spp.counts<-compMMC.raw %>% group_by(Species) %>% summarize (N=n()) 
MMC.spp.counts<-MMC.spp.counts %>% filter(N>=MinSize) #MinSize = minimum sample size that was adequate for all 4 metrics

for (i in 1:nrow(MMC.spp.counts)){
  #make dataset of new critter, substitute it in for mouse
  species.name<-MMC.spp.counts$Species[i]
  mouse<-compMMC.raw %>% filter(Species == species.name)
  #make table of summary stats for total mouse population, mouse.stats
  mouse.stats<-data.frame(Mu=rep(NA,nrow(ratio)),row.names=ratio[,1])
  mouse.stats$Mu<-mouse[,ratio[,1]] %>% apply(.,2,mean)
  mouse.stats$Sigma<-mouse[,ratio[,1]] %>% apply(.,2,sd)
  mouse.stats$CV<-mouse.stats$Sigma / mouse.stats$Mu * 100
  #confidence intervals, assuming normality, from Sokal and Rohlf
  mouse.stats$M.CI.U=mouse.stats$Mu + mouse.stats$Sigma * 1.96
  mouse.stats$M.CI.L=mouse.stats$Mu - mouse.stats$Sigma * 1.96
  mouse.stats$S.CI.L<-(((nrow(mouse)-1) * (mouse.stats$Sigma)^2) /
                         qchisq(c(0.975),df=nrow(mouse)-1)) %>% sqrt
  mouse.stats$S.CI.U<-(((nrow(mouse)-1) * (mouse.stats$Sigma)^2) /
                         qchisq(c(0.025),df=nrow(mouse)-1)) %>% sqrt
  
  #and take off running
  # source(paste(locateScripts,"molar_ratio_sample_size.R",sep="/"))
  source("molar_ratio_sample_size.R")
}

# repeat for areas ------------
ratio<-matrix(c("m3.m1A","m3m1.mean","m3m1.SD",
                "m2.m1A","m2m1.mean","m2m1.SD"),byrow = TRUE,ncol=3)
#reformat resampled.names to fit dataset
resampled.names<-c("N","m3m1.mean","m3m1.SD",#"m3m1.max","m3m1.min",
                   "m2m1.mean","m2m1.SD")#,"m2m1.max","m2m1.min")

#get sample size
compICM2.raw<-read_excel("input/comparison_ICM_roseman.xlsx", sheet = "Molar Size Data", na="NA")
#get ratios, other info for Roseman & Delezene data
compICM2.raw$m3.m1A<-compICM2.raw$M3.Area/compICM2.raw$M1.Area
compICM2.raw$m2.m1A<-compICM2.raw$M2.Area/compICM2.raw$M1.Area
compICM2<-compICM2.raw[!is.na(compICM2.raw$m3.m1A),]
compICM.pub<-compICM2[!is.na(compICM2$m2.m1A),]

ICM.spp.counts<-compICM.pub %>% group_by(Species, Sex) %>% summarize (N=n()) 
ICM.spp.counts<-ICM.spp.counts %>% filter(N>=MinSize) #MinSize = minimum sample size that was adequate for all 4 metrics

for (i in 1:nrow(ICM.spp.counts)){
  #make dataset of new critter, substitute it in for mouse
  species.name<-ICM.spp.counts$Species[i]
  species.sex<-ICM.spp.counts$Sex[i]
  mouse<-compICM.pub %>% filter(Species == species.name, Sex == species.sex)
  
  #make table of summary stats for total mouse population, mouse.stats
  mouse.stats<-data.frame(Mu=rep(NA,nrow(ratio)),row.names=ratio[,1])
  mouse.stats$Mu<-mouse[,ratio[,1]] %>% apply(.,2,mean)
  mouse.stats$Sigma<-mouse[,ratio[,1]] %>% apply(.,2,sd)
  mouse.stats$CV<-mouse.stats$Sigma / mouse.stats$Mu * 100
  #confidence intervals, assuming normality, from Sokal and Rohlf
  mouse.stats$M.CI.U=mouse.stats$Mu + mouse.stats$Sigma * 1.96
  mouse.stats$M.CI.L=mouse.stats$Mu - mouse.stats$Sigma * 1.96
  mouse.stats$S.CI.L<-(((nrow(mouse)-1) * (mouse.stats$Sigma)^2) /
                         qchisq(c(0.975),df=nrow(mouse)-1)) %>% sqrt
  mouse.stats$S.CI.U<-(((nrow(mouse)-1) * (mouse.stats$Sigma)^2) /
                         qchisq(c(0.025),df=nrow(mouse)-1)) %>% sqrt
  
  #and take off running
  # source(paste(locateScripts,"molar_ratio_sample_size.R",sep="/"))
  source("molar_ratio_sample_size.R")
}