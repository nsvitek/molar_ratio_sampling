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
MMC.spp.counts<-MMC.spp.counts %>% filter(N>=43) #43 = minimum sample size that was adequate for all 4 metrics

for (i in 1:nrow(MMC.spp.counts)){
  #make dataset of new critter, substitute it in for mouse
  species.name<-MMC.spp.counts$Species[i]
  mouse<-compMMC.raw %>% filter(Species == species.name)
  #and take off running
  source("molar_ratio_sample_model.R")
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
ICM.spp.counts<-compICM.pub %>% group_by(Species, Sex) %>% summarize (N=n()) 
ICM.spp.counts<-ICM.spp.counts %>% filter(N>=43) #43 = minimum sample size that was adequate for all 4 metrics

for (i in 1:nrow(ICM.spp.counts)){
  #make dataset of new critter, substitute it in for mouse
  species.name<-ICM.spp.counts$Species[i]
  species.sex<-ICM.spp.counts$Sex[i]
  mouse<-compICM.pub %>% filter(Species == species.name, Sex == species.sex)
  #and take off running
  source("molar_ratio_sample_model.R")
  # source(paste(locateScripts,"molar_ratio_sample_size.R",sep="/"))
  source("molar_ratio_sample_size.R")
}