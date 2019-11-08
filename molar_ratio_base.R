# !diagnostics off

# Use molar rows of Recent Peromyscus gossypinus & others as comparison to fossil record
# dependencies -----
library(dplyr) #organize
library(ggplot2) #plot
# library(gridExtra) #plot
library(ggforce) #allows drawing circles and ellipses
library(scales) #allows use of standard colors. May not be completely necessary
library(readxl) #read raw data
# library(reshape2) #for reformatting & comparing CV data
# library(ggthemes) #to get Paul Tol colors
# library(lmodel2)

locateScripts<-"C:/cygwin/home/N.S/scripts/molar-ratio-sampling"
# locateScripts<-"C:/scripts/scripts"
locateData<-"C:/Users/N.S/Dropbox/Documents/research/vitek-etal_inhibitory-cascade-isolated"
# locateData<-"D:/Dropbox/Documents/research/vitek-etal_inhibitory-cascade-isolated"

setwd(locateScripts)
source("function_bootstrap.R")

# functions ------
#calculate a 95% normal confidence interval from a sample, or resample
calc.ci<-function(sample){ #function to get elements of 95% CI circle
  #calculate mean of sample
  mean.xy<-c(mean(sample$m2.m1), mean(sample$m3.m1)) %>% 
    matrix(nrow=1,dimnames=list("mean",c("m2.m1","m3.m1"))) %>% data.frame
  #calculate distances between each subsample and mean, 
  distances<-((sample$m2.m1-mean.xy$m2.m1)^2 + (sample$m3.m1-mean.xy$m3.m1)^2) %>% 
    sqrt()
  #calculate  the 95% CI for those distances
  ci95<-quantile(distances,probs=c(0.95),na.rm=TRUE)
  #make into a data frame for ggplot
  ci<-cbind(mean.xy,ci95)
  return(ci)
}

# settings ------
replicates<-9999 #bootstrap number
n<-10 #resample number

# load data -----------
setwd(locateData)
#read in raw data
mouse.raw<-read.csv("inhibitory_results/input/gossypinus_measurements.csv")
#read in comparative data
compICM.raw<-read_excel("inhibitory_results/input/ICM_comparative_variance.xlsx")
compMMC.raw<-read_excel("inhibitory_results/input/ICM_macrocranion_test.xlsx")
#read in the data of macrocranion
fossil.raw<-read_excel("inhibitory_results/input/ICM_macrocranion_test.xlsx")
#read in data pertaining to CV of raw lengths
# comp.areaCV.raw
# comp.lengthCV.raw
## #format cotton mouse --------
#using crown area ICM:
mouse.raw$m1<-mouse.raw$m1.length*mouse.raw$m1.width
mouse.raw$m2<-mouse.raw$m2.length*mouse.raw$m2.width
mouse.raw$m3<-mouse.raw$m3.length*mouse.raw$m3.width
mouse.raw$m2.m1<-mouse.raw$m2/mouse.raw$m1
mouse.raw$m3.m1<-mouse.raw$m3/mouse.raw$m1

#using MMC:
mouse.raw$m2.m1L<-mouse.raw$m2.length/mouse.raw$m1.length
mouse.raw$m3.m1L<-mouse.raw$m3.length/mouse.raw$m1.length


## # test assumption of normality -----
#In Kavanagh et al. 2007, the area metrics are not transformed
#Does that protocol seems okay in terms of normality, according to tests?
shapiro.test(mouse.raw$m3.m1)
shapiro.test(mouse.raw$m2.m1)
#In Hlusko et al. 2016, metrics are "inverse-normalized before analysis to account for residual kurtosis"
#and ratios are not strictly normal (according to tests below)
shapiro.test(mouse.raw$m3.m1L)
shapiro.test(mouse.raw$m2.m1L)

hist(mouse.raw$m2.m1)
hist(mouse.raw$m2.m1L)

## #format ICM comparative data ------
#will work for ratios of area. 
comp.raw$m2.m1CV<-(comp.raw$m2.m1variance/comp.raw$m2.m1) *100
comp.raw$m3.m1CV<-(comp.raw$m3.m1variance/comp.raw$m3.m1) *100

#confidence intervals?
#sd<-sqrt(var) and 95% on each is side is 1.96*sd 
comp.raw$m2.m1CI<-sqrt(comp.raw$m2.m1variance)*1.96
comp.raw$m3.m1CI<-sqrt(comp.raw$m3.m1variance)*1.96

ggplot(data=comp.raw,aes(x=m2.m1, y=m3.m1, fill=species)) +
  geom_ellipse(aes(x0=m2.m1, y0=m3.m1, a=m2.m1CI, b=m3.m1CI, angle=0), alpha=0.1) +
  coord_fixed() 

ggplot(data=comp.raw,aes(x=m2.m1CV, y=m3.m1CV)) +
  geom_point(aes(color=species)) #+
# geom_hline(yintercept = mean(comp.raw$m3.m1CV))
#other thoughts are convex hulls. What do you want to communicate?


## #format MMC comparative data -----
## #format fossils -----
## #create simulated composite molar rows --------
#create resampled pseudoreplicates of molars sizes and ratios from composite cheek teeth
sim.isolate.area<-matrix(NA,ncol=5,nrow=replicates) %>% as.data.frame #empty holder frame
colnames(sim.isolate.area)<-c("m1","m2","m3","m2.m1","m3.m1")
for (i in 1:nrow(sim.isolate.area)){ #resample each tooth position and take mean
  sim.isolate.area$m1[i]<-mouse.raw$m1 %>% sample(.,size=n,replace=TRUE) %>% mean
  sim.isolate.area$m2[i]<-mouse.raw$m2 %>% sample(.,size=n,replace=TRUE) %>% mean
  sim.isolate.area$m3[i]<-mouse.raw$m3 %>% sample(.,size=n,replace=TRUE) %>% mean
}
sim.isolate.area$m2.m1<-sim.isolate.area$m2/sim.isolate.area$m1
sim.isolate.area$m3.m1<-sim.isolate.area$m3/sim.isolate.area$m1

# sim.isolate.length<-matrix(NA,ncol=5,nrow=replicates) %>% as.data.frame #empty holder frame
# colnames(sim.isolate.length)<-c("m1","m2","m3","m2.m1","m3.m1")
# for (i in 1:nrow(sim.isolate.length)){ #resample each tooth position and take mean
#   sim.isolate.length$m1[i]<-mouse.raw$m1.length %>% sample(.,size=n,replace=TRUE) %>% mean
#   sim.isolate.length$m2[i]<-mouse.raw$m2.length %>% sample(.,size=n,replace=TRUE) %>% mean
#   sim.isolate.length$m3[i]<-mouse.raw$m3.length %>% sample(.,size=n,replace=TRUE) %>% mean
# }
# sim.isolate.length$m2.m1<-sim.isolate.length$m2/sim.isolate.length$m1
# sim.isolate.length$m3.m1<-sim.isolate.length$m3/sim.isolate.length$m1

## #calculate confidence intervals ------
raw.ci<-calc.ci(mouse.raw)
sim.a.ci<-calc.ci(sim.isolate.area)
sim.l.ci<-calc.ci(sim.isolate.length)

#write CI to file for later use
write.csv(sim.ci$ci95,"inhibitory_results/CI_gossypinus.csv",row.names=FALSE)

## #plot data ------

#make the two polygons of "forbidden space" in Kavanagh model
poly.right<-data.frame(id=rep("one",4),x=c(1,1,2.2,2.2),y=c(0,1,2,0))
poly.left<-data.frame(id=rep("two",4),x=c(0,1,1,0),y=c(0,1,2.2,2.2))

#plot to look at distribution. Structuring by site?
ggplot() + theme(legend.position=c(.7,.1)) +
  geom_point(data=mouse.raw,aes(x=m2.m1, y=m3.m1,fill=Institution),size=4,pch=21) 

#plot complete molar rows vs. simulated composite molar rows
cairo_pdf("inhibitory_results/gossypinus_ICM.pdf",width=4,height=4)
ggplot()+ 
  geom_polygon(data=poly.right,aes(x=x,y=y)) +
  geom_polygon(data=poly.left,aes(x=x,y=y)) +
  geom_point(data=sim.isolate.area,aes(x=m2.m1,y=m3.m1),alpha=0.5,fill="gray20",color="gray20",pch=21) +
  geom_point(data=mouse.raw,aes(x=m2.m1, y=m3.m1, pch=Institution),size=4,fill=hue_pal()(7)[1]) +
  scale_shape_manual(values=c(21,22,23))+
  coord_cartesian(xlim = c(0.7,0.95),ylim = c (0.5,0.75))+ 
  xlab("M/2:M/1")+ylab("M/3:M/1") +
  theme_minimal() +theme(legend.position = "none")
dev.off()
embedFonts("inhibitory_results/gossypinus_ICM.pdf")

#repeat, but with lengths instead of areas
ggplot()+ 
  # geom_polygon(data=poly.right,aes(x=x,y=y)) +
  # geom_polygon(data=poly.left,aes(x=x,y=y)) +
  geom_point(data=sim.isolate.length,aes(x=m2.m1,y=m3.m1),alpha=0.5,fill="gray20",color="gray20",pch=21) +
  geom_point(data=mouse.raw,aes(x=m2.m1L, y=m3.m1L, pch=Institution),size=4,fill=hue_pal()(7)[1]) +
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
  geom_point(data=mouse.raw,aes(x=m2.m1, y=m3.m1),size=4,pch=21,fill=hue_pal()(7)[1]) +
  coord_cartesian(xlim = c(0.6,1),ylim = c (0.5,0.8))+
  theme_minimal()

## #descriptive statistics ----

comp.ic<-comparative.ic.raw %>% group_by(.,species) %>% #group by species
  summarize(.,mean.m2.m1=mean(m2.m1),sd.m2.m1=mean(m2.m1variance), #get mean and sd for each species
            mean.m3.m1=mean(m3.m1),sd.m3.m1=mean(m3.m1variance)) %>%
  mutate(.,cv.m2.m1=sd.m2.m1/mean.m2.m1*100, #calculate cv
       cv.m3.m1=sd.m3.m1/mean.m3.m1*100)

#need mean, standard deviation, and coefficient of variation for both ratios, both cotton mouse sets
icm.cv<-function(sample){
  cv.m2.m1<-sd(sample$m2.m1)/mean(sample$m2.m1)*100  
  cv.m3.m1<-sd(sample$m3.m1)/mean(sample$m3.m1)*100
  return(c(cv.m2.m1,cv.m3.m1))
}

raw.gossypinus<-c("Peromyscus_gossypinus_complete",mean(mouse.raw$m2.m1),sd(mouse.raw$m2.m1),
                  mean(mouse.raw$m3.m1),sd(mouse.raw$m3.m1),icm.cv(mouse.raw))
sim.gossypinus<-c("Peromyscus_gossypinus_simulated",mean(sim.isolate$m2.m1),sd(sim.isolate$m2.m1),
                  mean(sim.isolate$m3.m1),sd(sim.isolate$m3.m1),icm.cv(sim.isolate))

all.ic<-rbind(comp.ic,raw.gossypinus,sim.gossypinus) %>% as.matrix
for (col in 2:ncol(all.ic)){all.ic[,col]<-as.numeric(all.ic[,col]) %>% round(.,2)}
write.csv(all.ic,"inhibitory_results/ICM_gossypinus_descriptive.csv",quote=FALSE)

## #compare means -----
t.test(mouse.raw$m2.m1, sim.isolate.area$m2.m1,alternative="two.sided",paired=FALSE,var.equal=FALSE) %>% 
  capture.output(.,file="inhibitory_results/gossypinus_ttest_m2m1.txt")
t.test(mouse.raw$m3.m1, sim.isolate.area$m3.m1,alternative="two.sided",paired=FALSE,var.equal=FALSE) %>% 
  capture.output(.,file="inhibitory_results/gossypinus_ttest_m3m1.txt")

t.test(mouse.raw$m2.m1L, sim.isolate.length$m2.m1,alternative="two.sided",paired=FALSE,var.equal=FALSE) %>% 
  capture.output(.,file="inhibitory_results/gossypinus_ttest_m2m1L.txt")
t.test(mouse.raw$m3.m1L, sim.isolate.length$m3.m1,alternative="two.sided",paired=FALSE,var.equal=FALSE) %>% 
  capture.output(.,file="inhibitory_results/gossypinus_ttest_m3m1L.txt")
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
