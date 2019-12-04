###REDO THESE!


#show distribution of empirical CVs with regard to sample size. 3 part figure?

#for sample size recs for "true" populations, plot how many resamples are outside 95% CI of true mean & sd

#for composite simulations, plot mean of the 10k"fossil samples" and point estimates for each one 
#showing distribution with regards to the "true" value on a horizontal line for mean, plus 95%CI around true mean

#for composite simulations + "pseudopop", plot the SD of each "fossil sample" for each N
#but the horizontal line is the "true" SD of the N=70 complete molar row sample

#for covariance, center the "true" posterior distribution mode on zero and adjust the 2 HPD values accordingly
#For each of the 2 covariances for Area and Length [that's 6 graphs]
#and those are your horizontal lines (each measure) then plot the 10k  deviations from the "true" covariance
#per N

#plot
ggplot(data=ssi.long,aes(x=N,y= value)) +
  geom_point() +
  geom_point(data=ssi.long[which(ssi.long$value<0.05),],col="red")+
  facet_grid(facets = variable ~ ratio.name,scales="free_y") #+
# scale_x_discrete("Part of PETM") +
# scale_y_continuous("Thickness (mm)") +
# scale_color_ptol(guide=FALSE) + scale_shape(guide=FALSE)+
# ggtitle("M/2 Trigonid") +
# theme_few() + theme(axis.text.x = rotatedAxisElementText(0,"x"))

ggplot(data=resampled.stats, aes_string(x="N",y="m2m1.SD"))+
  geom_point(alpha=0.25,color="gray") +
  geom_smooth(color="black",method="loess") + #make smoother, prettier? or quantile?
  geom_hline(yintercept=var.CI[1],linetype="dashed",color="blue",size=1)+
  geom_hline(yintercept=var.CI[2],linetype="dashed",color="blue",size=1)+
  geom_point(data=resampled.stats[which(resampled.stats$m2m1.SD-var.CI[2]>0),],col="red")



#plot MMC CV, how gossypinus relates to sample size and to estimates from samples of other spp
cairo_pdf("output/MMC_species_CV.pdf")
ggplot(data=resampled.stats, aes(x=N,y=MMC.CV))+
  geom_point(alpha=0.25,color="gray") +
  geom_smooth(color="black",method="loess") + #make smoother, prettier? or quantile?
  geom_hline(yintercept=mean(MMC.CV$MMC.CV),linetype="dashed",color="black",size=1.5)+
  geom_point(data=MMC.CV,aes(x=N,y=MMC.CV,color=Order),size=1.5) +
  xlim(2,max(MMC.CV$N)) + theme_minimal()
dev.off()

#same but for MMC mean
cairo_pdf("output/MMC_species_mean.pdf")
ggplot(data=resampled.stats, aes(x=N,y=MMC.mean))+
  geom_point(alpha=0.25,color="gray") +
  geom_smooth(color="black",method="loess") + #make smoother, prettier?
  geom_point(data=MMC.CV,aes(x=N,y=MMC.mean,color=Order),size=1.5) +
  xlim(2,max(MMC.CV$N)) + theme_minimal()
dev.off()

#same but for MMC SD
cairo_pdf("output/MMC_species_SD.pdf")
ggplot(data=resampled.stats, aes(x=N,y=MMC.SD))+
  geom_point(alpha=0.25,color="gray") +
  geom_smooth(color="black",method="loess") + #make smoother, prettier?
  geom_point(data=MMC.CV,aes(x=N,y=MMC.SD,color=Order),size=1.5) +
  xlim(2,max(MMC.CV$N)) + theme_minimal()
dev.off()

#ICM: calculate mean by population/species in similar way
cairo_pdf("output/ICMm3m1_species_mean.pdf")
ggplot(data=resampled.stats, aes(x=N,y=m3m1.mean))+
  geom_point(alpha=0.25,color="gray") +
  geom_smooth(color="black",method="loess") + #make smoother, prettier?
  geom_point(data=compICM.pop,aes(x=N,y=m3m1.mean,color=Species),size=1.5) +
  xlim(2,max(compICM.pop$N)) + theme_minimal()
dev.off()

cairo_pdf("output/ICMm2m1_species_mean.pdf")
ggplot(data=resampled.stats, aes(x=N,y=m2m1.mean))+
  geom_point(alpha=0.25,color="gray") +
  geom_smooth(color="black",method="loess") + #make smoother, prettier?
  geom_point(data=compICM.pop,aes(x=N,y=m2m1.mean,color=Species),size=1.5) +
  xlim(2,max(compICM.pop$N)) + theme_minimal()
dev.off()

cairo_pdf("output/ICMm3m1_species_CV.pdf")
ggplot(data=resampled.stats, aes(x=N,y=m3m1.CV))+
  geom_point(alpha=0.25,color="gray") +
  geom_smooth(color="black",method="loess") + #make smoother, prettier?
  geom_point(data=compICM.pop,aes(x=N,y=m3m1.CV,color=Species),size=1.5) +
  xlim(2,max(compICM.pop$N)) + theme_minimal()
dev.off()

cairo_pdf("output/ICMm2m1_species_CV.pdf")
ggplot(data=resampled.stats, aes(x=N,y=m2m1.CV))+
  geom_point(alpha=0.25,color="gray") +
  geom_smooth(color="black",method="loess") + #make smoother, prettier?
  geom_point(data=compICM.pop,aes(x=N,y=m2m1.CV,color=Species),size=1.5) +
  xlim(2,max(compICM.pop$N)) + theme_minimal()
dev.off()

cairo_pdf("output/ICMm3m1_species_SD.pdf")
ggplot(data=resampled.stats, aes(x=N,y=m3m1.SD))+
  geom_point(alpha=0.25,color="gray") +
  geom_smooth(color="black",method="loess") + #make smoother, prettier?
  geom_point(data=compICM.pop,aes(x=N,y=m3m1.SD,color=Species),size=1.5) +
  xlim(2,max(compICM.pop$N)) + theme_minimal()
dev.off()

cairo_pdf("output/ICMm2m1_species_SD.pdf")
ggplot(data=resampled.stats, aes(x=N,y=m2m1.SD))+
  geom_point(alpha=0.25,color="gray") +
  geom_smooth(color="black",method="loess") + #make smoother, prettier?
  geom_point(data=compICM.pop,aes(x=N,y=m2m1.SD,color=Species),size=1.5) +
  xlim(2,max(compICM.pop$N)) + theme_minimal()
dev.off()

#Old
# #check to make sure things look right 
# ggplot(data=simulaton.stats[which(simulation.stats$N==50)],aes(x=m2.m1A, y=m3.m1A)) +
#   geom_ellipse(aes(x0=m2.m1A, y0=m3.m1A, a=m2.m1CI, b=m3.m1CI, angle=0), alpha=0.1) +
#   coord_fixed()
# 
# ggplot(data=comp.raw,aes(x=m2.m1CV, y=m3.m1CV)) +
#   geom_point(aes(color=species)) #+
# # geom_hline(yintercept = mean(comp.raw$m3.m1CV))
# #other thoughts are convex hulls. What do you want to communicate?
