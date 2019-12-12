###REDO THESE!
# CV in extant species -------
#CV of ratios by locality, Peromyscus gossypinus
####plot needs fixing up of axes, labels.
cairo_pdf("output/CV_geo_subsample.pdf")
ggplot(data=ratio.pop.CV.2, aes(x=variable,y=value)) +
  geom_point(data=filter(ratio.pop.CV.2,state=="Total"),color="black",size=3) +
  geom_point(data=filter(ratio.pop.CV.2,state!="Total"), aes(color=state),size=1.5) +
  theme_minimal()
dev.off()

#sense of ratio variation, Peromyscus gossypinus
####plot needs fixing up of axes, labels.
cairo_pdf("output/ICM_geo_gossy.pdf")
ggplot(mouse,aes(x=m2.m1A, y=m3.m1A, colour=state))+
  stat_chull(fill=NA)+
  geom_point(size=1.5)+
  theme_minimal()
dev.off()

#CV of ratios by sex, Primates
#####plots needs fixing up of axes, labels
cairo_pdf("output/CV_geo_subsample.pdf")
ggplot(data=ratio.sex.CV.2, aes(x=Species,y=value)) +
  facet_grid(rows = vars(variable),scales="free_y") +
  geom_point(data=filter(ratio.sex.CV.2,Sex=="Total"),color="black",size=3) +
  geom_point(data=filter(ratio.sex.CV.2,Sex!="Total"), aes(color=Sex),size=1.5) +
  theme_minimal()
dev.off()

#sense of how ratios vary by sex, Primates
#####plots needs fixing up of axes, labels
cairo_pdf("output/ICM_sex_primates.pdf")
ggplot(select(ICM.dimorph,Species,Sex,m3.m1A,m2.m1A) %>% melt(id=c("Species","Sex")),
  aes(x=Species, y=value, fill=Sex))+
  facet_grid(rows = vars(variable),scales="free_y") +
    geom_boxplot()
dev.off()

#show distribution of empirical CVs
#no longer a point in it being in regard to sample size. WHAT DO YOU WANT TO SHOW?
#####plots needs fixing up of axes, labels
cairo_pdf("output/CV_Species.pdf")
ggplot(data=CV.survey2, aes(x=M2,y=M3,color=Order,shape=Order))+
  facet_grid(rows = vars(dimension),scales="free_y") +
  geom_point()
dev.off()


# ICM expectations ------
pdf("output/mouse_ICM_RMA_area.pdf") 
plot(mouse.RMA.model.a,"SMA")
abline(a=mouse.RMA.model.a$regression.results$Intercept[3], 
       b=mouse.RMA.model.a$regression.results$Slope[3], col="red",lwd=2) #highlight RMA
abline(a=-1, b=2,col="black") #ICM slope
dev.off()

pdf("output/mouse_ICM_RMA_length.pdf")
plot(mouse.RMA.model.l,"SMA")
abline(a=mouse.RMA.model.l$regression.results$Intercept[3], 
       b=mouse.RMA.model.l$regression.results$Slope[3], col="red",lwd=2) #highlight RMA
abline(a=-1, b=2,col="black") #ICM slope
dev.off()

# composite molar rows -------
#plot complete molar rows vs. simulated composite molar rows

#make this at N=1,2,5,10?
ggplot()+ 
  geom_point(data=simulation.stats[which(simulation.stats$N==10),],
             aes(x=m2.m1A,y=m3.m1A),alpha=0.5,fill="gray20",color="gray20",pch=21) +
  geom_point(data=mouse,aes(x=m2.m1A, y=m3.m1A, pch=state),size=4,fill="darkred") +
  scale_shape_manual(values=c(21,22,23,24))+
  xlab("M2:M1")+ylab("M3:M1") +
  coord_cartesian(xlim = c(0.6,1.2),ylim = c (0.35,0.95))+
  theme_minimal() + theme(legend.position = "none")
