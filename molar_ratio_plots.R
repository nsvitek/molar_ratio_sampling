###REDO THESE!
library(ggthemes) #to get Paul Tol colors
single.column.width<-3.46
double.column.width<-7
color.mouse<-ptol_pal()(4)
color.order<-ptol_pal()(length(unique(CV.survey2$Order)))
# CV in extant species -------
#CV of ratios by locality, Peromyscus gossypinus
cairo_pdf("output/CV_geo_subsample.pdf",width=single.column.width,height=single.column.width)
ggplot(data=ratio.pop.CV.2, aes(x=variable,y=value)) +
  geom_point(data=filter(ratio.pop.CV.2,state=="Total"),color="black",size=3) +
  geom_point(data=filter(ratio.pop.CV.2,state!="Total"), aes(color=state),size=2) +
  scale_color_manual(values=color.mouse) +
  xlab("Ratio")+ylab("CV") +
  theme_minimal()
dev.off()

#CV of ratios by sex, Primates
cairo_pdf("output/CV_sex_subsample.pdf",width=single.column.width,height=single.column.width)
ggplot(data=ratio.sex.CV.2, aes(x=Species,y=value)) +
  facet_grid(rows = vars(variable),scales="free_y") +
  geom_point(data=filter(ratio.sex.CV.2,Sex=="Total"),color="black",size=3) +
  geom_point(data=filter(ratio.sex.CV.2,Sex!="Total"), aes(color=Sex),size=1.5) +
  theme_minimal()
dev.off()


# #sense of ratio variation, Peromyscus gossypinus
# ####plot needs fixing up of axes, labels.
# cairo_pdf("output/ICM_geo_gossy.pdf")
# ggplot(mouse,aes(x=m2.m1A, y=m3.m1A, colour=state))+
#   stat_chull(fill=NA)+
#   geom_point(size=1.5)+
#   theme_minimal()
# dev.off()

# cairo_pdf("output/MMC_geo_gossy.pdf")
# ggplot(mouse,aes(x=m2.m1L, y=m3.m1L, colour=state))+
#   stat_chull(fill=NA)+
#   geom_point(size=1.5)+
#   theme_minimal()
# dev.off()

# #sense of how ratios vary by sex, Primates
# #####plots needs fixing up of axes, labels
# cairo_pdf("output/ICM_sex_primates.pdf")
# ggplot(select(ICM.dimorph,Species,Sex,m3.m1A,m2.m1A) %>% melt(id=c("Species","Sex")),
#   aes(x=Species, y=value, fill=Sex))+
#   facet_grid(rows = vars(variable),scales="free_y") +
#     geom_boxplot()  + theme_minimal()
# dev.off()

#show distribution of empirical CVs
#no longer a point in it being in regard to sample size. WHAT DO YOU WANT TO SHOW?
#####plots needs fixing up of axes, labels
cairo_pdf("output/FigX_CV_Species.pdf",width=double.column.width,height=single.column.width)
ggplot(data=CV.survey2, aes(x=M2,y=M3,color=Order,shape=Order))+
  facet_grid(rows = vars(dimension),scales="free_y") +
  scale_color_manual(values=color.order)+
  geom_point() + theme_minimal() +
  xlab(expression(CV(M[2]/M[1]))) +
  ylab(expression(CV(M[3]/M[1])))
dev.off()


# ICM expectations ------
####THESE TWO RMA PLOTS STILL NEED WORK
pdf("output/mouse_ICM_RMA_area.pdf") 
plot(mouse.RMA.model.a,"SMA",main=NULL,xlab=expr(M[2]:M[1]),ylab=expr(M[3]:M[1]),
     pch=21,bg="black")
abline(a=mouse.RMA.model.a$regression.results$Intercept[3], 
       b=mouse.RMA.model.a$regression.results$Slope[3], col="red",lwd=2) #highlight RMA
abline(a=-1, b=2,col="black") #ICM slope
dev.off()

pdf("output/mouse_ICM_RMA_length.pdf")
plot(mouse.RMA.model.l,"SMA",main=NULL,xlab=expr(M[2]:M[1]),ylab=expr(M[3]:M[1]),
     pch=21,bg="black")
abline(a=mouse.RMA.model.l$regression.results$Intercept[3], 
       b=mouse.RMA.model.l$regression.results$Slope[3], col="red",lwd=2) #highlight RMA
abline(a=-1, b=2,col="black") #ICM slope
dev.off()

# composite molar rows -------
#plot complete molar rows vs. simulated composite molar rows

#make this at N=1,2,5,10?
simstats2plot<-simulation.stats %>% select(.,N,m2m1.mean,m3m1.mean,MMC.mean,MMC2.mean) %>% 
  filter(.,N %in% c(1,2,5,10)) %>% melt(.,id="N") 
simstats2plot$tooth<-"M3"
simstats2plot$tooth[grepl("2",simstats2plot$variable,perl=TRUE)]<-"M2"
simstats2plot$dimension<-"Area"
simstats2plot$dimension[grepl("MMC",simstats2plot$variable,perl=TRUE)]<-"Length"
simstats2plot$ID<-rep(c(1:replicates),16)
simstats2plot2<-dcast(simstats2plot, ID+ dimension + N~ tooth)

mouse2plot<-select(mouse,specimen,m2.m1A,m3.m1A,m2.m1L,m3.m1L) %>% melt(.,id="specimen")
mouse2plot$tooth<-"M3"
mouse2plot$tooth[grepl("2",mouse2plot$variable,perl=TRUE)]<-"M2"
mouse2plot$dimension<-"Area"
mouse2plot$dimension[grepl("L",mouse2plot$variable,perl=TRUE)]<-"Length"
mouse2plot2<-dcast(mouse2plot, specimen + dimension ~ tooth)
mouse2plot3<-mouse2plot2[rep(seq_len(nrow(mouse2plot2)), 4), ]
mouse2plot3$N<-rep(c(1,2,5,10),each=nrow(mouse)*2)

composite.vs.variation<-ggplot()+ 
  geom_point(data=simstats2plot2,
             aes(x=M2,y=M3),alpha=0.5,color=ptol_pal()(2)[1],size=0.33) +
  geom_point(data=mouse2plot3,aes(x=M2,y=M3),size=1,color=ptol_pal()(2)[2]) +
  facet_grid(facets = N ~ dimension, scales="fixed") +
  xlab(expr(M[2]:M[1]))+ylab(expr(M[3]:M[1])) +
  theme_minimal() 
ggsave(composite.vs.variation, filename = paste("output/","FigX_composite_vs_variation.pdf",sep=""),
       dpi = fig.dpi,width=col.width,height=col.width*2,units=fig.units)

#show relationship between sample size and true values of mean and SD
composite.size.plot<-ggplot(data=sim.melted, aes(x=N,y=value))+
  facet_grid(facets = statistic + tooth ~ dimension, scales="free_y", switch="y") +
  geom_point(alpha=0.25,aes(color=outside)) +
  geom_smooth(color="black",method="auto") + #make smoother, prettier?
  scale_color_manual(values=c("gray","red")) +
  ggtitle(species.name) + theme_minimal() +
  theme(legend.position="none",axis.title.y=element_blank(),strip.placement = "outside",
        plot.title=element_text(face="italic"))
ggsave(composite.size.plot, filename = paste("output/","composite_sample_size.png",sep=""),
       dpi = fig.dpi,width=col.width,height=col.width*2,units=fig.units)

