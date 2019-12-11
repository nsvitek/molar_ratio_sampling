# determine minimum sample size via resampling
# make pseudoreplicates, calculate summary stats from each -------

#Both: resample gossypinus dataset for all 4 ratios, check for sensitivity to sample size
#creates resampled.stats
resampled.base<-matrix(NA,nrow=replicates,ncol=length(resampled.names),
                       dimnames=list(c(NULL),resampled.names))
for (sample.n in 2:nrow(mouse)){
  resampled<-resampled.base
  resampled[,"N"]<-sample.n #sample size
  for(i in 1:replicates){
    mouse.resampled<-dplyr::sample_n(mouse,size=sample.n,replace=FALSE)
    if (any(str_detect(resampled.names,"MMC"))){
      resampled[i,"MMC.mean"]<-mean(mouse.resampled$m3.m1L)
      resampled[i,"MMC.SD"]<-sd(mouse.resampled$m3.m1L)
      # resampled[i,"MMC.max"]<-max(mouse.resampled$m3.m1L)
      # resampled[i,"MMC.min"]<-min(mouse.resampled$m3.m1L)
      resampled[i,"MMC2.mean"]<-mean(mouse.resampled$m2.m1L)
      resampled[i,"MMC2.SD"]<-sd(mouse.resampled$m2.m1L)
      # resampled[i,"MMC2.max"]<-max(mouse.resampled$m2.m1L)
      # resampled[i,"MMC2.min"]<-min(mouse.resampled$m2.m1L)
    }
    if (any(str_detect(resampled.names,"m1"))){
      resampled[i,"m3m1.mean"]<-mean(mouse.resampled$m3.m1A)
      resampled[i,"m3m1.SD"]<-sd(mouse.resampled$m3.m1A)
      # resampled[i,"m3m1.max"]<-max(mouse.resampled$m3.m1A)
      # resampled[i,"m3m1.min"]<-min(mouse.resampled$m3.m1A)
      resampled[i,"m2m1.mean"]<-mean(mouse.resampled$m2.m1A)
      resampled[i,"m2m1.SD"]<-sd(mouse.resampled$m2.m1A)
      # resampled[i,"m2m1.max"]<-max(mouse.resampled$m2.m1A)
      # resampled[i,"m2m1.min"]<-min(mouse.resampled$m2.m1A)
    }
  }
  if(sample.n==2){resampled.stats<-as.data.frame(resampled)}
  if(sample.n>=3){resampled.stats<-rbind(resampled.stats,resampled)}
}

# evaluate resamples --------
#get HPD interval for variables
HPD<-HPDinterval(as.mcmc(mouse.ratio.estimated.table))

#calculate summary stats for original mouse dataset to use as standar
for(ratio.choice in 1:nrow(ratio)){
  #create a vector of mean error
  vs.true<-data.frame(prop.mean.above.95ci=rep(NA,nrow(mouse)-1),
                      prop.mean.below.95ci=rep(NA,nrow(mouse)-1),
                      prop.sd.above.95ci=rep(NA,nrow(mouse)-1),
                      prop.sd.below.95ci=rep(NA,nrow(mouse)-1),
                      max.mean=rep(NA,nrow(mouse)-1),
                      min.mean=rep(NA,nrow(mouse)-1),
                      N=rep(2:nrow(mouse)))
  for (i in 2:nrow(mouse)){
    #little object of the sample of means
    resampled.means<-resampled.stats[which(resampled.stats$N==i),ratio[ratio.choice,2]] %>% unlist
    #calculate the proportion above and below the HPD
    vs.true$prop.mean.above.95ci[i-1]<-which(resampled.means > (ratio[ratio.choice,2] %>% HPD[.,2])) %>% 
      length(.)/replicates
    vs.true$prop.mean.below.95ci[i-1]<-which(resampled.means < (ratio[ratio.choice,2] %>% HPD[.,1])) %>% 
      length(.)/replicates
    #get range of means
    vs.true$max.mean[i-1]<-max(resampled.means)
    vs.true$min.mean[i-1]<-min(resampled.means)
    resampled.sd<-resampled.stats[which(resampled.stats$N==i),ratio[ratio.choice,3]] %>% unlist
    vs.true$prop.sd.above.95ci[i-1]<-which(resampled.sd > (ratio[ratio.choice,3] %>% HPD[.,2]) ) %>% 
      length(.)/replicates
    vs.true$prop.sd.below.95ci[i-1]<-which(resampled.sd < (ratio[ratio.choice,3] %>% HPD[.,1])) %>% 
      length(.)/replicates
  }
  vs.true$ratio.name<-ratio[ratio.choice,1]
  if(ratio.choice==1){ sample.size.indicators<-vs.true }
  if(ratio.choice>1){ sample.size.indicators<-rbind(sample.size.indicators,vs.true)}
}

#tally totals
sample.size.indicators$outside.mean<-(sample.size.indicators$prop.mean.above.95ci + 
                                        sample.size.indicators$prop.mean.below.95ci) %>%
  round(.,3)
sample.size.indicators$outside.var<-(sample.size.indicators$prop.sd.above.95ci + 
                                       sample.size.indicators$prop.sd.below.95ci) %>%
  round(.,3)

# reformat results for various purposes --------
#melt to to plot
resample.melted<-resampled.stats %>% melt(.,id=c("N"))
resample.melted$tooth<-"M3"
resample.melted$tooth[grepl("2",resample.melted$variable,perl=TRUE)]<-"M2"
resample.melted$dimension<-"Area"
resample.melted$dimension[grepl("MMC",resample.melted$variable,perl=TRUE)]<-"Length"
resample.melted$statistic<-"Mean"
resample.melted$statistic[grepl("SD",resample.melted$variable,perl=TRUE)]<-"Standard Deviation"
resample.melted$HPD.upper<-apply(resample.melted,1, function(x) HPD[x[2],2])
resample.melted$HPD.lower<-apply(resample.melted,1, function(x) HPD[x[2],1])
resample.melted$outside<-(resample.melted$HPD.lower > resample.melted$value) |
                            (resample.melted$HPD.upper < resample.melted$value)
 
cairo_pdf(paste("output/",species.name,"_resample.pdf",sep=""))
ggplot(data=resample.melted, aes(x=N,y=value))+
  facet_grid(facets = statistic + tooth ~ dimension, scales="free_y", switch="y") +
  geom_point(alpha=0.25,aes(color=outside)) +
  geom_smooth(color="black",method="auto") + #make smoother, prettier?
  scale_color_manual(values=c("gray","red")) +
  theme_minimal() + theme(legend.position="none",axis.title.y=element_blank(),strip.placement = "outside")
dev.off()

#melt sample size evaluation, this may not get used
ssi.long<-sample.size.indicators %>% select(N,ratio.name,outside.mean,outside.var) %>%
  melt(.,id=c("N","ratio.name"))

#recast to write evaluation to a table
sample.size.indicators.to.write<-sample.size.indicators %>%
  select(N,ratio.name,outside.mean,outside.var) %>%
  melt(.,id=c("N","ratio.name")) %>%
  dcast(.,N ~ ratio.name + variable)
write.csv(sample.size.indicators.to.write,
          paste("output/",species.name,"_sample_size_metrics.csv",sep=""),
          row.names = FALSE)
