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

# reformat master mouse confidence intervals -----
mouse.stats$ID<-row.names(mouse.stats)
mouse.CI<-mouse.stats %>% select(-Mu,-Sigma,-CV) %>% melt(.,id="ID")
mouse.CI$direction<-"upper"
mouse.CI$direction[grepl("L",mouse.CI$variable)]<-"lower"
mouse.CI$metric<-"MMC2"
mouse.CI$metric[grepl("m3.m1L",mouse.CI$ID)]<-"MMC"
mouse.CI$metric[grepl("m3.m1A",mouse.CI$ID)]<-"m3m1"
mouse.CI$metric[grepl("m2.m1A",mouse.CI$ID)]<-"m2m1"
mouse.CI$type<-"mean"
mouse.CI$type[grepl("S",mouse.CI$variable)]<-"SD"
mouse.CI$label<-paste(mouse.CI$metric,mouse.CI$type,sep=".")
mouse.CI<-mouse.CI %>% select(direction,label,value) %>%
  dcast(.,label ~ direction)
row.names(mouse.CI)<-mouse.CI$label
CI<-select(mouse.CI,-label)

# evaluate resamples --------
#calculate summary stats for original mouse dataset to use as standar
for(ratio.choice in 1:nrow(ratio)){
  #create a vector of mean error
  vs.true<-data.frame(prop.mean.above.95ci=rep(NA,nrow(mouse)-1),
                      prop.mean.below.95ci=rep(NA,nrow(mouse)-1),
                      prop.sd.above.95ci=rep(NA,nrow(mouse)-1),
                      prop.sd.below.95ci=rep(NA,nrow(mouse)-1),
                      max.mean=rep(NA,nrow(mouse)-1),
                      min.mean=rep(NA,nrow(mouse)-1),
                      avg.sd=rep(NA,nrow(mouse)-1),
                      N=rep(2:nrow(mouse)))
  for (i in 2:nrow(mouse)){
    #little object of the sample of means
    resampled.means<-resampled.stats[which(resampled.stats$N==i),ratio[ratio.choice,2]] %>% unlist
    #calculate the proportion above and below the CI
    vs.true$prop.mean.above.95ci[i-1]<-which(resampled.means > (ratio[ratio.choice,2] %>% CI[.,2])) %>% 
      length(.)/replicates
    vs.true$prop.mean.below.95ci[i-1]<-which(resampled.means < (ratio[ratio.choice,2] %>% CI[.,1])) %>% 
      length(.)/replicates
    #get range of means
    vs.true$max.mean[i-1]<-max(resampled.means)
    vs.true$min.mean[i-1]<-min(resampled.means)
    resampled.sd<-resampled.stats[which(resampled.stats$N==i),ratio[ratio.choice,3]] %>% unlist
    vs.true$prop.sd.above.95ci[i-1]<-which(resampled.sd > (ratio[ratio.choice,3] %>% CI[.,2]) ) %>% 
      length(.)/replicates
    vs.true$prop.sd.below.95ci[i-1]<-which(resampled.sd < (ratio[ratio.choice,3] %>% CI[.,1])) %>% 
      length(.)/replicates
    vs.true$avg.sd[i-1]<-mean(resampled.sd)
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

#this tally hard to calculate in a single line because of changing ratio choices. Resorting to for-loop.
sample.size.indicators$mean.var.wrong<-NULL
for (i in 1:nrow(sample.size.indicators)){
  ratio.pick<-which(ratio[,1]==sample.size.indicators$ratio.name[i])
  sample.size.indicators$mean.var.wrong[i]<-(sample.size.indicators$avg.sd[i] > (ratio[ratio.pick,3] %>% CI[.,2]) )| 
    (sample.size.indicators$avg.sd[i] < (ratio[ratio.pick,3] %>% CI[.,1]) )
}


# reformat results for various purposes --------
#melt to to plot
resample.melted<-resampled.stats %>% melt(.,id=c("N"))
resample.melted$tooth<-"M3"
resample.melted$tooth[grepl("2",resample.melted$variable,perl=TRUE)]<-"M2"
resample.melted$dimension<-"Area"
resample.melted$dimension[grepl("MMC",resample.melted$variable,perl=TRUE)]<-"Length"
resample.melted$statistic<-"Mean"
resample.melted$statistic[grepl("SD",resample.melted$variable,perl=TRUE)]<-"Standard Deviation"
resample.melted$CI.upper<-apply(resample.melted,1, function(x) CI[x[2],"upper"])
resample.melted$CI.lower<-apply(resample.melted,1, function(x) CI[x[2],"lower"])
resample.melted$outside<-(resample.melted$CI.lower > resample.melted$value) |
                            (resample.melted$CI.upper < resample.melted$value)
 
#recast to write evaluation to a table
sample.size.indicators.to.write<-sample.size.indicators %>%
  select(N,ratio.name,outside.mean,outside.var,mean.var.wrong) %>%
  melt(.,id=c("N","ratio.name")) %>%
  dcast(.,N ~ ratio.name + variable)
write.csv(sample.size.indicators.to.write,
          paste("output/",species.name,"_sample_size_metrics.csv",sep=""),
          row.names = FALSE)

# plot --------
sample.size.plot<-ggplot(data=resample.melted, aes(x=N,y=value))+
  facet_grid(facets = statistic + tooth ~ dimension, scales="free_y", switch="y") +
  geom_point(alpha=0.25,aes(color=outside)) +
  geom_smooth(color="black",method="auto") + #make smoother, prettier?
  scale_color_manual(values=c("gray","red")) +
  ggtitle(species.name) + theme_minimal() + 
  theme(legend.position="none",axis.title.y=element_blank(),strip.placement = "outside",
        plot.title=element_text(face="italic"))
ggsave(sample.size.plot, filename = paste("output/",species.name,"_resample.pdf",sep=""), dpi = 300)