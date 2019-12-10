# determine minimum sample size via resampling
# make pseudoreplicates, calculate summary stats from each -------
resampled.names<-c("N","MMC.mean","MMC.SD","MMC.max","MMC.min",
                   "MMC2.mean","MMC2.SD","MMC2.max","MMC2.min",
                   "m3m1.mean","m3m1.SD","m3m1.max","m3m1.min",
                   "m2m1.mean","m2m1.SD","m2m1.max","m2m1.min")

#Both: resample gossypinus dataset for all 4 ratios, check for sensitivity to sample size
#creates resampled.stats
resampled.base<-matrix(NA,nrow=replicates,ncol=length(resampled.names),
                       dimnames=list(c(NULL),resampled.names))
for (sample.n in 2:nrow(mouse)){
  resampled<-resampled.base
  resampled[,"N"]<-sample.n #sample size
  for(i in 1:replicates){
    mouse.resampled<-dplyr::sample_n(mouse,size=sample.n,replace=FALSE)
    resampled[i,"MMC.mean"]<-mean(mouse.resampled$m3.m1L)
    resampled[i,"MMC.SD"]<-sd(mouse.resampled$m3.m1L)
    resampled[i,"MMC.max"]<-max(mouse.resampled$m3.m1L)
    resampled[i,"MMC.min"]<-min(mouse.resampled$m3.m1L)
    
    resampled[i,"MMC2.mean"]<-mean(mouse.resampled$m2.m1L)
    resampled[i,"MMC2.SD"]<-sd(mouse.resampled$m2.m1L)
    resampled[i,"MMC2.max"]<-max(mouse.resampled$m2.m1L)
    resampled[i,"MMC2.min"]<-min(mouse.resampled$m2.m1L)
    
    resampled[i,"m3m1.mean"]<-mean(mouse.resampled$m3.m1A)
    resampled[i,"m3m1.SD"]<-sd(mouse.resampled$m3.m1A)
    resampled[i,"m3m1.max"]<-max(mouse.resampled$m3.m1A)
    resampled[i,"m3m1.min"]<-min(mouse.resampled$m3.m1A)
    
    resampled[i,"m2m1.mean"]<-mean(mouse.resampled$m2.m1A)
    resampled[i,"m2m1.SD"]<-sd(mouse.resampled$m2.m1A)
    resampled[i,"m2m1.max"]<-max(mouse.resampled$m2.m1A)
    resampled[i,"m2m1.min"]<-min(mouse.resampled$m2.m1A)
  }
  if(sample.n==2){resampled.stats<-as.data.frame(resampled)}
  if(sample.n>=3){resampled.stats<-rbind(resampled.stats,resampled)}
}

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
    #######COMPARE AGAINST HPD, NOT CONFIDENCE INTERVAL!!!!!!!
    #little object of the sample of means
    resampled.means<-resampled.stats[which(resampled.stats$N==i),ratio[ratio.choice,2]] %>% unlist
    #calculate the proportion above and below the HPD
    
    vs.true$prop.mean.above.95ci[i-1]<-which(resampled.means >= (ratio[1,2] %>% HPD[.,2])) %>% 
      length(.)/replicates
    vs.true$prop.mean.below.95ci[i-1]<-which(abs(resampled.means-mouse.stats$Mu[ratio.choice])>=(mouse.stats$Sigma[ratio.choice]*1.96)) %>% 
      length(.)/replicates
    
    #get range of means
    vs.true$max.mean[i-1]<-max(resampled.means)
    vs.true$min.mean[i-1]<-min(resampled.means)
    resampled.standarddevs<-resampled.stats[which(resampled.stats$N==i),ratio[ratio.choice,3]] %>% unlist
    vs.true$prop.sd.above.95ci[i-1]<-which((resampled.standarddevs-mouse.stats$var.CI.U[ratio.choice])>0) %>% 
      length(.)/replicates
    vs.true$prop.sd.below.95ci[i-1]<-which((resampled.standarddevs-mouse.stats$var.CI.L[ratio.choice])<0) %>% 
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
