#functions for molar ratio
#calculate a 95% normal confidence interval from a sample, or resample
calc.ci<-function(sample,X.variable,Y.variable){ #function to get elements of 95% CI circle
  #calculate mean of sample
  mean.xy<-c(mean(sample[,X.variable]), mean(sample[,Y.variable])) %>% 
    matrix(nrow=1,dimnames=list("mean",c(X.variable,Y.variable))) %>% data.frame
  #calculate distances between each subsample and mean, 
  distances<-((sample[,X.variable]-mean.xy[,X.variable])^2 +
                (sample[,Y.variable]-mean.xy[,Y.variable])^2) %>% 
    sqrt()
  #calculate  the 95% CI for those distances
  ci95<-quantile(distances,probs=c(0.95),na.rm=TRUE)
  #make into a data frame for ggplot
  ci<-cbind(mean.xy,ci95)
  return(ci)
}

#function is C. Roseman's function from Roseman and Delezene 2019
IA.dist.new <- function(V.m1, V.m2 ,mu.m1, mu.m2,C.m1.m2){
  #calculate expected m3 area based on size of m1 and m2 and the ICM equation of Kavanagh
  mu.m3 <- 2 * mu.m2 - mu.m1
  #expected variance of m3 based on Roseman and Delezen 2019 eq. 1, involving rule:
  #cov(X+Y) = var(X) + var(Y) + 2*cov(X,Y)
  V.m3 <- 4*V.m2 + V.m1 - 4*C.m1.m2 
  #Roseman and Delezene 2019 eq. 2
  C.m1.m3 <- 2*C.m1.m2 - V.m1
  #Roseman and Delezene 2019 eq. 3
  C.m2.m3 <- 2*V.m2 - C.m1.m2
  #put everything in one table
  out <- data.frame(mu.m2,mu.m3,V.m2,V.m3,C.m1.m2,C.m1.m3,C.m2.m3)
  #label the table
  names(out) <- 	c("mu.m2","mu.m3","V.m2","V.m3","C.m1.m2","C.m1.m3","C.m2.m3")
  return(out)
}

#transforms results to a relative scale. 
Transform.Expectations<-function(estimated.table,model.table){
  relative.delta <- cbind(posterior.mode(as.mcmc((estimated.table)/model.table)),
                          HPDinterval(as.mcmc((estimated.table)/model.table)))
  #Mean standardized results broken into individual vectors
  M3.means.relative.delta <-relative.delta["mu.m3",]
  M3.var.relative.delta <- relative.delta["V.m3",]
  M1.M3.cov.relative.delta <- relative.delta["C.m1.m3",]
  M2.M3.cov.relative.delta <- relative.delta["C.m2.m3",]
  #calculate this scaled-variance variable.
  #as label in plot explains, it's the variance [v]
  #over the mean-squared [mu, squared],
  #but this is all scaled by estimated/model results
  mstd.M3.relative.delta <-((estimated.table$V.m3/estimated.table$mu.m3^2)/
                              (model.table$V.m3/model.table$mu.m3^2))
  #this is the fourth variable to be plotted, mode and interval
  mstd.M3.relative.delta.summary <- cbind(posterior.mode(as.mcmc(mstd.M3.relative.delta)),
                                          (HPDinterval(as.mcmc(mstd.M3.relative.delta))))*100
  #add in the relative proportion expectations
  M2.prop.relative.delta<-relative.delta["mu.m2.prop",]
  M3.var.Young.relative.delta<-relative.delta["V.m3.prop",]
  return(list(M3.means.relative.delta=M3.means.relative.delta,
              M3.var.relative.delta=M3.var.relative.delta,
              M1.M3.cov.relative.delta=M1.M3.cov.relative.delta,
              M2.M3.cov.relative.delta=M2.M3.cov.relative.delta,
              mstd.M3.relative.delta.summary=mstd.M3.relative.delta.summary,
              M2.prop.relative.delta=M2.prop.relative.delta,
              M3.var.Young.relative.delta=M3.var.Young.relative.delta))
}

##pairwise U-test
Pairwise.U<-function(measurement,group.var){
  pairs.state<-combn(unique(group.var),m=2) %>% t %>% as.data.frame
  pairs.state$U.p.vals<-NULL
  for (i in 1:ncol(pairs.state)){
    pairs.state$U.p.vals[i]<-wilcox.test(measurement[which(group.var==pairs.state[i,1])],
                                         measurement[which(group.var==pairs.state[i,2])])$p.value
  }
  return(pairs.state)
}

#bootstrap ANOVA approach, these are mostly vectors. pop.sizes is a tibble from dplyr, replicates an integer
BootstrapANOVA<-function(measurement,group.var,pop.sizes,replicates){
  observed.value<-aov(lm(measurement~group.var)) %>% summary
  resam.distribution<-observed.value[[1]][1,4]
  resam.var<-rep(NA,length(group.var))
  group.var<-as.data.frame(group.var) #transform it so you can add columns
  for(j in 1:replicates){
    for(i in 1:nrow(pop.sizes)){
      resam.var[which(group.var==as.character(pop.sizes[i,1]))]<-sample(measurement, 
                                                                        size=as.integer(pop.sizes[i,2]), replace=TRUE)
    }
    group.var$dependent<-resam.var
    anova.resample<-aov(lm(dependent~.,data=group.var)) %>% summary
    statistic<-anova.resample[[1]]$`F value`[1]
    resam.distribution<-c(resam.distribution,statistic)
  }
  return(resam.distribution)
}

#convex hulls in ggplot
#from: https://cran.r-project.org/web/packages/ggplot2/vignettes/extending-ggplot2.html
StatChull <- ggproto("StatChull", Stat,
                     compute_group = function(data, scales) {
                       data[chull(data$x, data$y), , drop = FALSE]
                     },
                     
                     required_aes = c("x", "y")
)

stat_chull <- function(mapping = NULL, data = NULL, geom = "polygon",
                       position = "identity", na.rm = FALSE, show.legend = NA, 
                       inherit.aes = TRUE, ...) {
  ggplot2::layer(
    stat = StatChull, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}