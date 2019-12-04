# ICM expectations ################################
#do length or area measures more closely match ICM predictions, implying that one is
#a preferred metric?
#ICM E1: The slope and intercept? ------------
#fit the gossypinus data to RMA linear regression model (1 spp, so no phylogenetic correction)
#actually, there are 2 questions:
#1) Is intraspecific variation linear and does that line match ICM? [cf. Bernal et al. 2013]
#Reduced Major Axis (RMA) aka Standard Major Axis regression:
#looked at Sokal and Rohlf (1995) section 14.13 for info about Model II regression
#reasoning for MII use is that when body X & Y have error, the slope estimates are biased
mouse.RMA.model.a<-lmodel2(m3.m1A ~ m2.m1A + 1, data=mouse,nperm=replicates)
mouse.RMA.model.l<-lmodel2(m3.m1L ~ m2.m1L + 1, data=mouse,nperm=replicates)

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


#2) Can you model m3 = 2*m2-m1 and get accurate m3 sizes? 
#Do this in E2 as part of larger Bayesian model

#also model variance of these ratios.
mouse.RMA.A.P<- cov(mouse[,c("m3.m1A","m2.m1A")],
                    use = "pairwise.complete.obs") %>% diag %>% diag
mouse.RMA.A.prior <- list(R = list(V = mouse.RMA.A.P,n=nrow(mouse.RMA.A.P)))

mouse.RMA.L.P<- cov(mouse[,c("m3.m1L","m2.m1L")],
                    use = "pairwise.complete.obs") %>% diag %>% diag
mouse.RMA.L.prior <- list(R = list(V = mouse.RMA.A.P,n=nrow(mouse.RMA.A.P)))


mouse.RMA.model.a.var<-MCMCglmm(cbind(m3.m1A,m2.m1A)~ trait - 1 ,
                                family = c("gaussian","gaussian"),
                                prior= mouse.RMA.A.prior, data = mouse,
                                rcov = ~us(trait):units,nitt = 1500000,burnin = 500000, 
                                thin = 1000)
mouse.RMA.model.l.var<-MCMCglmm(cbind(m3.m1L,m2.m1L)~ trait - 1 ,
                                family = c("gaussian","gaussian"),
                                prior= mouse.RMA.L.prior, data = mouse,
                                rcov = ~us(trait):units,nitt = 1500000,burnin = 500000, 
                                thin = 1000)

#put together the corresponding numbers that you actually get from the "true population"
mouse.ratio.estimated.table <- data.frame(mouse.RMA.model.a.var$VCV[,"traitm3.m1A:traitm3.m1A.units"]*avg.repeat["area"],
                                          mouse.RMA.model.a.var$VCV[,"traitm2.m1A:traitm2.m1A.units"]*avg.repeat["area"],
                                          mouse.RMA.model.l.var$VCV[,"traitm3.m1L:traitm3.m1L.units"]*avg.repeat["length"],
                                          mouse.RMA.model.l.var$VCV[,"traitm2.m1L:traitm2.m1L.units"]*avg.repeat["length"])
names(mouse.ratio.estimated.table) <- 	c("V.m3.m1A","V.m2.m1A","V.m3.m1L","V.m2.m1L")


#E2a: Roseman & Delezene's variance-covariance: Area -------------
#Delezene's variance-covariance prediction and Bayesian modelling that produced their fig. 3
#do via MCMCglmm following Roseman code

#make variance matrix for prior
mouse.A.P<- cov(mouse[,c("m1.area","m2.area","m3.area")],
                use = "pairwise.complete.obs") %>% diag %>% diag
mouse.A.prior <- list(R = list(V = mouse.A.P,n=nrow(mouse.A.P)))

#make a Bayesian model of the true underlying population and relationships among its variables
#most of this seems pretty standard. Note that "rcov" matches the help guide:
#"formula for residual covariance structure. This has to be set up so that 
#each data point is associated with a unique residual. For example a multi-response 
#model might have the R-structure defined by ~us(trait):units"
mouseA.model.test <- MCMCglmm(cbind(m1.area,m2.area,m3.area)~ trait - 1 ,
                              family = c("gaussian","gaussian","gaussian"),
                              prior= mouse.A.prior, data = mouse,
                              rcov = ~us(trait):units,nitt = 1500000,burnin = 500000, 
                              thin = 1000)

#only run next line if you want to pull it back up later. 
save(mouseA.model.test,file = "output/mouseA.model.test.RData")
#It uses the IA.dist.new() function.

#calculate the means, variances, covariances from the "true population"
#that the Bayesian model fit to the observed data
#note that variance is corrected for measurement error ("avg.repeat" object)

#calculate the numbers you chould get from applying ICM model expectations 
mouseA.model.table <- IA.dist.new(mouseA.model.test$VCV[,"traitm1.area:traitm1.area.units"]*avg.repeat["area"],
                                  mouseA.model.test$VCV[,"traitm2.area:traitm2.area.units"]*avg.repeat["area"],
                                  mouseA.model.test$Sol[,"traitm1.area"],
                                  mouseA.model.test$Sol[,"traitm2.area"], 
                                  mouseA.model.test$VCV[,"traitm1.area:traitm2.area.units"])

#put together the corresponding numbers that you actually get from the "true population"
mouseA.estimated.table <- data.frame(mouseA.model.test$Sol[,"traitm2.area"],
                                     mouseA.model.test$Sol[,"traitm3.area"],
                                     mouseA.model.test$VCV[,"traitm2.area:traitm2.area.units"]*avg.repeat["area"],
                                     mouseA.model.test$VCV[,"traitm3.area:traitm3.area.units"]*avg.repeat["area"],
                                     mouseA.model.test$VCV[,"traitm1.area:traitm2.area.units"],
                                     mouseA.model.test$VCV[,"traitm1.area:traitm3.area.units"],
                                     mouseA.model.test$VCV[,"traitm2.area:traitm3.area.units"])
names(mouseA.estimated.table) <- 	c("mu.m2","mu.m3","V.m2","V.m3",
                                    "C.m1.m2","C.m1.m3","C.m2.m3")
#E2b: Roseman & Delezene's variance-covariance: Length -------------
#Delezene's variance-covariance prediction and Bayesian modelling that produced their fig. 3
#do via MCMCglmm following Roseman code

#make variance matrix for prior
mouse.L.P<- cov(mouse[,c("m1.length","m2.length","m3.length")],
                use = "pairwise.complete.obs") %>% diag %>% diag
mouse.L.prior <- list(R = list(V = mouse.L.P,n=nrow(mouse.L.P)))

#make a Bayesian model of the true underlying population and relationships among its variables
#most of this seems pretty standard. Note that "rcov" matches the help guide:
#"formula for residual covariance structure. This has to be set up so that 
#each data point is associated with a unique residual. For example a multi-response 
#model might have the R-structure defined by ~us(trait):units"
mouseL.model.test <- MCMCglmm(cbind(m1.length,m2.length,m3.length)~ trait - 1 ,
                              family = c("gaussian","gaussian","gaussian"),
                              prior= mouse.L.prior, data = mouse,
                              rcov = ~us(trait):units,nitt = 1500000,burnin = 500000, 
                              thin = 1000)

#only run next line if you want to pull it back up later. 
save(mouseL.model.test,file = "output/mouseL.model.test.RData")

#calculate the means, variances, covariances from the "true population"
#that the Bayesian model fit to the observed data
#note that variance is corrected for measurement error ("avg.repeat" object)

#calculate the numbers you chould get from applying ICM model expectations 
mouseL.model.table <- IA.dist.new(mouseL.model.test$VCV[,"traitm1.length:traitm1.length.units"]*avg.repeat["length"],
                                  mouseL.model.test$VCV[,"traitm2.length:traitm2.length.units"]*avg.repeat["length"],
                                  mouseL.model.test$Sol[,"traitm1.length"],
                                  mouseL.model.test$Sol[,"traitm2.length"], 
                                  mouseL.model.test$VCV[,"traitm1.length:traitm2.length.units"])

#put together the corresponding numbers that you actually get from the "true population"
mouseL.estimated.table <- data.frame(mouseL.model.test$Sol[,"traitm2.length"],
                                     mouseL.model.test$Sol[,"traitm3.length"],
                                     mouseL.model.test$VCV[,"traitm2.length:traitm2.length.units"]*avg.repeat["length"],
                                     mouseL.model.test$VCV[,"traitm3.length:traitm3.length.units"]*avg.repeat["length"],
                                     mouseL.model.test$VCV[,"traitm1.length:traitm2.length.units"],
                                     mouseL.model.test$VCV[,"traitm1.length:traitm3.length.units"],
                                     mouseL.model.test$VCV[,"traitm2.length:traitm3.length.units"])
names(mouseL.estimated.table) <- 	c("mu.m2","mu.m3","V.m2","V.m3",
                                    "C.m1.m2","C.m1.m3","C.m2.m3")
#E3a: M2 is 1/3: area ------
#calculate mean and CI. is 1/3 in CI for proportional m2 size?
#can you also get Young's prediction that varaince of proportional m1 and 3 should be equal?
#do via MCMCglmm following Roseman code
mouse.A.R.P<- cov(mouse[,c("prop.m1A","prop.m2A","prop.m3A")],
                  use = "pairwise.complete.obs") %>% diag %>% diag
mouse.A.rel.prior <- list(R = list(V = mouse.A.R.P,n=nrow(mouse.A.R.P)))

#this model: (1) looks for sexual dimorphism in m2 proportions and
#(2) predicts mean of relative m2 size (the "+1")
proportionA.model.test <- MCMCglmm(cbind(prop.m1A,prop.m2A,prop.m3A)~ trait - 1 ,
                                   family = c("gaussian","gaussian","gaussian"),
                                   prior= mouse.A.rel.prior,data = mouse,
                                   rcov = ~us(trait):units, nitt = 1500000, burnin = 500000, 
                                   thin = 1000)

#only run next line if you want to pull it back up later. 
save(proportionA.model.test,file = "output/proportionA.model.test.RData")


#create expected and estimated relative m2 sizes
mouseA.model.table$mu.m2.prop<-1/3
mouseA.estimated.table$mu.m2.prop<-proportionA.model.test$Sol[,"traitprop.m2A"] %>% as.numeric

#create expected and estimated m3 model variances according to Young's prediction
mouseA.model.table$V.m3.prop<-(proportionA.model.test$VCV[,"traitprop.m1A:traitprop.m1A.units"]*
                                 avg.repeat["area"]) %>% as.numeric
mouseA.estimated.table$V.m3.prop<-(proportionA.model.test$VCV[,"traitprop.m3A:traitprop.m3A.units"]*
                                     avg.repeat["area"]) %>% as.numeric

# #compare to the kind of calculation you learned in intro stats:
# mean(mouse$prop.m2A)
# error <- qnorm(0.975)*sd(mouse$prop.m2A)/sqrt(length(mouse$prop.m2A))
# c(mean(mouse$prop.m2A)-error,mean(mouse$prop.m2A)+error) #doesn't provide exactly the same CI as the model
# hist(mouse$prop.m2A)
# abline(v=0.33333,col="red")
# abline(v=mean(mouse$prop.m2A)-error,col="gray",lty=2)
# abline(v=mean(mouse$prop.m2A)+error,col="gray",lty=2)
# #following Hlusko et al. 2016:
# plot(mouse$m2.area,mouse$total.area)
# abline(a=0,b=3)

#E3b: M2 is 1/3: length ------
#calculate mean and CI. is 1/3 in CI for proportional m2 size?
#can you also get Young's prediction that varaince of proportional m1 and 3 should be equal?
#do via MCMCglmm following Roseman code
mouse.L.R.P<- cov(mouse[,c("prop.m1L","prop.m2L","prop.m3L")],
                  use = "pairwise.complete.obs") %>% diag %>% diag
mouse.L.rel.prior <- list(R = list(V = mouse.L.R.P,n=nrow(mouse.L.R.P)))

#this model: (1) looks for sexual dimorphism in m2 proportions and
#(2) predicts mean of relative m2 size (the "+1")
proportionL.model.test <- MCMCglmm(cbind(prop.m1L,prop.m2L,prop.m3L)~ trait - 1 ,
                                   family = c("gaussian","gaussian","gaussian"),
                                   prior= mouse.L.rel.prior,data = mouse,
                                   rcov = ~us(trait):units, nitt = 1500000, burnin = 500000, 
                                   thin = 1000)
#only run next line if you want to pull it back up later. 
save(proportionL.model.test,file = "output/proportionL.model.test.RData")


#create expected and estimated relative m2 sizes
mouseL.model.table$mu.m2.prop<-1/3
mouseL.estimated.table$mu.m2.prop<-proportionL.model.test$Sol[,"traitprop.m2L"] %>% as.numeric

#create expected and estimated m3 model variances according to Young's prediction
mouseL.model.table$V.m3.prop<-(proportionL.model.test$VCV[,"traitprop.m1L:traitprop.m1L.units"]*
                                 avg.repeat["area"]) %>% as.numeric
mouseL.estimated.table$V.m3.prop<-(proportionL.model.test$VCV[,"traitprop.m3L:traitprop.m3L.units"]*
                                     avg.repeat["area"]) %>% as.numeric

#Put calculations on standard scales for plot:area  -------
#look at summary statistics. Is zero part of the probability density interval (= model and data match)?
HPDinterval(as.mcmc(mouseA.model.table-mouseA.estimated.table)) %>% round(.,5)

p1<-Transform.Expectations(mouseA.estimated.table,mouseA.model.table)

#make a Roseman & Delezene-style plot
pdf("output/mouse_ICM_predictions_area.pdf")
plot(1:7,NULL, type = "n", ylim = c(0,190),xlim = c(0.5,7.5),xlab = "",
     ylab = "Observed/Expected X 100", xaxt="n")

points(1:7,c(p1$M1.M3.cov.relative.delta[1],p1$M2.M3.cov.relative.delta[1], 
             p1$mstd.M3.relative.delta.summary[1],p1$M3.var.relative.delta[1],
             p1$M3.var.Young.relative.delta[1],p1$M3.means.relative.delta[1],
             p1$M2.prop.relative.delta[1]),pch = 19, cex = 3)

segments(1:7,c(p1$M1.M3.cov.relative.delta[2],p1$M2.M3.cov.relative.delta[2], 
               p1$mstd.M3.relative.delta.summary[2],p1$M3.var.relative.delta[2],
               p1$M3.var.Young.relative.delta[2],p1$M3.means.relative.delta[2],
               p1$M2.prop.relative.delta[2]),
         1:7,c(p1$M1.M3.cov.relative.delta[3],p1$M2.M3.cov.relative.delta[3], 
               p1$mstd.M3.relative.delta.summary[3],p1$M3.var.relative.delta[3],
               p1$M3.var.Young.relative.delta[3],p1$M3.means.relative.delta[3],
               p1$M2.prop.relative.delta[3]),lwd = 3)

#add expectations line
segments(0.175,100,7.3,100,lty = 3, lwd= 3, col = "darkgrey")

#add labels
text(2.5,175,"area")
text(1,10,expression(sigma(M[1],M[3])),cex=0.8)
text(2,10,expression(sigma(M[2],M[3])),cex=0.8)
text(3,10,expression(frac(sigma^2*(M[3]),mu^2)),cex=0.8)
text(4,10,expression(sigma^2*(M[3])),cex=0.8)
text(5,10,expression(sigma^2*(M[3-Rel])),cex=0.8)
text(6,10,expression(mu*(M[3])),cex=0.8)
text(7,10,expression(mu*M[2-Rel]),cex=0.8)
dev.off()
#plot: length -----------
#same but for lengths
HPDinterval(as.mcmc(mouseL.model.table-mouseL.estimated.table)) %>% round(.,5)

p2<-Transform.Expectations(mouseL.estimated.table,mouseL.model.table)

#make a Roseman & Delezene-style plot
pdf("output/mouse_ICM_predictions_length.pdf")
plot(1:7,NULL, type = "n", ylim = c(0,190),xlim = c(0.5,7.5),xlab = "",
     ylab = "Observed/Expected X 100", xaxt="n")

points(1:7,c(p2$M1.M3.cov.relative.delta[1],p2$M2.M3.cov.relative.delta[1], 
             p2$mstd.M3.relative.delta.summary[1],p2$M3.var.relative.delta[1],
             p2$M3.var.Young.relative.delta[1],p2$M3.means.relative.delta[1],
             p2$M2.prop.relative.delta[1]),pch = 19, cex = 3)

segments(1:7,c(p2$M1.M3.cov.relative.delta[2],p2$M2.M3.cov.relative.delta[2], 
               p2$mstd.M3.relative.delta.summary[2],p2$M3.var.relative.delta[2],
               p2$M3.var.Young.relative.delta[2],p2$M3.means.relative.delta[2],
               p2$M2.prop.relative.delta[2]),
         1:7,c(p2$M1.M3.cov.relative.delta[3],p2$M2.M3.cov.relative.delta[3], 
               p2$mstd.M3.relative.delta.summary[3],p2$M3.var.relative.delta[3],
               p2$M3.var.Young.relative.delta[3],p2$M3.means.relative.delta[3],
               p2$M2.prop.relative.delta[3]),lwd = 3)

#add expectations line
segments(0.175,100,7.3,100,lty = 3, lwd= 3, col = "darkgrey")

#add labels
text(2.5,175,"length")
text(1,10,expression(sigma(M[1],M[3])),cex=0.8)
text(2,10,expression(sigma(M[2],M[3])),cex=0.8)
text(3,10,expression(frac(sigma^2*(M[3]),mu^2)),cex=0.8)
text(4,10,expression(sigma^2*(M[3])),cex=0.8)
text(5,10,expression(sigma^2*(M[3-Rel])),cex=0.8)
text(6,10,expression(mu*(M[3])),cex=0.8)
text(7,10,expression(mu*(M[2-Rel])),cex=0.8)
dev.off()
# write tables ---------
tables.made<-list(mouseA.model.table=mouseA.model.table,
                  mouseA.estimated.table=mouseA.estimated.table,
                  mouseL.model.table=mouseL.model.table,
                  mouseL.estimated.table=mouseL.estimated.table)
for(i in 1:length(tables.made)){
  summaries<-cbind(posterior.mode(as.mcmc(tables.made[[i]])),HPDinterval(as.mcmc(tables.made[[i]])))
  write.csv(summaries,file=paste("output/",names(tables.made)[i],".csv",sep=""))
}
