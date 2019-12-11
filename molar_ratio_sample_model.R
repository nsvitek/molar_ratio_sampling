# make Bayesian model -----------
avg.repeat<-mean(c(0.9313, 0.9586)) #repeatabilities for Peromyscus
#also model variance of these ratios.
mouse.RMA.L.P<- cov(mouse[,c("m3.m1L","m2.m1L")],
                    use = "pairwise.complete.obs") %>% diag %>% diag
mouse.RMA.L.prior <- list(R = list(V = mouse.RMA.L.P,n=nrow(mouse.RMA.L.P)))

mouse.RMA.model.l.var<-MCMCglmm(cbind(m3.m1L,m2.m1L)~ trait - 1 ,
                                family = c("gaussian","gaussian"),
                                prior= mouse.RMA.L.prior, data = mouse,
                                rcov = ~us(trait):units,nitt = 1500000,burnin = 500000, 
                                thin = 1000)

mouse.ratio.estimated.table <- data.frame(mouse.RMA.model.l.var$Sol[,"traitm3.m1L"],
                                          mouse.RMA.model.l.var$Sol[,"traitm2.m1L"],
                                          sqrt(mouse.RMA.model.l.var$VCV[,"traitm3.m1L:traitm3.m1L.units"]*avg.repeat),
                                          sqrt(mouse.RMA.model.l.var$VCV[,"traitm2.m1L:traitm2.m1L.units"]*avg.repeat))
names(mouse.ratio.estimated.table) <- 	c("MMC.mean","MMC2.mean",
                                         "MMC.SD","MMC2.SD")
