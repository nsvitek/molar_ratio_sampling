#some basic lambda and K to estimate phylogenetic signal of a trait on a basic, easily accessible mammal tree. 

library(phytools)
library(treeman)

data(mammals)
str(mammals)
MMC.CV$Species
MMC.CV$MatchSpp<-gsub("([A-Za-z]*) ([A-Za-z]*)","\\1_\\2",MMC.CV$Species,perl=TRUE)

#prune the tree
mammalsMMC<-rmTips(mammals, tids = mammals@tips[-which(mammals@tips %in% MMC.CV$MatchSpp)]) %>%
  as(.,"phylo")
plot(mammalsMMC) #chec: a reasonable tree?
rm(rm="mammals") #is an 8mb dataset, remove when finished with it
#create variable for testing phylogenetic signal
phylo.MMC<-MMC.CV$MMC.CV[which(MMC.CV$MatchSpp %in% mammalsMMC$tip.label)]
names(phylo.MMC)<-MMC.CV$MatchSpp[which(MMC.CV$MatchSpp %in% mammalsMMC$tip.label)]
#test
CV.K<-phylosig(mammalsMMC, phylo.MMC, method="K", test=TRUE, nsim=1000)
CV.L<-phylosig(mammalsMMC, phylo.MMC, method="lambda", test=TRUE, nsim=1000)

#likelihood ratio test
fitBrownian<-brownie.lite(paintSubTree(mammalsMMC,mammalsMMC$edge[1,1],state="1"),phylo.MMC)
lam1<-phylosig(mammalsMMC,x1,method="lambda")
LR<-2*(CV.L$logL-fitBrownian$logL1)
pchisq(LR,df=1,lower.tail=FALSE) #significantly different from "phylogenetic signal", or pure Brownian motion
