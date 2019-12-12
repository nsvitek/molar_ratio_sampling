###################################################################
#Show that the ascending-descending (AD) model of Young et al. 2015 
#matches the 'ICM consistent' space of Polly 2007
##################################################################

# AD Model #######################################################
#visually match Young et al. 2015 Null Model 3 (ascending-descending) space
K<-seq(-2,-0.5,length.out=10) #b
J<-seq(1,0.5,length.out=10) #a
plot(x=0.5,xlim=c(0,1),ylim=c(0,1),type="n")
for(i in 1:length(K)){abline(a=J[i],b=K[i])}
J.pred<-(1/(-3))*K+(1/3)
round(J,3)==round(J.pred,3) #check

#So:
# m3 = K * m1 + J
#where
# J = (-1/3) * K + (1/3)
# -2 < K < -0.5
# 0.5 < J < 1

# ICM-C Model ###################################################
#visualy match Polly 2007 ICM-consistent space
B1<-seq(0,-9,length.out=10) #b, substituting 9 for infinity for visual purposes
Q1<-seq(1,9,length.out=10) #a
B2<-seq(0,-99999999999,length.out=10) #b, substituting 99999999 for infinity for visual purposes
Q2<-seq(1,99999999999,length.out=10) #a
plot(x=1,y=1,xlim=c(0,2),ylim=c(0,2))
for(i in 1:length(B)){
  abline(a=B1[i],b=Q1[i])
  abline(a=B2[i],b=Q2[i],col="blue")
  }
Q.pred<-1-B2
round(Q2/10,0)==round(Q.pred/10,0) #check, some fuzziness for infinity

#So:
# m3/m1 = Q * m2/m1 + B
#where
# Q = 1 - B
# 1 < Q < Inf
# -Inf < B < 0

#Note, the strict ICM line is where Q = 2 and B = -1. 
#according to covariance variance and covariance rules, matching Roseman and Delezene,
# var(m3) = Q^2 * var(m2) + B^2 * var(m1) - 2 * Q * B * Cov(m2, m1)

ICM_var<-function(m1,m3){
  # m3<-3*m2-(2*m1)
  m2<-(m3+(2*m1))/3
  total<-m1+m2+m3
  m2_relative<-m2/total
  return(m2_relative)
}
ICM<-function(m1,m3){
  # m3<-1.9*m2-(0.9*m1)
  m2<-(m3+(1.0*m1))/2
  total<-m1+m2+m3
  m2_relative<-m2/total
  return(m2_relative)
}
ICM_var(1,0.1)
ICM_var(1,0.4)
ICM_var(1,0.9)
ICM(1,0.4)
ICM(1,0.9)

# match variables ###############################################
#Axes of AD space are:
# X = m1/(m1 + m2 + m3)
# Y = m3/(m1 + m2 + m3)
#and
# 1 = m1 + m2 + m3
#simplifying to
# X = m1
# Y = m3

#Axes of ICM-C space are:
# X = m2/m1
# Y = m3/m1
#and if you continue to require that
# 1 = m1 + m2 + m3
#then you can change the X axis to
# X = (1-m1-m3)/m1
#which can be subsituted in so that
# m3/m1 = Q * (1-m1-m3)/m1 + B
#and rearranged so that
# m3 = m1 * (B-Q)/(1+Q) + Q/(1+Q)
#if AD space and ICM-C space are the same, and the molars have been set to be the same, then
#this rearranged form of ICM-C molar relationships should be equivalent to the
#equation relating m1 and m3 size in AD space, which would make:
# m1 * K  + J = m3 = m1 * (B-Q)/(1+Q) + Q/(1+Q)
# K = (B-Q)/(1+Q)
# J = Q/(1+Q)
#substitute those values into the relationship between J and K:
# J = (-1/3) * K + (1/3)
# Q/(1+Q) = (-1/3) * (B-Q)/(1+Q) + (1/3)
# and simplify to see if the ICM-C relationship between Q and B holds:
# 3 * Q/(1+Q) = (-1) * (B-Q)/(1+Q) + 1
# 3Q = -B + Q + 1 + Q
# Q = 1 - B
#Relationship holds!
