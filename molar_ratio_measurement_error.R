# functions, dependencies -------
library(dplyr)
library(reshape2)

PME<-function(ANOVA,r){
  # Yezerinac et al. 1992 p. 474 % measurement error
  s.within<-ANOVA$`Mean Sq`[2]
  s.among<-(ANOVA$`Mean Sq`[1]-s.within)/r
  percent_measurement_error<-(s.within/(s.within+s.among))
  return(percent_measurement_error)
}
# format error data ----
#mouse.raw from molar_ratio_base.R
in_df<-dplyr::select(mouse.raw,-c(institution, catalog_number, genus, species, filestart, level_fromfloor,
                           voxels, state, county, locality, sex))

#make some sample areas for each molar
in_df$m1.a1<-in_df$m1.l1*in_df$m1.w1
in_df$m1.a2<-in_df$m1.l2*in_df$m1.w2
in_df$m1.a3<-in_df$m1.l3*in_df$m1.w3

in_df$m2.a1<-in_df$m2.l1*in_df$m2.w1
in_df$m2.a2<-in_df$m2.l2*in_df$m2.w2
in_df$m2.a3<-in_df$m2.l3*in_df$m2.w3

in_df$m3.a1<-in_df$m3.l1*in_df$m3.w1
in_df$m3.a2<-in_df$m3.l2*in_df$m3.w2
in_df$m3.a3<-in_df$m3.l3*in_df$m3.w3

#step 1 in trying to make data "tidy"? 
in_df_long<-melt(in_df,id=c("specimen")) 

#unpack variables from label

#break up variable for faceting
in_df_long$tooth<-gsub("(m[1-3]).*","\\1",in_df_long$variable)
in_df_long$orientation<-in_df_long$variable
in_df_long$orientation<-gsub("m[1-3]\\.l[1-3]","length",in_df_long$orientation)
in_df_long$orientation<-gsub("m[1-3]\\.w[1-3]","width",in_df_long$orientation)
in_df_long$orientation<-gsub("m[1-3]\\.a[1-3]","area",in_df_long$orientation)
in_df_long$replicate<-gsub("m[1-3]\\.[lwa]([1-3])","\\1",in_df_long$variable)

# calculate repeatibility ----
#The Roseman & Delezene method for evaluating "within-species repeatibility" (p.233)
#is identical to Yezerinac's equation for percent measurement error, contained in function.

#r is number of repeated measurements per variable. Here, each measurement taken 3 times.
r<-3                                  

#make empty matrix to hold percents, rows are tooth positions, columns 
levels.row<-unique(in_df_long$tooth)
levels.col<-unique(in_df_long$orientation)
percents<-matrix(NA,nrow=length(levels.row),
                 ncol=length(levels.col),
                 dimnames=list(levels.row,levels.col))

#fill in each cell in the matrix with a nest of for-loops
for (row in levels.row){
  for (col in levels.col){
    ANOVA<-anova(lm(value~specimen, data = filter(in_df_long,
                                                  tooth==row,orientation==col)))
    percents[row,col]<-PME(ANOVA,r)
    
  }
}

#take the average repeatibility for each metric across the tooth row (avg of m1 thru m3)
avg.repeat<-apply(percents,2,function(x) 1-mean(x))

write.csv(avg.repeat,"output/measurement_repeatibility.csv")
#clean workspace
rm(list=c("ANOVA","in_df","in_df_long","col","levels.col","levels.row","r","row"))