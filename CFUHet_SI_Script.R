#                Code for Supplementary Information
#
######################################################################
# Function wormbootGMM
wormbootGMM<-function(reps, batches, meanD, varD, meanU, varU, fUP, fDIFF){
  # Script for generating log-scale fake data based on a Gaussian mixture model
  # Prints p-values for t-tests and Mann-Whitney U tests in order of batch size
  # Returns an example data frame of simulated data
  # where "batch" is the number of individual measurements in the batch digest
  # "Count" is CFU/worm, and "logCount" is the original log-scale simulated data
  #"fUP" should be the larger of the fractions [0,1] of high-mode individuals
  # and "fDIFF" is the difference in modal fraction across the two populations
  fUPb<-fUP-fDIFF
  mySDD<-sqrt(varD)
  mySDU<-sqrt(varU)
  pv1<-numeric(reps)
  pv5<-numeric(reps)
  pv10<-numeric(reps)
  pv20<-numeric(reps)
  pv50<-numeric(reps)
  wpv1<-numeric(reps)
  wpv5<-numeric(reps)
  wpv10<-numeric(reps)
  wpv20<-numeric(reps)
  wpv50<-numeric(reps)
  
  for (j in 1:reps){
    up<-sum(rbinom(batches,1, fUP))
    down=batches-up
    tempA<-c(rnorm(down, mean=meanD, sd=mySDD), rnorm(up, mean=meanU, sd=mySDU))
    tempA[tempA<1.3]<-0 
    up<-sum(rbinom(batches,1, fUPb))
    down=batches-up
    tempB<-c(rnorm(down, mean=meanD, sd=mySDD), rnorm(up, mean=meanU, sd=mySDU))
    tempB[tempB<1.3]<-0 
    # batch 5
    temp5A<-numeric(batches)
    temp5B<-numeric(batches)
    temp10A<-numeric(batches)
    temp10B<-numeric(batches)
    temp20A<-numeric(batches)
    temp20B<-numeric(batches)
    temp50A<-numeric(batches)
    temp50B<-numeric(batches)
    
    for(k in 1:batches){
      up<-sum(rbinom(5, 1, fUP))
      down=5-up
      temp1<-c(rnorm(down, mean=meanD, sd=mySDD), rnorm(up, mean=meanU, sd=mySDU))
      up<-sum(rbinom(5, 1, fUPb))
      down=5-up
      temp2<-c(rnorm(down, mean=meanD, sd=mySDD), rnorm(up, mean=meanU, sd=mySDU))
      temp1[temp1<1.3]<-0
      temp2[temp2<1.3]<-0
      temp5A[k]<-log10(mean(10^temp1))
      temp5B[k]<-log10(mean(10^temp2))
      #now batch 10
      up<-sum(rbinom(10, 1, fUP))
      down=10-up
      temp1<-c(rnorm(down, mean=meanD, sd=mySDD), rnorm(up, mean=meanU, sd=mySDU))
      up<-sum(rbinom(10, 1, fUPb))
      down=10-up
      temp2<-c(rnorm(down, mean=meanD, sd=mySDD), rnorm(up, mean=meanU, sd=mySDU))
      temp1[temp1<1.3]<-0
      temp2[temp2<1.3]<-0
      temp10A[k]<-log10(mean(10^temp1))
      temp10B[k]<-log10(mean(10^temp2))
      # batch 20
      up<-sum(rbinom(20, 1, fUP))
      down=20-up
      temp1<-c(rnorm(down, mean=meanD, sd=mySDD), rnorm(up, mean=meanU, sd=mySDU))
      up<-sum(rbinom(20, 1, fUPb))
      down=20-up
      temp2<-c(rnorm(down, mean=meanD, sd=mySDD), rnorm(up, mean=meanU, sd=mySDU))
      temp1[temp1<1.3]<-0
      temp2[temp2<1.3]<-0
      temp20A[k]<-log10(mean(10^temp1))
      temp20B[k]<-log10(mean(10^temp2))
      #batch50
      up<-sum(rbinom(50, 1, fUP))
      down=50-up
      temp1<-c(rnorm(down, mean=meanD, sd=mySDD), rnorm(up, mean=meanU, sd=mySDU))
      up<-sum(rbinom(50, 1, fUPb))
      down=50-up
      temp2<-c(rnorm(down, mean=meanD, sd=mySDD), rnorm(up, mean=meanU, sd=mySDU))
      temp1[temp1<1.3]<-0
      temp2[temp2<1.3]<-0
      temp50A[k]<-log10(mean(10^temp1))
      temp50B[k]<-log10(mean(10^temp2))
    }
    
    #t tests on log data
    test1<-t.test(tempA, tempB)
    test5<-t.test(temp5A, temp5B)
    test10<-t.test(temp10A, temp10B)
    test20<-t.test(temp20A, temp20B)
    test50<-t.test(temp50A, temp50B)
    
    #Wilcoxon rank sum tests
    wtest1<-wilcox.test(tempA,tempB)
    wtest5<-wilcox.test(temp5A,temp5B)
    wtest10<-wilcox.test(temp10A,temp10B)
    wtest20<-wilcox.test(temp20A,temp20B)
    wtest50<-wilcox.test(temp50A,temp50B)
    
    pv1[j]<-test1$p.value
    pv5[j]<-test5$p.value
    pv10[j]<-test10$p.value
    pv20[j]<-test20$p.value
    pv50[j]<-test50$p.value
    
    wpv1[j]<-wtest1$p.value
    wpv5[j]<-wtest5$p.value
    wpv10[j]<-wtest10$p.value
    wpv20[j]<-wtest20$p.value
    wpv50[j]<-wtest50$p.value
  }
  
  pt<-c(sum(pv1<0.05)/reps, sum(pv5<0.05)/reps, sum(pv10<0.05)/reps, sum(pv20<0.05)/reps, sum(pv50<0.05)/reps)
  pw<-c(sum(wpv1<0.05)/reps, sum(wpv5<0.05)/reps, sum(wpv10<0.05)/reps, sum(wpv20<0.05)/reps, sum(wpv50<0.05)/reps)
  print(pt)
  print(pw)
  batch<-c(rep(1,batches), rep(5, batches), rep(10, batches), rep(20, batches), rep(50, batches))
  dataA<-c(tempA, temp5A, temp10A, temp20A, temp50A)
  dataB<-c(tempB, temp5B, temp10B, temp20B, temp50B)
  frameA<-data.frame(batch, logCount=dataA)
  frameA$set<-"1"
  frameB<-data.frame(batch, logCount=dataB)
  frameB$set<-"2"
  mydata<-rbind(frameA, frameB)
  mydata$Count<-10^mydata$logCount
  return(mydata)
}

#########################################################################################
#########################################################################################
# function simBetaCFU()
simBetaCFU<-function(a, b, nreps, K=100000){
  # Function for generating simulated CFU-like data from a family of beta distributions
  # The input parameters (a, b) are used for the parent distribution
  # In daughter distributions, both parameters are multiplied by a constant 
  # from myfactors=(2.5, 5, 10, 20, 50) by default
  # such that the mean of all distributions is the same, but higher moments change.
  # 
  # Parameters:
  # a: alpha for the initial beta(alpha, beta)
  # b: beta for the initial beta(alpha, beta)
  # nreps: Number of values to create from each distribution (same for all distributions),
  # K: Rescaling constant. Used to bring the products of rbeta into CFU/worm range (default 10e5)
  myfactors<-c(1, 2.5, 5, 10, 20, 50)
  mydata<-data.frame(matrix(vector(), nreps*length(myfactors), 2,
                            dimnames=list(c(), c("Factor", "Data"))),
                     stringsAsFactors=F)
  mydata$Factor<-c(rep(1, nreps), 
                   rep(2.5, nreps), 
                   rep(5, nreps), 
                   rep(10, nreps), 
                   rep(20, nreps), 
                   rep(50,nreps))
  temp<-vector()
  for (i in 1:length(myfactors)){
    ai<-a*myfactors[i]
    bi<-b*myfactors[i]
    temp<-c(temp,rbeta(nreps, ai, bi)*K)
  }
  mydata$Data<-temp
  mydata$Data[mydata$Data<1]<-0
  mydata$LogData<-log10(mydata$Data+1)
  mydata$alpha<-a*mydata$Factor
  mydata$beta<-b*mydata$Factor
  mydata$Factor<-as.factor(mydata$Factor)
  return(mydata)
}


#########################################################################################
#########################################################################################
# function wormSimBatchBetaCompare
wormSimBatchBetaCompare<-function(a1, b1, a2, b2, runs, batches, maxCFU, prange=0){
  #Function that simulates worm CFU data based on draws from beta distribution
  # and compares simulated distributions using t-tests & Wilcoxon
  #
  # If prange is specified:
  # Parameters for data sets A and B are re-drawn separately for each run with errors U(-prange, prange)
  #
  #a and b are underlying parameters of the beta distribution (double)
  #runs is the number of runs in each simulation (integer)
  #batches is the number of data points to generate for each run (integer)
  #maxCFU is the maximum value you want for the CFU/worm synthetic data (number)
  #return a data frame containing moments of data from each rep
  #column "batch" is the batch size (1,5,10,20,50) as a factor
  
  #depends on
  library(e1071)
  
  mydata<-data.frame(meanA=double(), 
                     meanB=double(), 
                     varA=double(), 
                     varB=double(), 
                     skewA=double(), 
                     skewB=double(), 
                     kurtA=double(), 
                     kurtB=double(),
                     stringsAsFactors = TRUE)
  
  pv1<-numeric(runs)
  pv5<-numeric(runs)
  pv10<-numeric(runs)
  pv20<-numeric(runs)
  pv50<-numeric(runs)
  
  wpv1<-numeric(runs)
  wpv5<-numeric(runs)
  wpv10<-numeric(runs)
  wpv20<-numeric(runs)
  wpv50<-numeric(runs)
  
  j<-1
  while (j <= runs){
    if (prange!=0){ #Has user specified randomized parameters?
      prange1 <- -1*prange
      a11<-a1*(1++runif(1, min=prange1, max=prange))
      a21<-a2*(1++runif(1, min=prange1, max=prange))
      b11<-b1*(1++runif(1, min=prange1, max=prange))
      b21<-b2*(1++runif(1, min=prange1, max=prange))
      #print(a11)
      #print(a21)
      #print(b11)
      #print(b21)
    } 
    tempA<-rbeta(batches, a11, b11)*maxCFU
    tempB<-rbeta(batches, a21, b21)*maxCFU
    temp5A<-numeric(batches)
    temp5B<-numeric(batches)
    temp10A<-numeric(batches)
    temp10B<-numeric(batches)
    temp20A<-numeric(batches)
    temp20B<-numeric(batches)
    temp50A<-numeric(batches)
    temp50B<-numeric(batches)
    
    for(k in 1:batches){
      temp1<-rbeta(5, a11, b11)*maxCFU
      temp5A[k]<-mean(temp1)
      temp2<-rbeta(5, a21, b21)*maxCFU
      temp5B[k]<-mean(temp2)
      temp1<-rbeta(10, a11, b11)*maxCFU
      temp10A[k]<-mean(temp1)
      temp2<-rbeta(10, a21, b21)*maxCFU
      temp10B[k]<-mean(temp2)
      temp1<-rbeta(20, a11, b11)*maxCFU
      temp20A[k]<-mean(temp1)
      temp2<-rbeta(20, a21, b21)*maxCFU
      temp20B[k]<-mean(temp2)
      temp1<-rbeta(50, a11, b11)*maxCFU
      temp50A[k]<-mean(temp1)
      temp2<-rbeta(50, a21, b21)*maxCFU
      temp50B[k]<-mean(temp2)
    }  
    #commit statistics to data frame
    datacheck<-c(tempA, tempB, temp5A, temp5B, temp10A, temp10B, temp20A, temp20B, temp50A, temp50B)
    if(sum(is.na(datacheck))==0){
      mydata<-rbind(mydata, 
                    data.frame(
                      meanA=c(mean(tempA),mean(temp5A),mean(temp10A),mean(temp20A),mean(temp50A)),
                      varA=c(var(tempA),var(temp5A),var(temp10A),var(temp20A),var(temp50A)),
                      skewA=c(skewness(tempA),skewness(temp5A),skewness(temp10A),skewness(temp20A),skewness(temp50A)),
                      kurtA=c(kurtosis(tempA),kurtosis(temp5A),kurtosis(temp5A),kurtosis(temp5A),kurtosis(temp5A)),
                      meanB=c(mean(tempB),mean(temp5B),mean(temp10B),mean(temp20B),mean(temp50B)),
                      varB=c(var(tempB),var(temp5B), var(temp10B), var(temp20B), var(temp50B)),
                      skewB=c(skewness(tempB),skewness(temp5B), skewness(temp10B),skewness(temp20B),skewness(temp50B)),
                      kurtB=c(kurtosis(tempB), kurtosis(temp5B), kurtosis(temp10B),kurtosis(temp20B), kurtosis(temp50B))
                    ))
      #t tests on data
      test1<-t.test(tempA, tempB)
      test5<-t.test(temp5A, temp5B)
      test10<-t.test(temp10A, temp10B)
      test20<-t.test(temp20A, temp20B)
      test50<-t.test(temp50A, temp50B)
      
      #Mann-Whitney U tests
      wtest1<-wilcox.test(tempA,tempB)
      wtest5<-wilcox.test(temp5A, temp5B)
      wtest10<-wilcox.test(temp10A, temp10B)
      wtest20<-wilcox.test(temp20A, temp20B)
      wtest50<-wilcox.test(temp50A, temp50B)
      
      pv1[j]<-test1$p.value
      pv5[j]<-test5$p.value
      pv10[j]<-test10$p.value
      pv20[j]<-test20$p.value
      pv50[j]<-test50$p.value
      
      wpv1[j]<-wtest1$p.value
      wpv5[j]<-wtest5$p.value
      wpv10[j]<-wtest10$p.value
      wpv20[j]<-wtest20$p.value
      wpv50[j]<-wtest50$p.value
      
      j<-j+1
    }
  }
  t.pvals<-c(sum(pv1<0.05)/runs, sum(pv5<0.05)/runs, sum(pv10<0.05)/runs, sum(pv20<0.05)/runs, sum(pv50<0.05)/runs)
  w.pvals<-c(sum(wpv1<0.05)/runs, sum(wpv5<0.05)/runs, sum(wpv10<0.05)/runs, sum(wpv20<0.05)/runs, sum(wpv50<0.05)/runs)
  print(t.pvals)
  print(w.pvals)
  batchlist=c(1,5,10,20,50)
  batch<-rep(batchlist,runs)
  mydata$batch<-as.factor(batch)
  return(mydata)
}
######################################################################################
######################################################################################


######################################################################################
######################################################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                       Section 1: 
#            Multi-Modality in CFU/Worm Data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# S. enterica and S. aureus bleach experiment data and simulations thereof
# The objects called here are generated in the main text script
# and the code calls the function "wormbootGMM()"
#library(tidyverse)
#library(cowplot)

library(mclust, quietly=TRUE)
library(sBIC)

# Using data for SE and SA with all experimental runs, zeros removed, run IDs as integers 1-3
# Combined data as "Pooled"
# Data are loaded from file on line 92 of main script WormAveraging.R
# The filtered object is created on line 124

mylen<-dim(SaSeCount2)[1]
SaSeCount2$Rep<-as.factor(SaSeCount2$Rep)
SaSeCountPool<-rbind(SaSeCount2, data.frame(Condition=SaSeCount2$Condition,
                                            Rep=rep("Pooled", mylen), # add pooled condition
                                            Count=SaSeCount2$Count,
                                            D=SaSeCount2$D,
                                            CFU=SaSeCount2$CFU,
                                            logCFU=SaSeCount2$logCFU
))
pSeCountAll<-SaSeCountPool %>%
  filter(Condition=="SE") %>%
  ggplot(aes(x=Rep, y=logCFU, color=Rep)) +
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + 
  theme_classic() + 
  scale_color_viridis_d(end=0.9)+
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=14),
        ) + 
  # facet_wrap(vars(Condition), scales="free_x") +
  labs(title=expression(paste(italic("S. enterica"), " LT2")), y=expression(log[10](CFU/Worm)), color="Replicate")
#ggsave("pSECountAll.png", width=4, height=3, units="in", dpi=300)
pSeCountAll

#ok let's try some gaussian mixture models (code from vignette)
X<-SaSeCountPool$logCFU[SaSeCountPool$Condition=="SE" & SaSeCountPool$Rep=="Pooled"]
mean(X)
var(X)
mean(SaSeCountPool$CFU[SaSeCountPool$Condition=="SE"& SaSeCountPool$Rep=="Pooled"])
sd(SaSeCountPool$CFU[SaSeCountPool$Condition=="SE"& SaSeCountPool$Rep=="Pooled"])

fit<- mclustBIC(X)
fit
#Top 3 models based on the BIC criterion: 
#  E,2       V,2       E,1 
#-174.8787 -177.0372 -180.1553 

plot(fit)
# build a tibble for ggplotting
fitBIC.df<-tibble(Components=rep(1:9),
                  E=fit[,1],
                  V=fit[,2],
               )

# plot
pSePooledBICMclust<-fitBIC.df %>%
  pivot_longer(
    cols=E:V,
    names_to="Var",
    values_to="BIC"
  ) %>%
  ggplot(aes(x=Components, y=BIC, linetype=factor(Var), pch=factor(Var)))+
  geom_line()+
  scale_x_continuous(breaks=seq(1, 9, by=1))+
  theme_classic() + 
  theme(text=element_text(size=14), 
        #axis.title.x = element_blank(), 
        #axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=14),
        legend.title = element_blank()) + 
  labs(title=expression(paste(italic("S. enterica"), " LT2, Pooled")), x="# Components", y="BIC")
pSePooledBICMclust


# Here we are still fitting to the combined Salmonella colonization data.
# Fitting a model with two components, allowing variances to be unequal ("V") or equal ("E")
fit2 <- Mclust(X, G=2, model="E")
summary(fit2)
plot(fit2, what="density", main="G2", xlab="logCFU")
rug(X)
fit2$parameters
#(mean1, var1)=(3.225159, 0.1877358) and (mean2, var2)=(4.662897, 0.1877358), containing 39% and 61% of the mass respectively

# let's put these plots together with the Se data plot
names(fit2)

# easier to ggplot if we call to density? This will be the same model as above
fit2d<-densityMclust(X, G=2)
summary(fit2d)
names(fit2d)
glimpse(fit2d)
rug(X)

# ok let's try ggplotting
as.numeric(fit2d$data)
fit2d_toplot<-tibble(logCFU=as.numeric(fit2d$data), density=as.numeric(fit2d$density))
glimpse(fit2d_toplot)
attach(fit2d_toplot)
newdata<-fit2d_toplot[order(logCFU),]
glimpse(newdata)
detach(fit2d_toplot)

pSeAllDensityMclust<-newdata %>%
  ggplot(aes(x=logCFU, y=density))+
  geom_line()+
  geom_rug(sides="b")+
  theme_classic() + 
  theme(text=element_text(size=14), 
        #axis.title.x = element_blank(), 
        #axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=14),
        legend.title = element_blank()) + 
  labs(title=expression(paste(italic("S. enterica"), " LT2, Pooled")), x=expression(log[10](CFU/Worm)), y="Density")
pSeAllDensityMclust

plot_grid(pSeCountAll, pSeAllDensityMclust, pSePooledBICMclust, ncol=3, labels="AUTO", align="h")
ggsave("FigS1_SentericaMclust.png", width=10, height=3, dpi=400, units="in")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# We can also look at each replicate individually
X1<-SaSeCountPool$logCFU[SaSeCountPool$Condition=="SE" & SaSeCountPool$Rep=="1"]
fit1<-mclustBIC(X1)
fit1 # Best support for G=2, model E
fit1 <- Mclust(X1, G=2, model="E")
summary(fit1)
plot(fit1, what="density", main="G2", xlab="logCFU")
rug(X1)
fit1$parameters
#(mean1, var1)=(3.3866, 0.05795) and (mean2, var2)=(4.74495, 0.05795), containing 34.6% and 65.4% of the mass respectively

X2<-SaSeCountPool$logCFU[SaSeCountPool$Condition=="SE" & SaSeCountPool$Rep=="2"]
fit2<-mclustBIC(X2)
fit2 #best support for E1, then E2
fit2 <- Mclust(X2, G=2, model="V")
fit2$parameters
#(mean1, var1)=(2.774017, 0.09134093) and (mean2, var2)=(4.319905, 0.32749870), containing 42.6% and 57.3% of the mass respectively

X3<-SaSeCountPool$logCFU[SaSeCountPool$Condition=="SE" & SaSeCountPool$Rep=="3"]
fit3<-mclustBIC(X3)
fit3 #Best support for E2
fit3 <- Mclust(X3, G=2, model="E")
summary(fit3)
plot(fit3, what="density", main="G2", xlab="logCFU")
rug(X3)
fit3$parameters
#(mean1, var1)=(3.7839, 0.06167) and (mean2, var2)=(4.9194, 0.06167), containing 42.1% and 57.9% of the mass respectively

#library(sBIC)
gMix = GaussianMixtures(maxNumComponents=10, phi=1, restarts=100)
set.seed(1234)
m = sBIC(X, gMix)
print(m)
matplot(
  cbind(m$BIC - m$BIC[1], m$sBIC - m$sBIC[1]),
  pch = c(1, 3),
  col = "black",
  xlab = "Number of components",
  ylab = expression(BIC - BIC(M[1])),
  las=1, xaxt="n"
)
axis(1, at = 1:10)
legend("bottomleft",
       c(expression(BIC), expression(bar(sBIC)[1])),
       pch = c(1, 3),
       y.intersp = 1.2)
#ggsave("pSE_GMM_BICvsBIC.png", width=4, height=3, units="in", dpi=300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                    FIGURE S2
# Using the all-SE 2-mode GMM, we can generate simulated data
# Note that we are allowing the weight in modes to be stochastic (binomial) 

# fit 1
#(mean1, var1)=(3.3866, 0.05795) and (mean2, var2)=(4.74495, 0.05795), containing 34.6% and 65.4% of the mass respectively
# fit 3
#(mean1, var1)=(3.7839, 0.06167) and (mean2, var2)=(4.9194, 0.06167), containing 42.1% and 57.9% of the mass respectively
# fit pooled
#(mean1, var1)=(3.225159, 0.1877358) and (mean2, var2)=(4.662897, 0.1877358), containing 39% and 61% of the mass respectively

SeBootGMM<-wormbootGMM(100, 50, 3.225159, 0.1877358, 4.662897, 0.1877358, 0.61, 0)
SeBootGMM$batch<-as.factor(SeBootGMM$batch)
SeBootGMM$set<-as.factor(SeBootGMM$set)

SeBootGMM.1.2<-wormbootGMM(100, 50, 3.225159, 0.1877358, 4.662897, 0.1877358, 0.61, 0.1)
SeBootGMM.1.2$batch<-as.factor(SeBootGMM.1.2$batch)
SeBootGMM.1.2$set<-as.factor(SeBootGMM.1.2$set)

SeBootGMM.1.2.05<-wormbootGMM(100, 50, 3.225159, 0.1877358, 4.662897, 0.1877358, 0.61, 0.05)
SeBootGMM.1.2.05$batch<-as.factor(SeBootGMM.1.2.05$batch)
SeBootGMM.1.2.05$set<-as.factor(SeBootGMM.1.2.05$set)


# Plot out each set of simulations and store as plot objects, to combine in figure
pSeBootGMM<-SeBootGMM %>%
  ggplot(aes(x=set, y=logCount, color=set))+
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + 
  ylim(-0.1,6) +
  theme_classic() + 
  scale_color_viridis_d(option="magma", begin=0.3, end=0.85)+
  theme(text=element_text(size=16), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=16),
        legend.position="none"
  ) + 
  labs(title="Pooled data parameterization", y=expression(log[10](CFU/worm)))+
  facet_wrap(~batch, nrow=1)+
  stat_compare_means(label.y=0.2)+
  stat_compare_means(method="t.test", label.y = 0.9)
pSeBootGMM

pSeBootGMM.1.2<-SeBootGMM.1.2 %>%
  ggplot(aes(x=set, y=logCount, color=set))+
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + 
  ylim(-0.1,6) +
  theme_classic() + 
  scale_color_viridis_d(option="magma", begin=0.3, end=0.85)+
  theme(text=element_text(size=16), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=16),
        legend.position="none"
  ) + 
  labs(title="Pooled data parameterization, 10% difference in mode proportions", y=expression(log[10](CFU/worm)))+
  facet_wrap(~batch, nrow=1)+
  stat_compare_means(label.y=0.2)+
  stat_compare_means(method="t.test", label.y = 0.9)
pSeBootGMM.1.2

pSeBootGMM.1.2.05<-SeBootGMM.1.2.05 %>%
  ggplot(aes(x=set, y=logCount, color=set))+
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + 
  ylim(-0.1,6) +
  theme_classic() + 
  scale_color_viridis_d(option="magma", begin=0.3, end=0.85)+
  theme(text=element_text(size=16), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=16),
        legend.position="none"
  ) + 
  labs(title="Pooled data parameterization, 5% difference in mode proportions", y=expression(log[10](CFU/worm)))+
  facet_wrap(~batch, nrow=1)+
  stat_compare_means(label.y=0.2)+
  stat_compare_means(method="t.test", label.y = 0.9)
pSeBootGMM.1.2.05

plot_grid(pSeBootGMM, pSeBootGMM.1.2, pSeBootGMM.1.2.05, ncol=1, labels="AUTO")
ggsave("FigS2_pSimSeGMM_run1v3.png", width=12, height=10, units="in", dpi=400)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                     Section 2:
# Simulating "CFU" data: effects of variation and skew on false-negative rates 
#
# Calls function simBetaCFU()
# to generate simulated data based on the beta distribution
# with a max of 5 logs CFU/individual (default)
# Start with right hand skew:
a<-0.1
b<-0.5
nreps<-50

simBetaRHS<-simBetaCFU(a,b,nreps)
glimpse(simBetaRHS)

pSimBetaRHS<-simBetaRHS %>%
	ggplot(aes(x=Factor, y=LogData, color=Factor))+
	geom_jitter(shape=16, position=position_jitter(0.05)) +
	geom_violin(fill=NA) + 
	theme_classic() + 
	theme(text=element_text(size=16), 
	axis.title.x = element_blank(), 
	plot.title=element_text(hjust=0.5, size=14),
	legend.position = "none") + 
	labs(title="Beta-distributed data", y=expression(log[10](CFU/worm)))
pSimBetaRHS + scale_x_discrete(labels=c("1"="B(0.1, 0.5)",
                                        "2.5"="B(0.25, 1.25)",
                                        "5"="B(0.5, 2.5)",
                                        "10"="B(1, 5)",
                                        "20"="B(2, 10)",
                                        "50"="B(5, 25)"))
ggsave("FigS3_pSimBetaRHS_1_5.png", width=8, height=4, units="in", dpi=400)
#ggsave("FigS3_pSimBetaRHS_1_5.pdf", width=8, height=4, units="in")

# note that the mean of these is all 100000*0.16667 = 16667 expected
# and all the actual arithmetic means are close to this, which is nice
simBetaRHS %>%
  group_by(Factor) %>%
  summarize(means=mean(Data),
            medians=median(Data)) %>%
  ungroup()


# comparisons
SimBetaRHS.ttests<-compare_means(Data~Factor, data=simBetaRHS, method="t.test")
SimBetaRHS.wtests<-compare_means(Data~Factor, data=simBetaRHS)
write.table(SimBetaRHS.ttests, "SimBetaRHS_ttests.txt", sep="\t")
write.table(SimBetaRHS.wtests, "SimBetaRHS_wtests.txt", sep="\t")

# we will focus on X5 Beta(0.5, 2.5) and X10 Beta(1,5)
# let's create batch data using wormSimBoot-based functions
simBatchBetaR5<-wormSimBatchBetaCompare(0.5, 2.5, 0.5, 2.5, 1000, 25, 100000)
simBatchBetaR10<-wormSimBatchBetaCompare(1, 5, 1, 5, 1000, 25, 100000)

simBatchBetaR5v10<-wormSimBatchBetaCompare(0.5, 2.5, 1, 5, 1000, 25, 100000)
simBatchBetaR5v20<-wormSimBatchBetaCompare(0.5, 2.5, 2, 10, 1000, 25, 100000)


#######
# now left hand skew
a<-0.5
b<-0.1
simBetaLHS<-simBetaCFU(a,b,nreps)
pSimBetaLHS<- simBetaLHS %>%
	ggplot(aes(x=Factor, y=LogData, color=Factor))+
	geom_jitter(shape=16, position=position_jitter(0.05)) +
	geom_violin(fill=NA) + 
	theme_classic() + 
	theme(text=element_text(size=16), 
	axis.title.x = element_blank(), 
	plot.title=element_text(hjust=0.5, size=16)) + 
	labs(title="Beta-distributed data", y="log10(CFU/worm)")
pSimBetaLHS
#ggsave("pSimBetaLHS_5_1.png", width=6, height=4, units="in", dpi=300)

# mean should be 100000*(5/6) = 83,333
# and it's close across the board
simBetaLHS %>%
  group_by(Factor) %>%
  summarize(means=mean(Data),
            medians=median(Data))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#        FIGURE S4  MOMENTS AND DISTANCES OF SIMULATED DATA; 
#                 RUN-TO-RUN VARIATION     
#
# ok back to beta starting at B(1,5) why not
# without much loss of generality let's draw parameter shifts from U(-0.1, 0.1) each run

wormSimBetaRand.1.5<-wormSimBatchBetaCompare(1,5,1,5,1000,24,100000, 0.1)
wormSimBetaRand.p5.2p5<-wormSimBatchBetaCompare(0.5,2.5,0.5,2.5,1000,24,100000, 0.1)
wormSimBetaRand.5.1<-wormSimBatchBetaCompare(5,1,5,1,1000,24,100000, 0.1)
wormSimBetaRand.2p5.p5<-wormSimBatchBetaCompare(2.5,0.5,2.5,0.5,1000,24,100000, 0.1)
wormSimBetaRand.5.5<-wormSimBatchBetaCompare(5,5,5,5,1000,24,100000, 0.1)
wormSimBetaRand.2p5.2p5<-wormSimBatchBetaCompare(2.5,2.5,2.5,2.5,1000,24,100000, 0.1)
wormSimBetaRand.1.1<-wormSimBatchBetaCompare(1,1,1,1,1000,24,100000, 0.1)

# calculate CV; include the distances
wormSimBetaRand.1.5$meandist<-abs(wormSimBetaRand.1.5$meanA-wormSimBetaRand.1.5$meanB)/(wormSimBetaRand.1.5$meanA+wormSimBetaRand.1.5$meanB)
wormSimBetaRand.1.5$cvA<-sqrt(wormSimBetaRand.1.5$varA)/wormSimBetaRand.1.5$meanA
wormSimBetaRand.1.5$cvB<-sqrt(wormSimBetaRand.1.5$varB)/wormSimBetaRand.1.5$meanB
wormSimBetaRand.1.5$cvdist<-abs(wormSimBetaRand.1.5$cvA-wormSimBetaRand.1.5$cvB)
wormSimBetaRand.1.5$skewdist<-abs(wormSimBetaRand.1.5$skewA-wormSimBetaRand.1.5$skewB)
wormSimBetaRand.1.5$kurtdist<-abs(wormSimBetaRand.1.5$kurtA-wormSimBetaRand.1.5$kurtB)

# Flip the data set around so it can get factor-gridded
wormSimBetaRand.1.5_A<-as_tibble(wormSimBetaRand.1.5) %>%
  select(batch, meanA, cvA, skewA, kurtA) %>%
  rename(mean=meanA,
         cv=cvA,
         skew=skewA,
         kurt=kurtA)
wormSimBetaRand.1.5_A$set<-"A"

wormSimBetaRand.1.5_B<-as_tibble(wormSimBetaRand.1.5) %>%
  select(batch, meanB, cvB, skewB, kurtB) %>%
  rename(mean=meanB,
         cv=cvB,
         skew=skewB,
         kurt=kurtB)
wormSimBetaRand.1.5_B$set<-"B"

wormSimBetaRand.1.5_dist<-as_tibble(wormSimBetaRand.1.5) %>%
  select(batch, meandist, cvdist, skewdist, kurtdist) %>%
  rename(mean=meandist,
         cv=cvdist,
         skew=skewdist,
         kurt=kurtdist)
wormSimBetaRand.1.5_dist$set<-"Distance"

# and combine
wormSimBetaRand.1.5_long<-rbind(wormSimBetaRand.1.5_A, wormSimBetaRand.1.5_B)
wormSimBetaRand.1.5_long$logMean<-log10(wormSimBetaRand.1.5_long$mean)

#plot out
pSimBatchBetaRand.1.5.AB<-wormSimBetaRand.1.5_long %>%
  relocate(set) %>%
  relocate(logMean, .after=mean) %>%
  select(-mean) %>%
  rename(average=logMean) %>%
  pivot_longer(average:kurt, names_to="moment", values_to="value") %>%
  ggplot(aes(x=batch, y=value, color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  scale_color_viridis_d(begin=0.3, end=0.8) +  
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(),  legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  #  labs(title=paste("\u03b2","(1,5) A", sep=""), y="Mean")
  facet_grid(rows=vars(moment), cols=vars(set), scales="free_y")

pSimBatchBetaRand.1.5.dist<-wormSimBetaRand.1.5_dist %>%
  relocate(set) %>%
  rename(average=mean) %>%
  pivot_longer(average:kurt, names_to="moment", values_to="value") %>%
  ggplot(aes(x=batch, y=value, color=batch))+
  geom_violin(fill=NA) + 
  theme_classic() + 
  scale_color_viridis_d(begin=0.3, end=0.8) +
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(),  legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  #  labs(title=paste("\u03b2","(1,5) A", sep=""), y="Mean")
  facet_grid(rows=vars(moment), cols=vars(set), scales="free_y") 

plot_grid(pSimBatchBetaRand.1.5.AB, pSimBatchBetaRand.1.5.dist, ncol=2, rel_widths = c(2,1), labels="AUTO")
ggsave("FigS4_pSimBatchBetaRand.1.5.moments.png", width=12, height=10, units="in", dpi=400)      

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################     FIGURE S5     #############################
#### add data from symmetric distribution for comparison
#note these runs use the flat U(-0.1,0.1) parameter error
wormSimBetaRand.5.5<-wormSimBatchBetaCompare(5,5,5,5,1000,24,100000, 0.1)
wormSimBetaRand.5.5$cvA<-sqrt(wormSimBetaRand.5.5$varA)/wormSimBetaRand.1.5$meanA
wormSimBetaRand.5.5$cvB<-sqrt(wormSimBetaRand.5.5$varB)/wormSimBetaRand.1.5$meanB
wormSimBetaRand.5.5$cvdist<-abs(wormSimBetaRand.5.5$cvA-wormSimBetaRand.1.5$cvB)
wormSimBetaRand.5.5$meandist<-abs(wormSimBetaRand.5.5$meanA-wormSimBetaRand.5.5$meanB)/(wormSimBetaRand.5.5$meanA+wormSimBetaRand.5.5$meanB)

pSimBatchBetaRand.1.5.meanA<-wormSimBetaRand.1.5 %>%
  ggplot(aes(x=batch, y=log10(meanA), color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
 # ylim(3.8, 5)+
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  scale_color_viridis_d(begin=0.3, end=0.8) +
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(),  
        legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title=paste("\u03b2","(1,5) A", sep=""), y="Mean")
#pSimBatchBetaRand.1.5.meanA

pSimBatchBetaRand.1.5.meandist<-wormSimBetaRand.1.5 %>%
  ggplot(aes(x=batch, y=meandist, color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  scale_color_viridis_d(begin=0.3, end=0.8) +
  ylim(-0.01,0.45)+
  theme(text=element_text(size=14), 
        #axis.title.x = element_blank(), 
        legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title="Distance", y="Dist(mean)")

pSimBatchBetaRand.5.5.meanA<-wormSimBetaRand.5.5 %>%
  ggplot(aes(x=batch, y=log10(meanA), color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
 # ylim(3.8, 5)+
  scale_color_viridis_d(begin=0.3, end=0.8) +
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(),  
        legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title=paste("\u03b2","(5,5) A", sep=""), y="Mean")

pSimBatchBetaRand.5.5.meandist<-wormSimBetaRand.5.5 %>%
  ggplot(aes(x=batch, y=meandist, color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  ylim(-0.01,0.45)+
  scale_color_viridis_d(begin=0.3, end=0.8) +
  theme(text=element_text(size=14), 
        #axis.title.x = element_blank(), 
        legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  #labs(title="Distance", y="Dist(mean)")
  labs(title=paste("\u03b2","(5,5) Distance"), y="Dist(Mean)")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Assembling Figure S5
plot_grid(pSimBatchBetaRand.1.5.meanA, pSimBatchBetaRand.5.5.meanA,
          pSimBatchBetaRand.1.5.meandist, pSimBatchBetaRand.5.5.meandist, 
          ncol=2, nrow=2, labels="AUTO", align = "h")
ggsave("FigS5_pSimBatchBetaRand.1.5.vs.5.5.means.png", width=8, height=6, units="in", dpi=400)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#               Figure S6
# Now more realistic - we allow triplicates, 
# and cap the number of both total worms and digests we will perform
#"small" is 100 worms and max 12 digests
#"medium" 200 worms max 24 digests
#"large" 500 worms max 48 digests

a<-1
b<-5
(a*b)/((a+b+1)*(a+b)^2)
2*(b-a)*sqrt(a+b+1)/(sqrt(a*b)*(a+b+2))

#nWorms<-100
#maxSamples<-12
#maxCFU<-100000

wormSim3.1.5.s<-wormSimBatchBeta3(1,5,100000,100,12,1000)
#wormSim3.1.5.s$set<-as.factor(wormSim3.1.5.s$set)
#wormSim3.1.5.s$logCFU<-log10(wormSim3.1.5.s$CFU)
wormSim3.1.5.m<-wormSimBatchBeta3(1,5,100000,200,24,1000)
wormSim3.1.5.l<-wormSimBatchBeta3(1,5,100000,500,48,1000)

wormSim3.5.5.s<-wormSimBatchBeta3(5,5,100000,100,12,1000)
wormSim3.5.5.m<-wormSimBatchBeta3(5,5,100000,200,24,1000)
wormSim3.5.5.l<-wormSimBatchBeta3(5,5,100000,500,48,1000)

wormSim3.3.7.s<-wormSimBatchBeta3(3,7,100000,100,12,1000)
wormSim3.3.7.m<-wormSimBatchBeta3(3,7,100000,200,24,1000)
wormSim3.3.7.l<-wormSimBatchBeta3(3,7,100000,500,48,1000)

# alternate version of function that returns means within run across three days of data
# wormSimBatchBeta3<-function(a, b, maxCFU, nWorms, maxSamples, reps, runs, sameParams=FALSE, returnMEANS=TRUE)
# here we use enough worms to allow 24 batches even at 50 worms

wormSim1.5.5.m.means<-wormSimBatchBeta3(5,5,100000,1200,24,1,1000)
wormSim1.1.5.m.means<-wormSimBatchBeta3(1,5,100000,1200,24,1,1000)
wormSim1.3.7.m.means<-wormSimBatchBeta3(3,7,100000,1200,24,1,1000)

wormSim3.5.5.m.means<-wormSimBatchBeta3(5,5,100000,1200,24,3,1000)
wormSim3.1.5.m.means<-wormSimBatchBeta3(1,5,100000,1200,24,3,1000)
wormSim3.5.5.m.means$batch<-as.factor(wormSim3.5.5.m.means$batch)
wormSim3.1.5.m.means$batch<-as.factor(wormSim3.1.5.m.means$batch)

wormSim3.1.5.m<-wormSimBatchBeta3(1,5,100000,1200,24,3,1000, returnMEANS=FALSE)
wormSim3.1.5.m$batch<-as.factor(wormSim3.1.5.m$batch)
wormSim3.1.5.m$day<-as.factor(wormSim3.1.5.m$day)

pwormSim3.1.5.m.A<-subset(wormSim3.1.5.m, set=="A") %>%
  ggplot(aes(x=batch, y=logCFU, color=day)) + 
  geom_violin(fill=NA) + 
  geom_point(shape=16, position=position_jitterdodge(0.2)) +
  ylim(1,5) + theme_classic() + 
  theme(text=element_text(size=14), 
        #axis.title.x = element_blank(), 
        plot.title=element_text(hjust=0.5, size=14),
        legend.position="none") + 
  labs(title="beta(1,5) set A", y="log10(CFU/worm)")
#pwormSim3.1.5.m.A
pwormSim3.1.5.m.B<-subset(wormSim3.1.5.m, set=="B") %>%
  ggplot(aes(x=batch, y=logCFU, color=day)) + 
  geom_violin(fill=NA) + 
  geom_point(shape=16, position=position_jitterdodge(0.2)) +
  ylim(1,5) + theme_classic() + 
  theme(text=element_text(size=14), 
        #axis.title.x = element_blank(), 
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title="beta(1,5) set B", y="log10(CFU/worm)")
plot_grid(pwormSim3.1.5.m.A, pwormSim3.1.5.m.B, rel_widths = c(1,1.15))
ggsave("pwormSim3.1.5.examplerun.png", width=8, height=4, units="in")

## plot differences in convergence of means for one-day vs. 3-day data
wormSim3.5.5.m.means$meandist<-abs(wormSim3.5.5.m.means$meanA-wormSim3.5.5.m.means$meanB)/(wormSim3.5.5.m.means$meanA+wormSim3.5.5.m.means$meanB)
pSimBatchBetaRand.5.5.mean3days<-wormSim3.5.5.m.means %>%
  ggplot(aes(x=batch, y=log10(meanA), color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  ylim(4.1,5)+
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(), legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title=paste("\u03b2","(5,5)x3"), y="Mean")
pSimBatchBetaRand.5.5.mean3days

pSimBatchBetaRand.5.5.meandist3days<-wormSim3.5.5.m.means %>%
  ggplot(aes(x=batch, y=meandist, color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  ylim(-0.01,0.2)+
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(), legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title=paste("\u03b2","(5,5)x3 Distance"), y="Dist(Mean)")
pSimBatchBetaRand.5.5.meandist3days
plot_grid(pSimBatchBetaRand.5.5.meanA,pSimBatchBetaRand.5.5.mean3days,
          pSimBatchBetaRand.5.5.meandist,pSimBatchBetaRand.5.5.meandist3days, 
          align="h")
ggsave("FigS6_pSimBatchBetaRand.5.5.meancompare3days.png", width=6, height=6, units="in")

## The mean converges nicely, and the distance between means is low, in the 3-day data

wormSim3.5.5.m.means$cvA<-sqrt(wormSim3.5.5.m.means$varA)/wormSim3.5.5.m.means$meanA
wormSim3.5.5.m.means$cvB<-sqrt(wormSim3.5.5.m.means$varB)/wormSim3.5.5.m.means$meanB
pSimBatchBetaRand.5.5.cv3days<-wormSim3.5.5.m.means %>%
  ggplot(aes(x=batch, y=cvA, color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(), legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title="Beta(5,5)x3", y="CV")
pSimBatchBetaRand.5.5.cv3days
plot_grid(pSimBatchBetaRand.5.5.cvA,pSimBatchBetaRand.5.5.cv3days)

pSimBatchBetaRand.1.5.mean3days<-wormSim3.1.5.m.means %>%
  ggplot(aes(x=batch, y=log10(meanA), color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  ylim(3.7,4.4)+
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(), legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title="Beta(1,5)x3", y="Mean")
plot_grid(pSimBatchBetaRand.1.5.meanA, pSimBatchBetaRand.1.5.mean3days, pSimBatchBetaRand.5.5.meanA, pSimBatchBetaRand.5.5.mean3days, nrow=2, ncol=2, labels="AUTO")
ggsave("pSimBatchBetaRand.1.5.vs.5.5.means.png", width=8, height=8, units="in")


###########  What if we increase the number of replicates?

wormSim10.5.5.m.means<-wormSimBatchBeta3(5,5,100000,1200,24,10,1000)
wormSim10.1.5.m.means<-wormSimBatchBeta3(1,5,100000,1200,24,10,1000)

wormSim25.5.5.m.means<-wormSimBatchBeta3(5,5,100000,1200,24,25,1000)
wormSim25.1.5.m.means<-wormSimBatchBeta3(1,5,100000,1200,24,25,1000)


################# NOT CURRENTLY USING ##############################
pwormSim3.1.5.s.batch1<-subset(wormSim3.1.5.s, batch==1) %>%
  ggplot(aes(x=set, y=logCFU, color=set)) + theme_classic() +
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  #    geom_violin(fill=NA) +  
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(), 
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title="1", y="log10(CFU/worm)")
pwormSim3.1.5.s.batch1

pwormSim3.1.5.l.batch1<-subset(wormSim3.1.5.l, batch==1) %>%
  ggplot(aes(x=set, y=logCFU, color=set)) + theme_classic() +
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  #    geom_violin(fill=NA) +  
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(), 
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title="1", y="log10(CFU/worm)")
pwormSim3.1.5.l.batch1

