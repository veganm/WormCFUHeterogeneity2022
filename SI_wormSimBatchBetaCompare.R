#########################################################################################
#########################################################################################
# function wormSimBatchBetaCompare
wormSimBatchBetaCompare<-function(a1, b1, a2, b2, runs=10, batches=24, batch_sizes=c(1,5,10,20,50),
                                  maxCFU=1e5, prange=0){
  #Function that simulates worm CFU data based on draws from beta distribution
  # and compares simulated distributions using t-tests & Wilcoxon
  #
  # If prange is specified:
  # Parameters for data sets 1(A) and 2(B) are re-drawn separately for each run 
  #   with errors U(-prange, prange)
  #
  #a and b are underlying parameters of the beta distribution (double)
  #runs is the number of runs in each simulation (integer)
  #batches is the number of data points to generate for each run (integer)
  #maxCFU is the maximum value you want for the CFU/worm synthetic data (number)
  #return a data frame containing moments of data from each rep
  #column "batch" is the batch size (1,5,10,20,50) as a factor
  
  #depends on
  library(e1071)
  
  temp<-vector("list", length=length(batch_sizes)*runs*batches)

  idx<-1
  for (i in seq_len(runs)){
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
    for (j in seq_along(batch_sizes)){ # at each batch size
      tempA<-numeric(length=batches) # create new data sets A and B
      tempB<-numeric(length=batches)
      for (k in seq_len(batches)){ # for each data point
        tempA[k]<-mean(rbeta(batch_sizes[j], a11, b11)*maxCFU, na.rm=TRUE)
        tempB[k]<-mean(rbeta(batch_sizes[j], a21, b21)*maxCFU, na.rm=TRUE)
      } # close data point
      
      # stats comparing data sets A and B
      test1<-t.test(tempA, tempB)
      wtest1<-wilcox.test(tempA,tempB)
      
      # store
      temp[[idx]]<-tibble( # store summaries
        dist_A=paste("Beta(", a1, ",", b1, ")", sep=""),
        dist_B=paste("Beta(", a2, ",", b2, ")", sep=""),
        Run=i,
        batch=batch_sizes[j],
        n=batches,
        mean_A=mean(tempA, na.rm=TRUE),
        logmean_A=log10(mean(tempA, na.rm=TRUE)),
        var_A=var(tempA,na.rm=TRUE),
        logvar_A=log10(var(tempA,na.rm=TRUE)),
        cv_A=sd(tempA,na.rm=TRUE)/mean(tempA,na.rm=TRUE),
        skew_A=skewness(tempA,na.rm=TRUE),
        kurt_A=kurtosis(tempA,na.rm=TRUE),
        mean_B=mean(tempB,na.rm=TRUE),
        logmean_B=log10(mean(tempB, na.rm=TRUE)),
        var_B=var(tempB,na.rm=TRUE),
        logvar_B=log10(var(tempB,na.rm=TRUE)),
        cv_B=sd(tempB,na.rm=TRUE)/mean(tempB,na.rm=TRUE),
        skew_B=skewness(tempB,na.rm=TRUE),
        kurt_B=kurtosis(tempB,na.rm=TRUE), 
        p.t=test1$p.value,
        p.w=wtest1$p.value
      ) # close storage
      idx<-idx+1
    } # close batch size
  } # close run loop
  
  dataSet<-dplyr::bind_rows(temp)  # unpack
  dataSet<-dataSet %>%
    mutate(
      comparison=paste(dist_A, "v", dist_B, sep=" "),
      meandist=abs(mean_A-mean_B)/(mean_A+mean_B),
      cvdist=abs((sqrt(var_A)/mean_A) - (sqrt(var_B)/mean_B)),
      skewdist=abs(skew_A-skew_B),
      kurtdist=abs(kurt_A-kurt_B)
    )
  return(dataSet)
}
