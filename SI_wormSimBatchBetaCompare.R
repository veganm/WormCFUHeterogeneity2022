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

