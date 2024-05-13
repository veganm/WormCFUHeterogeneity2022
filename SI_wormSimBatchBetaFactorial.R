wormSimBatchBetaFactorial<-function(a, b, maxCFU, runs = 4, reps=10, sameParams=FALSE){
  #Function that generates two sets (A and B) of beta-distributed "worm CFU counts" with factorial design
  # Simplified version of wormSimBatchBeta3() without embedded tests
  # and with a constant number of measurements per experiment at all batch sizes
  #a and b are the parameters of the beta distribution
  #maxCFU is the scaling factor used to turn beta values to CFUs
  #reps is the number of measurements to be performed in each "experiment" (default 10); 
  #  all batch sizes have the same number of replicate measurements 
  #runs is the number of independent "experiments" (default 4)
  #When sameParams=FALSE, Beta parameters are given different noise in data sets A and B
  #Returns a simulated biologically-averaged data set for analysis
  
  # make places to put the data
  mydata<-data.frame(CFU=double(), 
                     batch=double(), 
                     run=double(), 
                     set=character(), 
                     stringsAsFactors = TRUE)
  
  # as usual, the range of parameter noise
  prange<-0.1
  
  j<-1 # set index on runs
  while (j <= runs){
      temp5A<-numeric(reps)
      temp5B<-numeric(reps)
      temp10A<-numeric(reps)
      temp10B<-numeric(reps)
      temp20A<-numeric(reps)
      temp20B<-numeric(reps)
      temp50A<-numeric(reps)
      temp50B<-numeric(reps)
      
      a1<-a*(1+runif(1, min=-prange, max=prange))
      b1<-b*(1+runif(1, min=-prange, max=prange))
      if(sameParams){
        a2<-a1
        b2<-b1
      }else {
        a2<-a*(1+runif(1, min=-prange, max=prange))
        b2<-b*(1+runif(1, min=-prange, max=prange))
      }
      tempA<-rbeta(reps, a1, b1)*maxCFU
      tempB<-rbeta(reps, a2, b2)*maxCFU
      for(k in 1:reps){
        temp1<-rbeta(5, a1, b1)*maxCFU
        temp5A[k]<-mean(temp1, na.rm=TRUE)
        temp2<-rbeta(5, a2, b2)*maxCFU
        temp5B[k]<-mean(temp2, na.rm=TRUE)
        temp1<-rbeta(10, a1, b1)*maxCFU
        temp10A[k]<-mean(temp1, na.rm=TRUE)
        temp2<-rbeta(10, a2, b2)*maxCFU
        temp10B[k]<-mean(temp2, na.rm=TRUE)
        temp1<-rbeta(20, a1, b1)*maxCFU
        temp20A[k]<-mean(temp1, na.rm=TRUE)
        temp2<-rbeta(20, a2, b2)*maxCFU
        temp20B[k]<-mean(temp2, na.rm=TRUE)
        temp1<-rbeta(50, a1, b1)*maxCFU
        temp50A[k]<-mean(temp1, na.rm=TRUE)
        temp2<-rbeta(50, a2, b2)*maxCFU
        temp50B[k]<-mean(temp2, na.rm=TRUE)
        }
      #commit statistics to data frame
      datacheck<-c(tempA, tempB, temp5A, temp5B, temp10A, temp10B, temp20A, temp20B, temp50A, temp50B)
      if(sum(is.na(datacheck))==0){
        mydata<-rbind(mydata, 
                      data.frame(
                        CFU=datacheck,
                        batch=c(rep(1,(reps*2)), 
                                rep(5,(reps*2)), 
                                rep(10,(reps*2)), 
                                rep(20,(reps*2)), 
                                rep(50,(reps*2))),
                        run=rep(j, length(datacheck)),
                        set=c(rep("A", reps), rep("B",reps),
                              rep("A", reps), rep("B", reps),
                              rep("A", reps), rep("B", reps),
                              rep("A", reps), rep("B", reps),
                              rep("A", reps), rep("B", reps)
                        )))
      j<-j+1
      }
  }
  mydata$logCFU<-log10(mydata$CFU)
  mydata$logCFU[!is.finite(mydata$logCFU)]<-0
  
  return(mydata)
}
