wormSimBatchBetaFactorial<-function(a, b, runs = 4, nbatches=10, batch_sizes=c(1,5,10,20,50),
                                    maxCFU=1e5, prange=0.1, sameParams=FALSE){
  #Function that generates two sets (A and B) of beta-distributed "worm CFU counts" with factorial design
  # Simplified version of wormSimBatchBeta3() without embedded tests
  # and with a constant number of measurements per experiment at all batch sizes
  #a and b are the parameters of the beta distribution
  #maxCFU is the scaling factor used to turn beta values to CFUs
  #nbatches is the number of measurements to be performed in each "experiment" (default 10); 
  #  all batch sizes have the same number of replicate measurements 
  #runs is the number of independent "experiments" (default 4)
  # prange is the uniform range of parameter noise
  #When sameParams=FALSE, Beta parameters are re-drawn in data sets A and B
  #Returns a simulated biologically-averaged data set for analysis
  
  # needs
  pacman::p_load(tidyverse)

  # Someplace to put data
  temp<-vector("list", length=length(batch_sizes)*runs)
  
  # Iterate
  idx<-1
  for (i in seq_len(runs)){
    # new parameters each run
    a1<-a*(1+runif(1, min=-prange, max=prange))
    b1<-b*(1+runif(1, min=-prange, max=prange))
    if(sameParams){
      a2<-a1
      b2<-b1
    }else {
      a2<-a*(1+runif(1, min=-prange, max=prange))
      b2<-b*(1+runif(1, min=-prange, max=prange))
    }
    
    for(k in seq_along(batch_sizes)){
      tempA<-numeric()  # someplace to keep all the data within an "experiment"
      tempB<-numeric()
        for (m in seq_len(nbatches)){ # make individual data points
          tempA[m]<-mean(rbeta(batch_sizes[k], a1, b1)*maxCFU, na.rm=TRUE)
          tempB[m]<-mean(rbeta(batch_sizes[k], a2, b2)*maxCFU, na.rm=TRUE)
        } # end point by point generation
      # store
      temp[[idx]]<-tibble(
        dist=paste("Beta(", a, ",", b, ")", sep=""),
        batch=batch_sizes[k],
        run=i,
        set=c(rep("A", nbatches), rep("B", nbatches)),
        FinalCount=c(tempA, tempB)
        )
      print(c(i,k, idx)) #debug
      idx<-idx+1 # shift position
      } # end loop over batch sizes
    } # end loop over runs
  mydata<-dplyr::bind_rows(temp)  # unpack
  mydata$logCFU<-log10(mydata$FinalCount)    # create & correct log final counts
  mydata$logCFU[!is.finite(mydata$logCFU)]<-1
  mydata$logCFU[mydata$logCFU<0]<-1
  #print(min(mydata$logCFU)) # debug
  
  return(mydata)
}
