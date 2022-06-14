wormSimBatchBeta3<-function(a, b, maxCFU, nWorms, maxSamples, reps, runs, sameParams=FALSE, returnMEANS=TRUE){
  #Function that generates beta-distributed "worm CFU counts" with triplicate days
  #a and b are the parameters of the beta distribution
  #maxCFU is the scaling factor used to turn beta values to CFUs
  #nWorms is the total number of worms that we have available for each condition per day (>=12)
  #maxSamples is the maximum number of "digests" that we will perform, until nWorms limit is reached
  #reps is the number of independent replicates to be performed in each "experiment"
  #runs is the number of times to run the simulation
  #can return meanCFU or one run of data
  meanCFU<-data.frame(meanA=double(), 
                      meanB=double(),
                      varA=double(),
                      varB=double(),
                     batch=double(), 
                     stringsAsFactors = TRUE)
  
  wpv1<-numeric(runs)
  wpv5<-numeric(runs)
  wpv10<-numeric(runs)
  wpv20<-numeric(runs)
  wpv50<-numeric(runs)

  #determine how many batches we can create at each size
  batches5<-min(maxSamples, nWorms/5)
  batches10<-min(maxSamples, nWorms/10)
  batches20<-min(maxSamples, nWorms/20)
  batches50<-min(maxSamples, nWorms/50)
  
  #prange<-min(a,b)/10
  prange<-0.1
  
  j<-1
  while (j <= runs){
    i<-1
    my1A<-numeric()
    my5A<-numeric()
    my10A<-numeric()
    my20A<-numeric()
    my50A<-numeric()
    my1B<-numeric()
    my5B<-numeric()
    my10B<-numeric()
    my20B<-numeric()
    my50B<-numeric()
    mydata<-data.frame(CFU=double(), 
                       batch=double(), 
                       day=double(), 
                       set=character(), 
                       stringsAsFactors = TRUE)
    while(i <= reps){
     temp5A<-numeric(batches5)
     temp5B<-numeric(batches5)
     temp10A<-numeric(batches10)
     temp10B<-numeric(batches10)
     temp20A<-numeric(batches20)
     temp20B<-numeric(batches20)
     temp50A<-numeric(batches50)
     temp50B<-numeric(batches50)

      a1<-a*(1+runif(1, min=-prange, max=prange))
      b1<-b*(1+runif(1, min=-prange, max=prange))
      if(sameParams){
       a2<-a1
       b2<-b1
      }else {
        a2<-a*(1+runif(1, min=-prange, max=prange))
        b2<-b*(1+runif(1, min=-prange, max=prange))
      }
       tempA<-rbeta(maxSamples, a1, b1)*maxCFU
       tempB<-rbeta(maxSamples, a2, b2)*maxCFU
    for(k in 1:batches5){
      temp1<-rbeta(5, a1, b1)*maxCFU
      temp5A[k]<-mean(temp1)
      temp2<-rbeta(5, a2, b2)*maxCFU
      temp5B[k]<-mean(temp2)}
    for (k in 1:batches10){
      temp1<-rbeta(10, a1, b1)*maxCFU
      temp10A[k]<-mean(temp1)
      temp2<-rbeta(10, a2, b2)*maxCFU
      temp10B[k]<-mean(temp2)}
    for (k in 1:batches20){
      temp1<-rbeta(20, a1, b1)*maxCFU
      temp20A[k]<-mean(temp1)
      temp2<-rbeta(20, a2, b2)*maxCFU
      temp20B[k]<-mean(temp2)}
    for (k in 1:batches50){
      temp1<-rbeta(50, a1, b1)*maxCFU
      temp50A[k]<-mean(temp1)
      temp2<-rbeta(50, a2, b2)*maxCFU
      temp50B[k]<-mean(temp2)}
    #commit statistics to data frame
    datacheck<-c(tempA, tempB, temp5A, temp5B, temp10A, temp10B, temp20A, temp20B, temp50A, temp50B)
    if(sum(is.na(datacheck))==0){
      mydata<-rbind(mydata, 
                    data.frame(
                      CFU=datacheck,
                      batch=c(rep(1,(maxSamples*2)), rep(5,(batches5*2)), rep(10,(batches10*2)), rep(20,(batches20*2)), rep(50,(batches50*2))),
                      day=rep(i, length(datacheck)),
                      set=c(rep("A",maxSamples), rep("B",maxSamples),
                            rep("A", batches5), rep("B", batches5),
                            rep("A", batches10), rep("B", batches10),
                            rep("A", batches20), rep("B", batches20),
                            rep("A", batches50), rep("B", batches50)
                            )
                    ))}
    my1A<-c(my1A,tempA)
    my5A<-c(my5A,temp5A)
    my10A<-c(my10A,temp10A)
    my20A<-c(my20A,temp20A)
    my50A<-c(my50A,temp50A)
    my1B<-c(my1B,tempB)
    my5B<-c(my5B,temp5B)
    my10B<-c(my10B,temp10B)
    my20B<-c(my20B,temp20B)
    my50B<-c(my50B,temp50B)
    
    i<-i+1
    }
#    print(my1)
    meanCFU<-rbind(meanCFU, data.frame(
      meanA=c(mean(my1A), mean(my5A), mean(my10A), mean(my20A), mean(my50A)),
      meanB=c(mean(my1B), mean(my5B), mean(my10B), mean(my20B), mean(my50B)),
      varA=c(var(my1A), var(my5A), var(my10A), var(my20A), var(my50A)),
      varB=c(var(my1B), var(my5B), var(my10B), var(my20B), var(my50B)),
      batch=c(1, 5, 10, 20, 50)
    ))
    
      #Mann-Whitney U tests
      wtest1<-wilcox.test(mydata$CFU[mydata$batch==1 & mydata$set=="A"],mydata$CFU[mydata$batch==1 & mydata$set=="B"])
      wtest5<-wilcox.test(mydata$CFU[mydata$batch==5 & mydata$set=="A"],mydata$CFU[mydata$batch==5 & mydata$set=="B"])
      wtest10<-wilcox.test(mydata$CFU[mydata$batch==10 & mydata$set=="A"],mydata$CFU[mydata$batch==10 & mydata$set=="B"])
      wtest20<-wilcox.test(mydata$CFU[mydata$batch==20 & mydata$set=="A"],mydata$CFU[mydata$batch==20 & mydata$set=="B"])
      wtest50<-wilcox.test(mydata$CFU[mydata$batch==50 & mydata$set=="A"],mydata$CFU[mydata$batch==50 & mydata$set=="B"])
      
      wpv1[j]<-wtest1$p.value
      wpv5[j]<-wtest5$p.value
      wpv10[j]<-wtest10$p.value
      wpv20[j]<-wtest20$p.value
      wpv50[j]<-wtest50$p.value
      
      j<-j+1
  }
  w.pvals<-c(sum(wpv1<0.05)/runs, sum(wpv5<0.05)/runs, sum(wpv10<0.05)/runs, sum(wpv20<0.05)/runs, sum(wpv50<0.05)/runs)
  print(w.pvals)
  mydata$logCFU<-log10(mydata$CFU)
  mydata$logCFU[!is.finite(mydata$logCFU)]<-0
  
  if(returnMEANS){return(meanCFU)}
  else {return(mydata)}
  
}