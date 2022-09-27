wormbootOnCounts<-function(reps, mydata, Dcorrect){
  # Expects a number of reps for the bootstrap (reps)
  # and a data frame of worm CFU data for individuals (mydata)
  # where the number of colonies counted is in column "Count"
  # the dilution at which these colonies were measured is in column "D"
  # and the CFU/worm is in column "CFU"
  #
  ## The dilution correction factor (numeric) is given as "Dcorrect"
  # 
  # Returns a data frame of simulated batch digests
  # with batch sizes 1, 5, 10, 20, 50 worms/batch
  # values reported as inferred CFU/worm and log10(CFU/worm)
  capp<-dim(mydata)[1]
  batch5<-rep(0,reps)
  batch10<-rep(0,reps)
  batch20<-rep(0,reps)
  batch50<-rep(0, reps)
  for(i in 1:reps){
    idx5<-sample(1:capp,5,replace=TRUE)
    idx10<-sample(1:capp,10,replace=TRUE)
    idx20<-sample(1:capp,20,replace=TRUE)
    idx50<-sample(1:capp,50,replace=TRUE)
    
    # Assume Poisson count error
    temp_count<-rpois(5, mydata$Count[idx5])
    temp<-Dcorrect*temp_count*10^mydata$D[idx5]
    batch5[i]<-mean(temp)
    
    temp_count<-rpois(10, mydata$Count[idx10])
    temp<-Dcorrect*temp_count*10^mydata$D[idx10]
    batch10[i]<-mean(temp)
    
    temp_count<-rpois(20, mydata$Count[idx20])
    temp<-Dcorrect*temp_count*10^mydata$D[idx20]
    batch20[i]<-mean(temp)

    temp_count<-rpois(50, mydata$Count[idx50])
    temp<-Dcorrect*temp_count*10^mydata$D[idx50]
    batch50[i]<-mean(temp)
    }
  batch5log<-log10(batch5+1)
  batch10log<-log10(batch10+1)
  batch20log<-log10(batch20+1) 
  batch50log<-log10(batch50+1) 
  Batch<-c(rep(1,times=capp), rep(5, times=reps), rep(10, times=reps), rep(20, times=reps), rep(50, times=reps))
  mylogCFU<-log10(mydata$CFU+1)
  logCFU<-c(mylogCFU, batch5log, batch10log, batch20log, batch50log)
  CFU<-c(mydata$CFU, batch5, batch10, batch20, batch50)
  dataSet<-data.frame(Batch, CFU, logCFU)
  return(dataSet)
}
