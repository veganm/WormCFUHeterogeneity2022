wormSimBatchBeta3<-function(a, b, a_alt=1, b_alt=1, nWorms, maxSamples, runs, reps=3, batch_sizes=c(1,5,10,20,50),
                            maxCFU=1e5, prange=0.1, sameParams=TRUE){
  #Function that generates beta-distributed "worm CFU counts" with replicate days (default reps=3)
  #a and b are the parameters of the beta distribution
  #maxCFU is the scaling factor used to turn beta values to CFUs
  #nWorms is the total number of worms that we have available for each condition per replicate (>=12)
  #maxSamples is the maximum number of "digests" that we will perform, until nWorms limit is reached
  #reps is the number of independent replicates to be performed in each "experiment"
  #runs is the number of times to run the simulation
  #prange is the width of the uniform noise on beta parameters
  #   If sameparams=FALSE (default), two data sets A and B will be instantiated
  #   with data set B using base parameters (a_alt, b_alt)
  #   Otherwise, A and B have the same base parameters (a, b) within a replicate
  # Returns a list containing summaries of all runs and the most recent run of simulated data
  
  # needs
  pacman::p_load(e1071, tidyverse)
  
  # someplace to put data
  temp<-vector("list", length=length(batch_sizes)*runs)
  single_run<-vector("list", length=length(batch_sizes))
  
  idx<-1
  #runs is the number of times to run the simulation
  for (i in seq_len(runs)) {
    #reps is the number of independent replicates to be performed in each "experiment"
    for(k in seq_along(batch_sizes)){
      tempA<-numeric()  # someplace to keep all the data within an "experiment"
      tempB<-numeric()
      repA<-numeric()
      repB<-numeric()
      nbatches<-min(maxSamples, floor(nWorms/batch_sizes[k])) # how many data points can we make with this many worms?

      for (j in seq_len(reps)){
        temp1<-numeric(length=nbatches)  # someplace to put the data sets
        temp2<-numeric(length=nbatches)  
        a1<-a*(1+runif(1, min=-prange, max=prange))  # note there are always new parameters across replicates
        b1<-b*(1+runif(1, min=-prange, max=prange))
        if(sameParams){ # instantiate new parameters for set B?
          a_alt<-a
          b_alt<-b
        }
        a2<-a_alt*(1+runif(1, min=-prange, max=prange))
        b2<-b_alt*(1+runif(1, min=-prange, max=prange))
      
        for (m in seq_len(nbatches)){ # make individual data points
          temp1[m]<-mean(rbeta(batch_sizes[k], a1, b1)*maxCFU, na.rm=TRUE)
          temp2[m]<-mean(rbeta(batch_sizes[k], a2, b2)*maxCFU, na.rm=TRUE)
        } # end point by point generation
        tempA<-c(tempA, temp1) # concatenate data
        tempB<-c(tempB, temp2)
        repA<-c(repA, rep(j, nbatches)) # and batch IDs
        repB<-c(repB, rep(j, nbatches))
      } # end loop over replicates
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # construct raw data for this run
      mydata<-tibble(
          dist=c(rep(paste("Beta(", a, ",", b, ")", sep=""), nbatches*reps),
                 rep(paste("Beta(", a_alt, ",", b_alt, ")", sep=""),nbatches*reps)),
          batch=batch_sizes[k],
          set=c(rep("A", nbatches*reps), rep("B", nbatches*reps)),
          replicate=c(repA, repB),
          FinalCount=c(tempA, tempB)
        )
      mydata$logCFU<-log10(mydata$FinalCount)    # create & correct log final counts
      mydata$logCFU[!is.finite(mydata$logCFU)]<-1
      mydata$logCFU[mydata$logCFU<0]<-1
      #print(min(mydata$logCFU)) # debug
      
      mydata$replicate<-as.character(mydata$replicate) # make sure replicate is an index, not a number
      
      if(i==runs){#if this is the last run
        #print(i) # debug
        single_run[[k]]<-mydata
      }
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      # tests
      ttest1<-t.test(tempA, tempB)
      wtest1<-wilcox.test(tempA, tempB, exact=FALSE)
      
      # glm significance for "set"
      myglm<-mydata %>%
        glm(logCFU~replicate+set, family=Gamma, data=.)
      myglm.nested<-mydata %>%
        glm(logCFU~set, family=Gamma, data=.)
      myglm.unique<-mydata %>%
        mutate(unique_replicate=paste(set, replicate, sep="")) %>%
        glm(logCFU~unique_replicate, family=Gamma, data=.)
     myglm.lrt<-lrtest(myglm, myglm.nested)
     #print(AIC(myglm)) # debug
     #print(AIC(myglm.unique)) # debug
      
      # ANOVA
      my.aov<-aov(logCFU~replicate+set, data=mydata)
      
      # store summary stats
      temp[[idx]]<-tibble(
        dist=paste("Beta(", a, ",", b, ")", sep=""),
        distB=paste("Beta(", a_alt, ",", b_alt, ")", sep=""),
        Run=i,
        nreps=reps,
        batch=batch_sizes[k],
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
        p.t=ttest1$p.value,
        p.w=wtest1$p.value,
        p.glm=coef(summary(myglm))["setB",4], # obtain p values
        p.glm.lrt=myglm.lrt$`Pr(>Chisq)`[2],
        p.glm.AIC=exp((AIC(myglm)-AIC(myglm.nested))/2),
        p.glm.AICu=exp((AIC(myglm.unique)-AIC(myglm))/2),
        p.aov=summary(my.aov)[[1]][2,5]
      )
      idx<-idx+1
  } # end loop over batch sizes (k)

  } # end loop over runs (i)
  
  dataSet<-dplyr::bind_rows(temp)  # unpack
  dataSet<-dataSet %>%
    mutate(
      meandist=abs(mean_A-mean_B)/(mean_A+mean_B),
      cvdist=abs((sqrt(var_A)/mean_A) - (sqrt(var_B)/mean_B)),
      skewdist=abs(skew_A-skew_B),
      kurtdist=abs(kurt_A-kurt_B)
    )
  
  singleRunData<-dplyr::bind_rows(single_run) # unpack most recent run
  my_list<-list("data_summary" = dataSet, # returns summary data
                "single_run" = singleRunData) # and the most recent single run

  return(my_list)
}

