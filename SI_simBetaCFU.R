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
