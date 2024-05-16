#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sim_sample_noise<-function(raw_data, foldD=10, correction_constant=1, n_countable=50){
  # Function to simulate sampling noise on count data  
  # Called in sim_means_named()
  # Inputs:
  # raw_data (num vector): Numeric vector of final counts
  # FoldD (num): Fold dilution in dilution series. Default is 10X dilutions.
  # correction_constant (num): a multiplier to correct the fraction of a sample used for one measurement
  # n_countable (num): Upper bound for "countable" number of events (threshold of "too many to count")
  
  # make places to put stuff
  n_points<-length(raw_data) # how many data points are there?

  N_d<-rep(NA, times=n_points) # store noised counts
  
  for(i in seq_len(n_points)){ # Iterate over data
    d<-0 # start at dilution 0
    n_new<-raw_data[i]
    # print(n_new) # for debugging
    n_counted<-round(rbinom(n=1, size=n_new, prob=1/correction_constant))
    while(n_counted>n_countable){
      d<-d+1
      n_new<-rbinom(n=1, size=n_new, prob=1/foldD) # number in diluted sample
      n_counted<-round(rbinom(n=1, size=n_new, prob=1/correction_constant))
    }
    N_d[i]<-correction_constant * (n_counted*foldD^d)
  }
  return(N_d)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sim_means_named<-function(input_data, n_reps, batch_sizes, Batch=1, foldD=10, alpha_start=1, beta_start=2, correction_constant=1, n_countable=50){
  # A function to generate simulated distributions of means based on raw data
  # using the normal, lognormal, and Beta distributions (latter is non-unique) 
  # Data generated are average values for simulated "samples"
  # n_reps: number of values to generate at each batch size
  # batch_sizes: vector of batch sizes for simulations (number of individuals per batch)
  # Batch (num): Number of individuals/units in the batch for each measurement in input_data 
  #               Default (Batch==1) for individual (unbatched) measurements.
  # alpha_start and beta_start: initialization for beta distribution fits
  # correction_constant: multiplier to correct the fraction of a sample used for one measurement
  #   (e.g. for CFU data, measurements from 10 uL spots need a correction of (1000/10)=100 for an original volume of 1 mL)
  # n_countable: upper bound for "countable" number of events (threshold of "too many to count")
  #   (for use in simulating sampling noise)
  # Returns a tibble where observed mean values are stored for each "sample"

  # input_data should be a tibble where each row is one data point, 
  #  and all data points are taken from the same sampling frame (one experiment, sampling run,etc)
  #  and represent the same batch size.
  #  The data must have columns:
  #   - Count (num): raw counts (e.g. number of colonies counted for CFU data)
  # The data may have columns:
  #   - D (num):  Serial dilution factor. A value of 0 (default) is undiluted; 1 indicates the first N-fold dilution; etc.
  #               If this column is not present, the code uses the raw counts * correction_constant for all calculations.
  #               If this column is present, the code will check for a column named "CFU" and calculate it if missing.
  #   - FoldD(num): Fold dilution in dilution series. default is 10X dilutions.
  #   - FinalCount (num): Corrected counts. ***Will be calculated from raw counts if absent.***
  
    
  # load in packages
  pacman::p_load(MASS, tidyverse) 

    # Create an empty temporary vector of defined size to hold the data
  temp<-vector("list", length=length(batch_sizes)*n_reps)
  
  # Do we need to correct counts for dilution?
  # Divide by batch size (may be 1) to get per-unit or per-individual numbers
  my_names<-names(input_data) 
  if(length(which(my_names=="FinalCount"))>0){
    if(length(which(my_names=="D"))!=0){ # if we are given a dilution factor
      input_data$FinalCount<-(correction_constant *(Count*foldD^D))/Batch
    } else {
      input_data$FinalCount<-(correction_constant * Count)/Batch
    }
  } 
  # now we should have final counts
  
  # For simulations, we'll use the corrected count data
  mydata<-input_data %>%
    pull(FinalCount)
  
  # Get summary stats for data
  shared_mean<-mean(mydata, na.rm=TRUE)
  shared_var<-var(mydata, na.rm=TRUE)
  shared_sd<-sqrt(shared_var)
  shared_mean_log<-mean(log(mydata+1), na.rm=TRUE)
  shared_var_log<-var(log(mydata+1), na.rm=TRUE)
  shared_sd_log<-sqrt(shared_var_log)
  
  ####  Begin parameter generation
  # If there are individual-sample data
  if(Batch==1){
    ### Parameterize the beta distribution from data
    # First, range-normalize to [0,1]
    mydata_rn<-(mydata - min(mydata))/(max(mydata)-min(mydata))
    # then fit shape parameters
    mydata_rn[mydata_rn==0]<-0.001 # Beta(0,1)
    mydata_rn[mydata_rn==1]<-0.999 # not inclusive
    beta_fit<-fitdistr(mydata_rn, dbeta, list(shape1=alpha_start, shape2=beta_start))
    a<-beta_fit$estimate[[1]]
    b<-beta_fit$estimate[[2]]
    ##############################################
  } else { ############## if we have batched data, have to do this differently
    # First, the variance of the normal needs to be corrected for batch size
    shared_var<-shared_var*Batch
    shared_sd<-sqrt(shared_var)
    
    # same for the lognormal (from Caudill 2010)
    shared_cv<-shared_sd/shared_mean
    shared_var_log<-log(Batch*shared_cv^2 + 1)
    shared_sd_log<-sqrt(shared_var_log)
    
    # The beta is a problem, because this is not closed-form. 
    # So we'll have to hack it.
    # First, we'll simulate some data from a lognormal with the right moments
    mydata<-round(rlnorm(25, meanlog=shared_mean_log, sdlog=shared_sd_log))
    
    # now range-normalize to [0,1]
    mydata_rn<-(mydata - min(mydata))/(max(mydata)-min(mydata))
    # then fit shape parameters
    mydata_rn[mydata_rn==0]<-0.001 # Beta(0,1)
    mydata_rn[mydata_rn==1]<-0.999 # not inclusive
    beta_fit<-fitdistr(mydata_rn, dbeta, list(shape1=alpha_start, shape2=beta_start))
    # and flip the parameters to change the direction of skew
    b<-beta_fit$estimate[[1]]
    a<-beta_fit$estimate[[2]]
    # I said this was going to be a hack.
  } # end of parameter generation
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # now we have parameters and can create simulated data.
    
    # Establish an index which will be used to store data
    idx<-1
    
    # Iterate over each of the values in the vector of batch sizes
    for (i in 1:length(batch_sizes)){
      #print(i) # for debugging
      # Within each sample size, execute the indicated number of replicates
      for (j in 1:n_reps){
        #print(j) # for debugging
        # Within each replicate, create a tibble with the mean from each sample
        # and staple it into the temporary vector
        
        # first generate the raw data without heterogeneity
        sim_normal<-round(rnorm(batch_sizes[i], mean=shared_mean, sd=shared_sd))
        sim_zeros<-which(sim_normal<0)
        sim_normal[sim_zeros]<-0
        sim_lnorm<-round(rlnorm(batch_sizes[i], meanlog=shared_mean_log, sdlog=shared_sd_log))
        sim_beta<-rbeta(batch_sizes[i], shape1 = a, shape2 = b)
        # de-normalize the beta
        sim_beta<-round((sim_beta * (max(mydata)-min(mydata))) + min(mydata))
        
        # then add sampling noise
        # (function defined above in this script)
        sim_normal_noise<-sim_sample_noise(round(sum(sim_normal)), 
                                           foldD=foldD, correction_constant = correction_constant, n_countable = n_countable)/batch_sizes[i]
        sim_lnormal_noise<-sim_sample_noise(round(sum(sim_lnorm)), 
                                            foldD=foldD, correction_constant = correction_constant, n_countable = n_countable)/batch_sizes[i]
        sim_beta_noise<-sim_sample_noise(round(sum(sim_beta)), 
                                         foldD=foldD, correction_constant = correction_constant, n_countable = n_countable)/batch_sizes[i]
        
        temp[[idx]]<-tibble(dist_name=c("Normal", "Lognormal", "Beta",
                                        "Normal Sampling",
                                        "Lognormal Sampling",
                                        "Beta Sampling"),
                            sample_size=batch_sizes[i],
                            FinalCount=c(
                                    round(mean(sim_normal, na.rm=TRUE)),
                                    round(mean(sim_lnorm, na.rm=TRUE)),
                                    round(mean(sim_beta, na.rm=TRUE)),
                                    sim_normal_noise,
                                    sim_lnormal_noise,
                                    sim_beta_noise
                            )
                            ) # close assignment to temp[[]]
        idx<-idx+1
      }
    }
  
  # Rearrange the finished data from temp into a proper tibble
  # with columns dist_name, sample_size, and FinalCount
  distribution_means<-dplyr::bind_rows(temp)
  return(distribution_means)
}


