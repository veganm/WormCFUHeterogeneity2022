sim_means_named<-function(input_data, n_reps, batch_sizes, alpha_start=1, beta_start=2, correction_constant=1){
  # A function to generate simulated distributions of means
  # using the normal, lognormal, and Beta distributions (latter is non-unique) 
  # Data generated are average values for simulated "samples"
  # n_reps is the number of values to generate at each batch size
  # batch_sizes is a vector of batch sizes
  # alpha_start and beta_start are initialization for beta distribution fits
  # correction_constant is a multiplier to correct the fraction of a sample used for one measurement
  #   (e.g. for CFU data, measurements from 10 uL spots need a correction of (1000/10)=100 for an original volume of 1 mL)
  # Returns a tibble where observed mean values are stored for each "sample"

  # input_data should be a tibble where each row is one data point, 
  #  and all data points are taken from the same sampling frame (one experiment, sampling run,etc)
  #  The data must have columns:
  #   - Batch (num): Number of individuals/units in a batch. 
  #                 (Batch==1) for individual (unbatched) measurements.
  #   - Count (num): raw counts (e.g. number of colonies counted for CFU data)
  # The data may have columns:
  #   - D (num):  Ten-fold serial dilution. A value of 0 is undiluted; 1 indicates the first ten-fold dilution; etc.
  #               If this column is not present, the code uses the raw counts * correction_constant for all calculations.
  #               If this column is present, the code will check for a column named "CFU" and calculate it if missing.
  #   - CFU (num): CFU counts. Can be calculated from raw counts if dilution factor is given.
  
    
  # load in packages
  pacman::p_load(MASS, tidyverse) 
  
  # If there are individual-sample data
  
  
  # Get summary stats for data
  shared_mean<-mean(mydata, na.rm=TRUE)
  shared_var<-var(mydata, na.rm=TRUE)
  shared_sd<-sqrt(shared_var)
  
  # Parameterize the beta distribution
  # First, range-normalize to [0,1]
  mydata_rn<-(mydata - min(mydata))/(max(mydata)-min(mydata))
  # then fit shape parameters
  beta_fit<-fitdistr(mydata_rn, dbeta, list(shape1=alpha_start, shape2=beta_start))
  
  # Create an empty temporary vector of defined size to hold the data
  temp<-vector("list", length=length(sample_sizes)*n_reps)
  
  # and populate it
  # Establish an index which will be used to store data
  idx<-1
  # Iterate over each of the values in the vector of sample sizes
  for (i in 1:length(sample_sizes)){
    #print(i) # for debugging
    # Within each sample size, execute the indicated number of replicates
    for (j in 1:n_reps){
      #print(j) # for debugging
      # Within each replicate, create a tibble with the mean from each sample
      # and staple it into the temporary vector
      
      # first generate the raw data
      sim_normal<-rnorm(sample_sizes[i], mean=shared_mean, sd=shared_sd)
      sim_lnorm<-rlnorm(sample_sizes[i], meanlog=shared_mean, sdlog=shared_sd)
      sim_beta<-rbeta(sample_sizes[i], shape1 = a, shape2 = b)
      
      temp[[idx]]<-tibble(dist_name=c("Normal", "Lognormal", "Beta"),
                          sample_size=sample_sizes[i],
                          means=c(mean(sim_normal),
                                  mean(sim_lnorm),
                                  mean(sim_beta))
                          )
      )
      idx<-idx+1
    }
  }
  # Rearrange the finished data from temp into a proper tibble
  # with columns dist_name, sample_size, and means
  distribution_means<-dplyr::bind_rows(temp)
  return(distribution_means)
}
