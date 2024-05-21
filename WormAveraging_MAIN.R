pacman::p_load(ggplot2, tidyverse, cowplot, mclust, e1071, ggpubr, readxl, patchwork)
xTextSize<-14

# This code was produced and run with
# R 4.4.0 running in Rstudio 2024.04.0 build 735
# using Tidyverse 2.2.0
# and other packages 
# e1071 1.7-14
# ggplot2 3.5.1
# ggpubr 0.6.0
# mclust 6.1.1
# patchwork 1.2.0
# readxl 1.4.3
# sBIC 0.2.0 (SI)

#############################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FUNCTIONS
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bootOnCounts<-function(n_reps, mydata, batch_sizes=c(1,5,10,20,50), FoldD=10, correction_constant=20){
  # Expects a number of replicates for the bootstrap (n_reps)
  # and a data frame where each row represents one individual (mydata)
  # where the number of colonies counted is in column "Count"
  # the dilution at which these colonies were measured is in column "D",
  #   e.g. D=0 for undiluted sample, D=1 for the first FoldD-fold dilution, etc
  # FoldD(num) is fold dilution in dilution series 
  #  (Default is 10X dilutions, e.g. each step in serial dilution is 1:10 volume)
  #
  ## The dilution correction factor (numeric) is given as "correction_constant"
  # and is the ratio of the original volume and the volume plated/measured
  # e.g. for 10 uL spots and an original volume of 1 mL, correction_constant = 100
  # 
  # Returns a data frame of simulated batch digests
  # with batch sizes in vector batch_sizes, default (1, 5, 10, 20, 50) individuals/batch
  # values reported as inferred CFU/individual and log10(CFU/individual)
  
  # get size of data
  capp<-dim(mydata)[1]
  
  # make someplace to put stuff
  temp<-vector("list", length=length(batch_sizes)*n_reps)
  
  for(i in seq_len(n_reps)){
    batch_data<-rep(0, length(batch_sizes))
    for (j in seq_along(batch_sizes)){
      # Randomly pull samples from data for batching
      idx<-sample(1:capp,batch_sizes[j],replace=TRUE)
      # Assume Poisson count error and generate new counts
      temp_count<-rpois(batch_sizes[j], mydata$Count[idx])
      # Calculate CFU/worm
      batch_data[j]<-mean(correction_constant*temp_count*FoldD^mydata$D[idx], na.rm=TRUE)
    }
    temp[[i]]<-tibble(Batch=batch_sizes,
                      FinalCount=batch_data) 
  }
  
  # unfold data
  dataSet<-dplyr::bind_rows(temp)
  
  dataSet<-dataSet %>%
    mutate(logCFU=log10(FinalCount))
  
  return(dataSet)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bootOnCountsStats<-function(input_data, batch_sizes=c(1,5,10,20,50), nboot=1000, 
                            FoldD=10, correction_constant=20){
  # Function that performs many iterations of bootstrapping on counts
  # Assumes all data sets (indicated by unique values of Run) are independent replicates of one experiment
  # If this is true, results will indicate run-to-run variation and estimated false-positive rates
  # Calls bootOnCounts() for basic functionality.
  # 
  # IF A DATA SET OF INDIVIDUAL-BASED DATA IS PROVIDED
  #   - The data set must contain at least two runs of data, ideally with 10+ data points in each
  #   - Function assumes no subsampling (each individual is essentially a unit-1 subsample)
  #   - Generates resampled data plus Poisson noise for individual and batch-based measurements
  #   - Returns data statistics and t-test and Wilcoxon results of nboot resamples for each pair of samples at each batch size
  #
  # IF A DATA SET OF BATCH-BASED DATA IS PROVIDED
  #   - The data set must contain at least two runs of data, ideally with 10+ data points in each
  #   - Generates resampled data plus Poisson noise using the indicated weighting from input data
  #   - Returns data statistics and t-test and Wilcoxon results of nboot resamples for each pair of samples at each batch size
  # 
  # Takes a data set with columns
  #   Run (num): Each unique index represents one data set
  #   Count (num): Number of raw counts associated with each sample
  #   Batch (int): batch size or weight of each sample.
  #         Batch=1 is the minimum and indicates individual-based sampling.
  #   D (num, optional): Fold dilution at which counts were taken
  #         (e.g. D=0 indicates undiluted sample; D=1 indicates the first FoldD-fold dilution; etc)
  #         ***IF D IS NOT GIVEN, FinalCount MUST ALREADY EXIST IN THE DATA SET***
  #   FinalCount (num, optional): Total inferred counts, based on raw counts and adjusted for dilution and sampling.
  #         Will be calculated if not present.
  # Each row of the data set represents one measurement.
  # 
  # Other inputs:
  #   batch_sizes (int): If suppling individual measurements, a vector of batch sizes is needed for
  #   FoldD (num): Fold dilution in dilution series. Default is 10X dilutions.
  #   correction_constant (num): a multiplier to correct the fraction of a sample used for one measurement to the original volume
  #        (e.g. counts from 10 uL spots and an initial volume of 1 mL require correction (1000/10)=100)
  
  # load the necessary
  pacman::p_load(tidyverse)
  
  ## START
  # Figure out what is in the data set
  # Do we need to correct counts for dilution?
  # Divide by batch size (may be 1) to get per-unit or per-individual numbers
  my_names<-names(input_data) 
  if(length(which(my_names=="FinalCount"))==0){  # if FinalCount does not exist
    if(length(which(my_names=="D"))!=0){ # if we are given a dilution factor
      input_data$FinalCount<-(correction_constant *(input_data$Count*FoldD^input_data$D))/Batch
    } else if(length(which(my_names=="D"))==0){
      stop("Either fold dilution D or FinalCount must be a column in input_data!")
    }
      else{
      input_data$FinalCount<-(correction_constant * input_data$Count)/Batch
    }
  } 
  # now we should have final counts
  
  # Are data individual or batch-based?
  isIndividual<-max(input_data$Batch == 1)
  
  # How many runs of data are in the input data set?
  run_ids<-unique(input_data$Run)
  n_runs<-length(run_ids)
  if(n_runs<2){
    stop("At least two runs of data are required for bootstrapping.")
  }
  
  # initiate looping index
  idx<-1
  
  if (isIndividual){  ## For individual-based data
    # Someplace to put the data generated
    temp<-vector("list", length=length(batch_sizes)*n_runs*n_runs*nboot)
    
    for (m in seq_len(nboot)){  
      for (i in seq_len(n_runs)){
        for (j in seq_len(n_runs)){ # For each pair of runs
            # Obtain two sets of data from the input
            temp1<-input_data %>%
              dplyr::filter(Run==run_ids[i])
            temp2<-input_data %>%
              dplyr::filter(Run==run_ids[j])
            
            # Bootstrap both sets
            Boot1<-bootOnCounts(n_reps=dim(temp1)[1], mydata=temp1, 
                                batch_sizes=batch_sizes, FoldD=FoldD, correction_constant = correction_constant)
            Boot2<-bootOnCounts(n_reps=dim(temp2)[1], mydata=temp2, 
                                batch_sizes=batch_sizes, FoldD=FoldD, correction_constant = correction_constant)
            
            # fix any zeros thrown by resample
            if (sum(is.infinite(Boot1$logCFU))>0){Boot1$logCFU[is.infinite(Boot1$logCFU)]<-0}
            if (sum(is.infinite(Boot2$logCFU))>0){Boot2$logCFU[is.infinite(Boot2$logCFU)]<-0}
    
            # generate stats and carry out tests
            for (k in seq_along(batch_sizes)){
              Boot1s<-Boot1 %>%
                filter(Batch==batch_sizes[k])
              Boot2s<-Boot2 %>%
                filter(Batch==batch_sizes[k])
              Boot.t<-t.test(Boot1s$logCFU, Boot2s$logCFU)
              Boot.w<-wilcox.test(Boot1s$logCFU, Boot2s$logCFU)
              sw1<-shapiro.test(Boot1s$logCFU)
              sw2<-shapiro.test(Boot2s$logCFU)
              
              # store the results
              temp[[idx]]<-tibble(
                Boot=m,
                Batch=batch_sizes[k],
                Run1=i,
                Run2=j,
                meanCFU1=mean(Boot1s$FinalCount, na.rm=TRUE), 
                meanCFU2=mean(Boot2s$FinalCount, na.rm=TRUE),
                varCFU1=var(Boot1s$FinalCount, na.rm=TRUE), 
                varCFU2=var(Boot2s$FinalCount, na.rm=TRUE),
                cvCFU1=sd(Boot1s$FinalCount, na.rm=TRUE)/mean(Boot1s$FinalCount, na.rm=TRUE),
                cvCFU2=sd(Boot2s$FinalCount, na.rm=TRUE)/mean(Boot2s$FinalCount, na.rm=TRUE),
                p_sw1=sw1$p.value, 
                p_sw2=sw2$p.value,
                p_t=Boot.t$p.value,
                p_w=Boot.w$p.value
              )
              idx<-idx+1
              } # finish loop over batch sizes
          }
        } # finish loops over run pairs
    } # finish loop over bootstraps
  } else {# finish individual conditional
    # Someplace to put the data generated
    temp<-vector("list", length=n_runs*n_runs*nboot)
    
    # If we have batched data
    # Can't extrapolate to other batch sizes, but we can generate summary stats for the batching in data
    for (m in seq_len(nboot)){  
      for (i in seq_len(n_runs)){
        for (j in seq_len(n_runs)){ # For each pair of runs
          # Obtain two sets of data from the input
          temp1<-input_data %>%
            dplyr::filter(Run==run_ids[i])
          temp2<-input_data %>%
            dplyr::filter(Run==run_ids[j])
          n1<-dim(temp1)[1] # get sizes
          n2<-dim(temp2)[1]
          # Resamples with batch weights from data
          idx1<-sample(1:n1,n1,replace=TRUE)
          idx2<-sample(1:n2,n2,replace=TRUE)
          # Assume Poisson count error and generate new counts
          temp_count_1<-rpois(n1, temp1$Count[idx1])
          temp_count_2<-rpois(n2, temp2$Count[idx2])
          # Calculate biologically averaged total counts for resampled data
          Boot1s<-(correction_constant*temp_count_1*FoldD^temp1$D[idx1])/temp1$Batch
          Boot2s<-(correction_constant*temp_count_2*FoldD^temp2$D[idx2])/temp2$Batch
          
          Boot.t<-t.test(log10(Boot1s), log10(Boot2s))
          Boot.w<-wilcox.test(Boot1s, Boot2s)
          sw1<-shapiro.test(log10(Boot1s))
          sw2<-shapiro.test(log10(Boot2s))
          
          # store the results
          temp[[idx]]<-tibble(
            Boot=m,
            Run1=i,
            Run2=j,
            meanCFU1=mean(Boot1s, na.rm=TRUE), 
            meanCFU2=mean(Boot2s, na.rm=TRUE),
            varCFU1=var(Boot1s, na.rm=TRUE), 
            varCFU2=var(Boot2s, na.rm=TRUE),
            cvCFU1=sd(Boot1s, na.rm=TRUE)/mean(Boot1s, na.rm=TRUE),
            cvCFU2=sd(Boot2s, na.rm=TRUE)/mean(Boot2s, na.rm=TRUE),
            p_sw1=sw1$p.value, 
            p_sw2=sw2$p.value,
            p_t=Boot.t$p.value,
            p_w=Boot.w$p.value
          )
          idx<-idx+1
        }
      } # finish loop over run pairs
    } # finish loop over bootstraps
  } # end else loop
  
  # unfold and return results
  dataSet<-dplyr::bind_rows(temp)
  return(dataSet)
}

###############################################################################
#  DATA ANALYSIS AND FIGURES
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~   Figure 1A
# Real data from batch digests 2022-2-18,
# N2 worms with Staphylococcus aureus

SaSeCountAll<-read_xlsx("SaSeCount.xlsx")  # full data file

BatchDigests<-SaSeCountAll %>%
  dplyr::filter(Condition=="SA" & Run==6) # run with single worms and batch digests

BatchDigests$Batch<-as.factor(BatchDigests$Batch)
pBatchSA<-BatchDigests %>%
  ggplot(aes(x=Batch, y=logCFU, color=Batch)) + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + theme_classic() + 
  scale_color_viridis_d(option="magma", end=0.8)+
  theme(text=element_text(size=xTextSize), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=xTextSize),
        axis.title.y = element_text(size=xTextSize),
        axis.text.y = element_text(size=xTextSize-1),
        legend.text = element_text(size=xTextSize-1),
        legend.position.inside=c(0.9,0.3)) + 
  labs(title=expression(paste(italic("S. aureus"), " Newman")), y=expression(log[10](CFU/Worm)))
pBatchSA

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Now Salmonella enterica LT2-GFP and S. aureus-GFP in single worm digests (N2)

SaSeCount<-SaSeCountAll %>%
  dplyr::filter(Batch==1) # single worms only

# Note that runs 1-3 of S. aureus have a different TOD
# with dilutions 1-3 plated; 
# runs 4-6 start at dilution 0
SaSeCount %>%
  ggplot(aes(x=factor(Run), y=logCFU, color=factor(Run))) + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + 
  ylim(-0.1,6)+ theme_classic() + 
  theme(
    text=element_text(size=14), 
    axis.title.x = element_blank(), 
    plot.title=element_text(hjust=0.5,size=14),
    legend.position = "none") +
  facet_wrap(vars(Condition), scales="free_x")

# summary statistics
SaSeCount %>%
  group_by(Condition, Run) %>%
  summarise(count=n(),
            medianFinalCount=median(FinalCount, na.rm=TRUE),
            meanFinalCount=mean(FinalCount, na.rm=TRUE),
            sdFinalCount=sd(FinalCount, na.rm=TRUE),
            medianLogCFU=median(logCFU, na.rm=TRUE),
            meanLogCFU=mean(logCFU, na.rm=TRUE),
            sdLogCFU=sd(logCFU, na.rm = TRUE)
            )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#take out the zeros
SaSeCount2<-SaSeCount
SaSeCount2$Count[SaSeCount2$Count==0]<-NA
SaSeCount2<-SaSeCount2[complete.cases(SaSeCount2),]

SaSeCount2 %>%
	ggplot(aes(x=factor(Run), y=logCFU, color=factor(Run))) + 
	geom_jitter(shape=16, position=position_jitter(0.05)) +
	geom_violin(fill=NA) + 
	ylim(-0.1,6)+ theme_classic() + 
	theme(
  	text=element_text(size=14), 
  	#axis.title.x = element_blank(), 
	  plot.title=element_text(hjust=0.5,size=14),
  	legend.title = element_blank()) +
  facet_wrap(vars(Condition), scales="free_x")+
  labs(y=expression(log[10](CFU/worm)), x="Replicate")
	
# summary statistics
SaSeCount2 %>%
  group_by(Condition, Run) %>%
  summarise(count=n(),
            medianFinalCount=median(FinalCount, na.rm=TRUE),
            meanFinalCount=mean(FinalCount, na.rm=TRUE),
            sdFinalCount=sd(FinalCount, na.rm=TRUE),
            medianLogCFU=median(logCFU, na.rm=TRUE),
            meanLogCFU=mean(logCFU, na.rm=TRUE),
            sdLogCFU=sd(logCFU, na.rm = TRUE)
  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#        GENERATE BOOTSTRAPPED DATA FROM S. ENTERICA COLONIZATION
#
# here we will call the function bootOnCounts(n_reps, mydata, batch_sizes=c(1,5,10,20,50), FoldD=10, correction_constant=20)
# to generate simulated data and make comparisons

temp1<-SaSeCount %>%
  dplyr::filter(Condition=="SE" & Run=="3")
temp2<-SaSeCount %>%
  dplyr::filter(Condition=="SE" & Run=="2")

SeBoot1<-bootOnCounts(n_reps=dim(temp1)[1], mydata=temp1)
SeBoot2<-bootOnCounts(n_reps=dim(temp2)[1], mydata=temp2)
if (sum(is.infinite(SeBoot1$logCFU))>0){SeBoot1$logCFU[is.infinite(SeBoot1$logCFU)]<-0}
if (sum(is.infinite(SeBoot2$logCFU))>0){SeBoot2$logCFU[is.infinite(SeBoot2$logCFU)]<-0}

SeBoot1$Run<-as.factor("Run3")
SeBoot2$Run<-as.factor("Run2")
jointSeBoot<-rbind(SeBoot1, SeBoot2)

pJointSEBoot<-jointSeBoot %>%
	ggplot(aes(x=factor(Run), y=logCFU, color=factor(Run))) + 
	geom_jitter(shape=16, position=position_jitter(0.05)) +
	geom_violin(fill=NA) + 
  geom_hline(yintercept=log10(20), color="black", lty="dashed")+
	ylim(-0.1,6.2) + 
  theme_classic() +
  scale_color_viridis_d(begin=0.3, end=0.8) +
	theme(
	  text=element_text(size=xTextSize), 
	  axis.title.x = element_blank(),
	  axis.text.x = element_blank(),
	  axis.text.y = element_text(size=xTextSize-1),
	  axis.title.y = element_text(size=xTextSize),
	  plot.title=element_text(hjust=0.5,size=xTextSize),
	  #legend.position.inside=c(0.9,0.3),
	  legend.title = element_blank(),
	  legend.text = element_text(size=xTextSize-1)
	  )+
  facet_wrap(vars(Batch), ncol=5)+
  labs(title=expression(paste(italic("S. enterica"), " Simulated Batch Digests")), 
       y=expression(log[10](CFU/Worm)), x="Run")+
  stat_compare_means(method="t.test", label.y = 6.2, size=3.2)+
  stat_compare_means(method="wilcox.test", label.y=5.9, size=3.2)
pJointSEBoot

#rm(SeBoot1, SeBoot2, temp1, temp2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# same SE bootstrap with zeros set to TOD
temp<-dplyr::filter(SaSeCount, Condition=="SE")
idx<-which(temp$Count==0)
temp$Count[idx]<-1
temp$FinalCount[idx]<-20
temp$CFUPerWorm[idx]<-20
temp$logCFU[idx]<-log10(20)

temp1A<-temp %>%
  dplyr::filter(Run=="3")
temp2A<-temp %>%
  dplyr::filter(Run=="2")

SeBoot1A<-bootOnCounts(n_reps=dim(temp1A)[1], mydata=temp1A)
SeBoot2A<-bootOnCounts(n_reps=dim(temp2A)[1], mydata=temp2A)
SeBoot1A$Run<-as.factor("Run3")
SeBoot2A$Run<-as.factor("Run2")
jointSeBoot_nozeros<-rbind(SeBoot1A, SeBoot2A)
idx<-which(!is.finite(jointSeBoot_nozeros$logCFU))
jointSeBoot_nozeros$logCFU[idx]<-0

pJointSEBoot_nozeros<-jointSeBoot_nozeros %>%
  dplyr::filter(Batch<20) %>% # sufficient to make the point
  ggplot(aes(x=factor(Run), y=logCFU, color=factor(Run))) + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + 
  ylim(-0.1,6.2) + 
  theme_classic() +
  scale_color_viridis_d(begin=0.3, end=0.8) +
  theme(
    text=element_text(size=xTextSize), 
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=xTextSize-1),
    axis.title.y = element_text(size=xTextSize),
    plot.title=element_text(hjust=0.5,size=xTextSize),
    #legend.position.inside=c(0.9,0.3),
    legend.title = element_blank(),
    legend.text = element_text(size=xTextSize-1)
  )+
  facet_wrap(vars(Batch), ncol=5)+
  labs(title="Zeros removed", x="Run", y=expression(log[10](CFU/Worm)))+
  stat_compare_means(method="t.test", label.y = 6.2, size=3.2)+
  stat_compare_means(method="wilcox.test", label.y=5.9, size=3.2)
pJointSEBoot_nozeros

#rm(SeBoot1A, SeBoot2A, temp1A, temp2A)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# and with extra zeros
# implemented as higher TOD, call it 10^3
temp<-dplyr::filter(SaSeCount, Condition=="SE")
idx<-which(temp$logCFU<3)
temp$Count[idx]<-0
temp$D[idx]<-0
temp$FinalCount[idx]<-0
temp$CFUPerWorm[idx]<-0
temp$logCFU[idx]<-0

temp1B<-temp %>%
  dplyr::filter(Run=="3")
temp2B<-temp %>%
  dplyr::filter(Run=="2")

SeBoot1B<-bootOnCounts(n_reps=dim(temp1B)[1], mydata=temp1B)
SeBoot2B<-bootOnCounts(n_reps=dim(temp2B)[1], mydata=temp2B)
SeBoot1B$Run<-as.factor("Run3")
SeBoot2B$Run<-as.factor("Run2")
jointSeBoot_zeros3<-rbind(SeBoot1B, SeBoot2B)

pJointSEBoot_zeros3<-jointSeBoot_zeros3 %>%
  dplyr::filter(Batch<20) %>% # sufficient to make the point
  ggplot(aes(x=factor(Run), y=logCFU, color=factor(Run))) + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + 
  geom_hline(yintercept=3, color="black", lty="dashed")+
  ylim(-0.1,6.2) + 
  theme_classic() +
  scale_color_viridis_d(begin=0.3, end=0.8) +
  theme(
    text=element_text(size=xTextSize), 
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=xTextSize-1),
    axis.title.y = element_text(size=xTextSize),
    plot.title=element_text(hjust=0.5,size=xTextSize),
    #legend.position.inside=c(0.9,0.3),
    legend.title = element_blank(),
    legend.text = element_text(size=xTextSize-1)
  )+
  facet_wrap(vars(Batch), ncol=5)+
  labs(title="Zero-enriched", x="Run", y=expression(log[10](CFU/Worm)))+
  stat_compare_means(method="t.test", label.y = 6.2, size=3.2)+
  stat_compare_means(method="wilcox.test", label.y=5.9, size=3.2)
pJointSEBoot_zeros3

#rm(SeBoot1B, SeBoot2B, temp1B, temp2B, temp)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~    Build Figure 1
# Plot out SE bootstraps together with real data from S. aureus batch size experiment
#Fig1AB<-plot_grid(pBatchSA, pJointSEBoot, labels="AUTO", ncol=2, rel_widths=c(1,3))
#Fig1CD<-plot_grid(pJointSEBoot_nozeros, pJointSEBoot_zeros3, labels=c("C", "D"), ncol=2)
#plot_grid(Fig1AB, Fig1CD, ncol=1)

(pBatchSA + pJointSEBoot + plot_layout(ncol=2, widths=c(1,3))) /
  (pJointSEBoot_nozeros + pJointSEBoot_zeros3 + plot_layout(ncol=2)) +
  plot_layout(guides="collect")+
  plot_annotation(tag_levels = 'A')

ggsave("Fig1_pSaBatchJointSeBoot.png", width=12, height=8, units="in", dpi=400)

# and some summary statistics for the record
mean(BatchDigests$CFU[BatchDigests$Species=="SA" & BatchDigests$Batch==1 ])
median(BatchDigests$CFU[BatchDigests$Species=="SA" & BatchDigests$Batch==1 ])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~                 FIGURE 2
#~~~~    WITHIN AND BETWEEN REPLICATE HETEROGENEITY OF CFUS
#
# First, pull in the community colonization data set (10.1128/mSystems.00608-20)
# and the single worm data from NRRL1 evolution paper (10.1093/femsec/fiac115)
# so we can describe run to run variation
# just need the total counts for this


WormReps<-read_csv("WormCFUReplicates.csv")
names(WormReps)
WormReps$Condition[WormReps$Condition=="decc1"]<-"dec1" # fixing a strain name

WormReps<-WormReps[complete.cases(WormReps),]

# filter the original data set to remove anything with <18 measurements per strain on a given day
conditions<-unique(WormReps$Condition)
WormReps.temp<-WormReps
for (i in 1:length(conditions)){
	tempdata<-subset(WormReps.temp, Condition==conditions[i])
	temp2<-unique(tempdata$Date)
	for (j in 1:length(temp2)){
		mycount<-(length(WormReps.temp$logCFU[WormReps.temp$Condition==conditions[i] & 
		                                       WormReps.temp$Date==temp2[j]])) #check # worms in data set
		if(mycount<18) { #if too few worms
		  WormReps.temp$logCFU[WormReps.temp$Condition==conditions[i] & 
		                        WormReps.temp$Date==temp2[j]]<-NA #flag for removal
		}
	}
}
WormReps.temp$logCFU[is.infinite(WormReps.temp$logCFU)]<-NA
WormReps.big<-WormReps.temp[complete.cases(WormReps.temp),]
rm(WormReps.temp)

### let's see what's left
conditions.big<-unique(WormReps.big$Condition)
temp<-integer()
for (i in 1:length(conditions.big)){
	temp[i]<-length(unique(WormReps.big$Rep[WormReps.big$Condition==conditions.big[i]]))}
temp # number of replicates each experimental condition
temp2<-unique(WormReps.big$Condition) #lost one host type

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~  Building the filtered data set
# Let's filter out anything that doesn't have 2 or more data sets left
# these have enough runs to be worth it
condition_list<-temp2[which(temp>1)]
WormReps.big<-filter(WormReps.big, Condition %in% condition_list)
  
# Plot out individual replicates by host genotype
# first create new facet labels
condition.labs<-c("AU37", "daf-16", "daf-2", "dbl-1", "dec-1", "exp-1", "N2", "MO", "vhp-1")
names(condition.labs)<-c("AU37", "daf16", "daf2", "dbl1", "dec1", "exp1", "N2", "MO.I2.9", "vhp1")

pWormReps.big.byCondition<-WormReps.big %>%
  ggplot(aes(x=factor(Rep), y=logCFU)) + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + 	
  ylim(1,6) +	
  theme_classic() + 
  theme( 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size=xTextSize),
        axis.text.y = element_text(size=xTextSize-1),
        legend.position="none", 
        plot.title=element_text(hjust=0.5, size=xTextSize),
        strip.text.x = element_text(size=xTextSize-2, face="italic")) +
  labs(title="Minimal community colonization, by replicate", y=expression(log[10](CFU/Worm))) +
  facet_wrap(~Condition, scales="free_x", ncol=4, labeller=labeller(Condition=condition.labs))
  
pWormReps.big.byCondition

#ggsave("pAll8BigLogCFUHostByDay.png", width=12, height=9, units="in", dpi=300)

# Plot out the log-transformed CFU data by host genotype, replicates combined
WormReps.big %>%
	ggplot(aes(x=as.factor(Condition), y=logCFU))+
	geom_jitter(shape=16, position=position_jitter(0.05)) +
	geom_violin(fill=NA) + 	ylim(1,6) +	theme_classic() + 
	theme(text=element_text(size=14), axis.title.x = element_blank(), legend.position="none", plot.title=element_text(hjust=0.5, size=16)) +
	labs(title="Totals", y=expression(log[10](CFU/Worm)))

rm(temp, temp2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# OK now we are going to take a WHOLE BUNCH of summary stats and see what we get

AllCounts<-as_tibble(subset(WormReps.big, select=c("Condition", "Rep", "CFU", "logCFU")))
glimpse(AllCounts)
AllCounts$Rep<-as.factor(AllCounts$Rep)

# Incorporate the pathogen colonization data: prep for rbind

SaSeCount2_subset<- SaSeCount2 %>%
  subset(select=c("Condition", "Run", "CFU", "logCFU")) %>%
  dplyr::rename(Rep=Run)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Merge all the single-worm CFU data
AllCounts<-rbind(AllCounts, SaSeCount2_subset)
glimpse(AllCounts)

# Get a list of all experiments
conditions.big<-unique(AllCounts$Condition)

# set up the data frame to hold our calculations
AllCountsStats.r2r <- data.frame(Condition=character(),
                 Rep=character(),
                 maxCFU=double(),
                 minCFU=double(),
                 q90=double(),
                 q75=double(),
                 q50=double(),
                 q25=double(),
                 q10=double(),
                 meanCFU=double(),
                 varCFU=double(),
                 skewCFU=double(), # from package e1071
		            	kurtCFU=double(),
		            	Q1=double(), # Hogg's Q1, a measure of skewness
		            	Q2=double(), # Hogg's Q2, a measure of tail weight
                 stringsAsFactors=TRUE)

# loop over each condition (host genotype or colonizing bacteria) in the data set
# then each replicate (day) for that condition
# and record

for (i in seq_along(conditions.big)){
	temp<-unique(AllCounts$Rep[AllCounts$Condition==conditions.big[i]])
#	print(conditions.big[i]) # for debugging
	for (j in seq_along(temp)){
#		print(temp[j]) # for debugging
		tempCFU<-sort(AllCounts$CFU[AllCounts$Condition==conditions.big[i] & AllCounts$Rep==temp[j]])
		myq<-as.double(quantile(tempCFU, probs=c(0.1, 0.25, 0.5, 0.75, 0.9)))
		# get values for Hogg's calculations
		mylen<-length(tempCFU)
		idx05<-ceiling(mylen/20)
		if(idx05>1){L05<-mean(tempCFU[1:idx05])} else {L05<-tempCFU[1]}
		idx25<-ceiling(mylen/4)
		idx75<-floor(mylen*0.75)
		idx95<-floor(mylen*0.95)
		if(idx95<mylen){U05<-mean(tempCFU[- (1:(idx95-1))])} else {U05<-tail(tempCFU,1)}
		M5<-mean(tempCFU[idx25:idx75])
		idx50<-round(mylen/2)
		L5<-mean(tempCFU[1:idx50])
		U5<-mean(tempCFU[- (1:idx50)])
		#build the new entry
		AllCountsStats.r2r<-rbind(AllCountsStats.r2r, data.frame(Condition=conditions.big[i],
	  	Rep=temp[j],
	  	maxCFU=max(tempCFU),
	  	minCFU=min(tempCFU),
  		q90=myq[5],	
  		q75=myq[4],
	  	q50=myq[3],
	  	q25=myq[2],
	  	q10=myq[1],
	  	meanCFU=mean(tempCFU),
	  	varCFU=var(tempCFU),
	  	skewCFU=skewness(tempCFU),
	  	kurtCFU=kurtosis(tempCFU),
	  	Q1=(U05-M5)/(M5-L05),
	  	Q2=(U05-L05)/(U5-L5)))
	}	
}
dim(AllCountsStats.r2r)

# add in CV after the fact
AllCountsStats.r2r$cvCFU<-sqrt(AllCountsStats.r2r$varCFU)/AllCountsStats.r2r$meanCFU

#write summary stats to table
write.table(AllCountsStats.r2r, "WormDataRunToRun.txt", sep="\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's plot out the data and statistics by experimental condition
# first create new facet labels
condition.labs<-c("N2","AU37", "daf-16", "daf-2", "dbl-1", "dec-1", "exp-1", "vhp-1", "M. oxydans",  "S. aureus", "S. enterica")
names(condition.labs)<-c("N2","AU37", "daf16", "daf2", "dbl1", "dec1", "exp1", "vhp1", "MO.I2.9",  "SA", "SE")

# we lost a lot of replicates - rescale for plotting
unique(AllCounts$Rep)
unique(AllCounts$Rep[AllCounts$Condition=="N2"]) #largest number of reps
# fix N2
AllCounts$Rep[AllCounts$Condition=="N2" & AllCounts$Rep=="10"]<-"4"
AllCounts$Rep[AllCounts$Condition=="N2" & AllCounts$Rep=="11"]<-"5"

pAllCounts<-AllCounts %>%
  mutate(Condition=factor(Condition, levels=c("N2","AU37", "daf16", "daf2", "dbl1", "dec1", "exp1", "vhp1", "MO.I2.9", "SA", "SE"))) %>%
  ggplot(aes(x=Rep, y=logCFU, color=Rep)) + 
  geom_point(shape=16, position=position_jitterdodge(0.05)) +
  geom_violin(fill=NA) + 
  theme_classic() + 
  scale_color_viridis_d(option="inferno", end=0.85)+
  facet_wrap(~Condition, scales="free_x", labeller=labeller(Condition=condition.labs))+
  theme( 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=xTextSize), 
        axis.text.y = element_text(size=xTextSize-1), 
        #plot.title=element_text(hjust=0.5, size=14),
        legend.position = "none",
        strip.text.x = element_text(size=xTextSize-2, face="italic"))+
  labs(y=expression(log[10](CFU/Worm)))
pAllCounts
#~~~~~ this will be FIGURE 2A
#ggsave("pAll8BigSaSeCFUHostByDay.png", width=12, height=10, units="in", dpi=300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~ Plot out the moments by condition
names(AllCountsStats.r2r)

#some of these stats are hard to see on linear axes; mutate to plot
AllCountStats.r2r<-AllCountsStats.r2r %>%
  mutate(logmeanCFU=log10(meanCFU),
         logq50=log10(q50))

pAllCountStats<- AllCountStats.r2r %>% 
  select(Condition, logmeanCFU, skewCFU, kurtCFU, Q1, Q2, cvCFU, logq50) %>%
  pivot_longer(., cols = c(logmeanCFU, skewCFU, kurtCFU, Q1, Q2, cvCFU, logq50), names_to = "Var", values_to = "Value") %>%
  ggplot(aes(x=Condition, y=Value, color=Condition))+
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_boxplot(fill=NA) + 	
  theme_classic() + 
  scale_color_viridis_d(option="magma", begin=0.3, end=0.85)+
  theme(text=element_text(size=14), 
        axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x = element_blank(), 
        legend.position="none", 
        plot.title=element_text(hjust=0.5, size=16)) +
  facet_wrap(vars(Var), scales="free_y", ncol=4)
#ggsave("pRawMomentsRunToRunWormCFU.png", width=14, height=8, units="in", dpi=300)
pAllCountStats

### cool. Now let's do the pairwise comparisons
npairs<-dim(AllCountStats.r2r)[1]

# we have an alphabetical order problem let's fix it
AllCountStats.r2r.t<-as_tibble(AllCountStats.r2r)
AllCountStats.r2r.t<-arrange(AllCountStats.r2r, Condition)
hostlist<-unique(AllCountStats.r2r.t$Condition)
hostcount<-as.data.frame(table(AllCountStats.r2r.t$Condition))

host1<-rep(as.character(AllCountStats.r2r.t$Condition[1]), npairs-1)
host2<-AllCountStats.r2r.t$Condition[-1]
for (i in 2:(npairs-1)){
	host1<-c(host1, rep(as.character(AllCountStats.r2r.t$Condition[i]), npairs-i))
	idx<-seq(1,i)
	host2<-c(host2, AllCountStats.r2r.t$Condition[-idx])
}

names(AllCountStats.r2r.t)

#~~~~~~~~~~~~~~~~~~~
# now let's look at the differences in statistics

AllCountStats.r2r.dist <- tibble(meanCFUdist=numeric(),
                                 meanCFUCohenD=numeric(),
                                 mediandist=numeric(),
                                 cvdist=numeric(),
                                 skewdist=numeric(),
                                 kurtdist=numeric(),
                                 Q1dist=numeric(), # a measure of skewness
                                 Q2dist=numeric()) # a measure of tail weight

#~~~~~~~~~   Calculate all pairwise distances and store
#npairs<-3 #for testing
for(i in 1:(npairs-1)){
	#print(i)
	for (j in (i+1):npairs){
		#print(j)
	  AllCountStats.r2r.dist<-bind_rows(AllCountStats.r2r.dist, tibble(
		meanCFUdist=abs(AllCountStats.r2r.t$meanCFU[i]-AllCountStats.r2r.t$meanCFU[j])/(AllCountStats.r2r.t$meanCFU[i]+AllCountStats.r2r.t$meanCFU[j]),
	  meanCFUCohenD=abs(AllCountStats.r2r.t$meanCFU[i]-AllCountStats.r2r.t$meanCFU[j])/sqrt((AllCountStats.r2r.t$varCFU[i]+AllCountStats.r2r.t$varCFU[j])/2),
		mediandist=abs(AllCountStats.r2r.t$q50[i]-AllCountStats.r2r.t$q50[j])/(AllCountStats.r2r.t$q50[i]+AllCountStats.r2r.t$q50[j]),
		cvdist=abs(AllCountStats.r2r.t$cvCFU[i]-AllCountStats.r2r.t$cvCFU[j]),
		skewdist=abs(AllCountStats.r2r.t$skewCFU[i]-AllCountStats.r2r.t$skewCFU[j]),
		kurtdist=abs(AllCountStats.r2r.t$kurtCFU[i]-AllCountStats.r2r.t$kurtCFU[j]),
		Q1dist=abs(AllCountStats.r2r.t$Q1[i]-AllCountStats.r2r.t$Q1[j]),
		Q2dist=abs(AllCountStats.r2r.t$Q2[i]-AllCountStats.r2r.t$Q2[j])
		)
	  )
	}
}

# add the condition labels
AllCountStats.r2r.dist$host1<-host1 
AllCountStats.r2r.dist$host2<-host2
names(AllCountStats.r2r.dist)
idx<-which(host1==host2)
idx
AllCountStats.r2r.dist$same<-"Different"
AllCountStats.r2r.dist$same[idx]<-"Same"
AllCountStats.r2r.dist$same<-as.factor(AllCountStats.r2r.dist$same)

# sanity check 
AllCountStats.r2r.dist
AllCountStats.r2r.dist[idx,]

# and plot out comparisons of distances, same vs. different condition
my.labs<-c("Mean CFU", "Cohen's D", "Median", "CV", "Q1", "Skewness", "Q2", "Kurtosis")
names(my.labs)=c("meanCFUdist", "meanCFUCohenD", "mediandist", "cvdist", "Q1dist","skewdist", "Q2dist", "kurtdist")
pAllCountStats.r2r<-AllCountStats.r2r.dist %>%
  pivot_longer(., cols = c(meanCFUCohenD, meanCFUdist, mediandist, cvdist, Q1dist, skewdist, Q2dist, kurtdist), 
               names_to = "Var", values_to = "Value") %>%
  mutate(Var=factor(Var, 
                    levels=c("meanCFUdist", "meanCFUCohenD", "mediandist", "cvdist", "Q1dist","skewdist", "Q2dist", "kurtdist")) ) %>%
  ggplot(aes(x=same, y=Value, color=same))+
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_boxplot(fill=NA) +	
  theme_classic() + 
  scale_color_viridis_d(option="magma", begin=0.3, end=0.7)+
  theme( 
    axis.title.x = element_blank(), 
    axis.text.x = element_text(size=xTextSize-1),
    axis.title.y = element_text(size=xTextSize), 
    axis.text.y = element_text(size=xTextSize-1), 
    legend.position="none", 
      plot.title=element_text(hjust=0.5, size=xTextSize)) +
  labs(title="", y="Distance between pairs")+
  facet_wrap(vars(Var), scales="free_y", ncol=2, labeller=labeller(Var=my.labs))+
  stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),
                         label.x.npc = 0.5, label.y.npc=0.9, size=4)
pAllCountStats.r2r
#ggsave("pAll8SelectMomentsSummaryDistance.png", width=9, height=12, units="in", dpi=400)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~   FIGURE 2 ASSEMBLY
plot_grid(pAllCounts, pAllCountStats.r2r, rel_widths=c(1, 1.2), labels="AUTO")
ggsave("Fig2_AllCounts_logCFU_MomentDist.png", width=15, height=9, units="in", dpi=400)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#########       FIGURE 3          ##################################
#####        FALSE NEGATIVES
#####        with real data 
# Pairs 
#(N2 6/21/19 vs dbl-1 6/26/19)
#(SE replicate 2 vs. N2 adults 2/13/19)
# First have to import the raw data for the multispecies colonization

jointBootRawData<-read.table("jointBootFalseNegativeRawData.txt", header=TRUE)

#########  First set (N2 6/21/19 vs dbl-1 6/26/19 multispp) ######
unique(jointBootRawData$Date)
temp2A<-jointBootRawData %>%
  filter(Host=="N2" & Date=="6/21/2019")
temp1A<-jointBootRawData %>%
  filter(Host=="dbl1")
mean(temp2A$CFU)
sd(temp2A$CFU)/mean(temp2A$CFU)
skewness(temp2A$CFU)
mean(temp1A$CFU)
sd(temp1A$CFU)/mean(temp1A$CFU)
skewness(temp1A$CFU)

# here we will call the boot function (bootOnCounts)
# to generate simulated data and make comparisons
# 22 data points to match the smallest data set (N2)
# batch 50 is pretty large for these data so we'll stop at 20

tempBoot1A<-bootOnCounts(n_reps=22, mydata=temp1A, correction_constant=2)
tempBoot2A<-bootOnCounts(n_reps=22, mydata=temp2A, correction_constant=2)
tempBoot1A$Host<-as.factor("dbl1")
tempBoot2A$Host<-as.factor("N2")
jointTempBootA<-rbind(tempBoot1A, tempBoot2A)
#names(jointTempBootA)

pjointBootA<-jointTempBootA  %>%
  ggplot(aes(x=Host, y=logCFU, color=Host)) + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + 
  ylim(2,6)+ 
  theme_classic() + 
  scale_color_viridis_d(begin=0.3, end=0.7)+
  theme( 
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size=xTextSize), 
    axis.text.y = element_text(size=xTextSize-1), 
    axis.text.x = element_blank(),
    plot.title=element_text(hjust=0.5,size=14),
    legend.position="none") +
  labs(title=expression(paste("N2 vs ", italic(dbl-1))), y=expression(log[10](CFU/worm)))+
  facet_wrap(vars(Batch), nrow = 1)+
  stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),
                     label.x = 1, size=4)
pjointBootA

# some tests for normality just b/c
shapiro.test(tempBoot1A$logCFU[tempBoot1A$Batch==1])
shapiro.test(tempBoot1A$logCFU[tempBoot1A$Batch==5])
shapiro.test(tempBoot1A$logCFU[tempBoot1A$Batch==10])
shapiro.test(tempBoot1A$logCFU[tempBoot1A$Batch==20])
shapiro.test(tempBoot1A$logCFU[tempBoot1A$Batch==50])
shapiro.test(tempBoot2A$logCFU[tempBoot2A$Batch==1])
shapiro.test(tempBoot2A$logCFU[tempBoot2A$Batch==5])
shapiro.test(tempBoot2A$logCFU[tempBoot2A$Batch==10])
shapiro.test(tempBoot2A$logCFU[tempBoot2A$Batch==20])
shapiro.test(tempBoot2A$logCFU[tempBoot2A$Batch==50])

#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######   Second set (SE vs N2 Multi 2/13/19)  ########

temp1<-jointBootRawData %>%
  subset(Host=="N2" & Date=="2/13/2019") %>%
  select(Host, Count, D, CFU) %>%
  rename(Condition=Host)
glimpse(temp1)
temp2<-SaSeCount2 %>%
  filter(Condition =="SE") %>%
  subset(Rep=="2", select=c("Condition", "Count", "D", "CFU"))
glimpse(temp2)

mean(temp1$CFU)
sd(temp1$CFU)/mean(temp1$CFU)
skewness(temp1$CFU)
mean(temp2$CFU)
sd(temp2$CFU)/mean(temp2$CFU)
skewness(temp2$CFU)

# Boostrap data
tempBoot1B<-bootOnCounts(n_reps=24, mydata=temp1, correction_constant=2)
tempBoot2B<-bootOnCounts(n_reps=24, mydata=temp2, correction_constant=20)
tempBoot1B$Condition<-as.factor("Multi")
tempBoot2B$Condition<-as.factor("SE")
jointTempBootB<-rbind(tempBoot1B, tempBoot2B)
#glimpse(jointTempBootB)

pJointBootB<- jointTempBootB  %>%
  ggplot(aes(x=Condition, y=logCFU, color=Condition)) + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + 
  ylim(2,6)+ 
  theme_classic() + 
  scale_color_viridis_d(begin=0.3, end=0.7)+
  theme( 
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size=xTextSize), 
    axis.text.y = element_text(size=xTextSize-1), 
    axis.text.x = element_blank(),
    plot.title=element_text(hjust=0.5,size=14),
    legend.position="none") +
  labs(title=expression(paste("Multi-species vs ", italic("S. enterica"), "in N2 hosts")), y=expression(log[10](CFU/worm)))+
  facet_wrap(vars(Batch), nrow = 1)+
  stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),
                     label.x = 1, size=4)
pJointBootB

plot_grid(pjointBootA, pJointBootB, labels="AUTO", ncol=1)

shapiro.test(tempBoot1B$logCFU[tempBoot1B$Batch==1])
shapiro.test(tempBoot1B$logCFU[tempBoot1B$Batch==5])
shapiro.test(tempBoot1B$logCFU[tempBoot1B$Batch==10])
shapiro.test(tempBoot1B$logCFU[tempBoot1B$Batch==20])
shapiro.test(tempBoot1B$logCFU[tempBoot1B$Batch==50])
shapiro.test(tempBoot2B$logCFU[tempBoot2B$Batch==5])
shapiro.test(tempBoot2B$logCFU[tempBoot2B$Batch==10])

#ggsave("pjointSeBoot.svg", width=9, height=4, units="in", dpi=300)
ggsave("Fig3_pjointFalseNegativeBoot.png", width=10, height=7, units="in", dpi=300)

# let's do this a bunch of times
# pair 1
reps<-10000
wpv5<-numeric(reps)
wpv10<-numeric(reps)
wpv20<-numeric(reps)
wpv50<-numeric(reps)

for (j in 1:reps){
  tempBoot1<-bootOnCounts(n_reps=22, mydata=temp1A, correction_constant=2)
  tempBoot2<-bootOnCounts(n_reps=22, mydata=temp2A, correction_constant=2)
  wtest5<-wilcox.test(tempBoot1$CFU[tempBoot1$Batch==5], tempBoot2$CFU[tempBoot2$Batch==5])
  wtest10<-wilcox.test(tempBoot1$CFU[tempBoot1$Batch==10], tempBoot2$CFU[tempBoot2$Batch==10])
  wtest20<-wilcox.test(tempBoot1$CFU[tempBoot1$Batch==20], tempBoot2$CFU[tempBoot2$Batch==20])
  wtest50<-wilcox.test(tempBoot1$CFU[tempBoot1$Batch==50], tempBoot2$CFU[tempBoot2$Batch==50])
  wpv5[j]<-wtest5$p.value
  wpv10[j]<-wtest10$p.value
  wpv20[j]<-wtest20$p.value
  wpv50[j]<-wtest50$p.value
}

w.pvals.A<-c(sum(wpv5<0.05)/reps, sum(wpv10<0.05)/reps, 
             sum(wpv20<0.05)/reps, sum(wpv50<0.05)/reps)
w.pvals.A

#~~~~~~~~~~~~~~~~~~~~~~~~
# pair 2
wpv5<-numeric(reps)
wpv10<-numeric(reps)
wpv20<-numeric(reps)
wpv50<-numeric(reps)

for (j in 1:reps){
  tempBoot1<-bootOnCounts(n_reps=24, mydata=temp1, correction_constant=2)
  tempBoot2<-bootOnCounts(n_reps=24, mydata=temp2, correction_constant=20)
  wtest5<-wilcox.test(tempBoot1$CFU[tempBoot1$Batch==5], tempBoot2$CFU[tempBoot2$Batch==5])
  wtest10<-wilcox.test(tempBoot1$CFU[tempBoot1$Batch==10], tempBoot2$CFU[tempBoot2$Batch==10])
  wtest20<-wilcox.test(tempBoot1$CFU[tempBoot1$Batch==20], tempBoot2$CFU[tempBoot2$Batch==20])
  wtest50<-wilcox.test(tempBoot1$CFU[tempBoot1$Batch==50], tempBoot2$CFU[tempBoot2$Batch==50])
  wpv5[j]<-wtest5$p.value
  wpv10[j]<-wtest10$p.value
  wpv20[j]<-wtest20$p.value
  wpv50[j]<-wtest50$p.value
}

w.pvals.B<-c(sum(wpv5<0.05)/reps, sum(wpv10<0.05)/reps, 
             sum(wpv20<0.05)/reps, sum(wpv50<0.05)/reps)
w.pvals.B
