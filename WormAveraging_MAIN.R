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
wormbootOnCounts<-function(n_reps, mydata, batch_sizes=c(1,5,10,20,50), foldD=10, correction_constant=1){
  # Expects a number of replicates for the bootstrap (n_reps)
  # and a data frame of worm CFU data for individuals (mydata)
  # where the number of colonies counted is in column "Count"
  # the dilution at which these colonies were measured is in column "D",
  # FoldD(num) is fold dilution in dilution series 
  #  (Default is 10X dilutions, e.g. each step in serial dilution is 1:10 volume)
  #
  ## The dilution correction factor (numeric) is given as "correction_constant"
  # and is the ratio of the original volume and the volume plated
  # e.g. for 10 uL spots and an original volume of 1 mL, correction_constant = 100
  # 
  # Returns a data frame of simulated batch digests
  # with batch sizes in vector batch_sizes, default (1, 5, 10, 20, 50) worms/batch
  # values reported as inferred CFU/worm and log10(CFU/worm)
  
  # get size of data
  capp<-dim(mydata)[1]
  
  # make someplace to put stuff
  temp<-vector("list", length=length(batch_sizes)*n_reps)
  
  #batch5<-rep(0,n_reps)
  #batch10<-rep(0,n_reps)
  #batch20<-rep(0,n_reps)
  #batch50<-rep(0, n_reps)
  for(i in 1:n_reps){
    #idx5<-sample(1:capp,5,replace=TRUE)
    #idx10<-sample(1:capp,10,replace=TRUE)
    #idx20<-sample(1:capp,20,replace=TRUE)
    #idx50<-sample(1:capp,50,replace=TRUE)
    
    temp<-rep(0, length(batch_sizes))
    for (j in seq_len(batch_sizes)){
      # Randomly pull samples from data for batching
      idx<-sample(1:capp,batch_sizes[j],replace=TRUE)
      # Assume Poisson count error and generate new counts
      temp_count<-rpois(batch_sizes[j], mydata$Count[idx])
      # Calculate CFU/worm
      temp_final<-correction_constant*temp_count*foldD^mydata$D[idx]
      temp[j]<-mean(temp_final)
    }
    #temp_count<-rpois(5, mydata$Count[idx5])
    #temp<-correction_constant*temp_count*10^mydata$D[idx5]
    #batch5[i]<-mean(temp)
    
    #temp_count<-rpois(10, mydata$Count[idx10])
    #temp<-correction_constant*temp_count*10^mydata$D[idx10]
    #batch10[i]<-mean(temp)
    
    #temp_count<-rpois(20, mydata$Count[idx20])
    #temp<-correction_constant*temp_count*10^mydata$D[idx20]
    #batch20[i]<-mean(temp)
    
    #temp_count<-rpois(50, mydata$Count[idx50])
    #temp<-correction_constant*temp_count*10^mydata$D[idx50]
    #batch50[i]<-mean(temp)
    
    temp[[i]]<-tibble(Batch=batch_sizes,
                      FinalCount=temp) 
  }
  
  # unfold data
  dataSet<-dplyr::bind_rows(temp)
  
  dataSet<-dataSet %>%
    mutate(logCFU=log10(FinalCount))
  
  #batch5log<-log10(batch5+1)
  #batch10log<-log10(batch10+1)
  #batch20log<-log10(batch20+1) 
  #batch50log<-log10(batch50+1) 
  #Batch<-c(rep(1,times=capp), rep(5, times=n_reps), rep(10, times=n_reps), rep(20, times=n_reps), rep(50, times=n_reps))
  #mylogCFU<-log10(mydata$CFU+1)
  #logCFU<-c(mylogCFU, batch5log, batch10log, batch20log, batch50log)
  #CFU<-c(mydata$CFU, batch5, batch10, batch20, batch50)
  #dataSet<-data.frame(Batch, CFU, logCFU)
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

BatchDigests<-read.table("batchdigests.txt", header=TRUE)
BatchDigests$Batch<-as.factor(BatchDigests$Batch)
pBatchSA<-BatchDigests %>%
  ggplot(aes(x=Batch, y=logCFU, color=Batch)) + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + theme_classic() + 
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
#ggsave("pBatchSA.png", width=4, height=3, units="in", dpi=300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Now Salmonella enterica LT2-GFP single worms (fed on plates)
# and data for S. aureus-GFP

SaSeCount<-read_xlsx("SaSeCount.xlsx")  

SaSeCount %>%
  ggplot(aes(x=factor(Rep), y=logCFU, color=factor(Rep))) + 
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
  group_by(Condition, Rep) %>%
  summarise(count=n(),
            medianCFU=median(CFU, na.rm=TRUE),
            meanCFU=mean(CFU, na.rm=TRUE),
            sdCFU=sd(CFU, na.rm=TRUE),
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
	ggplot(aes(x=factor(Rep), y=logCFU, color=factor(Rep))) + 
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
  group_by(Condition, Rep) %>%
  summarise(count=n(),
            medianCFU=median(CFU, na.rm=TRUE),
            meanCFU=mean(CFU, na.rm=TRUE),
            sdCFU=sd(CFU, na.rm=TRUE),
            medianLogCFU=median(logCFU, na.rm=TRUE),
            meanLogCFU=mean(logCFU, na.rm=TRUE),
            sdLogCFU=sd(logCFU, na.rm = TRUE)
  )



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#        GENERATE BOOTSTRAPPED DATA FROM S. ENTERICA COLONIZATION
#
# here we will call the wormboot function (wormbootOnCounts)
# to generate simulated data and make comparisons

temp1<-SaSeCount %>%
  dplyr::filter(Condition=="SE" & Rep=="1")
temp2<-SaSeCount %>%
  dplyr::filter(Condition=="SE" & Rep=="2")

# summary statistics
mean(temp1$CFU)
median(temp1$CFU)
mean(temp2$CFU)
median(temp2$CFU)

SeBoot1<-wormbootOnCounts(dim(temp1)[1], temp1, 20)
SeBoot2<-wormbootOnCounts(dim(temp2)[1], temp2, 20)
SeBoot1$Rep<-as.factor("Rep1")
SeBoot2$Rep<-as.factor("Rep2")
jointSeBoot<-rbind(SeBoot1, SeBoot2)

pJointSEBoot<-jointSeBoot %>%
	ggplot(aes(x=factor(Rep), y=logCFU, color=factor(Rep))) + 
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
       y=expression(log[10](CFU/Worm)), x="Replicate")+
  stat_compare_means(method="t.test", label.y = 6.2, size=3.2)+
  stat_compare_means(method="wilcox.test", label.y=5.9, size=3.2)
pJointSEBoot

#and some statistical tests on the bootstrapped data
t.test(SeBoot1$logCFU[SeBoot1$Batch==1], SeBoot2$logCFU[SeBoot2$Batch==1])
t.test(SeBoot1$logCFU[SeBoot1$Batch==5], SeBoot2$logCFU[SeBoot2$Batch==5])
t.test(SeBoot1$logCFU[SeBoot1$Batch==10], SeBoot2$logCFU[SeBoot2$Batch==10])
t.test(SeBoot1$logCFU[SeBoot1$Batch==20], SeBoot2$logCFU[SeBoot2$Batch==20])
t.test(SeBoot1$logCFU[SeBoot1$Batch==50], SeBoot2$logCFU[SeBoot2$Batch==50])

# and tests for normality
# Shapiro-Wilk
shapiro.test(SeBoot1$logCFU[SeBoot1$Batch==1])
shapiro.test(SeBoot1$logCFU[SeBoot1$Batch==5])
shapiro.test(SeBoot1$logCFU[SeBoot1$Batch==10])
shapiro.test(SeBoot1$logCFU[SeBoot1$Batch==20])
shapiro.test(SeBoot1$logCFU[SeBoot1$Batch==50])
shapiro.test(SeBoot2$logCFU[SeBoot2$Batch==1])
shapiro.test(SeBoot2$logCFU[SeBoot2$Batch==5])
shapiro.test(SeBoot2$logCFU[SeBoot2$Batch==10])
shapiro.test(SeBoot2$logCFU[SeBoot2$Batch==20])
shapiro.test(SeBoot2$logCFU[SeBoot2$Batch==50])
rm(SeBoot1, SeBoot2, temp1, temp2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# same SE bootstrap with zeros set to TOD
temp<-dplyr::filter(SaSeCount, Condition=="SE")
idx<-which(temp$Count==0)
temp$Count[idx]<-1
temp$CFU[idx]<-20
temp$logCFU[idx]<-log10[20]

temp1A<-temp %>%
  dplyr::filter(Rep=="1")
temp2A<-temp %>%
  dplyr::filter(Rep=="2")

SeBoot1A<-wormbootOnCounts(dim(temp1A)[1], temp1A, 20)
SeBoot2A<-wormbootOnCounts(dim(temp2A)[1], temp2A, 20)
SeBoot1A$Rep<-as.factor("Rep1")
SeBoot2A$Rep<-as.factor("Rep2")
jointSeBoot_nozeros<-rbind(SeBoot1A, SeBoot2A)

pJointSEBoot_nozeros<-jointSeBoot_nozeros %>%
  dplyr::filter(Batch<20) %>% # sufficient to make the point
  ggplot(aes(x=factor(Rep), y=logCFU, color=factor(Rep))) + 
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
  labs(title="Zeros removed", x="Replicate", y=expression(log[10](CFU/Worm)))+
  stat_compare_means(method="t.test", label.y = 6.2, size=3.2)+
  stat_compare_means(method="wilcox.test", label.y=5.9, size=3.2)
pJointSEBoot_nozeros

rm(SeBoot1A, SeBoot2A, temp1A, temp2A)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# and with extra zeros

# implemented as higher TOD, call it 10^3
temp<-dplyr::filter(SaSeCount, Condition=="SE")
idx<-which(temp$logCFU<3)
temp$Count[idx]<-0
temp$D[idx]<-0
temp$CFU[idx]<-0
temp$logCFU[idx]<-0

temp1B<-temp %>%
  dplyr::filter(Rep=="1")
temp2B<-temp %>%
  dplyr::filter(Rep=="2")

SeBoot1B<-wormbootOnCounts(dim(temp1B)[1], temp1B, 20)
SeBoot2B<-wormbootOnCounts(dim(temp2B)[1], temp2B, 20)
SeBoot1B$Rep<-as.factor("Rep1")
SeBoot2B$Rep<-as.factor("Rep2")
jointSeBoot_zeros3<-rbind(SeBoot1B, SeBoot2B)

pJointSEBoot_zeros3<-jointSeBoot_zeros3 %>%
  dplyr::filter(Batch<20) %>% # sufficient to make the point
  ggplot(aes(x=factor(Rep), y=logCFU, color=factor(Rep))) + 
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
  labs(title="Zero-enriched", x="Replicate", y=expression(log[10](CFU/Worm)))+
  stat_compare_means(method="t.test", label.y = 6.2, size=3.2)+
  stat_compare_means(method="wilcox.test", label.y=5.9, size=3.2)
pJointSEBoot_zeros3

rm(SeBoot1B, SeBoot2B, temp1B, temp2B, temp)


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
SaSeCount2_subset<- subset(SaSeCount2, select=c("Condition", "Rep", "CFU", "logCFU")) 

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

# here we will call the wormboot function (wormbootOnCounts)
# to generate simulated data and make comparisons
# 22 data points to match the smallest data set (N2)
# batch 50 is pretty large for these data so we'll stop at 20

tempBoot1A<-wormbootOnCounts(22, temp1A, 2)
tempBoot2A<-wormbootOnCounts(22, temp2A, 2)
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
tempBoot1B<-wormbootOnCounts(24, temp1, 2)
tempBoot2B<-wormbootOnCounts(24, temp2, 20)
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
  tempBoot1<-WormbootOnCounts(22, temp1A, 2)
  tempBoot2<-wormbootOnCounts(22, temp2A, 2)
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
  tempBoot1<-wormbootOnCounts(24, temp1, 2)
  tempBoot2<-wormbootOnCounts(24, temp2, 20)
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
