library(ggplot2)
library(tidyverse)
library(cowplot)
library(mclust)
library(e1071)
library(ggpubr)

# Code variously stolen from
# https://www.datanovia.com/en/blog/how-to-change-ggplot-facet-labels/


save.image("WormDigestHetSims.Rdata")

load("WormDigestHetSims.Rdata")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~   Figure 1A
# Real data from batch digests 2022-2-18, SA and SE
# SA looks fine, but didn't get enough worms for SE, and needed to plate higher dilutions
BatchDigests<-read.table("batchdigests.txt", header=TRUE)
BatchDigests$Batch<-as.factor(BatchDigests$Batch)
pBatchSA<-subset(BatchDigests, Species=="SA") %>%
  ggplot(aes(x=Batch, y=logCFU, color=Batch)) + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + theme_classic() + 
  theme(text=element_text(size=16), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=14),
        legend.position=c(0.9,0.3)) + 
  labs(title=expression(paste(italic("S. aureus"), " Newman")), y=expression(log[10](CFU/Worm)))
pBatchSA
ggsave("pBatchSA.png", width=4, height=3, units="in", dpi=300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Now Salmonella enterica LT2-GFP single worms (fed on plates), same as in bleach protocol
# and corresponding data for S. aureus-GFP
# including single worm data from the batch digest run on 2022-2-18

SaSeCount<-read.table("SaSeCount.txt", header=TRUE)  
#SaSeCount$Condition<-as.factor(SaSeCount$Condition)
#SaSeCount$Date<-as.factor(SaSeCount$Date)
#SaSeCount$Rep<-as.factor(SaSeCount$Rep)

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
  facet_wrap(vars(Condition), scales="free_x") +
  labs(y="log10(CFU/worm)", x="Replicate")
ggsave("pSaSeLogCFUByReplicate.png", width=8, height=5, units="in", dpi=400)

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
SaSeCount2[SaSeCount2==0]<-NA
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
  labs(y="log10(CFU/worm)", x="Replicate")
	
#labs(title=expression(paste(italic("S. enterica"), " LT2")), y="log10(CFU/worm)")
#ggsave("pSeByDay.png",width=4, height=3, units="in", dpi=300)

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


#another version of this plot with the combined data (all 3 days as rep 0)
mylen<-dim(SaSeCount2)[1]
SaSeTemp<-rbind(SaSeCount2, data.frame(Treatment=SaSeCount2$Treatment, Condition=SaSeCount2$Condition,
                                   Date=rep("All", mylen),
                                  Rep=rep("0", mylen),
                                   Count=SaSeCount2$Count,
                                   D=SaSeCount2$D,
                                  CFU=SaSeCount2$CFU,
                                   logCFU=SaSeCount2$logCFU
))
SaSeTemp %>%
  ggplot(aes(x=Rep, y=logCFU, color=Rep)) +
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + 
  theme_classic() + 
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=14)) + 
  facet_wrap(vars(Condition), scales="free_x")
  #  labs(title=expression(paste(italic("S. enterica"), " LT2")), y="log10(CFU/worm)"
#ggsave("pSECountAll.png", width=4, height=3, units="in", dpi=300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#        GENERATE BOOTSTRAPPED DATA FROM S. ENTERICA COLONIZATION
#
# here we will call the wormboot function (wormboot.R)
# to generate simulated data and make comparisons
# 25 reps as usual

temp1<-SaSeCount %>%
  filter(Condition=="SE" & Rep=="1")
temp2<-SaSeCount %>%
  filter(Condition=="SE" & Rep=="2")

SeBoot1<-wormbootOnCounts(25, temp1, 20)
SeBoot2<-wormbootOnCounts(25, temp2, 20)
SeBoot1$Rep<-as.factor("Rep1")
SeBoot2$Rep<-as.factor("Rep2")
jointSeBoot<-rbind(SeBoot1, SeBoot2)
#SeBleachCount$Date<-as.factor(SeBleachCount$Date)

pJointSEBoot<-jointSeBoot %>%
	ggplot(aes(x=factor(Rep), y=logCFU, color=factor(Rep))) + 
	geom_jitter(shape=16, position=position_jitter(0.05)) +
	geom_violin(fill=NA) + 
	ylim(-0.1,6)+ theme_classic() + 
	theme(
	text=element_text(size=14), 
	axis.title.x = element_blank(),
	axis.text.x = element_blank(),
	plot.title=element_text(hjust=0.5,size=14),
	legend.position=c(0.9,0.3),
	legend.title = element_blank())+
	#legend.position="none") +
  facet_wrap(vars(Batch), ncol=5)+
	labs(x="Replicate", y=expression(log[10](CFU/Worm)))
#ggsave("pjointSeBoot.svg", width=9, height=4, units="in", dpi=300)
pJointSEBoot
#ggsave("pSeBootTwoRepsByColonies_unpruned_ByReplicate.png", width=8, height=5, units="in", dpi=400)

#and some statistical tests on the bootstrapped data
t.test(SeBoot1$CFU[SeBoot1$Batch==1], SeBoot2$CFU[SeBoot2$Batch==1])
t.test(SeBoot1$CFU[SeBoot1$Batch==5], SeBoot2$CFU[SeBoot2$Batch==5])
t.test(SeBoot1$CFU[SeBoot1$Batch==10], SeBoot2$CFU[SeBoot2$Batch==10])
t.test(SeBoot1$CFU[SeBoot1$Batch==20], SeBoot2$CFU[SeBoot2$Batch==20])
t.test(SeBoot1$CFU[SeBoot1$Batch==50], SeBoot2$CFU[SeBoot2$Batch==50])

wilcox.test(SeBoot1$CFU[SeBoot1$Batch==1], SeBoot2$CFU[SeBoot2$Batch==1])
wilcox.test(SeBoot1$CFU[SeBoot1$Batch==5], SeBoot2$CFU[SeBoot2$Batch==5])
wilcox.test(SeBoot1$CFU[SeBoot1$Batch==10], SeBoot2$CFU[SeBoot2$Batch==10])
wilcox.test(SeBoot1$CFU[SeBoot1$Batch==20], SeBoot2$CFU[SeBoot2$Batch==20])
wilcox.test(SeBoot1$CFU[SeBoot1$Batch==50], SeBoot2$CFU[SeBoot2$Batch==50])

t.test(SeBoot1$logCFU[SeBoot1$Batch==1], SeBoot2$logCFU[SeBoot2$Batch==1])
t.test(SeBoot1$logCFU[SeBoot1$Batch==5], SeBoot2$logCFU[SeBoot2$Batch==5])
t.test(SeBoot1$logCFU[SeBoot1$Batch==10], SeBoot2$logCFU[SeBoot2$Batch==10])
t.test(SeBoot1$logCFU[SeBoot1$Batch==20], SeBoot2$logCFU[SeBoot2$Batch==20])
t.test(SeBoot1$logCFU[SeBoot1$Batch==50], SeBoot2$logCFU[SeBoot2$Batch==50])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~    Build Figure 1
# Plot out SE bootstraps together with real data from S. aureus batch size experiment
plot_grid(pBatchSA, pJointSEBoot, labels="AUTO", ncol=2, rel_widths=c(1,3))

#ggsave("pjointSeBoot.svg", width=9, height=4, units="in", dpi=300)
ggsave("Fig1_pSaBatchJointSeBoot.png", width=10, height=4, units="in", dpi=400)

# and some summary statistics for the record
mean(BatchDigests$CFU[BatchDigests$Species=="SA" & BatchDigests$Batch==1 ])
median(BatchDigests$CFU[BatchDigests$Species=="SA" & BatchDigests$Batch==1 ])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~                 FIGURE 2
#~~~~    WITHIN AND BETWEEN REPLICATE HETEROGENEITY OF CFUS
#
# First, pull in the all-8 data set so we can describe run to run variation
# just need the total counts for this

All8<-read.table("All8.txt", header=TRUE)
names(All8)
hosts<-unique(All8$Host)
All8$logTotal<-log10(All8$Total)
All8[All8=="decc1"]<-"dec1"
All8$logTotal[is.infinite(All8$logTotal)]<-NA
All8<-All8[complete.cases(All8),]

# filter the original data set to remove anything with <18 measurements per strain on a given day
temp<-unique(All8$Host)
All8.temp<-All8
for (i in 1:length(temp)){
	#print(i)
	tempdata<-subset(All8.temp, Host==temp[i])
	temp2<-unique(tempdata$Date)
	for (j in 1:length(temp2)){
		#print(temp[i])
		mycount<-(length(All8.temp$logTotal[All8.temp$Host==temp[i] & All8.temp$Date==temp2[j]]))
		#print(mycount<18)
		if(mycount<18) {
			All8.temp$logTotal[All8.temp$Host==temp[i] & All8.temp$Date==temp2[j]]<-NA
		}
	}
}
All8.temp$logTotal[is.infinite(All8.temp$logTotal)]<-NA
All8.big<-All8.temp[complete.cases(All8.temp),]
rm(All8.temp, temp)

### let's see what's left
hosts.big<-unique(All8.big$Host)
temp<-integer()
for (i in 1:length(hosts.big)){
	temp[i]<-length(unique(All8.big$Date[All8.big$Host==hosts.big[i]]))}
temp
#[1] 2 1 3 3 2 2 2 1 1 1 4 1 2
temp2<-unique(All8.big$Host)
#[1] "AU37"   "ctls40" "daf16"  "daf2"   "dbl1"   "dec1"   "exp1"   "glp1"   "glp4"   "hsf1"   "N2"    
#[12] "phm2"   "vhp1"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~  Building the filtered All8 data set
# Let's filter out anything that doesn't have 2 or more data sets left
# these have enough runs to be worth it
temp2[which(temp>1)]
# [1] "AU37"  "daf16" "daf2"  "dbl1"  "dec1"  "exp1"  "N2"    "vhp1"
All8.big<-subset(All8.big, Host=="AU37" | Host=="daf16"| Host=="daf2" | Host=="dbl1" | Host=="dec1" | Host=="exp1" | Host=="N2" | Host=="vhp1")

# And get dates into a consistent format
#unique(All8.big$Date)
# 21319    22519  7122019  7292019  8162019  5102019  9112019  9162019  6262019  6142019  6192019
#[12]  6212019 10252019 10202021

All8.big$Date<-as.character(All8.big$Date)
All8.big[All8.big=="21319"]<-"2/13/19"
All8.big[All8.big=="22519"]<-"2/25/19"
All8.big[All8.big=="5102019"]<-"5/10/19"
All8.big[All8.big=="9112019"]<-"9/11/19"
All8.big[All8.big=="9162019"]<-"9/16/19"
All8.big[All8.big=="6142019"]<-"6/14/19"
All8.big[All8.big=="6192019"]<-"6/19/19"
All8.big[All8.big=="6212019"]<-"6/21/19"
All8.big[All8.big=="6262019"]<-"6/26/19"
All8.big[All8.big=="7122019"]<-"7/12/19"
All8.big[All8.big=="7292019"]<-"7/29/19"
All8.big[All8.big=="8162019"]<-"8/16/19"
All8.big[All8.big=="10252019"]<-"10/25/19"
All8.big[All8.big=="10202021"]<-"10/20/21"
All8.big$Date<-as.factor(All8.big$Date)

# Plot out individual replicates by host genotype
# first create new facet labels
host.labs<-c("AU37", "daf-16", "daf-2", "dbl-1", "dec-1", "exp-1", "N2", "vhp-1")
names(host.labs)<-c("AU37", "daf16", "daf2", "dbl1", "dec1", "exp1", "N2", "vhp1")
pAll8.big.byHost<-All8.big%>%
  ggplot(aes(x=Date, y=logTotal, color=Date)) + geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + 	ylim(1,6) +	theme_classic() + 
  labs(title="Minimal community colonization, by replicate", y=expression(log[10](CFU/Worm))) +
  facet_wrap(~Host, scales="free_x", ncol=4, labeller=labeller(Host=host.labs))+
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        legend.position="none", 
        plot.title=element_text(hjust=0.5, size=16),
        strip.text.x = element_text(size=12, face="italic"))
pAll8.big.byHost
ggsave("pAll8BigLogCFUHostByDay.png", width=12, height=9, units="in", dpi=300)

# Plot out the log-transformed CFU data by host genotype, replicates combined
pAll8.big<-All8.big %>%
	ggplot(aes(x=as.factor(Host), y=logTotal, color=as.factor(Host)))+
	geom_jitter(shape=16, position=position_jitter(0.05)) +
	geom_violin(fill=NA) + 	ylim(1,6) +	theme_classic() + 
	theme(text=element_text(size=14), axis.title.x = element_blank(), legend.position="none", plot.title=element_text(hjust=0.5, size=16)) +
	labs(title="Minimal Microbiome Totals", y=expression(log[10](CFU/Worm)))
pAll8.big

# 2022-06-03: bringing in single worm data from NRRL1 evolution paper
NRRL1counts<-read_csv("NRRL1counts.csv")
names(NRRL1counts)
#[1] "Condition" "CFU"       "logCFU"    "Date"
NRRL1counts$Date<-as.character(NRRL1counts$Date)

#as before loop over experiment (Condition) and rep (Date) and record
hosts2<-unique(NRRL1counts$Condition) #length 24

# Only mono-colonization by MO.I2.9 (hosts[10]) has 2 reps of >18 (both n=24)
# subset out & prep for rbind
hosts2[10]
MOcount<-subset(NRRL1counts, Condition==hosts2[10], select=c("Condition", "Date", "CFU", "logCFU")) %>%
  rename("Host"="Condition",
         "Total"="CFU",
         "logTotal"="logCFU") %>%
  as_tibble()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# OK now we are going to take a WHOLE BUNCH of summary stats and see what we get

AllCounts<-as_tibble(subset(All8.big, select=c("Host", "Date", "Total", "logTotal")))
glimpse(AllCounts)

# Incorporate the pathogen colonization data: prep for rbind
#SeCount<-SeCount[complete.cases(SeCount),]
#SeCount$Host<-"SE"
#SaBleachCount2$Host<-"SA"
SaSeCount2_subset<- subset(SaSeCount2, select=c("Condition", "Date", "CFU", "logCFU")) %>%
  rename("Host" = "Condition", 
        "Total"="CFU",
         "logTotal"="logCFU") %>%
  as_tibble()
SaSeCount2_subset$Date<-as.character(SaSeCount2_subset$Date)
glimpse(SaSeCount2_subset)
glimpse(MOcount)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Merge all the single-worm CFU data
AllCounts2<-rbind(AllCounts, SaSeCount2_subset, MOcount) %>%
  rename("Condition"="Host")
glimpse(AllCounts2)

# Get a list of all experiments
conditions.big<-unique(AllCounts2$Condition)

# set up the data frame to hold our calculations
AllCountsStats.r2r <- data.frame(Condition=character(),
                 Date=character(),
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
	temp<-unique(AllCounts2$Date[AllCounts2$Condition==conditions.big[i]])
#	print(conditions.big[i]) # for debugging
	for (j in seq_along(temp)){
#		print(temp[j]) # for debugging
		tempCFU<-sort(AllCounts2$Total[AllCounts2$Condition==conditions.big[i] & AllCounts2$Date==temp[j]])
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
		Date=temp[j],
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

pAllCounts<-AllCounts2 %>%
  mutate(Condition=factor(Condition, levels=c("N2","AU37", "daf16", "daf2", "dbl1", "dec1", "exp1", "vhp1", "MO.I2.9", "SA", "SE"))) %>%
  ggplot(aes(x=Date, y=logTotal, color=Date)) + 
  geom_point(shape=16, position=position_jitterdodge(0.05)) +
  geom_violin(fill=NA) + 
  theme_classic() + 
  facet_wrap(~Condition, scales="free_x", labeller=labeller(Condition=condition.labs))+
  theme(text=element_text(size=16), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        #plot.title=element_text(hjust=0.5, size=14),
        legend.position = "none",
        strip.text.x = element_text(size=12, face="italic"))+
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
  geom_boxplot(fill=NA) +	theme_classic() + 
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(), 
        legend.position="none", 
        plot.title=element_text(hjust=0.5, size=16)) +
  labs(title="", y="Distance between pairs")+
  facet_wrap(vars(Var), scales="free_y", ncol=2, labeller=labeller(Var=my.labs))+
  stat_compare_means(label.x.npc = 0.5, label.y.npc=0.9, size=4)
pAllCountStats.r2r
ggsave("pAll8SelectMomentsSummaryDistance.png", width=9, height=12, units="in", dpi=400)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~   FIGURE 2 ASSEMBLY
plot_grid(pAllCounts, pAllCountStats.r2r, rel_widths=c(1, 1.2), labels="AUTO")
ggsave("Fig2_AllCounts_logCFU_MomentDist.png", width=15, height=9, units="in")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# And some comparisons for the record
wilcox.test(meanCFUdist~same, data=All8.big.r2r.dist)
wilcox.test(cvdist~same, data=All8.big.r2r.dist)
wilcox.test(skewdist~same, data=All8.big.r2r.dist)
wilcox.test(kurtdist~same, data=All8.big.r2r.dist)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#########       FIGURE 3          ##################################
#####     TRYING FOR FALSE NEGATIVES
#####        with real data 
# Pairs 
#(N2 6/21/19 vs dbl-1 6/26/19?) in jointTempBootA 
#(SE day 2 vs. N2 adults 2/13/19) in jointTempBoot
# I could clean up this plotting but I won't

unique(All8.big$Date)
temp2A<-All8.big$Total[All8.big$Host=="N2" & All8.big$Date=="6/21/19"]
temp1A<-All8.big$Total[All8.big$Host=="dbl1" & All8.big$Date=="6/26/19"]

# here we will call the wormboot function (wormboot.R)
# to generate simulated data and make comparisons
# 21 reps to match the smallest data set (N2); other is 33
# batch 50 is pretty large for these data so we'll stop at 20
#note that wormboot() has changed to allow isRand=TRUE (default is FALSE)
#which creates U[up,down] errors in CFU counts

tempBoot1A<-wormboot(21, temp1A)
tempBoot2A<-wormboot(21, temp2A)
tempBoot1A$Host<-as.factor("dbl1")
tempBoot2A$Host<-as.factor("N2")
jointTempBootA<-rbind(tempBoot1A, tempBoot2A)
names(jointTempBootA)

temp1A<-subset(All8.big, Host=="N2" & Date=="6/21/19")
temp2A<-subset(All8.big, Host=="dbl1" & Date=="6/26/19")
jointTempA<-rbind(temp1A, temp2A)

#names(All8.big)
#names(SeCount2)

temp1<-subset(All8.big, Host=="N2" & Date=="2/13/19", select = c("Host", "logTotal"))
temp2<-subset(SeCount2, Date=="2", select=c("Species", "logCount"))
#mean(SeCount2$Count[SeCount2$Date=="2"])
#length(SeCount2$Count[SeCount2$Date=="2"])
#skewness(SeCount2$Count[SeCount2$Date=="2"])

temp2<-temp2 %>%
  rename(Host=Species, logTotal=logCount)
jointTemp<-rbind(temp1, temp2)
tempBoot1$Host<-as.factor("All8")
tempBoot2$Host<-as.factor("SE")
jointTempBoot<-rbind(tempBoot1, tempBoot2)


#########  First set (N2 6/21/19 vs dbl-1 6/26/19 multispp) ######

pJointTemp1A<- jointTempA  %>%
  ggplot(aes(x=Host, y=logTotal, color=Host)) + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + 
  ylim(-0.1,6)+ theme_classic() + 
  theme(
    text=element_text(size=14), 
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    plot.title=element_text(hjust=0.5,size=14),
    legend.position="none") +
  labs(title="Raw data", y="log10(CFU/worm)")+
  stat_compare_means(label.y = 0.4, label.x=1.2)
pJointTemp1A
pTempBoot2A<-subset(jointTempBootA, batch==5) %>%
  ggplot(aes(x=Host, y=logCount, color=Host)) + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + ylim(-0.1,6)+ theme_classic() + 
  theme(text=element_text(size=16), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=14),
        legend.position="none") + 
  labs(title="Batch 5", y="")+
  stat_compare_means(label.y = 0.4, label.x=1.2)
pTempBoot2A
pTempBoot3A<-subset(jointTempBootA, batch==10) %>%
  ggplot(aes(x=Host, y=logCount, color=Host)) + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + ylim(-0.1,6)+ theme_classic() + 
  theme(text=element_text(size=16), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=14),
        legend.position="none") + 
  labs(title="Batch 10", y="")+
  stat_compare_means(label.y = 0.4, label.x=1.2)
pTempBoot3A
pTempBoot4A<-subset(jointTempBootA, batch==20) %>%
  ggplot(aes(x=Host, y=logCount, color=Host)) + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + ylim(-0.1,6)+ theme_classic() + 
  theme(text=element_text(size=16), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=14),
        legend.position=c(0.7,0.3),
        legend.title = element_blank()) + 
  labs(title="Batch 20", y="")+
  stat_compare_means(label.y = 0.4, label.x=1.2)


######   Second set (SE vs N2 Multi)  ########

pJointTemp1<- jointTemp  %>%
  ggplot(aes(x=Host, y=logTotal, color=Host)) + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + 
  ylim(-0.1,6)+ theme_classic() + 
  theme(
    text=element_text(size=14), 
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    plot.title=element_text(hjust=0.5,size=14),
    legend.position="none") +
  labs(title="Raw data", y="log10(CFU/worm)")+
  stat_compare_means(label.y = 0.4, label.x=1.2)
pJointTemp1
pTempBoot2<-subset(jointTempBoot, batch==5) %>%
  ggplot(aes(x=Host, y=data, color=Host)) + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + ylim(-0.1,6)+ theme_classic() + 
  theme(text=element_text(size=16), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=14),
        legend.position="none") + 
  labs(title="Batch 5", y="")+
  stat_compare_means(label.y = 0.4, label.x=1.2)
pTempBoot2
pTempBoot3<-subset(jointTempBoot, batch==10) %>%
  ggplot(aes(x=Host, y=data, color=Host)) + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + ylim(-0.1,6)+ theme_classic() + 
  theme(text=element_text(size=16), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=14),
        legend.position="none") + 
  labs(title="Batch 10", y="")+
  stat_compare_means(label.y = 0.4, label.x=1.2)
pTempBoot3
pTempBoot4<-subset(jointTempBoot, batch==20) %>%
  ggplot(aes(x=Host, y=data, color=Host)) + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + ylim(-0.1,6)+ theme_classic() + 
  theme(text=element_text(size=16), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=14),
        legend.position=c(0.7,0.3),
        legend.title = element_blank()) + 
  labs(title="Batch 20", y="")+
  stat_compare_means(label.y = 0.4, label.x=1.2)
pTempBoot4
plot_grid(pJointTemp1A, pTempBoot2A, pTempBoot3A, pTempBoot4A,
          pJointTemp1, pTempBoot2, pTempBoot3, pTempBoot4, labels="AUTO", ncol=4, nrow=2, align="h")
#ggsave("pjointSeBoot.svg", width=9, height=4, units="in", dpi=300)
ggsave("Fig3_pjointFalseNegativeBoot.png", width=10, height=7, units="in", dpi=300)

tempBoot1$Count<-as.double(10^tempBoot1$data)
tempBoot2$Count<-as.double(10^tempBoot2$data)
wilcox.test(temp1$Total, temp2$Total)
wilcox.test(tempBoot1$Count[tempBoot1$batch==1], tempBoot2$Count[tempBoot2$batch==1])
wilcox.test(tempBoot1$Count[tempBoot1$batch==5], tempBoot2$Count[tempBoot2$batch==5])
wilcox.test(tempBoot1$Count[tempBoot1$batch==10], tempBoot2$Count[tempBoot2$batch==10])
wilcox.test(tempBoot1$Count[tempBoot1$batch==20], tempBoot2$Count[tempBoot2$batch==20])

# let's do this a bunch of times
#note that wormboot() has changed to allow isRand=TRUE (default is FALSE)
#which creates U[up,down] errors in CFU counts
reps<-10000
wpv1<-numeric(reps)
wpv5<-numeric(reps)
wpv10<-numeric(reps)
wpv20<-numeric(reps)
#temp2<-All8.big$Total[All8.big$Host=="N2" & All8.big$Date=="6/21/19"]
#temp1<-All8.big$Total[All8.big$Host=="dbl1" & All8.big$Date=="6/26/19"]
#names(SeCount2)
#unique(SeCount2$Date)
temp2<-SeCount2$Count[SeCount2$Date=="2"]
#mean(temp2)
temp1<-All8.big$Total[All8.big$Host=="N2" & All8.big$Date=="2/13/19"]
#mean(temp1)
#wilcox.test(temp1, temp2)

for (j in 1:reps){
  tempBoot1<-wormboot(24, temp1, isRand=TRUE)
  tempBoot2<-wormboot(24, temp2, isRand=TRUE)
  wtest1<-wilcox.test(tempBoot1$Count[tempBoot1$batch==1], tempBoot2$Count[tempBoot2$batch==1])
  wtest5<-wilcox.test(tempBoot1$Count[tempBoot1$batch==5], tempBoot2$Count[tempBoot2$batch==5])
  wtest10<-wilcox.test(tempBoot1$Count[tempBoot1$batch==10], tempBoot2$Count[tempBoot2$batch==10])
  wtest20<-wilcox.test(tempBoot1$Count[tempBoot1$batch==20], tempBoot2$Count[tempBoot2$batch==20])
  wpv1[j]<-wtest1$p.value
  wpv5[j]<-wtest5$p.value
  wpv10[j]<-wtest10$p.value
  wpv20[j]<-wtest20$p.value
}
w.pvals<-c(sum(wpv1<0.05)/reps, sum(wpv5<0.05)/reps, sum(wpv10<0.05)/reps, sum(wpv20<0.05)/reps)
w.pvals
#wtest5
#wtest10
#wtest20
