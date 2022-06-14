library(ggplot2)
library(tidyverse)
library(cowplot)
library(mclust)
library(e1071)
library(ggpubr)


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
  labs(title=expression(paste(italic("S. aureus"), " Newman")), y="log10(CFU/worm)")
pBatchSA
ggsave("pBatchSA.png", width=4, height=3, units="in", dpi=300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Now Salmonella enterica LT2-GFP single worms (fed on plates), same as in bleach protocol
# and corresponding data for S. aureus-GFP
SeBleachCount<-read.table("SeBleach.txt", header=TRUE)  
SeBleachCount$logCount<-log10(SeBleachCount$Count+1) #zero correct  
SeBleachCount[SeBleachCount=="B"]<-"Bleach"
SeBleachCount[SeBleachCount=="NB"]<-"NoBleach"
SeBleachCount$Condition<-as.factor(SeBleachCount$Condition)

#SeBleachCount[SeBleachCount==20210315]<-1 #replace experiment date with day ID
#SeBleachCount[SeBleachCount==20210308]<-2
SeBleachCount$Date<-as.factor(SeBleachCount$Date)

SaBleachCount<-read.table("SaBleach.txt", header=TRUE)  
SaBleachCount$logCount<-log10(SaBleachCount$Count+1) #zero correct  
SaBleachCount[SaBleachCount=="B"]<-"Bleach"
SaBleachCount[SaBleachCount=="NB"]<-"NoBleach"
SaBleachCount$Condition<-as.factor(SaBleachCount$Condition)

pSaBleachCount<-ggplot(SaBleachCount, aes(x=Condition, y=logCount, color=Condition)) + geom_jitter(shape=16, position=position_jitter(0.05)) 
p1<-pSaBleachCount + geom_violin(fill=NA) + ylim(-0.1,6)+ theme_classic() + theme(text=element_text(size=14), axis.title.x = element_blank(), plot.title=element_text(hjust=0.5)) +labs(title=expression(paste(italic("S. aureus"), " MSSA Newman")), y="log10(CFU/worm)")
pSeBleachCount<-ggplot(SeBleachCount, aes(x=Condition, y=logCount, color=Condition)) + geom_jitter(shape=16, position=position_jitter(0.05)) 
p2<-pSeBleachCount + geom_violin(fill=NA) + ylim(-0.1,6)+ theme_classic() + theme(text=element_text(size=14), axis.title.x = element_blank(), plot.title=element_text(hjust=0.5))+labs(title=expression(paste(italic("S. enterica"), " LT2")), y="log10(CFU/worm)")

pSeByDay<- SeBleachCount %>%
	ggplot(aes(x=Date, y=logCount, color=Date)) + 
	geom_jitter(shape=16, position=position_jitter(0.05)) +
	geom_violin(fill=NA) + 
	ylim(-0.1,6)+ theme_classic() + 
	theme(
	text=element_text(size=14), 
	axis.title.x = element_blank(), 
	plot.title=element_text(hjust=0.5,size=14)) +
	labs(title=expression(paste(italic("S. enterica"), " LT2")), y="log10(CFU/worm)")
pSeByDay
#ggsave("pSeByDay.png",width=4, height=3, units="in", dpi=300)

# here we will call the wormboot function (wormboot.R)
# to generate simulated data and make comparisons
# 25 reps as usual
SeBoot1<-wormboot(25, SeBleachCount$Count[SeBleachCount$Date=="1"])
SeBoot2<-wormboot(25, SeBleachCount$Count[SeBleachCount$Date=="2"])
SeBoot1$Day<-as.factor("Day1")
SeBoot2$Day<-as.factor("Day2")
jointSeBoot<-rbind(SeBoot1, SeBoot2)
SeBleachCount$Date<-as.factor(SeBleachCount$Date)

pBoot1<- SeBleachCount %>%
	ggplot(aes(x=Date, y=logCount, color=Date)) + 
	geom_jitter(shape=16, position=position_jitter(0.05)) +
	geom_violin(fill=NA) + 
	ylim(-0.1,6)+ theme_classic() + 
	theme(
	text=element_text(size=14), 
	axis.title.x = element_blank(),
	axis.text.x = element_blank(),
	plot.title=element_text(hjust=0.5,size=14),
	legend.position="none") +
	labs(title="Raw data", y="log10(CFU/worm)")
pBoot2<-subset(jointSeBoot, batch==5) %>%
	ggplot(aes(x=Day, y=data, color=Day)) + 
	geom_jitter(shape=16, position=position_jitter(0.05)) +
	geom_violin(fill=NA) + ylim(-0.1,6)+ theme_classic() + 
	theme(text=element_text(size=16), 
	axis.title.x = element_blank(), 
	axis.text.x = element_blank(),
	plot.title=element_text(hjust=0.5, size=14),
	legend.position="none") + 
	labs(title="Batch 5", y="")
#pBoot2
pBoot3<-subset(jointSeBoot, batch==10) %>%
	ggplot(aes(x=Day, y=data, color=Day)) + 
	geom_jitter(shape=16, position=position_jitter(0.05)) +
	geom_violin(fill=NA) + ylim(-0.1,6)+ theme_classic() + 
	theme(text=element_text(size=16), 
	axis.title.x = element_blank(), 
	axis.text.x = element_blank(),
	plot.title=element_text(hjust=0.5, size=14),
	legend.position="none") + 
	labs(title="Batch 10", y="")
#pBoot3
pBoot4<-subset(jointSeBoot, batch==20) %>%
	ggplot(aes(x=Day, y=data, color=Day)) + 
	geom_jitter(shape=16, position=position_jitter(0.05)) +
	geom_violin(fill=NA) + ylim(-0.1,6)+ theme_classic() + 
	theme(text=element_text(size=16), 
	axis.title.x = element_blank(), 
	axis.text.x = element_blank(),
	plot.title=element_text(hjust=0.5, size=14),
	legend.position=c(0.7,0.2)) + 
	labs(title="Batch 20", y="")
#pBoot4
pBoot5<-subset(jointSeBoot, batch==50) %>%
	ggplot(aes(x=Day, y=data, color=Day)) + 
	geom_jitter(shape=16, position=position_jitter(0.05)) +
	geom_violin(fill=NA) + ylim(-0.1,6)+ theme_classic() + 
	theme(text=element_text(size=16), 
	axis.title.x = element_blank(), 
	axis.text.x = element_blank(),
	plot.title=element_text(hjust=0.5, size=14),
	legend.position=c(0.7,0.2)) + 
	labs(title="Batch 50", y="")
#pBoot5
plot_grid(pBoot1, pBoot2, pBoot3, pBoot4, labels="AUTO", ncol=4, align="h")
#ggsave("pjointSeBoot.svg", width=9, height=4, units="in", dpi=300)
ggsave("pjointSeBoot.png", width=9, height=4, units="in", dpi=300)

SeBoot1$Count<-as.double(10^SeBoot1$data)
SeBoot2$Count<-as.double(10^SeBoot2$data)

#and some statistical tests on the bootstrapped data
t.test(SeBoot1$Count[SeBoot1$batch==5], SeBoot2$Count[SeBoot2$batch==5])
t.test(SeBoot1$Count[SeBoot1$batch==10], SeBoot2$Count[SeBoot2$batch==10])
t.test(SeBoot1$Count[SeBoot1$batch==20], SeBoot2$Count[SeBoot2$batch==20])
t.test(SeBoot1$Count[SeBoot1$batch==50], SeBoot2$Count[SeBoot2$batch==50])
wilcox.test(SeBoot1$Count[SeBoot1$batch==5], SeBoot2$Count[SeBoot2$batch==5])
wilcox.test(SeBoot1$Count[SeBoot1$batch==10], SeBoot2$Count[SeBoot2$batch==10])
wilcox.test(SeBoot1$Count[SeBoot1$batch==20], SeBoot2$Count[SeBoot2$batch==20])
wilcox.test(SeBoot1$Count[SeBoot1$batch==50], SeBoot2$Count[SeBoot2$batch==50])
t.test(SeBoot1$data[SeBoot1$batch==5], SeBoot2$data[SeBoot2$batch==5])
t.test(SeBoot1$data[SeBoot1$batch==10], SeBoot2$data[SeBoot2$batch==10])
t.test(SeBoot1$data[SeBoot1$batch==20], SeBoot2$data[SeBoot2$batch==20])
t.test(SeBoot1$data[SeBoot1$batch==50], SeBoot2$data[SeBoot2$batch==50])
wilcox.test(SeBoot1$data[SeBoot1$batch==5], SeBoot2$data[SeBoot2$batch==5])
wilcox.test(SeBoot1$data[SeBoot1$batch==10], SeBoot2$data[SeBoot2$batch==10])
wilcox.test(SeBoot1$data[SeBoot1$batch==20], SeBoot2$data[SeBoot2$batch==20])
wilcox.test(SeBoot1$data[SeBoot1$batch==50], SeBoot2$data[SeBoot2$batch==50])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~    Build Figure 1
# Plot out SE bootstraps together with real data from S. aureus batch size experiment
plot_grid(pBatchSA, pBoot1, pBoot2, pBoot3, pBoot4, labels="AUTO", ncol=5,  rel_widths=c(1.5, 1,1,1,1), align="h")
#ggsave("pjointSeBoot.svg", width=9, height=4, units="in", dpi=300)
ggsave("pSaBatchJointSeBoot.png", width=10, height=4, units="in", dpi=300)

#pBatchSE<-subset(BatchDigests, Species=="SE") %>%
#	ggplot(aes(x=Batch, y=logCFU, color=Batch)) + 
#	geom_jitter(shape=16, position=position_jitter(0.05)) +
#	geom_violin(fill=NA) + 
#	theme_classic() + 
#	theme(text=element_text(size=16), 
#	axis.title.x = element_blank(), 
#	plot.title=element_text(hjust=0.5, size=14),
#	legend.position=c(0.8,0.2)) + 
#	labs(title=expression(paste(italic("S. enterica"), " LT2")), y="log10(CFU/worm)")
#pBatchSE

mean(BatchDigests$CFU[BatchDigests$Species=="SA" & BatchDigests$Batch==1 ])
median(BatchDigests$CFU[BatchDigests$Species=="SA" & BatchDigests$Batch==1 ])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# S. enterica and S. aureus bleach experiment data and simulations thereof
library(mclust, quietly=TRUE)
library(sBIC)

# let's combine the SE data with the single worm data (n=22) from the batch digest run
# now we have three replicates on different days
temp<-subset(BatchDigests, Species=="SE" & Batch==1)
#length(temp$Batch)
tempSE<-data.frame(Condition=rep("Bleach", 22), Date=rep("20220218", 22), Count = temp$CFU, logCount=temp$logCFU, Species=rep("SE", 22))
SeCount<-rbind(SeBleachCount, tempSE)
rm(tempSE)
unique(SeCount$Date)
SeCount[SeCount=="20210308"]<-2
SeCount[SeCount=="20210315"]<-1
SeCount[SeCount=="20220218"]<-3
SeCount$Date<-as.factor(SeCount$Date)

#take out the zeros
SeCount[SeCount==0]<-NA
SeCount2<-SeCount[complete.cases(SeCount),]

# do the same for the S. aureus data
SaBleachCount[SaBleachCount=="20210212"]<-1
SaBleachCount[SaBleachCount=="20210207"]<-2
SaBleachCount[SaBleachCount=="20210623"]<-3
SaBleachCount$Date<-as.factor(SaBleachCount$Date)
# censor out zeros
SaBleachCount2<-SaBleachCount
SaBleachCount2[SaBleachCount2==0]<-NA
SaBleachCount2<-SaBleachCount[complete.cases(SaBleachCount2),]

pSEsingle<-SeCount2 %>%
	ggplot(aes(x=Date, y=logCount, color=Date)) + 
	geom_jitter(shape=16, position=position_jitter(0.05)) +
	geom_violin(fill=NA) + 
	theme_classic() + 
	theme(text=element_text(size=16), 
	axis.title.x = element_blank(), 
	axis.text.x = element_text(size=12),
	plot.title=element_text(hjust=0.5, size=14),
	legend.position = "none") + 
	labs(title=expression(paste(italic("S. enterica"), " LT2")), y="log10(CFU/worm)")
pSEsingle
ggsave("pSEsingle.png", width=4, height=3, units="in", dpi=300)

#another version of this plot with the combined data (all 3 days)
mylen<-dim(SeCount2)[1]
SeTemp<-rbind(SeCount2, data.frame(Condition=SeCount2$Condition,
                                     Date=rep("All", mylen),
                                     Count=SeCount2$Count,
                                     logCount=SeCount2$logCount,
                                     Species=SeCount2$Species
                                     ))
pSeCountAll<-ggplot(SeTemp, aes(x=Date, y=logCount, color=Date)) +
    geom_jitter(shape=16, position=position_jitter(0.05)) +
    geom_violin(fill=NA) + 
    theme_classic() + 
    theme(text=element_text(size=14), 
          axis.title.x = element_blank(), 
          axis.text.x = element_blank(),
          plot.title=element_text(hjust=0.5, size=14)) + 
    labs(title=expression(paste(italic("S. enterica"), " LT2")), y="log10(CFU/worm)")
pSeCountAll
ggsave("pSECountAll.png", width=4, height=3, units="in", dpi=300)

#ok let's try some gaussian mixture models (code from vignette)
X<-SeCount2$logCount
mean(X)
var(X)
mean(SeCount2$Count)
sd(SeCount2$Count)

fit<- mclustBIC(X)
# Best BIC values:
#           V,2        E,2         E,3
#BIC      -251.966 -252.38286 -258.396880
#BIC diff    0.000   -0.41682   -6.430838

plot(fit)

# fitting a model with two components, allowing variances to be unequal ("V")
fit2 <- Mclust(X, G=2, model="V")
summary(fit2)
plot(fit2, what="density", main="G2", xlab="logCFU")
rug(X)
fit2$parameters
#(mean1, var1)=(2.985569, 0.4568056) and (mean2, var2)=(4.742767, 0.1239928), containing 58% and 42% of the mass respectively
# let's put these plots together with the Se data plot

X1<-SeCount2$logCount[SeCount2$Date==1]
fit1 = Mclust(X1, G=2, model="V")
summary(fit1)
plot(fit1, what="density", main="G2", xlab="logCFU")
rug(X1)
fit1$parameters
#(mean1, var1)=(2.777394, 0.34993190) and (mean2, var2)=(4.753665, 0.05852001), containing 59.4% and 40.6% of the mass respectively

X2<-SeCount2$logCount[SeCount2$Date==2]
fit2 = Mclust(X2, G=2, model="V")
fit2$parameters
#(mean1, var1)=(2.774017, 0.09134093) and (mean2, var2)=(4.319905, 0.32749870), containing 42.6% and 57.3% of the mass respectively

X3<-SeCount2$logCount[SeCount2$Date==3]
fit3 = Mclust(X3, G=2, model="V")
#summary(fit3)
#plot(fit3, what="density", main="G2", xlab="logCFU")
#rug(X3)
fit3$parameters
#(mean1, var1)=(3.680638, 0.97808089) and (mean2, var2)=(4.905500, 0.03381337), containing 67.2% and 32.8% of the mass respectively

#library(sBIC)
gMix = GaussianMixtures(maxNumComponents=10, phi=1, restarts=100)
set.seed(1234)
m = sBIC(X, gMix)
print(m)
matplot(
	cbind(m$BIC - m$BIC[1], m$sBIC - m$sBIC[1]),
	pch = c(1, 3),
	col = "black",
	xlab = "Number of components",
	ylab = expression(BIC - BIC(M[1])),
	las=1, xaxt="n"
	)
axis(1, at = 1:10)
legend("bottomleft",
c(expression(BIC), expression(bar(sBIC)[1])),
pch = c(1, 3),
y.intersp = 1.2)
#ggsave("pSE_GMM_BICvsBIC.png", width=4, height=3, units="in", dpi=300)

# Using the all-SE 2-mode GMM, we can generate simulated data
# Note that we are allowing the weight in modes to be stochastic (binomial) 
SeBootGMM<-wormbootGMM(100, 50, 2.985569, 0.4568056, 4.742767, 0.1239928, 0.58, 0)
SeBootGMM$set<-as.factor(SeBootGMM$set)
SeBootGMM$batch<-as.factor(SeBootGMM$batch)

SeBootGMM.1.2<-wormbootGMM(100, 50, 2.985569, 0.4568, 4.742767, 0.12399, 0.57, 0.15)
SeBootGMM.1.2$set<-as.factor(SeBootGMM.1.2$set)
SeBootGMM.1.2$batch<-as.factor(SeBootGMM.1.2$batch)

SeBootGMM.1.2.05<-wormbootGMM(100, 50, 2.985569, 0.4568, 4.742767, 0.12399, 0.57, 0.05)

pSimSeGMM1<- subset(SeBootGMM.1.2, batch==1) %>%
	ggplot(aes(x=set, y=logCount, color=set))+
	geom_jitter(shape=16, position=position_jitter(0.05)) +
	geom_violin(fill=NA) + 
	ylim(-0.1,6) +
	theme_classic() + 
	theme(text=element_text(size=16), 
	axis.title.x = element_blank(), 
	plot.title=element_text(hjust=0.5, size=16),
	legend.position="none"
	) + 
	labs(title="1", y="log10(CFU/worm)")
#pSimSeGMM1
pSimSeGMM5<- subset(SeBootGMM.1.2, batch==5) %>%
	ggplot(aes(x=set, y=logCount, color=set))+
	geom_jitter(shape=16, position=position_jitter(0.05)) +
	geom_violin(fill=NA) + 
	ylim(-0.1,6) +
	theme_classic() + 
	theme(text=element_text(size=16), 
	axis.title.x = element_blank(), 
	plot.title=element_text(hjust=0.5, size=16),
	legend.position="none"
	) + 
	labs(title="5", y="")
pSimSeGMM10<- subset(SeBootGMM.1.2, batch==10) %>%
	ggplot(aes(x=set, y=logCount, color=set))+
	geom_jitter(shape=16, position=position_jitter(0.05)) +
	geom_violin(fill=NA) + 
	ylim(-0.1,6) +
	theme_classic() + 
	theme(text=element_text(size=16), 
	axis.title.x = element_blank(), 
	plot.title=element_text(hjust=0.5, size=16),
	legend.position="none"
	) + 
	labs(title="10", y="")
pSimSeGMM20<- subset(SeBootGMM.1.2, batch==20) %>%
	ggplot(aes(x=set, y=logCount, color=set))+
	geom_jitter(shape=16, position=position_jitter(0.05)) +
	geom_violin(fill=NA) + 
	ylim(-0.1,6) +
	theme_classic() + 
	theme(text=element_text(size=16), 
	axis.title.x = element_blank(), 
	plot.title=element_text(hjust=0.5, size=16),
	legend.position="none"
	) + 
	labs(title="20", y="")
pSimSeGMM50<- subset(SeBootGMM.1.2, batch==50) %>%
	ggplot(aes(x=set, y=logCount, color=set))+
	geom_jitter(shape=16, position=position_jitter(0.05)) +
	geom_violin(fill=NA) + 
	ylim(-0.1,6) +
	theme_classic() + 
	theme(text=element_text(size=16), 
	axis.title.x = element_blank(), 
	plot.title=element_text(hjust=0.5, size=16),
	legend.position=c(0.8,0.2)) + 
	labs(title="50", y="")
plot_grid(pSimSeGMM1, pSimSeGMM5, pSimSeGMM10, pSimSeGMM20, pSimSeGMM50, labels="AUTO", ncol=5, align="h")
ggsave("pSimSeGMM_day1v2.png", width=9, height=4, units="in", dpi=300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~    WITHIN AND BETWEEN REPLICATE HETEROGENEITY OF CFUS
#~~~~
# pull in the all-8 data set so we can describe run to run variation
# just need the total counts for this

All8<-read.table("All8.txt", header=TRUE)
names(All8)
hosts<-unique(All8$Host)
All8$logTotal<-log10(All8$Total)
All8[All8=="decc1"]<-"dec1"
All8$logTotal[is.infinite(All8$logTotal)]<-NA
All8<-All8[complete.cases(All8),]

# find hosts with sufficient runs
temp<-integer()
for (i in 1:length(hosts)){
	temp[i]<-length(unique(All8$Date[All8$Host==hosts[i]]))}
temp
#5  2  4  4  2  5  6  5  1  1  1 11  5  3
hosts
#"AU37"   "ctls40" "daf16"  "daf2"   "dbl1"   "decc1"  "eat14"  "exp1"   "glp1"   "glp4"   "hsf1"   "N2"  "phm2"   "vhp1" 
bighosts<-hosts[temp>=5]
#"AU37"  "daf2"  "dec1" "eat14" "exp1"  "N2"    "phm2" 
# these have at least five runs
All8b<-subset(All8, Host=="AU37" | Host=="daf2" | Host=="dec1" | Host=="eat14" | Host=="exp1" | Host=="N2" | Host=="phm2")
All8b$Date<-as.factor(All8b$Date)

pAll8b<-All8b %>%
	ggplot(aes(x=as.factor(Host), y=logTotal, color=as.factor(Host)))+
	geom_jitter(shape=16, position=position_jitter(0.05)) +
	geom_violin(fill=NA) + 	ylim(1,6) +	theme_classic() + 
	theme(text=element_text(size=14), axis.title.x = element_blank(), legend.position="none", plot.title=element_text(hjust=0.5, size=16)) +
	labs(title="Minimal Microbiome Totals", y="log10(CFU/Worm)")
pAll8b
#ggsave("pAll8LogTotal.png", width=6, height=4, units="in", dpi=300)

# Plot out individual replicates by host genotype
All8b%>%
  ggplot(aes(x=Date, y=logTotal, color=Date)) + geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + 	ylim(1,6) +	theme_classic() + 
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        legend.position="none", 
        plot.title=element_text(hjust=0.5, size=16)) +
  labs(title="Minimal community colonization, by replicate", y=expression(log[10](CFU/Worm))) +
  facet_wrap(~Host, scales="free_x", ncol=4)
#ggsave("pAll8LogCFUHostByDay.png", width=7, height=10, units="in", dpi=300)

# turns out a lot of the early runs were low n
# not great for quantifying heterogeneity
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
pAll8.big.byHost<-All8.big%>%
  ggplot(aes(x=Date, y=logTotal, color=Date)) + geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + 	ylim(1,6) +	theme_classic() + 
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        legend.position="none", 
        plot.title=element_text(hjust=0.5, size=16)) +
  labs(title="Minimal community colonization, by replicate", y=expression(log[10](CFU/Worm))) +
  facet_wrap(~Host, scales="free_x", ncol=4)

pAll8.big<-All8.big %>%
	ggplot(aes(x=as.factor(Host), y=logTotal, color=as.factor(Host)))+
	geom_jitter(shape=16, position=position_jitter(0.05)) +
	geom_violin(fill=NA) + 	ylim(1,6) +	theme_classic() + 
	theme(text=element_text(size=14), axis.title.x = element_blank(), legend.position="none", plot.title=element_text(hjust=0.5, size=16)) +
	labs(title="Minimal Microbiome Totals", y=expression(log[10](CFU/Worm)))

plot_grid()
ggsave("pAll8BigLogCFUHostByDay.png", width=12, height=9, units="in", dpi=300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########## OK now we are going to take a WHOLE BUNCH of summary stats and see what we get

AllCounts<-as.tibble(subset(All8.big, select=c("Host", "Date", "Total", "logTotal")))
names(AllCounts)

# Incorporate the pathogen colonization data
SeCount<-SeCount[complete.cases(SeCount),]
SeCount$Host<-"SE"
SaBleachCount2$Host<-"SA"
SeCount2<- subset(SeCount, select=c("Host", "Date", "Count", "logCount")) %>%
  rename("Total"="Count",
         "logTotal"="logCount") %>%
  as_tibble()
SaCount2<-subset(SaBleachCount2, select=c("Host", "Date", "Count", "logCount")) %>%
  rename("Total"="Count",
         "logTotal"="logCount") %>%
  as_tibble()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2022-06-03: bringing in single worm data from NRRL1
# note: took out CG mono-colonization data, too sparse
# and CG pairs, most had only one day & all had few points

NRRL1counts<-read_csv("NRRL1counts.csv")
names(NRRL1counts)
#[1] "Condition" "CFU"       "logCFU"    "Date"
NRRL1counts$Date<-as.character(NRRL1counts$Date)

#as before loop over experiment (Condition) and rep (Date) and record
hosts2<-unique(NRRL1counts$Condition) #length 24

#Only SSP colonization by MO.I2.9 (hosts[10]) has 2 reps of >18 (both n=24)
hosts2[10]
MOcount<-subset(NRRL1counts, Condition==hosts2[10], select=c("Condition", "Date", "CFU", "logCFU")) %>%
  rename("Host"="Condition",
         "Total"="CFU",
         "logTotal"="logCFU") %>%
  as_tibble()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Merge everything
AllCounts<-rbind(AllCounts, SaCount2, SeCount2, MOcount) %>%
  rename("Condition"="Host")

# Get a list of all experiments
conditions.big<-unique(AllCounts$Condition)

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
	temp<-unique(AllCounts$Date[AllCounts$Condition==conditions.big[i]])
#	print(conditions.big[i])
	for (j in seq_along(temp)){
#		print(temp[j])
		tempCFU<-sort(AllCounts$Total[AllCounts$Condition==conditions.big[i] & AllCounts$Date==temp[j]])
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's plot out

AllCounts %>%
  ggplot(aes(x=Date, y=logTotal, color=Date)) + 
  geom_point(shape=16, position=position_jitterdodge(0.05)) +
  geom_violin(fill=NA) + 
  theme_classic() + 
  theme(text=element_text(size=16), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=14),
        legend.position = "none") + 
  facet_wrap(~Condition, scales="free_x")
  
#~~~~~ FIGURE 2
ggsave("pAll8BigSaSeCFUHostByDay.png", width=12, height=10, units="in", dpi=300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~ Plot out the moments by condition
names(AllCountsStats.r2r)

#some of these stats are hard to see on linear axes; mutate to plot
AllCountStats.r2r %>%
  mutate(logmeanCFU=log10(meanCFU),
         logq50=log10(q50))

AllCountStats.r2r %>% 
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
  facet_wrap(vars(Var), scales="free_y")
  
ggsave("pComparingMomentsRunToRunWormCFU.png", width=14, height=10, units="in", dpi=300)

### cool. Now let's do the pairwise comparisons
# actually I think I want this as a data frame
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

#names(All8.big.r2r.t)

#~~~~~~~~~~~~~~~~~~~
# now let's look at the differences in statistics

AllCountStats.r2r.dist <- tibble(meanCFUdist=numeric(),
                                 meanCFUCohenD=numeric(),
                                 mediandist=numeric(),
                                 cvdist=numeric(),
                                 skewdist=numeric(),
                                 kurtdist=numeric(),
                                 Q1dist=numeric(),
                                 Q2dist=numeric())

#~~~~~~~~~   NOT WORKING YET
#npairs<-2 #for testing
for(i in 1:(npairs-1)){
	#print(i)
	for (j in (i+1):npairs){
		#print(j)
	  bind_rows(AllCountStats.r2r.dist, tibble(
		meanCFUdist=abs(AllCountStats.r2r.t$meanCFU[i]-AllCountStats.r2r.t$meanCFU[j])/(AllCountStats.r2r.t$meanCFU[i]+AllCountStats.r2r.t$meanCFU[j]),
	  meanCFUCohenD=abs(AllCountStats.r2r.t$meanCFU[i]-AllCountStats.r2r.t$meanCFU[j])/sqrt(min(AllCountStats.r2r.t$varCFU[i],AllCountStats.r2r.t$varCFU[i])),
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

# plot out 
AllCountStats.r2r.dist

pAll8MeanCFUdist<-AllCountStats.r2r.dist %>%
  ggplot(aes(x=same, y=meanCFUdist, color=same))+
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_boxplot(fill=NA) +	theme_classic() + 
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(), 
        legend.position="none", 
        plot.title=element_text(hjust=0.5, size=16)) +
  stat_compare_means(label.y = 0.99, label.x=1.2) +
  labs(title="", y="dist(Mean CFU)")
pAll8MeanCFUdist
ggsave("pAll8MeanCFUdist.png", width=6, height=4, units="in", dpi=300)


pAll8cvCFUdist<-All8.big.r2r.dist %>%
  ggplot(aes(x=same, y=cvdist, color=same))+
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_boxplot(fill=NA) +	theme_classic() + 
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(), 
        legend.position="none", 
        plot.title=element_text(hjust=0.5, size=16)) +
  stat_compare_means(label.y = 2.6, label.x=1.5) +
  labs(title="", y="dist(CV)")
pAll8cvCFUdist

pAll8SkewCFUdist<-All8.big.r2r.dist %>%
  ggplot(aes(x=same, y=skewdist, color=same))+
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_boxplot(fill=NA) +	theme_classic() + 
  theme(text=element_text(size=14), axis.title.x = element_blank(), legend.position="none", plot.title=element_text(hjust=0.5, size=16)) +
  stat_compare_means(label.y = 4.3, label.x=1.5) +
  labs(title="", y="dist(Skew)")
pAll8SkewCFUdist

pAll8KurtCFUdist<-All8.big.r2r.dist %>%
  ggplot(aes(x=same, y=kurtdist, color=same))+
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_boxplot(fill=NA) +	theme_classic() + 
  theme(text=element_text(size=14), axis.title.x = element_blank(), legend.position="none", plot.title=element_text(hjust=0.5, size=16)) +
  stat_compare_means(label.y = 17, label.x=1.5) +
  labs(title="", y="dist(Kurt)")
pAll8KurtCFUdist

#~~~~~~   FIGURE 3?
#~~~  Plot out only the comparisons of moments
plot_grid(pAll8MeanCFUdist, pAll8MeanCFUCohenD, pAll8q50CFUdist,
          pAll8cvCFUdist, pAll8SkewCFUdist, pAll8KurtCFUdist, 
          labels="AUTO", ncol=2)
ggsave("pMomentsCFUdist.png", width=7, height=9, units="in", dpi=300)
#ggsave("pAll8MeanCFUdist.png", width=6, height=4, units="in", dpi=300)

## Plotting summary statistics along with their comparisons
plot_grid(pmeanCFU, pcvCFU, pskewCFU, pkurtCFU, 
          pAll8MeanCFUdist, pAll8cvCFUdist, pAll8SkewCFUdist, pAll8KurtCFUdist, 
          ncol=4, nrow=2, labels="AUTO")
ggsave("pAll8SelectMomentsSummaryDistance.png", width=12, height=6, units="in", dpi=300)

wilcox.test(meanCFUdist~same, data=All8.big.r2r.dist)
wilcox.test(cvdist~same, data=All8.big.r2r.dist)
wilcox.test(skewdist~same, data=All8.big.r2r.dist)
wilcox.test(kurtdist~same, data=All8.big.r2r.dist)

#####################################################################
####################################################################
###########      FIGURE 5  MOMENTS      ##############################
# ok back to beta starting at B(1,5) why not
# without much loss of generality let's draw factors from a N(0,1)
# each run
# function is in wormSimBatchBetaRand.R

wormSimBetaRand.1.5<-wormSimBatchBetaRand(1,5,1000,24,100000)
wormSimBetaRand.p5.2p5<-wormSimBatchBetaRand(0.5,2.5,1000,24,100000)
wormSimBetaRand.5.1<-wormSimBatchBetaRand(5,1,1000,24,100000)
wormSimBetaRand.2p5.p5<-wormSimBatchBetaRand(2.5,0.5,1000,24,100000)
wormSimBetaRand.5.5<-wormSimBatchBetaRand(5,5,1000,24,100000)
wormSimBetaRand.2p5.2p5<-wormSimBatchBetaRand(2.5,2.5,1000,24,100000)
wormSimBetaRand.1.1<-wormSimBatchBetaRand(1,1,1000,24,100000)
#mgenbeta(3, 1, 5, 1)
#mgenbeta(3, 0.5, 2.5, 1)

#plot out
wormSimBetaRand.1.5$meandist<-abs(wormSimBetaRand.1.5$meanA-wormSimBetaRand.1.5$meanB)/(wormSimBetaRand.1.5$meanA+wormSimBetaRand.1.5$meanB)
pSimBatchBetaRand.1.5.meanA<-wormSimBetaRand.1.5 %>%
  ggplot(aes(x=batch, y=log10(meanA), color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  ylim(3.7,4.5)+
  geom_jitter(shape=16, position=position_jitter(0.05)) +
    theme(text=element_text(size=14), 
          axis.title.x = element_blank(),  legend.position="none",
          plot.title=element_text(hjust=0.5, size=14)) + 
    labs(title=paste("\u03b2","(1,5) A", sep=""), y="Mean")
#pSimBatchBetaRand.1.5.meanA
pSimBatchBetaRand.1.5.meanB<-wormSimBetaRand.1.5 %>%
  ggplot(aes(x=batch, y=log10(meanB), color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(),  legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title=paste("\u03b2","(1,5) B", sep=""), y="Mean")
pSimBatchBetaRand.1.5.meandist<-wormSimBetaRand.1.5 %>%
  ggplot(aes(x=batch, y=meandist, color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  theme(text=element_text(size=14), 
        #axis.title.x = element_blank(), 
        legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title="Distance", y="Dist(Mean)")
plot_grid(pSimBatchBetaRand.1.5.meanA, pSimBatchBetaRand.1.5.meanB, pSimBatchBetaRand.1.5.meandist, ncol=3)

wormSimBetaRand.1.5$cvA<-sqrt(wormSimBetaRand.1.5$varA)/wormSimBetaRand.1.5$meanA
wormSimBetaRand.1.5$cvB<-sqrt(wormSimBetaRand.1.5$varB)/wormSimBetaRand.1.5$meanB
wormSimBetaRand.1.5$cvdist<-abs(wormSimBetaRand.1.5$cvA-wormSimBetaRand.1.5$cvB)
pSimBatchBetaRand.1.5.cvA<-wormSimBetaRand.1.5 %>%
  ggplot(aes(x=batch, y=cvA, color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  ylim(0, 1.5)+
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(),  legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title="", y="CV")
pSimBatchBetaRand.1.5.cvB<-wormSimBetaRand.1.5 %>%
  ggplot(aes(x=batch, y=cvB, color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(),  legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title="", y="CV")
pSimBatchBetaRand.1.5.cvdist<-wormSimBetaRand.1.5 %>%
  ggplot(aes(x=batch, y=cvdist, color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  ylim(0,1)+
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(), legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title="", y="Dist(CV)")
plot_grid(pSimBatchBetaRand.1.5.cvA, pSimBatchBetaRand.1.5.cvB, pSimBatchBetaRand.1.5.cvdist, ncol=3)

wormSimBetaRand.1.5$skewdist<-abs(wormSimBetaRand.1.5$skewA-wormSimBetaRand.1.5$skewB)
pSimBatchBetaRand.1.5.skewA<-wormSimBetaRand.1.5 %>%
  ggplot(aes(x=batch, y=skewA, color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(),  legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title="", y="Skewness")
pSimBatchBetaRand.1.5.skewB<-wormSimBetaRand.1.5 %>%
  ggplot(aes(x=batch, y=skewB, color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(),  legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title="", y="Skewness")
pSimBatchBetaRand.1.5.skewdist<-wormSimBetaRand.1.5 %>%
  ggplot(aes(x=batch, y=skewdist, color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(), legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title="", y="Dist(skew)")
plot_grid(pSimBatchBetaRand.1.5.skewA, pSimBatchBetaRand.1.5.skewB, pSimBatchBetaRand.1.5.skewdist, ncol=3)

wormSimBetaRand.1.5$kurtdist<-abs(wormSimBetaRand.1.5$kurtA-wormSimBetaRand.1.5$kurtB)
pSimBatchBetaRand.1.5.kurtA<-wormSimBetaRand.1.5 %>%
  ggplot(aes(x=batch, y=kurtA, color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  theme(text=element_text(size=14), 
        #axis.title.x = element_blank(),  
        legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title="", y="Kurtosis")
pSimBatchBetaRand.1.5.kurtB<-wormSimBetaRand.1.5 %>%
  ggplot(aes(x=batch, y=kurtB, color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  theme(text=element_text(size=14), 
        #axis.title.x = element_blank(),  
        legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title="", y="Kurtosis")
pSimBatchBetaRand.1.5.kurtdist<-wormSimBetaRand.1.5 %>%
  ggplot(aes(x=batch, y=kurtdist, color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  theme(text=element_text(size=14), 
        #axis.title.x = element_blank(), 
        legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title="", y="Dist(kurtosis)")
plot_grid(pSimBatchBetaRand.1.5.kurtA, pSimBatchBetaRand.1.5.kurtB, pSimBatchBetaRand.1.5.kurtdist, ncol=3)

plot_grid(pSimBatchBetaRand.1.5.meanA, 
          pSimBatchBetaRand.1.5.meanB, 
          pSimBatchBetaRand.1.5.meandist, 
          pSimBatchBetaRand.1.5.cvA, 
          pSimBatchBetaRand.1.5.cvB, 
          pSimBatchBetaRand.1.5.cvdist, 
          pSimBatchBetaRand.1.5.skewA, 
          pSimBatchBetaRand.1.5.skewB, 
          pSimBatchBetaRand.1.5.skewdist, 
          pSimBatchBetaRand.1.5.kurtA, 
          pSimBatchBetaRand.1.5.kurtB, 
          pSimBatchBetaRand.1.5.kurtdist,
          nrow=4,
          ncol=3)
ggsave("pSimBatchBetaRand.1.5.moments.png", width=12, height=12, units="in")      


#################     FIGURE 6     #############################
#### add data from symmetric distribution for comparison
#note these runs use the flat U(-0.1,0.1) parameter error
wormSimBetaRand.5.5$cvA<-sqrt(wormSimBetaRand.5.5$varA)/wormSimBetaRand.1.5$meanA
wormSimBetaRand.5.5$cvB<-sqrt(wormSimBetaRand.5.5$varB)/wormSimBetaRand.1.5$meanB
wormSimBetaRand.5.5$cvdist<-abs(wormSimBetaRand.5.5$cvA-wormSimBetaRand.1.5$cvB)

pSimBatchBetaRand.5.5.cvA<-wormSimBetaRand.5.5 %>%
  ggplot(aes(x=batch, y=cvA, color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  ylim(-0.01, 1.5)+
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(),  legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title=paste("\u03b2","(5,5)", sep=""), y="CV")
pSimBatchBetaRand.5.5.cvA
pSimBatchBetaRand.5.5.cvdist<-wormSimBetaRand.5.5 %>%
  ggplot(aes(x=batch, y=cvdist, color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  ylim(0, 1)+
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(), legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title=paste("\u03b2","(5,5)", sep=""), y="Dist(CV)")
pSimBatchBetaRand.5.5.cvdist
#plot_grid(pSimBatchBetaRand.1.5.cvA, pSimBatchBetaRand.5.5.cvA, pSimBatchBetaRand.1.5.cvdist, pSimBatchBetaRand.5.5.cvdist, nrow=2, ncol=2)

wormSimBetaRand.5.5$meandist<-abs(wormSimBetaRand.5.5$meanA-wormSimBetaRand.5.5$meanB)/(wormSimBetaRand.5.5$meanA+wormSimBetaRand.5.5$meanB)
pSimBatchBetaRand.5.5.meanA<-wormSimBetaRand.5.5 %>%
  ggplot(aes(x=batch, y=log10(meanA), color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  ylim(4.1, 5)+
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(),  
        legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title=paste("\u03b2","(5,5) A", sep=""), y="Mean")
pSimBatchBetaRand.5.5.meandist<-wormSimBetaRand.5.5 %>%
  ggplot(aes(x=batch, y=meandist, color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  ylim(-0.01,0.45)+
  theme(text=element_text(size=14), 
        #axis.title.x = element_blank(), 
        legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title="Distance", y="Dist(mean)")
plot_grid(pSimBatchBetaRand.1.5.meanA, pSimBatchBetaRand.5.5.meanA,pSimBatchBetaRand.1.5.meandist, pSimBatchBetaRand.5.5.meandist, ncol=2, nrow=2, labels="AUTO", align = "h")
ggsave("pSimBatchBetaRand.1.5.vs.5.5.means.png", width=8, height=6, units="in")

#simBatchBeta.12.10<-wormSimBatchBetaCompare(1, 11, 1, 9, 1000, 25, 100000)
#pSimBatchBeta.12.10.b1<- subset(simBatchBeta.12.10, batch==1) %>%
#  ggplot(aes(x=set, y=log10(data), color=set))+
#  geom_jitter(shape=16, position=position_jitter(0.05)) +
#  geom_violin(fill=NA) + theme_classic() + 
#  theme(text=element_text(size=16), 
#        axis.title.x = element_blank(), 
#        plot.title=element_text(hjust=0.5, size=16)) + 
#  labs(title="1", y="log10(CFU/worm)")
#pSimBatchBeta.12.10.b1
#ggsave("pSimBetaRHS_1_5.png", width=6, height=4, units="in", dpi=300)
#skewness(simBatchBeta.12.10$data[simBatchBeta.12.10$batch==1 & simBatchBeta.12.10$set=="A"])
#skewness(simBatchBeta.12.10$data[simBatchBeta.12.10$batch==1 & simBatchBeta.12.10$set=="B"])
#kurtosis(simBatchBeta.12.10$data[simBatchBeta.12.10$batch==1 & simBatchBeta.12.10$set=="A"])
#kurtosis(simBatchBeta.12.10$data[simBatchBeta.12.10$batch==1 & simBatchBeta.12.10$set=="B"])

##### Now more realistic - we allow triplicates, and cap the number of both total worms and digests we will perform
#"small" is 100 worms and max 12 digests
#"medium" 200 worms max 24 digests
#"large" 500 worms max 48 digests

a<-4
b<-4
(a*b)/((a+b+1)*(a+b)^2)
2*(b-a)*sqrt(a+b+1)/(sqrt(a*b)*(a+b+2))

#nWorms<-100
#maxSamples<-12
#maxCFU<-100000

wormSim3.1.5.s<-wormSimBatchBeta3(1,5,100000,100,12,1000)
#wormSim3.1.5.s$set<-as.factor(wormSim3.1.5.s$set)
#wormSim3.1.5.s$logCFU<-log10(wormSim3.1.5.s$CFU)
wormSim3.1.5.m<-wormSimBatchBeta3(1,5,100000,200,24,1000)
wormSim3.1.5.l<-wormSimBatchBeta3(1,5,100000,500,48,1000)

wormSim3.5.5.s<-wormSimBatchBeta3(5,5,100000,100,12,1000)
wormSim3.5.5.m<-wormSimBatchBeta3(5,5,100000,200,24,1000)
wormSim3.5.5.l<-wormSimBatchBeta3(5,5,100000,500,48,1000)

wormSim3.3.7.s<-wormSimBatchBeta3(3,7,100000,100,12,1000)
wormSim3.3.7.m<-wormSimBatchBeta3(3,7,100000,200,24,1000)
wormSim3.3.7.l<-wormSimBatchBeta3(3,7,100000,500,48,1000)

# alternate version of function that returns means within run across three days of data
# wormSimBatchBeta3<-function(a, b, maxCFU, nWorms, maxSamples, reps, runs, sameParams=FALSE, returnMEANS=TRUE)
# here we use enough worms to allow 24 batches even at 50 worms

wormSim1.5.5.m.means<-wormSimBatchBeta3(5,5,100000,1200,24,1,1000)
wormSim1.1.5.m.means<-wormSimBatchBeta3(1,5,100000,1200,24,1,1000)
wormSim1.3.7.m.means<-wormSimBatchBeta3(3,7,100000,1200,24,1,1000)

wormSim3.5.5.m.means<-wormSimBatchBeta3(5,5,100000,1200,24,3,1000)
wormSim3.1.5.m.means<-wormSimBatchBeta3(1,5,100000,1200,24,3,1000)
wormSim3.5.5.m.means$batch<-as.factor(wormSim3.5.5.m.means$batch)
wormSim3.1.5.m.means$batch<-as.factor(wormSim3.1.5.m.means$batch)

wormSim3.1.5.m<-wormSimBatchBeta3(1,5,100000,1200,24,3,1000, returnMEANS=FALSE)
wormSim3.1.5.m$batch<-as.factor(wormSim3.1.5.m$batch)
wormSim3.1.5.m$day<-as.factor(wormSim3.1.5.m$day)

pwormSim3.1.5.m.A<-subset(wormSim3.1.5.m, set=="A") %>%
  ggplot(aes(x=batch, y=logCFU, color=day)) + 
  geom_violin(fill=NA) + 
  geom_point(shape=16, position=position_jitterdodge(0.2)) +
  ylim(1,5) + theme_classic() + 
  theme(text=element_text(size=14), 
        #axis.title.x = element_blank(), 
        plot.title=element_text(hjust=0.5, size=14),
        legend.position="none") + 
  labs(title="beta(1,5) set A", y="log10(CFU/worm)")
#pwormSim3.1.5.m.A
pwormSim3.1.5.m.B<-subset(wormSim3.1.5.m, set=="B") %>%
  ggplot(aes(x=batch, y=logCFU, color=day)) + 
  geom_violin(fill=NA) + 
  geom_point(shape=16, position=position_jitterdodge(0.2)) +
  ylim(1,5) + theme_classic() + 
  theme(text=element_text(size=14), 
        #axis.title.x = element_blank(), 
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title="beta(1,5) set B", y="log10(CFU/worm)")
plot_grid(pwormSim3.1.5.m.A, pwormSim3.1.5.m.B, rel_widths = c(1,1.15))
ggsave("pwormSim3.1.5.examplerun.png", width=8, height=4, units="in")

## plot differences in convergence of means for one-day vs. 3-day data
wormSim3.5.5.m.means$meandist<-abs(wormSim3.5.5.m.means$meanA-wormSim3.5.5.m.means$meanB)/(wormSim3.5.5.m.means$meanA+wormSim3.5.5.m.means$meanB)
pSimBatchBetaRand.5.5.mean3days<-wormSim3.5.5.m.means %>%
  ggplot(aes(x=batch, y=log10(meanA), color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  ylim(4.1,5)+
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(), legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title="Beta(5,5)x3", y="Mean")
pSimBatchBetaRand.5.5.mean3days
pSimBatchBetaRand.5.5.meandist3days<-wormSim3.5.5.m.means %>%
  ggplot(aes(x=batch, y=meandist, color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  ylim(-0.01,0.2)+
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(), legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title="Beta(5,5)x3", y="Dist(Mean)")
pSimBatchBetaRand.5.5.meandist3days
plot_grid(pSimBatchBetaRand.5.5.meanA,pSimBatchBetaRand.5.5.mean3days,pSimBatchBetaRand.5.5.meandist,pSimBatchBetaRand.5.5.meandist3days)
ggsave("pSimBatchBetaRand.5.5.meancompare3days.png", width=6, height=6, units="in")

## The mean converges nicely, and the distance between means is low, in the 3-day data

wormSim3.5.5.m.means$cvA<-sqrt(wormSim3.5.5.m.means$varA)/wormSim3.5.5.m.means$meanA
wormSim3.5.5.m.means$cvB<-sqrt(wormSim3.5.5.m.means$varB)/wormSim3.5.5.m.means$meanB
pSimBatchBetaRand.5.5.cv3days<-wormSim3.5.5.m.means %>%
  ggplot(aes(x=batch, y=cvA, color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(), legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title="Beta(5,5)x3", y="CV")
pSimBatchBetaRand.5.5.cv3days
plot_grid(pSimBatchBetaRand.5.5.cvA,pSimBatchBetaRand.5.5.cv3days)

pSimBatchBetaRand.1.5.mean3days<-wormSim3.1.5.m.means %>%
  ggplot(aes(x=batch, y=log10(meanA), color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  ylim(3.7,4.4)+
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(), legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title="Beta(1,5)x3", y="Mean")
plot_grid(pSimBatchBetaRand.1.5.meanA, pSimBatchBetaRand.1.5.mean3days, pSimBatchBetaRand.5.5.meanA, pSimBatchBetaRand.5.5.mean3days, nrow=2, ncol=2, labels="AUTO")
ggsave("pSimBatchBetaRand.1.5.vs.5.5.means.png", width=8, height=8, units="in")


###########  What if we increase the number of replicates?

wormSim10.5.5.m.means<-wormSimBatchBeta3(5,5,100000,1200,24,10,1000)
wormSim10.1.5.m.means<-wormSimBatchBeta3(1,5,100000,1200,24,10,1000)

wormSim25.5.5.m.means<-wormSimBatchBeta3(5,5,100000,1200,24,25,1000)
wormSim25.1.5.m.means<-wormSimBatchBeta3(1,5,100000,1200,24,25,1000)


################# NOT CURRENTLY USING ##############################
pwormSim3.1.5.s.batch1<-subset(wormSim3.1.5.s, batch==1) %>%
  ggplot(aes(x=set, y=logCFU, color=set)) + theme_classic() +
    geom_jitter(shape=16, position=position_jitter(0.05)) +
#    geom_violin(fill=NA) +  
    theme(text=element_text(size=14), 
          axis.title.x = element_blank(), 
          plot.title=element_text(hjust=0.5, size=14)) + 
    labs(title="1", y="log10(CFU/worm)")
pwormSim3.1.5.s.batch1

pwormSim3.1.5.l.batch1<-subset(wormSim3.1.5.l, batch==1) %>%
  ggplot(aes(x=set, y=logCFU, color=set)) + theme_classic() +
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  #    geom_violin(fill=NA) +  
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(), 
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title="1", y="log10(CFU/worm)")
pwormSim3.1.5.l.batch1



####################################################################
####################################################################
#########       FIGURE 8          ##################################
#####     TRYING FOR FALSE NEGATIVES
#####        with real data 
# Pairs 
#(N2 6/21/19 vs dbl-1 6/26/19?) in jointTempBootA 
#(SE day 2 vs. N2 adults 2/13/19) in jointTempBoot

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
#pTempBoot4A + scale_fill_discrete(name = "", breaks=c("dbl1", "N2"), labels = c("dbl-1", "N2"))
# doesn't work idk why

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
ggsave("pjointFalseNegativeBoot.png", width=10, height=7, units="in", dpi=300)

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
