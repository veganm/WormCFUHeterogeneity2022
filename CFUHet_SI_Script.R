library(tidyverse)
library(cowplot)

#                Code for Supplementary Information
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                       Section 1: 
#            Multi-Modality in CFU/Worm Data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# S. enterica and S. aureus bleach experiment data and simulations thereof
# The objects called here are generated in the main text script
# and the code calls the function "wormbootGMM()" from the script of the same name
library(mclust, quietly=TRUE)
library(sBIC)

# Using data for SE and SA with all three experimental runs, zeros removed, dates as integers 1-3
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                    FIGURE S2
# Using the all-SE 2-mode GMM, we can generate simulated data
# Note that we are allowing the weight in modes to be stochastic (binomial) 

temp<-wormbootGMM(100, 50, 2.985569, 0.4568056, 4.742767, 0.1239928, 0.58, 0)

SeBootGMM<-wormbootGMM(100, 50, 2.985569, 0.4568056, 4.742767, 0.1239928, 0.58, 0)
SeBootGMM$set<-as.factor(SeBootGMM$set)
SeBootGMM$batch<-as.factor(SeBootGMM$batch)

SeBootGMM.1.2<-wormbootGMM(100, 50, 2.985569, 0.4568, 4.742767, 0.12399, 0.57, 0.15)
SeBootGMM.1.2$set<-as.factor(SeBootGMM.1.2$set)
SeBootGMM.1.2$batch<-as.factor(SeBootGMM.1.2$batch)

SeBootGMM.1.2.05<-wormbootGMM(100, 50, 2.985569, 0.4568, 4.742767, 0.12399, 0.57, 0.05)

# change the data set labels to be less confusing
SeBootGMM$set<-as.numeric(SeBootGMM$set)
SeBootGMM<-as.tibble(SeBootGMM) %>%
  mutate(set=replace(set, set=="A", 1)) %>%
  mutate(set=replace(set, set=="B", 2))
glimpse(SeBootGMM)
SeBootGMM$set<-as.factor(SeBootGMM$set)

SeBootGMM.1.2$set<-as.numeric(SeBootGMM.1.2$set)
SeBootGMM.1.2<-as.tibble(SeBootGMM.1.2) %>%
  mutate(set=replace(set, set=="A", 1)) %>%
  mutate(set=replace(set, set=="B", 2))
glimpse(SeBootGMM.1.2)
SeBootGMM.1.2$set<-as.factor(SeBootGMM.1.2$set)

SeBootGMM.1.2.05<-as.tibble(SeBootGMM.1.2.05) %>%
  mutate(set=replace(set, set=="A", 1)) %>%
  mutate(set=replace(set, set=="B", 2))
glimpse(SeBootGMM.1.2)

# Plot out each set of simulations and store as plot objects, to combine in figure
pSeBootGMM<-SeBootGMM %>%
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
  labs(title="Pooled data parameterization", y=expression(log[10](CFU/worm)))+
  facet_wrap(~batch, nrow=1)+
  stat_compare_means(label.y=0.2)+
  stat_compare_means(method="t.test", label.y = 0.9)
pSeBootGMM

pSeBootGMM.1.2<-SeBootGMM.1.2 %>%
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
  labs(title="Pooled data parameterization, 15% difference in mode proportions", y=expression(log[10](CFU/worm)))+
  facet_wrap(~batch, nrow=1)+
  stat_compare_means(label.y=0.2)+
  stat_compare_means(method="t.test", label.y = 0.9)
pSeBootGMM.1.2

pSeBootGMM.1.2.05<-SeBootGMM.1.2.05 %>%
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
  labs(title="Pooled data parameterization, 5% difference in mode proportions", y=expression(log[10](CFU/worm)))+
  facet_wrap(~batch, nrow=1)+
  stat_compare_means(label.y=0.2)+
  stat_compare_means(method="t.test", label.y = 0.9)
pSeBootGMM.1.2.05

plot_grid(pSeBootGMM, pSeBootGMM.1.2, pSeBootGMM.1.2.05, ncol=1, labels="AUTO")
ggsave("FigS2_pSimSeGMM_day1v2.png", width=12, height=10, units="in", dpi=400)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                     Section 2:
# Simulating "CFU" data: effects of variation and skew on false-negative rates 
#
# Calls function simBetaCFU() from the script of the same name
# to generate fake data based on the beta distribution
# with a max of 5 logs CFU/individual (default)
# Start with right hand skew:
a<-0.1
b<-0.5
nreps<-50

simBetaRHS<-simBetaCFU(a,b,nreps)
glimpse(simBetaRHS)

pSimBetaRHS<-simBetaRHS %>%
	ggplot(aes(x=Factor, y=LogData, color=Factor))+
	geom_jitter(shape=16, position=position_jitter(0.05)) +
	geom_violin(fill=NA) + 
	theme_classic() + 
	theme(text=element_text(size=16), 
	axis.title.x = element_blank(), 
	plot.title=element_text(hjust=0.5, size=14),
	legend.position = "none") + 
	labs(title="Beta-distributed data", y=expression(log[10](CFU/worm)))
pSimBetaRHS + scale_x_discrete(labels=c("1"="B(0.1, 0.5)",
                                        "2.5"="B(0.25, 1.25)",
                                        "5"="B(0.5, 2.5)",
                                        "10"="B(1, 5)",
                                        "20"="B(2, 10)",
                                        "50"="B(5, 25)"))
ggsave("FigS3_pSimBetaRHS_1_5.png", width=8, height=4, units="in", dpi=400)
#ggsave("FigS3_pSimBetaRHS_1_5.pdf", width=8, height=4, units="in")

# note that the mean of these is all 100000*0.16667 = 16667 expected
# and all the actual arithmetic means are close to this, which is nice
simBetaRHS %>%
  group_by(Factor) %>%
  summarize(means=mean(Data),
            medians=median(Data)) %>%
  ungroup()


# comparisons
SimBetaRHS.ttests<-compare_means(Data~Factor, data=simBetaRHS, method="t.test")
SimBetaRHS.wtests<-compare_means(Data~Factor, data=simBetaRHS)
write.table(SimBetaRHS.ttests, "SimBetaRHS_ttests.txt", sep="\t")
write.table(SimBetaRHS.wtests, "SimBetaRHS_wtests.txt", sep="\t")

# we will focus on X5 Beta(0.5, 2.5) and X10 Beta(1,5)
# let's create batch data using wormSimBoot-based functions
simBatchBetaR5<-wormSimBatchBeta(0.5, 2.5, 1000, 25, 100000)
simBatchBetaR10<-wormSimBatchBeta(1, 5, 1000, 25, 100000)
simBatchBetaR5v10<-wormSimBatchBetaCompare(0.5, 2.5, 1, 5, 1000, 25, 100000)
simBatchBetaR5v20<-wormSimBatchBetaCompare(0.5, 2.5, 2, 10, 1000, 25, 100000)


#######
# now left hand skew
a<-0.5
b<-0.1
simBetaLHS<-simBetaCFU(a,b,nreps)
pSimBetaLHS<- simBetaLHS %>%
	ggplot(aes(x=Factor, y=LogData, color=Factor))+
	geom_jitter(shape=16, position=position_jitter(0.05)) +
	geom_violin(fill=NA) + 
	theme_classic() + 
	theme(text=element_text(size=16), 
	axis.title.x = element_blank(), 
	plot.title=element_text(hjust=0.5, size=16)) + 
	labs(title="Beta-distributed data", y="log10(CFU/worm)")
pSimBetaLHS
ggsave("pSimBetaLHS_5_1.png", width=6, height=4, units="in", dpi=300)

# mean should be 100000*(5/6) = 83,333
# and it's close across the board
simBetaLHS %>%
  group_by(Factor) %>%
  summarize(means=mean(Data),
            medians=median(Data))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#        FIGURE S4  MOMENTS AND DISTANCES OF SIMULATED DATA; 
#                 RUN-TO-RUN VARIATION     
#
# ok back to beta starting at B(1,5) why not
# without much loss of generality let's draw factors from a U(0,1) each run
# function is in wormSimBatchBetaRand.R

wormSimBetaRand.1.5<-wormSimBatchBetaRand(1,5,1000,24,100000)
wormSimBetaRand.p5.2p5<-wormSimBatchBetaRand(0.5,2.5,1000,24,100000)
wormSimBetaRand.5.1<-wormSimBatchBetaRand(5,1,1000,24,100000)
wormSimBetaRand.2p5.p5<-wormSimBatchBetaRand(2.5,0.5,1000,24,100000)
wormSimBetaRand.2p5.2p5<-wormSimBatchBetaRand(2.5,2.5,1000,24,100000)
wormSimBetaRand.1.1<-wormSimBatchBetaRand(1,1,1000,24,100000)
#mgenbeta(3, 1, 5, 1)
#mgenbeta(3, 0.5, 2.5, 1)

# calculate CV; include the distances
wormSimBetaRand.1.5$meandist<-abs(wormSimBetaRand.1.5$meanA-wormSimBetaRand.1.5$meanB)/(wormSimBetaRand.1.5$meanA+wormSimBetaRand.1.5$meanB)
wormSimBetaRand.1.5$cvA<-sqrt(wormSimBetaRand.1.5$varA)/wormSimBetaRand.1.5$meanA
wormSimBetaRand.1.5$cvB<-sqrt(wormSimBetaRand.1.5$varB)/wormSimBetaRand.1.5$meanB
wormSimBetaRand.1.5$cvdist<-abs(wormSimBetaRand.1.5$cvA-wormSimBetaRand.1.5$cvB)
wormSimBetaRand.1.5$skewdist<-abs(wormSimBetaRand.1.5$skewA-wormSimBetaRand.1.5$skewB)
wormSimBetaRand.1.5$kurtdist<-abs(wormSimBetaRand.1.5$kurtA-wormSimBetaRand.1.5$kurtB)

# Flip the data set around so it can get factor-gridded
wormSimBetaRand.1.5_A<-as.tibble(wormSimBetaRand.1.5) %>%
  select(batch, meanA, cvA, skewA, kurtA) %>%
  rename(mean=meanA,
         cv=cvA,
         skew=skewA,
         kurt=kurtA)
wormSimBetaRand.1.5_A$set<-"A"

wormSimBetaRand.1.5_B<-as.tibble(wormSimBetaRand.1.5) %>%
  select(batch, meanB, cvB, skewB, kurtB) %>%
  rename(mean=meanB,
         cv=cvB,
         skew=skewB,
         kurt=kurtB)
wormSimBetaRand.1.5_B$set<-"B"

wormSimBetaRand.1.5_dist<-as.tibble(wormSimBetaRand.1.5) %>%
  select(batch, meandist, cvdist, skewdist, kurtdist) %>%
  rename(mean=meandist,
         cv=cvdist,
         skew=skewdist,
         kurt=kurtdist)
wormSimBetaRand.1.5_dist$set<-"Distance"

# and combine
wormSimBetaRand.1.5_long<-rbind(wormSimBetaRand.1.5_A, wormSimBetaRand.1.5_B)
wormSimBetaRand.1.5_long$logMean<-log10(wormSimBetaRand.1.5_long$mean)

#plot out
pSimBatchBetaRand.1.5.AB<-wormSimBetaRand.1.5_long %>%
  relocate(set) %>%
  relocate(logMean, .after=mean) %>%
  select(-mean) %>%
  rename(average=logMean) %>%
  pivot_longer(average:kurt, names_to="moment", values_to="value") %>%
  ggplot(aes(x=batch, y=value, color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(),  legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  #  labs(title=paste("\u03b2","(1,5) A", sep=""), y="Mean")
  facet_grid(rows=vars(moment), cols=vars(set), scales="free_y")

pSimBatchBetaRand.1.5.dist<-wormSimBetaRand.1.5_dist %>%
  relocate(set) %>%
  rename(average=mean) %>%
  pivot_longer(average:kurt, names_to="moment", values_to="value") %>%
  ggplot(aes(x=batch, y=value, color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(),  legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  #  labs(title=paste("\u03b2","(1,5) A", sep=""), y="Mean")
  facet_grid(rows=vars(moment), cols=vars(set), scales="free_y") 

plot_grid(pSimBatchBetaRand.1.5.AB, pSimBatchBetaRand.1.5.dist, ncol=2, rel_widths = c(2,1), labels="AUTO")
ggsave("Fig3_pSimBatchBetaRand.1.5.moments.png", width=12, height=10, units="in", dpi=400)      

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################     FIGURE S5     #############################
#### add data from symmetric distribution for comparison
#note these runs use the flat U(-0.1,0.1) parameter error
wormSimBetaRand.5.5<-wormSimBatchBetaRand(5,5,1000,24,100000)
wormSimBetaRand.5.5$cvA<-sqrt(wormSimBetaRand.5.5$varA)/wormSimBetaRand.1.5$meanA
wormSimBetaRand.5.5$cvB<-sqrt(wormSimBetaRand.5.5$varB)/wormSimBetaRand.1.5$meanB
wormSimBetaRand.5.5$cvdist<-abs(wormSimBetaRand.5.5$cvA-wormSimBetaRand.1.5$cvB)
wormSimBetaRand.5.5$meandist<-abs(wormSimBetaRand.5.5$meanA-wormSimBetaRand.5.5$meanB)/(wormSimBetaRand.5.5$meanA+wormSimBetaRand.5.5$meanB)

pSimBatchBetaRand.1.5.meanA<-wormSimBetaRand.1.5 %>%
  ggplot(aes(x=batch, y=log10(meanA), color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  ylim(4.1, 5)+
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(),  
        legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title=paste("\u03b2","(1,5) A", sep=""), y="Mean")

pSimBatchBetaRand.1.5.meandist<-wormSimBetaRand.1.5 %>%
  ggplot(aes(x=batch, y=meandist, color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  ylim(-0.01,0.45)+
  theme(text=element_text(size=14), 
        #axis.title.x = element_blank(), 
        legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title="Distance", y="Dist(mean)")

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
  #labs(title="Distance", y="Dist(mean)")
  labs(title=paste("\u03b2","(5,5) Distance"), y="Dist(Mean)")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Assembling Figure S5
plot_grid(pSimBatchBetaRand.1.5.meanA, pSimBatchBetaRand.5.5.meanA,
          pSimBatchBetaRand.1.5.meandist, pSimBatchBetaRand.5.5.meandist, 
          ncol=2, nrow=2, labels="AUTO", align = "h")
ggsave("FigS5_pSimBatchBetaRand.1.5.vs.5.5.means.png", width=8, height=6, units="in", dpi=400)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#               Figure S6
# Now more realistic - we allow triplicates, 
# and cap the number of both total worms and digests we will perform
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
  labs(title=paste("\u03b2","(5,5)x3"), y="Mean")
pSimBatchBetaRand.5.5.mean3days

pSimBatchBetaRand.5.5.meandist3days<-wormSim3.5.5.m.means %>%
  ggplot(aes(x=batch, y=meandist, color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  ylim(-0.01,0.2)+
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(), legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  labs(title=paste("\u03b2","(5,5)x3 Distance"), y="Dist(Mean)")
pSimBatchBetaRand.5.5.meandist3days
plot_grid(pSimBatchBetaRand.5.5.meanA,pSimBatchBetaRand.5.5.mean3days,
          pSimBatchBetaRand.5.5.meandist,pSimBatchBetaRand.5.5.meandist3days, 
          align="h")
ggsave("FigS6_pSimBatchBetaRand.5.5.meancompare3days.png", width=6, height=6, units="in")

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

