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
	plot.title=element_text(hjust=0.5, size=14)) + 
	labs(title="Beta-distributed data", y=expression(log[10](CFU/worm)))
pSimBetaRHS
ggsave("FigS1_pSimBetaRHS_1_5.png", width=6, height=4, units="in", dpi=300)

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

