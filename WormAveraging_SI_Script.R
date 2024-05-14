#                Code for Supplementary Information

pacman::p_load(ggplot2, tidyverse, cowplot, mclust, e1071, ggpubr, readxl, patchwork, sBIC)

# Note that functions are in separate files
# and that this code relies on data objects from the main text script (WormAveraging_MAIN.R)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                       Section 1: 
#            Multi-Modality in CFU/Worm Data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# S. enterica and S. aureus colonization data and simulations thereof
# The objects called here are generated in the main text script
# and the code calls the function "wormbootGMM()"
#library(tidyverse)
#library(cowplot)
#library(mclust, quietly=TRUE)
#library(sBIC)

# Using data for SE and SA with all experimental runs, zeros removed, run IDs as integers 1-3
# Combined data as "Pooled"
# Data are loaded from file on line 92 of main script WormAveraging.R
# The filtered object is created on line 124

mylen<-dim(SaSeCount2)[1]
SaSeCount2$Rep<-as.factor(SaSeCount2$Rep)
SaSeCountPool<-rbind(SaSeCount2, data.frame(Condition=SaSeCount2$Condition,
                                            Rep=rep("Pooled", mylen), # add pooled condition
                                            Count=SaSeCount2$Count,
                                            D=SaSeCount2$D,
                                            CFU=SaSeCount2$CFU,
                                            logCFU=SaSeCount2$logCFU
))
pSeCountAll<-SaSeCountPool %>%
  filter(Condition=="SE") %>%
  ggplot(aes(x=Rep, y=logCFU, color=Rep)) +
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + 
  theme_classic() + 
  scale_color_viridis_d(end=0.9)+
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=14),
        ) + 
  # facet_wrap(vars(Condition), scales="free_x") +
  labs(title=expression(paste(italic("S. enterica"), " LT2")), y=expression(log[10](CFU/Worm)), color="Replicate")
#ggsave("pSECountAll.png", width=4, height=3, units="in", dpi=300)
pSeCountAll

#ok let's try some gaussian mixture models (code from vignette)
X<-SaSeCountPool$logCFU[SaSeCountPool$Condition=="SE" & SaSeCountPool$Rep=="Pooled"]
mean(X)
var(X)
mean(SaSeCountPool$CFU[SaSeCountPool$Condition=="SE"& SaSeCountPool$Rep=="Pooled"])
sd(SaSeCountPool$CFU[SaSeCountPool$Condition=="SE"& SaSeCountPool$Rep=="Pooled"])

fit<- mclustBIC(X)
fit
#Top 3 models based on the BIC criterion: 
#  E,2       V,2       E,1 
#-174.8787 -177.0372 -180.1553 

plot(fit)
# build a tibble for ggplotting
fitBIC.df<-tibble(Components=rep(1:9),
                  E=fit[,1],
                  V=fit[,2],
               )

# plot
pSePooledBICMclust<-fitBIC.df %>%
  pivot_longer(
    cols=E:V,
    names_to="Var",
    values_to="BIC"
  ) %>%
  ggplot(aes(x=Components, y=BIC, linetype=factor(Var), pch=factor(Var)))+
  geom_line()+
  scale_x_continuous(breaks=seq(1, 9, by=1))+
  theme_classic() + 
  theme(text=element_text(size=14), 
        #axis.title.x = element_blank(), 
        #axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=14),
        legend.title = element_blank()) + 
  labs(title=expression(paste(italic("S. enterica"), " LT2, Pooled")), x="# Components", y="BIC")
pSePooledBICMclust


# Here we are still fitting to the combined Salmonella colonization data.
# Fitting a model with two components, allowing variances to be unequal ("V") or equal ("E")
fit2 <- Mclust(X, G=2, model="E")
summary(fit2)
plot(fit2, what="density", main="G2", xlab="logCFU")
rug(X)
fit2$parameters
#(mean1, var1)=(3.225159, 0.1877358) and (mean2, var2)=(4.662897, 0.1877358), containing 39% and 61% of the mass respectively

# let's put these plots together with the Se data plot
names(fit2)

# easier to ggplot if we call to density? This will be the same model as above
fit2d<-densityMclust(X, G=2)
summary(fit2d)
names(fit2d)
glimpse(fit2d)
rug(X)

# ok let's try ggplotting
as.numeric(fit2d$data)
fit2d_toplot<-tibble(logCFU=as.numeric(fit2d$data), density=as.numeric(fit2d$density))
glimpse(fit2d_toplot)
attach(fit2d_toplot)
newdata<-fit2d_toplot[order(logCFU),]
glimpse(newdata)
detach(fit2d_toplot)

pSeAllDensityMclust<-newdata %>%
  ggplot(aes(x=logCFU, y=density))+
  geom_line()+
  geom_rug(sides="b")+
  theme_classic() + 
  theme(text=element_text(size=14), 
        #axis.title.x = element_blank(), 
        #axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=14),
        legend.title = element_blank()) + 
  labs(title=expression(paste(italic("S. enterica"), " LT2, Pooled")), x=expression(log[10](CFU/Worm)), y="Density")
pSeAllDensityMclust

plot_grid(pSeCountAll, pSeAllDensityMclust, pSePooledBICMclust, ncol=3, labels="AUTO", align="h")
ggsave("FigS1_SentericaMclust.png", width=10, height=3, dpi=400, units="in")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# We can also look at each replicate individually
X1<-SaSeCountPool$logCFU[SaSeCountPool$Condition=="SE" & SaSeCountPool$Rep=="1"]
fit1<-mclustBIC(X1)
fit1 # Best support for G=2, model E
fit1 <- Mclust(X1, G=2, model="E")
summary(fit1)
plot(fit1, what="density", main="G2", xlab="logCFU")
rug(X1)
fit1$parameters
#(mean1, var1)=(3.3866, 0.05795) and (mean2, var2)=(4.74495, 0.05795), containing 34.6% and 65.4% of the mass respectively

X2<-SaSeCountPool$logCFU[SaSeCountPool$Condition=="SE" & SaSeCountPool$Rep=="2"]
fit2<-mclustBIC(X2)
fit2 #best support for E1, then E2
fit2 <- Mclust(X2, G=2, model="V")
fit2$parameters
#(mean1, var1)=(2.774017, 0.09134093) and (mean2, var2)=(4.319905, 0.32749870), containing 42.6% and 57.3% of the mass respectively

X3<-SaSeCountPool$logCFU[SaSeCountPool$Condition=="SE" & SaSeCountPool$Rep=="3"]
fit3<-mclustBIC(X3)
fit3 #Best support for E2
fit3 <- Mclust(X3, G=2, model="E")
summary(fit3)
plot(fit3, what="density", main="G2", xlab="logCFU")
rug(X3)
fit3$parameters
#(mean1, var1)=(3.7839, 0.06167) and (mean2, var2)=(4.9194, 0.06167), containing 42.1% and 57.9% of the mass respectively

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

# fit 1
#(mean1, var1)=(3.3866, 0.05795) and (mean2, var2)=(4.74495, 0.05795), containing 34.6% and 65.4% of the mass respectively
# fit 3
#(mean1, var1)=(3.7839, 0.06167) and (mean2, var2)=(4.9194, 0.06167), containing 42.1% and 57.9% of the mass respectively
# fit pooled
#(mean1, var1)=(3.225159, 0.1877358) and (mean2, var2)=(4.662897, 0.1877358), containing 39% and 61% of the mass respectively

SeBootGMM<-wormbootGMM(100, 50, 3.225159, 0.1877358, 4.662897, 0.1877358, 0.61, 0)
SeBootGMM$batch<-as.factor(SeBootGMM$batch)
SeBootGMM$set<-as.factor(SeBootGMM$set)

SeBootGMM.1.2<-wormbootGMM(100, 50, 3.225159, 0.1877358, 4.662897, 0.1877358, 0.61, 0.1)
SeBootGMM.1.2$batch<-as.factor(SeBootGMM.1.2$batch)
SeBootGMM.1.2$set<-as.factor(SeBootGMM.1.2$set)

SeBootGMM.1.2.05<-wormbootGMM(100, 50, 3.225159, 0.1877358, 4.662897, 0.1877358, 0.61, 0.05)
SeBootGMM.1.2.05$batch<-as.factor(SeBootGMM.1.2.05$batch)
SeBootGMM.1.2.05$set<-as.factor(SeBootGMM.1.2.05$set)


# Plot out each set of simulations and store as plot objects, to combine in figure
pSeBootGMM<-SeBootGMM %>%
  ggplot(aes(x=set, y=logCount, color=set))+
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  geom_violin(fill=NA) + 
  ylim(-0.1,6) +
  theme_classic() + 
  scale_color_viridis_d(option="magma", begin=0.3, end=0.85)+
  theme(text=element_text(size=16), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
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
  scale_color_viridis_d(option="magma", begin=0.3, end=0.85)+
  theme(text=element_text(size=16), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=16),
        legend.position="none"
  ) + 
  labs(title="Pooled data parameterization, 10% difference in mode proportions", y=expression(log[10](CFU/worm)))+
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
  scale_color_viridis_d(option="magma", begin=0.3, end=0.85)+
  theme(text=element_text(size=16), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title=element_text(hjust=0.5, size=16),
        legend.position="none"
  ) + 
  labs(title="Pooled data parameterization, 5% difference in mode proportions", y=expression(log[10](CFU/worm)))+
  facet_wrap(~batch, nrow=1)+
  stat_compare_means(label.y=0.2)+
  stat_compare_means(method="t.test", label.y = 0.9)
pSeBootGMM.1.2.05

plot_grid(pSeBootGMM, pSeBootGMM.1.2, pSeBootGMM.1.2.05, ncol=1, labels="AUTO")
ggsave("FigS2_pSimSeGMM_run1v3.png", width=12, height=10, units="in", dpi=400)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                     Section 2:
# Simulating "CFU" data: effects of variation and skew on false-negative rates 
#
# Calls function simBetaCFU()
# to generate simulated data based on the beta distribution
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
simBatchBetaR5<-wormSimBatchBetaCompare(0.5, 2.5, 0.5, 2.5, 1000, 25, 100000)
simBatchBetaR10<-wormSimBatchBetaCompare(1, 5, 1, 5, 1000, 25, 100000)

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
#ggsave("pSimBetaLHS_5_1.png", width=6, height=4, units="in", dpi=300)

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
# without much loss of generality let's draw parameter shifts from U(-0.1, 0.1) each run

wormSimBetaRand.1.5<-wormSimBatchBetaCompare(1,5,1,5,1000,24,100000, 0.1)
wormSimBetaRand.p5.2p5<-wormSimBatchBetaCompare(0.5,2.5,0.5,2.5,1000,24,100000, 0.1)
wormSimBetaRand.5.1<-wormSimBatchBetaCompare(5,1,5,1,1000,24,100000, 0.1)
wormSimBetaRand.2p5.p5<-wormSimBatchBetaCompare(2.5,0.5,2.5,0.5,1000,24,100000, 0.1)
wormSimBetaRand.5.5<-wormSimBatchBetaCompare(5,5,5,5,1000,24,100000, 0.1)
wormSimBetaRand.2p5.2p5<-wormSimBatchBetaCompare(2.5,2.5,2.5,2.5,1000,24,100000, 0.1)
wormSimBetaRand.1.1<-wormSimBatchBetaCompare(1,1,1,1,1000,24,100000, 0.1)

# calculate CV; include the distances
wormSimBetaRand.1.5$meandist<-abs(wormSimBetaRand.1.5$meanA-wormSimBetaRand.1.5$meanB)/(wormSimBetaRand.1.5$meanA+wormSimBetaRand.1.5$meanB)
wormSimBetaRand.1.5$cvA<-sqrt(wormSimBetaRand.1.5$varA)/wormSimBetaRand.1.5$meanA
wormSimBetaRand.1.5$cvB<-sqrt(wormSimBetaRand.1.5$varB)/wormSimBetaRand.1.5$meanB
wormSimBetaRand.1.5$cvdist<-abs(wormSimBetaRand.1.5$cvA-wormSimBetaRand.1.5$cvB)
wormSimBetaRand.1.5$skewdist<-abs(wormSimBetaRand.1.5$skewA-wormSimBetaRand.1.5$skewB)
wormSimBetaRand.1.5$kurtdist<-abs(wormSimBetaRand.1.5$kurtA-wormSimBetaRand.1.5$kurtB)

# Flip the data set around so it can get factor-gridded
wormSimBetaRand.1.5_A<-as_tibble(wormSimBetaRand.1.5) %>%
  select(batch, meanA, cvA, skewA, kurtA) %>%
  rename(mean=meanA,
         cv=cvA,
         skew=skewA,
         kurt=kurtA)
wormSimBetaRand.1.5_A$set<-"A"

wormSimBetaRand.1.5_B<-as_tibble(wormSimBetaRand.1.5) %>%
  select(batch, meanB, cvB, skewB, kurtB) %>%
  rename(mean=meanB,
         cv=cvB,
         skew=skewB,
         kurt=kurtB)
wormSimBetaRand.1.5_B$set<-"B"

wormSimBetaRand.1.5_dist<-as_tibble(wormSimBetaRand.1.5) %>%
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
  scale_color_viridis_d(begin=0.3, end=0.8) +  
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
  geom_violin(fill=NA) + 
  theme_classic() + 
  scale_color_viridis_d(begin=0.3, end=0.8) +
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(),  legend.position="none",
        plot.title=element_text(hjust=0.5, size=14)) + 
  #  labs(title=paste("\u03b2","(1,5) A", sep=""), y="Mean")
  facet_grid(rows=vars(moment), cols=vars(set), scales="free_y") 

plot_grid(pSimBatchBetaRand.1.5.AB, pSimBatchBetaRand.1.5.dist, ncol=2, rel_widths = c(2,1), labels="AUTO")
ggsave("FigS4_pSimBatchBetaRand.1.5.moments.png", width=12, height=10, units="in", dpi=400)      

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################     FIGURE S5     #############################
#### add data from symmetric distribution for comparison
#note these runs use the flat U(-0.1,0.1) parameter error
wormSimBetaRand.5.5<-wormSimBatchBetaCompare(5,5,5,5,1000,24,100000, 0.1)
wormSimBetaRand.5.5$cvA<-sqrt(wormSimBetaRand.5.5$varA)/wormSimBetaRand.1.5$meanA
wormSimBetaRand.5.5$cvB<-sqrt(wormSimBetaRand.5.5$varB)/wormSimBetaRand.1.5$meanB
wormSimBetaRand.5.5$cvdist<-abs(wormSimBetaRand.5.5$cvA-wormSimBetaRand.1.5$cvB)
wormSimBetaRand.5.5$meandist<-abs(wormSimBetaRand.5.5$meanA-wormSimBetaRand.5.5$meanB)/(wormSimBetaRand.5.5$meanA+wormSimBetaRand.5.5$meanB)

pSimBatchBetaRand.1.5.meanA<-wormSimBetaRand.1.5 %>%
  ggplot(aes(x=batch, y=log10(meanA), color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
 # ylim(3.8, 5)+
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  scale_color_viridis_d(begin=0.3, end=0.8) +
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(),  
        legend.position="none",
        plot.title=element_text(hjust=0.5)) + 
  labs(title=paste("\u03b2","(1,5) A", sep=""), y="Mean")
#pSimBatchBetaRand.1.5.meanA

pSimBatchBetaRand.1.5.meandist<-wormSimBetaRand.1.5 %>%
  ggplot(aes(x=batch, y=meandist, color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  scale_color_viridis_d(begin=0.3, end=0.8) +
  ylim(-0.01,0.45)+
  theme(text=element_text(size=14), 
        #axis.title.x = element_blank(), 
        legend.position="none",
        plot.title=element_text(hjust=0.5)) + 
  labs(title="Distance", y="Dist(mean)", x="Batch Size")

pSimBatchBetaRand.5.5.meanA<-wormSimBetaRand.5.5 %>%
  ggplot(aes(x=batch, y=log10(meanA), color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
 # ylim(3.8, 5)+
  scale_color_viridis_d(begin=0.3, end=0.8) +
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(),  
        legend.position="none",
        plot.title=element_text(hjust=0.5)) + 
  labs(title=paste("\u03b2","(5,5) A", sep=""), y="Mean")

pSimBatchBetaRand.5.5.meandist<-wormSimBetaRand.5.5 %>%
  ggplot(aes(x=batch, y=meandist, color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  ylim(-0.01,0.45)+
  scale_color_viridis_d(begin=0.3, end=0.8) +
  theme(text=element_text(size=14), 
        #axis.title.x = element_blank(), 
        legend.position="none",
        plot.title=element_text(hjust=0.5)) + 
  #labs(title="Distance", y="Dist(mean)")
  labs(title=paste("\u03b2","(5,5) Distance"), y="Dist(Mean)", x="Batch Size")

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

#a<-1
#b<-5
#(a*b)/((a+b+1)*(a+b)^2)
#2*(b-a)*sqrt(a+b+1)/(sqrt(a*b)*(a+b+2))

#nWorms<-100
#maxSamples<-12
#maxCFU<-100000

# calls function
# wormSimBatchBeta3(a, b, maxCFU, nWorms, maxSamples, runs, reps=3, sameParams=FALSE, returnMEANS=TRUE)

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

# here we use enough worms to allow 24 batches even at 50 worms
wormSim3.5.5.m.means<-wormSimBatchBeta3(5,5,100000,1200,24,1000)
wormSim3.1.5.m.means<-wormSimBatchBeta3(1,5,100000,1200,24,1000)
wormSim3.5.5.m.means$batch<-as.factor(wormSim3.5.5.m.means$batch)
wormSim3.1.5.m.means$batch<-as.factor(wormSim3.1.5.m.means$batch)

# simulated data for plotting
wormSim3.1.5.m<-wormSimBatchBeta3(1,5,100000,1200,24,1000, returnMEANS=FALSE)
wormSim3.1.5.m$batch<-as.factor(wormSim3.1.5.m$batch)
wormSim3.1.5.m$day<-as.factor(wormSim3.1.5.m$day)
wormSim3.1.5.m$replicate<-as.factor(wormSim3.1.5.m$replicate)

wormSim3.1.5.m %>%
  ggplot(aes(x=batch, y=logCFU, color=replicate)) + 
  geom_violin(fill=NA) + 
  geom_point(shape=16, position=position_jitterdodge(0.2)) +
  ylim(1,5) + theme_classic() + 
  scale_color_viridis_d(option="mako", end=0.8)+
  theme(text=element_text(size=14), 
        #axis.title.x = element_blank(), 
        plot.title=element_text(hjust=0.5),
        ) + 
  facet_wrap(~set)+
  labs(title="Simulated data, Beta(1,5)", 
       x="Batch size",
       y=expression(log[10](Bacteria)),
       color="Replicate")
ggsave("FigureS7_pwormSim3.1.5.examplerun.png", width=8, height=4, units="in", dpi=300)

## plot differences in convergence of means for one-day vs. 3-day data
wormSim3.5.5.m.means$meandist<-abs(wormSim3.5.5.m.means$meanA-wormSim3.5.5.m.means$meanB)/(wormSim3.5.5.m.means$meanA+wormSim3.5.5.m.means$meanB)

pSimBatchBetaRand.5.5.mean3days<-wormSim3.5.5.m.means %>%
  ggplot(aes(x=batch, y=log10(meanA), color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  ylim(4.1,5)+
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  scale_color_viridis_d(begin=0.3, end=0.8) +
  theme(text=element_text(size=14), 
        axis.title.x = element_blank(), legend.position="none",
        plot.title=element_text(hjust=0.5)) + 
  labs(title=paste("\u03b2","(5,5)x3"), y="Mean")
#pSimBatchBetaRand.5.5.mean3days

pSimBatchBetaRand.5.5.meandist3days<-wormSim3.5.5.m.means %>%
  ggplot(aes(x=batch, y=meandist, color=batch))+
  geom_violin(fill=NA) + theme_classic() + 
  ylim(-0.01,0.2)+
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  scale_color_viridis_d(begin=0.3, end=0.8) +
  theme(text=element_text(size=14), 
        #axis.title.x = element_blank(), 
        legend.position="none",
        plot.title=element_text(hjust=0.5)) + 
  labs(title=paste("\u03b2","(5,5)x3 Distance"), y="Dist(Mean)", x="Batch Size")
#pSimBatchBetaRand.5.5.meandist3days

# Figure S6
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


###########  What if we instead take the populations of means in regressions?
wormSimF.1.5<-wormSimBatchBetaFactorial(a=1, b=5, reps=24, maxCFU=100000)
glimpse(wormSimF.1.5)
wormSimF.1.5$run<-as.factor(wormSimF.1.5$run)
wormSimF.1.5 %>%
 ggplot(aes(x=factor(batch), y=logCFU, color=run)) + 
  geom_violin(fill=NA) + 
  geom_point(shape=16, position=position_jitterdodge(0.2)) +
  ylim(1,5) + theme_classic() + 
  scale_color_viridis_d(option="mako", end=0.8)+
  theme(text=element_text(size=14), 
        #axis.title.x = element_blank(), 
        plot.title=element_text(hjust=0.5),
  ) + 
  facet_wrap(~set)+
  labs(title="Simulated data, Beta(1,5)", 
       x="Batch size",
       y=expression(log[10](Bacteria)),
       color="Replicate")

# put through GLM
wormsimF.1.5.glm.single<-wormSimF.1.5 %>%
  dplyr::filter(batch==1) %>%
  glm(logCFU~run*set, family=Gamma, data=.)
summary(wormsimF.1.5.glm.single)

# checks
plot(density(resid(wormsimF.1.5.glm.single, type='deviance')))
scatter.smooth(1:length(rstandard(wormsimF.1.5.glm.single, type='deviance')), rstandard(wormsimF.1.5.glm.single, type='deviance'), col='gray')
scatter.smooth(predict(wormsimF.1.5.glm.single, type='response'), #residuals vs fitted
               rstandard(wormsimF.1.5.glm.single, type='deviance'), col='gray')

wormsimF.1.5.glm.batch5<-wormSimF.1.5 %>%
  dplyr::filter(batch==5) %>%
  glm(logCFU~run*set, family=Gamma, data=.)
summary(wormsimF.1.5.glm.batch5)

wormsimF.1.5.glm.batch10<-wormSimF.1.5 %>%
  dplyr::filter(batch==10) %>%
  glm(logCFU~run*set, family=Gamma, data=.)
summary(wormsimF.1.5.glm.batch10)

wormsimF.1.5.glm.batch20<-wormSimF.1.5 %>%
  dplyr::filter(batch==20) %>%
  glm(logCFU~run*set, family=Gamma, data=.)
summary(wormsimF.1.5.glm.batch20)

wormsimF.1.5.glm.batch50<-wormSimF.1.5 %>%
  dplyr::filter(batch==50) %>%
  glm(logCFU~run*set, family=Gamma, data=.)
summary(wormsimF.1.5.glm.batch50)

# and the symmetric data
wormSimF.5.5<-wormSimBatchBetaFactorial(a=5, b=5, reps=24, maxCFU=100000)
glimpse(wormSimF.5.5)
wormSimF.5.5$run<-as.factor(wormSimF.5.5$run)
wormSimF.5.5 %>%
  ggplot(aes(x=factor(batch), y=logCFU, color=run)) + 
  geom_violin(fill=NA) + 
  geom_point(shape=16, position=position_jitterdodge(0.2)) +
  ylim(1,5) + theme_classic() + 
  scale_color_viridis_d(option="mako", end=0.8)+
  theme(text=element_text(size=14), 
        #axis.title.x = element_blank(), 
        plot.title=element_text(hjust=0.5),
  ) + 
  facet_wrap(~set)+
  labs(title="Simulated data, Beta(5,5)", 
       x="Batch size",
       y=expression(log[10](Bacteria)),
       color="Replicate")


# Compare glm and Wilcoxon false positives
# with three runs as above
# here we use enough worms to allow 24 batches even at 50 worms
wormSim3.5.5.m.means<-wormSimBatchBeta3(5,5,100000,1200,24,1000)
wormSim3.5.5.m.means$batch<-as.factor(wormSim3.5.5.m.means$batch)
# Wilcoxon p-values
#[1] 0.096 0.232 0.372 0.461 0.592
# glm p-values for runB term
#[1] 0.108 0.282 0.368 0.500 0.641

wormSim3.1.5.m.means<-wormSimBatchBeta3(1,5,100000,1200,24,1000)
wormSim3.1.5.m.means$batch<-as.factor(wormSim3.1.5.m.means$batch)
# Wilcoxon
#0.061 0.120 0.199 0.310 0.449
# glm
# 0.044 0.114 0.230 0.335 0.492

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# summarize the simulated factorial data
wormSimF.1.5.summary<-wormSimF.1.5 %>%
  group_by(set, batch, run) %>%
  summarize(meanCFU=mean(CFU, na.rm=TRUE),
            varCFU=var(CFU, na.rm=TRUE),
            q25CFU=quantile(CFU,probs=0.25, na.rm=TRUE),
            medianCFU=quantile(CFU,probs=0.5, na.rm=TRUE),
            q75CFU=quantile(CFU,probs=0.75, na.rm=TRUE),
            n=n()
            )
view(wormSimF.1.5.summary)

wormSimF.1.5.summary %>%
  ggplot(aes(x=batch, y=varCFU))+
  geom_point()+
  scale_y_log10()

wormSimF.5.5.summary<-wormSimF.5.5 %>%
  group_by(set, batch, run) %>%
  summarize(meanCFU=mean(CFU, na.rm=TRUE),
            varCFU=var(CFU, na.rm=TRUE),
            q25CFU=quantile(CFU,probs=0.25, na.rm=TRUE),
            medianCFU=quantile(CFU,probs=0.5, na.rm=TRUE),
            q75CFU=quantile(CFU,probs=0.75, na.rm=TRUE),
            n=n()
  )
view(wormSimF.5.5.summary)
wormSimF.5.5.summary %>%
  ggplot(aes(x=batch, y=varCFU))+
  geom_point()+
  scale_y_log10()

wormSimF.5.5.summary$Distribution<-"Beta(5,5)"
wormSimF.1.5.summary$Distribution<-"Beta(1,5)"
wormSimF.summary<-rbind(wormSimF.1.5.summary, wormSimF.5.5.summary)
wormSimF.summary %>%
  ggplot(aes(x=batch, y=varCFU, color=factor(Distribution)))+
  geom_jitter(width=1)+
  scale_y_log10()

# messing with distributions of effort
tempF.1.5<-wormSimBatchBetaFactorial(a=1, b=5, runs=5, reps=10, maxCFU=100000)
tempF.1.5$run<-as.factor(tempF.1.5$run)
tempF.1.5.summary<-tempF.1.5 %>%
  group_by(set, batch, run) %>%
  summarize(meanCFU=mean(CFU, na.rm=TRUE),
            varCFU=var(CFU, na.rm=TRUE),
            q25CFU=quantile(CFU,probs=0.25, na.rm=TRUE),
            medianCFU=quantile(CFU,probs=0.5, na.rm=TRUE),
            q75CFU=quantile(CFU,probs=0.75, na.rm=TRUE),
            n=n()
  )
tempF.1.5.summary$Distribution<-"Beta(1,5)"

tempF.5.5<-wormSimBatchBetaFactorial(a=5, b=5, runs=5, reps=10, maxCFU=100000)
tempF.5.5$run<-as.factor(tempF.5.5$run)
tempF.5.5.summary<-tempF.5.5 %>%
  group_by(set, batch, run) %>%
  summarize(meanCFU=mean(CFU, na.rm=TRUE),
            varCFU=var(CFU, na.rm=TRUE),
            q25CFU=quantile(CFU,probs=0.25, na.rm=TRUE),
            medianCFU=quantile(CFU,probs=0.5, na.rm=TRUE),
            q75CFU=quantile(CFU,probs=0.75, na.rm=TRUE),
            n=n()
  )
tempF.5.5.summary$Distribution<-"Beta(5,5)"
tempF.summary<-rbind(tempF.1.5.summary, tempF.5.5.summary)
tempF.summary %>%
  ggplot(aes(x=batch, y=varCFU, color=factor(Distribution)))+
  geom_jitter(width=1)+
  scale_y_log10()

# Plot the average CFU within each experiment and set
tempF.summary %>%
  ggplot(aes(x=batch, y=meanCFU, color=factor(Distribution)))+
  geom_jitter(width=1)+
  scale_y_log10()

# expand. Simulate data from named distributions for comparison
