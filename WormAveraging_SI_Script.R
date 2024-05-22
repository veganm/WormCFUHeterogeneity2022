#                Code for Supplementary Information

pacman::p_load(ggplot2, tidyverse, cowplot, mclust, e1071, ggpubr, ggpmisc, readxl, patchwork, sBIC)

# Note that functions are in separate files
# and that this code relies on data objects from the main text script (WormAveraging_MAIN.R)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                   Section 1:
#       Bootstrapping Single-Worm Data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# S. enterica colonization data and bootstrap simulations thereof
# The objects called here are generated in the main text script

#~~~~~~~~~~~
# S. enterica colonization raw data
# Bootstrap many times, over all three runs

input_data<-SaSeCount %>%
  dplyr::filter(Condition=="SE")
SeBootCombinations<-bootOnCountsStats(input_data=input_data, batch_sizes=c(1,5,10,20,50), nboot=1000, 
                                      FoldD=10, correction_constant=20)
glimpse(SeBootCombinations)

SeBootCombinations$Run1ID<-paste("Run", SeBootCombinations$Run1, sep=" ")
SeBootCombinations$Run2ID<-paste("Run", SeBootCombinations$Run2, sep=" ")

pSeBootCombinations_ttest<-SeBootCombinations %>%
  ggplot(aes(x=factor(Batch), y=p_t, color=factor(Batch)))+
  geom_violin(fill=NA)+
  geom_jitter(width=0.05, alpha=0.1)+
  geom_hline(yintercept=0.05)+
  theme_bw()+
  scale_color_viridis_d(option="magma", end=0.85)+
  facet_grid(vars(Run1ID), vars(Run2ID))+
  theme(
    text=element_text(size=xTextSize), 
    plot.title=element_text(hjust=0.5)
  )+
  labs(x="Batch Size", y="t-test p-value", color="Batch",
       title=expression(paste(italic("S. enterica"), " Bootstrap t-tests")))
pSeBootCombinations_ttest

pSeBootCombinations_wilcoxon<-SeBootCombinations %>%
  ggplot(aes(x=factor(Batch), y=p_w, color=factor(Batch)))+
  geom_violin(fill=NA)+
  geom_jitter(width=0.05, alpha=0.1)+
  geom_hline(yintercept=0.05)+
  theme_bw()+
  scale_color_viridis_d(option="magma", end=0.85)+
  facet_grid(vars(Run1ID), vars(Run2ID))+
  theme(
    text=element_text(size=xTextSize), 
    plot.title=element_text(hjust=0.5)
  )+
  labs(x="Batch Size", y="Wilcoxon p-value", color="Batch",
       title=expression(paste(italic("S. enterica"), " Bootstrap Wilcoxon")))
pSeBootCombinations_wilcoxon

#~~~~~~~~~~~~~~~
# how many comparisons are significant for each test?
SeBootCombinations<-SeBootCombinations %>%
  mutate(isSig_t=as.integer(as.logical(p_t<0.05)),
         isSig_w=as.integer(as.logical(p_w<0.05)),
         RunPair=paste(Run1ID, "v", Run2ID, sep=" "))

# plot fraction significant tests
SeBootCombinations_testsummary<-SeBootCombinations %>%
  group_by(RunPair, Batch) %>%
  summarize(countSigT=sum(isSig_t),
            countSigW=sum(isSig_w),
            n=n()) %>%
  mutate(fracSigT=countSigT/n,
         fracSigW=countSigW/n) %>%
  ungroup() %>%
  dplyr::select(RunPair, Batch, fracSigT, fracSigW) %>%
  pivot_longer(cols=starts_with("fracSig"), names_to = "test", 
               names_prefix = "fracSig", values_to = "fracSig")

# remove redundant pairs
pair_list=c("Run 1 v Run 1", "Run 2 v Run 2", "Run 3 v Run 3",
            "Run 1 v Run 2","Run 1 v Run 3","Run 2 v Run 3")
# plot
pSeBootCombinations_testsummary<-SeBootCombinations_testsummary %>%
  dplyr::filter(RunPair %in% pair_list) %>%
  ggplot(aes(x=factor(Batch), y=fracSig, fill=factor(test)))+
  geom_col(position=position_dodge())+
  geom_hline(yintercept=0.05)+
  scale_fill_manual(labels = c("t-test", "Wilcoxon"), values = c("blue", "red")) +
  theme_bw()+
  theme(
    text=element_text(size=xTextSize), 
    plot.title=element_text(hjust=0.5)
  )+
  labs(x="Batch Size", y="Fraction Significant", fill="Test",
       title="Pairwise Tests")+
  facet_wrap(~RunPair)
pSeBootCombinations_testsummary

# how do mean and variance change with batching?
SeBootCombinations_stats1<-SeBootCombinations %>%
  dplyr::select(Batch, Run1, meanCFU1, varCFU1) %>%
  rename(Run=Run1, meanCFU=meanCFU1, varCFU=varCFU1)
SeBootCombinations_stats2<-SeBootCombinations %>%
  dplyr::select(Batch, Run2, meanCFU2, varCFU2) %>%
  rename(Run=Run2, meanCFU=meanCFU2, varCFU=varCFU2)
SeBootCombinations_stats<-rbind(SeBootCombinations_stats1, SeBootCombinations_stats2)
rm(SeBootCombinations_stats1, SeBootCombinations_stats2)
SeBootCombinations_stats <-SeBootCombinations_stats %>%
  mutate(logMeanCFU=log10(meanCFU),
         sdCFU=sqrt(varCFU)/meanCFU,
         RunID=paste("Run", Run, sep=" ")
         )
# plot
pSeBootCombinations_stats_mean<-SeBootCombinations_stats %>%
  ggplot(aes(x=factor(Batch), y=logMeanCFU, color=factor(Batch)))+
  geom_boxplot()+
  #geom_jitter(width=0.1, alpha=0.1)+
  theme_bw()+
  scale_color_viridis_d(option="magma", end=0.85)+
  facet_wrap(~RunID)+
  theme(
    text=element_text(size=xTextSize), 
    plot.title=element_text(hjust=0.5),
    legend.title = element_blank()
  )+
  labs(x="Batch Size", y=expression(log[10]("Mean CFU")), fill="Batch",
       title=expression(paste(italic("S. enterica"), " Bootstrap Means")))
pSeBootCombinations_stats_mean

## how many bootstrap replicates are indistinguishable from Gaussian?
SeBootCombinations_SWtest1<-SeBootCombinations %>%
  dplyr::select(Batch, Run1, p_sw1) %>%
  rename(Run=Run1, p_sw=p_sw1)
SeBootCombinations_SWtest2<-SeBootCombinations %>%
  dplyr::select(Batch, Run2, p_sw2) %>%
  rename(Run=Run2, p_sw=p_sw2)
SeBootCombinations_SWtest<-rbind(SeBootCombinations_SWtest1, SeBootCombinations_SWtest2)
rm(SeBootCombinations_SWtest1, SeBootCombinations_SWtest2)
SeBootCombinations_SWtest<-SeBootCombinations_SWtest %>%
  mutate(isGaussian=as.integer(as.logical(p_sw>0.05)))
glimpse(SeBootCombinations_SWtest)

SeBootCombinations_SWtest_summary<-SeBootCombinations_SWtest %>%
  group_by(Run, Batch) %>%
  summarize(numGaussian=sum(isGaussian),
            n=n()) %>%
  mutate(fracGaussian=numGaussian/n) %>%
  ungroup()

# plot
pSeBootCombinations_SWtest_summary<-SeBootCombinations_SWtest_summary %>%
  ggplot(aes(x=factor(Batch), y=fracGaussian, fill=factor(Run)))+
  geom_col(position=position_dodge())+
  geom_hline(yintercept=0.05)+
  theme_bw()+
  scale_fill_viridis_d(end=0.9)+
  theme(
    text=element_text(size=xTextSize), 
    plot.title=element_text(hjust=0.5)
  )+
  labs(x="Batch Size", y="Fraction Significant", fill="Run",
       title=expression(paste(italic("S. enterica"), " Bootstrap Shapiro-Wilk")))
pSeBootCombinations_SWtest_summary

# assemble
pFigS1AB<-(pSeBootCombinations_ttest/pSeBootCombinations_wilcoxon) +
  plot_annotation(tag_levels = 'A')
#pFigS1AB
pFigS1CDE<-plot_grid(pSeBootCombinations_testsummary,
                       pSeBootCombinations_stats_mean ,
                       pSeBootCombinations_SWtest_summary,
                     ncol=1, labels=c('C', 'D', 'E'))
pFigS1CDE
wrap_elements(pFigS1AB) | wrap_elements(pFigS1CDE)
ggsave("FigS1_SeSingleWormBootStats.png", width=12, height=7, units="in", dpi=300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Same bootstrap with zeros removed
temp0<-dplyr::filter(SaSeCount, Condition=="SE")
idx<-which(temp0$Count==0)
temp0[idx,] # quick look
temp0$Count[idx]<-1
temp0$FinalCount[idx]<-20
temp0$CFUPerWorm[idx]<-20
temp0$logCFU[idx]<-log10(20)
SeBootCombinations0<-bootOnCountsStats(input_data=temp0, batch_sizes=c(1,5,10,20,50), nboot=1000, 
                                      FoldD=10, correction_constant=20)
glimpse(SeBootCombinations0)

SeBootCombinations0$Run1ID<-paste("Run", SeBootCombinations0$Run1, sep=" ")
SeBootCombinations0$Run2ID<-paste("Run", SeBootCombinations0$Run2, sep=" ")

#~~~~~~~~~~~~~~~
# how many comparisons are significant for each test?
SeBootCombinations0<-SeBootCombinations0 %>%
  mutate(isSig_t=as.integer(as.logical(p_t<0.05)),
         isSig_w=as.integer(as.logical(p_w<0.05)),
         RunPair=paste(Run1ID, "v", Run2ID, sep=" "))

# plot fraction significant tests
SeBootCombinations0_testsummary<-SeBootCombinations0 %>%
  group_by(RunPair, Batch) %>%
  summarize(countSigT=sum(isSig_t),
            countSigW=sum(isSig_w),
            n=n()) %>%
  mutate(fracSigT=countSigT/n,
         fracSigW=countSigW/n) %>%
  ungroup() %>%
  dplyr::select(RunPair, Batch, fracSigT, fracSigW) %>%
  pivot_longer(cols=starts_with("fracSig"), names_to = "test", 
               names_prefix = "fracSig", values_to = "fracSig")

# remove redundant pairs
pair_list=c("Run 1 v Run 1", "Run 2 v Run 2", "Run 3 v Run 3",
            "Run 1 v Run 2","Run 1 v Run 3","Run 2 v Run 3")
# plot
pSeBootCombinations0_testsummary<-SeBootCombinations0_testsummary %>%
  dplyr::filter(RunPair %in% pair_list) %>%
  ggplot(aes(x=factor(Batch), y=fracSig, fill=factor(test)))+
  geom_col(position=position_dodge())+
  geom_hline(yintercept=0.05)+
  scale_fill_manual(labels = c("t-test", "Wilcoxon"), values = c("blue", "red")) +
  theme_bw()+
  theme(
    text=element_text(size=xTextSize), 
    plot.title=element_text(hjust=0.5)
  )+
  labs(x="Batch Size", y="Fraction Significant", fill="Test",
       title="No Zeros")+
  facet_wrap(~RunPair)
pSeBootCombinations0_testsummary

# how do mean and variance change with batching?
SeBootCombinations0_stats1<-SeBootCombinations0 %>%
  dplyr::select(Batch, Run1, meanCFU1, varCFU1) %>%
  rename(Run=Run1, meanCFU=meanCFU1, varCFU=varCFU1)
SeBootCombinations0_stats2<-SeBootCombinations0 %>%
  dplyr::select(Batch, Run2, meanCFU2, varCFU2) %>%
  rename(Run=Run2, meanCFU=meanCFU2, varCFU=varCFU2)
SeBootCombinations0_stats<-rbind(SeBootCombinations0_stats1, SeBootCombinations0_stats2)
rm(SeBootCombinations0_stats1, SeBootCombinations0_stats2)
SeBootCombinations0_stats <-SeBootCombinations0_stats %>%
  mutate(logMeanCFU=log10(meanCFU),
         sdCFU=sqrt(varCFU)/meanCFU,
         RunID=paste("Run", Run, sep=" ")
  )
# plot
pSeBootCombinations0_stats_cv<-SeBootCombinations0_stats %>%
  #ggplot(aes(x=factor(Batch), y=logMeanCFU, color=factor(Batch)))+
  mutate(cvCFU=sdCFU/meanCFU) %>%
  ggplot(aes(x=factor(Batch), y=cvCFU, color=factor(Batch)))+
  geom_boxplot()+
  theme_bw()+
  scale_color_viridis_d(option="magma", end=0.85)+
  facet_wrap(~RunID)+
  theme(
    text=element_text(size=xTextSize), 
    plot.title=element_text(hjust=0.5),
    legend.title = element_blank()
  )+
  labs(x="Batch Size", y="sd(CFU)", fill="Batch",
       title="CV, Zeros Removed")
pSeBootCombinations0_stats_cv

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Same bootstrap with zeros enriched
# implemented as higher TOD, call it 10^3
temp3<-dplyr::filter(SaSeCount, Condition=="SE")
idx<-which(temp3$logCFU<3)
temp3$Count[idx]<-0
temp3$D[idx]<-0
temp3$FinalCount[idx]<-0
temp3$CFUPerWorm[idx]<-0
temp3$logCFU[idx]<-0

SeBootCombinations3<-bootOnCountsStats(input_data=temp3, batch_sizes=c(1,5,10,20,50), nboot=1000, 
                                       FoldD=10, correction_constant=20)
glimpse(SeBootCombinations3)

SeBootCombinations3$Run1ID<-paste("Run", SeBootCombinations3$Run1, sep=" ")
SeBootCombinations3$Run2ID<-paste("Run", SeBootCombinations3$Run2, sep=" ")

#~~~~~~~~~~~~~~~
# how many comparisons are significant for each test?
SeBootCombinations3<-SeBootCombinations3 %>%
  mutate(isSig_t=as.integer(as.logical(p_t<0.05)),
         isSig_w=as.integer(as.logical(p_w<0.05)),
         RunPair=paste(Run1ID, "v", Run2ID, sep=" "))

# plot fraction significant tests
SeBootCombinations3_testsummary<-SeBootCombinations3 %>%
  group_by(RunPair, Batch) %>%
  summarize(countSigT=sum(isSig_t),
            countSigW=sum(isSig_w),
            n=n()) %>%
  mutate(fracSigT=countSigT/n,
         fracSigW=countSigW/n) %>%
  ungroup() %>%
  dplyr::select(RunPair, Batch, fracSigT, fracSigW) %>%
  pivot_longer(cols=starts_with("fracSig"), names_to = "test", 
               names_prefix = "fracSig", values_to = "fracSig")

# remove redundant pairs
pair_list=c("Run 1 v Run 1", "Run 2 v Run 2", "Run 3 v Run 3",
            "Run 1 v Run 2","Run 1 v Run 3","Run 2 v Run 3")
# plot
pSeBootCombinations3_testsummary<-SeBootCombinations3_testsummary %>%
  dplyr::filter(RunPair %in% pair_list) %>%
  ggplot(aes(x=factor(Batch), y=fracSig, fill=factor(test)))+
  geom_col(position=position_dodge())+
  geom_hline(yintercept=0.05)+
  scale_fill_manual(labels = c("t-test", "Wilcoxon"), values = c("blue", "red")) +
  theme_bw()+
  theme(
    text=element_text(size=xTextSize), 
    plot.title=element_text(hjust=0.5)
  )+
  labs(x="Batch Size", y="Fraction Significant", fill="Test",
       title="Zero-Enriched")+
  facet_wrap(~RunPair)
pSeBootCombinations3_testsummary

# how do mean and variance change with batching?
SeBootCombinations3_stats1<-SeBootCombinations3 %>%
  dplyr::select(Batch, Run1, meanCFU1, varCFU1) %>%
  rename(Run=Run1, meanCFU=meanCFU1, varCFU=varCFU1)
SeBootCombinations3_stats2<-SeBootCombinations3 %>%
  dplyr::select(Batch, Run2, meanCFU2, varCFU2) %>%
  rename(Run=Run2, meanCFU=meanCFU2, varCFU=varCFU2)
SeBootCombinations3_stats<-rbind(SeBootCombinations3_stats1, SeBootCombinations3_stats2)
rm(SeBootCombinations3_stats1, SeBootCombinations3_stats2)
SeBootCombinations3_stats <-SeBootCombinations3_stats %>%
  mutate(logMeanCFU=log10(meanCFU),
         sdCFU=sqrt(varCFU)/meanCFU,
         RunID=paste("Run", Run, sep=" ")
  )
# plot
pSeBootCombinations3_stats_mean<-SeBootCombinations3_stats %>%
  #mutate(cvCFU=sdCFU/meanCFU) %>%
  ggplot(aes(x=factor(Batch), y=logMeanCFU, color=factor(Batch)))+
  #ggplot(aes(x=factor(Batch), y=sdCFU, color=factor(Batch)))+
  geom_boxplot()+
  theme_bw()+
  scale_color_viridis_d(option="magma", end=0.85)+
  facet_wrap(~RunID)+
  theme(
    text=element_text(size=xTextSize), 
    plot.title=element_text(hjust=0.5),
    legend.title = element_blank()
  )+
  labs(x="Batch Size", y=expression(log[10]("Mean CFU")), fill="Batch",
       title="Means, Zero-Enriched")
pSeBootCombinations3_stats_mean

# merge?
SeBootCombinations_stats$DataSet<-"Data"
SeBootCombinations0_stats$DataSet<-"No Zeros"
SeBootCombinations3_stats$DataSet<-"Zero Enriched"
SeBootCombinationsAll_stats<-rbind(SeBootCombinations_stats,
                                   SeBootCombinations0_stats,
                                   SeBootCombinations3_stats)
SeBootCombinationsAll_stats<-SeBootCombinationsAll_stats %>%
  mutate(cvCFU=sdCFU/meanCFU)

pSeBootCombinationsAll_stats_mean<-SeBootCombinationsAll_stats %>%
  ggplot(aes(x=factor(Batch), y=meanCFU, color=factor(DataSet)))+
  geom_boxplot()+
  theme_bw()+
  scale_y_log10()+
  scale_color_viridis_d(option="cividis", end=0.85)+
  facet_wrap(~RunID)+
  theme(
    text=element_text(size=xTextSize), 
    plot.title=element_text(hjust=0.5),
    legend.title = element_blank()
  )+
  labs(x="Batch Size", y="Mean CFU", fill="",
       title="Means")
pSeBootCombinationsAll_stats_sd<-SeBootCombinationsAll_stats %>%
  ggplot(aes(x=factor(Batch), y=sdCFU, color=factor(DataSet)))+
  geom_boxplot()+
  theme_bw()+
  scale_y_log10()+
  scale_color_viridis_d(option="cividis", end=0.85)+
  facet_wrap(~RunID)+
  theme(
    text=element_text(size=xTextSize), 
    plot.title=element_text(hjust=0.5),
    legend.title = element_blank()
  )+
  labs(x="Batch Size", y="sd(CFU)", fill="",
       title="Standard Deviation")
pSeBootCombinationsAll_stats_CV<-SeBootCombinationsAll_stats %>%
  ggplot(aes(x=factor(Batch), y=cvCFU, color=factor(DataSet)))+
  geom_boxplot()+
  theme_bw()+
  scale_y_log10()+
  scale_color_viridis_d(option="cividis", end=0.85)+
  facet_wrap(~RunID)+
  theme(
    text=element_text(size=xTextSize), 
    plot.title=element_text(hjust=0.5),
    legend.title = element_blank()
  )+
  labs(x="Batch Size", y="cv(CFU)", fill="",
       title="Coefficient of Variation")
#pSeBootCombinationsAll_stats_CV

# combine 
#pSeBootCombinations_stats + pSeBootCombinations0_stats  + pSeBootCombinations3_stats +
pSeBootCombinations_testsummary + 
  pSeBootCombinations0_testsummary + 
  pSeBootCombinations3_testsummary + 
  pSeBootCombinationsAll_stats_mean +
  pSeBootCombinationsAll_stats_sd +
  pSeBootCombinationsAll_stats_CV+
  plot_layout(ncol=2, byrow=FALSE, guides='collect')+
  plot_annotation(tag_levels = 'A') &
  theme(legend.position='bottom')
ggsave("FigS2_SeZeros.png", height=9, width=14, units="in", dpi=300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                       Section 2: 
#            Multi-Modality in CFU/Worm Data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# S. enterica and S. aureus colonization data and simulations thereof
# The objects called here are generated in the main text script
# and the code calls the function "wormbootGMM()"

# Using data for SE and SA with all experimental runs, zeros removed, run IDs as integers 1-3
# Combined data as "Pooled"
# Data are loaded from file on line 92 of main script WormAveraging.R
# The filtered object is created on line 124

mylen<-dim(SaSeCount2)[1]
SaSeCount2$Run<-as.factor(SaSeCount2$Run)
SaSeCount2 <- SaSeCount2 %>%
  dplyr::select(!FinalCount)
SaSeCountPool<-rbind(SaSeCount2, tibble(Condition=SaSeCount2$Condition,
                                        Run="Pooled", # add pooled condition
                                        Batch=1,
                                            Count=SaSeCount2$Count,
                                            D=SaSeCount2$D,
                                            CFUPerWorm=SaSeCount2$CFUPerWorm,
                                            logCFU=SaSeCount2$logCFU
))
pSeCountAll<-SaSeCountPool %>%
  filter(Condition=="SE") %>%
  ggplot(aes(x=Run, y=logCFU, color=Run)) +
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
  labs(title=expression(paste(italic("S. enterica"), " LT2")), y=expression(log[10](CFU/Worm)), color="Run")
pSeCountAll

#ok let's try some gaussian mixture models (code from vignette)
X<-SaSeCountPool$logCFU[SaSeCountPool$Condition=="SE" & SaSeCountPool$Run=="Pooled"]
mean(X)
var(X)
mean(SaSeCountPool$CFU[SaSeCountPool$Condition=="SE"& SaSeCountPool$Run=="Pooled"])
sd(SaSeCountPool$CFU[SaSeCountPool$Condition=="SE"& SaSeCountPool$Run=="Pooled"])

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
ggsave("FigS3_SentericaMclust.png", width=10, height=3, dpi=400, units="in")

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
ggsave("FigS4_pSimSeGMM_run1v3.png", width=12, height=10, units="in", dpi=400)


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
ggsave("FigS5_pSimBetaRHS_1_5.png", width=8, height=4, units="in", dpi=400)
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
#        FIGURE S6  MOMENTS AND DISTANCES OF SIMULATED DATA; 
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
  dplyr::relocate(set) %>%
  dplyr::relocate(logMean, .after=mean) %>%
  dplyr::select(-mean) %>%
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
  dplyr::relocate(set) %>%
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
ggsave("FigS6_pSimBatchBetaRand.1.5.moments.png", width=12, height=10, units="in", dpi=400)      

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################     FIGURE S7     #############################
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
# Assembling Figure S7
plot_grid(pSimBatchBetaRand.1.5.meanA, pSimBatchBetaRand.5.5.meanA,
          pSimBatchBetaRand.1.5.meandist, pSimBatchBetaRand.5.5.meandist, 
          ncol=2, nrow=2, labels="AUTO", align = "h")
ggsave("FigS7_pSimBatchBetaRand.1.5.vs.5.5.means.png", width=8, height=6, units="in", dpi=400)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#               Figure S8
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
ggsave("FigureS8_pwormSim3.1.5.examplerun.png", width=8, height=4, units="in", dpi=300)

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

# Figure S8
plot_grid(pSimBatchBetaRand.5.5.meanA,pSimBatchBetaRand.5.5.mean3days,
          pSimBatchBetaRand.5.5.meandist,pSimBatchBetaRand.5.5.meandist3days, 
          align="h")
ggsave("FigS8_pSimBatchBetaRand.5.5.meancompare3days.png", width=6, height=6, units="in")

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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

wormSimF.1.5.summary %>%
  ggplot(aes(x=batch, y=meanCFU))+
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

#~~~~~~~~~~~~~~~~~~~~~~~
# plot together
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

#~~~~~~~~~~~~~~~~~~~~~~~
# some ANOVAs
wormSimF.1.5.aov10<-wormSimF.1.5 %>%
  dplyr::filter(batch==10)%>%
  aov(logCFU~run+set, data=.)
summary(wormSimF.1.5.aov10)
hist(wormSimF.1.5.aov10$residuals) # quick check - design is even, so should be ok

wormSimF.1.5.aov20<-wormSimF.1.5 %>%
  dplyr::filter(batch==20)%>%
  aov(logCFU~run+set, data=.)
summary(wormSimF.1.5.aov20)
hist(wormSimF.1.5.aov20$residuals) # quick check - design is even, so should be ok

wormSimF.1.5.aov50<-wormSimF.1.5 %>%
  dplyr::filter(batch==50)%>%
  aov(logCFU~run+set, data=.)
summary(wormSimF.1.5.aov50)
hist(wormSimF.1.5.aov50$residuals) # quick check - design is even, so should be ok

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
###  again with a different number of runs to check the biological variances
# SQUARE DESIGN
wormSimF.1.5.12square<-wormSimBatchBetaFactorial(a=1, b=5, reps=12, runs=12, maxCFU=100000)
wormSimF.1.5.12square$run<-as.factor(wormSimF.1.5.12square$run)
wormSimF.1.5.12square.aov10<-wormSimF.1.5.12square %>% dplyr::filter(batch==10) %>% aov(logCFU~run+set, data=.)
wormSimF.1.5.12square.aov50<-wormSimF.1.5.12square %>% dplyr::filter(batch==50) %>% aov(logCFU~run+set, data=.)
wormSimF.1.5.12square.aov5<-wormSimF.1.5.12square %>% dplyr::filter(batch==5) %>% aov(logCFU~run+set, data=.)
wormSimF.1.5.12square.aov20<-wormSimF.1.5.12square %>% dplyr::filter(batch==20) %>% aov(logCFU~run+set, data=.)
wormSimF.1.5.12square.aov1<-wormSimF.1.5.12square %>% dplyr::filter(batch==1) %>% aov(logCFU~run+set, data=.)
summary(wormSimF.1.5.12square.aov1)
summary(wormSimF.1.5.12square.aov5)
summary(wormSimF.1.5.12square.aov10)
summary(wormSimF.1.5.12square.aov20)
summary(wormSimF.1.5.12square.aov50)

wormSimF.1.5.12square.aov1.summary<-summary(wormSimF.1.5.12square.aov1)
wormSimF.1.5.12square.aov1.summary[[1]]
wormSimF.1.5.12square.aov1.summary[[1]][,2] # SSEs
wormSimF.1.5.12square.aov1.summary[[1]][,5] # P-vals

# interesting. MSEs are very similar in batched
# see if this holds
# variability is higher
wormSimF.1.5.5square<-wormSimBatchBetaFactorial(a=1, b=5, reps=5, runs=5, maxCFU=100000)
wormSimF.1.5.5square$run<-as.factor(wormSimF.1.5.5square$run)
wormSimF.1.5.5square.aov10<-wormSimF.1.5.5square %>% dplyr::filter(batch==10) %>% aov(logCFU~run+set, data=.)
wormSimF.1.5.5square.aov50<-wormSimF.1.5.5square %>% dplyr::filter(batch==50) %>% aov(logCFU~run+set, data=.)
wormSimF.1.5.5square.aov5<-wormSimF.1.5.5square %>% dplyr::filter(batch==5) %>% aov(logCFU~run+set, data=.)
wormSimF.1.5.5square.aov20<-wormSimF.1.5.5square %>% dplyr::filter(batch==20) %>% aov(logCFU~run+set, data=.)
wormSimF.1.5.5square.aov1<-wormSimF.1.5.5square %>% dplyr::filter(batch==1) %>% aov(logCFU~run+set, data=.)
summary(wormSimF.1.5.5square.aov1)
summary(wormSimF.1.5.5square.aov5)
summary(wormSimF.1.5.5square.aov10)
summary(wormSimF.1.5.5square.aov20)
summary(wormSimF.1.5.5square.aov50)

#~~~~~~~~~~~~~~~~~~~~~
# First, non-square designs

factor_dim<-c(3, 5, 8, 10, 12, 15, 20) # squares
batch_sizes<-c(1, 5, 10, 20, 50)
iter<-25 # number of times to iterate the simulation

# Create an empty temporary vector of defined size to hold the data
temp<-vector("list", length=length(batch_sizes)*length(factor_dim)* iter)

idx<-1
for (k in seq_len(iter)){
  for (i in seq_along(factor_dim)){ # loop along dimension of square factorial design
    for (m in seq_along(factor_dim)){ # again
      # this code allows n_reps and n_runs to differ
      temp.beta<-wormSimBatchBetaFactorial(a=1, b=5, reps=factor_dim[i], runs=factor_dim[m], maxCFU=100000)
      temp.beta$run<-as.factor(temp.beta$run)
      for (j in seq_along(batch_sizes)){ # loop over batch sizes
        temp.beta.aov<-temp.beta %>% dplyr::filter(batch==batch_sizes[j]) %>% aov(logCFU~run+set, data=.)
        temp.beta.aov.summary<-summary(temp.beta.aov)
        tempA<- temp.beta %>% dplyr::filter(batch==batch_sizes[j] & set=="A")
        tempB<- temp.beta %>% dplyr::filter(batch==batch_sizes[j] & set=="B")
        wtest<-wilcox.test(tempA$CFU, tempB$CFU)
        temp[[idx]]<-tibble(
          Iteration=k,
          n_reps=factor_dim[i],
          n_runs=factor_dim[m],
          batch=batch_sizes[j],
          Component=c("Run", "Set"),
          MSE=temp.beta.aov.summary[[1]][1:2,3],
          p_aov=temp.beta.aov.summary[[1]][1:2,5],
          p_w=wtest$p.value
        )
        idx<-idx+1
      } # close wrap over batch sizes
    } # again
  } # close wrap over factorial dimension
} # close wrap over iterations

# plot out
wormSimF.1.5.RunRep.summary<-dplyr::bind_rows(temp)

wormSimF.1.5.RunRep.summary_wide <- wormSimF.1.5.RunRep.summary %>%
  #dplyr::filter(batch>1)%>%
  pivot_wider(names_from=Component, values_from = c(MSE, p_aov))
wormSimF.1.5.RunRep.summary_wide %>%
  ggplot(aes(x=MSE_Run, y=MSE_Set))+
  geom_point(aes(color=factor(batch)))+
  #scale_x_log10()+
  #scale_y_log10()+
  stat_poly_line(formula=y~x+0)+
  stat_poly_eq(use_label("eq"), formula = y ~ x + 0)

wormSimF.1.5.RunRep.summary_wide<-wormSimF.1.5.RunRep.summary_wide%>%
  mutate(MSE_ratio=MSE_Set/MSE_Run)

wormSimF.1.5.RunRep.summary_wide %>%
  ggplot(aes(y=MSE_ratio, x=factor(n_reps), color=factor(batch)))+
  geom_boxplot()+
  scale_y_log10()+
  geom_hline(yintercept = 1)+
  facet_grid(~n_runs)

# p values?
wormSimF.1.5.RunRep.summary %>%
  dplyr::filter(Component=="Set") %>%
  ggplot(aes(x=factor(batch), y=p_aov, color=factor(n_reps)))+
  geom_boxplot()+
  geom_hline(yintercept = 0.05)+
  scale_y_log10()+
  facet_grid(~n_runs)

wormSimF.1.5.RunRep.summary %>%
  dplyr::filter(Component=="Set") %>%
  ggplot(aes(x=factor(batch), y=p_w, color=factor(n_reps)))+
  geom_boxplot()+
  geom_hline(yintercept = 0.05)+
  scale_y_log10()+
  facet_grid(~n_runs)


##############
# now the squares
# loop it
factor_dim<-c(3, 5, 8, 10, 12, 15, 20) # squares
batch_sizes<-c(1, 5, 10, 20, 50)
iter<-25 # number of times to iterate the simulation

# Create an empty temporary vector of defined size to hold the data
temp<-vector("list", length=length(batch_sizes)*length(factor_dim)* iter)

idx<-1
for (k in seq_len(iter)){
  for (i in seq_along(factor_dim)){ # loop along dimension of square factorial design
    temp.square<-wormSimBatchBetaFactorial(a=1, b=5, reps=factor_dim[i], runs=factor_dim[i], maxCFU=100000)
    temp.square$run<-as.factor(temp.square$run)
    for (j in seq_along(batch_sizes)){ # loop over batch sizes
      temp.square.aov<-temp.square %>% dplyr::filter(batch==batch_sizes[j]) %>% aov(logCFU~run+set, data=.)
      temp.square.aov.summary<-summary(temp.square.aov)
      tempA<- temp.square %>% dplyr::filter(batch==batch_sizes[j] & set=="A")
      tempB<- temp.square %>% dplyr::filter(batch==batch_sizes[j] & set=="B")
      wtest<-wilcox.test(tempA$CFU, tempB$CFU)
      temp[[idx]]<-tibble(
        Iteration=k,
        n_reps=factor_dim[i],
        n_runs=factor_dim[i],
        batch=batch_sizes[j],
        Component=c("Run", "Set"),
        MSE=temp.square.aov.summary[[1]][1:2,3],
        p_aov=temp.square.aov.summary[[1]][1:2,5],
        p_w=wtest$p.value
      )
      idx<-idx+1
    } # close wrap over batch sizes
  } # close wrap over factorial dimension
} # close wrap over iterations

# plot out
wormSimF.1.5.square.summary<-dplyr::bind_rows(temp)
glimpse(wormSimF.1.5.square.summary)

wormSimF.1.5.square.summary_wide <- wormSimF.1.5.square.summary %>%
  #dplyr::filter(batch>1)%>%
  pivot_wider(names_from=Component, values_from = c(MSE, p_aov))
wormSimF.1.5.square.summary_wide %>%
  ggplot(aes(x=MSE_Run, y=MSE_Set))+
  geom_point(aes(size=n_reps, color=factor(batch)))+
  scale_x_log10()+
  scale_y_log10()+
  stat_poly_line(formula=y~x+0)+
  stat_poly_eq(use_label("eq"), formula = y ~ x + 0)

wormSimF.1.5.square.summary_wide<-wormSimF.1.5.square.summary_wide%>%
  mutate(MSE_ratio=MSE_Set/MSE_Run)

wormSimF.1.5.square.summary_wide %>%
  ggplot(aes(y=MSE_ratio, x=factor(n_reps), color=factor(batch)))+
  geom_boxplot()+
  scale_y_log10()+
  geom_hline(yintercept = 1)

# p values?
wormSimF.1.5.square.summary %>%
  dplyr::filter(Component=="Set") %>%
  ggplot(aes(x=factor(batch), y=p_aov, color=factor(n_reps)))+
  geom_boxplot()+
  geom_hline(yintercept = 0.05)+
  scale_y_log10()

wormSimF.1.5.square.summary %>%
  dplyr::filter(Component=="Set") %>%
  ggplot(aes(x=factor(batch), y=p_w, color=factor(n_reps)))+
  geom_boxplot()+
  geom_hline(yintercept = 0.05)+
  scale_y_log10()

#~~~~~~~~~~~~
# Now results from ANOVA of Beta(1,5) simulated data 
# vs data from Beta(2,5)

# Create an empty temporary vector of defined size to hold the data
temp<-vector("list", length=length(batch_sizes)*length(factor_dim)* iter)

idx<-1
for (k in seq_len(iter)){
  for (i in seq_along(factor_dim)){ # loop along dimension of square factorial design
    temp.square.A<-wormSimBatchBetaFactorial(a=1, b=5, reps=factor_dim[i], runs=factor_dim[i], maxCFU=100000)
    temp.square.B<-wormSimBatchBetaFactorial(a=2, b=5, reps=factor_dim[i], runs=factor_dim[i], maxCFU=100000)
    temp.square.A <- temp.square.A %>%
      dplyr::filter(set=="A")
    temp.square.B <- temp.square.B %>%
      dplyr::filter(set=="B")
    temp.square<-rbind(temp.square.A, temp.square.B)
    temp.square$run<-as.factor(temp.square$run)
    for (j in seq_along(batch_sizes)){ # loop over batch sizes
      temp.square.aov<-temp.square %>% dplyr::filter(batch==batch_sizes[j]) %>% aov(logCFU~run+set, data=.)
      temp.square.aov.summary<-summary(temp.square.aov)
      temp.square.A1 <- temp.square.A %>%
        dplyr::filter(batch==batch_sizes[j])
      temp.square.B1 <- temp.square.B %>%
        dplyr::filter(batch==batch_sizes[j])
      wtest<-wilcox.test(temp.square.A1$CFU, temp.square.B1$CFU)
      temp[[idx]]<-tibble(
        Iteration=k,
        n_reps=factor_dim[i],
        n_runs=factor_dim[i],
        batch=batch_sizes[j],
        Component=c("Run", "Set"),
        MSE=temp.square.aov.summary[[1]][1:2,3],
        p_aov=temp.square.aov.summary[[1]][1:2,5],
        p_w=wtest$p.value
      )
      idx<-idx+1
    } # close wrap over batch sizes
  } # close wrap over factorial dimension
} # close wrap over iterations

# plot out
wormSimF.TwoBeta.square.summary<-dplyr::bind_rows(temp)

# p values
wormSimF.TwoBeta.square.summary %>%
  dplyr::filter(Component=="Set") %>%
  ggplot(aes(x=factor(batch), y=p_aov, color=factor(n_reps)))+
  geom_boxplot()+
  geom_hline(yintercept = 0.05)+
  scale_y_log10()

wormSimF.TwoBeta.square.summary %>%
  dplyr::filter(Component=="Set") %>%
  ggplot(aes(x=factor(batch), y=p_w, color=factor(n_reps)))+
  geom_boxplot()+
  geom_hline(yintercept = 0.05)+
  scale_y_log10()

# MSE
wormSimF.TwoBeta.square.summary_wide <- wormSimF.TwoBeta.square.summary %>%
  #dplyr::filter(batch>1)%>%
  pivot_wider(names_from=Component, values_from = c(MSE, p_aov))
wormSimF.TwoBeta.square.summary_wide<-wormSimF.TwoBeta.square.summary_wide%>%
  mutate(MSE_ratio=MSE_Set/MSE_Run)
glimpse(wormSimF.TwoBeta.square.summary_wide)

wormSimF.TwoBeta.square.summary_wide %>%
  ggplot(aes(y=MSE_ratio, x=factor(n_reps), color=factor(batch)))+
  geom_boxplot()+
  scale_y_log10()+
  geom_hline(yintercept = 1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# combine simulated data sets
wormSimF.1.5.RunRep.summary$Condition<-"Beta(1,5)"
wormSimF.1.5.square.summary$Condition<-"Beta(1,5) Square"
wormSimF.TwoBeta.square.summary$Condition<-"B(1,5) vs B(2,5)"

wormSimF.Beta.Factorial<-rbind(wormSimF.1.5.RunRep.summary,
                               wormSimF.1.5.square.summary,
                               wormSimF.TwoBeta.square.summary)
glimpse(wormSimF.Beta.Factorial)

wormSimF.1.5.RunRep.summary_wide$Condition<-"Beta(1,5)"
wormSimF.1.5.square.summary_wide$Condition<-"Beta(1,5) Square"
wormSimF.TwoBeta.square.summary_wide$Condition<-"B(1,5) vs B(2,5)"

wormSimF.Beta.Factorial_wide<-rbind(wormSimF.1.5.RunRep.summary_wide,
                                    wormSimF.1.5.square.summary_wide,
                                    wormSimF.TwoBeta.square.summary_wide)
glimpse(wormSimF.Beta.Factorial_wide)

# Generate summary table for p values for set to set comparisons
wormSimF.Beta.Factorial.pvals<-wormSimF.Beta.Factorial %>%
  dplyr::filter(Component=="Run") %>%
  mutate(isSig_aov = as.integer(p_aov<=0.05),
         isSig_w=as.integer(p_w<=0.05)) %>%
  group_by(Condition, n_reps, n_runs, batch) %>%
  summarize(psig_aov=sum(isSig_aov),
            psig_w=sum(isSig_w),
            n=n()
            )
wormSimF.Beta.Factorial.pvals <- wormSimF.Beta.Factorial.pvals %>%
  mutate(fsig_aov=psig_aov/n,
         fsig_w=psig_w/n)
view(wormSimF.Beta.Factorial.pvals)

# Fraction of significant p-values for run in ANOVA
wormSimF.Beta.Factorial.pvals %>%
  ggplot(aes(x=factor(batch), y=fsig_aov, color=factor(Condition)))+
  geom_boxplot()+
  geom_hline(yintercept = 0.05)+
  labs(x="Batch Size", y="Fraction Significant, Run", title="ANOVA",
       color="")

# Fraction of significant p-values for run in Wilcoxon
wormSimF.Beta.Factorial.pvals %>%
  ggplot(aes(x=factor(batch), y=fsig_w, color=factor(Condition)))+
  geom_boxplot()+
  geom_hline(yintercept = 0.05)+
  labs(x="Batch Size", y="Fraction Significant, Run", title="Wilcoxon",
       color="")

# MSE ratio
wormSimF.Beta.Factorial_wide %>%
  ggplot(aes(y=MSE_ratio, x=factor(batch), color=factor(n_reps)))+
  geom_boxplot()+
  scale_y_log10()+
  geom_hline(yintercept = 1)+
  facet_wrap(~Condition)+
  labs(x="Batch Size", y=expression(MSE[Run]*"/"*MSE[Replicate]), title="MSE Ratio",
       color="Replicates")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Simulate data from named distributions for comparison
# let's take the SA single-worm data set (three runs with D0 plated)

SaCount<- SaSeCount %>%
  dplyr::filter(Condition=="SA")

SaCount6<-SaCount %>%
  dplyr::filter(Run==6)

# simulate data from stats of individual worm real data
SaCount6_distributions<-sim_means_named(input_data=SaCount6, n_reps=12, n_runs=12, correction_constant = 20)

# Diagnostic plots
SaCount6_distributions %>%
  ggplot(aes(x=factor(Batch), y=logFinalCount, color=factor(Run)))+
  geom_jitter(width=0.1)+
  facet_wrap(~dist_name)

SaCount6_distributions_summary<-SaCount6_distributions %>%
  group_by(dist_name, Batch, Run) %>%
  summarize(meanTotal=mean(FinalCount, na.rm=TRUE),
            sdTotal=sd(FinalCount, na.rm=TRUE),
            n=n()
  )

SaCount6_distributions_summary<-SaCount6_distributions_summary %>%
  mutate(cvTotal=sdTotal/meanTotal)
view(SaCount6_distributions_summary)

SaCount6_distributions_summary %>%
  ggplot(aes(x=factor(Batch), y=cvTotal, color=factor(dist_name)))+
  #geom_jitter(width=1)+
  geom_boxplot()+
  scale_color_viridis_d(option="turbo")

SaCount6_distributions_summary %>%
  ggplot(aes(x=factor(Batch), y=meanTotal, color=factor(dist_name)))+
  geom_jitter(width=0.1)+
  scale_y_log10()+
  scale_color_viridis_d()+
  facet_wrap(~dist_name)

SaCount6_distributions_summary %>%
  ggplot(aes(x=factor(Batch), y=meanTotal, color=factor(dist_name)))+
  geom_boxplot()+
  scale_y_log10()+
  scale_color_viridis_d(option="turbo")
