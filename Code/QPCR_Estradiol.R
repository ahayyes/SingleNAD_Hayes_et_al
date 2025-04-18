########## Code to analyse the qPCR data for 17-beta-estradiol ###########

#### Libraries ####
library(tidyverse)
library(ggplot2)
library(MASS) 
library(reshape2)
library(lme4) 
library(emmeans) 
library(DHARMa) 
library(rstatix)
library(MetBrewer)
library(patchwork)
library(ggpubr)

##### Initial Data Read In #####

df <- read.csv("data/Estradiol_qPCR.csv", header=TRUE)

#### Data Wrangling and Averaging ####

#Give each specific sample tube a unique ID
df$Rep <- paste(df$Sample, df$Rep) 

#add in the prevalence data
df <- df %>%
  mutate(d0_previntI1 = df$d0_intI1quantity / df$d0_16Squantity) %>%
  mutate(d7_previntI1 = df$d7_intI1quantity / df$d7_16Squantity)

df <- df %>%
  mutate(d0_prevcintI1 = df$d0_cintI1quantity / df$d0_16Squantity) %>%
  mutate(d7_prevcintI1 = df$d7_cintI1quantity / df$d7_16Squantity)


#IntI1 Dataset
df_intI1 <- df %>% select(Rep, Sample, d0_previntI1, d7_previntI1)

#cintI1 dataset
df_cintI1 <- df %>% select(Rep, Sample, d0_prevcintI1, d7_prevcintI1)


#Average the technical replicates for each dataset
df_avg_intI1 <- aggregate(df_intI1[c(2, 3:4)], list(df_intI1$Rep), mean)

df_avg_cintI1 <- aggregate(df_cintI1[c(2, 3:4)], list(df_cintI1$Rep), mean)


#Melt for use
df_melt_intI1 <- melt(df_avg_intI1, id=c("Group.1", "Sample"))

df_melt_cintI1 <- melt(df_avg_cintI1, id=c("Group.1", "Sample"))

#need to make sure that the Sample.Name column is a factor so it plots correctly, rather than trying to make it into a numeric to fit
df_melt_intI1$Sample <- as.factor(df_melt_intI1$Sample)

df_melt_cintI1$Sample.Name <- as.factor(df_melt_cintI1$Sample)


############# intI1 Stats ####################


#make a linear mixed effect model
stats1 <- lmer(value ~ variable * Sample + (1|Group.1), data=df_melt_intI1)

summary(stats1)
#what are we looking at withthis one?

plot(stats1)
qqnorm(resid(stats1))
qqline(resid(stats1))

simulationOutput <- simulateResiduals(fittedModel = stats1, plot = F)
plot(simulationOutput)

#use logged

stats2 <- lmer(log(value) ~ variable * Sample + (1|Group.1), data=df_melt_intI1)
summary(stats2)

drop1(stats2, test="Chisq")
#the Chisq is 0.3283 so this is not a significant term

stats3 <- lmer(log(value) ~ variable + Sample + (1|Group.1), data=df_melt_intI1)
#stats4 <- update(stats3, ~.-variable:Sample.Name)
summary(stats3)

plot(stats3)
qqnorm(resid(stats3))
qqline(resid(stats3))


simulationOutput <- simulateResiduals(fittedModel = stats3, plot = F)
plot(simulationOutput)

Anova(stats3, test="Chisq", type="III") 

emmeans(stats3, list(pairwise ~ variable+Sample), adjust ="fdr")





#### cintI1 stats ####

stats3 <- lmer(value ~ variable * Sample + (1|Group.1), data=df_melt_cintI1) 

summary(stats3)

drop1(stats3, test="Chisq")
#interaction term is not significant p=0.129

stats3.2 <- lmer(value~variable + Sample + (1|Group.1), data=df_melt_cintI1)

summary(stats3.2)

drop1(stats3.2, test="Chisq")
#Sample is not significantly affecting cintI1 prevalence in the community

simulationOutput <- simulateResiduals(fittedModel = stats3.2, plot = T)
#very badly fitting


stats3.2.2 <- lmer(log(value)~variable + Sample + (1|Group.1), data=df_melt_cintI1)

summary(stats3.2.2)

drop1(stats3.2.2, test="Chisq")
#Sample is not significantly affecting cintI1 prevalence in the community

simulationOutput <- simulateResiduals(fittedModel = stats3.2.2, plot=T)

plot(stats3.2.2)
qqnorm(resid(stats3.2.2))
qqline(resid(stats3.2.2))

#logged data fits the assumptions

ggplot(data = df_melt_sul1, aes(x = resid(stats3.2.2))) +
  geom_histogram(bins = 10, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')


#This dit is kind of okay? qqnorm and histogram looks nice so we will use it

drop1(stats3.2.2, test="Chisq")
#sample conc is not significant in explaining the data p=0.1782






stats3.3 <- lmer(value ~ variable + (1|Group.1), data=df_melt_cintI1)

summary(stats3.3)
#accounts for 8.152e-17 of the variation (the microcosm)

simulationOutput <- simulateResiduals(fittedModel = stats3.3, plot = T)

stats3.4 <- lmer(log(value)~variable + (1|Group.1), data=df_melt_cintI1)

simulationOutput <- simulateResiduals(fittedModel = stats3.4, plot = T)

summary(stats3.4)

Anova(stats3.4, test="Chisq", type = "III")

emmeans(stats3.4, list(pairwise ~ variable*Sample), adjust ="fdr")


#### Plots ############

#Create table for your stats from the stuff in the stats section
stat.test1 <- tibble::tribble(
  ~group1, ~group2, ~variable, ~p.adj,
  "0", "2", "d7_previntI1", 'NS',
  "0", "7", "d7_previntI1", '***',
  "0", "22", "d7_previntI1", '**',
  "0", "67", "d7_previntI1",'**',
  "0", "201", "d7_previntI1",'***',
  "0", "604", "d7_previntI1", '**',
  "0", "1800", "d7_previntI1",'**',
  "0", "5400", "d7_previntI1",'***')

#intI1
ggplot(df_melt_intI1, aes(x=Sample, y=value, fill=variable))+
  geom_boxplot(outlier.shape =  NA, lwd = 0.8)+
  geom_point(position=position_jitterdodge(), alpha=0.7, size = 2, shape = 21, show.legend = FALSE)+
  theme_bw()+
  labs(x="17-β-estradiol concentration (µg/L)", y="intI1 prevalence")+
  ylab(expression(paste(italic("intI1"), "prevalence")))+
  stat_pvalue_manual(stat.test1, y.position = 0.24, label = "p.adj", xmin="group2", xmax=NULL)+
  scale_fill_manual(values = met.brewer("Egypt", 2), name="Day", labels=c("D0", "D7")) +
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size=20))+
  stat_pvalue_manual(stat.test1, label = "p.adj",y.position = 0.24,  xmin="group2", xmax=NULL)

ggsave("Figure 2_estradiolqPCR.tiff", dpi = 300, height = 8, width = 10)

#cintI1
ggplot(df_melt_cintI1, aes(x=Sample.Name, y=value, fill=variable))+
  geom_boxplot(outlier.shape =  NA, lwd = 0.8)+
  geom_point(position=position_jitterdodge(), alpha=0.7, size = 2, shape = 21, show.legend = FALSE)+
  theme_bw()+
  labs(x="17-β-estradiol concentration (µg/L)", y="cintI1 prevalence")+
  ylab(expression(paste(italic("cintI1"), "prevalence")))+
  scale_fill_manual(values = met.brewer("OKeeffe2", 2), name="Day", labels=c("D0", "D7")) +
  scale_x_discrete(labels=c('0', "8", "26", "79", "240", "710", "2100", "6400")) +
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size=20))+
  theme(legend.position = "none")


ggsave("Supplementary Figure 3_EstradiolQPCR.tiff", dpi = 300, height = 8, width = 10)