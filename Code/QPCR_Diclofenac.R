################## QPCR Diclofenac #######################

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


##### Initial Data Read In #####
#read in data

df <- read.csv("data/Diclofenac_qPCR.csv", header=TRUE)

#Give each specific sample tube a unique ID
df$Rep <- paste(df$Sample.Name, df$Rep) 

#add in the prevalence data
df <- df %>%
  mutate(d0_previntI1 = df$quantity_d0intI1 / df$quantity_d016S) %>%
  mutate(d7_previntI1 = df$quantity_d7_intI1 / df$quantity_d7_16s)

df <- df %>%
  mutate(d0_prevcintI1 = df$quantity_d0cintI1 / df$quantity_d016S_2) %>%
  mutate(d7_prevcintI1 = df$quantity_d7cintI1 / df$quantity_d716S_2)

#IntI1 Dataset
df_intI1 <- df %>% select(Rep, Sample.Name, d0_previntI1, d7_previntI1)

#cintI1 dataset
df_cintI1 <- df %>% select(Rep, Sample.Name, d0_prevcintI1, d7_prevcintI1)

#Average the technical replicates for each dataset
df_avg_intI1 <- aggregate(df_intI1[c(2, 3:4)], list(df_intI1$Rep), mean)

df_avg_cintI1 <- aggregate(df_cintI1[c(2, 3:4)], list(df_cintI1$Rep), mean)


#Melt for use
df_melt_intI1 <- melt(df_avg_intI1, id=c("Group.1", "Sample.Name"))

df_melt_cintI1 <- melt(df_avg_cintI1, id=c("Group.1", "Sample.Name"))

#need to make sure that the Sample.Name column is a factor so it plots correctly, rather than trying to make it into a numeric to fit
df_melt_intI1$Sample.Name <- as.factor(df_melt_intI1$Sample.Name)

df_melt_cintI1$Sample.Name <- as.factor(df_melt_cintI1$Sample.Name)



#### intI1 Stats ####
stats1 <- lmer(value ~ variable * Sample.Name + (1|Group.1), data=df_melt_intI1)

summary(stats1)

drop1(stats1, test="Chisq")
#the Chisq value is 0.03271 so this is  significant and should be kept in the model

#Post hoc
Anova(stats1, test="Chisq", type="III") 
#This shows that the intI1 prevalence changes with time p=0.03571
#But not treatment p=0.92035   

simulationOutput <- simulateResiduals(fittedModel = stats1, plot = T)
#The model fits well within the parameters


######## cinti1 stats #########

df_cintI1_stats <- na.omit(df_melt_cintI1)

m3 <- lmer(value ~ variable * Sample.Name + (1|Group.1), data=df_cintI1_stats)

plot(m3)
qqnorm(resid(m3))
qqline(resid(m3))

simulationOutput <- simulateResiduals(fittedModel = m3, plot = T)

ggplot(data = df_cintI1_stats, aes(x = resid(m3))) +
  geom_histogram(bins = 10, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

#Not great fitting on the qqplots
#Slightly positive skew on the histogram

m3.2 <- lmer(log(value) ~ variable * Sample.Name + (1|Group.1), data=df_cintI1_stats)

plot(m3.2)
qqnorm(resid(m3.2))
qqline(resid(m3.2))

simulationOutput <- simulateResiduals(fittedModel = m3.2, plot = T)

ggplot(data = df_cintI1_stats, aes(x = resid(m3.2))) +
  geom_histogram(bins = 10, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

#Great fit for qqplot and dharma
#More of a weird fit for the histogram
#Will use this data

drop1(m3.2, test="Chisq")

#Interaction term is significant p=0.0008076  

Anova(m3.2, test="Chisq", type="III") 


summary(m3.2)

emmeans(m3.2, list(pairwise ~ variable*Sample.Name), adjust ="fdr")












######## Plots #############
p_intI1 <- ggplot(df_melt_intI1, aes(x=Sample.Name, y=value, fill=variable))+
  geom_boxplot(outlier.shape =  NA, lwd = 0.8)+
  geom_point(position=position_jitterdodge(), alpha=0.7, size = 2, shape = 21, show.legend = FALSE)+
  theme_bw()+
  labs(x="Diclofenac concentration (µg/L)", y="intI1 prevalence")+
  ylab(expression(paste(italic("intI1"), "prevalence")))+
  scale_x_discrete(labels=c('0', "8", "26", "79", "240", "710", "2100", "6400"))+
  scale_fill_manual(values = met.brewer("VanGogh3", 2), name="Day", labels=c("D0", "D7")) +
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size=20))+
  theme(legend.position = "none")


p_cintI1 <- ggplot(df_melt_cintI1, aes(x=Sample.Name, y=value, fill=variable))+
  geom_boxplot(outlier.shape =  NA, lwd = 0.8)+
  geom_point(position=position_jitterdodge(), alpha=0.7, size = 2, shape = 21, show.legend = FALSE)+
  theme_bw()+
  labs(x="Diclofenac concentration (µg/L)", y="cintI1 prevalence")+
  ylab(expression(paste(italic("cintI1"), "prevalence")))+
  scale_fill_manual(values = met.brewer("VanGogh3", 2), name="Day", labels=c("D0", "D7")) +
  scale_x_discrete(labels=c('0', "8", "26", "79", "240", "710", "2100", "6400")) +
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size=20))+
  theme(legend.position = "none")


p_intI1 / p_cintI1 + plot_annotation(tag_levels = 'A')

ggsave("figures/Supplementary Figure 1.tiff", dpi = 300, height = 10, width = 8)




