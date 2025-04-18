############################## Metformin qPCR ###################################



#### Load in Libraries ####
library(tidyverse)
library(ggplot2)
library(MASS) 
library(reshape2)
library(lme4) #makes the model
library(emmeans) #post hoc tests
library(DHARMa) #test the model is happy
library(rstatix)
library(MetBrewer)
library(patchwork)


#### Load in data ####
df <- read.csv("data/Metformin_qPCR.csv", header=TRUE)

#### Wrangle Data ####

#Give each specific sample tube a unique ID
df$Rep <- paste(df$Sample, df$Rep) 

#add in the prevalence data
df <- df %>%
  mutate(d0_previntI1 = df$quantity_d0intI1 / df$quantity_d016S) %>%
  mutate(d7_previntI1 = df$quantity_d7intI1 / df$quantity_d716S)

df <- df %>%
  mutate(d0_prevcintI1 = df$quantity_d0_cintI1 / df$quantity_d016S) %>%
  mutate(d7_prevcintI1 = df$quantitity_d7_cintI1 / df$quantity_d716S)

#### Make three datasets for each gene target ####

#IntI1 Dataset
df_intI1 <- df %>% select(Rep, Sample, d0_previntI1, d7_previntI1)

#cintI1 dataset
df_cintI1 <- df %>% select(Rep, Sample, d0_prevcintI1, d7_prevcintI1)


#Average the technical replicates for each dataset
df_avg_intI1 <- aggregate(df_intI1[c(3:4)], list(df_intI1$Rep), mean)

df_avg_cintI1 <- aggregate(df_cintI1[c(3:4)], list(df_cintI1$Rep), mean)


#Melt for use
df_melt_intI1 <- melt(df_avg_intI1, id=c("Group.1"))

df_melt_cintI1 <- melt(df_avg_cintI1, id=c("Group.1"))


#Split Group.1 into Sample and Rep

df_melt_intI1 <- separate(df_melt_intI1, Group.1, c("Sample", "Rep"), sep = " ")

df_melt_cintI1 <- separate(df_melt_cintI1, Group.1, c("Sample", "Rep"), sep = " ")

#Order the levels 
df_melt_intI1$Sample <- as.factor(df_melt_intI1$Sample)

levels(df_melt_intI1$Sample) <- c("0", "4.5", "13.", "40.7", "122", "367", "1100", "3300")


df_melt_cintI1$Sample <- as.factor(df_melt_cintI1$Sample)

levels(df_melt_cintI1$Sample) <- c("0", "4.5", "13.", "40.7", "122", "367", "1100", "3300")


###################### intI1 Stats ###################################

#Make new stats dfs
df_intI1_stats <- df_melt_intI1

df_intI1_stats$concrep <- paste(df_intI1_stats$Sample, df_intI1_stats$Rep)


#intI1 stats 
m1 <- lmer(value ~ variable * Sample + (1|concrep), data=df_intI1_stats)

summary(m1)

plot(m1)
qqnorm(resid(m1))
qqline(resid(m1))

simulationOutput <- simulateResiduals(fittedModel = m1, plot = T)
#bad fit 

m1.2 <- lmer(log(value) ~ variable * Sample + (1|concrep), data=df_intI1_stats)

plot(m1.2)
qqnorm(resid(m1.2))
qqline(resid(m1.2))

simulationOutput <- simulateResiduals(fittedModel = m1.2, plot = T)

#The fit of this is very nice

summary(m1.2)

drop1(m1.2, test="Chisq")
#The interaction term is not significant so we take it out of the model p=0.08137 

m1.3 <- lmer(log(value) ~ variable + Sample + (1|concrep), data=df_intI1_stats)

plot(m1.3)
qqnorm(resid(m1.3))
qqline(resid(m1.3))

ggplot(data = df_intI1_stats, aes(x = resid(m1.3))) +
  geom_histogram(bins = 10, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

simulationOutput <- simulateResiduals(fittedModel = m1.3, plot = T)

anova(m1.2, m1.3)
#m1.3 is smaller AIC 
#Use m1.3, no interaction

summary(m1.3)

Anova(m1.3, test="Chisq", type="III") 

emmeans(m1.3, list(pairwise ~ variable+Sample), adjust ="fdr")


#### cintI1 stats ####

df_cintI1_stats <- df_melt_cintI1

df_cintI1_stats$concrep <- paste(df_cintI1_stats$Sample, df_cintI1_stats$Rep)


m2 <- lmer(value ~ variable * Sample + (1|concrep), data=df_cintI1_stats)

summary(m2)

plot(m2)
qqnorm(resid(m2))
qqline(resid(m2))

#Very badly fitting
#Try logged data

m2.2 <- lmer(log(value) ~ variable * Sample + (1|concrep), data=df_cintI1_stats)

summary(m2.2)

plot(m2.2)
qqnorm(resid(m2.2))
qqline(resid(m2.2))

ggplot(data = df_cintI1_stats, aes(x = resid(m2.2))) +
  geom_histogram(bins = 10, fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')


simulationOutput <- simulateResiduals(fittedModel = m2.2, plot = T)

drop1(m2.2, test="Chisq")
#The interaction term is significant

Anova(m2.2, test="Chisq", type="III") 

emmeans(m2.2, list(pairwise ~ variable*Sample), adjust ="fdr")




#### Plots ####
#Changes size of fonts
p_intI1 <- ggplot(df_melt_intI1, aes(x=Sample, y=value, fill=variable))+
  geom_boxplot(outlier.shape =  NA, lwd = 0.8)+
  geom_point(position=position_jitterdodge(), alpha=0.7, size = 2, shape = 21, show.legend = FALSE)+
  theme_bw()+
  labs(x="Metformin concentration (µg/L)", y="intI1  prevalence")+
  ylab(expression(paste(italic("intI1 "), "prevalence")))+
  scale_fill_manual(values = met.brewer("Hokusai2", 2), name="Day", labels=c("D0", "D7")) +
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size=20))



pcintI1 <- ggplot(df_melt_cintI1, aes(x=Sample, y=value, fill=variable))+
  geom_boxplot(outlier.shape =  NA, lwd = 0.8)+
  geom_point(position=position_jitterdodge(), alpha=0.7, size = 2, shape = 21, show.legend = FALSE)+
  theme_bw()+
  labs(x="Metformin concentration (µg/L)", y="cintI1 prevalence")+
  ylab(expression(paste(italic("cintI1"), "prevalence")))+
  scale_fill_manual(values = met.brewer("Hokusai2", 2), name="Day", labels=c("D0", "D7")) +
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size=20))


p_intI1 / pcintI1 + plot_annotation(tag_levels = 'A') + plot_layout(guides = "collect")

ggsave("Supplementary Figure 2.tiff", height = 10, width = 8, dpi = 300, path = "figures")

