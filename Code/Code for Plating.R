#### Phenotypic Plating Analysis ####

#### Load Libraries ####
library(tidyverse)
library(reshape2)
library(ggplot2)
library(MASS) 
library(MetBrewer)
library(lme4) #makes the model
library(DHARMa) #testst the model is happy
library(elliptic) #to set limits
library(emmeans)
library(flextable)
library(rstatix)
library(ggpubr)

#### Load Data in ####

df <- read.csv("data/Full Experiment Results v4.csv", header=TRUE)


#### Melting the Dataframe ####
df_melt <- melt(df, id=c("Plate"))

names(df_melt)[names(df_melt)=="Plate"] <- "Treatment"


#Now split the variable column into the separate things you want
df_melt <- separate(df_melt, variable, c("Plate", "Type", "Dilution"), sep="_")

df_melt <- separate(df_melt, Treatment, c("Treatment", "Rep"), sep=" ")

df_melt$Dilution <- as.numeric(as.character(df_melt$Dilution))

#Make CFU/mL column

#we have the 50 becuase we want to times up to 1000uL and we plated 20uL
df_melt$CFU_mL <- df_melt$value * (50*df_melt$Dilution)

#Remove NAs and create a new dataframe
df_use <- na.omit(df_melt)

#Average across the replicates
df_avg <- group_by(df_use, Treatment, Plate, Type) %>%
  summarise(avg=mean(CFU_mL))


#Make plain and antibiotic column
df_plain$Treatment_Type <- paste(df_plain$Treatment, df_plain$Type)

df_avg_noplain$Treatment_Type <- paste(df_avg_noplain$Treatment, df_avg_noplain$Type)


#Add them back in
df_avg2 <- merge(df_avg_noplain, df_plain, by="Treatment_Type")

df_avg2 <- separate(df_avg2, Treatment_Type, c("Treatment", "Bacteria_Type"), sep=" ")


## Make proportion resistant
df_avg2$Prop_Res <- df_avg2$avg / df_avg2$Plain_CFU

df_avg2$Percent_Res <- df_avg2$Prop_Res * 100

#Tidy dataframe - remove unnecessary columns
df_3 <- df_avg2[c('Treatment', 'Bacteria_Type', 'Plate.x', "avg", "Plain_CFU", "Prop_Res", "Percent_Res")]

df_3 <- rename(df_3, "Abx_CFU" = avg)

df_3 <- rename(df_3, "Abx_Plate" = Plate.x)

df_3 <- separate(df_3, Treatment, c("Treatment", "Rep"), sep="_")


#Make if_else statement to cap all values at 1 - cannot have a greater than equal to sensitive growth on antibiotic plate
limitedProp_Res <- limit(df_3$Prop_Res, upper=1)

df_4 <- cbind(df_3, data.frame(limitedProp_Res))



#Make df_ecoli and coliform
df_ecoli <- subset(df_4, Bacteria_Type == "E.coli")

df_coliform <- subset(df_4, Bacteria_Type == "Coliform")



#### Stats ####

#creating microcosm individual observation label
df_stats <- df_4

df_stats$Microcosm <- paste(df_stats$Treatment, df_stats$Rep, sep="_")

df_stats$Microcosm <- as.factor(as.character(df_stats$Microcosm))


#Now we need to limit the raw data so that we can feed in into the model
df_stats$limitedABx_CFU <- if_else(df_stats$Abx_CFU > df_stats$Plain_CFU, df_stats$Plain_CFU, df_stats$Abx_CFU)

df_stats$non_resCFU <- df_stats$Plain_CFU - df_stats$Abx_CFU

df_stats$non_resCFU <- if_else(df_stats$non_resCFU < 0, 0, df_stats$non_resCFU)

df_stats <- mutate(df_stats, obs = 1:n(),
                   total =limitedABx_CFU + non_resCFU)

df_all <- group_by(df_stats, Bacteria_Type, Treatment, Abx_Plate, Microcosm) %>%
  summarise(limitedABx_CFU = sum(limitedABx_CFU), 
            non_resCFU = sum(non_resCFU),
            .groups = "drop") %>%
  mutate(prop = limitedABx_CFU /(limitedABx_CFU+non_resCFU))

df_all$obs <- 1:nrow(df_all)

#Make Abx a factor
df_all$Abx_Plate <- as.factor(as.character(df_all$Abx_Plate))

df_stats$Microcosm <- (as.character(df_stats$Microcosm))



## Stats on E. coli and coliform data together ##

m1<- lmer(sqrt(prop) ~ Treatment + Abx_Plate + (1|Microcosm), data = df_all)

summary(m1)

plot(m1)
qqnorm(resid(m1))
qqline(resid(m1))

ggplot(data = df_all, aes(x = resid(m1))) +
  geom_histogram(fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')


simulationOutput <- simulateResiduals(fittedModel = m1, plot = T)

#Chose this model because although the DHARMa output suggested some deviation from the ideal fit, both the general heterogeneity
#and the qqplots looked good and within range of normal, and  the histogram of residuals showed a normal distribution

#Previous models fitting showed that the multiplicative interaction between NAD treatment and type of antibiotic in the plate
#was not significant, so was removed. 

Anova(m1, test="Chisq", type="III") 
#NAD Treatment p=0.07685 
#Abx plate p=<2e-16

#i.e. NAD treatment does not affect phenotypic resistance in total E. coli and coliform 


## Stats for E. coli data ##

#df_all subset for e.coli 
df_stats_ecoli <- subset(df_all, Bacteria_Type == "E.coli" )


m2 <- lmer((prop^(1/3)) ~ Treatment + Abx_Plate + (1|Microcosm), data = df_stats_ecoli)

summary(m2)

plot(m2)
qqnorm(resid(m2))
qqline(resid(m2))

simulationOutput <- simulateResiduals(fittedModel = m2, plot = T)

ggplot(data = df_stats_ecoli, aes(x = resid(m2))) +
  geom_histogram(fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')


#Use the cube root transformation since the residuals fit nicely

Anova(m2.cube2, test="Chisq", type="III") 

#Abx Plate = p<2e-16
#NAD Treatment (NAD) p=0.5829    

#For E. coli, there is no effect of NAD treatment on the phenotypic resistance



## Coliform Stats ##

m3 <- lmer(sqrt(prop) ~ Treatment + Abx_Plate +(1|Microcosm), data=df_stats_coliform)

summary(m3)

plot(m3)

qqnorm(resid(m3))
qqline(resid(m3))

ggplot(data = df_stats_coliform, aes(x = resid(m3))) +
  geom_histogram(fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')

simulationOutput <- simulateResiduals(fittedModel = m3, plot = T)

#The fit is not perfect, but works reasonably well. When the model with the multiplicative interaction term is present (i.e. 
#lmer(sqrt(prop) ~ Treatment * Abx_Plate +(1|Microcosm))) is used, the fit is very good, but the interaction term is not significant,
#so has to be removed

Anova(m3, test="Chisq", type="III") 

#Treatment p=0.01323
#Avx plate p=<2e-16

emmeans(m3, list(pairwise ~ Treatment+Abx_Plate), type = "response", adjust = "fdr")

output_coliform  <-  emmeans(m3.3, list(pairwise ~ Treatment*Abx_Plate)) %>% as.data.frame()

options(max.print = 10000)

emmeans(m3, list(pairwise ~ Treatment*Abx_Plate), type = "response", adjust = "none")

#Read out the file, and read back in to remove unnecessary comparisons (i.e. Diclofenac Azi - Metformin Azi) since these 
#are not sensible comparisons

#Reload file and then add adjusted p values
output_coliform1 <- read.csv("data/unadjusted p values coliform.csv", header=TRUE)

output_coliform1 <- mutate(output_coliform1, padj = p.adjust(p.value, method = 'fdr'))

#Add significance values to boxplot

stat.test_coliform <-  tibble::tribble(
  ~Abx_Plate, ~group1, ~group2, ~p.adj, ~y, ~Treatment,
  "Amp", "Control", "17-ß-Estradiol", '**', "limitedProp_Res", "Control",
  "Gen", "Control", "17-ß-Estradiol", "*", "limitedProp_Red", "Control")



#### Boxplots ####

#Add label vector for the antibiotic
antibiotic.label <- c("Ampicillin", "Azithromycin", "Cefotaxime", "Ciprofloxacin", "Gentamicin", "Tetracycline")
names(antibiotic.label) <- c("Amp", "Azi", "Cef", "Cip", "Gen", "Tet")



#Boxplot for Supplementary Figure XXX
#Plating results for all E. coli and coliform colonies combined


#Total Growth (E.coli and coliforms together)
ggplot(df_final, aes(x=Treatment, y=limitedProp_Res, fill=Treatment)) +
  geom_boxplot(outlier.shape =  NA, lwd = 0.8)+
  geom_point(position=position_jitterdodge(), alpha=0.7, size = 2, shape = 21, show.legend = FALSE)+
  facet_wrap("Abx_Plate", scales = "free", labeller = labeller(Abx_Plate = antibiotic.label))+
  scale_fill_manual(values = c("#ddc000", "#79ad41", "#34b6c6", "#d7aca1"), name="Treatment")+
  theme_bw()+
  labs(x="Treatment", y="Proportion Resistant CFU/mL")+
  theme(text=element_text(size=20))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  theme(axis.title = element_text(size=20))+
  theme(legend.position="none",
        strip.background=element_rect(colour="black",
                                      fill="white"))

ggsave("Total E.coli and Coliform Data.tiff", path = "figures", dpi = 300, height = 10, width = 12)


#E.coli Data
ggplot(df_ecoli, aes(x=Treatment, y=limitedProp_Res, fill=Treatment))+
  geom_boxplot(outlier.shape =  NA, lwd = 0.8)+
  geom_point(position=position_jitterdodge(), alpha=0.7, size = 3, shape = 21, show.legend = FALSE)+
  # geom_point(position="jitter", size = 2)+
  facet_wrap("Abx_Plate", scales = "free", labeller = labeller(Abx_Plate = antibiotic.label))+
  scale_fill_manual(values = c("#ddc000", "#79ad41", "#34b6c6", "#d7aca1"), name="Treatment")+
  theme_bw()+
  labs(x="Treatment", y="Proportion Resistant CFU/mL")+
  theme(text=element_text(size=20))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  theme(axis.title = element_text(size=20))+
  theme(legend.position="none",
        strip.background=element_rect(colour="black",
                                      fill="white"))

ggsave("E coli Data.tiff", path = "figures", dpi = 300, height = 10, width = 12)


## Coliform Data 
ggplot(df_coliform, aes(x=Treatment, y=limitedProp_Res))+
  geom_boxplot(outlier.shape =  NA, lwd = 0.8, aes(fill=Treatment))+
  geom_point(position=position_jitterdodge(), alpha=0.7, size = 3, shape = 21, show.legend = FALSE,  aes(fill=Treatment))+
  facet_wrap("Abx_Plate", scales = "free", labeller = labeller(Abx_Plate = antibiotic.label))+
  scale_fill_manual(values = c("#ddc000", "#79ad41", "#34b6c6", "#d7aca1"), name="Treatment")+
  theme_bw()+
  labs(x="Treatment", y="Proportion Resistant CFU/mL")+
  theme(text=element_text(size=20))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  theme(axis.title = element_text(size=20))+
  theme(legend.position = "none", 
        strip.background=element_rect(colour="black", fill="white"))+
  stat_pvalue_manual(stat.test_coliform, label = "p.adj", tip.length = 0.04, y.position = 1.4, size = 5, bracket.size = 0.9)


ggsave("Coliform Data", path = "figures", dpi = 300, height = 8, width = 12)





