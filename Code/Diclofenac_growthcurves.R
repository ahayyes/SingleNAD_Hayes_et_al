###### Growth Curves and SELECT for Diclofenac ######

####Load Libraries ####
library(tidyverse)
library(ggplot2)
library(reshape)
library(dplyr)
library(growthcurver)
library(MetBrewer)
library(dunn.test)
library(plotrix)
library(patchwork)
library(emmeans)
library(DHARMa)

#Read in data
total <- read.csv("data/diclofenac_total_matrix.csv", header=TRUE)


#### SELECT Working to find LOEC ####
#shapiro test - normal if the p value is above 0.05
#want to find the point in exponential phase with the largest dose response. this is hour five for this data
shapiro.test(total$t6)
cor.test(total$conc, total$t6, method = c("spearman"))

#Time pint with greatest difference = 5 
dunn.test(total$t6, total$conc)


#Melt the dataframe so we can plot it
reads <- melt(total, id="conc")

#Change the names of the columns for easier interpretation
names(reads)[names(reads)=="variable"] <- "time"

names(reads)[names(reads)=="value"] <- "od"

#Make an averaged dataframe for plotting, including SE
avg <- reads %>% 
  group_by(conc, time) %>%
  summarise(
    mean=mean(od), 
    sd=sd(od), 
    se=std.error(od))

#remove 't' from the time column cells
avg$time <- gsub("t", "", as.character(avg$time), n)

#set conc as discrete 
avg$conc <- ordered(avg$conc)

#Turn your 'treatment' column into a character vector
avg$time <- as.character(avg$time)

#Then turn it back into a factor with the levels in the correct order
avg$time <- factor(avg$time, levels=unique(avg$time))

#Make Plot
ggplot(avg, aes(x=time, y=mean, colour=conc)) +
  geom_line(aes(group=conc), size=1.2) +
  scale_color_manual(values=met.brewer("Hiroshige", 8)) +
  scale_fill_manual(values=met.brewer("Hiroshige", 8))+
  labs(y="Mean OD(600nm)", x="Time (hours)", color="Diclofenac \nConcentration (µg/L)", 
       fill="Diclofenac \nConcentration (µg/L)") +
  labs(y=expression(`Mean OD`[(`600nm`)]))+
  geom_ribbon(aes(fill=conc, ymin=mean-se, ymax=mean+se), alpha=0.1, colour=NA)+
  theme_classic()+
  theme(axis.text=element_text(size=16)) + theme(axis.title=element_text(size = 18)) +
  theme(legend.text = element_text(size = 14)) + theme(legend.title = element_text(size = 16))



