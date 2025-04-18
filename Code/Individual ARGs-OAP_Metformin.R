#### Individual gene analyses Metformin ####

#### Load Libraries ####

library(tidyverse)
library(ggplot2)
library(reshape2) #to melt the dataframe
library(broom) #For the significance testing of ARGs 
library(dunn.test)
library(vegan)
library(MetBrewer)
library(data.table) #To help reshape the dataframe
library(DHARMa)
library(drc)
library(flextable) #to make table for word
library(officer)
library(dunn.test)

#### Read in Code and Merge into one file ####

list_files <- list.files(path="data", full.names=TRUE)

for (file in list_files){
  #if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- read.delim(file, header=TRUE, sep="\t")
  }
  #if the merged dataset does exist, append to it 
  if (exists("dataset")){
    temp_dataset <- read.delim(file, header=TRUE, sep="\t")
    dataset <- merge(dataset, temp_dataset)
    rm(temp_dataset)
  }
}

#Change column names 
dataset <- rename(dataset, "0 3" = metconcHr3)
dataset <- rename(dataset, "0 2" = metconcHr2)
dataset <- rename(dataset, "0 1" = metconcHr1)
dataset <- rename(dataset, "3300 3" = metconcAr3)
dataset <- rename(dataset, "3300 2" = metconcAr2)
dataset <- rename(dataset, "3300 1" = metconcAr1)
dataset <- rename(dataset, "1100 1" = metconcBr1)
dataset <- rename(dataset, "1100 2" = metconcBr2)
dataset <- rename(dataset, "1100 3" = metconcBr3)
dataset <- rename(dataset, "367 3" = metconcCr3)
dataset <- rename(dataset, "367 2" = metconcCr2)
dataset <- rename(dataset, "367 1" = metconcCr1)
dataset <- rename(dataset, "122 3" = metconcDr3)
dataset <- rename(dataset, "122 2" = metconcDr2)
dataset <- rename(dataset, "122 1" = metconcDr1)
dataset <- rename(dataset, "40.7 3" = metconcEr3)
dataset <- rename(dataset, "40.7 2" = metconcEr2)
dataset <- rename(dataset, "40.7 1" = metconcEr1)
dataset <- rename(dataset, "13.6 3" = metconcFr3)
dataset <- rename(dataset, "13.6 2" = metconcFr2)
dataset <- rename(dataset, "13.6 1" = metconcFr1)
dataset <- rename(dataset, "4.5 3" = metconcGr3)
dataset <- rename(dataset, "4.5 2" = metconcGr2)
dataset <- rename(dataset, "4.5 1" = metconcGr1)
dataset <- rename(dataset, "inoc 3" = inocrep3)
dataset <- rename(dataset, "inoc 2" = inocrep2)
dataset <- rename(dataset, "inoc 1" = inocrep1)

#Melt dataset so can be used in ggplot
df_melt <- melt(dataset, id=c("genes"))

#Change the name of the sample column
df_melt <- rename(df_melt, Sample=variable)


#Need to create a new column to split Treatment and Replicate
df_melt <- separate(df_melt, Sample, c('Treatment', 'Replicate'), sep=" ")

#Now average replicates within treatments
df_melt_avg <- group_by(df_melt, Treatment, genes) %>% 
  summarise(Avg=mean(value), .groups="drop")

df_melt_avg_evolved <- subset(df_melt_avg, Treatment!= "inoc")

dataset_evolved <- subset(dataset, select= -c(26:28))

#Melt evolved dataset so can be used in ggplot
df_melt_evolved <- melt(dataset_evolved, id=c("genes"))

#Change the name of the sample column
df_melt_evolved <- rename(df_melt_evolved, Sample=variable)


#Need to create a new column to split Treatment and Replicate
df_melt_evolved <- separate(df_melt_evolved, Sample, c('Treatment', 'Replicate'), sep=" ")



#### Stats ####

# run a Kruskal wallis test on each gene ###

#To test if there are differences across treatments, make new dataframe
df_stats <- subset(df_melt, Treatment!="inoc")


# set up empty results dataframe
models_kruskal <- select(dataset, genes) %>%
  distinct() %>%
  mutate(data = list(NA),
         model = list(NA),
         summary = list(NA))

for(i in 1:nrow(models_kruskal)){
  
  # grab amr class
  temp_amr_class <- models_kruskal$genes[i]
  
  # filter just for that amr class then the columns we need
  temp_data <- filter(df_stats, genes==temp_amr_class) %>%
    select(Treatment, value)
  
  # run kruskal wallis test
  temp_model <- kruskal.test(value ~ Treatment, data=temp_data)
  
  # get tidy output
  temp_output <- tidy(temp_model)
  
  # assign each output to the correct section in the results dataframe
  models_kruskal$data[[i]] <- temp_data
  models_kruskal$model[[i]] <- temp_model
  models_kruskal$summary[[i]] <- temp_output
  
}

models_kruskal

# look at one of the models
models_kruskal$model[[1]]

# a function called unnest() allows us to access different bits of the results dataframe

# 1. get the data frame again
unnest(models_kruskal, data) %>% select(-model, -summary)

# 2. get the output of the model
model_output <- unnest(models_kruskal, summary) %>% select(-model, -data) 

# add in a column for p adjustments
# can look at this column to determine which are significant across treatments
model_output <- mutate(model_output, padj = p.adjust(p.value, method = 'fdr'))




#### ToLC analyses (since we expect TolC to be affected by Metformin ####
#TolC appears to be changed with metformin treatment, so we are testing this specifically

tolc <- subset(df_melt_evolved, genes=="multidrug__TolC")

#Make conc numeric 
tolc$Treatment<- as.numeric(tolc$Treatment)

#Do Box Cox transformation to improve this fit

tolc_bxcx <- boxcox(tolc_m1, method="anova")

summary(tolc_bxcx)

#Look at residuals
plot(residuals(tolc_bxcx) ~ fitted(tolc_bxcx), main="Residuals vs Fitted")

qqnorm(residuals(tolc_bxcx))
qqline(residuals(tolc_bxcx))


#Predict confidence intervals from model

#Create dataframe first
newdata <- expand.grid(conc=exp(seq(log(0.1), log(3300), length=300)))

# predictions and confidence intervals
pm <- predict(tolc_bxcx, newdata=newdata, interval="confidence")


# new data with predictions
newdata$p <- pm[,1]
newdata$pmin <- pm[,2]
newdata$pmax <- pm[,3]

#Change the  first concentration so that ggplot will plot it 
tolc$conc0 <- tolc$Treatment
tolc$conc0[tolc$conc0 == 0] <- 0.1

#Plot
ggplot(tolc, aes(x = conc0, y = value)) +
  geom_point(size = 2, colour="#004b77", fill="#004b77") +
  theme_classic()+
  geom_ribbon(data=newdata, aes(x=conc, y=p, ymin=pmin, ymax=pmax), alpha=0.2, fill="#c7dde0") +
  geom_line(data=newdata, aes(x=conc, y=p), linewidth = 1, colour="#5cb0b8") +
  coord_trans(x="log") +
  ylim(0.02, 0.05)+
  scale_x_continuous(breaks=c(1, 10, 100, 1000, 3000))+
  theme(text=element_text(size=18))+
  theme(axis.title = element_text(size=20))+
  xlab("Metformin Concentration (µg/mL)") + ylab("tolC abundance") + 
  ylab(expression(paste(italic("TolC"), " Abundance")))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

#Save out plot
ggsave("Figure 4.tiff", path = "figures", dpi = 300, heigh = 6, width=12)


#Which conc are significantly different
dunn.test(tolc$value, tolc$Treatment)


