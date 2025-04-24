############### R Analyses for Transcriptomics Data ########################

########################### Load Libraries ############################

library(tidyverse)
library(ggplot2)
library(reshape2)
library(broom) 
library(dunn.test)
library(MetBrewer)
library(vegan)

################### For all Antibiotic Classes ############################

################# Read in Code and Merge into one file ######################

list_files <- list.files(path="data_class", full.names=TRUE)

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
dataset <- rename(dataset, "Control 1" = ControlrepA)
dataset <- rename(dataset, "Control 2" = ControlrepB)
dataset <- rename(dataset, "Control 3" = ControlrepC)
dataset <- rename(dataset, "Control 3_2" = ControlrepCV2)
dataset <- rename(dataset, "Control 4" = ControlrepD)
dataset <- rename(dataset, "Control 5" = ControlrepE)
dataset <- rename(dataset, "Estradiol 1" = C17repA)
dataset <- rename(dataset, "Estradiol 2" = C17repB)
dataset <- rename(dataset, "Estradiol 3" = C17repC)
dataset <- rename(dataset, "Estradiol 4" = C17repD)
dataset <- rename(dataset, "Estradiol 5" = C17repE)
dataset <- rename(dataset, "Diclofenac 1" = DrepA)
dataset <- rename(dataset, "Diclofenac 2" = DrepBV2)
dataset <- rename(dataset, "Diclofenac 3" = DrepC)
dataset <- rename(dataset, "Diclofenac 4" = DrepD)
dataset <- rename(dataset, "Diclofenac 5" = DrepE)
dataset <- rename(dataset, "Metformin 1" = MrepA)
dataset <- rename(dataset, "Metformin 1_2" = MrepAV2)
dataset <- rename(dataset, "Metformin 2" = MrepB)
dataset <- rename(dataset, "Metformin 2_2" = MrepBV2)
dataset <- rename(dataset, "Metformin 3" = MrepC)
dataset <- rename(dataset, "Metformin 4" = MrepD)
dataset <- rename(dataset, "Metformin 5" = MrepE)

# Some columns were sequenced twice to account for low yield
# Need to add these columns together
# aka Metformin 1 and Metformin 1_2

dataset <- dataset %>% mutate("Metformin 1"=c(`Metformin 1`+ `Metformin 1_2`)) %>% select(-`Metformin 1_2`)

dataset <- dataset %>% mutate("Metformin 2"=c(`Metformin 2`+ `Metformin 2_2`)) %>% select(-`Metformin 2_2`)

dataset <- dataset %>% mutate("Control 3"=c(`Control 3`+ `Control 3_2`)) %>% select(-`Control 3_2`)



# Melt dataset so can be used in ggplot
df_melt <- melt(dataset, id=c("antibiotic.class"))

# Change the name of the sample column
df_melt <- rename(df_melt, Sample=variable)

#Need to create a new column to split Treatment and Replicate
df_melt <- separate(df_melt, Sample, c('Treatment', 'Replicate'), sep=" ")

#Now average replicates within treatments
df_melt_avg <- group_by(df_melt, Treatment, antibiotic.class) %>% 
  summarise(Avg=mean(value), .groups="drop")

#Make datasets without multidrug and unclassified just in case 

df_melt_avg_nounclassified <- subset(df_melt_avg, antibiotic.class!="unclassified")

df_melt_nomulti <- subset(df_melt, antibiotic.class!= "multidrug")

df_melt_avg_nounclass_multi <- subset(df_melt_avg_nounclassified, antibiotic.class!="multidrug")

#Make column with concRep
df_melt$concrep <- paste(df_melt$Treatment, df_melt$Replicate)

df_melt_small <- subset(df_melt, antibiotic.class!="multidrug")
df_melt_small <- subset(df_melt_small, antibiotic.class!="unclassified")


############################# Stats ############################################

#Kruskal Wallis Test on each gene, and then do a p-value adjustment

# set up empty results dataframe
models_kruskal <- select(dataset, antibiotic.class) %>%
  distinct() %>%
  mutate(data = list(NA),
         model = list(NA),
         summary = list(NA))

for(i in 1:nrow(models_kruskal)){
  
  # grab amr class
  temp_amr_class <- models_kruskal$antibiotic.class[i]
  
  # filter just for that amr class then the columns we need
  temp_data <- filter(df_melt, antibiotic.class==temp_amr_class) %>%
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
model_output <- mutate(model_output, padj = p.adjust(p.value, method = 'fdr'))

#No significant changes prior or post adjustment for multiple testing


####################### NMDS on Gene Class Diversity ###########################
write.csv(dataset, "To be transposed gene classes.csv")

d_nmds <- read.csv("fullclass_transposed.csv", header=TRUE)

#Make Treatment and Rep Columns
d_nmds <- separate(d_nmds, treatment, c("Treatment", "Rep"), sep=" ")

#Do the NMDS 
pc = d_nmds

#make community matrix - extract columns with abundance information
com = pc[ ,3:ncol(pc)]
m_com = as.matrix(com)


#Set seed gives you the same results each time you run this code
set.seed(123)
nmds = metaMDS(m_com, distance = "bray")
nmds

#Plot this 
plot(nmds)

#vegan changed package as it is 2.6 and above use this code to get the data scores

## THIS IS CODE FOR OLDER VERSIONS OF VEGAN ##
#Extract Info to Add to Plot
# data.scores = as.data.frame(scores(nmds))

## THIS IS CURRENT CODE FOR DATA SCORES ##

#extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
data.scores = as.data.frame(scores(nmds)$sites)

#add columns to data frame 
data.scores$Sample = pc$Treatment


data.scores$Sample <- factor(data.scores$Sample, levels=unique(data.scores$Sample))

data.scores$Sample <- factor(data.scores$Sample, levels=c("Control", "Diclofenac", "Metformin", "Estradiol"))


#ANOSIM
ano = anosim(m_com, pc$Treatment, distance = "bray", permutations = 9999)
ano


#################### For all ARGs within the classes ###########################


#### Read in Code and Merge into one file ####

list_files <- list.files(path="data_subclass", full.names=TRUE)

for (file in list_files){
  #if the merged dataset doesn't exist, create it
  if (!exists("dataset2")){
    dataset2 <- read.delim(file, header=TRUE, sep="\t")
  }
  #if the merged dataset does exist, append to it 
  if (exists("dataset2")){
    temp_dataset2 <- read.delim(file, header=TRUE, sep="\t")
    dataset2 <- merge(dataset2, temp_dataset2)
    rm(temp_dataset2)
  }
}



#Change column names 
dataset2 <- rename(dataset2, "Control 1" = ControlrepA)
dataset2 <- rename(dataset2, "Control 2" = ControlrepB)
dataset2 <- rename(dataset2, "Control 3" = ControlrepC)
dataset2 <- rename(dataset2, "Control 3_2" = ControlrepCV2)
dataset2 <- rename(dataset2, "Control 4" = ControlrepD)
dataset2 <- rename(dataset2, "Control 5" = ControlrepE)
dataset2 <- rename(dataset2, "Estradiol 1" = C17repA)
dataset2 <- rename(dataset2, "Estradiol 2" = C17repB)
dataset2 <- rename(dataset2, "Estradiol 3" = C17repC)
dataset2 <- rename(dataset2, "Estradiol 4" = C17repD)
dataset2 <- rename(dataset2, "Estradiol 5" = C17repE)
dataset2 <- rename(dataset2, "Diclofenac 1" = DrepA)
dataset2 <- rename(dataset2, "Diclofenac 2" = DrepBV2)
dataset2 <- rename(dataset2, "Diclofenac 3" = DrepC)
dataset2 <- rename(dataset2, "Diclofenac 4" = DrepD)
dataset2 <- rename(dataset2, "Diclofenac 5" = DrepE)
dataset2 <- rename(dataset2, "Metformin 1" = MrepA)
dataset2 <- rename(dataset2, "Metformin 1_2" = MrepAV2)
dataset2 <- rename(dataset2, "Metformin 2" = MrepB)
dataset2 <- rename(dataset2, "Metformin 2_2" = MrepBV2)
dataset2 <- rename(dataset2, "Metformin 3" = MrepC)
dataset2 <- rename(dataset2, "Metformin 4" = MrepD)
dataset2 <- rename(dataset2, "Metformin 5" = MrepE)


dataset2 <- dataset2 %>% mutate("Metformin 1"=c(`Metformin 1`+ `Metformin 1_2`)) %>% select(-`Metformin 1_2`)

dataset2 <- dataset2 %>% mutate("Metformin 2"=c(`Metformin 2`+ `Metformin 2_2`)) %>% select(-`Metformin 2_2`)

dataset2 <- dataset2 %>% mutate("Control 3"=c(`Control 3`+ `Control 3_2`)) %>% select(-`Control 3_2`)

#Melt dataset so can be used in ggplot
df_melt2 <- melt(dataset2, id=c("antibiotic.gene"))

#Change the name of the sample column
df_melt2 <- rename(df_melt2, Sample=variable)


#Need to create a new column to split Treatment and Replicate
df_melt2 <- separate(df_melt2, Sample, c('Treatment', 'Replicate'), sep=" ")


#Now average replicates within treatments
df_melt_avg2 <- group_by(df_melt2, Treatment, antibiotic.gene) %>% 
  summarise(Avg=mean(value), .groups="drop")


#Make datasets without multidrug and unclassified just in case 

df_melt_avg2_nounclassified <- subset(df_melt_avg2, antibiotic.gene!="unclassified")

df_melt_nomulti2 <- subset(df_melt2, antibiotic.gene!= "multidrug")

df_melt_avg_small2 <- subset(df_melt_avg2_nounclassified, antibiotic.gene!="multidrug")

############################## Stats on all Genes ########################

#Kruskal Wallis Test on each gene, and then do a p-value adjustment
# run a Kruskal wallis test on each gene ###

# set up empty results dataframe
models_kruskal2 <- select(dataset2, antibiotic.gene) %>%
  distinct() %>%
  mutate(data = list(NA),
         model = list(NA),
         summary = list(NA))

for(i in 1:nrow(models_kruskal2)){
  
  # grab amr class
  temp_amr_class <- models_kruskal2$antibiotic.gene[i]
  
  # filter just for that amr class then the columns we need
  temp_data <- filter(df_melt2, antibiotic.gene==temp_amr_class) %>%
    select(Treatment, value)
  
  # run kruskal wallis test
  temp_model <- kruskal.test(value ~ Treatment, data=temp_data)
  
  # get tidy output
  temp_output <- tidy(temp_model)
  
  # assign each output to the correct section in the results dataframe
  models_kruskal2$data[[i]] <- temp_data
  models_kruskal2$model[[i]] <- temp_model
  models_kruskal2$summary[[i]] <- temp_output
  
}

models_kruskal2

# look at one of the models
models_kruskal2$model[[1]]

# a function called unnest() allows us to access different bits of the results dataframe

# 1. get the data frame again
unnest(models_kruskal2, data) %>% select(-model, -summary)

# 2. get the output of the model
model_output2 <- unnest(models_kruskal2, summary) %>% select(-model, -data) 

# add in a column for p adjustments
model_output2 <- mutate(model_output2, padj = p.adjust(p.value, method = 'fdr'))

#No significant changes  post adjustment for multiple testing

#################### NMDS on all ARGs ######################################

write.csv(dataset2, "To be transposed all genes.csv")

d_nmds <- read.csv("transposed_allgenes.csv", header=TRUE)

#Make Treatment and Rep Columns
d_nmds <- separate(d_nmds, treatment, c("Treatment", "Rep"), sep=" ")

#Do the NMDS 
pc = d_nmds

#make community matrix - extract columns with abundance information
com = pc[ ,3:ncol(pc)]
m_com = as.matrix(com)

set.seed(123)
nmds = metaMDS(m_com, distance = "bray")
nmds

#Plot this 
plot(nmds)

#BC vegan changed package as it is 2.6 and above use this code to get the data scores

#extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
data.scores = as.data.frame(scores(nmds)$sites)

#add columns to data frame 
data.scores$Sample = pc$Treatment


data.scores$Sample <- factor(data.scores$Sample, levels=unique(data.scores$Sample))

data.scores$Sample <- factor(data.scores$Sample, levels=c("Control", "Diclofenac", "Metformin", "Estradiol"))


ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  theme_classic()+
  scale_color_manual(values = c("#eace17", "#6dbc18", "#34b6c6", "#ed6a48"), name="Treatment",
                     labels=c("Control", "Diclofenac", "Metformin", "17-β-estradiol"))+
  geom_point(size = 4, aes(colour = Sample))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12), 
        axis.text.x = element_text(colour = "black", size = 12),
        legend.text = element_text(size = 12, colour ="black"),
        legend.position = "right", axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.key=element_blank()) +
  labs(x = "NMDS1", colour = "Treatment", y = "NMDS2", shape = "Type")




#ANOSIM
ano = anosim(m_com, pc$Treatment, distance = "bray", permutations = 9999)
ano


