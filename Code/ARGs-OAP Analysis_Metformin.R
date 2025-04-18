#### Metformin ARGs-OAP ####


#### Load Libraries ####

library(tidyverse)
library(ggplot2)
library(reshape2) #to melt the dataframe
library(broom) #For the significance testing of ARGs 
library(dunn.test)
library(MetBrewer)
library(vegan)


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
df_melt <- melt(dataset, id=c("antibiotic.class"))

#Change the name of the sample column
df_melt <- rename(df_melt, Sample=variable)


#Need to create a new column to split Treatment and Replicate
df_melt <- separate(df_melt, Sample, c('Treatment', 'Replicate'), sep=" ")

#Now average replicates within treatments
df_melt_avg <- group_by(df_melt, Treatment, antibiotic.class) %>% 
  summarise(Avg=mean(value), .groups="drop")


#Without unclassified and multidrug
#Subset the data to have both of these removed - might be easier for later downstream stuff
df_melt_avg_small <- subset(df_melt_avg, antibiotic.class!="unclassified")
df_melt_avg_small <- subset(df_melt_avg_small, antibiotic.class!="multidrug")

df_melt_multidrug <- subset(df_melt, antibiotic.class == "multidrug")
df_melt_avg_multidrug <- subset(df_melt_avg, antibiotic.class=="multidrug")



# Make a heatmap of all ARG classes, by treatment, including inoculum
df_melt_avg_small$Treatment <- as.factor(df_melt_avg_small$Treatment)

df_melt_avg_small$Treatment <- factor(df_melt_avg_small$Treatment, levels=unique(df_melt_avg_small$Treatment))

df_melt_avg_small$Treatment <- factor(df_melt_avg_small$Treatment, levels=c("inoc", "0", "4.5", "13.6", "40.7", "122", "367", "1100", "3300"))


# Heatmaps without inoculum
# Subset Data for Dataframe w/o inoculum
df_melt_avg_noinoc <- subset(df_melt_avg_small, Treatment!= "inoc")

# For all antibiotic resistance genes
ggplot(df_melt_avg_noinoc, aes(x=Treatment, y=reorder(antibiotic.class, desc(antibiotic.class)), fill=Avg))+
  geom_tile()+
  scale_fill_gradientn(colors=met.brewer("Hokusai2"))+
  theme_classic()+
  theme(text=element_text(size=14))+
  labs(x="Metformin concentration (µg/L)", y="Antibiotic Resistance Gene Class", fill = "Abundance")


ggsave("Supplementary Figure 10.tiff", height = 6, width = 8, path= "figures")


#### Stats ####

# run a Kruskal wallis test on each gene ###

#To test if there are differences across treatments, make new dataframe
df_stats <- subset(df_melt, Treatment!="inoc")


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
  temp_data <- filter(df_stats, antibiotic.class==temp_amr_class) %>%
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



#### NMDS ####

#do this on only evolved treatments 

#read in transposed dataframe
d_nmds <- read.csv("Transposed Dataframe_total_metformin.csv", header=TRUE)
d_nmds <- rename(d_nmds, Sample=antibiotic.class)

#make Treatment and Rep Columns
d_nmds <- separate(d_nmds, Sample, c("Treatment", "Rep"), sep=" ")


#exclude inoculum
d_nmds <- subset(d_nmds, Treatment!="inoc")


#do the NMDS 
pc = d_nmds

#make community matrix - extract columns with abundance information
com = pc[ ,3:ncol(pc)]
m_com = as.matrix(com)#Set seed gives you the same results each time you run this code

set.seed(123)
nmds = metaMDS(m_com, distance = "bray")
nmds

data.scores = as.data.frame(scores(nmds)$sites)


#add columns to data frame 
data.scores$Sample = pc$Treatment


head(data.scores)

data.scores$Sample <- factor(data.scores$Sample, levels=unique(data.scores$Sample))

data.scores$Sample <- factor(data.scores$Sample, levels=unique(data.scores$Sample))

data.scores$Sample <- factor(data.scores$Sample, levels=c("0", "4.5", "13.6", "40.7", "122", "367", "1100", "3300"))

ano = anosim(m_com, pc$Treatment, distance = "bray", permutations = 9999)
ano


ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  theme_classic()+
  geom_point(size = 4, aes(colour = Sample))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12), 
        axis.text.x = element_text(colour = "black", size = 12),
        legend.text = element_text(size = 12, colour ="black"),
        legend.position = "right", axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        legend.key=element_blank()) +
  labs(x = "NMDS1", colour = "Metformin \nConcentration (µg/L)", y = "NMDS2", shape = "Type")+   
  scale_colour_manual(values = met.brewer("Hokusai2", 9)) 





