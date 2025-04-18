############### BacMet Analysis for 17-beta-estradiol #############################


##################### Read in code and libraries ###############################

library(tidyverse) 
library(reshape2) 
library(ggplot2) 
library(MetBrewer) 
library(vegan) 
library(data.table) 
library(broom) 
library(DHARMa)
library(plotrix)
library(patchwork)


#Read in the filenames 
# list files
files <- list.files(path="data", full.names=TRUE)


# Write a function and then use purrr

aprils_cleaning_function <- function(file){
  
  # read in file
  df <- read_tsv(file, col_names=TRUE) %>%
    subset(., Percent.identity>="85") %>%
    subset(., Match.length>="25") %>%
    group_by(Gene) %>%
    tally() %>%
    mutate(., file = basename(tools::file_path_sans_ext(file)))
  
  return(df)
  
}



## Need this later - added in here with admin things ##
#create a function to split each string by an identifier and grab the second index
quick_strsplit <- function(x, index, split){
  # split the string
  temp <- strsplit(x, split = split)
  # iterate through the list and grab each indexed element
  temp <- sapply(temp, "[", index)
  return(temp)
}



#Make dataframe
df <- purrr::map_df(files, aprils_cleaning_function)

df$file <- gsub("bacmetconc", "", df$file)

df$file <- gsub("bacmet_", "", df$file)

#Change names 
df$file[df$file=="Ar1"] <- "5400 1"
df$file[df$file=="Ar2"] <- "5400 2"
df$file[df$file=="Ar3"] <- "5400 3"
df$file[df$file=="Br1"] <- "1800 1"
df$file[df$file=="Br2"] <- "1800 2"
df$file[df$file=="Br3"] <- "1800 3"
df$file[df$file=="Cr1"] <- "607 1"
df$file[df$file=="Cr2"] <- "607 2"
df$file[df$file=="Cr3"] <- "607 3"
df$file[df$file=="Dr1"] <- "201 1"
df$file[df$file=="Dr2"] <- "201 2"
df$file[df$file=="Dr3"] <- "201 3"
df$file[df$file=="Er1"] <- "67 1"
df$file[df$file=="Er2"] <- "67 2"
df$file[df$file=="Er3"] <- "67 3"
df$file[df$file=="Fr1"] <- "22 1"
df$file[df$file=="Fr2"] <- "22 2"
df$file[df$file=="Fr3"] <- "22 3"
df$file[df$file=="Gr1"] <- "7 1"
df$file[df$file=="Gr2"] <- "7 2"
df$file[df$file=="Gr3"] <- "7 3"
df$file[df$file=="Hr1"] <- "2 1"
df$file[df$file=="Hr2"] <- "2 2"
df$file[df$file=="Hr3"] <- "2 3"
df$file[df$file=="Ir1"] <- "0 1"
df$file[df$file=="Ir2"] <- "0 2"
df$file[df$file=="Ir3"] <- "0 3"
df$file[df$file=="bacmet_inocR1"] <- "inoc 1"
df$file[df$file=="bacmet_inocR2"] <- "inoc 2"
df$file[df$file=="bacmet_inocR3"] <- "inoc 3"

colnames(df)[3] <- "Name"

#### Read in 16S Data ####

# cells <- read.table("16S Numbers.txt", header=TRUE, sep="\t")
cells <- read.csv("16S Numbers Excel.csv", header=TRUE)


#Sorting Dataframe to remove unnecessary columns
cells <- cells %>% select(-SampleID) %>% select(Name, X.of16Sreads, CellNumber)

#Rename filenames
cells$Name[cells$Name=="Ar1"] <- "5400 1"
cells$Name[cells$Name=="Ar2"] <- "5400 2"
cells$Name[cells$Name=="Ar3"] <- "5400 3"
cells$Name[cells$Name=="Br1"] <- "1800 1"
cells$Name[cells$Name=="Br2"] <- "1800 2"
cells$Name[cells$Name=="Br3"] <- "1800 3"
cells$Name[cells$Name=="Cr1"] <- "607 1"
cells$Name[cells$Name=="Cr2"] <- "607 2"
cells$Name[cells$Name=="Cr3"] <- "607 3"
cells$Name[cells$Name=="Dr1"] <- "201 1"
cells$Name[cells$Name=="Dr2"] <- "201 2"
cells$Name[cells$Name=="Dr3"] <- "201 3"
cells$Name[cells$Name=="Er1"] <- "67 1"
cells$Name[cells$Name=="Er2"] <- "67 2"
cells$Name[cells$Name=="Er3"] <- "67 3"
cells$Name[cells$Name=="Fr1"] <- "22 1"
cells$Name[cells$Name=="Fr2"] <- "22 2"
cells$Name[cells$Name=="Fr3"] <- "22 3"
cells$Name[cells$Name=="Gr1"] <- "7 1"
cells$Name[cells$Name=="Gr2"] <- "7 2"
cells$Name[cells$Name=="Gr3"] <- "7 3"
cells$Name[cells$Name=="Hr1"] <- "2 1"
cells$Name[cells$Name=="Hr2"] <- "2 2"
cells$Name[cells$Name=="Hr3"] <- "2 3"
cells$Name[cells$Name=="Ir1"] <- "0 1"
cells$Name[cells$Name=="Ir2"] <- "0 2"
cells$Name[cells$Name=="Ir3"] <- "0 3"
cells$Name[cells$Name=="inocR1"] <- "inoc 1"
cells$Name[cells$Name=="inocR2"] <- "inoc 2"
cells$Name[cells$Name=="inocR3"] <- "inoc 3"


# Merge the two dataframes so you can have relative abundance of genes
merged <- merge(df, cells, by = "Name")


#Separate out Name to Conc and Treatment
merged <- tidyr::separate(merged, Name, c("Concentration", "Replicate"), sep=" ")

merged <- merged %>% mutate(Gene_per_16S = merged$n / merged$X.of16Sreads) #This provides us with a dataframe that has relative abundance of BMRGs

#Average out across the biological replicates (for graph reasons)

d_avg <- group_by(merged, Concentration, Gene) %>% 
  summarise(Avg=mean(Gene_per_16S), .groups="drop")

d_avg_evolved <- subset(d_avg, Concentration!="inoc") #This has just the evolved populations



#Set the levels
d_avg_evolved$Concentration <- factor(d_avg_evolved$Concentration, levels=unique(d_avg_evolved$Concentration))
d_avg_evolved$Concentration <- factor(d_avg_evolved$Concentration, levels=c("0", "2", "7", "22", "67", "201", "607", "1800", "5400"))


d_avg$Concentration <- factor(d_avg$Concentration, levels=unique(d_avg$Concentration))
d_avg$Concentration <- factor(d_avg$Concentration, levels=c("inoc", "0","2", "7", "22", "67", "201", "607", "1800", "5400"))


#### Changing the dataframe so that we have zeros (instead of presence/absense) ####

working_merged <- merged

working_merged$concrep <- paste(working_merged$Concentration, working_merged$Replicate)

working_merged <- select(working_merged, Gene, concrep, Gene_per_16S)

working_merged_spread <- spread(working_merged, key=concrep, value=Gene_per_16S)

#Changes NA to 0
working_merged_spread[is.na(working_merged_spread)] <- 0

#Now gather the dataframe again
gathered_merged <- gather(working_merged_spread, key="Concentration", "Relative_Abundance", 2:28)

#Now separate concrep column
gathered_merged <- separate(gathered_merged, Concentration, c("Concentration", "Rep"), sep=" ")

#Make only evolved population dataframe
gathered_merged_evolved <- subset(gathered_merged, Concentration!="inoc")



#### Stats on all BMRGs - evolved populations ####

# set up empty results dataframe - to fill with the models we create (using antibiotic.class as AMR genes)
models_kruskal1 <- select(gathered_merged_evolved, Gene) %>%
  distinct() %>%
  mutate(data = list(NA),
         model = list(NA),
         summary = list(NA))

for(i in 1:nrow(models_kruskal1)){
  
  # grab species
  temp_gene <- models_kruskal1$Gene[i]
  
  # filter just for that amr class then the columns we need
  temp_data <- filter(gathered_merged_evolved, Gene==temp_gene) %>%
    select(Concentration, Relative_Abundance)
  
  # run kruskal wallis test
  temp_model <- kruskal.test(Relative_Abundance ~ Concentration, temp_data)
  
  # get tidy output
  temp_output <- tidy(temp_model)
  
  # assign each output to the correct section in the results dataframe
  models_kruskal1$data[[i]] <- temp_data
  models_kruskal1$model[[i]] <- temp_model
  models_kruskal1$summary[[i]] <- temp_output
  
}

models_kruskal1


# look at one of the models
models_kruskal1$model[[1]]

# a function called unnest() allows us to access different bits of the results dataframe

# 1. get the data frame again
unnest(models_kruskal1, data) %>% select(-model, -summary)

# 2. get the output of the model
model_output1 <- unnest(models_kruskal1, summary) %>% select(-model, -data) 

# add in a column for p adjustments
model_output1 <- mutate(model_output1, padj = p.adjust(p.value, method = 'fdr'))

filter(model_output1, padj < 0.05) # here are our significant ones
#= no significant results = no genes significantly change with treatment when we look at them all




#### Filter for Efflux Pump Genes ####

#Read in list of efflux pump genes gathered from BacMet
efflux_list <- read.csv("Efflux Pump Genes.csv", header=TRUE)

#Make this a vector
efflux_vector <- efflux_list$Gene

#Use a for loop to filter the initial subset (gathered_merged) by the efflux pump genes
#We are interested in all genes, but really the evolved population, so do with all of them, and then filter by concentration != inoc

for (i in efflux_vector){
  if(!exists("df_efflux")){
    df_efflux <- d_avg[d_avg$Gene %like% i, ]
  }
  #if the merged dataset exists append to it
  if(exists("df_efflux")){
    temp_dataset <- d_avg[d_avg$Gene %like% i, ]
    df_efflux <- rbind(df_efflux, temp_dataset)
    rm(temp_dataset)
  }
}


#### Stats on Efflux Pump Genes ####

#Need to run the for loop again so we can subset by the genes in gathered_merged

for (i in efflux_vector){
  if(!exists("df_efflux_stats")){
    df_efflux_stats <- gathered_merged[gathered_merged$Gene %like% i, ]
  }
  #if the merged dataset exists append to it
  if(exists("df_efflux_stats")){
    temp_dataset <- gathered_merged[gathered_merged$Gene %like% i, ]
    df_efflux_stats <- rbind(df_efflux_stats, temp_dataset)
    rm(temp_dataset)
  }
}

## Do in the evolved populations
#Dataset for stats with just the evolved populations
df_efflux_stats_evolved <- subset(df_efflux_stats, Concentration!="inoc")


#USe the larger for loop for stats

# set up empty results dataframe - to fill with the models we create (using antibiotic.class as AMR genes)
models_kruskal_efflux1 <- select(df_efflux_stats_evolved, Gene) %>%
  distinct() %>%
  mutate(data = list(NA),
         model = list(NA),
         summary = list(NA))

for(i in 1:nrow(models_kruskal_efflux1)){
  
  # grab species
  temp_gene <- models_kruskal_efflux1$Gene[i]
  
  # filter just for that amr class then the columns we need
  temp_data <- filter(df_efflux_stats_evolved, Gene==temp_gene) %>%
    select(Concentration, Relative_Abundance)
  
  # run kruskal wallis test
  temp_model <- kruskal.test(Relative_Abundance ~ Concentration, temp_data)
  
  # get tidy output
  temp_output <- tidy(temp_model)
  
  # assign each output to the correct section in the results dataframe
  models_kruskal_efflux1$data[[i]] <- temp_data
  models_kruskal_efflux1$model[[i]] <- temp_model
  models_kruskal_efflux1$summary[[i]] <- temp_output
  
}

models_kruskal_efflux1


# look at one of the models
models_kruskal_efflux1$model[[1]]

# a function called unnest() allows us to access different bits of the results dataframe

# 1. get the data frame again
unnest(models_kruskal_efflux1, data) %>% select(-model, -summary)

# 2. get the output of the model
model_output_efflux1 <- unnest(models_kruskal_efflux1, summary) %>% select(-model, -data) 

# add in a column for p adjustments
model_output_efflux1 <- mutate(model_output_efflux1, padj = p.adjust(p.value, method = 'fdr'))

filter(model_output_efflux1, padj < 0.05) # here are our significant ones




#### Plasmid Genes ####

#Read in the list of plasmid genes
plasmid_list <- read.csv("Plasmid Related Genes.csv", header=TRUE)

#Make this a vector
plasmid_vector <- plasmid_list$Gene


#Need to run the for loop again so we can subset by the genes in gathered_merged

for (i in plasmid_vector){
  if(!exists("df_plasmid_stats")){
    df_plasmid_stats <- gathered_merged[gathered_merged$Gene %like% i, ]
  }
  #if the merged dataset exists append to it
  if(exists("df_plasmid_stats")){
    temp_dataset <- gathered_merged[gathered_merged$Gene %like% i, ]
    df_plasmid_stats <- rbind(df_plasmid_stats, temp_dataset)
    rm(temp_dataset)
  }
}

#Dataset for stats with just the evolved populations
df_plasmid_stats_evolved <- subset(df_plasmid_stats, Concentration!="inoc")


#USe the larger for loop for stats

# set up empty results dataframe - to fill with the models we create (using antibiotic.class as AMR genes)
models_kruskal_plasmid1 <- select(df_plasmid_stats_evolved, Gene) %>%
  distinct() %>%
  mutate(data = list(NA),
         model = list(NA),
         summary = list(NA))

for(i in 1:nrow(models_kruskal_plasmid1)){
  
  # grab species
  temp_gene <- models_kruskal_plasmid1$Gene[i]
  
  # filter just for that amr class then the columns we need
  temp_data <- filter(df_plasmid_stats_evolved, Gene==temp_gene) %>%
    select(Concentration, Relative_Abundance)
  
  # run kruskal wallis test
  temp_model <- kruskal.test(Relative_Abundance ~ Concentration, temp_data)
  
  # get tidy output
  temp_output <- tidy(temp_model)
  
  # assign each output to the correct section in the results dataframe
  models_kruskal_plasmid1$data[[i]] <- temp_data
  models_kruskal_plasmid1$model[[i]] <- temp_model
  models_kruskal_plasmid1$summary[[i]] <- temp_output
  
}

models_kruskal_plasmid1


# look at one of the models
models_kruskal_plasmid1$model[[1]]

# a function called unnest() allows us to access different bits of the results dataframe

# 1. get the data frame again
unnest(models_kruskal_plasmid1, data) %>% select(-model, -summary)

# 2. get the output of the model
model_output_plasmid1 <- unnest(models_kruskal_plasmid1, summary) %>% select(-model, -data) 

# add in a column for p adjustments
model_output_plasmid1 <- mutate(model_output_plasmid1, padj = p.adjust(p.value, method = 'fdr'))

filter(model_output_plasmid1, padj < 0.05) # here are our significant ones


### Make a Df of the genes of interest
evolved_genes_plasmidinterest <- subset(model_output_plasmid1, padj<0.05) #Genes significantly change in one treatment prior to adjustment


#Write it out as a CSV so can copy out as a table
write.csv(evolved_genes_plasmidinterest, "List of Plasmid associated genes sig diff after adjustment KW.csv", row.names=FALSE)

#make a list of the genes of interest
list_genes_plasmid <- evolved_genes_plasmidinterest$Gene



#Subset gathered_merged_evolved to only those genes in list genes
#This is so we have all of the information about these genes in one dataframe
#This is only the evolved populations since the stats was only on the evolved pops
d_plot_plasmid <- filter(gathered_merged_evolved, Gene %in% list_genes_plasmid)


#Make sure the data are continuous variables and make sure that the factors are ordered
d_plot_plasmid$Concentration <- factor(d_plot_plasmid$Concentration, levels=c("0", "2", "7", "22", "67", "201", "607", "1800", "5400" ))


#This enables us to plot this properly
d_plot_plasmid$Concentration <- as.numeric(as.character(d_plot_plasmid$Concentration))



## Do correlation For loop to check the correlation between dose and BMRG prevalence ###
df_plasmid_stats_coruse <- d_plot_plasmid

df_plasmid_stats_coruse$Concentration <- as.numeric(as.character(df_plasmid_stats_coruse$Concentration))

models_kruskal_plasmid_cor <- select(df_plasmid_stats_coruse, Gene) %>%
  distinct() %>%
  mutate(data = list(NA),
         model = list(NA),
         summary = list(NA))

for(i in 1:nrow(models_kruskal_plasmid_cor)){
  
  # grab species
  temp_gene <- models_kruskal_plasmid_cor$Gene[i]
  
  # filter just for that amr class then the columns we need
  temp_data <- filter(df_plasmid_stats_coruse, Gene==temp_gene) %>%
    select(Concentration, Relative_Abundance)
  
  # run correlation test
  temp_model <- cor.test(temp_data$Relative_Abundance, temp_data$Concentration)
  
  # get tidy output
  temp_output <- tidy(temp_model)
  
  # assign each output to the correct section in the results dataframe
  models_kruskal_plasmid_cor$data[[i]] <- temp_data
  models_kruskal_plasmid_cor$model[[i]] <- temp_model
  models_kruskal_plasmid_cor$summary[[i]] <- temp_output
  
}

models_kruskal_plasmid_cor


# look at one of the models
models_kruskal_plasmid_cor$model[[1]]

# a function called unnest() allows us to access different bits of the results dataframe

# 1. get the data frame again
unnest(models_kruskal_plasmid_cor, data) %>% select(-model, -summary)

# 2. get the output of the model
model_output_plasmidcor <- unnest(models_kruskal_plasmid_cor, summary) %>% select(-model, -data) 

# add in a column for p adjustments
model_output_plasmidcor <- mutate(model_output_plasmidcor, padj = p.adjust(p.value, method = 'fdr'))

filter(model_output_plasmidcor, padj < 0.05) # here are our significant ones

#Remove Gene name stuff 
d_plot_plasmid$Gene<- gsub("BAC*|", "", d_plot_plasmid$Gene)
d_plot_plasmid$Gene <- gsub("0*|", "", d_plot_plasmid$Gene)


## Look at the genes that showed a significant dose response and were significant different in at least one treatment 
## also were significantly correlated prior to correction



############################### Genes dose response ##############################

## arsB

d_arsB <- subset(d_plot_plasmid, Gene=="576|arsB|sp|P52146|ARSB2_ECOLX")

cor.test(d_arsB$Relative_Abundance, d_arsB$Concentration)

model_arsB <- lm(Relative_Abundance ~ Concentration, data=d_arsB)

summary(model_arsB)

simulationOutput <- simulateResiduals(fittedModel = model_arsB, plot = T)

#Nice fit of the data

anova(model_arsB, test="F")

preds_df_arsB <- group_by(d_arsB, Concentration) %>%
  do(data.frame(expand.grid(Concentration = seq(0, 5400, length.out = 100)
  ))) %>%
  ungroup()

preds_df_arsB <- cbind(preds_df_arsB, predict(model_arsB, newdata = preds_df_arsB, interval = 'confidence'))

df_arsB_avg <- d_arsB %>%
  group_by(Concentration) %>%
  summarise(
    mean=mean(Relative_Abundance), 
    sd=sd(Relative_Abundance),
    se=std.error(Relative_Abundance)) %>%
  mutate(Relative_Abundance=sd)

df_arsB_avg$Concentration <- as.numeric(df_arsB_avg$Concentration)


plot_arsB <- ggplot() +
  theme_bw()+
  geom_ribbon(aes(ymin=lwr, ymax=upr, x=Concentration), fill = "#BA553C", preds_df_arsB, alpha=0.1) +
  geom_line(aes(x=Concentration, y=fit), color = "#913119", size=1.22, preds_df_arsB)+
  geom_point(data=d_arsB, mapping = aes(x=Concentration, y=Relative_Abundance), colour = "#BA553C", size=4, alpha=0.5)+
  geom_point(data=df_arsB_avg, mapping = aes(x=Concentration, y=mean), fill = "#913119",  size=6, shape=21)+
  labs(y="Prevalence of arsR", x="17-β-estradiol concentration (µg/L)")+
  theme_classic()+
  theme(axis.text=element_text(size=16))+ 
  theme(axis.title=element_text(size = 18))+
  theme(legend.text = element_text(size = 14)) + 
  theme(legend.title = element_text(size = 16))+
  ylab(expression(paste("Prevalence of", italic("arsB"))))



# arsB v2

d_arsB_v2 <- subset(d_plot_plasmid, Gene =="578|arsB|tr|O5594|O5594_ACIMU")

cor.test(d_arsB_v2$Relative_Abundance, d_arsB_v2$Concentration)
#Significant COrrelation, p=0.03266

m_arsB2 <- lm(Relative_Abundance ~ Concentration, data=d_arsB_v2)

summary(m_arsB2)

simulationOutput <- simulateResiduals(fittedModel = m_arsB2, plot = T)

anova(m_arsB2, test="F")

preds_df_arsBv2 <- group_by(d_arsB_v2, Concentration) %>%
  do(data.frame(expand.grid(Concentration = seq(0, 5400, length.out = 100)
  ))) %>%
  ungroup()

preds_df_arsBv2 <- cbind(preds_df_arsBv2, predict(m_arsB2, newdata = preds_df_arsBv2, interval = 'confidence'))


df_arsBv2_avg <- d_arsB_v2 %>%
  group_by(Concentration) %>%
  summarise(
    mean=mean(Relative_Abundance), 
    sd=sd(Relative_Abundance),
    se=std.error(Relative_Abundance)) %>%
  mutate(Relative_Abundance=sd)

df_arsBv2_avg$Concentration <- as.numeric(df_arsBv2_avg$Concentration)


plot_arsB_2 <- ggplot() +
  theme_bw()+
  geom_ribbon(aes(ymin=lwr, ymax=upr, x=Concentration), fill = "#BA553C", preds_df_arsBv2, alpha=0.1) +
  geom_line(aes(x=Concentration, y=fit), color = "#913119", size=1.22, preds_df_arsBv2)+
  geom_point(data=d_arsB_v2, mapping = aes(x=Concentration, y=Relative_Abundance), colour = "#BA553C", size=4, alpha=0.5)+
  geom_point(data=df_arsBv2_avg, mapping = aes(x=Concentration, y=mean), fill = "#913119",  size=6, shape=21)+
  labs(y="Prevalence of arsR", x="17-β-estradiol concentration (µg/L)")+
  theme_classic()+
  theme(axis.text=element_text(size=16))+ 
  theme(axis.title=element_text(size = 18))+
  theme(legend.text = element_text(size = 14)) + 
  theme(legend.title = element_text(size = 16))+
  ylab(expression(paste("Prevalence of", italic("arsB"))))




### arsR
d_594 <- subset(d_plot_plasmid, Gene == "594|arsR|sp|P3739|ARSR_ECOLI")

cor.test(d_594$Relative_Abundance, d_594$Concentration)

m_arsR <- lm(Relative_Abundance ~ Concentration, data=d_594)

summary(m_arsR)

simulationOutput <- simulateResiduals(fittedModel = m_arsR, plot = T)

anova(m_arsR, test="F")

preds_df_594 <- group_by(d_594, Concentration) %>%
  do(data.frame(expand.grid(Concentration = seq(0, 5400, length.out = 100)
  ))) %>%
  ungroup()

preds_df_594 <- cbind(preds_df_594, predict(m_arsR, newdata = preds_df_594, interval = 'confidence'))


df_594_avg <- d_594 %>%
  group_by(Concentration) %>%
  summarise(
    mean=mean(Relative_Abundance), 
    sd=sd(Relative_Abundance),
    se=std.error(Relative_Abundance)) %>%
  mutate(Relative_Abundance=sd)

df_594_avg$Concentration <- as.numeric(df_594_avg$Concentration)

plot_arsR <- ggplot() +
  theme_bw()+
  geom_ribbon(aes(ymin=lwr, ymax=upr, x=Concentration), fill = "#BA553C", preds_df_594, alpha=0.1) +
  geom_line(aes(x=Concentration, y=fit), color = "#913119", size=1.22, preds_df_594)+
  geom_point(data=d_594, mapping = aes(x=Concentration, y=Relative_Abundance), colour = "#BA553C", size=4, alpha=0.5)+
  geom_point(data=df_594_avg, mapping = aes(x=Concentration, y=mean), fill = "#913119",  size=6, shape=21)+
  labs(y="Prevalence of arsR", x="17-β-estradiol concentration (µg/L)")+
  theme_classic()+
  theme(axis.text=element_text(size=16))+ 
  theme(axis.title=element_text(size = 18))+
  theme(legend.text = element_text(size = 14)) + 
  theme(legend.title = element_text(size = 16))+
  ylab(expression(paste("Prevalence of", italic("arsR"))))



### ncrA
d_72 <- subset(d_plot_plasmid, Gene == "72|ncrA|tr|Q8VTR6|Q8VTR6_HAFAL")

cor.test(d_72$Relative_Abundance, d_72$Concentration)

m_ncrA <- lm(sqrt(Relative_Abundance) ~ Concentration, data=d_72)

simulateResiduals(fittedModel = m_ncrA, plot = T)
hist(resid(m_ncrA))

qqnorm(resid(m_ncrA))
qqline(resid(m_ncrA))

summary(m_ncrA)

anova(m_ncrA, test="F")


#Make the plot
preds_df_72 <- group_by(d_72, Concentration) %>%
  do(data.frame(expand.grid(Concentration = seq(0, 5400, length.out = 100)
  ))) %>%
  ungroup()

preds_df_72 <- cbind(preds_df_72, predict(m_ncrA, newdata = preds_df_72, interval = 'confidence'))


df_72_avg <- d_72 %>%
  group_by(Concentration) %>%
  summarise(
    mean=mean(Relative_Abundance), 
    sd=sd(Relative_Abundance),
    se=std.error(Relative_Abundance)) %>%
  mutate(Relative_Abundance=sd)

df_72_avg$Concentration <- as.numeric(df_72_avg$Concentration)



#Graph for ncrA
plot_ncrA <- ggplot() +
  theme_bw()+
  geom_ribbon(aes(ymin=lwr, ymax=upr, x=Concentration), fill = "#BA553C", preds_df_72, alpha=0.1) +
  geom_line(aes(x=Concentration, y=(fit)), color = "#913119", size=1.22, preds_df_72)+
  geom_point(data=d_72, mapping = aes(x=Concentration, y=sqrt(Relative_Abundance)), colour = "#BA553C", size=4, alpha=0.5)+
  geom_point(data=df_72_avg, mapping = aes(x=Concentration, y=sqrt(mean)), fill = "#913119",  size=6, shape=21)+
  labs(y="Prevalence of arsR", x="17-β-estradiol concentration (µg/L)")+
  theme_classic()+
  theme(axis.text=element_text(size=18))+ 
  theme(axis.title=element_text(size = 20))+
  theme(legend.text = element_text(size = 14)) + 
  theme(legend.title = element_text(size = 16))+
  ylab(expression(paste("Prevalence of ", italic("ncrA"))))

### plot together

plot_arsB + plot_arsB_2 + plot_arsR + plot_ncrA

ggsave("figures/Figure 4.tiff", height = 8, width = 6, dpi = 300)