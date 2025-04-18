##################### BacMet Analysis  for Metformin ##########################

#################### Read in code and libraries ###############################
library(tidyverse) 
library(reshape2) 
library(ggplot2)
library(MetBrewer) 
library(vegan) 
library(data.table) 
library(broom) 



## Need this later - added in here with admin things ##
#create a function to split each string by an identifier and grab the second index
quick_strsplit <- function(x, index, split){
  # split the string
  temp <- strsplit(x, split = split)
  # iterate through the list and grab each indexed element
  temp <- sapply(temp, "[", index)
  return(temp)
}


######## #Read in the files
# list files
files <- list.files(path="data", full.names=TRUE)

# Make a function to read in files
aprils_cleaning_function <- function(file){
  
  # read in file
  df <- read_tsv(file, col_names=TRUE) %>%
    subset(., `Percent identity`>="85") %>%
    subset(., `Match length`>="25") %>% #these marks make us read the thing as a whole, ot as two separate words
    group_by(Gene) %>%
    tally() %>%
    mutate(., file = basename(tools::file_path_sans_ext(file)))
  
  return(df)
  
}

# make a dataframe
df <- purrr::map_df(files, aprils_cleaning_function)

df$file <- gsub("bacmet_", "", df$file)

df$file[df$file=="metconcAr1"] <- "3300 1"
df$file[df$file=="metconcAr2"] <- "3300 2"
df$file[df$file=="metconcAr3"] <- "3300 3"
df$file[df$file=="metconcBr1"] <- "1100 1"
df$file[df$file=="metconcBr2"] <- "1100 2"
df$file[df$file=="metconcBr3"] <- "1100 3"
df$file[df$file=="metconcCr1"] <- "367 1"
df$file[df$file=="metconcCr2"] <- "367 2"
df$file[df$file=="metconcCr3"] <- "367 3"
df$file[df$file=="metconcDr1"] <- "122 1"
df$file[df$file=="metconcDr2"] <- "122 2"
df$file[df$file=="metconcDr3"] <- "122 3"
df$file[df$file=="metconcEr1"] <- "40.7 1"
df$file[df$file=="metconcEr2"] <- "40.7 2"
df$file[df$file=="metconcEr3"] <- "40.7 3"
df$file[df$file=="metconcFr1"] <- "13.6 1"
df$file[df$file=="metconcFr2"] <- "13.6 2"
df$file[df$file=="metconcFr3"] <- "13.6 3"
df$file[df$file=="metconcGr1"] <- "4.5 1"
df$file[df$file=="metconcGr2"] <- "4.5 2"
df$file[df$file=="metconcGr3"] <- "4.5 3"
df$file[df$file=="metconcHr1"] <- "0 1"
df$file[df$file=="metconcHr2"] <- "0 2"
df$file[df$file=="metconcHr3"] <- "0 3"
df$file[df$file=="bacmetinocrep1"] <- "inoc 1"
df$file[df$file=="bacmetinocrep2"] <- "inoc 2"
df$file[df$file=="bacmetinocrep3"] <- "inoc 3"

colnames(df)[3] <- "Name"


########### Read in 16S Data ###

cells <- read.table("16S Data.txt", header=TRUE, sep="\t")

cells <- cells %>% select(-SampleID)

cells <- cells %>% select(Name, X16Sreads, CellNumber)

cells$Name <- gsub("dicloconc", "", cells$Name)

cells$Name[cells$Name=="metconcAr1"] <- "3300 1"
cells$Name[cells$Name=="metconcAr2"] <- "3300 2"
cells$Name[cells$Name=="metconcAr3"] <- "3300 3"
cells$Name[cells$Name=="metconcBr1"] <- "1100 1"
cells$Name[cells$Name=="metconcBr2"] <- "1100 2"
cells$Name[cells$Name=="metconcBr3"] <- "1100 3"
cells$Name[cells$Name=="metconcCr1"] <- "367 1"
cells$Name[cells$Name=="metconcCr2"] <- "367 2"
cells$Name[cells$Name=="metconcCr3"] <- "367 3"
cells$Name[cells$Name=="metconcDr1"] <- "122 1"
cells$Name[cells$Name=="metconcDr2"] <- "122 2"
cells$Name[cells$Name=="metconcDr3"] <- "122 3"
cells$Name[cells$Name=="metconcEr1"] <- "40.7 1"
cells$Name[cells$Name=="metconcEr2"] <- "40.7 2"
cells$Name[cells$Name=="metconcEr3"] <- "40.7 3"
cells$Name[cells$Name=="metconcFr1"] <- "13.6 1"
cells$Name[cells$Name=="metconcFr2"] <- "13.6 2"
cells$Name[cells$Name=="metconcFr3"] <- "13.6 3"
cells$Name[cells$Name=="metconcGr1"] <- "4.5 1"
cells$Name[cells$Name=="metconcGr2"] <- "4.5 2"
cells$Name[cells$Name=="metconcGr3"] <- "4.5 3"
cells$Name[cells$Name=="metconcHr1"] <- "0 1"
cells$Name[cells$Name=="metconcHr2"] <- "0 2"
cells$Name[cells$Name=="metconcHr3"] <- "0 3"
cells$Name[cells$Name=="inocrep1"] <- "inoc 1"
cells$Name[cells$Name=="inocrep2"] <- "inoc 2"
cells$Name[cells$Name=="inocrep3"] <- "inoc 3"

#Join the two dataframes together
merged <- left_join(df, cells)

#Separate Sample into Conc and Rep
merged <- tidyr::separate(merged, Name, c("Concentration", "Replicate"), sep=" ")

#Make relative abundance of gene (i.e. bacmet gene per 16S)
merged <- merged %>% mutate(Gene_per_16S = merged$n / merged$X16Sreads)

#Average the dataframe
d_avg <- group_by(merged, Concentration, Gene) %>% 
  summarise(Avg=mean(Gene_per_16S), .groups="drop")

#Remove some of the gross gene name business
d_avg$Gene <- gsub("BAC*|", "", d_avg$Gene)
d_avg$Gene <- gsub("0*|", "", d_avg$Gene)


#Add levels accurately
d_avg$Concentration <- factor(d_avg$Concentration, levels=unique(d_avg$Concentration))

d_avg$Concentration <- factor(d_avg$Concentration, levels=c("inoc", "0", "4.5", "13.6", "40.7", "122", "367", "1100", "3300"))

#Make evolved dataframe
d_avg_evolved <- subset(d_avg, Concentration!="inoc") 

#Set the levels
d_avg_evolved$Concentration <- factor(d_avg_evolved$Concentration, levels=unique(d_avg_evolved$Concentration))
d_avg_evolved$Concentration <- factor(d_avg_evolved$Concentration, levels=c("0", "4.5", "13.6", "40.7", "122", "367", "1100", "3300"))



#### Changing the dataframe so that we have zeros (instead of presence/absense) ####

working_merged <- merged

working_merged$concrep <- paste(working_merged$Concentration, working_merged$Replicate) 

working_merged <- select(working_merged, Gene, concrep, Gene_per_16S)

working_merged <- spread(working_merged, concrep, Gene_per_16S)

working_merged[is.na(working_merged)] <- 0

gathered_merged <- gather(working_merged, key="Concentration", "Relative_Abundance", 2:28) 

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
#= no significant results



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
d_plot_plasmid$Concentration <- factor(d_plot_plasmid$Concentration, levels=c("inoc", "0", "4.5", "13.6", "40.7", "122", "367", "1100", "3300"))

d_plot_plasmid_evolved <- subset(d_plot_plasmid, Concentration!="inoc")

d_plot_plasmid_evolved$Concentration <- as.numeric(as.character(d_plot_plasmid_evolved$Concentration))



## Do correlation For loop to check the correlation between dose and BMRG prevalence ###

models_kruskal_plasmid_cor <- select(d_plot_plasmid_evolved, Gene) %>%
  distinct() %>%
  mutate(data = list(NA),
         model = list(NA),
         summary = list(NA))

for(i in 1:nrow(models_kruskal_plasmid_cor)){
  
  # grab species
  temp_gene <- models_kruskal_plasmid_cor$Gene[i]
  
  # filter just for that amr class then the columns we need
  temp_data <- filter(d_plot_plasmid_evolved, Gene==temp_gene) %>%
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
output_models_plasmid_cor <- unnest(models_kruskal_plasmid_cor, summary) %>% select(-model, -data) 

# add in a column for p adjustments
output_models_plasmid_cor <- mutate(output_models_plasmid_cor, padj = p.adjust(p.value, method = 'fdr'))

filter(output_models_plasmid_cor, padj < 0.05) 
#no genes after adjustment are significant



#### NMDS on evolved populations only ####
d_nmds <- gathered_merged_evolved
# d_nmds <- subset(gathered_merged, select = c(1:3, 7))
d_nmds <- unite(d_nmds, Sample, Concentration, Rep, sep = " ")

d_nmds <- reshape2::dcast(d_nmds, Sample~Gene, sum)

#Make Treatment and Rep Columns
d_nmds <- separate(d_nmds, Sample, c("Treatment", "Rep"), sep=" ")

d_nmds_evolved <- subset(d_nmds, Treatment!="inoc")

#Do the NMDS 
pc = d_nmds_evolved

#make community matrix - extract columns with abundance information
com = pc[ ,3:ncol(pc)]
m_com = as.matrix(com)


#Set seed gives you the same results each time you run this code
set.seed(123)
nmds = metaMDS(m_com, distance = "bray")
nmds

ano = anosim(m_com, pc$Treatment, distance = "bray", permutations = 9999)
ano


