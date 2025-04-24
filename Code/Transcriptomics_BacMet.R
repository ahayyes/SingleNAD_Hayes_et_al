################### BacMet Analysis for Transcriptomics Data #####################

################# Read in code and libraries ####################################

library(tidyverse) 
library(reshape2) 
library(ggplot2)
library(MetBrewer) 
library(vegan)
library(data.table) 
library(broom)
library(plotrix)
library(patchwork)
library(rstatix)
library(ggpubr)
library(dunn.test)
library(patchwork)

############################## Read in Files #################################

#Read in the filenames
# list files
files <- list.files(path="data", full.names=TRUE)

# Make a function to read and make a df

aprils_cleaning_function <- function(file){
  
  # read in file
  df <- read_tsv(file, col_names=TRUE) %>%
    subset(., `Percent identity`>="85") %>%
    subset(., `Match length`>="25") %>%
    group_by(Gene) %>%
    tally() %>%
    mutate(., file = basename(tools::file_path_sans_ext(file)))
  
  return(df)
  
}

# Make the df
df <- purrr::map_df(files, aprils_cleaning_function)

# Remove bacmet from the files
df$file <- gsub("bacmet_", "", df$file)


#Change names 
df$file[df$file=="C17repA"] <- "17-beta 1"
df$file[df$file=="C17repB"] <- "17-beta 2"
df$file[df$file=="C17repC"] <- "17-beta 3"
df$file[df$file=="C17repD"] <- "17-beta 4"
df$file[df$file=="C17repE"] <- "17-beta 5"
df$file[df$file=="ControlrepA"] <- "Control 1"
df$file[df$file=="ControlrepB"] <- "Control 2"
df$file[df$file=="ControlrepC"] <- "Control 3"
df$file[df$file=="ControlrepCV2"] <- "Control 3V2"
df$file[df$file=="ControlrepD"] <- "Control 4"
df$file[df$file=="ControlrepE"] <- "Control 5"
df$file[df$file=="DrepA"] <- "Diclofenac 1"
df$file[df$file=="DrepB"] <- "Diclofenac 2"
df$file[df$file=="DrepBV2"] <- "Diclofenac 2V2"
df$file[df$file=="DrepC"] <- "Diclofenac 3"
df$file[df$file=="DrepD"] <- "Diclofenac 4"
df$file[df$file=="DrepE"] <- "Diclofenac 5"
df$file[df$file=="MrepA"] <- "Metformin 1"
df$file[df$file=="MrepAV2"] <- "Metformin 1V2"
df$file[df$file=="MrepB"] <- "Metformin 2"
df$file[df$file=="MrepBV2"] <- "Metformin 2V2"
df$file[df$file=="MrepC"] <- "Metformin 3"
df$file[df$file=="MrepD"] <- "Metformin 4"
df$file[df$file=="MrepE"] <- "Metformin 5"


colnames(df)[3] <- "Name"


######################## Read in 16S Data #######################################

cells <- read.csv("16S Numbers Excel.csv", header=TRUE)

#Sorting Dataframe to remove unnecessary columns
cells <- cells %>% select(Name, X.of16Sreads, CellNumber)

#Change names 
cells$Name[cells$Name=="C17repA"] <- "17-beta 1"
cells$Name[cells$Name=="C17repB"] <- "17-beta 2"
cells$Name[cells$Name=="C17repC"] <- "17-beta 3"
cells$Name[cells$Name=="C17repD"] <- "17-beta 4"
cells$Name[cells$Name=="C17repE"] <- "17-beta 5"
cells$Name[cells$Name=="ControlrepA"] <- "Control 1"
cells$Name[cells$Name=="ControlrepB"] <- "Control 2"
cells$Name[cells$Name=="ControlrepC"] <- "Control 3"
cells$Name[cells$Name=="ControlrepCV2"] <- "Control 3V2"
cells$Name[cells$Name=="ControlrepD"] <- "Control 4"
cells$Name[cells$Name=="ControlrepE"] <- "Control 5"
cells$Name[cells$Name=="DrepA"] <- "Diclofenac 1"
cells$Name[cells$Name=="DrepB"] <- "Diclofenac 2"
cells$Name[cells$Name=="DrepBV2"] <- "Diclofenac 2V2"
cells$Name[cells$Name=="DrepC"] <- "Diclofenac 3"
cells$Name[cells$Name=="DrepD"] <- "Diclofenac 4"
cells$Name[cells$Name=="DrepE"] <- "Diclofenac 5"
cells$Name[cells$Name=="MrepA"] <- "Metformin 1"
cells$Name[cells$Name=="MrepAV2"] <- "Metformin 1V2"
cells$Name[cells$Name=="MrepB"] <- "Metformin 2"
cells$Name[cells$Name=="MrepBV2"] <- "Metformin 2V2"
cells$Name[cells$Name=="MrepC"] <- "Metformin 3"
cells$Name[cells$Name=="MrepD"] <- "Metformin 4"
cells$Name[cells$Name=="MrepE"] <- "Metformin 5"


#Need to remove the doubles from the dataframe here - spread to make rows columns, 
#Then add the columns together

# Merge the two dataframes so you can have relative abundance of genes
merged <- merge(df, cells, by = "Name")

#Need to average the repeated ones
#Metformin 1
#Metformin 2
#Control 3

merged <- merged %>% mutate(Gene_per_16S = merged$n / merged$X.of16Sreads) 
#This provides us with a dataframe that has relative abundance of BMRGs
#We have relative abundance
#Now we can remove the duplicates

#Make a new dataframe to deal with this
rm_dups <- merged

rm_dups <- select(rm_dups, Name, Gene, Gene_per_16S )

rm_dups <- spread(rm_dups, key=Name, value=Gene_per_16S)

#Changes NA to 0
rm_dups[is.na(rm_dups)] <- 0

# Remove the dups
rm_dups <- rm_dups %>% mutate("Control 3"=c(`Control 3`+ `Control 3V2`)) %>% select(-`Control 3V2`)

rm_dups <- rm_dups %>% mutate("Metformin 1"=c(`Metformin 1`+ `Metformin 1V2`)) %>% select(-`Metformin 1V2`)

rm_dups <- rm_dups %>% mutate("Metformin 2"=c(`Metformin 2`+ `Metformin 2V2`)) %>% select(-`Metformin 2V2`)

colnames(rm_dups)[colnames(rm_dups) == 'Diclofenac 2V2'] <- 'Diclofenac 2'

gathered_merged <- gather(rm_dups, key="Concentration", "Relative_Abundance", 2:21)

gathered_merged$concrep <- paste(gathered_merged$Treatment, gathered_merged$Replicate)

#Average
d_avg <- group_by(gathered_merged, Treatment, Gene) %>% 
  summarise(Avg=mean(Relative_Abundance), .groups="drop")

#Set the levels
d_avg$Treatment <- factor(d_avg$Treatment, levels=unique(d_avg$Treatment))
d_avg$Treatment <- factor(d_avg$Treatment, levels=c("Control", "Diclofenac", "Metformin", "17-beta"))


#################### Stats on everything #################################

# set up empty results dataframe -

models_kruskal1 <- select(gathered_merged, Gene) %>%
  distinct() %>%
  mutate(data = list(NA),
         model = list(NA),
         summary = list(NA))

for(i in 1:nrow(models_kruskal1)){
  
  # grab species
  temp_gene <- models_kruskal1$Gene[i]
  
  # filter just for that amr class then the columns we need
  temp_data <- filter(gathered_merged, Gene==temp_gene) %>%
    select(Treatment, Relative_Abundance)
  
  # run kruskal wallis test
  temp_model <- kruskal.test(Relative_Abundance ~ Treatment, temp_data)
  
  # get tidy output
  temp_output <- tidy(temp_model)
  
  # assign each output to the correct section in the results dataframe
  models_kruskal1$data[[i]] <- temp_data
  models_kruskal1$model[[i]] <- temp_model
  models_kruskal1$summary[[i]] <- temp_output
  
}

# models_kruskal1

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
#= no significant results

#### Filter for Efflux Pump Genes ####

#Read in list of efflux pump genes gathered from BacMet
efflux_list <- read.csv("Efflux Pump Genes.csv", header=TRUE)

#Make this a vector
efflux_vector <- efflux_list$Gene

#Use a for loop to filter the initial subset (gathered_merged) by the efflux pump genes
#We are interested in all genes, but really the evolved population, so do with all of them, and then filter by concentration != inoc

for (i in efflux_vector){
  if(!exists("df_efflux")){
    df_efflux <- gathered_merged[gathered_merged$Gene %like% i, ]
  }
  #if the merged dataset exists append to it
  if(exists("df_efflux")){
    temp_dataset <- gathered_merged[gathered_merged$Gene %like% i, ]
    df_efflux <- rbind(df_efflux, temp_dataset)
    rm(temp_dataset)
  }
}


# create a function to split each string by an identifier and grab the second index
quick_strsplit <- function(x, index, split){
  # split the string
  temp <- strsplit(x, split = split)
  # iterate through the list and grab each indexed element
  temp <- sapply(temp, "[", index)
  return(temp)
}


df_efflux <- df_efflux %>% mutate(gene_name = quick_strsplit(df_efflux$Gene, index = 2, split = "\\|"))

df_efflux$Treatment <- factor(df_efflux$Treatment, levels=unique(df_efflux$Treatment))
df_efflux$Treatment <- factor(df_efflux$Treatment, levels=c("Control", "Diclofenac", "Metformin", "17-beta" ))

df_efflux$concrep <- factor(df_efflux$concrep, levels=c("Control 1", "Control 2", "Control 3", "Control 4", "Control 5",
                                                        "Diclofenac 1", "Diclofenac 2", "Diclofenac 3", "Diclofenac 4", 
                                                        "Diclofenac 5", 
                                                        "Metformin 1", "Metformin 2", "Metformin 3", 
                                                        "Metformin 4", "Metformin 5", 
                                                        "17-beta 1",  "17-beta 2",  "17-beta 3",  "17-beta 4",  "17-beta 5"))


########################## Stats on Efflux Genes ###################################

# set up empty results dataframe - to fill with the models we create (using antibiotic.class as AMR genes)
models_kruskal_efflux <- select(df_efflux, Gene) %>%
  distinct() %>%
  mutate(data = list(NA),
         model = list(NA),
         summary = list(NA))


for(i in 1:nrow(models_kruskal_efflux)){
  
  # grab species
  temp_gene <- models_kruskal_efflux$Gene[i]
  
  # filter just for that amr class then the columns we need
  temp_data <- filter(df_efflux, Gene==temp_gene) %>%
    select(Treatment, Relative_Abundance)
  
  # run kruskal wallis test
  temp_model <- kruskal.test(Relative_Abundance ~ Treatment, temp_data)
  
  # get tidy output
  temp_output <- tidy(temp_model)
  
  # assign each output to the correct section in the results dataframe
  models_kruskal_efflux$data[[i]] <- temp_data
  models_kruskal_efflux$model[[i]] <- temp_model
  models_kruskal_efflux$summary[[i]] <- temp_output
  
}


# look at one of the models
models_kruskal_efflux$model[[1]]

# a function called unnest() allows us to access different bits of the results dataframe

# 1. get the data frame again
unnest(models_kruskal_efflux, data) %>% select(-model, -summary)

# 2. get the output of the model
model_output_efflux<- unnest(models_kruskal_efflux, summary) %>% select(-model, -data) 

# add in a column for p adjustments
model_output_efflux <- mutate(model_output_efflux, padj = p.adjust(p.value, method = 'fdr'))

filter(model_output_efflux, padj < 0.05) # here are our significant ones
#= no significant results



######################### Plasmid Genes #########################################

#Read in the list of plasmid genes
plasmid_list <- read.csv("Plasmid Related Genes.csv", header=TRUE)

#Make this a vector
plasmid_vector <- plasmid_list$Gene


#Use the for loop to create df with averaged replicates for heatmap for plasmid genes
for (i in plasmid_vector){
  if(!exists("df_plasmid")){
    df_plasmid <- gathered_merged[gathered_merged$Gene %like% i, ]
  }
  #if the merged dataset exists append to it
  if(exists("df_plasmid")){
    temp_dataset <- gathered_merged[gathered_merged$Gene %like% i, ]
    df_plasmid <- rbind(df_plasmid, temp_dataset)
    rm(temp_dataset)
  }
}



#Set levels of df_efflux
df_plasmid$Treatment <- factor(df_plasmid$Treatment, levels=unique(df_plasmid$Treatment))
df_plasmid$Treatment <- factor(df_plasmid$Treatment, levels=c("Control", "Diclofenac", "Metformin", "17-beta"))
df_plasmid$concrep <- factor(df_plasmid$concrep, levels=c("Control 1", "Control 2", "Control 3", "Control 4", "Control 5",
                                                          "Diclofenac 1", "Diclofenac 2", "Diclofenac 3", "Diclofenac 4", 
                                                          "Diclofenac 5", 
                                                          "Metformin 1", "Metformin 2", "Metformin 3", 
                                                          "Metformin 4", "Metformin 5", 
                                                          "17-beta 1",  "17-beta 2",  "17-beta 3",  "17-beta 4",  "17-beta 5"))

df_plasmid <- df_plasmid %>% mutate(gene_name = quick_strsplit(df_plasmid$Gene, index = 2, split = "\\|"))



# set up empty results dataframe - to fill with the models we create (using antibiotic.class as AMR genes)
models_kruskal_plasmid1 <- select(df_plasmid, Gene) %>%
  distinct() %>%
  mutate(data = list(NA),
         model = list(NA),
         summary = list(NA))

for(i in 1:nrow(models_kruskal_plasmid1)){
  
  # grab species
  temp_gene <- models_kruskal_plasmid1$Gene[i]
  
  # filter just for that amr class then the columns we need
  temp_data <- filter(df_plasmid, Gene==temp_gene) %>%
    select(Treatment, Relative_Abundance)
  
  # run kruskal wallis test
  temp_model <- kruskal.test(Relative_Abundance ~ Treatment, temp_data)
  
  # get tidy output
  temp_output <- tidy(temp_model)
  
  # assign each output to the correct section in the results dataframe
  models_kruskal_plasmid1$data[[i]] <- temp_data
  models_kruskal_plasmid1$model[[i]] <- temp_model
  models_kruskal_plasmid1$summary[[i]] <- temp_output
  
}


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

#Make a list of the plasmid genes that are of interest

list_genes_plasmid <- subset(model_output_plasmid1, padj<0.05)
list_genes_plasmid <- list_genes_plasmid$Gene


#Subset gathered_merged_evolved to only those genes in list genes
#Only includes evolved communities

d_plot_plasmid <- filter(gathered_merged, Gene %in% list_genes_plasmid )

#Make sure the data are continuous variables and make sure that the factors are ordered
d_plot_plasmid$Treatment <- factor(d_plot_plasmid$Treatment, levels=c("Control", "Diclofenac", "Metformin", "17-beta"))

#Shorten the gene name
d_plot_plasmid <- d_plot_plasmid %>% mutate(Gene_name = quick_strsplit(d_plot_plasmid$Gene, index = 2, split = "\\|"))

d_plot_plasmid <- d_plot_plasmid %>% mutate(Gen_num = quick_strsplit(d_plot_plasmid$Gene,
                                                                     index = 1, split = "\\|"))
d_plot_plasmid$Name <- paste(d_plot_plasmid$Gene_name, d_plot_plasmid$Gen_num, " ")


ggplot(d_plot_plasmid, aes(x=Treatment, y=Relative_Abundance, fill=Treatment))+
  geom_boxplot(outlier.shape =  NA, lwd = 0.8)+
  geom_point(position=position_jitterdodge(), alpha=0.7, size = 2, shape = 21, show.legend = FALSE)+
  theme_bw()+
  labs(x="NAD Treatment", y="Relative Abundance")+
  scale_fill_manual(values = c("#eace17", "#6dbc18", "#34b6c6", "#ed6a48"), name="Treatment",
                    labels=c("Control", "Diclofenac", "Metformin", "17-β-estradiol"))+  # theme(text=element_text(size=20))+
  theme(text=element_text(size=20))+
  theme(axis.title = element_text(size=20))+
  theme(strip.text = element_text(face = "italic"))+
  facet_wrap(~Name, scales = "free")+
  theme(strip.background=element_rect(colour="black",
                                      fill="white"))


#Filter to contain genes that are present in more than one sample 
#Since 6 genes here are only present in 1 sample 

#The genes we are interested in are 
#arsB BAC0031
#arsR BAC0589
#arsR BAC0594
#merA BAC0652
#merP BAC0231

#The genes we aren't are 
#arsB BAC0574
#arsB BAC0575
#merD BAC0227
#merT BAC0233
#merT BAC0690
#merT BAC0693

#Subset to get right of them 
d_plot_plasmid_true <- subset(d_plot_plasmid, Name != "arsB BAC0574  ") 
d_plot_plasmid_true <- subset(d_plot_plasmid_true, Gen_num != "BAC0575") 
d_plot_plasmid_true <- subset(d_plot_plasmid_true, Gen_num != "BAC0227") 
d_plot_plasmid_true <- subset(d_plot_plasmid_true, Gen_num != "BAC0233") 
d_plot_plasmid_true <- subset(d_plot_plasmid_true, Gen_num != "BAC0690") 
d_plot_plasmid_true <- subset(d_plot_plasmid_true, Gen_num != "BAC0693") 


#### Gene 1 - arsB BAC0031

df_arsB_BAC0031 <- subset(d_plot_plasmid_true, Gen_num =="BAC0031")

dunn.test(df_arsB_BAC0031$Relative_Abundance, df_arsB_BAC0031$Treatment)

#Create table for your stats from the stuff in the stats section
stat.test1 <- tibble::tribble(
  ~group1, ~group2, ~variable, ~p.adj, ~Name, ~Treatment,
  "Control", "Diclofenac", "Relative_Abundance", 'NS', "arsB BAC0031  ", "Control", 
  "Control", "Metformin", "Relative_Abundance", '*', "arsB BAC0031  ", "Control", 
  "Control", "17-beta", "Relative_Abundance", 'NS', "arsB BAC0031  ", "Control")

p1 <- ggplot(df_arsB_BAC0031, aes(x=Treatment, y=Relative_Abundance, fill=Treatment))+
  geom_boxplot(outlier.shape =  NA, lwd = 0.8)+
  geom_point(position=position_jitterdodge(), alpha=0.7, size = 2, shape = 21, show.legend = FALSE)+
  theme_bw()+
  labs(x="NAD Treatment", y="Relative Abundance")+
  scale_fill_manual(values = c("#ddc000", "#79ad41", "#34b6c6", "#d7aca1"), name="Treatment", 
                    labels=c("Control", "Diclofenac", "Metformin", "17-β-estradiol"))+
  # scale_fill_manual(values = c("#eace17", "#6dbc18", "#34b6c6", "#ed6a48"), name="Treatment",
  #                   labels=c("Control", "Diclofenac", "Metformin", "17-β-estradiol"))+
  theme(text=element_text(size=22))+
  scale_x_discrete(labels = c("Control", "Diclofenac", "Metformin", "17-β-estradiol"))+
  theme(axis.title = element_text(size=20))+
  theme(strip.text = element_text(face = "italic"))+
  theme(strip.background=element_rect(colour="black",
                                      fill="white"))+
  facet_wrap(~Name)+
  xlab(" ")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  stat_pvalue_manual(stat.test1, y.position = 0.032, label = "p.adj", xmin="group2", xmax=NULL, 
                     tip.length = 0.04, size = 5, bracket.size = 0.9)



#### Gene 2 - arsR BAC0589

df_arsR_BAC089 <- subset(d_plot_plasmid_true, Gen_num =="BAC0589")

dunn.test(df_arsR_BAC089$Relative_Abundance, df_arsR_BAC089$Treatment)

stat.test5 <- tibble::tribble(
  ~group1, ~group2, ~variable, ~p.adj, ~Name, ~Treatment,
  "Control", "Diclofenac", "Relative_Abundance", 'NS', "arsR BAC0589  ", "Control", 
  "Control", "Metformin", "Relative_Abundance", 'NS', "arsR BAC0589  ", "Control", 
  "Control", "17-beta", "Relative_Abundance", 'NS', "arsR BAC0589  ", "Control")


p2 <- ggplot(subset(d_plot_plasmid_true, Gen_num=="BAC0589"), aes(x=Treatment, y=Relative_Abundance, fill=Treatment))+
  geom_boxplot(outlier.shape =  NA, lwd = 0.8)+
  geom_point(position=position_jitterdodge(), alpha=0.7, size = 2, shape = 21, show.legend = FALSE)+
  theme_bw()+
  labs(x="NAD Treatment", y="Relative Abundance")+
  # scale_fill_manual(values = c("#eace17", "#6dbc18", "#34b6c6", "#ed6a48"), name="Treatment",
  #                   labels=c("Control", "Diclofenac", "Metformin", "17-β-estradiol"))+
  scale_x_discrete(labels = c("Control", "Diclofenac", "Metformin", "17-β-estradiol"))+
  scale_fill_manual(values = c("#ddc000", "#79ad41", "#34b6c6", "#d7aca1"), name="Treatment", 
                    labels=c("Control", "Diclofenac", "Metformin", "17-β-estradiol"))+ 
  theme(text=element_text(size=22))+
  theme(axis.title = element_text(size=20))+
  theme(strip.text = element_text(face = "italic"))+
  theme(strip.background=element_rect(colour="black",
                                      fill="white"))+
  stat_pvalue_manual(stat.test5, y.position = 0.015, label = "p.adj", xmin="group2", xmax=NULL, 
                     size = 5)+
  xlab(" ")+
  facet_wrap(~Name)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylab(" ")


#### Gene 3 - arsR BAC0594

df_arsR_BAC0594 <- subset(d_plot_plasmid_true, Gen_num =="BAC0594")

dunn.test(df_arsR_BAC0594$Relative_Abundance, df_arsR_BAC0594$Treatment)

stat.test2 <- tibble::tribble(
  ~group1, ~group2, ~variable, ~p.adj, ~Name, ~Treatment,
  "Control", "Diclofenac", "Relative_Abundance", '*', "arsR BAC0594  ", "Control", 
  "Control", "Metformin", "Relative_Abundance", 'NS', "arsR BAC0594  ", "Control", 
  "Control", "17-beta", "Relative_Abundance", '*', "arsR BAC0594  ", "Control")


p3 <-  ggplot(subset(d_plot_plasmid_true, Gen_num=="BAC0594"), aes(x=Treatment, y=Relative_Abundance, fill=Treatment))+
  geom_boxplot(outlier.shape =  NA, lwd = 0.8)+
  geom_point(position=position_jitterdodge(), alpha=0.7, size = 2, shape = 21, show.legend = FALSE)+
  theme_bw()+
  labs(x="NAD Treatment", y="Relative Abundance")+
  # scale_fill_manual(values = c("#eace17", "#6dbc18", "#34b6c6", "#ed6a48"), name="Treatment",
  #                   labels=c("Control", "Diclofenac", "Metformin", "17-β-estradiol"))+ 
  scale_fill_manual(values = c("#ddc000", "#79ad41", "#34b6c6", "#d7aca1"), name="Treatment", 
                    labels=c("Control", "Diclofenac", "Metformin", "17-β-estradiol"))+
  theme(text=element_text(size=22))+
  theme(axis.title = element_text(size=20))+
  scale_x_discrete(name = "Treatment", labels = c("Control", "Diclofenac", "Metformin", "17-β-estradiol"))+
  theme(strip.text = element_text(face = "italic"))+
  theme(strip.background=element_rect(colour="black",
                                      fill="white"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  facet_wrap(~Name)+
  ylab(" ")+
  stat_pvalue_manual(stat.test2, y.position = 0.015, label = "p.adj", xmin="group2", xmax=NULL, 
                     size = 5)



#### Gene 4 - merA BAC0652

df_merA_BAC0652 <- subset(d_plot_plasmid_true, Gen_num =="BAC0652")

dunn.test(df_merA_BAC0652$Relative_Abundance, df_merA_BAC0652$Treatment)

stat.test3 <- tibble::tribble(
  ~group1, ~group2, ~variable, ~p.adj, ~Name, ~Treatment,
  "Control", "Diclofenac", "Relative_Abundance", 'NS', "merA BAC0652  ", "Control", 
  "Control", "Metformin", "Relative_Abundance", 'NS', "merA BAC0652  ", "Control", 
  "Control", "17-beta", "Relative_Abundance", '*', "merA BAC0652  ", "Control")

p4 <- ggplot(subset(d_plot_plasmid_true, Gen_num=="BAC0652"), aes(x=Treatment, y=Relative_Abundance, fill=Treatment))+
  geom_boxplot(outlier.shape =  NA, lwd = 0.8)+
  geom_point(position=position_jitterdodge(), alpha=0.7, size = 2, shape = 21, show.legend = FALSE)+
  theme_bw()+
  labs(x="NAD Treatment", y="Relative Abundance")+
  scale_fill_manual(values = c("#ddc000", "#79ad41", "#34b6c6", "#d7aca1"), name="Treatment", 
                    labels=c("Control", "Diclofenac", "Metformin", "17-β-estradiol"))+
  theme(text=element_text(size=22))+
  theme(axis.title = element_text(size=20))+
  theme(strip.text = element_text(face = "italic"))+
  scale_x_discrete(name = "Treatment", labels = c("Control", "Diclofenac", "Metformin", "17-β-estradiol"))+
  theme(strip.background=element_rect(colour="black",
                                      fill="white"))+
  facet_wrap(~Name)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  stat_pvalue_manual(stat.test3, y.position = 0.029, label = "p.adj", xmin="group2", xmax=NULL, 
                     size = 5)




#### Gene 5 - merA BAC0652

df_merP_BAC0231 <- subset(d_plot_plasmid_true, Gen_num =="BAC0231")

dunn.test(df_merP_BAC0231$Relative_Abundance, df_merP_BAC0231$Treatment)


stat.test4 <- tibble::tribble(
  ~group1, ~group2, ~variable, ~p.adj, ~Name, ~Treatment,
  "Control", "Diclofenac", "Relative_Abundance", 'NS', "merP BAC0231  ", "Control", 
  "Control", "Metformin", "Relative_Abundance", '*', "merP BAC0231  ", "Control", 
  "Control", "17-beta", "Relative_Abundance", 'NS', "merP BAC0231  ", "Control")



p5 <- ggplot(subset(d_plot_plasmid_true, Gen_num=="BAC0231"), aes(x=Treatment, y=Relative_Abundance, fill=Treatment))+
  geom_boxplot(outlier.shape =  NA, lwd = 0.8)+
  geom_point(position=position_jitterdodge(), alpha=0.7, size = 2, shape = 21, show.legend = FALSE)+
  theme_bw()+
  labs(x="NAD Treatment", y="Relative Abundance")+
  scale_fill_manual(values = c("#ddc000", "#79ad41", "#34b6c6", "#d7aca1"), name="Treatment", 
                    labels=c("Control", "Diclofenac", "Metformin", "17-β-estradiol"))+
  # scale_fill_manual(values = c("#eace17", "#6dbc18", "#34b6c6", "#ed6a48"), name="Treatment",
  #                   labels=c("Control", "Diclofenac", "Metformin", "17-β-estradiol"))+  # theme(text=element_text(size=20))+
  theme(text=element_text(size=22))+
  scale_x_discrete(name = "Treatment", labels = c("Control", "Diclofenac", "Metformin", "17-β-estradiol"))+
  theme(axis.title = element_text(size=20))+
  theme(strip.text = element_text(face = "italic"))+
  ylab(" ")+
  facet_wrap(~Name)+
  theme(strip.background=element_rect(colour="black",
                                      fill="white"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  stat_pvalue_manual(stat.test4, y.position = 0.002, label = "p.adj", xmin="group2", xmax=NULL, 
                     size = 5)



#### Multiplot with sig values ####

#Multiplopt

p1 + p2 + p3+ p4 + p5 + guide_area() +
  plot_layout(guides = 'collect')



ggsave("Figure 5.tiff", height = 14, width = 16, dpi = 300, path = "figures")






#### NMDS on total BMRG Diversity ####
d_nmds <- gathered_merged
d_nmds <- subset(gathered_merged, select = c(1, 4:5))

d_nmds <- dcast(d_nmds, concrep ~ Gene, sum)

d_nmds <- spread(d_nmds, key=concrep, value=Relative_Abundance)

write.csv(d_nmds, "To be transposed BMRGs.csv")

d_nmds <- read.csv("Transposed BMRGs.csv", header=TRUE)

#Make Treatment and Rep Columns
d_nmds <- separate(d_nmds, Gene, c("Treatment", "Rep"), sep=" ")

#This d_nmds is what we want !


#Do the NMDS 
pc = d_nmds

#make community matrix - extract columns with abundance information
com = pc[ ,3:ncol(pc)]
m_com = as.matrix(com)


#Set seed gives you the same results each time you run this code
set.seed(123)
nmds = metaMDS(m_com, distance = "bray")
nmds


# Call:
#   metaMDS(comm = m_com, distance = "bray") 
# 
# global Multidimensional Scaling using monoMDS
# 
# Data:     wisconsin(m_com) 
# Distance: bray 
# 
# Dimensions: 2 
# Stress:     4.850083e-05 
# Stress type 1, weak ties
# Best solution was not repeated after 20 tries
# The best solution was from try 9 (random start)
# Scaling: centring, PC rotation, halfchange scaling 
# Species: expanded scores based on ‘wisconsin(m_com)

plot(nmds)


data.scores = as.data.frame(scores(nmds)$sites)

#add columns to data frame 
data.scores$Sample = pc$Treatment


head(data.scores)

data.scores$Sample <- factor(data.scores$Sample, levels=unique(data.scores$Sample))

data.scores$Sample <- factor(data.scores$Sample, levels=c("Control", "Diclofenac", "Metformin", "17-beta"))


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

ggsave("NMDS Genes Total BMRGs.tiff", height = 8, width = 10, path = "figures")

ano = anosim(m_com, pc$Treatment, distance = "bray", permutations = 9999)
ano

# Call:
#   anosim(x = m_com, grouping = pc$Treatment, permutations = 9999,      distance = "bray") 
# Dissimilarity: bray 
# 
# ANOSIM statistic R: 0.03167 
# Significance: 0.2799 
# 
# Permutation: free
# Number of permutations: 9999
