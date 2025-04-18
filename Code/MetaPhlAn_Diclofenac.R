########################### Diclofenac MetaPhlAn Script #################################


#### Load Libraries ####

library(tidyverse)
library(ggplot2)
library(reshape2) 
library(MetBrewer)
library(vegan)
library(broom)

#### Read in Code and Merge into one file ####

files <- list.files(path="data", full.names=TRUE)

#Make a function to read in the files
read_in_metaphlan <- function(filename){
  mydata <- read.table(filename)
  mydata <- subset(mydata, str_detect(V1, "s__"))
  mydata <- subset(mydata, !str_detect(V1, "t__"))
  mydata <- separate(mydata, V1, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "species"), sep = "\\|") 
  mydata$filename <- gsub(".extendedFrags_profile.txt"," ",basename(filename))
  
  return(mydata)
  
}

#Create a dataframe
d <- purrr::map_df(files, read_in_metaphlan) #combine all the tables and apply the function to all files 

#Replace white space
unique(d$filename)
d$filename <- trimws(d$filename)
unique(d$filename)

unique(d$species)
d$species <- trimws((d$species))
unique(d$species)

#Run the dataframe again to remove the s__ from the species names 
#s__
d$species <- gsub("s__", " ", d$species)

#Add Column with Treatment name to it so we can average
#Need to add one specific column properly


d <- d %>% mutate(d, Treatment = 
                    case_when(filename=="paired1r1" ~ "6400 1",
                              filename=="paired1r2" ~ "6400 2",
                              filename=="paired1r3" ~ "6400 3",
                              filename=="paired2r1" ~ "2100 1",
                              filename=="paired2r2" ~ "2100 2",
                              filename=="paired2r3" ~ "2100 3", 
                              filename=="paired3r1" ~ "710 1", 
                              filename=="paired3r2" ~ "710 2", 
                              filename=="paired3r3" ~ "710 3",
                              filename=="conc4_rep1.txt" ~ "240 1",
                              filename=="conc4_rep2.txt" ~ "240 2", 
                              filename=="conc4_rep3.txt" ~ "240 3", 
                              filename=="paired5r1" ~ "79 1", 
                              filename=="paired5r2" ~ "79 2", 
                              filename=="paired5r3" ~ "79 3", 
                              filename=="paired6r1" ~ "26 1",
                              filename=="paired6r2" ~ "26 2", 
                              filename=="paired6r3" ~ "26 3", 
                              filename=="paired7r1" ~ "8 1", 
                              filename=="paired7r2" ~ "8 2", 
                              filename=="paired7r3" ~ "8 3", 
                              filename=="paired0r1" ~ "0 1", 
                              filename=="paired0r2" ~ "0 2", 
                              filename=="paired0r3" ~ "0 3",
                              filename=="pairedinocr1" ~ "inoc 1",
                              filename=="pairedinocr2" ~ "inoc 2",
                              filename=="pairedinocr3" ~ "inoc 3"))


#Need to create a new column to split Treatment and Replicate within the current Treatment Column
d <- separate(d, Treatment, c('Concentration', 'Replicate'), sep=" ")

#Filter for only bacteria to be present
d <- subset(d, Kingdom=="k__Bacteria")


#Now average within treatments 
d_avg <- group_by(d, Concentration, species) %>% 
  summarise(Avg=mean(V2), .groups="drop")

 
d_avg_genus <- group_by(d, Concentration, Genus) %>%
  summarise(Avg=mean(V2), .groups="drop")

#Inoculum
d_noinoc <- subset(d_avg, Concentration!="inoc")
d_evolved <- subset(d, Concentration!="inoc")


d_avg_genus$Concentration <- factor(d_avg_genus$Concentration, levels=c("inoc", "0", "8", "26", "79", "240", "710", "2100", "6400"))

#Substitute _ for space in species names
d_noinoc$species <- gsub("_", " ", d_noinoc$species)




################ Alpha Diversity ###############################

richness_df <- d %>% 
  group_by(Concentration, Replicate, Genus) %>%
  count(Genus) 


richness_df$Concentration <- factor(richness_df$Concentration, levels=c("inoc", "0", "8", "26", "79", "240", "710", "2100", "6400"))

#richness per genus
richness_df_summ <- d %>%
  group_by(Concentration) %>%
  count(Genus)

richness_df_summ$Concentration <- factor(richness_df_summ$Concentration, levels=c("inoc", "0", "8", "26", "79", "240", "710", "2100", "6400"))

#This df is all genera per concentration w avg and sd
richness_summ_2 <- richness_df %>%
  group_by(Concentration, Replicate) %>%
  summarise(n=n())%>%
  summarise(avg = mean(n), sd = sd(n))


## Plot evolved samples alpha diversity
ggplot(subset(richness_summ_2, Concentration !="inoc"))+
  geom_point(aes(x=Concentration, y=avg, color=Concentration), size = 5) +
  geom_errorbar(aes(x = Concentration, ymin=avg-sd, ymax=avg+sd, colour = Concentration), width=.4, linewidth = 1) +
  scale_color_manual(values = met.brewer("VanGogh3", 8))+
  # scale_colour_manual(values = c("#F95639", "#B476FA", "#79ad41", "#34b6c6", "black"), name="Treatment")+
  theme_classic()  +
  theme(text=element_text(size=14))+
  labs(y= "Observed Genera", x="Diclofenac Concentration (µg/L)")+
  guides(colour = "none")

ggsave("figures/Supplementary Figure 4_Diclofenac Comm Richness.tiff", height = 6, width = 8, dpi = 300)


## Stats

m_rich <- lm(log(n) ~ Concentration, data = richness_df)
summary(m_rich)

plot(m_rich$resid)

qqnorm(m_rich$resid)
qqline(m_rich$resid)

anova(m_rich)


############################### Beta Diversity ##############################################

d_NMDS <- read.csv("NMDS Try Transposed.csv", header=TRUE)

d_NMDS <- separate(d_NMDS, species, c("Concentration", "Replicate"), sep="_")

d_NMDS <- subset(d_NMDS, Concentration != "inoc" )


# do nmds

#do the NMDS 
pc = d_NMDS

#make community matrix - extract columns with abundance information
com = pc[ ,3:ncol(pc)]
m_com = as.matrix(com)



#Set seed gives you the same results each time you run this code
set.seed(123)
nmds = metaMDS(m_com, distance = "bray")
nmds

data.scores = as.data.frame(scores(nmds)$sites)

#add columns to data frame 
data.scores$Sample = pc$Concentration

head(data.scores)

data.scores$Sample <- factor(data.scores$Sample, levels=unique(data.scores$Sample))

data.scores$Sample <- factor(data.scores$Sample, levels=c("0", "8", "26", "79", "240", "710", "2100", "6400"))






ano = anosim(m_com, pc$Concentration, distance = "bray", permutations = 9999)
ano



#Plot in ggplot if required

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
  labs(x = "NMDS1", colour = "Diclofenac \nConcentration (µg/L)", y = "NMDS2", shape = "Type")+   
  scale_colour_manual(values = met.brewer("VanGogh3", 8)) 

# ggsave("NMDS evolved.tiff", height = 8, width = 8, path = "figures")



############ Supplementary Figure  7 - Genera by Treatment #####################

ggplot(subset(d_avg_genus, Concentration!="inoc"), aes(x=Concentration, y=Genus, fill=Avg))+
  geom_tile()+
  scale_y_discrete(limits=rev)+
  scale_fill_gradientn(colours=met.brewer("VanGogh3"))+
  theme_classic()+
  theme(text=element_text(size=20))+
  labs(x="Diclofenac concentration (µg/L)", y="Genus", fill="Percentage of \ntotal community")+
  theme(axis.text.y = element_text(face = "italic"))


ggsave("Supplementary Figure 7 - MetaPhlAn Genera diclofenac, evolved samples.tiff", width=16, height=14, path="figures", dpi=300)




######################################### Stats #########################################

# no inoculum d_evolved 

## Need to create a dataframe with all samples containing all the species, even if0 - not just presence absence as this will not 
# be allowed in testing the way below

d2 <- unite(d, col="Treatment", c("Concentration", "Replicate"), sep=" ")

d2 <- group_by(d2, Treatment, species) %>%
  select(Treatment, species, V2)

d2 <- dcast(d2, species~Treatment, sum)

#gives dataframe with 0s instead of absences
d2_evolved <- melt(d2, id.vars = c("species"))

d2_evolved <- separate(d2_evolved, variable, c("Concentration", "Replicate"), sep=" ")

d2_evolved <- subset(d2_evolved, Concentration !="inoc")


#Kruskal-Wallis tests

#### Stats on only in all evolved treatments ####

models_kruskal <- select(d2_evolved, species) %>%
  distinct() %>%
  mutate(data = list(NA),
         model = list(NA),
         summary = list(NA))

for(i in 1:nrow(models_kruskal)){
  
  # grab species
  temp_species <- models_kruskal$species[i]
  
  # filter just for that amr class then the columns we need
  temp_data <- filter(d2_evolved, species==temp_species) %>%
    select(Concentration, value)
  
  # run kruskal wallis test
  temp_model <- kruskal.test(value ~ Concentration, temp_data)
  
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

#nothing is significantly different between treatments






