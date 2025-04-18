#### Code to analyse the MetaPhlAn data 17-beta-Estradiol ####

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
d$Genus <- gsub("g__", " ", d$Genus)

#Add Column with Treatment name to it so we can average?
#Need to add one specific column properly


d <- d %>% mutate(d, Treatment = 
                    case_when(filename=="pairedAr1" ~ "5400 1",
                              filename=="pairedAr2" ~ "5400 2",
                              filename=="pairedAr3" ~ "5400 3",
                              filename=="pairedBr1" ~ "1800 1",
                              filename=="pairedBr2" ~ "1800 2",
                              filename=="pairedBr3" ~ "1800 3", 
                              filename=="pairedCr1" ~ "607 1", 
                              filename=="pairedCr2" ~ "607 2", 
                              filename=="pairedCr3" ~ "607 3",
                              filename=="pairedDr1" ~ "201 1",
                              filename=="pairedDr2" ~ "201 2", 
                              filename=="pairedDr3" ~ "201 3", 
                              filename=="pairedEr1" ~ "67 1", 
                              filename=="pairedEr2" ~ "67 2", 
                              filename=="pairedEr3" ~ "67 3", 
                              filename=="pairedFr1" ~ "22 1",
                              filename=="pairedFr2" ~ "22 2", 
                              filename=="pairedFr3" ~ "22 3", 
                              filename=="pairedGr1" ~ "7 1", 
                              filename=="pairedGr2" ~ "7 2", 
                              filename=="pairedGr3" ~ "7 3", 
                              filename=="pairedHr1" ~ "2 1", 
                              filename=="pairedHr2" ~ "2 2", 
                              filename=="pairedHr3" ~ "2 3",
                              filename=="pairedIr1" ~ "0 1", 
                              filename=="pairedIr2" ~ "0 2", 
                              filename=="pairedIr3" ~ "0 3",
                              filename=="pairedinocR1" ~ "inoc 1",
                              filename=="pairedinocR2" ~ "inoc 2",
                              filename=="pairedinocR3" ~ "inoc 3"))


#Need to create a new column to split Treatment and Replicate within the current Treatment Column
d <- separate(d, Treatment, c('Concentration', 'Replicate'), sep=" ")

#Filter for only bacteria to be present
d <- subset(d, Kingdom=="k__Bacteria")

#Now average within treatments 
#I think this works, because we are left with Species name, and Concentration, and the V2 which is the percentage
d_avg <- group_by(d, Concentration, species) %>% 
  summarise(Avg=mean(V2), .groups="drop")

#Inoculum
d_noinoc <- subset(d_avg, Concentration!="inoc")

#Substitute _ for space in species names
d_noinoc$species <- gsub("_", " ", d_noinoc$species)

#Make the factors of the concentrations process properly
d_avg$Concentration <- factor(d_avg$Concentration, levels=c("inoc", "0", "2", "7", "22", "67", "201", "607", "1800", "5400"))
d_noinoc$Concentration <- factor(d_noinoc$Concentration, levels=c("0", "2", "7", "22", "67", "201", "607", "1800", "5400"))


d_avg_genus <- group_by(d, Concentration, Genus) %>% 
  summarise(Avg=mean(V2), .groups="drop")


d_avg_genus$Concentration <- factor(d_avg_genus$Concentration, levels=c("inoc", "0", "2", "7", "22", "67", "201", "607", "1800", "5400"))



#################################### Alpha Diversity ################################

richness_df <- d %>% 
  group_by(Concentration, Replicate, Genus) %>%
  count(Genus) 


richness_df$Concentration <- factor(richness_df$Concentration,
                                    levels=c("inoc","0", "2", "7", "22", "67", "201", "607", "1800", "5400"))

#richness per genus, probably not important
richness_df_summ <- d %>%
  group_by(Concentration) %>%
  count(Genus)

richness_df_summ$Concentration <- factor(richness_df_summ$Concentration,
                                         levels=c("inoc","0", "2", "7", "22", "67", "201", "607", "1800", "5400"))


#This one is all genera per concentration w avg and sd
richness_summ_2 <- richness_df %>%
  group_by(Concentration, Replicate) %>%
  summarise(n=n())%>%
  summarise(avg = mean(n), sd = sd(n))


## Plot this
ggplot(richness_summ_2, aes(x=Concentration, y=avg))+
  geom_point(size = 3)+
  theme_classic()+   
  labs(x="17-β-estradiol concentration concentration (µg/L)", y="Number of Genera")

ggplot(subset(richness_summ_2, Concentration !="inoc"))+
  geom_point(aes(x=Concentration, y=avg, color=Concentration), size = 5) +
  geom_errorbar(aes(x = Concentration, ymin=avg-sd, ymax=avg+sd, colour = Concentration), width=.4, linewidth = 1) +
  scale_color_manual(values = met.brewer("OKeeffe2", 9))+
  theme_classic()  +
  theme(text=element_text(size=14))+
  labs(y= "Observed Genera", x="17-β-estradiol concentration Concentration (µg/L)")+
  guides(colour = "none")

ggsave("figures/Supplementary Figure6_Estradiol Comm Richness.tiff", height = 6, width = 8, dpi = 300)


## Stats

m_rich <- lm(log(n) ~ Concentration, data = richness_df)
summary(m_rich)

plot(m_rich$resid)

qqnorm(m_rich$resid)
qqline(m_rich$resid)

anova(m_rich, test = "F")


############################### Beta Diversity ##############################################

d_NMDS <- read.csv("NMDS Transposed.csv", header=TRUE)

d_NMDS <- separate(d_NMDS, species, c("Concentration", "Replicate"), sep="_")

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

data.scores$Sample <- factor(data.scores$Sample, #
                             levels=c("0", "2", "7", "22", "67", "201", "607", "1800", "5400"))

#ANOSIM
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
  labs(x = "NMDS1", colour = "17-β-estradiol \nConcentration (µg/L)", y = "NMDS2", shape = "Type")+   
  scale_colour_manual(values = met.brewer("OKeeffe2", 9)) 


############ Supplementary Figure  9 - Genera by Treatment #####################

ggplot(subset(d_avg_genus, Concentration!="inoc"), aes(x=Concentration, y=Genus, fill=Avg))+
  geom_tile()+
  scale_y_discrete(limits=rev)+
  scale_fill_gradientn(colours=met.brewer("OKeeffe2"))+
  theme_classic()+
  theme(text=element_text(size=20))+
  labs(x="17-β-estradiol concentration (µg/L)", y="Genus", fill="Percentage of \ntotal community")+
  theme(axis.text.y = element_text(face = "italic"))


# ggsave("Supplementary Figure 9 - MetaPhlAn Genera diclofenac, evolved samples.tiff", width=16, height=14, path="figures", dpi=300)



