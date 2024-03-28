# Analyses and Data for the manuscript 'Antimicrobial effects, and selection for AMR by non-antibiotic drugs on bacterial communities' submitted to Microbiome
DOI of paper XXX

Sequence data is deposited at ENA XXX

Plating and qPCR data is deposited here as part of this repository.

##  Outline 
This repository contains the datasets, analyses, and figures for the paper 'Antimicrobial effects, and selection for AMR by non-antibiotic drugs on bacterial communities'. The code present can recreate all figures and tables. The raw sequencing data for the metagenomics analyses are available at XXX. The code for running the bioinformatics pipelines (ARGs-OAP, Bacmet, MetaPhlAn) are not present here, but are available upon request. 

This repository contains the qPCR and plating datasets, all analyses (including for sequence data) and figures of the above-mentioned paper. It can recreate all figures in both the main text and the supplement.

## Information about the files 
**Plating Results** - full results of the plating experiment. Community was exposed to diclofenac, metformin, and 17-beta-estradiol at 20uM for seven days, and at the end of the seven day experiment plated onto coliform selective agar containing one of six (or no antibioitc) antibioitcs. The data is split by E. coli and presumptive coliforms since these can be identified as they grow either purple or pink on coliform selective agar. The antibioitcs used were ampicillin (amp), azithromyicn (azi), cefotaxime (cef), ciprofloxacin (cip), gentamicin (gen), and tetracycline (tet). Each treatment is recorded as NAD_biological replicate technical replicate (i.e. D_1 R1 means this is diclofenac biological replicate 1, technical replicate 1). The columns of the dataset indicate the dilution plated. 

**Diclofenac_ARG Transposed Dataframe_total.csv, Metformin_ARG Transposed Dataframe_total.csv, 17-beta-estradiol_ARG Transposed Dataframe_total.csv** - These data contain the transposed dataframes that need to be read into R to be able to do NMDS plots. These dataframes were created by reading the dataframe created from the ARG analysis scripts (e.g. XXX), and transposing them in Excel. 
