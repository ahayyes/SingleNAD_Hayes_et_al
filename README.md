# Analyses and Data for the manuscript 'Antimicrobial effects, and selection for AMR by non-antibiotic drugs on bacterial communities' submitted to Microbiome
DOI of paper XXX

Sequence data is deposited at ENA XXX

Plating and qPCR data is deposited here as part of this repository.

##  Outline 
This repository contains the datasets, analyses, and figures for the paper 'Antimicrobial effects, and selection for AMR by non-antibiotic drugs on bacterial communities'. The code present can recreate all figures and tables. The raw sequencing data for the metagenomics analyses are available at XXX. The code for running the bioinformatics pipelines (ARGs-OAP, Bacmet, MetaPhlAn) are not present here, but are available upon request. 

This repository contains the qPCR and plating datasets, all analyses (including for sequence data) and figures of the above-mentioned paper. It can recreate all figures in both the main text and the supplement.

## Information about the files 

**Growth Data**
There are three dataframes of growth data (i.e. Diclofenac_growth). These dataframes contain the matrix read out from the plate reader. I.e. they contain a first column of well (labelled with treatment and replicate) and subsequent columns of time points and OD read at that time point.

**Plating Data** 
Full results of the plating experiment. Community was exposed to diclofenac, metformin, and 17-beta-estradiol at 20uM for seven days, and at the end of the seven day experiment plated onto coliform selective agar containing one of six (or no antibioitc) antibioitcs. The data is split by E. coli and presumptive coliforms since these can be identified as they grow either purple or pink on coliform selective agar. The antibioitcs used were ampicillin (amp), azithromyicn (azi), cefotaxime (cef), ciprofloxacin (cip), gentamicin (gen), and tetracycline (tet). Each treatment is recorded as NAD_biological replicate technical replicate (i.e. D_1 R1 means this is diclofenac biological replicate 1, technical replicate 1). The columns of the dataset indicate the dilution plated. 

**qPCR Data**
There are three sets of qPCR results (i.e. Diclofenac_qPCR etc.). These results hold the day 0 (d0) and day 7 (d7) qPCR quantities for the four genes targeted - 16S _rRNA_, _intI1_, _cintI1_, and _sul1_

**ARG Data**
There are three sets of datasets for these again, (i.e. Diclofenac_ARG Transposed Dataframe). These data contain the transposed dataframes that need to be read into R to be able to do NMDS plots. These dataframes were created by reading the dataframe created from the ARG analysis scripts (e.g. XXX), and transposing them in Excel. 

**BMRG Data**
Again, there are three sets of datasets of transposed BMRG datasets to be read into R to make the NMDS. These dataframes were created by reading the dataframe created in R from the BMRG analysis scripts, and then transposing them into Excel.
