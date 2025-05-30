# Analyses and Data for the manuscript 'Antimicrobial effects, and selection for AMR by non-antibiotic drugs in bacterial communities' in Environment International
DOI of paper at Environment International https://doi.org/10.1016/j.envint.2025.109490 

Sequence data is deposited at ENA - Accession Number PRJEB74464. 

All other data is deposited here as part of this repository ![image](https://github.com/user-attachments/assets/b35d7bf5-a275-41c4-8bd1-efbfc0a69b6b)


##  Outline 
This repository contains the datasets, analyses, and figures for the paper 'Antimicrobial effects, and selection for AMR by non-antibiotic drugs on bacterial communities'. The code present can recreate all figures and tables. The raw sequencing data for the metagenomics analyses are available at ENA Accession Number PRJEB74464, study number ERP159148. The code for running the bioinformatics pipelines (ARGs-OAP, Bacmet, MetaPhlAn) are not present here, but all pipelines were run with default parameters. The code for running these are available on request. [ARGs-OAP](https://github.com/xinehc/ARGs_OAP) version 2.0 was used to identify ARGs present (both at antibiotic class and individual gene level) in all samples. In these analyses, ARG hits normalised to 16S rRNA counts (as generated by the ARGs-OAP pipeline) were used. [BacMet](http://bacmet.biomedicine.gu.se/) was used to determine metal and biocide resistance genes present in the samples. The BacMet outputs were normalised to the 16S reads generated in ARGs-OAP and hits were filtered to match the default parameters of ARGs-OAP (i.e. 80% identity match of at least 25bp). [MetaPhlAn 2.0](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn2) was used to identify the taxonomy of the samples, FLASH2 was used to pair the forward and reverse reads, and these paired reads were used in the MetaPhlAn pipeline. 

This repository contains the qPCR and plating datasets, and will contain all code for the analyses (once accepted) for the above-mentioned paper. 

## Information about the files 

**Growth Data**

There are several dataframes of growth data (i.e. Diclofenac_growth). These dataframes contain the matrix read out from the plate reader. I.e. they contain a first column of well (labelled with treatment and replicate) and subsequent columns of time points and OD read at that time point. For diclofenac, metformin, and 17-β-estradiol both high and low datasets are included, however the lowest concentration dataset is used to generate LOECs. The highest concentration dataset contains OD readings for growth in the higher concentations of these NADs, however a LOEC was not able to be identified with these concentrations, hence why the lower concentrations were also used separately.

**Plating Data** 

Full results of the plating experiment. Community was exposed to diclofenac, metformin, and 17-β-estradiol at 20uM for seven days, and at the end of the seven day experiment plated onto coliform selective agar containing one of six (or no antibiotic) antibiotics. The data is split by *E. coli* and presumptive coliforms since these can be identified as they grow either purple or pink on coliform selective agar. The antibioitcs used were ampicillin (amp), azithromyicn (azi), cefotaxime (cef), ciprofloxacin (cip), gentamicin (gen), and tetracycline (tet). Each treatment is recorded as NAD_biological replicate technical replicate (i.e. D_1 R1 means this is diclofenac biological replicate 1, technical replicate 1). The columns of the dataset indicate the dilution plated. 

**qPCR Data**

There are three sets of qPCR results (i.e. Diclofenac_qPCR etc.). These results hold the day 0 (d0) and day 7 (d7) qPCR quantities for the four genes targeted - 16S _rRNA_, _intI1_, _cintI1_, and _sul1_

**Note on the sequencing data**

The diclofenac and 17-β-estradiol experiments were sequenced using NEB Ultra ll FS library prep, and the metformin experiment used a PCR free library prep. There is a set of inoculum sequenced using both library preps to account for any biases when comparing to the inoculum. The inoculum to be used with the diclofenac and 17-β-estradiol analyses are 'nad_inoculum_r1_R1, nad_inoculum_r1_R2' etc. The inoculum files to be used metformin anlayses are 'nad_inoc_r1_R1, nad_inoc_r1_R2' etc. 


**Code**

R Code for these analyses are found here. There is code to analyse the ARGs-OAP, MetaPhlAn, and BacMet pipeline outputs in R. There is also code to analyse the qPCR, plating assays, and example code to analyse the growth curves.
