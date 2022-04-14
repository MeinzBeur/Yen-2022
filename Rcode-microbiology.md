Rcode-microbiology
================
Sven Le Moine Bauer
2022-04-14

## Introduction

This is the R code for the microbial analysis of described in **Feed
characteristics and potential effects on ruminal bacteria of ensiled
Saccharina latis-sima and Alaria esculenta for dairy cows** from Yen *et
al.*, 2022. The code covers the following parts:

-   Preparation of the data
-   Barplots of the bacterial composition
-   Principal component analysis and associated statistical tests.
-   Rarefaction curves and Shannon diversity
-   Differential abundance analysis using ANCOM-BC and ALDEx2
-   Supplementary figures

## Preparation of the data

The data used here can be downloaded from the Github repository, and
contains 3 files:

-   “filtered_table.qza” is the OTU table produced through QIIME.
-   “taxonomy.qza” contains the taxonomic assignments to each OTU.
    Produced through QIIME.
-   The “ch2-sample-meta.tsv” file contains the metadata for each
    sample.

First we need load the libraries used here, and set up the working
directory

``` r
library(qiime2R) # to import QIIME output files into R



# Set the directory to the one the script is saved in
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
```

And now let’s import the data:

``` r
SVs <- read_qza(file="filtered_table.qza") 
Tx <- read_qza(file="taxonomy.qza") 
metadata <- read.csv("ch2-sample-meta.tsv", row.names = 1, sep = "\t")
```
