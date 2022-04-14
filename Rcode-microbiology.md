Rcode-microbiology
================
Sven Le Moine Bauer
2022-04-14

## Introduction

This is the R code for the microbial analysis of described in **Feed
characteristics and potential effects on ruminal bacteria of ensiled
Saccharina latis-sima and Alaria esculenta for dairy cows** from Yen *et
al.*, 2022. While we tried to make the code as clear as possible, we
have no pretension of being highly skilled R coders, and it is likely
that some sections could be rewritten in a more idiomatic manner. Any
comments are welcome. The code covers the following parts:

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
directory.

``` r
library(qiime2R) # to import QIIME output files into R
library(phyloseq) # To be able to play with phyloseq objects


# Set the directory to the one the script is saved in
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
```

And now let’s import the data. The code assumes that the files are in
the same directory than the script.

``` r
# Import the metadata to R
metadata <- read.csv("ch2-sample-meta.tsv", row.names = 1, sep = "\t")
# Make a phyloseq object out of the QIIME files
physeq1 <- qza_to_phyloseq( 
  features = "filtered_table.qza",
  taxonomy = "taxonomy.qza")
```

There is a certain amount of things that need to be modified on these
files before doing the analysis.

``` r
# It is generaly not advised to keep sample names as numbers, so let's add "S" in front
otu <- as.data.frame(otu_table(physeq1)) # export the OTU table from phyloseq as data frame
otu <- otu[,order(as.numeric(colnames(otu)))] # Order the samples
colnames(otu) <- paste("S", colnames(otu), sep = "") # Paste "S" in front of each sample name
rownames(metadata) <- paste("S", rownames(metadata), sep = "") # Make the change in the metadata table too

# Let's remove Archaea from the dataset as they are not supposed to be targeted by the primer pair.
tax <- as.data.frame(tax_table(physeq1)) # export the taxonomy table from phyloseq as data frame
tax <- subset(tax, tax$Kingdom %in% "Bacteria") # Keep only Bacteria
otu<- subset(otu, row.names(otu) %in% row.names(tax)) # Subset the OTU table too

# While not mandatorily needed, it is also nicer to have proper OTU names rather than a nonsensical string.
identical(rownames(otu), rownames(tax)) # Checks if the order is same in the otu and tax tables. Has to be true!
```

    ## [1] TRUE

``` r
OTU_list <- paste("OTU_", c(1:1865), sep="") # Makes a list of new names
rownames(otu) <- OTU_list # Change the names to the otu table
rownames(tax) <- OTU_list # Change the names to the tax table

# We need to replace the NAs in the tax table by "Unassigned_XXX", otherwise we may have some nonsensical taxonomic amalgams later on.
tax[is.na(tax)] <- "" # Replaces NA by nothing.
for (i in 1:nrow(tax)){tax[i,] <- as.character(tax[i,])} # Make sure that everything is a string in the table
for (i in 1:nrow(tax)){
  if (tax[i,2] == ""){
    phylum <- paste("Unclassified_", tax[i,1], sep = "")
    tax[i, 2:7] <- phylum
  } else if (tax[i,3] == ""){
    class <- paste("Unclassified_", tax[i,2], sep = "")
    tax[i, 3:7] <- class
  } else if (tax[i,4] == ""){
    order <- paste("Unclassified_", tax[i,3], sep = "")
    tax[i, 4:7] <- order
  } else if (tax[i,5] == ""){
    family <- paste("Unclassified_", tax[i,4], sep = "")
    tax[i, 5:7] <- family
  } else if (tax[i,6] == ""){
    genus <- paste("Unclassified_", tax[i,5], sep = "")
    tax[i, 6:7] <- genus
  } else if (tax[i,7] == ""){
    species <- paste("Unclassified_", tax[i,6], sep = "")
    tax[i, 7] <- species
  }
}

# Some cleanup to finsish
rm(physeq1, OTU_list, i, class, family, genus, order, phylum, species)
```
