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

The script is written so that the code of any section can be run once
the “Preparation of the data” section has been run.

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
library(plyr) # To use the ddply function
library(ggplot2) # For plotting

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

There is a certain amount of things that needs to be modified on these
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

# Some cleanup to finish
rm(physeq1, OTU_list, i, class, family, genus, order, phylum, species)
```

## Barplots of the microbial composition

This is the code leading to figure 2. First, we want to use relative
abundance, so let’s convert the OTU table.

``` r
# Change counts by relative abundance per sample
otu_rel <- prop.table(as.matrix(otu), margin = 2)*100

# Make a new phyloseq object
OTU_REL <- otu_table(otu_rel, taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(tax))
samples = sample_data(metadata)
Phyloseq_rel <- phyloseq(OTU_REL, TAX, samples)
```

A we want to have a barplot at the phylum level and another one at the
genus level, we need to do some amalgams.

``` r
Phyloseq_genus <- tax_glom(Phyloseq_rel, taxrank = "Genus")
Phyloseq_phylum <- tax_glom(Phyloseq_rel, taxrank = "Phylum")
```

### Phylum barplot

Now we are all set. Let’s start with the barplot at the phylum level.
When building a barplot, one needs to choose which taxa should be
represented, and which ones should be pulled together in an “Other”
category as they are less abundant in the communities. While this
separation is arbitrary, it should balance the importance of information
carried by the taxa and the readability of the plot. For our phylum
barplot, we decided to pull together all phyla that are not present in
more than 1% in at least 1 sample. That makes that 7 phyla remained,
plus the “Other” group.

Let’s start by pulling together the low abundance phyla.

``` r
# For that we need to get the phyloseq object into the long format, and make sure that phyla are strings.
phylum_distribution <- psmelt(Phyloseq_phylum)
phylum_distribution$Phylum <- as.character(phylum_distribution$Phylum)

phylum_max <- ddply(phylum_distribution, ~Phylum, function(x) c(max=max(x$Abundance))) # What is the highest abundance for each phylum?
phylum_others <- phylum_max[phylum_max$max < 1,]$Phylum # Get the list of phyla that do not reach 1% in any sample.
phylum_others_OTU <- row.names(tax_table(Phyloseq_phylum))[apply(tax_table(Phyloseq_phylum), 1, function(u) any(u %in% phylum_others))] # What is the OTU name of these phyla?
Phyloseq_phylum <- merge_taxa(Phyloseq_phylum, phylum_others_OTU) # Now we can pool these OTUs together inside the phyloseq object
```

Now it still needs a bit of polishing to be able to plot the data

``` r
phylum_distribution <- psmelt(Phyloseq_phylum) # ggplot also needs data in the long format
phylum_distribution$Phylum <- as.character(phylum_distribution$Phylum) # Makes sure that phyla are strings
phylum_distribution$Phylum[is.na(phylum_distribution$Phylum)] <- "Others" # Gives the name "Others" to the pulled phyla

# We want to have the phyla on the plot orgasnised from most abundant to least abundant, with "Others" as last.
phylum_total <- (ddply(phylum_distribution, ~Phylum, function(x) c(sum=sum(x$Abundance)))) # Cumulative abundance of each OTU.
phylum_total <- phylum_total[order(phylum_total$sum),] # Order the phyla per cumulative abundance
order_barplot <- phylum_total$Phylum #Get the list of ordered phyla
order_barplot <- c(order_barplot[3], order_barplot[1:2], order_barplot[4:8]) # Move "Others" to the beginning

# ggplot has a nasty habit to shuffle samples and stuff, so let's factorise the samples, and the phylum order.
phylum_distribution$Sample <- factor(phylum_distribution$Sample, levels = colnames(otu))
phylum_distribution$Phylum <- factor(phylum_distribution$Phylum, levels = order_barplot)

 # Make a color palette
cols <- c("grey","hotpink","yellow","purple","powderblue","orange","green","steelblue")
```

And finally, we can plot figure 2a! Note that some minor modification
have been made to the code to fit the markdown output format (eg. the
size variable)

``` r
ggplot(phylum_distribution, aes(x = Sample,
                                y = Abundance, 
                                fill = Phylum)) +
  geom_bar(stat="identity", colour = "black") +
  labs(x = "", y = "Relative abundance")+ 
  scale_fill_manual(values = cols)+
  facet_grid(~factor(sample_Species, levels=c("Alaria-esculenta", "Saccharina-latissima","No-seaweed")), scales = 'free', space = 'free')+
  theme_bw(base_size = 12)+
  theme(axis.text.x = element_text(angle = 45, size = 8))+
  theme(legend.position = "right", legend.text=element_text(size=10))+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.border = element_blank()) + guides(fill=guide_legend(ncol = 1))
```

![](Rcode-microbiology_files/figure-gfm/plot-1.png)<!-- --> \### Genus
barplot The principle here is the same than for the phylum barplot, but
this time we decided to pull together all genera that are not present in
more than 1.6% in at least 1 sample. That makes that 19 genera kept,
plus the “Other” group.

Again, let’s start by pooling together the low abundance genera.

``` r
genus_distribution <- psmelt(Phyloseq_genus)
genus_distribution$Genus <- as.character(genus_distribution$Genus)
genus_max <- ddply(genus_distribution, ~Genus, function(x) c(max=max(x$Abundance)))
genus_others <- genus_max[genus_max$max < 1.6,]$Genus
genus_others_OTU <- row.names(tax_table(Phyloseq_genus))[apply(tax_table(Phyloseq_genus), 1, function(u) any(u %in% genus_others))]
Phyloseq_genus <- merge_taxa(Phyloseq_genus, genus_others_OTU)
```

The polishing part. Again, same as for the phylum plot.

``` r
genus_distribution <- psmelt(Phyloseq_genus)
genus_distribution$Genus <- as.character(genus_distribution$Genus)
genus_distribution$Genus[is.na(genus_distribution$Genus)] <- "Others"
genus_total <- (ddply(genus_distribution, ~Genus, function(x) c(sum=sum(x$Abundance))))
genus_total <- genus_total[order(genus_total$sum),]
order_barplot <- genus_total$Genus
order_barplot <- c(order_barplot[16], order_barplot[1:15], order_barplot[17:20])
genus_distribution$Sample <- factor(genus_distribution$Sample, levels = colnames(otu))
genus_distribution$Genus <- factor(genus_distribution$Genus, levels = order_barplot)
colg<-c("lightgrey","salmon","yellow","powderblue","green","darkred","tomato","indianred",
        "purple","orange","greenyellow","orangered","midnightblue","lightskyblue","slateblue",
        "yellowgreen", "seagreen","lightsalmon","lightseagreen","steelblue")
```

And the plot! This is figure 2b.

``` r
ggplot(genus_distribution, aes(x = Sample,
                               y = Abundance, 
                               fill = Genus)) +
  geom_bar(stat="identity", colour = "black") +
  labs(x = "", y = "Relative abundance")+ 
  scale_fill_manual(values = colg)+
  facet_grid(~factor(sample_Species, levels=c("Alaria-esculenta", "Saccharina-latissima","No-seaweed")), scales = 'free', space = 'free')+
  theme_bw(base_size = 12)+
  theme(axis.text.x = element_text(angle = 45, size = 8))+
  theme(legend.position = "right", legend.text=element_text(size=10), legend.key.size = unit(0.5, 'cm'))+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.border = element_blank()) + guides(fill=guide_legend(ncol = 1))#Figure2B
```

![](Rcode-microbiology_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->
