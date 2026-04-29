#R code used to generate OTU tables for Bioreactor Data. 
#Will export to python, and do further analyses.
library(tidyr)
library(dplyr)  # optional, but often useful with tidyr
library(readr)
library(phyloseq)
library(qiime2R)
library(microViz)
library(vegan)
library(ggplot2)
setwd("~/Juelab/BR2019_Data/BR2019_qiime2")
######make phyloseq object 2019
data.phyBR19 <- qza_to_phyloseq(
  features = "BR_2019_table_F260_R200.qza",
  tree = "BR_2019_rooted_tree_F260_R200.qza",
  taxonomy = "BR_2019_custom_taxonomy.qza",
  metadata = "BR_2019_metadata.tsv")
#metadata = "stats_data/metadata.tsv") for just the qiime2 metadata output
tax.table<-data.phyBR19@tax_table
head(tax.table)
# get rid of taxa not present in any samples, and of samples with no reads
data.phyBR19 <- prune_taxa(taxa_sums(data.phyBR19) > 0, data.phyBR19)
data.phyBR19 <- prune_samples(sample_sums(data.phyBR19) > 0, data.phyBR19)

sample_names(data.phyBR19)

# rarefy samples (get rid of any with too few reads)

OTUtab <- otu_table(data.phyBR19)
class(OTUtab) <- "matrix"
spiketab <- t(OTUtab) # transpose observations to rows
head(OTUtab)
head(spiketab)
write.csv(spiketab, "BR19_OTUs.csv")

rarecurve(spiketab, step=50, lwd=2, ylab="OTU")
set.seed(1234)
data_rare.phyBR19 <- rarefy_even_depth(data.phyBR19, min(sample.size=1000), 
                                       rngseed=FALSE, replace=FALSE, 
                                       trimOTUs=TRUE, verbose=TRUE)


#export rarefied OTUs
# Extract OTU table from rarefied data
OTUtab_rare <- otu_table(data_rare.phyBR19)
class(OTUtab_rare) <- "matrix"

# Transpose (same as before)
spiketab_rare <- t(OTUtab_rare)

# Check
head(spiketab_rare)

# Export
write.csv(spiketab_rare, "BR19_OTUs_rarefied.csv")


######2024
#########################################
#BR2024 data analysis using updated data
########################################
#set working direc to
setwd("~/Juelab/BR2024")
#####make phyloseq object 2024
data.phyBR24 <- qza_to_phyloseq(
  features = "BR_2024_table_F240_R220.qza",
  tree = "BR2024_rooted_tree_F240_R220.qza",
  taxonomy = "BR2024_taxonomy_F240_R220.qza",
  metadata = "BR2024_metadata_F240_R220.tsv")
#metadata = "stats_data/metadata.tsv") for just the qiime2 metadata output

# get rid of taxa not present in any samples, and of samples with no reads
data.phyBR24 <- prune_taxa(taxa_sums(data.phyBR24) > 0, data.phyBR24)
data.phyBR24 <- prune_samples(sample_sums(data.phyBR24) > 0, data.phyBR24)

sample_names(data.phyBR24)

# rarefy samples (get rid of any with too few reads)

OTUtab <- otu_table(data.phyBR24)
class(OTUtab) <- "matrix"
spiketab <- t(OTUtab) # transpose observations to rows

write.csv(spiketab, "BR24_OTUs.csv")

rarecurve(spiketab, step=50, lwd=2, ylab="OTU")
set.seed(1234)
data_rare.phyBR24 <- rarefy_even_depth(data.phyBR24, min(sample.size=1000), 
                                       rngseed=FALSE, replace=FALSE, 
                                       trimOTUs=TRUE, verbose=TRUE)

# Extract OTU table from rarefied data
OTUtab_rare <- otu_table(data_rare.phyBR24)
class(OTUtab_rare) <- "matrix"

# Transpose 
spiketab_rare <- t(OTUtab_rare)

head(spiketab_rare)

write.csv(spiketab_rare, "BR24_OTUs_rarefied.csv")