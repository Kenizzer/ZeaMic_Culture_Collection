# Data Processing for "Analysis of "Heterosis of leaf and rhizosphere microbiomes
# in field-grown maize"
# https://doi.org/10.1111/nph.16730
 
# Data retrieved from Zenodo: http://doi.org/10.5281/zenodo.3596465 


####### Packages #######
library(tidyverse)
library(phyloseq)
library(lmerTest)
library(genefilter)

# Load data (non vst transformed)
hyb.bact.asv <- readRDS('intermediate_data/phyloseq_b_asv.RDS')
# metadata
sample_data(hyb.bact.asv)

### Export fasta
# Subset by sample type
leaf.ps <- subset_samples(hyb.bact.asv, Type == "leaf")
rhizo.ps <- subset_samples(hyb.bact.asv, Type == "rhizosphere")
# Remove taxa with zero abundance
leaf.ps <- prune_taxa(taxa_sums(leaf.ps) > 0, leaf.ps)
rhizo.ps <- prune_taxa(taxa_sums(rhizo.ps) > 0, rhizo.ps)
# taxa dataframes
leaf.tax <- data.frame(leaf.ps@tax_table)
rhizo.tax <- data.frame(rhizo.ps@tax_table)
# class(leaf.tax)
# class(rhizo.tax)
# Create a taxonomy string for each row
leaf.tax$taxonomy <- apply(leaf.tax[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")], 1, function(x) paste(na.omit(x), collapse = "_"))
rhizo.tax$taxonomy <- apply(rhizo.tax[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")], 1, function(x) paste(na.omit(x), collapse = "_"))
# Combine ASV ID and taxonomy into a FASTA header
leaf.tax$header <- paste0(">", rownames(leaf.tax), "-", "PHYLLOSPHERE", "-", "USA_NC", " ", leaf.tax$taxonomy)
rhizo.tax$header <- paste0(">", rownames(rhizo.tax), "-", "RHIZOSPHERE", "-", "USA_NC", " ", rhizo.tax$taxonomy)
# Interleave header and write sequence lines
leaf_fasta_lines <- c(rbind(leaf.tax$header, leaf.tax$Sequence))
rhizo_fasta_lines <- c(rbind(rhizo.tax$header, rhizo.tax$Sequence))
writeLines(leaf_fasta_lines, "Phyllosphere_ASVs_with_taxonomy_Wagner2020.fasta")
writeLines(rhizo_fasta_lines, "Rhizosphere_ASVs_with_taxonomy_Wagner2020.fasta")