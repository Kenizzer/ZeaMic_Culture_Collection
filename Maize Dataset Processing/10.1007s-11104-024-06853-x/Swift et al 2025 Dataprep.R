# Swift, J.F., Kolp, M.R., Carmichael, A. et al. Drought stress homogenizes maize
# growth responses to diverse natural soil microbiomes. Plant Soil 509, 181–199
# (2025). https://doi.org/10.1007/s11104-024-06853-x
## Code by: Joel Swift 2026

####### Packages #######
library(tidyverse)
library(phyloseq)
library(lmerTest)
library(genefilter)

# Load data (non vst transformed)
hyb.bact.asv <- readRDS('intermediate_data/phyloseq_b_asv_clean_allRoots.RDS')
# metadata
hyb.bact.asv@sam_data
### Export fasta
# Subset by sample type
root.ps <- subset_samples(hyb.bact.asv, Type  == "root")
# Remove taxa with zero abundance
root.ps <- prune_taxa(taxa_sums(root.ps) > 0, root.ps)
# taxa dataframes
root.tax <- data.frame(root.ps@tax_table)
# Create a taxonomy string for each row
root.tax$taxonomy <- apply(root.tax[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")], 1, function(x) paste(na.omit(x), collapse = "_"))
# Combine ASV ID and taxonomy into a FASTA header
root.tax$header <- paste0(">", "Swift", "_", sub("^b", "", rownames(root.tax)), "-", "ROOT", "-", "USA_KS", " ", root.tax$taxonomy)
# Interleave header and write sequence lines
root_fasta_lines <- c(rbind(root.tax$header, root.tax$Sequence))
writeLines(root_fasta_lines, "Root_ASVs_with_taxonomy_Swift2024.fasta")