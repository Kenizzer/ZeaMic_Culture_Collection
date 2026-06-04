# Data processing for "Specific and conserved patterns of microbiota-structuring 
# by maize benzoxazinoids in the field"
# https://doi.org/10.1186/s40168-021-01049-2

# Note, Cadot et al. 2021 used the samples from the study below in their analysis
# we extracted these samples from the combine tables

# Data processing for "Root exudate metabolites drive plant-soil feedbacks on 
# growth and defense by shaping the rhizosphere microbiota"
#https://doi.org/10.1038/s41467-018-05122-7

####### Packages #######
library(tidyverse)
library(phyloseq)
library(lmerTest)
library(genefilter)

####### Set parameters #######
ASVprefix <- 'Cadot_' # for bacterial ASVs

####### Metadata file #######
smd <- read.csv("./metadata/metadata_Cadot2021.txt", sep = "\t")
smd$sample_name <- smd$SampleID
smd <- as.data.frame(smd)


# make the sample names match the ASV table
colnames(smd)
# Run + bac_f_bc + bac_r_bc
# Example R43_1193R_F70_799F
smd$SampleID <- smd$sample_ID_bac
rownames(smd) <- smd$SampleID

####### Load ASV table and taxa info #######
## ASV table is a sintax output that has the taxonomy attached to the end of each line
ASV <- read.delim('zh_ch_us_trimmed_qfiltered_renamed_ZOTU_Count_Sintax.txt', sep = "\t",header=TRUE)

# split out ASV Table
tax_col <- "Consensus.Lineage"
asv_table <- ASV %>%
  select(-all_of(tax_col)) %>%
  column_to_rownames("X.OTU.ID")
asv_table <- as.data.frame(t(asv_table))

# Get Taxonomy Table
tax <- ASV %>%
  select(`X.OTU.ID`, all_of(tax_col)) %>%
  dplyr::rename(OTU = X.OTU.ID) %>%
  separate(`Consensus.Lineage`,
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";", fill = "right", remove = TRUE) %>%
  mutate(across(-`OTU`, ~str_trim(str_remove_all(., "[a-z]__")))) %>%
  column_to_rownames("OTU")

### Sequences to ASV/taxonomy table
seqs <- Biostrings::readDNAStringSet("zh_ch_us_trimmed_qfiltered_renamed_ZOTU.fa")
# Add to taxonomy table and ASV table
tax$Sequence <- as.character(seqs[rownames(tax)])
rownames(tax) <- tax$Sequence
colnames(asv_table) <- tax$Sequence
head(tax)


####### Make Phyloseq object #######
bphy <- phyloseq(otu_table(asv_table,taxa_are_rows=FALSE),tax_table(as.matrix(tax)),sample_data(smd))
ntaxa(bphy) # 17653 ASVs

# Save original phyloseq object:
saveRDS(bphy,'intermediate_data/phyloseq_b_original.RDS')

####### Remove bad ASVs and record usable reads for each sample #######
# Remove ASVs unclassified at Kingdom level + mitochondria + chloroplast
bphy.nobadASVs<-subset_taxa(bphy,Kingdom%in%c('Bacteria','Archaea') & Family!='mitochondria' & Class!='Chloroplast')
ntaxa(bphy.nobadASVs) # 17653 ASVs remain
sample_data(bphy.nobadASVs)$UsableReads <-sample_sums(bphy.nobadASVs)

####### Calculate alpha diversity before thresholding: #######
# Calculate alpha diversity and add to data frame:
smd <- estimate_richness(bphy.nobadASVs,measures=c('Shannon','ACE','Observed')) %>%
  as.data.frame() %>% merge(.,sample_data(bphy.nobadASVs),by='row.names') %>%
  column_to_rownames('SampleID') %>% dplyr::rename(SampleID=Row.names)
sample_data(bphy.nobadASVs) <- smd  # update in phyloseq object

####### Across-sample thresholding: 5x25 (throw out "non-reproducible" ASVs) #######
threshold<-kOverA(5,A=25) # set threshold values (require k samples with A reads)
bphy.nobadASVs.thresholded <- filter_taxa(bphy.nobadASVs,threshold,TRUE)
ntaxa(bphy.nobadASVs.thresholded) # 2048 ASVs remain
ntaxa(bphy.nobadASVs) # 17653 ASVs originally
# What proportion of original reads remain?
sum(taxa_sums(bphy.nobadASVs.thresholded))/(sum(taxa_sums(bphy.nobadASVs)))
# 76.5 % of reads remain after thresholding

####### Remove samples with <500 usable reads #######
bphy.nobadASVs.thresholded.highcoverage <- subset_samples(bphy.nobadASVs.thresholded,UsableReads>=500) %>% prune_taxa(taxa_sums(.)>0,.)
nsamples(bphy.nobadASVs.thresholded.highcoverage) # 258 samples remain

####### Save total number of observations in sample metadata #######
# This is the nuisance variable we will include in statistical models to control for sequencing depth # 
sample_data(bphy.nobadASVs.thresholded.highcoverage)$Obs <- sample_sums(bphy.nobadASVs.thresholded.highcoverage) # add up "good" ASV observations in each sample
sample_data(bphy.nobadASVs.thresholded.highcoverage)$logObs <- log(sample_data(bphy.nobadASVs.thresholded.highcoverage)$Obs) # store natural log of "good" ASV observations

# how many observations in main dataset?
mean(sample_data(bphy.nobadASVs.thresholded.highcoverage)$Obs) # mean = 33003 observations
sd(sample_data(bphy.nobadASVs.thresholded.highcoverage)$Obs) # sd = 13714

# how many ASVs?
ntaxa(bphy.nobadASVs.thresholded.highcoverage) # 2048 bASVs

####### Extract sample data for "final" dataset #######
hyb.bact <- bphy.nobadASVs.thresholded.highcoverage # rename final phyloseq object
smd.bact <- as(sample_data(hyb.bact),'data.frame')
smd.bact$SampleID <- row.names(smd.bact) # Store Sample IDs as a column
summary(smd.bact) # Make sure everything looks right

####### Relabel ASVs for convenience #######
tax <- as(tax_table(hyb.bact),'matrix')
ASV <- as(otu_table(hyb.bact),'matrix')

# Make sure ASV sequences match in taxonomy & ASV tables
setdiff(rownames(tax),colnames(ASV)) # no differences
setdiff(colnames(ASV),rownames(tax)) # no differences

# For security, make a named list of sequences and ASV_IDs to use in relabeling:
ASVlist <- paste0(ASVprefix,'ASV_',1:ncol(ASV)) # list of ASV_IDs
names(ASVlist) <- colnames(ASV) # name ASV_IDs by their sequences

# Tax table: add column to hold sequences, THEN relabel rows
tax <- transform(tax, Sequence = row.names(tax))
for (i in 1:nrow(tax)) {
  seq <- tax[['Sequence']][i]
  row.names(tax)[i] <- ASVlist[[seq]] # rename rows by matching ASV_IDs from named list
  print(paste0('relabeling ',ASVprefix,'ASV_',i,'...')) }

# ASV table: relabel columns
for (i in 1:ncol(ASV)) { 
  seq <- colnames(ASV)[i]
  colnames(ASV)[i] <- ASVlist[[seq]] # rename columns by matching ASV_IDs from named list
  print(paste0('relabeling ',ASVprefix,'ASV_',i,' out of ',ncol(ASV),'...'))
}

# Add relabeled tax & OTU tables back into phyloseq object:
hyb.bact <- phyloseq(otu_table(ASV,taxa_are_rows=FALSE),tax_table(as(tax,'matrix')))
sample_data(hyb.bact) <- smd.bact

# Save cleaned dataset to file:
saveRDS(hyb.bact,'intermediate_data/Cadot_phyloseq_b_asv.RDS')



#### Generate ASV Lists ####
# Load data (non vst transformed)
hyb.bact.asv <- readRDS('intermediate_data/Cadot_phyloseq_b_asv.RDS')
# metadata
summary(sample_data(hyb.bact.asv))

### Hu 2018 
#rhizosphere
rhizo.ps <- subset_samples(hyb.bact.asv, location  == "Changins" & compartment == "rhizo")
# Remove taxa with zero abundance
rhizo.ps <- prune_taxa(taxa_sums(rhizo.ps) > 0, rhizo.ps)
# taxa dataframes
rhizo.tax <- data.frame(rhizo.ps@tax_table)
# Create a taxonomy string for each row
rhizo.tax$taxonomy <- apply(rhizo.tax[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")], 1, function(x) paste(na.omit(x), collapse = "_"))
# Combine ASV ID and taxonomy into a FASTA header
rhizo.tax$header <- paste0(">", rownames(rhizo.tax), "-", "RHIZOSPHERE", "-", "Switzerland_Changins", " ", rhizo.tax$taxonomy)
# Replace "Cadot" with "Hu" as the samples orginate from a different study
rhizo.tax$header <- gsub("^>Cadot", ">Hu", rhizo.tax$header)
# Interleave header and write sequence lines
rhizo_fasta_lines <- c(rbind(rhizo.tax$header, rhizo.tax$Sequence))
writeLines(rhizo_fasta_lines, "Rhizosphere_ASVs_with_taxonomy_Hu2018.fasta")

# soil
soil.ps <- subset_samples(hyb.bact.asv, location  == "Changins" & compartment == "soil")
# Remove taxa with zero abundance
soil.ps <- prune_taxa(taxa_sums(soil.ps) > 0, soil.ps)
# taxa dataframes
soil.tax <- data.frame(soil.ps@tax_table)
# Create a taxonomy string for each row
soil.tax$taxonomy <- apply(soil.tax[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")], 1, function(x) paste(na.omit(x), collapse = "_"))
# Combine ASV ID and taxonomy into a FASTA header
soil.tax$header <- paste0(">", rownames(soil.tax), "-", "SOIL", "-", "Switzerland_Changins", " ", soil.tax$taxonomy)
# Replace "Cadot" with "Hu" as the samples orginate from a different study
soil.tax$header <- gsub("^>Cadot", ">Hu", soil.tax$header)
# Interleave header and write sequence lines
soil_fasta_lines <- c(rbind(soil.tax$header, soil.tax$Sequence))
writeLines(soil_fasta_lines, "Soil_ASVs_with_taxonomy_Hu2018.fasta")

# root endosphere
root.ps <- subset_samples(hyb.bact.asv, location  == "Changins" & compartment == "root")
# Remove taxa with zero abundance
root.ps <- prune_taxa(taxa_sums(root.ps) > 0, root.ps)
# taxa dataframes
root.tax <- data.frame(root.ps@tax_table)
# Create a taxonomy string for each row
root.tax$taxonomy <- apply(root.tax[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")], 1, function(x) paste(na.omit(x), collapse = "_"))
# Combine ASV ID and taxonomy into a FASTA header
root.tax$header <- paste0(">", rownames(root.tax), "-", "ROOT", "-", "Switzerland_Changins", " ", root.tax$taxonomy)
# Replace "Cadot" with "Hu" as the samples orginate from a different study
root.tax$header <- gsub("^>Cadot", ">Hu", root.tax$header)
# Interleave header and write sequence lines
root_fasta_lines <- c(rbind(root.tax$header, root.tax$Sequence))
writeLines(root_fasta_lines, "Root_ASVs_with_taxonomy_Hu2018.fasta")


## Cadot 2021 - Ithaca Feild Site USA-NY

#rhizosphere
rhizo.ps <- subset_samples(hyb.bact.asv, location  == "Aurora" & compartment == "rhizo")
# Remove taxa with zero abundance
rhizo.ps <- prune_taxa(taxa_sums(rhizo.ps) > 0, rhizo.ps)
# taxa dataframes
rhizo.tax <- data.frame(rhizo.ps@tax_table)
# Create a taxonomy string for each row
rhizo.tax$taxonomy <- apply(rhizo.tax[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")], 1, function(x) paste(na.omit(x), collapse = "_"))
# Combine ASV ID and taxonomy into a FASTA header
rhizo.tax$header <- paste0(">", rownames(rhizo.tax), "-", "RHIZOSPHERE", "-", "USA_NY", " ", rhizo.tax$taxonomy)
# Interleave header and write sequence lines
rhizo_fasta_lines <- c(rbind(rhizo.tax$header, rhizo.tax$Sequence))
writeLines(rhizo_fasta_lines, "Rhizosphere_ASVs_with_taxonomy_Cadot2021_Aurora.fasta")

# soil
soil.ps <- subset_samples(hyb.bact.asv, location  == "Aurora" & compartment == "soil")
# Remove taxa with zero abundance
soil.ps <- prune_taxa(taxa_sums(soil.ps) > 0, soil.ps)
# taxa dataframes
soil.tax <- data.frame(soil.ps@tax_table)
# Create a taxonomy string for each row
soil.tax$taxonomy <- apply(soil.tax[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")], 1, function(x) paste(na.omit(x), collapse = "_"))
# Combine ASV ID and taxonomy into a FASTA header
soil.tax$header <- paste0(">", rownames(soil.tax), "-", "SOIL", "-", "USA_NY", " ", soil.tax$taxonomy)
# Interleave header and write sequence lines
soil_fasta_lines <- c(rbind(soil.tax$header, soil.tax$Sequence))
writeLines(soil_fasta_lines, "Soil_ASVs_with_taxonomy_Cadot2021_Aurora.fasta")

# root endosphere
root.ps <- subset_samples(hyb.bact.asv, location  == "Aurora" & compartment == "root")
# Remove taxa with zero abundance
root.ps <- prune_taxa(taxa_sums(root.ps) > 0, root.ps)
# taxa dataframes
root.tax <- data.frame(root.ps@tax_table)
# Create a taxonomy string for each row
root.tax$taxonomy <- apply(root.tax[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")], 1, function(x) paste(na.omit(x), collapse = "_"))
# Combine ASV ID and taxonomy into a FASTA header
root.tax$header <- paste0(">", rownames(root.tax), "-", "ROOT", "-", "USA_NY", " ", root.tax$taxonomy)
# Interleave header and write sequence lines
root_fasta_lines <- c(rbind(root.tax$header, root.tax$Sequence))
writeLines(root_fasta_lines, "Root_ASVs_with_taxonomy_Cadot2021_Aurora.fasta")


## Cadot 2021 - Reckenholz Feild Site Switzerland_Nyon

#rhizosphere
rhizo.ps <- subset_samples(hyb.bact.asv, location  == "Reckenholz" & compartment == "rhizo")
# Remove taxa with zero abundance
rhizo.ps <- prune_taxa(taxa_sums(rhizo.ps) > 0, rhizo.ps)
# taxa dataframes
rhizo.tax <- data.frame(rhizo.ps@tax_table)
# Create a taxonomy string for each row
rhizo.tax$taxonomy <- apply(rhizo.tax[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")], 1, function(x) paste(na.omit(x), collapse = "_"))
# Combine ASV ID and taxonomy into a FASTA header
rhizo.tax$header <- paste0(">", rownames(rhizo.tax), "-", "RHIZOSPHERE", "-", "Switzerland_Nyon", " ", rhizo.tax$taxonomy)
# Interleave header and write sequence lines
rhizo_fasta_lines <- c(rbind(rhizo.tax$header, rhizo.tax$Sequence))
writeLines(rhizo_fasta_lines, "Rhizosphere_ASVs_with_taxonomy_Cadot2021_Reckenholz.fasta")

# soil
soil.ps <- subset_samples(hyb.bact.asv, location  == "Reckenholz" & compartment == "soil")
# Remove taxa with zero abundance
soil.ps <- prune_taxa(taxa_sums(soil.ps) > 0, soil.ps)
# taxa dataframes
soil.tax <- data.frame(soil.ps@tax_table)
# Create a taxonomy string for each row
soil.tax$taxonomy <- apply(soil.tax[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")], 1, function(x) paste(na.omit(x), collapse = "_"))
# Combine ASV ID and taxonomy into a FASTA header
soil.tax$header <- paste0(">", rownames(soil.tax), "-", "SOIL", "-", "Switzerland_Nyon", " ", soil.tax$taxonomy)
# Interleave header and write sequence lines
soil_fasta_lines <- c(rbind(soil.tax$header, soil.tax$Sequence))
writeLines(soil_fasta_lines, "Soil_ASVs_with_taxonomy_Cadot2021_Reckenholz.fasta")

# root endosphere
root.ps <- subset_samples(hyb.bact.asv, location  == "Reckenholz" & compartment == "root")
# Remove taxa with zero abundance
root.ps <- prune_taxa(taxa_sums(root.ps) > 0, root.ps)
# taxa dataframes
root.tax <- data.frame(root.ps@tax_table)
# Create a taxonomy string for each row
root.tax$taxonomy <- apply(root.tax[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")], 1, function(x) paste(na.omit(x), collapse = "_"))
# Combine ASV ID and taxonomy into a FASTA header
root.tax$header <- paste0(">", rownames(root.tax), "-", "ROOT", "-", "Switzerland_Nyon", " ", root.tax$taxonomy)
# Interleave header and write sequence lines
root_fasta_lines <- c(rbind(root.tax$header, root.tax$Sequence))
writeLines(root_fasta_lines, "Root_ASVs_with_taxonomy_Cadot2021_Reckenholz.fasta")