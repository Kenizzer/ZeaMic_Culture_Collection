# Data processing for "Plant secondary metabolite-dependent plant-soil feedbacks can 
# improve crop yield in the field"

# https://doi.org/10.7554/eLife.84988


####### Packages #######
library(tidyverse)
library(phyloseq)
library(lmerTest)
library(genefilter)


####### Set parameters #######
ASVprefix <- 'Gfeller_' # for bacterial ASVs

####### Metadata file #######
smd <- readxl::read_xlsx("./metadata/metadata_Gfeller2023.xlsx", sheet = 2)
smd$sample_name <- smd$SampleID
smd <- as.data.frame(smd)


# make the sample names match the ASV table
colnames(smd)
# Run + bac_f_bc + bac_r_bc
# Example runBE10_F66_R43
smd$SampleID <- paste0(smd$Run, "_", smd$bac_f_bc, "_", smd$bac_r_bc)
rownames(smd) <- smd$SampleID

####### Load ASV table and taxa info #######
ASV <- read.delim('bacteria_DAT100.tab', sep = "\t",header=TRUE)
rownames(ASV) <- gsub("-", "_", rownames(ASV))
tax <- read.delim('bacteria_TAXA100.tab',sep='\t',header=TRUE)

### Sequences to ASV/taxonomy table
seqs <- Biostrings::readDNAStringSet("bacteria_SEQ100.fasta")
# Add to taxonomy table and ASV table
tax$Sequence <- as.character(seqs[rownames(tax)])
rownames(tax) <- tax$Sequence
colnames(ASV) <- tax$Sequence
head(tax)


####### Make Phyloseq object #######
bphy <- phyloseq(otu_table(ASV,taxa_are_rows=FALSE),tax_table(as.matrix(tax)),sample_data(smd))
ntaxa(bphy) # 12748 ASVs

# Save original phyloseq object:
saveRDS(bphy,'intermediate_data/phyloseq_b_original.RDS')

####### Remove bad ASVs and record usable reads for each sample #######
# Remove ASVs unclassified at Kingdom level + mitochondria + chloroplast
bphy.nobadASVs<-subset_taxa(bphy,Kingdom%in%c('Bacteria','Archaea') & Family!='mitochondria' & Class!='Chloroplast')
ntaxa(bphy.nobadASVs) # 8765 ASVs remain
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
ntaxa(bphy.nobadASVs.thresholded) # 1166 ASVs remain
ntaxa(bphy.nobadASVs) # 8765 ASVs originally
# What proportion of original reads remain?
sum(taxa_sums(bphy.nobadASVs.thresholded))/(sum(taxa_sums(bphy.nobadASVs)))
# 87.9 % of reads remain after thresholding

####### Remove samples with <500 usable reads #######
bphy.nobadASVs.thresholded.highcoverage <- subset_samples(bphy.nobadASVs.thresholded,UsableReads>=500) %>% prune_taxa(taxa_sums(.)>0,.)
nsamples(bphy.nobadASVs.thresholded.highcoverage) # 248 samples remain

####### Save total number of observations in sample metadata #######
# This is the nuisance variable we will include in statistical models to control for sequencing depth # 
sample_data(bphy.nobadASVs.thresholded.highcoverage)$Obs <- sample_sums(bphy.nobadASVs.thresholded.highcoverage) # add up "good" ASV observations in each sample
sample_data(bphy.nobadASVs.thresholded.highcoverage)$logObs <- log(sample_data(bphy.nobadASVs.thresholded.highcoverage)$Obs) # store natural log of "good" ASV observations

# how many observations in main dataset?
mean(sample_data(bphy.nobadASVs.thresholded.highcoverage)$Obs) # mean = 11689 observations
sd(sample_data(bphy.nobadASVs.thresholded.highcoverage)$Obs) # sd = 5263

# how many ASVs?
ntaxa(bphy.nobadASVs.thresholded.highcoverage) # 1166 bASVs

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
saveRDS(hyb.bact,'intermediate_data/Gfeller2023_phyloseq_b_asv.RDS')






#### Generate ASV Lists ####
# Load data (non vst transformed)
hyb.bact.asv <- readRDS('intermediate_data/Gfeller2023_phyloseq_b_asv.RDS')
# metadata
summary(sample_data(hyb.bact.asv))

### Export fasta
#rhizosphere
# Subset by species
rhizo.ps <- subset_samples(hyb.bact.asv, Species  == "maize" & SampleType == "rhizosphere")
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
writeLines(rhizo_fasta_lines, "Rhizosphere_ASVs_with_taxonomy_Gfeller2023.fasta")

# soil
soil.ps <- subset_samples(hyb.bact.asv, Species  == "maize" & SampleType == "soil")
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
writeLines(soil_fasta_lines, "Soil_ASVs_with_taxonomy_Gfeller2023.fasta")

# root endosphere
root.ps <- subset_samples(hyb.bact.asv, Species  == "maize" & SampleType == "root")
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
writeLines(root_fasta_lines, "Root_ASVs_with_taxonomy_Gfeller2023.fasta")