# Analysis of "Heritable microbiome variation is correlated with source environment
# in locally adapted maize varieties"
# https://doi.org/10.1038/s41477-024-01654-7
## Code by: Joel Swift 2026


####### Packages #######
library(tidyverse)
library(phyloseq)
library(lmerTest)
library(genefilter)

####### Set parameters #######
ASVprefix <- 'He_' # for bacterial ASVs

####### Metadata file #######
smd <- readRDS("./metadata/10.1038-s41477-024-01654-7_metadata.RDS")

####### Load ASV table and taxa info #######
ASV <- readRDS('merged.16S.seqtab.nochim.RDS')
tax <- read.delim('taxa_bacteria_He_etal.txt',sep='\t',header=TRUE)

####### Make Phyloseq object #######
bphy <- phyloseq(otu_table(ASV,taxa_are_rows=FALSE),tax_table(as.matrix(tax)),sample_data(smd))
ntaxa(bphy) # 475646 ASVs

# Save original phyloseq object:
saveRDS(bphy,'intermediate_data/phyloseq_b_original.RDS')
####### Look for plant contamination: #######
plantASVs <- subset_taxa(bphy,Family=='mitochondria' | str_detect(Class,'Chloroplast')) 
sum(taxa_sums(plantASVs)) # 127665188 observations
sum(taxa_sums(bphy)) # out of 687809128

####### Remove bad ASVs and record usable reads for each sample #######
# Remove ASVs unclassified at Kingdom level + mitochondria + chloroplast
bphy.nobadASVs<-subset_taxa(bphy,Kingdom%in%c('Bacteria','Archaea') & Family!='mitochondria' & Class!='Chloroplast')
ntaxa(bphy.nobadASVs) # 181161 ASVs remain
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
ntaxa(bphy.nobadASVs.thresholded) # 3742 ASVs remain
ntaxa(bphy.nobadASVs) # 181161 ASVs originally
# What proportion of original reads remain?
sum(taxa_sums(bphy.nobadASVs.thresholded))/(sum(taxa_sums(bphy.nobadASVs)))
# 96.0 % of reads remain after thresholding

####### Remove samples with <500 usable reads #######
bphy.nobadASVs.thresholded.highcoverage <- subset_samples(bphy.nobadASVs.thresholded,UsableReads>=500) %>% prune_taxa(taxa_sums(.)>0,.)
nsamples(bphy.nobadASVs.thresholded.highcoverage) # 1045 samples remain

####### Save total number of observations in sample metadata #######
# This is the nuisance variable we will include in statistical models to control for sequencing depth # 
sample_data(bphy.nobadASVs.thresholded.highcoverage)$Obs <- sample_sums(bphy.nobadASVs.thresholded.highcoverage) # add up "good" ASV observations in each sample
sample_data(bphy.nobadASVs.thresholded.highcoverage)$logObs <- log(sample_data(bphy.nobadASVs.thresholded.highcoverage)$Obs) # store natural log of "good" ASV observations

# how many observations in main dataset?
mean(sample_data(bphy.nobadASVs.thresholded.highcoverage)$Obs) # mean = 147562 observations
sd(sample_data(bphy.nobadASVs.thresholded.highcoverage)$Obs) # sd = 92629

# how many ASVs?
ntaxa(bphy.nobadASVs.thresholded.highcoverage) # 3742 bASVs

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
saveRDS(hyb.bact,'intermediate_data/He2024_phyloseq_b_asv.RDS')


#### Generate ASV Lists ####
# Load data (non vst transformed)
hyb.bact.asv <- readRDS('intermediate_data/He2024_phyloseq_b_asv.RDS')
# metadata
summary(sample_data(hyb.bact.asv))

### Export fasta
# Subset by sample type
rhizo.ps <- subset_samples(hyb.bact.asv, tissue  == "Rhizosphere")
# Remove taxa with zero abundance
rhizo.ps <- prune_taxa(taxa_sums(rhizo.ps) > 0, rhizo.ps)
# taxa dataframes
rhizo.tax <- data.frame(rhizo.ps@tax_table)
# Create a taxonomy string for each row
rhizo.tax$taxonomy <- apply(rhizo.tax[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")], 1, function(x) paste(na.omit(x), collapse = "_"))
# Combine ASV ID and taxonomy into a FASTA header
rhizo.tax$header <- paste0(">", rownames(rhizo.tax), "-", "RHIZOSPHERE", "-", "Germany_Cologne", " ", rhizo.tax$taxonomy)
# Interleave header and write sequence lines
rhizo_fasta_lines <- c(rbind(rhizo.tax$header, rhizo.tax$Sequence))
writeLines(rhizo_fasta_lines, "Rhizosphere_ASVs_with_taxonomy_He2024.fasta")

# Subset by sample type
Root.ps <- subset_samples(hyb.bact.asv, tissue  == "Root")
# Remove taxa with zero abundance
Root.ps <- prune_taxa(taxa_sums(Root.ps) > 0, Root.ps)
# taxa dataframes
Root.tax <- data.frame(Root.ps@tax_table)
# Create a taxonomy string for each row
Root.tax$taxonomy <- apply(Root.tax[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")], 1, function(x) paste(na.omit(x), collapse = "_"))
# Combine ASV ID and taxonomy into a FASTA header
Root.tax$header <- paste0(">", rownames(Root.tax), "-", "ROOT", "-", "Germany_Cologne", " ", Root.tax$taxonomy)
# Interleave header and write sequence lines
Root_fasta_lines <- c(rbind(Root.tax$header, Root.tax$Sequence))
writeLines(Root_fasta_lines, "Root_ASVs_with_taxonomy_He2024.fasta")