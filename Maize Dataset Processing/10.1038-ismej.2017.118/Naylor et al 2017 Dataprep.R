# Data processing for "Drought and host selection influence bacterial community dynamics
# in the grass root microbiome"
# https://doi.org/10.1038/ismej.2017.118
## Code by: Joel Swift 2026

####### Packages #######
library(tidyverse)
library(phyloseq)
library(lmerTest)
library(genefilter)

####### Set parameters #######
ASVprefix <- 'Naylor' # for bacterial ASVs

####### Metadata file #######
smd <- readRDS("./metadata/10.1038-ismej.2017.118_metadata.RDS")

####### Load ASV table and taxa info #######
ASV <- readRDS('Naylor_merged.16S.seqtab.nochim.RDS')
tax <- read.delim('Naylor_taxa_bacteria.txt',sep='\t',header=TRUE)

####### Make Phyloseq object #######
bphy <- phyloseq(otu_table(ASV,taxa_are_rows=FALSE),tax_table(as.matrix(tax)),sample_data(smd))
ntaxa(bphy) # 72831 ASVs

# Save original phyloseq object:
saveRDS(bphy,'intermediate_data/phyloseq_b_original.RDS')
####### Look for plant contamination: #######
plantASVs <- subset_taxa(bphy,Family=='mitochondria' | str_detect(Class,'Chloroplast')) 
sum(taxa_sums(plantASVs)) # 72652 observations
sum(taxa_sums(bphy)) # out of 16484002

####### Remove bad ASVs and record usable reads for each sample #######
# Remove ASVs unclassified at Kingdom level + mitochondria + chloroplast
bphy.nobadASVs<-subset_taxa(bphy,Kingdom%in%c('Bacteria','Archaea') & Family!='mitochondria' & Class!='Chloroplast')
ntaxa(bphy.nobadASVs) # 56071 ASVs remain
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
ntaxa(bphy.nobadASVs.thresholded) # 2978 ASVs remain
ntaxa(bphy.nobadASVs) # 56071 ASVs originally
# What proportion of original reads remain?
sum(taxa_sums(bphy.nobadASVs.thresholded))/(sum(taxa_sums(bphy.nobadASVs)))
# 86.0 % of reads remain after thresholding

####### Remove samples with <500 usable reads #######
bphy.nobadASVs.thresholded.highcoverage <- subset_samples(bphy.nobadASVs.thresholded,UsableReads>=500) %>% prune_taxa(taxa_sums(.)>0,.)
nsamples(bphy.nobadASVs.thresholded.highcoverage) # 495 samples remain

####### Save total number of observations in sample metadata #######
# This is the nuisance variable we will include in statistical models to control for sequencing depth # 
sample_data(bphy.nobadASVs.thresholded.highcoverage)$Obs <- sample_sums(bphy.nobadASVs.thresholded.highcoverage) # add up "good" ASV observations in each sample
sample_data(bphy.nobadASVs.thresholded.highcoverage)$logObs <- log(sample_data(bphy.nobadASVs.thresholded.highcoverage)$Obs) # store natural log of "good" ASV observations

# how many observations in main dataset?
mean(sample_data(bphy.nobadASVs.thresholded.highcoverage)$Obs) # mean = 23899 observations
sd(sample_data(bphy.nobadASVs.thresholded.highcoverage)$Obs) # sd = 15305

# how many ASVs?
ntaxa(bphy.nobadASVs.thresholded.highcoverage) # 2978 bASVs

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
ASVlist <- paste0(ASVprefix,'_ASV_',1:ncol(ASV)) # list of ASV_IDs
names(ASVlist) <- colnames(ASV) # name ASV_IDs by their sequences

# Tax table: add column to hold sequences, THEN relabel rows
tax <- transform(tax, Sequence = row.names(tax))
for (i in 1:nrow(tax)) {
  seq <- tax[['Sequence']][i]
  row.names(tax)[i] <- ASVlist[[seq]] # rename rows by matching ASV_IDs from named list
  print(paste0('relabeling ',ASVprefix,'_ASV_',i,'...')) }

# ASV table: relabel columns
for (i in 1:ncol(ASV)) { 
  seq <- colnames(ASV)[i]
  colnames(ASV)[i] <- ASVlist[[seq]] # rename columns by matching ASV_IDs from named list
  print(paste0('relabeling ',ASVprefix,'_ASV_',i,' out of ',ncol(ASV),'...'))
}

# Add relabeled tax & OTU tables back into phyloseq object:
hyb.bact <- phyloseq(otu_table(ASV,taxa_are_rows=FALSE),tax_table(as(tax,'matrix')))
sample_data(hyb.bact) <- smd.bact

# Save cleaned dataset to file:
saveRDS(hyb.bact,'intermediate_data/phyloseq_b_asv.RDS')


#### Generate ASV Lists ####
# Load data (non vst transformed)
hyb.bact.asv <- readRDS('intermediate_data/phyloseq_b_asv.RDS')

### Export fasta
# Subset by plant_body_site and isolation source
zea.ps <- subset_samples(hyb.bact.asv, isolation_source  == "Zea mays")
rhizo.ps <- subset_samples(zea.ps, plant_body_site  == "Rhizosphere")
# Remove taxa with zero abundance
rhizo.ps <- prune_taxa(taxa_sums(rhizo.ps) > 0, rhizo.ps)
# taxa dataframes
rhizo.ps <- data.frame(rhizo.ps@tax_table)
# Create a taxonomy string for each row
rhizo.ps$taxonomy <- apply(rhizo.ps[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")], 1, function(x) paste(na.omit(x), collapse = "_"))
# Combine ASV ID and taxonomy into a FASTA header
rhizo.ps$header <- paste0(">", rownames(rhizo.ps), "-", "RHIZOSPHERE", "-", "USA_CA", " ", rhizo.ps$taxonomy)
# Interleave header and write sequence lines
rhizo_fasta_lines <- c(rbind(rhizo.ps$header, rhizo.ps$Sequence))
writeLines(rhizo_fasta_lines, "Rhizosphere_ASVs_with_taxonomy_Naylor2017.fasta")

# Subset by plant_body_site and isolation source
zea.ps <- subset_samples(hyb.bact.asv, isolation_source  == "Zea mays")
root.ps <- subset_samples(zea.ps, plant_body_site  == "Root")
# Remove taxa with zero abundance
root.ps <- prune_taxa(taxa_sums(root.ps) > 0, root.ps)
# taxa dataframes
root.ps <- data.frame(root.ps@tax_table)
# Create a taxonomy string for each row
root.ps$taxonomy <- apply(root.ps[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")], 1, function(x) paste(na.omit(x), collapse = "_"))
# Combine ASV ID and taxonomy into a FASTA header
root.ps$header <- paste0(">", rownames(root.ps), "-", "ROOT", "-", "USA_CA", " ", root.ps$taxonomy)
# Interleave header and write sequence lines
root_fasta_lines <- c(rbind(root.ps$header, root.ps$Sequence))
writeLines(root_fasta_lines, "Root_ASVs_with_taxonomy_Naylor2017.fasta")