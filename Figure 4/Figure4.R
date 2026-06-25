# Abundance Stats from ASV lists
## Merging datasets together to identify new additions for the collection
## Code by: Joel Swift 2026

# Load packages and set ggplot theme
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggpubr)
theme_set(theme_classic())

# Function to extract occupancy, relative abundance, and tag asv if in ZeaMiC collection
Extract_Occ_RA <- function(phyloseq_obj, Blast_ASV_list){
  # Extract STUDY name
  sd <- as.data.frame(sample_data(phyloseq_obj))
  study_label <- unique(sd$STUDY)
  # Extract OTU table, transpose, re-covert to df, add ASV column for later merging
  otutab <- data.frame(phyloseq_obj@otu_table) 
  otutab <- data.frame(t(otutab)) #Taxa are rows now!
  #otutab$ASV <- rownames(otutab)
  # Extract taxonomy info
  taxonomy <- data.frame(tax_table(phyloseq_obj))
  taxonomy$ASV <- rownames(taxonomy)
  
  # Generate P/A, RA values per taxon
  otu_PA <- 1 * ((otutab > 0) == 1)                                                           # presence-absence data
  otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                                     # occupancy calculation
  otu_rel <- apply(decostand(otutab, method="total", MARGIN=2),1, mean)                       # mean relative abundance
  occ_abun <- data.frame(otu_occ=otu_occ, otu_rel=otu_rel) %>% rownames_to_column('ASV')      # combining occupancy and abundance data frame
  
  # Join them together by ASV col
  occ_rel_tax <- left_join(taxonomy,occ_abun, by = "ASV")
  occ_rel_tax_tagged <- occ_rel_tax %>%
    mutate(Study = study_label,
           InCollection = ifelse(ASV %in% Blast_ASV_list$ASV_id, "Yes", "No"))
  return(occ_rel_tax_tagged)
} 

# Wrapper function to return global and per compartment results for OCC-RA 
Extract_Occ_RA_by_compartment <- function(phyloseq_obj, Blast_ASV_list, compartment_var = "compartment") {
  # Extract metadata and compartment uniq list
  sd <- as.data.frame(sample_data(phyloseq_obj))
  comps <- unique(sd[[compartment_var]])
  # global OCC-RA dataframe, Just wrap the orginal function
  res_global <- Extract_Occ_RA(phyloseq_obj, Blast_ASV_list)
  # Run the function on each comaprtment seperately
  ## Also prune taxa that are not found within the compartment 
  res_comp <- lapply(comps, function(comp) {
    samples_keep <- rownames(sd)[sd[[compartment_var]] == comp]
    ps_sub <- phyloseq_obj %>%
      prune_samples(samples_keep, .) %>%
      prune_taxa(taxa_sums(.) > 0, .)
    df <- Extract_Occ_RA(ps_sub, Blast_ASV_list)
    df$Compartment <- comp
    df})
  # Add global and compartment OCC-RA dataframes to a list, tag name with compartment
  names(res_comp) <- comps
  return(c(list(Global = res_global), res_comp))
}


###### 1) LOAD ALL STUDIES 99 ######
##### 1.1) Meier #####
Meier_bacterial_ps <- readRDS("Phyloseq_objs/Meier_bacterial_phyloseq.RDS")
# Load ASV Lists
Meier_rhizosphere_asvs <- read.csv("ASV_LISTS/99/Meier_rhizosphere_ASVs.csv")
# Add ASV prefix for study (first author first name capitalized)
Meier_rhizosphere_asvs$ASV_id <- paste0("Meier_", Meier_rhizosphere_asvs$ASV_id)

##### 1.2) Wagner ##### 
Wagner_bacterial_ps <- readRDS("Phyloseq_objs/Wagner_bacterial_phyloseq.RDS")
# Load ASV Lists
Wagner_list1 <- read.csv("ASV_LISTS/99/Wagner_phyllosphere_ASVs.csv")
Wagner_list2 <- read.csv("ASV_LISTS/99/Wagner_rhizosphere_ASVs.csv")
# concatenate them together
Wagner_asvs <- bind_rows(Wagner_list1, Wagner_list2) %>% distinct(ASV_id, .keep_all = TRUE)
# Add ASV prefix for study (first author first name capitalized)
Wagner_asvs$ASV_id <- paste0("Wagner_", Wagner_asvs$ASV_id)

##### 1.3) He ##### 
He_bacterial_ps <- readRDS("Phyloseq_objs/He2024_phyloseq_b_asv.RDS")
He_bacterial_ps <- phyloseq::subset_samples(He_bacterial_ps, tissue != "Soil")
# Load ASV Lists
He_list1 <- read.csv("ASV_LISTS/99/He_rhizosphere_ASVs.csv")
He_list2 <- read.csv("ASV_LISTS/99/He_root_ASVs.csv")
# concatenate them together
He_asvs <- bind_rows(He_list1, He_list2) %>% distinct(ASV_id, .keep_all = TRUE)
# Add ASV prefix for study (first author first name capitalized)
He_asvs$ASV_id <- paste0("He_", He_asvs$ASV_id)

##### 1.4) Swift ##### 
Swift_bacterial_ps <- readRDS("Phyloseq_objs/Swift2024_phyloseq_b_asv.RDS")
# change out ASV names to be like other studies
old_names <- taxa_names(Swift_bacterial_ps)
new_names <- paste0("Swift_ASV_", seq_along(old_names))
taxa_names(Swift_bacterial_ps) <- new_names
# Load ASV Lists
Swift_list1 <- read.csv("ASV_LISTS/99/Swift_root_ASVs.csv")
Swift_asvs <- Swift_list1
# Add ASV prefix for study (first author first name capitalized)
Swift_asvs$ASV_id <- paste0("Swift_", Swift_list1$ASV_id)

##### 1.5) Walters ##### 
Walters_bacterial_ps <- readRDS("Phyloseq_objs/Walters2018_phyloseq_b_asv.RDS")
# Load ASV Lists
Walters_list1 <- read.csv("ASV_LISTS/99/Walters_rhizosphere_ASVs.csv")
Walters_asvs <- Walters_list1
# Add ASV prefix for study (first author first name capitalized)
Walters_asvs$ASV_id <- paste0("Walters_", Walters_list1$ASV_id)

##### 1.6) Naylor ##### 
Naylor_bacterial_ps <- readRDS("Phyloseq_objs/Naylor2017_phyloseq_b_asv.RDS")
Naylor_bacterial_ps <- subset_samples(Naylor_bacterial_ps, isolation_source == "Zea_mays")
# Load ASV Lists
Naylor_list1 <- read.csv("ASV_LISTS/99/Naylor_rhizosphere_ASVs.csv")
Naylor_list2 <- read.csv("ASV_LISTS/99/Naylor_root_ASVs.csv")
# concatenate them together
Naylor_asvs <- bind_rows(Naylor_list1, Naylor_list2) %>% distinct(ASV_id, .keep_all = TRUE)
# Add ASV prefix for study (first author first name capitalized)
Naylor_asvs$ASV_id <- paste0("Naylor_", Naylor_asvs$ASV_id)

##### 1.7) Gfeller Cadot ##### 
GfellerCadot_bacterial_ps <- readRDS("Phyloseq_objs/GfellerCadot2023_phyloseq_b_asv.RDS")
GfellerCadot_bacterial_ps <- subset_samples(GfellerCadot_bacterial_ps, Species=="Maize")
# Load ASV Lists
GfellerCadot_list1 <- read.csv("ASV_LISTS/99/GfellerCadot_rhizosphere_ASVs.csv")
GfellerCadot_list2 <- read.csv("ASV_LISTS/99/GfellerCadot_root_ASVs.csv")
# concatenate them together
GfellerCadot_asvs <- bind_rows(GfellerCadot_list1, GfellerCadot_list2) %>% distinct(ASV_id, .keep_all = TRUE)
# Add ASV prefix for study (first author first name capitalized)
GfellerCadot_asvs$ASV_id <- paste0("GfellerCadot_", GfellerCadot_asvs$ASV_id)

##### 1.8) Gfeller ##### 
Gfeller_bacterial_ps <- readRDS("Phyloseq_objs/Gfeller2023_phyloseq_b_asv.RDS")
Gfeller_bacterial_ps <- subset_samples(Gfeller_bacterial_ps, Species=="maize")
# Load ASV Lists
Gfeller_list1 <- read.csv("ASV_LISTS/99/Gfeller_rhizosphere_ASVs.csv")
Gfeller_list2 <- read.csv("ASV_LISTS/99/Gfeller_root_ASVs.csv")
Gfeller_list3 <- read.csv("ASV_LISTS/99/Gfeller_soil_ASVs.csv")
# concatenate them together
Gfeller_asvs <- bind_rows(Gfeller_list1, Gfeller_list2, Gfeller_list3) %>% distinct(ASV_id, .keep_all = TRUE)
# Add ASV prefix for study (first author first name capitalized)
Gfeller_asvs$ASV_id <- paste0("Gfeller_", Gfeller_asvs$ASV_id)

##### 1.9) Cadot #####
Cadot_bacterial_ps <- readRDS("Phyloseq_objs/Cadot_phyloseq_b_asv.RDS")
Cadot_bacterial_ps <- subset_samples(Cadot_bacterial_ps, location  != "Changins")
Cadot_bacterial_ps <- prune_taxa(taxa_sums(Cadot_bacterial_ps) > 0, Cadot_bacterial_ps)
# Load ASV Lists
Cadot_list1 <- read.csv("ASV_LISTS/99/Cadot_rhizosphere_ASVs.csv")
Cadot_list2 <- read.csv("ASV_LISTS/99/Cadot_root_ASVs.csv")
Cadot_list3 <- read.csv("ASV_LISTS/99/Cadot_soil_ASVs.csv")
# concatenate them together
Cadot_asvs <- bind_rows(Cadot_list1, Cadot_list2, Cadot_list3) %>% distinct(ASV_id, .keep_all = TRUE)
# Add ASV prefix for study (first author first name capitalized)
Cadot_asvs$ASV_id <- paste0("Cadot_", Cadot_asvs$ASV_id)

##### 1.10) Hu #####
Cadot_bacterial_ps <- readRDS("Phyloseq_objs/Cadot_phyloseq_b_asv.RDS")
Hu_bacterial_ps <- subset_samples(Cadot_bacterial_ps, location  == "Changins")
Hu_bacterial_ps <- prune_taxa(taxa_sums(Hu_bacterial_ps) > 0, Hu_bacterial_ps)
# Fix ASV names, as they are subset of the Cadot study
colnames(Hu_bacterial_ps@otu_table) <- str_replace(colnames(Hu_bacterial_ps@otu_table), "^Cadot_", "Hu_")
rownames(Hu_bacterial_ps@tax_table) <- str_replace(rownames(Hu_bacterial_ps@tax_table), "^Cadot_", "Hu_")
# Load ASV Lists
Hu_list1 <- read.csv("ASV_LISTS/99/Hu_rhizosphere_ASVs.csv")
Hu_list2 <- read.csv("ASV_LISTS/99/Hu_root_ASVs.csv")
Hu_list3 <- read.csv("ASV_LISTS/99/Hu_soil_ASVs.csv")
# concatenate them together
Hu_asvs <- bind_rows(Hu_list1, Hu_list2, Hu_list3) %>% distinct(ASV_id, .keep_all = TRUE)
# Add ASV prefix for study (first author first name capitalized)
Hu_asvs$ASV_id <- paste0("Hu_", Hu_asvs$ASV_id)

##### 1.11) Cleanup datasets #####
###### CLEAN DATASETS #####
Meier_bacterial_ps@sam_data$STUDY <- "Meier (2022)"
Wagner_bacterial_ps@sam_data$STUDY <- "Wagner (2020)"
He_bacterial_ps@sam_data$STUDY <- "He (2024)"
Swift_bacterial_ps@sam_data$STUDY <- "Swift (2025)"
Walters_bacterial_ps@sam_data$STUDY <- "Walters (2018)"
Naylor_bacterial_ps@sam_data$STUDY <- "Naylor (2017)"
GfellerCadot_bacterial_ps@sam_data$STUDY <- "Gfeller & Cadot (2023)"
Gfeller_bacterial_ps@sam_data$STUDY <- "Gfeller (2023)"
Cadot_bacterial_ps@sam_data$STUDY <- "Cadot (2021)"
Hu_bacterial_ps@sam_data$STUDY <- "Hu (2018)"

sample_data(Meier_bacterial_ps)$compartment <- Meier_bacterial_ps@sam_data$isolation_source
sample_data(Meier_bacterial_ps)$compartment <- str_to_title(sample_data(Meier_bacterial_ps)$compartment)
unique(sample_data(Meier_bacterial_ps)$compartment)

sample_data(Wagner_bacterial_ps)$compartment <- Wagner_bacterial_ps@sam_data$Type
sample_data(Wagner_bacterial_ps)$compartment <- str_to_title(sample_data(Wagner_bacterial_ps)$compartment)
unique(sample_data(Wagner_bacterial_ps)$compartment)

sample_data(He_bacterial_ps)$compartment <- He_bacterial_ps@sam_data$tissue
sample_data(He_bacterial_ps)$compartment <- str_to_title(sample_data(He_bacterial_ps)$compartment)
unique(sample_data(He_bacterial_ps)$compartment)

sample_data(Swift_bacterial_ps)$compartment <- Swift_bacterial_ps@sam_data$Type
sample_data(Swift_bacterial_ps)$compartment <- str_to_title(sample_data(Swift_bacterial_ps)$compartment)
unique(sample_data(Swift_bacterial_ps)$compartment)

sample_data(Walters_bacterial_ps)$compartment <- Walters_bacterial_ps@sam_data$sample_type
sample_data(Walters_bacterial_ps)$compartment <- str_to_title(sample_data(Walters_bacterial_ps)$compartment)
unique(sample_data(Walters_bacterial_ps)$compartment)

sample_data(Naylor_bacterial_ps)$compartment <- Naylor_bacterial_ps@sam_data$plant_body_site
sample_data(Naylor_bacterial_ps)$compartment <- str_to_title(sample_data(Naylor_bacterial_ps)$compartment)
unique(sample_data(Naylor_bacterial_ps)$compartment)

sample_data(GfellerCadot_bacterial_ps)$compartment <- GfellerCadot_bacterial_ps@sam_data$Sample_type
sample_data(GfellerCadot_bacterial_ps)$compartment <- str_to_title(sample_data(GfellerCadot_bacterial_ps)$compartment)
unique(sample_data(GfellerCadot_bacterial_ps)$compartment)

sample_data(Gfeller_bacterial_ps)$compartment <- Gfeller_bacterial_ps@sam_data$SampleType
sample_data(Gfeller_bacterial_ps)$compartment <- str_to_title(sample_data(Gfeller_bacterial_ps)$compartment)
unique(sample_data(Gfeller_bacterial_ps)$compartment)

sample_data(Cadot_bacterial_ps)$compartment <- Cadot_bacterial_ps@sam_data$compartment
sample_data(Cadot_bacterial_ps)$compartment <- sample_data(Cadot_bacterial_ps)$compartment %>% recode("rhizo" = "rhizosphere")
sample_data(Cadot_bacterial_ps)$compartment <- str_to_title(sample_data(Cadot_bacterial_ps)$compartment)
unique(sample_data(Cadot_bacterial_ps)$compartment)

sample_data(Hu_bacterial_ps)$compartment <- Hu_bacterial_ps@sam_data$compartment
sample_data(Hu_bacterial_ps)$compartment <- sample_data(Hu_bacterial_ps)$compartment %>% recode("rhizo" = "rhizosphere")
sample_data(Hu_bacterial_ps)$compartment <- str_to_title(sample_data(Hu_bacterial_ps)$compartment)
unique(sample_data(Hu_bacterial_ps)$compartment)

#### 2) Function Call ####
Meier        <- Extract_Occ_RA_by_compartment(Meier_bacterial_ps, Meier_rhizosphere_asvs)
Wagner       <- Extract_Occ_RA_by_compartment(Wagner_bacterial_ps, Wagner_asvs)
He           <- Extract_Occ_RA_by_compartment(He_bacterial_ps, He_asvs)
Swift        <- Extract_Occ_RA_by_compartment(Swift_bacterial_ps, Swift_asvs)
Walters      <- Extract_Occ_RA_by_compartment(Walters_bacterial_ps, Walters_asvs)
Naylor       <- Extract_Occ_RA_by_compartment(Naylor_bacterial_ps, Naylor_asvs)
GfellerCadot <- Extract_Occ_RA_by_compartment(GfellerCadot_bacterial_ps, GfellerCadot_asvs)
Gfeller      <- Extract_Occ_RA_by_compartment(Gfeller_bacterial_ps, Gfeller_asvs)
Cadot        <- Extract_Occ_RA_by_compartment(Cadot_bacterial_ps, Cadot_asvs)
Hu           <- Extract_Occ_RA_by_compartment(Hu_bacterial_ps, Hu_asvs)

### Function to check results
check_fun <- function(df) {
  # Checks
  rel_sum <- sum(df$otu_rel, na.rm = TRUE) # sum to ~1
  occ_ok <- all(df$otu_occ <= 1, na.rm = TRUE)
  # Create a list report
  list(total_otu_rel   = rel_sum,
       occupancy_valid = occ_ok,
       rows_with_occ_over_1 = df %>% filter(otu_occ > 1))}

check_fun(Meier$Global)
check_fun(Wagner$Global)
check_fun(He$Global)
check_fun(Swift$Global)
check_fun(Walters$Global)
check_fun(Naylor$Global)
check_fun(GfellerCadot$Global)
check_fun(Gfeller$Global)
check_fun(Cadot$Global)
check_fun(Hu$Global)
# ALL GOOD

### Remove species column (Empty in all studies)
Meier$Global$Species <- NULL
He$Global$Species <- NULL
Walters$Global$Species <- NULL
Naylor$Global$Species <- NULL
Hu$Global$Species <- NULL
Cadot$Global$Species <- NULL

# check all column names
all.equal(colnames(Meier$Global),
          colnames(Wagner$Global),
          colnames(He$Global),
          colnames(Swift$Global),
          colnames(Walters$Global),
          colnames(Naylor$Global),
          colnames(GfellerCadot$Global),
          colnames(Gfeller$Global),
          colnames(Cadot$Global),
          colnames(Hu$Global))

### Clean Env
rm(list = c(ls(pattern = "*list*"), ls(pattern = "*_ps*"), ls(pattern = "*asvs*")), old_names, new_names, check_fun)
ls() # Cleaned

#### 3) Grouping study together by region, collapsing identical ASVs/Genera ####
All_studies_99 <- rbind(Meier$Global, Wagner$Global, He$Global, Swift$Global, 
                     Walters$Global, Naylor$Global, GfellerCadot$Global, 
                     Gfeller$Global, Cadot$Global, Hu$Global)

# Rhizo
## Have to remove speices columns from dataframes to get everything to the same dimensions
Meier$Rhizosphere$Species <- NULL
He$Rhizosphere$Species <- NULL
Walters$Rhizosphere$Species <- NULL
Naylor$Rhizosphere$Species <- NULL
Hu$Rhizosphere$Species <- NULL
Cadot$Rhizosphere$Species <- NULL
Rhizo_studies_99 <- rbind(Meier$Rhizosphere, Wagner$Rhizosphere, He$Rhizosphere,
                       Walters$Rhizosphere, Naylor$Rhizosphere, GfellerCadot$Rhizosphere, 
                       Gfeller$Rhizosphere, Cadot$Rhizosphere, Hu$Rhizosphere)
# Root
## Have to remove speices columns from dataframes to get everything to the same dimensions
He$Root$Species <- NULL
Naylor$Root$Species <- NULL
GfellerCadot$Root$Species <- NULL
Swift$Root$Species <- NULL
Gfeller$Root$Species <- NULL
Cadot$Root$Species <- NULL
Hu$Root$Species <- NULL
Root_studies_99 <- rbind(He$Root, Naylor$Root, GfellerCadot$Root, Swift$Root, 
                      Gfeller$Root, Cadot$Root, Hu$Root)
# Soil
Naylor$Soil$Species <- NULL
GfellerCadot$Soil$Species <- NULL
Gfeller$Soil$Species <- NULL
Cadot$Soil$Species <- NULL
Hu$Soil$Species <- NULL
Walters$Soil$Species <- NULL
Soil_studies_99 <- rbind(Naylor$Soil, GfellerCadot$Soil, Gfeller$Soil, Cadot$Soil,
                      Hu$Soil, Walters$Soil)

#leaf
Leaf_studies_99 <- Wagner$Leaf






###### 1) LOAD ALL STUDIES 97 ######
##### 1.1) Meier #####
Meier_bacterial_ps <- readRDS("Phyloseq_objs/Meier_bacterial_phyloseq.RDS")
# Load ASV Lists
Meier_rhizosphere_asvs <- read.csv("ASV_LISTS/97/Meier_rhizosphere_ASVs.csv")
# Add ASV prefix for study (first author first name capitalized)
Meier_rhizosphere_asvs$ASV_id <- paste0("Meier_", Meier_rhizosphere_asvs$ASV_id)

##### 1.2) Wagner ##### 
Wagner_bacterial_ps <- readRDS("Phyloseq_objs/Wagner_bacterial_phyloseq.RDS")
# Load ASV Lists
Wagner_list1 <- read.csv("ASV_LISTS/97/Wagner_phyllosphere_ASVs.csv")
Wagner_list2 <- read.csv("ASV_LISTS/97/Wagner_rhizosphere_ASVs.csv")
# concatenate them together
Wagner_asvs <- bind_rows(Wagner_list1, Wagner_list2) %>% distinct(ASV_id, .keep_all = TRUE)
# Add ASV prefix for study (first author first name capitalized)
Wagner_asvs$ASV_id <- paste0("Wagner_", Wagner_asvs$ASV_id)

##### 1.3) He ##### 
He_bacterial_ps <- readRDS("Phyloseq_objs/He2024_phyloseq_b_asv.RDS")
He_bacterial_ps <- phyloseq::subset_samples(He_bacterial_ps, tissue != "Soil")
# Load ASV Lists
He_list1 <- read.csv("ASV_LISTS/97/He_rhizosphere_ASVs.csv")
He_list2 <- read.csv("ASV_LISTS/97/He_root_ASVs.csv")
# concatenate them together
He_asvs <- bind_rows(He_list1, He_list2) %>% distinct(ASV_id, .keep_all = TRUE)
# Add ASV prefix for study (first author first name capitalized)
He_asvs$ASV_id <- paste0("He_", He_asvs$ASV_id)

##### 1.4) Swift ##### 
Swift_bacterial_ps <- readRDS("Phyloseq_objs/Swift2024_phyloseq_b_asv.RDS")
# change out ASV names to be like other studies
old_names <- taxa_names(Swift_bacterial_ps)
new_names <- paste0("Swift_ASV_", seq_along(old_names))
taxa_names(Swift_bacterial_ps) <- new_names
# Load ASV Lists
Swift_list1 <- read.csv("ASV_LISTS/97/Swift_root_ASVs.csv")
Swift_asvs <- Swift_list1
# Add ASV prefix for study (first author first name capitalized)
Swift_asvs$ASV_id <- paste0("Swift_", Swift_list1$ASV_id)

##### 1.5) Walters ##### 
Walters_bacterial_ps <- readRDS("Phyloseq_objs/Walters2018_phyloseq_b_asv.RDS")
# Load ASV Lists
Walters_list1 <- read.csv("ASV_LISTS/97/Walters_rhizosphere_ASVs.csv")
Walters_asvs <- Walters_list1
# Add ASV prefix for study (first author first name capitalized)
Walters_asvs$ASV_id <- paste0("Walters_", Walters_list1$ASV_id)

##### 1.6) Naylor ##### 
Naylor_bacterial_ps <- readRDS("Phyloseq_objs/Naylor2017_phyloseq_b_asv.RDS")
Naylor_bacterial_ps <- subset_samples(Naylor_bacterial_ps, isolation_source == "Zea_mays")
# Load ASV Lists
Naylor_list1 <- read.csv("ASV_LISTS/97/Naylor_rhizosphere_ASVs.csv")
Naylor_list2 <- read.csv("ASV_LISTS/97/Naylor_root_ASVs.csv")
# concatenate them together
Naylor_asvs <- bind_rows(Naylor_list1, Naylor_list2) %>% distinct(ASV_id, .keep_all = TRUE)
# Add ASV prefix for study (first author first name capitalized)
Naylor_asvs$ASV_id <- paste0("Naylor_", Naylor_asvs$ASV_id)

##### 1.7) Gfeller Cadot ##### 
GfellerCadot_bacterial_ps <- readRDS("Phyloseq_objs/GfellerCadot2023_phyloseq_b_asv.RDS")
GfellerCadot_bacterial_ps <- subset_samples(GfellerCadot_bacterial_ps, Species=="Maize")
# Load ASV Lists
GfellerCadot_list1 <- read.csv("ASV_LISTS/97/GfellerCadot_rhizosphere_ASVs.csv")
GfellerCadot_list2 <- read.csv("ASV_LISTS/97/GfellerCadot_root_ASVs.csv")
# concatenate them together
GfellerCadot_asvs <- bind_rows(GfellerCadot_list1, GfellerCadot_list2) %>% distinct(ASV_id, .keep_all = TRUE)
# Add ASV prefix for study (first author first name capitalized)
GfellerCadot_asvs$ASV_id <- paste0("GfellerCadot_", GfellerCadot_asvs$ASV_id)

##### 1.8) Gfeller ##### 
Gfeller_bacterial_ps <- readRDS("Phyloseq_objs/Gfeller2023_phyloseq_b_asv.RDS")
Gfeller_bacterial_ps <- subset_samples(Gfeller_bacterial_ps, Species=="maize")
# Load ASV Lists
Gfeller_list1 <- read.csv("ASV_LISTS/97/Gfeller_rhizosphere_ASVs.csv")
Gfeller_list2 <- read.csv("ASV_LISTS/97/Gfeller_root_ASVs.csv")
Gfeller_list3 <- read.csv("ASV_LISTS/97/Gfeller_soil_ASVs.csv")
# concatenate them together
Gfeller_asvs <- bind_rows(Gfeller_list1, Gfeller_list2, Gfeller_list3) %>% distinct(ASV_id, .keep_all = TRUE)
# Add ASV prefix for study (first author first name capitalized)
Gfeller_asvs$ASV_id <- paste0("Gfeller_", Gfeller_asvs$ASV_id)

##### 1.9) Cadot #####
Cadot_bacterial_ps <- readRDS("Phyloseq_objs/Cadot_phyloseq_b_asv.RDS")
Cadot_bacterial_ps <- subset_samples(Cadot_bacterial_ps, location  != "Changins")
Cadot_bacterial_ps <- prune_taxa(taxa_sums(Cadot_bacterial_ps) > 0, Cadot_bacterial_ps)
# Load ASV Lists
Cadot_list1 <- read.csv("ASV_LISTS/97/Cadot_rhizosphere_ASVs.csv")
Cadot_list2 <- read.csv("ASV_LISTS/97/Cadot_root_ASVs.csv")
Cadot_list3 <- read.csv("ASV_LISTS/97/Cadot_soil_ASVs.csv")
# concatenate them together
Cadot_asvs <- bind_rows(Cadot_list1, Cadot_list2, Cadot_list3) %>% distinct(ASV_id, .keep_all = TRUE)
# Add ASV prefix for study (first author first name capitalized)
Cadot_asvs$ASV_id <- paste0("Cadot_", Cadot_asvs$ASV_id)

##### 1.10) Hu #####
Cadot_bacterial_ps <- readRDS("Phyloseq_objs/Cadot_phyloseq_b_asv.RDS")
Hu_bacterial_ps <- subset_samples(Cadot_bacterial_ps, location  == "Changins")
Hu_bacterial_ps <- prune_taxa(taxa_sums(Hu_bacterial_ps) > 0, Hu_bacterial_ps)
# Fix ASV names, as they are subset of the Cadot study
colnames(Hu_bacterial_ps@otu_table) <- str_replace(colnames(Hu_bacterial_ps@otu_table), "^Cadot_", "Hu_")
rownames(Hu_bacterial_ps@tax_table) <- str_replace(rownames(Hu_bacterial_ps@tax_table), "^Cadot_", "Hu_")
# Load ASV Lists
Hu_list1 <- read.csv("ASV_LISTS/97/Hu_rhizosphere_ASVs.csv")
Hu_list2 <- read.csv("ASV_LISTS/97/Hu_root_ASVs.csv")
Hu_list3 <- read.csv("ASV_LISTS/97/Hu_soil_ASVs.csv")
# concatenate them together
Hu_asvs <- bind_rows(Hu_list1, Hu_list2, Hu_list3) %>% distinct(ASV_id, .keep_all = TRUE)
# Add ASV prefix for study (first author first name capitalized)
Hu_asvs$ASV_id <- paste0("Hu_", Hu_asvs$ASV_id)

##### 1.11) Cleanup datasets #####
###### CLEAN DATASETS #####
Meier_bacterial_ps@sam_data$STUDY <- "Meier (2022)"
Wagner_bacterial_ps@sam_data$STUDY <- "Wagner (2020)"
He_bacterial_ps@sam_data$STUDY <- "He (2024)"
Swift_bacterial_ps@sam_data$STUDY <- "Swift (2025)"
Walters_bacterial_ps@sam_data$STUDY <- "Walters (2018)"
Naylor_bacterial_ps@sam_data$STUDY <- "Naylor (2017)"
GfellerCadot_bacterial_ps@sam_data$STUDY <- "Gfeller & Cadot (2023)"
Gfeller_bacterial_ps@sam_data$STUDY <- "Gfeller (2023)"
Cadot_bacterial_ps@sam_data$STUDY <- "Cadot (2021)"
Hu_bacterial_ps@sam_data$STUDY <- "Hu (2018)"

sample_data(Meier_bacterial_ps)$compartment <- Meier_bacterial_ps@sam_data$isolation_source
sample_data(Meier_bacterial_ps)$compartment <- str_to_title(sample_data(Meier_bacterial_ps)$compartment)
unique(sample_data(Meier_bacterial_ps)$compartment)

sample_data(Wagner_bacterial_ps)$compartment <- Wagner_bacterial_ps@sam_data$Type
sample_data(Wagner_bacterial_ps)$compartment <- str_to_title(sample_data(Wagner_bacterial_ps)$compartment)
unique(sample_data(Wagner_bacterial_ps)$compartment)

sample_data(He_bacterial_ps)$compartment <- He_bacterial_ps@sam_data$tissue
sample_data(He_bacterial_ps)$compartment <- str_to_title(sample_data(He_bacterial_ps)$compartment)
unique(sample_data(He_bacterial_ps)$compartment)

sample_data(Swift_bacterial_ps)$compartment <- Swift_bacterial_ps@sam_data$Type
sample_data(Swift_bacterial_ps)$compartment <- str_to_title(sample_data(Swift_bacterial_ps)$compartment)
unique(sample_data(Swift_bacterial_ps)$compartment)

sample_data(Walters_bacterial_ps)$compartment <- Walters_bacterial_ps@sam_data$sample_type
sample_data(Walters_bacterial_ps)$compartment <- str_to_title(sample_data(Walters_bacterial_ps)$compartment)
unique(sample_data(Walters_bacterial_ps)$compartment)

sample_data(Naylor_bacterial_ps)$compartment <- Naylor_bacterial_ps@sam_data$plant_body_site
sample_data(Naylor_bacterial_ps)$compartment <- str_to_title(sample_data(Naylor_bacterial_ps)$compartment)
unique(sample_data(Naylor_bacterial_ps)$compartment)

sample_data(GfellerCadot_bacterial_ps)$compartment <- GfellerCadot_bacterial_ps@sam_data$Sample_type
sample_data(GfellerCadot_bacterial_ps)$compartment <- str_to_title(sample_data(GfellerCadot_bacterial_ps)$compartment)
unique(sample_data(GfellerCadot_bacterial_ps)$compartment)

sample_data(Gfeller_bacterial_ps)$compartment <- Gfeller_bacterial_ps@sam_data$SampleType
sample_data(Gfeller_bacterial_ps)$compartment <- str_to_title(sample_data(Gfeller_bacterial_ps)$compartment)
unique(sample_data(Gfeller_bacterial_ps)$compartment)

sample_data(Cadot_bacterial_ps)$compartment <- Cadot_bacterial_ps@sam_data$compartment
sample_data(Cadot_bacterial_ps)$compartment <- sample_data(Cadot_bacterial_ps)$compartment %>% recode("rhizo" = "rhizosphere")
sample_data(Cadot_bacterial_ps)$compartment <- str_to_title(sample_data(Cadot_bacterial_ps)$compartment)
unique(sample_data(Cadot_bacterial_ps)$compartment)

sample_data(Hu_bacterial_ps)$compartment <- Hu_bacterial_ps@sam_data$compartment
sample_data(Hu_bacterial_ps)$compartment <- sample_data(Hu_bacterial_ps)$compartment %>% recode("rhizo" = "rhizosphere")
sample_data(Hu_bacterial_ps)$compartment <- str_to_title(sample_data(Hu_bacterial_ps)$compartment)
unique(sample_data(Hu_bacterial_ps)$compartment)

#### 2) Function Call ####
Meier        <- Extract_Occ_RA_by_compartment(Meier_bacterial_ps, Meier_rhizosphere_asvs)
Wagner       <- Extract_Occ_RA_by_compartment(Wagner_bacterial_ps, Wagner_asvs)
He           <- Extract_Occ_RA_by_compartment(He_bacterial_ps, He_asvs)
Swift        <- Extract_Occ_RA_by_compartment(Swift_bacterial_ps, Swift_asvs)
Walters      <- Extract_Occ_RA_by_compartment(Walters_bacterial_ps, Walters_asvs)
Naylor       <- Extract_Occ_RA_by_compartment(Naylor_bacterial_ps, Naylor_asvs)
GfellerCadot <- Extract_Occ_RA_by_compartment(GfellerCadot_bacterial_ps, GfellerCadot_asvs)
Gfeller      <- Extract_Occ_RA_by_compartment(Gfeller_bacterial_ps, Gfeller_asvs)
Cadot        <- Extract_Occ_RA_by_compartment(Cadot_bacterial_ps, Cadot_asvs)
Hu           <- Extract_Occ_RA_by_compartment(Hu_bacterial_ps, Hu_asvs)

### Function to check results
check_fun <- function(df) {
  # Checks
  rel_sum <- sum(df$otu_rel, na.rm = TRUE) # sum to ~1
  occ_ok <- all(df$otu_occ <= 1, na.rm = TRUE)
  # Create a list report
  list(total_otu_rel   = rel_sum,
       occupancy_valid = occ_ok,
       rows_with_occ_over_1 = df %>% filter(otu_occ > 1))}

check_fun(Meier$Global)
check_fun(Wagner$Global)
check_fun(He$Global)
check_fun(Swift$Global)
check_fun(Walters$Global)
check_fun(Naylor$Global)
check_fun(GfellerCadot$Global)
check_fun(Gfeller$Global)
check_fun(Cadot$Global)
check_fun(Hu$Global)
# ALL GOOD

### Remove species column (Empty in all studies)
Meier$Global$Species <- NULL
He$Global$Species <- NULL
Walters$Global$Species <- NULL
Naylor$Global$Species <- NULL
Hu$Global$Species <- NULL
Cadot$Global$Species <- NULL

# check all column names
all.equal(colnames(Meier$Global),
          colnames(Wagner$Global),
          colnames(He$Global),
          colnames(Swift$Global),
          colnames(Walters$Global),
          colnames(Naylor$Global),
          colnames(GfellerCadot$Global),
          colnames(Gfeller$Global),
          colnames(Cadot$Global),
          colnames(Hu$Global))

### Clean Env
rm(list = c(ls(pattern = "*list*"), ls(pattern = "*_ps*"), ls(pattern = "*asvs*")), old_names, new_names, check_fun, Extract_Occ_RA, Extract_Occ_RA_by_compartment)
ls() # Cleaned

#### 3) Grouping study together by region, collapsing identical ASVs/Genera ####
All_studies_97 <- rbind(Meier$Global, Wagner$Global, He$Global, Swift$Global, 
                        Walters$Global, Naylor$Global, GfellerCadot$Global, 
                        Gfeller$Global, Cadot$Global, Hu$Global)

# Rhizo
## Have to remove speices columns from dataframes to get everything to the same dimensions
Meier$Rhizosphere$Species <- NULL
He$Rhizosphere$Species <- NULL
Walters$Rhizosphere$Species <- NULL
Naylor$Rhizosphere$Species <- NULL
Hu$Rhizosphere$Species <- NULL
Cadot$Rhizosphere$Species <- NULL
Rhizo_studies_97 <- rbind(Meier$Rhizosphere, Wagner$Rhizosphere, He$Rhizosphere,
                          Walters$Rhizosphere, Naylor$Rhizosphere, GfellerCadot$Rhizosphere, 
                          Gfeller$Rhizosphere, Cadot$Rhizosphere, Hu$Rhizosphere)
# Root
## Have to remove speices columns from dataframes to get everything to the same dimensions
He$Root$Species <- NULL
Naylor$Root$Species <- NULL
GfellerCadot$Root$Species <- NULL
Swift$Root$Species <- NULL
Gfeller$Root$Species <- NULL
Cadot$Root$Species <- NULL
Hu$Root$Species <- NULL
Root_studies_97 <- rbind(He$Root, Naylor$Root, GfellerCadot$Root, Swift$Root, 
                         Gfeller$Root, Cadot$Root, Hu$Root)
# Soil
Naylor$Soil$Species <- NULL
GfellerCadot$Soil$Species <- NULL
Gfeller$Soil$Species <- NULL
Cadot$Soil$Species <- NULL
Hu$Soil$Species <- NULL
Walters$Soil$Species <- NULL
Soil_studies_97 <- rbind(Naylor$Soil, GfellerCadot$Soil, Gfeller$Soil, Cadot$Soil,
                         Hu$Soil, Walters$Soil)

#leaf
Leaf_studies_97 <- Wagner$Leaf









#### PLOTTING 99+97 ####

# Collapse the information from All_studies_99 and All_studies_97 into a single dataframe
# capture the ASVs that are found at A) 99 B) at 99+97 and C) not at all
All_studies_combined <- bind_rows(All_studies_99 %>% mutate(Level = "99%"),
                                  All_studies_97 %>% mutate(Level = "97%")) %>%
  group_by(ASV) %>%
  summarise(otu_rel = max(otu_rel, na.rm = TRUE),
            otu_occ = max(otu_occ, na.rm = TRUE),
            Study = first(Study),
            InCollection_99 = ifelse(any(Level == "99%" & InCollection == "Yes"), "Yes", "No"),
            InCollection_97 = ifelse(any(Level == "97%" & InCollection == "Yes"), "Yes", "No"),
            .groups = "drop") %>%
  mutate(Status = case_when(InCollection_99 == "Yes" ~ "Matched at 99%",
                            InCollection_97 == "Yes" ~ "Additional at 97%",
                            TRUE ~ "Not in collection"),
         Status = factor(Status, levels = c("Not in collection",
                                            "Additional at 97%",
                                            "Matched at 99%")))

# OA Scatter plot with nos, 97, 99s plotted by color
OA_99_97_plot <- ggplot(All_studies_combined, aes(x = log10(otu_rel), y = otu_occ, fill = Status, alpha = Status)) + 
  geom_point(data = . %>% filter(Status == "Not in collection"), aes(fill = "Not in collection", alpha = "Not in collection"), shape = 21, size = 1, color = "black") +    ## nos
  geom_point(data = . %>% filter(Status == "Additional at 97%"), aes(fill = "Additional at 97%", alpha = "Additional at 97%"), shape = 21, size = 2, color = "black") + ## 97s
  geom_point(data = . %>% filter(Status == "Matched at 99%"), aes(fill = "Matched at 99%", alpha = "Matched at 99%"), shape = 21, size = 2, color = "black") + ## 99
  scale_fill_manual(values = c("Not in collection" = "white", "Additional at 97%" = "#7850C8", "Matched at 99%" = "#50C878"), breaks = c("Not in collection", "Additional at 97%", "Matched at 99%"), name = "In ZeaMiC Collection") +
  scale_alpha_manual(values = c("Not in collection" = 0.3, "Additional at 97%" = 1, "Matched at 99%"=1), breaks = c("Not in collection", "Additional at 97%", "Matched at 99%"), name = "In ZeaMiC Collection") +
  labs(x = expression(log[10]("RA")), y = "Occupancy") +
  facet_wrap(~Study, nrow = 2) + 
  theme(legend.position = "top", strip.text = element_text(size = 7.5))

# Collapse the information from All_studies_99 and All_studies_97 into a single dataframe
# Here we summarize the relative abundance captured by A) 99 B) 97 - 99 and C) not at all
# So we see the additional RA captured when we have 99 vs 97% sequence similarity

RA_combine <- bind_rows(All_studies_99 %>% mutate(Similarity = "99%"), All_studies_97 %>% mutate(Similarity = "97%")) %>%
  group_by(Study, Similarity, InCollection) %>%
  summarise(total = sum(otu_rel, na.rm = TRUE), .groups = "drop") %>%
  # Collect the RA at 99/97 and get the delta between for plotting
  pivot_wider(names_from = c(Similarity, InCollection), values_from = total, values_fill = 0) %>%
  mutate(in_99    = `99%_Yes`,
         in_97    = `97%_Yes`,
         not_97   = `97%_No`,
         added_97 = in_97 - in_99) %>%
  select(Study, in_99, added_97, not_97) %>%
  pivot_longer(cols = -Study,
               names_to = "Component",
               values_to = "value") %>%
  mutate(Component = factor(Component, levels = c("in_99", "added_97", "not_97"), labels = c("Matched at 99%", "Additional at 97%","Not in collection")))

RA_plot <- ggplot(RA_combine, aes(x = "Collection coverage", y = value, fill = Component)) +
  geom_col(color = "black", width = 0.6) +
  facet_wrap(~Study, nrow = 2) +
  scale_fill_manual(values = c("Matched at 99%" = "#50C878", "Additional at 97%" = "#7850C8", "Not in collection" = "white"), name = NULL) +
  labs(x = NULL, y = "Relative Abundance") +
  theme(legend.position = "top", strip.text = element_text(size = 7.5))

# SAVE
ggsave("./Sequence_99-97_OCC_RA_plot_part1.svg", OA_99_97_plot, width = 7, height = 4.5, unit = "in")
ggsave("./Sequence_99-97_RA_plot_part2.svg", RA_plot, width = 7, height = 4.5, unit = "in")

# Further edits and combination was done in inkscape

