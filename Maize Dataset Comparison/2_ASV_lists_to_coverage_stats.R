# Abundance Stats from ASV lists
# Map ASV lists back to studies and pull the stats
## Code by: Joel Swift 2026

# Load pacakges and set ggplot theme
library(tidyverse)
library(phyloseq)
library(vegan)
theme_set(theme_classic())


# Function to take phyloseq object and isolate ASV list and return a occupancy/abundance dataframe
OCC_Abund_func <- function(ps_object, isolate_asvs_df) {
  # Rarefy the dataset based on the 10% decile (use decile as seed)
  decile_value <- quantile(sample_sums(ps_object), probs = 0.1)
  ps_rare <- rarefy_even_depth(
    ps_object,
    sample.size = decile_value,
    rngseed = decile_value,
    replace = TRUE,
    trimOTUs = TRUE,
    verbose = TRUE)
  otu_tab <- data.frame(otu_table(ps_rare))
  # Occupancy and abundance
  otu_PA <- 1 * ((t(otu_tab) > 0) == 1)
  otu_occ <- rowSums(otu_PA) / ncol(otu_PA)
  otu_rel <- apply(decostand(t(otu_tab), method = "total", MARGIN = 2), 1, mean)
  occ_abun <- data.frame(otu_occ = otu_occ, otu_rel = otu_rel) %>% rownames_to_column("otu")
  # Add taxonomy
  tax_table_df <- data.frame(tax_table(ps_rare)) %>% tibble::rownames_to_column("otu")
  occ_abun_tax <- occ_abun %>% left_join(tax_table_df %>% select(otu, Kingdom:Genus), by = "otu")
  # Tag ASVs that match isolates
  occ_abun_tax <- occ_abun_tax %>% mutate(InCollection = ifelse(otu %in% isolate_asvs_df$ASV_id, "Yes", "No"))
  # Create a summary list for output
  result <- list(occ_abun_tax = occ_abun_tax)
  return(result)
}

# Function to summarize the output in terms of relative abundance captured by ZeaMic collection
summarize_isolate_output <- function(model_output, summarised = FALSE) {
  # Extract the main dataframe
  occ_abun_tax <- model_output$occ_abun_tax
  # Create the base formatted table
  formatted_df <- occ_abun_tax %>%
    select(Genus, otu_occ, otu_rel, InCollection) %>%
    mutate(Abundpercent = round(otu_rel * 100, 2), Occpercent = round(otu_occ * 100, 2)) %>%
    select(Genus, Abundpercent, Occpercent, InCollection) %>%
    arrange(desc(Abundpercent))
  # Optionally summarise by InCollection
  if (summarised) {
    summarised_df <- formatted_df %>%
      group_by(InCollection) %>%
      summarise(TotalAbundance = sum(Abundpercent, na.rm = TRUE)) %>%
      ungroup()
    return(summarised_df)
  } else {
    return(formatted_df)
  }
}


#### 99% ####
###### 1) Meier et al. 2022 ######
# Load phyloseq datasets
Meier_bacterial_ps <- readRDS("Phyloseq_objs/Meier_bacterial_phyloseq.RDS")
# Load ASV Lists
Meier_rhizosphere_asvs <- read.csv("ASV_LISTS/99/Meier_rhizosphere_ASVs.csv")
# Add ASV prefix for study (first author first name capitalized)
Meier_rhizosphere_asvs$ASV_id <- paste0("Meier_", Meier_rhizosphere_asvs$ASV_id)
# run function for O/A table
Meier_df <- OCC_Abund_func(Meier_bacterial_ps, Meier_rhizosphere_asvs)

###### 2) Wagner et al. 2020 ######
# Load phyloseq datasets
Wagner_bacterial_ps <- readRDS("Phyloseq_objs/Wagner_bacterial_phyloseq.RDS")
# Load ASV Lists
Wagner_list1 <- read.csv("ASV_LISTS/99/Wagner_phyllosphere_ASVs.csv")
Wagner_list2 <- read.csv("ASV_LISTS/99/Wagner_rhizosphere_ASVs.csv")
# concatenate them togher
Wagner_asvs <- bind_rows(Wagner_list1, Wagner_list2) %>% distinct(ASV_id, .keep_all = TRUE)
# Add ASV prefix for study (first author first name capitalized)
Wagner_asvs$ASV_id <- paste0("Wagner_", Wagner_asvs$ASV_id)
# run function for O/A table
Wagner_df <- OCC_Abund_func(Wagner_bacterial_ps, Wagner_asvs)

###### 3) He et al. 2020 ######
# Load phyloseq datasets
He_bacterial_ps <- readRDS("Phyloseq_objs/He2024_phyloseq_b_asv.RDS")
He_bacterial_ps <- phyloseq::subset_samples(He_bacterial_ps, tissue != "Soil")
# Load ASV Lists
He_list1 <- read.csv("ASV_LISTS/99/He_rhizosphere_ASVs.csv")
He_list2 <- read.csv("ASV_LISTS/99/He_root_ASVs.csv")
# concatenate them together
He_asvs <- bind_rows(He_list1, He_list2) %>% distinct(ASV_id, .keep_all = TRUE)
# Add ASV prefix for study (first author first name capitalized)
He_asvs$ASV_id <- paste0("He_", He_asvs$ASV_id)
# run function for O/A table
He_df <- OCC_Abund_func(He_bacterial_ps, He_asvs)

###### 4) Naylor et al. 2017 ######
# Load phyloseq datasets
Naylor_bacterial_ps <- readRDS("Phyloseq_objs/Naylor2017_phyloseq_b_asv.RDS")
Naylor_bacterial_ps <- subset_samples(Naylor_bacterial_ps, isolation_source == "Zea_mays")
# Load ASV Lists
Naylor_list1 <- read.csv("ASV_LISTS/99/Naylor_rhizosphere_ASVs.csv")
Naylor_list2 <- read.csv("ASV_LISTS/99/Naylor_root_ASVs.csv")
# concatenate them together
Naylor_asvs <- bind_rows(Naylor_list1, Naylor_list2) %>% distinct(ASV_id, .keep_all = TRUE)
# Add ASV prefix for study (first author first name capitalized)
Naylor_asvs$ASV_id <- paste0("Naylor_", Naylor_asvs$ASV_id)
# run function for O/A table
Naylor_df <- OCC_Abund_func(Naylor_bacterial_ps, Naylor_asvs)

###### 5) GfellerCadot et al. 2023 ######
# Load phyloseq datasets
GfellerCadot_bacterial_ps <- readRDS("Phyloseq_objs/GfellerCadot2023_phyloseq_b_asv.RDS")
GfellerCadot_bacterial_ps <- subset_samples(GfellerCadot_bacterial_ps, Species=="Maize")
GfellerCadot_bacterial_ps <- prune_taxa(taxa_sums(GfellerCadot_bacterial_ps) > 0, GfellerCadot_bacterial_ps)
# Load ASV Lists
GfellerCadot_list1 <- read.csv("ASV_LISTS/99/GfellerCadot_rhizosphere_ASVs.csv")
GfellerCadot_list2 <- read.csv("ASV_LISTS/99/GfellerCadot_root_ASVs.csv")
GfellerCadot_list3 <- read.csv("ASV_LISTS/99/GfellerCadot_soil_ASVs.csv")
# concatenate them together
GfellerCadot_asvs <- bind_rows(GfellerCadot_list1, GfellerCadot_list2, GfellerCadot_list3) %>% distinct(ASV_id, .keep_all = TRUE)
# Add ASV prefix for study (first author first name capitalized)
GfellerCadot_asvs$ASV_id <- paste0("GfellerCadot_", GfellerCadot_asvs$ASV_id)
# run function for O/A table
GfellerCadot_df <- OCC_Abund_func(GfellerCadot_bacterial_ps, GfellerCadot_asvs)

###### 6) Gfeller et al. 2023 ######
# Load phyloseq datasets
Gfeller_bacterial_ps <- readRDS("Phyloseq_objs/Gfeller2023_phyloseq_b_asv.RDS")
Gfeller_bacterial_ps <- subset_samples(Gfeller_bacterial_ps, Species=="maize")
Gfeller_bacterial_ps <- prune_taxa(taxa_sums(Gfeller_bacterial_ps) > 0, Gfeller_bacterial_ps)
# Load ASV Lists
Gfeller_list1 <- read.csv("ASV_LISTS/99/Gfeller_rhizosphere_ASVs.csv")
Gfeller_list2 <- read.csv("ASV_LISTS/99/Gfeller_root_ASVs.csv")
Gfeller_list3 <- read.csv("ASV_LISTS/99/Gfeller_soil_ASVs.csv")
# concatenate them together
Gfeller_asvs <- bind_rows(Gfeller_list1, Gfeller_list2, Gfeller_list3) %>% distinct(ASV_id, .keep_all = TRUE)
# Add ASV prefix for study (first author first name capitalized)
Gfeller_asvs$ASV_id <- paste0("Gfeller_", Gfeller_asvs$ASV_id)
# run function for O/A table
Gfeller_df <- OCC_Abund_func(Gfeller_bacterial_ps, Gfeller_asvs)

###### 7) Cadot et al. 2021 ######
# Load phyloseq datasets
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
# run function for O/A table
Cadot_df <- OCC_Abund_func(Cadot_bacterial_ps, Cadot_asvs)

###### 8) Hu et al. 2018 ######
# Load phyloseq datasets
Cadot_bacterial_ps <- readRDS("Phyloseq_objs/Cadot_phyloseq_b_asv.RDS")
Hu_bacterial_ps <- subset_samples(Cadot_bacterial_ps, location  == "Changins")
Hu_bacterial_ps <- prune_taxa(taxa_sums(Hu_bacterial_ps) > 0, Hu_bacterial_ps)
taxa_names(Hu_bacterial_ps) <- sub("^Cadot_", "Hu_", taxa_names(Hu_bacterial_ps))
# Load ASV Lists
Hu_list1 <- read.csv("ASV_LISTS/99/Hu_rhizosphere_ASVs.csv")
Hu_list2 <- read.csv("ASV_LISTS/99/Hu_root_ASVs.csv")
Hu_list3 <- read.csv("ASV_LISTS/99/Hu_soil_ASVs.csv")
# concatenate them together
Hu_asvs <- bind_rows(Hu_list1, Hu_list2, Hu_list3) %>% distinct(ASV_id, .keep_all = TRUE)
# Add ASV prefix for study (first author first name capitalized)
Hu_asvs$ASV_id <- paste0("Hu_", Hu_asvs$ASV_id)
# run function for O/A table
Hu_df <- OCC_Abund_func(Hu_bacterial_ps, Hu_asvs)

###### 9) Walters et al. 2018 ######
# Load phyloseq datasets
Walters_bacterial_ps <- readRDS("Phyloseq_objs/Walters2018_phyloseq_b_asv.RDS")
# Load ASV Lists
Walters_list1 <- read.csv("ASV_LISTS/99/Walters_rhizosphere_ASVs.csv")
Walters_asvs <- Walters_list1
# Add ASV prefix for study (first author first name capitalized)
Walters_asvs$ASV_id <- paste0("Walters_", Walters_list1$ASV_id)
# run function for O/A table
Walters_df <- OCC_Abund_func(Walters_bacterial_ps, Walters_asvs)

###### 10) Swift et al. 2024 ######
# Load phyloseq datasets
Swift_bacterial_ps <- readRDS("Phyloseq_objs/Swift2024_phyloseq_b_asv.RDS")
Swift_bacterial_ps@otu_table
# Load ASV Lists
Swift_list1 <- read.csv("ASV_LISTS/99/Swift_root_ASVs.csv")
Swift_asvs <- Swift_list1
# Add ASV prefix for study (first author first name capitalized)
Swift_asvs$ASV_id <- paste0("b", Swift_list1$ASV_id)
# run function for O/A table
Swift_df <- OCC_Abund_func(Swift_bacterial_ps, Swift_asvs)

#### Collate Results 99 ####
objs <- list(
  Swift        = Swift_df,
  Naylor       = Naylor_df,
  He           = He_df,
  Wagner       = Wagner_df,
  Meier        = Meier_df,
  Walters      = Walters_df,
  GfellerCadot = GfellerCadot_df,
  Gfeller      = Gfeller_df,
  Cadot        = Cadot_df,
  Hu           = Hu_df
)


# Extract the "Yes" percentage from each
result99 <- map_dfr(names(objs), function(nm) {
  
  out <- summarize_isolate_output(objs[[nm]], summarised = TRUE)
  
  tibble(
    object = nm,
    Yes = out %>%
      filter(InCollection == "Yes") %>%
      pull(TotalAbundance)
  )
})

result99

#### 97% ####
###### 1) Meier et al. 2022 ######
# Load phyloseq datasets
Meier_bacterial_ps <- readRDS("Phyloseq_objs/Meier_bacterial_phyloseq.RDS")
# Load ASV Lists
Meier_rhizosphere_asvs <- read.csv("ASV_LISTS/97/Meier_rhizosphere_ASVs.csv")
# Add ASV prefix for study (first author first name capitalized)
Meier_rhizosphere_asvs$ASV_id <- paste0("Meier_", Meier_rhizosphere_asvs$ASV_id)
# run function for O/A table
Meier_df <- OCC_Abund_func(Meier_bacterial_ps, Meier_rhizosphere_asvs)

###### 2) Wagner et al. 2020 ######
# Load phyloseq datasets
Wagner_bacterial_ps <- readRDS("Phyloseq_objs/Wagner_bacterial_phyloseq.RDS")
# Load ASV Lists
Wagner_list1 <- read.csv("ASV_LISTS/97/Wagner_phyllosphere_ASVs.csv")
Wagner_list2 <- read.csv("ASV_LISTS/97/Wagner_rhizosphere_ASVs.csv")
# concatenate them togher
Wagner_asvs <- bind_rows(Wagner_list1, Wagner_list2) %>% distinct(ASV_id, .keep_all = TRUE)
# Add ASV prefix for study (first author first name capitalized)
Wagner_asvs$ASV_id <- paste0("Wagner_", Wagner_asvs$ASV_id)
# run function for O/A table
Wagner_df <- OCC_Abund_func(Wagner_bacterial_ps, Wagner_asvs)

###### 3) He et al. 2020 ######
# Load phyloseq datasets
He_bacterial_ps <- readRDS("Phyloseq_objs/He2024_phyloseq_b_asv.RDS")
He_bacterial_ps <- phyloseq::subset_samples(He_bacterial_ps, tissue != "Soil")
# Load ASV Lists
He_list1 <- read.csv("ASV_LISTS/97/He_rhizosphere_ASVs.csv")
He_list2 <- read.csv("ASV_LISTS/97/He_root_ASVs.csv")
# concatenate them together
He_asvs <- bind_rows(He_list1, He_list2) %>% distinct(ASV_id, .keep_all = TRUE)
# Add ASV prefix for study (first author first name capitalized)
He_asvs$ASV_id <- paste0("He_", He_asvs$ASV_id)
# run function for O/A table
He_df <- OCC_Abund_func(He_bacterial_ps, He_asvs)

###### 4) Naylor et al. 2017 ######
# Load phyloseq datasets
Naylor_bacterial_ps <- readRDS("Phyloseq_objs/Naylor2017_phyloseq_b_asv.RDS")
Naylor_bacterial_ps <- subset_samples(Naylor_bacterial_ps, isolation_source == "Zea_mays")
# Load ASV Lists
Naylor_list1 <- read.csv("ASV_LISTS/97/Naylor_rhizosphere_ASVs.csv")
Naylor_list2 <- read.csv("ASV_LISTS/97/Naylor_root_ASVs.csv")
# concatenate them together
Naylor_asvs <- bind_rows(Naylor_list1, Naylor_list2) %>% distinct(ASV_id, .keep_all = TRUE)
# Add ASV prefix for study (first author first name capitalized)
Naylor_asvs$ASV_id <- paste0("Naylor_", Naylor_asvs$ASV_id)
# run function for O/A table
Naylor_df <- OCC_Abund_func(Naylor_bacterial_ps, Naylor_asvs)

###### 5) GfellerCadot et al. 2023 ######
# Load phyloseq datasets
GfellerCadot_bacterial_ps <- readRDS("Phyloseq_objs/GfellerCadot2023_phyloseq_b_asv.RDS")
GfellerCadot_bacterial_ps <- subset_samples(GfellerCadot_bacterial_ps, Species=="Maize")
GfellerCadot_bacterial_ps <- prune_taxa(taxa_sums(GfellerCadot_bacterial_ps) > 0, GfellerCadot_bacterial_ps)
# Load ASV Lists
GfellerCadot_list1 <- read.csv("ASV_LISTS/97/GfellerCadot_rhizosphere_ASVs.csv")
GfellerCadot_list2 <- read.csv("ASV_LISTS/97/GfellerCadot_root_ASVs.csv")
GfellerCadot_list3 <- read.csv("ASV_LISTS/97/GfellerCadot_soil_ASVs.csv")
# concatenate them together
GfellerCadot_asvs <- bind_rows(GfellerCadot_list1, GfellerCadot_list2, GfellerCadot_list3) %>% distinct(ASV_id, .keep_all = TRUE)
# Add ASV prefix for study (first author first name capitalized)
GfellerCadot_asvs$ASV_id <- paste0("GfellerCadot_", GfellerCadot_asvs$ASV_id)
# run function for O/A table
GfellerCadot_df <- OCC_Abund_func(GfellerCadot_bacterial_ps, GfellerCadot_asvs)

###### 6) Gfeller et al. 2023 ######
# Load phyloseq datasets
Gfeller_bacterial_ps <- readRDS("Phyloseq_objs/Gfeller2023_phyloseq_b_asv.RDS")
Gfeller_bacterial_ps <- subset_samples(Gfeller_bacterial_ps, Species=="maize")
Gfeller_bacterial_ps <- prune_taxa(taxa_sums(Gfeller_bacterial_ps) > 0, Gfeller_bacterial_ps)
# Load ASV Lists
Gfeller_list1 <- read.csv("ASV_LISTS/97/Gfeller_rhizosphere_ASVs.csv")
Gfeller_list2 <- read.csv("ASV_LISTS/97/Gfeller_root_ASVs.csv")
Gfeller_list3 <- read.csv("ASV_LISTS/97/Gfeller_soil_ASVs.csv")
# concatenate them together
Gfeller_asvs <- bind_rows(Gfeller_list1, Gfeller_list2, Gfeller_list3) %>% distinct(ASV_id, .keep_all = TRUE)
# Add ASV prefix for study (first author first name capitalized)
Gfeller_asvs$ASV_id <- paste0("Gfeller_", Gfeller_asvs$ASV_id)
# run function for O/A table
Gfeller_df <- OCC_Abund_func(Gfeller_bacterial_ps, Gfeller_asvs)

###### 7) Cadot et al. 2021 ######
# Load phyloseq datasets
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
# run function for O/A table
Cadot_df <- OCC_Abund_func(Cadot_bacterial_ps, Cadot_asvs)

###### 8) Hu et al. 2018 ######
# Load phyloseq datasets
Cadot_bacterial_ps <- readRDS("Phyloseq_objs/Cadot_phyloseq_b_asv.RDS")
Hu_bacterial_ps <- subset_samples(Cadot_bacterial_ps, location  == "Changins")
Hu_bacterial_ps <- prune_taxa(taxa_sums(Hu_bacterial_ps) > 0, Hu_bacterial_ps)
taxa_names(Hu_bacterial_ps) <- sub("^Cadot_", "Hu_", taxa_names(Hu_bacterial_ps))
# Load ASV Lists
Hu_list1 <- read.csv("ASV_LISTS/97/Hu_rhizosphere_ASVs.csv")
Hu_list2 <- read.csv("ASV_LISTS/97/Hu_root_ASVs.csv")
Hu_list3 <- read.csv("ASV_LISTS/97/Hu_soil_ASVs.csv")
# concatenate them together
Hu_asvs <- bind_rows(Hu_list1, Hu_list2, Hu_list3) %>% distinct(ASV_id, .keep_all = TRUE)
# Add ASV prefix for study (first author first name capitalized)
Hu_asvs$ASV_id <- paste0("Hu_", Hu_asvs$ASV_id)
# run function for O/A table
Hu_df <- OCC_Abund_func(Hu_bacterial_ps, Hu_asvs)

###### 9) Walters et al. 2018 ######
# Load phyloseq datasets
Walters_bacterial_ps <- readRDS("Phyloseq_objs/Walters2018_phyloseq_b_asv.RDS")
# Load ASV Lists
Walters_list1 <- read.csv("ASV_LISTS/97/Walters_rhizosphere_ASVs.csv")
Walters_asvs <- Walters_list1
# Add ASV prefix for study (first author first name capitalized)
Walters_asvs$ASV_id <- paste0("Walters_", Walters_list1$ASV_id)
# run function for O/A table
Walters_df <- OCC_Abund_func(Walters_bacterial_ps, Walters_asvs)

###### 10) Swift et al. 2024 ######
# Load phyloseq datasets
Swift_bacterial_ps <- readRDS("Phyloseq_objs/Swift2024_phyloseq_b_asv.RDS")
Swift_bacterial_ps@otu_table
# Load ASV Lists
Swift_list1 <- read.csv("ASV_LISTS/97/Swift_root_ASVs.csv")
Swift_asvs <- Swift_list1
# Add ASV prefix for study (first author first name capitalized)
Swift_asvs$ASV_id <- paste0("b", Swift_list1$ASV_id)
# run function for O/A table
Swift_df <- OCC_Abund_func(Swift_bacterial_ps, Swift_asvs)

#### Collate Results 97 ####
objs <- list(
  Swift        = Swift_df,
  Naylor       = Naylor_df,
  He           = He_df,
  Wagner       = Wagner_df,
  Meier        = Meier_df,
  Walters      = Walters_df,
  GfellerCadot = GfellerCadot_df,
  Gfeller      = Gfeller_df,
  Cadot        = Cadot_df,
  Hu           = Hu_df
)


# Extract the "Yes" percentage from each
result97 <- map_dfr(names(objs), function(nm) {
  
  out <- summarize_isolate_output(objs[[nm]], summarised = TRUE)
  
  tibble(
    object = nm,
    Yes = out %>%
      filter(InCollection == "Yes") %>%
      pull(TotalAbundance)
  )
})

result99
result97


# ==================================================
# Metadata table
# ==================================================
datasets <- tribble(
  ~study,            ~compartment, ~level, ~file,                                             ~ps_object,
  
  "Wagner",          "phy",        "99",   "ASV_LISTS/99/Wagner_phyllosphere_ASVs.csv",       Wagner_bacterial_ps,
  "Wagner",          "phy",        "97",   "ASV_LISTS/97/Wagner_phyllosphere_ASVs.csv",       Wagner_bacterial_ps,
  "Wagner",          "rhz",        "99",   "ASV_LISTS/99/Wagner_rhizosphere_ASVs.csv",        Wagner_bacterial_ps,
  "Wagner",          "rhz",        "97",   "ASV_LISTS/97/Wagner_rhizosphere_ASVs.csv",        Wagner_bacterial_ps,
  
  "He",              "rhz",        "99",   "ASV_LISTS/99/He_rhizosphere_ASVs.csv",            He_bacterial_ps,
  "He",              "rhz",        "97",   "ASV_LISTS/97/He_rhizosphere_ASVs.csv",            He_bacterial_ps,
  "He",              "rot",        "99",   "ASV_LISTS/99/He_root_ASVs.csv",                   He_bacterial_ps,
  "He",              "rot",        "97",   "ASV_LISTS/97/He_root_ASVs.csv",                   He_bacterial_ps,
  
  "Naylor",          "rhz",        "99",   "ASV_LISTS/99/Naylor_rhizosphere_ASVs.csv",        Naylor_bacterial_ps,
  "Naylor",          "rhz",        "97",   "ASV_LISTS/97/Naylor_rhizosphere_ASVs.csv",        Naylor_bacterial_ps,
  "Naylor",          "rot",        "99",   "ASV_LISTS/99/Naylor_root_ASVs.csv",               Naylor_bacterial_ps,
  "Naylor",          "rot",        "97",   "ASV_LISTS/97/Naylor_root_ASVs.csv",               Naylor_bacterial_ps,
  
  "GfellerCadot",    "soi",        "99",   "ASV_LISTS/99/GfellerCadot_soil_ASVs.csv",         GfellerCadot_bacterial_ps,
  "GfellerCadot",    "soi",        "97",   "ASV_LISTS/97/GfellerCadot_soil_ASVs.csv",         GfellerCadot_bacterial_ps,
  "GfellerCadot",    "rhz",        "99",   "ASV_LISTS/99/GfellerCadot_rhizosphere_ASVs.csv",  GfellerCadot_bacterial_ps,
  "GfellerCadot",    "rhz",        "97",   "ASV_LISTS/97/GfellerCadot_rhizosphere_ASVs.csv",  GfellerCadot_bacterial_ps,
  "GfellerCadot",    "rot",        "99",   "ASV_LISTS/99/GfellerCadot_root_ASVs.csv",         GfellerCadot_bacterial_ps,
  "GfellerCadot",    "rot",        "97",   "ASV_LISTS/97/GfellerCadot_root_ASVs.csv",         GfellerCadot_bacterial_ps,
  
  "Gfeller",         "soi",        "99",   "ASV_LISTS/99/Gfeller_soil_ASVs.csv",              Gfeller_bacterial_ps,
  "Gfeller",         "soi",        "97",   "ASV_LISTS/97/Gfeller_soil_ASVs.csv",              Gfeller_bacterial_ps,
  "Gfeller",         "rhz",        "99",   "ASV_LISTS/99/Gfeller_rhizosphere_ASVs.csv",       Gfeller_bacterial_ps,
  "Gfeller",         "rhz",        "97",   "ASV_LISTS/97/Gfeller_rhizosphere_ASVs.csv",       Gfeller_bacterial_ps,
  "Gfeller",         "rot",        "99",   "ASV_LISTS/99/Gfeller_root_ASVs.csv",              Gfeller_bacterial_ps,
  "Gfeller",         "rot",        "97",   "ASV_LISTS/97/Gfeller_root_ASVs.csv",              Gfeller_bacterial_ps,
  
  "Cadot",           "soi",        "99",   "ASV_LISTS/99/Cadot_soil_ASVs.csv",                Cadot_bacterial_ps,
  "Cadot",           "soi",        "97",   "ASV_LISTS/97/Cadot_soil_ASVs.csv",                Cadot_bacterial_ps,
  "Cadot",           "rhz",        "99",   "ASV_LISTS/99/Cadot_rhizosphere_ASVs.csv",         Cadot_bacterial_ps,
  "Cadot",           "rhz",        "97",   "ASV_LISTS/97/Cadot_rhizosphere_ASVs.csv",         Cadot_bacterial_ps,
  "Cadot",           "rot",        "99",   "ASV_LISTS/99/Cadot_root_ASVs.csv",                Cadot_bacterial_ps,
  "Cadot",           "rot",        "97",   "ASV_LISTS/97/Cadot_root_ASVs.csv",                Cadot_bacterial_ps,
  
  "Hu",              "soi",        "99",   "ASV_LISTS/99/Hu_soil_ASVs.csv",                   Hu_bacterial_ps,
  "Hu",              "soi",        "97",   "ASV_LISTS/97/Hu_soil_ASVs.csv",                   Hu_bacterial_ps,
  "Hu",              "rhz",        "99",   "ASV_LISTS/99/Hu_rhizosphere_ASVs.csv",            Hu_bacterial_ps,
  "Hu",              "rhz",        "97",   "ASV_LISTS/97/Hu_rhizosphere_ASVs.csv",            Hu_bacterial_ps,
  "Hu",              "rot",        "99",   "ASV_LISTS/99/Hu_root_ASVs.csv",                   Hu_bacterial_ps,
  "Hu",              "rot",        "97",   "ASV_LISTS/97/Hu_root_ASVs.csv",                   Hu_bacterial_ps,
  
  "Walters",         "rhz",        "99",   "ASV_LISTS/99/Walters_rhizosphere_ASVs.csv",       Walters_bacterial_ps,
  "Walters",         "rhz",        "97",   "ASV_LISTS/97/Walters_rhizosphere_ASVs.csv",       Walters_bacterial_ps,
  "Walters",         "soi",        "99",   "ASV_LISTS/99/Walters_soil_ASVs.csv",              Walters_bacterial_ps,
  "Walters",         "soi",        "97",   "ASV_LISTS/97/Walters_soil_ASVs.csv",              Walters_bacterial_ps
)

# ==================================================
# Run all analyses and summarize
# ==================================================
results_stratified <- pmap_dfr(
  datasets,
  function(study, compartment, level, file, ps_object) {
    
    # read ASV table
    df <- read.csv(file)
    
    # add study prefix
    df$ASV_id <- paste0(study, "_", df$ASV_id)
    
    # run neutral model
    model_out <- OCC_Abund_func(ps_object, df)
    
    # extract Yes abundance
    yes_abundance <- summarize_isolate_output(
      model_out,
      summarised = TRUE
    ) %>%
      filter(InCollection == "Yes") %>%
      pull(TotalAbundance)
    
    tibble(
      Study = study,
      Compartment = compartment,
      Identity = level,
      Yes = yes_abundance
    )
  }
)

# optional prettier labels
results_stratified <- results_stratified %>%
  mutate(
    Compartment = recode(
      Compartment,
      phy = "Phyllosphere",
      rot = "Root endosphere",
      rhz = "Rhizosphere",
      soi = "Soil"
    )
  )

results_stratified

results_all <- result99 %>%
  rename(Yes_99 = Yes) %>%
  left_join(result97 %>% rename(Yes_97 = Yes), by = "object")

write.csv(x = results_stratified, file = "stats/RA_stats_stratified.csv")
write.csv(x = results_all, file = "stats/RA_stats.csv")