# Abundance Stats from ASV lists
# Map ASV lists back to studies and pull the stats

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

# Output stats
summarize_isolate_output(Meier_df, summarised = FALSE)
summarize_isolate_output(Meier_df, summarised = TRUE)

# InCollection `sum(Abundpercent)`
# 1 No                          72.5
# 2 Yes                         23.0



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

# Occupancy-Abundance with normal model curve
wagner_plot <- ggplot() +
  geom_point(data = subset(Wagner_df$occ_abun_tax, InCollection == "No"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 0.5, fill = "white") +
  geom_point(data = subset(Wagner_df$occ_abun_tax, InCollection == "Yes"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 1, fill = "#50C878") +
  geom_line(color='black', data=Wagner_df$obs_np, size=1, aes(y=freq.pred, x=log10(p)), alpha=.25) +
  geom_line(color='black', lty='twodash', size=1, data=Wagner_df$obs_np, aes(y=pred.upr, x=log10(p)), alpha=.25)+
  geom_line(color='black', lty='twodash', size=1, data=Wagner_df$obs_np, aes(y=pred.lwr, x=log10(p)), alpha=.25)+
  labs(x="log10(mean relative abundance)", y="Occupancy", title = "Wagner et al. (2020) - North Carlonia Field Grown Maize Rhizosphere/Phyllosphere")

ggsave("OA_curves/wagner_OA_curve.svg", wagner_plot, height = 19, width = 19, unit = "cm")

# Output stats
summarize_isolate_output(Wagner_df, summarised = FALSE)
summarize_isolate_output(Wagner_df, summarised = TRUE)

# InCollection TotalAbundance
# 1 No                     58.1
# 2 Yes                    41.7



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

# Occupancy-Abundance with normal model curve
he_plot <- ggplot() +
  geom_point(data = subset(He_df$occ_abun_tax, InCollection == "No"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 0.5, fill = "white") +
  geom_point(data = subset(He_df$occ_abun_tax, InCollection == "Yes"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 1, fill = "#50C878") +
  geom_line(color='black', data=He_df$obs_np, size=1, aes(y=freq.pred, x=log10(p)), alpha=.25) +
  geom_line(color='black', lty='twodash', size=1, data=He_df$obs_np, aes(y=pred.upr, x=log10(p)), alpha=.25)+
  geom_line(color='black', lty='twodash', size=1, data=He_df$obs_np, aes(y=pred.lwr, x=log10(p)), alpha=.25)+
  labs(x="log10(mean relative abundance)", y="Occupancy", title = "He et al. (2024) - German Soil, pot grown Maize Rhizosphere/Root Endosphere")

ggsave("OA_curves/He_OA_curve.svg", he_plot, height = 19, width = 19, unit = "cm")

# Output stats
summarize_isolate_output(He_df, summarised = FALSE)
summarize_isolate_output(He_df, summarised = TRUE)

# InCollection TotalAbundance
# 1 No                     82.6
# 2 Yes                    14.9



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

# Occupancy-Abundance with normal model curve
naylor_plot <- ggplot() +
  geom_point(data = subset(Naylor_df$occ_abun_tax, InCollection == "No"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 0.5, fill = "white") +
  geom_point(data = subset(Naylor_df$occ_abun_tax, InCollection == "Yes"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 1, fill = "#50C878") +
  geom_line(color='black', data=Naylor_df$obs_np, size=1, aes(y=freq.pred, x=log10(p)), alpha=.25) +
  geom_line(color='black', lty='twodash', size=1, data=Naylor_df$obs_np, aes(y=pred.upr, x=log10(p)), alpha=.25)+
  geom_line(color='black', lty='twodash', size=1, data=Naylor_df$obs_np, aes(y=pred.lwr, x=log10(p)), alpha=.25)+
  labs(x="log10(mean relative abundance)", y="Occupancy", title = "Naylor et al. (2017) - California Field Grown Maize Rhizosphere/Root Endosphere")

ggsave("OA_curves/Naylor_OA_curve.svg", naylor_plot, height = 19, width = 19, unit = "cm")

# Output stats
summarize_isolate_output(Naylor_df, summarised = FALSE)
summarize_isolate_output(Naylor_df, summarised = TRUE)

# InCollection TotalAbundance
# 1 No                    91.4 
# 2 Yes                    8.76



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

# Occupancy-Abundance with normal model curve
gfellercadot_plot <- ggplot() +
  geom_point(data = subset(GfellerCadot_df$occ_abun_tax, InCollection == "No"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 0.5, fill = "white") +
  geom_point(data = subset(GfellerCadot_df$occ_abun_tax, InCollection == "Yes"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 1, fill = "#50C878") +
  geom_line(color='black', data=GfellerCadot_df$obs_np, size=1, aes(y=freq.pred, x=log10(p)), alpha=.25) +
  geom_line(color='black', lty='twodash', size=1, data=GfellerCadot_df$obs_np, aes(y=pred.upr, x=log10(p)), alpha=.25)+
  geom_line(color='black', lty='twodash', size=1, data=GfellerCadot_df$obs_np, aes(y=pred.lwr, x=log10(p)), alpha=.25)+
  labs(x="log10(mean relative abundance)", y="Occupancy", title = "Gfeller Cadot et al. (2023) - Switzerland Field Maize Soil/Rhizosphere/Root Endosphere")

ggsave("OA_curves/gfellercadot_OA_curve.svg", gfellercadot_plot, height = 19, width = 19, unit = "cm")


# Output stats
summarize_isolate_output(GfellerCadot_df, summarised = FALSE)
summarize_isolate_output(GfellerCadot_df, summarised = TRUE)

# InCollection TotalAbundance
# <chr>                 <dbl>
# 1 No                     81.2
# 2 Yes                    18.6



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

# Occupancy-Abundance with normal model curve
gfeller_plot <- ggplot() +
  geom_point(data = subset(Gfeller_df$occ_abun_tax, InCollection == "No"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 0.5, fill = "white") +
  geom_point(data = subset(Gfeller_df$occ_abun_tax, InCollection == "Yes"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 1, fill = "#50C878") +
  geom_line(color='black', data=Gfeller_df$obs_np, size=1, aes(y=freq.pred, x=log10(p)), alpha=.25) +
  geom_line(color='black', lty='twodash', size=1, data=Gfeller_df$obs_np, aes(y=pred.upr, x=log10(p)), alpha=.25)+
  geom_line(color='black', lty='twodash', size=1, data=Gfeller_df$obs_np, aes(y=pred.lwr, x=log10(p)), alpha=.25)+
  labs(x="log10(mean relative abundance)", y="Occupancy", title = "Gfeller et al. (2023) - Switzerland Field Maize Soil/Rhizosphere/Root Endosphere")

ggsave("OA_curves/gfeller_OA_curve.svg", gfeller_plot, height = 19, width = 19, unit = "cm")

# Output stats
summarize_isolate_output(Gfeller_df, summarised = FALSE)
summarize_isolate_output(Gfeller_df, summarised = TRUE)

# InCollection TotalAbundance
# 1 No                     86.1
# 2 Yes                    13.8



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

# Occupancy-Abundance with normal model curve
cadot_plot <- ggplot() +
  geom_point(data = subset(Cadot_df$occ_abun_tax, InCollection == "No"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 0.5, fill = "white") +
  geom_point(data = subset(Cadot_df$occ_abun_tax, InCollection == "Yes"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 1, fill = "#50C878") +
  geom_line(color='black', data=Cadot_df$obs_np, size=1, aes(y=freq.pred, x=log10(p)), alpha=.25) +
  geom_line(color='black', lty='twodash', size=1, data=Cadot_df$obs_np, aes(y=pred.upr, x=log10(p)), alpha=.25)+
  geom_line(color='black', lty='twodash', size=1, data=Cadot_df$obs_np, aes(y=pred.lwr, x=log10(p)), alpha=.25)+
  labs(x="log10(mean relative abundance)", y="Occupancy", title = "Cadot et al. (2021) - Switzerland/New York Field Maize Soil/Rhizosphere/Root Endosphere")

ggsave("OA_curves/cadot_OA_curve.svg", cadot_plot, height = 19, width = 19, unit = "cm")

# Output stats
summarize_isolate_output(Cadot_df, summarised = FALSE)
summarize_isolate_output(Cadot_df, summarised = TRUE)

# InCollection TotalAbundance
# 1 No                     72.2
# 2 Yes                    27.8



###### 8) Hu et al. 2018 ######
# Load phyloseq datasets
Cadot_bacterial_ps <- readRDS("Phyloseq_objs/Cadot_phyloseq_b_asv.RDS")
Hu_bacterial_ps <- subset_samples(Cadot_bacterial_ps, location  == "Changins")
Hu_bacterial_ps <- prune_taxa(taxa_sums(Hu_bacterial_ps) > 0, Hu_bacterial_ps)
# Load ASV Lists
Hu_list1 <- read.csv("ASV_LISTS/99/Hu_rhizosphere_ASVs.csv")
Hu_list2 <- read.csv("ASV_LISTS/99/Hu_root_ASVs.csv")
Hu_list3 <- read.csv("ASV_LISTS/99/Hu_soil_ASVs.csv")
# concatenate them together
Hu_asvs <- bind_rows(Hu_list1, Hu_list2, Hu_list3) %>% distinct(ASV_id, .keep_all = TRUE)
# Add ASV prefix for study (first author first name capitalized)
Hu_asvs$ASV_id <- paste0("Cadot_", Hu_asvs$ASV_id)

# run function for O/A table
Hu_df <- OCC_Abund_func(Hu_bacterial_ps, Hu_asvs)

# Occupancy-Abundance with normal model curve
Hu_plot <- ggplot() +
  geom_point(data = subset(Hu_df$occ_abun_tax, InCollection == "No"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 0.5, fill = "white") +
  geom_point(data = subset(Hu_df$occ_abun_tax, InCollection == "Yes"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 1, fill = "#50C878") +
  geom_line(color='black', data=Hu_df$obs_np, size=1, aes(y=freq.pred, x=log10(p)), alpha=.25) +
  geom_line(color='black', lty='twodash', size=1, data=Hu_df$obs_np, aes(y=pred.upr, x=log10(p)), alpha=.25)+
  geom_line(color='black', lty='twodash', size=1, data=Hu_df$obs_np, aes(y=pred.lwr, x=log10(p)), alpha=.25)+
  labs(x="log10(mean relative abundance)", y="Occupancy", title = "Hu et al. (2023) - Switzerland Maize Rhizosphere/Root Endosphere")

ggsave("OA_curves/Hu_OA_curve.svg", Hu_plot, height = 19, width = 19, unit = "cm")

# Output stats
summarize_isolate_output(Hu_df, summarised = FALSE)
summarize_isolate_output(Hu_df, summarised = TRUE)

# InCollection TotalAbundance
# 1 No                     85.4
# 2 Yes                    13.6



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

# Occupancy-Abundance with normal model curve
walters_plot <- ggplot() +
  geom_point(data = subset(Walters_df$occ_abun_tax, InCollection == "No"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 0.5, fill = "white") +
  geom_point(data = subset(Walters_df$occ_abun_tax, InCollection == "Yes"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 1, fill = "#50C878") +
  geom_line(color='black', data=Walters_df$obs_np, size=1, aes(y=freq.pred, x=log10(p)), alpha=.25) +
  geom_line(color='black', lty='twodash', size=1, data=Walters_df$obs_np, aes(y=pred.upr, x=log10(p)), alpha=.25)+
  geom_line(color='black', lty='twodash', size=1, data=Walters_df$obs_np, aes(y=pred.lwr, x=log10(p)), alpha=.25)+
  labs(x="log10(mean relative abundance)", y="Occupancy", title = "Walters et al. (2018) - NY,MO,IL Maize Rhizosphere")

ggsave("OA_curves/walters_OA_curve.svg", walters_plot, height = 19, width = 19, unit = "cm")

# Output stats
summarize_isolate_output(Walters_df, summarised = FALSE)
summarize_isolate_output(Walters_df, summarised = TRUE)

# InCollection TotalAbundance
# 1 No                     74.7
# 2 Yes                    18.4

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

# Occupancy-Abundance with normal model curve
Swift_plot <- ggplot() +
  geom_point(data = subset(Swift_df$occ_abun_tax, InCollection == "No"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 0.5, fill = "white") +
  geom_point(data = subset(Swift_df$occ_abun_tax, InCollection == "Yes"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 1, fill = "#50C878") +
  geom_line(color='black', data=Swift_df$obs_np, size=1, aes(y=freq.pred, x=log10(p)), alpha=.25) +
  geom_line(color='black', lty='twodash', size=1, data=Swift_df$obs_np, aes(y=pred.upr, x=log10(p)), alpha=.25)+
  geom_line(color='black', lty='twodash', size=1, data=Swift_df$obs_np, aes(y=pred.lwr, x=log10(p)), alpha=.25)+
  labs(x="log10(mean relative abundance)", y="Occupancy", title = "Swift et al. (2024) - KS Maize Rhizosphere")

ggsave("OA_curves/Swift_OA_curve.svg", Swift_plot, height = 19, width = 19, unit = "cm")

# Output stats
summarize_isolate_output(Swift_df, summarised = FALSE)
summarize_isolate_output(Swift_df, summarised = TRUE)

# InCollection TotalAbundance
# 1 No                     34.2
# 2 Yes                    65.6




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

# Occupancy-Abundance with normal model curve
meier_plot <- ggplot() +
  geom_point(data = subset(Meier_df$occ_abun_tax, InCollection == "No"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 0.5, fill = "white") +
  geom_point(data = subset(Meier_df$occ_abun_tax, InCollection == "Yes"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 1, fill = "#50C878") +
  geom_line(color='black', data=Meier_df$obs_np, size=1, aes(y=freq.pred, x=log10(p)), alpha=.25) +
  geom_line(color='black', lty='twodash', size=1, data=Meier_df$obs_np, aes(y=pred.upr, x=log10(p)), alpha=.25)+
  geom_line(color='black', lty='twodash', size=1, data=Meier_df$obs_np, aes(y=pred.lwr, x=log10(p)), alpha=.25)+
  labs(x="log10(mean relative abundance)", y="Occupancy", title = "Meier et al. (2022) - Nebraska Field Grown Maize Rhizosphere")

ggsave("OA_curves/97/meier_OA_curve.svg", meier_plot, height = 19, width = 18, unit = "cm")

# Output stats
summarize_isolate_output(Meier_df, summarised = FALSE)
summarize_isolate_output(Meier_df, summarised = TRUE)

# InCollection `sum(Abundpercent)`
# 1 No                     60.9
# 2 Yes                    34.6



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

# Occupancy-Abundance with normal model curve
wagner_plot <- ggplot() +
  geom_point(data = subset(Wagner_df$occ_abun_tax, InCollection == "No"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 0.5, fill = "white") +
  geom_point(data = subset(Wagner_df$occ_abun_tax, InCollection == "Yes"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 1, fill = "#50C878") +
  geom_line(color='black', data=Wagner_df$obs_np, size=1, aes(y=freq.pred, x=log10(p)), alpha=.25) +
  geom_line(color='black', lty='twodash', size=1, data=Wagner_df$obs_np, aes(y=pred.upr, x=log10(p)), alpha=.25)+
  geom_line(color='black', lty='twodash', size=1, data=Wagner_df$obs_np, aes(y=pred.lwr, x=log10(p)), alpha=.25)+
  labs(x="log10(mean relative abundance)", y="Occupancy", title = "Wagner et al. (2020) - North Carlonia Field Grown Maize Rhizosphere/Phyllosphere")

ggsave("OA_curves/97/wagner_OA_curve.svg", wagner_plot, height = 19, width = 19, unit = "cm")

# Output stats
summarize_isolate_output(Wagner_df, summarised = FALSE)
summarize_isolate_output(Wagner_df, summarised = TRUE)

# InCollection TotalAbundance
# 1 No                     49.6
# 2 Yes                    50.2



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

# Occupancy-Abundance with normal model curve
he_plot <- ggplot() +
  geom_point(data = subset(He_df$occ_abun_tax, InCollection == "No"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 0.5, fill = "white") +
  geom_point(data = subset(He_df$occ_abun_tax, InCollection == "Yes"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 1, fill = "#50C878") +
  geom_line(color='black', data=He_df$obs_np, size=1, aes(y=freq.pred, x=log10(p)), alpha=.25) +
  geom_line(color='black', lty='twodash', size=1, data=He_df$obs_np, aes(y=pred.upr, x=log10(p)), alpha=.25)+
  geom_line(color='black', lty='twodash', size=1, data=He_df$obs_np, aes(y=pred.lwr, x=log10(p)), alpha=.25)+
  labs(x="log10(mean relative abundance)", y="Occupancy", title = "He et al. (2024) - German Soil, pot grown Maize Rhizosphere/Root Endosphere")

ggsave("OA_curves/97/He_OA_curve.svg", he_plot, height = 19, width = 19, unit = "cm")

# Output stats
summarize_isolate_output(He_df, summarised = FALSE)
summarize_isolate_output(He_df, summarised = TRUE)

# InCollection TotalAbundance
# 1 No                     60.4
# 2 Yes                    37.1



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

# Occupancy-Abundance with normal model curve
naylor_plot <- ggplot() +
  geom_point(data = subset(Naylor_df$occ_abun_tax, InCollection == "No"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 0.5, fill = "white") +
  geom_point(data = subset(Naylor_df$occ_abun_tax, InCollection == "Yes"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 1, fill = "#50C878") +
  geom_line(color='black', data=Naylor_df$obs_np, size=1, aes(y=freq.pred, x=log10(p)), alpha=.25) +
  geom_line(color='black', lty='twodash', size=1, data=Naylor_df$obs_np, aes(y=pred.upr, x=log10(p)), alpha=.25)+
  geom_line(color='black', lty='twodash', size=1, data=Naylor_df$obs_np, aes(y=pred.lwr, x=log10(p)), alpha=.25)+
  labs(x="log10(mean relative abundance)", y="Occupancy", title = "Naylor et al. (2017) - California Field Grown Maize Rhizosphere/Root Endosphere")

ggsave("OA_curves/97/Naylor_OA_curve.svg", naylor_plot, height = 19, width = 19, unit = "cm")

# Output stats
summarize_isolate_output(Naylor_df, summarised = FALSE)
summarize_isolate_output(Naylor_df, summarised = TRUE)

# InCollection TotalAbundance
# 1 No                    87.4 
# 2 Yes                   12.1



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

# Occupancy-Abundance with normal model curve
gfellercadot_plot <- ggplot() +
  geom_point(data = subset(GfellerCadot_df$occ_abun_tax, InCollection == "No"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 0.5, fill = "white") +
  geom_point(data = subset(GfellerCadot_df$occ_abun_tax, InCollection == "Yes"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 1, fill = "#50C878") +
  geom_line(color='black', data=GfellerCadot_df$obs_np, size=1, aes(y=freq.pred, x=log10(p)), alpha=.25) +
  geom_line(color='black', lty='twodash', size=1, data=GfellerCadot_df$obs_np, aes(y=pred.upr, x=log10(p)), alpha=.25)+
  geom_line(color='black', lty='twodash', size=1, data=GfellerCadot_df$obs_np, aes(y=pred.lwr, x=log10(p)), alpha=.25)+
  labs(x="log10(mean relative abundance)", y="Occupancy", title = "Gfeller Cadot et al. (2023) - Switzerland Field Maize Soil/Rhizosphere/Root Endosphere")

ggsave("OA_curves/97/gfellercadot_OA_curve.svg", gfellercadot_plot, height = 19, width = 19, unit = "cm")


# Output stats
summarize_isolate_output(GfellerCadot_df, summarised = FALSE)
summarize_isolate_output(GfellerCadot_df, summarised = TRUE)

# InCollection TotalAbundance
# <chr>                 <dbl>
# 1 No                     74.9
# 2 Yes                    25.2



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

# Occupancy-Abundance with normal model curve
gfeller_plot <- ggplot() +
  geom_point(data = subset(Gfeller_df$occ_abun_tax, InCollection == "No"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 0.5, fill = "white") +
  geom_point(data = subset(Gfeller_df$occ_abun_tax, InCollection == "Yes"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 1, fill = "#50C878") +
  geom_line(color='black', data=Gfeller_df$obs_np, size=1, aes(y=freq.pred, x=log10(p)), alpha=.25) +
  geom_line(color='black', lty='twodash', size=1, data=Gfeller_df$obs_np, aes(y=pred.upr, x=log10(p)), alpha=.25)+
  geom_line(color='black', lty='twodash', size=1, data=Gfeller_df$obs_np, aes(y=pred.lwr, x=log10(p)), alpha=.25)+
  labs(x="log10(mean relative abundance)", y="Occupancy", title = "Gfeller et al. (2023) - Switzerland Field Maize Soil/Rhizosphere/Root Endosphere")

ggsave("OA_curves/97/gfeller_OA_curve.svg", gfeller_plot, height = 19, width = 19, unit = "cm")

# Output stats
summarize_isolate_output(Gfeller_df, summarised = FALSE)
summarize_isolate_output(Gfeller_df, summarised = TRUE)

# InCollection TotalAbundance
# 1 No                     65.9
# 2 Yes                    34.1



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

# Occupancy-Abundance with normal model curve
cadot_plot <- ggplot() +
  geom_point(data = subset(Cadot_df$occ_abun_tax, InCollection == "No"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 0.5, fill = "white") +
  geom_point(data = subset(Cadot_df$occ_abun_tax, InCollection == "Yes"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 1, fill = "#50C878") +
  geom_line(color='black', data=Cadot_df$obs_np, size=1, aes(y=freq.pred, x=log10(p)), alpha=.25) +
  geom_line(color='black', lty='twodash', size=1, data=Cadot_df$obs_np, aes(y=pred.upr, x=log10(p)), alpha=.25)+
  geom_line(color='black', lty='twodash', size=1, data=Cadot_df$obs_np, aes(y=pred.lwr, x=log10(p)), alpha=.25)+
  labs(x="log10(mean relative abundance)", y="Occupancy", title = "Cadot et al. (2021) - Switzerland/New York Field Maize Soil/Rhizosphere/Root Endosphere")

ggsave("OA_curves/97/cadot_OA_curve.svg", cadot_plot, height = 19, width = 19, unit = "cm")

# Output stats
summarize_isolate_output(Cadot_df, summarised = FALSE)
summarize_isolate_output(Cadot_df, summarised = TRUE)

# InCollection TotalAbundance
# 1 No                     56.1
# 2 Yes                    43.9



###### 8) Hu et al. 2018 ######
# Load phyloseq datasets
Cadot_bacterial_ps <- readRDS("Phyloseq_objs/Cadot_phyloseq_b_asv.RDS")
Hu_bacterial_ps <- subset_samples(Cadot_bacterial_ps, location  == "Changins")
Hu_bacterial_ps <- prune_taxa(taxa_sums(Hu_bacterial_ps) > 0, Hu_bacterial_ps)
# Load ASV Lists
Hu_list1 <- read.csv("ASV_LISTS/97/Hu_rhizosphere_ASVs.csv")
Hu_list2 <- read.csv("ASV_LISTS/97/Hu_root_ASVs.csv")
Hu_list3 <- read.csv("ASV_LISTS/97/Hu_soil_ASVs.csv")
# concatenate them together
Hu_asvs <- bind_rows(Hu_list1, Hu_list2, Hu_list3) %>% distinct(ASV_id, .keep_all = TRUE)
# Add ASV prefix for study (first author first name capitalized)
Hu_asvs$ASV_id <- paste0("Cadot_", Hu_asvs$ASV_id)

# run function for O/A table
Hu_df <- OCC_Abund_func(Hu_bacterial_ps, Hu_asvs)

# Occupancy-Abundance with normal model curve
Hu_plot <- ggplot() +
  geom_point(data = subset(Hu_df$occ_abun_tax, InCollection == "No"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 0.5, fill = "white") +
  geom_point(data = subset(Hu_df$occ_abun_tax, InCollection == "Yes"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 1, fill = "#50C878") +
  geom_line(color='black', data=Hu_df$obs_np, size=1, aes(y=freq.pred, x=log10(p)), alpha=.25) +
  geom_line(color='black', lty='twodash', size=1, data=Hu_df$obs_np, aes(y=pred.upr, x=log10(p)), alpha=.25)+
  geom_line(color='black', lty='twodash', size=1, data=Hu_df$obs_np, aes(y=pred.lwr, x=log10(p)), alpha=.25)+
  labs(x="log10(mean relative abundance)", y="Occupancy", title = "Hu et al. (2023) - Switzerland Maize Rhizosphere/Root Endosphere")

ggsave("OA_curves/97/Hu_OA_curve.svg", Hu_plot, height = 19, width = 19, unit = "cm")

# Output stats
summarize_isolate_output(Hu_df, summarised = FALSE)
summarize_isolate_output(Hu_df, summarised = TRUE)

# InCollection TotalAbundance
# 1 No                     61.8
# 2 Yes                    37.5



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

# Occupancy-Abundance with normal model curve
walters_plot <- ggplot() +
  geom_point(data = subset(Walters_df$occ_abun_tax, InCollection == "No"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 0.5, fill = "white") +
  geom_point(data = subset(Walters_df$occ_abun_tax, InCollection == "Yes"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 1, fill = "#50C878") +
  geom_line(color='black', data=Walters_df$obs_np, size=1, aes(y=freq.pred, x=log10(p)), alpha=.25) +
  geom_line(color='black', lty='twodash', size=1, data=Walters_df$obs_np, aes(y=pred.upr, x=log10(p)), alpha=.25)+
  geom_line(color='black', lty='twodash', size=1, data=Walters_df$obs_np, aes(y=pred.lwr, x=log10(p)), alpha=.25)+
  labs(x="log10(mean relative abundance)", y="Occupancy", title = "Walters et al. (2018) - NY,MO,IL Maize Rhizosphere")

ggsave("OA_curves/97/walters_OA_curve.svg", walters_plot, height = 19, width = 19, unit = "cm")

# Output stats
summarize_isolate_output(Walters_df, summarised = FALSE)
summarize_isolate_output(Walters_df, summarised = TRUE)

# InCollection TotalAbundance
# 1 No                     46.1
# 2 Yes                    47.0


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

# Occupancy-Abundance with normal model curve
Swift_plot <- ggplot() +
  geom_point(data = subset(Swift_df$occ_abun_tax, InCollection == "No"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 0.5, fill = "white") +
  geom_point(data = subset(Swift_df$occ_abun_tax, InCollection == "Yes"), aes(x=log10(otu_rel), y=otu_occ), shape = 21, size = 3, alpha = 1, fill = "#50C878") +
  geom_line(color='black', data=Swift_df$obs_np, size=1, aes(y=freq.pred, x=log10(p)), alpha=.25) +
  geom_line(color='black', lty='twodash', size=1, data=Swift_df$obs_np, aes(y=pred.upr, x=log10(p)), alpha=.25)+
  geom_line(color='black', lty='twodash', size=1, data=Swift_df$obs_np, aes(y=pred.lwr, x=log10(p)), alpha=.25)+
  labs(x="log10(mean relative abundance)", y="Occupancy", title = "Swift et al. (2024) - KS Maize Rhizosphere")

ggsave("OA_curves/97/Swift_OA_curve.svg", Swift_plot, height = 19, width = 19, unit = "cm")

# Output stats
summarize_isolate_output(Swift_df, summarised = FALSE)
summarize_isolate_output(Swift_df, summarised = TRUE)

# InCollection TotalAbundance
# 1 No                     10.9
# 2 Yes                    88.9





#### Stratified by Compartment ####
###### 2) Wagner et al. 2020 ######
Wagner_phy99 <- read.csv("ASV_LISTS/99/Wagner_phyllosphere_ASVs.csv")
Wagner_phy97 <- read.csv("ASV_LISTS/97/Wagner_phyllosphere_ASVs.csv")
Wagner_rhz99 <- read.csv("ASV_LISTS/99/Wagner_rhizosphere_ASVs.csv")
Wagner_rhz97 <- read.csv("ASV_LISTS/97/Wagner_rhizosphere_ASVs.csv")

Wagner_phy99$ASV_id <- paste0("Wagner_", Wagner_phy99$ASV_id)
Wagner_phy97$ASV_id <- paste0("Wagner_", Wagner_phy97$ASV_id)
Wagner_rhz99$ASV_id <- paste0("Wagner_", Wagner_rhz99$ASV_id)
Wagner_rhz97$ASV_id <- paste0("Wagner_", Wagner_rhz97$ASV_id)

# run function for O/A table
Wagner_phy99 <- OCC_Abund_func(Wagner_bacterial_ps, Wagner_phy99)
Wagner_phy97 <- OCC_Abund_func(Wagner_bacterial_ps, Wagner_phy97)
Wagner_rhz99 <- OCC_Abund_func(Wagner_bacterial_ps, Wagner_rhz99)
Wagner_rhz97 <- OCC_Abund_func(Wagner_bacterial_ps, Wagner_rhz97)


# Output stats
summarize_isolate_output(Wagner_phy99, summarised = TRUE)
summarize_isolate_output(Wagner_phy97, summarised = TRUE)
summarize_isolate_output(Wagner_rhz99, summarised = TRUE)
summarize_isolate_output(Wagner_rhz97, summarised = TRUE)

# > summarize_isolate_output(Wagner_phy99, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     64.8
# 2 Yes                    35.0
# > summarize_isolate_output(Wagner_phy97, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     61.2
# 2 Yes                    38.6
# > summarize_isolate_output(Wagner_rhz99, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     59.0
# 2 Yes                    40.8
# > summarize_isolate_output(Wagner_rhz97, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     51.2
# 2 Yes                    48.6


###### 3) He et al. 2020 ######
# Load ASV lists and add prefix
He_rhz99 <- read.csv("ASV_LISTS/99/He_rhizosphere_ASVs.csv")
He_rhz97 <- read.csv("ASV_LISTS/97/He_rhizosphere_ASVs.csv")
He_rot99 <- read.csv("ASV_LISTS/99/He_root_ASVs.csv")
He_rot97 <- read.csv("ASV_LISTS/97/He_root_ASVs.csv")

He_rhz99$ASV_id <- paste0("He_", He_rhz99$ASV_id)
He_rhz97$ASV_id <- paste0("He_", He_rhz97$ASV_id)
He_rot99$ASV_id <- paste0("He_", He_rot99$ASV_id)
He_rot97$ASV_id <- paste0("He_", He_rot97$ASV_id)

# run function for O/A table
He_rhz99 <- OCC_Abund_func(He_bacterial_ps, He_rhz99)
He_rhz97 <- OCC_Abund_func(He_bacterial_ps, He_rhz97)
He_rot99 <- OCC_Abund_func(He_bacterial_ps, He_rot99)
He_rot97 <- OCC_Abund_func(He_bacterial_ps, He_rot97)

# Output stats
summarize_isolate_output(He_rhz99, summarised = TRUE)
summarize_isolate_output(He_rhz97, summarised = TRUE)
summarize_isolate_output(He_rot99, summarised = TRUE)
summarize_isolate_output(He_rot97, summarised = TRUE)

# > summarize_isolate_output(He_rhz99, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
# 1 No                     83.1
# 2 Yes                    14.4
# > summarize_isolate_output(He_rhz97, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
# 1 No                     61.1
# 2 Yes                    36.4
# > summarize_isolate_output(He_rot99, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
# 1 No                     82.5
# 2 Yes                    15.0
# > summarize_isolate_output(He_rot97, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
# 1 No                     60.4
# 2 Yes                    37.1


###### 4) Naylor et al. 2017 ######
# Load ASV lists and add prefix
Naylor_rhz99 <- read.csv("ASV_LISTS/99/Naylor_rhizosphere_ASVs.csv")
Naylor_rhz97 <- read.csv("ASV_LISTS/97/Naylor_rhizosphere_ASVs.csv")
Naylor_rot99 <- read.csv("ASV_LISTS/99/Naylor_root_ASVs.csv")
Naylor_rot97 <- read.csv("ASV_LISTS/97/Naylor_root_ASVs.csv")

Naylor_rhz99$ASV_id <- paste0("Naylor_", Naylor_rhz99$ASV_id)
Naylor_rhz97$ASV_id <- paste0("Naylor_", Naylor_rhz97$ASV_id)
Naylor_rot99$ASV_id <- paste0("Naylor_", Naylor_rot99$ASV_id)
Naylor_rot97$ASV_id <- paste0("Naylor_", Naylor_rot97$ASV_id)

# run function for O/A table
Naylor_rhz99 <- OCC_Abund_func(Naylor_bacterial_ps, Naylor_rhz99)
Naylor_rhz97 <- OCC_Abund_func(Naylor_bacterial_ps, Naylor_rhz97)
Naylor_rot99 <- OCC_Abund_func(Naylor_bacterial_ps, Naylor_rot99)
Naylor_rot97 <- OCC_Abund_func(Naylor_bacterial_ps, Naylor_rot97)

# Output stats
summarize_isolate_output(Naylor_rhz99, summarised = TRUE)
summarize_isolate_output(Naylor_rhz97, summarised = TRUE)
summarize_isolate_output(Naylor_rot99, summarised = TRUE)
summarize_isolate_output(Naylor_rot97, summarised = TRUE)

# > # Output stats
#   > summarize_isolate_output(Naylor_rhz99, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                    90.8 
# 2 Yes                    8.85
# > summarize_isolate_output(Naylor_rhz97, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     88.3
# 2 Yes                    11.4
# > summarize_isolate_output(Naylor_rot99, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                    90.5 
# 2 Yes                    9.12
# > summarize_isolate_output(Naylor_rot97, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     88.0
# 2 Yes                    11.7

###### 5) GfellerCadot et al. 2023 ######
# Load ASV Lists
GfellerCadot_soi99 <- read.csv("ASV_LISTS/99/GfellerCadot_soil_ASVs.csv")
GfellerCadot_soi97 <- read.csv("ASV_LISTS/97/GfellerCadot_soil_ASVs.csv")
GfellerCadot_rhz99 <- read.csv("ASV_LISTS/99/GfellerCadot_rhizosphere_ASVs.csv")
GfellerCadot_rhz97 <- read.csv("ASV_LISTS/97/GfellerCadot_rhizosphere_ASVs.csv")
GfellerCadot_rot99 <- read.csv("ASV_LISTS/99/GfellerCadot_root_ASVs.csv")
GfellerCadot_rot97 <- read.csv("ASV_LISTS/97/GfellerCadot_root_ASVs.csv")

GfellerCadot_soi99$ASV_id <- paste0("GfellerCadot_", GfellerCadot_soi99$ASV_id)
GfellerCadot_soi97$ASV_id <- paste0("GfellerCadot_", GfellerCadot_soi97$ASV_id)
GfellerCadot_rhz99$ASV_id <- paste0("GfellerCadot_", GfellerCadot_rhz99$ASV_id)
GfellerCadot_rhz97$ASV_id <- paste0("GfellerCadot_", GfellerCadot_rhz97$ASV_id)
GfellerCadot_rot99$ASV_id <- paste0("GfellerCadot_", GfellerCadot_rot99$ASV_id)
GfellerCadot_rot97$ASV_id <- paste0("GfellerCadot_", GfellerCadot_rot97$ASV_id)

# run function for O/A table
GfellerCadot_soi99 <- OCC_Abund_func(GfellerCadot_bacterial_ps, GfellerCadot_soi99)
GfellerCadot_soi97 <- OCC_Abund_func(GfellerCadot_bacterial_ps, GfellerCadot_soi97)
GfellerCadot_rhz99 <- OCC_Abund_func(GfellerCadot_bacterial_ps, GfellerCadot_rhz99)
GfellerCadot_rhz97 <- OCC_Abund_func(GfellerCadot_bacterial_ps, GfellerCadot_rhz97)
GfellerCadot_rot99 <- OCC_Abund_func(GfellerCadot_bacterial_ps, GfellerCadot_rot99)
GfellerCadot_rot97 <- OCC_Abund_func(GfellerCadot_bacterial_ps, GfellerCadot_rot97)

# Output stats
summarize_isolate_output(GfellerCadot_soi99, summarised = TRUE)
summarize_isolate_output(GfellerCadot_soi97, summarised = TRUE)
summarize_isolate_output(GfellerCadot_rhz99, summarised = TRUE)
summarize_isolate_output(GfellerCadot_rhz97, summarised = TRUE)
summarize_isolate_output(GfellerCadot_rot99, summarised = TRUE)
summarize_isolate_output(GfellerCadot_rot97, summarised = TRUE)

# > # Output stats
#   > summarize_isolate_output(GfellerCadot_soi99, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     82.7
# 2 Yes                    17.1
# > summarize_isolate_output(GfellerCadot_soi97, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     68.2
# 2 Yes                    31.6
# > summarize_isolate_output(GfellerCadot_rhz99, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     81.4
# 2 Yes                    18.4
# > summarize_isolate_output(GfellerCadot_rhz97, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     66.0
# 2 Yes                    33.8
# > summarize_isolate_output(GfellerCadot_rot99, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     81.4
# 2 Yes                    18.4
# > summarize_isolate_output(GfellerCadot_rot97, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     65.6
# 2 Yes                    34.2


###### 6) Gfeller et al. 2023 ######
# Load ASV lists and add prefix
Gfeller_soi99 <- read.csv("ASV_LISTS/99/Gfeller_soil_ASVs.csv")
Gfeller_soi97 <- read.csv("ASV_LISTS/97/Gfeller_soil_ASVs.csv")
Gfeller_rhz99 <- read.csv("ASV_LISTS/99/Gfeller_rhizosphere_ASVs.csv")
Gfeller_rhz97 <- read.csv("ASV_LISTS/97/Gfeller_rhizosphere_ASVs.csv")
Gfeller_rot99 <- read.csv("ASV_LISTS/99/Gfeller_root_ASVs.csv")
Gfeller_rot97 <- read.csv("ASV_LISTS/97/Gfeller_root_ASVs.csv")

Gfeller_soi99$ASV_id <- paste0("Gfeller_", Gfeller_soi99$ASV_id)
Gfeller_soi97$ASV_id <- paste0("Gfeller_", Gfeller_soi97$ASV_id)
Gfeller_rhz99$ASV_id <- paste0("Gfeller_", Gfeller_rhz99$ASV_id)
Gfeller_rhz97$ASV_id <- paste0("Gfeller_", Gfeller_rhz97$ASV_id)
Gfeller_rot99$ASV_id <- paste0("Gfeller_", Gfeller_rot99$ASV_id)
Gfeller_rot97$ASV_id <- paste0("Gfeller_", Gfeller_rot97$ASV_id)

# run function for O/A table
Gfeller_soi99 <- OCC_Abund_func(Gfeller_bacterial_ps, Gfeller_soi99)
Gfeller_soi97 <- OCC_Abund_func(Gfeller_bacterial_ps, Gfeller_soi97)
Gfeller_rhz99 <- OCC_Abund_func(Gfeller_bacterial_ps, Gfeller_rhz99)
Gfeller_rhz97 <- OCC_Abund_func(Gfeller_bacterial_ps, Gfeller_rhz97)
Gfeller_rot99 <- OCC_Abund_func(Gfeller_bacterial_ps, Gfeller_rot99)
Gfeller_rot97 <- OCC_Abund_func(Gfeller_bacterial_ps, Gfeller_rot97)

# Output stats
summarize_isolate_output(Gfeller_soi99, summarised = TRUE)
summarize_isolate_output(Gfeller_soi97, summarised = TRUE)
summarize_isolate_output(Gfeller_rhz99, summarised = TRUE)
summarize_isolate_output(Gfeller_rhz97, summarised = TRUE)
summarize_isolate_output(Gfeller_rot99, summarised = TRUE)
summarize_isolate_output(Gfeller_rot97, summarised = TRUE)

# > # Output stats
#   > summarize_isolate_output(Gfeller_soi99, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     89.1
# 2 Yes                    10.9
# > summarize_isolate_output(Gfeller_soi97, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     71.2
# 2 Yes                    28.7
# > summarize_isolate_output(Gfeller_rhz99, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     86.3
# 2 Yes                    13.6
# > summarize_isolate_output(Gfeller_rhz97, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     65.9
# 2 Yes                    34.1
# > summarize_isolate_output(Gfeller_rot99, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     86.6
# 2 Yes                    13.4
# > summarize_isolate_output(Gfeller_rot97, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     66.5
# 2 Yes                    33.4

###### 7) Cadot et al. 2021 ######
# Load ASV Lists
Cadot_soi99 <- read.csv("ASV_LISTS/99/Cadot_soil_ASVs.csv")
Cadot_soi97 <- read.csv("ASV_LISTS/97/Cadot_soil_ASVs.csv")
Cadot_rhz99 <- read.csv("ASV_LISTS/99/Cadot_rhizosphere_ASVs.csv")
Cadot_rhz97 <- read.csv("ASV_LISTS/97/Cadot_rhizosphere_ASVs.csv")
Cadot_rot99 <- read.csv("ASV_LISTS/99/Cadot_root_ASVs.csv")
Cadot_rot97 <- read.csv("ASV_LISTS/97/Cadot_root_ASVs.csv")

Cadot_soi99$ASV_id <- paste0("Cadot_", Cadot_soi99$ASV_id)
Cadot_soi97$ASV_id <- paste0("Cadot_", Cadot_soi97$ASV_id)
Cadot_rhz99$ASV_id <- paste0("Cadot_", Cadot_rhz99$ASV_id)
Cadot_rhz97$ASV_id <- paste0("Cadot_", Cadot_rhz97$ASV_id)
Cadot_rot99$ASV_id <- paste0("Cadot_", Cadot_rot99$ASV_id)
Cadot_rot97$ASV_id <- paste0("Cadot_", Cadot_rot97$ASV_id)

# Run function for O/A table
Cadot_soi99 <- OCC_Abund_func(Cadot_bacterial_ps, Cadot_soi99)
Cadot_soi97 <- OCC_Abund_func(Cadot_bacterial_ps, Cadot_soi97)
Cadot_rhz99 <- OCC_Abund_func(Cadot_bacterial_ps, Cadot_rhz99)
Cadot_rhz97 <- OCC_Abund_func(Cadot_bacterial_ps, Cadot_rhz97)
Cadot_rot99 <- OCC_Abund_func(Cadot_bacterial_ps, Cadot_rot99)
Cadot_rot97 <- OCC_Abund_func(Cadot_bacterial_ps, Cadot_rot97)

# Output stats
summarize_isolate_output(Cadot_soi99, summarised = TRUE)
summarize_isolate_output(Cadot_soi97, summarised = TRUE)
summarize_isolate_output(Cadot_rhz99, summarised = TRUE)
summarize_isolate_output(Cadot_rhz97, summarised = TRUE)
summarize_isolate_output(Cadot_rot99, summarised = TRUE)
summarize_isolate_output(Cadot_rot97, summarised = TRUE)

# > # Output stats
#   > summarize_isolate_output(Cadot_soi99, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     74.2
# 2 Yes                    25.8
# > summarize_isolate_output(Cadot_soi97, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     57.7
# 2 Yes                    42.4
# > summarize_isolate_output(Cadot_rhz99, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     74.2
# 2 Yes                    25.8
# > summarize_isolate_output(Cadot_rhz97, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     57.7
# 2 Yes                    42.4
# > summarize_isolate_output(Cadot_rot99, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     74.2
# 2 Yes                    25.8
# > summarize_isolate_output(Cadot_rot97, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     57.7
# 2 Yes                    42.4


###### 8) Hu et al. 2018 ######
# Load ASV Lists
Hu_soi99 <- read.csv("ASV_LISTS/99/Hu_soil_ASVs.csv")
Hu_soi97 <- read.csv("ASV_LISTS/97/Hu_soil_ASVs.csv")
Hu_rhz99 <- read.csv("ASV_LISTS/99/Hu_rhizosphere_ASVs.csv")
Hu_rhz97 <- read.csv("ASV_LISTS/97/Hu_rhizosphere_ASVs.csv")
Hu_rot99 <- read.csv("ASV_LISTS/99/Hu_root_ASVs.csv")
Hu_rot97 <- read.csv("ASV_LISTS/97/Hu_root_ASVs.csv")

Hu_soi99$ASV_id <- paste0("Cadot_", Hu_soi99$ASV_id)
Hu_soi97$ASV_id <- paste0("Cadot_", Hu_soi97$ASV_id)
Hu_rhz99$ASV_id <- paste0("Cadot_", Hu_rhz99$ASV_id)
Hu_rhz97$ASV_id <- paste0("Cadot_", Hu_rhz97$ASV_id)
Hu_rot99$ASV_id <- paste0("Cadot_", Hu_rot99$ASV_id)
Hu_rot97$ASV_id <- paste0("Cadot_", Hu_rot97$ASV_id)

# run function for O/A table
Hu_soi99 <- OCC_Abund_func(Hu_bacterial_ps, Hu_soi99)
Hu_soi97 <- OCC_Abund_func(Hu_bacterial_ps, Hu_soi97)
Hu_rhz99 <- OCC_Abund_func(Hu_bacterial_ps, Hu_rhz99)
Hu_rhz97 <- OCC_Abund_func(Hu_bacterial_ps, Hu_rhz97)
Hu_rot99 <- OCC_Abund_func(Hu_bacterial_ps, Hu_rot99)
Hu_rot97 <- OCC_Abund_func(Hu_bacterial_ps, Hu_rot97)

# Output stats
summarize_isolate_output(Hu_soi99, summarised = TRUE)
summarize_isolate_output(Hu_soi97, summarised = TRUE)
summarize_isolate_output(Hu_rhz99, summarised = TRUE)
summarize_isolate_output(Hu_rhz97, summarised = TRUE)
summarize_isolate_output(Hu_rot99, summarised = TRUE)
summarize_isolate_output(Hu_rot97, summarised = TRUE)

# > # Output stats
#   > summarize_isolate_output(Hu_soi99, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     85.4
# 2 Yes                    13.6
# > summarize_isolate_output(Hu_soi97, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     61.7
# 2 Yes                    37.3
# > summarize_isolate_output(Hu_rhz99, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     85.4
# 2 Yes                    13.6
# > summarize_isolate_output(Hu_rhz97, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     61.7
# 2 Yes                    37.4
# > summarize_isolate_output(Hu_rot99, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     85.4
# 2 Yes                    13.6
# > summarize_isolate_output(Hu_rot97, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     61.7
# 2 Yes                    37.4


###### 9) Walters et al. 2018 ######
# Load ASV Lists
Walters_rhz99 <- read.csv("ASV_LISTS/99/Walters_rhizosphere_ASVs.csv")
Walters_rhz97 <- read.csv("ASV_LISTS/97/Walters_rhizosphere_ASVs.csv")
Walters_soi99 <- read.csv("ASV_LISTS/99/Walters_soil_ASVs.csv")
Walters_soi97 <- read.csv("ASV_LISTS/97/Walters_soil_ASVs.csv")

Walters_rhz99$ASV_id <- paste0("Walters_", Walters_rhz99$ASV_id)
Walters_rhz97$ASV_id <- paste0("Walters_", Walters_rhz97$ASV_id)
Walters_soi99$ASV_id <- paste0("Walters_", Walters_soi99$ASV_id)
Walters_soi97$ASV_id <- paste0("Walters_", Walters_soi97$ASV_id)


# run function for O/A table
Walters_rhz99 <- OCC_Abund_func(Walters_bacterial_ps, Walters_rhz99)
Walters_rhz97 <- OCC_Abund_func(Walters_bacterial_ps, Walters_rhz97)
Walters_soi99 <- OCC_Abund_func(Walters_bacterial_ps, Walters_soi99)
Walters_soi97 <- OCC_Abund_func(Walters_bacterial_ps, Walters_soi97)

# Output stats
summarize_isolate_output(Walters_rhz99, summarised = TRUE)
summarize_isolate_output(Walters_rhz97, summarised = TRUE)
summarize_isolate_output(Walters_soi99, summarised = TRUE)
summarize_isolate_output(Walters_soi97, summarised = TRUE)

#> summarize_isolate_output(Walters_rhz99, summarised = TRUE)
## A tibble: 2 × 2
#InCollection TotalAbundance
#<chr>                 <dbl>
#  1 No                     74.7
#2 Yes                    18.4
#> summarize_isolate_output(Walters_rhz97, summarised = TRUE)
## A tibble: 2 × 2
#InCollection TotalAbundance
#<chr>                 <dbl>
#  1 No                     46.1
#2 Yes                    47.0
# > summarize_isolate_output(Walters_soi99, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     76.7
# 2 Yes                    16.4
# > summarize_isolate_output(Walters_soi97, summarised = TRUE)
# # A tibble: 2 × 2
# InCollection TotalAbundance
# <chr>                 <dbl>
#   1 No                     53.0
# 2 Yes                    40.1
