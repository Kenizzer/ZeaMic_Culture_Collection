# Parse Blast alignment results
## Generate ASV lists of taxa with alignment to zeamic for other scripts
## Presence/absence plots for supplement

# Load pacakges and set ggplot theme
library(tidyverse)
library(patchwork)
theme_set(theme_classic())

## load ZeaMic Isolates
all_isolates_df <- read.csv("../ZeaMic Isolates with Full Taxonomy.csv")
all_isolates_df$X <- NULL

#### load blast results 99 ####
blast_df <- read.delim("../../BLAST_ASVs-from-studies/99results.txt", header = FALSE, sep = "\t")

# Cleanup blast dataframe
colnames(blast_df) <- c("ZeaMic_Isolate", "ASV", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# Extract info from ASV column 
blast_df <- blast_df %>% extract(ASV, into = c("study", "ASV_id", "sample_type", "country", "region"),
                     regex = "(.+)_(ASV_\\d+)-([^-]+)-([^_]+)_(.+)",
                     remove = TRUE)
# Rename studies 
blast_df <- blast_df %>% mutate(study = dplyr::recode(study, 
    Naylor = "Naylor et al. (2017)",
    Wagner = "Wagner et al. (2020)", 
    Meier  = "Meier et al. (2022)", 
    He = "He et al. (2024)", 
    GfellerCadot = "Gfeller and Cadot et al. (2023)", 
    Gfeller = "Gfeller et al. (2023)",
    Cadot = "Cadot et al. (2021)",
    Hu = "Hu et al. (2018)",
    Walters = "Walters et al. (2018)",
    Swift = "Swift et al. (2025)"))

## Presence and absence table + isolates not in blast output
presence_absence <- blast_df %>%
  # Get unique combinations of isolate, study, and sample_type
  distinct(ZeaMic_Isolate, study, sample_type) %>%
  # Add a presence indicator
  mutate(present = 1) %>%
  # Create a combined column name for study and sample_type
  unite("study_sample", study, sample_type, sep = "_") %>%
  # Pivot to wide format
  pivot_wider(names_from = study_sample, 
              values_from = present,
              values_fill = 0)

presence_absence_complete <- all_isolates_df %>%
  left_join(presence_absence, by = "ZeaMic_Isolate") %>%
  # Replace NA with 0 for isolates not in BLAST results
  mutate(across(-ZeaMic_Isolate, ~replace_na(., 0)))

# Print table
print(presence_absence_complete)

# Convert the wide format back to long format for ggplot2
heatmap_data <- presence_absence_complete %>%
  pivot_longer(cols = -c(ZeaMic_Isolate, phylum, class, order, family, genus),
               names_to = "study_sample",
               values_to = "presence") %>%
  # Separate study and sample type for better visualization
  separate(study_sample, 
           into = c("study", "sample_type"), 
           sep = "_(?=[^_]+$)")  # Split at the last underscore

# Get number of isolate with significant alignments per study
summary(as.factor(heatmap_data$study))
heatmap_data %>% filter(study == "Naylor et al. (2017)") %>% group_by(sample_type) %>% summarise(total_isolates = sum(presence))
heatmap_data %>% filter(study == "He et al. (2024)") %>% group_by(sample_type) %>% summarise(total_isolates = sum(presence))
heatmap_data %>% filter(study == "Wagner et al. (2020)") %>% group_by(sample_type) %>% summarise(total_isolates = sum(presence))
heatmap_data %>% filter(study == "Meier et al. (2022)") %>% group_by(sample_type) %>% summarise(total_isolates = sum(presence))
heatmap_data %>% filter(study == "Walters et al. (2018)") %>% group_by(sample_type) %>% summarise(total_isolates = sum(presence))
heatmap_data %>% filter(study == "Gfeller et al. (2023)") %>% group_by(sample_type) %>% summarise(total_isolates = sum(presence))
heatmap_data %>% filter(study == "Gfeller and Cadot et al. (2023)") %>% group_by(sample_type) %>% summarise(total_isolates = sum(presence))
heatmap_data %>% filter(study == "Cadot et al. (2021)") %>% group_by(sample_type) %>% summarise(total_isolates = sum(presence))
heatmap_data %>% filter(study == "Hu et al. (2018)") %>% group_by(sample_type) %>% summarise(total_isolates = sum(presence))
heatmap_data %>% filter(study == "Swift et al. (2025)") %>% group_by(sample_type) %>% summarise(total_isolates = sum(presence))

#### get list of ASVs to map back to individual study ####
ASV_presence <- blast_df %>%
  # Get unique combinations of isolate, study, sample_type, and ASV
  distinct(ZeaMic_Isolate, study, sample_type, ASV_id) %>%
  # Add a presence indicator
  mutate(present = 1) %>%
  # Create a combined column name for study and sample_type
  unite("study_sample", study, sample_type, sep = "_") %>%
  # Pivot to wide format
  pivot_wider(names_from = study_sample, 
              values_from = present,
              values_fill = 0)

# 99

#Naylor
Naylor_root_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Naylor et al. (2017)_ROOT`) %>%
  filter(`Naylor et al. (2017)_ROOT` == 1) %>%
  select(ZeaMic_Isolate, ASV_id)
Naylor_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Naylor et al. (2017)_RHIZOSPHERE`) %>%
  filter(`Naylor et al. (2017)_RHIZOSPHERE` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
write.csv(Naylor_root_ASVs, "99/Naylor_root_ASVs.csv", row.names = FALSE)
write.csv(Naylor_rhizosphere_ASVs, "99/Naylor_rhizosphere_ASVs.csv", row.names = FALSE)


#Meier
Meier_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Meier et al. (2022)_RHIZOSPHERE`) %>%
  filter(`Meier et al. (2022)_RHIZOSPHERE` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
write.csv(Meier_rhizosphere_ASVs, "99/Meier_rhizosphere_ASVs.csv", row.names = FALSE)


#Wagner
Wagner_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Wagner et al. (2020)_RHIZOSPHERE`) %>%
  filter(`Wagner et al. (2020)_RHIZOSPHERE` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
Wagner_phyllosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Wagner et al. (2020)_PHYLLOSPHERE`) %>%
  filter(`Wagner et al. (2020)_PHYLLOSPHERE` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
write.csv(Wagner_rhizosphere_ASVs, "99/Wagner_rhizosphere_ASVs.csv", row.names = FALSE)
write.csv(Wagner_phyllosphere_ASVs, "99/Wagner_phyllosphere_ASVs.csv", row.names = FALSE)

#He et al. (2024)
He_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `He et al. (2024)_RHIZOSPHERE`) %>%
  filter(`He et al. (2024)_RHIZOSPHERE` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
He_root_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `He et al. (2024)_ROOT`) %>%
  filter(`He et al. (2024)_ROOT` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
write.csv(He_rhizosphere_ASVs, "99/He_rhizosphere_ASVs.csv", row.names = FALSE)
write.csv(He_root_ASVs, "99/He_root_ASVs.csv", row.names = FALSE)

#GfellerCadot
colnames(ASV_presence)
GfellerCadot_soil_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Gfeller and Cadot et al. (2023)_SOIL`) %>%
  filter(`Gfeller and Cadot et al. (2023)_SOIL` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
GfellerCadot_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Gfeller and Cadot et al. (2023)_RHIZOSPHERE`) %>%
  filter(`Gfeller and Cadot et al. (2023)_RHIZOSPHERE` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
GfellerCadot_root_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Gfeller and Cadot et al. (2023)_ROOT`) %>%
  filter(`Gfeller and Cadot et al. (2023)_ROOT` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
write.csv(GfellerCadot_soil_ASVs, "99/GfellerCadot_soil_ASVs.csv", row.names = FALSE)
write.csv(GfellerCadot_rhizosphere_ASVs, "99/GfellerCadot_rhizosphere_ASVs.csv", row.names = FALSE)
write.csv(GfellerCadot_root_ASVs, "99/GfellerCadot_root_ASVs.csv", row.names = FALSE)

#Gfeller
colnames(ASV_presence)
Gfeller_soil_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Gfeller et al. (2023)_SOIL`) %>%
  filter(`Gfeller et al. (2023)_SOIL` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
Gfeller_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Gfeller et al. (2023)_RHIZOSPHERE`) %>%
  filter(`Gfeller et al. (2023)_RHIZOSPHERE` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
Gfeller_root_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Gfeller et al. (2023)_ROOT`) %>%
  filter(`Gfeller et al. (2023)_ROOT` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
write.csv(Gfeller_soil_ASVs, "99/Gfeller_soil_ASVs.csv", row.names = FALSE)
write.csv(Gfeller_rhizosphere_ASVs, "99/Gfeller_rhizosphere_ASVs.csv", row.names = FALSE)
write.csv(Gfeller_root_ASVs, "99/Gfeller_root_ASVs.csv", row.names = FALSE)


#Cadot
colnames(ASV_presence)
Cadot_soil_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Cadot et al. (2021)_SOIL`) %>%
  filter(`Cadot et al. (2021)_SOIL` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
Cadot_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Cadot et al. (2021)_RHIZOSPHERE`) %>%
  filter(`Cadot et al. (2021)_RHIZOSPHERE` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
Cadot_root_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Cadot et al. (2021)_ROOT`) %>%
  filter(`Cadot et al. (2021)_ROOT` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
write.csv(Cadot_soil_ASVs, "99/Cadot_soil_ASVs.csv", row.names = FALSE)
write.csv(Cadot_rhizosphere_ASVs, "99/Cadot_rhizosphere_ASVs.csv", row.names = FALSE)
write.csv(Cadot_root_ASVs, "99/Cadot_root_ASVs.csv", row.names = FALSE)

#Hu
colnames(ASV_presence)
Hu_soil_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Hu et al. (2018)_SOIL`) %>%
  filter(`Hu et al. (2018)_SOIL` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
Hu_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Hu et al. (2018)_RHIZOSPHERE`) %>%
  filter(`Hu et al. (2018)_RHIZOSPHERE` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
Hu_root_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Hu et al. (2018)_ROOT`) %>%
  filter(`Hu et al. (2018)_ROOT` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
write.csv(Hu_soil_ASVs, "99/Hu_soil_ASVs.csv", row.names = FALSE)
write.csv(Hu_rhizosphere_ASVs, "99/Hu_rhizosphere_ASVs.csv", row.names = FALSE)
write.csv(Hu_root_ASVs, "99/Hu_root_ASVs.csv", row.names = FALSE)

#walters
colnames(ASV_presence)
Walters_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Walters et al. (2018)_RHIZOSPHERE`) %>%
  filter(`Walters et al. (2018)_RHIZOSPHERE` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
Walters_soil_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Walters et al. (2018)_SOIL`) %>%
  filter(`Walters et al. (2018)_SOIL` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
write.csv(Walters_rhizosphere_ASVs, "99/Walters_rhizosphere_ASVs.csv", row.names = FALSE)
write.csv(Walters_soil_ASVs, "99/Walters_soil_ASVs.csv", row.names = FALSE)

#Swift
colnames(ASV_presence)
Swift_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Swift et al. (2025)_ROOT`) %>%
  filter(`Swift et al. (2025)_ROOT` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
write.csv(Swift_rhizosphere_ASVs, "99/Swift_root_ASVs.csv", row.names = FALSE)




#### 97 ####
blast_df <- read.delim("../../BLAST_ASVs-from-studies/97results.txt", header = FALSE, sep = "\t")

# Cleanup blast dataframe
colnames(blast_df) <- c("ZeaMic_Isolate", "ASV", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# Extract info from ASV column 
blast_df <- blast_df %>% extract(ASV, into = c("study", "ASV_id", "sample_type", "country", "region"),
                                 regex = "(.+)_(ASV_\\d+)-([^-]+)-([^_]+)_(.+)",
                                 remove = TRUE)
# Rename studies 
blast_df <- blast_df %>% mutate(study = dplyr::recode(study, 
                                                      Naylor = "Naylor et al. (2017)",
                                                      Wagner = "Wagner et al. (2020)", 
                                                      Meier  = "Meier et al. (2022)", 
                                                      He = "He et al. (2024)", 
                                                      GfellerCadot = "Gfeller and Cadot et al. (2023)", 
                                                      Gfeller = "Gfeller et al. (2023)",
                                                      Cadot = "Cadot et al. (2021)",
                                                      Hu = "Hu et al. (2018)",
                                                      Walters = "Walters et al. (2018)",
                                                      Swift = "Swift et al. (2025)"))

## Presence and absence table + isolates not in blast output
presence_absence <- blast_df %>%
  # Get unique combinations of isolate, study, and sample_type
  distinct(ZeaMic_Isolate, study, sample_type) %>%
  # Add a presence indicator
  mutate(present = 1) %>%
  # Create a combined column name for study and sample_type
  unite("study_sample", study, sample_type, sep = "_") %>%
  # Pivot to wide format
  pivot_wider(names_from = study_sample, 
              values_from = present,
              values_fill = 0)

presence_absence_complete <- all_isolates_df %>%
  left_join(presence_absence, by = "ZeaMic_Isolate") %>%
  # Replace NA with 0 for isolates not in BLAST results
  mutate(across(-ZeaMic_Isolate, ~replace_na(., 0)))

# Print table
print(presence_absence_complete)

# Convert the wide format back to long format for ggplot2
heatmap_data <- presence_absence_complete %>%
  pivot_longer(cols = -c(ZeaMic_Isolate, genus, phylum, class, order, family),
               names_to = "study_sample",
               values_to = "presence") %>%
  # Separate study and sample type for better visualization
  separate(study_sample, 
           into = c("study", "sample_type"), 
           sep = "_(?=[^_]+$)")  # Split at the last underscore

# Get number of isolate with significant alignments per study
summary(as.factor(heatmap_data$study))
heatmap_data %>% filter(study == "Naylor et al. (2017)") %>% group_by(sample_type) %>% summarise(total_isolates = sum(presence))
heatmap_data %>% filter(study == "He et al. (2024)") %>% group_by(sample_type) %>% summarise(total_isolates = sum(presence))
heatmap_data %>% filter(study == "Wagner et al. (2020)") %>% group_by(sample_type) %>% summarise(total_isolates = sum(presence))
heatmap_data %>% filter(study == "Meier et al. (2022)") %>% group_by(sample_type) %>% summarise(total_isolates = sum(presence))
heatmap_data %>% filter(study == "Walters et al. (2018)") %>% group_by(sample_type) %>% summarise(total_isolates = sum(presence))
heatmap_data %>% filter(study == "Gfeller et al. (2023)") %>% group_by(sample_type) %>% summarise(total_isolates = sum(presence))
heatmap_data %>% filter(study == "Gfeller and Cadot et al. (2023)") %>% group_by(sample_type) %>% summarise(total_isolates = sum(presence))
heatmap_data %>% filter(study == "Cadot et al. (2021)") %>% group_by(sample_type) %>% summarise(total_isolates = sum(presence))
heatmap_data %>% filter(study == "Hu et al. (2018)") %>% group_by(sample_type) %>% summarise(total_isolates = sum(presence))
heatmap_data %>% filter(study == "Swift et al. (2025)") %>% group_by(sample_type) %>% summarise(total_isolates = sum(presence))

### get list of ASVs to map back to individual study ###
ASV_presence <- blast_df %>%
  # Get unique combinations of isolate, study, sample_type, and ASV
  distinct(ZeaMic_Isolate, study, sample_type, ASV_id) %>%
  # Add a presence indicator
  mutate(present = 1) %>%
  # Create a combined column name for study and sample_type
  unite("study_sample", study, sample_type, sep = "_") %>%
  # Pivot to wide format
  pivot_wider(names_from = study_sample, 
              values_from = present,
              values_fill = 0)

#Naylor
Naylor_root_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Naylor et al. (2017)_ROOT`) %>%
  filter(`Naylor et al. (2017)_ROOT` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
Naylor_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Naylor et al. (2017)_RHIZOSPHERE`) %>%
  filter(`Naylor et al. (2017)_RHIZOSPHERE` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
write.csv(Naylor_root_ASVs, "97/Naylor_root_ASVs.csv", row.names = FALSE)
write.csv(Naylor_rhizosphere_ASVs, "97/Naylor_rhizosphere_ASVs.csv", row.names = FALSE)


#Meier
Meier_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Meier et al. (2022)_RHIZOSPHERE`) %>%
  filter(`Meier et al. (2022)_RHIZOSPHERE` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
write.csv(Meier_rhizosphere_ASVs, "97/Meier_rhizosphere_ASVs.csv", row.names = FALSE)


#Wagner
Wagner_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Wagner et al. (2020)_RHIZOSPHERE`) %>%
  filter(`Wagner et al. (2020)_RHIZOSPHERE` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
Wagner_phyllosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Wagner et al. (2020)_PHYLLOSPHERE`) %>%
  filter(`Wagner et al. (2020)_PHYLLOSPHERE` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
write.csv(Wagner_rhizosphere_ASVs, "97/Wagner_rhizosphere_ASVs.csv", row.names = FALSE)
write.csv(Wagner_phyllosphere_ASVs, "97/Wagner_phyllosphere_ASVs.csv", row.names = FALSE)

#He et al. (2024)"
He_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `He et al. (2024)_RHIZOSPHERE`) %>%
  filter(`He et al. (2024)_RHIZOSPHERE` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
He_root_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `He et al. (2024)_ROOT`) %>%
  filter(`He et al. (2024)_ROOT` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
write.csv(He_rhizosphere_ASVs, "97/He_rhizosphere_ASVs.csv", row.names = FALSE)
write.csv(He_root_ASVs, "97/He_root_ASVs.csv", row.names = FALSE)

#GfellerCadot
colnames(ASV_presence)
GfellerCadot_soil_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Gfeller and Cadot et al. (2023)_SOIL`) %>%
  filter(`Gfeller and Cadot et al. (2023)_SOIL` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
GfellerCadot_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Gfeller and Cadot et al. (2023)_RHIZOSPHERE`) %>%
  filter(`Gfeller and Cadot et al. (2023)_RHIZOSPHERE` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
GfellerCadot_root_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Gfeller and Cadot et al. (2023)_ROOT`) %>%
  filter(`Gfeller and Cadot et al. (2023)_ROOT` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
write.csv(GfellerCadot_rhizosphere_ASVs, "97/GfellerCadot_rhizosphere_ASVs.csv", row.names = FALSE)
write.csv(GfellerCadot_root_ASVs, "97/GfellerCadot_root_ASVs.csv", row.names = FALSE)
write.csv(GfellerCadot_soil_ASVs, "97/GfellerCadot_soil_ASVs.csv", row.names = FALSE)

#Gfeller
colnames(ASV_presence)
Gfeller_soil_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Gfeller et al. (2023)_SOIL`) %>%
  filter(`Gfeller et al. (2023)_SOIL` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
Gfeller_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Gfeller et al. (2023)_RHIZOSPHERE`) %>%
  filter(`Gfeller et al. (2023)_RHIZOSPHERE` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
Gfeller_root_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Gfeller et al. (2023)_ROOT`) %>%
  filter(`Gfeller et al. (2023)_ROOT` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
write.csv(Gfeller_soil_ASVs, "97/Gfeller_soil_ASVs.csv", row.names = FALSE)
write.csv(Gfeller_rhizosphere_ASVs, "97/Gfeller_rhizosphere_ASVs.csv", row.names = FALSE)
write.csv(Gfeller_root_ASVs, "97/Gfeller_root_ASVs.csv", row.names = FALSE)


#Cadot
colnames(ASV_presence)
Cadot_soil_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Cadot et al. (2021)_SOIL`) %>%
  filter(`Cadot et al. (2021)_SOIL` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
Cadot_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Cadot et al. (2021)_RHIZOSPHERE`) %>%
  filter(`Cadot et al. (2021)_RHIZOSPHERE` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
Cadot_root_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Cadot et al. (2021)_ROOT`) %>%
  filter(`Cadot et al. (2021)_ROOT` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
write.csv(Cadot_soil_ASVs, "97/Cadot_soil_ASVs.csv", row.names = FALSE)
write.csv(Cadot_rhizosphere_ASVs, "97/Cadot_rhizosphere_ASVs.csv", row.names = FALSE)
write.csv(Cadot_root_ASVs, "97/Cadot_root_ASVs.csv", row.names = FALSE)

#Hu
colnames(ASV_presence)
Hu_soil_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Hu et al. (2018)_SOIL`) %>%
  filter(`Hu et al. (2018)_SOIL` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
Hu_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Hu et al. (2018)_RHIZOSPHERE`) %>%
  filter(`Hu et al. (2018)_RHIZOSPHERE` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
Hu_root_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Hu et al. (2018)_ROOT`) %>%
  filter(`Hu et al. (2018)_ROOT` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
write.csv(Hu_soil_ASVs, "97/Hu_soil_ASVs.csv", row.names = FALSE)
write.csv(Hu_rhizosphere_ASVs, "97/Hu_rhizosphere_ASVs.csv", row.names = FALSE)
write.csv(Hu_root_ASVs, "97/Hu_root_ASVs.csv", row.names = FALSE)

#walters
colnames(ASV_presence)
Walters_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Walters et al. (2018)_RHIZOSPHERE`) %>%
  filter(`Walters et al. (2018)_RHIZOSPHERE` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
Walters_soil_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Walters et al. (2018)_SOIL`) %>%
  filter(`Walters et al. (2018)_SOIL` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
write.csv(Walters_rhizosphere_ASVs, "97/Walters_rhizosphere_ASVs.csv", row.names = FALSE)
write.csv(Walters_soil_ASVs, "97/Walters_soil_ASVs.csv", row.names = FALSE)

#Swift
colnames(ASV_presence)
Swift_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Swift et al. (2025)_ROOT`) %>%
  filter(`Swift et al. (2025)_ROOT` == 1) %>%
  select(ZeaMic_Isolate ,ASV_id)
write.csv(Swift_rhizosphere_ASVs, "97/Swift_root_ASVs.csv", row.names = FALSE)
