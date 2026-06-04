# Parse Blast alignment results
## Generate ASV lists of taxa with alignment to zeamic for other scripts
## Presence/absence plots for supplement

# Load pacakges and set ggplot theme
library(tidyverse)
library(patchwork)
theme_set(theme_classic())

## load ZeaMic Isolates
all_isolates_df <- read.csv("ZeaMic Isolates with Full Taxonomy.csv")
all_isolates_df$X <- NULL

#### load blast results 99 ####
blast_df <- read.delim("../BLAST_ASVs-from-studies/99results.txt", header = FALSE, sep = "\t")

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

##### Plotting isolate and genera P/A clustered by taxonomy #####
# Arrange the isolates by class
isolate_99_data <- heatmap_data %>%
  arrange(class, ZeaMic_Isolate) %>%
  mutate(ZeaMic_Isolate = factor(ZeaMic_Isolate, levels = unique(ZeaMic_Isolate)))
# Arrange the genera by class
genera_99_data <- heatmap_data %>%
  arrange(class, genus) %>%
  mutate(genus = factor(genus, levels = unique(genus)))
# Class barplot levels
# Class barplots for both plots, as we organize the genus and isolate by the factor level class
class_df_iso <- isolate_99_data %>% distinct(ZeaMic_Isolate, class)
class_df_gen <- genera_99_data %>% distinct(genus, class)

### Plots
## Very minimal class barplot, will be attached to other plots
class_bar_iso <- ggplot(class_df_iso, aes(x = 1, y = ZeaMic_Isolate, fill = class)) +
  geom_tile() +
  scale_fill_manual(values = c("Actinomycetes"= "#23733a",
                              "Bacilli" = "#b8d437",
                              "Deinococci" = "#ffc323",
                              "Chitinophagia" = "#f0852f",
                              "Flavobacteriia" = "#e51f41",
                              "Alphaproteobacteria" = "#c6529f",
                              "Betaproteobacteria" = "#405caa",
                              "Gammaproteobacteria" = "#46bbe7"), name = "Class") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = levels(heatmap_data$ZeaMic_Isolate), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid = element_blank(),
        panel.border = element_blank(), 
        plot.margin = margin(t = 5, r = 2, b = 5, l = 5), 
        legend.position = "right")

class_bar_gen <- ggplot(class_df_gen, aes(x = 1, y = genus, fill = class)) +
  geom_tile() +
  scale_fill_manual(values = c("Actinomycetes"= "#23733a",
                               "Bacilli" = "#b8d437",
                               "Deinococci" = "#ffc323",
                               "Chitinophagia" = "#f0852f",
                               "Flavobacteriia" = "#e51f41",
                               "Alphaproteobacteria" = "#c6529f",
                               "Betaproteobacteria" = "#405caa",
                               "Gammaproteobacteria" = "#46bbe7"), name = "Class") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = levels(heatmap_data$genus), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid = element_blank(),
        panel.border = element_blank(), 
        plot.margin = margin(t = 5, r = 2, b = 5, l = 5), 
        legend.position = "right")

## Isolate heatmap
isolate_plot <- ggplot(isolate_99_data, aes(x = sample_type, y = ZeaMic_Isolate, fill = factor(presence))) +
  geom_tile(color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("0" = "white", "1" = "#50C878"),
                    labels = c("Absent", "Present"),
                    name = "Status") +
  facet_wrap(~study, scales = "free_x", nrow = 1, labeller = label_wrap_gen(width = 12)) +
  scale_y_discrete(limits = levels(heatmap_data$ZeaMic_Isolate), expand = c(0, 0)) +
  labs(x = "Sample Type", y = "Isolate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(face = "bold", size = 6),
        panel.grid = element_blank(), 
        axis.text.y = element_text(size = 5),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 2))

## Genus heatmap
genus_plot <- ggplot(genera_99_data, aes(x = sample_type, y = genus, fill = factor(presence))) +
  geom_tile(color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("0" = "white", "1" = "#50C878"),
                    labels = c("Absent", "Present"),
                    name = "Status") +
  facet_wrap(~study, scales = "free_x", nrow = 1, labeller = label_wrap_gen(width = 12)) +
  scale_y_discrete(limits = levels(heatmap_data$ZeaMic_Isolate), expand = c(0, 0)) +
  labs(x = "Sample Type", y = "Genus") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(face = "bold", size = 6),
        panel.grid = element_blank(), 
        axis.text.y = element_text(size = 8),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 2))

### Arranging
isolate_plot_save <- class_bar_iso + isolate_plot + plot_layout(widths = c(0.2, 10), guides = "collect") &
  theme(legend.position = "top", legend.title = element_text(size = 8), legend.text = element_text(size = 6), legend.key.size = unit(0.5, "cm"))

ggsave("Presence_absence_plots/99perc Isolate Plot.svg", isolate_plot_save, height = 23, width = 18.4, units = "cm")
ggsave("Presence_absence_plots/99perc Isolate Plot.png", isolate_plot_save, height = 23, width = 18.4, units = "cm")

genus_plot_save <- class_bar_gen + genus_plot + plot_layout(widths = c(0.2, 10), guides = "collect") &
  theme(legend.position = "top", legend.title = element_text(size = 8), legend.text = element_text(size = 6), legend.key.size = unit(0.5, "cm"))

ggsave("Presence_absence_plots/99perc Genus Plot.svg", genus_plot_save, height = 23, width = 18.4, units = "cm")
ggsave("Presence_absence_plots/99perc Genus Plot.png", genus_plot_save, height = 23, width = 18.4, units = "cm")

# Clean env
rm(class_df_gen, class_df_iso, class_bar_iso, class_bar_gen, genus_plot,
   genus_plot_save, isolate_plot, isolate_plot_save,isolate_99_data, genera_99_data)


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
  select(ASV_id)
Naylor_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Naylor et al. (2017)_RHIZOSPHERE`) %>%
  filter(`Naylor et al. (2017)_RHIZOSPHERE` == 1) %>%
  select(ASV_id)
write.csv(Naylor_root_ASVs, "ASV_LISTS/99/Naylor_root_ASVs.csv", row.names = FALSE)
write.csv(Naylor_rhizosphere_ASVs, "ASV_LISTS/99/Naylor_rhizosphere_ASVs.csv", row.names = FALSE)


#Meier
Meier_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Meier et al. (2022)_RHIZOSPHERE`) %>%
  filter(`Meier et al. (2022)_RHIZOSPHERE` == 1) %>%
  select(ASV_id)
write.csv(Meier_rhizosphere_ASVs, "ASV_LISTS/99/Meier_rhizosphere_ASVs.csv", row.names = FALSE)


#Wagner
Wagner_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Wagner et al. (2020)_RHIZOSPHERE`) %>%
  filter(`Wagner et al. (2020)_RHIZOSPHERE` == 1) %>%
  select(ASV_id)
Wagner_phyllosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Wagner et al. (2020)_PHYLLOSPHERE`) %>%
  filter(`Wagner et al. (2020)_PHYLLOSPHERE` == 1) %>%
  select(ASV_id)
write.csv(Wagner_rhizosphere_ASVs, "ASV_LISTS/99/Wagner_rhizosphere_ASVs.csv", row.names = FALSE)
write.csv(Wagner_phyllosphere_ASVs, "ASV_LISTS/99/Wagner_phyllosphere_ASVs.csv", row.names = FALSE)

#He et al. (2024)
He_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `He et al. (2024)_RHIZOSPHERE`) %>%
  filter(`He et al. (2024)_RHIZOSPHERE` == 1) %>%
  select(ASV_id)
He_root_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `He et al. (2024)_ROOT`) %>%
  filter(`He et al. (2024)_ROOT` == 1) %>%
  select(ASV_id)
write.csv(He_rhizosphere_ASVs, "ASV_LISTS/99/He_rhizosphere_ASVs.csv", row.names = FALSE)
write.csv(He_root_ASVs, "ASV_LISTS/99/He_root_ASVs.csv", row.names = FALSE)

#GfellerCadot
colnames(ASV_presence)
GfellerCadot_soil_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Gfeller and Cadot et al. (2023)_SOIL`) %>%
  filter(`Gfeller and Cadot et al. (2023)_SOIL` == 1) %>%
  select(ASV_id)
GfellerCadot_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Gfeller and Cadot et al. (2023)_RHIZOSPHERE`) %>%
  filter(`Gfeller and Cadot et al. (2023)_RHIZOSPHERE` == 1) %>%
  select(ASV_id)
GfellerCadot_root_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Gfeller and Cadot et al. (2023)_ROOT`) %>%
  filter(`Gfeller and Cadot et al. (2023)_ROOT` == 1) %>%
  select(ASV_id)
write.csv(GfellerCadot_soil_ASVs, "ASV_LISTS/99/GfellerCadot_soil_ASVs.csv", row.names = FALSE)
write.csv(GfellerCadot_rhizosphere_ASVs, "ASV_LISTS/99/GfellerCadot_rhizosphere_ASVs.csv", row.names = FALSE)
write.csv(GfellerCadot_root_ASVs, "ASV_LISTS/99/GfellerCadot_root_ASVs.csv", row.names = FALSE)

#Gfeller
colnames(ASV_presence)
Gfeller_soil_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Gfeller et al. (2023)_SOIL`) %>%
  filter(`Gfeller et al. (2023)_SOIL` == 1) %>%
  select(ASV_id)
Gfeller_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Gfeller et al. (2023)_RHIZOSPHERE`) %>%
  filter(`Gfeller et al. (2023)_RHIZOSPHERE` == 1) %>%
  select(ASV_id)
Gfeller_root_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Gfeller et al. (2023)_ROOT`) %>%
  filter(`Gfeller et al. (2023)_ROOT` == 1) %>%
  select(ASV_id)
write.csv(Gfeller_soil_ASVs, "ASV_LISTS/99/Gfeller_soil_ASVs.csv", row.names = FALSE)
write.csv(Gfeller_rhizosphere_ASVs, "ASV_LISTS/99/Gfeller_rhizosphere_ASVs.csv", row.names = FALSE)
write.csv(Gfeller_root_ASVs, "ASV_LISTS/99/Gfeller_root_ASVs.csv", row.names = FALSE)


#Cadot
colnames(ASV_presence)
Cadot_soil_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Cadot et al. (2021)_SOIL`) %>%
  filter(`Cadot et al. (2021)_SOIL` == 1) %>%
  select(ASV_id)
Cadot_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Cadot et al. (2021)_RHIZOSPHERE`) %>%
  filter(`Cadot et al. (2021)_RHIZOSPHERE` == 1) %>%
  select(ASV_id)
Cadot_root_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Cadot et al. (2021)_ROOT`) %>%
  filter(`Cadot et al. (2021)_ROOT` == 1) %>%
  select(ASV_id)
write.csv(Cadot_soil_ASVs, "ASV_LISTS/99/Cadot_soil_ASVs.csv", row.names = FALSE)
write.csv(Cadot_rhizosphere_ASVs, "ASV_LISTS/99/Cadot_rhizosphere_ASVs.csv", row.names = FALSE)
write.csv(Cadot_root_ASVs, "ASV_LISTS/99/Cadot_root_ASVs.csv", row.names = FALSE)

#Hu
colnames(ASV_presence)
Hu_soil_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Hu et al. (2018)_SOIL`) %>%
  filter(`Hu et al. (2018)_SOIL` == 1) %>%
  select(ASV_id)
Hu_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Hu et al. (2018)_RHIZOSPHERE`) %>%
  filter(`Hu et al. (2018)_RHIZOSPHERE` == 1) %>%
  select(ASV_id)
Hu_root_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Hu et al. (2018)_ROOT`) %>%
  filter(`Hu et al. (2018)_ROOT` == 1) %>%
  select(ASV_id)
write.csv(Hu_soil_ASVs, "ASV_LISTS/99/Hu_soil_ASVs.csv", row.names = FALSE)
write.csv(Hu_rhizosphere_ASVs, "ASV_LISTS/99/Hu_rhizosphere_ASVs.csv", row.names = FALSE)
write.csv(Hu_root_ASVs, "ASV_LISTS/99/Hu_root_ASVs.csv", row.names = FALSE)

#walters
colnames(ASV_presence)
Walters_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Walters et al. (2018)_RHIZOSPHERE`) %>%
  filter(`Walters et al. (2018)_RHIZOSPHERE` == 1) %>%
  select(ASV_id)
Walters_soil_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Walters et al. (2018)_SOIL`) %>%
  filter(`Walters et al. (2018)_SOIL` == 1) %>%
  select(ASV_id)
write.csv(Walters_rhizosphere_ASVs, "ASV_LISTS/99/Walters_rhizosphere_ASVs.csv", row.names = FALSE)
write.csv(Walters_soil_ASVs, "ASV_LISTS/99/Walters_soil_ASVs.csv", row.names = FALSE)

#Swift
colnames(ASV_presence)
Swift_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Swift et al. (2025)_ROOT`) %>%
  filter(`Swift et al. (2025)_ROOT` == 1) %>%
  select(ASV_id)
write.csv(Swift_rhizosphere_ASVs, "ASV_LISTS/99/Swift_root_ASVs.csv", row.names = FALSE)




#### 97 ####
blast_df <- read.delim("../BLAST_ASVs-from-studies/97results.txt", header = FALSE, sep = "\t")

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

##### Plotting isolate and genera P/A clustered by taxonomy #####
# Arrange the isolates by class
isolate_97_data <- heatmap_data %>%
  arrange(class, ZeaMic_Isolate) %>%
  mutate(ZeaMic_Isolate = factor(ZeaMic_Isolate, levels = unique(ZeaMic_Isolate)))
# Arrange the genera by class
genera_97_data <- heatmap_data %>%
  arrange(class, genus) %>%
  mutate(genus = factor(genus, levels = unique(genus)))
# Class barplot levels
# Class barplots the same for both plots, as we organize the genus and isolate
## by the factor level class
## by the factor level class
class_df_iso <- isolate_97_data %>% distinct(ZeaMic_Isolate, class)
class_df_gen <- genera_97_data %>% distinct(genus, class)

### Plots
## Very minimal class barplot, will be attached to other plots
class_bar_iso <- ggplot(class_df_iso, aes(x = 1, y = ZeaMic_Isolate, fill = class)) +
  geom_tile() +
  scale_fill_manual(values = c("Actinomycetes"= "#23733a",
                               "Bacilli" = "#b8d437",
                               "Deinococci" = "#ffc323",
                               "Chitinophagia" = "#f0852f",
                               "Flavobacteriia" = "#e51f41",
                               "Alphaproteobacteria" = "#c6529f",
                               "Betaproteobacteria" = "#405caa",
                               "Gammaproteobacteria" = "#46bbe7"), name = "Class") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = levels(heatmap_data$ZeaMic_Isolate), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid = element_blank(),
        panel.border = element_blank(), 
        plot.margin = margin(t = 5, r = 2, b = 5, l = 5), 
        legend.position = "right")

class_bar_gen <- ggplot(class_df_gen, aes(x = 1, y = genus, fill = class)) +
  geom_tile() +
  scale_fill_manual(values = c("Actinomycetes"= "#23733a",
                               "Bacilli" = "#b8d437",
                               "Deinococci" = "#ffc323",
                               "Chitinophagia" = "#f0852f",
                               "Flavobacteriia" = "#e51f41",
                               "Alphaproteobacteria" = "#c6529f",
                               "Betaproteobacteria" = "#405caa",
                               "Gammaproteobacteria" = "#46bbe7"), name = "Class") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = levels(heatmap_data$genus), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid = element_blank(),
        panel.border = element_blank(), 
        plot.margin = margin(t = 5, r = 2, b = 5, l = 5), 
        legend.position = "right")

## Isolate heatmap
isolate_plot <- ggplot(isolate_97_data, aes(x = sample_type, y = ZeaMic_Isolate, fill = factor(presence))) +
  geom_tile(color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("0" = "white", "1" = "#50C878"),
                    labels = c("Absent", "Present"),
                    name = "Status") +
  facet_wrap(~study, scales = "free_x", nrow = 1, labeller = label_wrap_gen(width = 12)) +
  scale_y_discrete(limits = levels(heatmap_data$ZeaMic_Isolate), expand = c(0, 0)) +
  labs(x = "Sample Type", y = "Isolate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(face = "bold", size = 6),
        panel.grid = element_blank(), 
        axis.text.y = element_text(size = 5),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 2))

## Genus heatmap
genus_plot <- ggplot(genera_97_data, aes(x = sample_type, y = genus, fill = factor(presence))) +
  geom_tile(color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("0" = "white", "1" = "#50C878"),
                    labels = c("Absent", "Present"),
                    name = "Status") +
  facet_wrap(~study, scales = "free_x", nrow = 1, labeller = label_wrap_gen(width = 12)) +
  scale_y_discrete(limits = levels(heatmap_data$ZeaMic_Isolate), expand = c(0, 0)) +
  labs(x = "Sample Type", y = "Genus") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(face = "bold", size = 6),
        panel.grid = element_blank(), 
        axis.text.y = element_text(size = 8),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 2))

### Arranging
isolate_plot_save <- class_bar_iso + isolate_plot + plot_layout(widths = c(0.2, 10), guides = "collect") &
  theme(legend.position = "top", legend.title = element_text(size = 8), legend.text = element_text(size = 6), legend.key.size = unit(0.5, "cm"))

ggsave("Presence_absence_plots/97perc Isolate Plot.svg", isolate_plot_save, height = 23, width = 18.4, units = "cm")
ggsave("Presence_absence_plots/97perc Isolate Plot.png", isolate_plot_save, height = 23, width = 18.4, units = "cm")

genus_plot_save <- class_bar_gen + genus_plot + plot_layout(widths = c(0.2, 10), guides = "collect") &
  theme(legend.position = "top", legend.title = element_text(size = 8), legend.text = element_text(size = 6), legend.key.size = unit(0.5, "cm"))

ggsave("Presence_absence_plots/97perc Genus Plot.svg", genus_plot_save, height = 23, width = 18.4, units = "cm")
ggsave("Presence_absence_plots/97perc Genus Plot.png", genus_plot_save, height = 23, width = 18.4, units = "cm")

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
  select(ASV_id)
Naylor_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Naylor et al. (2017)_RHIZOSPHERE`) %>%
  filter(`Naylor et al. (2017)_RHIZOSPHERE` == 1) %>%
  select(ASV_id)
write.csv(Naylor_root_ASVs, "ASV_LISTS/97/Naylor_root_ASVs.csv", row.names = FALSE)
write.csv(Naylor_rhizosphere_ASVs, "ASV_LISTS/97/Naylor_rhizosphere_ASVs.csv", row.names = FALSE)


#Meier
Meier_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Meier et al. (2022)_RHIZOSPHERE`) %>%
  filter(`Meier et al. (2022)_RHIZOSPHERE` == 1) %>%
  select(ASV_id)
write.csv(Meier_rhizosphere_ASVs, "ASV_LISTS/97/Meier_rhizosphere_ASVs.csv", row.names = FALSE)


#Wagner
Wagner_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Wagner et al. (2020)_RHIZOSPHERE`) %>%
  filter(`Wagner et al. (2020)_RHIZOSPHERE` == 1) %>%
  select(ASV_id)
Wagner_phyllosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Wagner et al. (2020)_PHYLLOSPHERE`) %>%
  filter(`Wagner et al. (2020)_PHYLLOSPHERE` == 1) %>%
  select(ASV_id)
write.csv(Wagner_rhizosphere_ASVs, "ASV_LISTS/97/Wagner_rhizosphere_ASVs.csv", row.names = FALSE)
write.csv(Wagner_phyllosphere_ASVs, "ASV_LISTS/97/Wagner_phyllosphere_ASVs.csv", row.names = FALSE)

#He et al. (2024)"
He_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `He et al. (2024)_RHIZOSPHERE`) %>%
  filter(`He et al. (2024)_RHIZOSPHERE` == 1) %>%
  select(ASV_id)
He_root_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `He et al. (2024)_ROOT`) %>%
  filter(`He et al. (2024)_ROOT` == 1) %>%
  select(ASV_id)
write.csv(He_rhizosphere_ASVs, "ASV_LISTS/97/He_rhizosphere_ASVs.csv", row.names = FALSE)
write.csv(He_root_ASVs, "ASV_LISTS/97/He_root_ASVs.csv", row.names = FALSE)

#GfellerCadot
colnames(ASV_presence)
GfellerCadot_soil_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Gfeller and Cadot et al. (2023)_SOIL`) %>%
  filter(`Gfeller and Cadot et al. (2023)_SOIL` == 1) %>%
  select(ASV_id)
GfellerCadot_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Gfeller and Cadot et al. (2023)_RHIZOSPHERE`) %>%
  filter(`Gfeller and Cadot et al. (2023)_RHIZOSPHERE` == 1) %>%
  select(ASV_id)
GfellerCadot_root_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Gfeller and Cadot et al. (2023)_ROOT`) %>%
  filter(`Gfeller and Cadot et al. (2023)_ROOT` == 1) %>%
  select(ASV_id)
write.csv(GfellerCadot_rhizosphere_ASVs, "ASV_LISTS/97/GfellerCadot_rhizosphere_ASVs.csv", row.names = FALSE)
write.csv(GfellerCadot_root_ASVs, "ASV_LISTS/97/GfellerCadot_root_ASVs.csv", row.names = FALSE)
write.csv(GfellerCadot_soil_ASVs, "ASV_LISTS/97/GfellerCadot_soil_ASVs.csv", row.names = FALSE)

#Gfeller
colnames(ASV_presence)
Gfeller_soil_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Gfeller et al. (2023)_SOIL`) %>%
  filter(`Gfeller et al. (2023)_SOIL` == 1) %>%
  select(ASV_id)
Gfeller_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Gfeller et al. (2023)_RHIZOSPHERE`) %>%
  filter(`Gfeller et al. (2023)_RHIZOSPHERE` == 1) %>%
  select(ASV_id)
Gfeller_root_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Gfeller et al. (2023)_ROOT`) %>%
  filter(`Gfeller et al. (2023)_ROOT` == 1) %>%
  select(ASV_id)
write.csv(Gfeller_soil_ASVs, "ASV_LISTS/97/Gfeller_soil_ASVs.csv", row.names = FALSE)
write.csv(Gfeller_rhizosphere_ASVs, "ASV_LISTS/97/Gfeller_rhizosphere_ASVs.csv", row.names = FALSE)
write.csv(Gfeller_root_ASVs, "ASV_LISTS/97/Gfeller_root_ASVs.csv", row.names = FALSE)


#Cadot
colnames(ASV_presence)
Cadot_soil_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Cadot et al. (2021)_SOIL`) %>%
  filter(`Cadot et al. (2021)_SOIL` == 1) %>%
  select(ASV_id)
Cadot_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Cadot et al. (2021)_RHIZOSPHERE`) %>%
  filter(`Cadot et al. (2021)_RHIZOSPHERE` == 1) %>%
  select(ASV_id)
Cadot_root_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Cadot et al. (2021)_ROOT`) %>%
  filter(`Cadot et al. (2021)_ROOT` == 1) %>%
  select(ASV_id)
write.csv(Cadot_soil_ASVs, "ASV_LISTS/97/Cadot_soil_ASVs.csv", row.names = FALSE)
write.csv(Cadot_rhizosphere_ASVs, "ASV_LISTS/97/Cadot_rhizosphere_ASVs.csv", row.names = FALSE)
write.csv(Cadot_root_ASVs, "ASV_LISTS/97/Cadot_root_ASVs.csv", row.names = FALSE)

#Hu
colnames(ASV_presence)
Hu_soil_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Hu et al. (2018)_SOIL`) %>%
  filter(`Hu et al. (2018)_SOIL` == 1) %>%
  select(ASV_id)
Hu_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Hu et al. (2018)_RHIZOSPHERE`) %>%
  filter(`Hu et al. (2018)_RHIZOSPHERE` == 1) %>%
  select(ASV_id)
Hu_root_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Hu et al. (2018)_ROOT`) %>%
  filter(`Hu et al. (2018)_ROOT` == 1) %>%
  select(ASV_id)
write.csv(Hu_soil_ASVs, "ASV_LISTS/97/Hu_soil_ASVs.csv", row.names = FALSE)
write.csv(Hu_rhizosphere_ASVs, "ASV_LISTS/97/Hu_rhizosphere_ASVs.csv", row.names = FALSE)
write.csv(Hu_root_ASVs, "ASV_LISTS/97/Hu_root_ASVs.csv", row.names = FALSE)

#walters
colnames(ASV_presence)
Walters_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Walters et al. (2018)_RHIZOSPHERE`) %>%
  filter(`Walters et al. (2018)_RHIZOSPHERE` == 1) %>%
  select(ASV_id)
Walters_soil_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Walters et al. (2018)_SOIL`) %>%
  filter(`Walters et al. (2018)_SOIL` == 1) %>%
  select(ASV_id)
write.csv(Walters_rhizosphere_ASVs, "ASV_LISTS/97/Walters_rhizosphere_ASVs.csv", row.names = FALSE)
write.csv(Walters_soil_ASVs, "ASV_LISTS/97/Walters_soil_ASVs.csv", row.names = FALSE)

#Swift
colnames(ASV_presence)
Swift_rhizosphere_ASVs <- ASV_presence %>%
  select(ZeaMic_Isolate, ASV_id, `Swift et al. (2025)_ROOT`) %>%
  filter(`Swift et al. (2025)_ROOT` == 1) %>%
  select(ASV_id)
write.csv(Swift_rhizosphere_ASVs, "ASV_LISTS/97/Swift_root_ASVs.csv", row.names = FALSE)
