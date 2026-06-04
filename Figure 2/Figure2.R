# Barplot and OA curve for Swift et al. 2025

# Packages w/ version numbers.
library(tidyverse); packageVersion('tidyverse')
library(phyloseq); packageVersion('phyloseq')
library(ggpubr); packageVersion('ggpubr')
library(ggtext); packageVersion('ggtext')
library(fantaxtic); packageVersion('fantaxtic')
library(vegan); packageVersion('vegan')


# Theme set
theme_set(theme_classic())

#### Load dataset for use throughout ####
bact.asv <- readRDS('phyloseq_b_asv_clean_allRoots.RDS')

#### Top 35 Genera Taxonomic Barplot ####
# genus glom 
bact.gen <- tax_glom(bact.asv, "Genus") # 133 taxa
# Relative Abundance Transform
bact.gen.relab <- transform_sample_counts(bact.gen, function(x) x/sum(x))
# Extract top 35 Genera
bact.gen.relab.top35 <- fantaxtic::get_top_taxa(physeq_obj = bact.gen.relab, n = 35, relative = TRUE, discard_other = FALSE, other_label = "Low Abundance Taxa")
# Melt data frame with phyloseq function for plotting.
bact.gen.relab.top35 <- psmelt(bact.gen.relab.top35)
# change Massilia to new genera name
bact.gen.relab.top35 <- bact.gen.relab.top35 %>% mutate(Genus = if_else(Genus == "Massilia", "Telluria", Genus))
rm(bact.gen,bact.gen.relab) # Dump intermediate objects

# Load ZeaMiC isolate genera list to assist in coloring plot
ZeaMic_isolate_genera<- data.frame(Genus = c("Enterobacter", "Enterobacter", "Enterobacter", "Monashia", "Burkholderia",
  "Burkholderia","Monashia", "Trinickia", "Luteibacter", "Leifsonia", "Monashia",
  "Burkholderia", "Mesorhizobium", "Burkholderia", "Labrys", "Microbacterium",
  "Luteibacter", "Paraburkholderia", "Leifsonia", "Leifsonia", "Burkholderia",
  "Paraburkholderia", "Burkholderia", "Rhizobium", "Paraburkholderia", "Rhizobium",
  "Pseudomonas", "Pseudomonas", "Chryseobacterium", "Cellulomonas", "Chryseobacterium",
  "Microbacterium", "Telluria", "Telluria", "Rahnella", "Enterobacter", "Arthrobacter",
  "Flavobacterium", "Brevibacillus", "Mitsuaria", "Rhodococcus", "Stenotrophomonas",
  "Pseudomonas", "Pseudomonas", "Stenotrophomonas", "Brucella", "Curtobacterium",
  "Enterobacter", "Chryseobacterium", "Herbaspirillum", "Pseudomonas", "Methylobacterium",
  "Arthrobacter", "Franconibacter", "Pantoea", "Kocuria", "Chitinophaga", "Achromobacter",
  "Achromobacter", "Luteibacter", "Bradyrhizobium", "Achromobacter", "Rhizobium",
  "Enterobacter", "Paraburkholderia", "Rhizobium", "Sinomonas", "Arthrobacter",
  "Agrobacterium", "Rhizobium", "Dyella", "Luteibacter", "Burkholderia",
  "Variovorax", "Burkholderia", "Rhizobium", "Burkholderia", "Dyella", "Burkholderia",
  "Rhizobium", "Variovorax", "Peribacillus", "Bacillus", "Deinococcus", "Streptomyces",
  "Nocardioides", "Rhizobium", "Sphingobium"))
# unique list
ZeaMic_isolate_genera <- ZeaMic_isolate_genera %>% distinct()
# Assign a new column to color plot if genera in collection
bact.gen.relab.top35 <- bact.gen.relab.top35 %>% 
  mutate(InCollection = ifelse(Genus %in% ZeaMic_isolate_genera$Genus, "Yes", "No"))
# cut to only necessary columns
bact.gen.relab.top35 <- bact.gen.relab.top35 %>%
  select(OTU, Drought.or.Watered, Sample,Abundance, Kingdom:InCollection)

# Set the plotting order levels (By in collection column, then alphabetical, with Low Abundance Taxa at the end)
genus_order <- bact.gen.relab.top35 %>%
  distinct(Genus, InCollection) %>%
  mutate(is_low = grepl("low", Genus, ignore.case = TRUE)) %>%
  arrange(is_low, desc(InCollection), Genus) %>%
  pull(Genus)

bact.gen.relab.top35 <- bact.gen.relab.top35 %>%
  mutate(Genus = factor(Genus, levels = genus_order))

# Set the colors and order
incollection_colors <- c("Yes" = "#7fbf7b","No"  = "#af8dc3", "Low" = "#000000")
bact.gen.relab.top35 <- bact.gen.relab.top35 %>%
  mutate(InCollection = case_when( grepl("low", Genus, ignore.case = TRUE) ~ "Low",TRUE ~ InCollection))  

bact.gen.relab.top35$InCollection <- factor(bact.gen.relab.top35$InCollection, levels = c("Yes", "No", "Low") )

# Genera color plot
library(Polychrome)
set.seed(101239811)
P35 <- createPalette(35,  c("#d53e4f", "#fc8d59", "#fee08b",
                            "#e6f598", "#99d594", "#3288bd"))
names(P35) <- NULL
P36 <- c(P35, "#000000")
rm(P35)

telluria_order <- bact.gen.relab.top35 %>%
  filter(Genus == "Telluria") %>%
  arrange(desc(Abundance))%>%
  pull(Sample)

bact.gen.relab.top35 <- bact.gen.relab.top35 %>% mutate(Sample = factor(Sample, levels = telluria_order))

barplot_Genera <- bact.gen.relab.top35 %>%
  ggplot(., aes(x=Sample, y=Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance") +
  xlab("sample") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values=P36) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right')


#### Occupancy Abundance Plot ######
# Load phyloseq dataset
Swift_bacterial_ps <- readRDS("phyloseq_b_asv_clean_allRoots.RDS")
# Load ASV Lists
Swift_list1 <- read.csv("Swift_root_ASVs97.csv")
Swift_asvs <- Swift_list1
# Add ASV prefix for study (first author first name capitalized)
Swift_asvs$ASV_id <- paste0("b", Swift_list1$ASV_id)

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
  # Occupancy and abundance
  otu_tab <- data.frame(otu_table(ps_rare))
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

# run function for O/A table
Swift_df <- OCC_Abund_func(Swift_bacterial_ps, Swift_asvs)

# Occupancy-Abundance
OCC_Abund_Gen <- ggplot(data = Swift_df$occ_abun_tax, aes(x = log10(otu_rel), y = otu_occ)) +
  geom_point(aes(fill = Genus_for_color), alpha = 0.8, shape = 21, size = 5, color = "black") +
  scale_fill_manual(values=P36, name = "Genus") +
  labs(x = expression(log[10]("mean Relative Abundance")), y = "Occupancy")


Tosave <- ggarrange(OCC_Abund_Gen, barplot_Genera, align = "hv", common.legend = TRUE, legend = "bottom")

ggsave("Figure2_draft.svg", Tosave, height = 18, width = 24, unit = "cm")

## Further edits were made manually in Inkscape.