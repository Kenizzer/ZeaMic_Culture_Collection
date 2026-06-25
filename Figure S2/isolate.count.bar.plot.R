## Figure S2 ##
# Code by: Nichole Ginnan 2026

#Load necessary libraries
library(ggplot2)
library(readxl)
library(dplyr)
library(RColorBrewer)
library(ggbreak)

# Read the updated Excel file (adjust path as needed)
df <- read_excel("./Crude_collection_table.xlsx")
bubble_df <- df %>% count(Class, Genus, name = "Genus_count", sort = TRUE)
# Arrange by Class and count
bubble_df <- bubble_df %>%
  arrange(Class, desc(Genus_count)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))

# Plot
plot <- ggplot(bubble_df, aes(x = Genus, y = Genus_count, fill = Class)) +
  geom_col(color = "black") +
  scale_fill_manual(values = c(
    "Actinomycetes" ="#23733a",
    "Alphaproteobacteria"= "#c6529e",
    "Bacilli" ="#b8d437",
    "Betaproteobacteria"= "#3f5cab",
    "Chitinophagia"="#f0852d",
    "Gammaproteobacteria" ="#45bbe8",
    "Flavobacteriia" ="#e61f42")) +
  scale_y_break(c(15, 55), scales=1) +
  scale_y_break(c(125, 540), scales=0.2) +
  geom_text(
    aes(label = Genus_count),
    vjust = -0.2,
    size = 3
  ) +
  labs(
    x = "Genus",
    y = "Number of Isolates",
    fill = "Class") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "right",
        axis.text.y = element_text(face = "italic"))

ggsave("Figures2.svg", plot, dpi=300, width=250, height=100, unit = "mm")