# Distributions of subsamples of MRB Isolates

library(tidyverse)
library(patchwork)
theme_set(theme_classic())

# load result file
MAZ <- read.csv("99-97_results_MRBcollection.csv")
## Combine with a column to denote source species
MAZ$species <- as.factor("MRB")
## test plot
collection_palette <- c(MRB = "black")


#### Both 99 and 97 Plot ####
## Swift 
a <- MAZ %>%
  filter(study == "Swift" & version == 99) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 66.8, color = "ZeaMic"), linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red"))

b <- MAZ %>%
  filter(study == "Swift" & version == 97) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 92.0, color = "ZeaMic"), linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red"))

swift <- (a + b)
swift


## Naylor 
a <- MAZ %>%
  filter(study == "Naylor" & version == 99) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 9.4, color = "ZeaMic"), linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red"))

b <- MAZ %>%
  filter(study == "Naylor" & version == 97) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 12.2, color = "ZeaMic"), linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red"))

naylor <- (a + b)
naylor

## He 
a <- MAZ %>%
  filter(study == "He" & version == 99) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 15.2, color = "ZeaMic"), linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red"))

b <- MAZ %>%
  filter(study == "He" & version == 97) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 37.4, color = "ZeaMic"), linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red"))

he <- (a + b)
he

### Wagner
a <- MAZ %>% filter(study == "Wagner" & version == 99) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 65.4, color = "ZeaMic"), linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red"))

b <- MAZ %>% filter(study == "Wagner" & version == 97) %>% 
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 77.5, color = "ZeaMic"), linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red"))

wagner <- (a + b)
wagner


## Meier 
a <- MAZ %>%
  filter(study == "Meier" & version == 99) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 29.7, color = "ZeaMic"), linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red"))

b <- MAZ %>%
  filter(study == "Meier" & version == 97) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 43.8, color = "ZeaMic"), linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red"))

meier <- (a + b)
meier


## Walters 
a <- MAZ %>%
  filter(study == "Walters" & version == 99) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 19.9, color = "ZeaMic"), linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red"))

b <- MAZ %>%
  filter(study == "Walters" & version == 97) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 49.0, color = "ZeaMic"), linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red"))

walters <- (a + b)
walters


###Gfeller and cadot 2023
a <- MAZ %>%
  filter(study == "GfellerCadot" & version == 99) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 20.6,  color = "ZeaMic"), linetype = "dashed", linewidth = 1) + 
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red"))

b <- MAZ %>% 
  filter(study == "GfellerCadot" & version == 97) %>% 
  ggplot(., aes(x = TotalAbundance, fill = species)) + 
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") + 
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") + 
  geom_vline(aes(xintercept = 37.0, color = "ZeaMic"), linetype = "dashed", linewidth = 1) + 
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red"))

gfellercadot <- (a + b)
gfellercadot


###Gfeller 2023
a <- MAZ %>%
  filter(study == "Gfeller" & version == 99) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 15.0,  color = "ZeaMic"), linetype = "dashed", linewidth = 1) + 
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red"))

b <- MAZ %>% 
  filter(study == "Gfeller" & version == 97) %>% 
  ggplot(., aes(x = TotalAbundance, fill = species)) + 
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") + 
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") + 
  geom_vline(aes(xintercept = 36.1, color = "ZeaMic"), linetype = "dashed", linewidth = 1) + 
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red"))

gfeller <- (a + b)
gfeller


###Cadot
a <- MAZ %>%
  filter(study == "Cadot" & version == 99) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 28.9,  color = "ZeaMic"), linetype = "dashed", linewidth = 1) + 
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red"))

b <- MAZ %>% 
  filter(study == "Cadot" & version == 97) %>% 
  ggplot(., aes(x = TotalAbundance, fill = species)) + 
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") + 
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") + 
  geom_vline(aes(xintercept = 44.8, color = "ZeaMic"), linetype = "dashed", linewidth = 1) + 
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red"))

cadot <- (a + b)
cadot


###Hu
a <- MAZ %>%
  filter(study == "Hu" & version == 99) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 14.5,  color = "ZeaMic"), linetype = "dashed", linewidth = 1) + 
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red"))

b <- MAZ %>% 
  filter(study == "Hu" & version == 97) %>% 
  ggplot(., aes(x = TotalAbundance, fill = species)) + 
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") + 
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") + 
  geom_vline(aes(xintercept = 37.9, color = "ZeaMic"), linetype = "dashed", linewidth = 1) + 
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red"))

hu <- (a + b)
hu



# Final Combination 
final_plot <- wrap_plots(
  swift,
  naylor,
  he,
  wagner,
  meier,
  walters,
  gfellercadot,
  gfeller,
  cadot,
  hu,
  ncol = 2,
  guides = "collect")&
  theme(legend.position = "top",
        legend.key = element_rect(fill = "#b3b3b3ff",color = NA),
        legend.background = element_rect(fill = "#b3b3b3ff", color = "black"))


final_plot
ggsave("Distrutions_99-97_MRB.svg", final_plot, width = 30, height = 30, unit = "cm")



#### 99 only plot ####
## Swift
swift99 <- MAZ %>%
  filter(study == "Swift" & version == 99) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 66.8, color = "ZeaMic"), linetype = "solid", linewidth = 1) +
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red")) +
  scale_x_continuous(n.breaks = 4, limits = c(0, 75), labels = scales::label_number(accuracy = 1))
## Naylor 
naylor99 <- MAZ %>%
  filter(study == "Naylor" & version == 99) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 9.4, color = "ZeaMic"), linetype = "solid", linewidth = 1) +
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red")) +
  scale_x_continuous(n.breaks = 3, limits = c(0, 25), labels = scales::label_number(accuracy = 1))
## He 
he99 <- MAZ %>%
  filter(study == "He" & version == 99) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 15.2, color = "ZeaMic"), linetype = "solid", linewidth = 1) +
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red")) +
  scale_x_continuous(n.breaks = 3, limits = c(0, 25), labels = scales::label_number(accuracy = 1))
### Wagner
wagner99 <- MAZ %>% filter(study == "Wagner" & version == 99) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 65.4, color = "ZeaMic"), linetype = "solid", linewidth = 1) +
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red")) +
  scale_x_continuous(n.breaks = 4, limits = c(0, 75), labels = scales::label_number(accuracy = 1))
## Meier 
meier99 <- MAZ %>%
  filter(study == "Meier" & version == 99) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 29.7, color = "ZeaMic"), linetype = "solid", linewidth = 1) +
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red")) +
  scale_x_continuous(n.breaks = 3, limits = c(0, 50), labels = scales::label_number(accuracy = 1))
## Walters 
walters99 <- MAZ %>%
  filter(study == "Walters" & version == 99) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 19.9, color = "ZeaMic"), linetype = "solid", linewidth = 1) +
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red")) +
  scale_x_continuous(n.breaks = 3, limits = c(0, 50), labels = scales::label_number(accuracy = 1))
###Gfeller and cadot 2023
gfellercadot99 <- MAZ %>%
  filter(study == "GfellerCadot" & version == 99) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 20.6,  color = "ZeaMic"), linetype = "solid", linewidth = 1) + 
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red")) +
  scale_x_continuous(n.breaks = 3, limits = c(0, 25), labels = scales::label_number(accuracy = 1))
###Gfeller 2023
gfeller99 <- MAZ %>%
  filter(study == "Gfeller" & version == 99) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 15.0,  color = "ZeaMic"), linetype = "solid", linewidth = 1) + 
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red")) +
  scale_x_continuous(n.breaks = 3, limits = c(0, 25), labels = scales::label_number(accuracy = 1))
###Cadot
cadot99 <- MAZ %>%
  filter(study == "Cadot" & version == 99) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 28.9,  color = "ZeaMic"), linetype = "solid", linewidth = 1) + 
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red")) +
  scale_x_continuous(n.breaks = 3, limits = c(0, 50), labels = scales::label_number(accuracy = 1))
###Hu
hu99 <- MAZ %>%
  filter(study == "Hu" & version == 99) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 14.5,  color = "ZeaMic"), linetype = "solid", linewidth = 1) + 
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red")) +
  scale_x_continuous(n.breaks = 3, limits = c(0, 25), labels = scales::label_number(accuracy = 1))
# Final Combination 
final_plot99 <- wrap_plots(
  cadot99,
  gfellercadot99,
  gfeller99,
  he99,
  hu99,
  meier99,
  naylor99,
  swift99,
  wagner99,
  walters99,
  nrow = 2,
  guides = "collect") &
  theme(legend.position = "top",
        legend.key = element_rect(fill = "#b3b3b3ff",color = NA),
        legend.background = element_rect(fill = "#b3b3b3ff", color = "black"))
final_plot99
ggsave("Figure panels/Figure7A.svg", final_plot99, width = 24, height = 12, unit = "cm")


#### 97 only plot ####
## Swift 
swift97 <- MAZ %>%
  filter(study == "Swift" & version == 97) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 92, color = "ZeaMic"), linetype = "solid", linewidth = 1) +
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red")) +
  scale_x_continuous(n.breaks = 4, limits = c(0, 100), labels = scales::label_number(accuracy = 1))
## Naylor 
naylor97 <- MAZ %>%
  filter(study == "Naylor" & version == 97) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 12.2, color = "ZeaMic"), linetype = "solid", linewidth = 1) +
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red")) +
  scale_x_continuous(n.breaks = 3, limits = c(0, 25), labels = scales::label_number(accuracy = 1))
## He 
he97 <- MAZ %>%
  filter(study == "He" & version == 97) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 37.4, color = "ZeaMic"), linetype = "solid", linewidth = 1) +
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red")) +
  scale_x_continuous(n.breaks = 3, limits = c(0, 50), labels = scales::label_number(accuracy = 1))
### Wagner
wagner97 <- MAZ %>% filter(study == "Wagner" & version == 97) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 77.5, color = "ZeaMic"), linetype = "solid", linewidth = 1) +
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red")) +
  scale_x_continuous(n.breaks = 4, limits = c(0, 100), labels = scales::label_number(accuracy = 1))
## Meier 
meier97 <- MAZ %>%
  filter(study == "Meier" & version == 97) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 43.8, color = "ZeaMic"), linetype = "solid", linewidth = 1) +
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red")) +
  scale_x_continuous(n.breaks = 3, limits = c(0, 50), labels = scales::label_number(accuracy = 1))
## Walters 
walters97 <- MAZ %>%
  filter(study == "Walters" & version == 97) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 49.0, color = "ZeaMic"), linetype = "solid", linewidth = 1) +
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red")) +
  scale_x_continuous(n.breaks = 3, limits = c(0, 50), labels = scales::label_number(accuracy = 1))
###Gfeller and cadot 2023
gfellercadot97 <- MAZ %>%
  filter(study == "GfellerCadot" & version == 97) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 37.0,  color = "ZeaMic"), linetype = "solid", linewidth = 1) + 
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red")) +
  scale_x_continuous(n.breaks = 3, limits = c(0, 50), labels = scales::label_number(accuracy = 1))
###Gfeller 2023
gfeller97 <- MAZ %>%
  filter(study == "Gfeller" & version == 97) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 36.1,  color = "ZeaMic"), linetype = "solid", linewidth = 1) + 
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red")) +
  scale_x_continuous(n.breaks = 3, limits = c(0, 50), labels = scales::label_number(accuracy = 1))
###Cadot
cadot97 <- MAZ %>%
  filter(study == "Cadot" & version == 97) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 44.8,  color = "ZeaMic"), linetype = "solid", linewidth = 1) + 
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red")) +
  scale_x_continuous(n.breaks = 3, limits = c(0, 50), labels = scales::label_number(accuracy = 1))
###Hu
hu97 <- MAZ %>%
  filter(study == "Hu" & version == 97) %>%
  ggplot(., aes(x = TotalAbundance, fill = species)) +
  geom_histogram(position = "identity", binwidth = 0.5, color = "black") +
  facet_wrap(~interaction(study, version, sep = " - ")) +
  labs(x = "Relative Abundance", y = "Frequency") +
  scale_fill_manual(values = collection_palette, name = "Collection") +
  geom_vline(aes(xintercept = 37.9,  color = "ZeaMic"), linetype = "solid", linewidth = 1) + 
  scale_color_manual(name = "Full Collection", values = c("ZeaMic" = "red")) +
  scale_x_continuous(n.breaks = 3, limits = c(0, 50), labels = scales::label_number(accuracy = 1))
# Final Combination 
final_plot97 <- wrap_plots(
  cadot97,
  gfellercadot97,
  gfeller97,
  he97,
  hu97,
  meier97,
  naylor97,
  swift97,
  wagner97,
  walters97,
  nrow = 2,
  guides = "collect") &
  theme(legend.position = "top",
        legend.key = element_rect(fill = "#b3b3b3ff",color = NA),
        legend.background = element_rect(fill = "#b3b3b3ff", color = "black"))



final_plot97
ggsave("Supplement/FigureS8.svg", final_plot97, width = 24, height = 12, unit = "cm")
ggsave("Supplement/FigureS8.png", final_plot97, width = 24, height = 12, unit = "cm")

