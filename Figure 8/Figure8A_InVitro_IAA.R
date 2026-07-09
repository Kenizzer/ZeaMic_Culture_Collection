##### In Vivo ZeaMiC analysis (IAA Production) #####
## Code by: Joel Swift 2026

library(tidyverse)
library(patchwork)

theme_set(theme_classic())

## load data
IAA_data <- readxl::read_xlsx("./InVitro_IAA_production.xlsx")

## Convert any negative values to zeros
IAA_data <- IAA_data %>% mutate(`IAA (µg mL-1)` = if_else(`IAA (µg mL-1)` < 0, 0, `IAA (µg mL-1)`))


## Order isolates by mean IAA, but force NiuAA isolates to the end
iso_order <- IAA_data %>%
  group_by(`Isolate ID`) %>%
  summarize(mean_IAA = mean(`IAA (µg mL-1)`), .groups = "drop") %>%
  arrange(desc(mean_IAA)) %>%
  pull(`Isolate ID`) %>%
  {\(x) c(x[!str_detect(x, "^NiuAA")],
          x[str_detect(x, "^NiuAA")])}()

## Number of non-NiuAA isolates (used for shading)
n_other <- sum(!str_detect(iso_order, "^NiuAA"))

## Plot
IAA_PLOT <- IAA_data %>%
  mutate(Group = if_else(str_detect(`Isolate ID`, "^NiuAA"), "NiuAA", "Other"),
         `Isolate ID` = factor(`Isolate ID`, levels = iso_order)) %>%
  ggplot(aes(x = `Isolate ID`, y = `IAA (µg mL-1)`)) +
  geom_boxplot(aes(color = Group, fill = Group), outlier.shape = NA) +
  geom_point(aes(color = Group), position = position_jitter(width = 0.15), size = 1) +
  scale_color_manual(values = c("#C77B2C", "#1C3D5A")) +
  scale_fill_manual(values = c("#EAC79E", "#C7D5E0")) +
  annotate("rect", xmin = n_other + 0.5, xmax = length(iso_order) + 0.5, ymin = -Inf, ymax = Inf, alpha = 0.12, fill = "#C77B2C") +
  labs(x = "ZeaMiC Isolate", y = expression(italic("In vitro")~IAA~production~(mu*g~mL^{-1}))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), legend.position = "none", 
        axis.text.y  = element_text(size = 10))

ggsave("Figure panels/Figure8A.svg", IAA_PLOT, width = 26, height = 10, unit = "cm")