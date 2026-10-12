# Code by Rob McMinn/Joel Swift 2026
library(tidyverse)
library(readxl)
library(patchwork)
theme_set(theme_bw())

##### In Vivo DATA #####
# =============================================================================
# Root surface area  —  per-trial standardization (Hedge's g vs control) + figure
#
#   Treatment is "Control" or "Inoculated"; each trial has its own controls.
#
# Pipeline:  read raw plant values  ->  standardize each inoculated plant to its
#            trial's control (Hedge's g)  ->  box-and-whisker per isolate, with
#            the NiuAA isolates pinned right.
# =============================================================================

raw <- read_excel("./InVivo_RSA.xlsx") %>% rename(value = `Surface.Area.mm2`)
# Per-trial control statistics 
ctrl <- raw %>%
  filter(Treatment == "Control") %>%
  group_by(Trial) %>%
  summarise(m_c = mean(value), v_c = var(value), n_c = n(), .groups = "drop")
# Per-replicate Hedge's g, standardized within each trial 
Vivo_df <- raw %>%
  filter(Treatment == "Inoculated") %>%
  group_by(Isolate, Trial) %>%
  mutate(n_i = n(), v_i = var(value)) %>%
  ungroup() %>%
  left_join(ctrl, by = "Trial") %>%
  mutate(
    pooled_sd = sqrt(((n_i - 1) * v_i + (n_c - 1) * v_c) / (n_i + n_c - 2)),
    J         = 1 - 3 / (4 * (n_i + n_c) - 9),     # remove "* J" below for plain Cohen's d
    g         = (value - m_c) / pooled_sd * J
  ) %>%
  select(Isolate, Trial, g)

## Order isolates by mean IAA, but force NiuAA isolates to the end
iso_order_vivo <- Vivo_df %>%
  group_by(Isolate) %>%
  summarize(mean_g = mean(g), .groups = "drop") %>%
  arrange(desc(mean_g)) %>%
  pull(Isolate) %>%
  {\(x) c(x[!str_detect(x, "^NiuAA")], x[str_detect(x, "^NiuAA")])}()

###### In Vitro Data ######
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
  {\(x) c(x[!str_detect(x, "^NiuAA")], x[str_detect(x, "^NiuAA")])}()
## Number of non-NiuAA isolates (used for shading)
n_other <- sum(!str_detect(iso_order, "^NiuAA"))


### Add in genus names
Zeamic_taxonomy <- read.csv("ZeaMic Isolates with Full Taxonomy.csv")
Zeamic_taxonomy$X <- NULL
Zeamic_taxonomy <- Zeamic_taxonomy %>% rename(Isolate = ZeaMic_Isolate)
IAA_data <- IAA_data %>% rename(Isolate = `Isolate ID`)

# add taxonomy
IAA_data <- left_join(IAA_data, Zeamic_taxonomy, by = "Isolate")
Vivo_df <- left_join(Vivo_df, Zeamic_taxonomy, by = "Isolate")

##### Plots #####
## Very minimal class barplot, will be attached to other plots
class_bar_iso <- Zeamic_taxonomy %>% 
  mutate(Isolate = factor(Isolate, levels = iso_order_vivo)) %>%
  ggplot(., aes(x = 1, y = Isolate, fill = class)) +
  geom_tile() +
  scale_fill_manual(values = c(
    "Actinomycetes" = "#23733a",
    "Bacilli" = "#b8d437",
    "Deinococci" = "#ffc323",
    "Chitinophagia" = "#f0852f",
    "Cytophagia" = "#c68642",
    "Flavobacteriia" = "#e51f41",
    "Sphingobacteriia" = "#a63d40",
    "Alphaproteobacteria" = "#c6529f",
    "Betaproteobacteria" = "#405caa",
    "Gammaproteobacteria" = "#46bbe7"), name = "Class") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid = element_blank(),
        panel.border = element_blank(), 
        plot.margin = margin(t = 5, r = 2, b = 5, l = 5), 
        legend.position = "right")

## invivo
Invivo_PLOT <- Vivo_df %>%
  mutate(Group = if_else(str_detect(Isolate, "^NiuAA"), "NiuAA", "Other"),
         Isolate = factor(Isolate, levels = iso_order_vivo)) %>%
  ggplot(., aes(y = Isolate, x = g)) +
  geom_boxplot(aes(color = Group, fill = Group), outlier.shape = NA) +
  geom_point(aes(color = Group), size = 1) +
  scale_color_manual(values = c("#C77B2C", "#1C3D5A")) +
  scale_fill_manual(values = c("#EAC79E", "#C7D5E0")) +
  annotate("rect", ymin = n_other + 0.5, ymax = length(iso_order) + 0.5, xmin = -Inf, xmax = Inf, alpha = 0.12, fill = "#C77B2C") +
  labs(y = "ZeaMiC Isolate", x = "Standardized Root Surface Area") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "#555555", linewidth = 0.5) +
  theme(axis.text.y = element_text(size = 7, face = "plain"))

## in vitro
IAA_PLOT <- IAA_data %>%
  mutate(Group = if_else(str_detect(Isolate, "^NiuAA"), "NiuAA", "Other"),
         Isolate = factor(Isolate, levels = iso_order_vivo)) %>%
  ggplot(aes(y = Isolate, x = `IAA (µg mL-1)`)) +
  geom_boxplot(aes(color = Group, fill = Group), outlier.shape = NA) +
  geom_point(aes(color = Group), position = position_jitter(width = 0.15), size = 1) +
  scale_color_manual(values = c("#C77B2C", "#1C3D5A")) +
  scale_fill_manual(values = c("#EAC79E", "#C7D5E0")) +
  annotate("rect", ymin = n_other + 0.5, ymax = length(iso_order) + 0.5, xmin = -Inf, xmax = Inf, alpha = 0.12, fill = "#C77B2C") +
  labs(y = "ZeaMiC Isolate", x = expression(italic("In vitro")~IAA~production~(mu*g~mL^{-1}))) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "#555555", linewidth = 0.5) +
  theme(axis.text.y = element_text(size = 7, face = "plain"))

##### Patchwork combination #####
Figure7 <- (class_bar_iso + Invivo_PLOT + IAA_PLOT) +
  plot_layout(widths = c(0.2, 5, 5), guides = "collect", axes = "collect") & 
  theme(legend.position = "none")

ggsave("Figure7.svg", Figure7, width = 18.5, height = 24, unit = "cm")
ggsave("Figure7.png", Figure7, width = 18.5, height = 24, unit = "cm")
