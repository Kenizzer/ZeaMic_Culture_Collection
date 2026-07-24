# =============================================================================
# Root surface area  —  per-trial standardization (Hedge's g vs control) + figure
#
#   Treatment is "Control" or "Inoculated"; each trial has its own controls.
#
# Pipeline:  read raw plant values  ->  standardize each inoculated plant to its
#            trial's control (Hedge's g)  ->  box-and-whisker per isolate, with
#            the NiuAA isolates pinned right.
# =============================================================================

# Code by Rob McMinn 2026

library(readxl)
library(dplyr)
library(ggplot2)

raw <- read_excel("./InVivo_RSA.xlsx") %>% rename(value = `Surface.Area.mm2`)

# Per-trial control statistics 
ctrl <- raw %>%
  filter(Treatment == "Control") %>%
  group_by(Trial) %>%
  summarise(m_c = mean(value), v_c = var(value), n_c = n(), .groups = "drop")

# Per-replicate Hedge's g, standardized within each trial 
g_df <- raw %>%
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
 
# FIGURE
df <- g_df  # Isolate + per-replicate g

BOX_FILL <- "#C7D5E0"; BOX_EDGE <- "#1C3D5A"; PT     <- "#1C3D5A"   # main isolates
NIU_FILL <- "#EAC79E"; NIU_EDGE <- "#C77B2C"; NIU_PT <- "#C77B2C"   # Niu isolates

# Niu detector
is_niu <- function(x) grepl("^Niu", x)
df$group <- ifelse(is_niu(df$Isolate), "Niu", "main")

ord <- df %>%
  group_by(Isolate, group) %>%
  summarise(med = median(g), .groups = "drop") %>%
  arrange(group == "Niu", desc(med))
lev        <- ord$Isolate
df$Isolate <- factor(df$Isolate, levels = lev)

n_main  <- sum(ord$group == "main")
n_all   <- length(lev)
divider <- n_main + 0.5

# identity columns by color
df$box_fill <- ifelse(df$group == "Niu", NIU_FILL, BOX_FILL)
df$box_edge <- ifelse(df$group == "Niu", NIU_EDGE, BOX_EDGE)
df$pt_col   <- ifelse(df$group == "Niu", NIU_PT,   PT)

lab_col  <- ifelse(is_niu(lev), NIU_PT, BOX_EDGE)
lab_face <- ifelse(is_niu(lev), "bold", "plain")



df %>%
  group_by(Isolate) %>%
  summarize(n = n(), .groups = "drop") %>%
  summarize(
    min_n = min(n),
    max_n = max(n))



set.seed(7)
p <- ggplot(df, aes(x = Isolate, y = g)) +
  annotate("rect", xmin = divider, xmax = n_all + 0.5, ymin = -Inf, ymax = Inf,
           fill = NIU_EDGE, alpha = 0.12) +                                    # Niu band
  geom_hline(yintercept = 0, linetype = "dotted", colour = "#555555",
             linewidth = 0.5) +                                               # zero line
  geom_vline(xintercept = divider, linetype = "dashed", colour = NIU_EDGE,
             linewidth = 0.5) +                                               # divider
  geom_boxplot(aes(group = Isolate, fill = box_fill, colour = box_edge),
               width = 0.62, outlier.shape = NA, linewidth = 0.4, fatten = 2) +
  geom_jitter(aes(colour = pt_col), width = 0.05, height = 0,
              size = 1.1, alpha = 0.7, stroke = 0.2) +
  scale_fill_identity() +
  scale_colour_identity() +
  labs(x = "ZeaMiC Isolate", y = "Standardized Root Surface Area") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1,
                                size = 6, colour = lab_col, face = lab_face),
    axis.text.y  = element_text(size = 10))

# Save 
ggsave("Figure panels/Figure7B.svg", p, width = 26, height = 10, unit = "cm")

