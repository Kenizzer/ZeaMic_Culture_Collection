# Metadata processing for https://doi.org/10.1038/s41477-024-01654-7

# Heritable microbiome variation is correlated with source 
# environment in locally adapted maize varieties

library(tidyverse)


# load in metadata files
meta.tab <- read.csv("He2024_metadata_table1.csv")

# identify columns to keep
colnames(meta.tab)

# columns to keep
meta.tab <- meta.tab %>%
  select(Run, Sample.Name, Cultivar, geo_loc_name_country, geo_loc_name_country_continent,
         Organism, tissue, treatment)

# fix column names
meta.tab <- janitor::clean_names(meta.tab)
meta.tab <- meta.tab %>% dplyr::rename(SampleID = run)

# phyloseq specific formatting
rownames(meta.tab) <- meta.tab$SampleID # make the sample IDs the row names for import into Phyloseq

# save copies
write.csv(meta.tab, "10.1038-s41477-024-01654-7_metadata.csv")
saveRDS(meta.tab, "10.1038-s41477-024-01654-7_metadata.RDS")
