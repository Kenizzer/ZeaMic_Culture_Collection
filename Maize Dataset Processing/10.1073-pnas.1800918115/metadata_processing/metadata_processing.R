# Metadata concatenation for 10.7554-eLife.75790

# Association analyses of host genetics, root-colonizing microbes, and plant 
# phenotypes under different nitrogen conditions in maize

library(tidyverse)

# load in metadata files
Walters1 <- read.csv("Walters2018_metadata_table1.csv")
Walters2 <- read.csv("Walters2018_metadata_table2.csv")

# identify columns to keep
colnames(Walters1)
colnames(Walters2)

head(Walters1)
head(Walters2)

# columns to keep
Walters1 <- Walters1 %>%
  select(Run, sample_name, geo_loc_name_country)

Walters2 <- Walters2 %>%
  select(Run, geo_loc_name_country , host_common_name, env, inbreds, sample_type, AGE, age_unit)

ALL_meta <- bind_rows(Walters1, Walters2)

# fix column names
ALL_meta <- janitor::clean_names(ALL_meta)
ALL_meta <- ALL_meta %>% dplyr::rename(SampleID = run)

# phyloseq specific formatting
rownames(ALL_meta) <- ALL_meta$SampleID # make the sample IDs the row names for import into Phyloseq

# Recode compartment names to be useful
# soil -> rhizosphere 
# bulk soil -> soil
ALL_meta <- ALL_meta %>% 
  mutate(sample_type = recode(sample_type,
                              "soil" = "rhizosphere",
                              "bulk soil" = "soil", 
                              .missing = "rhizosphere"))

# save copies
write.csv(ALL_meta, "10.1073-pnas.1800918115_metadata.csv")
saveRDS(ALL_meta, "10.1073-pnas.1800918115_metadata.RDS")
