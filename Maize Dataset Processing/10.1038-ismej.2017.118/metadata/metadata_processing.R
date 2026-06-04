# Metadata concatenation for 10.1038/ismej.2017.118

# Drought and host selection influence bacterial community 
# dynamics in the grass root microbiome

library(tidyverse)

# load in metadata files
Naylor_meta <- read.csv("Naylor_metadata_raw.csv")

# identify columns to keep
colnames(Naylor_meta)

# columns to keep
Naylor_meta <- Naylor_meta %>%
  select(Run, Sample.Name, geo_loc_name_country_continent, plant_body_site, watering_regm, isolation_source, 
         host_life_stage, SRA.Study, Collection_Date)

# fix column names
clean_meta <- janitor::clean_names(Naylor_meta)
clean_meta <- clean_meta %>% dplyr::rename(SampleID = run)

# phyloseq specific formatting
rownames(clean_meta) <- clean_meta$SampleID # make the sample IDs the row names for import into Phyloseq

# save copies
write.csv(clean_meta, "10.1038-ismej.2017.118_metadata.csv")
saveRDS(clean_meta, "10.1038-ismej.2017.118_metadata.RDS")
