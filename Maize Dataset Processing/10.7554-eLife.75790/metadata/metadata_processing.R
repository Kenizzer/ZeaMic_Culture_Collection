# Metadata concatenation for 10.7554-eLife.75790

# Association analyses of host genetics, root-colonizing microbes, and plant 
# phenotypes under different nitrogen conditions in maize

library(tidyverse)


# load in metadata files
SRP298238 <- read.csv("SRP298238_meta.csv")
SRP298247 <- read.csv("SRP298247_meta.csv")
SRP341638 <- read.csv("SRP341638_meta.csv")
SRP341639 <- read.csv("SRP341639_meta.csv")
SRP341643 <- read.csv("SRP341643_meta.csv")

# identify columns to keep
colnames(SRP298238)
colnames(SRP298247)
colnames(SRP341638)
colnames(SRP341639)
colnames(SRP341643)

# columns to keep
SRP298238 <- SRP298238 %>%
  select(Run, Block, genotype, geo_loc_name_country_continent, geo_loc_name,
       isolation_source, Pedigree, raw_seq_count, Sample.Name,
       split_plot_block, split_plot, Subpopulation,Subsample, year)

SRP298247 <- SRP298247 %>%
  select(Run, Block, genotype, geo_loc_name_country_continent, geo_loc_name,
         isolation_source, Pedigree, raw_seq_count, Sample.Name,
         Subpopulation, Subsample, year)

SRP341638 <- SRP341638 %>%
  dplyr::rename(genotype = host_genotype) %>%
  select(Run, Block, genotype, geo_loc_name_country_continent, geo_loc_name,
         isolation_source, Pedigree, Sample.Name,
         split_plot_block, split_plot, Subpopulation,Subsample, year)

SRP341639 <- SRP341639 %>%
  dplyr::rename(genotype = host_genotype) %>%
  select(Run, Block, genotype, geo_loc_name_country_continent, geo_loc_name,
         isolation_source, Pedigree, Sample.Name,
         split_plot_block, split_plot, Subpopulation,Subsample, year)

SRP341643 <- SRP341643 %>%
  dplyr::rename(genotype = host_genotype) %>%
  select(Run, Block, genotype, geo_loc_name_country_continent, geo_loc_name,
         isolation_source, Pedigree, Sample.Name,
         split_plot_block, split_plot, Subpopulation,Subsample, year)

ALL_meta <- bind_rows(SRP298238, SRP298247, SRP341638, SRP341639, SRP341643)

# fix column names
ALL_meta <- janitor::clean_names(ALL_meta)
ALL_meta <- ALL_meta %>% dplyr::rename(SampleID = run)

# phyloseq specific formatting
rownames(ALL_meta) <- ALL_meta$SampleID # make the sample IDs the row names for import into Phyloseq

# save copies
write.csv(ALL_meta, "10.7554-eLife.75790_metadata.csv")
saveRDS(ALL_meta, "10.7554-eLife.75790_metadata.RDS")