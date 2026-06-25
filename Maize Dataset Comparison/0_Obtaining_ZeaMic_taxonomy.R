## Code by: Joel Swift 2026
# ZeaMic isolate taxonomy

library(tidyverse)

## All Isolate names 
all_isolates_df <-  data.frame(ZeaMic_Isolate = c("HAY10-12D","HAY10-2E","HAY10-8C","HAY1-11D","HAY1-12B","HAY1-1A",
                                                  "HAY1-2A","HAY1-2D","HAY1-4A","HAY1-4F","HAY1-5C","HAY1-7F","HAY1-8B",
                                                  "HAY1-8D","HAY1-8F","KNZ12-10D","KNZ12-1A","KNZ12-1H","KNZ12-8B",
                                                  "KNZ13-11D","KNZ13-1D","KNZ14-4A","KNZ14-4F","KNZ14-9A","KNZ1-6G",
                                                  "KNZ2-6A","KU2C","KU3C","KU4A","KU4C","KU4D","KU4F","MI12","MI13",
                                                  "MI2","MI27","MI9","MSU2B","MSU2C","MSU3D","MSU4B","MSU4C","MSUFA",
                                                  "MSUFE","NiuAA1","NiuAA2","NiuAA3","NiuAA4","NiuAA5","NiuAA6","NiuAA7",
                                                  "SB1","SD1074B","SD173","SD811","SVR3-10G","SVR4-12A","SVR4-3F","SVR4-6G",
                                                  "SVR4-7G","SVR5-10H","SVR5-12B","SVR5-2C","SVR5-6F","SVR5-6H","SVR5-9E",
                                                  "T428B","T572","TLI6-9E","TLI7-10B","TLI7-1G","TLI7-4A","TLI7-4D",
                                                  "TLI7-8E","TLI7-8G","TLI8-12F","TLI8-1H","TLI8-3E","TLI8-6A","TLI9-4D",
                                                  "TLI9-5A","LBA21","LBA71","LDE1","LMG1","LML1","LRC7.O","LSP13",
                                                  "E-01", "R-06", "R-17", "R-28", "R-29", "R-33"), 
                               genus = c("Enterobacter","Enterobacter","Enterobacter","Monashia","Burkholderia",
                                         "Burkholderia","Monashia","Trinickia","Luteibacter","Leifsonia","Monashia",
                                         "Burkholderia","Mesorhizobium","Burkholderia","Labrys","Microbacterium",
                                         "Luteibacter","Paraburkholderia","Leifsonia","Leifsonia","Burkholderia",
                                         "Paraburkholderia","Burkholderia","Rhizobium","Paraburkholderia","Rhizobium",
                                         "Pseudomonas","Pseudomonas","Chryseobacterium","Cellulomonas","Chryseobacterium",
                                         "Microbacterium","Telluria","Telluria","Rahnella","Enterobacter","Arthrobacter",
                                         "Flavobacterium","Brevibacillus","Mitsuaria","Rhodococcus","Stenotrophomonas",
                                         "Pseudomonas","Pseudomonas","Stenotrophomonas","Brucella","Curtobacterium",
                                         "Enterobacter","Chryseobacterium","Herbaspirillum","Pseudomonas","Methylobacterium",
                                         "Arthrobacter","Franconibacter","Pantoea","Kocuria","Chitinophaga","Achromobacter",
                                         "Achromobacter","Luteibacter","Bradyrhizobium","Achromobacter","Rhizobium","Enterobacter",
                                         "Paraburkholderia","Rhizobium","Sinomonas","Arthrobacter","Agrobacterium","Rhizobium",
                                         "Dyella","Luteibacter","Burkholderia","Variovorax","Burkholderia","Rhizobium",
                                         "Burkholderia","Dyella","Burkholderia","Rhizobium","Variovorax","Peribacillus",
                                         "Bacillus","Deinococcus","Streptomyces","Nocardioides","Rhizobium","Sphingobium",
                                         "Dyadobacter", "Pedobacter", "Exiguobacterium", "Acidovorax", "Ralstonia", "Mucilaginibacter"))

# Get full taxonomy for isolates
## Getting taxonomy for genera
library(taxize)
# Query NCBI taxonomy
tax_info <- tax_name(query = unique(all_isolates_df$genus), get = c("phylum", "class", "order", "family"), db = "ncbi")

# Convert to dataframe
tax_info_df <- as.data.frame(tax_info)
tax_info_df <- tax_info_df %>% rename("genus"= "query")

##merge back to isolates
Fulltax_tab <- all_isolates_df %>% left_join(tax_info_df, by = "genus") %>% select(!db)
Fulltax_tab$X <- NULL # remove column 

write.csv(Fulltax_tab, "ZeaMic Isolates with Full Taxonomy.csv")