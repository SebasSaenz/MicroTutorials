#diversity

df_taxa <- read_tsv("rawdata/taxa.tsv")

df_taxa %>% 
  select(Species, contains("Intensity")) %>% 
  pivot_longer(-Species, names_to = "rawfile", values_to = "intensity") %>% 
  mutate(rawfile = str_remove(rawfile, "Intensity "),
         Species = if_else(is.na(Species), "Unclassified", Species)) %>%
  group_by(Species, rawfile) %>% 
  summarise(sum_intensity = sum(intensity), .groups = "drop") %>% 
  inner_join(metadata, by = "rawfile") %>%
  mutate(presence = if_else(sum_intensity >0, 1, 0)) %>% 
  group_by(Species, section) %>% 
  summarise(sum_diversity = sum(presence)) %>% 
  mutate(diversity = if_else(sum_diversity >5, 1,0)) %>% 
  group_by(section) %>% 
  summarise(sum_diversity = sum(diversity))
  


