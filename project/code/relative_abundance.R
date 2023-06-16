# Analysis of microbial abundance
# Script written by Johan S. Saenz
# Data: 13-06-23

# Set-WD -------------------------------------------------------
setwd("/Users/sebastiansaenz/Documents/github/MicroTutorials/project")

# Install-libraries ------------------------------------------------------
# install.packages("Tidyverse")
#install.packages("patchwork") 

# Load packages and files ---------------------------------------------
library(tidyverse)
library(patchwork)

# color-palette ---------------------------------------------

simple_color <- c("#d8b365", "#5ab4ac")

# load df ------------------- 

df_taxa <- read_tsv("rawdata/taxa.tsv")

metadata <- read.table(file = "rawdata/metadata_metapro.txt",
                       header = TRUE,
                       sep = "\t",
                       check.names = FALSE) %>% 
  rename(rawfile=Probenname, section=Abschnitt) %>% 
  mutate(rawfile = str_replace(rawfile, "\\.", "_")) %>% 
  separate(section,
           into = c("section", "origin"),
           sep = "\\s") 

# relative abundance ----------------------------

df_rel_abund <- df_taxa %>% 
  select(Phylum, contains("Intensity")) %>% 
  pivot_longer(-Phylum, names_to = "rawfile", values_to = "intensity") %>% 
  mutate(rawfile = str_remove(rawfile, "Intensity "),
         Phylum = if_else(is.na(Phylum), "Unclassified", Phylum)) %>% 
  group_by(Phylum, rawfile) %>% 
  summarise(sum_intensity = sum(intensity), .groups = "drop") %>% 
  group_by(rawfile) %>% 
  mutate(rel_abundance = 100*(sum_intensity/sum(sum_intensity)))

pooled_treshold <- 1

df_pool <- df_rel_abund %>% 
  group_by(Phylum) %>% 
  summarise(pool = mean(rel_abundance) < pooled_treshold )
bar_col <-  c('#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598','#abdda4','#66c2a5','#3288bd')

df_rel_abund %>%
  inner_join(df_pool, by ="Phylum") %>% 
  mutate(Phylum_new = if_else(pool == TRUE, "Others", Phylum)) %>%
  group_by(Phylum_new, rawfile) %>% 
  summarise(sum_rel = sum(rel_abundance), .groups = "drop") %>%
  inner_join(metadata, by = "rawfile") %>% 
  group_by(Phylum_new, section) %>% 
  summarise(mean_rel = mean(sum_rel), .groups = "drop") %>% 
  ggplot(aes(x = section,
             y = mean_rel,
             fill = Phylum_new)) +
  geom_col() +
  scale_y_continuous(limits = c(0, 100),
                     expand = c(0, 0)) +
  scale_fill_manual(values = bar_col) +
  labs(y = "Relative abundance %",
       x =NULL) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(size =15))


  
#statistical test
df_rel_abund %>% 
  inner_join(metadata, by = "rawfile") %>% 
  filter(Phylum == "Actinobacteriota") %>% 
  kruskal.test(rel_abundance~section, data =.)

df_rel_abund %>% 
  inner_join(metadata, by = "rawfile") %>% 
  filter(Phylum == "Actinobacteriota") %>% 
  ggplot(aes(x = section,
             y = rel_abundance,
             color = section)) +
  geom_boxplot(width = 0.5, 
               color = "black", 
               outlier.colour = "white") +
  geom_jitter(width = 0.3) +
  scale_color_manual(values = simple_color) +
  scale_y_continuous(limits = c(0, 7),
                     breaks = seq(0, 7, 1)) +
  labs(y = "Relative abundance (%)",
       title = "Kruskal-Wallis chi-squared = 9.1429, df = 1, p-value = 0.002497") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(size = 8))

