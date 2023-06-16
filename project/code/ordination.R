setwd("/Users/sebastiansaenz/Documents/github/MicroTutorials/project")

# libraries ----------------------
library(tidyverse)
library(vegan)

# data------------------
df_proteins <- read_tsv("rawdata/final_proteins.tsv")

metadata <- read.table(file = "rawdata/metadata_metapro.txt",
                       header = TRUE,
                       sep = "\t",
                       check.names = FALSE) %>% 
  rename(rawfile=Probenname, section=Abschnitt) %>% 
  mutate(rawfile = str_replace(rawfile, "\\.", "_")) %>% 
  separate(section,
           into = c("section", "origin"),
           sep = "\\s") 

# Calculate relative abundance ---------------------

rel_abundance_df <- df_proteins %>%
  rename(proteinid = `Protein IDs`) %>% 
  select(proteinid, contains("Intensity"), -Intensity) %>% 
  pivot_longer(-proteinid, names_to = "rawfile", values_to = "intensity") %>%
  mutate(intensity = as.numeric(intensity),
         rawfile = str_remove(rawfile, "Intensity ")) %>% 
  group_by(rawfile) %>% 
  mutate(rel_abundance = 100*(intensity/sum(intensity))) %>%
  select(proteinid, rawfile, rel_abundance) %>% 
  pivot_wider(names_from = rawfile, values_from = rel_abundance)

# Create a matrix and transpose data--------------
matrix <- rel_abundance_df[2:21] %>% 
  t()

#Calculate distance -----------------
nmds1 <- metaMDS(matrix,  #perform nmds
                 distance = "bray", 
                 try = 20, 
                 trymax = 100, 
                 maxit = 1000, 
                 k = 3)


# qc checking
nmds1$stress  
stressplot(nmds1)


#Extracting scores
data_scores <- as.data.frame(scores(nmds1, display=c("sites")))

#Addd metadata to dataframe
data_scores$rawfile <- as.character(row.names(data_scores)) 

#joing metadata with nmds scorres
data_nmds <- left_join(data_scores, metadata, by = "rawfile")

data_nmds$section <- factor(data_nmds$section,
                             levels = c("Caecum", "Jejunum"))
#make plot

data_nmds %>%
  ggplot() +
  geom_point(aes(x = NMDS1,
                 y = NMDS2,
                 colour =section),
             size = 4
             #alpha = 0.7
  ) +
  scale_color_manual(values = c('#66C2A5','#FC8D62')) +
  theme_classic() +
  theme(panel.background = element_blank(), #remove background
        panel.grid.major = element_blank(), #remove grid
        panel.grid.minor = element_blank(),#remove grid
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        legend.position = "top") #remove legend title

