# Analysis of metaproteomics run
# Script written by Johan S. Saenz
# Data: 12-06-23

# Set-WD -------------------------------------------------------
setwd("/Users/sebastiansaenz/Documents/github/MicroTutorials/project")

# Install-libraries ------------------------------------------------------
# install.packages("Tidyverse")
install.packages("patchwork") 

# load packages and files ---------------------------------------------
library(tidyverse)
library(patchwork)

# color-palette ---------------------------------------------

simple_color <- c("#d8b365", "#5ab4ac")
# --------------
df_summary <-  read.table(file = "rawdata/final_summary.tsv",
           header = TRUE,
           sep = "\t",
           check.names = FALSE) %>% 
  rename(rawfile = 'Raw file',
         experiment = "Experiment",
         msms = 'MS/MS',
         msms_identified = 'MS/MS Identified',
        msms_identifed_perc = 'MS/MS Identified [%]',
        peptide_sequences = "Peptide Sequences Identified") %>% 
  select(rawfile, msms_identifed_perc, peptide_sequences) %>% 
  filter(rawfile != "Total")

metadata <- read.table(file = "rawdata/metadata_metapro.txt",
                       header = TRUE,
                       sep = "\t",
                       check.names = FALSE) %>% 
  rename(rawfile=Probenname, section=Abschnitt) %>% 
  mutate(rawfile = str_replace(rawfile, "\\.", "_")) %>% 
  separate(section,
           into = c("section", "origin"),
           sep = "\\s") 
# Create a summary table


sapply(df_summary[3:6], min)
sapply(df_summary[2:3], max)
sapply(df_summary[2:3], median)
sapply(df_summary[2:3], mean)

mytable <- data.frame(Stats = c("Min", "Max", "Median", "Mean"),
           msms_identifed_perc =c(1.96, 58.35, 30.295, 26.159),
           peptide_sequences = c(642.00,27841.00, 10622.000, 10646.400))

write_tsv(mytable, file = "rawdata/my_summarytable.tsv")

#########################

peptides_plot <- df_summary %>% # peptides plot
  pivot_longer(-rawfile, names_to = "qc", values_to = "values") %>%
  inner_join(metadata, by = "rawfile") %>% 
  filter(qc == "peptide_sequences") %>% 
  ggplot(aes(x = section,
             y = values,
             color = section)) +
  geom_boxplot(width = 0.2) +
  geom_jitter(size = 3, shape = 8) +
  labs(y = "Number peptide sequences",
       x = NULL) +
  scale_y_continuous(limits = c(0, 30000), breaks = seq(0, 30000 ,3000)) +
  scale_color_manual(values = simple_color) +
  theme_classic() +
  theme(axis.title.y = element_text(size=10),
        legend.position = "nonen") 


ms_plot <- df_summary %>% # msms plot
  pivot_longer(-rawfile, names_to = "qc", values_to = "values") %>%
  inner_join(metadata, by = "rawfile") %>% 
  filter(qc == "msms_identifed_perc") %>% 
  ggplot(aes(x = section,
             y = values,
             color = section)) +
  geom_boxplot(width = 0.2) +
  geom_jitter(size = 3, shape = 8) +
  labs(y = "MS/MS identified %",
       x = NULL) +
  scale_y_continuous(limits = c(0, 70), breaks = seq(0, 70 ,10)) +
  scale_color_manual(values = simple_color) +
  theme_classic() +
  theme(axis.title.y = element_text(size=10),
        legend.position = "nonen") 


Fig_1 <- peptides_plot + ms_plot + plot_annotation(tag_levels = 'A')

ggsave(Fig_1, file = "figures/Fig_1.png", height = 5, width = 10, dpi = 400)

# statistical-test ------------------------------------------------

df_summary %>% 
  pivot_longer(-rawfile, names_to = "qc", values_to = "values") %>%
  inner_join(metadata, by = "rawfile") %>% 
  filter(qc == "msms_identifed_perc") %>% 
  kruskal.test(values~section, data = .)



