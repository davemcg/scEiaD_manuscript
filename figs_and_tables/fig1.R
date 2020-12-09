library(tidyverse)
library(tidygraph)
library(ggdag)
library(ggraph)
library(cowplot)
library(patchwork)
a <- meta_filter %>% 
  mutate(Platform = gsub('C1', 'SMARTSeq_v4', Platform)) %>%
  group_by(organism, Platform) %>% 
  summarise(accession = length(unique(study_accession)), paper = length(unique(PMID)), batch = length(unique(batch))) %>% 
  pivot_longer(cols = accession:batch) %>% 
  mutate(Count = value,
         Type = name,
         organism = gsub(' ','\n',organism)) %>%  
  ggplot(aes(x=Platform, y=Count, fill = Type, group = Count)) + 
  facet_wrap(~organism) + 
  geom_bar(stat = 'identity', position = position_dodge2()) + 
  cowplot::theme_cowplot() + 
  coord_flip()

dagified <- dagify('counts' ~ 'fastq',
                   'SnakePOP' ~ 'CellType',
                   'SnakePOP' ~ 'counts',
                   'LatentDims' ~ 'SnakePOP',
                   'cluster' ~ 'LatentDims',
                   'umap' ~ 'LatentDims' ,
                   'ML' ~ 'CellType',
                   'projectedCellType' ~ 'ML',
                   'ML' ~ 'LatentDims',
                   'DiffTesting' ~ 'counts',
                   'DiffTesting' ~ 'projectedCellType',
                   'DiffTesting' ~ 'CellType',
                   'DiffTesting' ~ 'cluster',
                   'DiffTesting' ~ 'LatentDims',
                   'projectedCellType' ~ 'LatentDims',
                   'projectedCellType' ~ 'CellType')

tidy_dagitty(dagified)
set.seed(2534)
b <- as_tbl_graph(dagified) %>% mutate(bop = name) %>% ggraph(layout = 'graphopt') + 
  geom_edge_link(arrow = arrow(length = unit(3, 'mm')), 
                 end_cap = circle(8, 'mm'), start_cap = circle(2, 'mm')) + 
  #geom_node_point(size = 10) +
  geom_node_label(aes(label = name), color = 'black', alpha = 1)  + 
    theme_nothing()

  

# CellType
ct_order <- qc %>% mutate(CellType = case_when(is.na(CellType) ~ 'Unlabelled', grepl('RPC', CellType) ~ 'RPCs', TRUE ~ CellType)) %>% group_by(CellType) %>% count() %>% arrange(-n) %>% filter(n > 10000 | CellType %in% c('Astrocytes','RPE', 'Horizontal Cells'), CellType != 'Unlabelled') %>% pull(CellType)

c <- qc %>% 
  mutate(CellType = case_when(is.na(CellType) ~ 'Unlabelled', grepl('RPC', CellType) ~ 'RPCs', TRUE ~ CellType)) %>% 
  filter(CellType %in% ct_order) %>% 
  filter(!is.na(CellType), 
         !is.na(study_accession), 
         !CellType %in% c('Doublet', 'Doublets'),
         !grepl('RPE|Vascul', CellType)) %>% 
  mutate(organism = case_when(grepl('Homo', organism) ~ 'HS',
                              grepl('Mus', organism) ~ 'MM',
                              TRUE ~ 'MF')) %>% 
  group_by(CellType, organism, study_accession) %>% 
  summarise(Count = n()) %>% 
  summarise(Studies = length(study_accession), Count = sum(Count)) %>%   
  mutate(CellType = factor(CellType, levels = ct_order)) %>% 
  ggplot(aes(label = Studies, x=CellType, y = Count, fill = organism,)) + 
  geom_bar(stat = 'identity') + 
  geom_text(position = position_stack(vjust = 0.5)) + 
  coord_flip() + 
  theme_cowplot()

# (a / c) | b
