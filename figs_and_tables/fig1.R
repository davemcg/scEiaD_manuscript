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
                   'LatentDims' ~ 'scVIrefquery',
                   'cluster' ~ 'LatentDims',
                   'umap' ~ 'LatentDims' ,
                   'ML' ~ 'CellType',
                   'projectedCellType' ~ 'ML',
                   'ML' ~ 'LatentDims',
                   'scVIrefquery' ~ 'SnakePOP',
                   'DiffTesting' ~ 'counts',
                   'DiffTesting' ~ 'projectedCellType',
                   'DiffTesting' ~ 'CellType',
                   'DiffTesting' ~ 'cluster',
                   'SnakePOP' ~ 'BatchInfo',
                   'scVIrefquery' ~ 'BatchInfo',
                   'projectedCellType' ~ 'LatentDims')

tidy_dagitty(dagified)
set.seed(51345)
b <- tidy_dagitty(dagified) %>% as_tbl_graph() %>% 
  mutate(type = case_when(name %in% c('fastq', 'CellType', 'BatchInfo') ~ 'Input',
                          name %in% c('DiffTesting', 'Velocity', 'Trajectory','umap', 'projectedCellType','cluster') ~ 'Outputs'),
         name = case_when(name == 'CellType' ~ 'Published\nCell Types',
                          name == 'projectedCellType' ~ 'Learned\nCell Types',
                          name == 'DiffTesting' ~ 'Diff\nTesting',
                          name == 'scVIrefquery' ~ 'scVI Ref +\nQuery',
                          name == 'LatentDims' ~ 'Latent Dims',
                          name == 'BatchInfo' ~ 'Batch\nInfo',
                          TRUE ~ name)) %>% 
  ggraph(layout = 'gem') + 
  geom_edge_link(arrow = arrow(length = unit(3, 'mm')), 
                 end_cap = circle(9, 'mm'), start_cap = circle(2, 'mm')) + 
  #geom_node_point(size = 10) +
  geom_node_label(aes(label = name, color = type), alpha = 1)  + 
  scale_color_manual(values = c('Red','Blue'), na.value = 'Black') +
  theme_nothing()



# CellType
ct_order <- qc %>% mutate(CellType = case_when(is.na(CellType) ~ 'Unlabelled', grepl('RPC', CellType) ~ 'RPCs', TRUE ~ CellType)) %>% group_by(CellType) %>% dplyr::count() %>% arrange(-n) %>% filter(n > 10000 | CellType %in% c('Astrocytes','RPE', 'Horizontal Cells'), CellType != 'Unlabelled') %>% pull(CellType)

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
