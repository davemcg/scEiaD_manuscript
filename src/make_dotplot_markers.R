library(pool)
library(RSQLite)
library(tidyverse)
library(glue)
scEiaD_2020_v01 <- dbPool(drv = SQLite(), dbname = "~/data/scEiaD/2021_02_05_MOARTABLES__anthology_limmaFALSE___5000-transform-counts-universe-batch-scVIprojectionSO-8-0.1-50-0.6.sqlite", idleTimeout = 3600000)

meta_filter <- fst::read_fst('~/git/plaeApp/inst/app/www/meta_filter.fst') %>% as_tibble()
meta_filter <- meta_filter %>% mutate(CellType_predict = case_when(CellType_predict == 'Photoreceptor Precursors' ~ 'PR Precursors',
                                                                   CellType_predict == 'AC/HC_Precurs' ~ 'AC/HC Precursors',
                                                                   CellType_predict == 'RPC' ~ 'RPCs',
                                                                   CellType_predict == 'Mesenchymal/RPE/Endothelial' ~ 'Endothelial',
                                                                   TRUE ~ CellType_predict)) %>%
  mutate(CellType_predict = case_when(!is.na(TabulaMurisCellType_predict) ~ 'Tabula Muris',
                                      TRUE ~ CellType_predict))

markers_summary <- scEiaD_2020_v01 %>% 
  tbl('wilcox_diff_testing') %>% 
  group_by(Base, Gene) %>% 
  summarise(auc_count = sum(AUC > 0.5), 
            pval = min(`p.value`),
            FDR = min(FDR),
            mean_auc = mean(AUC)) %>% 
  mutate(cluster = Base) %>% 
  as_tibble()

marker_info <-  markers_summary %>% 
  #filter(pval < 1) %>% 
  #group_by(cluster) %>% 
  #top_n(4, meauc) %>% 
  #filter(grepl('Amacrine|Rod|Cone|Retinal|Muller|Horizon|Bipol|RPE|Astro', cluster)) %>% 
  as_tibble() %>% 
  left_join(scEiaD_2020_v01 %>% 
              tbl('genes') %>% 
              collect() %>% 
              mutate(gene = str_extract(Gene, 'ENSG\\d+'))) %>% 
  # left_join(marker_list %>% 
  #             bind_rows(.id = 'cluster') %>% 
  #             group_by(cluster)) %>% 
  #             #top_n(5, L4) %>% 
  #             #filter(grepl('Amacrine|Rod|Cone|Retinal|Muller|Horizon|Bipol|Astro|RPE', cluster))) #%>% 
  left_join(scEiaD_2020_v01 %>% tbl('haystack') %>% collect())

genes <- marker_info$Gene %>% unique()

grouping_features <- 'CellType_predict'
exp_stats <- scEiaD_2020_v01 %>% tbl('grouped_stats') %>%
  filter(Gene %in% genes) %>%
  group_by_at(vars(one_of(c('Gene', grouping_features)))) %>%
  summarise(cpm = sum(cpm * cell_exp_ct) / sum(cell_exp_ct),
            cell_exp_ct = sum(cell_exp_ct, na.rm = TRUE)) %>%
  collect() %>%
  tidyr::drop_na() %>%
  left_join(., meta_filter %>%
              group_by_at(vars(one_of(grouping_features))) %>%
              summarise(Count = n())) %>%
  mutate(cell_exp_ct = ifelse(is.na(cell_exp_ct), 0, cell_exp_ct)) %>%
  mutate(`%` = round((cell_exp_ct / Count) * 100, 2),
         Expression = round(cpm * (`%` / 100), 2)) %>%
  select_at(vars(one_of(c('Gene', grouping_features, 'cell_exp_ct', 'Count', '%', 'Expression')))) %>%
  arrange(-Expression) 
# D_KL > 0.15
top_markers <- marker_info %>% 
  left_join(exp_stats %>% dplyr::rename(cluster = CellType_predict)) %>% 
  mutate(pval = as.numeric(pval), FDR = as.numeric(FDR)) %>% 
  filter((FDR < 1 |  mean_auc > 0.10), D_KL > 0.15, `%` > 10) %>% 
  group_by(Gene) %>% 
  slice_max(order_by = `%`) %>% 
  filter(grepl('Amacrine|Rod|Cone|Retinal|Muller|Horizon|Mast|Bipol|Astro|RPE|Micro', cluster), !grepl('^NA', Gene))

#save(top_markers, marker_info, exp_stats, file = '~/data/massive_integrated_eye_scRNA/top_markers.Rdata')
source('src/pubmed_query.R')
