library(pool)
library(RSQLite)
library(tidyverse)
library(glue)
scEiaD_2020_v01 <- dbPool(drv = SQLite(), dbname = "~/data/massive_integrated_eye_scRNA/MOARTABLES__anthology_limmaFALSE___5000-transform-counts-universe-batch-scVIprojectionSO-8-0.1-500-0.6.sqlite", idleTimeout = 3600000)

meta_filter <- fst::read_fst('~/git/plaeApp/inst/app/www/meta_filter.fst') %>% as_tibble()
meta_filter <- meta_filter %>% mutate(CellType_predict = case_when(CellType_predict == 'Photoreceptor Precursors' ~ 'PR Precursors',
                                                                   CellType_predict == 'AC/HC_Precurs' ~ 'AC/HC Precursors',
                                                                   CellType_predict == 'RPC' ~ 'RPCs',
                                                                   CellType_predict == 'Mesenchymal/RPE/Endothelial' ~ 'Endothelial',
                                                                   TRUE ~ CellType_predict)) %>%
  mutate(CellType_predict = case_when(!is.na(TabulaMurisCellType_predict) ~ 'Tabula Muris',
                                      TRUE ~ CellType_predict))
load('~/data/massive_integrated_eye_scRNA/n_features-5000__transform-counts__partition-universe__covariate-batch__method-scVIprojectionSO__dims-8__knn-0.6__neighbors-500__dist-0.1__CellType_predict.sceWilcox.Rdata')
#genes <- c('ENSG00000144834','ENSG00000139970','ENSG00000109846','ENSG00000026025','ENSG00000067715','ENSG00000123307','ENSG00000139053','ENSG00000143774','ENSG00000105372','ENSG00000130066','ENSG00000104435','ENSG00000139970','ENSG00000050165','ENSG00000135821','ENSG00000104435','ENSG00000131711','ENSG00000016082','ENSG00000198668','ENSG00000149489','ENSG00000129535','ENSG00000105372','ENSG00000140988')

#markers <- scEiaD_2020_v01 %>% tbl('genes') %>% collect() %>% mutate(ENS = str_extract(Gene, 'ENSG\\d+')) %>% filter(ENS %in% genes) %>% pull(Gene)

markers_summary <- list()
for (i in names(markers_wilcox)){
  print(i)
  temp <- markers_wilcox[[i]][,4:ncol(markers_wilcox[[i]])] %>% as.matrix()
  count <-  apply(temp, 1, function(x) sum(x > 0.7, na.rm = TRUE))
  pval = markers_wilcox[[i]][,1]
  FDR = markers_wilcox[[i]][,2]
  med_auc <- apply(temp, 1, function(x) median(x, na.rm = TRUE))
  mean_auc <- apply(temp, 1, function(x) mean(x, na.rm = TRUE))
  gene <- row.names(temp)
  cluster <- i
  markers_summary[[i]] <- cbind(gene, pval, FDR, count, med_auc, mean_auc, cluster) %>% as_tibble()
}
markers_summary <- markers_summary %>% bind_rows() %>% mutate(count = as.numeric(count), med_auc = as.numeric(med_auc), mean_auc = as.numeric(mean_auc))


# marker_info_DT <- marker_list %>% 
#   bind_rows(.id = 'cluster') %>% 
#   group_by(cluster) %>% 
#   #top_n(5, L4) %>% 
#   filter(grepl('Amacrine|Rod|Cone|Retinal|Muller|Horizon|Bipol|Astro', cluster)) %>% 
#   filter(!Gene %in%
#            (marker_list %>% 
#            bind_rows(.id = 'cluster') %>% 
#            group_by(cluster) %>% 
#            top_n(5, L4) %>% 
#            filter(grepl('Amacrine|Rod|Cone|Retinal|Muller|Horizon|Bipol|Astro', cluster)) %>% 
#            group_by(Gene) %>% 
#            summarise(Count = n()) %>% 
#            arrange(-Count) %>% 
#            filter(Count > 1) %>% 
#            pull(Gene))) %>% 
#   ungroup() %>% 
#   group_by(cluster) %>% 
#   top_n(3, L4)  


marker_list <- list()
for (i in (meta_filter$CellType_predict %>% unique())){
  if (is.na(i)){
    next()
  } else {
    var <- i
    print(var)
    test <- glue("%-1*{var}%")
    ct <- glue("%{var}%")
    markers <- scEiaD_2020_v01 %>% 
      tbl('PB_results') %>% 
      filter(PB_Test == 'Pairwise CellType (Predict) against CellType (Predict)', 
             test %like% (!!ct), 
             logCPM > 5)  %>% 
      mutate(logFC = case_when(comparison %like% (!!test) ~ abs(logFC), 
                               TRUE ~ logFC)) %>%  group_by(Gene) %>% 
      summarise(L1 = sum(logFC > 1), L4 = sum(logFC > 4), 
                Diff = mean(logFC)) %>% 
      arrange(-Diff, -L4) %>% 
      collect()
    marker_list[[var]] <- markers
  }
}

marker_info <-  markers_summary %>% 
  #filter(pval < 1) %>% 
  group_by(cluster) %>% 
  #top_n(4, meauc) %>% 
  #filter(grepl('Amacrine|Rod|Cone|Retinal|Muller|Horizon|Bipol|RPE|Astro', cluster)) %>% 
  data.frame() %>% 
  left_join(scEiaD_2020_v01 %>% 
              tbl('genes') %>% 
              collect() %>% 
              mutate(gene = str_extract(Gene, 'ENSG\\d+'))) %>% 
  left_join(marker_list %>% 
              bind_rows(.id = 'cluster') %>% 
              group_by(cluster)) %>% 
              #top_n(5, L4) %>% 
              #filter(grepl('Amacrine|Rod|Cone|Retinal|Muller|Horizon|Bipol|Astro|RPE', cluster))) #%>% 
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
  filter((FDR < 1 |  med_auc > 0.20 | L1 > 50), D_KL > 0.15, `%` > 20) %>% 
  group_by(Gene) %>% 
  slice_max(order_by = med_auc) %>% 
  filter(grepl('Amacrine|Rod|Cone|Retinal|Muller|Horizon|Bipol|Astro|RPE|Micro', cluster))

#save(top_markers, marker_info, exp_stats, file = '~/data/massive_integrated_eye_scRNA/top_markers.Rdata')
source('src/pubmed_query.R')
