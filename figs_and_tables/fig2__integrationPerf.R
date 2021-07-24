# celltype prediction accuracy
load('~/data/scEiaD/merged_xgboost_2021-03-17.Rdata')

accuracy <- accuracy %>% 
  filter(score != 0) %>% 
  mutate(
    norm = str_extract(file, 'transform-[a-zA-Z]+') %>% gsub('transform-','',.),
    nf = str_extract(file, 'n_features-\\d+') %>% gsub('n_features-', '', .) %>% as.numeric(),
    dims = str_extract(file, 'dims-\\d+') %>% gsub('dims-', '', .) %>% as.numeric(),
    method = str_extract(file, 'method-[a-zA-Z]+') %>% gsub('method-','',.),
    knn = str_extract(file, 'knn-\\d+\\.\\d+|knn-\\d+') %>% gsub('knn-', '', .) %>% as.numeric(),
    epochs = str_extract(file, 'epochs-\\d+') %>% gsub('epochs-','',.) %>% as.numeric())


load('~/data/scEiaD/giga_all_methods/gigascience/2021_07_22.perf.allmethod.Rdata')
perfO <- perf %>% filter(dims %in% c(8,30), knn == 20, method != 'scVIprojection')

load('~/data/scEiaD/merged_stats_2021-03-17.Rdata')

perf <- perf %>% filter(grepl('scVI',method)) 
cluster_stats <- perf %>% select(clusterN, clusterMedian, dims, nf, knn, method, normalization, set) %>% unique() 

perfO <- perfO %>% 
  mutate(normalization = case_when(method == 'scArches' ~ 'standard',
                                   TRUE ~ normalization))

#cluster_stats <- perf %>% select(clusterN, clusterMedian, dims, nf, knn, method, normalization, set) %>% unique() 
# 
# perf_well <- perf %>% unique() %>% filter(set == 'onlyWELL') %>% 
#   select(Score, Group, Value, set, dims:normalization) %>% 
#   filter(Score %in% c('LISI','Silhouette')) 

perf_all_methods <- perfO %>% unique() %>% 
  select(Score, Group, Value, set, dims:normalization) %>% 
  filter(Score %in% c('LISI','Silhouette', 'ARI', 'PCR', 'NMI')) 

perf_scVI <- perf %>% unique() %>%
  select(Score, Group, Value, set, dims:normalization) %>% 
  filter(Score %in% c('LISI','Silhouette', 'ARI', 'PCR', 'NMI')) 

lisi <- perf_all_methods %>% 
  filter(Score == 'LISI', nf == 2000, dims %in% c(8,30), knn == 20) %>%
  pivot_wider(values_from = Value, names_from = c('Group')) %>% 
  ggplot(aes(y=Batch, x=-Cluster, shape = normalization, alpha = dims)) + geom_point(aes(color=method), size = 5) +
  cowplot::theme_cowplot() + scale_color_manual(values = pals::alphabet() %>% unname()) +
  ggtitle('LISI') +
  scale_alpha_continuous(range = c(0.6, 1), breaks =c(8,30))

ari <- perf_all_methods %>% 
  filter(Score == 'ARI', Group == 'CellType-Cluster', nf == 2000, dims %in% c(8,30), knn == 20) %>%
  ggplot(aes(y=Value, x=normalization)) + geom_jitter(aes(color=method, shape = normalization, alpha = dims), size = 5, width = 0.2) +
  cowplot::theme_cowplot() + scale_color_manual(values = pals::alphabet() %>% unname()) + ylab('ARI') +
  ggtitle('ARI') + theme(legend.position="none") +
  scale_alpha_continuous(range = c(0.6, 1), breaks =c(8,30))

nmi <- perf_all_methods %>% 
  filter(Score == 'NMI', Group == 'CellType-Cluster', nf == 2000, dims %in% c(8,30), knn == 20) %>%
  ggplot(aes(y=Value, x=normalization)) + geom_jitter(aes(color=method, shape = normalization, alpha = dims), size = 5, width = 0.2) +
  cowplot::theme_cowplot() + scale_color_manual(values = pals::alphabet() %>% unname()) + ylab('NMI') +
  ggtitle('NMI') + theme(legend.position="none") +
  scale_alpha_continuous(range = c(0.6, 1), breaks =c(8,30))

pcr <- perf_all_methods %>%
  filter(Score == 'PCR', nf == 2000, dims %in% c(30), knn == 20) %>%
  ggplot(aes(y=Value, x=normalization)) + geom_jitter(aes(color=method, shape = normalization, alpha = dims), size = 5, width = 0.2) +
  cowplot::theme_cowplot() + scale_color_manual(values = pals::alphabet() %>% unname()) + ylab('PCR') +
  ggtitle('PCR') + theme(legend.position="none") +
  scale_alpha_continuous(range = c(0.6, 1), breaks =c(8,30))

silhouette <- perf_all_methods %>% filter(Score == 'Silhouette', nf == 2000, dims %in% c(6,8,10, 30), knn == 20) %>% 
  pivot_wider(values_from = Value, names_from = c('Group')) %>% 
  ggplot(aes(y=-Batch, x=Cluster, shape = normalization, alpha = dims)) + geom_point(aes(color=method), size = 5) +
  cowplot::theme_cowplot() + scale_color_manual(values = pals::alphabet() %>% unname()) +
  ggtitle("Silhouette") + theme(legend.position="none") +
  scale_alpha_continuous(range = c(0.6, 1), breaks =c(8,30))

legend <- get_legend(
  # create some space to the left of the legend
  lisi + theme(legend.box.margin = margin(0, 0, 0, 12))
)

zscore_allMethods <- perf_all_methods %>% 
  filter(nf == 2000, dims %in% c(8, 30), knn == 20) %>% 
  pivot_wider(names_from = c('Score','Group'), values_from = Value) %>% 
  mutate(sumZScale = 
           -scale(LISI_CellType)[,1]  +
           #-scale(LISI_SubCellType)[,1]  +
           scale(Silhouette_CellType)[,1] +
           #scale(Silhouette_SubCellType)[,1] +
           
           scale(`NMI_CellType-Cluster`)[,1] +
           #scale(`NMI_SubCellType-Cluster`)[,1] +
           scale(`ARI_CellType-Cluster`)[,1] +
           #scale(`ARI_SubCellType-Cluster`)[,1] +
           
           scale(LISI_Batch)[,1] + # Z score
           -scale(LISI_Cluster)[,1] +
           -scale(Silhouette_Batch)[,1] + # Z score
           scale(Silhouette_Cluster)[,1]) 


# scale zscore sum
purity <- 3; balance <- 1; mixing <- 1
zscore_allMethods_highPurity_plot <- perf_all_methods %>%
  filter(nf == 2000, dims %in% c(8, 30), knn == 20) %>%
  pivot_wider(names_from = c('Score','Group'), values_from = Value) %>%
  mutate(sumZScale =
           -scale(LISI_CellType)[,1] * purity +
           scale(Silhouette_CellType)[,1] * purity  +
           scale(Silhouette_Cluster)[,1] * purity  +
           -scale(LISI_Cluster)[,1] * purity  +
           
           scale(`NMI_CellType-Cluster`)[,1] * balance  +
           scale(`ARI_CellType-Cluster`)[,1]* balance  +
           
           scale(LISI_Batch)[,1] * mixing  + # Z score
           -scale(Silhouette_Batch)[,1] * mixing) %>% 
  ggplot(aes(x=method, y = sumZScale, shape = normalization, group = dims, alpha = dims)) + 
  geom_point(aes(color=method), size = 4, position = position_dodge(width = 1)) +
  cowplot::theme_cowplot() + 
  scale_color_manual(values = pals::alphabet() %>% unname(), guide = FALSE) + 
  coord_flip() + 
  scale_alpha_continuous(range = c(0.4, 1), breaks =c(8,30))

purity <- 1; balance <- 1; mixing <- 3
zscore_allMethods_highMixing_plot <- perf_all_methods %>%
  filter(nf == 2000, dims %in% c(8, 30), knn == 20) %>%
  pivot_wider(names_from = c('Score','Group'), values_from = Value) %>%
  mutate(sumZScale =
           -scale(LISI_CellType)[,1] * purity +
           scale(Silhouette_CellType)[,1] * purity  +
           scale(Silhouette_Cluster)[,1] * purity  +
           -scale(LISI_Cluster)[,1] * purity  +
           
           scale(`NMI_CellType-Cluster`)[,1] * balance  +
           scale(`ARI_CellType-Cluster`)[,1]* balance  +
           
           scale(LISI_Batch)[,1] * mixing  + # Z score
           -scale(Silhouette_Batch)[,1] * mixing) %>% 
  ggplot(aes(x=method, y = sumZScale, shape = normalization, group = dims, alpha = dims)) + 
  geom_point(aes(color=method), size = 4, position = position_dodge(width = 1)) +
  cowplot::theme_cowplot() + 
  scale_color_manual(values = pals::alphabet() %>% unname(), guide = FALSE) + 
  coord_flip() + 
  scale_alpha_continuous(range = c(0.4, 1), breaks =c(8,30))


# full random
zscore_allMethods_random <- list()
for (i in seq(1:1000)){
  zscore_allMethods_random[[i]] <- perf_all_methods %>%
  filter(nf == 2000, dims %in% c(8, 30), knn == 20) %>%
  pivot_wider(names_from = c('Score','Group'), values_from = Value) %>%
  mutate(sumZScale =
           -scale(LISI_CellType)[,1] * (sample(1:100, 1) / 10) +
           scale(Silhouette_CellType)[,1] * (sample(1:100, 1) / 10)  +
           scale(Silhouette_Cluster)[,1] * (sample(1:100, 1) / 10)  +
           -scale(LISI_Cluster)[,1] * (sample(1:100, 1) / 10)  +
           
           scale(`NMI_CellType-Cluster`)[,1] * (sample(1:100, 1) / 10)  +
           scale(`ARI_CellType-Cluster`)[,1]* (sample(1:100, 1) / 10)  +
           
           scale(LISI_Batch)[,1] * (sample(1:100, 1) / 10)  + # Z score
           -scale(Silhouette_Batch)[,1] * (sample(1:100, 1) / 10))  %>% 
    mutate(Rank = rank(-sumZScale),
           Bootstrap = i)
}
mean_rank <- zscore_allMethods_random %>% bind_rows() %>% group_by(dims, nf, knn, method, normalization) %>% summarise(Rank = mean(Rank)) %>% mutate(params = paste(method, normalization, dims, nf)) %>% arrange(Rank)
zscore_ridges_plot <- zscore_allMethods_random %>% bind_rows() %>%  mutate(params = paste(method, normalization, dims, nf)) %>% ungroup() %>% select(params, Rank) %>% mutate(params = factor(params, levels = mean_rank$params)) %>% ggplot(aes(x=Rank,y=params)) + ggridges::geom_density_ridges() + cowplot::theme_cowplot() + scale_y_discrete(limits=rev)

zscore_sum_all_methods <- zscore_allMethods  %>% 
  filter(!grepl('proj', method)) %>% 
  mutate(normalization = case_when(method == 'scArches' ~ 'standard',
                                   TRUE ~ normalization)) %>% 
  arrange(-sumZScale) %>% 
  ggplot(aes(x=method, y = sumZScale, shape = normalization, group = dims, alpha = dims)) + 
  geom_point(aes(color=method), size = 4, position = position_dodge(width = 1)) +
  cowplot::theme_cowplot() + 
  scale_color_manual(values = pals::alphabet() %>% unname(), guide = FALSE) + 
  coord_flip() + 
  scale_alpha_continuous(range = c(0.4, 1), breaks =c(8,30))

zscore_droplet_scVI_optimize <- perf_scVI %>% 
  #filter(method == 'scVIprojection') %>% 
  filter(grepl('scVI', method), epochs == 5) %>% 
  pivot_wider(names_from = c('Score','Group'), values_from = Value) %>% 
  left_join(accuracy %>% group_by(norm, nf, dims, method, knn, epochs) %>% summarise(xgboost_score = mean(score))) %>% 
  mutate( sumZScale = 
            -scale(LISI_CellType)[,1] +
            scale(Silhouette_CellType)[,1] +
            scale(`NMI_CellType-Cluster`)[,1] +
            scale(`ARI_CellType-Cluster`)[,1] +
            
            scale(`xgboost_score`)[,1] +
            
            scale(LISI_Batch)[,1] + 
            -scale(LISI_Cluster)[,1] +
            -scale(Silhouette_Batch)[,1]  + 
            scale(Silhouette_Cluster)[,1] #+
          
          #scale(`PCR_After-Before`)[,1]
  ) %>% 
  #left_join(cluster_stats) %>% 
  #filter(!is.na(clusterN)) %>% 
  #filter(knn == 0.6, method == 'scVI') %>% 
  #filter(!is.na(clusterN)) %>% 
  mutate(nf = as.factor(nf),
         `scVI latent dims` = as.factor(dims),
         method = gsub('scVIprojection', 'scVI-projection',method),
         method = gsub('scVI$', 'scVI-standard', method),
         method = factor(method, levels = c('scVI-standard','scVI-projection'))) %>% 
  ggplot(aes(x=sumZScale, y = nf, color = `scVI latent dims`, fill = method)) + 
  #ggbeeswarm::geom_quasirandom(groupOnX=FALSE) +
  geom_boxplot() +
  scale_color_manual(values = pals::brewer.set1(n = 10) %>% unname()) +
  scale_fill_manual(values = c('white','gray')) +
  cowplot::theme_cowplot() +
  ylab('Number of\nHVGs')

# zscore_onlyWell <- perf_well %>% 
#   pivot_wider(names_from = c('Score','Group'), values_from = Value) %>% 
#   mutate(sumZScale = 
#            scale(LISI_Batch)[,1] + # Z score
#            2 * -scale(LISI_Cluster)[,1] + 
#            -scale(Silhouette_Batch)[,1] + # Z score
#            2 * scale(Silhouette_Cluster)[,1] 
#   ) %>% 
#   arrange(-sumZScale) %>% DT::datatable()


# cowplot::plot_grid(cowplot::plot_grid(lisi + theme(legend.position="none"), silhouette, nrow = 1), cowplot::plot_grid(ari, legend, rel_widths = c(1,0.2), ncol = 2), nrow = 2)