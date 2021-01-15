load('~/data/scEiaD//metrics_2020_11_28.Rdata')
perfO <- perf

load('~/data/scEiaD//merged_stats_2021_01_01.Rdata')

perf <- perf %>% 
  mutate(knn = case_when(method == 'scArches' ~ 7, TRUE ~ knn))
cluster_stats <- perf %>% select(clusterN, clusterMedian, dims, nf, knn, method, normalization, set) %>% unique() 

perfO <- perfO %>% 
  mutate(knn = case_when(method == 'scArches' ~ 7, TRUE ~ knn),
         normalization = case_when(method == 'scArches' ~ 'standard',
                                   TRUE ~ normalization))

#cluster_stats <- perf %>% select(clusterN, clusterMedian, dims, nf, knn, method, normalization, set) %>% unique() 
# 
# perf_well <- perf %>% unique() %>% filter(set == 'onlyWELL') %>% 
#   select(Score, Group, Value, set, dims:normalization) %>% 
#   filter(Score %in% c('LISI','Silhouette')) 

perf_tabula <- perfO %>% unique() %>% filter(set %in% c('TabulaDroplet')) %>% 
  select(Score, Group, Value, set, dims:normalization) %>% 
  filter(Score %in% c('LISI','Silhouette', 'ARI', 'PCR', 'NMI')) 

perf_scVI <- perf %>% unique() %>% filter(set %in% c('TabulaDroplet', 'universe')) %>% 
  select(Score, Group, Value, set, dims:normalization) %>% 
  filter(Score %in% c('LISI','Silhouette', 'ARI', 'PCR', 'NMI')) 

lisi <- perf_tabula %>% 
  filter(Score == 'LISI', nf == 2000, dims %in% c(8,30), knn == 7) %>%
  pivot_wider(values_from = Value, names_from = c('Group')) %>% 
  ggplot(aes(y=Batch, x=-Cluster, shape = normalization, alpha = dims)) + geom_point(aes(color=method), size = 5) +
  cowplot::theme_cowplot() + scale_color_manual(values = pals::alphabet() %>% unname()) +
  ggtitle('LISI') +
  scale_alpha_continuous(range = c(0.6, 1), breaks =c(8,30))

ari <- perf_tabula %>% 
  filter(Score == 'ARI', Group == 'CellType-Cluster', nf == 2000, dims %in% c(8,30), knn == 7) %>%
  ggplot(aes(y=Value, x=normalization)) + geom_jitter(aes(color=method, shape = normalization, alpha = dims), size = 5, width = 0.2) +
  cowplot::theme_cowplot() + scale_color_manual(values = pals::alphabet() %>% unname()) + ylab('ARI') +
  ggtitle('ARI') + theme(legend.position="none") +
  scale_alpha_continuous(range = c(0.6, 1), breaks =c(8,30))

nmi <- perf_tabula %>% 
  filter(Score == 'NMI', Group == 'CellType-Cluster', nf == 2000, dims %in% c(8,30), knn == 7) %>%
  ggplot(aes(y=Value, x=normalization)) + geom_jitter(aes(color=method, shape = normalization, alpha = dims), size = 5, width = 0.2) +
  cowplot::theme_cowplot() + scale_color_manual(values = pals::alphabet() %>% unname()) + ylab('NMI') +
  ggtitle('NMI') + theme(legend.position="none") +
  scale_alpha_continuous(range = c(0.6, 1), breaks =c(8,30))

pcr <- perf_tabula %>%
  filter(Score == 'PCR', nf == 2000, dims %in% c(30), knn == 7) %>%
  ggplot(aes(y=Value, x=normalization)) + geom_jitter(aes(color=method, shape = normalization, alpha = dims), size = 5, width = 0.2) +
  cowplot::theme_cowplot() + scale_color_manual(values = pals::alphabet() %>% unname()) + ylab('PCR') +
  ggtitle('PCR') + theme(legend.position="none") +
  scale_alpha_continuous(range = c(0.6, 1), breaks =c(8,30))

silhouette <- perf_tabula %>% filter(Score == 'Silhouette', nf == 2000, dims %in% c(8,30), knn == 7) %>% 
  pivot_wider(values_from = Value, names_from = c('Group')) %>% 
  ggplot(aes(y=-Batch, x=Cluster, shape = normalization, alpha = dims)) + geom_point(aes(color=method), size = 5) +
  cowplot::theme_cowplot() + scale_color_manual(values = pals::alphabet() %>% unname()) +
  ggtitle("Silhouette") + theme(legend.position="none") +
  scale_alpha_continuous(range = c(0.6, 1), breaks =c(8,30))

legend <- get_legend(
  # create some space to the left of the legend
  lisi + theme(legend.box.margin = margin(0, 0, 0, 12))
)

zscore_tabula <- perf_tabula %>% 
  filter(nf == 2000, dims %in% c(8, 30), knn == 7) %>% 
  pivot_wider(names_from = c('Score','Group'), values_from = Value) %>% 
  mutate(sumZScale = 
           -scale(LISI_CellType)[,1] +
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

zscore_sum_all_methods <- zscore_tabula  %>% 
  filter(!grepl('proj', method)) %>% 
  mutate(normalization = case_when(method == 'scArches' ~ 'standard',
                                   TRUE ~ normalization)) %>% 
  arrange(-sumZScale) %>% 
  ggplot(aes(x=method, y = sumZScale, shape = normalization, group = dims, alpha = dims)) + 
  geom_point(aes(color=method), size = 4, position = position_dodge(width = 1)) +
  cowplot::theme_cowplot() + 
  scale_color_manual(values = pals::alphabet() %>% unname(), guide = FALSE) + 
  coord_flip() + 
  scale_alpha_continuous(range = c(0.6, 1), breaks =c(8,30))


zscore_droplet_scVI_optimize <- perf_scVI%>% 
  filter(knn > 0.4, knn < 10) %>% 
  filter(grepl('scVI', method)) %>% 
  pivot_wider(names_from = c('Score','Group'), values_from = Value) %>% 
  mutate(sumZScale = 
           -scale(LISI_CellType)[,1] +
           scale(Silhouette_CellType)[,1] +
           scale(`NMI_CellType-Cluster`)[,1] +
           scale(`ARI_CellType-Cluster`)[,1] +
           
           #scale(LISI_Batch)[,1] * 0.5 + 
           -scale(LISI_Cluster)[,1] +
           #-scale(Silhouette_Batch)[,1] * 0.5  + 
           scale(Silhouette_Cluster)[,1] #+
           
           #scale(`PCR_After-Before`)[,1]
  ) %>% 
  left_join(cluster_stats) %>% 
  filter(!is.na(clusterN)) %>% 
  #filter(knn == 0.6, method == 'scVI') %>% 
  filter(!is.na(clusterN)) %>% 
  mutate(nf = as.factor(nf),
         `scVI latent dims` = as.factor(dims)) %>% 
  ggplot(aes(x=sumZScale, y = nf, color = `scVI latent dims`, shape = method)) + 
  ggbeeswarm::geom_quasirandom(size = 3, groupOnX=FALSE) +
  scale_color_manual(values = pals::brewer.set1(n = 10) %>% unname()) +
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