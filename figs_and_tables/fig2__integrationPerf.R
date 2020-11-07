load('~/data/massive_integrated_eye_scRNA/metrics_2020_10_07.Rdata')

perf <- perf %>% 
  mutate(knn = case_when(method == 'scArches' ~ 7, TRUE ~ knn))
cluster_stats <- perf %>% select(clusterN, clusterMedian, dims, nf, knn, method, normalization, set) %>% unique() 

perf_well <- perf %>% unique() %>% filter(set == 'onlyWELL') %>% 
  select(Score, Group, Value, set, dims:normalization) %>% 
  filter(Score %in% c('LISI','Silhouette')) 

perf_tabula <- perf %>% unique() %>% filter(set == 'TabulaDroplet') %>% 
  select(Score, Group, Value, set, dims:normalization) %>% 
  filter(Score %in% c('LISI','Silhouette', 'ARI', 'PCR', 'NMI')) %>% 
  filter(method != 'CCA')

lisi <- perf_tabula %>% 
  filter(Score == 'LISI', nf == 2000, dims %in% c(8,30), knn == 7) %>%
  pivot_wider(values_from = Value, names_from = c('Group')) %>% 
  ggplot(aes(y=Batch, x=-Cluster, shape = normalization)) + geom_point(aes(color=method), size = 5) +
  cowplot::theme_cowplot() + scale_color_manual(values = pals::alphabet() %>% unname()) +
  ggtitle('LISI')

ari <- perf_tabula %>% 
  filter(Score == 'ARI', Group == 'CellType-Cluster', nf == 2000, dims %in% c(8, 30), knn == 7) %>%
  ggplot(aes(y=Value, x=normalization)) + geom_point(aes(color=method, shape = normalization), size = 5) +
  cowplot::theme_cowplot() + scale_color_manual(values = pals::alphabet() %>% unname()) + ylab('ARI (higher is better)') +
  ggtitle('ARI') + theme(legend.position="none")

silhouette <- perf_tabula %>% filter(Score == 'Silhouette', nf == 2000, dims %in% c(8,30), knn == 7) %>% 
  pivot_wider(values_from = Value, names_from = c('Group')) %>% 
  ggplot(aes(y=-Batch, x=Cluster, shape = normalization)) + geom_point(aes(color=method), size = 5) +
  cowplot::theme_cowplot() + scale_color_manual(values = pals::alphabet() %>% unname()) +
  ggtitle("Silhouette") + theme(legend.position="none")

legend <- get_legend(
  # create some space to the left of the legend
  lisi + theme(legend.box.margin = margin(0, 0, 0, 12))
)

zscore_tabula <- perf_tabula %>% 
  filter(nf == 2000, dims %in% c(8,30), knn == 7) %>% 
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
           scale(Silhouette_Cluster)[,1] +
           scale(`PCR_After-Before`)[,1]) %>% 
  left_join(cluster_stats) %>% 
  filter(!is.na(clusterN))

zscore_sum_all_methods <- zscore_tabula  %>% 
  arrange(-sumZScale) %>% 
  ggplot(aes(x=method, y = sumZScale, shape = normalization, group = dims)) + 
  geom_point(aes(color=method), size = 4, position = position_dodge(width = 1)) +
  cowplot::theme_cowplot() + 
  scale_color_manual(values = pals::alphabet() %>% unname()) + 
  coord_flip() 


zscore_droplet_scVI_optimize <- perf_tabula %>% 
  #filter(knn > 0.4) %>% 
  #filter(method == 'scVI') %>% 
  pivot_wider(names_from = c('Score','Group'), values_from = Value) %>% 
  mutate(sumZScale = 
           -scale(LISI_CellType)[,1] +
           scale(Silhouette_CellType)[,1] +
           scale(`NMI_CellType-Cluster`)[,1] +
           scale(`ARI_CellType-Cluster`)[,1] +
           
           scale(LISI_Batch)[,1] + 
           -scale(LISI_Cluster)[,1] +
           -scale(Silhouette_Batch)[,1] + 
           scale(Silhouette_Cluster)[,1] +
           
           scale(`PCR_After-Before`)[,1]
  ) %>% 
  left_join(cluster_stats) %>% 
  filter(!is.na(clusterN)) %>% 
  filter(clusterMedian < 5000) %>% arrange(-sumZScale)

zscore_onlyWell <- perf_well %>% 
  pivot_wider(names_from = c('Score','Group'), values_from = Value) %>% 
  mutate(sumZScale = 
           
           scale(LISI_Batch)[,1] + # Z score
           2 * -scale(LISI_Cluster)[,1] + 
           -scale(Silhouette_Batch)[,1] + # Z score
           2 * scale(Silhouette_Cluster)[,1] 
  ) %>% 
  arrange(-sumZScale) %>% DT::datatable()


# cowplot::plot_grid(cowplot::plot_grid(lisi + theme(legend.position="none"), silhouette, nrow = 1), cowplot::plot_grid(ari, legend, rel_widths = c(1,0.2), ncol = 2), nrow = 2)