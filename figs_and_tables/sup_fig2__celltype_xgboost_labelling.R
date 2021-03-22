ct_processing <-  meta_filter %>% 
  #filter(CellType %in% common) %>% 
  filter(!is.na(CellType), CellType_predict != 'Tabula Muris') %>% 
  filter(!(!Platform %in% c('10xv2','10xv3','DropSeq') & CellType_predict == 'Unlabelled'  )) %>% 
  filter(!(!Platform %in% c('10xv2','10xv3','DropSeq') & CellType_predict == 'Muller Glia Progenitor'  )) %>% 
  mutate(CellType = case_when(is.na(CellType) ~ 'Unlabelled', 
                              CellType == 'AC/HC_Precurs' ~ 'AC/HC Precursors',
                              TRUE ~ CellType)) %>% 
  mutate(CellType_predict = case_when(CellType_predict == 'RPC' ~ 'RPCs',
                                      CellType_predict == 'AC/HC_Precurs' ~ 'AC/HC Precursors',
                                      grepl('Mesenchymal', CellType_predict) ~ 'Endothelial',
                                      is.na(CellType_predict) ~ 'Unlabelled', 
                                      TRUE ~ CellType_predict)) 


ct_processing_NF <-  meta_filter %>% 
  mutate(CellType = case_when(is.na(CellType) ~ 'Unlabelled', 
                              CellType == 'AC/HC_Precurs' ~ 'AC/HC Precursors',
                              TRUE ~ CellType)) %>% 
  mutate(CellType_predict = case_when(CellType_predict == 'RPC' ~ 'RPCs',
                                      CellType_predict == 'AC/HC_Precurs' ~ 'AC/HC Precursors',
                                      grepl('Mesenchymal', CellType_predict) ~ 'Endothelial',
                                      is.na(CellType_predict) ~ 'Unlabelled', 
                                      TRUE ~ CellType_predict)) 

##########
# Cell Counts
predictedCT <- ct_processing_NF %>% 
  group_by(organism, CellType_predict) %>% 
  summarise(`Published Count` = n()) %>% pivot_wider(values_from = `Published Count`, names_from = c(organism))
colnames(predictedCT) <- c('CellType','HS Transferred','MF Transferred','MM Transferred')
prelabelledCT <- ct_processing_NF %>% 
  group_by(organism, CellType) %>% 
  summarise(`Published Count` = n()) %>% pivot_wider(values_from = `Published Count`, names_from = c(organism))
colnames(prelabelledCT) <- c('CellType','HS Published','MF Published','MM Published')
joinedCT <- left_join(prelabelledCT, predictedCT)
joinedCT[is.na(joinedCT)] <- 0
ctTable <- joinedCT %>% flextable()
##########

##########
# Study Counts
predictedSA <- ct_processing_NF %>% 
  filter(!is.na(CellType_predict)) %>% 
  select(CellType_predict, organism, study_accession) %>% unique() %>% 
  group_by(organism, CellType_predict) %>% 
  summarise(`Study Count` = n()) %>% pivot_wider(values_from = `Study Count`, names_from = c(organism))
colnames(predictedSA) <- c('CellType','HS Studies (transferred)','MF Studies (transferred)','MM Studies (transferred)')
prelabelledSA <- ct_processing_NF %>% 
  filter(!is.na(CellType)) %>% 
  select(CellType, organism, study_accession) %>% unique() %>% 
  group_by(organism, CellType) %>% 
  summarise(`Study Count` = n()) %>% pivot_wider(values_from = `Study Count`, names_from = c(organism))
colnames(prelabelledSA) <- c('CellType','HS Studies (published)','MF Studies (published)','MM Studies (published)')
joinedSA <- left_join(prelabelledSA, predictedSA)
joinedSA[is.na(joinedSA)] <- 0
saTable <- joinedSA %>% flextable()
##########
# 
# ct_alluvial <- ct_processing %>%
#   select(CellType, CellType_predict) %>%
#   group_by(CellType, CellType_predict) %>%
#   summarise(Count = n()) %>%
#   ggplot(aes(y = sqrt(Count), axis1 = CellType, axis2 = CellType_predict)) +
#   geom_alluvium(aes(fill = `CellType`), alpha = 0.8) +
#   geom_stratum(alpha = 0) +
#   ggrepel::geom_label_repel(stat = "stratum",
#                             direction = 'x',
#                             fill = alpha(c("white"),0.5),
#                             size = 2,
#                             aes(label = after_stat(stratum))) +
#   coord_cartesian(xlim = c(0.5,2.5)) +
#   theme_void() +
#   scale_fill_manual(values = c(pals::polychrome() %>% unname(),
#                                pals::alphabet() %>% unname())) +
#   theme(legend.position = "none")

#########
#  confusion matrix
########
ct_mat <- ct_processing %>%
  select(CellType, CellType_predict) %>%
  group_by(CellType, CellType_predict) %>%
  summarise(Count = n()) %>% 
  mutate(Ratio = Count / sum(Count)) %>% 
  select(-Count) %>% 
  #filter(!CellType_predict %in% c('Astrocytes', 'Artery',  'Choriocapillaris', 'Tabula Muris','B-Cell','T-Cell', 'Smooth Muscle Cell', 'Unlabelled'), !CellType %in% c('Astrocytes','Artery','Choriocapillaris','Tabula Muris','B-Cell','T-Cell', 'Smooth Muscle Cell')) %>% 
  #filter(!CellType_predict %in% c('Astrocytes', 'Artery',  'Choriocapillaris', 'Tabula Muris','Smooth Muscle Cell', 'Unlabelled'), !CellType %in% c('Astrocytes','Artery','Choriocapillaris', 'Smooth Muscle Cell')) %>% 
  pivot_wider(values_from = Ratio, names_from = CellType_predict) 
ct_mat[is.na(ct_mat)] <- 0
#ct_mat <- data.frame(ct_mat)


rowN <- ct_mat$CellType
ct_mat <- ct_mat[,-1] 
ct_mat <- ct_mat[, colnames(ct_mat) %>% sort()]
ct_mat <- ct_mat %>% as.matrix()
row.names(ct_mat) <- rowN
ct_mat <- ct_mat[row.names(ct_mat) %>% rev(), ]

ct_confusion <- Heatmap(ct_mat, cluster_rows = FALSE, cluster_columns = FALSE, col=viridis(20), name =  'Recall')
