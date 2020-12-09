library(tidyverse)

library(multiROC)


# # load('Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_spec_genes-0__n_features2000__counts__TabulaDroplet__batch__scVI__dims8__preFilter__mindist0.001__nneighbors500.umap.Rdata')
# # umapRef <- umap


pr_calculator <- function(test_predictions, partition){
  if (partition != 'All'){
    test_predictions <- test_predictions %>% 
      left_join(umapRef %>% select(Barcode, study_accession), by = 'Barcode') %>% 
      filter(study_accession == partition)
  } 
  cell_type2id <- test_predictions %>% select(cell_type_id, CellType) %>% distinct()
  target_cell_types <-unique(cell_type2id$CellType)
  pr_df <- test_predictions %>% 
    select(-cell_type_id, -CellType) %>% 
    rename(cell_type_id = true_cell_id) %>% 
    inner_join(cell_type2id) %>% 
    filter(CellType%in% target_cell_types) %>% 
    select(all_of(target_cell_types), cell_type_id, CellType) 
  # test_predictions <- predictionsC %>%  
  #   mutate(max_pred_prob = rowMaxs(.[,-(31:ncol(predictionsC))] %>% as.matrix )) %>% 
  #   select(-cell_type_id) %>% 
  #   rename(PredCellType = CellType) %>% 
  #   inner_join(umapRef %>% select(Barcode, TrueCellType = CellType)) %>% 
  #   filter(!is.na(TrueCellType)) %>% 
  #   mutate(pred_correct  = ifelse(PredCellType == TrueCellType, 'Correct', 'incorrect')) 
  # 
  # 
  # sum(test_predictions$pred_correct == 'Correct') / nrow(test_predictions) 
  # 
  # cell_type2id <- predictions %>% select(cell_type_id, CellType) %>% distinct %>% filter(CellType!='None')
  
  # pr_df <- test_predictions %>% 
  #   select(all_of(target_cell_types), PredCellType, TrueCellType) 
  one_hot_labled <- model.matrix(~0+as.factor(pr_df$CellType)) %>% as.data.frame
  colnames(one_hot_labled) <- strsplit(colnames(one_hot_labled), ')', fixed = T ) %>% 
    sapply(function(x)x[length(x)]) %>% paste0( '_true')
  pr_df <- pr_df %>% select(all_of(target_cell_types))
  colnames(pr_df) <- paste0(colnames(pr_df), '_pred_xgb')
  pr_df <- bind_cols( one_hot_labled, pr_df)
  colnames(pr_df) <- str_replace_all(colnames(pr_df), ' ','-')
  pr_res <- multi_pr(pr_df %>% as.data.frame, force_diag = T)
  
  out <- plot_pr_data(pr_res) %>% mutate(Study = partition)
  out
}

pr_list <- list()
for (i in c('All', umap$study_accession %>% unique())){
  #print(i)
  suppressWarnings({
    try({pr_list[[i]] <- pr_calculator(test_predictions, i)}, silent = TRUE)
    })
}

xgboost_pr_table <- pr_list %>% bind_rows() %>%
  filter(!Group%in% c('Macro', 'Micro'), Study == 'All') %>% 
  mutate(AUC = round(AUC, 2),
         `Cell Type` = str_replace_all(Group, '-', ' ')) %>%
  select(AUC, `Cell Type`) %>%
  unique() %>%
  arrange(-AUC) 

xgboost_pr_table_full <- pr_list %>% bind_rows() %>%
  filter(!Group%in% c('Macro', 'Micro')) %>% 
  mutate(AUC = round(AUC, 2),
         `Cell Type` = str_replace_all(Group, '-', ' ')) %>%
  select(AUC, `Cell Type`, Study) %>%
  unique() %>%
  arrange(-AUC) 
  
  
# test <- pr_calculator(predictions, 'SRP157927')
pr_split_curves <- pr_list %>% bind_rows() %>%
  filter(!Group%in% c('Macro', 'Micro'), Study != 'All') %>%
  mutate(AUC = round(AUC, 2),
         `Cell Type` = str_replace_all(Group, '-', ' ')) %>%
  left_join(., xgboost_pr_table %>% select(AUC_T = AUC, `Cell Type`), by = c('Cell Type')) %>% 
  #filter(`Cell Type` %in% cells_to_plot) %>%
  mutate(`Cell Type` = paste0(`Cell Type`, ' (', AUC_T, ')')) %>%     
  ggplot() +
  geom_line(aes(x = Recall,y= Precision, color = Study))+
  scale_color_manual(values = pals::cols25() %>% unname())+
  xlim(c(0, 1)) +
  ylim(c(0, 1))+
  facet_wrap(~`Cell Type`) +
  cowplot::theme_cowplot()
# 
# pr_curves <-pr_list %>%
#   filter(!Group%in% c('Macro', 'Micro'), Study == 'All') %>%
#   mutate(AUC = round(AUC, 2),
#          `Cell Type` = str_replace_all(Group, '-', ' ')) %>%
#   #filter(`Cell Type` %in% cells_to_plot) %>%
#   #mutate(`Cell Type` = paste0(`Cell Type`, ' (', AUC, ')')) %>%     
#   ggplot() +
#   geom_line(aes(x = Recall,y= Precision, color = Study))+
#   scale_color_manual(values = pals::cols25() %>% unname())+
#   xlim(c(0, 1)) +
#   ylim(c(0, 1))+
#   facet_wrap(~`Cell Type`) +
#   cowplot::theme_cowplot()
