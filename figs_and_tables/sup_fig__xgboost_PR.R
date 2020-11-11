library(tidyverse)
library(matrixStats)
library(multiROC)
load('~/PREDICTIONSMus_musculus_Macaca_fascicularis_Homo_sapiens__n_spec_genes-0__n_features2000__counts__TabulaDroplet__batch__scVI__dims8__preFilter__mindist0.001__nneighbors500.Rdata')
load('Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_spec_genes-0__n_features2000__counts__TabulaDroplet__batch__scVI__dims8__preFilter__mindist0.001__nneighbors500.umap.Rdata')
umapRef <- umap
test_predictions <- predictions %>%  
  mutate(max_pred_prob = rowMaxs(.[,-(31:33)] %>% as.matrix )) %>% 
  select(-cell_type_id) %>% 
  rename(PredCellType = CellType) %>% 
  inner_join(umapRef %>% select(Barcode, TrueCellType = CellType)) %>% 
  filter(!is.na(TrueCellType)) %>% 
  mutate(pred_correct  = ifelse(PredCellType == TrueCellType, 'Correct', 'incorrect')) 
sum(test_predictions$pred_correct == 'Correct') / nrow(test_predictions) 

cell_type2id <- predictions %>% select(cell_type_id, CellType) %>% distinct %>% filter(CellType!='None')
target_cell_types <-unique(cell_type2id$CellType)
pr_df <- test_predictions %>% 
  select(all_of(target_cell_types), PredCellType, TrueCellType) 
pr_df <- filter(pr_df, TrueCellType  %in% colnames(pr_df)) 

one_hot_labled <- model.matrix(~0+as.factor(pr_df$TrueCellType)) %>% as.data.frame
colnames(one_hot_labled) <- strsplit(colnames(one_hot_labled), ')', fixed = T ) %>% 
  sapply(function(x)x[length(x)]) %>% paste0( '_true')
pr_df <- pr_df %>% select(all_of(target_cell_types))
colnames(pr_df) <- paste0(colnames(pr_df), '_pred_xgb') 

pr_df <- bind_cols( one_hot_labled, pr_df)
colnames(pr_df) <- str_replace_all(colnames(pr_df), ' ','-')
pr_res <- multi_pr(pr_df %>% as.data.frame, force_diag = T)
cells_to_plot <-c("Rods","Retinal Ganglion Cells","Bipolar Cells","Amacrine Cells","Muller Glia","Cones","Horizontal Cells")
plot_pr_df <- plot_pr_data(pr_res) %>% 
  filter(!Group%in% c('Macro', 'Micro')) %>% 
  mutate(AUC = round(AUC, 2), 
         `Cell Type` = str_replace_all(Group, '-', ' ')) %>% 
  #filter(`Cell Type` %in% cells_to_plot) %>%            
  mutate(`Cell Type` = paste0(`Cell Type`, ' (', AUC, ')'))

xgboost_pr <- ggplot(plot_pr_df) + 
  geom_line(aes(x = Recall,y= Precision))+
  scale_color_manual(values = pals::cols25() %>% unname())+
  xlim(c(.7, 1)) +
  ylim(c(.7, 1))+
  facet_wrap(~`Cell Type`) +
  cowplot::theme_cowplot()



