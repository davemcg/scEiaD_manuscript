library(ggtext)
library(glue)
make_dotplot <- function(input, db, meta_filter, cat_to_color_df){
  # ### this makes it a little easier to test
  # input <- list()
  # input[['dotplot_Gene']] <- markers
  # input[['dotplot_groups']] <- c('CellType_predict')
  # #input[['dotplot_filter_cat']] <- 'CellType_predict'
  # #input[['dotplot_filter_on']] <- (meta_filter %>%
  # #                                   filter(grepl('Amacrine|Rod|Cone|Retinal|Muller|Horizon|Bipol|RPE|Astro', CellType_predict)) %>% pull(CellType_predict) %>% unique())
  # db <- scEiaD_2020_v01
  # 
  gene <- input$dotplot_Gene
  grouping_features <- input$dotplot_groups
  
  # ct order
  ct_order <- data.frame(CT = c('Cones','Rods','RPE', 'Retinal Ganglion Cells', 'Astrocytes','Amacrine Cells','Bipolar Cells','Rod Bipolar Cells','Horizontal Cells','Muller Glia', 'Microglia'),
                         Order = c(1,2,3,4,5,6,7,8,9,10, 11))
  if (input$dotplot_filter_cat != ''){
    dotplot_data <- db %>% tbl('grouped_stats') %>%
      filter(Gene %in% gene) %>%
      collect() %>%
      filter(!!as.symbol(input$dotplot_filter_cat) %in% input$dotplot_filter_on)
    
  } else {
    dotplot_data <- db %>% tbl('grouped_stats') %>%
      filter(Gene %in% gene) %>%
      collect()
  }
  
  dotplot_data <- dotplot_data %>%
    filter(Gene %in% gene) %>%
    group_by_at(vars(one_of(c('Gene', grouping_features)))) %>%
    summarise(cpm = sum(cpm * cell_exp_ct) / sum(cell_exp_ct),
              cell_exp_ct = sum(cell_exp_ct, na.rm = TRUE)) %>%
    collect() %>%
    tidyr::drop_na() %>%
    full_join(., meta_filter %>%
                group_by_at(vars(one_of(grouping_features))) %>%
                summarise(Count = n())) %>%
    mutate(cell_exp_ct = ifelse(is.na(cell_exp_ct), 0, cell_exp_ct)) %>%
    mutate(`%` = round((cell_exp_ct / Count) * 100, 2),
           Expression = cpm * (`%` / 100)) %>%
    select_at(vars(one_of(c('Gene', grouping_features, 'cell_exp_ct', 'Count', '%', 'Expression')))) %>%
    filter(!is.na(Gene))
  #meanE <- dotplot_data %>% group_by(Gene) %>% summarise(meanE = mean(Expression))
  #dotplot_data <- dotplot_data %>% left_join(meanE) %>% mutate(Expression = Expression - meanE)
  if (length(grouping_features) == 2){
    colnames(dotplot_data)[c(2,3)] <- c("Group1","Group2")
    dotplot_data <- dotplot_data %>%
      mutate(Column = paste(Group1,Group2, sep = ':'),
             Column = case_when(grepl('^\\d', Column) ~ paste0('X', Column),
                                TRUE ~ Column))
  } else {
    colnames(dotplot_data)[c(2)] <- c("Group1")
    dotplot_data <- dotplot_data %>%
      mutate(Column = Group1,
             Column = case_when(grepl('^\\d', Column) ~ paste0('X', Column),
                                TRUE ~ Column))
  }
  
  
  dotplot <- dotplot_data %>%
    left_join(top_markers %>% 
                filter(Gene %in% (!!gene)) %>% 
                left_join(cat_to_color_df %>% 
                            filter(meta_category == 'CellType_predict') %>% 
                            select(cluster = value, gcolor = color)) %>% 
                select(Gene, gcolor, cluster)) %>% 
    mutate(Gene = gsub('\\(.*)','', Gene) %>% gsub(' ','', .)) %>% 
    left_join(cat_to_color_df %>% 
                filter(meta_category == 'CellType_predict') %>% 
                select(Group1 = value, color)) %>% 
    left_join(ct_order, by = c('Column' = 'CT')) %>% 
    left_join(ct_order %>% mutate(Order2 =rev(Order)) %>% select(CT, Order2), by = c('cluster' = 'CT')) %>% 
    mutate(Column = glue("<b style='color:{color}'>{Group1}</b>" %>% as.character), 
           Gene = glue("<b style='color:{gcolor}'>{Gene}</b>")) %>% 
    ungroup() %>% 
    mutate(cluster = fct_reorder(cluster, Order2),
           Column = fct_reorder(Column, Order)) %>% 
    # Gene = fct_reorder(Gene, as.factor(cluster))) %>%     
    ggplot(aes(x=Column, y = Gene, size = `%`, color = Expression)) +
    geom_point() + scale_radius(range=c(0, 10)) +
    cowplot::theme_cowplot() +
    scale_color_viridis_c(option = 'magma', guide = guide_legend(override.aes = list(size=10), label.position = "bottom",
                                                                 label.theme = element_text(angle = 90, hjust = 0.5, vjust = 0.5))) +
    scale_size_area(guide = guide_legend(label.position = "bottom",
                                         label.theme = element_text(angle = 90, hjust = 0.5, vjust = 0.5))) +
    theme(axis.line  = element_blank(),
          axis.ticks = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab('') + xlab('') +
    scale_y_discrete(position = 'right')+
    #coord_flip() + 
    theme(axis.text.y.right =  element_markdown(),
          axis.text.x = element_markdown(angle = 45, vjust = 1)) +
    theme(plot.margin=unit(c(0.5,0,0,0.5),"cm"))
  
  return(dotplot)
}