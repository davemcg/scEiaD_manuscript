library(ggplot2)
library(Cairo)
library(scattermore)
library(pals)
library(ggrepel)
library(patchwork)
library(purrr)
library(pool)
library(RSQLite)
library(dplyr)
library(magick)
library(stringr)
library(fst)

tabulamuris_predict_labels <-scEiaD_2020_v01 %>% tbl('tabulamuris_predict_labels') %>%  collect
celltype_predict_labels <-scEiaD_2020_v01 %>% tbl('celltype_predict_labels') %>% 
  mutate(CellType_predict = case_when(CellType_predict == 'RPC' ~ 'RPCs',
                                      CellType_predict == 'Mesenchymal/RPE/Endothelial' ~ 'Endothelial',
                                      TRUE ~ CellType_predict)) %>% collect() %>% 
  mutate( UMAP_a = UMAP_2 * x_dir,
          UMAP_b = UMAP_1 * y_dir) %>%
  mutate(UMAP_1 = UMAP_a, UMAP_2 = UMAP_b)
celltype_labels <-scEiaD_2020_v01 %>% tbl('celltype_labels') %>% collect() %>% 
  mutate(  UMAP_a = UMAP_2 * x_dir,
           UMAP_b = UMAP_1 * y_dir) %>%
  mutate(UMAP_1 = UMAP_a, UMAP_2 = UMAP_b)
cluster_labels <-scEiaD_2020_v01 %>% tbl('cluster_labels') %>% collect() %>% 
  mutate(  UMAP_a = UMAP_2 * x_dir,
           UMAP_b = UMAP_1 * y_dir) %>%
  mutate(UMAP_1 = UMAP_a, UMAP_2 = UMAP_b)
mf <- meta_filter 

# generate color_mappings
categorical_columns <- c("Phase","batch","study_accession","library_layout","organism","Platform",
                         "Covariate","CellType","CellType_predict","TabulaMurisCellType","TabulaMurisCellType_predict",
                         "GSE","Summary","Design","Citation","PMID","Stage","cluster",
                         "Doublet","TechType", "SubCellType", 'subcluster' )
#"SubCellType" and subcluster are problems
meta_filter <- meta_filter %>% mutate(SubCellType = tidyr::replace_na(SubCellType, 'None'),
                                      subcluster = as.character(subcluster)) 


map_color <- function(column, meta_filter){
  # master_colorlist <- c(pals::polychrome()[3:length(pals::polychrome())], pals::alphabet2())
  master_colorlist <- c(pals::cols25()[1:23],pals::alphabet())
  values <- meta_filter %>% pull(!!column) %>% unique %>% sort
  if(length(values) > length(master_colorlist) ){
    r= round(length(values) / length(master_colorlist)) +1
    master_colorlist <- rep(master_colorlist, r)
  }
  
  colors <- master_colorlist[1:length(values)]
  return(tibble(meta_category = column,value = values, color=colors))
  
}


cat_to_color_df <- lapply(categorical_columns, function(col) map_color(col, meta_filter)) %>% bind_rows()
# %>%
#   bind_rows(sub_cluster_map)

make_meta_scatter_umap_plot <- function(input, mf, meta_filter,
                                        celltype_predict_labels,
                                        celltype_labels,
                                        tabulamuris_predict_labels,
                                        cluster_labels,
                                        cat_to_color_df,
                                        velocity = FALSE
){
  
  meta_column <- input$meta_column
  transform <- input$meta_column_transform
  
  pt_size <- input$pt_size_meta %>% as.numeric()
  filter_column <- input$meta_column
  # cut down to match tech selected
  #tech <- input$gene_and_meta_scatter_tech
  #mf <- mf %>% filter(TechType == tech)
  #meta_filter <- meta_filter %>% filter(TechType == tech)
  #celltype_predict_labels <- celltype_predict_labels %>% filter(TechType == tech)
  #celltype_labels <- celltype_labels %>% filter(TechType == tech)
  #tabulamuris_predict_labels <- tabulamuris_predict_labels %>% filter(TechType == tech)
  #cluster_labels <- cluster_labels %>% filter(TechType == tech)
  
  if (transform == 'log2' && is.numeric(meta_filter[,meta_column] %>% pull(1))){
    cat('log2 time')
    meta_filter[,meta_column] <- log2(meta_filter[,meta_column] + 1)
  }
  if (!is.null(input$meta_filter_on)){
    p_data <- meta_filter %>%
      #filter(!grepl('Doub|\\/Margin\\/Periocular', CellType)) %>%
      filter_at(vars(all_of(input$meta_filter_cat)), all_vars(. %in% input$meta_filter_on))
  } else {
    p_data <- meta_filter %>% filter(!is.na(!!as.symbol(input$meta_column)))
  }
  
  # metadata NUMERIC plot --------------
  if (is.numeric(meta_filter[,meta_column] %>% pull(1)) ){
    color_range <- range(p_data[,meta_column] %>% pull(1))
    suppressWarnings(plot <- ggplot() +
                       geom_scattermost(cbind(mf %>% pull(UMAP_1),
                                              mf %>% pull(UMAP_2)),
                                        pointsize = pt_size - 1, color = '#D3D3D333',
                                        pixels = c(1750,750)) +
                       geom_scattermost(cbind(p_data$UMAP_1, p_data$UMAP_2),
                                        color = viridis::viridis(100, alpha=0.3)
                                        [1+99*((p_data[,meta_column] %>% pull(1))-color_range[1])/diff(color_range)],
                                        pointsize= pt_size - 1,
                                        pixels=c(1750,750),
                                        interpolate=FALSE) +
                       geom_point(data=data.frame(x=double(0)), aes(x,x,color=x))  +
                       scale_color_gradientn(  #add the manual guide for the empty aes
                         limits=c(min(p_data[,meta_column] %>% pull(1)),
                                  max(p_data[,meta_column] %>% pull(1))),
                         colors=viridis::viridis(100),
                         name=meta_column %>% gsub('_',' (', .) %>% gsub('$',')',.)) +
                       theme_cowplot() +
                       guides(colour = guide_legend(override.aes = list(size=7, alpha = 1))) +
                       theme(axis.line = element_blank(),
                             axis.title = element_blank(),
                             axis.ticks = element_blank(),
                             axis.text = element_blank()) +
                       annotate("text", -Inf, Inf, label = "Metadata", hjust = 0, vjust = 1, size = 6))
    # metadata CATEGORICAL plot --------------
  } else {
    cur_color_df <- cat_to_color_df %>%
      filter( meta_category %in%  meta_column,
              value %in% p_data[[meta_column]]) %>% distinct
    
    color_list <- cur_color_df %>% pull(color)
    # replaced join with this for speed
    k <- cur_color_df$color
    names(k) <- cur_color_df$value
    np_color <- {k[p_data[[meta_column]] ]} %>% paste0(., '33')
    names(np_color) <- NULL
    color_data <- cur_color_df  %>%
      select(value) %>%
      mutate( x=0)
    names(color_list) <- color_data$value
    
    suppressWarnings(plot <- ggplot() +
                       geom_scattermost(cbind(mf  %>% pull(UMAP_1),
                                              mf %>% pull(UMAP_2)),
                                        pointsize = input$pt_size_back, color = '#D3D3D333',
                                        pixels = c(1750,1750)) +
                       geom_scattermore(data = p_data, alpha = 0.2,
                                        aes(x= UMAP_1, y=UMAP_2),
                                        color = np_color,
                                        pointsize= pt_size,
                                        pixels=c(1750,1750), 
                                        interpolate=FALSE) +
                       #geom_point(data=data.frame(x=double(0)), aes(x,x,color=x))  +
                       geom_point(data=color_data, aes(x,x,color=value), alpha = 0) +
                       scale_colour_manual(name= meta_column %>% gsub('_',' (', .) %>% gsub('$',')',.),
                                           values = color_list) +
                       theme_cowplot() +
                       theme(axis.line = element_blank(),
                             axis.title = element_blank(),
                             axis.ticks = element_blank(),
                             axis.text = element_blank()) +
                       guides(colour = guide_legend(override.aes = list(alpha = 1, size = 7)))
    )
    
  }
  
  
  more <- NULL
  if ('1' %in% input$label_toggle){
    more <- geom_text_repel(data = celltype_labels, bg.color = 'white',
                            aes(x = UMAP_1, y = UMAP_2, label = CellType))
  }
  if ('2' %in% input$label_toggle){
    more <- geom_text_repel(data = celltype_predict_labels, bg.color = 'white',
                            aes(x = UMAP_1, y = UMAP_2, label = CellType_predict))
  }
  if ('3' %in% input$label_toggle){
    more <- geom_text_repel(data = cluster_labels, bg.color = 'white',
                            aes(x = UMAP_1, y = UMAP_2, label = cluster),
                            max.iter = 20)
  }
  if ('4' %in% input$label_toggle){
    more <- geom_text_repel(data = tabulamuris_predict_labels, bg.color = 'white',
                            aes(x = UMAP_1, y = UMAP_2, label = TabulaMurisCellType_predict),
                            max.iter = 20)
  }
  col_size <- {meta_filter[[meta_column]]} %>% n_distinct
  out <- list()
  out$plot <- plot + more
  if (velocity){
    out$plot <- plot + geom_segment(data=grid.df, alpha = 0.6, mapping=aes(y=start.1, x=-start.2, 
                                                             yend=end.1, xend=-end.2), arrow=arrow(length=unit(0.05, "inches"))) + more
  }
  out$col_size <- col_size
  out
}


celltype_predict_labels <- celltype_predict_labels %>% 
  mutate(CellType_predict = 
           case_when(CellType_predict %in% 
                       c('RPE', 'Endothelial','Choriocapillaris','Artery', 'Vein', 'Melanocytes') ~ 'Connective', 
                     TRUE ~ CellType_predict)) %>% 
  group_by(CellType_predict) %>% 
  summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2)) %>% 
  filter(!CellType_predict %in% c('B-Cell','T-Cell','Smooth Muscle Cell', 'Red Blood Cells')) %>% 
  mutate(CellType_predict = case_when(CellType_predict == 'Photoreceptor Precursors' ~ 'PR Precursors',
                                      CellType_predict == 'AC/HC_Precurs' ~ 'AC/HC Precursors',
                                      CellType_predict == 'RPC' ~ 'RPCs',
                                      CellType_predict == 'Mesenchymal/RPE/Endothelial' ~ 'Endothelial',
                                      TRUE ~ CellType_predict))

