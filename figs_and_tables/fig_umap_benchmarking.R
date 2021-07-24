library(scattermore)
source('~/git/scEiaD_manuscript/figs_and_tables/fig2__integrationPerf.R')
load('~/data/scEiaD/giga_all_methods/gigascience/umap/n_features-2000__transform-standard__partition-universe__covariate-batch__method-none__dims-30__epochs-5__dist-0.1__neighbors-50__knn-20__umap.Rdata')
umapNone <- umap

load('~/data/scEiaD/giga_all_methods/gigascience/umap/n_features-2000__transform-libSize__partition-universe__covariate-batch__method-combat__dims-8__epochs-5__dist-0.1__neighbors-50__knn-20__umap.Rdata')
umapCombat <- umap

umap_plot_maker <- function(umap, color_against = 'CellType', red = 'UMAP', ptsize = 2){
  celltype_col <- 'CellType' 
  if (color_against == 'CellType_predict' ){
    umap$CellType = umap$CellType_predict
  }
  # filter
  umapFig <- umap %>% 
    #rename(Stage = integration_group) %>% 
    mutate(CellType = gsub('Rod Bipolar Cells', 'Bipolar Cells', !!as.symbol(celltype_col))) %>% 
    
    filter(!CellType %in% c('RPE', 'Vein', 'Smooth Muscle Cell', 'Schwann','Macrophage','Mast')) %>% 
    filter(!is.na(CellType), 
           !is.na(study_accession), 
           !CellType %in% c('Doublet', 'Doublets'),
           !grepl('Mesenchyme/Lens', CellType))  %>% 
    mutate(Size = case_when(organism == 'Homo sapiens' ~ 0.015,
                            TRUE ~ 0.01)) 
  
  # attach colors to cell types
  cell_types <- umapFig %>% 
    pull(CellType) %>% unique() %>% sort()
  type_val <- setNames(c(pals::cols25(),pals::alphabet())[1:length(cell_types)], cell_types)
  type_col <- scale_colour_manual(values = type_val)
  type_fill <- scale_fill_manual(values = type_val)
  # cell type known
  if (color_against == 'CellType' | color_against == 'CellType_predict'){
    umapFig %>% 
      ggplot() + 
      geom_scattermore(aes(x=umapFig[,paste0(red,'_1')] %>% pull(1), 
                           y = umapFig[,paste0(red,'_2')] %>% pull(1), 
                           colour = CellType), pointsize = (ptsize/3), alpha = 0.1) + 
      guides(colour = guide_legend(override.aes = list(size=2, alpha = 1))) + 
      theme_void() +
      type_col + 
      xlab(paste(red, '1')) + ylab(paste(red, '2')) 
  } else if (color_against == 'study_accession'){
    umapFig %>% 
      ggplot() + 
      geom_scattermore(aes(x=umapFig[,paste0(red,'_1')] %>% pull(1), 
                           y = umapFig[,paste0(red,'_2')] %>% pull(1), 
                           colour = study_accession), pointsize = (ptsize/3), alpha = 0.1) + 
      guides(colour = guide_legend(override.aes = list(size=2, alpha = 1))) + 
      theme_void() +
      scale_colour_manual(values = pals::glasbey()) +
      xlab(paste(red, '1')) + ylab(paste(red, '2')) 
  } else if (color_against == 'organism'){
    umapFig %>% 
      ggplot() + 
      geom_scattermore(aes(x=umapFig[,paste0(red,'_1')] %>% pull(1), 
                           y = umapFig[,paste0(red,'_2')] %>% pull(1), 
                           colour = organism), pointsize = (ptsize/3), alpha = 0.1) + 
      guides(colour = guide_legend(override.aes = list(size=2, alpha = 1))) + 
      theme_void() +
      scale_colour_manual(values = pals::brewer.set3(n=3)) +
      xlab(paste(red, '1')) + ylab(paste(red, '2')) 
  }
}

