library(scattermore)
source('figs_and_tables/fig2__integrationPerf.R')
load('~/data/massive_integrated_eye_scRNA/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_spec_genes-0__n_features2000__libSize__TabulaDroplet__batch__none__dims30__preFilter__mindist0.3__nneighbors30.umap.Rdata')
umapNone <- umap

load('~/data/massive_integrated_eye_scRNA/Mus_musculus_Macaca_fascicularis_Homo_sapiens__n_spec_genes-0__n_features2000__libSize__TabulaDroplet__batch__combat__dims30__preFilter__mindist0.3__nneighbors30.umap.Rdata')
umapCombat <- umap

umap_plot_maker <- function(umap, color_against = 'CellType', red = 'UMAP', ptsize = 2){
  celltype_col <- 'CellType' 
  # filter
  umapFig <- umap %>% 
    #rename(Stage = integration_group) %>% 
    mutate(CellType = gsub('Rod Bipolar Cells', 'Bipolar Cells', !!as.symbol(celltype_col))) %>% 
    filter(!is.na(CellType), 
           !is.na(study_accession), 
           !CellType %in% c('Doublet', 'Doublets'),
           !grepl('Mesenchyme/Lens', CellType))  %>% 
    mutate(Size = case_when(organism == 'Homo sapiens' ~ 0.015,
                            TRUE ~ 0.01)) 
  
  # attach colors to cell types
  cell_types <- umapFig %>% 
    pull(CellType) %>% unique() %>% sort()
  type_val <- setNames(c(pals::alphabet(), pals::alphabet2())[1:length(cell_types)], cell_types)
  type_col <- scale_colour_manual(values = type_val)
  type_fill <- scale_fill_manual(values = type_val)
  # cell type known
  if (color_against == 'CellType'){
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
  }
}

