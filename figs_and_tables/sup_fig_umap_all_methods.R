umap_files <- list.files('~/data/scEiaD/giga_all_methods/gigascience/umap/', full.names = TRUE)

umap_files <- umap_files %>% enframe() %>%  mutate(method = str_extract(value, 'method-[a-zA-Z]+') %>% gsub('method-','',.)) %>% arrange(method) %>% pull(value)

umap_data <- list()
for (i in umap_files){
  print(i)
  dims = str_extract(i, 'dims-\\d+') %>% gsub('dims-','',.) %>% as.integer()
  norm = str_extract(i, 'transform-[a-zA-Z]+') %>% gsub('transform-','',.)
  load(i)
  umap$dims <- dims
  umap$norm <- norm
  umap_data[[i]] <- umap
}


# plot by celltype
plots_celltype <- list()
for (i in umap_files){
  plots_celltype[[i]] <- umap_plot_maker(umap_data[[i]], color_against <- 'CellType', red = 'UMAP') + 
    theme(legend.position = 'none') +
    ggtitle(paste(umap_data[[i]]$Method %>% unique(), 'dims:', umap_data[[i]]$dims%>% unique(), 'normalization:', umap_data[[i]]$norm %>% unique()))
}
legend <- get_legend(
  umap_plot_maker(umap_data[[1]], color_against <- 'CellType', red = 'UMAP') + 
    ggtitle(paste(umap_data[[1]]$Method %>% unique(), 'dims:', umap_data[[1]]$dims%>% unique(), 'normalization:', umap_data[[1]]$norm %>% unique())) +
    theme(legend.box.margin = margin(0, 0, 0, 12))
)

pdf('all_methods_celltype.pdf', width = 12, height = 50)
plot_grid(plot_grid(plotlist = plots_celltype, ncol = 3), plot_grid(legend), ncol = 1, rel_heights = c(20,4))
dev.off()

# plot by study accession
plots_batch <- list()
for (i in umap_files){
  plots_batch[[i]] <- umap_plot_maker(umap_data[[i]], color_against <- 'study_accession', red = 'UMAP') + 
    theme(legend.position = 'none') +
    ggtitle(paste(umap_data[[i]]$Method %>% unique(), 'dims:', umap_data[[i]]$dims%>% unique(), 'normalization:', umap_data[[i]]$norm %>% unique()))
}
legend_batch <- get_legend(
  umap_plot_maker(umap_data[[1]], color_against <- 'study_accession', red = 'UMAP') + 
    ggtitle(paste(umap_data[[1]]$Method %>% unique(), 'dims:', umap_data[[1]]$dims%>% unique(), 'normalization:', umap_data[[1]]$norm %>% unique())) +
    theme(legend.box.margin = margin(0, 0, 0, 12))
)
pdf('all_methods_batch.pdf', width = 12, height = 50)
plot_grid(plot_grid(plotlist = plots_batch, ncol = 3), plot_grid(legend_batch), ncol = 1, rel_heights = c(20,4))
dev.off()

# plot by organism
plots_organism <- list()
for (i in umap_files){
  plots_organism[[i]] <- umap_plot_maker(umap_data[[i]], color_against <- 'organism', red = 'UMAP') + 
    theme(legend.position = 'none') +
    ggtitle(paste(umap_data[[i]]$Method %>% unique(), 'dims:', umap_data[[i]]$dims%>% unique(), 'normalization:', umap_data[[i]]$norm %>% unique()))
}
legend_org <- get_legend(
  umap_plot_maker(umap_data[[1]], color_against <- 'organism', red = 'UMAP') + 
    ggtitle(paste(umap_data[[1]]$Method %>% unique(), 'dims:', umap_data[[1]]$dims%>% unique(), 'normalization:', umap_data[[1]]$norm %>% unique())) +
    theme(legend.box.margin = margin(0, 0, 0, 12))
)
pdf('all_methods_organism.pdf', width = 12, height = 50)
plot_grid(plot_grid(plotlist = plots_organism, ncol = 3), plot_grid(legend_org), ncol = 1, rel_heights = c(20,4))
dev.off()

