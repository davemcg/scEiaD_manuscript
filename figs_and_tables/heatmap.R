genes <- c('ENSG00000175130','ENSG00000104435','ENSG00000067715','ENSG00000144834','ENSG00000139970','ENSG00000109846','ENSG00000010278','ENSG00000172216','ENSG00000026025','ENSG00000177606','ENSG00000067715','ENSG00000123307','ENSG00000154975','ENSG00000188133','ENSG00000129244','ENSG00000139053','ENSG00000143774','ENSG00000021300','ENSG00000090097','ENSG00000112041','ENSG00000087086','ENSG00000105372','ENSG00000010278','ENSG00000026025','ENSG00000130066','ENSG00000104435','ENSG00000137285','ENSG00000119547','ENSG00000144834','ENSG00000139970','ENSG00000115310','ENSG00000050165','ENSG00000140522','ENSG00000162998','ENSG00000135821','ENSG00000104435','ENSG00000137285','ENSG00000131711','ENSG00000172020','ENSG00000139970','ENSG00000016082','ENSG00000253368','ENSG00000198668','ENSG00000183036','ENSG00000084234','ENSG00000149489','ENSG00000129535','ENSG00000185527','ENSG00000112619','ENSG00000130561','ENSG00000105372','ENSG00000231500','ENSG00000177600','ENSG00000083845','ENSG00000140988')

markers <- scEiaD_2020_v01 %>% tbl('genes') %>% collect() %>% mutate(ENS = str_extract(Gene, 'ENSG\\d+')) %>% filter(ENS %in% genes) %>% pull(Gene)
dat <- scEiaD_2020_v01 %>% tbl('cpm') %>%
  filter(Gene %in% markers) %>% collect()

dat_w <- dat %>% pivot_wider(values_from = cpm, names_from = Barcode)

# only keep main retina cell types
dat_w <- dat_w[,colnames(dat_w) %in% 
                 (meta_filter %>% 
                    filter(grepl('Amacrine|Rod|Cone|Retinal|Muller|Horizon|Bipol|RPE|Endo|Astro', CellType_predict)) %>% 
                    pull(Barcode)) ]

dat_m <- meta_filter$Barcode[!(meta_filter %>% 
                                 filter(grepl('Amacrine|Rod|Cone|Retinal|Muller|Horizon|Bipol|RPE|Endo|Astro', CellType_predict)) %>% 
                                 pull(Barcode)) %in%
                               dat$Barcode ]

mat_empty <- matrix(0L, nrow = nrow(dat_w), ncol = length(dat_m))
colnames(mat_empty) <- dat_m




dat_all <- cbind(dat_w, mat_empty)
sample_down <- dat_all[,sample(ncol(dat_all), 1e4)]
sample_down[is.na(sample_down)] <- 0

labels <- colnames(sample_down) %>% enframe(value = 'Barcode') %>% left_join(meta_filter %>% select(Barcode, CellType_predict)) %>% pull(CellType_predict)
ha_column = HeatmapAnnotation(df = data.frame(CT  = labels), 
                              col = list(CT = c("AC/HC Precursors" = "#F0A0FF" ,
                                                "Amacrine Cells" = "#0075DC" ,
                                                "Artery" = "#993F00" ,
                                                "Astrocytes" = "#4C005C" ,
                                                "B-Cell" = "#191919" ,
                                                "Bipolar Cells" = "#005C31" ,
                                                "Choriocapillaris" = "#2BCE48" ,
                                                "Cones" = "#FFCC99" ,
                                                "Early RPCs" = "#808080" ,
                                                "Endothelial" = "#94FFB5" ,
                                                "Fibroblasts" = "#8F7C00" ,
                                                "Horizontal Cells" = "#9DCC00" ,
                                                "Late RPCs" = "#C20088" ,
                                                "Macrophage" = "#003380" ,
                                                "Melanocytes" = "#FFA405" ,
                                                "Microglia" = "#FFA8BB" ,
                                                "Muller Glia" = "#426600" ,
                                                "Neurogenic Cells" = "#FF0010" ,
                                                "Pericytes" = "#5EF1F2" ,
                                                "PR Precursors" = "#00998F" ,
                                                "Red Blood Cells" = "#E0FF66" ,
                                                "Retinal Ganglion Cells" = "#740AFF" ,
                                                "Rod Bipolar Cells" = "#990000" ,
                                                "Rods" = "#FFFF80" ,
                                                "RPCs" = "#FFE100" ,
                                                "RPE" = "#FF5005" ,
                                                "Schwann" = "#AA0DFE" ,
                                                "Smooth Muscle Cell" = "#3283FE" ,
                                                "T-Cell" = "#85660D" ,
                                                "Vein" = "#782AB6" )))

#colnames(sample_down) <- NA
Heatmap((as.matrix(sample_down)), 
        show_column_dend = FALSE, 
        show_row_dend = FALSE, 
        show_row_names = FALSE,
        show_column_names = FALSE,
        jitter = TRUE,
        col=viridis::viridis(100),
        top_annotation = ha_column)
#Heatmap(sample_down)



labels <- colnames(sample_down) %>% enframe(value = 'Barcode') %>% left_join(meta_filter %>% select(Barcode, CellType_predict)) %>% pull(CellType_predict)
ha_column = HeatmapAnnotation(df = data.frame(CT  = labels))
Line = meta  %>% 
  filter(col %in% colnames(sample_dists)) %>% 
  dplyr::select(col, Line) %>% 
  mutate(sample=factor(col, levels=colnames(sample_dists))) %>% 
  arrange(sample) %>% 
  unique() %>% 
  pull(Line)), 
col = list(Line = c('W0255-132' = pals::kelly(n=22)[6], 
                    'W0255-142' = pals::kelly(n=22)[10],
                    'W0118-151' = pals::kelly(n=22)[11],
                    'W0118-254' = pals::kelly(n=22)[12],
                    'D3C' = pals::kelly(n=22)[15]),
           Info = c("Day 2 iPS" =  viridis::magma(20)[3], "Day 11 RPE" = viridis::magma(20)[8], "Day 25 RPE" = viridis::magma(20)[13])))

Heatmap(sample_dists %>% as.matrix(), 
        col=viridis::viridis(100),
        name = 'Sample\nDistances',
        top_annotation = ha_column)