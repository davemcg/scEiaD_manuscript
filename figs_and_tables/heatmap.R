dar <- scEiaD_2020_v01 %>% tbl('cpm') %>%
  filter(Gene %in% c('CRX (ENSG00000105392)', 'TFAP2A (ENSG00000137203)')) %>% collect()

dar_w <- dar %>% pivot_wider(values_from = cpm, names_from = Barcode)

dar_m <- meta_filter$Barcode[!meta_filter$Barcode %in%  dar$Barcode ]

mat_empty <- matrix(0L, nrow = nrow(dar_w), ncol = length(dar_m))
colnames(mat_empty) <- dar_m


sample_down <- dar_all[,sample(ncol(dar_all), 1e4)]
sample_down[is.na(sample_down)] <- 0
