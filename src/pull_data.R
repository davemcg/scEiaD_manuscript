# Extract data for manuscript figures

library(pool)
library(RSQLite)
library(tidyverse)

scEiaD <- dbPool(drv = SQLite(), dbname = "~/data/massive_integrated_eye_scRNA/MOARTABLES__anthology_limmaFALSE___Mus_musculus_Macaca_fascicularis_Homo_sapiens-5000-counts-TabulaDroplet-batch-scVI-8-0.1-15-7.sqlite", idleTimeout = 3600000)

meta_filter <- scEiaD %>% tbl('metadata_filter')
