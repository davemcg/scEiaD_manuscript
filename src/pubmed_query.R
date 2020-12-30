# return PMID for query
library(easyPubMed)
library(tidyverse)
#load('~/data/massive_integrated_eye_scRNA/top_markers.Rdata')
library(pool)
library(RSQLite)
# scEiaD_2020_v01 <- dbPool(drv = SQLite(), dbname = "~/data/massive_integrated_eye_scRNA/MOARTABLES__anthology_limmaFALSE___5000-transform-counts-universe-batch-scVIprojectionSO-8-0.1-500-0.6.sqlite", idleTimeout = 3600000)

pubmed_counter <- function(query, sleep_time = 1){
  Sys.sleep(sleep_time)
  pub_query <- query
  entrez_id <- get_pubmed_ids(pub_query)
  if (entrez_id$Count == 0){
    out <- NA
  } else {
    abstracts_txt <- fetch_pubmed_data(entrez_id, format = "abstract") 
    out <- grep('PMID', abstracts_txt, value = TRUE) %>% str_extract(., '\\d\\d\\d+') %>% unique()}
  out
}

# x <- marker_info %>% 
#   left_join(exp_stats %>% dplyr::rename(cluster = CellType_predict)) %>% 
#   #filter((FDR < 1 |  med_auc > 0.20 | L1 > 50), D_KL > 0.15,`%` > 20) %>%
#   filter((FDR < 1 |  med_auc > 0.20 | L1 > 50), `%` > 20) %>% 
#   #group_by(cluster) %>% 
#   #top_n(100, med_auc) %>% 
#   ungroup() %>% 
#   group_by(Gene) %>% 
#   slice_max(order_by = `%`) %>% 
#   mutate(pm_query = paste0(gsub(' \\(.*','', Gene), ' AND ', cluster) %>% 
#            gsub('Cells','',.)) %>% 
#   mutate(pm_query2 = paste0(gsub(' \\(.*','', Gene), ' AND Retina'))

x <- top_markers %>% 
    mutate(pm_query = paste0(gsub(' \\(.*','', Gene), ' AND ', cluster) %>%
             gsub('Cells','',.)) %>%
    mutate(pm_query2 = paste0(gsub(' \\(.*','', Gene), ' AND Retina'))
pmid <- list()
for (i in x$pm_query){
  print(i)
  pmid[[i]] <- pubmed_counter(i)
}

pmid2 <- list()
for (i in x$pm_query2){
  print(i)
  #Sys.sleep(1)
  pmid2[[i]] <- pubmed_counter(i)
}


set.seed(1423)
all_genes <- scEiaD_2020_v01 %>% tbl('Genes') %>% pull(Gene)
rand_not_marker <- all_genes[!all_genes %in% x$Gene] %>% gsub(' \\(.*','', .)
rand_not_marker <- rand_not_marker[sample(1:length(rand_not_marker), 100)]
x2 <- rand_not_marker %>% enframe() %>% mutate(pm_query2 = paste0(gsub(' \\(.*','', value), ' AND Retina'))
pmid3 <- list()
for (i in x2$pm_query2){
  print(i)
  #Sys.sleep(1)
  pmid3[[i]] <- pubmed_counter(i)
}
# pmid2 <- list()
# for (i in x$pm_query2){
#   print(i)
#   #Sys.sleep(1)
#   pmid2[[i]] <- pubmed_counter(i)
# }

save(pmid, pmid2, pmid3, file = '~/data/massive_integrated_eye_scRNA/top_marker_pmid.Rdata')
top_markers$c1 <- pmid %>% map(function(x) sum(!is.na(x))) %>% unlist() %>% enframe() %>% pull(value)
top_markers$c2 <- pmid2 %>% map(function(x) sum(!is.na(x))) %>% unlist() %>% enframe() %>% pull(value)
save(top_markers, marker_info, exp_stats, file = '~/data/massive_integrated_eye_scRNA/top_markers.Rdata')
#pmid %>% map(function(x) sum(!is.na(x))) %>% unlist() %>% enframe()
#pmid %>% map(function(x) paste %>% unlist() %>% enframe()