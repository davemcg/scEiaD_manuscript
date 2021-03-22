library(ggalluvial)
counts <- c(qc %>% 
              mutate(organism = case_when(is.na(organism) ~ 'Homo sapiens',
                                          TRUE ~ organism)) %>% 
              mutate(QC = gsub(' QC','', QC),
                     QC = gsub('Uni', '\nUni', QC), 
                     QC = gsub('Doubl', '\nDoubl', QC), 
                     TechType = case_when(Platform %in% c('10xv2','10xv3','DropSeq') ~ 'Droplet',
                                          TRUE ~ 'Well')) %>% 
              mutate(OrgTech = paste0(organism, '\n', TechType)) %>% 
              mutate(OrgTech = gsub('fasci','\nfasci', OrgTech)) %>% 
              group_by(OrgTech) %>% 
              summarise(Count = n()) %>% pull(Count) %>% rev(),
            qc %>% 
              mutate(QC = gsub(' QC','', QC),
                     QC = gsub('Uni', '\nUni', QC), 
                     QC = gsub('Doubl', '\nDoubl', QC), 
                     TechType = case_when(Platform %in% c('10xv2','10xv3','DropSeq') ~ 'Droplet',
                                          TRUE ~ 'Well')) %>% 
              mutate(OrgTech = paste0(organism, '\n', TechType)) %>% 
              mutate(OrgTech = gsub('fasci','\nfasci', OrgTech)) %>% 
              group_by(QC) %>% 
              summarise(Count = n()) %>% pull(Count) %>% rev()
) 

supFig1_plot <- qc %>% 
  mutate(organism = case_when(is.na(organism) ~ 'Homo sapiens',
                                                   TRUE ~ organism)) %>% 
  mutate(QC = gsub(' QC','', QC),
         QC = gsub('Uni', '\nUni', QC), 
         QC = gsub('Doubl', '\nDoubl', QC), 
         TechType = case_when(Platform %in% c('10xv2','10xv3','DropSeq') ~ 'Droplet',
                              TRUE ~ 'Well')) %>% 
  mutate(OrgTech = paste0(organism, '\n', TechType)) %>% 
  mutate(OrgTech = gsub('fasci','\nfasci', OrgTech)) %>% 
  group_by(OrgTech, QC) %>% 
  summarise(Count = n()) %>% 
  #ggplot(aes(y = sqrt(Count), axis1 = organism, axis2 = TechType, axis3 = QC)) +
  ggplot(aes(y = sqrt(Count), axis1 = OrgTech, axis2 = QC, label = Count)) +
  geom_alluvium(aes(fill = `QC`), alpha = 0.4) +
  geom_stratum(aes(fill = `OrgTech`), alpha = 0.6) +
  geom_stratum(aes(fill = `QC`), alpha = 0.6) + 
  ggrepel::geom_text_repel(stat = "stratum", segment.color = NA,
                           nudge_x = c(rep(c(-20),5), rep(c(20), 4)),
                           direction = "y", size = 4,
                           aes(label = after_stat(stratum))) + theme_void() +
  shadowtext::geom_shadowtext(stat = "stratum", label=counts, bg.color = 'white', color = 'black') +
  scale_fill_manual(values = c(pals::glasbey() %>% unname())) +
  coord_cartesian(xlim = c(0,2.8)) +
  theme(legend.position = "none")
