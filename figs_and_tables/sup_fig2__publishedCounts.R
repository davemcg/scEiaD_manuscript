nxn_sc <- read_tsv('http://www.nxn.se/single-cell-studies/data.tsv')
sum <- nxn_sc %>% 
  mutate(Date = as.Date(x = as.character(Date), format = '%Y%m%d'), 
         Year = format(Date, '%Y')) %>% 
  group_by(Year) %>% 
  summarise(Count = sum(`Reported cells total`, na.rm = TRUE)) %>% 
  mutate(Sum = cumsum(Count)) %>% 
  ggplot(aes(x=Year, y = Sum)) +
  geom_bar(stat = 'identity') +
  ylab("Sum of SC Counts") +
  cowplot::theme_cowplot() + 
  coord_flip()
cowplot::plot_grid(sum, sum + scale_y_log10() + ylab('Log(Sum of SC Counts)'))


