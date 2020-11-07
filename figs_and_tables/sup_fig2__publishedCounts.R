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

median <- nxn_sc %>% 
  mutate(Date = as.Date(x = as.character(Date), format = '%Y%m%d'), 
         Year = format(Date, '%Y')) %>% 
  group_by(Year) %>% 
  summarise(Count = sum(`Reported cells total`, na.rm = TRUE), Median = median(`Reported cells total`, na.rm = TRUE)) %>% 
  ggplot(aes(x=Year, y = Median)) +
  geom_bar(stat = 'identity') +
  ylab("Median of SC Counts") +
  cowplot::theme_cowplot() + 
  coord_flip()

mean <- nxn_sc %>% 
  mutate(Date = as.Date(x = as.character(Date), format = '%Y%m%d'), 
         Year = format(Date, '%Y')) %>% 
  group_by(Year) %>% 
  summarise(Count = sum(`Reported cells total`, na.rm = TRUE), Mean = mean(`Reported cells total`, na.rm = TRUE)) %>% 
  ggplot(aes(x=Year, y = Mean)) +
  geom_bar(stat = 'identity') +
  ylab("Mean of SC Counts / Dataset") +
  cowplot::theme_cowplot() + 
  coord_flip() +
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))

