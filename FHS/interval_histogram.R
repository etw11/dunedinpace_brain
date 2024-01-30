################
# CPR, 2023-06-29
################

######
# Load packages
######
pacman::p_load(tidyverse, haven, sjlabelled, readr)

######
# Load Data
######
whitman_new <- read_csv(here::here("Data", "DaysfromDNAmtoMRIv32data.csv")) %>% 
  mutate(day_difference = mri_date - date8)


library(ggplot2)

# i am plotting the difference in days between blood draw used for DNA methylation and the MRI scan. Specifically, 
# I am subtracting DNA_meth_date - MRI_date, so negative days represents an MRI scan before blood draw and positive
# means MRI scan after blood draw

fhs_scanint_hist <- ggplot(data = whitman_new, aes(x = day_difference)) + 
  geom_histogram(fill='#00BA38') + # switch color to 
  labs(title = "FHS-OC", #switch to FHS-OC, of course
       x = "Days between blood draw and scan",
       y = 'Frequency') +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20, vjust = 1),
        axis.title.x = element_text(size = 20),
        title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.position = 'none')

# change as needed to save as an .Rdata object
save(fhs_scanint_hist, file = here::here("Output", "fhs_scanint_hist_mri_minus_dnam.Rdata"))
