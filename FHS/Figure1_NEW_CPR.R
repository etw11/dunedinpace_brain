################
# CPR, 2023-06-29
################

######
# Load packages
######
pacman::p_load(tidyverse, haven, sjlabelled, readr)

######
# Load data
######
whitman_new <- read_dta("Data/WhitmanFhamBrainAnalysis.dta")

# Clean
whitman_new_short <-whitman_new %>% 
  select(age8, dunedinpace, Left_Lateral_Ventricle) %>% 
  na.omit()
# 948 (OG had 903)

fhs <-whitman_new_short %>% 
  select(age8) %>%
  rename(Age = age8) %>%  
  mutate(dataset = 'FHS')


# Load ADNI
load("Data/adni_age.RData")

adni_age <-as.data.frame(adni_age)
names(adni_age) <-"Age"

# Clean
adni <-adni_age %>%
  mutate(Age = Age/12) %>% 
  mutate(dataset = "ADNI")

# Bind
full_fhs_adni <-bind_rows(fhs, as_tibble(adni))

###########
# Fig 1
############ 
figure1 <-ggplot(data = full_fhs_adni, aes(x=Age, group = dataset,  fill = dataset)) +
  geom_density(alpha = .5) +
  geom_vline(xintercept = 45, color = '#F8766D', linetype = 'dashed', size = 1.5) +
  # feel free to adjust if any participants are outside of this age range
  xlim(40, 100)+
  labs(x = 'Age in years',
       y = 'Density',
       title = "Dataset Age Distributions") +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15, vjust = 1),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 20))+
  scale_fill_manual(values = c("#619CFF", "#00BA38"))


save(figure1, file = here::here("Output/Figures", "Framingham_Figure1Densities_NewData.RData"))

