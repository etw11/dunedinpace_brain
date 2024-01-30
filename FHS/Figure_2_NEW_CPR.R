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
whitman_new <- read_dta("Data/WhitmanFhamBrainAnalysis.dta")

whitman_new_clean <-whitman_new %>%
  filter(!is.na(age8) & !is.na(dunedinpace) & !is.na(Left_Lateral_Ventricle)) 
# 948 vs. 903 in OG

######
# Check variables
######
# Total Brain Volumn (TBV):
whitman_new_clean %>% select(contains("TBV"))
# Hippocampal volume (HC):
whitman_new_clean %>% select(contains("HV")) 
# White matter hypointensities (WMH):
whitman_new_clean %>% select(contains("WMH"))  
# ** this variable should be log-transformed
# Cortical thickness (CT):
whitman_new_clean %>% select(contains("MCT"))  
# Cortical surface area (SA):
whitman_new_clean %>% select(contains("TSA"))

######
# Scale and residulize on age and sex
######
tbv <-lm(scale(TBV) ~ sex + age8 + age8*sex + age8^2 + age8^2*sex + ICV, data = whitman_new_clean)
hv <-update(tbv, scale(HV) ~ .)
wmh <-update(tbv, scale(log(WMH)) ~ . - ICV)
mct <-update(tbv, scale(MCT) ~ . - ICV)
tsa <-update(tbv, scale(TSA) ~ . - ICV)

whitman_new_clean$tbv_resid <-resid(tbv)
whitman_new_clean$hv_resid <-resid(hv)
whitman_new_clean$wmh_resid <-resid(wmh)
whitman_new_clean$mct_resid <-resid(mct)
whitman_new_clean$tsa_resid <-resid(tsa)


whitman_new_clean <-whitman_new_clean %>% select(dunedinpace, contains('resid'))
######
# Create figures
######

# Function to create scatter plot with different y-axis variables
create_scatter_plot <- function(y_variable) {
  # Map y_variable to y_axis_label
  y_axis_label <- switch(
    y_variable,
    "tbv_resid" = "TBV",
    "hv_resid" = "HV",
    "wmh_resid" = "WMH",
    "mct_resid" = "MCT",
    "tsa_resid" = "TSA"
  )
  
  ggplot(data = whitman_new_clean, aes(x = scale(dunedinpace), y = !!sym(y_variable))) +
    geom_point(color = '#00BA38') +
    geom_smooth(method = 'lm', color = 'black', se = FALSE) +
    ylim(-3, 2.4) +
    xlim(-3.4, 5.1) +
    labs(y = y_axis_label,
         x = "DunedinPACE",
         title = 'Dunedin') +
    theme_classic() +
    theme(
      axis.text.x = element_text(vjust = 0.5, size = 10),
      axis.text.y = element_text(vjust = 0.5, size = 10),
      axis.title.y = element_text(size = 15, vjust = 0),
      axis.title.x = element_text(size = 15),
      title = element_text(size = 0)
    )
}

# List of y-axis variables
y_variables <- c("tbv_resid", "hv_resid", "wmh_resid", "mct_resid", "tsa_resid")

# Create a list to store the plots
plots <- list()

# Iterate over y-axis variables and create plots
for (variable in y_variables) {
  plot <- create_scatter_plot(variable)
  plots[[variable]] <- plot
}


# if you could save each plot as an .R object, then I could load the scatterplots on my machine and take care of sizing. 
# let me know if this works for you.

save(plots, file = here::here("Output/Figures", "Framingham_Figure2Associations_NewData.RData"))



