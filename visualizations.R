# NBA visualizations
# ethan whitman
# 1/30/24


# NOTE - this script must be run after running the ADNI and Dunedin analysis scripts - it depends on lots of objects from those analyses.



######################################################
#################### MAKE FIGURES ####################
######################################################

# FIGURE 2
# ADNI 
adni_age_at_meth <- merge(adni_methclocks, adni_baseline_dates, by = "RID")

adni_age_at_meth$age_at_meth <- adni_age_at_meth$age_at_baseline + interval(ymd(adni_age_at_meth$baseline_date),  mdy(adni_age_at_meth$EXAMDATE)) %/% months(1)

# Dunedin
load('/Users/ew198/Documents/methylation/data/dunedin/Dunedin_Age.rdata')
dunedin_chronage <- merge(image_dunpa_subset, Dunedin_Age, by = 'snum')

# FHS
load("/Users/ew198/Documents/methylation/from_dan/figures/OriginalDataFigures/Framingham_Figure1Densities_OriginalData.RData")
load("/Users/ew198/Documents/methylation/from_dan/figures/NewDataFigures/Framingham_Figure1Densities_NewData.RData")

ggplot() +
  geom_histogram(data = figure1$data[figure1$data$dataset=='ADNI',], aes(x=Age),  fill = '#619CFF', alpha = .5) +
  geom_histogram(data = figure1$data[figure1$data$dataset=='FHS',], aes(x=Age),  fill = '#00BA38', alpha = .5) +
  geom_histogram(data = data.frame(Age=round(dunedin_chronage$AgeatInt45,2)), aes(x=Age),  fill = '#F8766D', alpha = .5) +
  xlim(39, 100)+
  labs(x = 'Age in years',
       y = 'Frequency') +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15, vjust = 1),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 20))





# FIGURE 3
# loading FHS scatterplots from Dan

load("/Users/ew198/Documents/methylation/from_dan/figures/OriginalDataFigures/Framingham_Figure2Associations_OriginalData.RData")
load("/Users/ew198/Documents/methylation/from_dan/figures/NewDataFigures/Framingham_Figure2Associations_NewData.RData")

# mri forest plot

# ADNI

imaging_meth_summary_acceleration <- data.frame(imaging = rep(c('TBV', 'HC', 'WMHyper', 'CT', 'SA', 'WMHypo'), 5),
                                                clock = c(rep('dunedinpace', 6), rep('horvath', 6), rep('hannum',6), rep('phenoage',6), rep('grimage',6)),
                                                beta = c(tbvrel.dp.plm.res[2,1],  hcrel.dp.plm.res[2,1],        wmh.dp.plm.res[2,1],   ct.dp.plm.res[2,1],       sa.dp.plm.res[2,1], wmhypo.dp.plm.res[2,1],  
                                                         tbvrel.hth.plm.res[2,1], hcrel.hth.plm.res[2,1],       wmh.hth.plm.res[2,1],   ct.hth.plm.res[2,1],    sa.hth.plm.res[2,1],  wmhypo.hth.plm.res[2,1],  
                                                         tbvrel.hnm.plm.res[2,1], hcrel.hnm.plm.res[2,1],       wmh.hnm.plm.res[2,1],   ct.hnm.plm.res[2,1],    sa.hnm.plm.res[2,1],  wmhypo.hnm.plm.res[2,1],  
                                                         tbvrel.phe.plm.res[2,1], hcrel.phe.plm.res[2,1],       wmh.phe.plm.res[2,1],   ct.phe.plm.res[2,1],    sa.phe.plm.res[2,1],  wmhypo.phe.plm.res[2,1],  
                                                         tbvrel.grim.plm.res[2,1],hcrel.grim.plm.res[2,1],      wmh.grim.plm.res[2,1],   ct.grim.plm.res[2,1], sa.grim.plm.res[2,1],  wmhypo.grim.plm.res[2,1]),
                                                stderror = c(tbvrel.dp.plm.res[2,2],  hcrel.dp.plm.res[2,2],   wmh.dp.plm.res[2,1],    ct.dp.plm.res[2,2],       sa.dp.plm.res[2,2], wmhypo.dp.plm.res[2,2],  
                                                             tbvrel.hth.plm.res[2,2], hcrel.hth.plm.res[2,2],  wmh.hth.plm.res[2,1],    ct.hth.plm.res[2,2],    sa.hth.plm.res[2,2],  wmhypo.hth.plm.res[2,2],
                                                             tbvrel.hnm.plm.res[2,2], hcrel.hnm.plm.res[2,2],  wmh.hnm.plm.res[2,1],    ct.hnm.plm.res[2,2],    sa.hnm.plm.res[2,2],  wmhypo.hnm.plm.res[2,2],
                                                             tbvrel.phe.plm.res[2,2], hcrel.phe.plm.res[2,2],  wmh.phe.plm.res[2,1],    ct.phe.plm.res[2,2],    sa.phe.plm.res[2,2],  wmhypo.phe.plm.res[2,2],
                                                             tbvrel.grim.plm.res[2,2],hcrel.grim.plm.res[2,2], wmh.grim.plm.res[2,1],   ct.grim.plm.res[2,2], sa.grim.plm.res[2,2],  wmhypo.grim.plm.res[2,2]))



# Dunedin

dunedin_meth_imaging_pointrange <- data.frame(imaging = rep(c('TBV', 'CT', 'SA', 'HC', 'WMHyper', 'WMHypo'),5),
                                              clock = c(rep('dunedinpace', 6), rep('horvath', 6), rep('hannum',6), rep('phenoage',6), rep('grimage',6)),
                                              beta = c(tbv.dp.res[2,1],   ct.dp.res[2,1],  sa.dp.res[2,1],  hc.dp.res[2,1],       wmh.dp.res[2,1],   wmhypo.dp.res[2,1], 
                                                       tbv.hth.res[2,1],  ct.hth.res[2,1], sa.hth.res[2,1], hc.hth.res[2,1],      wmh.hth.res[2,1],  wmhypo.hth.res[2,1],
                                                       tbv.hnm.res[2,1],  ct.hnm.res[2,1], sa.hnm.res[2,1], hc.hnm.res[2,1],      wmh.hnm.res[2,1],  wmhypo.hnm.res[2,1],
                                                       tbv.phe.res[2,1],  ct.phe.res[2,1], sa.phe.res[2,1], hc.phe.res[2,1],      wmh.phe.res[2,1],  wmhypo.phe.res[2,1],
                                                       tbv.grim.res[2,1],  ct.grim.res[2,1], sa.grim.res[2,1], hc.grim.res[2,1],  wmh.grim.res[2,1], wmhypo.grim.res[2,1]),
                                              
                                              
                                              std_error = c(tbv.dp.res[2,2], ct.dp.res[2,2], ct.dp.res[2,2],        hc.dp.res[2,2],    wmh.dp.res[2,2], wmhypo.dp.res[2,2],
                                                            tbv.hth.res[2,2],    ct.hth.res[2,2], ct.hth.res[2,2],  hc.hth.res[2,2],   wmh.hth.res[2,2], wmhypo.hth.res[2,2],
                                                            tbv.hnm.res[2,2],    ct.hnm.res[2,2], ct.hnm.res[2,2],  hc.hnm.res[2,2],   wmh.hnm.res[2,2], wmhypo.hnm.res[2,2],
                                                            tbv.phe.res[2,2],    ct.phe.res[2,2], ct.phe.res[2,2],  hc.phe.res[2,2],   wmh.phe.res[2,2], wmhypo.phe.res[2,2],
                                                            tbv.grim.res[2,2],   ct.grim.res[2,2], ct.grim.res[2,2], hc.grim.res[2,2], wmh.grim.res[2,2], wmhypo.grim.res[2,2])
)



mri_forestplot_table <- data.frame(dataset = c(rep('ADNI', 5), rep('Dunedin', 5), rep('FHS-OC', 5)),
                                   imaging = rep(c('TBV', 'CT', 'SA', 'HC', 'WMHypo'),3),
                                   beta = c(tbvrel.dp.plm.res[2,1], ct.dp.plm.res[2,1], sa.dp.plm.res[2,1], hcrel.dp.plm.res[2,1], wmhypo.dp.plm.res[2,1],
                                            tbv.dp.res[2,1], ct.dp.res[2,1], sa.dp.res[2,1], hc.dp.res[2,1], wmh.dp.res[2,1],
                                            -0.03752815, -0.092464541, -0.009306709, -0.071963606, 0.093908079),
                                   std_error = c(tbvrel.dp.plm.res[2,2], ct.dp.plm.res[2,2], sa.dp.plm.res[2,2], hcrel.dp.plm.res[2,2], wmhypo.dp.plm.res[2,2],
                                                 tbv.dp.res[2,2], ct.dp.res[2,2], sa.dp.res[2,2], hc.dp.res[2,2], wmh.dp.res[2,2],
                                                 0.0172009, 0.0325889, 0.02985215, 0.03030376, 0.03555245)
)


dataset_colors <- rev(c('#F8766D', '#00BA38', '#619CFF'))

mri_levels <- c('TBV', 'HC', 'WMHypo', 'CT', 'SA')


mri_pointrange <- ggplot(mri_forestplot_table, 
                         aes(x=factor(imaging, levels = rev(mri_levels)), y=beta, 
                             color = factor(dataset, levels = rev(c('Dunedin', 'FHS-OC', 'ADNI'))), 
                             group=factor(dataset, levels = rev(c('Dunedin', 'FHS-OC', 'ADNI'))))) + 
  geom_pointrange(aes(ymin=beta-(1.96*std_error), ymax=beta+(1.96*std_error)),position=position_dodge(width=0.6)) +
  geom_hline(yintercept = 0)+
  labs(title = 'DunedinPACE and Brain Anatomy', 
       x = 'MRI phenotype') +
  coord_flip() +
  theme_classic()+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15, hjust = .5),
        title = element_text(size = 0),
        legend.text = element_text(size = 15))+
  guides(color = guide_legend(reverse = TRUE, title = 'Dataset'))+
  scale_color_manual(values=c(rep(dataset_colors,3)))



# scatterplots

# ADNI

tbv.p.dat.plot <- tbv.p.dat
tbv.p.dat.plot$tbv_resid <- resid(lm(TBV~sex+mean_age+mean_age^2+sex*mean_age+sex*mean_age^2+ICV, data = tbv.p.dat.plot))

summary(lm(tbv_resid~zPoAm45, data = tbv.p.dat.plot))

adni_tbv_scatter <- ggplot(data=as.data.frame(tbv.p.dat.plot), aes(y=tbv_resid, x=zPoAm45)) +
  geom_point(color = '#619CFF') +
  geom_smooth(method = 'lm', color = 'black', se=FALSE) +
  ylim(-3,2.5)+
  xlim(-4.4,5.7)+
  labs(y = "TBV",
       x = "DunedinPACE", 
       title = 'ADNI') +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(vjust = 0.5, size = 10),
        axis.title.y = element_text(size = 15, vjust = 0),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 0))


hc.p.dat.plot <- hc.p.dat
hc.p.dat.plot$hc_resid <- resid(lm(HC~sex+mean_age+mean_age^2+sex*mean_age+sex*mean_age^2+ICV, data = hc.p.dat.plot))

summary(lm(hc_resid~zPoAm45, data = hc.p.dat.plot))

adni_hc_scatter <- ggplot(data=as.data.frame(hc.p.dat.plot), aes(y=hc_resid, x=zPoAm45)) +
  geom_point(color = '#619CFF') +
  geom_smooth(method = 'lm', color = 'black', se = FALSE) +
  ylim(-3.1,3.1)+
  xlim(-4.4,5.7)+
  labs(y = "HC",
       x = "DunedinPACE", 
       title = 'ADNI') +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(vjust = 0.5, size = 10),
        axis.title.y = element_text(size = 15, vjust = 0),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 0))

wmhypo.p.dat.plot <- fs.wmhypo.p.dat
wmhypo.p.dat.plot$wmhypo_resid <- resid(lm(wmhypo.combat.log~sex+mean_age+mean_age^2+sex*mean_age+sex*mean_age^2, data = wmhypo.p.dat.plot))


adni_wmh_scatter <- ggplot(data=as.data.frame(wmhypo.p.dat.plot), aes(y=wmhypo_resid, x=zPoAm45)) +
  geom_point(color = '#619CFF') +
  geom_smooth(method = 'lm', color = 'black', se=FALSE) +
  ylim(-5.3,5.1)+
  xlim(-4.4,5.7)+
  labs(y = "WMHypo",
       x = "DunedinPACE", 
       title = 'ADNI') +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(vjust = 0.5, size = 10),
        axis.title.y = element_text(size = 15, vjust = 0),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 0))


ct.p.dat.plot <- fs.p.dat
ct.p.dat.plot$mean_thickness <- as.numeric(ct.p.dat.plot$mean_thickness)
ct.p.dat.plot$ct_resid <- resid(lm(mean_thickness~sex+mean_age+mean_age^2+sex*mean_age+sex*mean_age^2, data = ct.p.dat.plot))

summary(lm(ct_resid~zPoAm45, data = ct.p.dat.plot))

adni_ct_scatter <- ggplot(data=as.data.frame(ct.p.dat.plot), aes(y=ct_resid, x=zPoAm45)) +
  geom_point(color = '#619CFF') +
  geom_smooth(method = 'lm', color = 'black', se=FALSE) +
  ylim(-5.1,3.2)+
  xlim(-4.4,5.7)+
  labs(y = "CT",
       x = "DunedinPACE", 
       title = 'ADNI') +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(vjust = 0.5, size = 10),
        axis.title.y = element_text(size = 15, vjust = 0),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 0))


sa.p.dat.plot <- fs.p.dat
sa.p.dat.plot$sa_resid <- resid(lm(sa~sex+mean_age+mean_age^2+sex*mean_age+sex*mean_age^2, data = sa.p.dat.plot))

summary(lm(sa_resid~zPoAm45, data = sa.p.dat.plot))

adni_sa_scatter <- ggplot(data=as.data.frame(sa.p.dat.plot), aes(y=sa_resid, x=zPoAm45)) +
  geom_point(color = '#619CFF') +
  geom_smooth(method = 'lm', color = 'black', se=FALSE) +
  ylim(-2.7,3)+
  xlim(-4.4,5.7)+
  labs(y = "SA",
       x = "DunedinPACE", 
       title = 'ADNI') +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(vjust = 0.5, size = 10),
        axis.title.y = element_text(size = 15, vjust = 0),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 0))



# Dunedin

image_dunpa_subset_plot <- image_dunpa_subset
image_dunpa_subset_plot$tbv_resid <- resid(lm(scale(img_BVtot45)~sex+scale(ICV), data = image_dunpa_subset_plot))

summary(lm(tbv_resid~dunedinpace, data = image_dunpa_subset_plot))

dunedin_tbv_scatter <- ggplot(data=image_dunpa_subset_plot, aes(y=tbv_resid, x=dunedinpace)) +
  geom_point(color = '#F8766D') +
  geom_smooth(method = 'lm', color = 'black', se=FALSE) +
  ylim(-3,2.5)+
  xlim(-4.4,5.7)+
  labs(y = "TBV",
       x = "DunedinPACE", 
       title = 'Dunedin') +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(vjust = 0.5, size = 10),
        axis.title.y = element_text(size = 15, vjust = 0),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 0)) 


image_dunpa_subset_plot$hc_resid <- resid(lm(scale(img_HippocampVol_meanBilat45)~sex+scale(ICV), data = image_dunpa_subset_plot))

summary(lm(hc_resid~dunedinpace, data = image_dunpa_subset_plot))

dunedin_hc_scatter <- ggplot(data=image_dunpa_subset_plot, aes(y=hc_resid, x=dunedinpace)) +
  geom_point(color = '#F8766D') +
  geom_smooth(method = 'lm', color = 'black', se = FALSE) +
  ylim(-3.1,3.1)+
  xlim(-4.4,5.7)+
  labs(y = "HC",
       x = "DunedinPACE", 
       title = 'Dunedin') +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(vjust = 0.5, size = 10),
        axis.title.y = element_text(size = 15, vjust = 0),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 0)) 

image_dunpa_subset_plot$wmhypo_resid <- resid(lm(scale(WM.hypointensities_log)~sex, data = image_dunpa_subset_plot))

summary(lm(wmhypo_resid~dunedinpace, data = image_dunpa_subset_plot))

dunedin_wmh_scatter <- ggplot(data=image_dunpa_subset_plot, aes(y=wmhypo_resid, x=dunedinpace)) +
  geom_point(color = '#F8766D') +
  geom_smooth(method = 'lm', color = 'black', se=FALSE) +
  ylim(-5.3,5.1)+
  xlim(-4.4,5.7)+
  labs(y = "WMHypo",
       x = "DunedinPACE", 
       title = 'Dunedin') +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(vjust = 0.5, size = 10),
        axis.title.y = element_text(size = 15, vjust = 0),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 0)) 


image_dunpa_subset_plot$ct_resid <- resid(lm(scale(img_CT_AVG45)~sex, data = image_dunpa_subset_plot))

summary(lm(ct_resid~dunedinpace, data = image_dunpa_subset_plot))

dunedin_ct_scatter <- ggplot(data=image_dunpa_subset_plot, aes(y=ct_resid, x=dunedinpace)) +
  geom_point(color = '#F8766D') +
  geom_smooth(method = 'lm', color = 'black', se=FALSE) +
  ylim(-5.1,3.2)+
  xlim(-4.4,5.7)+
  labs(y = "CT",
       x = "DunedinPACE", 
       title = 'Dunedin') +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(vjust = 0.5, size = 10),
        axis.title.y = element_text(size = 15, vjust = 0),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 0)) 

image_dunpa_subset_plot$sa_resid <- resid(lm(scale(total_sa)~sex, data = image_dunpa_subset_plot))

summary(lm(sa_resid~dunedinpace, data = image_dunpa_subset_plot))

dunedin_sa_scatter <- ggplot(data=image_dunpa_subset_plot, aes(y=sa_resid, x=dunedinpace)) +
  geom_point(color = '#F8766D') +
  geom_smooth(method = 'lm', color = 'black', se=FALSE) +
  ylim(-2.7,3)+
  xlim(-4.4,5.7)+
  labs(y = "SA",
       x = "DunedinPACE", 
       title = 'Dunedin') +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(vjust = 0.5, size = 10),
        axis.title.y = element_text(size = 15, vjust = 0),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 0)) 


### FHS-OC

# this is the green shade i want:
#619CFF

## "scatterplots" - can only plot slope right now because Dan has raw data

fhs_tbv_scatter <- ggplot(data=plots$tbv_resid$data, aes(y=tbv_resid, x=scale(dunedinpace)))+
  geom_point(color = '#00BA38') +
  geom_smooth(method = 'lm', color = 'black', se=FALSE) +
  ylim(-3,2.5)+
  xlim(-4.4,5.7)+
  labs(y = "TBV",
       x = "DunedinPACE", 
       title = 'FHS-OC') +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(vjust = 0.5, size = 10),
        axis.title.y = element_text(size = 15, vjust = 0),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 0)) 


fhs_hc_scatter <- ggplot(data=plots$hv_resid$data, aes(y=hv_resid, x=scale(dunedinpace)))+
  geom_point(color = '#00BA38') +
  geom_smooth(method = 'lm', color = 'black', se=FALSE)+
  ylim(-3.1,3.1)+
  xlim(-4.4,5.7)+
  labs(y = "HC",
       x = "DunedinPACE", 
       title = 'FHS-OC')  +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(vjust = 0.5, size = 10),
        axis.title.y = element_text(size = 15, vjust = 0),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 0)) 

fhs_wmh_scatter <- ggplot(data=plots$wmh_resid$data, aes(y=wmh_resid, x=scale(dunedinpace)))+
  geom_point(color = '#00BA38') +
  geom_smooth(method = 'lm', color = 'black', se=FALSE) +
  ylim(-5.3,5.1)+
  xlim(-4.4,5.7)+
  labs(y = "WMHypo",
       x = "DunedinPACE", 
       title = 'FHS-OC') +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(vjust = 0.5, size = 10),
        axis.title.y = element_text(size = 15, vjust = 0),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 0))


fhs_sa_scatter <- ggplot(data=plots$tsa_resid$data, aes(y=tsa_resid, x=scale(dunedinpace)))+
  geom_point(color = '#00BA38') +
  geom_smooth(method = 'lm', color = 'black', se=FALSE) +
  ylim(-2.7,3)+
  xlim(-4.4,5.7)+
  labs(y = "SA",
       x = "DunedinPACE", 
       title = 'Dunedin') +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(vjust = 0.5, size = 10),
        axis.title.y = element_text(size = 15, vjust = 0),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 0)) 

fhs_ct_scatter <- ggplot(data=plots$mct_resid$data, aes(y=mct_resid, x=scale(dunedinpace)))+
  geom_point(color = '#00BA38') +
  geom_smooth(method = 'lm', color = 'black', se=FALSE) +
  ylim(-5.1,3.2)+
  xlim(-4.4,5.7)+
  labs(y = "CT",
       x = "DunedinPACE", 
       title = 'Dunedin') +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(vjust = 0.5, size = 10),
        axis.title.y = element_text(size = 15, vjust = 0),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 0)) 


grid.arrange(mri_pointrange, 
             dunedin_tbv_scatter, fhs_tbv_scatter, adni_tbv_scatter,
             dunedin_hc_scatter,  fhs_hc_scatter,  adni_hc_scatter,
             dunedin_wmh_scatter, fhs_wmh_scatter, adni_wmh_scatter,
             dunedin_ct_scatter, fhs_ct_scatter, adni_ct_scatter,
             dunedin_sa_scatter, fhs_sa_scatter, adni_sa_scatter,
             nrow = 6, ncol = 3, layout_matrix = cbind(c(1,2,5,8,11,14), c(1,3,6,9,12,15), c(1,4,7,10,13,16)),
             heights = c(3,1,1,1,1,1))


# just pointrange
mri_pointrange

# just scatterplots - 800 x 800
grid.arrange(dunedin_tbv_scatter, fhs_tbv_scatter, adni_tbv_scatter,
             dunedin_hc_scatter,  fhs_hc_scatter,  adni_hc_scatter,
             dunedin_wmh_scatter, fhs_wmh_scatter, adni_wmh_scatter,
             dunedin_ct_scatter, fhs_ct_scatter, adni_ct_scatter,
             dunedin_sa_scatter, fhs_sa_scatter, adni_sa_scatter,
             nrow = 5, ncol = 3)



#### WMHyper figures for supplement


wmhyper_forestplot_table <- data.frame(dataset = c('Dunedin', 'ADNI'),
                                       beta = c(wmh.dp.res[2,1], wmh.dp.plm.res[2,1]),
                                       std_error = c(wmh.dp.res[2,2], wmh.dp.plm.res[2,2])
)


wmh_dataset_colors <- rev(c('#F8766D','#619CFF'))


wmhyper_pointrange <- ggplot(wmhyper_forestplot_table, 
                             aes(x=dataset, y=beta, 
                                 color = factor(dataset, levels = rev(c('Dunedin', 'ADNI'))), 
                                 group=factor(dataset, levels = rev(c('Dunedin', 'ADNI'))))) + 
  geom_pointrange(aes(ymin=beta-(1.96*std_error), ymax=beta+(1.96*std_error)),position=position_dodge(width=0.6)) +
  geom_hline(yintercept = 0)+
  ylim(-.1, .21)+
  coord_flip() +
  theme_classic()+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 0),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 0),
        legend.text = element_text(size = 15),
        legend.position='none')+
  guides(color = guide_legend(reverse = TRUE, title = 'Dataset'))+
  scale_color_manual(values=c(wmh_dataset_colors))

image_dunpa_subset_plot$wmhyper_resid <- resid(lm(scale(wmh)~sex+scale(ICV), data = image_dunpa_subset_plot))
dunedin_wmhyper_scatter <- ggplot(data=image_dunpa_subset_plot, aes(y=wmhyper_resid, x=dunedinpace)) +
  geom_point(color = '#F8766D') +
  geom_smooth(method = 'lm', color = 'black', se=FALSE) +
  ylim(-5.3,5.1)+
  xlim(-4.4,5.7)+
  labs(y = "WMHyper",
       x = "DunedinPACE", 
       title = 'Dunedin') +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(vjust = 0.5, size = 10),
        axis.title.y = element_text(size = 15, vjust = 0),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 0)) 


wmh.p.dat.plot <- wmh.p.dat
wmh.p.dat.plot$wmhyper_resid <- resid(lm(wmh~sex+mean_age+mean_age^2+sex*mean_age+sex*mean_age^2, data = wmh.p.dat.plot))

summary(lm(wmh_resid~zPoAm45, data = wmh.p.dat.plot))

adni_wmhyper_scatter <- ggplot(data=as.data.frame(wmh.p.dat.plot), aes(y=wmhyper_resid, x=zPoAm45)) +
  geom_point(color = '#619CFF') +
  geom_smooth(method = 'lm', color = 'black', se=FALSE) +
  ylim(-5.3,5.1)+
  xlim(-4.4,5.7)+
  labs(y = "WMHyper",
       x = "DunedinPACE", 
       title = 'ADNI') +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(vjust = 0.5, size = 10),
        axis.title.y = element_text(size = 15, vjust = 0),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 0))

grid.arrange(dunedin_wmhyper_scatter, adni_wmhyper_scatter, nrow = 1)


###### FIGURE 4 ######

framingham_meth_imaging_pointrange <- data.frame(imaging = rep(c('TBV', 'HC', 'SA', 'WMH', 'CT'),5),
                                                 clock = c(rep('dunedinpace', 5), rep('horvath', 5), rep('hannum',5), rep('phenoage',5), rep('grimage',5)),
                                                 beta = c(-0.0375281497447516,  -0.071963606,  -.00931068, .09390808, -0.092464541,
                                                          -.01124434, -.01852601, .01476129,  .03215491, 0.000884573,
                                                          .01196693,  -.01386389, .01766917,  .04898107, -0.003534288,
                                                          -.01633458, -.04337284, .03100742,  .07708464, -0.046354683,
                                                          -.09036906, -.0922754,  .00927577,  .0580993, -0.133374675),
                                                 
                                                 
                                                 ci_l = c(-0.071241917, -0.131358964,  -.06677337,  .0341427 ,-0.156338781,
                                                          -.05895267, -.07547289, -.04091806,  -.02602951,-0.054003595,
                                                          -.03792766, -.07342726, -.04058882,  -.01185608,-0.06101394,
                                                          -.06400632, -.10022542, -.024434,    .01930126 ,-0.101974641,
                                                          -.14063551, -.15239417, -.04955799,  -.00329772,-0.191710234),
                                                 
                                                 ci_h = c(-0.003814383, -0.0125682483423219,  .04815201,  .15367346,-0.028590301,
                                                          .03646399, .03842087,  .07044065,  .09033933, 0.055772741,
                                                          .06186153, .04569949,  .07592715,  .10981821, 0.053945363,
                                                          .03133716, .01347973,  .08644883,  .13486802, 0.009265276,
                                                          -.0401026, -.03215663, .06810954,  .11949632, -0.075039116)
)


## comparing clocks figure


framingham_meth_imaging_pointrange$std_error <- abs(framingham_meth_imaging_pointrange$beta - framingham_meth_imaging_pointrange$ci_l) / 1.96

# adni
imaging_meth_summary_acceleration_group_plot <- imaging_meth_summary_acceleration
colnames(imaging_meth_summary_acceleration_group_plot)[4] <- 'std_error'

## tbv panel


clock_colors <- rev(c('red', 'blue', 'navyblue','royalblue3', 'skyblue2'))

clock_levels <- rev(c('dunedinpace', 'horvath', 'hannum', 'phenoage', 'grimage'))

# TBV
tbv_comparing_clocks_pointrange <- data.frame(
  dataset = c(rep('ADNI', 5), rep('Dunedin',5), rep('FHS-OC', 5)),
  rbind(imaging_meth_summary_acceleration_group_plot[imaging_meth_summary_acceleration_group_plot$imaging == 'TBV',],
        dunedin_meth_imaging_pointrange[dunedin_meth_imaging_pointrange$imaging == 'TBV',],
        framingham_meth_imaging_pointrange[framingham_meth_imaging_pointrange$imaging == 'TBV',c('imaging', 'clock', 'beta', 'std_error')])
)


tbv_compclocks_pr <- ggplot(tbv_comparing_clocks_pointrange, 
                            aes(x=factor(dataset, levels = c('ADNI', 'FHS-OC', 'Dunedin')),
                                y=beta, group = factor(clock, levels = clock_levels), color = factor(clock, levels = clock_levels))) + 
  geom_pointrange(aes(ymin=beta-(1.96*std_error), ymax=beta+(1.96*std_error)), position=position_dodge(width=0.6)) +
  geom_hline(yintercept = 0)+
  ylim(-.27, .27)+
  labs(title = 'TBV') +
  coord_flip() +
  theme_classic()+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 0),
        axis.title.x = element_text(size = 15),
        legend.text = element_text(size=15),
        legend.position = 'none',
        title = element_text(size = 20))+
  guides(color = guide_legend(reverse = TRUE, title = 'Clock')) +
  scale_color_manual(values=c(rep(clock_colors,3)))


# HC
hc_comparing_clocks_pointrange <- data.frame(
  dataset = c(rep('ADNI', 5), rep('Dunedin',5), rep('FHS-OC', 5)),
  rbind(imaging_meth_summary_acceleration_group_plot[imaging_meth_summary_acceleration_group_plot$imaging == 'HC',],
        dunedin_meth_imaging_pointrange[dunedin_meth_imaging_pointrange$imaging == 'HC',],
        framingham_meth_imaging_pointrange[framingham_meth_imaging_pointrange$imaging == 'HC',c('imaging', 'clock', 'beta', 'std_error')])
)

hc_compclocks_pr <- ggplot(hc_comparing_clocks_pointrange, 
                           aes(x=factor(dataset, levels = c('ADNI', 'FHS-OC', 'Dunedin')),
                               y=beta, group = factor(clock, levels = clock_levels), color = factor(clock, levels = clock_levels))) + 
  geom_pointrange(aes(ymin=beta-(1.96*std_error), ymax=beta+(1.96*std_error)), position=position_dodge(width=0.6)) +
  geom_hline(yintercept = 0)+
  ylim(-.27, .27)+
  labs(title = 'HC') +
  coord_flip() +
  theme_classic()+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 0),
        axis.title.x = element_text(size = 15),
        legend.text = element_text(size=15),
        legend.position = 'none',
        title = element_text(size = 20))+
  guides(color = guide_legend(reverse = TRUE, title = 'Clock')) +
  scale_color_manual(values=c(rep(clock_colors,3)))



# WMHypo
wmhypo_comparing_clocks_pointrange <- data.frame(
  dataset = c(rep('ADNI', 5), rep('Dunedin',5), rep('FHS-OC', 5)),
  rbind(imaging_meth_summary_acceleration_group_plot[imaging_meth_summary_acceleration_group_plot$imaging == 'WMHypo',],
        dunedin_meth_imaging_pointrange[dunedin_meth_imaging_pointrange$imaging == 'WMHypo',],
        framingham_meth_imaging_pointrange[framingham_meth_imaging_pointrange$imaging == 'WMH',c('imaging', 'clock', 'beta', 'std_error')])
)

wmhypo_compclocks_pr <- ggplot(wmhypo_comparing_clocks_pointrange, 
                               aes(x=factor(dataset, levels = c('ADNI', 'FHS-OC', 'Dunedin')),
                                   y=beta, group = factor(clock, levels = clock_levels), color = factor(clock, levels = clock_levels))) + 
  geom_pointrange(aes(ymin=beta-(1.96*std_error), ymax=beta+(1.96*std_error)), position=position_dodge(width=0.6)) +
  geom_hline(yintercept = 0)+
  ylim(-.27, .27)+
  labs(title = 'WMHypo') +
  coord_flip() +
  theme_classic()+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 0),
        axis.title.x = element_text(size = 15),
        legend.text = element_text(size=15),
        legend.position = 'none',
        title = element_text(size = 20))+
  guides(color = guide_legend(reverse = TRUE, title = 'Clock')) +
  scale_color_manual(values=c(rep(clock_colors,3)))




# CT
ct_comparing_clocks_pointrange <- data.frame(
  dataset = c(rep('ADNI', 5), rep('Dunedin',5), rep('FHS-OC', 5)),
  rbind(imaging_meth_summary_acceleration_group_plot[imaging_meth_summary_acceleration_group_plot$imaging == 'CT',],
        dunedin_meth_imaging_pointrange[dunedin_meth_imaging_pointrange$imaging == 'CT',],
        framingham_meth_imaging_pointrange[framingham_meth_imaging_pointrange$imaging == 'CT',c('imaging', 'clock', 'beta', 'std_error')])
)


ct_compclocks_pr <- ggplot(ct_comparing_clocks_pointrange, 
                           aes(x=factor(dataset, levels = c('ADNI', 'FHS-OC', 'Dunedin')),
                               y=beta, group = factor(clock, levels = clock_levels), color = factor(clock, levels = clock_levels))) + 
  geom_pointrange(aes(ymin=beta-(1.96*std_error), ymax=beta+(1.96*std_error)), position=position_dodge(width=0.6)) +
  geom_hline(yintercept = 0)+
  ylim(-.27, .27)+
  labs(title = 'CT') +
  coord_flip() +
  theme_classic()+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 0),
        axis.title.x = element_text(size = 15),
        legend.text = element_text(size=15),
        legend.position = 'none',
        title = element_text(size = 20))+
  guides(color = guide_legend(reverse = TRUE, title = 'Clock')) +
  scale_color_manual(values=c(rep(clock_colors,3)))



# SA
sa_comparing_clocks_pointrange <- data.frame(
  dataset = c(rep('ADNI', 5), rep('Dunedin',5), rep('FHS-OC', 5)),
  rbind(imaging_meth_summary_acceleration_group_plot[imaging_meth_summary_acceleration_group_plot$imaging == 'SA',],
        dunedin_meth_imaging_pointrange[dunedin_meth_imaging_pointrange$imaging == 'SA',],
        framingham_meth_imaging_pointrange[framingham_meth_imaging_pointrange$imaging == 'SA',c('imaging', 'clock', 'beta', 'std_error')])
)


sa_compclocks_pr <- ggplot(sa_comparing_clocks_pointrange, 
                           aes(x=factor(dataset, levels = c('ADNI', 'FHS-OC', 'Dunedin')),
                               y=beta, group = factor(clock, levels = clock_levels), color = factor(clock, levels = clock_levels))) + 
  geom_pointrange(aes(ymin=beta-(1.96*std_error), ymax=beta+(1.96*std_error)), position=position_dodge(width=0.6)) +
  geom_hline(yintercept = 0)+
  ylim(-.27, .27)+
  labs(title = 'SA') +
  coord_flip() +
  theme_classic()+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 0),
        axis.title.x = element_text(size = 15),
        legend.text = element_text(size=15),
        legend.position = 'none',
        title = element_text(size = 20))+
  guides(color = guide_legend(reverse = TRUE, title = 'Clock')) +
  scale_color_manual(values=c(rep(clock_colors,3)), 
                     labels=rev(c('DunedinPACE', 'Horvath', 'Hannum', 'PhenoAge', 'GrimAge')))




# 1000 x 700
grid.arrange(tbv_compclocks_pr, hc_compclocks_pr, wmhypo_compclocks_pr, ct_compclocks_pr, sa_compclocks_pr, nrow = 2)


summary_cc <- rbind(tbv_comparing_clocks_pointrange, hc_comparing_clocks_pointrange, wmhypo_comparing_clocks_pointrange, ct_comparing_clocks_pointrange, sa_comparing_clocks_pointrange)

summary_cc$ci_l <- summary_cc$beta - (1.96*summary_cc$std_error)
summary_cc$ci_h <- summary_cc$beta + (1.96*summary_cc$std_error)


write.csv(summary_cc, file = '/Users/ew198/Documents/methylation/results/comparing_clocks_table_11_6_23.csv')

