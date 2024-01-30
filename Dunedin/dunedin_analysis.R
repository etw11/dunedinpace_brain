# responding to NBA reviewers
# ethan whitman
# 10/27/23

library(haven)
library(superheat)


# mode function
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

###########################################
############ LOAD DUNEDIN DATA ############ 
###########################################

load('/Users/ew198/Documents/methylation/DunedinClocks45.rdata')
dunedinpace <- DunedinClocks45

load('/Users/ew198/Documents/methylation/MoreClocks_2022Sep22.rdata')
all_clocks <- data.frame(DunedinClocks45, MoreClocks_2022Sep22)

### load imaging data
load("/Users/ew198/Documents/individual_fns/behavior/Ethan_2022_0406.rdata")
behavior <- Ethan_2022_0406
behavior_all_clocks <- data.frame(behavior, dunedinpace = all_clocks$zrDunedinPACE45, horvath = all_clocks$zrHorvath45, hannum = all_clocks$zrHannum45, pheno = all_clocks$zrPheno45, grim = all_clocks$zrGrim45)

gid_clocks <- all_clocks$snum[complete.cases(all_clocks)]

fs <- read.csv('/Users/ew198/Documents/methylation/data/dunedin/aseg_FreeSurfer6.0.csv')
wmh <- read.csv('/Users/ew198/Documents/methylation/WMHvolume_UBO_052419.csv')
wmh$wmhVol_log <- log(as.numeric(wmh$wmhVol))
wmh$wmhVol_log[wmh$wmhVol_log == '.'] <- NA
brainage <- read.csv('/Users/ew198/Documents/methylation/BrainAge_liem_fs5.3_062419.csv') 
brainage$brainAgeGap45_ctrd[brainage$brainAgeGap45_ctrd == '.'] <- NA
surface_area <- read.csv('/Users/ew198/Documents/methylation/SA_HCPMPP.csv')
icv <- read.csv('/Users/ew198/Documents/methylation/Structural_wholebrain_HCPMPP.csv')

load("/Users/ew198/Documents/individual_fns/behavior/Ethan_2022_0406.rdata")
behavior <- Ethan_2022_0406
behavior_dunedinpace <- data.frame(behavior, pace = dunedinpace$zrDunedinPACE45, acc = dunedinpace$zAADunedinPACE45)
meth_subset <- data.frame(behavior[!is.na(dunedinpace$zrDunedinPACE45), ], dunedinpace[!is.na(dunedinpace$zrDunedinPACE45),])
gid_dunedinpace <- dunedinpace$snum[!is.na(dunedinpace$zrDunedinPACE45)]

structural_imaging <- data.frame(behavior_dunedinpace[,c(1,2,6:17)], wmh$wmhVol_log, brainage$brainAgeGap45_ctrd, as.numeric(surface_area$SA_TOT), as.numeric(icv$eTIV_aseg))
colnames(structural_imaging) <- c(colnames(structural_imaging)[1:14], 'brainAGE', 'wmh', 'total_sa', 'ICV')
structural_imaging_compete <- structural_imaging[complete.cases(structural_imaging),]
gid_imaging <- structural_imaging_compete$snum

gid_imaging_dunedinpace <- intersect(gid_imaging, gid_dunedinpace)

index <- behavior_dunedinpace$snum %in% gid_imaging_dunedinpace

image_dunpa_subset <- data.frame(behavior, 
                                 dunedinpace = all_clocks$zrDunedinPACE45, 
                                 horvath = all_clocks$zrHorvath45,
                                 hannum = all_clocks$zrHannum45,
                                 phenoage = all_clocks$zrPheno45,
                                 grimmage = all_clocks$zrGrim45,
                                 wmh = as.numeric(wmh$wmhVol_log), 
                                 brainage_gap = as.numeric(brainage$brainAgeGap45_ctrd),
                                 total_sa = as.numeric(surface_area$SA_TOT),
                                 ICV = as.numeric(icv$eTIV_aseg))[index,]
image_dunpa_subset$relative_hc <- resid(lm(img_HippocampVol_meanBilat45~ICV, data = image_dunpa_subset))

image_dunpa_subset <- merge(image_dunpa_subset, fs[,c('snum', 'WM.hypointensities')], by = 'snum')
image_dunpa_subset$WM.hypointensities_log <- log(as.numeric(image_dunpa_subset$WM.hypointensities))

# object is called WBCs45v2
load('/Users/ew198/Documents/methylation/WBCs45v2.rdata')

image_dunpa_subset <- merge(image_dunpa_subset, WBCs45v2[,-7], by = 'snum')
image_dunpa_subset$DNA_APOE_4count <- image_dunpa_subset$DNA_APOE_4count

####################################################
################# ORIGINAL RESULTS ################# 
####################################################

# dunedinpace

tbv.dp.res <-summary(lm(scale(img_BVtot45)~scale(dunedinpace)+sex+ICV, data = image_dunpa_subset))$coefficients
ct.dp.res <-summary(lm(scale(img_CT_AVG45)~scale(dunedinpace)+sex, data = image_dunpa_subset))$coefficients
sa.dp.res <-summary(lm(scale(total_sa)~scale(dunedinpace)+sex, data = image_dunpa_subset))$coefficients
hc.dp.res <-summary(lm(scale(img_HippocampVol_meanBilat45)~scale(dunedinpace)+sex+ICV, data = image_dunpa_subset))$coefficients
wmh.dp.res <-summary(lm(scale(wmh)~scale(dunedinpace)+sex, data = image_dunpa_subset))$coefficients
bag.dp.res <-summary(lm(scale(brainage_gap)~scale(dunedinpace)+sex, data = image_dunpa_subset))$coefficients


# horvath

tbv.hth.res <-summary(lm(scale(img_BVtot45)~scale(horvath)+sex+ICV, data = image_dunpa_subset))$coefficients
ct.hth.res <-summary(lm(scale(img_CT_AVG45)~scale(horvath)+sex, data = image_dunpa_subset))$coefficients
sa.hth.res <-summary(lm(scale(total_sa)~scale(horvath)+sex, data = image_dunpa_subset))$coefficients
hc.hth.res <-summary(lm(scale(img_HippocampVol_meanBilat45)~scale(horvath)+sex+ICV, data = image_dunpa_subset))$coefficients
wmh.hth.res <-summary(lm(scale(wmh)~scale(horvath)+sex, data = image_dunpa_subset))$coefficients
bag.hth.res <-summary(lm(scale(brainage_gap)~scale(horvath)+sex, data = image_dunpa_subset))$coefficients


# hannum

tbv.hnm.res <-summary(lm(scale(img_BVtot45)~scale(hannum)+sex+ICV, data = image_dunpa_subset))$coefficients
ct.hnm.res <-summary(lm(scale(img_CT_AVG45)~scale(hannum)+sex, data = image_dunpa_subset))$coefficients
sa.hnm.res <-summary(lm(scale(total_sa)~scale(hannum)+sex, data = image_dunpa_subset))$coefficients
hc.hnm.res <-summary(lm(scale(img_HippocampVol_meanBilat45)~scale(hannum)+sex+ICV, data = image_dunpa_subset))$coefficients
wmh.hnm.res <-summary(lm(scale(wmh)~scale(hannum)+sex, data = image_dunpa_subset))$coefficients
bag.hnm.res <-summary(lm(scale(brainage_gap)~scale(hannum)+sex, data = image_dunpa_subset))$coefficients

# phenoage

tbv.phe.res <-summary(lm(scale(img_BVtot45)~scale(phenoage)+sex+ICV, data = image_dunpa_subset))$coefficients
ct.phe.res <-summary(lm(scale(img_CT_AVG45)~scale(phenoage)+sex, data = image_dunpa_subset))$coefficients
sa.phe.res <-summary(lm(scale(total_sa)~scale(phenoage)+sex, data = image_dunpa_subset))$coefficients
hc.phe.res <-summary(lm(scale(img_HippocampVol_meanBilat45)~scale(phenoage)+sex+ICV, data = image_dunpa_subset))$coefficients
rel_hc.phe.res <-summary(lm(scale(relative_hc)~scale(phenoage)+sex, data = image_dunpa_subset))$coefficients
wmh.phe.res <-summary(lm(scale(wmh)~scale(phenoage)+sex, data = image_dunpa_subset))$coefficients
bag.phe.res <-summary(lm(scale(brainage_gap)~scale(phenoage)+sex, data = image_dunpa_subset))$coefficients

#grimmage

tbv.grim.res <-summary(lm(scale(img_BVtot45)~scale(grimmage)+sex+ICV, data = image_dunpa_subset))$coefficients
ct.grim.res <-summary(lm(scale(img_CT_AVG45)~scale(grimmage)+sex, data = image_dunpa_subset))$coefficients
sa.grim.res <-summary(lm(scale(total_sa)~scale(grimmage)+sex, data = image_dunpa_subset))$coefficients
hc.grim.res <-summary(lm(scale(img_HippocampVol_meanBilat45)~scale(grimmage)+sex+ICV, data = image_dunpa_subset))$coefficients
rel_hc.grim.res <-summary(lm(scale(relative_hc)~scale(grimmage)+sex, data = image_dunpa_subset))$coefficients
wmh.grim.res <-summary(lm(scale(wmh)~scale(grimmage)+sex, data = image_dunpa_subset))$coefficients
bag.grim.res <-summary(lm(scale(brainage_gap)~scale(grimmage)+sex, data = image_dunpa_subset))$coefficients




####################################################
############ WMHypo WMHyper correlation ############ 
####################################################

### Dunedin WMHypo-WMHyper correlation

hist(image_dunpa_subset$WM.hypointensities_log)

cor.test(image_dunpa_subset$wmh, image_dunpa_subset$WM.hypointensities_log)


# WMHypo

wmhypo.dp.res <-summary(lm(scale(WM.hypointensities_log)~scale(dunedinpace)+sex, data = image_dunpa_subset))$coefficients
wmhypo.hth.res <-summary(lm(scale(WM.hypointensities_log)~scale(horvath)+sex, data = image_dunpa_subset))$coefficients
wmhypo.hnm.res <-summary(lm(scale(WM.hypointensities_log)~scale(hannum)+sex, data = image_dunpa_subset))$coefficients
wmhypo.phe.res <-summary(lm(scale(WM.hypointensities_log)~scale(phenoage)+sex, data = image_dunpa_subset))$coefficients
wmhypo.grim.res <-summary(lm(scale(WM.hypointensities_log)~scale(grimmage)+sex, data = image_dunpa_subset))$coefficients

# WBC control

wmhypo.dp.res.wca <-summary(lm(scale(WM.hypointensities_log)~scale(dunedinpace)+sex+ CD4T45 + NK45 + Mono45 + Gran45 + PlasmaBlast45 + CD8pCD28nCD45RAn_45 + CD8naive45, data = image_dunpa_subset_wbc))$coefficients
wmhypo.hth.res.wca <-summary(lm(scale(WM.hypointensities_log)~scale(horvath)+sex+ CD4T45 + NK45 + Mono45 + Gran45 + PlasmaBlast45 + CD8pCD28nCD45RAn_45 + CD8naive45, data = image_dunpa_subset_wbc))$coefficients
wmhypo.hnm.res.wca <-summary(lm(scale(WM.hypointensities_log)~scale(hannum)+sex+ CD4T45 + NK45 + Mono45 + Gran45 + PlasmaBlast45 + CD8pCD28nCD45RAn_45 + CD8naive45, data = image_dunpa_subset_wbc))$coefficients
wmhypo.pheno.res.wca <-summary(lm(scale(WM.hypointensities_log)~scale(phenoage)+sex+ CD4T45 + NK45 + Mono45 + Gran45 + PlasmaBlast45 + CD8pCD28nCD45RAn_45 + CD8naive45, data = image_dunpa_subset_wbc))$coefficients
wmhypo.grim.res.wca <-summary(lm(scale(WM.hypointensities_log)~scale(grimmage)+sex+ CD4T45 + NK45 + Mono45 + Gran45 + PlasmaBlast45 + CD8pCD28nCD45RAn_45 + CD8naive45, data = image_dunpa_subset_wbc))$coefficients


# APOE control
load('/Users/ew198/Documents/methylation/WBCs45v2.rdata')

image_dunpa_subset_apoe <- merge(image_dunpa_subset, WBCs45v2[,-7], by = 'snum')

wmhypo.dp.res.apoe <-summary(lm(scale(WM.hypointensities_log)~scale(dunedinpace)+sex+DNA_APOE_4count.y, data = image_dunpa_subset_apoe))$coefficients
wmhypo.hth.res.apoe <-summary(lm(scale(WM.hypointensities_log)~scale(horvath)+sex+DNA_APOE_4count.y, data = image_dunpa_subset_apoe))$coefficients
wmhypo.hnm.res.apoe <-summary(lm(scale(WM.hypointensities_log)~scale(hannum)+sex+DNA_APOE_4count.y, data = image_dunpa_subset_apoe))$coefficients
wmhypo.pheno.res.apoe <-summary(lm(scale(WM.hypointensities_log)~scale(phenoage)+sex+DNA_APOE_4count.y, data = image_dunpa_subset_apoe))$coefficients
wmhypo.grim.res.apoe <-summary(lm(scale(WM.hypointensities_log)~scale(grimmage)+sex+DNA_APOE_4count.y, data = image_dunpa_subset_apoe))$coefficients



##### visualization ######


# WMHyper vs WMHypo correlation

cor.test(image_dunpa_subset$wmh, image_dunpa_subset$WM.hypointensities_log)


ggplot(data=data.frame(WMHyper = image_dunpa_subset$wmh, 
                       WMHypo = image_dunpa_subset$WM.hypointensities_log), 
       aes(x=WMHyper, y=WMHypo)) +
  geom_point() +
  labs(title='Dunedin')+
  #ylim(-5.5, 5.5)+
  #xlim(-5.5,5.5) +
  #geom_abline(intercept=0,slope=1, size = 1) +
  geom_smooth(method = 'lm', color = 'red') +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(vjust = 0.5, size = 10),
        axis.title.y = element_text(size = 15, vjust = 0),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 15)) 

dunedin_wmhypo_dist <- ggplot(data = image_dunpa_subset, aes(x = as.numeric(WM.hypointensities))) + 
  geom_histogram(fill='#F8766D') +
  xlim(-2000, 53000)+
  ylim(0, 600)+
  labs(title = "Dunedin",
       x = "WMHypo volume (mm3)",
       y = 'Frequency') +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20, vjust = 1),
        axis.title.x = element_text(size = 20),
        title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.position = 'none')

# WMHypo forest plot

wmhypo_wmhyper_table_dunedin <- data.frame(imaging = rep(c('WMHyper', 'WMHypo'), 5),
                                   clock = c(rep('dunedinpace', 2), rep('horvath', 2), rep('hannum',2), rep('phenoage',2), rep('grimage',2)),
                                   beta = c(wmh.dp.res[2,1],  wmhypo.dp.res[2,1], 
                                            wmh.hth.res[2,1],  wmhypo.hth.res[2,1],
                                            wmh.hnm.res[2,1],  wmhypo.hnm.res[2,1],
                                            wmh.phe.res[2,1],  wmhypo.phe.res[2,1],
                                            wmh.grim.res[2,1],  wmhypo.grim.res[2,1]),
                                   stderror = c(wmh.dp.res[2,2],  wmhypo.dp.res[2,2], 
                                                wmh.hth.res[2,2],  wmhypo.hth.res[2,2],
                                                wmh.hnm.res[2,2],  wmhypo.hnm.res[2,2],
                                                wmh.phe.res[2,2],  wmhypo.phe.res[2,2],
                                                wmh.grim.res[2,2],  wmhypo.grim.res[2,2]))

# point range plot

level_order <- rev(c('horvath', 'hannum', 'phenoage', 'grimage', 'dunedinpace'))

ggplot(wmhypo_wmhyper_table_dunedin, aes(x=clock, y=beta, color = imaging, group=imaging)) + 
  geom_pointrange(aes(ymin=beta-(1.96*stderror), ymax=beta+(1.96*stderror)),position=position_dodge(width=0.6)) +
  geom_hline(yintercept = 0)+
  labs(title = 'WMHyper vs WMHypo') +
  coord_flip() +
  theme_classic()+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 0),
        axis.title.x = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))




####################################################
################# APOE stratification ############## 
####################################################

image_dunpa_subset_noncarriers <- image_dunpa_subset[image_dunpa_subset$DNA_APOE_4present == 0,]
image_dunpa_subset_carriers <- image_dunpa_subset[image_dunpa_subset$DNA_APOE_4present == 1,]

# dunedinpace noncarriers

tbv.dp.res.nc <-summary(lm(scale(img_BVtot45)~scale(dunedinpace)+sex+ICV, data = image_dunpa_subset_noncarriers))$coefficients
hc.dp.res.nc <-summary(lm(scale(img_HippocampVol_meanBilat45)~scale(dunedinpace)+sex+ICV, data = image_dunpa_subset_noncarriers))$coefficients
wmh.dp.res.nc <-summary(lm(scale(wmh)~scale(dunedinpace)+sex, data = image_dunpa_subset_noncarriers))$coefficients
wmhypo.dp.res.nc <-summary(lm(scale(WM.hypointensities_log)~scale(dunedinpace)+sex, data = image_dunpa_subset_noncarriers))$coefficients
ct.dp.res.nc <-summary(lm(scale(img_CT_AVG45)~scale(dunedinpace)+sex, data = image_dunpa_subset_noncarriers))$coefficients
sa.dp.res.nc <-summary(lm(scale(total_sa)~scale(dunedinpace)+sex, data = image_dunpa_subset_noncarriers))$coefficients

# dunedinpace carriers

tbv.dp.res.c <-summary(lm(scale(img_BVtot45)~scale(dunedinpace)+sex+ICV, data = image_dunpa_subset_carriers))$coefficients
hc.dp.res.c <-summary(lm(scale(img_HippocampVol_meanBilat45)~scale(dunedinpace)+sex+ICV, data = image_dunpa_subset_carriers))$coefficients
wmh.dp.res.c <-summary(lm(scale(wmh)~scale(dunedinpace)+sex, data = image_dunpa_subset_carriers))$coefficients
wmhypo.dp.res.c <-summary(lm(scale(WM.hypointensities_log)~scale(dunedinpace)+sex, data = image_dunpa_subset_carriers))$coefficients
ct.dp.res.c <-summary(lm(scale(img_CT_AVG45)~scale(dunedinpace)+sex, data = image_dunpa_subset_carriers))$coefficients
sa.dp.res.c <-summary(lm(scale(total_sa)~scale(dunedinpace)+sex, data = image_dunpa_subset_carriers))$coefficients

# apoe stratification forest plot

dunedin_apoe_strat_table <- data.frame(imaging = rep(c('TBV', 'HC', 'WMHyper', 'CT', 'SA', 'WMHypo'), 3),
                                       genotype = c(rep('full sample', 6), rep('APOE4 carriers',6), rep('APOE4 noncarriers', 6)),
                                       beta = c(tbv.dp.res[2,1],  hc.dp.res[2,1], wmh.dp.res[2,1],  ct.dp.res[2,1], sa.dp.res[2,1],  wmhypo.dp.res[2,1],
                                                tbv.dp.res.c[2,1],  hc.dp.res.c[2,1], wmh.dp.res.c[2,1],  ct.dp.res.c[2,1], sa.dp.res.c[2,1],  wmhypo.dp.res.c[2,1],
                                                tbv.dp.res.nc[2,1],  hc.dp.res.nc[2,1], wmh.dp.res.nc[2,1],  ct.dp.res.nc[2,1], sa.dp.res.nc[2,1],  wmhypo.dp.res.nc[2,1]),
                                       stderror = c(tbv.dp.res[2,2],  hc.dp.res[2,2], wmh.dp.res[2,2],  ct.dp.res[2,2], sa.dp.res[2,2],  wmhypo.dp.res[2,2],
                                                    tbv.dp.res.c[2,2],  hc.dp.res.c[2,2], wmh.dp.res.c[2,2],  ct.dp.res.c[2,2], sa.dp.res.c[2,2],  wmhypo.dp.res.c[2,2],
                                                    tbv.dp.res.nc[2,2],  hc.dp.res.nc[2,2], wmh.dp.res.nc[2,2],  ct.dp.res.nc[2,2], sa.dp.res.nc[2,2],  wmhypo.dp.res.nc[2,2]))

# point range plot

imaging_level_order_strat <- rev(c('TBV', 'HC', 'WMHypo', 'CT', 'SA', 'WMHyper'))
genotype_levels <- c('full sample', 'APOE4 noncarriers', 'APOE4 carriers')
dunedin_apoe_colors <- c('#F8766D', 'grey70', 'yellow3')

dunedin_apoe_strat <- ggplot(dunedin_apoe_strat_table, 
                             aes(x=factor(imaging, imaging_level_order_strat), 
                                 y=beta, color = factor(genotype, levels=genotype_levels), group=factor(genotype, levels = genotype_levels))) + 
  geom_pointrange(aes(ymin=beta-(1.96*stderror), ymax=beta+(1.96*stderror)),position=position_dodge(width=0.6)) +
  ylim(-.35, .35)+
  geom_hline(yintercept = 0)+
  labs(title = 'Dunedin') +
  coord_flip() +
  theme_classic()+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 0),
        axis.title.x = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        title = element_text(size=15))+
  guides(color = guide_legend(reverse = TRUE, title = 'Genotype')) +
  scale_color_manual(values=c(rep(dunedin_apoe_colors,3)))

# save this
save(dunedin_apoe_strat, file = '/Users/ew198/Documents/methylation/results/figures/dunedin_apoe_strat.Rdata')

grid.arrange(dunedin_apoe_strat, adni_apoe_strat, fhs_apoe_strat, nrow = 3)


################################################################# 
################# time between blood draw and scan ############## 
################################################################# 

load('/Users/ew198/Documents/methylation/data/dunedin/Dunedin_Age.rdata')

image_dunpa_subset_td <- merge(image_dunpa_subset, Dunedin_Age[,c('snum', 'MRIChronAge45', 'AgeatInt45', 'ScanUnit_IntvDays', 'ScanUnit_absdays')], by = 'snum')

tbv.dp.res.td <-summary(lm(scale(img_BVtot45)~scale(dunedinpace)+sex+ICV+ScanUnit_IntvDays, data = image_dunpa_subset_td))$coefficients
hc.dp.res.td <-summary(lm(scale(img_HippocampVol_meanBilat45)~scale(dunedinpace)+sex+ICV+ScanUnit_IntvDays, data = image_dunpa_subset_td))$coefficients
wmh.dp.res.td <-summary(lm(scale(wmh)~scale(dunedinpace)+sex+ScanUnit_IntvDays, data = image_dunpa_subset_td))$coefficients
wmhypo.dp.res.td <-summary(lm(scale(WM.hypointensities_log)~scale(dunedinpace)+sex+ScanUnit_IntvDays, data = image_dunpa_subset_td))$coefficients
ct.dp.res.td <-summary(lm(scale(img_CT_AVG45)~scale(dunedinpace)+sex+ScanUnit_IntvDays, data = image_dunpa_subset_td))$coefficients
sa.dp.res.td <-summary(lm(scale(total_sa)~scale(dunedinpace)+sex+ScanUnit_IntvDays, data = image_dunpa_subset_td))$coefficients


hist(image_dunpa_subset_td$ScanUnit_absdays)


dunedin_scanint_hist <- ggplot(data = image_dunpa_subset_td, aes(x = ScanUnit_IntvDays)) + 
  geom_histogram(fill='#F8766D') +
  labs(title = "Dunedin",
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

time_diff_table <- data.frame(dataset = c( rep('Dunedin',6), rep('FHS-OC',6), rep('ADNI',6) ),
                              imaging = rep(c('TBV', 'HC', 'WMHyper', 'CT', 'SA', 'WMHypo'),3),
                              beta_td_control = c(tbv.dp.res.td[2,1], hc.dp.res.td[2,1], wmh.dp.res.td[2,1], ct.dp.res.td[2,1], sa.dp.res.td[2,1], wmhypo.dp.res.td[2,1],
                                                 rep(NA, 6), 
                                                  tbvrel.dp.plm.res.td[2,1], hcrel.dp.plm.res.td[2,1], wmh.dp.plm.res.td[2,1], ct.dp.plm.res.td[2,1], sa.dp.plm.res.td[2,1], wmhypo.dp.plm.res.td[2,1]),
                              std_error_td_control = c(tbv.dp.res.td[2,2], hc.dp.res.td[2,2], wmh.dp.res.td[2,2], ct.dp.res.td[2,2], sa.dp.res.td[2,2], wmhypo.dp.res.td[2,2],
                                                    rep(NA, 6), 
                                                     tbvrel.dp.plm.res.td[2,2], hcrel.dp.plm.res.td[2,2], wmh.dp.plm.res.td[2,2], ct.dp.plm.res.td[2,2], sa.dp.plm.res.td[2,2], wmhypo.dp.plm.res.td[2,2]),
                              beta = c(tbv.dp.res[2,1], hc.dp.res[2,1], wmh.dp.res[2,1], ct.dp.res[2,1], sa.dp.res[2,1], wmhypo.dp.res[2,1],
                                        c(-.04, -.07, NA, .09, -.09, .01), 
                                        tbvrel.dp.plm.res[2,1], hcrel.dp.plm.res[2,1], wmh.dp.plm.res[2,1], ct.dp.plm.res[2,1], sa.dp.plm.res[2,1], wmhypo.dp.plm.res[2,1]),
                              std_error = c(tbv.dp.res[2,2], hc.dp.res[2,2], wmh.dp.res[2,2], ct.dp.res[2,2], sa.dp.res[2,2], wmhypo.dp.res[2,2],
                                           rep(NA, 6), 
                                            tbvrel.dp.plm.res[2,2], hcrel.dp.plm.res[2,2], wmh.dp.plm.res[2,2], ct.dp.plm.res[2,2], sa.dp.plm.res[2,2], wmhypo.dp.plm.res[2,2])
                                       )
                              
write.csv(time_diff_table, file = '/Users/ew198/Documents/methylation/results/time_diff_table.csv')



