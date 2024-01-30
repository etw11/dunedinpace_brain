# making NBA revisions with the ADNI data
# 10/30/23
# ethan whitman

library(haven)
library(superheat)
library(raveio)
library(tidyverse)
install.packages('devtools')
devtools::install_github("jcbeer/longCombat")
library(longCombat)
library(invgamma)
library(lme4)
library(viridis)
library(nlme)

# mode function
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

########################################
############ LOAD ADNI DATA ############ 
########################################

# initialization

adni_methclocks <- as.data.frame(read.delim('/Users/ew198/Documents/methylation/from_karen/ADNI_methclocks_pluspheno.txt'))

adni_meth <- read.delim("/Users/ew198/Documents/methylation/adni_methIDlist.txt", header = TRUE)
adni_meth$AGE_MONTH <- adni_meth$AGE*12
adnimerge <- read_csv("/Users/ew198/Documents/methylation/forEthan_070622/Data/ADNIMERGE.csv")

adni_sex <- data.frame(RID = adnimerge$RID, sex = adnimerge$PTGENDER)
adni_sex <- adni_sex[!duplicated(adni_sex),]

adni_baseline_dates <- data.frame(RID = unique(adni_meth$RID), 
                                  baseline_date = rep(NA, length(unique(adni_meth$RID))),
                                  age_at_baseline = rep(NA, length(unique(adni_meth$RID))))
x <- 1
for (i in unique(adni_meth$RID)){
  temp_df <- as.data.frame(adnimerge[adnimerge$RID == i,])
  temp_row <- temp_df[temp_df$VISCODE == 'bl',]
  adni_baseline_dates$baseline_date[x] <- as.character(temp_row$EXAMDATE)
  adni_baseline_dates$age_at_baseline[x] <- temp_row$AGE*12
  x <- x + 1
}

adni_meth <- merge(adni_meth, adni_baseline_dates, by = 'RID')

# ADNI methylation dates

adni_meth_datedrawn <- data.frame(RID = unique(adni_meth$RID), 
                                  DateDrawn1 = rep(NA, length(unique(adni_meth$RID))),
                                  DateDrawn2 = rep(NA, length(unique(adni_meth$RID))),
                                  DateDrawn3 = rep(NA, length(unique(adni_meth$RID))),
                                  DateDrawn4 = rep(NA, length(unique(adni_meth$RID))),
                                  DateDrawn5 = rep(NA, length(unique(adni_meth$RID))),
                                  age_at_timepoint1 = rep(NA, length(unique(adni_meth$RID))),
                                  age_at_timepoint2 = rep(NA, length(unique(adni_meth$RID))),
                                  age_at_timepoint3 = rep(NA, length(unique(adni_meth$RID))),
                                  age_at_timepoint4 = rep(NA, length(unique(adni_meth$RID))),
                                  age_at_timepoint5 = rep(NA, length(unique(adni_meth$RID))),
                                  first_timepoint = rep(NA, length(unique(adni_meth$RID))),
                                  last_timepoint = rep(NA, length(unique(adni_meth$RID))),
                                  age_at_last_timepoint = rep(NA, length(unique(adni_meth$RID))))

x <- 1
for (i in unique(adni_meth$RID)){
  temp_df <- adni_meth[adni_meth$RID == i,]
  timepoints <- nrow(temp_df)
  for (t in 1:timepoints){
    adni_meth_datedrawn[x, (t+1)] <- temp_df$DateDrawn[t]
    adni_meth_datedrawn[x, (t+6)] <- temp_df$AGE_MONTH[1] + interval(ymd(temp_df$baseline_date[1]), mdy(temp_df$DateDrawn[t])) %/% months(1)
  }
  adni_meth_datedrawn[x, 12] <- temp_df$DateDrawn[1]
  adni_meth_datedrawn[x, 13] <- temp_df$DateDrawn[timepoints]
  adni_meth_datedrawn[x, 14] <- temp_df$AGE_MONTH[1] + interval(ymd(temp_df$baseline_date[1]), mdy(temp_df$DateDrawn[timepoints])) %/% months(1)
  x <- x+1
}

# age plot age at methylation collection
adni_meth_datedrawn <- merge(adni_meth_datedrawn, adni_baseline_dates, by = 'RID')
adni_meth_datedrawn <- adni_meth_datedrawn[order(adni_meth_datedrawn$age_at_baseline),]
adni_meth_datedrawn <- data.frame(adni_meth_datedrawn, plot_id = c(1:649))




# imaging initialization
adnimerge <- read_csv("/Users/ew198/Documents/methylation/forEthan_070622/Data/ADNIMERGE.csv")
MRImeta<-read_csv("/Users/ew198/Documents/methylation/forEthan_070622/Data/MRI3META.csv")
methIDs<-read_table("/Users/ew198/Documents/methylation/forEthan_070622/Data/adni_methIDlist.txt") %>% dplyr::select(RID,VISCODE,Edate)
methIDs<-methIDs %>% mutate(methID = paste(RID,VISCODE,sep="_"))

UCD_WMH <-read_csv("/Users/ew198/Documents/methylation/forEthan_070622/Data/ADNI_UCD_WMH_05_02_22.csv") %>% dplyr::select(RID,COLPROT,VISCODE,VISCODE2,EXAMDATE,MAGNETICFIELDSTRENGTH,TOTAL_WMH,TOTAL_BRAIN)

#Select all rows from ADNI data that belong to a participant with methylation
methWMH <- semi_join(UCD_WMH, methIDs,by=c("RID")) 
methMergeStrict<-semi_join(adnimerge, methIDs,by=c("RID","VISCODE")) 

methMerge<-dplyr::semi_join(adnimerge, methIDs,by=c("RID")) 
methMerge<-semi_join(adnimerge, methIDs,by=c("RID")) 
methMerge <- merge(methMerge, adni_baseline_dates, by = 'RID')




###### freesurfer initialization ###### 



# FreeSurfer version 4
UCSF_FS4 <- read.csv("/Users/ew198/Documents/methylation/forEthan_070622/Data/UCSFFSX_11_02_15.csv")

methIDs <- read_table("/Users/ew198/Documents/methylation/forEthan_070622/Data/adni_methIDlist.txt")[,c('RID', 'VISCODE', 'DateDrawn')]
methIDs<-methIDs %>% mutate(methID = paste(RID,VISCODE,sep="_"))

methFS4 <- semi_join(UCSF_FS4, methIDs,by=c("RID")) 

methFS4 <- merge(methFS4, adni_baseline_dates, by='RID')

# read in ADNI freesurfer dictionary

adni_fs5_dict <- read.csv('/Users/ew198/Documents/methylation/forEthan_070622/Data/UCSFFSL51_DICT_03_01_22.csv')

sa_index <-  grepl("SA", adni_fs5_dict$FLDNAME)
sa_index[24:25] <- FALSE
ta_index <-  grepl("TA", adni_fs5_dict$FLDNAME)
ta_index[11] <- FALSE

cortex_ct_sa_df <- data.frame(sa_index = adni_fs5_dict$FLDNAME[sa_index], sa_label = adni_fs5_dict$TEXT[sa_index],
                              ta_index = adni_fs5_dict$FLDNAME[ta_index], ts_label = adni_fs5_dict$TEXT[ta_index])


mean_thickness_fs4 <- rep(NA, nrow(methFS4))
for (i in 1:nrow(methFS4)){
  mean_thickness_fs4[i] <- sum((methFS4[i,c(cortex_ct_sa_df$sa_index)] * methFS4[i,c(cortex_ct_sa_df$ta_index)])) / sum(methFS4[i,c(cortex_ct_sa_df$sa_index)])
}

total_sa_fs4 <- rep(NA, nrow(methFS4))
for (i in 1:nrow(methFS4)){
  total_sa_fs4[i] <- sum(methFS4[i,c(cortex_ct_sa_df$sa_index)])
}

methFS4$mean_thickness <- mean_thickness_fs4
methFS4$total_sa <- total_sa_fs4

# filter out failing QC scans - added 8/9

methFS4 <- methFS4[methFS4$OVERALLQC == 'Pass',]



# FreeSurfer version 6

UCSF_FS6 <- read.csv('/Users/ew198/Documents/methylation/forEthan_070622/Data/UCSFFSX6_06_07_22.csv')
methFS6 <- semi_join(UCSF_FS6, methIDs,by=c("RID")) 
methFS6 <- merge(methFS6, adni_baseline_dates, by='RID')

# filter bad QC
methFS6 <- methFS6[methFS6$OVERALLQC=='Pass',]

# calculating mean thickness and total sa in FS6

adni_fs6_dict <- read.csv('/Users/ew198/Documents/methylation/forEthan_070622/Data/UCSFFSX6_DICT_08_27_19.csv')

sa_index <-  grepl("SA", adni_fs6_dict$FLDNAME)
sa_index[c(144,299)] <- FALSE
ta_index <-  grepl("TA", adni_fs6_dict$FLDNAME)
ta_index[10] <- FALSE


cortex_ct_sa_df <- data.frame(sa_index = adni_fs6_dict$FLDNAME[sa_index], sa_label = adni_fs6_dict$TEXT[sa_index],
                              ta_index = adni_fs6_dict$FLDNAME[ta_index], ts_label = adni_fs6_dict$TEXT[ta_index])


mean_thickness <- rep(NA, nrow(methFS6))
for (i in 1:nrow(methFS6)){
  mean_thickness[i] <- sum((methFS6[i,c(cortex_ct_sa_df$sa_index)] * methFS6[i,c(cortex_ct_sa_df$ta_index)])) / sum(methFS6[i,c(cortex_ct_sa_df$sa_index)])
}

total_sa <- rep(NA, nrow(methFS6))
for (i in 1:nrow(methFS6)){
  total_sa[i] <- sum(methFS6[i,c(cortex_ct_sa_df$sa_index)])
}

methFS6$mean_thickness <- mean_thickness
methFS6$total_sa <- total_sa




##### freesurfer version 5.3

UCSF_FS5long <-read_csv("/Users/ew198/Documents/methylation/forEthan_070622/Data/UCSFFSL51_03_01_22.csv")
methFS5long <- semi_join(UCSF_FS5long, methIDs,by=c("RID"))

### calculating mean thickness in ADNI
# also do a sanity check on the total surface area value

# read in ADNI freesurfer dictionary

adni_fs5_dict <- read.csv('/Users/ew198/Documents/methylation/forEthan_070622/Data/UCSFFSL51_DICT_03_01_22.csv')

sa_index <-  grepl("SA", adni_fs5_dict$FLDNAME)
sa_index[24:25] <- FALSE
ta_index <-  grepl("TA", adni_fs5_dict$FLDNAME)
ta_index[11] <- FALSE

cortex_ct_sa_df <- data.frame(sa_index = adni_fs5_dict$FLDNAME[sa_index], sa_label = adni_fs5_dict$TEXT[sa_index],
                              ta_index = adni_fs5_dict$FLDNAME[ta_index], ts_label = adni_fs5_dict$TEXT[ta_index])

mean_thickness <- rep(NA, nrow(methFS5long))
for (i in 1:nrow(methFS5long)){
  mean_thickness[i] <- sum((methFS5long[i,c(cortex_ct_sa_df$sa_index)] * methFS5long[i,c(cortex_ct_sa_df$ta_index)])) / sum(methFS5long[i,c(cortex_ct_sa_df$sa_index)])
}

total_sa <- rep(NA, nrow(methFS5long))
for (i in 1:nrow(methFS5long)){
  total_sa[i] <- sum(methFS5long[i,c(cortex_ct_sa_df$sa_index)])
}

methFS5long$mean_thickness <- mean_thickness
methFS5long$total_sa <- total_sa

# filtering out QC fails - added 8/9
methFS5long <- methFS5long[!is.na(methFS5long$OVERALLQC),]
methFS5long <- methFS5long[(methFS5long$OVERALLQC == 'Pass'),]

methFS5<- merge(methFS5long, adni_baseline_dates, by='RID')






### ComBat on the freesurver observations


# freesurfer measures

# take highest UID for duplicate scans
methFS4_new <- methFS4
for (i in unique(methFS4$RID)){
  temp_df <- methFS4[methFS4$RID==i,]
  if (sum(duplicated(temp_df$EXAMDATE)) == 1){
    print(paste0('single duplicate in ', i))
    dup_scans <- temp_df[temp_df$EXAMDATE==(Mode(temp_df$EXAMDATE)),]
    scan_loser <- min(dup_scans$IMAGEUID)
    methFS4_new <- methFS4_new[methFS4_new$IMAGEUID != scan_loser,]
  } else if (sum(duplicated(temp_df$EXAMDATE)) > 1){
    print(paste0('double duplicate in ', i))
    duplicated_dates <- data.frame(table(temp_df$EXAMDATE))[data.frame(table(temp_df$EXAMDATE))$Freq>1,]$Var1
    for (d in 1:length(duplicated_dates)){
      dup_scans <- temp_df[ymd(temp_df$EXAMDATE) == ymd(duplicated_dates[d]),]
      scan_loser <- min(dup_scans$IMAGEUID)
      methFS4_new <- methFS4_new[methFS4_new$IMAGEUID != scan_loser,]
    }
  }
}
methFS4 <- methFS4_new


methFS5_new <- methFS5
for (i in unique(methFS5$RID)){
  temp_df <- methFS5[methFS5$RID==i,]
  if (sum(duplicated(temp_df$EXAMDATE)) == 1){
    print(paste0('duplicate in ', i))
    dup_scans <- temp_df[temp_df$EXAMDATE==(Mode(temp_df$EXAMDATE)),]
    scan_loser <- min(dup_scans$IMAGEUID)
    methFS5_new <- methFS5_new[methFS5_new$IMAGEUID != scan_loser,]
  } else if (sum(duplicated(temp_df$EXAMDATE)) > 1){
    print(paste0('double duplicate in ', i))
    duplicated_dates <- data.frame(table(temp_df$EXAMDATE))[data.frame(table(temp_df$EXAMDATE))$Freq>1,]$Var1
    for (d in 1:length(duplicated_dates)){
      dup_scans <- temp_df[ymd(temp_df$EXAMDATE) == ymd(duplicated_dates[d]),]
      scan_loser <- min(dup_scans$IMAGEUID)
      methFS5_new <- methFS5_new[methFS5_new$IMAGEUID != scan_loser,]
    }
  }
}
methFS5 <- methFS5_new

methFS6_new <- methFS6
for (i in unique(methFS6$RID)){
  temp_df <- methFS6[methFS6$RID==i,]
  if (sum(duplicated(temp_df$EXAMDATE)) == 1){
    print(paste0('duplicate in ', i))
    dup_scans <- temp_df[temp_df$EXAMDATE==(Mode(temp_df$EXAMDATE)),]
    scan_loser <- min(dup_scans$IMAGEUID)
    methFS6_new <- methFS6_new[methFS6_new$IMAGEUID != scan_loser,]
  } else if (sum(duplicated(temp_df$EXAMDATE)) > 1){
    print(paste0('double duplicate in ', i))
    duplicated_dates <- data.frame(table(temp_df$EXAMDATE))[data.frame(table(temp_df$EXAMDATE))$Freq>1,]$Var1
    for (d in 1:length(duplicated_dates)){
      dup_scans <- temp_df[ymd(temp_df$EXAMDATE) == ymd(duplicated_dates[d]),]
      scan_loser <- min(dup_scans$IMAGEUID)
      methFS6_new <- methFS6_new[methFS6_new$IMAGEUID != scan_loser,]
    }
  }
}
methFS6 <- methFS6_new

# fix date formatting
methFS4$EXAMDATE <- as.Date(methFS4$EXAMDATE)
methFS4$VERSION <- as.Date(methFS4$VERSION)
methFS4$RUNDATE <- as.Date(methFS4$RUNDATE)

methFS5$EXAMDATE <- as.Date(methFS5$EXAMDATE)
methFS5$VERSION <- as.Date(methFS5$VERSION)
methFS5$RUNDATE <- as.Date(methFS5$RUNDATE)

methFS6$EXAMDATE <- as.Date(methFS6$EXAMDATE)
methFS6$VERSION <- as.Date(methFS6$VERSION)
methFS6$RUNDATE <- as.Date(methFS6$RUNDATE)


methFScomb_temp <- rbind(methFS4[,c('RID', 'EXAMDATE', 'VERSION', 'RUNDATE', 'OVERALLQC', 'baseline_date', 'age_at_baseline', 'mean_thickness', 'total_sa')],
                         methFS5[,c('RID', 'EXAMDATE', 'VERSION', 'RUNDATE','OVERALLQC',  'baseline_date', 'age_at_baseline', 'mean_thickness', 'total_sa')],
                         methFS6[,c('RID', 'EXAMDATE', 'VERSION', 'RUNDATE','OVERALLQC',  'baseline_date', 'age_at_baseline', 'mean_thickness', 'total_sa')])

methFScomb_temp <- data.frame(methFScomb_temp, FSVERSION=c(rep(4, nrow(methFS4)), rep(5, nrow(methFS5)), rep(6, nrow(methFS6)) ) )
methFScomb_temp$age_at_scan <- methFScomb_temp$age_at_baseline + (interval(methFScomb_temp$baseline_date, methFScomb_temp$EXAMDATE) %/% months(1))
methFScomb_temp <- merge(methFScomb_temp, adni_sex, by = 'RID')

# run longCombat_adni_fs.R to get harmonized values for CT and SA
methFScomb_temp$scanage_scale <- scale(methFScomb_temp$age_at_scan)
methFScomb_forharm <- methFScomb_temp
methFScomb_forharm$age_years <- methFScomb_forharm$age_at_scan/12
adni_combat <- longCombat(idvar='RID', 
                          timevar='age_years',
                          batchvar='FSVERSION', 
                          features=c('mean_thickness', 'total_sa'), 
                          formula='sex',
                          ranef='(1+age_years|RID)',
                          data=methFScomb_forharm[,c('RID', 'mean_thickness', 'total_sa', 'FSVERSION', 'sex', 'age_at_scan', 'age_years')])

adni_fsharm <- adni_combat$data_combat
adni_fsharm$sex <- methFScomb_forharm$sex
adni_fsharm$scanage_scale <- scale(adni_fsharm$age_years)

if (identical(methFScomb_temp$RID, adni_fsharm$RID) && sum(abs(methFScomb_temp$scanage_scale - adni_fsharm$scanage_scale)) < 0.0000000001 ) {
  print('IDENTICAL')
  methFScomb_harmed <- data.frame(methFScomb_temp, adni_fsharm[,-6])
} else{
  print('NOT IDENTICAL')
}





methFS4_harmed <- methFScomb_harmed[methFScomb_harmed$FSVERSION==4,]
methFS5_harmed <- methFScomb_harmed[methFScomb_harmed$FSVERSION==5,]




# PAIRING IMAGING OBSERVATIONS WITH METHYLATION

adni_tbv_datescanned <- data.frame(RID = unique(methMerge$RID), 
                                   DateScanned1 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned2 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned3 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned4 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned5 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned6 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned7 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned8 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned9 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned10 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned11 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned12 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned13 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned14 = rep(NA, length(unique(methMerge$RID))),
                                   
                                   age_at_scan1 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan2 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan3 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan4 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan5 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan6 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan7 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan8 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan9 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan10 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan11 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan12 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan13 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan14 = rep(NA, length(unique(methMerge$RID))),
                                   
                                   first_scan = rep(NA, length(unique(methMerge$RID))),
                                   last_scan = rep(NA, length(unique(methMerge$RID))),
                                   age_at_last_scan = rep(NA, length(unique(methMerge$RID))))

x <- 1
for (i in unique(methMerge$RID)){
  temp_df <- as.data.frame(methMerge[methMerge$RID == i,])
  temp_df <- temp_df[!is.na(temp_df$WholeBrain),]
  temp_df <- temp_df[order(temp_df$EXAMDATE),]
  temp_df$AGE_MONTH <- 12*temp_df$AGE
  timepoints <- nrow(temp_df)
  for (t in 1:timepoints){
    adni_tbv_datescanned[x, (t+1)] <- as.character(temp_df$EXAMDATE[t])
    adni_tbv_datescanned[x, (t+15)] <- temp_df$AGE_MONTH[1] + interval(ymd(temp_df$baseline_date[1]), ymd(temp_df$EXAMDATE[t])) %/% months(1)
  }
  adni_tbv_datescanned[x, 30] <- as.character(temp_df$EXAMDATE[1])
  adni_tbv_datescanned[x, 31] <- as.character(temp_df$EXAMDATE[timepoints])
  adni_tbv_datescanned[x, 32] <- temp_df$AGE_MONTH[1] + interval(ymd(temp_df$baseline_date[1]), ymd(temp_df$EXAMDATE[timepoints])) %/% months(1)
  x <- x+1
}

# age plot age at whole brain volume collection
adni_tbv_datescanned <- merge(adni_tbv_datescanned, 
                              data.frame(RID = adni_meth_datedrawn$RID, meth_timepoint1 = adni_meth_datedrawn$age_at_timepoint1),
                              by = 'RID')
adni_tbv_datescanned <- merge(adni_tbv_datescanned, adni_baseline_dates, by = 'RID')
adni_tbv_datescanned <- adni_tbv_datescanned[order(adni_tbv_datescanned$age_at_baseline),]
adni_tbv_datescanned <- data.frame(adni_tbv_datescanned, plot_id = c(1:649))


# hippocampus volume

adni_hc_datescanned <- data.frame(RID = unique(methMerge$RID), 
                                  DateScanned1 = rep(NA, length(unique(methMerge$RID))),
                                  DateScanned2 = rep(NA, length(unique(methMerge$RID))),
                                  DateScanned3 = rep(NA, length(unique(methMerge$RID))),
                                  DateScanned4 = rep(NA, length(unique(methMerge$RID))),
                                  DateScanned5 = rep(NA, length(unique(methMerge$RID))),
                                  DateScanned6 = rep(NA, length(unique(methMerge$RID))),
                                  DateScanned7 = rep(NA, length(unique(methMerge$RID))),
                                  DateScanned8 = rep(NA, length(unique(methMerge$RID))),
                                  DateScanned9 = rep(NA, length(unique(methMerge$RID))),
                                  DateScanned10 = rep(NA, length(unique(methMerge$RID))),
                                  DateScanned11 = rep(NA, length(unique(methMerge$RID))),
                                  DateScanned12 = rep(NA, length(unique(methMerge$RID))),
                                  DateScanned13 = rep(NA, length(unique(methMerge$RID))),
                                  
                                  age_at_scan1 = rep(NA, length(unique(methMerge$RID))),
                                  age_at_scan2 = rep(NA, length(unique(methMerge$RID))),
                                  age_at_scan3 = rep(NA, length(unique(methMerge$RID))),
                                  age_at_scan4 = rep(NA, length(unique(methMerge$RID))),
                                  age_at_scan5 = rep(NA, length(unique(methMerge$RID))),
                                  age_at_scan6 = rep(NA, length(unique(methMerge$RID))),
                                  age_at_scan7 = rep(NA, length(unique(methMerge$RID))),
                                  age_at_scan8 = rep(NA, length(unique(methMerge$RID))),
                                  age_at_scan9 = rep(NA, length(unique(methMerge$RID))),
                                  age_at_scan10 = rep(NA, length(unique(methMerge$RID))),
                                  age_at_scan11 = rep(NA, length(unique(methMerge$RID))),
                                  age_at_scan12 = rep(NA, length(unique(methMerge$RID))),
                                  age_at_scan13 = rep(NA, length(unique(methMerge$RID))),
                                  
                                  first_scan = rep(NA, length(unique(methMerge$RID))),
                                  last_scan = rep(NA, length(unique(methMerge$RID))),
                                  age_at_last_scan = rep(NA, length(unique(methMerge$RID))))


x <- 1
for (i in unique(methMerge$RID)){
  temp_df <- as.data.frame(methMerge[methMerge$RID == i,])
  temp_df <- temp_df[!is.na(temp_df$Hippocampus),]
  temp_df <- temp_df[order(temp_df$EXAMDATE),]
  temp_df$AGE_MONTH <- 12*temp_df$AGE
  timepoints <- nrow(temp_df)
  if (timepoints > 0){
    for (t in 1:timepoints){
      adni_hc_datescanned[x, (t+1)] <- as.character(temp_df$EXAMDATE[t])
      adni_hc_datescanned[x, (t+14)] <- temp_df$AGE_MONTH[1] + interval(ymd(temp_df$baseline_date[1]), ymd(temp_df$EXAMDATE[t])) %/% months(1)
    }
    adni_hc_datescanned[x, 29] <- as.character(temp_df$EXAMDATE[1])
    adni_hc_datescanned[x, 30] <- as.character(temp_df$EXAMDATE[timepoints])
    adni_hc_datescanned[x, 31] <- temp_df$AGE_MONTH[1] + interval(ymd(temp_df$baseline_date[1]), ymd(temp_df$EXAMDATE[timepoints])) %/% months(1)
  }
  x <- x+1
}


#hippocampus
adni_hc_datescanned <- merge(adni_hc_datescanned, 
                             data.frame(RID = adni_meth_datedrawn$RID, meth_timepoint1 = adni_meth_datedrawn$age_at_timepoint1),
                             by = 'RID')
adni_hc_datescanned <- merge(adni_hc_datescanned, adni_baseline_dates, by = 'RID')
adni_hc_datescanned <- adni_hc_datescanned[order(adni_hc_datescanned$age_at_baseline),]
adni_hc_datescanned <- data.frame(adni_hc_datescanned, plot_id = c(1:649))



# WMH

wmh_merge <- merge(adni_meth_datedrawn, methWMH,  by = 'RID', all = TRUE)
wmh_merge$wmhVol_log <- log(wmh_merge$TOTAL_WMH)

adni_wmh_datescanned <- data.frame(RID = unique(methMerge$RID), 
                                   DateScanned1 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned2 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned3 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned4 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned5 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned6 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned7 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned8 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned9 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned10 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned11 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned12 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned13 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned14 = rep(NA, length(unique(methMerge$RID))),
                                   
                                   age_at_scan1 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan2 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan3 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan4 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan5 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan6 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan7 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan8 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan9 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan10 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan11 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan12 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan13 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan14 = rep(NA, length(unique(methMerge$RID))),
                                   
                                   first_scan = rep(NA, length(unique(methMerge$RID))),
                                   last_scan = rep(NA, length(unique(methMerge$RID))),
                                   age_at_last_scan = rep(NA, length(unique(methMerge$RID))))

x <- 1
for (i in unique(wmh_merge$RID)){
  temp_df <- as.data.frame(wmh_merge[wmh_merge$RID == i,])
  temp_df <- temp_df[!is.na(temp_df$TOTAL_WMH),]
  temp_df <- temp_df[order(temp_df$EXAMDATE),]
  timepoints <- nrow(temp_df)
  if (timepoints > 0){
    for (t in 1:timepoints){
      adni_wmh_datescanned[x, (t+1)] <- as.character(temp_df$EXAMDATE[t])
      adni_wmh_datescanned[x, (t+15)] <- temp_df$age_at_baseline[1] + interval(ymd(temp_df$EXAMDATE[1]), ymd(temp_df$EXAMDATE[t])) %/% months(1)
    }
    adni_wmh_datescanned[x, 30] <- as.character(temp_df$EXAMDATE[1])
    adni_wmh_datescanned[x, 31] <- as.character(temp_df$EXAMDATE[timepoints])
    adni_wmh_datescanned[x, 32] <- temp_df$age_at_baseline[1] + interval(ymd(temp_df$EXAMDATE[1]), ymd(temp_df$EXAMDATE[timepoints])) %/% months(1)
  }
  x <- x+1
}

# age plot age at whole brain volume collection
adni_wmh_datescanned <- merge(adni_wmh_datescanned, 
                              data.frame(RID = adni_meth_datedrawn$RID, meth_timepoint1 = adni_meth_datedrawn$age_at_timepoint1),
                              by = 'RID')
adni_wmh_datescanned <- merge(adni_wmh_datescanned, adni_baseline_dates, by = 'RID')
adni_wmh_datescanned <- adni_wmh_datescanned[order(adni_wmh_datescanned$age_at_baseline),]
adni_wmh_datescanned <- data.frame(adni_wmh_datescanned, plot_id = c(1:649))






# date of diagnosis

adni_dx_datediagnosed <- data.frame(RID = unique(methMerge$RID), 
                                    DateDiagnosed1 = rep(NA, length(unique(methMerge$RID))),
                                    DateDiagnosed2 = rep(NA, length(unique(methMerge$RID))),
                                    DateDiagnosed3 = rep(NA, length(unique(methMerge$RID))),
                                    DateDiagnosed4 = rep(NA, length(unique(methMerge$RID))),
                                    DateDiagnosed5 = rep(NA, length(unique(methMerge$RID))),
                                    DateDiagnosed6 = rep(NA, length(unique(methMerge$RID))),
                                    DateDiagnosed7 = rep(NA, length(unique(methMerge$RID))),
                                    DateDiagnosed8 = rep(NA, length(unique(methMerge$RID))),
                                    DateDiagnosed9 = rep(NA, length(unique(methMerge$RID))),
                                    DateDiagnosed10 = rep(NA, length(unique(methMerge$RID))),
                                    DateDiagnosed11 = rep(NA, length(unique(methMerge$RID))),
                                    DateDiagnosed12 = rep(NA, length(unique(methMerge$RID))),
                                    DateDiagnosed13 = rep(NA, length(unique(methMerge$RID))),
                                    DateDiagnosed14 = rep(NA, length(unique(methMerge$RID))),
                                    DateDiagnosed15 = rep(NA, length(unique(methMerge$RID))),
                                    DateDiagnosed16 = rep(NA, length(unique(methMerge$RID))),
                                    DateDiagnosed17 = rep(NA, length(unique(methMerge$RID))),
                                    
                                    age_at_dx1 = rep(NA, length(unique(methMerge$RID))),
                                    age_at_dx2 = rep(NA, length(unique(methMerge$RID))),
                                    age_at_dx3 = rep(NA, length(unique(methMerge$RID))),
                                    age_at_dx4 = rep(NA, length(unique(methMerge$RID))),
                                    age_at_dx5 = rep(NA, length(unique(methMerge$RID))),
                                    age_at_dx6 = rep(NA, length(unique(methMerge$RID))),
                                    age_at_dx7 = rep(NA, length(unique(methMerge$RID))),
                                    age_at_dx8 = rep(NA, length(unique(methMerge$RID))),
                                    age_at_dx9 = rep(NA, length(unique(methMerge$RID))),
                                    age_at_dx10 = rep(NA, length(unique(methMerge$RID))),
                                    age_at_dx11 = rep(NA, length(unique(methMerge$RID))),
                                    age_at_dx12 = rep(NA, length(unique(methMerge$RID))),
                                    age_at_dx13 = rep(NA, length(unique(methMerge$RID))),
                                    age_at_dx14 = rep(NA, length(unique(methMerge$RID))),
                                    age_at_dx15 = rep(NA, length(unique(methMerge$RID))),
                                    age_at_dx16 = rep(NA, length(unique(methMerge$RID))),
                                    age_at_dx17 = rep(NA, length(unique(methMerge$RID))),
                                    
                                    first_dx = rep(NA, length(unique(methMerge$RID))),
                                    last_dx = rep(NA, length(unique(methMerge$RID))),
                                    age_at_last_dx = rep(NA, length(unique(methMerge$RID))))

x <- 1
for (i in unique(methMerge$RID)){
  temp_df <- as.data.frame(methMerge[methMerge$RID == i,])
  temp_df$DX[temp_df$DX == ""] <- NA
  temp_df <- temp_df[!is.na(temp_df$DX),]
  temp_df <- temp_df[order(temp_df$EXAMDATE),]
  temp_df$AGE_MONTH <- 12*temp_df$AGE
  timepoints <- nrow(temp_df)
  for (t in 1:timepoints){
    adni_dx_datediagnosed[x, (t+1)] <- as.character(temp_df$EXAMDATE[t])
    adni_dx_datediagnosed[x, (t+18)] <- temp_df$AGE_MONTH[1] + interval(ymd(temp_df$baseline_date[1]), ymd(temp_df$EXAMDATE[t])) %/% months(1)
  }
  adni_dx_datediagnosed[x, 38] <- as.character(temp_df$EXAMDATE[1])
  adni_dx_datediagnosed[x, 39] <- as.character(temp_df$EXAMDATE[timepoints])
  adni_dx_datediagnosed[x, 40] <- temp_df$AGE_MONTH[1] + interval(ymd(temp_df$baseline_date[1]), ymd(temp_df$EXAMDATE[timepoints])) %/% months(1)
  x <- x+1
}

# age plot age at diagnosis
adni_dx_datediagnosed <- merge(adni_dx_datediagnosed, 
                               data.frame(RID = adni_meth_datedrawn$RID, meth_timepoint1 = adni_meth_datedrawn$age_at_timepoint1),
                               by = 'RID')
adni_dx_datediagnosed <- merge(adni_dx_datediagnosed, adni_baseline_dates, by = 'RID')
adni_dx_datediagnosed <- adni_dx_datediagnosed[order(adni_dx_datediagnosed$age_at_baseline),]
adni_dx_datediagnosed <- data.frame(adni_dx_datediagnosed, plot_id = c(1:649))



















adni_fs6_datescanned <- data.frame(RID = unique(methMerge$RID), 
                                   DateScanned1 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned2 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned3 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned4 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned5 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned6 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned7 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned8 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned9 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned10 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned11 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned12 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned13 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned14 = rep(NA, length(unique(methMerge$RID))),
                                   
                                   age_at_scan1 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan2 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan3 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan4 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan5 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan6 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan7 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan8 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan9 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan10 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan11 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan12 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan13 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan14 = rep(NA, length(unique(methMerge$RID))),
                                   
                                   first_scan = rep(NA, length(unique(methMerge$RID))),
                                   last_scan = rep(NA, length(unique(methMerge$RID))),
                                   age_at_last_scan = rep(NA, length(unique(methMerge$RID))))

x <- 1
for (i in unique(methMerge$RID)){
  temp_df <- as.data.frame(methFS6[methFS6$RID == i,])
  temp_df <- temp_df[order(temp_df$EXAMDATE),]
  timepoints <- nrow(temp_df)
  if (timepoints > 0){
    for (t in 1:timepoints){
      adni_fs6_datescanned[x, (t+1)] <- as.character(temp_df$EXAMDATE[t])
      adni_fs6_datescanned[x, (t+15)] <- temp_df$age_at_baseline[1] + interval(ymd(temp_df$baseline_date[1]), ymd(temp_df$EXAMDATE[t])) %/% months(1)
    }
    adni_fs6_datescanned[x, 30] <- as.character(temp_df$EXAMDATE[1])
    adni_fs6_datescanned[x, 31] <- as.character(temp_df$EXAMDATE[timepoints])
    adni_fs6_datescanned[x, 32] <- temp_df$age_at_baseline[1] + interval(ymd(temp_df$baseline_date[1]), ymd(temp_df$EXAMDATE[timepoints])) %/% months(1)
  }
  x <- x+1
}

# age plot age at whole brain volume collection
adni_fs6_datescanned <- merge(adni_fs6_datescanned, 
                              data.frame(RID = adni_meth_datedrawn$RID, meth_timepoint1 = adni_meth_datedrawn$age_at_timepoint1),
                              by = 'RID')
adni_fs6_datescanned <- merge(adni_fs6_datescanned, adni_baseline_dates, by = 'RID')
adni_fs6_datescanned <- adni_fs6_datescanned[order(adni_fs6_datescanned$age_at_baseline),]




# pair with methylation

adni_fs4_datescanned <- data.frame(RID = unique(methMerge$RID), 
                                   DateScanned1 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned2 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned3 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned4 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned5 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned6 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned7 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned8 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned9 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned10 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned11 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned12 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned13 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned14 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned15 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned16 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned17 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned18 = rep(NA, length(unique(methMerge$RID))),
                                   
                                   
                                   age_at_scan1 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan2 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan3 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan4 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan5 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan6 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan7 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan8 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan9 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan10 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan11 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan12 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan13 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan14 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan15 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan16 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan17 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan18 = rep(NA, length(unique(methMerge$RID))),
                                   
                                   first_scan = rep(NA, length(unique(methMerge$RID))),
                                   last_scan = rep(NA, length(unique(methMerge$RID))),
                                   age_at_last_scan = rep(NA, length(unique(methMerge$RID))))

x <- 1
for (i in unique(methMerge$RID)){
  temp_df <- as.data.frame(methFS4_harmed[methFS4_harmed$RID == i,])
  temp_df <- temp_df[order(temp_df$EXAMDATE),]
  timepoints <- nrow(temp_df)
  if (timepoints > 0){
    for (t in 1:timepoints){
      adni_fs4_datescanned[x, (t+1)] <- as.character(temp_df$EXAMDATE[t])
      adni_fs4_datescanned[x, (t+19)] <- temp_df$age_at_baseline[1] + interval(ymd(temp_df$baseline_date[1]), ymd(temp_df$EXAMDATE[t])) %/% months(1)
    }
    adni_fs4_datescanned[x, 37] <- as.character(temp_df$EXAMDATE[1])
    adni_fs4_datescanned[x, 38] <- as.character(temp_df$EXAMDATE[timepoints])
    adni_fs4_datescanned[x, 39] <- temp_df$age_at_baseline[1] + interval(ymd(temp_df$baseline_date[1]), ymd(temp_df$EXAMDATE[timepoints])) %/% months(1)
  }
  x <- x+1
}




# pair fs5 with methylation

adni_fs5_datescanned <- data.frame(RID = unique(methMerge$RID), 
                                   DateScanned1 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned2 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned3 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned4 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned5 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned6 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned7 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned8 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned9 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned10 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned11 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned12 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned13 = rep(NA, length(unique(methMerge$RID))),
                                   DateScanned14 = rep(NA, length(unique(methMerge$RID))),
                                   
                                   age_at_scan1 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan2 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan3 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan4 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan5 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan6 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan7 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan8 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan9 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan10 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan11 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan12 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan13 = rep(NA, length(unique(methMerge$RID))),
                                   age_at_scan14 = rep(NA, length(unique(methMerge$RID))),
                                   
                                   first_scan = rep(NA, length(unique(methMerge$RID))),
                                   last_scan = rep(NA, length(unique(methMerge$RID))),
                                   age_at_last_scan = rep(NA, length(unique(methMerge$RID))))

x <- 1
for (i in unique(methMerge$RID)){
  temp_df <- as.data.frame(methFS5_harmed[methFS5_harmed$RID == i,])
  temp_df <- temp_df[order(temp_df$EXAMDATE),]
  timepoints <- nrow(temp_df)
  if (timepoints > 0){
    for (t in 1:timepoints){
      adni_fs5_datescanned[x, (t+1)] <- as.character(temp_df$EXAMDATE[t])
      adni_fs5_datescanned[x, (t+15)] <- temp_df$age_at_baseline[1] + interval(ymd(temp_df$baseline_date[1]), ymd(temp_df$EXAMDATE[t])) %/% months(1)
    }
    adni_fs5_datescanned[x, 30] <- as.character(temp_df$EXAMDATE[1])
    adni_fs5_datescanned[x, 31] <- as.character(temp_df$EXAMDATE[timepoints])
    adni_fs5_datescanned[x, 32] <- temp_df$age_at_baseline[1] + interval(ymd(temp_df$baseline_date[1]), ymd(temp_df$EXAMDATE[timepoints])) %/% months(1)
  }
  x <- x+1
}



#### pairing to meth-imaging observations within 6 months of eachother


imaging_observations <- list(tbv = adni_tbv_datescanned,
                             hc = adni_hc_datescanned,
                             wmh = adni_wmh_datescanned,
                             dx = adni_dx_datediagnosed,
                             fs5 = adni_fs5_datescanned, 
                             fs6 = adni_fs6_datescanned,
                             fs4 = adni_fs4_datescanned)


all_time_differences_mats <- list()
for (o in 1:length(imaging_observations)){
  temp_merged <- merge(adni_meth_datedrawn, imaging_observations[[o]], by = 'RID')
  all_time_differences_mats[[o]] <- array(rep(NA, 8*17*649), dim = c(8, 17, 649))
  for (m in 2:6){
    for (i in 18:34)
      if(grepl('Date', colnames(temp_merged[,i, drop = FALSE]))){
        all_time_differences_mats[[o]][(m-1), (i-17),] <- interval(mdy(temp_merged[,m]), ymd(temp_merged[,i])) %/% days(1)
      }
  }
}


# all imaging/meth timepoint pairs within and without 6mo

imaging_meth_obs_pairs <- list()
imaging_meth_obs_pairs_all <- list()
for (o in 1:length(all_time_differences_mats)){
  pairs <- data.frame(RID = 0, meth = 0, imaging = 0, time_diff = 0)
  all_pairs <- data.frame(RID = 0, meth = 0, imaging = 0, time_diff = 0)
  x <- 1
  y <- 1
  for (p in 1:649){
    person_mat <- all_time_differences_mats[[o]][,,p]
    rid_temp <- temp_merged$RID[p]
    for (m in 1:5){
      for (i in 1:17){
        if (abs(person_mat[m,i]) < 183 && !is.na(person_mat[m,i])){
          new_pair <- c(rid_temp, m, i, abs(person_mat[m,i]))
          pairs[x,] <- t(new_pair)
          x <- x+1
        }
        if (!is.na(person_mat[m,i])){
          new_pair <- c(rid_temp, m, i, abs(person_mat[m,i]))
          all_pairs[y,] <- t(new_pair)
          y <- y+1
        }
      }
    }
  }
  imaging_meth_obs_pairs[[o]] <- pairs
  imaging_meth_obs_pairs_all[[o]] <- all_pairs
}



# PAIRS WITHIN 6 MONTHS

length_temp <- rep(1, 7)
length_prev <- rep(0, 7)
count <- 1
while(sum(length_temp == length_prev) < length(imaging_meth_obs_pairs)){
  length_prev <- c(nrow(imaging_meth_obs_pairs[[1]]), 
                   nrow(imaging_meth_obs_pairs[[2]]), 
                   nrow(imaging_meth_obs_pairs[[3]]),
                   nrow(imaging_meth_obs_pairs[[4]]), 
                   nrow(imaging_meth_obs_pairs[[5]]),
                   nrow(imaging_meth_obs_pairs[[6]]),
                   nrow(imaging_meth_obs_pairs[[7]]))
  for (o in 1:length(imaging_meth_obs_pairs)){
    pairs_no_repeats <- data.frame(RID = 0, meth = 0, imaging = 0, time_diff = 0)
    pairs_temp <- imaging_meth_obs_pairs[[o]]
    for (s in unique(pairs_temp$RID)){
      person_temp <- pairs_temp[pairs_temp$RID == s,]
      ### duplicate meth
      if (length(unique(person_temp$meth)) != nrow(person_temp) && length(unique(person_temp$imaging)) == nrow(person_temp)){
        print(paste0('duplicate meth for ', s))
        duplicate <- person_temp[person_temp$meth == Mode(person_temp$meth),]
        competing_imaging <- duplicate$imaging
        
        temp_index <- abs(all_time_differences_mats[[o]][Mode(person_temp$meth),c(competing_imaging),which(sort(adni_meth_datedrawn$RID) == s)])
        
        if (length(unique(temp_index)) > 1){
          imaging_loser_timediff <- temp_index[which(temp_index != min(temp_index) )] 
          person_temp <- person_temp[-c(which(person_temp$time_diff %in% imaging_loser_timediff)),]
          pairs_no_repeats <- rbind(pairs_no_repeats, person_temp)
        } else if (length(unique(temp_index)) == 1){
          print(paste0('exact time diff in ', o, ' ', s))
          imaging_loser_dup_random <- competing_imaging[rep(c(2,1), 999)[count]]
          count <- count+1
          person_temp <- person_temp[-c(which(person_temp$imaging %in% imaging_loser_dup_random)),]
          pairs_no_repeats <- rbind(pairs_no_repeats, person_temp)
        }
        
        ## duplicate imaging
      } else if (length(unique(person_temp$imaging)) != nrow(person_temp) && length(unique(person_temp$meth)) == nrow(person_temp)){
        print(paste0('duplicate imaging for ', s))
        duplicate <- person_temp[person_temp$imaging == Mode(person_temp$imaging),]
        competing_meth <- duplicate$meth
        
        temp_index <- abs(all_time_differences_mats[[o]][c(competing_meth),Mode(person_temp$imaging),which(sort(adni_meth_datedrawn$RID) == s)])
        
        meth_loser_timediffs <- temp_index[which(temp_index != min(temp_index) )] 
        
        
        person_temp <- person_temp[-c(which(person_temp$time_diff %in% meth_loser_timediffs)),]
        pairs_no_repeats <- rbind(pairs_no_repeats, person_temp)
      } else if (length(unique(person_temp$meth)) != nrow(person_temp) && length(unique(person_temp$imaging)) != nrow(person_temp) ){
        # double duplicate
        print(c(o, s))
        print('double duplicate in imaging observation c(imaging observation, RID)')
        person_temp_new <- data.frame(RID = 0, meth = 0, imaging = 0, time_diff = 0)
        i <- 1
        for (m in unique(person_temp$meth)){
          person_temp_sub <- person_temp[person_temp$meth == m,]
          person_temp_sub <- person_temp_sub[order(person_temp_sub$time_diff),]
          person_temp_new[i,] <- person_temp_sub[1,]
          i <- i+1
        }
        pairs_no_repeats <- rbind(pairs_no_repeats, person_temp_new)
      } else if (length(unique(person_temp$meth)) == nrow(person_temp) && length(unique(person_temp$imaging)) == nrow(person_temp)){
        pairs_no_repeats <- rbind(pairs_no_repeats, person_temp)
      }
    }
    pairs_no_repeats <- pairs_no_repeats[2:nrow(pairs_no_repeats),]
    imaging_meth_obs_pairs[[o]] <- pairs_no_repeats
    length_temp[o] <- nrow(imaging_meth_obs_pairs[[o]])
  }
}





# comparing freesurfer versions
imaging_meth_obs_pairs[[5]] # fs5
imaging_meth_obs_pairs[[6]] # fs6
imaging_meth_obs_pairs[[7]] # fs4




#### pulling data from across sources and running associations




############# TBV ############# 

obs_temp <- imaging_meth_obs_pairs[[1]]
tbv_paired_values <- data.frame(RID = 0, 
                                zPoAm45 = 0,
                                zHorvath = 0,
                                zHannum = 0,
                                zPheno = 0,
                                zGrim = 0,
                                datedrawn = 0, age_at_meth = 0, TBV = 0, datescanned = 0, age_at_scan = 0, ICV = 0, DX = 0,
                                CD4T = 0,
                                NK = 0,
                                Mono = 0,
                                Gran = 0,
                                PlasmaBlast = 0,
                                CD8pCD28nCD45RAn = 0,
                                CD8.naive = 0,
                                APOE4 = 0)
for (i in 1:nrow(imaging_meth_obs_pairs[[1]])){
  obs_pair <- data.frame(RID = 0, zPoAm45 = 0, zHorvath = 0, zHannum = 0, zPheno = 0, zGrim = 0, datedrawn = 0, age_at_meth = 0, TBV = 0, datescanned = 0, age_at_scan = 0, ICV = 0, DX = 0,
                         CD4T = 0, NK = 0, Mono = 0, Gran = 0, PlasmaBlast = 0, CD8pCD28nCD45RAn = 0, CD8.naive = 0, APOE4 = 0)
  obs_pair_index_temp <- obs_temp[i,]
  
  meth_temp <- adni_methclocks[adni_methclocks$RID == obs_pair_index_temp$RID,]
  meth_temp <- meth_temp[order(mdy(meth_temp$Edate)),]
  
  age_at_meth <- adni_meth[adni_meth$RID == obs_pair_index_temp$RID,]
  age_at_meth <- age_at_meth[order(mdy(age_at_meth$Edate)),]
  
  meth_index_temp <- obs_pair_index_temp$meth
  
  obs_pair[,c('RID', 'zPoAm45', 'zHorvath', 'zHannum', 'zPheno', 'zGrim', 'datedrawn', 'age_at_meth',
              'CD4T', 'NK', 'Mono', 'Gran', 'PlasmaBlast', 'CD8pCD28nCD45RAn', 'CD8.naive')] <- c(meth_temp$RID[meth_index_temp], 
                                                                                                  meth_temp$zPoAm45[meth_index_temp], 
                                                                                                  meth_temp$zHorvath[meth_index_temp],
                                                                                                  meth_temp$zHannum[meth_index_temp],
                                                                                                  meth_temp$zPheno[meth_index_temp],
                                                                                                  meth_temp$zGrim[meth_index_temp],
                                                                                                  meth_temp$EXAMDATE[meth_index_temp], 
                                                                                                  age_at_meth$AGE_MONTH[meth_index_temp] + interval(ymd(age_at_meth$baseline_date[meth_index_temp]),  mdy(meth_temp$Edate[meth_index_temp])) %/% months(1),
                                                                                                  
                                                                                                  meth_temp$CD4T[meth_index_temp],
                                                                                                  meth_temp$NK[meth_index_temp],
                                                                                                  meth_temp$Mono[meth_index_temp],
                                                                                                  meth_temp$Gran[meth_index_temp],
                                                                                                  meth_temp$PlasmaBlast[meth_index_temp],
                                                                                                  meth_temp$CD8pCD28nCD45RAn[meth_index_temp],
                                                                                                  meth_temp$CD8.naive[meth_index_temp]
              )
  
  
  imaging_temp <- methMerge[methMerge$RID == obs_pair_index_temp$RID,]
  imaging_temp <- imaging_temp[!is.na(imaging_temp$WholeBrain),]
  imaging_temp <- imaging_temp[order(ymd(imaging_temp$EXAMDATE)),]
  
  imaging_index_temp <- obs_pair_index_temp$imaging
  
  obs_pair[,c('TBV', 'datescanned', 'age_at_scan', 'ICV', 'DX', 'APOE4')] <- c(imaging_temp$WholeBrain[imaging_index_temp], 
                                                                               as.character(imaging_temp$EXAMDATE[imaging_index_temp]), 
                                                                               imaging_temp$age_at_baseline[imaging_index_temp] + interval(ymd(imaging_temp$baseline_date[imaging_index_temp]),  ymd(imaging_temp$EXAMDATE[imaging_index_temp])) %/% months(1),
                                                                               imaging_temp$ICV[imaging_index_temp],
                                                                               imaging_temp$DX[imaging_index_temp],
                                                                               imaging_temp$APOE4[imaging_index_temp])
  
  
  
  tbv_paired_values <- rbind(tbv_paired_values, obs_pair)
}
tbv_paired_values <- tbv_paired_values[2:nrow(tbv_paired_values),]

tbv_paired_values$mean_age <- scale(rowMeans(data.frame(as.numeric(tbv_paired_values$age_at_meth), as.numeric(tbv_paired_values$age_at_scan))))
tbv_paired_values$zPoAm45 <- as.numeric(tbv_paired_values$zPoAm45)
tbv_paired_values$zHorvath <- as.numeric(tbv_paired_values$zHorvath)
tbv_paired_values$zHannum <- as.numeric(tbv_paired_values$zHannum)
tbv_paired_values$zPheno <- as.numeric(tbv_paired_values$zPheno)
tbv_paired_values$zGrim <- as.numeric(tbv_paired_values$zGrim)

tbv_paired_values$TBV <- as.numeric(tbv_paired_values$TBV)
tbv_paired_values$ICV <- as.numeric(tbv_paired_values$ICV)
adni_sex <- data.frame(RID = adnimerge$RID, sex = adnimerge$PTGENDER)
adni_sex <- adni_sex[!duplicated(adni_sex),]
tbv_paired_values <- merge(adni_sex, tbv_paired_values, by = 'RID')
tbv_paired_values$sex <- as.factor(tbv_paired_values$sex)
tbv_paired_values$RID <- as.factor(tbv_paired_values$RID)


# ultimately decided that I am going to use Annchen's versions because there is some weird thing with the different date labels that isn't lining up exactly

## checking against Annchen's versions

load('/Users/ew198/Documents/methylation/data/adni/annchen/meth_merged.Rdata')
TBV_subset_ak <- meth_merged %>% filter(!is.na(WholeBrain_val)) # discrepancy is 4136
tbv_paired_values <- TBV_subset_ak

tbv.p.dat <- pdata.frame(data.frame(RID = as.character(TBV_subset_ak$RID),
                                    TBV = scale(TBV_subset_ak$WholeBrain_val), 
                                    zPoAm45=scale(TBV_subset_ak$zPoAm45), 
                                    zHorvath=scale(TBV_subset_ak$zHorvath),
                                    zHannum=scale(TBV_subset_ak$zHannum),
                                    zPheno=scale(TBV_subset_ak$zPheno),
                                    zGrim=scale(TBV_subset_ak$zGrim),
                                    sex = TBV_subset_ak$PTGENDER.x,
                                    mean_age = TBV_subset_ak$WholeBrain_AGE.avg,
                                    mean_age_sq = TBV_subset_ak$WholeBrain_AGE.avg^2,
                                    ICV = scale(TBV_subset_ak$WholeBrainICV),
                                    DX = TBV_subset_ak$DX_val,
                                    datescanned = TBV_subset_ak$WholeBrain_date, 
                                    datedrawn = TBV_subset_ak$DateDrawn_asDate),
                         index = c("RID"), drop.index = F, row.names = T)


### TBV controlling for ICV

#dunedinpace
pm1 <- plm(TBV~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex+ICV, model = "pooling", data = tbv.p.dat, na.action = na.omit)
G <- length(unique(tbv.p.dat$RID))
N <- length(tbv.p.dat$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
tbvrel.dp.plm.res <- coeftest(pm1, vcov = firm_c_vcov)

#horvath
pm1 <- plm(TBV~zHorvath+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex+ICV, model = "pooling", data = tbv.p.dat, na.action = na.omit)
G <- length(unique(tbv_paired_values$RID))
N <- length(tbv_paired_values$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
tbvrel.hth.plm.res <- coeftest(pm1, vcov = firm_c_vcov)

#hannum
pm1 <- plm(TBV~zHannum+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex+ICV, model = "pooling", data = tbv.p.dat, na.action = na.omit)
G <- length(unique(tbv_paired_values$RID))
N <- length(tbv_paired_values$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
tbvrel.hnm.plm.res <- coeftest(pm1, vcov = firm_c_vcov)

#phenoage
pm1 <- plm(TBV~zPheno+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex+ICV, model = "pooling", data = tbv.p.dat, na.action = na.omit)
G <- length(unique(tbv_paired_values$RID))
N <- length(tbv_paired_values$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
tbvrel.phe.plm.res <- coeftest(pm1, vcov = firm_c_vcov)

#grimage
pm1 <- plm(TBV~zGrim+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex+ICV, model = "pooling", data = tbv.p.dat, na.action = na.omit)
G <- length(unique(tbv_paired_values$RID))
N <- length(tbv_paired_values$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
tbvrel.grim.plm.res <- coeftest(pm1, vcov = firm_c_vcov)


############# HC ############# 


obs_temp <- imaging_meth_obs_pairs[[2]]
hc_paired_values <- data.frame(RID = 0, 
                               zPoAm45 = 0,
                               zHorvath = 0,
                               zHannum = 0,
                               zPheno = 0,
                               zGrim = 0,
                               datedrawn = 0, age_at_meth = 0, HC = 0, datescanned = 0, age_at_scan = 0, ICV = 0, TBV = 0, DX = 0,
                               CD4T = 0,
                               NK = 0,
                               Mono = 0,
                               Gran = 0,
                               PlasmaBlast = 0,
                               CD8pCD28nCD45RAn = 0,
                               CD8.naive = 0,
                               APOE4 = 0
)
for (i in 1:nrow(imaging_meth_obs_pairs[[2]])){
  obs_pair <- data.frame(RID = 0, zPoAm45 = 0, zHorvath = 0, zHannum = 0, zPheno = 0, zGrim = 0, datedrawn = 0, age_at_meth = 0, HC = 0, datescanned = 0, age_at_scan = 0, ICV = 0, TBV = 0, DX = 0,
                         CD4T = 0, NK = 0, Mono = 0, Gran = 0, PlasmaBlast = 0, CD8pCD28nCD45RAn = 0, CD8.naive = 0, APOE4 = 0)
  obs_pair_index_temp <- obs_temp[i,]
  
  meth_temp <- adni_methclocks[adni_methclocks$RID == obs_pair_index_temp$RID,]
  meth_temp <- meth_temp[order(mdy(meth_temp$Edate)),]
  
  age_at_meth <- adni_meth[adni_meth$RID == obs_pair_index_temp$RID,]
  age_at_meth <- age_at_meth[order(mdy(age_at_meth$Edate)),]
  
  meth_index_temp <- obs_pair_index_temp$meth
  
  obs_pair[,c('RID', 'zPoAm45', 'zHorvath', 'zHannum', 'zPheno', 'zGrim', 'datedrawn', 'age_at_meth',
              'CD4T', 'NK', 'Mono', 'Gran', 'PlasmaBlast', 'CD8pCD28nCD45RAn', 'CD8.naive')] <- c(meth_temp$RID[meth_index_temp], 
                                                                                                  meth_temp$zPoAm45[meth_index_temp], 
                                                                                                  meth_temp$zHorvath[meth_index_temp],
                                                                                                  meth_temp$zHannum[meth_index_temp],
                                                                                                  meth_temp$zPheno[meth_index_temp],
                                                                                                  meth_temp$zGrim[meth_index_temp],
                                                                                                  meth_temp$EXAMDATE[meth_index_temp], 
                                                                                                  age_at_meth$AGE_MONTH[meth_index_temp] + interval(ymd(age_at_meth$baseline_date[meth_index_temp]),  mdy(meth_temp$Edate[meth_index_temp])) %/% months(1),
                                                                                                  
                                                                                                  meth_temp$CD4T[meth_index_temp],
                                                                                                  meth_temp$NK[meth_index_temp],
                                                                                                  meth_temp$Mono[meth_index_temp],
                                                                                                  meth_temp$Gran[meth_index_temp],
                                                                                                  meth_temp$PlasmaBlast[meth_index_temp],
                                                                                                  meth_temp$CD8pCD28nCD45RAn[meth_index_temp],
                                                                                                  meth_temp$CD8.naive[meth_index_temp]
              )
  
  
  
  
  imaging_temp <- methMerge[methMerge$RID == obs_pair_index_temp$RID,]
  imaging_temp <- imaging_temp[!is.na(imaging_temp$Hippocampus),]
  imaging_temp <- imaging_temp[order(ymd(imaging_temp$EXAMDATE)),]
  
  imaging_index_temp <- obs_pair_index_temp$imaging
  
  obs_pair[,c('HC', 'datescanned', 'age_at_scan', 'ICV', 'TBV', 'DX', 'APOE4')] <- c(imaging_temp$Hippocampus[imaging_index_temp], 
                                                                                     as.character(imaging_temp$EXAMDATE[imaging_index_temp]), 
                                                                                     imaging_temp$age_at_baseline[imaging_index_temp] + interval(ymd(imaging_temp$baseline_date[imaging_index_temp]),  ymd(imaging_temp$EXAMDATE[imaging_index_temp])) %/% months(1),
                                                                                     imaging_temp$ICV[imaging_index_temp],
                                                                                     imaging_temp$WholeBrain[imaging_index_temp],
                                                                                     imaging_temp$DX[imaging_index_temp],
                                                                                     imaging_temp$APOE4[imaging_index_temp]
  )
  
  
  
  hc_paired_values <- rbind(hc_paired_values, obs_pair)
}
hc_paired_values <- hc_paired_values[2:nrow(hc_paired_values),]

hc_paired_values$mean_age <- scale(rowMeans(data.frame(as.numeric(hc_paired_values$age_at_meth), as.numeric(hc_paired_values$age_at_scan))))
hc_paired_values$zPoAm45 <- as.numeric(hc_paired_values$zPoAm45)
hc_paired_values$zHorvath <- as.numeric(hc_paired_values$zHorvath)
hc_paired_values$zHannum <- as.numeric(hc_paired_values$zHannum)
hc_paired_values$zPheno <- as.numeric(hc_paired_values$zPheno)
hc_paired_values$zGrim <- as.numeric(hc_paired_values$zGrim)

hc_paired_values$HC <- as.numeric(hc_paired_values$HC)
hc_paired_values$ICV <- as.numeric(hc_paired_values$ICV)
hc_paired_values$TBV <- as.numeric(hc_paired_values$TBV)
adni_sex <- data.frame(RID = adnimerge$RID, sex = adnimerge$PTGENDER)
adni_sex <- adni_sex[!duplicated(adni_sex),]
hc_paired_values <- merge(adni_sex, hc_paired_values, by = 'RID')
hc_paired_values$sex <- as.factor(hc_paired_values$sex)
hc_paired_values$RID <- as.factor(hc_paired_values$RID)


load('/Users/ew198/Documents/methylation/data/adni/annchen/meth_merged.Rdata')
HC_subset_ak <- meth_merged %>% filter(!is.na(Hippocampus_val)) # discrepancy is 4136
hc_paired_values <- HC_subset_ak

hc.p.dat <- pdata.frame(data.frame(RID = as.character(HC_subset_ak$RID),
                                    HC = scale(HC_subset_ak$Hippocampus_val), 
                                    zPoAm45=scale(HC_subset_ak$zPoAm45), 
                                    zHorvath=scale(HC_subset_ak$zHorvath),
                                    zHannum=scale(HC_subset_ak$zHannum),
                                    zPheno=scale(HC_subset_ak$zPheno),
                                    zGrim=scale(HC_subset_ak$zGrim),
                                    sex = HC_subset_ak$PTGENDER.x,
                                    mean_age = HC_subset_ak$Hippocampus_AGE.avg,
                                    mean_age_sq = HC_subset_ak$Hippocampus_AGE.avg^2,
                                    ICV = scale(HC_subset_ak$HippocampusICV),
                                    DX = HC_subset_ak$DX_val,
                                    datescanned = HC_subset_ak$Hippocampus_date, 
                                    datedrawn = HC_subset_ak$DateDrawn_asDate),
                         index = c("RID"), drop.index = F, row.names = T)


### relative HC (controlling for ICV)

#dunedinpace
pm1 <- plm(HC~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex+ICV, model = "pooling", data = hc.p.dat, na.action = na.omit)
G <- length(unique(hc.p.dat$RID))
N <- length(hc.p.dat$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
hcrel.dp.plm.res <- coeftest(pm1, vcov = firm_c_vcov)

#horvath
pm1 <- plm(HC~zHorvath+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex+ICV, model = "pooling", data = hc.p.dat, na.action = na.omit)
G <- length(unique(hc_paired_values$RID))
N <- length(hc_paired_values$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
hcrel.hth.plm.res <- coeftest(pm1, vcov = firm_c_vcov)

#hannum
pm1 <- plm(HC~zHannum+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex+ICV, model = "pooling", data = hc.p.dat, na.action = na.omit)
G <- length(unique(hc_paired_values$RID))
N <- length(hc_paired_values$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
hcrel.hnm.plm.res <- coeftest(pm1, vcov = firm_c_vcov)

#phenoage
pm1 <- plm(HC~zPheno+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex+ICV, model = "pooling", data = hc.p.dat, na.action = na.omit)
G <- length(unique(hc_paired_values$RID))
N <- length(hc_paired_values$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
hcrel.phe.plm.res <- coeftest(pm1, vcov = firm_c_vcov)

#grimage
pm1 <- plm(HC~zGrim+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex+ICV, model = "pooling", data = hc.p.dat, na.action = na.omit)
G <- length(unique(hc_paired_values$RID))
N <- length(hc_paired_values$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
hcrel.grim.plm.res <- coeftest(pm1, vcov = firm_c_vcov)




############# WMH ############# 

obs_temp <- imaging_meth_obs_pairs[[3]]
wmh_paired_values <- data.frame(RID = 0, 
                                zPoAm45 = 0,
                                zHorvath = 0,
                                zHannum = 0,
                                zPheno = 0,
                                zGrim = 0,
                                datedrawn = 0, age_at_meth = 0, wmh = 0, datescanned = 0, age_at_scan = 0, TBV = 0, WMV = 0, DX = 0,
                                CD4T = 0,
                                NK = 0,
                                Mono = 0,
                                Gran = 0,
                                PlasmaBlast = 0,
                                CD8pCD28nCD45RAn = 0,
                                CD8.naive = 0)
for (i in 1:nrow(imaging_meth_obs_pairs[[3]])){
  obs_pair <- data.frame(RID = 0, zPoAm45 = 0, zHorvath = 0, zHannum = 0, zPheno = 0, zGrim = 0, datedrawn = 0, age_at_meth = 0, wmh = 0, datescanned = 0, age_at_scan = 0, TBV = 0, WMV = 0, DX = 0,
                         CD4T = 0, NK = 0, Mono = 0, Gran = 0, PlasmaBlast = 0, CD8pCD28nCD45RAn = 0, CD8.naive = 0)
  obs_pair_index_temp <- obs_temp[i,]
  
  meth_temp <- adni_methclocks[adni_methclocks$RID == obs_pair_index_temp$RID,]
  meth_temp <- meth_temp[order(mdy(meth_temp$Edate)),]
  
  age_at_meth <- adni_meth[adni_meth$RID == obs_pair_index_temp$RID,]
  age_at_meth <- age_at_meth[order(mdy(age_at_meth$Edate)),]
  
  meth_index_temp <- obs_pair_index_temp$meth
  
  obs_pair[,c('RID', 'zPoAm45', 'zHorvath', 'zHannum', 'zPheno', 'zGrim', 'datedrawn', 'age_at_meth',
              'CD4T', 'NK', 'Mono', 'Gran', 'PlasmaBlast', 'CD8pCD28nCD45RAn', 'CD8.naive')] <- c(meth_temp$RID[meth_index_temp], 
                                                                                                  meth_temp$zPoAm45[meth_index_temp], 
                                                                                                  meth_temp$zHorvath[meth_index_temp],
                                                                                                  meth_temp$zHannum[meth_index_temp],
                                                                                                  meth_temp$zPheno[meth_index_temp],
                                                                                                  meth_temp$zGrim[meth_index_temp],
                                                                                                  meth_temp$EXAMDATE[meth_index_temp], 
                                                                                                  age_at_meth$AGE_MONTH[meth_index_temp] + interval(ymd(age_at_meth$baseline_date[meth_index_temp]),  mdy(meth_temp$Edate[meth_index_temp])) %/% months(1),
                                                                                                  
                                                                                                  meth_temp$CD4T[meth_index_temp],
                                                                                                  meth_temp$NK[meth_index_temp],
                                                                                                  meth_temp$Mono[meth_index_temp],
                                                                                                  meth_temp$Gran[meth_index_temp],
                                                                                                  meth_temp$PlasmaBlast[meth_index_temp],
                                                                                                  meth_temp$CD8pCD28nCD45RAn[meth_index_temp],
                                                                                                  meth_temp$CD8.naive[meth_index_temp]
              )
  
  
  
  
  imaging_temp <- wmh_merge[wmh_merge$RID == obs_pair_index_temp$RID,]
  imaging_temp <- imaging_temp[!is.na(imaging_temp$wmhVol_log),]
  imaging_temp <- imaging_temp[order(ymd(imaging_temp$EXAMDATE)),]
  
  imaging_index_temp <- obs_pair_index_temp$imaging
  
  obs_pair[,c('wmh', 'datescanned', 'age_at_scan', 'TBV')] <- c(imaging_temp$wmhVol_log[imaging_index_temp], 
                                                                as.character(imaging_temp$EXAMDATE[imaging_index_temp]), 
                                                                imaging_temp$age_at_baseline[imaging_index_temp] + interval(ymd(imaging_temp$baseline_date[imaging_index_temp]),  ymd(imaging_temp$EXAMDATE[imaging_index_temp])) %/% months(1),
                                                                imaging_temp$TOTAL_BRAIN[imaging_index_temp]
  )
  
  
  
  wmh_paired_values <- rbind(wmh_paired_values, obs_pair)
}
wmh_paired_values <- wmh_paired_values[2:nrow(wmh_paired_values),]

wmh_paired_values$mean_age <- scale(rowMeans(data.frame(as.numeric(wmh_paired_values$age_at_meth), as.numeric(wmh_paired_values$age_at_scan))))
wmh_paired_values$zPoAm45 <- as.numeric(wmh_paired_values$zPoAm45)
wmh_paired_values$zHorvath <- as.numeric(wmh_paired_values$zHorvath)
wmh_paired_values$zHannum <- as.numeric(wmh_paired_values$zHannum)
wmh_paired_values$zPheno <- as.numeric(wmh_paired_values$zPheno)
wmh_paired_values$zGrim <- as.numeric(wmh_paired_values$zGrim)

wmh_paired_values$wmh <- as.numeric(wmh_paired_values$wmh)
wmh_paired_values$TBV <- as.numeric(wmh_paired_values$TBV)
wmh_paired_values$WMV <- as.numeric(wmh_paired_values$WMV)
adni_sex <- data.frame(RID = adnimerge$RID, sex = adnimerge$PTGENDER)
adni_sex <- adni_sex[!duplicated(adni_sex),]
wmh_paired_values <- merge(adni_sex, wmh_paired_values, by = 'RID')
wmh_paired_values$sex <- as.factor(wmh_paired_values$sex)
wmh_paired_values$RID <- as.factor(wmh_paired_values$RID)


load('/Users/ew198/Documents/methylation/data/adni/annchen/meth_merged.Rdata')
WMH_subset_ak <- meth_merged %>% filter(!is.na(logWMH_val)) # discrepancy is 4136
wmh_paired_values <- WMH_subset_ak

wmh.p.dat <- pdata.frame(data.frame(RID = as.character(WMH_subset_ak$RID),
                                   wmh = scale(WMH_subset_ak$logWMH_val), 
                                   zPoAm45=scale(WMH_subset_ak$zPoAm45), 
                                   zHorvath=scale(WMH_subset_ak$zHorvath),
                                   zHannum=scale(WMH_subset_ak$zHannum),
                                   zPheno=scale(WMH_subset_ak$zPheno),
                                   zGrim=scale(WMH_subset_ak$zGrim),
                                   sex = WMH_subset_ak$PTGENDER.x,
                                   mean_age = WMH_subset_ak$logWMH_AGE.avg,
                                   mean_age_sq = WMH_subset_ak$logWMH_AGE.avg^2,
                                   DX = WMH_subset_ak$DX_val,
                                   datescanned = WMH_subset_ak$logWMH_date, 
                                   datedrawn = WMH_subset_ak$DateDrawn_asDate),
                        index = c("RID"), drop.index = F, row.names = T)

#dunedinpace
pm1 <- plm(wmh~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = wmh.p.dat, na.action = na.omit)
G <- length(unique(wmh.p.dat$RID))
N <- length(wmh.p.dat$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmh.dp.plm.res <- coeftest(pm1, vcov = firm_c_vcov)

#horvath
pm1 <- plm(wmh~zHorvath+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = wmh.p.dat, na.action = na.omit)
G <- length(unique(wmh_paired_values$RID))
N <- length(wmh_paired_values$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmh.hth.plm.res <- coeftest(pm1, vcov = firm_c_vcov)

#hannum
pm1 <- plm(wmh~zHannum+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = wmh.p.dat, na.action = na.omit)
G <- length(unique(wmh_paired_values$RID))
N <- length(wmh_paired_values$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmh.hnm.plm.res <- coeftest(pm1, vcov = firm_c_vcov)

#phenoage
pm1 <- plm(wmh~zPheno+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = wmh.p.dat, na.action = na.omit)
G <- length(unique(wmh_paired_values$RID))
N <- length(wmh_paired_values$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmh.phe.plm.res <- coeftest(pm1, vcov = firm_c_vcov)

#grimage
pm1 <- plm(wmh~zGrim+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = wmh.p.dat, na.action = na.omit)
G <- length(unique(wmh_paired_values$RID))
N <- length(wmh_paired_values$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmh.grim.plm.res <- coeftest(pm1, vcov = firm_c_vcov)






# FS5

obs_temp <- imaging_meth_obs_pairs[[5]]
fs5_paired_values <- data.frame(RID = 0, 
                                zPoAm45 = 0,
                                zHorvath = 0,
                                zHannum = 0,
                                zPheno = 0,
                                zGrim = 0,
                                datedrawn = 0, age_at_meth = 0, QC = 0, total_sa = 0, mean_thickness = 0, total_sa.combat=0, mean_thickness.combat=0, 
                                datescanned = 0, age_at_scan = 0, DX = 0,
                                CD4T = 0,
                                NK = 0,
                                Mono = 0,
                                Gran = 0,
                                PlasmaBlast = 0,
                                CD8pCD28nCD45RAn = 0,
                                CD8.naive = 0)
for (i in 1:nrow(imaging_meth_obs_pairs[[5]])){
  obs_pair <- data.frame(RID = 0, zPoAm45 = 0, zHorvath = 0, zHannum = 0, zPheno = 0, zGrim = 0, datedrawn = 0, age_at_meth = 0, QC = 0, total_sa = 0, mean_thickness = 0, datescanned = 0, age_at_scan = 0, DX = 0,
                         CD4T = 0, NK = 0, Mono = 0, Gran = 0, PlasmaBlast = 0, CD8pCD28nCD45RAn = 0, CD8.naive = 0)
  obs_pair_index_temp <- obs_temp[i,]
  
  meth_temp <- adni_methclocks[adni_methclocks$RID == obs_pair_index_temp$RID,]
  meth_temp <- meth_temp[order(mdy(meth_temp$Edate)),]
  
  #adni_meth <- merge(adni_meth, adni_baseline_dates, by = 'RID')
  age_at_meth <- adni_meth[adni_meth$RID == obs_pair_index_temp$RID,]
  age_at_meth <- age_at_meth[order(mdy(age_at_meth$Edate)),]
  
  meth_index_temp <- obs_pair_index_temp$meth
  
  obs_pair[,c('RID', 'zPoAm45', 'zHorvath', 'zHannum', 'zPheno', 'zGrim', 'datedrawn', 'age_at_meth',
              'CD4T', 'NK', 'Mono', 'Gran', 'PlasmaBlast', 'CD8pCD28nCD45RAn', 'CD8.naive')] <- c(meth_temp$RID[meth_index_temp], 
                                                                                                  meth_temp$zPoAm45[meth_index_temp], 
                                                                                                  meth_temp$zHorvath[meth_index_temp],
                                                                                                  meth_temp$zHannum[meth_index_temp],
                                                                                                  meth_temp$zPheno[meth_index_temp],
                                                                                                  meth_temp$zGrim[meth_index_temp],
                                                                                                  meth_temp$EXAMDATE[meth_index_temp], 
                                                                                                  age_at_meth$AGE_MONTH[meth_index_temp] + interval(ymd(age_at_meth$baseline_date[meth_index_temp]),  mdy(meth_temp$Edate[meth_index_temp])) %/% months(1),
                                                                                                  
                                                                                                  meth_temp$CD4T[meth_index_temp],
                                                                                                  meth_temp$NK[meth_index_temp],
                                                                                                  meth_temp$Mono[meth_index_temp],
                                                                                                  meth_temp$Gran[meth_index_temp],
                                                                                                  meth_temp$PlasmaBlast[meth_index_temp],
                                                                                                  meth_temp$CD8pCD28nCD45RAn[meth_index_temp],
                                                                                                  meth_temp$CD8.naive[meth_index_temp]
              )
  
  
  imaging_temp <- methFS5_harmed[methFS5_harmed$RID == obs_pair_index_temp$RID,]
  imaging_temp <- imaging_temp[!is.na(imaging_temp$total_sa),]
  imaging_temp <- imaging_temp[order(ymd(imaging_temp$EXAMDATE)),]
  
  imaging_index_temp <- obs_pair_index_temp$imaging
  
  obs_pair[,c('QC', 'total_sa', 'mean_thickness', 'total_sa.combat', 'mean_thickness.combat', 'datescanned', 'age_at_scan')] <- c(imaging_temp$OVERALLQC[imaging_index_temp],
                                                                                                                                  imaging_temp$total_sa[imaging_index_temp], 
                                                                                                                                  imaging_temp$mean_thickness[imaging_index_temp],
                                                                                                                                  imaging_temp$total_sa.combat[imaging_index_temp], 
                                                                                                                                  imaging_temp$mean_thickness.combat[imaging_index_temp],
                                                                                                                                  as.character(imaging_temp$EXAMDATE[imaging_index_temp]), 
                                                                                                                                  imaging_temp$age_at_baseline[imaging_index_temp] + interval(ymd(imaging_temp$baseline_date[imaging_index_temp]),  ymd(imaging_temp$EXAMDATE[imaging_index_temp])) %/% months(1)
                                                                                                                                  # imaging_temp$DX[imaging_index_temp],
                                                                                                                                  # imaging_temp$APOE4[imaging_index_temp]
  )
  
  
  
  fs5_paired_values <- rbind(fs5_paired_values, obs_pair)
}
fs5_paired_values <- fs5_paired_values[2:nrow(fs5_paired_values),]

fs5_paired_values$mean_age <- rowMeans(data.frame(as.numeric(fs5_paired_values$age_at_meth), as.numeric(fs5_paired_values$age_at_scan)))
fs5_paired_values$zPoAm45 <- as.numeric(fs5_paired_values$zPoAm45)
fs5_paired_values$zHorvath <- as.numeric(fs5_paired_values$zHorvath)
fs5_paired_values$zHannum <- as.numeric(fs5_paired_values$zHannum)
fs5_paired_values$zPheno <- as.numeric(fs5_paired_values$zPheno)
fs5_paired_values$zGrim <- as.numeric(fs5_paired_values$zGrim)
fs5_paired_values$mean_thickness <- as.numeric(fs5_paired_values$mean_thickness)
fs5_paired_values$total_sa <- as.numeric(fs5_paired_values$total_sa)
fs5_paired_values$mean_thickness.combat <- as.numeric(fs5_paired_values$mean_thickness.combat)
fs5_paired_values$total_sa.combat <- as.numeric(fs5_paired_values$total_sa.combat)

adni_sex <- data.frame(RID = adnimerge$RID, sex = adnimerge$PTGENDER)
adni_sex <- adni_sex[!duplicated(adni_sex),]
fs5_paired_values <- merge(adni_sex, fs5_paired_values, by = 'RID')
fs5_paired_values$sex <- as.factor(fs5_paired_values$sex)
fs5_paired_values$RID <- as.factor(fs5_paired_values$RID)



# FS4

obs_temp <- imaging_meth_obs_pairs[[7]]
fs4_paired_values <- data.frame(RID = 0, 
                                zPoAm45 = 0,
                                zHorvath = 0,
                                zHannum = 0,
                                zPheno = 0,
                                zGrim = 0,
                                datedrawn = 0, age_at_meth = 0, QC = 0, total_sa = 0, mean_thickness = 0, total_sa.combat=0, mean_thickness.combat=0,
                                datescanned = 0, age_at_scan = 0, DX = 0,
                                CD4T = 0,
                                NK = 0,
                                Mono = 0,
                                Gran = 0,
                                PlasmaBlast = 0,
                                CD8pCD28nCD45RAn = 0,
                                CD8.naive = 0)
for (i in 1:nrow(imaging_meth_obs_pairs[[7]])){
  obs_pair <- data.frame(RID = 0, zPoAm45 = 0, zHorvath = 0, zHannum = 0, zPheno = 0, zGrim = 0, datedrawn = 0, age_at_meth = 0, QC = 0, total_sa = 0, mean_thickness = 0, datescanned = 0, age_at_scan = 0, DX = 0,
                         CD4T = 0, NK = 0, Mono = 0, Gran = 0, PlasmaBlast = 0, CD8pCD28nCD45RAn = 0, CD8.naive = 0)
  obs_pair_index_temp <- obs_temp[i,]
  
  meth_temp <- adni_methclocks[adni_methclocks$RID == obs_pair_index_temp$RID,]
  meth_temp <- meth_temp[order(mdy(meth_temp$Edate)),]
  
  age_at_meth <- adni_meth[adni_meth$RID == obs_pair_index_temp$RID,]
  age_at_meth <- age_at_meth[order(mdy(age_at_meth$Edate)),]
  
  meth_index_temp <- obs_pair_index_temp$meth
  
  obs_pair[,c('RID', 'zPoAm45', 'zHorvath', 'zHannum', 'zPheno', 'zGrim', 'datedrawn', 'age_at_meth',
              'CD4T', 'NK', 'Mono', 'Gran', 'PlasmaBlast', 'CD8pCD28nCD45RAn', 'CD8.naive')] <- c(meth_temp$RID[meth_index_temp], 
                                                                                                  meth_temp$zPoAm45[meth_index_temp], 
                                                                                                  meth_temp$zHorvath[meth_index_temp],
                                                                                                  meth_temp$zHannum[meth_index_temp],
                                                                                                  meth_temp$zPheno[meth_index_temp],
                                                                                                  meth_temp$zGrim[meth_index_temp],
                                                                                                  meth_temp$EXAMDATE[meth_index_temp], 
                                                                                                  age_at_meth$AGE_MONTH[meth_index_temp] + interval(ymd(age_at_meth$baseline_date[meth_index_temp]),  mdy(meth_temp$Edate[meth_index_temp])) %/% months(1),
                                                                                                  
                                                                                                  meth_temp$CD4T[meth_index_temp],
                                                                                                  meth_temp$NK[meth_index_temp],
                                                                                                  meth_temp$Mono[meth_index_temp],
                                                                                                  meth_temp$Gran[meth_index_temp],
                                                                                                  meth_temp$PlasmaBlast[meth_index_temp],
                                                                                                  meth_temp$CD8pCD28nCD45RAn[meth_index_temp],
                                                                                                  meth_temp$CD8.naive[meth_index_temp]
              )
  
  
  imaging_temp <- methFS4_harmed[methFS4_harmed$RID == obs_pair_index_temp$RID,]
  imaging_temp <- imaging_temp[order(ymd(imaging_temp$EXAMDATE)),]
  
  imaging_index_temp <- obs_pair_index_temp$imaging
  
  obs_pair[,c('QC', 'total_sa', 'mean_thickness', 'total_sa.combat', 'mean_thickness.combat', 'datescanned', 'age_at_scan')] <- c(imaging_temp$OVERALLQC[imaging_index_temp],
                                                                                                                                  imaging_temp$total_sa[imaging_index_temp], 
                                                                                                                                  imaging_temp$mean_thickness[imaging_index_temp],
                                                                                                                                  imaging_temp$total_sa.combat[imaging_index_temp], 
                                                                                                                                  imaging_temp$mean_thickness.combat[imaging_index_temp],
                                                                                                                                  as.character(imaging_temp$EXAMDATE[imaging_index_temp]), 
                                                                                                                                  imaging_temp$age_at_baseline[imaging_index_temp] + interval(ymd(imaging_temp$baseline_date[imaging_index_temp]),  ymd(imaging_temp$EXAMDATE[imaging_index_temp])) %/% months(1)
                                                                                                                                  # imaging_temp$DX[imaging_index_temp],
                                                                                                                                  # imaging_temp$APOE4[imaging_index_temp]
  )
  
  
  
  fs4_paired_values <- rbind(fs4_paired_values, obs_pair)
}
fs4_paired_values <- fs4_paired_values[2:nrow(fs4_paired_values),]

fs4_paired_values$mean_age <- rowMeans(data.frame(as.numeric(fs4_paired_values$age_at_meth), as.numeric(fs4_paired_values$age_at_scan)))
fs4_paired_values$zPoAm45 <- as.numeric(fs4_paired_values$zPoAm45)
fs4_paired_values$zHorvath <- as.numeric(fs4_paired_values$zHorvath)
fs4_paired_values$zHannum <- as.numeric(fs4_paired_values$zHannum)
fs4_paired_values$zPheno <- as.numeric(fs4_paired_values$zPheno)
fs4_paired_values$zGrim <- as.numeric(fs4_paired_values$zGrim)
fs4_paired_values$mean_thickness <- as.numeric(fs4_paired_values$mean_thickness)
fs4_paired_values$total_sa <- as.numeric(fs4_paired_values$total_sa)
fs4_paired_values$mean_thickness.combat <- as.numeric(fs4_paired_values$mean_thickness.combat)
fs4_paired_values$total_sa.combat <- as.numeric(fs4_paired_values$total_sa.combat)

adni_sex <- data.frame(RID = adnimerge$RID, sex = adnimerge$PTGENDER)
adni_sex <- adni_sex[!duplicated(adni_sex),]
fs4_paired_values <- merge(adni_sex, fs4_paired_values, by = 'RID')
fs4_paired_values$sex <- as.factor(fs4_paired_values$sex)
fs4_paired_values$RID <- as.factor(fs4_paired_values$RID)



## combining fs4 and fs5

fs_paired_values <- rbind(fs5_paired_values, fs4_paired_values)
sum(fs_paired_values$QC != 'Pass')


for (i in unique(fs_paired_values$RID)){
  df_temp <- fs_paired_values[fs_paired_values$RID == i,]
  if (length(unique(df_temp$datescanned)) != nrow(df_temp)){
    print(i)
  }
}

fs.p.dat <- pdata.frame(data.frame(RID = as.character(fs_paired_values$RID),
                                   sa = scale(fs_paired_values$total_sa),
                                   mean_thickness = scale(fs_paired_values$mean_thickness),
                                   sa.combat = scale(fs_paired_values$total_sa.combat),
                                   mean_thickness.combat = scale(fs_paired_values$mean_thickness.combat),
                                   zPoAm45=scale(fs_paired_values$zPoAm45), 
                                   zHorvath=scale(fs_paired_values$zHorvath),
                                   zHannum=scale(fs_paired_values$zHannum),
                                   zPheno=scale(fs_paired_values$zPheno),
                                   zGrim=scale(fs_paired_values$zGrim),
                                   sex = fs_paired_values$sex,
                                   mean_age = scale(fs_paired_values$mean_age),
                                   mean_age_sq = scale(fs_paired_values$mean_age^2),
                                   datescanned = fs_paired_values$datescanned,
                                   datedrawn = fs_paired_values$datedrawn),
                        index = c("RID"), drop.index = F, row.names = T)


# CT

#dunedinpace
pm1 <- plm(mean_thickness.combat~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.p.dat, na.action = na.omit)
G <- length(unique(fs.p.dat$RID))
N <- length(fs.p.dat$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
ct.dp.plm.res <- coeftest(pm1, vcov = firm_c_vcov)

#horvath
pm1 <- plm(mean_thickness.combat~zHorvath+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.p.dat, na.action = na.omit)
G <- length(unique(fs_paired_values$RID))
N <- length(fs_paired_values$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
ct.hth.plm.res <- coeftest(pm1, vcov = firm_c_vcov)

#hannum
pm1 <- plm(mean_thickness.combat~zHannum+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.p.dat, na.action = na.omit)
G <- length(unique(fs_paired_values$RID))
N <- length(fs_paired_values$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
ct.hnm.plm.res <- coeftest(pm1, vcov = firm_c_vcov)

#phenoage
pm1 <- plm(mean_thickness.combat~zPheno+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.p.dat, na.action = na.omit)
G <- length(unique(fs_paired_values$RID))
N <- length(fs_paired_values$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
ct.phe.plm.res <- coeftest(pm1, vcov = firm_c_vcov)

#grimage
pm1 <- plm(mean_thickness.combat~zGrim+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.p.dat, na.action = na.omit)
G <- length(unique(fs_paired_values$RID))
N <- length(fs_paired_values$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
ct.grim.plm.res <- coeftest(pm1, vcov = firm_c_vcov)


# SA

#dunedinpace
pm1 <- plm(sa.combat~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.p.dat, na.action = na.omit)
G <- length(unique(fs.p.dat$RID))
N <- length(fs.p.dat$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
sa.dp.plm.res <- coeftest(pm1, vcov = firm_c_vcov)

#horvath
pm1 <- plm(sa.combat~zHorvath+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.p.dat, na.action = na.omit)
G <- length(unique(fs_paired_values$RID))
N <- length(fs_paired_values$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
sa.hth.plm.res <- coeftest(pm1, vcov = firm_c_vcov)

#hannum
pm1 <- plm(sa.combat~zHannum+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.p.dat, na.action = na.omit)
G <- length(unique(fs_paired_values$RID))
N <- length(fs_paired_values$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
sa.hnm.plm.res <- coeftest(pm1, vcov = firm_c_vcov)

#phenoage
pm1 <- plm(sa.combat~zPheno+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.p.dat, na.action = na.omit)
G <- length(unique(fs_paired_values$RID))
N <- length(fs_paired_values$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
sa.phe.plm.res <- coeftest(pm1, vcov = firm_c_vcov)

#grimage
pm1 <- plm(sa.combat~zGrim+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.p.dat, na.action = na.omit)
G <- length(unique(fs_paired_values$RID))
N <- length(fs_paired_values$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
sa.grim.plm.res <- coeftest(pm1, vcov = firm_c_vcov)




##### saving out files so I dont have to run all of this again

# save(tbv.p.dat, file='/Users/ew198/Documents/methylation/data/adni/annchen/tbv.p.dat.Rdata')
# save(hc.p.dat, file='/Users/ew198/Documents/methylation/data/adni/annchen/hc.p.dat.Rdata')
# save(wmh.p.dat, file='/Users/ew198/Documents/methylation/data/adni/annchen/wmh.p.dat.Rdata')
# save(fs.p.dat, file='/Users/ew198/Documents/methylation/data/adni/annchen/fs.p.dat.Rdata')

load('/Users/ew198/Documents/methylation/data/adni/annchen/tbv.p.dat.Rdata')
load('/Users/ew198/Documents/methylation/data/adni/annchen/hc.p.dat.Rdata')
load('/Users/ew198/Documents/methylation/data/adni/annchen/wmh.p.dat.Rdata')
load('/Users/ew198/Documents/methylation/data/adni/annchen/fs.p.dat.Rdata')
load('/Users/ew198/Documents/methylation/data/adni/annchen/fs.wmhypo.p.dat.Rdata')



##################################################
############ STRAFITY RESULTS BY APOE ############ 
##################################################

# get apoe carriership
adni_meth_apoe <- methMerge[!duplicated(methMerge$RID),c('RID', 'APOE4')]

adni_meth_apoe_carriers <- adni_meth_apoe$RID[adni_meth_apoe$APOE4 != 0]
adni_meth_apoe_noncarriers <- adni_meth_apoe$RID[adni_meth_apoe$APOE4 == 0]

# APOE4 carriers

# TBV
pm1 <- plm(TBV~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex+ICV, model = "pooling", data = tbv.p.dat[tbv.p.dat$RID %in% adni_meth_apoe_carriers,], na.action = na.omit)
G <- length(unique(tbv.p.dat[tbv.p.dat$RID %in% adni_meth_apoe_carriers,]$RID))
N <- length(tbv.p.dat[tbv.p.dat$RID %in% adni_meth_apoe_carriers,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
tbvrel.dp.plm.res.ac <- coeftest(pm1, vcov = firm_c_vcov)

# HC
pm1 <- plm(HC~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex+ICV, model = "pooling", data = hc.p.dat[hc.p.dat$RID %in% adni_meth_apoe_carriers,], na.action = na.omit)
G <- length(unique(hc.p.dat[hc.p.dat$RID %in% adni_meth_apoe_carriers,]$RID))
N <- length(hc.p.dat[hc.p.dat$RID %in% adni_meth_apoe_carriers,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
hcrel.dp.plm.res.ac <- coeftest(pm1, vcov = firm_c_vcov)

# WMH
pm1 <- plm(wmh~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = wmh.p.dat[wmh.p.dat$RID %in% adni_meth_apoe_carriers,], na.action = na.omit)
G <- length(unique(wmh.p.dat[wmh.p.dat$RID %in% adni_meth_apoe_carriers,]$RID))
N <- length(wmh.p.dat[wmh.p.dat$RID %in% adni_meth_apoe_carriers,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmh.dp.plm.res.ac <- coeftest(pm1, vcov = firm_c_vcov)

# CT
pm1 <- plm(mean_thickness.combat~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.p.dat[fs.p.dat$RID %in% adni_meth_apoe_carriers,], na.action = na.omit)
G <- length(unique(fs.p.dat[fs.p.dat$RID %in% adni_meth_apoe_carriers,]$RID))
N <- length(fs.p.dat[fs.p.dat$RID %in% adni_meth_apoe_carriers,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
ct.dp.plm.res.ac <- coeftest(pm1, vcov = firm_c_vcov)

# SA
pm1 <- plm(sa.combat~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.p.dat[fs.p.dat$RID %in% adni_meth_apoe_carriers,], na.action = na.omit)
G <- length(unique(fs.p.dat[fs.p.dat$RID %in% adni_meth_apoe_carriers,]$RID))
N <- length(fs.p.dat[fs.p.dat$RID %in% adni_meth_apoe_carriers,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
sa.dp.plm.res.ac <- coeftest(pm1, vcov = firm_c_vcov)

# WMHypo - added after next section
pm1 <- plm(wmhypo.combat.log~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.wmhypo.p.dat[fs.wmhypo.p.dat$RID %in% adni_meth_apoe_carriers,], na.action = na.omit)
G <- length(unique(fs.wmhypo.p.dat[fs.wmhypo.p.dat$RID %in% adni_meth_apoe_carriers,]$RID))
N <- length(fs.wmhypo.p.dat[fs.wmhypo.p.dat$RID %in% adni_meth_apoe_carriers,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmhypo.dp.plm.res.ac <- coeftest(pm1, vcov = firm_c_vcov)


# APOE4 noncarriers

# TBV
pm1 <- plm(TBV~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex+ICV, model = "pooling", data = tbv.p.dat[tbv.p.dat$RID %in% adni_meth_apoe_noncarriers,], na.action = na.omit)
G <- length(unique(tbv.p.dat[tbv.p.dat$RID %in% adni_meth_apoe_noncarriers,]$RID))
N <- length(tbv.p.dat[tbv.p.dat$RID %in% adni_meth_apoe_noncarriers,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
tbvrel.dp.plm.res.nc <- coeftest(pm1, vcov = firm_c_vcov)

# HC
pm1 <- plm(HC~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex+ICV, model = "pooling", data = hc.p.dat[hc.p.dat$RID %in% adni_meth_apoe_noncarriers,], na.action = na.omit)
G <- length(unique(hc.p.dat[hc.p.dat$RID %in% adni_meth_apoe_noncarriers,]$RID))
N <- length(hc.p.dat[hc.p.dat$RID %in% adni_meth_apoe_noncarriers,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
hcrel.dp.plm.res.nc <- coeftest(pm1, vcov = firm_c_vcov)

# WMH
pm1 <- plm(wmh~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = wmh.p.dat[wmh.p.dat$RID %in% adni_meth_apoe_noncarriers,], na.action = na.omit)
G <- length(unique(wmh.p.dat[wmh.p.dat$RID %in% adni_meth_apoe_noncarriers,]$RID))
N <- length(wmh.p.dat[wmh.p.dat$RID %in% adni_meth_apoe_noncarriers,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmh.dp.plm.res.nc <- coeftest(pm1, vcov = firm_c_vcov)

# CT
pm1 <- plm(mean_thickness.combat~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.p.dat[fs.p.dat$RID %in% adni_meth_apoe_noncarriers,], na.action = na.omit)
G <- length(unique(fs.p.dat[fs.p.dat$RID %in% adni_meth_apoe_noncarriers,]$RID))
N <- length(fs.p.dat[fs.p.dat$RID %in% adni_meth_apoe_noncarriers,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
ct.dp.plm.res.nc <- coeftest(pm1, vcov = firm_c_vcov)

# SA
pm1 <- plm(sa.combat~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.p.dat[fs.p.dat$RID %in% adni_meth_apoe_noncarriers,], na.action = na.omit)
G <- length(unique(fs.p.dat[fs.p.dat$RID %in% adni_meth_apoe_noncarriers,]$RID))
N <- length(fs.p.dat[fs.p.dat$RID %in% adni_meth_apoe_noncarriers,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
sa.dp.plm.res.nc <- coeftest(pm1, vcov = firm_c_vcov)

# WMHypo - added after next section 
pm1 <- plm(wmhypo.combat.log~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.wmhypo.p.dat[fs.wmhypo.p.dat$RID %in% adni_meth_apoe_noncarriers,], na.action = na.omit)
G <- length(unique(fs.wmhypo.p.dat[fs.wmhypo.p.dat$RID %in% adni_meth_apoe_noncarriers,]$RID))
N <- length(fs.wmhypo.p.dat[fs.wmhypo.p.dat$RID %in% adni_meth_apoe_noncarriers,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmhypo.dp.plm.res.nc <- coeftest(pm1, vcov = firm_c_vcov)


# apoe stratification forest plot

adni_apoe_strat_table <- data.frame(imaging = rep(c('TBV', 'HC', 'WMHyper', 'CT', 'SA', 'WMHypo'), 3),
                               genotype = c(rep('full sample', 6), rep('APOE4 carriers',6), rep('APOE4 noncarriers', 6)),
                               beta = c(tbvrel.dp.plm.res[2,1],  hcrel.dp.plm.res[2,1], wmh.dp.plm.res[2,1],  ct.dp.plm.res[2,1], sa.dp.plm.res[2,1],  wmhypo.dp.plm.res[2,1],
                                        tbvrel.dp.plm.res.ac[2,1],  hcrel.dp.plm.res.ac[2,1], wmh.dp.plm.res.ac[2,1],  ct.dp.plm.res.ac[2,1], sa.dp.plm.res.ac[2,1],  wmhypo.dp.plm.res.ac[2,1],
                                        tbvrel.dp.plm.res.nc[2,1],  hcrel.dp.plm.res.nc[2,1], wmh.dp.plm.res.nc[2,1],  ct.dp.plm.res.nc[2,1], sa.dp.plm.res.nc[2,1],  wmhypo.dp.plm.res.nc[2,1]),
                               stderror = c(tbvrel.dp.plm.res[2,2],  hcrel.dp.plm.res[2,2], wmh.dp.plm.res[2,2],  ct.dp.plm.res[2,2], sa.dp.plm.res[2,2],  wmhypo.dp.plm.res[2,2],
                                            tbvrel.dp.plm.res.ac[2,2],  hcrel.dp.plm.res.ac[2,2], wmh.dp.plm.res.ac[2,2],  ct.dp.plm.res.ac[2,2], sa.dp.plm.res.ac[2,2],  wmhypo.dp.plm.res.ac[2,2],
                                            tbvrel.dp.plm.res.nc[2,2],  hcrel.dp.plm.res.nc[2,2], wmh.dp.plm.res.nc[2,2],  ct.dp.plm.res.nc[2,2], sa.dp.plm.res.nc[2,2],  wmhypo.dp.plm.res.nc[2,2]))

# point range plot

imaging_level_order_strat <- rev(c('TBV', 'HC', 'WMHypo', 'CT', 'SA', 'WMHyper'))
genotype_levels <- c('full sample', 'APOE4 noncarriers', 'APOE4 carriers')
adni_apoe_colors <- c('#619CFF', 'grey70', 'yellow3')

adni_apoe_strat <- ggplot(adni_apoe_strat_table, 
                          aes(x=factor(imaging, imaging_level_order_strat), 
                              y=beta, color = factor(genotype, levels=genotype_levels), group=factor(genotype, levels = genotype_levels))) + 
  geom_pointrange(aes(ymin=beta-(1.96*stderror), ymax=beta+(1.96*stderror)),position=position_dodge(width=0.6)) +
  ylim(-.35, .35)+
  geom_hline(yintercept = 0)+
  labs(title = 'ADNI') +
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
  scale_color_manual(values=c(rep(adni_apoe_colors,3)))
#save(adni_apoe_strat, file = '/Users/ew198/Documents/methylation/results/figures/adni_apoe_strat.Rdata')

grid.arrange(dunedin_apoe_strat, adni_apoe_strat, nrow = 2)



##################################################
############ WMHyper WMHypo correlation ########## 
##################################################

# first run combat on WMHypo estimates because they are processed with different versions of freesurfer
# WMHypo is ST128SV according to the FS5 dict
# WMHypo is also ST128SV in the FS6 dict

methWMHypo_comb_temp <- rbind(methFS4[,c('RID', 'EXAMDATE', 'VERSION', 'RUNDATE', 'OVERALLQC', 'baseline_date', 'age_at_baseline', 'ST128SV', 'mean_thickness', 'total_sa')],
                         methFS5[,c('RID', 'EXAMDATE', 'VERSION', 'RUNDATE','OVERALLQC',  'baseline_date', 'age_at_baseline', 'ST128SV', 'mean_thickness', 'total_sa')],
                         methFS6[,c('RID', 'EXAMDATE', 'VERSION', 'RUNDATE','OVERALLQC',  'baseline_date', 'age_at_baseline', 'ST128SV', 'mean_thickness', 'total_sa')])

methWMHypo_comb_temp <- data.frame(methWMHypo_comb_temp, FSVERSION=c(rep(4, nrow(methFS4)), rep(5, nrow(methFS5)), rep(6, nrow(methFS6)) ) )
methWMHypo_comb_temp$age_at_scan <- methWMHypo_comb_temp$age_at_baseline + (interval(methWMHypo_comb_temp$baseline_date, methWMHypo_comb_temp$EXAMDATE) %/% months(1))
methWMHypo_comb_temp <- merge(methWMHypo_comb_temp, adni_sex, by = 'RID')

# run longCombat_adni_fs.R to get harmonized values for CT and SA
methWMHypo_comb_temp$scanage_scale <- scale(methWMHypo_comb_temp$age_at_scan)
methWMHypo_comb_temp_forharm <- methWMHypo_comb_temp
methWMHypo_comb_temp_forharm$age_years <- methWMHypo_comb_temp_forharm$age_at_scan/12
adni_WMHypo_combat <- longCombat(idvar='RID', 
                          timevar='age_years',
                          batchvar='FSVERSION', 
                          features=c('ST128SV', 'mean_thickness', 'total_sa'), 
                          formula='sex',
                          ranef='(1+age_years|RID)',
                          data=methWMHypo_comb_temp_forharm[,c('RID', 'ST128SV', 'mean_thickness', 'total_sa', 'FSVERSION', 'sex', 'age_at_scan', 'age_years')])

adni_WMHypo_harm <- adni_WMHypo_combat$data_combat
adni_WMHypo_harm$sex <- methWMHypo_comb_temp_forharm$sex
adni_WMHypo_harm$scanage_scale <- scale(adni_WMHypo_harm$age_years)

if (identical(methWMHypo_comb_temp_forharm$RID, adni_WMHypo_harm$RID) && sum(abs(methWMHypo_comb_temp_forharm$scanage_scale - adni_WMHypo_harm$scanage_scale)) < 0.0000000001 ) {
  print('IDENTICAL')
  WMHypo_comb_harmed <- data.frame(methWMHypo_comb_temp_forharm, adni_WMHypo_harm[,-6])
} else{
  print('NOT IDENTICAL')
}

WMHypo_comb_harmed$WMHypo.combat <- WMHypo_comb_harmed$ST128SV.combat
WMHypo_comb_harmed$WMHypo.combat_log <- log(WMHypo_comb_harmed$WMHypo.combat)

wmh_paired_values$age_at_scan_formerge <- as.character(round(wmh_paired_values$logWMH_AGE.img*12,1))
WMHypo_comb_harmed$age_at_scan_formerge <- as.character(round(WMHypo_comb_harmed$age_at_scan,1))

# save(wmh_paired_values_WMHypo, file = '/Users/ew198/Documents/methylation/data/adni/wmh_paired_values_WMHypo.Rdata')

# correlation between WMHyper and WMHypo
cor.test(wmh_paired_values_WMHypo$wmh, wmh_paired_values_WMHypo$WMHypo.combat_log)


# find WMHypo and methylation pairs - luckily already did this with fs.p.dat
fs_paired_values$mean_thickness_merge <- as.character(round(fs_paired_values$mean_thickness, 6))
WMHypo_comb_harmed$mean_thickness_merge <- as.character(round(WMHypo_comb_harmed$mean_thickness,6))

View(fs_paired_values[fs_paired_values$RID==123,])
View(WMHypo_comb_harmed[WMHypo_comb_harmed$RID==123,])

WMHypo_comb_harmed$RID <- factor(WMHypo_comb_harmed$RID)
fs_paired_values_WMHypo <- left_join(fs_paired_values, WMHypo_comb_harmed[,c('RID', 'EXAMDATE', 'baseline_date', 'age_at_baseline', 'age_at_scan', 'WMHypo.combat', 'WMHypo.combat_log', 'mean_thickness_merge')], 
                                  by = c('RID', 'mean_thickness_merge'))

fs.wmhypo.p.dat <- pdata.frame(data.frame(RID = as.character(fs_paired_values_WMHypo$RID),
                                   sa = scale(fs_paired_values_WMHypo$total_sa),
                                   wmhypo.combat = scale(fs_paired_values_WMHypo$WMHypo.combat),
                                   wmhypo.combat.log = scale(fs_paired_values_WMHypo$WMHypo.combat_log),
                                   zPoAm45=scale(fs_paired_values_WMHypo$zPoAm45), 
                                   zHorvath=scale(fs_paired_values_WMHypo$zHorvath),
                                   zHannum=scale(fs_paired_values_WMHypo$zHannum),
                                   zPheno=scale(fs_paired_values_WMHypo$zPheno),
                                   zGrim=scale(fs_paired_values_WMHypo$zGrim),
                                   sex = fs_paired_values_WMHypo$sex,
                                   mean_age = scale(fs_paired_values_WMHypo$mean_age),
                                   mean_age_sq = scale(fs_paired_values_WMHypo$mean_age^2),
                                   datescanned = fs_paired_values_WMHypo$datescanned,
                                   datedrawn = fs_paired_values_WMHypo$datedrawn),
                        index = c("RID"), drop.index = F, row.names = T)

# save out!!
#save(fs.wmhypo.p.dat, file='/Users/ew198/Documents/methylation/data/adni/annchen/fs.wmhypo.p.dat.Rdata')

# run analyses with WMHypo

# WMHypo

# dunedinPACE
pm1 <- plm(wmhypo.combat.log~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.wmhypo.p.dat, na.action = na.omit)
G <- length(unique(fs_paired_values_WMHypo$RID))
N <- length(fs_paired_values_WMHypo$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmhypo.dp.plm.res <- coeftest(pm1, vcov = firm_c_vcov)

# horvath
pm1 <- plm(wmhypo.combat.log~zHorvath+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.wmhypo.p.dat, na.action = na.omit)
G <- length(unique(fs_paired_values_WMHypo$RID))
N <- length(fs_paired_values_WMHypo$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmhypo.hth.plm.res <- coeftest(pm1, vcov = firm_c_vcov)

# hannum
pm1 <- plm(wmhypo.combat.log~zHannum+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.wmhypo.p.dat, na.action = na.omit)
G <- length(unique(fs_paired_values_WMHypo$RID))
N <- length(fs_paired_values_WMHypo$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmhypo.hnm.plm.res <- coeftest(pm1, vcov = firm_c_vcov)

# phenoage
pm1 <- plm(wmhypo.combat.log~zPheno+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.wmhypo.p.dat, na.action = na.omit)
G <- length(unique(fs_paired_values_WMHypo$RID))
N <- length(fs_paired_values_WMHypo$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmhypo.phe.plm.res <- coeftest(pm1, vcov = firm_c_vcov)

# grimage
pm1 <- plm(wmhypo.combat.log~zGrim+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.wmhypo.p.dat, na.action = na.omit)
G <- length(unique(fs_paired_values_WMHypo$RID))
N <- length(fs_paired_values_WMHypo$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmhypo.grim.plm.res <- coeftest(pm1, vcov = firm_c_vcov)


# just WMHyper sample for effect size comparison

# dunedinPACE
pm1 <- plm(WMHypo.combat_log~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = pdata.frame(wmh_paired_values_WMHypo), na.action = na.omit)
G <- length(unique(pdata.frame(wmh_paired_values_WMHypo)$RID))
N <- length(pdata.frame(wmh_paired_values_WMHypo)$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmhypo.dp.plm.res.ss <- coeftest(pm1, vcov = firm_c_vcov)

# horvath
pm1 <- plm(WMHypo.combat_log~zHorvath+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = pdata.frame(wmh_paired_values_WMHypo), na.action = na.omit)
G <- length(unique(pdata.frame(wmh_paired_values_WMHypo)$RID))
N <- length(pdata.frame(wmh_paired_values_WMHypo)$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmhypo.hth.plm.res.ss <- coeftest(pm1, vcov = firm_c_vcov)

# hannum
pm1 <- plm(WMHypo.combat_log~zHannum+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = pdata.frame(wmh_paired_values_WMHypo), na.action = na.omit)
G <- length(unique(pdata.frame(wmh_paired_values_WMHypo)$RID))
N <- length(pdata.frame(wmh_paired_values_WMHypo)$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmhypo.hnm.plm.res.ss <- coeftest(pm1, vcov = firm_c_vcov)

# phenoage
pm1 <- plm(WMHypo.combat_log~zPheno+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = pdata.frame(wmh_paired_values_WMHypo), na.action = na.omit)
G <- length(unique(pdata.frame(wmh_paired_values_WMHypo)$RID))
N <- length(pdata.frame(wmh_paired_values_WMHypo)$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmhypo.phe.plm.res.ss <- coeftest(pm1, vcov = firm_c_vcov)

# grimage
pm1 <- plm(WMHypo.combat_log~zGrim+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = pdata.frame(wmh_paired_values_WMHypo), na.action = na.omit)
G <- length(unique(pdata.frame(wmh_paired_values_WMHypo)$RID))
N <- length(pdata.frame(wmh_paired_values_WMHypo)$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmhypo.grim.plm.res.ss <- coeftest(pm1, vcov = firm_c_vcov)







### sensitivity tests

# APOE4 covariate 

fs.wmhypo.p.dat.apoe <- merge(fs.wmhypo.p.dat, adni_meth_apoe, by = 'RID')

# dunedinPACE
pm1 <- plm(wmhypo.combat.log~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex + APOE4, model = "pooling", data = fs.wmhypo.p.dat.apoe, na.action = na.omit)
G <- length(unique(fs_paired_values_WMHypo$RID))
N <- length(fs_paired_values_WMHypo$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmhypo.dp.plm.res.apoe <- coeftest(pm1, vcov = firm_c_vcov)

# horvath
pm1 <- plm(wmhypo.combat.log~zHorvath+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex + APOE4, model = "pooling", data = fs.wmhypo.p.dat.apoe, na.action = na.omit)
G <- length(unique(fs_paired_values_WMHypo$RID))
N <- length(fs_paired_values_WMHypo$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmhypo.hth.plm.res.apoe <- coeftest(pm1, vcov = firm_c_vcov)

# hannum
pm1 <- plm(wmhypo.combat.log~zHannum+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex + APOE4, model = "pooling", data = fs.wmhypo.p.dat.apoe, na.action = na.omit)
G <- length(unique(fs_paired_values_WMHypo$RID))
N <- length(fs_paired_values_WMHypo$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmhypo.hnm.plm.res.apoe <- coeftest(pm1, vcov = firm_c_vcov)

# phenoage
pm1 <- plm(wmhypo.combat.log~zPheno+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex + APOE4, model = "pooling", data = fs.wmhypo.p.dat.apoe, na.action = na.omit)
G <- length(unique(fs_paired_values_WMHypo$RID))
N <- length(fs_paired_values_WMHypo$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmhypo.phe.plm.res.apoe <- coeftest(pm1, vcov = firm_c_vcov)

# grimage
pm1 <- plm(wmhypo.combat.log~zGrim+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex + APOE4, model = "pooling", data = fs.wmhypo.p.dat.apoe, na.action = na.omit)
G <- length(unique(fs_paired_values_WMHypo$RID))
N <- length(fs_paired_values_WMHypo$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmhypo.grim.plm.res.apoe <- coeftest(pm1, vcov = firm_c_vcov)

# WBC abundance

fs.wmhypo.p.dat.wca <- pdata.frame(data.frame(RID = as.character(fs_paired_values_WMHypo$RID),
                                          sa = scale(fs_paired_values_WMHypo$total_sa),
                                          wmhypo.combat = scale(fs_paired_values_WMHypo$WMHypo.combat),
                                          wmhypo.combat.log = scale(fs_paired_values_WMHypo$WMHypo.combat_log),
                                          zPoAm45=scale(fs_paired_values_WMHypo$zPoAm45), 
                                          zHorvath=scale(fs_paired_values_WMHypo$zHorvath),
                                          zHannum=scale(fs_paired_values_WMHypo$zHannum),
                                          zPheno=scale(fs_paired_values_WMHypo$zPheno),
                                          zGrim=scale(fs_paired_values_WMHypo$zGrim),
                                          sex = fs_paired_values_WMHypo$sex,
                                          mean_age = scale(fs_paired_values_WMHypo$mean_age),
                                          mean_age_sq = scale(fs_paired_values_WMHypo$mean_age^2),
                                          CD4T = as.numeric(fs_paired_values_WMHypo$CD4T),
                                          NK = as.numeric(fs_paired_values_WMHypo$NK),
                                          Mono = as.numeric(fs_paired_values_WMHypo$Mono),
                                          Gran = as.numeric(fs_paired_values_WMHypo$Gran),
                                          PlasmaBlast = as.numeric(fs_paired_values_WMHypo$PlasmaBlast),
                                          CD8pCD28nCD45RAn = as.numeric(fs_paired_values_WMHypo$CD8pCD28nCD45RAn),
                                          CD8.naive = as.numeric(fs_paired_values_WMHypo$CD8.naive)),
                               index = c("RID"), drop.index = F, row.names = T)



#dunedinpace
pm1 <- plm(wmhypo.combat.log~zPoAm45+sex+mean_age+ sex*mean_age+mean_age_sq+mean_age_sq*sex+ CD4T + NK + Mono + Gran + PlasmaBlast + CD8pCD28nCD45RAn + CD8.naive, model = "pooling", data = fs.wmhypo.p.dat.wca, na.action = na.omit)
G <- length(unique(fs_paired_values_WMHypo$RID))
N <- length(fs_paired_values_WMHypo$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmhypo.dp.plm.res.wca <- coeftest(pm1, vcov = firm_c_vcov)

#horvath
pm1 <- plm(wmhypo.combat.log~zHorvath+sex+mean_age+ sex*mean_age+mean_age_sq+mean_age_sq*sex+ CD4T + NK + Mono + Gran + PlasmaBlast + CD8pCD28nCD45RAn + CD8.naive, model = "pooling", data = fs.wmhypo.p.dat.wca, na.action = na.omit)
G <- length(unique(fs_paired_values_WMHypo$RID))
N <- length(fs_paired_values_WMHypo$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmhypo.hth.plm.res.wca <- coeftest(pm1, vcov = firm_c_vcov)

#hannum
pm1 <- plm(wmhypo.combat.log~zHannum+sex+mean_age+ sex*mean_age+mean_age_sq+mean_age_sq*sex+ CD4T + NK + Mono + Gran + PlasmaBlast + CD8pCD28nCD45RAn + CD8.naive, model = "pooling", data = fs.wmhypo.p.dat.wca, na.action = na.omit)
G <- length(unique(fs_paired_values_WMHypo$RID))
N <- length(fs_paired_values_WMHypo$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmhypo.hnm.plm.res.wca <- coeftest(pm1, vcov = firm_c_vcov)

#phenoage
pm1 <- plm(wmhypo.combat.log~zPheno+sex+mean_age+ sex*mean_age+mean_age_sq+mean_age_sq*sex+ CD4T + NK + Mono + Gran + PlasmaBlast + CD8pCD28nCD45RAn + CD8.naive, model = "pooling", data = fs.wmhypo.p.dat.wca, na.action = na.omit)
G <- length(unique(fs_paired_values_WMHypo$RID))
N <- length(fs_paired_values_WMHypo$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmhypo.phe.plm.res.wca <- coeftest(pm1, vcov = firm_c_vcov)

#grimage
pm1 <- plm(wmhypo.combat.log~zGrim+sex+mean_age+ sex*mean_age+mean_age_sq+mean_age_sq*sex+ CD4T + NK + Mono + Gran + PlasmaBlast + CD8pCD28nCD45RAn + CD8.naive, model = "pooling", data = fs.wmhypo.p.dat.wca, na.action = na.omit)
G <- length(unique(fs_paired_values_WMHypo$RID))
N <- length(fs_paired_values_WMHypo$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmhypo.grim.plm.res.wca <- coeftest(pm1, vcov = firm_c_vcov)



##### visualization ######


# WMHyper vs WMHypo correlation

cor.test(wmh_paired_values_WMHypo$wmh, wmh_paired_values_WMHypo$WMHypo.combat_log)


ggplot(data=data.frame(WMHyper = wmh_paired_values_WMHypo$wmh, 
                       WMHypo = wmh_paired_values_WMHypo$WMHypo.combat_log), 
       aes(x=WMHyper, y=WMHypo)) +
  geom_point() +
  labs(title = 'ADNI')+
 # ylim(-3, 3.5)+
 # xlim(-3,3) +
 # geom_abline(intercept=0,slope=1, size = 2) +
  geom_smooth(method = 'lm', color = 'red') +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(vjust = 0.5, size = 10),
        axis.title.y = element_text(size = 15, vjust = 0),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 15))

ggplot(data=data.frame(WMHyper = wmh_paired_values_WMHypo$wmh, 
                       WMHypo = wmh_paired_values_WMHypo$WMHypo.combat_log), 
       aes(x=scale(WMHyper), y=scale(WMHypo))) +
  geom_point() +
  labs(title = 'ADNI')+
  ylim(-3, 3.5)+
  xlim(-3,3.5) +
  geom_abline(intercept=0,slope=1, size = 2) +
  geom_smooth(method = 'lm', color = 'red') +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(vjust = 0.5, size = 10),
        axis.title.y = element_text(size = 15, vjust = 0),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 15))

adni_wmhypo_dist <- ggplot(data = wmh_paired_values_WMHypo, aes(x = WMHypo.combat)) + 
  geom_histogram(fill='#619CFF') +
  labs(title = "ADNI",
       x = "WMHypo volume (mm3)",
       y = 'Frequency') +
  xlim(-2000, 53000)+
  ylim(0, 400)+
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20, vjust = 1),
        axis.title.x = element_text(size = 20),
        title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.position = 'none')

adni_wmhyper_dist <- ggplot(data = wmh_paired_values_WMHypo, aes(x = exp(wmh))) + 
  geom_histogram(fill='#619CFF') +
  labs(title = "ADNI",
       x = "WMHyper volume (mm3)",
       y = 'Frequency') +
  xlim(-5, 85)+
  ylim(0, 400)+
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20, vjust = 1),
        axis.title.x = element_text(size = 20),
        title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.position = 'none')

grid.arrange(adni_wmhypo_dist, dunedin_wmhypo_dist, nrow =2)

grid.arrange(adni_wmhypo_dist, dunedin_wmhypo_dist, nrow =2)

# WMHypo forest plot

wmhypo_wmhyper_table <- data.frame(imaging = rep(c('WMHyper', 'WMHypo'), 5),
                                   clock = c(rep('dunedinpace', 2), rep('horvath', 2), rep('hannum',2), rep('phenoage',2), rep('grimage',2)),
                                   beta = c(wmh.dp.plm.res[2,1],  wmhypo.dp.plm.res[2,1], 
                                            wmh.hth.plm.res[2,1],  wmhypo.hth.plm.res[2,1],
                                            wmh.hnm.plm.res[2,1],  wmhypo.hnm.plm.res[2,1],
                                            wmh.phe.plm.res[2,1],  wmhypo.phe.plm.res[2,1],
                                            wmh.grim.plm.res[2,1],  wmhypo.grim.plm.res[2,1]),
                                   stderror = c(wmh.dp.plm.res[2,2],  wmhypo.dp.plm.res[2,2], 
                                                wmh.hth.plm.res[2,2],  wmhypo.hth.plm.res[2,2],
                                                wmh.hnm.plm.res[2,2],  wmhypo.hnm.plm.res[2,2],
                                                wmh.phe.plm.res[2,2],  wmhypo.phe.plm.res[2,2],
                                                wmh.grim.plm.res[2,2],  wmhypo.grim.plm.res[2,2]))

# point range plot

level_order <- rev(c('horvath', 'hannum', 'phenoage', 'grimage', 'dunedinpace'))

unequal_ss <- ggplot(wmhypo_wmhyper_table, aes(x=clock, y=beta, color = imaging, group=imaging)) + 
  geom_pointrange(aes(ymin=beta-(1.96*stderror), ymax=beta+(1.96*stderror)),position=position_dodge(width=0.6)) +
  geom_hline(yintercept = 0)+
  labs(title = 'WMHyper vs WMHypo') +
  ylim(-.075,.2)+
  coord_flip() +
  theme_classic()+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 0),
        axis.title.x = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))


# equal sample size

wmhypo_wmhyper_table_ss <- data.frame(imaging = rep(c('WMHyper', 'WMHypo'), 5),
                                   clock = c(rep('dunedinpace', 2), rep('horvath', 2), rep('hannum',2), rep('phenoage',2), rep('grimage',2)),
                                   beta = c(wmh.dp.plm.res[2,1],  wmhypo.dp.plm.res.ss[2,1], 
                                            wmh.hth.plm.res[2,1],  wmhypo.hth.plm.res.ss[2,1],
                                            wmh.hnm.plm.res[2,1],  wmhypo.hnm.plm.res.ss[2,1],
                                            wmh.phe.plm.res[2,1],  wmhypo.phe.plm.res.ss[2,1],
                                            wmh.grim.plm.res[2,1],  wmhypo.grim.plm.res.ss[2,1]),
                                   stderror = c(wmh.dp.plm.res[2,2],  wmhypo.dp.plm.res.ss[2,2], 
                                                wmh.hth.plm.res[2,2],  wmhypo.hth.plm.res.ss[2,2],
                                                wmh.hnm.plm.res[2,2],  wmhypo.hnm.plm.res.ss[2,2],
                                                wmh.phe.plm.res[2,2],  wmhypo.phe.plm.res.ss[2,2],
                                                wmh.grim.plm.res[2,2],  wmhypo.grim.plm.res.ss[2,2]))

# point range plot

level_order <- rev(c('horvath', 'hannum', 'phenoage', 'grimage', 'dunedinpace'))

equal_ss <- ggplot(wmhypo_wmhyper_table_ss, aes(x=clock, y=beta, color = imaging, group=imaging)) + 
  geom_pointrange(aes(ymin=beta-(1.96*stderror), ymax=beta+(1.96*stderror)),position=position_dodge(width=0.6)) +
  geom_hline(yintercept = 0)+
  labs(title = 'WMHyper vs WMHypo - equiv samples') +
  ylim(-.075,.2)+
  coord_flip() +
  theme_classic()+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 0),
        axis.title.x = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))

grid.arrange(unequal_ss, equal_ss, nrow =2)


##################################################
############ STRAFITY RESULTS BY DX ##############
##################################################


# this loads the longitudinal diagnostic groups - could look at this another time
# load("/Users/ew198/Documents/methylation/all_dx.R")

### pairing with cross-sectional diagnostic data


# imaging_observations <- list(tbv = adni_tbv_datescanned,
#                              hc = adni_hc_datescanned,
#                              wmh = adni_wmh_datescanned,
#                              dx = adni_dx_datediagnosed,
#                              fs5 = adni_fs5_datescanned, 
#                              fs6 = adni_fs6_datescanned,
#                              fs4 = adni_fs4_datescanned)
# 
# dx_time_differences_mats <- list()
# for (o in 1:length(imaging_observations)){
#   temp_merged <- merge(adni_dx_datediagnosed, imaging_observations[[o]], by = 'RID')
#   dx_time_differences_mats[[o]] <- array(rep(NA, 17*17*649), dim = c(17, 17, 649))
#   for (d in 2:18){
#     for (i in 45:58)
#       if(typeof(temp_merged[,i]) == 'character'){
#         dx_time_differences_mats[[o]][(d-1), (i-44),] <- interval(ymd(temp_merged[,d]), ymd(temp_merged[,i])) %/% days(1)
#       }
#   }
# }
# 
# 
# # min imaging/dx timepoint pairs
# 
# imaging_dx_obs_pairs <- list()
# imaging_dx_obs_pairs_pos <- list()
# for (o in 1:length(dx_time_differences_mats)){
#   pairs <- data.frame(RID = 0, dx = 0, imaging = 0, time_diff = 0)
#   pos_pairs <- data.frame(RID = 0, dx = 0, imaging = 0, time_diff = 0)
#   x <- 1
#   y <- 1
#   for (p in 1:649){
#     person_mat <- dx_time_differences_mats[[o]][,,p]
#     rid_temp <- temp_merged$RID[p]
#     for (i in 1:17){
#       for (d in 1:17){
#         if (abs(person_mat[d,i]) == min(abs(person_mat[,i]), na.rm = TRUE) && !is.na(person_mat[d,i])){
#           new_pair <- c(rid_temp, d, i, abs(person_mat[d,i]))
#           pairs[x,] <- t(new_pair)
#           x <- x+1
#         }
#         if (person_mat[d,i] == min(person_mat[,i][person_mat[,i] >= 0], na.rm = TRUE) && !is.na(person_mat[d,i]) && person_mat[d,i] >= 0){
#           new_pair <- c(rid_temp, d, i, abs(person_mat[d,i]))
#           pos_pairs[y,] <- t(new_pair)
#           y <- y+1
#         }
#       }
#     }
#   }
#   imaging_dx_obs_pairs[[o]] <- pairs
#   imaging_dx_obs_pairs_pos[[o]] <- pos_pairs
# }
# 
# 
# 
# # filter by imaging observations included because of their temporal proximity to DNA methylation observations
# 
# 
# imaging_dx_obs_pairs_filtered <- list()
# for (o in 1:length(imaging_dx_obs_pairs)){
#   imaging_dx_obs_pairs_filtered[[o]] <- left_join(imaging_meth_obs_pairs[[o]], imaging_dx_obs_pairs[[o]], by = c('RID', 'imaging'))
# }
# 
# 
# 
# # pulling DX info for each variable
# 
# # TBV 
# 
# obs_temp <- imaging_dx_obs_pairs_filtered[[1]]
# tbv_paired_dx <- data.frame(RID = 0, datescanned = 0, TBV = 0, datediagnosed = 0, DX = 0, datedrawn = 0)
# for (i in 1:nrow(imaging_dx_obs_pairs_filtered[[1]])){
#   obs_pair <- data.frame(RID = 0, datescanned = 0, TBV = 0, datediagnosed = 0,  DX = 0, datedrawn = 0)
#   obs_pair_index_temp <- obs_temp[i,]
#   dx_temp <- methMerge[methMerge$RID == obs_pair_index_temp$RID,]
#   dx_temp <- dx_temp[!is.na(dx_temp$DX),]
#   dx_temp <- dx_temp[order(ymd(dx_temp$EXAMDATE)),]
#   dx_index_temp <- obs_pair_index_temp$dx
#   obs_pair[,c('RID', 'datediagnosed', 'DX')] <- c(dx_temp$RID[dx_index_temp],
#                                                   as.character(dx_temp$EXAMDATE[dx_index_temp]),
#                                                   dx_temp$DX[dx_index_temp])
#   
#   imaging_temp <- methMerge[methMerge$RID == obs_pair_index_temp$RID,]
#   imaging_temp <- imaging_temp[!is.na(imaging_temp$WholeBrain),]
#   imaging_temp <- imaging_temp[order(ymd(imaging_temp$EXAMDATE)),]
#   imaging_index_temp <- obs_pair_index_temp$imaging
#   obs_pair[,c('datescanned', 'TBV')] <- c(as.character(imaging_temp$EXAMDATE[imaging_index_temp]),
#                                           imaging_temp$WholeBrain[imaging_index_temp])
#   meth_temp <- adni_methclocks[adni_methclocks$RID == obs_pair_index_temp$RID,]
#   meth_temp <- meth_temp[order(mdy(meth_temp$DateDrawn)),]
#   meth_index_temp <- obs_pair_index_temp$meth
#   obs_pair[,c('datedrawn')] <- c(meth_temp$DateDrawn[meth_index_temp])
#   tbv_paired_dx <- rbind(tbv_paired_dx, obs_pair)
# }
# tbv_paired_dx <- tbv_paired_dx[2:nrow(tbv_paired_dx),]
# tbv_paired_dx$RID <- as.factor(tbv_paired_dx$RID)
# 
# 
# # HC 
# 
# obs_temp <- imaging_dx_obs_pairs_filtered[[2]]
# hc_paired_dx <- data.frame(RID = 0, datescanned = 0, HC = 0, datediagnosed = 0, DX = 0, datedrawn = 0)
# for (i in 1:nrow(imaging_dx_obs_pairs_filtered[[2]])){
#   obs_pair <- data.frame(RID = 0, datescanned = 0, HC = 0, datediagnosed = 0,  DX = 0, datedrawn = 0)
#   obs_pair_index_temp <- obs_temp[i,]
#   dx_temp <- methMerge[methMerge$RID == obs_pair_index_temp$RID,]
#   dx_temp <- dx_temp[!is.na(dx_temp$DX),]
#   dx_temp <- dx_temp[order(ymd(dx_temp$EXAMDATE)),]
#   dx_index_temp <- obs_pair_index_temp$dx
#   obs_pair[,c('RID', 'datediagnosed', 'DX')] <- c(dx_temp$RID[dx_index_temp],
#                                                   as.character(dx_temp$EXAMDATE[dx_index_temp]),
#                                                   dx_temp$DX[dx_index_temp])
#   imaging_temp <- methMerge[methMerge$RID == obs_pair_index_temp$RID,]
#   imaging_temp <- imaging_temp[!is.na(imaging_temp$Hippocampus),]
#   imaging_temp <- imaging_temp[order(ymd(imaging_temp$EXAMDATE)),]
#   imaging_index_temp <- obs_pair_index_temp$imaging
#   obs_pair[,c('datescanned', 'HC')] <- c(as.character(imaging_temp$EXAMDATE[imaging_index_temp]),
#                                          imaging_temp$Hippocampus[imaging_index_temp])
#   meth_temp <- adni_methclocks[adni_methclocks$RID == obs_pair_index_temp$RID,]
#   meth_temp <- meth_temp[order(mdy(meth_temp$DateDrawn)),]
#   meth_index_temp <- obs_pair_index_temp$meth
#   obs_pair[,c('datedrawn')] <- c(meth_temp$DateDrawn[meth_index_temp])
#   hc_paired_dx <- rbind(hc_paired_dx, obs_pair)
# }
# hc_paired_dx <- hc_paired_dx[2:nrow(hc_paired_dx),]
# hc_paired_dx$RID <- as.factor(hc_paired_dx$RID)
# 
# 
# 
# # WMH
# 
# obs_temp <- imaging_dx_obs_pairs_filtered[[3]]
# wmh_paired_dx <- data.frame(RID = 0, datescanned = 0, wmh = 0, datediagnosed = 0, DX = 0, datedrawn = 0)
# for (i in 1:nrow(imaging_dx_obs_pairs_filtered[[3]])){
#   obs_pair <- data.frame(RID = 0, datescanned = 0, wmh = 0, datediagnosed = 0,  DX = 0, datedrawn = 0)
#   obs_pair_index_temp <- obs_temp[i,]
#   dx_temp <- methMerge[methMerge$RID == obs_pair_index_temp$RID,]
#   dx_temp <- dx_temp[!is.na(dx_temp$DX),]
#   dx_temp <- dx_temp[order(ymd(dx_temp$EXAMDATE)),]
#   dx_index_temp <- obs_pair_index_temp$dx
#   obs_pair[,c('RID', 'datediagnosed', 'DX')] <- c(dx_temp$RID[dx_index_temp],
#                                                   as.character(dx_temp$EXAMDATE[dx_index_temp]),
#                                                   dx_temp$DX[dx_index_temp])
#   imaging_temp <- wmh_merge[wmh_merge$RID == obs_pair_index_temp$RID,]
#   imaging_temp <- imaging_temp[!is.na(imaging_temp$wmhVol_log),]
#   imaging_temp <- imaging_temp[order(ymd(imaging_temp$EXAMDATE)),]
#   imaging_index_temp <- obs_pair_index_temp$imaging
#   obs_pair[,c('datescanned', 'wmh')] <- c(as.character(imaging_temp$EXAMDATE[imaging_index_temp]),
#                                           imaging_temp$wmhVol_log[imaging_index_temp])
#   meth_temp <- adni_methclocks[adni_methclocks$RID == obs_pair_index_temp$RID,]
#   meth_temp <- meth_temp[order(mdy(meth_temp$DateDrawn)),]
#   meth_index_temp <- obs_pair_index_temp$meth
#   obs_pair[,c('datedrawn')] <- c(meth_temp$DateDrawn[meth_index_temp])
#   wmh_paired_dx <- rbind(wmh_paired_dx, obs_pair)
# }
# wmh_paired_dx <- wmh_paired_dx[2:nrow(wmh_paired_dx),]
# wmh_paired_dx$RID <- as.factor(wmh_paired_dx$RID)



# FS5

obs_temp <- imaging_dx_obs_pairs_filtered[[5]]
fs5_paired_dx <- data.frame(RID = 0, datescanned = 0, mean_thickness = 0, total_sa = 0, datediagnosed = 0, DX = 0, datedrawn = 0)
for (i in 1:nrow(imaging_dx_obs_pairs_filtered[[5]])){
  obs_pair <- data.frame(RID = 0, datescanned = 0, mean_thickness = 0, total_sa = 0, datediagnosed = 0,  DX = 0, datedrawn = 0)
  obs_pair_index_temp <- obs_temp[i,]
  dx_temp <- methMerge[methMerge$RID == obs_pair_index_temp$RID,]
  dx_temp <- dx_temp[!is.na(dx_temp$DX),]
  dx_temp <- dx_temp[order(ymd(dx_temp$EXAMDATE)),]
  dx_index_temp <- obs_pair_index_temp$dx
  obs_pair[,c('RID', 'datediagnosed', 'DX')] <- c(dx_temp$RID[dx_index_temp],
                                                  as.character(dx_temp$EXAMDATE[dx_index_temp]),
                                                  dx_temp$DX[dx_index_temp])
  imaging_temp <- methFS5[methFS5$RID == obs_pair_index_temp$RID,]
  imaging_temp <- imaging_temp[!is.na(imaging_temp$total_sa),]
  imaging_temp <- imaging_temp[order(ymd(imaging_temp$EXAMDATE)),]
  imaging_index_temp <- obs_pair_index_temp$imaging
  obs_pair[,c('datescanned', 'total_sa', 'mean_thickness')] <- c(as.character(imaging_temp$EXAMDATE[imaging_index_temp]),
                                                                 imaging_temp$total_sa[imaging_index_temp], 
                                                                 imaging_temp$mean_thickness[imaging_index_temp])
  meth_temp <- adni_methclocks[adni_methclocks$RID == obs_pair_index_temp$RID,]
  meth_temp <- meth_temp[order(mdy(meth_temp$DateDrawn)),]
  meth_index_temp <- obs_pair_index_temp$meth
  obs_pair[,c('datedrawn')] <- c(meth_temp$DateDrawn[meth_index_temp])
  fs5_paired_dx <- rbind(fs5_paired_dx, obs_pair)
}
fs5_paired_dx <- fs5_paired_dx[2:nrow(fs5_paired_dx),]
fs5_paired_dx$RID <- as.factor(fs5_paired_dx$RID)


# FS4

obs_temp <- imaging_dx_obs_pairs_filtered[[7]]
fs4_paired_dx <- data.frame(RID = 0, datescanned = 0, mean_thickness = 0, total_sa = 0, datediagnosed = 0, DX = 0, datedrawn = 0)
for (i in 1:nrow(imaging_dx_obs_pairs_filtered[[7]])){
  obs_pair <- data.frame(RID = 0, datescanned = 0, mean_thickness = 0, total_sa = 0, datediagnosed = 0,  DX = 0, datedrawn = 0)
  obs_pair_index_temp <- obs_temp[i,]
  dx_temp <- methMerge[methMerge$RID == obs_pair_index_temp$RID,]
  dx_temp <- dx_temp[!is.na(dx_temp$DX),]
  dx_temp <- dx_temp[order(ymd(dx_temp$EXAMDATE)),]
  dx_index_temp <- obs_pair_index_temp$dx
  
  obs_pair[,c('RID', 'datediagnosed', 'DX')] <- c(dx_temp$RID[dx_index_temp],
                                                  as.character(dx_temp$EXAMDATE[dx_index_temp]),
                                                  dx_temp$DX[dx_index_temp])
  imaging_temp <- methFS4[methFS4$RID == obs_pair_index_temp$RID,]
  imaging_temp <- imaging_temp[!is.na(imaging_temp$total_sa),]
  imaging_temp <- imaging_temp[order(ymd(imaging_temp$EXAMDATE)),]
  imaging_index_temp <- obs_pair_index_temp$imaging
  obs_pair[,c('datescanned', 'total_sa', 'mean_thickness')] <- c(as.character(imaging_temp$EXAMDATE[imaging_index_temp]),
                                                                 imaging_temp$total_sa[imaging_index_temp], 
                                                                 imaging_temp$mean_thickness[imaging_index_temp])
  meth_temp <- adni_methclocks[adni_methclocks$RID == obs_pair_index_temp$RID,]
  meth_temp <- meth_temp[order(mdy(meth_temp$DateDrawn)),]
  meth_index_temp <- obs_pair_index_temp$meth
  obs_pair[,c('datedrawn')] <- c(meth_temp$DateDrawn[meth_index_temp])
  fs4_paired_dx <- rbind(fs4_paired_dx, obs_pair)
}
fs4_paired_dx <- fs4_paired_dx[2:nrow(fs4_paired_dx),]
fs4_paired_dx$RID <- as.factor(fs4_paired_dx$RID)
#combine fs4 and fs5
fs_paired_dx <- rbind(fs5_paired_dx, fs4_paired_dx)


# merge in DX data for WMH and FreeSurfer

wmh_paired_values <- wmh_paired_values[, !(names(wmh_paired_values) %in% "DX")]
wmh_paired_values_dx <- merge(wmh_paired_values, wmh_paired_dx[,c('RID', 'datescanned', 'datediagnosed', 'DX', 'wmh')], by = c('RID', 'datescanned', 'wmh'))
fs_paired_values <- fs_paired_values[, !(names(fs_paired_values) %in% "DX")]
fs_paired_values_dx <- merge(fs_paired_values, fs_paired_dx[,c('RID', 'datescanned', 'datediagnosed', 'DX', 'mean_thickness')], by = c('RID', 'datescanned', 'mean_thickness'))
fs_paired_values_WMHypo <- fs_paired_values_WMHypo[, !(names(fs_paired_values_WMHypo) %in% "DX")]
fs_paired_values_WMHypo_dx <- merge(fs_paired_values_WMHypo, fs_paired_dx[,c('RID', 'datescanned', 'datediagnosed', 'DX')], by = c('RID', 'datescanned'))


# DEMENTIA
tbv_dem <- tbv.p.dat$RID[tbv.p.dat$DX=='Dementia']
hc_dem <- hc.p.dat$RID[hc.p.dat$DX=='Dementia']
wmh_dem <- wmh.p.dat$RID[wmh_paired_values_dx$DX=='Dementia']
fs_dem <- fs.p.dat$RID[fs_paired_values_dx$DX=='Dementia']

# TBV
pm1 <- plm(TBV~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex+ICV, model = "pooling", data = tbv.p.dat[tbv.p.dat$RID %in% tbv_dem,], na.action = na.omit)
G <- length(unique(tbv.p.dat[tbv.p.dat$RID %in% tbv_dem,]$RID))
N <- length(tbv.p.dat[tbv.p.dat$RID %in% tbv_dem,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
tbvrel.dp.plm.res.dem <- coeftest(pm1, vcov = firm_c_vcov)

# HC
pm1 <- plm(HC~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex+ICV, model = "pooling", data = hc.p.dat[hc.p.dat$RID %in% hc_dem,], na.action = na.omit)
G <- length(unique(hc.p.dat[hc.p.dat$RID %in% hc_dem,]$RID))
N <- length(hc.p.dat[hc.p.dat$RID %in% hc_dem,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
hcrel.dp.plm.res.dem <- coeftest(pm1, vcov = firm_c_vcov)

# WMH
pm1 <- plm(wmh~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = wmh.p.dat[wmh.p.dat$RID %in% wmh_dem,], na.action = na.omit)
G <- length(unique(wmh.p.dat[wmh.p.dat$RID %in% wmh_dem,]$RID))
N <- length(wmh.p.dat[wmh.p.dat$RID %in% wmh_dem,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmh.dp.plm.res.dem <- coeftest(pm1, vcov = firm_c_vcov)

# CT
pm1 <- plm(mean_thickness.combat~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.p.dat[fs.p.dat$RID %in% fs_dem,], na.action = na.omit)
G <- length(unique(fs.p.dat[fs.p.dat$RID %in% fs_dem,]$RID))
N <- length(fs.p.dat[fs.p.dat$RID %in% fs_dem,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
ct.dp.plm.res.dem <- coeftest(pm1, vcov = firm_c_vcov)

# SA
pm1 <- plm(sa.combat~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.p.dat[fs.p.dat$RID %in% fs_dem,], na.action = na.omit)
G <- length(unique(fs.p.dat[fs.p.dat$RID %in% fs_dem,]$RID))
N <- length(fs.p.dat[fs.p.dat$RID %in% fs_dem,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
sa.dp.plm.res.dem <- coeftest(pm1, vcov = firm_c_vcov)

# WMHypo - added after next section
pm1 <- plm(wmhypo.combat.log~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.wmhypo.p.dat[fs.wmhypo.p.dat$RID %in% fs_dem,], na.action = na.omit)
G <- length(unique(fs.wmhypo.p.dat[fs.wmhypo.p.dat$RID %in% fs_dem,]$RID))
N <- length(fs.wmhypo.p.dat[fs.wmhypo.p.dat$RID %in% fs_dem,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmhypo.dp.plm.res.dem <- coeftest(pm1, vcov = firm_c_vcov)



# MCI
tbv_mci <- tbv.p.dat$RID[tbv.p.dat$DX=='MCI']
hc_mci <- hc.p.dat$RID[hc.p.dat$DX=='MCI']
wmh_mci <- wmh.p.dat$RID[wmh_paired_values_dx$DX=='MCI']
fs_mci <- fs.p.dat$RID[fs_paired_values_dx$DX=='MCI']

# TBV
pm1 <- plm(TBV~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex+ICV, model = "pooling", data = tbv.p.dat[tbv.p.dat$RID %in% tbv_mci,], na.action = na.omit)
G <- length(unique(tbv.p.dat[tbv.p.dat$RID %in% tbv_mci,]$RID))
N <- length(tbv.p.dat[tbv.p.dat$RID %in% tbv_mci,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
tbvrel.dp.plm.res.mci <- coeftest(pm1, vcov = firm_c_vcov)

# HC
pm1 <- plm(HC~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex+ICV, model = "pooling", data = hc.p.dat[hc.p.dat$RID %in% hc_mci,], na.action = na.omit)
G <- length(unique(hc.p.dat[hc.p.dat$RID %in% hc_mci,]$RID))
N <- length(hc.p.dat[hc.p.dat$RID %in% hc_mci,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
hcrel.dp.plm.res.mci <- coeftest(pm1, vcov = firm_c_vcov)

# WMH
pm1 <- plm(wmh~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = wmh.p.dat[wmh.p.dat$RID %in% wmh_mci,], na.action = na.omit)
G <- length(unique(wmh.p.dat[wmh.p.dat$RID %in% wmh_mci,]$RID))
N <- length(wmh.p.dat[wmh.p.dat$RID %in% wmh_mci,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmh.dp.plm.res.mci <- coeftest(pm1, vcov = firm_c_vcov)

# CT
pm1 <- plm(mean_thickness.combat~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.p.dat[fs.p.dat$RID %in% fs_mci,], na.action = na.omit)
G <- length(unique(fs.p.dat[fs.p.dat$RID %in% fs_mci,]$RID))
N <- length(fs.p.dat[fs.p.dat$RID %in% fs_mci,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
ct.dp.plm.res.mci <- coeftest(pm1, vcov = firm_c_vcov)

# SA
pm1 <- plm(sa.combat~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.p.dat[fs.p.dat$RID %in% fs_mci,], na.action = na.omit)
G <- length(unique(fs.p.dat[fs.p.dat$RID %in% fs_mci,]$RID))
N <- length(fs.p.dat[fs.p.dat$RID %in% fs_mci,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
sa.dp.plm.res.mci <- coeftest(pm1, vcov = firm_c_vcov)

# WMHypo - added after next section
pm1 <- plm(wmhypo.combat.log~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.wmhypo.p.dat[fs.wmhypo.p.dat$RID %in% fs_mci,], na.action = na.omit)
G <- length(unique(fs.wmhypo.p.dat[fs.wmhypo.p.dat$RID %in% fs_mci,]$RID))
N <- length(fs.wmhypo.p.dat[fs.wmhypo.p.dat$RID %in% fs_mci,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmhypo.dp.plm.res.mci <- coeftest(pm1, vcov = firm_c_vcov)




# CN
tbv_cn <- tbv.p.dat$RID[tbv.p.dat$DX=='CN']
hc_cn <- hc.p.dat$RID[hc.p.dat$DX=='CN']
wmh_cn <- wmh.p.dat$RID[wmh_paired_values_dx$DX=='CN']
fs_cn <- fs.p.dat$RID[fs_paired_values_dx$DX=='CN']

# TBV
pm1 <- plm(TBV~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex+ICV, model = "pooling", data = tbv.p.dat[tbv.p.dat$RID %in% tbv_cn,], na.action = na.omit)
G <- length(unique(tbv.p.dat[tbv.p.dat$RID %in% tbv_cn,]$RID))
N <- length(tbv.p.dat[tbv.p.dat$RID %in% tbv_cn,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
tbvrel.dp.plm.res.cn <- coeftest(pm1, vcov = firm_c_vcov)

# HC
pm1 <- plm(HC~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex+ICV, model = "pooling", data = hc.p.dat[hc.p.dat$RID %in% hc_cn,], na.action = na.omit)
G <- length(unique(hc.p.dat[hc.p.dat$RID %in% hc_cn,]$RID))
N <- length(hc.p.dat[hc.p.dat$RID %in% hc_cn,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
hcrel.dp.plm.res.cn <- coeftest(pm1, vcov = firm_c_vcov)

# WMH
pm1 <- plm(wmh~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = wmh.p.dat[wmh.p.dat$RID %in% wmh_cn,], na.action = na.omit)
G <- length(unique(wmh.p.dat[wmh.p.dat$RID %in% wmh_cn,]$RID))
N <- length(wmh.p.dat[wmh.p.dat$RID %in% wmh_cn,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmh.dp.plm.res.cn <- coeftest(pm1, vcov = firm_c_vcov)

# CT
pm1 <- plm(mean_thickness.combat~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.p.dat[fs.p.dat$RID %in% fs_cn,], na.action = na.omit)
G <- length(unique(fs.p.dat[fs.p.dat$RID %in% fs_cn,]$RID))
N <- length(fs.p.dat[fs.p.dat$RID %in% fs_cn,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
ct.dp.plm.res.cn <- coeftest(pm1, vcov = firm_c_vcov)

# SA
pm1 <- plm(sa.combat~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.p.dat[fs.p.dat$RID %in% fs_cn,], na.action = na.omit)
G <- length(unique(fs.p.dat[fs.p.dat$RID %in% fs_cn,]$RID))
N <- length(fs.p.dat[fs.p.dat$RID %in% fs_cn,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
sa.dp.plm.res.cn <- coeftest(pm1, vcov = firm_c_vcov)

# WMHypo - added after next section
pm1 <- plm(wmhypo.combat.log~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex, model = "pooling", data = fs.wmhypo.p.dat[fs.wmhypo.p.dat$RID %in% fs_cn,], na.action = na.omit)
G <- length(unique(fs.wmhypo.p.dat[fs.wmhypo.p.dat$RID %in% fs_cn,]$RID))
N <- length(fs.wmhypo.p.dat[fs.wmhypo.p.dat$RID %in% fs_cn,]$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmhypo.dp.plm.res.cn <- coeftest(pm1, vcov = firm_c_vcov)


# visualization stratified by DX

# why is WMH flipping signs - i guess the covariates have a slightly different effect depending on the overall sample?

pm1_cn <- plm(wmh~zPoAm45+mean_age+sex, model = "pooling", data = wmh.p.dat[wmh.p.dat$RID %in% wmh_cn,], na.action = na.omit)
pm1_mci <- plm(wmh~zPoAm45+mean_age+sex, model = "pooling", data = wmh.p.dat[wmh.p.dat$RID %in% wmh_mci,], na.action = na.omit)
pm1_dem <- plm(wmh~zPoAm45+mean_age+sex, model = "pooling", data = wmh.p.dat[wmh.p.dat$RID %in% wmh_dem,], na.action = na.omit)
pm1_all <- plm(wmh~zPoAm45+mean_age+sex, model = "pooling", data = wmh.p.dat, na.action = na.omit)


adni_dx_strat_table <- data.frame(imaging = rep(c('TBV', 'HC', 'WMHyper', 'CT', 'SA', 'WMHypo'), 3),
                                    dx = c(rep('CN', 6), rep('MCI',6), rep('Dementia', 6)),
                                    beta = c(tbvrel.dp.plm.res.cn[2,1],  hcrel.dp.plm.res.cn[2,1], wmh.dp.plm.res.cn[2,1],  ct.dp.plm.res.cn[2,1], sa.dp.plm.res.cn[2,1],  wmhypo.dp.plm.res.cn[2,1],
                                             tbvrel.dp.plm.res.mci[2,1],  hcrel.dp.plm.res.mci[2,1], wmh.dp.plm.res.mci[2,1],  ct.dp.plm.res.mci[2,1], sa.dp.plm.res.mci[2,1],  wmhypo.dp.plm.res.mci[2,1],
                                             tbvrel.dp.plm.res.dem[2,1],  hcrel.dp.plm.res.dem[2,1], wmh.dp.plm.res.dem[2,1],  ct.dp.plm.res.dem[2,1], sa.dp.plm.res.dem[2,1],  wmhypo.dp.plm.res.dem[2,1]),
                                    stderror = c(tbvrel.dp.plm.res.cn[2,2],  hcrel.dp.plm.res.cn[2,2], wmh.dp.plm.res.cn[2,2],  ct.dp.plm.res.cn[2,2], sa.dp.plm.res.cn[2,2],  wmhypo.dp.plm.res.cn[2,2],
                                                 tbvrel.dp.plm.res.mci[2,2],  hcrel.dp.plm.res.mci[2,2], wmh.dp.plm.res.mci[2,2],  ct.dp.plm.res.mci[2,2], sa.dp.plm.res.mci[2,2],  wmhypo.dp.plm.res.mci[2,2],
                                                 tbvrel.dp.plm.res.dem[2,2],  hcrel.dp.plm.res.dem[2,2], wmh.dp.plm.res.dem[2,2],  ct.dp.plm.res.dem[2,2], sa.dp.plm.res.dem[2,2],  wmhypo.dp.plm.res.dem[2,2]))

# point range plot

imaging_level_order_strat <- rev(c('TBV', 'HC', 'WMHypo', 'CT', 'SA', 'WMHyper'))
dx_levels <- c('CN', 'MCI', 'Dementia')
adni_dx_colors <- c('chartreuse3', 'orange2', 'brown3')

adni_dx_strat <- ggplot(adni_dx_strat_table, 
                             aes(x=factor(imaging, imaging_level_order_strat), 
                                 y=beta, color = factor(dx, levels=dx_levels), group=factor(dx, levels = dx_levels))) + 
  geom_pointrange(aes(ymin=beta-(1.96*stderror), ymax=beta+(1.96*stderror)),position=position_dodge(width=0.6)) +
  geom_hline(yintercept = 0)+
  labs(title = 'ADNI') +
  coord_flip() +
  theme_classic()+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 0),
        axis.title.x = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        title = element_text(size=15))+
  guides(color = guide_legend(reverse = TRUE, title = 'Diagnosis')) +
  scale_color_manual(values=c(rep(adni_dx_colors,3)))


######################################################
############ METHYLATION SCAN TIME DIFF ##############
######################################################

tbv.p.dat$time_diff <- interval(ymd(tbv.p.dat$datedrawn), ymd(tbv.p.dat$datescanned)) %/% days()
hc.p.dat$time_diff <- interval(ymd(hc.p.dat$datedrawn), ymd(hc.p.dat$datescanned)) %/% days()
wmh.p.dat$time_diff <- interval(ymd(wmh.p.dat$datedrawn), ymd(wmh.p.dat$datescanned)) %/% days()
fs.p.dat$time_diff <- interval(mdy(fs.p.dat$datedrawn), ymd(fs.p.dat$datescanned)) %/% days()
fs.wmhypo.p.dat$time_diff <- interval(mdy(fs.wmhypo.p.dat$datedrawn), ymd(fs.wmhypo.p.dat$datescanned)) %/% days()

# TBV
pm1 <- plm(TBV~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex+ICV+time_diff, model = "pooling", data = tbv.p.dat, na.action = na.omit)
G <- length(unique(tbv.p.dat$RID))
N <- length(tbv.p.dat$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
tbvrel.dp.plm.res.td <- coeftest(pm1, vcov = firm_c_vcov)

# HC
pm1 <- plm(HC~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex+ICV+time_diff, model = "pooling", data = hc.p.dat, na.action = na.omit)
G <- length(unique(hc.p.dat$RID))
N <- length(hc.p.dat$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
hcrel.dp.plm.res.td <- coeftest(pm1, vcov = firm_c_vcov)

# WMH
pm1 <- plm(wmh~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex+time_diff, model = "pooling", data = wmh.p.dat, na.action = na.omit)
G <- length(unique(wmh.p.dat$RID))
N <- length(wmh.p.dat$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmh.dp.plm.res.td <- coeftest(pm1, vcov = firm_c_vcov)

# CT
pm1 <- plm(mean_thickness.combat~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex+time_diff, model = "pooling", data = fs.p.dat, na.action = na.omit)
G <- length(unique(fs.p.dat$RID))
N <- length(fs.p.dat$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
ct.dp.plm.res.td <- coeftest(pm1, vcov = firm_c_vcov)

# SA
pm1 <- plm(sa.combat~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex+time_diff, model = "pooling", data = fs.p.dat, na.action = na.omit)
G <- length(unique(fs.p.dat$RID))
N <- length(fs.p.dat$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
sa.dp.plm.res.td <- coeftest(pm1, vcov = firm_c_vcov)

# WMHypo - added after next section
pm1 <- plm(wmhypo.combat.log~zPoAm45+sex+mean_age+sex*mean_age+mean_age_sq+mean_age_sq*sex+time_diff, model = "pooling", data = fs.wmhypo.p.dat, na.action = na.omit)
G <- length(unique(fs.wmhypo.p.dat$RID))
N <- length(fs.wmhypo.p.dat$RID)
dfa <- (G/(G - 1)) * (N - 1)/pm1$df.residual
firm_c_vcov <- dfa * vcovHC(pm1, type = "HC0", cluster = "group", adjust = T)
wmhypo.dp.plm.res.td <- coeftest(pm1, vcov = firm_c_vcov)


adni_td_table <- data.frame(imaging = rep(c('TBV', 'HC', 'WMHyper', 'CT', 'SA', 'WMHypo'), 2),
                       clock = c(rep('Original', 6), rep('Controlling for time difference', 6)),
                       beta = c(tbvrel.dp.plm.res[2,1],  hcrel.dp.plm.res[2,1],        wmh.dp.plm.res[2,1],   ct.dp.plm.res[2,1],       sa.dp.plm.res[2,1], wmhypo.dp.plm.res[2,1],  
                                tbvrel.dp.plm.res.td[2,1], hcrel.dp.plm.res.td[2,1],       wmh.dp.plm.res.td[2,1],   ct.dp.plm.res.td[2,1],    sa.dp.plm.res.td[2,1],  wmhypo.dp.plm.res.td[2,1]),
                       stderror = c(tbvrel.dp.plm.res[2,2],  hcrel.dp.plm.res[2,2],        wmh.dp.plm.res[2,2],   ct.dp.plm.res[2,2],       sa.dp.plm.res[2,2], wmhypo.dp.plm.res[2,2],  
                                    tbvrel.dp.plm.res.td[2,2], hcrel.dp.plm.res.td[2,2],       wmh.dp.plm.res.td[2,2],   ct.dp.plm.res.td[2,2],    sa.dp.plm.res.td[2,2],  wmhypo.dp.plm.res.td[2,2])
)



# visualization

adni_tbv_scanint_hist <- ggplot(data = data.frame(tbv.p.dat), aes(x = time_diff)) + 
  geom_histogram(fill='#619CFF') +
  ylim(0,1500)+
  xlim(-200, 200)+
  labs(title = "ADNI - TBV",
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

adni_hc_scanint_hist <- ggplot(data = data.frame(hc.p.dat), aes(x = time_diff)) + 
  geom_histogram(fill='#619CFF') +
  ylim(0,1500)+
  xlim(-200, 200)+
  labs(title = "ADNI - HC",
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

adni_wmh_scanint_hist <- ggplot(data = data.frame(wmh.p.dat), aes(x = time_diff)) + 
  geom_histogram(fill='#619CFF') +
  ylim(0,1500)+
  xlim(-200, 200)+
  labs(title = "ADNI - WMHypo",
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

adni_fs_scanint_hist <- ggplot(data = data.frame(fs.p.dat), aes(x = time_diff)) + 
  geom_histogram(fill='#619CFF') +
  ylim(0,1500)+
  xlim(-200, 200)+
  labs(title = "ADNI - FreeSurfer",
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

grid.arrange(adni_tbv_scanint_hist, adni_hc_scanint_hist, adni_fs_scanint_hist, nrow =1)

######################################################
################# RACE DEMOGRAPHICS ##################
######################################################

adni_race <- adnimerge[,c('RID', 'PTETHCAT', 'PTRACCAT')]
adni_race <- adni_race[!duplicated(adni_race$RID),]

meth_race <- adni_race[adni_race$RID %in% adni_meth$RID,]
tbv_race <- adni_race[adni_race$RID %in% tbv_paired_values$RID,]
hc_race <- adni_race[adni_race$RID %in% hc_paired_values$RID,]
wmh_race <- adni_race[adni_race$RID %in% wmh_paired_values$RID,]
fs_race <- adni_race[adni_race$RID %in% fs_paired_values$RID,]

table(wmh_race$PTRACCAT)
table(wmh_race$PTETHCAT)

sum(as.numeric(meth_race$PTETHCAT == 'Not Hisp/Latino')*as.numeric(meth_race$PTRACCAT == 'White'))





###############################################################
#################### LONGITUDINAL ANALYSIS ####################
###############################################################

# merge in WMHypo with other FS variables
if (identical(methFScomb_harmed$RID, adni_WMHypo_harm$RID) == TRUE){
  print('IDENTICAL - COMBINED')
  methFScomb_harmed$WMHypo.combat <- adni_WMHypo_harm$ST128SV.combat
  methFScomb_harmed$WMHypo.combat_log <- log(adni_WMHypo_harm$ST128SV.combat)
} else {
  print('MISMATCH')
}

# longitudinal function

adni_long <- function(var, clock_pace, clock_acc) {
  
  
  adni_methclocks <- as.data.frame(read.delim('/Users/ew198/Documents/methylation/from_karen/ADNI_methclocks_pluspheno.txt'))
  
  adni_meth <- read.delim("/Users/ew198/Documents/methylation/adni_methIDlist.txt", header = TRUE)
  adni_meth$AGE_MONTH <- adni_meth$AGE*12
  adnimerge <- read.csv("/Users/ew198/Documents/methylation/forEthan_070622/Data/ADNIMERGE.csv")
  
  adni_baseline_dates <- data.frame(RID = unique(adni_meth$RID), 
                                    baseline_date = rep(NA, length(unique(adni_meth$RID))),
                                    age_at_baseline = rep(NA, length(unique(adni_meth$RID))))
  
  x <- 1
  for (i in unique(adni_meth$RID)){
    temp_df <- as.data.frame(adnimerge[adnimerge$RID == i,])
    temp_row <- temp_df[temp_df$VISCODE == 'bl',]
    adni_baseline_dates$baseline_date[x] <- as.character(temp_row$EXAMDATE)
    adni_baseline_dates$age_at_baseline[x] <- temp_row$AGE*12
    x <- x + 1
  }
  
  if (var == 'wmhVol_log'){
    methMerge_temp <- wmh_merge
    methMerge_temp$AGE_MONTH <- methMerge_temp$age_at_baseline
  } else if (var == 'mean_thickness.combat'){
    methMerge_temp <- methFScomb_harmed
  } else if (var == 'total_sa.combat'){
    methMerge_temp <- methFScomb_harmed
  } else if (var == 'WMHypo.combat_log'){
    methMerge_temp <- methFScomb_harmed
  } else {
    methMerge_temp <- methMerge
  }
  
  
  adni_long_dates <- data.frame(RID = unique(methMerge_temp$RID), 
                                DateScanned1 = rep(NA, length(unique(methMerge_temp$RID))),
                                DateScanned2 = rep(NA, length(unique(methMerge_temp$RID))),
                                DateScanned3 = rep(NA, length(unique(methMerge_temp$RID))),
                                DateScanned4 = rep(NA, length(unique(methMerge_temp$RID))),
                                DateScanned5 = rep(NA, length(unique(methMerge_temp$RID))),
                                DateScanned6 = rep(NA, length(unique(methMerge_temp$RID))),
                                DateScanned7 = rep(NA, length(unique(methMerge_temp$RID))),
                                DateScanned8 = rep(NA, length(unique(methMerge_temp$RID))),
                                DateScanned9 = rep(NA, length(unique(methMerge_temp$RID))),
                                DateScanned10 = rep(NA, length(unique(methMerge_temp$RID))),
                                DateScanned11 = rep(NA, length(unique(methMerge_temp$RID))),
                                DateScanned12 = rep(NA, length(unique(methMerge_temp$RID))),
                                DateScanned13 = rep(NA, length(unique(methMerge_temp$RID))),
                                DateScanned14 = rep(NA, length(unique(methMerge_temp$RID))),
                                DateScanned15 = rep(NA, length(unique(methMerge_temp$RID))),
                                DateScanned16 = rep(NA, length(unique(methMerge_temp$RID))),
                                DateScanned17 = rep(NA, length(unique(methMerge_temp$RID))),
                                DateScanned18 = rep(NA, length(unique(methMerge_temp$RID))),
                                
                                age_at_scan1 = rep(NA, length(unique(methMerge_temp$RID))),
                                age_at_scan2 = rep(NA, length(unique(methMerge_temp$RID))),
                                age_at_scan3 = rep(NA, length(unique(methMerge_temp$RID))),
                                age_at_scan4 = rep(NA, length(unique(methMerge_temp$RID))),
                                age_at_scan5 = rep(NA, length(unique(methMerge_temp$RID))),
                                age_at_scan6 = rep(NA, length(unique(methMerge_temp$RID))),
                                age_at_scan7 = rep(NA, length(unique(methMerge_temp$RID))),
                                age_at_scan8 = rep(NA, length(unique(methMerge_temp$RID))),
                                age_at_scan9 = rep(NA, length(unique(methMerge_temp$RID))),
                                age_at_scan10 = rep(NA, length(unique(methMerge_temp$RID))),
                                age_at_scan11 = rep(NA, length(unique(methMerge_temp$RID))),
                                age_at_scan12 = rep(NA, length(unique(methMerge_temp$RID))),
                                age_at_scan13 = rep(NA, length(unique(methMerge_temp$RID))),
                                age_at_scan14 = rep(NA, length(unique(methMerge_temp$RID))),
                                age_at_scan15 = rep(NA, length(unique(methMerge_temp$RID))),
                                age_at_scan16 = rep(NA, length(unique(methMerge_temp$RID))),
                                age_at_scan17 = rep(NA, length(unique(methMerge_temp$RID))),
                                age_at_scan18 = rep(NA, length(unique(methMerge_temp$RID))),
                                
                                var_at_scan1 = rep(NA, length(unique(methMerge_temp$RID))),
                                var_at_scan2 = rep(NA, length(unique(methMerge_temp$RID))),
                                var_at_scan3 = rep(NA, length(unique(methMerge_temp$RID))),
                                var_at_scan4 = rep(NA, length(unique(methMerge_temp$RID))),
                                var_at_scan5 = rep(NA, length(unique(methMerge_temp$RID))),
                                var_at_scan6 = rep(NA, length(unique(methMerge_temp$RID))),
                                var_at_scan7 = rep(NA, length(unique(methMerge_temp$RID))),
                                var_at_scan8 = rep(NA, length(unique(methMerge_temp$RID))),
                                var_at_scan9 = rep(NA, length(unique(methMerge_temp$RID))),
                                var_at_scan10 = rep(NA, length(unique(methMerge_temp$RID))),
                                var_at_scan11 = rep(NA, length(unique(methMerge_temp$RID))),
                                var_at_scan12 = rep(NA, length(unique(methMerge_temp$RID))),
                                var_at_scan13 = rep(NA, length(unique(methMerge_temp$RID))),
                                var_at_scan14 = rep(NA, length(unique(methMerge_temp$RID))),
                                var_at_scan15 = rep(NA, length(unique(methMerge_temp$RID))),
                                var_at_scan16 = rep(NA, length(unique(methMerge_temp$RID))),
                                var_at_scan17 = rep(NA, length(unique(methMerge_temp$RID))),
                                var_at_scan18 = rep(NA, length(unique(methMerge_temp$RID))),
                                
                                icv_at_scan1 = rep(NA, length(unique(methMerge_temp$RID))),
                                icv_at_scan2 = rep(NA, length(unique(methMerge_temp$RID))),
                                icv_at_scan3 = rep(NA, length(unique(methMerge_temp$RID))),
                                icv_at_scan4 = rep(NA, length(unique(methMerge_temp$RID))),
                                icv_at_scan5 = rep(NA, length(unique(methMerge_temp$RID))),
                                icv_at_scan6 = rep(NA, length(unique(methMerge_temp$RID))),
                                icv_at_scan7 = rep(NA, length(unique(methMerge_temp$RID))),
                                icv_at_scan8 = rep(NA, length(unique(methMerge_temp$RID))),
                                icv_at_scan9 = rep(NA, length(unique(methMerge_temp$RID))),
                                icv_at_scan10 = rep(NA, length(unique(methMerge_temp$RID))),
                                icv_at_scan11 = rep(NA, length(unique(methMerge_temp$RID))),
                                icv_at_scan12 = rep(NA, length(unique(methMerge_temp$RID))),
                                icv_at_scan13 = rep(NA, length(unique(methMerge_temp$RID))),
                                icv_at_scan14 = rep(NA, length(unique(methMerge_temp$RID))),
                                icv_at_scan15 = rep(NA, length(unique(methMerge_temp$RID))),
                                icv_at_scan16 = rep(NA, length(unique(methMerge_temp$RID))),
                                icv_at_scan17 = rep(NA, length(unique(methMerge_temp$RID))),
                                icv_at_scan18 = rep(NA, length(unique(methMerge_temp$RID))),
                                
                                
                                first_scan = rep(NA, length(unique(methMerge_temp$RID))),
                                last_scan = rep(NA, length(unique(methMerge_temp$RID))),
                                age_at_last_scan = rep(NA, length(unique(methMerge_temp$RID))), 
                                fs_version = rep(NA, length(unique(methMerge_temp$RID))) )
  
  methMerge_temp[,match(print(var), colnames(methMerge_temp))] <- as.numeric(methMerge_temp[,match(print(var), colnames(methMerge_temp))])
  
  x <- 1
  for (i in unique(methMerge_temp$RID)){
    temp_df <- as.data.frame(methMerge_temp[methMerge_temp$RID == i,])
    temp_df <- temp_df[!is.na(temp_df[,match(print(var), colnames(temp_df))]),]
    temp_df <- temp_df[order(temp_df$EXAMDATE),]
    if (var != 'mean_thickness.combat' && var != 'total_sa.combat' && var != 'wmhVol_log' && var != 'WMHypo.combat_log'){
      temp_df$AGE_MONTH <- 12*temp_df$AGE
    } else if (var == 'mean_thickness.combat' | var == 'total_sa.combat' | var == 'wmhVol_log' | var == 'WMHypo.combat_log'){
      temp_df$AGE_MONTH <- temp_df$age_at_baseline
      temp_df$EXAMDATE_bl <- temp_df$baseline_date
    }
    timepoints <- nrow(temp_df)
    if (timepoints > 0){
      if (var != 'wmhVol_log'){
        for (t in 1:timepoints){
          adni_long_dates[x, (t+1)] <- as.character(temp_df$EXAMDATE[t])
          adni_long_dates[x, (t+19)] <- temp_df$AGE_MONTH[1] + interval(ymd(temp_df$EXAMDATE_bl[1]), ymd(temp_df$EXAMDATE[t])) %/% months(1)
          adni_long_dates[x, (t+37)] <- temp_df[,match(print(var), colnames(temp_df))][t]
          if (var == "WholeBrain" | var == "Hippocampus"){
            adni_long_dates[x, (t+55)] <- temp_df[,'ICV'][t]
          }
        }
        adni_long_dates[x, 74] <- as.character(temp_df$EXAMDATE[1])
        adni_long_dates[x, 75] <- as.character(temp_df$EXAMDATE[timepoints])
        adni_long_dates[x, 76] <- temp_df$AGE_MONTH[1] + interval(ymd(temp_df$EXAMDATE_bl[1]), ymd(temp_df$EXAMDATE[timepoints])) %/% months(1)
        if (var == 'mean_thickness.combat' | var == 'total_sa.combat' | var == 'WMHypo.combat_log'){
          adni_long_dates[x, 77] <- temp_df$FSVERSION[1]
        }
        x <- x+1
      } else if (var == 'wmhVol_log'){
        for (t in 1:timepoints){
          adni_long_dates[x, (t+1)] <- as.character(temp_df$EXAMDATE[t])
          adni_long_dates[x, (t+19)] <- temp_df$AGE_MONTH[1] + interval(ymd(temp_df$baseline_date[1]), ymd(temp_df$EXAMDATE[t])) %/% months(1)
          adni_long_dates[x, (t+37)] <- temp_df[,match(print(var), colnames(temp_df))][t]
          adni_long_dates[x, (t+55)] <- temp_df[,'TOTAL_BRAIN'][t]
        }
        adni_long_dates[x, 74] <- as.character(temp_df$EXAMDATE[1])
        adni_long_dates[x, 75] <- as.character(temp_df$EXAMDATE[timepoints])
        adni_long_dates[x, 76] <- temp_df$AGE_MONTH[1] + interval(ymd(temp_df$baseline_date[1]), ymd(temp_df$EXAMDATE[timepoints])) %/% months(1)
        x <- x+1
      }
    } else if (timepoints == 0){
      x <- x + 1
    }
  }
  
  # age plot age at whole brain volume collection
  adni_long_dates <- merge(adni_long_dates, 
                           data.frame(RID = adni_meth_datedrawn$RID, meth_timepoint1 = adni_meth_datedrawn$age_at_timepoint1),
                           by = 'RID')
  adni_long_dates <- merge(adni_long_dates, adni_baseline_dates, by = 'RID')
  adni_long_dates <- adni_long_dates[order(adni_long_dates$age_at_baseline),]
  # adni_long_dates <- data.frame(adni_long_dates, plot_id = c(1:649))
  
  adni_var_long_date <- melt(adni_long_dates, id.vars = 'RID', measure.vars = c('DateScanned1', 
                                                                                'DateScanned2', 
                                                                                'DateScanned3', 
                                                                                'DateScanned4', 
                                                                                'DateScanned5', 
                                                                                'DateScanned6', 
                                                                                'DateScanned7', 
                                                                                'DateScanned8', 
                                                                                'DateScanned9', 
                                                                                'DateScanned10', 
                                                                                'DateScanned11', 
                                                                                'DateScanned12', 
                                                                                'DateScanned13', 
                                                                                'DateScanned14',
                                                                                'DateScanned15',
                                                                                'DateScanned16',
                                                                                'DateScanned17',
                                                                                'DateScanned18'))
  
  adni_var_long_age <- melt(adni_long_dates, id.vars = 'RID', measure.vars = c('age_at_scan1', 
                                                                               'age_at_scan2', 
                                                                               'age_at_scan3', 
                                                                               'age_at_scan4', 
                                                                               'age_at_scan5', 
                                                                               'age_at_scan6', 
                                                                               'age_at_scan7', 
                                                                               'age_at_scan8', 
                                                                               'age_at_scan9', 
                                                                               'age_at_scan10', 
                                                                               'age_at_scan11', 
                                                                               'age_at_scan12', 
                                                                               'age_at_scan13', 
                                                                               'age_at_scan14',
                                                                               'age_at_scan15',
                                                                               'age_at_scan16',
                                                                               'age_at_scan17',
                                                                               'age_at_scan18'))
  
  adni_var_long_var <- melt(adni_long_dates, id.vars = 'RID', measure.vars = c('var_at_scan1', 
                                                                               'var_at_scan2', 
                                                                               'var_at_scan3', 
                                                                               'var_at_scan4', 
                                                                               'var_at_scan5', 
                                                                               'var_at_scan6', 
                                                                               'var_at_scan7', 
                                                                               'var_at_scan8', 
                                                                               'var_at_scan9', 
                                                                               'var_at_scan10', 
                                                                               'var_at_scan11', 
                                                                               'var_at_scan12', 
                                                                               'var_at_scan13', 
                                                                               'var_at_scan14',
                                                                               'var_at_scan15',
                                                                               'var_at_scan16',
                                                                               'var_at_scan17',
                                                                               'var_at_scan18'))
  
  adni_var_long_icv <- melt(adni_long_dates, id.vars = 'RID', measure.vars = c('icv_at_scan1', 
                                                                               'icv_at_scan2', 
                                                                               'icv_at_scan3', 
                                                                               'icv_at_scan4', 
                                                                               'icv_at_scan5', 
                                                                               'icv_at_scan6', 
                                                                               'icv_at_scan7', 
                                                                               'icv_at_scan8', 
                                                                               'icv_at_scan9', 
                                                                               'icv_at_scan10', 
                                                                               'icv_at_scan11', 
                                                                               'icv_at_scan12', 
                                                                               'icv_at_scan13', 
                                                                               'icv_at_scan14',
                                                                               'icv_at_scan15',
                                                                               'icv_at_scan16',
                                                                               'icv_at_scan17',
                                                                               'icv_at_scan18'))
  
  if (var == "WholeBrain" | var == "Hippocampus" | var == "wmhVol_log"){
    adni_meth_var_long <- data.frame(adni_var_long_date, adni_var_long_age, adni_var_long_var, adni_var_long_icv)
  } else {
    adni_meth_var_long <- data.frame(adni_var_long_date, adni_var_long_age, adni_var_long_var)
  }
  
  adni_meth_baseline_values <- data.frame(RID = rep(0, length(unique(adni_methclocks$RID))),
                                          sex = rep(0, length(unique(adni_methclocks$RID))),
                                          baseline_pace = rep(0, length(unique(adni_methclocks$RID))),
                                          baseline_acc = rep(0, length(unique(adni_methclocks$RID))))
  
  x <- 1
  for (i in unique(adni_methclocks$RID)){
    temp_df <- as.data.frame(adni_methclocks[adni_methclocks$RID == i,])
    temp_df <- temp_df[order(temp_df$AGE),]
    adni_meth_baseline_values$RID[x] <- temp_df$RID[1]
    adni_meth_baseline_values$sex[x] <- temp_df$PTGENDER[1]
    adni_meth_baseline_values$baseline_pace[x] <- temp_df[1,match(print(clock_pace), colnames(temp_df))]
    adni_meth_baseline_values$baseline_acc[x] <- temp_df[1,match(print(clock_acc), colnames(temp_df))]
    x <- x+1
  }
  
  adni_meth_baseline_values_dates <- merge(adni_meth_baseline_values, adni_baseline_dates, by = 'RID')
  
  meth_predict_var <- merge(adni_meth_var_long, adni_meth_baseline_values, by = 'RID')
  names(meth_predict_var)[names(meth_predict_var) == 'value'] <- 'scandate'
  names(meth_predict_var)[names(meth_predict_var) == 'value.1'] <- 'scanage'
  names(meth_predict_var)[names(meth_predict_var) == 'value.2'] <- 'var'
  meth_predict_var <- meth_predict_var[complete.cases(meth_predict_var),]
  
  # filtering for pre-methylation scans
  meth_predict_var_post <- meth_predict_var[1,]
  for (i in unique(adni_methclocks$RID)){
    temp_baseline <- adni_methclocks[adni_methclocks$RID==i,]
    temp_baseline <- temp_baseline[order(mdy(temp_baseline$Edate)),]
    temp_df <- meth_predict_var[meth_predict_var$RID==i,]
    temp_df <- temp_df[(ymd(temp_df$scandate) >= mdy(temp_baseline$EXAMDATE[1])),]
    meth_predict_var_post <- rbind(meth_predict_var_post, temp_df)
  }
  meth_predict_var_post <- meth_predict_var_post[2:nrow(meth_predict_var_post),]
  
  
  meth_predict_var_post <- meth_predict_var_post[c(order(meth_predict_var_post$baseline_pace)),]
  color_pace <- data.frame(RID = unique(meth_predict_var_post$RID), plot_color_pace = plasma(length(unique(meth_predict_var_post$RID)))[1:length(unique(meth_predict_var_post$RID))])
  meth_predict_var_post <- merge(meth_predict_var_post, color_pace, by = 'RID')
  
  meth_predict_var_post <- meth_predict_var_post[c(order(meth_predict_var_post$baseline_acc)),]
  color_acc <- data.frame(RID = unique(meth_predict_var_post$RID), plot_color_acc = plasma(length(unique(meth_predict_var_post$RID)))[1:length(unique(meth_predict_var_post$RID))])
  meth_predict_var_post <- merge(meth_predict_var_post, color_acc, by = 'RID')
  
  meth_predict_var_post$scanage_scale <- scale(as.numeric(meth_predict_var_post$scanage))
  meth_predict_var_post$var_scale <- scale(meth_predict_var_post$var)
  if(var == "WholeBrain" | var == "Hippocampus"){
    meth_predict_var_post$icv_scale <- scale(as.numeric(meth_predict_var_post$value.3))
  } else if (var == "wmhVol_log"){
    meth_predict_var_post$tbv_scale <- scale(as.numeric(meth_predict_var_post$value.3))
  }
  
  meth_predict_var_post <- meth_predict_var_post[order(meth_predict_var_post$RID, meth_predict_var_post$scandate),]
  
  return(meth_predict_var_post)
  
}




# dunedinpace

meth_predict_tbv_dunedinpace <- adni_long('WholeBrain', 'DunedinPoAm_45', 'zPoAm45')
meth_predict_hc_dunedinpace <- adni_long('Hippocampus', 'DunedinPoAm_45', 'zPoAm45')
meth_predict_wmh_dunedinpace <- adni_long('wmhVol_log', 'DunedinPoAm_45', 'zPoAm45')
meth_predict_ct_dunedinpace <- adni_long('mean_thickness.combat', 'DunedinPoAm_45', 'zPoAm45')
meth_predict_sa_dunedinpace <- adni_long('total_sa.combat', 'DunedinPoAm_45', 'zPoAm45')
meth_predict_wmhypo_dunedinpace <- adni_long('WMHypo.combat_log', 'DunedinPoAm_45', 'zPoAm45')

# # horvath long
# 
# meth_predict_tbv_horvath <- adni_long('WholeBrain', 'DNAmAge', 'zHorvath')
# meth_predict_hc_horvath <- adni_long('Hippocampus', 'DNAmAge', 'zHorvath')
# meth_predict_wmh_horvath <- adni_long('wmhVol_log', 'DNAmAge', 'zHorvath')
# meth_predict_ct_horvath <- adni_long('mean_thickness.combat', 'DunedinPoAm_45', 'zHorvath')
# meth_predict_sa_horvath <- adni_long('total_sa.combat', 'DunedinPoAm_45', 'zHorvath')
# meth_predict_wmhypo_horvath <- adni_long('WMHypo.combat_log', 'DunedinPoAm_45', 'zHorvath')
# 
# # hannum long
# 
# meth_predict_tbv_hannum <- adni_long('WholeBrain', 'DNAmAgeHannum', 'zHannum')
# meth_predict_hc_hannum <- adni_long('Hippocampus', 'DNAmAgeHannum', 'zHannum')
# meth_predict_wmh_hannum <- adni_long('wmhVol_log', 'DNAmAgeHannum', 'zHannum')
# meth_predict_ct_hannum <- adni_long('mean_thickness.combat', 'DunedinPoAm_45', 'zHannum')
# meth_predict_sa_hannum <- adni_long('total_sa.combat', 'DunedinPoAm_45', 'zHannum')
# meth_predict_wmhypo_hannum <- adni_long('WMHypo.combat_log', 'DunedinPoAm_45', 'zHannum')
# 
# # phenoage long
# 
# meth_predict_tbv_phenoage <- adni_long('WholeBrain', 'DNAmPhenoAge', 'zPheno')
# meth_predict_hc_phenoage <- adni_long('Hippocampus', 'DNAmPhenoAge', 'zPheno')
# meth_predict_wmh_phenoage <- adni_long('wmhVol_log', 'DNAmPhenoAge', 'zPheno')
# meth_predict_ct_phenoage <- adni_long('mean_thickness.combat', 'DunedinPoAm_45', 'zPheno')
# meth_predict_sa_phenoage <- adni_long('total_sa.combat', 'DunedinPoAm_45', 'zPheno')
# meth_predict_wmhypo_phenoage <- adni_long('WMHypo.combat_log', 'DunedinPoAm_45', 'zPheno')
# 
# # grimage long
# 
# meth_predict_tbv_grimage <- adni_long('WholeBrain', 'DNAmGrimAge', 'zGrim')
# meth_predict_hc_grimage <- adni_long('Hippocampus', 'DNAmGrimAge', 'zGrim')
# meth_predict_wmh_grimage <- adni_long('wmhVol_log', 'DNAmGrimAge', 'zGrim')
# meth_predict_ct_grimage <- adni_long('mean_thickness.combat', 'DunedinPoAm_45', 'zGrim')
# meth_predict_sa_grimage <- adni_long('total_sa.combat', 'DunedinPoAm_45', 'zGrim')
# meth_predict_wmhypo_grimage <- adni_long('WMHypo.combat_log', 'DunedinPoAm_45', 'zGrim')


# getting longitudinal diagnostic groups

dx_df <- data.frame(RID = methMerge$RID, EXAMDATE = methMerge$EXAMDATE, Dx = methMerge$DX, Dx_bl = methMerge$DX_bl)
dx_df$Dx[dx_df$Dx == ""] <- NA

# dividing into CN stable, CN -> MCI, CN/MCI -> AD

long_dx <- data.frame(RID = unique(dx_df$RID), group = c(rep(NA, 649)))

x <- 1
for (i in unique(dx_df$RID)){
  temp_df <- as.data.frame(dx_df[dx_df$RID == i,])
  temp_df <- temp_df[!is.na(temp_df$Dx),]
  temp_df <- temp_df[order(temp_df$EXAMDATE),]
  timepoints <- nrow(temp_df)
  if (timepoints > 0){
    if (temp_df$Dx[1] == 'CN' && temp_df$Dx[timepoints] == 'CN'){
      long_dx$group[x] <- 'sCN'
    } else if (temp_df$Dx[1] == 'CN' && temp_df$Dx[timepoints] == 'MCI'){
      long_dx$group[x] <- 'CN_to_MCI'
    } else if (temp_df$Dx[1] == 'MCI' && temp_df$Dx[timepoints] == 'CN'){
      long_dx$group[x] <- 'MCI_to_CN'
    } else if (temp_df$Dx[1] == 'CN' && temp_df$Dx[timepoints] == 'Dementia'){
      long_dx$group[x] <- 'CN_to_Dem'
    } else if (temp_df$Dx[1] == 'MCI' && temp_df$Dx[timepoints] == 'MCI'){
      long_dx$group[x] <- 'sMCI'
    } else if (temp_df$Dx[1] == 'MCI' && temp_df$Dx[timepoints] == 'Dementia'){
      long_dx$group[x] <- 'MCI_to_Dem'
    } else if (temp_df$Dx[1] == 'Dementia' && temp_df$Dx[timepoints] == 'Dementia'){
      long_dx$group[x] <- 'sDem'
    }
  }
  x <- x+1
}


############## longitudinal analyses ############## 

control = lmeControl(msMaxIter = 10000, msMaxEval = 10000)


# TBV

# restrict to three timepoitns
tbv_three_index <- meth_predict_tbv_dunedinpace$RID %in% data.frame(table(meth_predict_tbv_dunedinpace$RID)[table(meth_predict_tbv_dunedinpace$RID) >= 3])$Var1
meth_predict_tbv_dunedinpace_three <- meth_predict_tbv_dunedinpace[tbv_three_index,]

# restrict by diagnosis - step moved up on 8/17
meth_predict_tbv_dunedinpace_three_dx <- left_join(meth_predict_tbv_dunedinpace_three, long_dx, by= 'RID')
meth_predict_tbv_dunedinpace_three_healthy <- meth_predict_tbv_dunedinpace_three_dx[meth_predict_tbv_dunedinpace_three_dx$group == 'sCN',]
meth_predict_tbv_dunedinpace_three_healthy$var_resid <- resid(lm(var_scale~icv_scale, data = meth_predict_tbv_dunedinpace_three_healthy))


# three
# meth_predict_tbv_dunedinpace_three$var_resid <- resid(lm(var_scale~icv_scale, data = meth_predict_tbv_dunedinpace_three))
# # lme() doesn't like this formula for some reason. lme() and lmer() seem to give comparable results though, so using lmer() for this one
tbv_randomslopes_three <- coef(lmer(var_scale~1 + (scanage_scale|RID), data = meth_predict_tbv_dunedinpace_three_healthy))$RID[,1]
tbv_randomslopes_icv_three <- coef(lme(var_resid~1, random=~scanage_scale|RID, control = control, data = meth_predict_tbv_dunedinpace_three_healthy))$scanage_scale


tbv_obs_length_three <- rep(NA, length(unique(meth_predict_tbv_dunedinpace_three_healthy$RID)))
baseline_age_tbv_three <- rep(NA, length(unique(meth_predict_tbv_dunedinpace_three_healthy$RID)))
x <- 1
for (i in unique(meth_predict_tbv_dunedinpace_three_healthy$RID)){
  temp_df <- meth_predict_tbv_dunedinpace_three_healthy[meth_predict_tbv_dunedinpace_three_healthy$RID == i,]
  tbv_obs_length_three[x] <- max(temp_df$scanage) - min(temp_df$scanage)
  baseline_age_tbv_three[x] <- min(temp_df$scanage)
  x <- x + 1
}


random_slopes_tbv_three_df <- data.frame(RID = unique(meth_predict_tbv_dunedinpace_three_healthy$RID), 
                                         sex = meth_predict_tbv_dunedinpace_three_healthy$sex[!duplicated(meth_predict_tbv_dunedinpace_three_healthy$RID)],
                                         mean_age = meth_predict_tbv_dunedinpace_three_healthy$scanage_scale[!duplicated(meth_predict_tbv_dunedinpace_three_healthy$RID)],
                                         baseline_age_tbv = baseline_age_tbv_three,
                                         obs_length = tbv_obs_length_three,
                                         growth_curve = tbv_randomslopes_three, 
                                         growth_curve_icv = tbv_randomslopes_icv_three, 
                                         baseline_acc = unique(meth_predict_tbv_dunedinpace_three_healthy$baseline_acc),
                                         baseline_pace = unique(meth_predict_tbv_dunedinpace_three_healthy$baseline_pace))


# this is now the model we want 8/17
stable_tbv_dunedinpace <- summary(lm(scale(growth_curve_icv)~scale(baseline_acc)+sex+scale(obs_length)+scale(baseline_age_tbv), data = random_slopes_tbv_three_df))






# HC

# restrict to three timepoitns
hc_three_index <- meth_predict_hc_dunedinpace$RID %in% data.frame(table(meth_predict_hc_dunedinpace$RID)[table(meth_predict_hc_dunedinpace$RID) >= 3])$Var1
meth_predict_hc_dunedinpace_three <- meth_predict_hc_dunedinpace[hc_three_index,]

# restrict by diagnosis - step moved up on 8/17
meth_predict_hc_dunedinpace_three_dx <- left_join(meth_predict_hc_dunedinpace_three, long_dx, by= 'RID')
meth_predict_hc_dunedinpace_three_healthy <- meth_predict_hc_dunedinpace_three_dx[meth_predict_hc_dunedinpace_three_dx$group == 'sCN',]
meth_predict_hc_dunedinpace_three_healthy$var_resid <- resid(lm(var_scale~icv_scale, data = meth_predict_hc_dunedinpace_three_healthy))


# three
# meth_predict_hc_dunedinpace_three$var_resid <- resid(lm(var_scale~icv_scale, data = meth_predict_hc_dunedinpace_three))
# # lme() doesn't like this formula for some reason. lme() and lmer() seem to give comparable results though, so using lmer() for this one
hc_randomslopes_three <- coef(lmer(var_scale~1 + (scanage_scale|RID), data = meth_predict_hc_dunedinpace_three_healthy))$RID[,1]
hc_randomslopes_icv_three <- coef(lme(var_resid~1, random=~scanage_scale|RID, control = control, data = meth_predict_hc_dunedinpace_three_healthy))$scanage_scale


hc_obs_length_three <- rep(NA, length(unique(meth_predict_hc_dunedinpace_three_healthy$RID)))
baseline_age_hc_three <- rep(NA, length(unique(meth_predict_hc_dunedinpace_three_healthy$RID)))
x <- 1
for (i in unique(meth_predict_hc_dunedinpace_three_healthy$RID)){
  temp_df <- meth_predict_hc_dunedinpace_three_healthy[meth_predict_hc_dunedinpace_three_healthy$RID == i,]
  hc_obs_length_three[x] <- max(temp_df$scanage) - min(temp_df$scanage)
  baseline_age_hc_three[x] <- min(temp_df$scanage)
  x <- x + 1
}


random_slopes_hc_three_df <- data.frame(RID = unique(meth_predict_hc_dunedinpace_three_healthy$RID), 
                                        sex = meth_predict_hc_dunedinpace_three_healthy$sex[!duplicated(meth_predict_hc_dunedinpace_three_healthy$RID)],
                                        mean_age = meth_predict_hc_dunedinpace_three_healthy$scanage_scale[!duplicated(meth_predict_hc_dunedinpace_three_healthy$RID)],
                                        baseline_age_hc = baseline_age_hc_three,
                                        obs_length = hc_obs_length_three,
                                        growth_curve = hc_randomslopes_three, 
                                        growth_curve_icv = hc_randomslopes_icv_three, 
                                        baseline_acc = unique(meth_predict_hc_dunedinpace_three_healthy$baseline_acc),
                                        baseline_pace = unique(meth_predict_hc_dunedinpace_three_healthy$baseline_pace))


# this is now the model we want 8/17
stable_hc_dunedinpace <- summary(lm(scale(growth_curve_icv)~scale(baseline_acc)+sex+scale(obs_length)+scale(baseline_age_hc), data = random_slopes_hc_three_df))


# WMH

# restrict to three timepoitns
wmh_three_index <- meth_predict_wmh_dunedinpace$RID %in% data.frame(table(meth_predict_wmh_dunedinpace$RID)[table(meth_predict_wmh_dunedinpace$RID) >= 3])$Var1
meth_predict_wmh_dunedinpace_three <- meth_predict_wmh_dunedinpace[wmh_three_index,]

# restrict by diagnosis - step moved up on 8/17
meth_predict_wmh_dunedinpace_three_dx <- left_join(meth_predict_wmh_dunedinpace_three, long_dx, by= 'RID')
meth_predict_wmh_dunedinpace_three_healthy <- meth_predict_wmh_dunedinpace_three_dx[meth_predict_wmh_dunedinpace_three_dx$group == 'sCN',]


# three
# meth_predict_wmh_dunedinpace_three$var_resid <- resid(lm(var_scale~icv_scale, data = meth_predict_wmh_dunedinpace_three))
# # lme() doesn't like this formula for some reason. lme() and lmer() seem to give comparable results though, so using lmer() for this one
wmh_randomslopes_three <- coef(lmer(var_scale~1 + (scanage_scale|RID), data = meth_predict_wmh_dunedinpace_three_healthy))$RID[,1]


wmh_obs_length_three <- rep(NA, length(unique(meth_predict_wmh_dunedinpace_three_healthy$RID)))
baseline_age_wmh_three <- rep(NA, length(unique(meth_predict_wmh_dunedinpace_three_healthy$RID)))
x <- 1
for (i in unique(meth_predict_wmh_dunedinpace_three_healthy$RID)){
  temp_df <- meth_predict_wmh_dunedinpace_three_healthy[meth_predict_wmh_dunedinpace_three_healthy$RID == i,]
  wmh_obs_length_three[x] <- max(temp_df$scanage) - min(temp_df$scanage)
  baseline_age_wmh_three[x] <- min(temp_df$scanage)
  x <- x + 1
}


random_slopes_wmh_three_df <- data.frame(RID = unique(meth_predict_wmh_dunedinpace_three_healthy$RID), 
                                         sex = meth_predict_wmh_dunedinpace_three_healthy$sex[!duplicated(meth_predict_wmh_dunedinpace_three_healthy$RID)],
                                         mean_age = meth_predict_wmh_dunedinpace_three_healthy$scanage_scale[!duplicated(meth_predict_wmh_dunedinpace_three_healthy$RID)],
                                         baseline_age_wmh = baseline_age_wmh_three,
                                         obs_length = wmh_obs_length_three,
                                         growth_curve = wmh_randomslopes_three, 
                                         baseline_acc = unique(meth_predict_wmh_dunedinpace_three_healthy$baseline_acc),
                                         baseline_pace = unique(meth_predict_wmh_dunedinpace_three_healthy$baseline_pace))

save(random_slopes_wmh_three_df, file = '/Volumes/data_commons-hariri-long/Scripts/random_slopes_wmh_three_df_8_17_23.Rdata')


# this is now the model we want 8/17
stable_wmh_dunedinpace <- summary(lm(scale(growth_curve)~scale(baseline_acc)+sex+scale(obs_length)+scale(baseline_age_wmh), data = random_slopes_wmh_three_df))


# SA

# restrict to three timepoints
sa_three_index <- meth_predict_sa_dunedinpace$RID %in% data.frame(table(meth_predict_sa_dunedinpace$RID)[table(meth_predict_sa_dunedinpace$RID) >= 3])$Var1
meth_predict_sa_dunedinpace_three <- meth_predict_sa_dunedinpace[sa_three_index,]

# restrict by diagnosis - step moved up on 8/17
meth_predict_sa_dunedinpace_three_dx <- left_join(meth_predict_sa_dunedinpace_three, long_dx, by= 'RID')
meth_predict_sa_dunedinpace_three_healthy <- meth_predict_sa_dunedinpace_three_dx[meth_predict_sa_dunedinpace_three_dx$group == 'sCN',]


# three
# meth_predict_sa_dunedinpace_three$var_resid <- resid(lm(var_scale~icv_scale, data = meth_predict_sa_dunedinpace_three))
# # lme() doesn't like this formula for some reason. lme() and lmer() seem to give comparable results though, so using lmer() for this one
sa_randomslopes_three <- coef(lmer(var_scale~1 + (scanage_scale|RID), data = meth_predict_sa_dunedinpace_three_healthy))$RID[,1]


sa_obs_length_three <- rep(NA, length(unique(meth_predict_sa_dunedinpace_three_healthy$RID)))
baseline_age_sa_three <- rep(NA, length(unique(meth_predict_sa_dunedinpace_three_healthy$RID)))
x <- 1
for (i in unique(meth_predict_sa_dunedinpace_three_healthy$RID)){
  temp_df <- meth_predict_sa_dunedinpace_three_healthy[meth_predict_sa_dunedinpace_three_healthy$RID == i,]
  sa_obs_length_three[x] <- max(temp_df$scanage) - min(temp_df$scanage)
  baseline_age_sa_three[x] <- min(temp_df$scanage)
  x <- x + 1
}


random_slopes_sa_three_df <- data.frame(RID = unique(meth_predict_sa_dunedinpace_three_healthy$RID), 
                                        sex = meth_predict_sa_dunedinpace_three_healthy$sex[!duplicated(meth_predict_sa_dunedinpace_three_healthy$RID)],
                                        mean_age = meth_predict_sa_dunedinpace_three_healthy$scanage_scale[!duplicated(meth_predict_sa_dunedinpace_three_healthy$RID)],
                                        baseline_age_sa = baseline_age_sa_three,
                                        obs_length = sa_obs_length_three,
                                        growth_curve = sa_randomslopes_three, 
                                        baseline_acc = unique(meth_predict_sa_dunedinpace_three_healthy$baseline_acc),
                                        baseline_pace = unique(meth_predict_sa_dunedinpace_three_healthy$baseline_pace))


# this is now the model we want 8/17
stable_sa_dunedinpace <- summary(lm(scale(growth_curve)~scale(baseline_acc)+sex+scale(obs_length)+scale(baseline_age_sa), data = random_slopes_sa_three_df))

save(random_slopes_sa_three_df, file = '/Volumes/data_commons-hariri-long/Scripts/random_slopes_sa_three_df_8_17_23.Rdata')

# CT

# restrict to three timepoitns
ct_three_index <- meth_predict_ct_dunedinpace$RID %in% data.frame(table(meth_predict_ct_dunedinpace$RID)[table(meth_predict_ct_dunedinpace$RID) >= 3])$Var1
meth_predict_ct_dunedinpace_three <- meth_predict_ct_dunedinpace[ct_three_index,]

# restrict by diagnosis - step moved up on 8/17
meth_predict_ct_dunedinpace_three_dx <- left_join(meth_predict_ct_dunedinpace_three, long_dx, by= 'RID')
meth_predict_ct_dunedinpace_three_healthy <- meth_predict_ct_dunedinpace_three_dx[meth_predict_ct_dunedinpace_three_dx$group == 'sCN',]


# three
# meth_predict_ct_dunedinpace_three$var_resid <- resid(lm(var_scale~icv_scale, data = meth_predict_ct_dunedinpace_three))
# # lme() doesn't like this formula for some reason. lme() and lmer() seem to give comparable results though, so using lmer() for this one
ct_randomslopes_three <- coef(lmer(var_scale~1 + (scanage_scale|RID), data = meth_predict_ct_dunedinpace_three_healthy))$RID[,1]


ct_obs_length_three <- rep(NA, length(unique(meth_predict_ct_dunedinpace_three_healthy$RID)))
baseline_age_ct_three <- rep(NA, length(unique(meth_predict_ct_dunedinpace_three_healthy$RID)))
x <- 1
for (i in unique(meth_predict_ct_dunedinpace_three_healthy$RID)){
  temp_df <- meth_predict_ct_dunedinpace_three_healthy[meth_predict_ct_dunedinpace_three_healthy$RID == i,]
  ct_obs_length_three[x] <- max(temp_df$scanage) - min(temp_df$scanage)
  baseline_age_ct_three[x] <- min(temp_df$scanage)
  x <- x + 1
}


random_slopes_ct_three_df <- data.frame(RID = unique(meth_predict_ct_dunedinpace_three_healthy$RID), 
                                        sex = meth_predict_ct_dunedinpace_three_healthy$sex[!duplicated(meth_predict_ct_dunedinpace_three_healthy$RID)],
                                        mean_age = meth_predict_ct_dunedinpace_three_healthy$scanage_scale[!duplicated(meth_predict_ct_dunedinpace_three_healthy$RID)],
                                        baseline_age_ct = baseline_age_ct_three,
                                        obs_length = ct_obs_length_three,
                                        growth_curve = ct_randomslopes_three, 
                                        baseline_acc = unique(meth_predict_ct_dunedinpace_three_healthy$baseline_acc),
                                        baseline_pace = unique(meth_predict_ct_dunedinpace_three_healthy$baseline_pace))


# this is now the model we want 8/17
stable_ct_dunedinpace <- summary(lm(scale(growth_curve)~scale(baseline_acc)+sex+scale(obs_length)+scale(baseline_age_ct), data = random_slopes_ct_three_df))


# WMHypo

# restrict to three timepoitns
wmhypo_three_index <- meth_predict_wmhypo_dunedinpace$RID %in% data.frame(table(meth_predict_wmhypo_dunedinpace$RID)[table(meth_predict_wmhypo_dunedinpace$RID) >= 3])$Var1
meth_predict_wmhypo_dunedinpace_three <- meth_predict_wmhypo_dunedinpace[wmhypo_three_index,]

# restrict by diagnosis - step moved up on 8/17
meth_predict_wmhypo_dunedinpace_three_dx <- left_join(meth_predict_wmhypo_dunedinpace_three, long_dx, by= 'RID')
meth_predict_wmhypo_dunedinpace_three_healthy <- meth_predict_wmhypo_dunedinpace_three_dx[meth_predict_wmhypo_dunedinpace_three_dx$group == 'sCN',]


# three
# meth_predict_wmhypo_dunedinpace_three$var_resid <- resid(lm(var_scale~icv_scale, data = meth_predict_wmhypo_dunedinpace_three))
# # lme() doesn't like this formula for some reason. lme() and lmer() seem to give comparable results though, so using lmer() for this one
wmhypo_randomslopes_three <- coef(lmer(var_scale~1 + (scanage_scale|RID), data = meth_predict_wmhypo_dunedinpace_three_healthy))$RID[,1]


wmhypo_obs_length_three <- rep(NA, length(unique(meth_predict_wmhypo_dunedinpace_three_healthy$RID)))
baseline_age_wmhypo_three <- rep(NA, length(unique(meth_predict_wmhypo_dunedinpace_three_healthy$RID)))
x <- 1
for (i in unique(meth_predict_wmhypo_dunedinpace_three_healthy$RID)){
  temp_df <- meth_predict_wmhypo_dunedinpace_three_healthy[meth_predict_wmhypo_dunedinpace_three_healthy$RID == i,]
  wmhypo_obs_length_three[x] <- max(temp_df$scanage) - min(temp_df$scanage)
  baseline_age_wmhypo_three[x] <- min(temp_df$scanage)
  x <- x + 1
}


random_slopes_wmhypo_three_df <- data.frame(RID = unique(meth_predict_wmhypo_dunedinpace_three_healthy$RID), 
                                        sex = meth_predict_wmhypo_dunedinpace_three_healthy$sex[!duplicated(meth_predict_wmhypo_dunedinpace_three_healthy$RID)],
                                        mean_age = meth_predict_wmhypo_dunedinpace_three_healthy$scanage_scale[!duplicated(meth_predict_wmhypo_dunedinpace_three_healthy$RID)],
                                        baseline_age_wmhypo = baseline_age_wmhypo_three,
                                        obs_length = wmhypo_obs_length_three,
                                        growth_curve = wmhypo_randomslopes_three, 
                                        baseline_acc = unique(meth_predict_wmhypo_dunedinpace_three_healthy$baseline_acc),
                                        baseline_pace = unique(meth_predict_wmhypo_dunedinpace_three_healthy$baseline_pace))


# this is now the model we want 8/17
stable_wmhypo_dunedinpace <- summary(lm(scale(growth_curve)~scale(baseline_acc)+sex+scale(obs_length)+scale(baseline_age_wmhypo), data = random_slopes_wmhypo_three_df))




# save out so i dont need to do this again
# save(random_slopes_tbv_three_df, file = '/Users/ew198/Documents/methylation/data/adni/longitudinal/random_slopes_tbv_three_df.Rdata')
# save(random_slopes_hc_three_df, file = '/Users/ew198/Documents/methylation/data/adni/longitudinal/random_slopes_hc_three_df.Rdata')
# save(random_slopes_wmh_three_df, file = '/Users/ew198/Documents/methylation/data/adni/longitudinal/random_slopes_wmh_three_df.Rdata')
# save(random_slopes_ct_three_df, file = '/Users/ew198/Documents/methylation/data/adni/longitudinal/random_slopes_ct_three_df.Rdata')
# save(random_slopes_sa_three_df, file = '/Users/ew198/Documents/methylation/data/adni/longitudinal/random_slopes_sa_three_df.Rdata')
# save(random_slopes_wmhypo_three_df, file = '/Users/ew198/Documents/methylation/data/adni/longitudinal/random_slopes_wmhypo_three_df.Rdata')


# TBV

tbv_rs_plot <- ggplot(meth_predict_tbv_dunedinpace_three_healthy, 
                      aes(x=scanage_scale, y=var_scale)) +
  geom_point(color = '#619CFF', size = .75) +
  geom_line(aes(y=predict(lme(var_scale~1, random=~scanage_scale|RID, data = meth_predict_tbv_dunedinpace_three_healthy)), group=RID), 
            color = '#619CFF') +
  theme_classic() +
  labs(x = 'Age (standardized)', 
       y = 'TBV (standardized)',
       title = 'TBV')+
  theme(axis.text.x = element_text(vjust = 0.5, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15, vjust = 1),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 17), 
        legend.position = 'none')


# HC

hc_rs_plot <- ggplot(meth_predict_hc_dunedinpace_three_healthy, 
                     aes(x=scanage_scale, y=var_scale)) +
  geom_point(color = '#619CFF', size = .75) +
  geom_line(aes(y=predict(lme(var_scale~1, random=~scanage_scale|RID, data = meth_predict_hc_dunedinpace_three_healthy)), group=RID), 
            color = '#619CFF') +
  theme_classic() +
  labs(x = 'Age (standardized)', 
       y = 'HC (standardized)',
       title = 'HC')+
  theme(axis.text.x = element_text(vjust = 0.5, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15, vjust = 1),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 17), 
        legend.position = 'none')

# WMHyper

wmh_rs_plot <- ggplot(meth_predict_wmh_dunedinpace_three_healthy, 
                      aes(x=scanage_scale, y=var_scale)) +
  geom_point(color = '#619CFF', size = .75) +
  geom_line(aes(y=predict(lme(var_scale~1, random=~scanage_scale|RID, data = meth_predict_wmh_dunedinpace_three_healthy)), group=RID), 
            color = '#619CFF') +
  theme_classic() +
  labs(x = 'Age (standardized)', 
       y = 'WMHyper (standardized)',
       title = 'WMHyper')+
  theme(axis.text.x = element_text(vjust = 0.5, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15, vjust = 1),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 17), 
        legend.position = 'none')


# CT

ct_rs_plot <- ggplot(meth_predict_ct_dunedinpace_three_healthy, 
                     aes(x=scanage_scale, y=var_scale)) +
  geom_point(color = '#619CFF', size = .75) +
  geom_line(aes(y=predict(lme(var_scale~1, random=~scanage_scale|RID, data = meth_predict_ct_dunedinpace_three_healthy)), group=RID), 
            color = '#619CFF') +
  theme_classic() +
  labs(x = 'Age (standardized)', 
       y = 'CT (standardized)',
       title = 'CT')+
  theme(axis.text.x = element_text(vjust = 0.5, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15, vjust = 1),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 17), 
        legend.position = 'none')

# SA

sa_rs_plot <- ggplot(meth_predict_sa_dunedinpace_three_healthy, 
                     aes(x=scanage_scale, y=var_scale)) +
  geom_point(color = '#619CFF', size = .75) +
  geom_line(aes(y=predict(lme(var_scale~1, random=~scanage_scale|RID, data = meth_predict_sa_dunedinpace_three_healthy)), group=RID), 
            color = '#619CFF') +
  theme_classic() +
  labs(x = 'Age (standardized)', 
       y = 'SA (standardized)',
       title = 'SA')+
  theme(axis.text.x = element_text(vjust = 0.5, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15, vjust = 1),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 17), 
        legend.position = 'none')

# WMHypo

wmhypo_rs_plot <- ggplot(meth_predict_wmhypo_dunedinpace_three_healthy, 
                     aes(x=scanage_scale, y=var_scale)) +
  geom_point(color = '#619CFF', size = .75) +
  geom_line(aes(y=predict(lme(var_scale~1, random=~scanage_scale|RID, data = meth_predict_wmhypo_dunedinpace_three_healthy)), group=RID), 
            color = '#619CFF') +
  theme_classic() +
  labs(x = 'Age (standardized)', 
       y = 'WMHypo volume (standardized)',
       title = 'WMHypo')+
  theme(axis.text.x = element_text(vjust = 0.5, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15, vjust = 1),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 17), 
        legend.position = 'none')



# random slopes historgrams

tbv_gc_hist <- ggplot(data = random_slopes_tbv_three_df, aes(x = growth_curve_icv)) + 
  geom_histogram(fill='#619CFF') +
  ylim(0,33)+
  labs(title = "TBV",
       x = "TBV growth curve slope",
       y = 'Frequency') +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15, vjust = 1),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.position = 'none')

hc_gc_hist <- ggplot(data = random_slopes_hc_three_df, aes(x = growth_curve_icv)) + 
  geom_histogram(fill='#619CFF') +
  ylim(0,33)+
  labs(title = "HC",
       x = "HC growth curve slope",
       y = 'Frequency') +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15, vjust = 1),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.position = 'none')

wmh_gc_hist <- ggplot(data = random_slopes_wmh_three_df, aes(x = growth_curve)) + 
  geom_histogram(fill='#619CFF') +
  ylim(0,33)+
  labs(title = "WMHyper",
       x = "WMHyper growth curve slope",
       y = 'Frequency') +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15, vjust = 1),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.position = 'none')

ct_gc_hist <- ggplot(data = random_slopes_ct_three_df, aes(x = growth_curve)) + 
  geom_histogram(fill='#619CFF') +
  ylim(0,33)+
  labs(title = "CT",
       x = "CT growth curve slope",
       y = 'Frequency') +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15, vjust = 1),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.position = 'none')

sa_gc_hist <- ggplot(data = random_slopes_sa_three_df, aes(x = growth_curve)) + 
  geom_histogram(fill='#619CFF') +
  ylim(0,33)+
  labs(title = "SA",
       x = "SA growth curve slope",
       y = 'Frequency') +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15, vjust = 1),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.position = 'none')

wmhypo_gc_hist <- ggplot(data = random_slopes_wmhypo_three_df, aes(x = growth_curve)) + 
  geom_histogram(fill='#619CFF') +
  ylim(0,33)+
  labs(title = "WMHypo",
       x = "WMHypo growth curve slope",
       y = 'Frequency') +
  theme_classic()+
  theme(axis.text.x = element_text(vjust = 0.5, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15, vjust = 1),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.position = 'none')

# 1100 x 600 - for scatterplots
grid.arrange(tbv_rs_plot, hc_rs_plot, wmhypo_rs_plot, ct_rs_plot, sa_rs_plot, wmh_rs_plot, nrow = 2)


# 1100 x 400 for histograms
grid.arrange(tbv_gc_hist, hc_gc_hist, wmh_gc_hist, ct_gc_hist, sa_gc_hist, wmhypo_gc_hist, nrow = 2)

# sample size descriptive table


three_sample_desc <- data.frame(measure = c('TBV', 'HC', 'WMHypo', 'CT', 'SA', 'WMHyper'),
                                n = c(nrow(random_slopes_tbv_three_df), nrow(random_slopes_hc_three_df), nrow(random_slopes_wmhypo_three_df), nrow(random_slopes_ct_three_df), 
                                      nrow(random_slopes_sa_three_df), nrow(random_slopes_wmh_three_df)),
                                mean_age_baseline = c(mean(random_slopes_tbv_three_df$baseline_age_tbv), mean(random_slopes_hc_three_df$baseline_age), mean(random_slopes_wmhypo_three_df$baseline_age), 
                                                      mean(random_slopes_ct_three_df$baseline_age), mean(random_slopes_sa_three_df$baseline_age), mean(random_slopes_wmh_three_df$baseline_age))/12,
                                
                                mean_obs_length = c(mean(random_slopes_tbv_three_df$obs_length), mean(random_slopes_hc_three_df$obs_length), mean(random_slopes_wmhypo_three_df$obs_length), 
                                                    mean(random_slopes_ct_three_df$obs_length), mean(random_slopes_sa_three_df$obs_length), mean(random_slopes_wmh_three_df$obs_length))/12,
                                
                                malesex = c(sum(random_slopes_tbv_three_df$sex == 'Male')/nrow(random_slopes_tbv_three_df),
                                            sum(random_slopes_hc_three_df$sex == 'Male')/nrow(random_slopes_hc_three_df),
                                            sum(random_slopes_wmhypo_three_df$sex == 'Male')/nrow(random_slopes_wmhypo_three_df),
                                            sum(random_slopes_ct_three_df$sex == 'Male')/nrow(random_slopes_ct_three_df),
                                            sum(random_slopes_sa_three_df$sex == 'Male')/nrow(random_slopes_sa_three_df),
                                            sum(random_slopes_wmh_three_df$sex == 'Male')/nrow(random_slopes_wmh_three_df))
)



write.csv(three_sample_desc, file = '/Users/ew198/Documents/methylation/results/three_sample_desc_11_6_23.csv')


# long effects forest plot


long_effects_stable_foresttable <- data.frame(imaging = c('TBV', 'HC', 'WMHypo', 'CT', 'SA', 'WMHyper'),
                                          beta = c(stable_tbv_dunedinpace$coefficients[2,1], stable_hc_dunedinpace$coefficients[2,1], stable_wmhypo_dunedinpace$coefficients[2,1], 
                                                   stable_ct_dunedinpace$coefficients[2,1], stable_sa_dunedinpace$coefficients[2,1], stable_wmh_dunedinpace$coefficients[2,1]),
                                          
                                          
                                          std_error = c(stable_tbv_dunedinpace$coefficients[2,2], stable_hc_dunedinpace$coefficients[2,2], stable_wmhypo_dunedinpace$coefficients[2,2], 
                                                        stable_ct_dunedinpace$coefficients[2,2], stable_sa_dunedinpace$coefficients[2,2], stable_wmh_dunedinpace$coefficients[2,1])
)

mri_levels_all <- c('TBV', 'HC', 'WMHypo', 'CT', 'SA', 'WMHyper')

long_healthy_pointrange <- ggplot(long_effects_stable_foresttable, aes(x=factor(imaging,levels = rev(mri_levels_all)), y=beta)) + 
  geom_pointrange(aes(ymin=beta-(1.96*std_error), ymax=beta+(1.96*std_error)), color='#619CFF', position=position_dodge(width=0.6)) +
  geom_hline(yintercept = 0)+
  labs(x = 'MRI phenotype') +
  coord_flip() +
  theme_classic() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        title = element_text(size = 20),
        legend.text = element_text(size = 20))



################################################
########### SENSITIVITY TEST TABLE #############
################################################


sensitivity_test_table <- data.frame(imaging = c( rep('TBV', 3), rep('HC', 3), rep('WMHypo', 3), rep('CT', 3), rep('SA', 3)),
                                     test = rep(c('Dunedin', 'FHS-OC', 'ADNI'), 5),
                                     
                                     # MAIN ANALYSES
                                     main_beta = c(tbv.dp.res[2,1], -0.03752815, tbvrel.dp.plm.res[2,1],
                                                   hc.dp.res[2,1], -0.071963606, hcrel.dp.plm.res[2,1],
                                                   wmhypo.dp.res[2,1], 0.093908079, wmhypo.dp.plm.res[2,1],
                                                   ct.dp.res[2,1], -0.092464541, ct.dp.plm.res[2,1],
                                                   sa.dp.res[2,1], -0.009306709, sa.dp.plm.res[2,1]),
                                     
                                     
                                     main_p = c(tbv.dp.res[2,4], 0.02917189, tbvrel.dp.plm.res[2,4],
                                                hc.dp.res[2,4], 0.017619982, hcrel.dp.plm.res[2,4],
                                                wmhypo.dp.res[2,4], 0.008314051, wmhypo.dp.plm.res[2,4],
                                                ct.dp.res[2,4], 0.004597945, ct.dp.plm.res[2,4],
                                                sa.dp.res[2,4], 0.754979375, sa.dp.plm.res[2,4]),
                                     
                                     main_lowbound = c(tbv.dp.res[2,1]-1.96*tbv.dp.res[2,2], -0.071241917, tbvrel.dp.plm.res[2,1]-1.96*tbvrel.dp.plm.res[2,2],
                                                       hc.dp.res[2,1]-1.96*hc.dp.res[2,2], -0.131358964, hcrel.dp.plm.res[2,1]-1.96*hcrel.dp.plm.res[2,2],
                                                       wmhypo.dp.res[2,1]-1.96*wmhypo.dp.res[2,2], 0.024225273, wmhypo.dp.plm.res[2,1]-1.96*wmhypo.dp.plm.res[2,2],
                                                       ct.dp.res[2,1]-1.96*ct.dp.res[2,2], -0.156338781, ct.dp.plm.res[2,1]-1.96*ct.dp.plm.res[2,2],
                                                       sa.dp.res[2,1]-1.96*sa.dp.res[2,2], -0.067816925, sa.dp.plm.res[2,1]-1.96*sa.dp.plm.res[2,2]),
                                     
                                     main_highbound = c(tbv.dp.res[2,1]+1.96*tbv.dp.res[2,2], -0.003814383, tbvrel.dp.plm.res[2,1]+1.96*tbvrel.dp.plm.res[2,2],
                                                        hc.dp.res[2,1]+1.96*hc.dp.res[2,2], -0.012568248, hcrel.dp.plm.res[2,1]+1.96*hcrel.dp.plm.res[2,2],
                                                        wmhypo.dp.res[2,1]+1.96*wmhypo.dp.res[2,2], 0.163590885, wmhypo.dp.plm.res[2,1]+1.96*wmhypo.dp.plm.res[2,2],
                                                        ct.dp.res[2,1]+1.96*ct.dp.res[2,2], -0.028590301, ct.dp.plm.res[2,1]+1.96*ct.dp.plm.res[2,2],
                                                        sa.dp.res[2,1]+1.96*sa.dp.res[2,2], 0.049203508, sa.dp.plm.res[2,1]+1.96*sa.dp.plm.res[2,2]),
                                     
                                     
                                     # WBC controls
                                     wbc_beta = c(tbvrel.dp.res.wca[2,1], -0.073081905, tbvrel.dp.plm.res.wca[2,1],
                                                  hcrel.dp.res.wca[2,1], -0.087006019, hcrel.dp.plm.res.wca[2,1],
                                                  wmhypo.dp.res.wca[2,1], 0.093584634, wmhypo.dp.plm.res.wca[2,1],
                                                  ct.dp.res.wca[2,1], -0.091602146, ct.dp.plm.res.wca[2,1],
                                                  sa.dp.res.wca[2,1], -0.01457126, sa.dp.plm.res.wca[2,1]),
                                     
                                     
                                     wbc_p = c(tbvrel.dp.res.wca[2,4], 0.019958554, tbvrel.dp.plm.res.wca[2,4],
                                               hcrel.dp.res.wca[2,4], 0.009969498, hcrel.dp.plm.res.wca[2,4],
                                               wmhypo.dp.res.wca[2,4], 0.007192243, wmhypo.dp.plm.res.wca[2,4],
                                               ct.dp.res.wca[2,4], 0.005818308, ct.dp.plm.res.wca[2,4],
                                               sa.dp.res.wca[2,4], 0.662513509, sa.dp.plm.res.wca[2,4]),
                                     
                                     wbc_lowbound = c(tbvrel.dp.res.wca[2,1]-1.96*tbvrel.dp.res.wca[2,2], -0.134605951, tbvrel.dp.plm.res.wca[2,1]-1.96*tbvrel.dp.plm.res.wca[2,2],
                                                      hcrel.dp.res.wca[2,1]-1.96*hcrel.dp.res.wca[2,2], -0.15312997, hcrel.dp.plm.res.wca[2,1]-1.96*hcrel.dp.plm.res.wca[2,2],
                                                      wmhypo.dp.res.wca[2,1]-1.96*wmhypo.dp.res.wca[2,2], 0.025407461, wmhypo.dp.plm.res.wca[2,1]-1.96*wmhypo.dp.plm.res.wca[2,2],
                                                      ct.dp.res.wca[2,1]-1.96*ct.dp.res.wca[2,2], -0.156633271, ct.dp.plm.res.wca[2,1]-1.96*ct.dp.plm.res.wca[2,2],
                                                      sa.dp.res.wca[2,1]-1.96*sa.dp.res.wca[2,2], -0.08007452, sa.dp.plm.res.wca[2,1]-1.96*sa.dp.plm.res.wca[2,2]),
                                     
                                     wbc_highbound = c(tbvrel.dp.res.wca[2,1]+1.96*tbvrel.dp.res.wca[2,2], -0.01155786, tbvrel.dp.plm.res.wca[2,1]+1.96*tbvrel.dp.plm.res.wca[2,2],
                                                       hcrel.dp.res.wca[2,1]+1.96*hcrel.dp.res.wca[2,2], -0.020882069, hcrel.dp.plm.res.wca[2,1]+1.96*hcrel.dp.plm.res.wca[2,2],
                                                       wmhypo.dp.res.wca[2,1]+1.96*wmhypo.dp.res.wca[2,2], 0.161761807, wmhypo.dp.plm.res.wca[2,1]+1.96*wmhypo.dp.plm.res.wca[2,2],
                                                       ct.dp.res.wca[2,1]+1.96*ct.dp.res.wca[2,2], -0.026571022, ct.dp.plm.res.wca[2,1]+1.96*ct.dp.plm.res.wca[2,2],
                                                       sa.dp.res.wca[2,1]+1.96*sa.dp.res.wca[2,2], 0.050932, sa.dp.plm.res.wca[2,1]+1.96*sa.dp.plm.res.wca[2,2]),
                                     
                                     
                                     # APOE E4 controls
                                     apoe_beta = c(tbvrel.dp.res.apoe[2,1], -0.083887164, tbvrel.dp.plm.res.apoe[2,1],
                                                   hcrel.dp.res.apoe[2,1], -0.100351872, hcrel.dp.plm.res.apoe[2,1],
                                                   wmhypo.dp.res.apoe[2,1], 0.075762139, wmhypo.dp.plm.res.apoe[2,1],
                                                   ct.dp.res.apoe[2,1], -0.091728743, ct.dp.plm.res.apoe[2,1],
                                                   sa.dp.res.apoe[2,1], -0.023891654, sa.dp.plm.res.apoe[2,1]),
                                     
                                     
                                     apoe_p = c(tbvrel.dp.res.apoe[2,4], 0.003700659, tbvrel.dp.plm.res.apoe[2,4],
                                                hcrel.dp.res.apoe[2,4], 0.001630896, hcrel.dp.plm.res.apoe[2,4],
                                                wmhypo.dp.res.apoe[2,4], 0.017662179, wmhypo.dp.plm.res.apoe[2,4],
                                                ct.dp.res.apoe[2,4], 0.002922081, ct.dp.plm.res.apoe[2,4],
                                                sa.dp.res.apoe[2,4], 0.432701521, sa.dp.plm.res.apoe[2,4]),
                                     
                                     apoe_lowbound = c(tbvrel.dp.res.apoe[2,1]-1.96*tbvrel.dp.res.apoe[2,2], -0.140451249, tbvrel.dp.plm.res.apoe[2,1]-1.96*tbvrel.dp.plm.res.apoe[2,2],
                                                       hcrel.dp.res.apoe[2,1]-1.96*hcrel.dp.res.apoe[2,2], -0.162668425, hcrel.dp.plm.res.apoe[2,1]-1.96*hcrel.dp.plm.res.apoe[2,2],
                                                       wmhypo.dp.res.apoe[2,1]-1.96*wmhypo.dp.res.apoe[2,2], 0.013211496, wmhypo.dp.plm.res.apoe[2,1]-1.96*wmhypo.dp.plm.res.apoe[2,2],
                                                       ct.dp.res.apoe[2,1]-1.96*ct.dp.res.apoe[2,2], -0.15205071, ct.dp.plm.res.apoe[2,1]-1.96*ct.dp.plm.res.apoe[2,2],
                                                       sa.dp.res.apoe[2,1]-1.96*sa.dp.res.apoe[2,2], -0.08363512, sa.dp.plm.res.apoe[2,1]-1.96*sa.dp.plm.res.apoe[2,2]),
                                     
                                     apoe_highbound = c(tbvrel.dp.res.apoe[2,1]+1.96*tbvrel.dp.res.apoe[2,2], -0.027323078, tbvrel.dp.plm.res.apoe[2,1]+1.96*tbvrel.dp.plm.res.apoe[2,2],
                                                        hcrel.dp.res.apoe[2,1]+1.96*hcrel.dp.res.apoe[2,2], -0.038035319, hcrel.dp.plm.res.apoe[2,1]+1.96*hcrel.dp.plm.res.apoe[2,2],
                                                        wmhypo.dp.res.apoe[2,1]+1.96*wmhypo.dp.res.apoe[2,2], 0.138312781, wmhypo.dp.plm.res.apoe[2,1]+1.96*wmhypo.dp.plm.res.apoe[2,2],
                                                        ct.dp.res.apoe[2,1]+1.96*ct.dp.res.apoe[2,2], -0.031406776, ct.dp.plm.res.apoe[2,1]+1.96*ct.dp.plm.res.apoe[2,2],
                                                        sa.dp.res.apoe[2,1]+1.96*sa.dp.res.apoe[2,2], 0.035851811, sa.dp.plm.res.apoe[2,1]+1.96*sa.dp.plm.res.apoe[2,2])
)


write.csv(sensitivity_test_table, file = '/Users/ew198/Documents/methylation/drafts/sensitivity_test_table_8_17_23.csv')



### longitudinal race table


tbv_long_race <- adni_race[adni_race$RID %in% random_slopes_tbv_three_df$RID,]
hc_long_race <- adni_race[adni_race$RID %in% random_slopes_hc_three_df$RID,]
wmh_long_race <- adni_race[adni_race$RID %in% random_slopes_wmh_three_df$RID,]
fs_long_race <- adni_race[adni_race$RID %in% random_slopes_ct_three_df$RID,]

table(wmh_long_race$PTRACCAT)
table(wmh_long_race$PTETHCAT)



