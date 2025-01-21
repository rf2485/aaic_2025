source("1_data_preparation.R")
library(gtsummary)

failed_qc <- c('sub-CC510255', #SCD abnormality in left temporal pole
               'sub-CC510438', #CTL abnormality in left frontal lobe
               'sub-CC620821', #SCD segmentation errors from large ventricles
               'sub-CC621011', #CTL segmentation errors from large ventricles
               'sub-CC621080', #SCD segmentation errors
               'sub-CC710551', #CTL motion artifacts in DWI
               'sub-CC711027', #SCD severe motion artifacts in T1
               'sub-CC721434' #CTL segmentation errors from large ventricles
)

amy_filter <- c("SCD", "Left.Amygdala", "Right.Amygdala")

#define remove_outliers function
remove_outliers <- function(x, na.rm = TRUE)
{
  ## Find 25% and 75% Quantiles using inbuild function
  quant <- quantile(x, probs=c(.25, .75), na.rm = na.rm)

  ## Find Interquantile range and multiply it by 1.5
  ## to derive factor for range calculation
  H <- 1.5 * IQR(x, na.rm = na.rm)

  y <- x

  ## fill the outlier elements with NA
  y[x < (quant[1] - H)] <- NA
  y[x > (quant[2] + H)] <- NA

  y
}

#extract SCD status from dwi_over_55 table
scd_status <- dwi_over_55 %>% select(participant_id, SCD) %>%
  filter(!participant_id %in% failed_qc)
scd_status$SCD <- factor(scd_status$SCD,
                         levels = c(1,0),
                         labels = c('SCD', 'Control'))

#import aseg stats table (subcortical volumes)
volumes <- read_tsv("freesurfer/asegtable.tsv") %>%
  rename(participant_id=`Measure:volume`) %>%
  mutate(across(c(3:ncol(.)), .fns = ~.*1000/EstimatedTotalIntraCranialVol)) %>% #normalize by intracranial volume
  # mutate(across(where(is.double), remove_outliers)) %>% #remove outliers (change to NA)
  left_join(scd_status, .)
names(volumes) <- make.names(names(volumes))

#read in diffusion and MTR tables
aparc2meas_files <- list.files(path = "freesurfer", 
                               pattern = "aparc.*aseg.*\\.*tsv", 
                               full.names = T)
for (i in 1:length(aparc2meas_files)) {
  assign(gsub(".tsv", "", 
              gsub("freesurfer.aparc.", "", make.names(aparc2meas_files[i]))), 
         read.delim(aparc2meas_files[i]))
}

fit_FWF <- aseg2fit_FWF %>%
  rename(participant_id = Measure.mean) %>%
  # mutate(across(where(is.double), remove_outliers)) %>% #remove outliers (change to NA)
  left_join(scd_status, .)
  

fit_NDI <- aseg2fit_NDI %>%
  rename(participant_id = Measure.mean) %>%
  # mutate(across(where(is.double), remove_outliers)) %>% #remove outliers (change to NA)
  left_join(scd_status, .)

fit_ODI <- aseg2fit_ODI %>%
  rename(participant_id = Measure.mean) %>%
  # mutate(across(where(is.double), remove_outliers)) %>% #remove outliers (change to NA)
  left_join(scd_status, .)

dti_fa <- aseg2dti_fa %>%
  rename(participant_id = Measure.mean) %>%
  # mutate(across(where(is.double), remove_outliers)) %>% #remove outliers (change to NA)
  left_join(scd_status, .)

dti_md <- aseg2dti_md %>%
  rename(participant_id = Measure.mean) %>%
  # mutate(across(where(is.double), remove_outliers)) %>% #remove outliers (change to NA)
  left_join(scd_status, .)

dki_kfa <- aseg2dki_kfa %>%
  rename(participant_id = Measure.mean) %>%
  # mutate(across(where(is.double), remove_outliers)) %>% #remove outliers (change to NA)
  left_join(scd_status, .)

dki_mk <- aseg2dki_mk %>%
  rename(participant_id = Measure.mean) %>%
  # mutate(across(where(is.double), remove_outliers)) %>% #remove outliers (change to NA)
  left_join(scd_status, .)

#extract SCD status from mti_over_55 table for each TR
scd_status_tr30 <- mti_over_55_tr30 %>% select(participant_id, SCD) %>%
  filter(!participant_id %in% failed_qc)
scd_status_tr30$SCD <- factor(scd_status_tr30$SCD,
                              levels = c(1,0),
                              labels = c('SCD', 'Control'))

scd_status_tr50 <- mti_over_55_tr50 %>% select(participant_id, SCD) %>%
  filter(!participant_id %in% failed_qc)
scd_status_tr50$SCD <- factor(scd_status_tr50$SCD,
                              levels = c(1,0),
                              labels = c('SCD', 'Control'))

#mtr stats for each TR
mtr_tr30 <- aseg2mtr %>%
  rename(participant_id = Measure.mean) %>%
  # mutate(across(where(is.double), remove_outliers)) %>% #remove outliers (change to NA)
  left_join(scd_status_tr30, .)

mtr_tr50 <- aseg2mtr %>%
  rename(participant_id = Measure.mean) %>%
  # mutate(across(where(is.double), remove_outliers)) %>% #remove outliers (change to NA)
  left_join(scd_status_tr50, .)

mtr_wm <- read_tsv("freesurfer/mtr_wm.tsv") %>%
  rename(participant_id=`Measure:mean`, mtr_wm=Seg0001)
mtr_gm <- read_tsv("freesurfer/mtr_gm.tsv") %>%
  rename(participant_id=`Measure:mean`, mtr_gm=Seg0001)
mtr_wm_gm_ratio <- mti_over_55 %>%
  select(participant_id, coil, mt_tr) %>%
  left_join(., mtr_wm) %>%
  left_join(., mtr_gm)
mtr_wm_gm_ratio$mtr_wm_gm_ratio <- mtr_wm_gm_ratio$mtr_wm / mtr_wm_gm_ratio$mtr_gm
mtr_wm_gm_ratio$coil <- as.factor(mtr_wm_gm_ratio$coil)
mtr_wm_gm_ratio$mt_tr <- factor(mtr_wm_gm_ratio$mt_tr,
                                levels = c(30, 50),
                                labels = c("TR=30ms", "TR=50ms"))
mtr_wm_gm_ratio %>% select(mt_tr, mtr_wm_gm_ratio) %>% tbl_summary(by = mt_tr, missing = "no") %>% add_p()

left_amygdala <- volumes %>% dplyr::select(c("participant_id", "SCD", "Left.Amygdala")) %>%
  rename(volume = "Left.Amygdala")
left_amygdala <- dti_fa %>% dplyr::select(c("participant_id", "SCD", "Left.Amygdala")) %>%
  rename(dti_fa = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- dti_md %>% dplyr::select(c("participant_id", "SCD", "Left.Amygdala")) %>%
  rename(dti_md = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- dki_kfa %>% dplyr::select(c("participant_id", "SCD", "Left.Amygdala")) %>%
  rename(dki_kfa = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- dki_mk %>% dplyr::select(c("participant_id", "SCD", "Left.Amygdala")) %>%
  rename(dki_mk = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- fit_FWF %>% dplyr::select(c("participant_id", "SCD", "Left.Amygdala")) %>%
  rename(fit_FWF = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- fit_NDI %>% dplyr::select(c("participant_id", "SCD", "Left.Amygdala")) %>%
  rename(fit_NDI = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- fit_ODI %>% dplyr::select(c("participant_id", "SCD", "Left.Amygdala")) %>%
  rename(fit_ODI = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- mtr_tr30 %>% dplyr::select(c("participant_id", "SCD", "Left.Amygdala")) %>%
  rename(mtr_tr30 = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- mtr_tr50 %>% dplyr::select(c("participant_id", "SCD", "Left.Amygdala")) %>%
  rename(mtr_tr50 = "Left.Amygdala") %>% full_join(left_amygdala, .)

left_amygdala %>% 
  select(!participant_id) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              # missing = "no"
              missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% 
  # filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

right_amygdala <- volumes %>% dplyr::select(c("participant_id", "SCD", "Right.Amygdala")) %>%
  rename(volume = "Right.Amygdala")
right_amygdala <- dti_fa %>% dplyr::select(c("participant_id", "SCD", "Right.Amygdala")) %>%
  rename(dti_fa = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- dti_md %>% dplyr::select(c("participant_id", "SCD", "Right.Amygdala")) %>%
  rename(dti_md = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- dki_kfa %>% dplyr::select(c("participant_id", "SCD", "Right.Amygdala")) %>%
  rename(dki_kfa = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- dki_mk %>% dplyr::select(c("participant_id", "SCD", "Right.Amygdala")) %>%
  rename(dki_mk = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- fit_FWF %>% dplyr::select(c("participant_id", "SCD", "Right.Amygdala")) %>%
  rename(fit_FWF = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- fit_NDI %>% dplyr::select(c("participant_id", "SCD", "Right.Amygdala")) %>%
  rename(fit_NDI = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- fit_ODI %>% dplyr::select(c("participant_id", "SCD", "Right.Amygdala")) %>%
  rename(fit_ODI = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- mtr_tr30 %>% dplyr::select(c("participant_id", "SCD", "Right.Amygdala")) %>%
  rename(mtr_tr30 = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- mtr_tr50 %>% dplyr::select(c("participant_id", "SCD", "Right.Amygdala")) %>%
  rename(mtr_tr50 = "Right.Amygdala") %>% full_join(right_amygdala, .)

right_amygdala %>% 
  select(!participant_id) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              # missing = "no"
              missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% 
  # filter_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Region of Interest**",
                estimate ~ "**Effect Size**")

