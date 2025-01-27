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


#extract SCD status and demographics from dwi_mti_over_55 table
scd_status <- dwi_mti_over_55 %>% select(participant_id, SCD, mt_tr, 
                                         Income, Ethnicity, Sex, age, age_education_completed) %>%
  filter(!participant_id %in% failed_qc)


#import aseg stats table (subcortical volumes)
volumes <- read_tsv("freesurfer/asegtable.tsv") %>%
  rename(participant_id=`Measure:volume`) %>%
  mutate(across(c(3:ncol(.)), .fns = ~.*1000/EstimatedTotalIntraCranialVol)) %>% #normalize by intracranial volume
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
  left_join(scd_status, .)
  

fit_NDI <- aseg2fit_NDI %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .)

fit_ODI <- aseg2fit_ODI %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .)

dti_fa <- aseg2dti_fa %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .)

dti_md <- aseg2dti_md %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .)

dki_kfa <- aseg2dki_kfa %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .)

dki_mk <- aseg2dki_mk %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .)

mtr <- aseg2mtr %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .)

left_amygdala <- volumes %>% dplyr::select(participant_id:age_education_completed, Left.Amygdala) %>%
  rename(volume = "Left.Amygdala")
left_amygdala <- fit_NDI %>% dplyr::select(participant_id:age_education_completed, Left.Amygdala) %>%
  rename(fit_NDI = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- dti_fa %>% dplyr::select(participant_id:age_education_completed, Left.Amygdala) %>%
  rename(dti_fa = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- dki_kfa %>% dplyr::select(participant_id:age_education_completed, Left.Amygdala) %>%
  rename(dki_kfa = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- dti_md %>% dplyr::select(participant_id:age_education_completed, Left.Amygdala) %>%
  rename(dti_md = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- dki_mk %>% dplyr::select(participant_id:age_education_completed, Left.Amygdala) %>%
  rename(dki_mk = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- mtr %>% dplyr::select(participant_id:age_education_completed, Left.Amygdala) %>%
  rename(mtr = "Left.Amygdala") %>% full_join(left_amygdala, .) %>%
  filter(!participant_id %in% failed_qc) %>%
  filter(mt_tr=="TR=30ms")

left_amygdala <- set_label(left_amygdala,
  volume = "Anatomical Volume",
  fit_NDI = "NODDI Neurite Density Index (NDI)",
  dti_fa = "DTI Fractional Anisotropy (FA)",
  dki_kfa = "DKI Kurtosis Fractional Anisotropy (KFA)",
  dti_md = "DTI Mean Diffusivity (MD)",
  dki_mk = "DKI Mean Kurtosis (MK)",
  mtr = "Magnetization Transfer Ratio"
)

right_amygdala <- volumes %>% dplyr::select(participant_id:age_education_completed, Right.Amygdala) %>%
  rename(volume = "Right.Amygdala")
right_amygdala <- fit_NDI %>% dplyr::select(participant_id:age_education_completed, Right.Amygdala) %>%
  rename(fit_NDI = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- dti_fa %>% dplyr::select(participant_id:age_education_completed, Right.Amygdala) %>%
  rename(dti_fa = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- dki_kfa %>% dplyr::select(participant_id:age_education_completed, Right.Amygdala) %>%
  rename(dki_kfa = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- dti_md %>% dplyr::select(participant_id:age_education_completed, Right.Amygdala) %>%
  rename(dti_md = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- dki_mk %>% dplyr::select(participant_id:age_education_completed, Right.Amygdala) %>%
  rename(dki_mk = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- mtr %>% dplyr::select(participant_id:age_education_completed, Right.Amygdala) %>%
  rename(mtr = "Right.Amygdala") %>% full_join(right_amygdala, .) %>%
  filter(!participant_id %in% failed_qc) %>%
  filter(mt_tr=="TR=30ms")

right_amygdala <- set_label(right_amygdala,
                           volume = "Anatomical Volume",
                           fit_NDI = "NODDI Neurite Density Index (NDI)",
                           dti_fa = "DTI Fractional Anisotropy (FA)",
                           dki_kfa = "DKI Kurtosis Fractional Anisotropy (KFA)",
                           dti_md = "DTI Mean Diffusivity (MD)",
                           dki_mk = "DKI Mean Kurtosis (MK)",
                           mtr = "Magnetization Transfer Ratio"
)

#group means of mtr tr30 subset
left_amygdala %>%
  select(SCD, Sex, Income, Ethnicity, age, age_education_completed) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})") %>% add_p()

left_amygdala_table <- left_amygdala %>% 
  select(SCD, volume, fit_NDI, dti_fa, dki_kfa, dti_md, dki_mk, mtr) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              # missing = "no"
              missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% 
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Metric**",
                estimate ~ "**Effect Size**")

right_amygdala_table <- right_amygdala %>% 
  select(SCD, volume, fit_NDI, dti_fa, dki_kfa, dti_md, dki_mk, mtr) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              # missing = "no"
              missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>%  bold_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Metric**",
                estimate ~ "**Effect Size**")

tbl_merge(list(left_amygdala_table, right_amygdala_table),
          tab_spanner = c("**Left Amygdala**", "**Right Amygdala**"))

#scatterplots: neurodegen and demyelination
NDI_by_KFA <- lm(fit_NDI ~ dki_kfa + SCD + Sex + Income + age_education_completed, right_amygdala)
pr <- predict_response(NDI_by_KFA, c( "dki_kfa [all]", "SCD[all]"))
plot(pr, show_residuals = T) + theme(legend.title = element_blank())
NDI_by_MK <- lm(fit_NDI ~ dki_mk + SCD + Sex + Income + age_education_completed, right_amygdala)
pr <- predict_response(NDI_by_MK, c( "dki_mk [all]", "SCD[all]"))
plot(pr, show_residuals = T) + theme(legend.title = element_blank())
MTR_by_KFA <- lm(mtr ~ dki_kfa + SCD + Sex + Income + age_education_completed, right_amygdala)
pr <- predict_response(MTR_by_KFA, c( "dki_kfa [all]", "SCD[all]"))
plot(pr, show_residuals = T) + theme(legend.title = element_blank())
MTR_by_MK <- lm(mtr ~ dki_mk + SCD + Sex + Income + age_education_completed, right_amygdala)
pr <- predict_response(MTR_by_MK, c( "dki_mk [all]", "SCD[all]"))
plot(pr, show_residuals = T) + theme(legend.title = element_blank())

