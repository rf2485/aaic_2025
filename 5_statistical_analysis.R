source("1_data_preparation.R")
library(gtsummary)
library(ggeffects)
library(ggtext)

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

#read in diffusion and MTR tables
aparc2meas_files <- list.files(path = "freesurfer", 
                               pattern = "aparc.*aseg.*\\.*tsv", 
                               full.names = T)
for (i in 1:length(aparc2meas_files)) {
  assign(gsub(".tsv", "", 
              gsub("freesurfer.aparc.", "", make.names(aparc2meas_files[i]))), 
         read.delim(aparc2meas_files[i]))
}

fit_NDI <- aseg2fit_NDI %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .)

fit_ODI <- aseg2fit_ODI %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .)

dki_kfa <- aseg2dki_kfa %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .)

dki_mk <- aseg2dki_mk %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .)

dki_ak <- aseg2dki_ak %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .)

dki_rk <- aseg2dki_rk %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .)

mtr <- aseg2mtr %>%
  rename(participant_id = Measure.mean) %>%
  left_join(scd_status, .)

left_amygdala <- fit_NDI %>% dplyr::select(participant_id:age_education_completed, Left.Amygdala) %>%
  rename(fit_NDI = "Left.Amygdala")
left_amygdala <- fit_ODI %>% dplyr::select(participant_id:age_education_completed, Left.Amygdala) %>%
  rename(fit_ODI = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- dki_kfa %>% dplyr::select(participant_id:age_education_completed, Left.Amygdala) %>%
  rename(dki_kfa = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- dki_mk %>% dplyr::select(participant_id:age_education_completed, Left.Amygdala) %>%
  rename(dki_mk = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- dki_ak %>% dplyr::select(participant_id:age_education_completed, Left.Amygdala) %>%
  rename(dki_ak = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- dki_rk %>% dplyr::select(participant_id:age_education_completed, Left.Amygdala) %>%
  rename(dki_rk = "Left.Amygdala") %>% full_join(left_amygdala, .)
left_amygdala <- mtr %>% dplyr::select(participant_id:age_education_completed, Left.Amygdala) %>%
  rename(mtr = "Left.Amygdala") %>% full_join(left_amygdala, .) %>%
  filter(!participant_id %in% failed_qc) %>%
  filter(mt_tr=="TR=30ms")

left_amygdala <- set_label(left_amygdala,
  fit_NDI = "NODDI Neurite Density (ND)",
  fit_ODI = "NODDI Orientation Dispersion (OD)",
  dki_kfa = "DKI Kurtosis Fractional Anisotropy (KFA)",
  dki_mk = "DKI Mean Kurtosis (MK)",
  dki_ak = "DKI Axial Kurtosis (AK)",
  dki_rk = "DKI Radial Kurtosis (RK)",
  mtr = "Magnetization Transfer Ratio (MTR)"
)

right_amygdala <- fit_NDI %>% dplyr::select(participant_id:age_education_completed, Right.Amygdala) %>%
  rename(fit_NDI = "Right.Amygdala")
right_amygdala <- fit_ODI %>% dplyr::select(participant_id:age_education_completed, Right.Amygdala) %>%
  rename(fit_ODI = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- dki_kfa %>% dplyr::select(participant_id:age_education_completed, Right.Amygdala) %>%
  rename(dki_kfa = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- dki_mk %>% dplyr::select(participant_id:age_education_completed, Right.Amygdala) %>%
  rename(dki_mk = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- dki_ak %>% dplyr::select(participant_id:age_education_completed, Right.Amygdala) %>%
  rename(dki_ak = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- dki_rk %>% dplyr::select(participant_id:age_education_completed, Right.Amygdala) %>%
  rename(dki_rk = "Right.Amygdala") %>% full_join(right_amygdala, .)
right_amygdala <- mtr %>% dplyr::select(participant_id:age_education_completed, Right.Amygdala) %>%
  rename(mtr = "Right.Amygdala") %>% full_join(right_amygdala, .) %>%
  filter(!participant_id %in% failed_qc) %>%
  filter(mt_tr=="TR=30ms")

right_amygdala <- set_label(right_amygdala,
                           fit_NDI = "NODDI Neurite Density (ND)",
                           fit_ODI = "NODDI Orientation Dispersion (OD)",
                           dki_kfa = "DKI Kurtosis Fractional Anisotropy (KFA)",
                           dki_mk = "DKI Mean Kurtosis (MK)",
                           dki_ak = "DKI Axial Kurtosis (AK)",
                           dki_rk = "DKI Radial Kurtosis (RK)",
                           mtr = "Magnetization Transfer Ratio (MTR)"
)

#group means of mtr tr30 subset
left_amygdala %>%
  select(SCD, Sex, Income, Ethnicity, age, age_education_completed) %>%
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})") %>% add_p()

left_amygdala_table <- left_amygdala %>% 
  select(SCD, fit_NDI, fit_ODI, dki_kfa, dki_mk, dki_ak, dki_rk, mtr) %>% 
  tbl_summary(by = SCD, statistic = all_continuous() ~ "{mean} ({sd})",
              # missing = "no"
              missing_text = "Excluded Outliers"
  ) %>%
  add_difference(test = list(everything() ~ 'cohens_d')) %>%
  modify_column_hide(conf.low) %>%
  add_p() %>% add_q() %>% bold_p(q=T) %>% bold_p() %>%
  modify_header(statistic ~ "**Test Statistic**", 
                label ~ "**Metric**",
                estimate ~ "**Effect Size**")

right_amygdala_table <- right_amygdala %>% 
  select(SCD, fit_NDI, fit_ODI, dki_kfa, dki_mk, dki_ak, dki_rk, mtr) %>% 
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

#so education can be compared to other covariates:
right_amygdala_na_omit <- mutate(right_amygdala,
                      across(where(is.factor),~fct_explicit_na(.,'Unknown'))) %>%
  na.omit(.)

#scatterplots: neurodegen and demyelination/protein aggregation
#first assess if each covariate improves the model. Discard if ANOVA is not significant.
KFA_by_NDI <- lm(dki_kfa ~ fit_NDI, right_amygdala_na_omit)
KFA_by_NDI_scd <- lm(dki_kfa ~ fit_NDI + SCD, right_amygdala_na_omit)
anova(KFA_by_NDI, KFA_by_NDI_scd)
KFA_by_NDI_age <- lm(dki_kfa ~ fit_NDI + age, right_amygdala_na_omit)
anova(KFA_by_NDI, KFA_by_NDI_age) #improvement
KFA_by_NDI_age_sex <- lm(dki_kfa ~ fit_NDI + age + Sex, right_amygdala_na_omit)
anova(KFA_by_NDI_age, KFA_by_NDI_age_sex) #improvement
KFA_by_NDI_age_sex_education <- lm(dki_kfa ~ fit_NDI + age + Sex + age_education_completed, right_amygdala_na_omit)
anova(KFA_by_NDI_age_sex, KFA_by_NDI_age_sex_education)
KFA_by_NDI_age_sex_income <- lm(dki_kfa ~ fit_NDI + age + Sex + Income, right_amygdala_na_omit)
anova(KFA_by_NDI_age_sex, KFA_by_NDI_age_sex_income)
KFA_by_NDI_age_sex_ethnicity <- lm(dki_kfa ~ fit_NDI + age + Sex + Ethnicity, right_amygdala_na_omit)
anova(KFA_by_NDI_age_sex, KFA_by_NDI_age_sex_ethnicity)

#image width 570 height 481
#only education has NAs and we are not using it, so switch back to full dataset
KFA_by_NDI_age_sex <- lm(dki_kfa ~ fit_NDI + age + Sex, right_amygdala)
adj_r_squared <- summary(KFA_by_NDI_age_sex)$adj.r.squared
model_p_value <- pf(summary(KFA_by_NDI_age_sex)$fstatistic[1], 
                    summary(KFA_by_NDI_age_sex)$fstatistic[2], 
                    summary(KFA_by_NDI_age_sex)$fstatistic[3], 
                    lower.tail=F)
pr <- predict_response(KFA_by_NDI_age_sex, c( "fit_NDI [all]"))
plot(pr, show_data = T, dot_alpha = 1, colors = "darkorange") + 
  labs(
    title = "A. Mean Right Amygdala KFA and ND, Corrected by Age and Sex",
    subtitle = paste0("** p = ", signif(model_p_value, 2),
                      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5, face = "bold"))

MK_by_NDI_age_sex <- lm(dki_mk ~ fit_NDI + age + Sex, right_amygdala)
adj_r_squared <- summary(MK_by_NDI_age_sex)$adj.r.squared
model_p_value <- pf(summary(MK_by_NDI_age_sex)$fstatistic[1], 
                    summary(MK_by_NDI_age_sex)$fstatistic[2], 
                    summary(MK_by_NDI_age_sex)$fstatistic[3], 
                    lower.tail=F)
pr <- predict_response(MK_by_NDI_age_sex, c( "fit_NDI [all]"))
plot(pr, show_data = T, dot_alpha = 1, colors = "darkorange") + 
  labs(
    title = "B. Mean Right Amygdala MK and ND, Corrected by Age and Sex",
    subtitle = paste0("*** p < 0.001",
                      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5, face = "bold"))

RK_by_NDI_age_sex <- lm(dki_rk ~ fit_NDI + age + Sex, right_amygdala)
adj_r_squared <- summary(RK_by_NDI_age_sex)$adj.r.squared
model_p_value <- pf(summary(RK_by_NDI_age_sex)$fstatistic[1], 
                    summary(RK_by_NDI_age_sex)$fstatistic[2], 
                    summary(RK_by_NDI_age_sex)$fstatistic[3], 
                    lower.tail=F)
pr <- predict_response(RK_by_NDI_age_sex, c( "fit_NDI [all]"))
plot(pr, show_data = T, dot_alpha = 1, colors = "darkorange") + 
  labs(
    title = "C. Mean Right Amygdala RK and ND, Corrected by Age and Sex",
    subtitle = paste0("*** p < 0.001",
                      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5, face = "bold"))

KFA_by_ODI_age_sex <- lm(dki_kfa ~ fit_ODI + age + Sex, right_amygdala)
adj_r_squared <- summary(KFA_by_ODI_age_sex)$adj.r.squared
model_p_value <- pf(summary(KFA_by_ODI_age_sex)$fstatistic[1], 
                    summary(KFA_by_ODI_age_sex)$fstatistic[2], 
                    summary(KFA_by_ODI_age_sex)$fstatistic[3], 
                    lower.tail=F)
pr <- predict_response(KFA_by_ODI_age_sex, c( "fit_ODI [all]"))
plot(pr, show_data = T, dot_alpha = 1, colors = "purple") + 
  labs(
    title = "D. Mean Right Amygdala KFA and OD, Corrected by Age and Sex",
    subtitle = paste0("p = ", signif(model_p_value, 2),
                      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5))

MK_by_ODI_age_sex <- lm(dki_mk ~ fit_ODI + age + Sex, right_amygdala)
adj_r_squared <- summary(MK_by_ODI_age_sex)$adj.r.squared
model_p_value <- pf(summary(MK_by_ODI_age_sex)$fstatistic[1], 
                    summary(MK_by_ODI_age_sex)$fstatistic[2], 
                    summary(MK_by_ODI_age_sex)$fstatistic[3], 
                    lower.tail=F)
pr <- predict_response(MK_by_ODI_age_sex, c( "fit_ODI [all]"))
plot(pr, show_data = T, dot_alpha = 1, colors = "purple") + 
  labs(
    title = "E. Mean Right Amygdala MK and OD, Corrected by Age and Sex",
    subtitle = paste0("\\* p = ", signif(model_p_value, 2),
                      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5, face = "bold"))

RK_by_ODI_age_sex <- lm(dki_rk ~ fit_ODI + age + Sex, right_amygdala)
adj_r_squared <- summary(RK_by_ODI_age_sex)$adj.r.squared
model_p_value <- pf(summary(RK_by_ODI_age_sex)$fstatistic[1], 
                    summary(RK_by_ODI_age_sex)$fstatistic[2], 
                    summary(RK_by_ODI_age_sex)$fstatistic[3], 
                    lower.tail=F)
pr <- predict_response(RK_by_ODI_age_sex, c( "fit_ODI [all]"))
plot(pr, show_data = T, dot_alpha = 1, colors = "purple") + 
  labs(
    title = "F. Mean Right Amygdala RK and OD, Corrected by Age and Sex",
    subtitle = paste0("p = ", signif(model_p_value, 2),
                      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5))

KFA_by_MTR_age_sex <- lm(dki_kfa ~ mtr + age + Sex, right_amygdala)
adj_r_squared <- summary(KFA_by_MTR_age_sex)$adj.r.squared
model_p_value <- pf(summary(KFA_by_MTR_age_sex)$fstatistic[1], 
                    summary(KFA_by_MTR_age_sex)$fstatistic[2], 
                    summary(KFA_by_MTR_age_sex)$fstatistic[3], 
                    lower.tail=F)
pr <- predict_response(KFA_by_MTR_age_sex, c( "mtr [all]"))
plot(pr, show_data = T, dot_alpha = 1, colors = "turquoise3") + 
  labs(
    title = "G. Mean Right Amygdala KFA and MTR, Corrected by Age and Sex",
    subtitle = paste0("p = ", signif(model_p_value, 2),
                      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5))

MK_by_MTR_age_sex <- lm(dki_mk ~ mtr + age + Sex, right_amygdala)
adj_r_squared <- summary(MK_by_MTR_age_sex)$adj.r.squared
model_p_value <- pf(summary(MK_by_MTR_age_sex)$fstatistic[1], 
                    summary(MK_by_MTR_age_sex)$fstatistic[2], 
                    summary(MK_by_MTR_age_sex)$fstatistic[3], 
                    lower.tail=F)
pr <- predict_response(MK_by_MTR_age_sex, c( "mtr [all]"))
plot(pr, show_data = T, dot_alpha = 1, colors = "turquoise3") + 
  labs(
    title = "H. Mean Right Amygdala MK and MTR, Corrected by Age and Sex",
    subtitle = paste0("p = ", signif(model_p_value, 2),
                      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5))


RK_by_MTR_age_sex <- lm(dki_rk ~ mtr + age + Sex, right_amygdala)
adj_r_squared <- summary(RK_by_MTR_age_sex)$adj.r.squared
model_p_value <- pf(summary(RK_by_MTR_age_sex)$fstatistic[1], 
                    summary(RK_by_MTR_age_sex)$fstatistic[2], 
                    summary(RK_by_MTR_age_sex)$fstatistic[3], 
                    lower.tail=F)
pr <- predict_response(RK_by_MTR_age_sex, c( "mtr [all]"))
plot(pr, show_data = T, dot_alpha = 1, colors = "turquoise3") + 
  labs(
    title = "I. Mean Right Amygdala RK and MTR, Corrected by Age and Sex",
    subtitle = paste0("p = ", signif(model_p_value, 2),
                      ", adj-R<sup>2</sup> = ", signif(adj_r_squared, 2))) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5))
