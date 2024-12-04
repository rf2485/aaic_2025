source("1_data_preparation.R")
library(arsenal)

failed_qc <- c('sub-CC510255', #abnormality in left temporal pole
               'sub-CC510438', #abnormality in left frontal lobe
               # 'sub-CC610308', #parietal lobe cutoff
               # 'sub-CC610469', #parietal lobe cutoff
               # 'sub-CC620466', #parietal lobe cutoff
               'sub-CC620821', #segmentation errors from large ventricles
               'sub-CC621011', #segmentation errors from large ventricles
               'sub-CC621080', #segmentation errors
               # 'sub-CC710214', #parietal lobe cutoff
               'sub-CC710551', #motion artifacts in DWI
               'sub-CC711027', #severe motion artifacts in T1
               # 'sub-CC712027', #parietal lobe cutoff
               'sub-CC721434' #segmentation errors from large ventricles
)

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

##wm-gm ratio of MTR for each TR
mtr_wm <- read_tsv("freesurfer/mtr_wm.tsv") %>%
  rename(participant_id=`Measure:mean`, mtr_wm=Seg0001) 
mtr_gm <- read_tsv("freesurfer/mtr_gm.tsv") %>%
  rename(participant_id=`Measure:mean`, mtr_gm=Seg0001)
mtr_wm_gm_ratio <- mti_over_55 %>%
  select(participant_id, coil, mt_tr) %>%
  left_join(., mtr_wm) %>%
  left_join(., mtr_gm) %>%
  filter(!participant_id %in% failed_qc) %>%
  mutate(across(where(is.double), remove_outliers))
mtr_wm_gm_ratio$mtr_wm_gm_ratio <- mtr_wm_gm_ratio$mtr_wm / mtr_wm_gm_ratio$mtr_gm
mtr_wm_gm_ratio$coil <- as.factor(mtr_wm_gm_ratio$coil)
mtr_wm_gm_ratio$mt_tr <- as.factor(mtr_wm_gm_ratio$mt_tr)
mtr_wm_gm_ratio_table <- tableby(mt_tr ~ mtr_wm_gm_ratio,
                                 data = mtr_wm_gm_ratio, numeric.test="kwt", total = F)
summary(mtr_wm_gm_ratio_table, text = T)
#TR=30ms has best contrast, using only subjects with TR=30ms mti images

#extract SCD status from mti_over_55 table for subjects with mt_tr==30
scd_status <- mti_over_55 %>% 
  filter(!participant_id %in% failed_qc & mt_tr==30) %>%
  select(participant_id, SCD)
scd_status$SCD <- factor(scd_status$SCD,
                         levels = c(1,0),
                         labels = c('SCD', 'Control'))

#import aseg stats table (subcortical volumes)
aseg = read_tsv("freesurfer/asegtable.tsv") %>%
  rename(participant_id=`Measure:volume`) #%>%
names(aseg) <- make.names(names(aseg))
# aseg_table <- tableby(formulize('SCD', names(aseg)[3:54]), 
#                          data = aseg, numeric.test="wt", total = F) 
# summary(aseg_table, text = TRUE)

#import JHU stats table (WM volumes)
jhu_volume <- read_tsv("freesurfer/jhu_volume.tsv") %>%
  rename(participant_id=`Measure:volume`,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049
         ) %>%
  select(!starts_with("Seg00"))
# jhu_volume[jhu_volume == 0] <- NA
volumes <- left_join(scd_status, jhu_volume) %>% #join JHU volumes with SCD status
  left_join(., aseg) %>%
  mutate(across(c(3:ncol(.)), .fns = ~.*1000/EstimatedTotalIntraCranialVol)) %>% #normalize by intracranial volume
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
volumes_table <- tableby(formulize('SCD', names(volumes)[3:88]), 
                      data = volumes, numeric.test="wt", total = F) 
summary(volumes_table, text = TRUE)

####diffusion and MTR group means
#read in diffusion and MTR tables
aparc2meas_files <- list.files(path = "freesurfer", 
                               pattern = "aparc.*aseg.*\\.*tsv", 
                               full.names = T)
for (i in 1:length(aparc2meas_files)) {
  assign(gsub(".tsv", "", 
              gsub("freesurfer.aparc.", "", make.names(aparc2meas_files[i]))), 
         read.delim(aparc2meas_files[i]))
}
jhu2meas_files <- list.files(path = "freesurfer",
                             pattern = "jhu2.*\\.*tsv",
                             full.names = T)
for (i in 1:length(jhu2meas_files)) {
  assign(gsub(".tsv", "",
              gsub("freesurfer.", "", make.names(jhu2meas_files[i]))),
         read.delim(jhu2meas_files[i]))
}

dti_fa <- aseg2dti_fa %>%
  left_join(., jhu2dti_fa) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049) %>%
  left_join(scd_status, .) %>%
  select(!CSF & !ends_with("Ventricle") & !starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
fa_table <- tableby(formulize('SCD', names(dti_fa)[3:ncol(dti_fa)]),
                     data = dti_fa, numeric.test="wt", total = FALSE)
summary(fa_table, text = T)

dti_md <- aseg2dti_md %>%
  left_join(., jhu2dti_md) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049
  ) %>%
  left_join(scd_status, .) %>%
  select(!CSF & !ends_with("Ventricle") & !starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
md_table <- tableby(formulize('SCD', names(dti_md)[3:ncol(dti_md)]),
                    data = dti_md, numeric.test="wt", total = FALSE)
summary(md_table, text = T)

dti_rd <- aseg2dti_rd %>%
  left_join(., jhu2dti_rd) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049
  ) %>%
  left_join(scd_status, .) %>%
  select(!CSF & !ends_with("Ventricle") & !starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
rd_table <- tableby(formulize('SCD', names(dti_rd)[3:ncol(dti_rd)]),
                    data = dti_rd, numeric.test="wt", total = FALSE)
summary(rd_table, text = T)

dti_ad <- aseg2dti_ad %>%
  left_join(., jhu2dti_ad) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049
  ) %>%
  left_join(scd_status, .) %>%
  select(!CSF & !ends_with("Ventricle") & !starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
ad_table <- tableby(formulize('SCD', names(dti_ad)[3:ncol(dti_ad)]),
                    data = dti_ad, numeric.test="wt", total = FALSE)
summary(ad_table, text = T)

dki_kfa <- aseg2dki_kfa %>%
  left_join(., jhu2dki_kfa) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049
  ) %>%
  left_join(scd_status, .) %>%
  select(!CSF & !ends_with("Ventricle") & !starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
kfa_table <- tableby(formulize('SCD', names(dki_kfa)[3:ncol(dki_kfa)]),
                    data = dki_kfa, numeric.test="wt", total = FALSE)
summary(kfa_table, text = T)

dki_mk <- aseg2dki_mk %>%
  left_join(., jhu2dki_mk) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049
  ) %>%
  left_join(scd_status, .) %>%
  select(!CSF & !ends_with("Ventricle") & !starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
mk_table <- tableby(formulize('SCD', names(dki_mk)[3:ncol(dki_mk)]),
                     data = dki_mk, numeric.test="wt", total = FALSE)
summary(mk_table, text = T)

dki_rk <- aseg2dki_rk %>%
  left_join(., jhu2dki_rk) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049
  ) %>%
  left_join(scd_status, .) %>%
  select(!CSF & !ends_with("Ventricle") & !starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
rk_table <- tableby(formulize('SCD', names(dki_rk)[3:ncol(dki_rk)]),
                     data = dki_rk, numeric.test="wt", total = FALSE)
summary(rk_table, text = T)

dki_ak <- aseg2dki_ak %>%
  left_join(., jhu2dki_ak) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049
  ) %>%
  left_join(scd_status, .) %>%
  select(!CSF & !ends_with("Ventricle") & !starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
ak_table <- tableby(formulize('SCD', names(dki_ak)[3:ncol(dki_ak)]),
                     data = dki_ak, numeric.test="wt", total = FALSE)
summary(ak_table, text = T)

smi_matlab_Da <- aseg2smi_matlab_Da %>%
  left_join(., jhu2smi_matlab_Da) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049
  ) %>%
  left_join(scd_status, .) %>%
  select(!starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
Da_table <- tableby(formulize('SCD', names(smi_matlab_Da)[3:ncol(smi_matlab_Da)]),
                     data = smi_matlab_Da, numeric.test="wt", total = FALSE)
summary(Da_table, text = T)

smi_matlab_DePar <- aseg2smi_matlab_DePar %>%
  left_join(., jhu2smi_matlab_DePar) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049
  ) %>%
  left_join(scd_status, .) %>%
  select(!starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
DePar_table <- tableby(formulize('SCD', names(smi_matlab_DePar)[3:ncol(smi_matlab_DePar)]),
                    data = smi_matlab_DePar, numeric.test="wt", total = FALSE)
summary(DePar_table, text = T)

smi_matlab_DePerp <- aseg2smi_matlab_DePerp %>%
  left_join(., jhu2smi_matlab_DePerp) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049
  ) %>%
  left_join(scd_status, .) %>%
  select(!starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
DePerp_table <- tableby(formulize('SCD', names(smi_matlab_DePerp)[3:ncol(smi_matlab_DePerp)]),
                       data = smi_matlab_DePerp, numeric.test="wt", total = FALSE)
summary(DePerp_table, text = T)

smi_matlab_f <- aseg2smi_matlab_f %>%
  left_join(., jhu2smi_matlab_f) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049
  ) %>%
  left_join(scd_status, .) %>%
  select(!starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
f_table <- tableby(formulize('SCD', names(smi_matlab_f)[3:ncol(smi_matlab_f)]),
                       data = smi_matlab_f, numeric.test="wt", total = FALSE)
summary(f_table, text = T)

smi_matlab_p2 <- aseg2smi_matlab_p2 %>%
  left_join(., jhu2smi_matlab_p2) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049
  ) %>%
  left_join(scd_status, .) %>%
  select(!starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
p2_table <- tableby(formulize('SCD', names(smi_matlab_p2)[3:ncol(smi_matlab_p2)]),
                       data = smi_matlab_p2, numeric.test="wt", total = FALSE)
summary(p2_table, text = T)

mtr_tr30 <- aseg2mtr %>%
  left_join(., jhu2mtr) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049
  ) %>%
  left_join(scd_status, .) %>%
  select(!CSF & !ends_with("Ventricle") & !ends_with("Vent") & !starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
mtr_tr30_table <- tableby(formulize('SCD', names(mtr_tr30)[3:ncol(mtr_tr30)]),
                          data = mtr_tr30, numeric.test="wt", total = FALSE)
summary(mtr_tr30_table, text = T)

g_ratio_tr30 <- aseg2g_ratio %>%
  left_join(., jhu2g_ratio) %>%
  rename(participant_id = Measure.mean,
         middle_cerebellar_peduncle=Seg0001,
         genu_corpus_callosum=Seg0003,
         body_corpus_callosum=Seg0004,
         splenium_corpus_callosum=Seg0005,
         inferior_cerebellar_peduncle_R=Seg0011,
         inferior_cerebellar_peduncle_L=Seg0012,
         superior_cerebellar_peduncle_R=Seg0013,
         superior_cerebellar_peduncle_L=Seg0014,
         cerebral_peduncle_R=Seg0015,
         cerebral_peduncle_L=Seg0016,
         anterior_limb_internal_capsule_R=Seg0017,
         anterior_limb_internal_capsule_L=Seg0018,
         posterior_limb_internal_capsule_R=Seg0019,
         posterior_limb_internal_capsule_L=Seg0020,
         retrolenticular_part_internal_capsule_R=Seg0021,
         retrolenticular_part_internal_capsule_L=Seg0022,
         anterior_corona_radiata_R=Seg0023,
         anterior_corona_radiata_L=Seg0024,
         superior_corona_radiata_R=Seg0025,
         superior_corona_radiata_L=Seg0026,
         posterior_corona_radiata_R=Seg0027,
         posterior_corona_radiata_L=Seg0028,
         posterior_thalamic_radiation_R=Seg0029,
         posterior_thalamic_radiation_L=Seg0030,
         sagittal_stratum_R=Seg0031,
         sagittal_stratum_L=Seg0032,
         external_capsule_R=Seg0033,
         external_capsule_L=Seg0034,
         upper_cingulum_R=Seg0035,
         upper_cingulum_L=Seg0036,
         lower_cingulum_R=Seg0037,
         lower_cingulum_L=Seg0038,
         fornix_cres_R=Seg0039,
         fornix_cres_L=Seg0040,
         superior_longitudinal_fasciculus_R=Seg0041,
         superior_longitudinal_fasciculus_L=Seg0042,
         superior_fronto_occipital_fasciculus_R=Seg0043,
         superior_fronto_occipital_fasciculus_L=Seg0044,
         inferior_fronto_occipital_fasciculus_R=Seg0045,
         inferior_fronto_occipital_fasciculus_L=Seg0046,
         uncinate_fasciculus_R=Seg0047,
         uncinate_fasciculus_L=Seg0048,
         tapetum_R=Seg0049
  ) %>%
  left_join(scd_status, .) %>%
  select(!starts_with("Seg00")) %>%
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
g_ratio_tr30_table <- tableby(formulize('SCD', names(g_ratio_tr30)[3:ncol(g_ratio_tr30)]),
                          data = g_ratio_tr30, numeric.test="wt", total = FALSE
                          )
summary(g_ratio_tr30_table, text = T)
