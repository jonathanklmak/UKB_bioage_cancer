#==============================================================================
# FILENAME:   3_analysis_UKB_BioAge_cancer.R
# PROJECT:    UKB_bioage_cancer
# PURPOSE:    To assess the association between BioAge and cancers in UKB
# AUTHOR:     Jonathan Mak
# CREATED:    2022-07-14
# UPDATED:    2023-01-25
# R VERSION:  4.1.3
#==============================================================================

setwd("Z:/UKB_bioage_cancer")

### Required packages
library(haven)
library(dplyr)
library(ggplot2)
library(cowplot)      # For combining ggplots
library(ggforestplot)
library(psych)
library(broom)        
library(gtsummary)
library(htmlTable)
library(tableone)     # For creating Table 1
library(survival)     # Survival analysis
#install.packages(c("Hmisc","rms"), dependencies=TRUE)
library(Hmisc); library(rms); library(splines) # For spline models


#==================== LOAD UKB BIOAGE AND CANCER DATA =========================

### Import pre-processed bioage-cancer data (performed in Stata, "Program/2.2_data_prep_UKB_BioAge")
ukb_bioage_cancer <- read_dta("Data/Cleaned_data/UKB_BioAge_cancer.dta")
ukb_bioage_cancer <- haven::as_factor(ukb_bioage_cancer, levels="label") # Change variables to factors using labels in Stata file
colnames(ukb_bioage_cancer); dim(ukb_bioage_cancer) # Total N = 308,156

### Create a "complete-case" dataset for sensitivity analysis (convert all values with "missing" indicator to NA)
ukb_bioage_cancer_complete <- ukb_bioage_cancer %>%  mutate(across(where(is.factor), ~na_if(., "missing")))

### Define covariates to be included in Cox models
cov_ukb_age_sex <- c("birth_year_cat","sex") # Age- and sex-adjusted model
cov_ukb_multivariable <- c("birth_year_cat","sex","assessment_centre","ethnicity_4group","bmi_cat","smoking","alcohol_3group","PA_cat","education","di_5group") # Multivariable model
cov_ukb_bc <- c("birth_year_cat","sex","assessment_centre","ethnicity_4group","bmi_cat","smoking","alcohol_3group","PA_cat","education","di_5group","family_history_breast_cancer","bc_screening","menopause","hormone_therapy","oral_contraceptive","parity") # Breast-cancer specific model
cov_ukb_pc <- c("birth_year_cat","sex","assessment_centre","ethnicity_4group","bmi_cat","smoking","alcohol_3group","PA_cat","education","di_5group","family_history_prostate_cancer","psa_test","diabetes") # Prostate-cancer specific model
cov_ukb_lc <- c("birth_year_cat","sex","assessment_centre","ethnicity_4group","bmi_cat","smoking","alcohol_3group","PA_cat","education","di_5group","family_history_lung_cancer") # Lung-cancer specific model
cov_ukb_cc <- c("birth_year_cat","sex","assessment_centre","ethnicity_4group","bmi_cat","smoking","alcohol_3group","PA_cat","education","di_5group","family_history_bowel_cancer","cc_screening","fresh_veg_cat","red_meat_cat","processed_meat_cat") # Colorectal-cancer specific model
cov_ukb_mel <- c("birth_year_cat","sex","assessment_centre","ethnicity_4group","bmi_cat","smoking","alcohol_3group","PA_cat","education","di_5group","time_outdoors_summer_4group","childhood_sunburn","solarium_use","skin_tan_ease","skin_color_4group","hair_color_3group","sun_protection") # Melanoma specific model


#====================== DESCRIPTIVE STATISTICS ================================

### Sample characteristics by sex in UKB
table_var <- c("age_baseline","kdm","kdm_res","phenoage","phenoage_res","hd_log","fev_1000","sbp","bun","hba1c","totchol","creat_umol","glucose_mmol","waist","rdw","albumin_gl","alp","trig","mcv","uap","lymph","rbc","crp","dbp","died","any_cancer_diag","bc_diag","pc_diag","lc_diag","cc_diag","mel_diag")
table_catvar <- c("died","any_cancer_diag","bc_diag","pc_diag","lc_diag","cc_diag","mel_diag")
print(CreateTableOne(vars = table_var, data = ukb_bioage_cancer, factorVars = table_catvar, strata="sex", addOverall=T, includeNA = T),
      formatOptions = list(big.mark = ","), quote = TRUE, noSpaces = TRUE)
rm(list=c("table_catvar","table_var"))

### Detailed sample characteristics
table_var <- c("birth_year_cat","assessment_centre","ethnicity_4group","bmi_cat",
               "smoking","alcohol_3group","fresh_veg_cat","red_meat_cat","processed_meat_cat",
               "PA_cat","education","di_5group","diabetes",
               "menopause","hormone_therapy","oral_contraceptive","parity",
               "bc_screening","psa_test","cc_screening",
               "family_history_breast_cancer","family_history_prostate_cancer",
               "family_history_lung_cancer","family_history_bowel_cancer",
               "time_outdoors_summer_4group","childhood_sunburn","solarium_use",
               "skin_tan_ease","skin_color_4group","hair_color_3group","sun_protection")                 
table_catvar <- c("birth_year_cat","assessment_centre","ethnicity_4group","bmi_cat",
                  "smoking","alcohol_3group","fresh_veg_cat","red_meat_cat","processed_meat_cat",
                  "PA_cat","education","di_5group","diabetes",
                  "menopause","hormone_therapy","oral_contraceptive","parity",
                  "bc_screening","psa_test","cc_screening",
                  "family_history_breast_cancer","family_history_prostate_cancer",
                  "family_history_lung_cancer","family_history_bowel_cancer",
                  "time_outdoors_summer_4group","childhood_sunburn","solarium_use",
                  "skin_tan_ease","skin_color_4group","hair_color_3group","sun_protection")
print(CreateTableOne(vars = table_var, data = ukb_bioage_cancer, factorVars = table_catvar, strata="sex", addOverall=T, includeNA = T),
      formatOptions = list(big.mark = ","), quote = TRUE, noSpaces = TRUE)
rm(list=c("table_catvar","table_var"))

### Correlations among the 18 biomarkers
round(cor(ukb_bioage_cancer[,c("fev_1000","sbp","bun","hba1c","totchol","creat_umol","glucose_mmol","waist","rdw","albumin_gl","alp","trig","mcv","uap","lymph","rbc","crp","dbp")],use="pairwise.complete.obs"), 2)
round(rcorr(as.matrix(ukb_bioage_cancer[,c("fev_1000","sbp","bun","hba1c","totchol","creat_umol","glucose_mmol","waist","rdw","albumin_gl","alp","trig","mcv","uap","lymph","rbc","crp","dbp")]), type="pearson")$P,3)



#================= COX MODELS FOR BA MEASURES AND CANCERS =====================

# BA measures to be tested
BA_measures <- c("kdm_res_sd", # KDM based on 18 biomarkers
                 "kdm_original_res_sd", # Original Levine KDM
                 "kdm_noglu_res_sd", # KDM, excluding glucose and HbA1c (sensitivity analysis)
                 "phenoage_res_sd", # PhenoAge based on 18 biomarkers
                 "phenoage_original_res_sd", # Original Levine ?PhenoAge
                 "phenoage_noglu_res_sd", # PhenoAge, excluding glucose and HbA1c (sensitivity analysis)
                 "hd_log_sd", # HD based on 18 biomarkers
                 "hd_log_noglu_sd") # HD, excluding glucose and HbA1c (sensitivity analysis)

# Cox models of BA measures on cancers
HR_BA_anycancer <- data.frame(BA=rep(BA_measures,each=2), 
                              outcome="Any cancer", 
                              model=rep(c("Age- and sex-adjusted model","Multivariable  model"), length(BA_measures)),
                              beta=NA, se=NA, HR=NA, lb=NA, ub=NA, p=NA, p_PH_BA=NA)
HR_BA_bc <- data.frame(BA=rep(BA_measures,each=3), 
                       outcome="Breast cancer in women",
                       model=rep(c("Age- and sex-adjusted model","Multivariable  model", "Breast cancer-specific model"), length(BA_measures)),
                       beta=NA, se=NA, HR=NA, lb=NA, ub=NA, p=NA, p_PH_BA=NA)
HR_BA_pc <- data.frame(BA=rep(BA_measures,each=3), 
                       outcome="Prostate cancer in men",
                       model=rep(c("Age- and sex-adjusted model","Multivariable  model", "Prostate cancer-specific model"), length(BA_measures)),
                       beta=NA, se=NA, HR=NA, lb=NA, ub=NA, p=NA, p_PH_BA=NA)
HR_BA_lc <- data.frame(BA=rep(BA_measures,each=3), 
                       outcome="Lung cancer",
                       model=rep(c("Age- and sex-adjusted model","Multivariable  model", "Lung cancer-specific model"), length(BA_measures)),
                       beta=NA, se=NA, HR=NA, lb=NA, ub=NA, p=NA, p_PH_BA=NA)
HR_BA_cc <- data.frame(BA=rep(BA_measures,each=3),
                       outcome="Colorectal cancer", 
                       model=rep(c("Age- and sex-adjusted model","Multivariable  model", "Colorectal cancer-specific model"), length(BA_measures)),
                       beta=NA, se=NA, HR=NA, lb=NA, ub=NA, p=NA, p_PH_BA=NA)
HR_BA_mel <- data.frame(BA=rep(BA_measures,each=3),
                        outcome="Melanoma",
                        model=rep(c("Age- and sex-adjusted model","Multivariable  model", "Melanoma-specific model"), length(BA_measures)),
                        beta=NA, se=NA, HR=NA, lb=NA, ub=NA, p=NA, p_PH_BA=NA)

for (j in 1:length(BA_measures)) {
  
  ###### Any cancer ######
  
  # Age- & sex-adjusted model
  fit <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ . , 
               ukb_bioage_cancer[,c("interview_date","birth_date","censor_date_any_cancer","any_cancer_diag",BA_measures[j],cov_ukb_age_sex)])
  coef <- summary(fit)
  HR_BA_anycancer[j*2-1,c(4,5,9)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_anycancer[j*2-1,6:8] <- coef$conf.int[1,c(1,3,4)]
  HR_BA_anycancer[j*2-1,10] <- cox.zph(fit)$table[1,3] # Test for proportional hazards assumption
  
  # Multivariable model
  fit <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ . , 
               ukb_bioage_cancer[,c("interview_date","birth_date","censor_date_any_cancer","any_cancer_diag",BA_measures[j],cov_ukb_multivariable)])
  coef <- summary(fit)
  HR_BA_anycancer[j*2,c(4,5,9)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_anycancer[j*2,6:8] <- coef$conf.int[1,c(1,3,4)]
  HR_BA_anycancer[j*2,10] <- cox.zph(fit)$table[1,3] # Test for proportional hazards assumption
  
  
  ###### Breast cancer in women ######
  
  # Age- & sex-adjusted model
  fit <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ . , 
               ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women",c("interview_date","birth_date","censor_date_bc","bc_diag",BA_measures[j],cov_ukb_age_sex)])
  coef <- summary(fit)
  HR_BA_bc[j*3-2,c(4,5,9)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_bc[j*3-2,6:8] <- coef$conf.int[1,c(1,3,4)]
  HR_BA_bc[j*3-2,10] <- cox.zph(fit)$table[1,3] # Test for proportional hazards assumption
  
  # Multivariable model
  fit <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ . , 
               ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women",c("interview_date","birth_date","censor_date_bc","bc_diag",BA_measures[j],cov_ukb_multivariable)])
  coef <- summary(fit)
  HR_BA_bc[j*3-1,c(4,5,9)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_bc[j*3-1,6:8] <- coef$conf.int[1,c(1,3,4)]
  HR_BA_bc[j*3-1,10] <- cox.zph(fit)$table[1,3] # Test for proportional hazards assumption
  
  # Breast-cancer specific model
  fit <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ . , 
               ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women",c("interview_date","birth_date","censor_date_bc","bc_diag",BA_measures[j],cov_ukb_bc)])
  coef <- summary(fit)
  HR_BA_bc[j*3,c(4,5,9)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_bc[j*3,6:8] <- coef$conf.int[1,c(1,3,4)]
  HR_BA_bc[j*3,10] <- cox.zph(fit)$table[1,3] # Test for proportional hazards assumption
  
  
  ###### Prostate cancer in men ######
  
  # Age- & sex-adjusted model
  fit <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_pc-birth_date)/365.25, pc_diag) ~ . , 
               ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men",c("interview_date","birth_date","censor_date_pc","pc_diag",BA_measures[j],cov_ukb_age_sex)])
  coef <- summary(fit)
  HR_BA_pc[j*3-2,c(4,5,9)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_pc[j*3-2,6:8] <- coef$conf.int[1,c(1,3,4)]
  HR_BA_pc[j*3-2,10] <- cox.zph(fit)$table[1,3] # Test for proportional hazards assumption
  
  # Multivariable model
  fit <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_pc-birth_date)/365.25, pc_diag) ~ . , 
               ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men",c("interview_date","birth_date","censor_date_pc","pc_diag",BA_measures[j],cov_ukb_multivariable)])
  coef <- summary(fit)
  HR_BA_pc[j*3-1,c(4,5,9)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_pc[j*3-1,6:8] <- coef$conf.int[1,c(1,3,4)]
  HR_BA_pc[j*3-1,10] <- cox.zph(fit)$table[1,3] # Test for proportional hazards assumption
  
  # Breast-cancer specific model
  fit <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_pc-birth_date)/365.25, pc_diag) ~ . , 
               ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men",c("interview_date","birth_date","censor_date_pc","pc_diag",BA_measures[j],cov_ukb_pc)])
  coef <- summary(fit)
  HR_BA_pc[j*3,c(4,5,9)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_pc[j*3,6:8] <- coef$conf.int[1,c(1,3,4)]
  HR_BA_pc[j*3,10] <- cox.zph(fit)$table[1,3] # Test for proportional hazards assumption
  
  
  ###### Lung cancer ###### 
  
  # Age- & sex-adjusted model
  fit <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
               ukb_bioage_cancer[,c("interview_date","birth_date","censor_date_lc","lc_diag",BA_measures[j],cov_ukb_age_sex)])
  coef <- summary(fit)
  HR_BA_lc[j*3-2,c(4,5,9)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_lc[j*3-2,6:8] <- coef$conf.int[1,c(1,3,4)]
  HR_BA_lc[j*3-2,10] <- cox.zph(fit)$table[1,3] # Test for proportional hazards assumption
  
  # Multivariable model
  fit <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
               ukb_bioage_cancer[,c("interview_date","birth_date","censor_date_lc","lc_diag",BA_measures[j],cov_ukb_multivariable)])
  coef <- summary(fit)
  HR_BA_lc[j*3-1,c(4,5,9)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_lc[j*3-1,6:8] <- coef$conf.int[1,c(1,3,4)]
  HR_BA_lc[j*3-1,10] <- cox.zph(fit)$table[1,3] # Test for proportional hazards assumption
  
  # Lung-cancer specific model
  fit <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
               ukb_bioage_cancer[,c("interview_date","birth_date","censor_date_lc","lc_diag",BA_measures[j],cov_ukb_lc)])
  coef <- summary(fit)
  HR_BA_lc[j*3,c(4,5,9)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_lc[j*3,6:8] <- coef$conf.int[1,c(1,3,4)]
  HR_BA_lc[j*3,10] <- cox.zph(fit)$table[1,3] # Test for proportional hazards assumption
  
  
  ###### Colorectal cancer ###### 
  
  # Age- & sex-adjusted model
  fit <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ . , 
               ukb_bioage_cancer[,c("interview_date","birth_date","censor_date_cc","cc_diag",BA_measures[j],cov_ukb_age_sex)])
  coef <- summary(fit)
  HR_BA_cc[j*3-2,c(4,5,9)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_cc[j*3-2,6:8] <- coef$conf.int[1,c(1,3,4)]
  HR_BA_cc[j*3-2,10] <- cox.zph(fit)$table[1,3] # Test for proportional hazards assumption
  
  # Multivariable model
  fit <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ . , 
               ukb_bioage_cancer[,c("interview_date","birth_date","censor_date_cc","cc_diag",BA_measures[j],cov_ukb_multivariable)])
  coef <- summary(fit)
  HR_BA_cc[j*3-1,c(4,5,9)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_cc[j*3-1,6:8] <- coef$conf.int[1,c(1,3,4)]
  HR_BA_cc[j*3-1,10] <- cox.zph(fit)$table[1,3] # Test for proportional hazards assumption
  
  # Colorectal-cancer specific model
  fit <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ . , 
               ukb_bioage_cancer[,c("interview_date","birth_date","censor_date_cc","cc_diag",BA_measures[j],cov_ukb_cc)])
  coef <- summary(fit)
  HR_BA_cc[j*3,c(4,5,9)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_cc[j*3,6:8] <- coef$conf.int[1,c(1,3,4)]
  HR_BA_cc[j*3,10] <- cox.zph(fit)$table[1,3] # Test for proportional hazards assumption
  
  
  ###### Melanoma ###### 
  # Age- & sex-adjusted model
  fit <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ . , 
               ukb_bioage_cancer[ukb_bioage_cancer$sun_protection!="missing"&ukb_bioage_cancer$ethnicity_4group!="missing",c("interview_date","birth_date","censor_date_mel","mel_diag",BA_measures[j],cov_ukb_age_sex)])
  coef <- summary(fit)
  HR_BA_mel[j*3-2,c(4,5,9)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_mel[j*3-2,6:8] <- coef$conf.int[1,c(1,3,4)]
  HR_BA_mel[j*3-2,10] <- cox.zph(fit)$table[1,3] # Test for proportional hazards assumption
  
  # Multivariable model
  fit <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ . , 
               ukb_bioage_cancer[ukb_bioage_cancer$sun_protection!="missing"&ukb_bioage_cancer$ethnicity_4group!="missing",c("interview_date","birth_date","censor_date_mel","mel_diag",BA_measures[j],cov_ukb_multivariable)])
  coef <- summary(fit)
  HR_BA_mel[j*3-1,c(4,5,9)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_mel[j*3-1,6:8] <- coef$conf.int[1,c(1,3,4)]
  HR_BA_mel[j*3-1,10] <- cox.zph(fit)$table[1,3] # Test for proportional hazards assumption
  
  # Melanoma specific model
  fit <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ . , 
               ukb_bioage_cancer[ukb_bioage_cancer$sun_protection!="missing"&ukb_bioage_cancer$ethnicity_4group!="missing",c("interview_date","birth_date","censor_date_mel","mel_diag",BA_measures[j],cov_ukb_mel)])
  coef <- summary(fit)
  HR_BA_mel[j*3,c(4,5,9)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_mel[j*3,6:8] <- coef$conf.int[1,c(1,3,4)]
  HR_BA_mel[j*3,10] <- cox.zph(fit)$table[1,3] # Test for proportional hazards assumption
  
}				
rm(list=c("fit","coef","j"))

### Output results
HR_BA_cancer <- rbind(HR_BA_anycancer,HR_BA_bc,HR_BA_pc,HR_BA_lc,HR_BA_cc,HR_BA_mel)
write.table(HR_BA_cancer, file = "Output/Cox_BA_cancers.txt", sep = "\t", row.names = F, na="")




#============== SUBGROUP ANALYSES FOR BA MEASURES AND CANCERS =================

###### Any cancer ######

### Subgroup analysis, by age

# KDM
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ 
        kdm_res_sd * age60 + birth_year_cat + sex + assessment_centre + bmi_cat + ethnicity_4group + smoking + alcohol_3group + PA_cat + education + di_5group, 
      ukb_bioage_cancer) # p interaction = 0.407010
fit_cox_anycancer_kdm_60minus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ . , 
                                       ukb_bioage_cancer[ukb_bioage_cancer$age60=="Age <60",c("interview_date","birth_date","censor_date_any_cancer","any_cancer_diag","kdm_res_sd",cov_ukb_multivariable)])
fit_cox_anycancer_kdm_60plus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ . , 
                                      ukb_bioage_cancer[ukb_bioage_cancer$age60=="Age 60+",c("interview_date","birth_date","censor_date_any_cancer","any_cancer_diag","kdm_res_sd",cov_ukb_multivariable)])
fit_cox_anycancer_kdm_60minus %>% summary()
fit_cox_anycancer_kdm_60plus %>% summary()

# PhenoAge
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ 
        phenoage_res_sd * age60 + birth_year_cat + sex + assessment_centre + bmi_cat + ethnicity_4group + smoking + alcohol_3group + PA_cat + education + di_5group, 
      ukb_bioage_cancer) # p interaction = 0.11603
fit_cox_anycancer_phenoage_60minus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ . , 
                                            ukb_bioage_cancer[ukb_bioage_cancer$age60=="Age <60",c("interview_date","birth_date","censor_date_any_cancer","any_cancer_diag","phenoage_res_sd",cov_ukb_multivariable)])
fit_cox_anycancer_phenoage_60plus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ . , 
                                           ukb_bioage_cancer[ukb_bioage_cancer$age60=="Age 60+",c("interview_date","birth_date","censor_date_any_cancer","any_cancer_diag","phenoage_res_sd",cov_ukb_multivariable)])
fit_cox_anycancer_phenoage_60minus %>% summary()
fit_cox_anycancer_phenoage_60plus %>% summary()

# HD
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ 
        hd_log_sd * age60 + birth_year_cat + sex + assessment_centre + bmi_cat + ethnicity_4group + smoking + alcohol_3group + PA_cat + education + di_5group, 
      ukb_bioage_cancer) # p interaction = 0.143812
fit_cox_anycancer_hd_60minus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ . , 
                                      ukb_bioage_cancer[ukb_bioage_cancer$age60=="Age <60",c("interview_date","birth_date","censor_date_any_cancer","any_cancer_diag","hd_log_sd",cov_ukb_multivariable)])
fit_cox_anycancer_hd_60plus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ . , 
                                     ukb_bioage_cancer[ukb_bioage_cancer$age60=="Age 60+",c("interview_date","birth_date","censor_date_any_cancer","any_cancer_diag","hd_log_sd",cov_ukb_multivariable)])
fit_cox_anycancer_hd_60minus %>% summary()
fit_cox_anycancer_hd_60plus %>% summary()


### Subgroup analysis, by sex

# KDM
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ 
        kdm_res_sd * sex + birth_year_cat + assessment_centre + bmi_cat + ethnicity_4group + smoking + alcohol_3group + PA_cat + education + di_5group, 
      ukb_bioage_cancer) # p interaction = 0.130674
fit_cox_anycancer_kdm_women <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ . , 
                                     ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women",c("interview_date","birth_date","censor_date_any_cancer","any_cancer_diag","kdm_res_sd",cov_ukb_multivariable)])
fit_cox_anycancer_kdm_men <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ . , 
                                   ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men",c("interview_date","birth_date","censor_date_any_cancer","any_cancer_diag","kdm_res_sd",cov_ukb_multivariable)])
fit_cox_anycancer_kdm_women %>% summary()
fit_cox_anycancer_kdm_men %>% summary()

# PhenoAge
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ 
        phenoage_res_sd * sex + birth_year_cat + assessment_centre + bmi_cat + ethnicity_4group + smoking + alcohol_3group + PA_cat + education + di_5group, 
      ukb_bioage_cancer) # p interaction = 0.01225
fit_cox_anycancer_phenoage_women <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ . , 
                                          ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women",c("interview_date","birth_date","censor_date_any_cancer","any_cancer_diag","phenoage_res_sd",cov_ukb_multivariable)])
fit_cox_anycancer_phenoage_men <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ . , 
                                        ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men",c("interview_date","birth_date","censor_date_any_cancer","any_cancer_diag","phenoage_res_sd",cov_ukb_multivariable)])
fit_cox_anycancer_phenoage_women %>% summary()
fit_cox_anycancer_phenoage_men %>% summary()

# HD
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ 
        hd_log_sd * sex + birth_year_cat + assessment_centre + bmi_cat + ethnicity_4group + smoking + alcohol_3group + PA_cat + education + di_5group, 
      ukb_bioage_cancer) # p interaction = 1.76e-07
fit_cox_anycancer_hd_women <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ . , 
                                    ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women",c("interview_date","birth_date","censor_date_any_cancer","any_cancer_diag","hd_log_sd",cov_ukb_multivariable)])
fit_cox_anycancer_hd_men <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ . , 
                                  ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men",c("interview_date","birth_date","censor_date_any_cancer","any_cancer_diag","hd_log_sd",cov_ukb_multivariable)])
fit_cox_anycancer_hd_women %>% summary()
fit_cox_anycancer_hd_men %>% summary()


### Subgroup analysis, by ethnicity

# KDM
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ 
        kdm_res_sd * ethnicity_2group + birth_year_cat + sex + assessment_centre + bmi_cat + smoking + alcohol_3group + PA_cat + education + di_5group, 
      ukb_bioage_cancer) # p interaction = 0.104831
fit_cox_anycancer_kdm_white <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ . , 
                                     ukb_bioage_cancer[ukb_bioage_cancer$ethnicity_2group=="White",c("interview_date","birth_date","censor_date_any_cancer","any_cancer_diag","kdm_res_sd",cov_ukb_multivariable)])
fit_cox_anycancer_kdm_nonwhite <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ . , 
                                        ukb_bioage_cancer[ukb_bioage_cancer$ethnicity_2group=="Non-white",c("interview_date","birth_date","censor_date_any_cancer","any_cancer_diag","kdm_res_sd",cov_ukb_multivariable)])
fit_cox_anycancer_kdm_white %>% summary()
fit_cox_anycancer_kdm_nonwhite %>% summary()

# PhenoAge
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ 
        phenoage_res_sd * ethnicity_2group + birth_year_cat + sex + assessment_centre + bmi_cat + smoking + alcohol_3group + PA_cat + education + di_5group, 
      ukb_bioage_cancer) # p interaction = 0.114867
fit_cox_anycancer_phenoage_white <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ . , 
                                          ukb_bioage_cancer[ukb_bioage_cancer$ethnicity_2group=="White",c("interview_date","birth_date","censor_date_any_cancer","any_cancer_diag","phenoage_res_sd",cov_ukb_multivariable)])
fit_cox_anycancer_phenoage_nonwhite <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ . , 
                                             ukb_bioage_cancer[ukb_bioage_cancer$ethnicity_2group=="Non-white",c("interview_date","birth_date","censor_date_any_cancer","any_cancer_diag","phenoage_res_sd",cov_ukb_multivariable)])
fit_cox_anycancer_phenoage_white %>% summary()
fit_cox_anycancer_phenoage_nonwhite %>% summary()

# HD
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ 
        hd_log_sd * ethnicity_2group + birth_year_cat + sex + assessment_centre + bmi_cat + smoking + alcohol_3group + PA_cat + education + di_5group, 
      ukb_bioage_cancer) # p interaction = 0.046443
fit_cox_anycancer_hd_white <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ . , 
                                    ukb_bioage_cancer[ukb_bioage_cancer$ethnicity_2group=="White",c("interview_date","birth_date","censor_date_any_cancer","any_cancer_diag","hd_log_sd",cov_ukb_multivariable)])
fit_cox_anycancer_hd_nonwhite <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ . , 
                                       ukb_bioage_cancer[ukb_bioage_cancer$ethnicity_2group=="Non-white",c("interview_date","birth_date","censor_date_any_cancer","any_cancer_diag","hd_log_sd",cov_ukb_multivariable)])
fit_cox_anycancer_hd_white %>% summary()
fit_cox_anycancer_hd_nonwhite %>% summary()



###### Breast cancer in women ######

### Subgroup analysis, by age

# KDM
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ 
        kdm_res_sd*age60+birth_year_cat+assessment_centre+ethnicity_4group+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group+family_history_breast_cancer+bc_screening+menopause+hormone_therapy+oral_contraceptive+parity, 
      ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women",]) # p interaction = 0.168648
fit_cox_bc_kdm_60minus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ . , 
                                ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women" & ukb_bioage_cancer$age60=="Age <60",c("interview_date","birth_date","censor_date_bc","bc_diag","kdm_res_sd",cov_ukb_bc)])
fit_cox_bc_kdm_60plus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ . , 
                               ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women" & ukb_bioage_cancer$age60=="Age 60+",c("interview_date","birth_date","censor_date_bc","bc_diag","kdm_res_sd",cov_ukb_bc)])
fit_cox_bc_kdm_60minus %>% summary()
fit_cox_bc_kdm_60plus %>% summary()

# PhenoAge
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ 
        phenoage_res_sd*age60+birth_year_cat+assessment_centre+ethnicity_4group+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group+family_history_breast_cancer+bc_screening+menopause+hormone_therapy+oral_contraceptive+parity, 
      ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women",]) # p interaction = 0.184078

fit_cox_bc_phenoage_60minus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ . , 
                                     ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women" & ukb_bioage_cancer$age60=="Age <60",c("interview_date","birth_date","censor_date_bc","bc_diag","phenoage_res_sd",cov_ukb_bc)])
fit_cox_bc_phenoage_60plus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ . , 
                                    ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women" & ukb_bioage_cancer$age60=="Age 60+",c("interview_date","birth_date","censor_date_bc","bc_diag","phenoage_res_sd",cov_ukb_bc)])
fit_cox_bc_phenoage_60minus %>% summary()
fit_cox_bc_phenoage_60plus %>% summary()

# HD
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ 
        hd_log_sd*age60+birth_year_cat+assessment_centre+ethnicity_4group+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group+family_history_breast_cancer+bc_screening+menopause+hormone_therapy+oral_contraceptive+parity, 
      ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women",]) # p interaction = 0.963072
fit_cox_bc_hd_60minus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ . , 
                               ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women" & ukb_bioage_cancer$age60=="Age <60",c("interview_date","birth_date","censor_date_bc","bc_diag","hd_log_sd",cov_ukb_bc)])
fit_cox_bc_hd_60plus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ . , 
                              ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women" & ukb_bioage_cancer$age60=="Age 60+",c("interview_date","birth_date","censor_date_bc","bc_diag","hd_log_sd",cov_ukb_bc)])
fit_cox_bc_hd_60minus %>% summary()
fit_cox_bc_hd_60plus %>% summary()


### Subgroup analysis, by menopause

# KDM
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ 
        kdm_res_sd*menopause+birth_year_cat+assessment_centre+ethnicity_4group+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group+family_history_breast_cancer+bc_screening+menopause+hormone_therapy+oral_contraceptive+parity, 
      ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women"& ukb_bioage_cancer$menopause!="missing",]) # p interaction = 0.000923

fit_cox_bc_kdm_premenopause <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ . , 
                                     ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women" & ukb_bioage_cancer$menopause=="Premenopausal",c("interview_date","birth_date","censor_date_bc","bc_diag","kdm_res_sd",cov_ukb_bc[!cov_ukb_bc%in%"menopause"])])
fit_cox_bc_kdm_postmenopause <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ . , 
                                      ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women" & ukb_bioage_cancer$menopause=="Postmenopausal",c("interview_date","birth_date","censor_date_bc","bc_diag","kdm_res_sd",cov_ukb_bc[!cov_ukb_bc%in%"menopause"])])
fit_cox_bc_kdm_premenopause %>% summary()
fit_cox_bc_kdm_postmenopause %>% summary()

# PhenoAge
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ 
        phenoage_res_sd*menopause+birth_year_cat+assessment_centre+ethnicity_4group+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group+family_history_breast_cancer+bc_screening+menopause+hormone_therapy+oral_contraceptive+parity, 
      ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women"& ukb_bioage_cancer$menopause!="missing",]) # p interaction = 0.001881
fit_cox_bc_phenoage_premenopause <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ . , 
                                          ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women" & ukb_bioage_cancer$menopause=="Premenopausal",c("interview_date","birth_date","censor_date_bc","bc_diag","phenoage_res_sd",cov_ukb_bc[!cov_ukb_bc%in%"menopause"])])
fit_cox_bc_phenoage_postmenopause <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ . , 
                                           ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women" & ukb_bioage_cancer$menopause=="Postmenopausal",c("interview_date","birth_date","censor_date_bc","bc_diag","phenoage_res_sd",cov_ukb_bc[!cov_ukb_bc%in%"menopause"])])
fit_cox_bc_phenoage_premenopause %>% summary()
fit_cox_bc_phenoage_postmenopause %>% summary()

# HD
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ 
        hd_log_sd*menopause+birth_year_cat+assessment_centre+ethnicity_4group+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group+family_history_breast_cancer+bc_screening+menopause+hormone_therapy+oral_contraceptive+parity, 
      ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women"& ukb_bioage_cancer$menopause!="missing",]) # p interaction = 0.092332
fit_cox_bc_hd_premenopause <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ . , 
                                    ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women" & ukb_bioage_cancer$menopause=="Premenopausal",c("interview_date","birth_date","censor_date_bc","bc_diag","hd_log_sd",cov_ukb_bc[!cov_ukb_bc%in%"menopause"])])
fit_cox_bc_hd_postmenopause <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ . , 
                                     ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women" & ukb_bioage_cancer$menopause=="Postmenopausal",c("interview_date","birth_date","censor_date_bc","bc_diag","hd_log_sd",cov_ukb_bc[!cov_ukb_bc%in%"menopause"])])
fit_cox_bc_hd_premenopause %>% summary()
fit_cox_bc_hd_postmenopause %>% summary()


### Subgroup analysis, by ethnicity

# KDM
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ 
        kdm_res_sd*ethnicity_2group+birth_year_cat+assessment_centre+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group+family_history_breast_cancer+bc_screening+menopause+hormone_therapy+oral_contraceptive+parity, 
      ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women",]) # p interaction = 0.808441
fit_cox_bc_kdm_white <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ . , 
                              ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women" & ukb_bioage_cancer$ethnicity_2group=="White",c("interview_date","birth_date","censor_date_bc","bc_diag","kdm_res_sd",cov_ukb_bc)])
fit_cox_bc_kdm_nonwhite <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ . , 
                                 ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women" & ukb_bioage_cancer$ethnicity_2group=="Non-white",c("interview_date","birth_date","censor_date_bc","bc_diag","kdm_res_sd",cov_ukb_bc)])
fit_cox_bc_kdm_white %>% summary()
fit_cox_bc_kdm_nonwhite %>% summary()

# PhenoAge
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ 
        phenoage_res_sd*ethnicity_2group+birth_year_cat+assessment_centre+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group+family_history_breast_cancer+bc_screening+menopause+hormone_therapy+oral_contraceptive+parity, 
      ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women",]) # p interaction = 0.317720
fit_cox_bc_phenoage_white <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ . , 
                                   ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women" & ukb_bioage_cancer$ethnicity_2group=="White",c("interview_date","birth_date","censor_date_bc","bc_diag","phenoage_res_sd",cov_ukb_bc)])
fit_cox_bc_phenoage_nonwhite <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ . , 
                                      ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women" & ukb_bioage_cancer$ethnicity_2group=="Non-white",c("interview_date","birth_date","censor_date_bc","bc_diag","phenoage_res_sd",cov_ukb_bc)])
fit_cox_bc_phenoage_white %>% summary()
fit_cox_bc_phenoage_nonwhite %>% summary()

# HD
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ 
        hd_log_sd*ethnicity_2group+birth_year_cat+assessment_centre+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group+family_history_breast_cancer+bc_screening+menopause+hormone_therapy+oral_contraceptive+parity, 
      ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women",]) # p interaction = 0.712201
fit_cox_bc_hd_white <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ . , 
                             ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women" & ukb_bioage_cancer$ethnicity_2group=="White",c("interview_date","birth_date","censor_date_bc","bc_diag","hd_log_sd",cov_ukb_bc)])
fit_cox_bc_hd_nonwhite <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ . , 
                                ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women" & ukb_bioage_cancer$ethnicity_2group=="Non-white",c("interview_date","birth_date","censor_date_bc","bc_diag","hd_log_sd",cov_ukb_bc)])
fit_cox_bc_hd_white %>% summary()
fit_cox_bc_hd_nonwhite %>% summary()



###### Prostate cancer in men ######

### Subgroup analysis, by age

# KDM
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_pc-birth_date)/365.25, pc_diag) ~ 
        kdm_res_sd*age60+birth_year_cat+sex+assessment_centre+bmi_cat+ethnicity_4group+smoking+alcohol_3group+PA_cat+education+di_5group+family_history_prostate_cancer+psa_test+diabetes, 
      ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men",]) # p interaction = 0.39750
fit_cox_pc_kdm_60minus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_pc-birth_date)/365.25, pc_diag) ~ . , 
                                ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men" & ukb_bioage_cancer$age60=="Age <60",c("interview_date","birth_date","censor_date_pc","pc_diag","kdm_res_sd",cov_ukb_pc)])
fit_cox_pc_kdm_60plus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_pc-birth_date)/365.25, pc_diag) ~ . , 
                               ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men" & ukb_bioage_cancer$age60=="Age 60+",c("interview_date","birth_date","censor_date_pc","pc_diag","kdm_res_sd",cov_ukb_pc)])
fit_cox_pc_kdm_60minus %>% summary()
fit_cox_pc_kdm_60plus %>% summary()

# PhenoAge
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_pc-birth_date)/365.25, pc_diag) ~ 
        phenoage_res_sd*age60+birth_year_cat+sex+assessment_centre+bmi_cat+ethnicity_4group+smoking+alcohol_3group+PA_cat+education+di_5group+family_history_prostate_cancer+psa_test+diabetes, 
      ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men",]) # p interaction = 0.20733
fit_cox_phenoage_60minus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_pc-birth_date)/365.25, pc_diag) ~ . , 
                                  ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men" & ukb_bioage_cancer$age60=="Age <60",c("interview_date","birth_date","censor_date_pc","pc_diag","phenoage_res_sd",cov_ukb_pc)])
fit_cox_phenoage_60plus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_pc-birth_date)/365.25, pc_diag) ~ . , 
                                 ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men" & ukb_bioage_cancer$age60=="Age 60+",c("interview_date","birth_date","censor_date_pc","pc_diag","phenoage_res_sd",cov_ukb_pc)])
fit_cox_phenoage_kdm_60minus %>% summary()
fit_cox_phenoage_kdm_60plus %>% summary()

# HD
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_pc-birth_date)/365.25, pc_diag) ~ 
        hd_log_sd*age60+birth_year_cat+sex+assessment_centre+bmi_cat+ethnicity_4group+smoking+alcohol_3group+PA_cat+education+di_5group+family_history_prostate_cancer+psa_test+diabetes, 
      ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men",]) # p interaction = 0.73969 
fit_cox_phenoage_hd_60minus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_pc-birth_date)/365.25, pc_diag) ~ . , 
                                     ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men" & ukb_bioage_cancer$age60=="Age <60",c("interview_date","birth_date","censor_date_pc","pc_diag","hd_log_sd",cov_ukb_pc)])
fit_cox_phenoage_hd_60plus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_pc-birth_date)/365.25, pc_diag) ~ . , 
                                    ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men" & ukb_bioage_cancer$age60=="Age 60+",c("interview_date","birth_date","censor_date_pc","pc_diag","hd_log_sd",cov_ukb_pc)])
fit_cox_phenoage_hd_60minus %>% summary()
fit_cox_phenoage_hd_60plus %>% summary()


### Subgroup analysis, by ethnicity

# KDM
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_pc-birth_date)/365.25, pc_diag) ~ 
        kdm_res_sd*ethnicity_2group+birth_year_cat+sex+assessment_centre+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group+family_history_prostate_cancer+psa_test+diabetes, 
      ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men",]) # p interaction = 0.069912
fit_cox_pc_kdm_white <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_pc-birth_date)/365.25, pc_diag) ~ . , 
                              ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men" & ukb_bioage_cancer$ethnicity_2group=="White",c("interview_date","birth_date","censor_date_pc","pc_diag","kdm_res_sd",cov_ukb_pc)])
fit_cox_pc_kdm_nonwhite <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_pc-birth_date)/365.25, pc_diag) ~ . , 
                                 ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men" & ukb_bioage_cancer$ethnicity_2group=="Non-white",c("interview_date","birth_date","censor_date_pc","pc_diag","kdm_res_sd",cov_ukb_pc)])
fit_cox_pc_kdm_white %>% summary()
fit_cox_pc_kdm_nonwhite %>% summary()

# PhenoAge
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_pc-birth_date)/365.25, pc_diag) ~ 
        phenoage_res_sd*ethnicity_2group+birth_year_cat+sex+assessment_centre+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group+family_history_prostate_cancer+psa_test+diabetes, 
      ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men",]) # p interaction = 0.000336
fit_cox_phenoage_white <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_pc-birth_date)/365.25, pc_diag) ~ . , 
                                ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men" & ukb_bioage_cancer$ethnicity_2group=="White",c("interview_date","birth_date","censor_date_pc","pc_diag","phenoage_res_sd",cov_ukb_pc)])
fit_cox_phenoage_nonwhite <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_pc-birth_date)/365.25, pc_diag) ~ . , 
                                   ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men" & ukb_bioage_cancer$ethnicity_2group=="Non-white",c("interview_date","birth_date","censor_date_pc","pc_diag","phenoage_res_sd",cov_ukb_pc)])
fit_cox_phenoage_white %>% summary()
fit_cox_phenoage_nonwhite %>% summary()

# HD
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_pc-birth_date)/365.25, pc_diag) ~ 
        hd_log_sd*ethnicity_2group+birth_year_cat+sex+assessment_centre+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group+family_history_prostate_cancer+psa_test+diabetes, 
      ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men",]) # p interaction = 0.20413 
fit_cox_phenoage_hd_white <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_pc-birth_date)/365.25, pc_diag) ~ . , 
                                   ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men" & ukb_bioage_cancer$ethnicity_2group=="White",c("interview_date","birth_date","censor_date_pc","pc_diag","hd_log_sd",cov_ukb_pc)])
fit_cox_phenoage_hd_nonwhite <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_pc-birth_date)/365.25, pc_diag) ~ . , 
                                      ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men" & ukb_bioage_cancer$ethnicity_2group=="Non-white",c("interview_date","birth_date","censor_date_pc","pc_diag","hd_log_sd",cov_ukb_pc)])
fit_cox_phenoage_hd_white %>% summary()
fit_cox_phenoage_hd_nonwhite %>% summary()



###### Lung cancer ######

### Subgroup analysis, by age

# KDM
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ 
        kdm_res_sd * age60 + birth_year_cat + sex + assessment_centre + bmi_cat + ethnicity_4group + smoking + alcohol_3group + PA_cat + education + di_5group + family_history_lung_cancer, 
      ukb_bioage_cancer) # p interaction = 0.883868
fit_cox_lc_kdm_60minus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
                                ukb_bioage_cancer[ukb_bioage_cancer$age60=="Age <60",c("interview_date","birth_date","censor_date_lc","lc_diag","kdm_res_sd",cov_ukb_lc)])
fit_cox_lc_kdm_60plus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
                               ukb_bioage_cancer[ukb_bioage_cancer$age60=="Age 60+",c("interview_date","birth_date","censor_date_lc","lc_diag","kdm_res_sd",cov_ukb_lc)])
fit_cox_lc_kdm_60minus %>% summary()
fit_cox_lc_kdm_60plus %>% summary()

# PhenoAge
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ 
        phenoage_res_sd * age60 + birth_year_cat + sex + assessment_centre + bmi_cat + ethnicity_4group + smoking + alcohol_3group + PA_cat + education + di_5group + family_history_lung_cancer, 
      ukb_bioage_cancer) # p interaction = 0.637857
fit_cox_lc_phenoage_60minus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
                                     ukb_bioage_cancer[ukb_bioage_cancer$age60=="Age <60",c("interview_date","birth_date","censor_date_lc","lc_diag","phenoage_res_sd",cov_ukb_lc)])
fit_cox_lc_phenoage_60plus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
                                    ukb_bioage_cancer[ukb_bioage_cancer$age60=="Age 60+",c("interview_date","birth_date","censor_date_lc","lc_diag","phenoage_res_sd",cov_ukb_lc)])
fit_cox_lc_phenoage_60minus %>% summary()
fit_cox_lc_phenoage_60plus %>% summary()

# HD
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ 
        hd_log_sd * age60 + birth_year_cat + sex + assessment_centre + bmi_cat + ethnicity_4group + smoking + alcohol_3group + PA_cat + education + di_5group + family_history_lung_cancer, 
      ukb_bioage_cancer) # p interaction = 0.941558
fit_cox_lc_hd_60minus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
                               ukb_bioage_cancer[ukb_bioage_cancer$age60=="Age <60",c("interview_date","birth_date","censor_date_lc","lc_diag","hd_log_sd",cov_ukb_lc)])
fit_cox_lc_hd_60plus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
                              ukb_bioage_cancer[ukb_bioage_cancer$age60=="Age 60+",c("interview_date","birth_date","censor_date_lc","lc_diag","hd_log_sd",cov_ukb_lc)])
fit_cox_lc_hd_60minus %>% summary()
fit_cox_lc_hd_60plus %>% summary()


### Subgroup analysis, by sex

# KDM
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ 
        kdm_res_sd * sex + birth_year_cat + assessment_centre + bmi_cat + ethnicity_4group + smoking + alcohol_3group + PA_cat + education + di_5group + family_history_lung_cancer, 
      ukb_bioage_cancer) # p interaction = 0.050391
fit_cox_lc_kdm_women <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
                              ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women",c("interview_date","birth_date","censor_date_lc","lc_diag","kdm_res_sd",cov_ukb_lc)])
fit_cox_lc_kdm_men <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
                            ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men",c("interview_date","birth_date","censor_date_lc","lc_diag","kdm_res_sd",cov_ukb_lc)])
fit_cox_lc_kdm_women %>% summary()
fit_cox_lc_kdm_men %>% summary()

# PhenoAge
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ 
        phenoage_res_sd * sex + birth_year_cat + assessment_centre + bmi_cat + ethnicity_4group + smoking + alcohol_3group + PA_cat + education + di_5group + family_history_lung_cancer, 
      ukb_bioage_cancer) # p interaction = 0.290766
fit_cox_lc_phenoage_women <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
                                   ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women",c("interview_date","birth_date","censor_date_lc","lc_diag","phenoage_res_sd",cov_ukb_lc)])
fit_cox_lc_phenoage_men <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
                                 ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men",c("interview_date","birth_date","censor_date_lc","lc_diag","phenoage_res_sd",cov_ukb_lc)])
fit_cox_lc_phenoage_women %>% summary()
fit_cox_lc_phenoage_men %>% summary()

# HD
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ 
        hd_log_sd * sex + birth_year_cat + assessment_centre + bmi_cat + ethnicity_4group + smoking + alcohol_3group + PA_cat + education + di_5group + family_history_lung_cancer, 
      ukb_bioage_cancer) # p interaction = 0.812911
fit_cox_lc_hd_women <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
                             ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women",c("interview_date","birth_date","censor_date_lc","lc_diag","hd_log_sd",cov_ukb_lc)])
fit_cox_lc_hd_men <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
                           ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men",c("interview_date","birth_date","censor_date_lc","lc_diag","hd_log_sd",cov_ukb_lc)])
fit_cox_lc_hd_women %>% summary()
fit_cox_lc_hd_men %>% summary()


### Subgroup analysis, by smoking

# KDM
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ 
        kdm_res_sd * eversmoke + birth_year_cat + assessment_centre + bmi_cat + ethnicity_4group + alcohol_3group + PA_cat + education + di_5group + family_history_lung_cancer, 
      ukb_bioage_cancer %>% mutate(eversmoke=case_when(smoking=="Never"~0,smoking=="Previous"|smoking=="Current"~1))) # p interaction = 2.64e-11
fit_cox_lc_kdm_nonsmoker <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
                                  ukb_bioage_cancer[ukb_bioage_cancer$smoking=="Never",c("interview_date","birth_date","censor_date_lc","lc_diag","kdm_res_sd",cov_ukb_lc[!cov_ukb_lc%in%"smoking"])])
fit_cox_lc_kdm_smoker <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
                               ukb_bioage_cancer[ukb_bioage_cancer$smoking=="Previous"|ukb_bioage_cancer$smoking=="Current",c("interview_date","birth_date","censor_date_lc","lc_diag","kdm_res_sd",cov_ukb_lc[!cov_ukb_lc%in%"smoking"])])
fit_cox_lc_kdm_nonsmoker %>% summary()
fit_cox_lc_kdm_smoker %>% summary()

# PhenoAge
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ 
        phenoage_res_sd * eversmoke + birth_year_cat + assessment_centre + bmi_cat + ethnicity_4group + alcohol_3group + PA_cat + education + di_5group + family_history_lung_cancer, 
      ukb_bioage_cancer %>% mutate(eversmoke=case_when(smoking=="Never"~0,smoking=="Previous"|smoking=="Current"~1))) # p interaction = 3.70e-13
fit_cox_lc_phenoage_nonsmoker <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
                                       ukb_bioage_cancer[ukb_bioage_cancer$smoking=="Never",c("interview_date","birth_date","censor_date_lc","lc_diag","phenoage_res_sd",cov_ukb_lc[!cov_ukb_lc%in%"smoking"])])
fit_cox_lc_phenoage_smoker <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
                                    ukb_bioage_cancer[ukb_bioage_cancer$smoking=="Previous"|ukb_bioage_cancer$smoking=="Current",c("interview_date","birth_date","censor_date_lc","lc_diag","phenoage_res_sd",cov_ukb_lc[!cov_ukb_lc%in%"smoking"])])
fit_cox_lc_phenoage_nonsmoker %>% summary()
fit_cox_lc_phenoage_smoker %>% summary()

# HD
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ 
        hd_log_sd * eversmoke + birth_year_cat + assessment_centre + bmi_cat + ethnicity_4group + alcohol_3group + PA_cat + education + di_5group + family_history_lung_cancer, 
      ukb_bioage_cancer %>% mutate(eversmoke=case_when(smoking=="Never"~0,smoking=="Previous"|smoking=="Current"~1))) # p interaction = 7.29e-05
fit_cox_lc_hd_nonsmoker <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
                                 ukb_bioage_cancer[ukb_bioage_cancer$smoking=="Never",c("interview_date","birth_date","censor_date_lc","lc_diag","hd_log_sd",cov_ukb_lc[!cov_ukb_lc%in%"smoking"])])
fit_cox_lc_hd_smoker <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
                              ukb_bioage_cancer[ukb_bioage_cancer$smoking=="Previous"|ukb_bioage_cancer$smoking=="Current",c("interview_date","birth_date","censor_date_lc","lc_diag","hd_log_sd",cov_ukb_lc[!cov_ukb_lc%in%"smoking"])])
fit_cox_lc_hd_nonsmoker %>% summary()
fit_cox_lc_hd_smoker %>% summary()


### Subgroup analysis, by ethnicity

# KDM
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ 
        kdm_res_sd * ethnicity_2group + birth_year_cat + sex + assessment_centre + bmi_cat + smoking + alcohol_3group + PA_cat + education + di_5group + family_history_lung_cancer, 
      ukb_bioage_cancer) # p interaction = 0.001053
fit_cox_lc_kdm_white <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
                              ukb_bioage_cancer[ukb_bioage_cancer$ethnicity_2group=="White",c("interview_date","birth_date","censor_date_lc","lc_diag","kdm_res_sd",cov_ukb_lc)])
fit_cox_lc_kdm_nonwhite <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
                                 ukb_bioage_cancer[ukb_bioage_cancer$ethnicity_2group=="Non-white",c("interview_date","birth_date","censor_date_lc","lc_diag","kdm_res_sd",cov_ukb_lc)])
fit_cox_lc_kdm_white %>% summary()
fit_cox_lc_kdm_nonwhite %>% summary()

# PhenoAge
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ 
        phenoage_res_sd * ethnicity_2group + birth_year_cat + sex + assessment_centre + bmi_cat + smoking + alcohol_3group + PA_cat + education + di_5group + family_history_lung_cancer, 
      ukb_bioage_cancer) # p interaction = 0.002584
fit_cox_lc_phenoage_white <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
                                   ukb_bioage_cancer[ukb_bioage_cancer$ethnicity_2group=="White",c("interview_date","birth_date","censor_date_lc","lc_diag","phenoage_res_sd",cov_ukb_lc)])
fit_cox_lc_phenoage_nonwhite <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
                                      ukb_bioage_cancer[ukb_bioage_cancer$ethnicity_2group=="Non-white",c("interview_date","birth_date","censor_date_lc","lc_diag","phenoage_res_sd",cov_ukb_lc)])
fit_cox_lc_phenoage_white %>% summary()
fit_cox_lc_phenoage_nonwhite %>% summary()

# HD
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ 
        hd_log_sd * ethnicity_2group + birth_year_cat + sex + assessment_centre + bmi_cat + smoking + alcohol_3group + PA_cat + education + di_5group + family_history_lung_cancer, 
      ukb_bioage_cancer) # p interaction = 0.010988
fit_cox_lc_hd_white <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
                             ukb_bioage_cancer[ukb_bioage_cancer$ethnicity_2group=="White",c("interview_date","birth_date","censor_date_lc","lc_diag","hd_log_sd",cov_ukb_lc)])
fit_cox_lc_hd_nonwhite <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
                                ukb_bioage_cancer[ukb_bioage_cancer$ethnicity_2group=="Non-white",c("interview_date","birth_date","censor_date_lc","lc_diag","hd_log_sd",cov_ukb_lc)])
fit_cox_lc_hd_white %>% summary()
fit_cox_lc_hd_nonwhite %>% summary()


###### Colorectal cancer ######

### Subgroup analysis, by age

# KDM
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ 
        kdm_res_sd * age60 + birth_year_cat + sex + assessment_centre + bmi_cat + ethnicity_4group + smoking + alcohol_3group + PA_cat + education +di_5group+family_history_bowel_cancer+cc_screening+fresh_veg_cat+red_meat_cat+processed_meat_cat, 
      ukb_bioage_cancer) # p interaction = 0.109533
fit_cox_cc_kdm_60minus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ . , 
                                ukb_bioage_cancer[ukb_bioage_cancer$age60=="Age <60",c("interview_date","birth_date","censor_date_cc","cc_diag","kdm_res_sd",cov_ukb_cc)])
fit_cox_cc_kdm_60plus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ . , 
                               ukb_bioage_cancer[ukb_bioage_cancer$age60=="Age 60+",c("interview_date","birth_date","censor_date_cc","cc_diag","kdm_res_sd",cov_ukb_cc)])
fit_cox_cc_kdm_60minus %>% summary()
fit_cox_cc_kdm_60plus %>% summary()

# PhenoAge
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ 
        phenoage_res_sd * age60 + birth_year_cat + sex + assessment_centre + bmi_cat + ethnicity_4group + smoking + alcohol_3group + PA_cat + education +di_5group+family_history_bowel_cancer+cc_screening+fresh_veg_cat+red_meat_cat+processed_meat_cat, 
      ukb_bioage_cancer) # p interaction = 0.472748
fit_cox_cc_phenoage_60minus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ . , 
                                     ukb_bioage_cancer[ukb_bioage_cancer$age60=="Age <60",c("interview_date","birth_date","censor_date_cc","cc_diag","phenoage_res_sd",cov_ukb_cc)])
fit_cox_cc_phenoage_60plus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ . , 
                                    ukb_bioage_cancer[ukb_bioage_cancer$age60=="Age 60+",c("interview_date","birth_date","censor_date_cc","cc_diag","phenoage_res_sd",cov_ukb_cc)])
fit_cox_cc_phenoage_60minus %>% summary()
fit_cox_cc_phenoage_60plus %>% summary()

# HD
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ 
        hd_log_sd * age60 + birth_year_cat + sex + assessment_centre + bmi_cat + ethnicity_4group + smoking + alcohol_3group + PA_cat + education +di_5group+family_history_bowel_cancer+cc_screening+fresh_veg_cat+red_meat_cat+processed_meat_cat, 
      ukb_bioage_cancer) # p interaction = 0.292623
fit_cox_cc_hd_60minus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ . , 
                               ukb_bioage_cancer[ukb_bioage_cancer$age60=="Age <60",c("interview_date","birth_date","censor_date_cc","cc_diag","hd_log_sd",cov_ukb_cc)])
fit_cox_cc_hd_60plus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ . , 
                              ukb_bioage_cancer[ukb_bioage_cancer$age60=="Age 60+",c("interview_date","birth_date","censor_date_cc","cc_diag","hd_log_sd",cov_ukb_cc)])
fit_cox_cc_hd_60minus %>% summary()
fit_cox_cc_hd_60plus %>% summary()


### Subgroup analysis, by sex

# KDM
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ 
        kdm_res_sd * sex + birth_year_cat + assessment_centre + bmi_cat + ethnicity_4group + smoking + alcohol_3group + PA_cat + education +di_5group+family_history_bowel_cancer+cc_screening+fresh_veg_cat+red_meat_cat+processed_meat_cat, 
      ukb_bioage_cancer) # p interaction = 0.001386
fit_cox_cc_kdm_women <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ . , 
                              ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women",c("interview_date","birth_date","censor_date_cc","cc_diag","kdm_res_sd",cov_ukb_cc)])
fit_cox_cc_kdm_men <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ . , 
                            ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men",c("interview_date","birth_date","censor_date_cc","cc_diag","kdm_res_sd",cov_ukb_cc)])
fit_cox_cc_kdm_women %>% summary()
fit_cox_cc_kdm_men %>% summary()

# PhenoAge
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ 
        phenoage_res_sd * sex + birth_year_cat + assessment_centre + bmi_cat + ethnicity_4group + smoking + alcohol_3group + PA_cat + education +di_5group+family_history_bowel_cancer+cc_screening+fresh_veg_cat+red_meat_cat+processed_meat_cat, 
      ukb_bioage_cancer) # p interaction = 0.03179
fit_cox_cc_phenoage_women <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ . , 
                                   ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women",c("interview_date","birth_date","censor_date_cc","cc_diag","phenoage_res_sd",cov_ukb_cc)])
fit_cox_cc_phenoage_men <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ . , 
                                 ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men",c("interview_date","birth_date","censor_date_cc","cc_diag","phenoage_res_sd",cov_ukb_cc)])
fit_cox_cc_phenoage_women %>% summary()
fit_cox_cc_phenoage_men %>% summary()

# HD
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ 
        hd_log_sd * sex + birth_year_cat + assessment_centre + bmi_cat + ethnicity_4group + smoking + alcohol_3group + PA_cat + education +di_5group+family_history_bowel_cancer+cc_screening+fresh_veg_cat+red_meat_cat+processed_meat_cat, 
      ukb_bioage_cancer) # p interaction = 0.394094
fit_cox_cc_hd_women <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ . , 
                             ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women",c("interview_date","birth_date","censor_date_cc","cc_diag","hd_log_sd",cov_ukb_cc)])
fit_cox_cc_hd_men <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ . , 
                           ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men",c("interview_date","birth_date","censor_date_cc","cc_diag","hd_log_sd",cov_ukb_cc)])
fit_cox_cc_hd_women %>% summary()
fit_cox_cc_hd_men %>% summary()


### Subgroup analysis, by ethnicity

# KDM
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ 
        kdm_res_sd * ethnicity_2group + birth_year_cat + sex + assessment_centre + bmi_cat + smoking + alcohol_3group + PA_cat + education +di_5group+family_history_bowel_cancer+cc_screening+fresh_veg_cat+red_meat_cat+processed_meat_cat, 
      ukb_bioage_cancer) # p interaction = 0.480912
fit_cox_cc_kdm_white <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ . , 
                              ukb_bioage_cancer[ukb_bioage_cancer$ethnicity_2group=="White",c("interview_date","birth_date","censor_date_cc","cc_diag","kdm_res_sd",cov_ukb_cc)])
fit_cox_cc_kdm_nonwhite <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ . , 
                                 ukb_bioage_cancer[ukb_bioage_cancer$ethnicity_2group=="Non-white",c("interview_date","birth_date","censor_date_cc","cc_diag","kdm_res_sd",cov_ukb_cc)])
fit_cox_cc_kdm_white %>% summary()
fit_cox_cc_kdm_nonwhite %>% summary()

# PhenoAge
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ 
        phenoage_res_sd * ethnicity_2group + birth_year_cat + sex + assessment_centre + bmi_cat + smoking + alcohol_3group + PA_cat + education +di_5group+family_history_bowel_cancer+cc_screening+fresh_veg_cat+red_meat_cat+processed_meat_cat, 
      ukb_bioage_cancer) # p interaction = 0.708330
fit_cox_cc_phenoage_white <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ . , 
                                   ukb_bioage_cancer[ukb_bioage_cancer$ethnicity_2group=="White",c("interview_date","birth_date","censor_date_cc","cc_diag","phenoage_res_sd",cov_ukb_cc)])
fit_cox_cc_phenoage_nonwhite <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ . , 
                                      ukb_bioage_cancer[ukb_bioage_cancer$ethnicity_2group=="Non-white",c("interview_date","birth_date","censor_date_cc","cc_diag","phenoage_res_sd",cov_ukb_cc)])
fit_cox_cc_phenoage_white %>% summary()
fit_cox_cc_phenoage_nonwhite %>% summary()

# HD
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ 
        hd_log_sd * ethnicity_2group + birth_year_cat + sex + assessment_centre + bmi_cat + smoking + alcohol_3group + PA_cat + education +di_5group+family_history_bowel_cancer+cc_screening+fresh_veg_cat+red_meat_cat+processed_meat_cat, 
      ukb_bioage_cancer) # p interaction = 0.617469
fit_cox_cc_hd_white <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ . , 
                             ukb_bioage_cancer[ukb_bioage_cancer$ethnicity_2group=="White",c("interview_date","birth_date","censor_date_cc","cc_diag","hd_log_sd",cov_ukb_cc)])
fit_cox_cc_hd_nonwhite <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ . , 
                                ukb_bioage_cancer[ukb_bioage_cancer$ethnicity_2group=="Non-white",c("interview_date","birth_date","censor_date_cc","cc_diag","hd_log_sd",cov_ukb_cc)])
fit_cox_cc_hd_white %>% summary()
fit_cox_cc_hd_nonwhite %>% summary()


###### Melanoma ######

### Subgroup analysis, by age

# KDM
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ 
        kdm_res_sd * age60 + birth_year_cat + sex + assessment_centre + bmi_cat + ethnicity_4group + smoking + alcohol_3group + PA_cat + education + di_5group + 
        PA_cat+time_outdoors_summer_4group+childhood_sunburn+solarium_use+skin_tan_ease+skin_color_4group+hair_color_3group+sun_protection, 
      ukb_bioage_cancer) # p interaction = 0.000892
fit_cox_mel_kdm_60minus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ . , 
                                 ukb_bioage_cancer[ukb_bioage_cancer$age60=="Age <60",c("interview_date","birth_date","censor_date_mel","mel_diag","kdm_res_sd",cov_ukb_mel)])
fit_cox_mel_kdm_60plus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ . , 
                                ukb_bioage_cancer[ukb_bioage_cancer$age60=="Age 60+",c("interview_date","birth_date","censor_date_mel","mel_diag","kdm_res_sd",cov_ukb_mel)])
fit_cox_mel_kdm_60minus %>% summary()
fit_cox_mel_kdm_60plus %>% summary()

# PhenoAge
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ 
        phenoage_res_sd * age60 + birth_year_cat + sex + assessment_centre + bmi_cat + ethnicity_4group + smoking + alcohol_3group + PA_cat + education + di_5group + 
        PA_cat+time_outdoors_summer_4group+childhood_sunburn+solarium_use+skin_tan_ease+skin_color_4group+hair_color_3group+sun_protection, 
      ukb_bioage_cancer) # p interaction = 0.00516
fit_cox_mel_phenoage_60minus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ . , 
                                      ukb_bioage_cancer[ukb_bioage_cancer$age60=="Age <60",c("interview_date","birth_date","censor_date_mel","mel_diag","phenoage_res_sd",cov_ukb_mel)])
fit_cox_mel_phenoage_60plus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ . , 
                                     ukb_bioage_cancer[ukb_bioage_cancer$age60=="Age 60+",c("interview_date","birth_date","censor_date_mel","mel_diag","phenoage_res_sd",cov_ukb_mel)])
fit_cox_mel_phenoage_60minus %>% summary()
fit_cox_mel_phenoage_60plus %>% summary()

# HD
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ 
        hd_log_sd * age60 + birth_year_cat + sex + assessment_centre + bmi_cat + ethnicity_4group + smoking + alcohol_3group + PA_cat + education + di_5group + 
        PA_cat+time_outdoors_summer_4group+childhood_sunburn+solarium_use+skin_tan_ease+skin_color_4group+hair_color_3group+sun_protection, 
      ukb_bioage_cancer) # p interaction = 0.172445
fit_cox_mel_hd_60minus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ . , 
                                ukb_bioage_cancer[ukb_bioage_cancer$age60=="Age <60",c("interview_date","birth_date","censor_date_mel","mel_diag","hd_log_sd",cov_ukb_mel)])
fit_cox_mel_hd_60plus <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ . , 
                               ukb_bioage_cancer[ukb_bioage_cancer$age60=="Age 60+",c("interview_date","birth_date","censor_date_mel","mel_diag","hd_log_sd",cov_ukb_mel)])
fit_cox_mel_hd_60minus %>% summary()
fit_cox_mel_hd_60plus %>% summary()


### Subgroup analysis, by sex

# KDM
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ 
        kdm_res_sd * sex + birth_year_cat + assessment_centre + bmi_cat + ethnicity_4group + smoking + alcohol_3group + PA_cat + education + di_5group + 
        PA_cat+time_outdoors_summer_4group+childhood_sunburn+solarium_use+skin_tan_ease+skin_color_4group+hair_color_3group+sun_protection, 
      ukb_bioage_cancer) # p interaction = 0.09018
fit_cox_mel_kdm_women <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ . , 
                               ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women",c("interview_date","birth_date","censor_date_mel","mel_diag","kdm_res_sd",cov_ukb_mel)])
fit_cox_mel_kdm_men <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ . , 
                             ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men",c("interview_date","birth_date","censor_date_mel","mel_diag","kdm_res_sd",cov_ukb_mel)])
fit_cox_mel_kdm_women %>% summary()
fit_cox_mel_kdm_men %>% summary()


# PhenoAge
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ 
        phenoage_res_sd * sex + birth_year_cat + assessment_centre + bmi_cat + ethnicity_4group + smoking + alcohol_3group + PA_cat + education + di_5group + 
        PA_cat+time_outdoors_summer_4group+childhood_sunburn+solarium_use+skin_tan_ease+skin_color_4group+hair_color_3group+sun_protection, 
      ukb_bioage_cancer) # p interaction = 0.81496
fit_cox_mel_phenoage_women <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ . , 
                                    ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women",c("interview_date","birth_date","censor_date_mel","mel_diag","phenoage_res_sd",cov_ukb_mel)])
fit_cox_mel_phenoage_men <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ . , 
                                  ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men",c("interview_date","birth_date","censor_date_mel","mel_diag","phenoage_res_sd",cov_ukb_mel)])
fit_cox_mel_phenoage_women %>% summary()
fit_cox_mel_phenoage_men %>% summary()

# HD
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ 
        hd_log_sd * sex + birth_year_cat + assessment_centre + bmi_cat + ethnicity_4group + smoking + alcohol_3group + PA_cat + education + di_5group + 
        PA_cat+time_outdoors_summer_4group+childhood_sunburn+solarium_use+skin_tan_ease+skin_color_4group+hair_color_3group+sun_protection, 
      ukb_bioage_cancer) # p interaction = 0.01721
fit_cox_mel_hd_women <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ . , 
                              ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women",c("interview_date","birth_date","censor_date_mel","mel_diag","hd_log_sd",cov_ukb_mel)])
fit_cox_mel_hd_men <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ . , 
                            ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men",c("interview_date","birth_date","censor_date_mel","mel_diag","hd_log_sd",cov_ukb_mel)])
fit_cox_mel_hd_women %>% summary()
fit_cox_mel_hd_men %>% summary()


### Subgroup analysis, by ethnicity

# KDM
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ 
        kdm_res_sd * ethnicity_2group + birth_year_cat + sex + assessment_centre + bmi_cat + smoking + alcohol_3group + PA_cat + education + di_5group + 
        PA_cat+time_outdoors_summer_4group+childhood_sunburn+solarium_use+skin_tan_ease+skin_color_4group+hair_color_3group+sun_protection, 
      ukb_bioage_cancer) # p interaction = 0.902966
fit_cox_mel_kdm_white <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ . , 
                               ukb_bioage_cancer[ukb_bioage_cancer$ethnicity_2group=="White",c("interview_date","birth_date","censor_date_mel","mel_diag","kdm_res_sd",cov_ukb_mel)])
fit_cox_mel_kdm_nonwhite <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ . , 
                                  ukb_bioage_cancer[ukb_bioage_cancer$ethnicity_2group=="Non-white",c("interview_date","birth_date","censor_date_mel","mel_diag","kdm_res_sd",cov_ukb_mel)])
fit_cox_mel_kdm_white %>% summary()
fit_cox_mel_kdm_nonwhite %>% summary()

# PhenoAge
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ 
        phenoage_res_sd * ethnicity_2group + birth_year_cat + sex + assessment_centre + bmi_cat + smoking + alcohol_3group + PA_cat + education + di_5group + 
        PA_cat+time_outdoors_summer_4group+childhood_sunburn+solarium_use+skin_tan_ease+skin_color_4group+hair_color_3group+sun_protection, 
      ukb_bioage_cancer) # p interaction = 0.99071
fit_cox_mel_phenoage_white <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ . , 
                                    ukb_bioage_cancer[ukb_bioage_cancer$ethnicity_2group=="White",c("interview_date","birth_date","censor_date_mel","mel_diag","phenoage_res_sd",cov_ukb_mel)])
fit_cox_mel_phenoage_nonwhite <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ . , 
                                       ukb_bioage_cancer[ukb_bioage_cancer$ethnicity_2group=="Non-white",c("interview_date","birth_date","censor_date_mel","mel_diag","phenoage_res_sd",cov_ukb_mel)])
fit_cox_mel_phenoage_white %>% summary()
fit_cox_mel_phenoage_nonwhite %>% summary()

# HD
coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ 
        hd_log_sd * ethnicity_2group + birth_year_cat + sex + assessment_centre + bmi_cat + smoking + alcohol_3group + PA_cat + education + di_5group + 
        PA_cat+time_outdoors_summer_4group+childhood_sunburn+solarium_use+skin_tan_ease+skin_color_4group+hair_color_3group+sun_protection, 
      ukb_bioage_cancer) # p interaction = 0.56681
fit_cox_mel_hd_white <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ . , 
                              ukb_bioage_cancer[ukb_bioage_cancer$ethnicity_2group=="White",c("interview_date","birth_date","censor_date_mel","mel_diag","hd_log_sd",cov_ukb_mel)])
fit_cox_mel_hd_nonwhite <- coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ . , 
                                 ukb_bioage_cancer[ukb_bioage_cancer$ethnicity_2group=="Non-white",c("interview_date","birth_date","censor_date_mel","mel_diag","hd_log_sd",cov_ukb_mel)])
fit_cox_mel_hd_white %>% summary()
fit_cox_mel_hd_nonwhite %>% summary()



#======================== TIME-VARYING HAZARD RATIOS ==========================

# BA measures to be tested
BA_measures <- c("kdm_res_sd", # KDM based on 18 biomarkers
                 "kdm_original_res_sd", # Original Levine KDM
                 "kdm_noglu_res_sd", # KDM, excluding glucose and HbA1c (sensitivity analysis)
                 "phenoage_res_sd", # PhenoAge based on 18 biomarkers
                 "phenoage_original_res_sd", # Original Levine ?PhenoAge
                 "phenoage_noglu_res_sd", # PhenoAge, excluding glucose and HbA1c (sensitivity analysis)
                 "hd_log_sd", # HD based on 18 biomarkers
                 "hd_log_noglu_sd") # HD, excluding glucose and HbA1c (sensitivity analysis)

# Test for proportional hazards
#HR_BA_cancer <- read.delim("Output/Cox_BA_cancers.txt")
nonPH <- HR_BA_cancer %>% mutate(nonPH = p_PH_BA< .05) %>% # Non-proportional hazards if p<.05
  filter(BA %in% BA_measures[c(1,4,7)] & 
           (outcome=="Any cancer" & model=="Multivariable  model" |
              outcome=="Breast cancer in women" & model=="Breast cancer-specific model" |
              outcome=="Prostate cancer in men" & model=="Prostate cancer-specific model" |
              outcome=="Lung cancer" & model=="Lung cancer-specific model" |
              outcome=="Colorectal cancer" & model=="Colorectal cancer-specific model" |
              outcome=="Melanoma" & model=="Melanoma-specific model" ) &
           nonPH)
nonPH # There is evidence of non-proportional-hazards in 7 models


####### Time-varying HR for PhenoAge and any cancer #######

dat_timesplit <- survSplit(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ . , 
                           ukb_bioage_cancer[,c("interview_date","birth_date","censor_date_any_cancer","any_cancer_diag","phenoage_res_sd",cov_ukb_multivariable)],
                           cut=seq(45,80,5), episode= "agegroup", id="eid")
dat_timesplit[1:40,]
dat_timesplit$agegroup <- as.factor( dat_timesplit$agegroup )
str( dat_timesplit )
fit_cox_tv_phenoage_anycancer <- coxph( Surv(tstart,tstop,any_cancer_diag) ~ phenoage_res_sd:agegroup +
                                          birth_year_cat + sex+assessment_centre + ethnicity_4group + bmi_cat + smoking + alcohol_3group + PA_cat + education + di_5group, data=dat_timesplit )
summary( fit_cox_tv_phenoage_anycancer )

plot_cox_tv_phenoage_anycancer <- ggplot(data.frame(
  hr=summary(fit_cox_tv_phenoage_anycancer)$coef[paste0("phenoage_res_sd:agegroup", 1:9), 2],
  lb=summary(fit_cox_tv_phenoage_anycancer)$conf.int[paste0("phenoage_res_sd:agegroup", 1:9), 3],
  ub=summary(fit_cox_tv_phenoage_anycancer)$conf.int[paste0("phenoage_res_sd:agegroup", 1:9), 4],
  startage=seq(40,80,5),
  endage=seq(45,85,5))) + 
  geom_hline(aes(yintercept=1, linetype="Reference line (HR=1)", color = "Reference line (HR=1)")) +
  geom_hline(aes(yintercept=nonPH[nonPH$BA=="phenoage_res_sd" & nonPH$outcome=="Any cancer","HR"], linetype="Time-fixed HR", color = "Time-fixed HR"), size=.6) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, 
                ymin = nonPH[nonPH$BA=="phenoage_res_sd" & nonPH$outcome=="Any cancer","lb"], 
                ymax = nonPH[nonPH$BA=="phenoage_res_sd" & nonPH$outcome=="Any cancer","ub"]), 
            fill = "#D55E00", alpha = .01) +
  geom_rect(aes(xmin = startage, xmax = endage, 
                ymin = lb, ymax = ub),
            fill = "#0072B2", alpha = .15) + 
  geom_step(aes(x=startage, y=hr), size = .6, color="#0072B2") + 
  geom_segment(aes(x = tail(startage,1), y = tail(hr,1), xend = tail(endage,1), yend = tail(hr,1), color="Time-varying HR", linetype="Time-varying HR"), size = .6) +
  scale_y_continuous(breaks = seq(.6,1.6,.2), trans=scales::log_trans()) +
  coord_cartesian(ylim=c(.6,1.6)) + 
  labs(title="Any cancer", x="Age (years)", y="HR per 1 SD increase") + 
  scale_color_manual(name = "PhenoAge residual", values = c("Time-varying HR"="#0072B2", "Time-fixed HR"="#D55E00", "Reference line (HR=1)"="grey60"), guide=guide_legend(reverse = T)) + 
  scale_linetype_manual(name = "PhenoAge residual", values = c("Time-varying HR"=1, "Time-fixed HR"=6, "Reference line (HR=1)"=1), guide=guide_legend(reverse = T)) +
  theme_bw() +
  theme(plot.title = element_text(face="bold", size = 11),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.position = c(.75,.18),
        legend.title = element_text(face="bold", size = 8),
        legend.text = element_text(size = 7),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.key.height = unit(4, "mm"),
        panel.grid.major=element_line(colour="grey96"),
        panel.grid.minor=element_line(colour="grey96"))

plot_cox_tv_phenoage_anycancer


####### Time-varying HR for PhenoAge and breast cancer #######

dat_timesplit <- survSplit(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ . , 
                           ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women",c("interview_date","birth_date","censor_date_bc","bc_diag","phenoage_res_sd",cov_ukb_bc)],
                           cut=seq(45,80,5), episode= "agegroup", id="eid")
dat_timesplit$agegroup <- as.factor( dat_timesplit$agegroup )
fit_cox_tv_phenoage_bc <- coxph( Surv(tstart,tstop,bc_diag) ~ phenoage_res_sd:agegroup +
                                          birth_year_cat + sex+assessment_centre + ethnicity_4group + bmi_cat + smoking + alcohol_3group + PA_cat + education + di_5group + family_history_breast_cancer + bc_screening + menopause + hormone_therapy + oral_contraceptive + parity, 
                                 data=dat_timesplit )
summary( fit_cox_tv_phenoage_bc )

plot_cox_tv_phenoage_bc <- ggplot(data.frame(
  hr=summary(fit_cox_tv_phenoage_bc)$coef[paste0("phenoage_res_sd:agegroup", 1:9), 2],
  lb=summary(fit_cox_tv_phenoage_bc)$conf.int[paste0("phenoage_res_sd:agegroup", 1:9), 3],
  ub=summary(fit_cox_tv_phenoage_bc)$conf.int[paste0("phenoage_res_sd:agegroup", 1:9), 4],
  startage=seq(40,80,5),
  endage=seq(45,85,5))) + 
  geom_hline(aes(yintercept=1, linetype="Reference line (HR=1)", color = "Reference line (HR=1)")) +
  geom_hline(aes(yintercept=nonPH[nonPH$BA=="phenoage_res_sd" & nonPH$outcome=="Breast cancer in women","HR"], linetype="Time-fixed HR", color = "Time-fixed HR"), size=.6) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, 
                ymin = nonPH[nonPH$BA=="phenoage_res_sd" & nonPH$outcome=="Breast cancer in women","lb"], 
                ymax = nonPH[nonPH$BA=="phenoage_res_sd" & nonPH$outcome=="Breast cancer in women","ub"]), 
            fill = "#D55E00", alpha = .01) +
  geom_rect(aes(xmin = startage, xmax = endage, 
                ymin = lb, ymax = ub),
            fill = "#0072B2", alpha = .15) + 
  geom_step(aes(x=startage, y=hr), size = .6, color="#0072B2") + 
  geom_segment(aes(x = tail(startage,1), y = tail(hr,1), xend = tail(endage,1), yend = tail(hr,1), color="Time-varying HR", linetype="Time-varying HR"), size = .6) +
  scale_y_continuous(breaks = seq(.6,1.6,.2), trans=scales::log_trans()) +
  coord_cartesian(ylim=c(.6,1.6)) + 
  labs(title="Breast cancer in women", x="Age (years)", y="HR per 1 SD increase") + 
  scale_color_manual(name = "PhenoAge residual", values = c("Time-varying HR"="#0072B2", "Time-fixed HR"="#D55E00", "Reference line (HR=1)"="grey60"), guide=guide_legend(reverse = T)) + 
  scale_linetype_manual(name = "PhenoAge residual", values = c("Time-varying HR"=1, "Time-fixed HR"=6, "Reference line (HR=1)"=1), guide=guide_legend(reverse = T)) +
  theme_bw() +
  theme(plot.title = element_text(face="bold", size = 11),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.position = c(.75,.18),
        legend.title = element_text(face="bold", size = 8),
        legend.text = element_text(size = 7),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.key.height = unit(4, "mm"),
        panel.grid.major=element_line(colour="grey96"),
        panel.grid.minor=element_line(colour="grey96"))

plot_cox_tv_phenoage_bc


####### Time-varying HR for KDM and colorectal cancer #######

dat_timesplit <- survSplit(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ . , 
                           ukb_bioage_cancer[,c("interview_date","birth_date","censor_date_cc","cc_diag","kdm_res_sd",cov_ukb_cc)],
                           cut=seq(45,80,5), episode= "agegroup", id="eid")
dat_timesplit$agegroup <- as.factor( dat_timesplit$agegroup )
fit_cox_tv_kdm_cc <- coxph( Surv(tstart,tstop,cc_diag) ~ kdm_res_sd:agegroup +
                                          birth_year_cat + sex+assessment_centre + ethnicity_4group + bmi_cat + smoking + alcohol_3group + PA_cat + education + di_5group + family_history_bowel_cancer + cc_screening + fresh_veg_cat + red_meat_cat + processed_meat_cat, data=dat_timesplit )
summary( fit_cox_tv_kdm_cc )

plot_cox_tv_kdm_cc <- ggplot(data.frame(
  hr=summary(fit_cox_tv_kdm_cc)$coef[paste0("kdm_res_sd:agegroup", 1:9), 2],
  lb=summary(fit_cox_tv_kdm_cc)$conf.int[paste0("kdm_res_sd:agegroup", 1:9), 3],
  ub=summary(fit_cox_tv_kdm_cc)$conf.int[paste0("kdm_res_sd:agegroup", 1:9), 4],
  startage=seq(40,80,5),
  endage=seq(45,85,5))) + 
  geom_hline(aes(yintercept=1, linetype="Reference line (HR=1)", color = "Reference line (HR=1)")) +
  geom_hline(aes(yintercept=nonPH[nonPH$BA=="kdm_res_sd" & nonPH$outcome=="Colorectal cancer","HR"], linetype="Time-fixed HR", color = "Time-fixed HR"), size=.6) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, 
                ymin = nonPH[nonPH$BA=="kdm_res_sd" & nonPH$outcome=="Colorectal cancer","lb"], 
                ymax = nonPH[nonPH$BA=="kdm_res_sd" & nonPH$outcome=="Colorectal cancer","ub"]), 
            fill = "#D55E00", alpha = .01) +
  geom_rect(aes(xmin = startage, xmax = endage, 
                ymin = lb, ymax = ub),
            fill = "#0072B2", alpha = .15) + 
  geom_step(aes(x=startage, y=hr), size = .6, color="#0072B2") + 
  geom_segment(aes(x = tail(startage,1), y = tail(hr,1), xend = tail(endage,1), yend = tail(hr,1), color="Time-varying HR", linetype="Time-varying HR"), size = .6) +
  scale_y_continuous(breaks = seq(.6,1.6,.2), trans=scales::log_trans()) +
  coord_cartesian(ylim=c(.6,1.6)) + 
  labs(title="Colorectal cancer", x="Age (years)", y="HR per 1 SD increase") + 
  scale_color_manual(name = "KDM residual", values = c("Time-varying HR"="#0072B2", "Time-fixed HR"="#D55E00", "Reference line (HR=1)"="grey60"), guide=guide_legend(reverse = T)) + 
  scale_linetype_manual(name = "KDM residual", values = c("Time-varying HR"=1, "Time-fixed HR"=6, "Reference line (HR=1)"=1), guide=guide_legend(reverse = T)) +
  theme_bw() +
  theme(plot.title = element_text(face="bold", size = 11),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.position = c(.75,.18),
        legend.title = element_text(face="bold", size = 8),
        legend.text = element_text(size = 7),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.key.height = unit(4, "mm"),
        panel.grid.major=element_line(colour="grey96"),
        panel.grid.minor=element_line(colour="grey96"))

plot_cox_tv_kdm_cc


####### Time-varying HR for PhenoAge and colorectal cancer #######

dat_timesplit <- survSplit(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ . , 
                           ukb_bioage_cancer[,c("interview_date","birth_date","censor_date_cc","cc_diag","phenoage_res_sd",cov_ukb_cc)],
                           cut=seq(45,80,5), episode= "agegroup", id="eid")
dat_timesplit$agegroup <- as.factor( dat_timesplit$agegroup )
fit_cox_tv_phenoage_cc <- coxph( Surv(tstart,tstop,cc_diag) ~ phenoage_res_sd:agegroup +
                                   birth_year_cat + sex+assessment_centre + ethnicity_4group + bmi_cat + smoking + alcohol_3group + PA_cat + education + di_5group + family_history_bowel_cancer + cc_screening + fresh_veg_cat + red_meat_cat + processed_meat_cat, data=dat_timesplit )
summary( fit_cox_tv_phenoage_cc )

plot_cox_tv_phenoage_cc <- ggplot(data.frame(
  hr=summary(fit_cox_tv_phenoage_cc)$coef[paste0("phenoage_res_sd:agegroup", 1:9), 2],
  lb=summary(fit_cox_tv_phenoage_cc)$conf.int[paste0("phenoage_res_sd:agegroup", 1:9), 3],
  ub=summary(fit_cox_tv_phenoage_cc)$conf.int[paste0("phenoage_res_sd:agegroup", 1:9), 4],
  startage=seq(40,80,5),
  endage=seq(45,85,5))) + 
  geom_hline(aes(yintercept=1, linetype="Reference line (HR=1)", color = "Reference line (HR=1)")) +
  geom_hline(aes(yintercept=nonPH[nonPH$BA=="phenoage_res_sd" & nonPH$outcome=="Colorectal cancer","HR"], linetype="Time-fixed HR", color = "Time-fixed HR"), size=.6) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, 
                ymin = nonPH[nonPH$BA=="phenoage_res_sd" & nonPH$outcome=="Colorectal cancer","lb"], 
                ymax = nonPH[nonPH$BA=="phenoage_res_sd" & nonPH$outcome=="Colorectal cancer","ub"]), 
            fill = "#D55E00", alpha = .01) +
  geom_rect(aes(xmin = startage, xmax = endage, 
                ymin = lb, ymax = ub),
            fill = "#0072B2", alpha = .15) + 
  geom_step(aes(x=startage, y=hr), size = .6, color="#0072B2") + 
  geom_segment(aes(x = tail(startage,1), y = tail(hr,1), xend = tail(endage,1), yend = tail(hr,1), color="Time-varying HR", linetype="Time-varying HR"), size = .6) +
  scale_y_continuous(breaks = seq(.6,1.6,.2), trans=scales::log_trans()) +
  coord_cartesian(ylim=c(.6,1.6)) + 
  labs(title="Colorectal cancer", x="Age (years)", y="HR per 1 SD increase") + 
  scale_color_manual(name = "PhenoAge residual", values = c("Time-varying HR"="#0072B2", "Time-fixed HR"="#D55E00", "Reference line (HR=1)"="grey60"), guide=guide_legend(reverse = T)) + 
  scale_linetype_manual(name = "PhenoAge residual", values = c("Time-varying HR"=1, "Time-fixed HR"=6, "Reference line (HR=1)"=1), guide=guide_legend(reverse = T)) +
  theme_bw() +
  theme(plot.title = element_text(face="bold", size = 11),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.position = c(.75,.18),
        legend.title = element_text(face="bold", size = 8),
        legend.text = element_text(size = 7),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.key.height = unit(4, "mm"),
        panel.grid.major=element_line(colour="grey96"),
        panel.grid.minor=element_line(colour="grey96"))

plot_cox_tv_phenoage_cc


####### Combine time-varying HR plots #######

jpeg("Output/Timevarying_HR_BioAge_cancer.jpg", width = 16, height = 13, units = "cm", res = 1000, quality = 100)
plot_grid(plot_cox_tv_phenoage_anycancer, plot_cox_tv_phenoage_bc, plot_cox_tv_kdm_cc, plot_cox_tv_phenoage_cc, labels = "AUTO",  label_size = 12)
dev.off()



#============= COX MODELS FOR INDIVIDUAL BIOMARKERS AND CANCERS ===============

# 18 biomarkers to be tested
biomarkers <- c("fev_1000","sbp","bun","hba1c","totchol","creat_umol","glucose_mmol",
                "waist","rdw","albumin_gl","alp","trig","mcv","uap","lymph","rbc","crp","dbp")
names(biomarkers) <-c("Forced expiratory volume",
                      "Systolic blood pressure",
                      "Blood urea nitrogen",
                      "Glycated hemoglobin",
                      "Total cholesterol",
                      "Creatinine",
                      "Serum glucose",
                      "Waist circumference",
                      "Red cell distribution width",
                      "Albumin",
                      "Alkaline phosphatase",
                      "Triglyceride",
                      "Mean cell volume",
                      "Uric acid",
                      "Lymphocyte",
                      "Red blood cell count",
                      "C-reactive protein",
                      "Diastolic blood pressure")

ukb_biomarker_cancer <- cbind(ukb_bioage_cancer[,!colnames(ukb_bioage_cancer) %in% biomarkers],
                              scale(ukb_bioage_cancer[,biomarkers], center=T, scale=T)) # Standardization of biomarkers

# Cox models of biomarkers on cancers, fully adjusted models
HR_biomarker_anycancer <- data.frame(biomarker=biomarkers, outcome="Any cancer", beta=NA, se=NA, HR=NA, lb=NA, ub=NA, p=NA)
HR_biomarker_bc <- data.frame(biomarker=biomarkers, outcome="Breast cancer in women", beta=NA, se=NA, HR=NA, lb=NA, ub=NA, p=NA)
HR_biomarker_pc <- data.frame(biomarker=biomarkers, outcome="Prostate cancer in men", beta=NA, se=NA, HR=NA, lb=NA, ub=NA, p=NA)
HR_biomarker_lc <- data.frame(biomarker=biomarkers, outcome="Lung cancer", beta=NA, se=NA, HR=NA, lb=NA, ub=NA, p=NA)
HR_biomarker_cc <- data.frame(biomarker=biomarkers, outcome="Colorectal cancer", beta=NA, se=NA, HR=NA, lb=NA, ub=NA, p=NA)
HR_biomarker_mel <- data.frame(biomarker=biomarkers, outcome="Melanoma", beta=NA, se=NA, HR=NA, lb=NA, ub=NA, p=NA)

for (j in 1:length(biomarkers)) {
  
  ###### Any cancer ######
  coef <- summary(coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ . , 
                        ukb_biomarker_cancer[,c("interview_date","birth_date","censor_date_any_cancer","any_cancer_diag",biomarkers[j],cov_ukb_multivariable)]))
  HR_biomarker_anycancer[j,c(3,4,8)] <- coef$coefficients[1,c(1,3,5)]
  HR_biomarker_anycancer[j,5:7] <- coef$conf.int[1,c(1,3,4)]
  
  ###### Breast cancer in women ######
  coef <- summary(coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ . , 
                        ukb_biomarker_cancer[ukb_biomarker_cancer$sex=="Women",c("interview_date","birth_date","censor_date_bc","bc_diag",biomarkers[j],cov_ukb_bc)]))
  HR_biomarker_bc[j,c(3,4,8)] <- coef$coefficients[1,c(1,3,5)]
  HR_biomarker_bc[j,5:7] <- coef$conf.int[1,c(1,3,4)]
  
  ###### Prostate cancer in men ###### 
  coef <- summary(coxph(Surv((interview_date-birth_date)/365.25, (censor_date_pc-birth_date)/365.25, pc_diag) ~ . , 
                        ukb_biomarker_cancer[ukb_biomarker_cancer$sex=="Men",c("interview_date","birth_date","censor_date_pc","pc_diag",biomarkers[j],cov_ukb_pc)]))
  HR_biomarker_pc[j,c(3,4,8)] <- coef$coefficients[1,c(1,3,5)]
  HR_biomarker_pc[j,5:7] <- coef$conf.int[1,c(1,3,4)]
  
  ###### Lung cancer ###### 
  coef <- summary(coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
                        ukb_biomarker_cancer[,c("interview_date","birth_date","censor_date_lc","lc_diag",biomarkers[j],cov_ukb_lc)]))
  HR_biomarker_lc[j,c(3,4,8)] <- coef$coefficients[1,c(1,3,5)]
  HR_biomarker_lc[j,5:7] <- coef$conf.int[1,c(1,3,4)]
  
  ###### Colorectal cancer ###### 
  coef <- summary(coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ . , 
                        ukb_biomarker_cancer[,c("interview_date","birth_date","censor_date_cc","cc_diag",biomarkers[j],cov_ukb_cc)]))
  HR_biomarker_cc[j,c(3,4,8)] <- coef$coefficients[1,c(1,3,5)]
  HR_biomarker_cc[j,5:7] <- coef$conf.int[1,c(1,3,4)]
  
  ###### Melanoma ###### 
  coef <- summary(coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ . , 
                        ukb_biomarker_cancer[ukb_biomarker_cancer$sun_protection!="missing"&ukb_biomarker_cancer$ethnicity_4group!="missing",c("interview_date","birth_date","censor_date_mel","mel_diag",biomarkers[j],cov_ukb_mel)]))
  HR_biomarker_mel[j,c(3,4,8)] <- coef$coefficients[1,c(1,3,5)]
  HR_biomarker_mel[j,5:7] <- coef$conf.int[1,c(1,3,4)]
}				

### Output results
HR_biomarker_cancer <- rbind(HR_biomarker_anycancer,HR_biomarker_bc,HR_biomarker_pc,HR_biomarker_lc,HR_biomarker_cc,HR_biomarker_mel)
write.table(HR_biomarker_cancer, 
            file = "Output/Cox_biomarkers_cancers.txt", sep = "\t", row.names = F, na="")

### Forest plot
#HR_biomarker_cancer <- read.delim("Output/Cox_biomarkers_cancers.txt")
#HR_BA_cancer <- read.delim("Output/Cox_BA_cancers.txt")
jpeg(file="Output/Cox_biomarkers_cancers.jpg", width = 18, height = 16, units = "cm", res = 600, quality=100)
ggforestplot::forestplot(
  df = rbind(HR_biomarker_cancer %>% 
               rename(name=biomarker) %>%
               mutate(group="Clinical biomarkers"), 
             HR_BA_cancer %>% 
               filter(BA %in% c("kdm_res_sd","phenoage_res_sd","hd_log_sd") & 
                        (outcome=="Any cancer" & model=="Multivariable  model" |
                           outcome=="Breast cancer in women" & model=="Breast cancer-specific model" |
                           outcome=="Prostate cancer in men" & model=="Prostate cancer-specific model" |
                           outcome=="Lung cancer" & model=="Lung cancer-specific model" |
                           outcome=="Colorectal cancer" & model=="Colorectal cancer-specific model" |
                           outcome=="Melanoma" & model=="Melanoma-specific model" )) %>% 
               select(BA,outcome,beta,se,HR,lb,ub,p) %>%
               rename(name=BA) %>%
               mutate(group="Biological age measures")) %>%
    mutate(name=c(rep(names(biomarkers),6),
                  rep(c("KDM residual","PhenoAge residual","HD (log)"), 6)),
           outcome=factor(outcome, levels = c("Melanoma", "Colorectal cancer", "Lung cancer", "Prostate cancer in men", "Breast cancer in women", "Any cancer")),
           group=factor(group, levels = c("Clinical biomarkers", "Biological age measures"))),
  name = name,
  estimate = beta,
  se = se,
  pvalue = p,
  logodds = T,
  psignif = .05/15,
  xlab = "Hazard ratio per 1-SD increase",
  xtickbreaks = c(.6,.7,.8,.9,1,1.1,1.2,1.3,1.4),
  colour = outcome,
  shape = outcome) +
  labs(colour = "Outcome", shape = "Outcome") +
  ggsci::scale_color_npg(palette = c("nrc"), alpha = 1) +
  scale_shape_manual(values = c(23L, 22L, 24L, 23L, 22L, 21L)) +
  guides(colour = guide_legend(override.aes = list(size=1.5), reverse = TRUE)) + 
  theme(axis.title = element_text(size = 11),
        axis.text = element_text(size = 8),
        legend.position="right",
        legend.title = element_text(face = "bold", size=10),
        legend.background = element_rect(colour = "black", fill = "white", linetype="solid", size=.1),
        legend.text = element_text(size = 8))  +
  ggforce::facet_col(
    facets = ~group,
    scales = "free_y",
    space = "free"
  )
dev.off()






#================== SPLINE MODELS FOR BA MEASURES AND CANCERS =================

#HR_BA_cancer <- read.delim("Output/Cox_BA_cancers.txt") # Results from linear models

###### Any cancer ######

### KDM
dd <- datadist(ukb_bioage_cancer); dd$limits$kdm_res_sd[2] <- 0; options(datadist="dd")
fit_cox_anycancer_kdm_spline <- cph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ 
                                      rcs(kdm_res_sd,4) +
                                      birth_year_cat+sex+assessment_centre+ethnicity_4group+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group,
                                    data=ukb_bioage_cancer)
anova(fit_cox_anycancer_kdm_spline) # p<.0001
plot_cox_anycancer_kdm_spline <- 
  ggplot(data = Predict(fit_cox_anycancer_kdm_spline, kdm_res_sd, ref.zero = T, fun=exp)) +
  geom_hline(aes(yintercept = 1), color="grey") +
  geom_abline(intercept = 0, 
              slope = HR_BA_cancer %>% filter(BA=="kdm_res_sd", outcome=="Any cancer", model=="Multivariable  model") %>% select(beta) %>% as.numeric(),
              linetype="dashed") +
  geom_line(aes(x=kdm_res_sd, y=yhat)) +
  annotate("text", x=1, y=.55, size=3.5, label="P for non-linearity<.001") + 
  scale_x_continuous(breaks=seq(-3,3,1), limits=c(-3,3), labels=c("-3SD","-2SD","-1SD","Mean","+1SD","+2SD","+3SD")) +
  scale_y_continuous(breaks=c(.5,1,1.5,2,2.5), trans=scales::log_trans()) +
  coord_cartesian(ylim=c(.5,2.5)) + 
  labs(title="Any cancer", x="KDM residual", y="Hazard Ratio") + 
  theme_bw() + 
  theme(plot.title = element_text(size=11, face="bold"), 
        plot.caption = element_blank(), 
        axis.title = element_text(face="bold", size=9),
        axis.text = element_text(size=7))

### PhenoAge residual
dd <- datadist(ukb_bioage_cancer); dd$limits$phenoage_res_sd[2] <- 0; options(datadist="dd")
fit_cox_anycancer_phenoage_spline <- cph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ 
                                           rcs(phenoage_res_sd,4) +
                                           birth_year_cat+sex+assessment_centre+ethnicity_4group+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group,
                                         data=ukb_bioage_cancer)
anova(fit_cox_anycancer_phenoage_spline) # p<.0001
plot_cox_anycancer_phenoage_spline <- 
  ggplot(data = Predict(fit_cox_anycancer_phenoage_spline, phenoage_res_sd, ref.zero = T, fun=exp)) +
  geom_hline(aes(yintercept = 1), color="grey") +
  geom_abline(intercept = 0, 
              slope = HR_BA_cancer %>% filter(BA=="phenoage_res_sd", outcome=="Any cancer", model=="Multivariable  model") %>% select(beta) %>% as.numeric(),
              linetype="dashed") +
  geom_line(aes(x=phenoage_res_sd, y=yhat)) +
  annotate("text", x=1, y=.55, size=3.5, label="P for non-linearity<.001") + 
  scale_x_continuous(breaks=seq(-3,3,1), limits=c(-3,3), labels=c("-3SD","-2SD","-1SD","Mean","+1SD","+2SD","+3SD")) +
  scale_y_continuous(breaks=c(.5,1,1.5,2,2.5), trans=scales::log_trans()) +
  coord_cartesian(ylim=c(.5,2.5)) + 
  labs(title="",x="PhenoAge residual", y="Hazard Ratio") + 
  theme_bw() + 
  theme(plot.title = element_text(size=11, face="bold"), 
        plot.caption = element_blank(), 
        axis.title = element_text(face="bold", size=9),
        axis.text = element_text(size=7))

### Homeostatic dysregulation
dd <- datadist(ukb_bioage_cancer); dd$limits$hd_log_sd[2] <- 0; options(datadist="dd")
fit_cox_anycancer_hd_spline <- cph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ 
                                     rcs(hd_log_sd,4) +
                                     birth_year_cat+sex+assessment_centre+ethnicity_4group+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group,
                                   data=ukb_bioage_cancer)
anova(fit_cox_anycancer_hd_spline) # p=.95
plot_cox_anycancer_hd_spline <- 
  ggplot(data = Predict(fit_cox_anycancer_hd_spline, hd_log_sd, ref.zero = T, fun=exp)) +
  geom_hline(aes(yintercept = 1), color="grey") +
  geom_abline(intercept = 0, 
              slope = HR_BA_cancer %>% filter(BA=="hd_log_sd", outcome=="Any cancer", model=="Multivariable  model") %>% select(beta) %>% as.numeric(),
              linetype="dashed") +
  geom_line(aes(x=hd_log_sd, y=yhat)) +
  annotate("text", x=1, y=.55, size=3.5, label="P for non-linearity=.95") + 
  scale_x_continuous(breaks=seq(-3,3,1), limits=c(-3,3), labels=c("-3SD","-2SD","-1SD","Mean","+1SD","+2SD","+3SD")) +
  scale_y_continuous(breaks=c(.5,1,1.5,2,2.5), trans=scales::log_trans()) +
  coord_cartesian(ylim=c(.5,2.5)) + 
  labs(title="",x="HD (log)", y="Hazard Ratio") + 
  theme_bw() + 
  theme(plot.title = element_text(size=11, face="bold"), 
        plot.caption = element_blank(), 
        axis.title = element_text(face="bold", size=9),
        axis.text = element_text(size=7))


###### Breast cancer in women ######

### KDM
dd <- datadist(ukb_bioage_cancer); dd$limits$kdm_res_sd[2] <- 0; options(datadist="dd")
fit_cox_bc_kdm_spline <- cph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ 
                               rcs(kdm_res_sd,4) +
                               birth_year_cat+assessment_centre+ethnicity_4group+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group+family_history_breast_cancer+bc_screening+menopause+hormone_therapy+oral_contraceptive+parity,
                             data=ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women",])
anova(fit_cox_bc_kdm_spline) # p=.88
plot_cox_bc_kdm_spline <- 
  ggplot(data = Predict(fit_cox_bc_kdm_spline, kdm_res_sd, ref.zero = T, fun=exp)) +
  geom_hline(aes(yintercept = 1), color="grey") +
  geom_abline(intercept = 0, 
              slope = HR_BA_cancer %>% filter(BA=="kdm_res_sd", outcome=="Breast cancer in women", model=="Breast cancer-specific model") %>% select(beta) %>% as.numeric(),
              linetype="dashed") +
  geom_line(aes(x=kdm_res_sd, y=yhat)) +
  annotate("text", x=1, y=.55, size=3.5, label="P for non-linearity=.88") + 
  scale_x_continuous(breaks=seq(-3,3,1), limits=c(-3,3), labels=c("-3SD","-2SD","-1SD","Mean","+1SD","+2SD","+3SD")) +
  scale_y_continuous(breaks=c(.5,1,1.5,2,2.5), trans=scales::log_trans()) +
  coord_cartesian(ylim=c(.5,2.5)) + 
  labs(title="Breast cancer in women",x="KDM residual", y="Hazard Ratio") + 
  theme_bw() + 
  theme(plot.title = element_text(size=11, face="bold"), 
        plot.caption = element_blank(), 
        axis.title = element_text(face="bold", size=9),
        axis.text = element_text(size=7))

### PhenoAge
dd <- datadist(ukb_bioage_cancer); dd$limits$phenoage_res_sd[2] <- 0; options(datadist="dd")
fit_cox_bc_phenoage_spline <- cph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ 
                                    rcs(phenoage_res_sd,4) +
                                    birth_year_cat+assessment_centre+ethnicity_4group+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group+family_history_breast_cancer+bc_screening+menopause+hormone_therapy+oral_contraceptive+parity,
                                  data=ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women",])
anova(fit_cox_bc_phenoage_spline) # p=.78
plot_cox_bc_phenoage_spline <- 
  ggplot(data = Predict(fit_cox_bc_phenoage_spline, phenoage_res_sd, ref.zero = T, fun=exp)) +
  geom_hline(aes(yintercept = 1), color="grey") +
  geom_abline(intercept = 0, 
              slope = HR_BA_cancer %>% filter(BA=="phenoage_res_sd", outcome=="Breast cancer in women", model=="Breast cancer-specific model") %>% select(beta) %>% as.numeric(),
              linetype="dashed") +
  geom_line(aes(x=phenoage_res_sd, y=yhat)) +
  annotate("text", x=1, y=.55, size=3.5, label="P for non-linearity=.78") + 
  scale_x_continuous(breaks=seq(-3,3,1), limits=c(-3,3), labels=c("-3SD","-2SD","-1SD","Mean","+1SD","+2SD","+3SD")) +
  scale_y_continuous(breaks=c(.5,1,1.5,2,2.5), trans=scales::log_trans()) +
  coord_cartesian(ylim=c(.5,2.5)) + 
  labs(title="",x="PhenoAge residual", y="Hazard Ratio") + 
  theme_bw() + 
  theme(plot.title = element_text(size=11, face="bold"), 
        plot.caption = element_blank(), 
        axis.title = element_text(face="bold", size=9),
        axis.text = element_text(size=7))

### HD
dd <- datadist(ukb_bioage_cancer); dd$limits$hd_log_sd[2] <- 0; options(datadist="dd")
fit_cox_bc_hd_spline <- cph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ 
                              rcs(hd_log_sd,4) +
                              birth_year_cat+assessment_centre+ethnicity_4group+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group+family_history_breast_cancer+bc_screening+menopause+hormone_therapy+oral_contraceptive+parity,
                            data=ukb_bioage_cancer[ukb_bioage_cancer$sex=="Women",])
anova(fit_cox_bc_hd_spline) # p=.055
plot_cox_bc_hd_spline <- 
  ggplot(data = Predict(fit_cox_bc_hd_spline, hd_log_sd, ref.zero = T, fun=exp)) +
  geom_hline(aes(yintercept = 1), color="grey") +
  geom_abline(intercept = 0, 
              slope = HR_BA_cancer %>% filter(BA=="hd_log_sd", outcome=="Breast cancer in women", model=="Breast cancer-specific model") %>% select(beta) %>% as.numeric(),
              linetype="dashed") +
  geom_line(aes(x=hd_log_sd, y=yhat)) +
  annotate("text", x=1, y=.55, size=3.5, label="P for non-linearity=.055") + 
  scale_x_continuous(breaks=seq(-3,3,1), limits=c(-3,3), labels=c("-3SD","-2SD","-1SD","Mean","+1SD","+2SD","+3SD")) +
  scale_y_continuous(breaks=c(.5,1,1.5,2,2.5), trans=scales::log_trans()) +
  coord_cartesian(ylim=c(.5,2.5)) + 
  labs(title="",x="HD (log)", y="Hazard Ratio") + 
  theme_bw() + 
  theme(plot.title = element_text(size=11, face="bold"), 
        plot.caption = element_blank(), 
        axis.title = element_text(face="bold", size=9),
        axis.text = element_text(size=7))


###### Prostate cancer in men ######

### KDM
dd <- datadist(ukb_bioage_cancer); dd$limits$kdm_res_sd[2] <- 0; options(datadist="dd")
fit_cox_pc_kdm_spline <- cph(Surv((interview_date-birth_date)/365.25, (censor_date_pc-birth_date)/365.25, pc_diag) ~ 
                               rcs(kdm_res_sd,4) +
                               birth_year_cat+assessment_centre+ethnicity_4group+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group+family_history_prostate_cancer+psa_test+diabetes,
                             data=ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men",])
anova(fit_cox_pc_kdm_spline) # p=.037
plot_cox_pc_kdm_spline <- 
  ggplot(data = Predict(fit_cox_pc_kdm_spline, kdm_res_sd, ref.zero = T, fun=exp)) +
  geom_hline(aes(yintercept = 1), color="grey") +
  geom_abline(intercept = 0, 
              slope = HR_BA_cancer %>% filter(BA=="kdm_res_sd", outcome=="Prostate cancer in men", model=="Prostate cancer-specific model") %>% select(beta) %>% as.numeric(),
              linetype="dashed") +
  geom_line(aes(x=kdm_res_sd, y=yhat)) +
  annotate("text", x=1, y=.55, size=3.5, label="P for non-linearity=.037") + 
  scale_x_continuous(breaks=seq(-3,3,1), limits=c(-3,3), labels=c("-3SD","-2SD","-1SD","Mean","+1SD","+2SD","+3SD")) +
  scale_y_continuous(breaks=c(.5,1,1.5,2,2.5), trans=scales::log_trans()) +
  coord_cartesian(ylim=c(.5,2.5)) + 
  labs(title="Prostate cancer in men",x="KDM residual", y="Hazard Ratio") + 
  theme_bw() + 
  theme(plot.title = element_text(size=11, face="bold"), 
        plot.caption = element_blank(), 
        axis.title = element_text(face="bold", size=9),
        axis.text = element_text(size=7))

### PhenoAge
dd <- datadist(ukb_bioage_cancer); dd$limits$phenoage_res_sd[2] <- 0; options(datadist="dd")
fit_cox_pc_phenoage_spline <- cph(Surv((interview_date-birth_date)/365.25, (censor_date_pc-birth_date)/365.25, pc_diag) ~ 
                                    rcs(phenoage_res_sd,4) +
                                    birth_year_cat+assessment_centre+ethnicity_4group+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group+family_history_prostate_cancer+psa_test+diabetes,
                                  data=ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men",])
anova(fit_cox_pc_phenoage_spline) # p=.17
plot_cox_pc_phenoage_spline <- 
  ggplot(data = Predict(fit_cox_pc_phenoage_spline, phenoage_res_sd, ref.zero = T, fun=exp)) +
  geom_hline(aes(yintercept = 1), color="grey") +
  geom_abline(intercept = 0, 
              slope = HR_BA_cancer %>% filter(BA=="phenoage_res_sd", outcome=="Prostate cancer in men", model=="Prostate cancer-specific model") %>% select(beta) %>% as.numeric(),
              linetype="dashed") +
  geom_line(aes(x=phenoage_res_sd, y=yhat)) +
  annotate("text", x=1, y=.55, size=3.5, label="P for non-linearity=.17") + 
  scale_x_continuous(breaks=seq(-3,3,1), limits=c(-3,3), labels=c("-3SD","-2SD","-1SD","Mean","+1SD","+2SD","+3SD")) +
  scale_y_continuous(breaks=c(.5,1,1.5,2,2.5), trans=scales::log_trans()) +
  coord_cartesian(ylim=c(.5,2.5)) + 
  labs(title="",x="PhenoAge residual", y="Hazard Ratio") + 
  theme_bw() + 
  theme(plot.title = element_text(size=11, face="bold"), 
        plot.caption = element_blank(), 
        axis.title = element_text(face="bold", size=9),
        axis.text = element_text(size=7))

### HD
dd <- datadist(ukb_bioage_cancer); dd$limits$hd_log_sd[2] <- 0; options(datadist="dd")
fit_cox_pc_hd_spline <- cph(Surv((interview_date-birth_date)/365.25, (censor_date_pc-birth_date)/365.25, pc_diag) ~ 
                              rcs(hd_log_sd,4) +
                              birth_year_cat+assessment_centre+ethnicity_4group+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group+family_history_prostate_cancer+psa_test+diabetes,
                            data=ukb_bioage_cancer[ukb_bioage_cancer$sex=="Men",])
anova(fit_cox_pc_hd_spline) # p=.024
plot_cox_pc_hd_spline <- 
  ggplot(data = Predict(fit_cox_pc_hd_spline, hd_log_sd, ref.zero = T, fun=exp)) +
  geom_hline(aes(yintercept = 1), color="grey") +
  geom_abline(intercept = 0, 
              slope = HR_BA_cancer %>% filter(BA=="hd_log_sd", outcome=="Prostate cancer in men", model=="Prostate cancer-specific model") %>% select(beta) %>% as.numeric(),
              linetype="dashed") +
  geom_line(aes(x=hd_log_sd, y=yhat)) +
  annotate("text", x=1, y=.55, size=3.5, label="P for non-linearity=.024") + 
  scale_x_continuous(breaks=seq(-3,3,1), limits=c(-3,3), labels=c("-3SD","-2SD","-1SD","Mean","+1SD","+2SD","+3SD")) +
  scale_y_continuous(breaks=c(.5,1,1.5,2,2.5), trans=scales::log_trans()) +
  coord_cartesian(ylim=c(.5,2.5)) + 
  labs(title="",x="HD (log)", y="Hazard Ratio") + 
  theme_bw() + 
  theme(plot.title = element_text(size=11, face="bold"), 
        plot.caption = element_blank(), 
        axis.title = element_text(face="bold", size=9),
        axis.text = element_text(size=7))


###### Lung cancer ######

### KDM
dd <- datadist(ukb_bioage_cancer); dd$limits$kdm_res_sd[2] <- 0; options(datadist="dd")
fit_cox_lc_kdm_spline <- cph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ 
                               rcs(kdm_res_sd,4) +
                               birth_year_cat+sex+assessment_centre+ethnicity_4group+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group+family_history_lung_cancer,
                             data=ukb_bioage_cancer)
anova(fit_cox_lc_kdm_spline) # p=.011
plot_cox_lc_kdm_spline <- 
  ggplot(data = Predict(fit_cox_lc_kdm_spline, kdm_res_sd, ref.zero = T, fun=exp)) +
  geom_hline(aes(yintercept = 1), color="grey") +
  geom_abline(intercept = 0, 
              slope = HR_BA_cancer %>% filter(BA=="kdm_res_sd", outcome=="Lung cancer", model=="Lung cancer-specific model") %>% select(beta) %>% as.numeric(),
              linetype="dashed") +
  geom_line(aes(x=kdm_res_sd, y=yhat)) +
  annotate("text", x=1, y=.55, size=3.5, label="P for non-linearity=.011") + 
  scale_x_continuous(breaks=seq(-3,3,1), limits=c(-3,3), labels=c("-3SD","-2SD","-1SD","Mean","+1SD","+2SD","+3SD")) +
  scale_y_continuous(breaks=c(.5,1,1.5,2,2.5), trans=scales::log_trans()) +
  coord_cartesian(ylim=c(.5,2.5)) + 
  labs(title="Lung cancer",x="KDM residual", y="Hazard Ratio") + 
  theme_bw() + 
  theme(plot.title = element_text(size=11, face="bold"), 
        plot.caption = element_blank(), 
        axis.title = element_text(face="bold", size=9),
        axis.text = element_text(size=7))

### PhenoAge
dd <- datadist(ukb_bioage_cancer); dd$limits$phenoage_res_sd[2] <- 0; options(datadist="dd")
fit_cox_lc_phenoage_spline <- cph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ 
                                    rcs(phenoage_res_sd,4) +
                                    birth_year_cat+sex+assessment_centre+ethnicity_4group+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group+family_history_lung_cancer,
                                  data=ukb_bioage_cancer)
anova(fit_cox_lc_phenoage_spline) # p=.011
plot_cox_lc_phenoage_spline <-
  ggplot(data = Predict(fit_cox_lc_phenoage_spline, phenoage_res_sd, ref.zero = T, fun=exp)) +
  geom_hline(aes(yintercept = 1), color="grey") +
  geom_abline(intercept = 0, 
              slope = HR_BA_cancer %>% filter(BA=="phenoage_res_sd", outcome=="Lung cancer", model=="Lung cancer-specific model") %>% select(beta) %>% as.numeric(),
              linetype="dashed") +
  geom_line(aes(x=phenoage_res_sd, y=yhat)) +
  annotate("text", x=1, y=.55, size=3.5, label="P for non-linearity=.011") + 
  scale_x_continuous(breaks=seq(-3,3,1), limits=c(-3,3), labels=c("-3SD","-2SD","-1SD","Mean","+1SD","+2SD","+3SD")) +
  scale_y_continuous(breaks=c(.5,1,1.5,2,2.5), trans=scales::log_trans()) +
  coord_cartesian(ylim=c(.5,2.5)) + 
  labs(title="",x="PhenoAge residual", y="Hazard Ratio") + 
  theme_bw() + 
  theme(plot.title = element_text(size=11, face="bold"), 
        plot.caption = element_blank(), 
        axis.title = element_text(face="bold", size=9),
        axis.text = element_text(size=7))

### HD
dd <- datadist(ukb_bioage_cancer); dd$limits$hd_log_sd[2] <- 0; options(datadist="dd")
fit_cox_lc_hd_spline <- cph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ 
                              rcs(hd_log_sd,4) +
                              birth_year_cat+sex+assessment_centre+ethnicity_4group+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group+family_history_lung_cancer,
                            data=ukb_bioage_cancer)
anova(fit_cox_lc_hd_spline) # p=.30
plot_cox_lc_hd_spline <- 
  ggplot(data = Predict(fit_cox_lc_hd_spline, hd_log_sd, ref.zero = T, fun=exp)) +
  geom_hline(aes(yintercept = 1), color="grey") +
  geom_abline(intercept = 0, 
              slope = HR_BA_cancer %>% filter(BA=="hd_log_sd", outcome=="Lung cancer", model=="Lung cancer-specific model") %>% select(beta) %>% as.numeric(),
              linetype="dashed") +
  geom_line(aes(x=hd_log_sd, y=yhat)) +
  annotate("text", x=1, y=.55, size=3.5, label="P for non-linearity=.30") + 
  scale_x_continuous(breaks=seq(-3,3,1), limits=c(-3,3), labels=c("-3SD","-2SD","-1SD","Mean","+1SD","+2SD","+3SD")) +
  scale_y_continuous(breaks=c(.5,1,1.5,2,2.5), trans=scales::log_trans()) +
  coord_cartesian(ylim=c(.5,2.5)) + 
  labs(title="",x="HD (log)", y="Hazard Ratio") + 
  theme_bw() + 
  theme(plot.title = element_text(size=11, face="bold"), 
        plot.caption = element_blank(), 
        axis.title = element_text(face="bold", size=9),
        axis.text = element_text(size=7))


###### Colorectal cancer ######

### KDM
dd <- datadist(ukb_bioage_cancer); dd$limits$kdm_res_sd[2] <- 0; options(datadist="dd")
fit_cox_cc_kdm_spline <- cph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ 
                               rcs(kdm_res_sd,4) +
                               birth_year_cat+sex+assessment_centre+ethnicity_4group+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group+family_history_bowel_cancer+cc_screening+fresh_veg_cat+red_meat_cat+processed_meat_cat,
                             data=ukb_bioage_cancer)
anova(fit_cox_cc_kdm_spline) # p=.025
plot_cox_cc_kdm_spline <-
  ggplot(data = Predict(fit_cox_cc_kdm_spline, kdm_res_sd, ref.zero = T, fun=exp)) +
  geom_hline(aes(yintercept = 1), color="grey") +
  geom_abline(intercept = 0, 
              slope = HR_BA_cancer %>% filter(BA=="kdm_res_sd", outcome=="Colorectal cancer", model=="Colorectal cancer-specific model") %>% select(beta) %>% as.numeric(),
              linetype="dashed") +
  geom_line(aes(x=kdm_res_sd, y=yhat)) +
  annotate("text", x=1, y=.55, size=3.5, label="P for non-linearity=.025") + 
  scale_x_continuous(breaks=seq(-3,3,1), limits=c(-3,3), labels=c("-3SD","-2SD","-1SD","Mean","+1SD","+2SD","+3SD")) +
  scale_y_continuous(breaks=c(.5,1,1.5,2,2.5), trans=scales::log_trans()) +
  coord_cartesian(ylim=c(.5,2.5)) + 
  labs(title="Colorectal cancer",x="KDM residual", y="Hazard Ratio") + 
  theme_bw() + 
  theme(plot.title = element_text(size=11, face="bold"), 
        plot.caption = element_blank(), 
        axis.title = element_text(face="bold", size=9),
        axis.text = element_text(size=7))

### PhenoAge
dd <- datadist(ukb_bioage_cancer); dd$limits$phenoage_res_sd[2] <- 0; options(datadist="dd")
fit_cox_cc_phenoage_spline <- cph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ 
                                    rcs(phenoage_res_sd,4) +
                                    birth_year_cat+sex+assessment_centre+ethnicity_4group+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group+family_history_bowel_cancer+cc_screening+fresh_veg_cat+red_meat_cat+processed_meat_cat,
                                  data=ukb_bioage_cancer)
anova(fit_cox_cc_phenoage_spline) # p=.011
plot_cox_cc_phenoage_spline <- 
  ggplot(data = Predict(fit_cox_cc_phenoage_spline, phenoage_res_sd, ref.zero = T, fun=exp)) +
  geom_hline(aes(yintercept = 1), color="grey") +
  geom_abline(intercept = 0, 
              slope = HR_BA_cancer %>% filter(BA=="phenoage_res_sd", outcome=="Colorectal cancer", model=="Colorectal cancer-specific model") %>% select(beta) %>% as.numeric(),
              linetype="dashed") +
  geom_line(aes(x=phenoage_res_sd, y=yhat)) +
  annotate("text", x=1, y=.55, size=3.5, label="P for non-linearity=.011") + 
  scale_x_continuous(breaks=seq(-3,3,1), limits=c(-3,3), labels=c("-3SD","-2SD","-1SD","Mean","+1SD","+2SD","+3SD")) +
  scale_y_continuous(breaks=c(.5,1,1.5,2,2.5), trans=scales::log_trans()) +
  coord_cartesian(ylim=c(.5,2.5)) + 
  labs(title="",x="PhenoAge residual", y="Hazard Ratio") + 
  theme_bw() + 
  theme(plot.title = element_text(size=11, face="bold"), 
        plot.caption = element_blank(), 
        axis.title = element_text(face="bold", size=9),
        axis.text = element_text(size=7))

### HD
dd <- datadist(ukb_bioage_cancer); dd$limits$hd_log_sd[2] <- 0; options(datadist="dd")
fit_cox_cc_hd_spline <- cph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ 
                              rcs(hd_log_sd,4) +
                              birth_year_cat+sex+assessment_centre+ethnicity_4group+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group+family_history_bowel_cancer+cc_screening+fresh_veg_cat+red_meat_cat+processed_meat_cat,
                            data=ukb_bioage_cancer)
anova(fit_cox_cc_hd_spline) # p=.61
plot_cox_cc_hd_spline <- 
  ggplot(data = Predict(fit_cox_cc_hd_spline, hd_log_sd, ref.zero = T, fun=exp)) +
  geom_hline(aes(yintercept = 1), color="grey") +
  geom_abline(intercept = 0, 
              slope = HR_BA_cancer %>% filter(BA=="hd_log_sd", outcome=="Colorectal cancer", model=="Colorectal cancer-specific model") %>% select(beta) %>% as.numeric(),
              linetype="dashed") +
  geom_line(aes(x=hd_log_sd, y=yhat)) +
  annotate("text", x=1, y=.55, size=3.5, label="P for non-linearity=.61") + 
  scale_x_continuous(breaks=seq(-3,3,1), limits=c(-3,3), labels=c("-3SD","-2SD","-1SD","Mean","+1SD","+2SD","+3SD")) +
  scale_y_continuous(breaks=c(.5,1,1.5,2,2.5), trans=scales::log_trans()) +
  coord_cartesian(ylim=c(.5,2.5)) + 
  labs(title="",x="HD (log)", y="Hazard Ratio") + 
  theme_bw() + 
  theme(plot.title = element_text(size=11, face="bold"), 
        plot.caption = element_blank(), 
        axis.title = element_text(face="bold", size=9),
        axis.text = element_text(size=7))


###### Melanoma ######

### KDM
dd <- datadist(ukb_bioage_cancer); dd$limits$kdm_res_sd[2] <- 0; options(datadist="dd")
fit_cox_mel_kdm_spline <- cph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ 
                                rcs(kdm_res_sd,4) +
                                birth_year_cat+sex+assessment_centre+ethnicity_4group+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group+PA_cat+time_outdoors_summer_4group+childhood_sunburn+solarium_use+skin_tan_ease+skin_color_4group+hair_color_3group+sun_protection,
                              data=ukb_bioage_cancer)
anova(fit_cox_mel_kdm_spline) # p=0.73
plot_cox_mel_kdm_spline <- 
  ggplot(data = Predict(fit_cox_mel_kdm_spline, kdm_res_sd, ref.zero = T, fun=exp)) +
  geom_hline(aes(yintercept = 1), color="grey") +
  geom_abline(intercept = 0, 
              slope = HR_BA_cancer %>% filter(BA=="kdm_res_sd", outcome=="Melanoma", model=="Melanoma-specific model") %>% select(beta) %>% as.numeric(),
              linetype="dashed") +
  geom_line(aes(x=kdm_res_sd, y=yhat)) +
  annotate("text", x=1, y=.55, size=3.5, label="P for non-linearity=.73") + 
  scale_x_continuous(breaks=seq(-3,3,1), limits=c(-3,3), labels=c("-3SD","-2SD","-1SD","Mean","+1SD","+2SD","+3SD")) +
  scale_y_continuous(breaks=c(.5,1,1.5,2,2.5), trans=scales::log_trans()) +
  coord_cartesian(ylim=c(.5,2.5)) + 
  labs(title="Melanoma",x="KDM residual", y="Hazard Ratio") + 
  theme_bw() + 
  theme(plot.title = element_text(size=11, face="bold"), 
        plot.caption = element_blank(), 
        axis.title = element_text(face="bold", size=9),
        axis.text = element_text(size=7))

### PhenoAge
dd <- datadist(ukb_bioage_cancer); dd$limits$phenoage_res_sd[2] <- 0; options(datadist="dd")
fit_cox_mel_phenoage_spline <- cph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ 
                                     rcs(phenoage_res_sd,4) +
                                     birth_year_cat+sex+assessment_centre+ethnicity_4group+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group+PA_cat+time_outdoors_summer_4group+childhood_sunburn+solarium_use+skin_tan_ease+skin_color_4group+hair_color_3group+sun_protection,
                                   data=ukb_bioage_cancer)
anova(fit_cox_mel_phenoage_spline) # p=.06
plot_cox_mel_phenoage_spline <- 
  ggplot(data = Predict(fit_cox_mel_phenoage_spline, phenoage_res_sd, ref.zero = T, fun=exp)) +
  geom_hline(aes(yintercept = 1), color="grey") +
  geom_abline(intercept = 0, 
              slope = HR_BA_cancer %>% filter(BA=="phenoage_res_sd", outcome=="Melanoma", model=="Melanoma-specific model") %>% select(beta) %>% as.numeric(),
              linetype="dashed") +
  geom_line(aes(x=phenoage_res_sd, y=yhat)) +
  annotate("text", x=1, y=.55, size=3.5, label="P for non-linearity=.06") + 
  scale_x_continuous(breaks=seq(-3,3,1), limits=c(-3,3), labels=c("-3SD","-2SD","-1SD","Mean","+1SD","+2SD","+3SD")) +
  scale_y_continuous(breaks=c(.5,1,1.5,2,2.5), trans=scales::log_trans()) +
  coord_cartesian(ylim=c(.5,2.5)) + 
  labs(title="",x="PhenoAge residual", y="Hazard Ratio") + 
  theme_bw() + 
  theme(plot.title = element_text(size=11, face="bold"), 
        plot.caption = element_blank(), 
        axis.title = element_text(face="bold", size=9),
        axis.text = element_text(size=7))

### HD
dd <- datadist(ukb_bioage_cancer); dd$limits$hd_log_sd[2] <- 0; options(datadist="dd")
fit_cox_mel_hd_spline <- cph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ 
                               rcs(hd_log_sd,4) +
                               birth_year_cat+sex+assessment_centre+ethnicity_4group+bmi_cat+smoking+alcohol_3group+PA_cat+education+di_5group+PA_cat+time_outdoors_summer_4group+childhood_sunburn+solarium_use+skin_tan_ease+skin_color_4group+hair_color_3group+sun_protection,
                             data=ukb_bioage_cancer)
anova(fit_cox_mel_hd_spline) # p=0.52
plot_cox_mel_hd_spline <- 
  ggplot(data = Predict(fit_cox_mel_hd_spline, hd_log_sd, ref.zero = T, fun=exp)) +
  geom_hline(aes(yintercept = 1), color="grey") +
  geom_abline(intercept = 0, 
              slope = HR_BA_cancer %>% filter(BA=="hd_log_sd", outcome=="Melanoma", model=="Melanoma-specific model") %>% select(beta) %>% as.numeric(),
              linetype="dashed") +
  geom_line(aes(x=hd_log_sd, y=yhat)) +
  annotate("text", x=1, y=.55, size=3.5, label="P for non-linearity=.52") + 
  scale_x_continuous(breaks=seq(-3,3,1), limits=c(-3,3), labels=c("-3SD","-2SD","-1SD","Mean","+1SD","+2SD","+3SD")) +
  scale_y_continuous(breaks=c(.5,1,1.5,2,2.5), trans=scales::log_trans()) +
  coord_cartesian(ylim=c(.5,2.5)) + 
  labs(title="",x="HD (log)", y="Hazard Ratio") + 
  theme_bw() + 
  theme(plot.title = element_text(size=11, face="bold"), 
        plot.caption = element_blank(), 
        axis.title = element_text(face="bold", size=9),
        axis.text = element_text(size=7))



###### Combine plots for non-linear associations ######

### Updated Jan 2023: keep only those with significant p for non-linearity
jpeg("Output/Spline_BioAge_cancer.jpg", width = 14, height = 17, units = "cm", res = 1000, quality = 100)
plot_grid(
  # Row 1: Any cancer
  plot_grid(plot_cox_anycancer_kdm_spline, plot_cox_anycancer_phenoage_spline, nrow = 1, labels = c("A","")) +
    theme(plot.background = element_rect(color = "black")),
  # Row 2: Breast cancer
  #plot_grid(plot_cox_bc_kdm_spline, plot_cox_bc_phenoage_spline, plot_cox_bc_hd_spline, nrow = 1, labels = c("B","","")) +
  #  theme(plot.background = element_rect(color = "black")),
  # Row 3: Prostate cancer
  plot_grid(plot_cox_pc_kdm_spline, plot_cox_pc_hd_spline, nrow = 1, labels = c("B","")) +
    theme(plot.background = element_rect(color = "black")),
  # Row 4: Lung cancer
  plot_grid(plot_cox_lc_kdm_spline, plot_cox_lc_phenoage_spline, nrow = 1, labels = c("C","")) +
    theme(plot.background = element_rect(color = "black")),
  # Row 5: Colorectal cancer
  plot_grid(plot_cox_cc_kdm_spline, plot_cox_cc_phenoage_spline, nrow = 1, labels = c("D","")) +
    theme(plot.background = element_rect(color = "black")),
  # Row 6: Melanoma
  #plot_grid(plot_cox_mel_kdm_spline, plot_cox_mel_phenoage_spline, plot_cox_mel_hd_spline, nrow = 1, labels = c("F","","")) +
  #  theme(plot.background = element_rect(color = "black")),
  nrow = 4, align = "h")
dev.off()



#=============== SENSITIVITY ANALYSIS: EXCLUDE MISSING DATA ===================

# Note: In this analysis, only complete data were used (i.e. individuals with
#       missing data on any covariate were excluded)

# BA measures to be tested
BA_measures_2 <- c("kdm_res_sd", # KDM based on 18 biomarkers
                   "phenoage_res_sd", # PhenoAge based on 18 biomarkers
                   "hd_log_sd") # HD based on 18 biomarkers

# Cox models of BA measures on cancers
HR_BA_anycancer_completecase <- data.frame(BA=BA_measures_2,
                                           outcome="Any cancer",
                                           beta=NA, se=NA, HR=NA, lb=NA, ub=NA, p=NA)
HR_BA_bc_completecase <- data.frame(BA=BA_measures_2,
                                    outcome="Breast cancer in women",
                                    beta=NA, se=NA, HR=NA, lb=NA, ub=NA, p=NA)
HR_BA_pc_completecase <- data.frame(BA=BA_measures_2,
                                    outcome="Prostate cancer in men",
                                    beta=NA, se=NA, HR=NA, lb=NA, ub=NA, p=NA)
HR_BA_lc_completecase <- data.frame(BA=BA_measures_2,
                                    outcome="Lung cancer",
                                    beta=NA, se=NA, HR=NA, lb=NA, ub=NA, p=NA)
HR_BA_cc_completecase <- data.frame(BA=BA_measures_2,
                                    outcome="Colorectal cancer", 
                                    beta=NA, se=NA, HR=NA, lb=NA, ub=NA, p=NA)
HR_BA_mel_completecase <- data.frame(BA=BA_measures_2,
                                     outcome="Melanoma",
                                     beta=NA, se=NA, HR=NA, lb=NA, ub=NA, p=NA)

for (j in 1:length(BA_measures_2)) {
  
  ###### Any cancer ######
  coef <- summary(coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ . , 
                        ukb_bioage_cancer_complete[,c("interview_date","birth_date","censor_date_any_cancer","any_cancer_diag",BA_measures_2[j],cov_ukb_multivariable)]))
  HR_BA_anycancer_completecase[j,c(3,4,8)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_anycancer_completecase[j,5:7] <- coef$conf.int[1,c(1,3,4)]
  
  ###### Breast cancer in women ######
  coef <- summary(coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ . , 
                        ukb_bioage_cancer_complete[ukb_bioage_cancer_complete$sex=="Women",c("interview_date","birth_date","censor_date_bc","bc_diag",BA_measures_2[j],cov_ukb_bc)]))
  HR_BA_bc_completecase[j,c(3,4,8)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_bc_completecase[j,5:7] <- coef$conf.int[1,c(1,3,4)]
  
  ###### Prostate cancer in men ###### 
  coef <- summary(coxph(Surv((interview_date-birth_date)/365.25, (censor_date_pc-birth_date)/365.25, pc_diag) ~ . , 
                        ukb_bioage_cancer_complete[ukb_bioage_cancer_complete$sex=="Men",c("interview_date","birth_date","censor_date_pc","pc_diag",BA_measures_2[j],cov_ukb_pc)]))
  HR_BA_pc_completecase[j,c(3,4,8)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_pc_completecase[j,5:7] <- coef$conf.int[1,c(1,3,4)]
  
  ###### Lung cancer ###### 
  coef <- summary(coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
                        ukb_bioage_cancer_complete[,c("interview_date","birth_date","censor_date_lc","lc_diag",BA_measures_2[j],cov_ukb_lc)]))
  HR_BA_lc_completecase[j,c(3,4,8)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_lc_completecase[j,5:7] <- coef$conf.int[1,c(1,3,4)]
  
  ###### Colorectal cancer ###### 
  coef <- summary(coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ . , 
                        ukb_bioage_cancer_complete[,c("interview_date","birth_date","censor_date_cc","cc_diag",BA_measures_2[j],cov_ukb_cc)]))
  HR_BA_cc_completecase[j,c(3,4,8)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_cc_completecase[j,5:7] <- coef$conf.int[1,c(1,3,4)]
  
  ###### Melanoma ###### 
  coef <- summary(coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ . , 
                        ukb_bioage_cancer_complete[,c("interview_date","birth_date","censor_date_mel","mel_diag",BA_measures_2[j],cov_ukb_mel)]))
  HR_BA_mel_completecase[j,c(3,4,8)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_mel_completecase[j,5:7] <- coef$conf.int[1,c(1,3,4)]
}				
rm(list=c("coef","j"))

### Output results
HR_BA_cancer_completecase <- rbind(HR_BA_anycancer_completecase,HR_BA_bc_completecase,HR_BA_pc_completecase,HR_BA_lc_completecase,HR_BA_cc_completecase,HR_BA_mel_completecase)
write.table(HR_BA_cancer_completecase, file = "Output/Cox_BA_cancers_completecase.txt", sep = "\t", row.names = F, na="")





#======== SENSITIVITY ANALYSIS: EXCLUDE THOSE WITH <2 YEARS FOLLOW-UP =========

# BA measures to be tested
BA_measures_2 <- c("kdm_res_sd", # KDM based on 18 biomarkers
                   "phenoage_res_sd", # PhenoAge based on 18 biomarkers
                   "hd_log_sd") # HD based on 18 biomarkers

# Cox models of BA measures on cancers
HR_BA_anycancer_excl_first_2y <- data.frame(BA=BA_measures_2,
                                            outcome="Any cancer",
                                            beta=NA, se=NA, HR=NA, lb=NA, ub=NA, p=NA)
HR_BA_bc_excl_first_2y <- data.frame(BA=BA_measures_2,
                                     outcome="Breast cancer in women",
                                     beta=NA, se=NA, HR=NA, lb=NA, ub=NA, p=NA)
HR_BA_pc_excl_first_2y <- data.frame(BA=BA_measures_2,
                                     outcome="Prostate cancer in men",
                                     beta=NA, se=NA, HR=NA, lb=NA, ub=NA, p=NA)
HR_BA_lc_excl_first_2y <- data.frame(BA=BA_measures_2,
                                     outcome="Lung cancer",
                                     beta=NA, se=NA, HR=NA, lb=NA, ub=NA, p=NA)
HR_BA_cc_excl_first_2y <- data.frame(BA=BA_measures_2,
                                     outcome="Colorectal cancer", 
                                     beta=NA, se=NA, HR=NA, lb=NA, ub=NA, p=NA)
HR_BA_mel_excl_first_2y <- data.frame(BA=BA_measures_2,
                                      outcome="Melanoma",
                                      beta=NA, se=NA, HR=NA, lb=NA, ub=NA, p=NA)

for (j in 1:length(BA_measures_2)) {
  
  ###### Any cancer ######
  coef <- summary(coxph(Surv((interview_date-birth_date)/365.25, (censor_date_any_cancer-birth_date)/365.25, any_cancer_diag) ~ . , 
                        ukb_bioage_cancer[(ukb_bioage_cancer$censor_date_any_cancer-ukb_bioage_cancer$interview_date)/365.25 >= 2 ,c("interview_date","birth_date","censor_date_any_cancer","any_cancer_diag",BA_measures_2[j],cov_ukb_multivariable)]))
  HR_BA_anycancer_excl_first_2y[j,c(3,4,8)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_anycancer_excl_first_2y[j,5:7] <- coef$conf.int[1,c(1,3,4)]
  
  ###### Breast cancer in women ######
  coef <- summary(coxph(Surv((interview_date-birth_date)/365.25, (censor_date_bc-birth_date)/365.25, bc_diag) ~ . , 
                        ukb_bioage_cancer[(ukb_bioage_cancer$censor_date_any_cancer-ukb_bioage_cancer$interview_date)/365.25 >= 2 & ukb_bioage_cancer$sex=="Women",c("interview_date","birth_date","censor_date_bc","bc_diag",BA_measures_2[j],cov_ukb_bc)]))
  HR_BA_bc_excl_first_2y[j,c(3,4,8)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_bc_excl_first_2y[j,5:7] <- coef$conf.int[1,c(1,3,4)]
  
  ###### Prostate cancer in men ###### 
  coef <- summary(coxph(Surv((interview_date-birth_date)/365.25, (censor_date_pc-birth_date)/365.25, pc_diag) ~ . , 
                        ukb_bioage_cancer[(ukb_bioage_cancer$censor_date_any_cancer-ukb_bioage_cancer$interview_date)/365.25 >= 2 & ukb_bioage_cancer$sex=="Men",c("interview_date","birth_date","censor_date_pc","pc_diag",BA_measures_2[j],cov_ukb_pc)]))
  HR_BA_pc_excl_first_2y[j,c(3,4,8)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_pc_excl_first_2y[j,5:7] <- coef$conf.int[1,c(1,3,4)]
  
  ###### Lung cancer ###### 
  coef <- summary(coxph(Surv((interview_date-birth_date)/365.25, (censor_date_lc-birth_date)/365.25, lc_diag) ~ . , 
                        ukb_bioage_cancer[(ukb_bioage_cancer$censor_date_any_cancer-ukb_bioage_cancer$interview_date)/365.25 >= 2 ,c("interview_date","birth_date","censor_date_lc","lc_diag",BA_measures_2[j],cov_ukb_lc)]))
  HR_BA_lc_excl_first_2y[j,c(3,4,8)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_lc_excl_first_2y[j,5:7] <- coef$conf.int[1,c(1,3,4)]
  
  ###### Colorectal cancer ###### 
  coef <- summary(coxph(Surv((interview_date-birth_date)/365.25, (censor_date_cc-birth_date)/365.25, cc_diag) ~ . , 
                        ukb_bioage_cancer[(ukb_bioage_cancer$censor_date_any_cancer-ukb_bioage_cancer$interview_date)/365.25 >= 2 ,c("interview_date","birth_date","censor_date_cc","cc_diag",BA_measures_2[j],cov_ukb_cc)]))
  HR_BA_cc_excl_first_2y[j,c(3,4,8)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_cc_excl_first_2y[j,5:7] <- coef$conf.int[1,c(1,3,4)]
  
  ###### Melanoma ###### 
  coef <- summary(coxph(Surv((interview_date-birth_date)/365.25, (censor_date_mel-birth_date)/365.25, mel_diag) ~ . , 
                        ukb_bioage_cancer[(ukb_bioage_cancer$censor_date_any_cancer-ukb_bioage_cancer$interview_date)/365.25 >= 2 & ukb_bioage_cancer$sun_protection!="missing"&ukb_bioage_cancer$ethnicity_4group!="missing",c("interview_date","birth_date","censor_date_mel","mel_diag",BA_measures_2[j],cov_ukb_mel)]))
  HR_BA_mel_excl_first_2y[j,c(3,4,8)] <- coef$coefficients[1,c(1,3,5)]
  HR_BA_mel_excl_first_2y[j,5:7] <- coef$conf.int[1,c(1,3,4)]
}				
rm(list=c("coef","j"))

### Output results
HR_BA_cancer_excl_first_2y <- rbind(HR_BA_anycancer_excl_first_2y,HR_BA_bc_excl_first_2y,HR_BA_pc_excl_first_2y,HR_BA_lc_excl_first_2y,HR_BA_cc_excl_first_2y,HR_BA_mel_excl_first_2y)
write.table(HR_BA_cancer_excl_first_2y, file = "Output/Cox_BA_cancers_excl_first_2y.txt", sep = "\t", row.names = F, na="")


# =============================== END OF FILE  ===============================