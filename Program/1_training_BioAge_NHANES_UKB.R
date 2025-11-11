#==============================================================================
# FILENAME:   1_training_BioAge_NHANES_UKB.R
# PROJECT:    UKB_bioage_cancer
# PURPOSE:    To create BioAge measures in NHANES (training) and project in UK Biobank
# AUTHOR:     Jonathan Mak & Christopher McMurran
# CREATED:    2022-05-23
# UPDATED:    2023-01-13
# R VERSION:  4.1.3
#==============================================================================
### Note:
# The codes for creating biological age (BA) algorithms are based on the BioAge
# package, more information can be found in: https://github.com/dayoonkwon/BioAge

setwd("Z:/UKB_bioage_cancer")

### Required packages
#devtools::install_github("dayoonkwon/BioAge", dependencies = TRUE)
library(BioAge)
library(haven)
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2); library(patchwork); library(grid); library(gridExtra)
library(psych)
library(splines) # for natural spline model to calculate BA residuals
library(ggeffects)


#================== IMPORT RELEVANT UK BIOBANK DATASETS =======================

### Import datasets
ukb_frailty <- read_dta("Data/Raw_UKB_data/UKB_frailty.dta") # Frailty and sociodemographic data
ukb_blood_biomarker <- read.delim("Data/Raw_UKB_data/ukb_blood_biomarkers.txt") # Blood biomarker data
ukb_biomarker_2 <- read.delim("Data/Raw_UKB_data/ukb_biomarkers_2.txt") # Other biomarker data (e.g., pulse, sbp, lymph...)

### Merge datasets
ukb <- merge(ukb_frailty,
             ukb_blood_biomarker %>% rename(., eid=f.eid), all = T) %>% 
  merge(., ukb_biomarker_2 %>% rename(., eid=f.eid), all = T) 
rm(list=c("ukb_frailty","ukb_blood_biomarker","ukb_biomarker_2"))

### Rename biomarkers (more information can be found in the Excel sheet: "Documents/BioAge_variables.xlsx")
biomarker_labels <- read_excel("Documents/BioAge_variables.xlsx", sheet = "UKB_variables") %>% as.data.frame() 
biomarker_labels$field_name <- ifelse(!is.na(biomarker_labels$`Field ID`), paste0("f.",biomarker_labels$`Field ID`,".0.0"), NA)
colnames(ukb)[colnames(ukb) %in% intersect(colnames(ukb), biomarker_labels$field_name)] <-
  as.vector(na.omit(biomarker_labels$`Variable Name`[match(colnames(ukb),biomarker_labels$field_name)]))

### Change units of biomarkers
ukb$lnalp <- log(ukb$alp)
ukb$waist <- ukb$waist_circumference
ukb$fev <- ukb$fev_1000*1000
ukb$bun <- ukb$bun_mmol*2.8
ukb$lnbun <- log(ukb$bun)
ukb$lncreat_umol <- log(ukb$creat_umol)
ukb$ttbl <- ukb$ttbl_umol*.0585
ukb$uap <- ukb$uap_umol*.0168
ukb$rbc <- ukb$rbc_L
ukb$wbc <- ukb$wbc_L
ukb$crp <- ukb$crp_mgL/10
ukb$lncrp <- log(ukb$crp)
ukb$hba1c <- ukb$hba1c_mmol*0.0915 + 2.15
ukb$lnhba1c <- log(ukb$hba1c)
ukb$hdl <- ukb$hdl_mmol*38.67
ukb$ldl <- ukb$ldl_mmol*38.67
ukb$trig <- ukb$trig_mmol*88.57
ukb$totchol <- ukb$totchol_mmol*38.67

### Keep important variables
ukb <- ukb %>% 
  select("eid","bmi","sex","age_baseline","assessment_centre",
         "ethnicity","smoking","alcohol","education","di",
         "waist_circumference","hip_circumference","whr",
         "weight_loss","exhaustion","low_activity","slow",
         "low_grip","fp_total","fp","fp_2group","tv_fi",
         "bmi","waist","fev_1000","albumin_gL","alp","bun","creat_umol",
         "glucose_mmol","ttbl","uap","lymph","mcv","monopa",
         "rbc","rdw","wbc","crp","dbp","sbp",
         "pulse","hba1c","hdl","trig","totchol"
         #"basopa","eosnpa","neut","cyst","ggt","ldl"
         ) %>% 
  rename(age=age_baseline, # Rename variables to be consistent with NHANES
         sampleID=eid,
         gender=sex) %>%
  mutate(gender=case_when(gender==1~1,gender==0~2)) %>% # Recode gender to be consistent with NHANES
  filter(!is.na(age) & !is.na(gender))

### Exclude individuals withdrawn from UKB (Updated 25 May 2022)
withdrawn_id <- read.delim("Documents/withdraw22224_233.txt", header=F) # List of participants withdrawn from UKB
ukb <- ukb[!(ukb$sampleID %in% withdrawn_id[,1]),] # Exclude 208 participants withdrawn from UKB
rm(withdrawn_id)



#==================== SELECTION OF BIOMARKERS IN NHANES3 ======================

### Notes:
# Biomarkers to be included in BA algorithms are selected based on:
# (1) Correlation with chronological age (usually r>.1 is used in the literature)
# (2) Data availability (we selected those with missingness <20%, which are mostly
#     biomarkers that are routinely collected in clinical practice)


### Summary of variables
describe(BioAge::NHANES3)
data.frame(missing_percent=sapply(NHANES3, function(x) round(mean(is.na(x))*100, 1))) # proportion of missing data
names(which(colMeans(is.na(NHANES3)) > 0.2)) # several biomarkers have high missing, some with even 100% missing -> those missing >20% were not selected in the algorithm

### Select biomarkers to use in NHANES3 (select those with absolute correlation >0.1)
candidate_biomarkers_nhanes3 <- c("bmi","waist","fev_1000",
                                  "albumin_gL","alp","bun","creat_umol","glucose_mmol",
                                  "ttbl","uap","lymph","mcv","monopa",
                                  "rbc","rdw","wbc","crp","dbp","sbp",
                                  "pulse","hba1c","hdl","trig","totchol"
                                  #"grip_r","grip_l","basopa","eosnpa","neut","cyst","ggt","ldl" # Exclude the ones with high missingness
)

# Make male and female datasets
NHANES3_m <- filter(NHANES3, gender == 1)
NHANES3_f <- filter(NHANES3, gender == 2)

# Calculate Pearson correlations for all, male, and female datasets
pearson_m_nhanes <- data.frame(BM=candidate_biomarkers_nhanes3, r="", Selected=FALSE)
pearson_f_nhanes <- data.frame(BM=candidate_biomarkers_nhanes3, r="", Selected=FALSE)
pearson_all_nhanes <- data.frame(BM=candidate_biomarkers_nhanes3, r="", Selected=FALSE)

for (i in 1:length(candidate_biomarkers_nhanes3)){
  BM <- candidate_biomarkers_nhanes3[i]
  
  # All
  r_all <- cor(NHANES3[,BM], NHANES3$age, use="pairwise.complete.obs")
  pearson_all_nhanes[pearson_all_nhanes$BM==BM, "r"] <- r_all %>% round(digits=2)
  if(abs(r_all)>0.1 & !is.na(r_all)){pearson_all_nhanes[pearson_all_nhanes$BM==BM, "Selected"] <- TRUE}
  
  # Male
  r_m <- cor(NHANES3_m[,BM], NHANES3_m$age, use="pairwise.complete.obs")
  pearson_m_nhanes[pearson_m_nhanes$BM==BM, "r"] <- r_m %>% round(digits=2)
  if(abs(r_m)>0.1 & !is.na(r_m)){pearson_m_nhanes[pearson_m_nhanes$BM==BM, "Selected"] <- TRUE}
  
  # Female  
  r_f <- cor(NHANES3_f[,BM], NHANES3_f$age, use="pairwise.complete.obs")
  pearson_f_nhanes[pearson_f_nhanes$BM==BM, "r"] <- r_f %>% round(digits=2)
  if(abs(r_f)>0.1 & !is.na(r_f)){pearson_f_nhanes[pearson_f_nhanes$BM==BM, "Selected"] <- TRUE}
}
rm(list=c("BM","r_all","r_m","r_f","i","NHANES3_m","NHANES3_f"))

biomarkers_all_nhanes <- data.frame(BM=candidate_biomarkers_nhanes3[pearson_all_nhanes$Selected],
                                    label=biomarker_labels$Description[match(candidate_biomarkers_nhanes3[pearson_all_nhanes$Selected],biomarker_labels$`Variable Name`)])
pearson_all_nhanes[pearson_all_nhanes$BM %in% biomarkers_all_nhanes[,1],] # 19 biomarkers have r>0.1 with age

### Correlations between the 19 biomarkers
round(cor(NHANES3[,biomarkers_all_nhanes[,1]],use="pairwise.complete.obs"), 2)

### Plots to compare biomarkers from UKB (red) with NHANES3 of age 37-73 (age range of UKB, but populations not exactly equivalent)
BM_plot_list <- list()
for (i in 1:dim(biomarkers_all_nhanes)[1]){
BM_plot_list[[i]] <- ggplot() + 
  theme_classic() + 
  geom_density(data=NHANES3[NHANES3$age>=37&NHANES3$age<=73,], aes_string(biomarkers_all_nhanes[i,1])) + 
  geom_density(data=ukb, aes_string(biomarkers_all_nhanes[i,1]), fill="red", alpha=0.3) + 
  labs(x = biomarkers_all_nhanes[i,2], y = "Density")
}
jpeg(file="Output/Densityplot_UKB_NHANES3_included_biomarkers.jpg", width = 32, height = 30, units="cm", res = 300, quality=100)
grid.arrange(grobs=BM_plot_list,ncol=4)
dev.off()
rm(list=c("BM_plot_list","i"))

# Note: Seems that the distribution of pulse is relatively lower in NHANES3 (and 
#       it also doesn't work when including pulse in the later model). Therefore,
#       pulse is excluded from the list


biomarkers_all_nhanes <- subset(biomarkers_all_nhanes, biomarkers_all_nhanes$BM!="pulse")

# Note: A total of 18 biomarkers are included to create BA measures:
# 1. "waist" - Waist circumference (cm)
# 2. "fev_1000" - Forced Expiratory Volume in the first 1.0 second (L)
# 3. "albumin_gL" - Albumin (g/dL)
# 4. "alp" - Alkaline phosphatase (U/L)
# 5. "bun" - Blood urea nitrogen (mg/dL)
# 6. "creat_umol" - Creatinine (umol/L)
# 7. "glucose_mmol" - Serum glucose (mmol/L)
# 8. "uap" - Uric acid (mg/dL)
# 9. "lymph" - Lymphocyte percent (%)
# 10. "mcv" - Mean cell volume (fL)
# 11. "rbc" - Red blood cell count (million cells/uL)
# 12. "rdw" - Red cell distribution width (%)
# 13. "crp" - C-reactive protein (mg/dL)
# 14. "dbp" - Diastolic blood pressure (mm Hg)
# 15. "sbp" - Systolic blood pressure (mm Hg)
# 16. "hba1c" - Glycohemoglobin (%)
# 17. "trig" - Triglyceride (mg/dL)
# 18. "totchol" - Total cholesterol (mg/dL)



#===================== TRANING AND TESTING BA MEASURES =======================

# Note: BA algorithms (KDM, PhenoAge, HD) were trained and tested in 3 steps:
#  (1) Training in NHANES III data
#  (2) Testing in NHANES IV data (i.e. compare new algorithms with published versions)
#  (3) Projecting the algorithms onto other datasets


######### Step 1: Training BA algorithms in NHANES III #########


### Train KDM bioage in NHANES III (separate training for men and women)

# Descriptive statistics of the dataset used for training of KDM
NHANES3 %>% 
  filter(age >= 30 & age <= 75 & pregnant == 0) %>% # Reference sample is NHANES III nonpregnant participants aged 30-75 years
  select(sampleID, year, wave, gender, age, biomarkers_all_nhanes[,1]) %>%
  na.omit() %>%
  describe()

# KDM algorithm based on 18 newly selected clinical biomarkers
kdm_nhanes_trained <- kdm_nhanes(biomarkers=biomarkers_all_nhanes[,1])

# KDM algorithm using original biomarkers in Levine et al (PMID: 23213031)
kdm_nhanes_trained_original <- kdm_nhanes(biomarkers=c("fev_1000","sbp","bun","hba1c","totchol","creat_umol","albumin_gL","alp","crp"))

# For sensitivity analysis for KDM, exclude HbA1c and serum glucose
kdm_nhanes_trained_noglu <- kdm_nhanes(biomarkers=biomarkers_all_nhanes[-c(7,16),1])


### Train PhenoAge in NHANES III

# Descriptive statistics of the dataset used for training of PhenoAge
NHANES3 %>% 
  filter(age >= 20 & age <= 84) %>%
  select(sampleID, year, wave, gender, age, biomarkers_all_nhanes[,1]) %>%
  na.omit() %>%
  describe()

# PhenoAge based on 18 newly selected clinical biomarkers
phenoage_nhanes_trained <- phenoage_nhanes(biomarkers=biomarkers_all_nhanes[,1])

# PhenoAge using original biomarkers in Levine et al
phenoage_nhanes_trained_original <- phenoage_nhanes(biomarkers=c("creat_umol","glucose_mmol","rdw","albumin_gL","alp","mcv","lymph","crp","wbc"))

# Sensitivity analysis for PhenoAge, exclude HbA1c and serum glucose
phenoage_nhanes_trained_noglu <- phenoage_nhanes(biomarkers=biomarkers_all_nhanes[-c(7,16),1])


### Train HD in NHANES III

# Descriptive statistics of the dataset used for training of KDM
NHANES3 %>% 
  filter(age >= 20 & age <= 30 & pregnant == 0 & bmi < 30) %>%
  mutate(albumin = ifelse(albumin >= 3.5 & albumin <= 5, albumin, NA),
         albumin_gL = ifelse(is.na(albumin), NA, albumin_gL),
         alp = ifelse(gender == 2, ifelse(alp >= 37 & alp <= 98, alp, NA), ifelse(alp >= 45 & alp <= 115, alp, NA)),
         lnalp = ifelse(is.na(alp), NA, lnalp),
         bap = ifelse(gender == 2, ifelse(bap <= 14, bap, NA), ifelse(bap <= 20, bap, NA)),
         bun = ifelse(gender == 2, ifelse(bun >= 6 & bun <= 21, bun, NA), ifelse(bun >= 8 & bun <= 24, bun, NA)),
         lnbun = ifelse(is.na(bun), NA, lnbun),
         creat = ifelse(gender == 2, ifelse(creat >= 0.6 & creat <= 1.1, creat, NA), ifelse(creat >= 0.8 & creat <= 1.3, creat, NA)),
         creat_umol = ifelse(is.na(creat), NA, creat_umol),
         lncreat = ifelse(is.na(creat), NA, lncreat),
         lncreat_umol = ifelse(is.na(creat), NA, lncreat_umol),
         glucose = ifelse(glucose >= 60 & glucose <= 100, glucose, NA),
         glucose_mmol = ifelse(is.na(glucose), NA, glucose_mmol),
         glucose_fasting = ifelse(glucose_fasting >= 65 & glucose_fasting <= 110, glucose_fasting, NA),
         ttbl = ifelse(ttbl >= 0.1 & ttbl <= 1.4, ttbl, NA),
         uap = ifelse(gender == 2, ifelse(uap >= 2 & uap <= 7, uap, NA), ifelse(uap >= 2.1 & uap <= 8.5, uap, NA)),
         lnuap = ifelse(is.na(uap), NA, lnuap),
         basopa = ifelse(basopa >= 0 & basopa <= 2, basopa, NA),
         eosnpa = ifelse(eosnpa >=1 & eosnpa <= 7, eosnpa, NA),
         mcv = ifelse(gender == 2, ifelse(mcv >= 78 & mcv <= 101, mcv, NA), ifelse(mcv >= 82 & mcv <= 102, mcv, NA)),
         monopa = ifelse(monopa >= 3 & monopa <= 10, monopa, NA),
         neut = ifelse(neut >= 45 & neut <= 74, neut, NA),
         rbc = ifelse(gender == 2, ifelse(rbc >= 3.5 & rbc <= 5.5, rbc, NA), ifelse(rbc >= 4.2 & rbc <= 6.9, rbc, NA)),
         rdw = ifelse(rdw >= 11.5 & rdw <= 14.5, rdw, NA),
         cadmium = ifelse(cadmium >= 2.7 & cadmium <= 10.7, cadmium, NA),
         crp = ifelse(crp < 2, crp, NA),
         crp_cat = ifelse(is.na(crp), NA, crp_cat),
         lncrp = ifelse(is.na(crp), NA, lncrp),
         cyst = ifelse(cyst >= 0.51 & cyst <= 0.98, cyst, NA),
         ggt = ifelse(gender == 2, ifelse(ggt <= 37.79, ggt, NA), ifelse(ggt <= 55.19, ggt, NA)),
         insulin = ifelse(insulin >= 2.52 & insulin <= 24.1, insulin, NA),
         hba1c = ifelse(hba1c >= 4 & hba1c <= 5.6, hba1c, NA),
         lnhba1c = ifelse(is.na(hba1c), NA, lnhba1c),
         hdl = ifelse(gender == 2, ifelse(hdl >= 40 & hdl <= 86, hdl, NA), ifelse(hdl >= 35 & hdl <= 80, hdl, NA)),
         ldl = ifelse(ldl >= 80 & ldl <= 130, ldl, NA),
         trig = ifelse(trig >= 54 & trig <= 110, trig, NA),
         lymph = ifelse(lymph >= 20 & lymph <= 40, lymph, NA),
         wbc = ifelse(wbc >= 4.5 & wbc <= 11, wbc, NA),
         uap = ifelse(gender == 2, ifelse(uap >= 2.7 & uap <= 6.3, uap, NA), ifelse(uap >= 3.7 & uap <= 8, uap, NA)),
         sbp = ifelse(sbp < 120, sbp, NA),
         dbp = ifelse(dbp < 80, dbp, NA),
         meanbp = ifelse(meanbp < 93.33, meanbp, NA),
         pulse = ifelse(pulse >= 60 & pulse <= 100, pulse, NA),
         totchol = ifelse(totchol < 200, totchol, NA),
         fev = ifelse(fev >= mean(fev, na.rm = TRUE) * 0.8, fev, NA),
         fev_1000 = ifelse(is.na(fev), NA, fev_1000),
         vitaminA = ifelse(vitaminA >= 1.05 & vitaminA <= 2.27, vitaminA, NA),
         vitaminE = ifelse(vitaminE <= 28, vitaminE, NA),
         vitaminB12 = ifelse(vitaminB12 >= 100 & vitaminB12 <= 700, vitaminB12, NA),
         vitaminC = ifelse(vitaminC >= 23 & vitaminC <= 85, vitaminC, NA)) %>%
  select(sampleID, year, wave, gender, age, biomarkers_all_nhanes[,1]) %>%
  na.omit() %>%
  describe()

# HD based on 18 newly selected clinical biomarkers
hd_nhanes_trained <- hd_nhanes(biomarkers=biomarkers_all_nhanes[,1])

# Sensitivity analysis for HD, exclude HbA1c and serum glucose
hd_nhanes_trained_noglu <- hd_nhanes(biomarkers=biomarkers_all_nhanes[-c(7,16),1])



######### Step 2: Testing BA algorithms in NHANES IV #########

### Assemble NHANES IV dataset with projected biological aging measures for analysis
BioAge_nhanes_trained_data <- merge(hd_nhanes_trained$data, kdm_nhanes_trained$data) %>% 
  merge(., phenoage_nhanes_trained$data) %>% 
  merge(., kdm_nhanes_trained_original$data %>% rename(kdm_original=kdm,kdm_advance_original=kdm_advance)) %>% 
  merge(., phenoage_nhanes_trained_original$data %>% rename(phenoage_original=phenoage,phenoage_advance_original=phenoage_advance)) %>%
  merge(., hd_nhanes_trained_noglu$data %>% rename(hd_noglu=hd,hd_log_noglu=hd_log)) %>% 
  merge(., kdm_nhanes_trained_noglu$data %>% rename(kdm_noglu=kdm,kdm_advance_noglu=kdm_advance)) %>% 
  merge(., phenoage_nhanes_trained_noglu$data %>% rename(phenoage_noglu=phenoage,phenoage_advance_noglu=phenoage_advance))

### Create BA residuals using regression model of a natural spline of CA with 3 degrees of freedom
get_BA_resids <- function(BA){
  data = BioAge_nhanes_trained_data %>% drop_na(BA)
  model <- parse(text = sprintf("lm(%s ~ ns(age, df = 3), data = data)", BA)) %>% eval()
  model_predict <- ggpredict(model, terms = c("age"))
  data[,"BA_res"] <- NA
  data[!is.na(data[BA]),"BA_res"] <- resid(model)
  return(residuals(model))
}
for(BA in c("kdm","kdm_original","kdm_noglu","phenoage","phenoage_original","phenoage_noglu")){
  BA_res <- paste0(BA, "_res")
  BioAge_nhanes_trained_data[,BA_res] = NA
  BioAge_nhanes_trained_data[!is.na(BioAge_nhanes_trained_data[BA]),BA_res] <- get_BA_resids(BA)
}
rm(list=c("BA","BA_res"))

summary(BioAge_nhanes_trained_data %>% 
  select("kdm_original","kdm_original_res","kdm","kdm_res","kdm_noglu","kdm_noglu_res",
         "phenoage_original","phenoage_original_res","phenoage","phenoage_res","phenoage_noglu","phenoage_noglu_res",
         "hd","hd_log","hd_noglu","hd_log_noglu"))


### Checking of the modified BioAge measures in NHANES4

# Correlations between BioAge measures and CA
agevar_1 <- c(
  "kdm"="KDM",
  "phenoage"="PhenoAge",
  "hd_log"="HD (log)",
  "kdm_original"="Levine original\nKDM",
  "phenoage_original"="Levine original\nPhenoAge")
jpeg(file="Output/NHANES4_BioAge_CA_correlations.jpg", width = 22, height = 15, units="cm", res = 300, quality=100)
plot_ba(data=BioAge_nhanes_trained_data, agevar=names(agevar_1), label=as.vector(agevar_1))
dev.off()

# Correlations between age residuals and CA
agevar_2 <- c(
  "kdm_res"="KDM residual",
  "phenoage_res"="PhenoAge residual",
  "hd_log" = "HD (log)",
  "kdm_original_res"="Levine original\nKDM residual",
  "phenoage_original_res"="Levine original\nPhenoAge residual",
  "age"="Chronological age")
get_axis_type <- function(labels){
  return(rep("float", length(labels)) %>% setNames(names(labels))) # Create function to generate axis_type variables for BAA plots
}
jpeg(file="Output/NHANES4_BioAgeadvance_CA_correlations.jpg", width = 25, height = 15, units="cm", res = 300, quality=100)
plot_baa(data=BioAge_nhanes_trained_data, agevar=names(agevar_2), label=agevar_2, axis_type=get_axis_type(agevar_2))
dev.off()

# Associations of biological aging measures with mortality
htmltools::save_html(html = table_surv(data=BioAge_nhanes_trained_data %>% drop_na(names(agevar_2[-6])), # Restrict data to be the same for all BA measures
                                       agevar=names(agevar_2[-6]), label=agevar_2[-6]), 
                     file = "Output/NHANES_BioAge_mortality_associations.html")



######### Step 3: Projecting BA algorithms onto UKB #########

#### Missing data pattern in UKB 
data.frame(missing_percent=sapply(ukb[,biomarkers_all_nhanes[,1]], function(x) round(mean(is.na(x))*100, 1))) # proportion of missing data
table(rowSums(is.na(ukb[,biomarkers_all_nhanes[,1]]))) # Number of missing biomarkers per individual; n=331,699 with all the 18 biomarkers available


### Projecting trained BA measures onto UKB (here we excluded UKB participants missing any of the 18 included biomarkers)

# KDM using 18 newly selected biomarkers (separate training for gender)
ukb_kdm_m <- kdm_calc(data = ukb[complete.cases( ukb[, biomarkers_all_nhanes[,1] ]), ] %>% 
                        filter(gender == 1), 
                      biomarkers=biomarkers_all_nhanes[,1],
                      fit = kdm_nhanes_trained$fit$male, 
                      s_ba2 = kdm_nhanes_trained$fit$male$s_b2)$data
ukb_kdm_f <- kdm_calc(data = ukb[complete.cases( ukb[, biomarkers_all_nhanes[,1] ]), ] %>% 
                        filter(gender == 2), 
                      biomarkers=biomarkers_all_nhanes[,1],
                      fit = kdm_nhanes_trained$fit$female, 
                      s_ba2 = kdm_nhanes_trained$fit$female$s_b2)$data
ukb_kdm_all <- rbind(ukb_kdm_m, ukb_kdm_f) # Combine the KDM datasets
rm(list=c("ukb_kdm_m","ukb_kdm_f"))

# Levine original KDM
ukb_kdm_original_m <- kdm_calc(data = ukb[complete.cases( ukb[, biomarkers_all_nhanes[,1] ]), ] %>%
                                 filter(gender == 1), 
                               biomarkers=c("fev_1000","sbp","bun","hba1c","totchol","creat_umol","albumin_gL","alp","crp"),
                               fit = kdm_nhanes_trained_original$fit$male, 
                               s_ba2 = kdm_nhanes_trained_original$fit$male$s_b2)$data
ukb_kdm_original_f <- kdm_calc(data = ukb[complete.cases( ukb[, biomarkers_all_nhanes[,1] ]), ] %>% 
                                 filter(gender == 2), 
                               biomarkers=c("fev_1000","sbp","bun","hba1c","totchol","creat_umol","albumin_gL","alp","crp"),
                               fit = kdm_nhanes_trained_original$fit$female, 
                               s_ba2 = kdm_nhanes_trained_original$fit$female$s_b2)$data
ukb_kdm_original_all <- rbind(ukb_kdm_original_m, ukb_kdm_original_f) %>% # Combine the KDM datasets
  rename(kdm_original=kdm, kdm_advance_original=kdm_advance)
rm(list=c("ukb_kdm_original_m","ukb_kdm_original_f"))

# Sensitivity analysis for KDM, remove HbA1c and serum glucose
ukb_kdm_noglu_m <- kdm_calc(data = ukb[complete.cases( ukb[, biomarkers_all_nhanes[,1] ]), ] %>%
                              filter(gender == 1), 
                            biomarkers=biomarkers_all_nhanes[-c(7,16),1],
                            fit = kdm_nhanes_trained_noglu$fit$male, 
                            s_ba2 = kdm_nhanes_trained_noglu$fit$male$s_b2)$data
ukb_kdm_noglu_f <- kdm_calc(data = ukb[complete.cases( ukb[, biomarkers_all_nhanes[,1] ]), ] %>% 
                              filter(gender == 2), 
                            biomarkers=biomarkers_all_nhanes[-c(7,16),1],
                            fit = kdm_nhanes_trained_noglu$fit$female, 
                            s_ba2 = kdm_nhanes_trained_noglu$fit$female$s_b2)$data
ukb_kdm_noglu_all <- rbind(ukb_kdm_noglu_m, ukb_kdm_noglu_f) %>% # Combine the KDM datasets
  rename(kdm_noglu=kdm, kdm_advance_noglu=kdm_advance)
rm(list=c("ukb_kdm_noglu_m","ukb_kdm_noglu_f"))

# PhenoAge using 18 newly selected biomarkers
ukb_pheno_all <- phenoage_calc(data = ukb[complete.cases( ukb[, biomarkers_all_nhanes[,1] ]), ], 
                               biomarkers=biomarkers_all_nhanes[,1], 
                               fit = phenoage_nhanes_trained$fit)$data
dim(subset(ukb_pheno_all, ukb_pheno_all$phenoage==Inf))[1] # 34 individuals with phenoage=inf
ukb_pheno_all[which(ukb_pheno_all$phenoage==Inf),c("phenoage", "phenoage_advance")] <- NA  # Exclude individuals with infinite PhenoAge

# Levine original PhenoAge
ukb_pheno_original_all <- phenoage_calc(data = ukb[complete.cases( ukb[, biomarkers_all_nhanes[,1] ]), ], 
                                        biomarkers=c("creat_umol","glucose_mmol","rdw","albumin_gL","alp","mcv","lymph","crp","wbc"), 
                                        fit = phenoage_nhanes_trained_original$fit)$data %>%
  rename(phenoage_original=phenoage, phenoage_advance_original=phenoage_advance)
dim(subset(ukb_pheno_original_all, ukb_pheno_original_all$phenoage_original==Inf))[1] # 40 individuals with phenoage=inf
ukb_pheno_original_all[which(ukb_pheno_original_all$phenoage_original==Inf),c("phenoage_original", "phenoage_advance_original")] <- NA  # Exclude individuals with infinite PhenoAge

# Sensitivity analysis for PhenoAge, remove HbA1c and serum glucose
ukb_pheno_noglu_all <- phenoage_calc(data = ukb[complete.cases( ukb[, biomarkers_all_nhanes[,1] ]), ], 
                                     biomarkers=biomarkers_all_nhanes[-c(7,16),1], 
                                     fit = phenoage_nhanes_trained_noglu$fit)$data %>%
  rename(phenoage_noglu=phenoage, phenoage_advance_noglu=phenoage_advance)
dim(subset(ukb_pheno_noglu_all, ukb_pheno_noglu_all$phenoage_noglu==Inf))[1] # 37 individuals with phenoage=inf
ukb_pheno_noglu_all[which(ukb_pheno_noglu_all$phenoage_noglu==Inf),c("phenoage_noglu", "phenoage_advance_noglu")] <- NA  # Exclude individuals with infinite PhenoAge

# HD using 18 newly selected biomarkers [Note: "NHANES3_CLEAN" = reference set for HD (age 20-30, non-pregnant, values in normal range)]
ukb_hd_m <- hd_calc(data = ukb[complete.cases( ukb[, biomarkers_all_nhanes[,1] ]), ] %>% 
                      filter(gender == 1),
                    reference = NHANES3_CLEAN %>% filter(gender == 1),
                    biomarkers=biomarkers_all_nhanes[,1])$data 
ukb_hd_f <- hd_calc(data = ukb[complete.cases( ukb[, biomarkers_all_nhanes[,1] ]), ] %>% 
                      filter(gender == 2),
                    reference = NHANES3_CLEAN %>% filter(gender == 2),
                    biomarkers=biomarkers_all_nhanes[,1])$data
ukb_hd_all <- rbind(ukb_hd_m, ukb_hd_f) # Combine the HD datasets
rm(list=c("ukb_hd_m","ukb_hd_f"))

# Sensitivity analysis for HD, remove HbA1c and serum glucose
ukb_hd_noglu_m <- hd_calc(data = ukb[complete.cases( ukb[, biomarkers_all_nhanes[,1] ]), ] %>% 
                            filter(gender == 1),
                          reference = NHANES3_CLEAN %>% filter(gender == 1),
                          biomarkers=biomarkers_all_nhanes[-c(7,16),1])$data 
ukb_hd_noglu_f <- hd_calc(data = ukb[complete.cases( ukb[, biomarkers_all_nhanes[,1] ]), ] %>% 
                            filter(gender == 2),
                          reference = NHANES3_CLEAN %>% filter(gender == 2),
                          biomarkers=biomarkers_all_nhanes[-c(7,16),1])$data
ukb_hd_noglu_all <- rbind(ukb_hd_noglu_m, ukb_hd_noglu_f) %>% # Combine the HD datasets
  rename(hd_noglu=hd, hd_log_noglu=hd_log)
rm(list=c("ukb_hd_noglu_m","ukb_hd_noglu_f"))



### Merge all biological aging measures
ukb_bioage <- left_join(ukb[complete.cases( ukb[, biomarkers_all_nhanes[,1] ]), ] %>% 
                          select("sampleID","gender","age","tv_fi","fp_total",biomarkers_all_nhanes[,1]), 
                        ukb_hd_all[, c("sampleID", "hd", "hd_log")], by = "sampleID") %>%
  left_join(., ukb_hd_noglu_all[, c("sampleID", "hd_log_noglu")], by = "sampleID") %>%
  left_join(., ukb_kdm_all[, c("sampleID", "kdm")], by = "sampleID") %>%
  left_join(., ukb_kdm_original_all[, c("sampleID", "kdm_original")], by = "sampleID") %>%
  left_join(., ukb_kdm_noglu_all[, c("sampleID", "kdm_noglu")], by = "sampleID") %>%
  left_join(., ukb_pheno_all[, c("sampleID","phenoage")], by = "sampleID") %>%
  left_join(., ukb_pheno_original_all[, c("sampleID", "phenoage_original")], by = "sampleID") %>%
  left_join(., ukb_pheno_noglu_all[, c("sampleID", "phenoage_noglu")], by = "sampleID")
summary(ukb_bioage %>% select(kdm, kdm_original, kdm_noglu, phenoage, phenoage_original, phenoage_noglu, hd_log, hd_log_noglu)) 

rm(list=c("ukb_hd_all","ukb_hd_noglu_all",
          "ukb_kdm_all","ukb_kdm_original_all","ukb_kdm_noglu_all",
          "ukb_pheno_all","ukb_pheno_original_all","ukb_pheno_noglu_all"))


### Exclude extreme outliers (those outside 5SD from mean)
table(abs(as.vector(scale(ukb_bioage$kdm)))>5) # n=55 with KDM outside 5 SD
table(abs(as.vector(scale(ukb_bioage$kdm_original)))>5) # n=32 with KDM (original) outside 5 SD
table(abs(as.vector(scale(ukb_bioage$kdm_noglu)))>5) # n=56 with KDM (sensitivity analysis) outside 5 SD
table(abs(as.vector(scale(ukb_bioage$phenoage)))>5) # n=45 with PhenoAge outside 5 SD
table(abs(as.vector(scale(ukb_bioage$phenoage_original)))>5) # n=79 with PhenoAge (original) outside 5 SD
table(abs(as.vector(scale(ukb_bioage$phenoage_noglu)))>5) # n=48 with PhenoAge (sensitivity analysis) outside 5 SD
table(abs(as.vector(scale(ukb_bioage$hd_log)))>5) # n=62 with HD outside 5 SD
table(abs(as.vector(scale(ukb_bioage$hd_log_noglu)))>5) # n=79 with HD (sensitivity analysis) outside 5 SD

ukb_bioage$kdm <- ifelse(abs(as.vector(scale(ukb_bioage$kdm)))>5, NA, ukb_bioage$kdm)
ukb_bioage$kdm_original <- ifelse(abs(as.vector(scale(ukb_bioage$kdm_original)))>5, NA, ukb_bioage$kdm_original)
ukb_bioage$kdm_noglu <- ifelse(abs(as.vector(scale(ukb_bioage$kdm_noglu)))>5, NA, ukb_bioage$kdm_noglu)
ukb_bioage$phenoage <- ifelse(abs(as.vector(scale(ukb_bioage$phenoage)))>5, NA, ukb_bioage$phenoage)
ukb_bioage$phenoage_original <- ifelse(abs(as.vector(scale(ukb_bioage$phenoage_original)))>5, NA, ukb_bioage$phenoage_original)
ukb_bioage$phenoage_noglu <- ifelse(abs(as.vector(scale(ukb_bioage$phenoage_noglu)))>5, NA, ukb_bioage$phenoage_noglu)
ukb_bioage$hd_log <- ifelse(abs(as.vector(scale(ukb_bioage$hd_log)))>5, NA, ukb_bioage$hd_log)
ukb_bioage$hd_log_noglu <- ifelse(abs(as.vector(scale(ukb_bioage$hd_log_noglu)))>5, NA, ukb_bioage$hd_log_noglu)

summary(ukb_bioage %>% select(kdm, kdm_original, kdm_noglu, phenoage, phenoage_original, phenoage_noglu, hd_log, hd_log_noglu)) 


### Regression model for BioAge, using a natural spline of CA with 3 degrees of freedom (to create age residuals)
# Function to generate BA residuals and plot model used
get_BA_resids <- function(BA){
  data = ukb_bioage %>% drop_na(BA)
  # Basic model = regress on age alone
  model <- parse(text = sprintf("lm(%s ~ ns(age, df = 3), data = data)", BA)) %>% eval()
  model_predict <- ggpredict(model, terms = c("age"))
  data[,"BA_res"] <- NA
  data[!is.na(data[BA]),"BA_res"] <- resid(model)
  return(residuals(model))
}
for(BA in c("kdm","kdm_original","kdm_noglu","phenoage","phenoage_original","phenoage_noglu")){
  BA_res <- paste0(BA, "_res")
  ukb_bioage[,BA_res] = NA
  ukb_bioage[!is.na(ukb_bioage[BA]),BA_res] <- get_BA_resids(BA)
}
rm(list=c("BA","BA_res"))


### Some plots for BAs in UKB

# Density plots of BAs
jpeg(file="Output/UKB_BioAge_density_plot.jpg", width = 22, height = 15, units="cm", res = 300, quality=100)
(ggplot(ukb_bioage, aes(x=kdm, fill=factor(gender, labels=c("Men","Women")))) +
  geom_density(alpha=0.4) + 
  scale_color_brewer(palette="Dark2") + 
  labs(x="KDM", y = "Density") + 
  xlim(0, 100) + 
  theme_bw() + 
  theme(legend.position="right",
        legend.title=element_blank())) + 
  (ggplot(ukb_bioage, aes(x=phenoage, fill=factor(gender, labels=c("Men","Women")))) +
  geom_density(alpha=0.4) + 
  scale_color_brewer(palette="Dark2") + 
  labs(x="PhenoAge", y = "Density") + 
  xlim(0, 100) + 
  theme_bw() + 
  theme(legend.position="right",
        legend.title=element_blank())) +
  (ggplot(ukb_bioage, aes(x=hd_log, fill=factor(gender, labels=c("Men","Women")))) +
     geom_density(alpha=0.4) + 
     scale_color_brewer(palette="Dark2") + 
     labs(x="HD (log)", y = "Density") + 
     theme_bw() + 
     theme(legend.position="right",
           legend.title=element_blank())) +
  (ggplot(ukb_bioage, aes(x=kdm_original, fill=factor(gender, labels=c("Men","Women")))) +
     geom_density(alpha=0.4) + 
     scale_color_brewer(palette="Dark2") + 
     labs(x="Levine original KDM", y = "Density") + 
     xlim(0, 100) + 
     theme_bw() + 
     theme(legend.position="right",
           legend.title=element_blank())) + 
  (ggplot(ukb_bioage, aes(x=phenoage_original, fill=factor(gender, labels=c("Men","Women")))) +
     geom_density(alpha=0.4) + 
     scale_color_brewer(palette="Dark2") + 
     labs(x="Levine original PhenoAge", y = "Density") + 
     xlim(0, 100) + 
     theme_bw() + 
     theme(legend.position="right",
           legend.title=element_blank())) +
  plot_layout(guides = "collect") &  theme(legend.position="right")
dev.off()


# Correlations between BioAge measures and CA

# 1. Correlations between BA and CA
jpeg(file="Output/UKB_BioAge_CA_correlations.jpg", width = 8, height = 15, units="cm", res = 300, quality=100)
plot_ba(data=ukb_bioage ,
        agevar=names(agevar_1[1:3]),
        label=as.vector(agevar_1[1:3])) + 
  scale_color_manual(values = c("#7294D4", "#E6A0C4"), name = "Sex", labels = c("Men", "Women")) +
  facet_wrap(~ method, ncol = 1, scales = "free") +
  theme(legend.position = "top", plot.margin=unit(c(-0.5, 1, 0.5, 0.5), units="line"),
        axis.text = element_text(size = 10), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 11, face = "bold"), axis.line = element_line(colour = "black"),
        legend.title = element_text(size=9, face = "bold"), legend.text = element_text(size = 9),
        strip.background = element_rect(fill = "gray89"), strip.text = element_text(size=10,face = "bold"))
dev.off()

# 2. Correlations between BA residuals
jpeg(file="Output/UKB_BioAgeResiduals_correlations.jpg", width = 18, height = 15, units="cm", res = 300, quality=100)
plot_baa(data=ukb_bioage, 
         agevar=names(agevar_2[c(1:3,6)]), 
         label=agevar_2[c(1:3,6)], 
         axis_type=get_axis_type(agevar_2[c(1:3,6)]))
dev.off()

# 3. Correlations between BA residuals, version 2
jpeg(file="Output/UKB_BioAgeResiduals_correlations_v2.jpg", width = 18, height = 12, units="cm", res = 300, quality=100)
plot_baa(data=ukb_bioage, 
         agevar=names(c("kdm"="KDM",
                        "kdm_res"="KDM\nresidual",
                        "phenoage"="PhenoAge",
                        "phenoage_res"="PhenoAge\nresidual",
                        "hd_log"="HD (log)",
                        "age"="Chronological\nage")), 
         label=c("kdm"="KDM",
                 "kdm_res"="KDM\nresidual",
                 "phenoage"="PhenoAge",
                 "phenoage_res"="PhenoAge\nresidual",
                 "hd_log"="HD (log)",
                 "age"="Chronological\nage"), 
         axis_type=get_axis_type(c("kdm"="KDM",
                                   "kdm_res"="KDM\nresidual",
                                   "phenoage"="PhenoAge",
                                   "phenoage_res"="PhenoAge\nresidual",
                                   "hd_log"="HD (log)",
                                   "age"="Chronological\nage")))
dev.off()



#===================== SAVE DATA FOR FURTHER ANALYSES =========================

write.table(ukb_bioage %>% 
              rename(eid=sampleID) %>% 
              # Select important variables
              select("eid",biomarkers_all_nhanes[,1],"kdm","kdm_res","kdm_original","kdm_original_res","kdm_noglu","kdm_noglu_res","phenoage","phenoage_res","phenoage_original","phenoage_original_res","phenoage_noglu","phenoage_noglu_res","hd","hd_log","hd_log_noglu"), 
            file = "Data/Cleaned_data/UKB_BioAge.txt", sep = "\t", row.names = F, na="")
#save.image("Data/R_data/UKB_BioAge.RData")


# =============================== END OF FILE  ================================