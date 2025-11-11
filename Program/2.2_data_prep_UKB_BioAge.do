/*==============================================================================
FILENAME:		2.2_data_prep_UKB_BioAge.do
PROJECT:		UKB_bioage_cancer
PURPOSE:		To prepare dataset for testing the association between BioAge measures 
			and incident cancers in UK Biobank
AUTHOR:			Jonathan Mak
CREATED:		2022-07-14
UPDATED:		2023-01-20
STATA VERSION:	16.0
==============================================================================*/

cd "Z:/UKB_BioAge/"

/*========================== Combine UKB datasets ============================*/

/// Import previously saved UKB BioAge variables (see "Program/1_training_BioAge_NHANES_UKB.R")
import delimited "Data/Cleaned_data/UKB_BioAge.txt", clear

/// Combine with the prepared UKB cancer & covariate data (codes for these covariates shown elsewhere; contact jonathan.mak@ki.se for more information)
merge 1:1 eid using "Data/Raw_UKB_data/UKB_cancer.dta", nogenerate keep(matched master) // File containing cancer variables
merge 1:1 eid using "Data/Raw_UKB_data/UKB_frailty.dta", nogenerate keep(matched master) // File containing several covariates, e.g., education, ethnicity, deprivation
merge 1:1 eid using "Data/Raw_UKB_data/UKB_physical_activity.dta", nogenerate keep(matched master) // File containing physical activity variables
merge 1:1 eid using "Data/Raw_UKB_data/UKB_family_history.dta", nogenerate keep(matched master) // File containing family history of cancer variables
merge 1:1 eid using "Data/Raw_UKB_data/UKB_smoking.dta", nogenerate keep(matched master) // File containing smoking variables
merge 1:1 eid using "Data/Raw_UKB_data/UKB_sun_exposure.dta", nogenerate keep(matched master) // File containing sun exposure variables
merge 1:1 eid using "Data/Raw_UKB_data/UKB_cancer_screening.dta", nogenerate keep(matched master) // File containing cancer screening variables
merge 1:1 eid using "Data/Raw_UKB_data/UKB_diet.dta", nogenerate keep(matched master)  // File containing diet variables
merge 1:1 eid using "Data/Raw_UKB_data/UKB_women.dta", nogenerate keep(matched master) // File containing women-specific factor variables, e.g., menopause, parity
merge 1:1 eid using "Data/Raw_UKB_data/UKB_death.dta", nogenerate keep(matched master) // File containing mortality data


/*=========================== Sample selection ===============================*/

/// Exclude those who requested to withdraw from the study prior (Updated 25 May 2022)

/// Exclude those with missing BioAge variables
drop if kdm==. | kdm_res==. | kdm_original==. | kdm_original_res==. | phenoage==. | phenoage_res==. | phenoage_original==. | phenoage_original_res==. | hd==. | hd_log==.

/// Exclude those with a cancer diagnosis before baseline (except non-melanoma skin cancer, ICD-10 C44)
drop if date_any_cancer_diag_first <= interview_date & date_any_cancer_diag_first!=.


/*====================== Recode and label variables ==========================*/

/// Create standardized BioAge measures
egen hd_log_sd = std(hd_log)
egen hd_log_noglu_sd = std(hd_log_noglu)
egen kdm_res_sd = std(kdm_res)
egen kdm_original_res_sd = std(kdm_original_res)
egen kdm_noglu_res_sd = std(kdm_noglu_res)
egen phenoage_res_sd = std(phenoage_res)
egen phenoage_original_res_sd = std(phenoage_original_res)
egen phenoage_noglu_res_sd = std(phenoage_noglu_res)

/// Recode covariate variables (for subsequent analyses)

* Recode BMI variable
gen bmi_cat = .
replace bmi_cat = 1 if bmi < 18.5
replace bmi_cat = 2 if bmi >= 18.5 & bmi < 25
replace bmi_cat = 3 if bmi >= 25 & bmi < 30
replace bmi_cat = 4 if bmi >= 30 & bmi < .
label define bmi 1 "<18.5" 2 "18.5- <25" 3 "25- <30" 4 "30+" 999 "missing", replace
label values bmi_cat bmi
label variable bmi_cat "BMI categories"

* Recode alcohol variable
gen alcohol_3group=1 if alcohol==4|alcohol==5|alcohol==6|alcohol==7
replace alcohol_3group=2 if alcohol==2|alcohol==3
replace alcohol_3group=3 if alcohol==1
label define alcohol_3group 1 "Less than 3 times a month" 2 "1–4 times a week" 3 "Daily or almost daily" 999 "missing", replace
label values alcohol_3group alcohol_3group
label variable alcohol_3group "Alcohol consumption (in 3 groups)"

gen alcohol_2group=0 if alcohol==4|alcohol==5|alcohol==6|alcohol==7
replace alcohol_2group=1 if alcohol==1|alcohol==2|alcohol==3
label define alcohol_2group 0 "Less than weekly" 1 "Weekly" 999 "missing", replace
label values alcohol_2group alcohol_2group
label variable alcohol_2group "Alcohol consumption (in 2 groups)"

* Recode ethnicity variable
gen ethnicity_4group=1 if ethnicity==1
replace ethnicity_4group=2 if ethnicity==2|ethnicity==5
replace ethnicity_4group=3 if ethnicity==3
replace ethnicity_4group=4 if ethnicity==4|ethnicity==6
label define ethnicity_4group 1 "White" 2 "Asian" 3 "Black" 4 "Others" 999 "missing", replace
label values ethnicity_4group ethnicity_4group
label variable ethnicity_4group "Ethnicity (in 4 groups)"

gen ethnicity_2group=1 if ethnicity==1
replace ethnicity_2group=2 if ethnicity>1 & ethnicity<.
label define ethnicity_2group 1 "White" 2 "Non-white", replace
label values ethnicity_2group ethnicity_2group
label variable ethnicity_2group "Ethnicity (in 2 groups)"

* Create Deprivation Index quintile
xtile di_5group = di, n(5)
label define di_5group 1 "1" 2 "2" 3 "3" 4 "4" 5 "5" 999 "missing", replace
label values di_5group di_5group
label variable di_5group "Deprivation index quintiles"

* Create an age indicator of < and >= 60 years
gen age60=.
replace age60=1 if age_baseline<60
replace age60=2 if age_baseline>=60 & age_baseline<.
label define age60 1 "Age <60" 2 "Age 60+", replace
label values age60 age60
label variable age60 "Age < or >= 60 years"

* Create year of birth variable with 10-year interval
gen birth_year_cat = .
replace birth_year_cat=1 if birth_year<1940
replace birth_year_cat=2 if birth_year>=1940 & birth_year<1950
replace birth_year_cat=3 if birth_year>=1950 & birth_year<1960
replace birth_year_cat=4 if birth_year>=1960 & birth_year<.
label define birth_year_cat 1 "1930–1939" 2 "1940–1949" 3 "1950-1959" 4 ">=1960", replace
label values birth_year_cat birth_year_cat
label variable birth_year_cat "Year of birth (in 10-year interval)"

* Self-reported diabetes
rename i19_diabetes diabetes

* Recode missing covariate variables as an "unknown" category
recode bmi_cat alcohol_2group alcohol_3group ethnicity_4group education smoking di_5group PA_cat time_outdoors_summer_4group skin_color_4group hair_color_3group childhood_sunburn solarium_use skin_tan_ease sun_protection psa_test cc_screening bc_screening diabetes fresh_veg_cat red_meat_cat processed_meat_cat menopause hormone_therapy oral_contraceptive parity (missing=999)
label define smoking 0 "Never" 1 "Previous" 2 "Current" 999 "missing", replace
label values smoking smoking
label define education 1 "High" 2 "Intermediate" 3 "Low" 999 "missing", replace
label values education education
label define time_outdoors_summer_4group 0 "<1h/day" 1 "1-2h/day" 2 "3-5h/day" 3 ">5h/day" 999 "missing", replace
label values time_outdoors_summer_4group time_outdoors_summer_4group
label define skin_color_4group 0 "Black, brown" 1 "Light, dark olive" 2 "Fair" 3 "Very fair" 999 "missing", replace
label values skin_color_4group skin_color_4group
label define skin_tan_ease 1 "Get very tanned" 2 "Get moderately tanned" 3 "Get mildly or occasionally tanned" 4 "Never tan, only burn" 999 "missing", replace
label values skin_tan_ease skin_tan_ease
label define hair_color_3group 0 "Black, dark brown, other" 1 "Light brown" 2 "Blonde, red" 999 "missing", replace
label values hair_color_3group hair_color_3group
label define PA_cat 1 "Low" 2 "Moderate" 3 "High" 999 "missing", replace
label values PA_cat PA_cat
label define sun_protection 0 "Never/rarely" 1 "Sometimes" 2 "Most of the time" 3 "Always" 4 "Do not go out in sunshine" 999 "missing", replace
label values sun_protection sun_protection
label define fresh_veg_cat 0 "<5 portions a day" 1 ">=5 portions a day" 999 "missing", replace
label values fresh_veg_cat fresh_veg_cat
label define meat 0 "Less than twice a week" 1 "Twice a week or more" 999 "missing", replace
label values processed_meat_cat red_meat_cat meat
label define menopause 0 "Premenopausal" 1 "Postmenopausal" 999 "missing", replace
label values menopause menopause
label define hormone_therapy 0 "Never" 1 "Ever" 999 "missing", replace
label values hormone_therapy hormone_therapy
label define oral_contraceptive 0 "Never" 1 "Ever" 999 "missing", replace
label values oral_contraceptive oral_contraceptive
label define parity 0 "0" 1 "1 to 2" 2 ">=3" 999 "missing", replace
label values parity parity
label define yesno 0 "No" 1 "Yes" 999 "missing", replace
label values childhood_sunburn solarium_use psa_test cc_screening bc_screening diabetes yesno


/// Variables for survival analysis

* Create a birth date variable
gen birth_date = mdy(birth_month, 15, birth_year)  // Assume all individuals were born on 15th of the month
format birth_date %d 

* Recode cancer variable for those diagnosed after the overall censoring date 2020-02-29
/* Note: (Updated Jul 2022) The overall censoring date for cancer register data were 29 February 2020 for England and Wales, and 31 January 2021 for Scotland. More information can be found in https://biobank.ndph.ox.ac.uk/showcase/exinfo.cgi?src=Data_providers_and_dates */

replace any_cancer_diag = 0 if date_any_cancer_diag_first > date("20200229","YMD") & date_any_cancer_diag_first != .
replace bc_diag = 0 if date_bc_diag_first > date("20200229","YMD") & date_bc_diag_first != .
replace pc_diag = 0 if date_pc_diag_first > date("20200229","YMD") & date_pc_diag_first != .
replace lc_diag = 0 if date_lc_diag_first > date("20200229","YMD") & date_lc_diag_first != .
replace lc_small_cell_carcinoma = 0 if lc_diag==0
replace lc_squmaous_cell_carcinoma = 0 if lc_diag==0
replace lc_large_cell_carcinoma = 0 if lc_diag==0
replace lc_adenocarcinoma = 0 if lc_diag==0
replace lc_bac = 0 if lc_diag==0
replace cc_diag = 0 if date_cc_diag_first > date("20200229","YMD") & date_cc_diag_first != .
replace mel_diag = 0 if date_mel_diag_first > date("20200229","YMD") & date_mel_diag_first != .

* Create a death variables to show whether the participant died during follow-up
gen died=0
replace died=1 if death_date!=. & death_date < date("20200229","YMD")
label variable died "Died during follow-up"

* Create censoring dates
gen censor_date_any_cancer = date("20200229","YMD")
replace censor_date_any_cancer = date_any_cancer_diag_first if any_cancer_diag==1 // Incident cancer
replace censor_date_any_cancer = death_date if death_date < censor_date_any_cancer & any_cancer_diag==0 // Those who died were censored
format censor_date_any_cancer %d 

gen censor_date_bc = date("20200229","YMD")
replace censor_date_bc = date_bc_diag_first if bc_diag==1 // Incident cancer
replace censor_date_bc = death_date if death_date < censor_date_bc & bc_diag==0 // Those who died were censored
format censor_date_bc %d 

gen censor_date_pc = date("20200229","YMD")
replace censor_date_pc = date_pc_diag_first if pc_diag==1 // Incident cancer
replace censor_date_pc = death_date if death_date < censor_date_pc & pc_diag==0 // Those who died were censored
format censor_date_pc %d 

gen censor_date_lc = date("20200229","YMD")
replace censor_date_lc = date_lc_diag_first if lc_diag==1 // Incident cancer
replace censor_date_lc = death_date if death_date < censor_date_lc & lc_diag==0 // Those who died were censored
format censor_date_lc %d 

gen censor_date_cc = date("20200229","YMD")
replace censor_date_cc = date_cc_diag_first if cc_diag==1 // Incident cancer
replace censor_date_cc = death_date if death_date < censor_date_cc & cc_diag==0 // Those who died were censored
format censor_date_cc %d 

gen censor_date_mel = date("20200229","YMD")
replace censor_date_mel = date_mel_diag_first if mel_diag==1 // Incident cancer
replace censor_date_mel = death_date if death_date < censor_date_mel & mel_diag==0 // Those who died were censored
format censor_date_mel %d 



/*================== Save dataset for further analysis =======================*/

keep eid waist fev_1000 albumin_gl alp bun creat_umol glucose_mmol uap lymph mcv rbc rdw crp dbp sbp hba1c trig totchol kdm kdm_res kdm_original kdm_original_res kdm_noglu kdm_noglu_res phenoage phenoage_res phenoage_original phenoage_original_res phenoage_noglu phenoage_noglu_res hd hd_log hd_log_noglu hd_log_sd hd_log_noglu_sd kdm_res_sd phenoage_res_sd kdm_original_res_sd kdm_noglu_res_sd phenoage_original_res_sd phenoage_noglu_res_sd interview_date death_date any_cancer_diag date_any_cancer_diag_first bc_diag date_bc_diag_first pc_diag date_pc_diag_first lc_diag date_lc_diag_first lc_small_cell_carcinoma lc_squmaous_cell_carcinoma lc_large_cell_carcinoma lc_adenocarcinoma lc_bac cc_diag date_cc_diag_first mel_diag date_mel_diag_first sex age_baseline assessment_centre smoking education PA_cat family_history_prostate_cancer family_history_breast_cancer family_history_bowel_cancer family_history_lung_cancer skin_tan_ease time_outdoors_summer_4group skin_color_4group childhood_sunburn hair_color_3group solarium_use sun_protection bmi_cat alcohol_2group alcohol_3group ethnicity_4group ethnicity_2group di di_5group age60 birth_year_cat birth_date died censor_date_any_cancer censor_date_bc censor_date_pc censor_date_lc censor_date_cc censor_date_mel psa_test cc_screening bc_screening diabetes fresh_veg_cat red_meat_cat processed_meat_cat menopause hormone_therapy oral_contraceptive parity 

save "Data/Cleaned_data/UKB_BioAge_cancer.dta", replace


/*=============================== END OF FILE  ===============================*/
