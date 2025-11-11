/*==============================================================================
FILENAME:		2.1_data_prep_UKB_cancer.do
PURPOSE:		To prepare data for the top 5 cancers in UK Biobank, including breast,
			prostate, lung, colorectal, and melanoma skin cancer
AUTHOR:			Jonathan Mak
CREATED			2021-04-20
UPDATED:		2023-01-13
STATA VERSION:		16.0
==============================================================================*/

cd "Z:/UKB_frailty/"

/*====================== Import UKB cancer data ==============================*/

/// Self-reported cancer variables
import delimited "Data/Raw_data/Cancer/ukb_cancer_self_reported.txt", clear
foreach var of varlist f2000100-f2000105 {
replace `var'="" if `var'=="NA"
}
rename feid eid
save "Data/Raw_data/Cancer/ukb_cancer_self_reported.dta", replace

/// Cancer registries variables
import delimited "Data/Raw_data/Cancer/ukb_cancer_registry_Jun2022.txt", clear
foreach var of varlist f4000500-f40021210 {
replace `var'="" if `var'=="NA"
}
rename feid eid
save "Data/Raw_data/Cancer/ukb_cancer_registry.dta", replace

/// Cancer screening variables
import delimited "Data/Raw_data/Cancer/ukb_cancer_screening.txt", clear
foreach var of varlist f234500-f380930 {
replace `var'="" if `var'=="NA"
}
rename feid eid
* Prostate cancer screening
rename f236500 psa_test
label variable psa_test "Ever had prostate specific antigen (PSA) test"
destring psa_test, replace
replace psa_test=. if psa_test<0
* Colorectal cancer screening
rename f234500 cc_screening
label variable cc_screening "Ever had bowel cancer screening"
destring cc_screening, replace
replace cc_screening=. if cc_screening<0
* Breast cancer screening
rename f267400 bc_screening
label variable bc_screening "Ever had breast cancer screening / mammogram"
destring bc_screening, replace
replace bc_screening=. if bc_screening<0
keep eid psa_test cc_screening bc_screening
save "Data/Raw_data/Cancer/ukb_cancer_screening.dta", replace


/*======================= Create cancer variables ============================*/

/// Import date of UK Biobank baseline interview variable
import delimited "Data/Raw_data/Covariates/ukb_interview_date.txt", clear
gen interview_date=date(f5300,"YMD")
format interview_date %td
drop f5300
rename feid eid
label variable eid "Participant ID"
label variable interview_date "Date of baseline assessment"

/// Combine self-reported and cancer registries data
merge 1:1 eid using "Data/Raw_data/Cancer/ukb_cancer_registry.dta", nogenerate
merge 1:1 eid using "Data/Raw_data/Cancer/ukb_cancer_self_reported.dta", nogenerate
drop if interview_date==.

/// Self-reported cancers at baseline

*** 1. Self-reported breast cancer: code 1002
gen bc_self = 0
foreach var in f2000100 f2000101 f2000102 f2000103 f2000104 f2000105 {
replace bc_self=1 if regexm(`var', "1002")
}

*** 2. Self-reported prostate cancer: code 1044
gen pc_self = 0
foreach var in f2000100 f2000101 f2000102 f2000103 f2000104 f2000105 {
replace pc_self=1 if regexm(`var', "1044")
}

*** 3. Self-reported lung cancer: code 1001, 1027, 1028, 1080
gen lc_self = 0
foreach var in f2000100 f2000101 f2000102 f2000103 f2000104 f2000105 {
replace lc_self=1 if regexm(`var', "1001|1027|1028|1080")
}

*** 4. Self-reported colorectal cancer: code 1020, 1022, 1023
gen cc_self = 0
foreach var in f2000100 f2000101 f2000102 f2000103 f2000104 f2000105 {
replace cc_self=1 if regexm(`var', "1020|1022|1023")
}

*** 5. Self-reported melanoma: code 1059
gen mel_self = 0
foreach var in f2000100 f2000101 f2000102 f2000103 f2000104 f2000105 {
replace mel_self=1 if regexm(`var', "1059")
}

drop f2000100 f2000101 f2000102 f2000103 f2000104 f2000105

/// Diagnosed cancers (from cancer registries)

* Reshape dataset from wide to long
rename (*0) (*)							  // Remove the last "0" digit in variable names
rename f400090 f40009
reshape long f40005 f40006 f40008 f40011 f40012 f40013 f40019 f40021, i(eid interview_date bc_self pc_self lc_self cc_self mel_self f40009) j(instance) 

* Rename variables and remove empty entries
gen date_cancer_diag = date(f40005,"YMD") // Date of cancer diagnosis
format date_cancer_diag %td
drop f40005
rename f40008 age_cancer_diag  			  // Age at cancer diagnosis
rename f40006 icd10_cancer	  			  // ICD-10 codes
rename f40013 icd9_cancer   		      // ICD-9 codes
rename f40011 cancer_histology			  // Histology of cancer tumour
rename f40009 occ_cancer			 	  // Reported occurrences of cancer
rename f40012 behavior_cancer			  // Behaviour of cancer tumour
rename f40019 cancer_format				  // Cancer record format
rename f40021 cancer_origin				  // Cancer record origin
destring age_cancer_diag occ_cancer behavior_cancer cancer_format, replace

bysort eid: drop if instance>0 & date_cancer_diag==. // Remove additional empty rows for each participant

*** 0. Any diagnosed cancer (except non-melanoma skin cancer, ICD-10 C44 / ICD-9 173)
gen any_cancer_diag=.
replace any_cancer_diag=1 if (icd10_cancer!="" | icd9_cancer !="") & !regexm(icd10_cancer, "C44") & !regexm(icd9_cancer, "173")
sort eid date_cancer_diag instance
by eid: gen any_cancer_diag_first = any_cancer_diag if sum(any_cancer_diag<.)==1 & sum(any_cancer_diag[_n-1]<.)==0 // First diagnosis of any cancer
gen date_any_cancer_diag_first = date_cancer_diag if any_cancer_diag_first==1 // Date of first diagnosis of any cancer
format date_any_cancer_diag_first %td
gen age_any_cancer_diag_first = age_cancer_diag if any_cancer_diag_first==1 // Age at first diagnosis of any cancer

*** 1. Diagnosed breast cancer: ICD-10 "C50" / ICD-9 "174"
gen bc_diag=.
replace bc_diag=1 if regexm(icd10_cancer, "C50") | regexm(icd9_cancer, "174")
sort eid date_cancer_diag instance
by eid: gen bc_diag_first = bc_diag if sum(bc_diag<.)==1 & sum(bc_diag[_n-1]<.)==0 // First diagnosis of breast cancer
gen date_bc_diag_first = date_cancer_diag if bc_diag_first==1 // Date of first diagnosis of breast cancer
format date_bc_diag_first %td
gen age_bc_diag_first = age_cancer_diag if bc_diag_first==1 // Age at first diagnosis of breast cancer

*** 2. Diagnosed prostate cancer: ICD-10 "C61" / ICD-9 "185"
gen pc_diag=.
replace pc_diag=1 if regexm(icd10_cancer, "C61") | regexm(icd9_cancer, "185")
sort eid date_cancer_diag instance
by eid: gen pc_diag_first = pc_diag if sum(pc_diag<.)==1 & sum(pc_diag[_n-1]<.)==0 // First diagnosis of prostate cancer
gen date_pc_diag_first = date_cancer_diag if pc_diag_first==1 // Date of first diagnosis of prostate cancer
format date_pc_diag_first %td
gen age_pc_diag_first = age_cancer_diag if pc_diag_first==1 // Age at first diagnosis of prostate cancer

*** 3. Diagnosed lung cancer: ICD-10 "C33, C34" / ICD-9 "162"
gen lc_diag=.
replace lc_diag=1 if regexm(icd10_cancer, "C33|C34") | regexm(icd9_cancer, "162")
sort eid date_cancer_diag instance
by eid: gen lc_diag_first = lc_diag if sum(lc_diag<.)==1 & sum(lc_diag[_n-1]<.)==0 // First diagnosis of lung cancer
gen date_lc_diag_first = date_cancer_diag if lc_diag_first==1 // Date of first diagnosis of lung cancer
format date_lc_diag_first %td
gen age_lc_diag_first = age_cancer_diag if lc_diag_first==1 // Age at first diagnosis of lung cancer

*********** Lung cancer subtypes
by eid: replace cancer_histology = cancer_histology[_n+1] if cancer_histology=="" & date_cancer_diag==date_cancer_diag[_n+1] & lc_diag==1 & lc_diag[_n+1]==1 // Fill in some missing lung cancer histology
* Small cell carcinoma (8041–8045, 8246)
gen lc_small_cell_carcinoma = 1 if lc_diag_first==1 & regexm(cancer_histology, "8041|8042|8043|8044|8045|8246")
* Squamous cell carcinoma (8050–8078, 8083,8084)
gen lc_squmaous_cell_carcinoma = 1 if lc_diag_first==1 & regexm(cancer_histology, "805|806|8070|8071|8072|8073|8074|8075|8076|8077|8078|8083|8084")
* Large cell carcinoma (8012–8031, 8035, 8310)
gen lc_large_cell_carcinoma = 1 if lc_diag_first==1 & regexm(cancer_histology, "8012|8013|8014|8015|8016|8017|8018|8019|802|8030|8031|8035|8310")
* Adenocarcinoma (8140, 8211, 8230,8231, 8255–8260, 8323, 8480–8490, 8550,8551, 8570–8574, 8576)
gen lc_adenocarcinoma = 1 if lc_diag_first==1 & regexm(cancer_histology, "8140|8211|8230|8231|8255|8256|8257|8258|8259|8260|8323|8490|8550|8551|8570|8571|8572|8573|8574|8576")
* Bronchioalveolar carcinoma (8250–8254)
gen lc_bac = 1 if lc_diag_first==1 & regexm(cancer_histology, "8250|8251|8252|8253|8254")

*** 4. Diagnosed colorectal cancer: ICD-10 "C18, C19, C20" / ICD-9 "153, 154"
gen cc_diag=.
replace cc_diag=1 if regexm(icd10_cancer, "C18|C19|C20") | regexm(icd9_cancer, "153|154")
sort eid date_cancer_diag instance
by eid: gen cc_diag_first = cc_diag if sum(cc_diag<.)==1 & sum(cc_diag[_n-1]<.)==0 // First diagnosis of colorectal cancer
gen date_cc_diag_first = date_cancer_diag if cc_diag_first==1 // Date of first diagnosis of colorectal cancer
format date_cc_diag_first %td
gen age_cc_diag_first = age_cancer_diag if cc_diag_first==1 // Age at first diagnosis of colorectal cancer

*** 5. Diagnosed melanoma: ICD-10 "C43" / ICD-9 "172"
gen mel_diag=.
replace mel_diag=1 if regexm(icd10_cancer, "C43") | regexm(icd9_cancer, "172")
sort eid date_cancer_diag instance
by eid: gen mel_diag_first = mel_diag if sum(mel_diag<.)==1 & sum(mel_diag[_n-1]<.)==0 // First diagnosis of melanoma
gen date_mel_diag_first = date_cancer_diag if mel_diag_first==1 // Date of first diagnosis of melanoma
format date_mel_diag_first %td
gen age_mel_diag_first = age_cancer_diag if mel_diag_first==1 // Age at first diagnosis of melanoma


/// Re-structure dataset to one row for one individual
collapse (max) interview_date bc_self pc_self lc_self cc_self mel_self any_cancer_diag_first date_any_cancer_diag_first age_any_cancer_diag_first bc_diag_first date_bc_diag_first age_bc_diag_first pc_diag_first date_pc_diag_first age_pc_diag_first lc_diag_first date_lc_diag_first age_lc_diag_first lc_small_cell_carcinoma lc_squmaous_cell_carcinoma lc_large_cell_carcinoma lc_adenocarcinoma lc_bac cc_diag_first date_cc_diag_first age_cc_diag_first mel_diag_first date_mel_diag_first age_mel_diag_first, by(eid)
rename any_cancer_diag_first any_cancer_diag
rename bc_diag_first bc_diag
rename pc_diag_first pc_diag
rename lc_diag_first lc_diag
rename cc_diag_first cc_diag
rename mel_diag_first mel_diag
label variable eid "Participant ID"
label variable interview_date "Date of UKB baseline interview"
label variable bc_self "Self-reported breast cancer"
label variable pc_self "Self-reported prostate cancer"
label variable lc_self "Self-reported lung cancer"
label variable cc_self "Self-reported colorectal cancer"
label variable mel_self "Self-reported melanoma"
label variable any_cancer_diag "Any diagnosed cancer (except non-melanoma skin cancer, ICD-10 C44)"
label variable date_any_cancer_diag_first "Date of first diagnosis of any cancer (except non-melanoma skin cancer, ICD-10 C44)"
label variable age_any_cancer_diag_first "Age at first diagnosis of any cancer (except non-melanoma skin cancer, ICD-10 C44)"
label variable bc_diag "Diagnosed breast cancer"
label variable date_bc_diag_first "Date of first diagnosis of breast cancer"
label variable age_bc_diag_first "Age at first diagnosis of breast cancer"
label variable pc_diag "Diagnosed prostate cancer"
label variable date_pc_diag_first "Date of first diagnosis of prostate cancer"
label variable age_pc_diag_first "Age at first diagnosis of prostate cancer"
label variable lc_diag "Diagnosed lung cancer"
label variable date_lc_diag_first "Date of first diagnosis of lung cancer"
label variable age_lc_diag_first "Age at first diagnosis of lung cancer"
label variable lc_small_cell_carcinoma "Lung cancer subtype: small cell carcinoma"
label variable lc_squmaous_cell_carcinoma "Lung cancer subtype: squamous cell carcinoma"
label variable lc_large_cell_carcinoma "Lung cancer subtype: large cell carcinoma"
label variable lc_adenocarcinoma "Lung cancer subtype: adenocarcinoma"
label variable lc_bac "Lung cancer subtype: bronchioalveolar carcinoma"
label variable cc_diag "Diagnosed colorectal cancer"
label variable date_cc_diag_first "Date of first diagnosis of colorectal cancer"
label variable age_cc_diag_first "Age at first diagnosis of colorectal cancer"
label variable mel_diag "Diagnosed melanoma"
label variable date_mel_diag_first "Date of first diagnosis of melanoma"
label variable age_mel_diag_first "Age at first diagnosis of melanoma"

replace any_cancer_diag = 0 if any_cancer_diag==.
replace bc_diag = 0 if bc_diag==.
replace pc_diag = 0 if pc_diag==.
replace lc_diag = 0 if lc_diag==.
replace cc_diag = 0 if cc_diag==.
replace mel_diag = 0 if mel_diag==.
replace lc_small_cell_carcinoma = 0 if lc_small_cell_carcinoma==.
replace lc_squmaous_cell_carcinoma = 0 if lc_squmaous_cell_carcinoma==.
replace lc_large_cell_carcinoma = 0 if lc_large_cell_carcinoma==.
replace lc_adenocarcinoma = 0 if lc_adenocarcinoma==.
replace lc_bac = 0 if lc_bac==.

/// Save dateset for further analysis
save "Data/Cleaned_data/UKB_cancer.dta", replace


/*=============================== END OF FILE  ===============================*/
