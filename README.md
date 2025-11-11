# UKB_bioage_cancer

This document provides explanations on the codes used in the paper "Clinical biomarker-based biological aging and risk of cancer in the UK Biobank".

Jonathan K. L. Mak, Department of Medical Epidemiology and Biostatistics, Karolinska Institutet, Sweden

Email: jonathan.mak@ki.se 


## Citation

Mak, J.K.L., McMurran, C.E., Kuja-Halkola, R. et al. Clinical biomarker-based biological aging and risk of cancer in the UK Biobank. Br J Cancer (2023). https://doi.org/10.1038/s41416-023-02288-w


## Abstract

**Background:** Despite a clear link between aging and cancer, there has been inconclusive evidence on how biological age (BA) may be associated with cancer incidence.

**Methods:** We studied 308,156 UK Biobank participants with no history of cancer at enrolment. Using 18 age-associated clinical biomarkers, we computed three BA measures (Klemera-Doubal method [KDM], PhenoAge, homeostatic dysregulation [HD]) and assessed their associations with incidence of any cancer and five common cancers (breast, prostate, lung, colorectal, and melanoma) using Cox proportional-hazards models.

**Results:** A total of 35,426 incident cancers were documented during a median follow-up of 10.9 years. Adjusting for common cancer risk factors, 1-standard deviation (SD) increment in the age-adjusted KDM (hazard ratio = 1.04, 95% confidence interval = 1.03–1.05), age-adjusted PhenoAge (1.09, 1.07–1.10), and HD (1.02, 1.01–1.03) was significantly associated with a higher risk of any cancer. All BA measures were also associated with increased risks of lung and colorectal cancers, but only PhenoAge was associated with breast cancer risk. Furthermore, we observed an inverse association between BA measures and prostate cancer, although it was attenuated after removing glycated hemoglobin and serum glucose from the BA algorithms.

**Conclusions:** Advanced BA quantified by clinical biomarkers is associated with increased risks of any cancer, lung cancer, and colorectal cancer.



## Codes and documents

#### Program/1_training_BioAge_NHANES_UKB.R

Codes used for training and testing BA algorithms in NHANES, and projecting the algorithms onto UK Biobank. All BA measures were constructed based on the "BioAge" R package (https://github.com/dayoonkwon/BioAge).

#### Program/2.1_data_prep_UKB_cancer.do

Codes used for creating cancer diagnosis variables in UK Biobank.

#### Program/2.2_data_prep_UKB_BioAge.do

Codes used for data preparation (combining biological age measures, cancer outcomes, and covariates) in UK Biobank.

#### Program/3_analysis_UKB_BioAge_cancer.R

Codes used for survival analyses.

#### BioAge_variables.xlsx

Excel containing the list of biomarkers used in UK Biobank (including the Field ID used) and NHANES III.
