#Compare TCGA downloaded data and the survival data from: 
#Almorox, L., Antequera, L., Rojas, I., Herrera, L. J., & Ortuño, F. M. (2024). 
#Gene Expression Analysis for Uterine Cervix and Corpus Cancer Characterization. 
#Genes, 15(3), 312. https://doi.org/10.3390/genes15030312
#Githhub: https://github.com/Almorox/MDPI_Journal_GENES_Uterine_Cancers-Characterization_through_Gene_Expression_Analysis/tree/main/Data_References

### Library ###

library(TCGAbiolinks)
library(survminer)
library(survival)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)
library(dplyr)
library(stringr)

#First load TCGA previous downloaded
#load("~/5_Global_Analisis/5_Survival/Clinical_cesc.RData")

### Load table from paper ###
url <- "https://raw.githubusercontent.com/Almorox/MDPI_Journal_GENES_Uterine_Cancers-Characterization_through_Gene_Expression_Analysis/main/Data_References/CESC_CGCI/CESC_CGCI_clinical.tsv"
clinical_almorox <- read_tsv(url)
#vroom::vroom_write(file = "~/5_Global_Analisis/5_Survival/CESC_FULL_almorox.tsv",clinical_almorox)

#Check only CESC project
clinical_almorox <- clinical_almorox %>% dplyr::filter(clinical_almorox$project_id == "TCGA-CESC")
#Check duplicates
#each column are row duplicated, the only difference are in  "treatment_or_therapy" and "treatment_type" columns
#in this way, I select unic IDs
clinical_almorox_unic <- clinical_almorox %>% distinct(case_submitter_id,.keep_all = TRUE)
dim(clinical_almorox_unic)
#[1] 304 158


### Clean data ###

#Days_to_last_follow_up

#class(clinical_almorox_unic$days_to_last_follow_up)
#[1] "character"
#print(clinical_almorox_unic$days_to_last_follow_up)
# [1] "4"    "2234" "149"  "555"  "5271" "986"  "12"   "'--"  "18"   "109"  "621" 
#[12] "0"    "834"  "828"  "596"  "'--"  "'--"  "'--"  "919"  "1582" "34"   "25"  
#[23] "21"   "'--"  "'--"  "788"  "'--"  "'--"  "1345" "1099" "552"  "'--"  "540" 
#[34] "607"  "473"  "186"  "9"    "241"  "6375" "1800" "'--"  "'--"  "'--"  "'--" 

clinical_almorox_unic <- clinical_almorox_unic %>%
  mutate(days_to_last_follow_up = as.integer(ifelse(grepl("^\\d+$", days_to_last_follow_up), days_to_last_follow_up, NA)))
# Verify
#class(clinical_almorox_unic$days_to_last_follow_up)
#[1] "integer"
#print(clinical_almorox_unic$days_to_last_follow_up)
#[1]    4 2234  149  555 5271  986   12   NA   18  109  621    0  834  828  596   NA
#[17]   NA   NA  919 1582   34   25   21   NA   NA  788   NA   NA 1345 1099  552   NA
#[33]  540  607  473  186    9  241 6375 1800   NA   NA   NA   NA   NA   NA   NA 1263

#Days_to_death

#class(clinical_almorox_unic$days_to_death)
#[1] "character"
#print(clinical_almorox_unic$days_to_death)
#[1] "'--"  "'--"  "'--"  "'--"  "'--"  "'--"  "'--"  "543"  "'--"  "'--"  "'--"  "'--"  "'--"  "'--"  "'--"  "144"  "355" 
#[18] "2052" "'--"  "'--"  "'--"  "'--"  "'--"  "2859" "348"  "'--"  "469"  "1394" "'--"  "'--"  "'--"  "4086" "'--"  "'--" 
#[35] "'--"  "253"  "'--"  "'--"  "'--"  "'--"  "14"   "506"  "570"  "951"  "1245" "861"  "2094" "'--"  "305"  "'--"  "'--" 

clinical_almorox_unic <- clinical_almorox_unic %>%
  mutate(days_to_death = as.integer(ifelse(grepl("^\\d+$", days_to_death), days_to_death, NA)))
# Verify
#class(clinical_almorox_unic$days_to_death)
#[1] "integer"
#print(clinical_almorox_unic$days_to_death)
#[1]   NA   NA   NA   NA   NA   NA   NA  543   NA   NA   NA   NA   NA   NA   NA  144  355 2052   NA   NA   NA   NA   NA 2859
#[25]  348   NA  469 1394   NA   NA   NA 4086   NA   NA   NA  253   NA   NA   NA   NA   14  506  570  951 1245  861 2094   NA
#[49]  305   NA   NA  607   NA   NA 1186   NA   NA   NA   NA   NA   NA  636   NA   NA   NA   NA   NA   NA  642   NA   NA   NA
#[73]  275  879   NA   NA   NA   NA  284   NA   NA   NA  370   NA   NA   NA   NA   NA   NA   NA   NA  523   NA   NA   NA   NA

#Vital status
#class(clinical_almorox_unic$vital_status)
#[1] "character"
clinical_almorox_unic$vital_status <- as.factor(clinical_almorox_unic$vital_status)
class(clinical_almorox_unic$vital_status)
#[1] "factor"

#Stage#
class(clinical_cesc$stage_event_clinical_stage)
#[1] "factor"
clinical_almorox_unic <- clinical_almorox_unic %>%
  mutate(figo_stage = ifelse(grepl("[A-Za-z]", figo_stage), figo_stage, NA))
clinical_almorox_unic$figo_stage <- as.factor(clinical_almorox_unic$figo_stage)
#vroom::vroom_write(file = "~/5_Global_Analisis/5_Survival/CESC_almorox_Cleaned.tsv",clinical_almorox_unic)


### Compare  ###

### Once verified that main columns are in the correct class and the data has corectly NAs
### I'll compare the two data sets: 
### Almorox's paper where accessed on 29 January 2024.
### I downloaded Clinical_cesc on      02 Feb 2025

#Make comparable data frames#
# Identificar IDs compartidos
shared_ids <- intersect(clinical_cesc$bcr_patient_barcode, clinical_almorox_unic$case_submitter_id)
#length(shared_ids) #[1] 304 it means that all in almorox is in clinical_cesc
# So, which IDs in 'clinical_cesc' are not in almorox
unique_cesc_ids <- setdiff(clinical_cesc$bcr_patient_barcode, clinical_almorox_unic$case_submitter_id) 
length(unique_cesc_ids) #[1] 3
unique_cesc_ids <- as.character(unique_cesc_ids)
#In this sence, i will create a clinical_cesc with the same IDS to be comparable
clinical_cesc_comparable <- clinical_cesc[!(clinical_cesc$bcr_patient_barcode %in% unique_cesc_ids), ]

#clinical_cesc_comparable has the same IDs to clinical_almorox_unic
#so, the stage of this shared IDs are the same?
#First, filter where stage is not NA
cesc_filtered_stage <- clinical_cesc_comparable[grepl("Stage", clinical_cesc_comparable$stage_event_clinical_stage), ] #297
almorox_filtered_stage <- clinical_almorox_unic[grepl("Stage", clinical_almorox_unic$figo_stage), ] #297
#Second, compare if ID is equal to stage in both dataframes
colnames(cesc_filtered_stage)[colnames(cesc_filtered_stage) == "bcr_patient_barcode"] <- "case_submitter_id"

merged_data <- merge(cesc_filtered_stage, almorox_filtered_stage, by = "case_submitter_id")
merged_data$stage_event_clinical_stage <- as.character(merged_data$stage_event_clinical_stage)
merged_data$figo_stage <- as.character(merged_data$figo_stage)

merged_data <- merged_data %>% select(case_submitter_id,
                                      stage_event_clinical_stage,
                                      figo_stage,
                                      days_to_last_followup,
                                      days_to_last_follow_up,
                                      days_to_death.x,
                                      days_to_death.y,
                                      vital_status.x,
                                      vital_status.y)
#FIGO
merged_data$figo_stage_match <- ifelse(is.na(merged_data$stage_event_clinical_stage) | 
                                         is.na(merged_data$figo_stage), 
                                       FALSE, 
                                       merged_data$stage_event_clinical_stage == merged_data$figo_stage)

# days_to_follow_match
merged_data$days_to_follow_match <- ifelse(is.na(merged_data$days_to_last_followup) | 
                                             is.na(merged_data$days_to_last_follow_up), 
                                           FALSE, 
                                           merged_data$days_to_last_followup == merged_data$days_to_last_follow_up)
#Death_merged
merged_data$death_merged <- ifelse(is.na(merged_data$days_to_death.x) | 
                                     is.na(merged_data$days_to_death.y), 
                                   FALSE, 
                                   merged_data$days_to_death.x == merged_data$days_to_death.y)
# Vital status
merged_data$status_merged <- ifelse(is.na(merged_data$vital_status.x) | 
                                      is.na(merged_data$vital_status.y), 
                                    FALSE, 
                                    merged_data$vital_status.x == merged_data$vital_status.y)


table(merged_data$figo_stage_match)
#TRUE 
#297                                    #Same samples with all same stage
table(merged_data$days_to_follow_match)
#FALSE  TRUE 
#102   195                              
table(merged_data$death_merged)
#FALSE  TRUE 
#237    60
table(merged_data$status_merged)
#FALSE  TRUE 
#290     7                        ###its suspicious that almost al NAs are Alive

###Conclusion: all data that i download of TCGA is the same that Almorox_data, 
#although not all Almorox_ data is available, in thise sence, il use Almorox data and after that ill use cesc ddata with dead status infered from days to death
##After that ill check perfectly the article and disccuss with my tutor

