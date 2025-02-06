# This script normaliced and preprocced data dor TCGA cervical cancer 
#(CC HPV negative) vs (CC HPV-16 positive)

### Load library  ###
library(dplyr)
library(stringr)
library(vroom)

###  Get metadata ###

HPV_IDsample <- vroom('~/0_HPV_Distribution/0_HPV_genotypes.tsv') 

TCGA_HPV <- HPV_IDsample %>% 
  dplyr::select(1,6) %>% 
  dplyr::slice(2:305) %>% 
  dplyr::mutate(
    tipos_HPV = if_else(
      str_detect(HPV_estimated, "HPV16"), "HPV16",
      if_else(HPV_estimated == "negative", "negative", NA_character_)
    )
  ) %>% 
  dplyr::filter(!is.na(tipos_HPV)) %>% 
  arrange(tipos_HPV)

TCGA_HPV <- TCGA_HPV %>% dplyr::select(1,3)
table(TCGA_HPV$tipos_HPV) 
#HPV16 negative 
#168       19 

### Get unstranded matrix (check script 1_Get_RNA_data.R) (is needed Outputquery and unstranded objets)
#load("~/1_Get_Data_TCGA/0_Image_data_RNA.RData")
Output_query_TCGA_unique <- Output_query_TCGA %>%
  distinct(cases.submitter_id, .keep_all = TRUE)

TCGA_HPV_RNAseq <- TCGA_HPV %>%
  left_join(Output_query_TCGA_unique, by = c("TCGA Case ID" = "cases.submitter_id"))
TCGA_HPV_RNAseq <- TCGA_HPV_RNAseq %>% filter(sample_type == "Primary Tumor")
#table(TCGA_HPV_RNAseq$tipos_HPV)
#HPV16 negative 
#165      19
vroom_write(file = "~/Survival_RAC/1_CC_negative_16.tsv", TCGA_HPV_RNAseq)
Cases_CCneg_HPV_16 <- TCGA_HPV_RNAseq %>% pull(cases)
unstranded_counts_hpv16 <- unstranded %>% dplyr::select(gene, Cases_CCneg_HPV_16)
head(unstranded_counts_hpv16)
# gene     TCGA-BI-A0VR-01A-11R…¹ TCGA-BI-A0VS-01A-11R…² TCGA-BI-A20A-01A-11R…³ TCGA-C5-A0TN-01A-21R…⁴
#<chr>                     <int>                  <int>                  <int>                  <int>
# 1 ENSG000…                   5825                   7522                  10173                   3095
#2 ENSG000…                      0                      1                      0                      0
#3 ENSG000…                   3850                   5529                   2067                   1997
#4 ENSG000…                   1076                    784                   1000                    455
#5 ENSG000…                   1246                    943                   1034                    499
#6 ENSG000…                    827                    820                    992                    162

#So i have unstranded counts with CC negative HPV and CC HPV16
#unstranded_counts_hpv16 <- as.tibble(unstranded_counts_hpv16)
#vroom_write(unstranded_counts_hpv16, "~/Survival_RAC/2_unstranded_counts_hpv16.tsv", delim = "\t")

