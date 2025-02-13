#Survival analisis first exploration, 
setwd("~/5_Global_Analisis/5_Survival/")
# Cargar librería
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
# Definir la consulta para obtener datos clínicos del proyecto TCGA-CESC
query <- GDCquery(
  project = "TCGA-CESC",
  data.category = "Clinical",
  data.format = "bcr xml"
)

# Descargar los datos
GDCdownload(query)
# Preparar los datos clínicos
clinical_cesc <- GDCprepare_clinic(query, clinical.info = "patient")
#Save (cause all take this results for serious studies)
save.image("~/5_Global_Analisis/5_Survival/Clinical_cesc.RData")

(colnames(clinical_cesc))
#[1] "bcr_patient_barcode"       "additional_studies"        "tissue_source_site"       
#[4] "patient_id"                "bcr_patient_uuid"          "informed_consent_verified"
ncol(clinical_cesc) #[1] 103

head(clinical_cesc)
#BIOINFOMAGICIAN
# script to run survival analysis using TCGA data
# setwd("~/Desktop/demo/survivalAnalysis")

library(TCGAbiolinks)
library(survminer)
library(survival)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)
library(dplyr)
library(stringr)

# getting clinical data for TCGA-BRCA cohort -------------------
colnames(clinical_cesc)
which(colnames(clinical_cesc) %in% c("vital_status", "days_to_last_followup", "days_to_death"))
clinical_cesc[,c(21,22,23)]

# looking at some variables associated with survival 
table(clinical_cesc$vital_status)
#Alive  Dead 
#7     1 
# Mantener los vivos ya identificados y solo inferir los muertos
# Asegurarnos de que los valores de 'vital_status' sean texto, no números
clinical_cesc$vital_status <- as.character(clinical_cesc$vital_status)
clinical_cesc$vital_status[clinical_cesc$vital_status == "1"] <- "Alive"
clinical_cesc$vital_status[clinical_cesc$vital_status == "2"] <- "Dead"

# Mantener los vivos ya identificados y solo inferir los muertos
clinical_cesc$vital_status <- ifelse(is.na(clinical_cesc$vital_status) & !is.na(clinical_cesc$days_to_death) & clinical_cesc$days_to_death > 0, "Dead", clinical_cesc$vital_status)

# Verificar los cambios
table(clinical_cesc$vital_status)

# change certain values the way they are encoded
clinical_cesc$deceased <- ifelse(clinical_cesc$vital_status == "Alive", FALSE, TRUE)
# Actualizar el estado vital según los días de muerte

# Verificar cambios
table(clinical_cesc$deceased)

# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
clinical_cesc$overall_survival <- ifelse(clinical_cesc$vital_status == "Alive",
                                         clinical_cesc$days_to_last_followup,
                                         clinical_cesc$days_to_death)
# get gene expression data -----------

# build a query to get gene expression data for entire cohort
#query_brca_all = GDCquery(
#  project = "TCGA-CESC",
#  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
#  experimental.strategy = "RNA-Seq",
# workflow.type = "STAR - Counts",
#  data.type = "Gene Expression Quantification",
#  sample.type = "Primary Tumor",
#  access = "open")

#output_brca <- getResults(query_brca_all)
# get 20 primary tissue sample barcodes
#tumor <- output_brca$cases[1:20]
# OR
#tumor <- output_brca[output_brca$sample_type == "Primary Tumor", "cases"][1:20]
#tumor

# # get gene expression data from 20 primary tumors 
#query_brca <- GDCquery(
#  project = "TCGA-BRCA",
#  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
#  experimental.strategy = "RNA-Seq",
#  workflow.type = "STAR - Counts",
#  data.type = "Gene Expression Quantification",
#  sample.type = c("Primary Tumor", "Solid Tissue Normal"),
#  access = "open",
#  barcode = tumor)

# download data
#GDCdownload(query_brca)

# get counts
#tcga_brca_data <- GDCprepare(query_brca, summarizedExperiment = TRUE)
#brca_matrix <- assay(tcga_brca_data, "unstranded")
#brca_matrix[1:10,1:10]


# extract gene and sample metadata from summarizedExperiment object
#gene_metadata <- as.data.frame(rowData(tcga_brca_data))
#coldata <- as.data.frame(colData(tcga_brca_data))


# vst transform counts to be used in survival analysis ---------------
# Setting up countData object   
#dds <- DESeqDataSetFromMatrix(countData = brca_matrix,
#                              colData = coldata,
#                              design = ~ 1)

# Removing genes with sum total of 10 reads across all samples
#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]


# vst 
#vsd <- vst(dds, blind=FALSE)
#brca_matrix_vst <- assay(vsd)
#brca_matrix_vst[1:10,1:10]

###Load data ###
Count_matrix_cesc <- vroom::vroom(file = "~/Survival_RAC/1_Get_prepro.R/Final_HPV16_Annot.tsv")
rownames(Count_matrix_cesc) <- Count_matrix_cesc$hgnc_symbol
head(Count_matrix_cesc) #This count matrix has been preprocessed: TPM normaliced and batch effect corrected
#hgnc_symbol `TCGA-C5-A1M5-01A-11R-A13Y-07` `TCGA-EK-A2R9-01A-11R-A18M-07` TCGA-EA-A5O9-01A-11R-A28H…¹
#<chr>                                <dbl>                          <dbl>                       <dbl>
#1 TSPAN6                               4645.                          7318.                      11490.
#2 DPM1                                 2466.                          3319.                       1372.
#3 SCYL3                                 727.                           563.                        537.
#4 FIRRM                                 683.                           523.                       1094.
#5 FGR                                   363.                           436.                        141.
#6 FUCA2                                2308.                          1930.                       2494.

metadata1 <- vroom::vroom(file = "~/Survival_RAC/1_Get_prepro.R/Factors_HPV16.tsv")
rownames(metadata1) <- metadata1$specimenID
library(vroom)

almorox <- vroom(file = "~/5_Global_Analisis/5_Survival/CESC_FULL_almorox.tsv")
metadata <- metadata1 %>%
  left_join(almorox %>% distinct(case_submitter_id, .keep_all = TRUE), 
            by = c("cases.submitter_id" = "case_submitter_id"))

#specimenID                   cases.submitter_id HPV_clade sample_type   HPV_type
#<chr>                        <chr>              <chr>     <chr>         <chr>   
#1 TCGA-C5-A1M5-01A-11R-A13Y-07 TCGA-C5-A1M5       A9        Primary Tumor HPV33   
#2 TCGA-EK-A2R9-01A-11R-A18M-07 TCGA-EK-A2R9       A9        Primary Tumor HPV33   
#3 TCGA-EA-A5O9-01A-11R-A28H-07 TCGA-EA-A5O9       A9        Primary Tumor HPV16   
#4 TCGA-C5-A902-01A-11R-A37O-07 TCGA-C5-A902       A9        Primary Tumor HPV16   
#5 TCGA-C5-A8XH-01A-11R-A37O-07 TCGA-C5-A8XH       A9        Primary Tumor HPV16   
#6 TCGA-C5-A3HL-01A-11R-A213-07 TCGA-C5-A3HL       A9        Primary Tumor HPV16 

#Filter specific genes for survival analisis Rac pathway
#RAC <- Count_matrix_cesc %>% filter(str_detect(hgnc_symbol, "^RAC\\d+$"))
#hgnc_symbol `TCGA-C5-A1M5-01A-11R-A13Y-07` `TCGA-EK-A2R9-01A-11R-A18M-07` TCGA-EA-A5O9-01A-11R-A28H…¹
#<chr>                                <dbl>                          <dbl>                       <dbl>
#1 RAC2                                 2519.                          1539.                       1732.
#2 RAC1                                38718.                         30282.                      41959.
# ℹ 265 more variables: `TCGA-C5-A902-01A-11R-A37O-07` <dbl>, `TCGA-C5-A8XH-01A-11R-A37O-07` <dbl>,
#TIAM <- Count_matrix_cesc %>% filter(str_detect(hgnc_symbol, "TIAM"))
#prex_genes <- Count_matrix_cesc %>%filter(str_detect(hgnc_symbol, "PREX"))
#vav_genes <- Count_matrix_cesc %>%filter(str_detect(hgnc_symbol, "VAV"))
#pak_genes <- Count_matrix_cesc %>%filter(str_detect(hgnc_symbol, "PAK"))
#cof_genes <- Count_matrix_cesc %>% filter(str_detect(hgnc_symbol, "CFL"))
#limk_genes <- Count_matrix_cesc %>% filter(str_detect(hgnc_symbol, "LIMK"))




# Get data for RAC gene and add gene metadata information to it -------------
#CASE SUBMITER  CASE_ID GENE  Count
#FOR RAC
RAC1_expression <- Count_matrix_cesc %>%
  filter(hgnc_symbol == "RAC3") %>%
  pivot_longer(-hgnc_symbol, names_to = "specimenID", values_to = "counts")

metadata_rac1 <- metadata %>%
  left_join(RAC1_expression, by = "specimenID")
head(metadata_rac1)
#specimenID              cases.submitter_id HPV_clade sample_type HPV_type hgnc_symbol RAC1_expression
#<chr>                   <chr>              <chr>     <chr>       <chr>    <chr>                 <dbl>
#1 TCGA-C5-A1M5-01A-11R-A… TCGA-C5-A1M5       A9        Primary Tu… HPV33    RAC1                 38718.
#2 TCGA-EK-A2R9-01A-11R-A… TCGA-EK-A2R9       A9        Primary Tu… HPV33    RAC1                 30282.
#3 TCGA-EA-A5O9-01A-11R-A… TCGA-EA-A5O9       A9        Primary Tu… HPV16    RAC1                 41959.
#4 TCGA-C5-A902-01A-11R-A… TCGA-C5-A902       A9        Primary Tu… HPV16    RAC1                 21488.
#5 TCGA-C5-A8XH-01A-11R-A… TCGA-C5-A8XH       A9        Primary Tu… HPV16    RAC1                 16873.
#6 TCGA-C5-A3HL-01A-11R-A… TCGA-C5-A3HL       A9        Primary Tu… HPV16    RAC1                 36503.


#brca_tp53 <- brca_matrix_vst %>% 
#  as.data.frame() %>% 
#  rownames_to_column(var = 'gene_id') %>% 
#  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
#  left_join(., gene_metadata, by = "gene_id") %>% 
#  filter(gene_name == "TP53")




# get median value
median_value <- median(metadata_rac1$counts)



# denote which cases have higher or lower expression than median count
metadata_rac1$strata <- ifelse(metadata_rac1$counts >= median_value, "HIGH", "LOW")
  
clinical_cesc <- metadata_rac1
# Mantener los vivos ya identificados y solo inferir los muertos
# Asegurarnos de que los valores de 'vital_status' sean texto, no números
clinical_cesc$vital_status <- as.character(clinical_cesc$vital_status)
clinical_cesc$vital_status[clinical_cesc$vital_status == "1"] <- "Alive"
clinical_cesc$vital_status[clinical_cesc$vital_status == "2"] <- "Dead"

# Mantener los vivos ya identificados y solo inferir los muertos
clinical_cesc$vital_status <- ifelse(is.na(clinical_cesc$vital_status) & !is.na(clinical_cesc$days_to_death) & clinical_cesc$days_to_death > 0, "Dead", clinical_cesc$vital_status)

# Verificar los cambios
table(clinical_cesc$vital_status)

# change certain values the way they are encoded
clinical_cesc$deceased <- ifelse(clinical_cesc$vital_status == "Alive", FALSE, TRUE)
# Actualizar el estado vital según los días de muerte

# Verificar cambios
table(clinical_cesc$deceased)

# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
clinical_cesc$overall_survival <- ifelse(clinical_cesc$vital_status == "Alive",
                                         clinical_cesc$days_to_last_follow_up,
                                         clinical_cesc$days_to_death)
# get gene expression d





# Add clinical information to brca_tp53
Full_rac1 <- clinical_cesc # %>% left_join(clinical_cesc, by = c("cases.submitter_id" = "bcr_patient_barcode"))
print(Full_rac1$overall_survival)
#[1] 2052   NA   NA   NA   NA   NA

Full_rac1$overall_survival <- as.numeric(Full_rac1$overall_survival)

# fitting survival curve -----------
fit <- survfit(Surv(overall_survival, deceased) ~ strata, data = Full_rac1)
fit
#261 observations deleted due to missingness 
#n events median 0.95LCL 0.95UCL
#strata=HIGH 4      1   2052      NA      NA
#strata=LOW  3      0     NA      NA      NA

ggsurvplot(fit,
           data = Full_rac1,
           pval = T,
           risk.table = T)



fit2 <- survdiff(Surv(overall_survival, deceased) ~ strata, data = Full_rac1)
#n=7, 261 observations deleted due to missingness.
#N Observed Expected
#strata=HIGH 4        1        1
#strata=LOW  3        0        0
#(O-E)^2/E (O-E)^2/V
#strata=HIGH         0       NaN
#strata=LOW        NaN       NaN
#Chisq= 0  on 0 degrees of freedom, p= 1 
