library(devtools)
library(openxlsx)

source_url("https://github.com/beyergroup/Global-analyses-of-the-human-structural-proteome-to-identify-a-new-type-of-disease-biomarker/tree/main/R/Functions_BuildingAndHandelingLinearModels.R?raw=TRUE")

# Loading data
## Peptide and protein quantities
LiPPep <- read.xlsx("Supplementary_Table_12_Data_for_GitHub_scripts.xlsx", "LiPPep_Preprocessed", startRow = 3, rowNames = T)
TrpPep <- read.xlsx("Supplementary_Table_12_Data_for_GitHub_scripts.xlsx", "TrpPep_Preprocessed", startRow = 3, rowNames = T)

TrpProt <- read.xlsx("Supplementary_Table_12_Data_for_GitHub_scripts.xlsx", "TrpProt_Preprocessed", startRow = 3, rowNames = T)
TrpProt2LiPPep <- read.xlsx("Supplementary_Table_12_Data_for_GitHub_scripts.xlsx", "TrpProt_map2LiPPep_Preprocessed", startRow = 3, rowNames = T)
TrpProt2TrpPep <- read.xlsx("Supplementary_Table_12_Data_for_GitHub_scripts.xlsx", "TrpProt_map2TrpPep_Preprocessed", startRow = 3, rowNames = T)

## Further peptide and sample information
PepProtInfo <- read.xlsx("Supplementary_Table_12_Data_for_GitHub_scripts.xlsx", "PepProtInfo", startRow = 3, rowNames = T)
SampleInfo <- read.xlsx("Supplementary_Table_12_Data_for_GitHub_scripts.xlsx", "SampleInfo", startRow = 3, rowNames = T)

# Getting data into right formats for linear models
## Getting peptides/proteins to fit in different models
### Model for structural variation, fitting LiPPep data using only full-tryptic peptides
PepsLiPFT <- Reduce(intersect, list(row.names(LiPPep), 
                                    row.names(TrpPep), 
                                    row.names(TrpProt2TrpPep),
                                    row.names(PepProtInfo)[PepProtInfo$IsTryptic == "full"]))

### Model for structural variation, fitting LiPPep data using only half-tryptic peptides
PepsLiPHT <- Reduce(intersect, list(row.names(LiPPep),  
                                    row.names(TrpProt2LiPPep),
                                    row.names(PepProtInfo)[PepProtInfo$IsTryptic == "half"]))

PepsLiPHT <- PepsLiPHT[!PepsLiPHT %in% row.names(TrpPep)] # removing half-tryptic peptides also measured in trypsin-only data

### Model for PK-independent variation, fitting TrpPep data using only full-tryptic peptides
PepsTrpFT <- Reduce(intersect, list(row.names(TrpPep), 
                                    row.names(TrpProt2TrpPep),
                                    row.names(PepProtInfo)[PepProtInfo$IsTryptic == "full"]))


## Setting sample information into right formats
SampleInfo <- SampleInfo[!is.na(SampleInfo$TPC), ] # removing sample with no corresponding TPC level

FactorCohort <- as.factor(SampleInfo$Cohort) # saving cohort info as factor
SampleInfo$Sex <- as.factor(SampleInfo$Sex) # setting sex info as factor
SampleInfo$Cohort <- as.factor(SampleInfo$Cohort) # setting cohort info as factor
contrasts(SampleInfo$Sex) <- contr.sum(2) # adding contrasts via deviation coding
contrasts(SampleInfo$Cohort) <- contr.sum(2) # adding contrasts via deviation coding

SampleInfo$AgeCohort <- SampleInfo$Age*(as.numeric(SampleInfo$Cohort)-1) # creating age-cohort interactions
SampleInfo$AgeSex <- SampleInfo$Age* sapply(SampleInfo$Sex, function(x) {
    ifelse(x==as.numeric(row.names(contrasts(SampleInfo$Sex)))[1], contrasts(SampleInfo$Sex)[1], contrasts(SampleInfo$Sex)[2])
    }) # creating age-sex interactions
SampleInfo$CohortSex <- sapply(SampleInfo$Sex, function(x){
    ifelse(x==as.numeric(row.names(contrasts(SampleInfo$Sex)))[1], contrasts(SampleInfo$Sex)[1], contrasts(SampleInfo$Sex)[2])
    })*(as.numeric(SampleInfo$Cohort)-1) # creating cohort-sex interactions

SampleInfo$Cohort <- FactorCohort # adding cohort info again, so that contrasts will be modelled using dummy coding

# Running and evaluating models
## Model for structural variation, fitting LiPPep data using only full-tryptic peptides
LM_LiPPepFT <- RunLM(data = LiPPep[PepsLiPFT, row.names(SampleInfo)],
                     dataPep = TrpPep[PepsLiPFT, row.names(SampleInfo)],
                     dataProt = TrpProt2TrpPep[PepsLiPFT, row.names(SampleInfo)],
                     info = SampleInfo, 
                     Age = T, Cohort = T, Sex = T, AgeCohort = T, AgeSex = T, CohortSex = T, Pep = T, Prot = T, TPC = T, Batch = T)

### Only including peptides which could be modelled for at least 80 samples
LM_LiPPepFT <- FilterLM(LM_LiPPepFT, 80)

### Extracting coefficients and p-values for the cohort coefficients
Cohort_LiPPepFT <- ExtractInfoLM(LM_LiPPepFT, "CohortT0")
Cohort_LiPPepFT$proteinwiseFDR <- DoProteinwiseFDR(setNames(Cohort_LiPPepFT$Pvalue, row.names(Cohort_LiPPepFT)), "fdr") ### correcting p-values protein-wise

## Model for structural variation, fitting LiPPep data using only half-tryptic peptides
LM_LiPPepHT <- RunLM(data = LiPPep[PepsLiPHT, row.names(SampleInfo)],
                     dataProt = TrpProt2LiPPep[PepsLiPHT, row.names(SampleInfo)],
                     info = SampleInfo, 
                     Age = T, Cohort = T, Sex = T, AgeCohort = T, AgeSex = T, CohortSex = T, Pep = NULL, Prot = T, TPC = T, Batch = T)

### Only including peptides which could be modelled for at least 80 samples
LM_LiPPepHT <- FilterLM(LM_LiPPepHT, 80)

### Extracting coefficients and p-values for the cohort coefficients
Cohort_LiPPepHT <- ExtractInfoLM(LM_LiPPepHT, "CohortT0")
Cohort_LiPPepHT$proteinwiseFDR <- DoProteinwiseFDR(setNames(Cohort_LiPPepHT$Pvalue, row.names(Cohort_LiPPepHT)), "fdr") ### correcting p-values protein-wise


## Model for Pk-independent variation, fitting TrpPep data using only full-tryptic peptides
LM_TrpPepFT <- RunLM(data = TrpPep[PepsTrpFT, row.names(SampleInfo)],
                     dataProt = TrpProt2TrpPep[PepsTrpFT, row.names(SampleInfo)],
                     info = SampleInfo, 
                     Age = T, Cohort = T, Sex = T, AgeCohort = T, AgeSex = T, CohortSex = T, Pep = NULL, Prot = T, TPC = T, Batch = T)

### Only including peptides which could be modelled for at least 80 samples
LM_TrpPepFT <- FilterLM(LM_TrpPepFT, 80)

### Extracting coefficients and p-values for the cohort coefficients
Cohort_TrpPepFT <- ExtractInfoLM(LM_TrpPepFT, "CohortT0")
Cohort_TrpPepFT$proteinwiseFDR <- DoProteinwiseFDR(setNames(Cohort_TrpPepFT$Pvalue, row.names(Cohort_TrpPepFT)), "fdr") ### correcting p-values protein-wise


## Model for protein abundance variation, fitting TrpProt data
LM_TrpProt <- RunLM(data = TrpProt[, row.names(SampleInfo)],
                    info = SampleInfo, 
                    Age = T, Cohort = T, Sex = T, AgeCohort = T, AgeSex = T, CohortSex = T, Pep = NULL, Prot = NULL, TPC = T, Batch = T)

### Only including peptides which could be modelled for at least 80 samples
LM_TrpProt <- FilterLM(LM_TrpProt, 80)

### Extracting coefficients and p-values for the cohort coefficients
Cohort_TrpProt<- ExtractInfoLM(LM_TrpProt, "CohortT0")

