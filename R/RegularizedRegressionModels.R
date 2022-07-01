library(devtools)
library(openxlsx)

source_url("https://github.com/LuiseNagel/Global-analyses-of-the-human-structural-proteome-to-identify-a-new-type-of-disease-biomarker-/blob/main/R/Functions_BuildingAndHandelingLinearModels.R?raw=TRUE")
source_url("https://github.com/LuiseNagel/Global-analyses-of-the-human-structural-proteome-to-identify-a-new-type-of-disease-biomarker-/blob/main/R/Functions_RegularizedRegression.R?raw=TRUE")


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


# Preparing data for modelling
## Getting peptides/proteins to fit in different models
PepsLiPFT <- Reduce(intersect, list(row.names(LiPPep[!is.na(rowSums(LiPPep)),]), 
                                    row.names(TrpPep), 
                                    row.names(TrpProt2TrpPep),
                                    row.names(PepProtInfo)[PepProtInfo$IsTryptic == "full"]))

PepsTrpFT <- PepsLiPFT

ProtsTrp <- row.names(TrpProt)[!is.na(rowSums(TrpProt))]

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

## Defining 5 runs of different folds
GetFolds <- function(seed){
    set.seed(seed)
    c(sample(rep(1:5, times = c(11, 10, 10, 10, 10))), 
           sample(rep(1:5, times = c(11, 10, 10, 10, 10))))
}
FoldRun1 <- GetFolds(57)
FoldRun2 <- GetFolds(99)
FoldRun3 <- GetFolds(42)
FoldRun4 <- GetFolds(121)
FoldRun5 <- GetFolds(1234)

# Rumnning linear models on all fold runs and estimating scaled residuals for training and test datasets
## Main Models
### Strucutral variation model: LiPPep ~ Age + Cohort + Sex + AgeCohort + AgeSex + CohortSex + TrpPep + TrpProt + TPC + Batch
SR_LiPPep_AllFoldRun <- RunLMandGetScaledResidualsAllFoldRuns(foldList = list(FoldRun1, FoldRun2, FoldRun3, FoldRun4, FoldRun5),
                                                              data = LiPPep[PepsLiPFT, row.names(SampleInfo)],
                                                              dataPep = TrpPep[PepsLiPFT, row.names(SampleInfo)],
                                                              dataProt = TrpProt2TrpPep[PepsLiPFT, row.names(SampleInfo)],
                                                              info = SampleInfo, Age = T, Cohort = T, Sex = T, AgeCohort = T, AgeSex = T, CohortSex = T, Pep = T, Prot = T, TPC = T, Batch = T,
                                                              pwFDR = T)

### PK-independent variation model: TrpPep ~ Age + Cohort + Sex + AgeCohort + AgeSex + CohortSex + TrpProt + TPC + Batch
SR_TrpPep_AllFoldRun <- RunLMandGetScaledResidualsAllFoldRuns(foldList = list(FoldRun1, FoldRun2, FoldRun3, FoldRun4, FoldRun5),
                                                              data = TrpPep[PepsTrpFT, row.names(SampleInfo)],
                                                              dataPep = NULL,
                                                              dataProt = TrpProt2TrpPep[PepsTrpFT, row.names(SampleInfo)],
                                                              info = SampleInfo,
                                                              Age = T, Cohort = T, Sex = T, AgeCohort = T, AgeSex = T, CohortSex = T, Pep = NULL, Prot = T, TPC = T, Batch = T,
                                                              pwFDR = T)

### Protein abundance model: LiPPep ~ Age + Cohort + Sex + AgeCohort + AgeSex + CohortSex + TPC + Batch
SR_TrpProt_AllFoldRun <- RunLMandGetScaledResidualsAllFoldRuns(foldList = list(FoldRun1, FoldRun2, FoldRun3, FoldRun4, FoldRun5),
                                                               data = TrpProt[ProtsTrp, row.names(SampleInfo)],
                                                               dataPep = NULL,
                                                               dataProt = NULL,
                                                               info = SampleInfo,
                                                               Age = T, Cohort = T, Sex = T, AgeCohort = T, AgeSex = T, CohortSex = T, Pep = NULL, Prot = NULL, TPC = T, Batch = T,
                                                               pwFDR = F)


## Variation Models
### Strucutral, PK-independent and abundance variation model: LiPPep ~ Age + Cohort + Sex + AgeCohort + AgeSex + CohortSex + TPC + Batch
SR_LiPPep_noPepProt_AllFoldRun <- RunLMandGetScaledResidualsAllFoldRuns(foldList = list(FoldRun1, FoldRun2, FoldRun3, FoldRun4, FoldRun5),
                                                                        data = LiPPep[PepsLiPFT, row.names(SampleInfo)],
                                                                        dataPep = NULL,
                                                                        dataProt = NULL,
                                                                        info = SampleInfo, Age = T, Cohort = T, Sex = T, AgeCohort = T, AgeSex = T, CohortSex = T, Pep = NULL, Prot = NULL, TPC = T, Batch = T,
                                                                        pwFDR = T)

### Strucutral and abundance variation model: LiPPep ~ Age + Cohort + Sex + AgeCohort + AgeSex + CohortSex + TrpPep + TPC + Batch
SR_LiPPep_noProt_AllFoldRun <- RunLMandGetScaledResidualsAllFoldRuns(foldList = list(FoldRun1, FoldRun2, FoldRun3, FoldRun4, FoldRun5),
                                                                     data = LiPPep[PepsLiPFT, row.names(SampleInfo)],
                                                                     dataPep = TrpPep[PepsLiPFT, row.names(SampleInfo)],
                                                                     dataProt = NULL,
                                                                     info = SampleInfo, Age = T, Cohort = T, Sex = T, AgeCohort = T, AgeSex = T, CohortSex = T, Pep = T, Prot = NULL, TPC = T, Batch = T,
                                                                     pwFDR = T)

### Strucutral and PK-independent variation model: LiPPep ~ Age + Cohort + Sex + AgeCohort + AgeSex + CohortSex + TrpProt + TPC + Batch
SR_LiPPep_noPep_AllFoldRun <- RunLMandGetScaledResidualsAllFoldRuns(foldList = list(FoldRun1, FoldRun2, FoldRun3, FoldRun4, FoldRun5),
                                                                    data = LiPPep[PepsLiPFT, row.names(SampleInfo)],
                                                                    dataPep = NULL,
                                                                    dataProt = TrpProt2TrpPep[PepsLiPFT, row.names(SampleInfo)],
                                                                    info = SampleInfo, Age = T, Cohort = T, Sex = T, AgeCohort = T, AgeSex = T, CohortSex = T, Pep = NULL, Prot = T, TPC = T, Batch = T,
                                                                    pwFDR = T)

# Building LASSO Models based on scaled residuals of training data, testing them on test data and summarizing predictions
## Get Lambdas
DFnLambda <- GetLambdas()

## Run models and extract summarized predictions for the test folds
### Strucutral variation model: LiPPep ~ Age + Cohort + Sex + AgeCohort + AgeSex + CohortSex + TrpPep + TrpProt + TPC + Batch
PredTest_LiPPep_AllFoldRun <- GetTestPredOnFoldRuns(SR_LiPPep_AllFoldRun, 
                                                    info = SampleInfo, 
                                                    nLambda = DFnLambda$LiPPep)
SumPred_LiPPep <- GetMeanPredOfFoldRuns(PredTest_LiPPep_AllFoldRun, SampleInfo)

### PK-independent variation model: TrpPep ~ Age + Cohort + Sex + AgeCohort + AgeSex + CohortSex + TrpProt + TPC + Batch
PredTest_TrpPep_AllFoldRun <- GetTestPredOnFoldRuns(SR_TrpPep_AllFoldRun, 
                                                    info = SampleInfo, 
                                                    nLambda = DFnLambda$TrpPep)
SumPred_TrpPep <- GetMeanPredOfFoldRuns(PredTest_TrpPep_AllFoldRun, SampleInfo)

### Protein abundance model: LiPPep ~ Age + Cohort + Sex + AgeCohort + AgeSex + CohortSex + TPC + Batch
PredTest_TrpProt_AllFoldRun <- GetTestPredOnFoldRuns(SR_TrpProt_AllFoldRun, 
                                                     info = SampleInfo, 
                                                     nLambda = DFnLambda$TrpProt)
SumPred_TrpProt <- GetMeanPredOfFoldRuns(PredTest_TrpProt_AllFoldRun, SampleInfo)

### Strucutral, PK-independent and abundance variation model: LiPPep ~ Age + Cohort + Sex + AgeCohort + AgeSex + CohortSex + TPC + Batch
PredTest_LiPPep_noPepProt_AllFoldRun <- GetTestPredOnFoldRuns(SR_LiPPep_noPepProt_AllFoldRun, 
                                                              info = SampleInfo, 
                                                              nLambda = DFnLambda$LiPPep_noPep)
SumPred_LiPPep_noPepProt <- GetMeanPredOfFoldRuns(PredTest_LiPPep_noPepProt_AllFoldRun, SampleInfo)

### Strucutral and abundance variation model: LiPPep ~ Age + Cohort + Sex + AgeCohort + AgeSex + CohortSex + TrpPep + TPC + Batch
PredTest_LiPPep_noProt_AllFoldRun <- GetTestPredOnFoldRuns(SR_LiPPep_noProt_AllFoldRun, 
                                                           info = SampleInfo, 
                                                           nLambda = DFnLambda$LiPPep_noProt)
SumPred_LiPPep_noProt <- GetMeanPredOfFoldRuns(PredTest_LiPPep_noProt_AllFoldRun, SampleInfo)

### Strucutral and PK-independent variation model: LiPPep ~ Age + Cohort + Sex + AgeCohort + AgeSex + CohortSex + TrpProt + TPC + Batch
PredTest_LiPPep_noPep_AllFoldRun <- GetTestPredOnFoldRuns(SR_LiPPep_noPep_AllFoldRun, 
                                                          info = SampleInfo, 
                                                          nLambda = DFnLambda$LiPPep_noPepProt)
SumPred_LiPPep_noPep <- GetMeanPredOfFoldRuns(PredTest_LiPPep_noPep_AllFoldRun, SampleInfo)
