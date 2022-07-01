library(devtools)
library(openxlsx)

source_url("https://github.com/LuiseNagel/Global-analyses-of-the-human-structural-proteome-to-identify-a-new-type-of-disease-biomarker/blob/main/R/Functions_StructuralVariabilityCSF.R?raw=TRUE")


# Loading data
## Peptide and protein quantities
LiPPep <- read.xlsx("Supplementary_Table_12_Data_for_GitHub_scripts.xlsx", "LiPPep_Preprocessed", startRow = 3, rowNames = T)
TrpPep <- read.xlsx("Supplementary_Table_12_Data_for_GitHub_scripts.xlsx", "TrpPep_Preprocessed", startRow = 3, rowNames = T)

TrpProt2TrpPep <- read.xlsx("Supplementary_Table_12_Data_for_GitHub_scripts.xlsx", "TrpProt_map2TrpPep_Preprocessed", startRow = 3, rowNames = T)

## Further peptide and sample information
PepProtInfo <- read.xlsx("Supplementary_Table_12_Data_for_GitHub_scripts.xlsx", "PepProtInfo", startRow = 3, rowNames = T)
SampleInfo <- read.xlsx("Supplementary_Table_12_Data_for_GitHub_scripts.xlsx", "SampleInfo", startRow = 3, rowNames = T)
SampleInfo <- SampleInfo[SampleInfo$Cohort == "H0",] # use only healthy donors

# Removing outliers using Cook's Distance
## Building linear models and estimating Cook's Distance from them
Peps2Use <- Reduce(intersect, list(row.names(LiPPep), 
                                   row.names(TrpPep),
                                   row.names(TrpProt2TrpPep),
                                   row.names(PepProtInfo)[PepProtInfo$IsTryptic == "full"]))

### LiP peptides
CD_LiPPep <- lapply(as.list(Peps2Use), function(i){
    lm(as.numeric(LiPPep[i, row.names(SampleInfo)]) ~ as.numeric(TrpProt2TrpPep[i, row.names(SampleInfo)]))
})
names(CD_LiPPep) <- Peps2Use
CD_LiPPep <- EstimateCooksDist(CD_LiPPep)
colnames(CD_LiPPep) <- row.names(SampleInfo)

### Trp peptides
CD_TrpPep <- lapply(as.list(Peps2Use), function(i){
    lm(as.numeric(TrpPep[i, row.names(SampleInfo)]) ~ as.numeric(TrpProt2TrpPep[i, row.names(SampleInfo)]))
})
names(CD_TrpPep) <- Peps2Use
CD_TrpPep <- EstimateCooksDist(CD_TrpPep)
colnames(CD_TrpPep) <- row.names(SampleInfo)

## Removing outliers according to Cook's Distance
LiPPep <- FilterCooksDist(LiPPep[Peps2Use, row.names(SampleInfo)], CD_LiPPep)
TrpPep <- FilterCooksDist(TrpPep[Peps2Use, row.names(SampleInfo)], CD_TrpPep)


# Estimating variability scores and corresponding values
## Estimating variability scores
VarPeps <- data.frame(SD_LipPep = apply(LiPPep, 1, sd, na.rm = T),
                      SD_TrpPep = apply(TrpPep, 1, sd, na.rm = T),
                      row.names = Peps2Use)
VarPeps$VarS <- VarPeps$SD_LipPep - VarPeps$SD_TrpPep

## Estimating corresponding p-value using   Levene's test                     
VarPeps$Pval <- sapply(Peps2Use, function(i){
    MyDf <- data.frame(Int = c(as.numeric(TrpPep[i,]), as.numeric(LiPPep[i,])),
                       Type = c(rep("Trp", ncol(TrpPep)), rep("LiP", ncol(LiPPep))))
    MyDf <- MyDf[apply(MyDf, 1, function(x) {sum(is.na(x))==0}),]
    LevTest <- levene.test(MyDf$Int, MyDf$Type)
    return(LevTest$p.value)
})

VarPeps$Classification <- sapply(Peps2Use, function(i){
    if(VarPeps[i, "VarS"] > 1 & -log10(VarPeps[i, "Pval"]) > 5){
        out <- "HighVar"
    }
    else if(VarPeps[i, "VarS"] > 0.5 & -log10(VarPeps[i, "Pval"]) > 5){
        out <- "MiddleVar"
    }
    else{
        out <- "NonVar"
    }
    return(out)
})

table(VarPeps$Classification)
## Estimating multimodel scores
VarPeps$MM_LiPPep <- GetMultiModal(LiPPep)
VarPeps$MM_TrpPep <- GetMultiModal(TrpPep)
VarPeps$IsMultiModal <- sapply(c(1:nrow(VarPeps)), function(i){
    ifelse(VarPeps[i, "MM_LiPPep"]<0.05&VarPeps[i, "MM_TrpPep"]>0.05, TRUE, FALSE)
})

