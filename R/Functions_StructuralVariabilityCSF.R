library(lawstat)
library(multimode)
library(plyr)

# Function for estimating Cook's Distance from linear models
EstimateCooksDist <- function(LM_List){
    CooksD <- do.call(rbind.fill, lapply(LM_List, function(x){
        as.data.frame(t(cooks.distance(x)))
    }))
    row.names(CooksD) <- names(LM_List)
    return(CooksD)
}

# Function for filtering peptide levels according to Cook's Distance
FilterCooksDist <- function(PepData, CooksD){
    CooksD <- t(apply(CooksD, 1, function(x){
        cutoff <- 4/sum(!is.na(x))
        x<cutoff
    }))
    PepData[CooksD == FALSE] <- NA
    PepData[is.na(CooksD)] <- NA
    return(PepData)
}

# Function for estimating if peptide is multimodal
GetMultiModal <- function(PepData){
    sapply(row.names(PepData), function(i){
        set.seed(42)
        x <- as.numeric(PepData[i,])
        modetest(na.omit(x))$p.value
    })
}

