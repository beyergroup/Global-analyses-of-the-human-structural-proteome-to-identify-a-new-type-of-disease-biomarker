# Function for running linear models with/without including different variables in the model
RunLM <- function(data, dataPep, dataProt, info, Age = NULL, Cohort = NULL, Sex = NULL, AgeCohort = NULL, AgeSex = NULL, CohortSex = NULL, Pep = NULL, Prot = NULL, TPC = NULL, Batch = NULL){
    ## Setting covariables for linear model
    if(!is.null(Age)){Age <- info$Age}
    if(!is.null(Cohort)){Cohort <- info$Cohort}
    if(!is.null(Sex)){Sex <- info$Sex}
    if(!is.null(AgeCohort)){AgeCohort <- info$AgeCohort}
    if(!is.null(AgeSex)){AgeSex <- info$AgeSex}
    if(!is.null(CohortSex)){CohortSex <- info$CohortSex}
    if(!is.null(TPC)){TPC <- info$TPC}
    if(!is.null(Batch)){Batch <- as.factor(info$Batch)}
    
    ## Running linear models
    LM <- lapply(as.list(1:nrow(data)), function(i){
        if(!is.null(Pep)){Peptide <-  as.numeric(dataPep[i,])}
        if(!is.null(Prot)){Protein <- as.numeric(dataProt[i,])}
        lm(formula = as.formula(paste0("as.numeric(data[i,]) ~ ",
                                       switch(is.null(Age)+1, "Age", NULL),
                                       switch(is.null(Cohort)+1, " + Cohort", NULL),
                                       switch(is.null(Sex)+1, " + Sex", NULL), 
                                       switch((is.null(AgeCohort))+1, "+ AgeCohort", NULL),
                                       switch((is.null(AgeSex))+1, "+ AgeSex", NULL),
                                       switch((is.null(CohortSex))+1, "+ CohortSex", NULL),
                                       switch(is.null(Pep)+1, " + Peptide", NULL),
                                       switch(is.null(Prot)+1, " + Protein", NULL),
                                       switch(is.null(TPC)+1, " + TPC", NULL),
                                       switch(is.null(Batch)+1, " + Batch", NULL))), 
           na.action = na.exclude)
    })
    names(LM) <- row.names(data)
    return(LM)
}


# Function for filtering list of linear models based on number of samples modelled in every peptide
FilterLM <- function(LM, minResid){
    NResid <- unlist(lapply(LM, function(x) length(x$residuals)))
    LM <- LM[NResid >= minResid]
    return(LM)
}

# Function for extracting information from linear model
ExtractInfoLM <- function(LM, variableName){
    
    ## Fetching coefficients for defined variable from linear model
    coeffVar <- unlist(lapply(LM, function(x){
        as.data.frame(t(summary.lm(x)$coefficients[variableName, 1]))
        }))
    
    ## Fetching p-values for defined variable from linear model
    pvalVar <- unlist(lapply(LM, function(x){
        as.data.frame(t(summary.lm(x)$coefficients[variableName, 4]))
    }))
    out <- data.frame(Coefficient = coeffVar, 
                      Pvalue = pvalVar, 
                      row.names = names(LM))
    return(out)
}

# Function for protein-wise FDR correction
DoProteinwiseFDR <- function(varPval, method){
    names <- names(varPval)
    varPval <- data.frame(varPval,
                       Prot = gsub("[-,;,_].*", "", names(varPval))) # neglecting isoform information for this step
    varPval <- split(varPval, varPval$Prot) # splitting data by protein
    #cat(length(varPval), "different proteins contained in the data.\n")
    
    ## Performing FDR on all peptides of each protein 
    protwFDR <- lapply(varPval, function(x){
        x[,1] <- p.adjust(x[,1], method = method)
        return(x)
    })
    
    ## Adding protein-wise FDRs in one data.frame and adjusting row.names
    protwFDR <- do.call(rbind, lapply(protwFDR, function(x) x[, 1, drop=F]))
    row.names(protwFDR) <- gsub("[.]", ";", row.names(protwFDR))
    outFDR <- protwFDR[,1]
    names(outFDR) <- row.names(protwFDR)
    names(outFDR) <- unname(sapply(names(outFDR), function(x){
        ifelse(grepl("_", x), x, names[grepl(x, names)])
    }))
    names(outFDR) <- gsub(".*;", "", names(outFDR))
    return(outFDR)
}


# Function for running RunLMandGetScaledResidualsOnFolds on all FoldRuns
RunLMandGetScaledResidualsAllFoldRuns <- function(foldList, data, dataPep, dataProt, info, Age = NULL, Cohort = NULL, Sex = NULL, AgeCohort = NULL, AgeSex = NULL, CohortSex = NULL, Pep = NULL, Prot = NULL, TPC = NULL, Batch = NULL, pwFDR = T){
    ListSR_FoldRuns <- list(FoldRun1 = RunLMandGetScaledResidualsOnFolds(folds = foldList[[1]], data, dataPep, dataProt, info, Age, Cohort, Sex, AgeCohort, AgeSex, CohortSex, Pep, Prot, TPC, Batch, FRun = 1),
                            FoldRun2 = RunLMandGetScaledResidualsOnFolds(folds = foldList[[2]], data, dataPep, dataProt, info, Age, Cohort, Sex, AgeCohort, AgeSex, CohortSex, Pep, Prot, TPC, Batch, FRun = 2),
                            FoldRun3 = RunLMandGetScaledResidualsOnFolds(folds = foldList[[3]], data, dataPep, dataProt, info, Age, Cohort, Sex, AgeCohort, AgeSex, CohortSex, Pep, Prot, TPC, Batch, FRun = 3),
                            FoldRun4 = RunLMandGetScaledResidualsOnFolds(folds = foldList[[4]], data, dataPep, dataProt, info, Age, Cohort, Sex, AgeCohort, AgeSex, CohortSex, Pep, Prot, TPC, Batch, FRun = 4),
                            FoldRun5 = RunLMandGetScaledResidualsOnFolds(folds = foldList[[5]], data, dataPep, dataProt, info, Age, Cohort, Sex, AgeCohort, AgeSex, CohortSex, Pep, Prot, TPC, Batch, FRun = 5))
    return(ListSR_FoldRuns)
}


# Function for running linear models on training folds and estimating scaled residuals for test and training dataset from them
RunLMandGetScaledResidualsOnFolds <- function(folds, data, dataPep, dataProt, info, Age = NULL, Cohort = NULL, Sex = NULL, AgeCohort = NULL, AgeSex = NULL, CohortSex = NULL, Pep = NULL, Prot = NULL, TPC = NULL, Batch = NULL, FRun = NULL, pwFDR = T){
    
    if(!is.null(FRun)){
        message("Running Fold Run ", FRun, ".")
    }
    
    ## Rnunning linear models and estimated scaled residuals for every fold per RunFold
    SR_AllFolds <- lapply(as.list(1:5), function(i){
        
     
        ## Splitting data into training and test subsets
        infoTrain <- info[folds!=i,]
        dataTrain <- data[, folds!=i]
        dataPepTrain <- dataPep[, folds!=i]
        dataProtTrain <- dataProt[, folds!=i]
        
        infoTest <- info[folds==i,]
        dataTest <- data[, folds==i]
        dataPepTest <- dataPep[, folds==i]
        dataProtTest <- dataProt[, folds==i]
        
        ## Running linear model on training data

        LMTrain <- RunLM(data = dataTrain,
                          dataPep = dataPepTrain,
                          dataProt = dataProtTrain,
                          info = infoTrain, 
                          Age, Cohort, Sex, AgeCohort, AgeSex, CohortSex, Pep, Prot, TPC, Batch)
        
        ### Extracting coefficients and p-values for the cohort coefficients
        Cohort_LMTrain <- ExtractInfoLM(LMTrain, "CohortT0")
        
        if(pwFDR){
            Cohort_LMTrain$proteinwiseFDR <- DoProteinwiseFDR(setNames(Cohort_LMTrain$Pvalue, row.names(Cohort_LMTrain)), "fdr") ### correcting p-values protein-wise
            SigHits <- row.names(Cohort_LMTrain)[Cohort_LMTrain$proteinwiseFDR < 0.05]
        }
        else{
            SigHits <- row.names(Cohort_LMTrain)[Cohort_LMTrain$Pval < 0.05]
        }
        
        ## Filtering linear models to only containing those, were the peptide shows a cohort effect
        LMTrain <- LMTrain[SigHits]

        ## Predicting quantities of training and test data assuming all samples are H0 or all are T0
        PredH0_TrainData <- PredictLM(LM = LMTrain, dataPep = dataPep, dataProt = dataProt, info = infoTrain, 
                                      Age = Age, Cohort = "H0", Sex = Sex, 
                                      AgeCohort = AgeCohort, AgeSex = AgeSex, CohortSex = CohortSex, 
                                      Pep = Pep, Prot = Prot, TPC = TPC, Batch = Batch) 
        PredT0_TrainData <- PredictLM(LM = LMTrain, dataPep = dataPep, dataProt = dataProt, info = infoTrain, 
                                      Age = Age, Cohort = "T0", Sex = Sex, 
                                      AgeCohort = AgeCohort, AgeSex = AgeSex, CohortSex = CohortSex, 
                                      Pep = Pep, Prot = Prot, TPC = TPC, Batch = Batch)
        
        PredH0_TestData <- PredictLM(LM = LMTrain, dataPep = dataPep, dataProt = dataProt,info = infoTest, 
                                     Age = Age, Cohort = "H0", Sex = Sex, 
                                     AgeCohort = AgeCohort, AgeSex = AgeSex, CohortSex = CohortSex, 
                                     Pep = Pep, Prot = Prot, TPC = TPC, Batch = Batch)
        PredT0_TestData <- PredictLM(LM = LMTrain, dataPep = dataPep, dataProt = dataProt, info = infoTest, 
                                     Age = Age, Cohort = "T0", Sex = Sex, 
                                     AgeCohort = AgeCohort, AgeSex = AgeSex, CohortSex = CohortSex, 
                                     Pep = Pep, Prot = Prot, TPC = TPC, Batch = Batch)
        
        # Estimating scalred residuals
        SR_TrainData <- EstimateScaledResid(data[row.names(PredH0_TrainData), colnames(PredH0_TrainData)], PredH0_TrainData, PredT0_TrainData)
        SR_TestData <- EstimateScaledResid(data[row.names(PredH0_TestData), colnames(PredH0_TestData)], PredH0_TestData, PredT0_TestData)
        
        message("Done with fold ", i, '!')
        return(list("SR_TrainD" = SR_TrainData,
                    "SR_TestD" = SR_TestData))
    })
    names(SR_AllFolds) <- c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5")
    return(SR_AllFolds)
}


# Function for predicting peptide quanitities from the linear models allowing for setting cohort varibale to fixed value (H0 or)
PredictLM <- function(LM, dataPep, dataProt, info, Age, Cohort, Sex, AgeCohort, AgeSex, CohortSex, Pep, Prot, TPC, Batch){
    Peps <- names(LM)
    Samples <- row.names(info)
    
    if(is.null(Pep)){
        Pep <- F
    }
    if(is.null(Prot)){
        Prot <- F
    }
    
    ## if cohort is set to H0 or T0 adjust this in the interaction terms
    if(Cohort != T){
        info$AgeCohort <- sapply(info$Age, function(x) ifelse(Cohort == "H0", 0, x))
        info$CohortSex <- sapply(row.names(info), function(x){info$CohortSex[info$Cohort == Cohort & as.vector(info$Sex) == as.vector(info[x, "Sex"])][1]})
    }
    
    ## set sex and cohort information into numeric data for estimating predicted values
    info$Sex <- sapply(info$Sex, function(x) ifelse(x==as.numeric(row.names(contrasts(info$Sex))[1]), contrasts(info$Sex)[1], contrasts(info$Sex)[2]))
    info$Cohort <- sapply(info$Cohort, function(x){
        if(Cohort != T){x <- Cohort}
        ifelse(x==row.names(contrasts(info$Cohort))[1], contrasts(info$Cohort)[1], contrasts(info$Cohort)[2])
    })
    
    ## Limiting data to correct peptides and samples
    info <- info[Samples, ]
    dataPep <- dataPep[Peps, Samples]
    dataProt <- dataProt[Peps, Samples]
    
    ## Predicting peptide quantities
    Pred <- as.data.frame(do.call(rbind, lapply(as.list(Peps), function(i){
      
        unlist(lapply(as.list(Samples), function(j){
            
            LM_Pep <- LM[[i]]
            if(info[j, "Batch"]%in% LM_Pep$xlevels$Batch){
                Pred <- unname(LM_Pep$coefficients["(Intercept)"] +
                                   ifelse(Age, LM_Pep$coefficients["Age"] * info[j, "Age"], 0) +
                                   unname(LM_Pep$coefficients["CohortT0"] * info[j, "Cohort"]) +
                                   ifelse(Sex, LM_Pep$coefficients["Sex1"] * info[j, "Sex"], 0) +
                                   ifelse(AgeCohort, LM_Pep$coefficients["AgeCohort"] * info[j, "AgeCohort"], 0) +
                                   ifelse(AgeSex, LM_Pep$coefficients["AgeSex"] * info[j, "AgeSex"], 0) +
                                   ifelse(CohortSex, LM_Pep$coefficients["CohortSex"] * info[j, "CohortSex"], 0) +
                                   ifelse(Pep, LM_Pep$coefficients["Peptide"] * dataPep[i,j], 0) +
                                   ifelse(Prot, LM_Pep$coefficients["Protein"] * dataProt[i,j], 0) +
                                   ifelse(TPC, LM_Pep$coefficients["TPC"] * info[j, "TPC"], 0) +
                                   ifelse(Batch, ifelse(paste0("Batch", info[j, "Batch"]) %in% names(LM_Pep$coefficients), LM_Pep$coefficients[paste0("Batch", info[j, "Batch"])], 0)),0) 
            }
            else{Pred <- NA} # taking care of option of a batch not being modelled due to fold seperation, set predictied value to NA then
            return(Pred)
        }))
    })))
    
    colnames(Pred) <- Samples
    row.names(Pred) <- Peps
    return(Pred)
}

# Function for estimating scaled residuals, limiting the denominator to not being below 1 SD
EstimateScaledResid <- function(LiP, PredH0, PredT0, LimitPredH0 = PredH0, LimitPredT0 = PredT0){
    
    ## Estimationg nominator and denominator
    Nom <- LiP[row.names(PredH0), colnames(PredH0)] - PredH0
    Denom <- PredT0 - PredH0
    
    ## Defining limit denominator can have
    DenomLim <- LimitPredT0 - LimitPredH0
    Limit <- apply(DenomLim, 1, sd, na.rm = T)
    
    ## Changing denominator if needed
    sapply(names(Limit), function(i){
        l <- Limit[i]
        Denom[row.names(Denom)==i & Denom>-l & Denom<0] <<- -l
        Denom[row.names(Denom)==i & Denom<l & Denom>0] <<- l
        return(NULL)
    })
    
    ## Estimating scaled residual
    SR <- Nom/Denom
    return(SR)
}


