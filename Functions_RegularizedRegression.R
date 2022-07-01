library(glmnet)

# Function for getting mean prediction of all models per number of predictors for each sample
GetMeanPredOfFoldRuns <- function(predTest, info){
    
    # Get dataframe with all predictions for 1 predictor
    df1Pred <- do.call(cbind, lapply(predTest, function(y){
        predFolds <- do.call(rbind, lapply(y, function(z){
            z$`1Pred`
        }))
        predFolds[row.names(info),]
    }))
    
    # Get dataframe with all predictions for 5 predicto2
    df5Pred <- do.call(cbind, lapply(predTest, function(y){
        predFolds <- do.call(rbind, lapply(y, function(z){
            z$`5Pred`
        }))
        predFolds[row.names(info),]
    }))
    
    # Get dataframe with all predictions for 10 predictors
    df10Pred <- do.call(cbind, lapply(predTest, function(y){
        predFolds <- do.call(rbind, lapply(y, function(z){
            z$`10Pred`
        }))
        predFolds[row.names(info),]
    }))
    
    AllPredDF <- data.frame(Pred1 = rowMeans(df1Pred, na.rm = T),
                            Pred5 = rowMeans(df5Pred, na.rm = T),
                            Pred10 = rowMeans(df10Pred, na.rm = T),
                            row.names = row.names(info))
    return(AllPredDF)
}


# Function for callling function GetTestPredOnFolds an all fold runs
GetTestPredOnFoldRuns <- function(SR_AllFoldRun, info, nPred = c(1, 5, 10), nLambda){
    
    FoldRuns <- list("FoldRun1", "FoldRun2", "FoldRun3", "FoldRun4", "FoldRun5")
    predTest_FoldRuns <- lapply(FoldRuns, function(i){
        message("Running Fold Run ", gsub("FoldRun", "", i), ".")
        GetTestPredOnFolds(SR_AllFoldRun[[i]], info, nPred, nLambda[[i]])
    })
    names(predTest_FoldRuns) <- as.character(FoldRuns)

    return(predTest_FoldRuns)
}

# Function for building LASSO models on training data and running it on test data for folds
GetTestPredOnFolds <- function(SR_Fold, info, nPred, nLambda){
    predTest <- lapply(as.list(1:5), function(i){
        #message("Running Fold ", i, ".")
        x <- SR_Fold[[i]]
        xnLambda <- nLambda[i]
        TrainAndEvalLASSO(x$SR_TrainD, x$SR_TestD, info, nPred, xnLambda)
    })
    names(predTest) <- c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5")
    return(predTest)
}

# Function for training LASSO models and running them on test fold vor individual folds
TrainAndEvalLASSO <- function(SRTrain, SRTest, info, nPred, nLambda){
    
    ## seperating info file into training and test samples
    infoTrain <- info[colnames(SRTrain),]
    infoTest <- info[colnames(SRTest),]
    
    # Sorting data
    ## Imputing NAs in scaled residuals
    SRTrain <- ImpNAsInSR(SRTrain)
    SRTest <- ImpNAsInSR(SRTest)
    
    LASSOModels <- TrainLASSO(dataTrain = SRTrain, 
                              infoTrain = infoTrain,
                              nPred = nPred,
                              nLambda = nLambda)
    
    predSRTrain <- lapply(LASSOModels, PredictOnLASSO, dataTest = SRTest)
    return(predSRTrain)
}

# Function for imputing mean scaled residual over all samples for NAs 
ImpNAsInSR <- function(SR){
    SR <- as.data.frame(t(apply(SR, 1, function(x){
        sapply(x, function(y){
            ifelse(is.na(y), mean(x, na.rm = T), y)
        })
    })))
    return(SR)
}

# Function for training LASSO models on training data
TrainLASSO <- function(dataTrain, infoTrain, nPred, nLambda){
    mLASSO <- lapply(as.list(nPred), function(x){
        lambda <- GetLamdaViaCrossVali(dataTrain, infoTrain, x, nLambda)
        set.seed(42)
        LASSO <- glmnet(x = t(dataTrain),
                        y = infoTrain$Cohort,
                        family = "binomial",
                        alpha = 1,
                        lambda =  lambda)
        return(LASSO)
    })
    names(mLASSO) <- c("1Pred", "5Pred", "10Pred")
    return(mLASSO)
}


# Function for getting lambda via cross validation
GetLamdaViaCrossVali <- function(dataTrain, infoTrain, dfmax, nLambda){
    set.seed(42)
    
    # if less peptides available than needed as predictors, set max number of peptides to number of available peptides
    if(dfmax > nrow(dataTrain)){
        dfmax <- nrow(dataTrain)
    }
    
    dfmax2 <- dfmax+1
    if(dfmax2 > nrow(dataTrain)){
        dfmax2 <- nrow(dataTrain)
    }
    
    ## run cross validation scheme
    cvLASSO <- cv.glmnet(x = t(dataTrain),
                         y = infoTrain$Cohort,
                         family = "binomial",
                         alpha = 1,
                         dfmax = dfmax2,
                         nlambda = nLambda)
    
    ## retrieve lamba for the training model for a specific number of predictors from the cross validation models
    lambdas <- cvLASSO$glmnet.fit$lambda[cvLASSO$glmnet.fit$df == dfmax]
    
    if(length(lambdas) == 0){
        dfmax <- dfmax-1
        lambdas <- cvLASSO$glmnet.fit$lambda[cvLASSO$glmnet.fit$df == dfmax]
    }
    
    MCV <- cvLASSO$cvm[cvLASSO$glmnet.fit$df == dfmax]
    lambda <- lambdas[which.min(MCV)]
    return(lambda)
}

# Function for predicting cohort based on LASSO models
PredictOnLASSO <- function(mLASSO, dataTest){
    predValues <- predict(mLASSO, 
                          s = mLASSO$lambda, 
                          newx = t(dataTest))
    colnames(predValues) <- "Pred"
    #print(sum(as.matrix(mLASSO$beta)!=0))
    return(predValues)
}

# Function for fetching lambas for running cross-validation models
GetLambdas <- function(){list(LiPPep = data.frame(FoldRun1 = c(5000, 5000, 5000, 5000, 5000),
                                      FoldRun2 = c(5000, 5000, 5000, 5000, 3000),
                                      FoldRun3 = c(5000, 5000, 5000, 5000, 5000),
                                      FoldRun4 = c(5000, 5000, 5000, 5000, 1000),
                                      FoldRun5 = c(5000, 5000, 5000, 5000, 5000),
                                      row.names = c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5")),
                  TrpPep = data.frame(FoldRun1 = c(5000, 5000, 5000, 5000, 5000),
                                      FoldRun2 = c(5000, 5000, 5000, 5000, 5000),
                                      FoldRun3 = c(5000, 5000, 5000, 5000, 5000),
                                      FoldRun4 = c(5000, 5000, 1000, 1000, 5000),
                                      FoldRun5 = c(5000, 5000, 5000, 5000, 5000),
                                      row.names = c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5")),
                  TrpProt = data.frame(FoldRun1 = c(5000, 5000, 5000, 5000, 5000),
                                      FoldRun2 = c(5000, 20000, 20000, 5000, 20000),
                                      FoldRun3 = c(5000, 5000, 5000, 5000, 5000),
                                      FoldRun4 = c(5000, 5000, 5000, 5000, 5000),
                                      FoldRun5 = c(5000, 5000, 5000, 5000, 5000),
                                      row.names = c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5")),
                  LiPPep_noPepProt = data.frame(FoldRun1 = c(5000, 5000, 5000, 5000, 5000),
                                                FoldRun2 = c(30000, 5000, 5000, 5000, 5000),
                                                FoldRun3 = c(5000, 3000, 5000, 5000, 5000),
                                                FoldRun4 = c(5000, 5000, 5000, 5000, 5000),
                                                FoldRun5 = c(5000, 5000, 5000, 5000, 5000),
                                                row.names = c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5")),
                  LiPPep_noProt = data.frame(FoldRun1 = c(5000, 5000, 5000, 5000, 5000),
                                             FoldRun2 = c(5000, 5000, 5000, 5000, 5000),
                                             FoldRun3 = c(5000, 5000, 5000, 1000, 10000),
                                             FoldRun4 = c(20000, 5000, 5000, 5000, 5000),
                                             FoldRun5 = c(5000, 5000, 5000, 5000, 5000),
                                             row.names = c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5")),
                  LiPPep_noPep = data.frame(FoldRun1 = c(5000, 5000, 5000, 5000, 5000),
                                             FoldRun2 = c(5000, 5000, 5000, 5000, 5000),
                                             FoldRun3 = c(5000, 5000, 5000, 5000, 5000),
                                             FoldRun4 = c(5000, 5000, 5000, 5000, 5000),
                                             FoldRun5 = c(5000, 5000, 5000, 5000, 5000),
                                             row.names = c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5")))
}


