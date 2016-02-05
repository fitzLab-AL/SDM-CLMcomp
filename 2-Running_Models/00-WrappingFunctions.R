################################################################################
##
## Wrapping functions to save code afterwards
##
################################################################################

# GLM - SDM
glmFit <- function(commM, xData, train, iterNum){
    spNames <- colnames(commM)
    glmList <- foreach(i=1:iterNum, .verbose=F, .errorhandling="pass") %dopar% {
        iTrain <- train[,i]
        sub.commM <- subset(commM, iTrain)
        sub.xData <- subset(xData, iTrain)
        forms <- lapply(spNames, FUN=function(x){formula(paste(x, "~ etr_year_ave + aet_year_ave + wdei_year_ave + tmax_high_quart + prcp_low_quart + prcp_high_quart + I(etr_year_ave^2) + I(aet_year_ave^2) + I(wdei_year_ave^2) + I(tmax_high_quart^2) + I(prcp_low_quart^2) + I(prcp_high_quart^2)"))})
        glmfit <- lapply(forms, FUN=function(m, d){glm(m, family=binomial, data=d)}, cbind(sub.commM, sub.xData))
        glmfit
    }
    return(glmList)
}

# CQO - CLM
cqoFit <- function(commM, xData, train, iterNum){
    require(VGAM)
    cqoList <- foreach(i=1:iterNum, .packages=c("VGAM"), .verbose=F, .errorhandling="pass") %dopar% {
        iTrain <- train[,i]
        sub.commM <- subset(commM, iTrain)
        sub.xData <- subset(xData, iTrain)
        form <- formula(sub.commM ~ etr_year_ave + aet_year_ave + wdei_year_ave + tmax_high_quart + prcp_low_quart + prcp_high_quart)
        #    cqofit <- cqo(form, family=binomialff(mv=TRUE), data=sub.xData, trace=TRUE, Rank=2)    
        cqofit <- cqo(form, family=binomialff(multiple.responses=TRUE), data=sub.xData, trace=TRUE)
        cqofit
    }
    return(cqoList)
}

# GAM - SDM
gamFit <- function(commM, xData, train, iterNum){
    require(gam)
    spNames <- colnames(commM)
    gamList <- foreach(i=1:iterNum, .packages=c("gam"), .verbose=F, .errorhandling="pass") %dopar% {
        iTrain <- train[,i]
        sub.commM <- subset(commM, iTrain)
        sub.xData <- subset(xData, iTrain)
        forms=list()
        forms[[1]]=formula(paste(spNames[1],"~s(etr_year_ave,df=1) + s(aet_year_ave, df=1) + s(wdei_year_ave, df=1) + s(tmax_high_quart, df=1) + s(prcp_low_quart, df=1) + s(prcp_high_quart, df=1)"))    
        forms[[2]]=formula(paste(spNames[2],"~s(etr_year_ave,df=1) + s(aet_year_ave, df=1) + s(wdei_year_ave, df=1) + s(tmax_high_quart, df=1) + s(prcp_low_quart, df=1) + s(prcp_high_quart, df=1)"))
        forms[[3]]=formula(paste(spNames[3],"~s(etr_year_ave,df=1) + s(aet_year_ave, df=1) + s(wdei_year_ave, df=1) + s(tmax_high_quart, df=1) + s(prcp_low_quart, df=1) + s(prcp_high_quart, df=1)"))
        forms[[4]]=formula(paste(spNames[4],"~s(etr_year_ave,df=1) + s(aet_year_ave, df=1) + s(wdei_year_ave, df=1) + s(tmax_high_quart, df=1) + s(prcp_low_quart, df=1) + s(prcp_high_quart, df=1)"))
        forms[[5]]=formula(paste(spNames[5],"~s(etr_year_ave,df=1) + s(aet_year_ave, df=1) + s(wdei_year_ave, df=1) + s(tmax_high_quart, df=1) + s(prcp_low_quart, df=1) + s(prcp_high_quart, df=1)"))
        forms[[6]]=formula(paste(spNames[6],"~s(etr_year_ave,df=1) + s(aet_year_ave, df=1) + s(wdei_year_ave, df=1) + s(tmax_high_quart, df=1) + s(prcp_low_quart, df=1) + s(prcp_high_quart, df=1)"))
        forms[[7]]=formula(paste(spNames[7],"~s(etr_year_ave,df=1) + s(aet_year_ave, df=1) + s(wdei_year_ave, df=1) + s(tmax_high_quart, df=1) + s(prcp_low_quart, df=1) + s(prcp_high_quart, df=1)"))
        forms[[8]]=formula(paste(spNames[8],"~s(etr_year_ave,df=1) + s(aet_year_ave, df=1) + s(wdei_year_ave, df=1) + s(tmax_high_quart, df=1) + s(prcp_low_quart, df=1) + s(prcp_high_quart, df=1)"))
        forms[[9]]=formula(paste(spNames[9],"~s(etr_year_ave,df=1) + s(aet_year_ave, df=1) + s(wdei_year_ave, df=1) + s(tmax_high_quart, df=1) + s(prcp_low_quart, df=1) + s(prcp_high_quart, df=1)"))
        forms[[10]]=formula(paste(spNames[10],"~s(etr_year_ave,df=1) + s(aet_year_ave, df=1) + s(wdei_year_ave, df=1) + s(tmax_high_quart, df=1) + s(prcp_low_quart, df=1) + s(prcp_high_quart, df=1)"))
        forms[[11]]=formula(paste(spNames[11],"~s(etr_year_ave,df=1) + s(aet_year_ave, df=1) + s(wdei_year_ave, df=1) + s(tmax_high_quart, df=1) + s(prcp_low_quart, df=1) + s(prcp_high_quart, df=1)"))
        forms[[12]]=formula(paste(spNames[12],"~s(etr_year_ave,df=1) + s(aet_year_ave, df=1) + s(wdei_year_ave, df=1) + s(tmax_high_quart, df=1) + s(prcp_low_quart, df=1) + s(prcp_high_quart, df=1)"))
        forms[[13]]=formula(paste(spNames[13],"~s(etr_year_ave,df=1) + s(aet_year_ave, df=1) + s(wdei_year_ave, df=1) + s(tmax_high_quart, df=1) + s(prcp_low_quart, df=1) + s(prcp_high_quart, df=1)"))
        forms[[14]]=formula(paste(spNames[14],"~s(etr_year_ave,df=1) + s(aet_year_ave, df=1) + s(wdei_year_ave, df=1) + s(tmax_high_quart, df=1) + s(prcp_low_quart, df=1) + s(prcp_high_quart, df=1)"))
        forms[[15]]=formula(paste(spNames[15],"~s(etr_year_ave,df=1) + s(aet_year_ave, df=1) + s(wdei_year_ave, df=1) + s(tmax_high_quart, df=1) + s(prcp_low_quart, df=1) + s(prcp_high_quart, df=1)"))
        forms[[16]]=formula(paste(spNames[16],"~s(etr_year_ave,df=1) + s(aet_year_ave, df=1) + s(wdei_year_ave, df=1) + s(tmax_high_quart, df=1) + s(prcp_low_quart, df=1) + s(prcp_high_quart, df=1)"))
        forms[[17]]=formula(paste(spNames[17],"~s(etr_year_ave,df=1) + s(aet_year_ave, df=1) + s(wdei_year_ave, df=1) + s(tmax_high_quart, df=1) + s(prcp_low_quart, df=1) + s(prcp_high_quart, df=1)"))
        forms[[18]]=formula(paste(spNames[18],"~s(etr_year_ave,df=1) + s(aet_year_ave, df=1) + s(wdei_year_ave, df=1) + s(tmax_high_quart, df=1) + s(prcp_low_quart, df=1) + s(prcp_high_quart, df=1)"))
        forms[[19]]=formula(paste(spNames[19],"~s(etr_year_ave,df=1) + s(aet_year_ave, df=1) + s(wdei_year_ave, df=1) + s(tmax_high_quart, df=1) + s(prcp_low_quart, df=1) + s(prcp_high_quart, df=1)"))
        gamfit <- lapply(forms, FUN=function(m, d){gam(m, family=binomial, data=d)}, cbind(sub.commM, sub.xData))
        gamfit
    }
    return(gamList)
}      

# CAO - CLM
caoFit <- function(commM, xData, train, iterNum){
    require(VGAM)
    caoList <- foreach(i=1:iterNum, .packages=c("VGAM"), .verbose=F, .errorhandling="pass") %dopar% {
        iTrain <- train[,i]
        sub.commM <- subset(commM, iTrain)
        sub.xData <- subset(xData, iTrain)
        form <- formula(sub.commM ~ etr_year_ave + aet_year_ave + wdei_year_ave + tmax_high_quart + prcp_low_quart + prcp_high_quart)
        caofit <- cao(form, family=binomialff(multiple.responses=TRUE), data=sub.xData, trace=TRUE, df1.nl=1)
        caofit
    }
    return(caoList)
}


# MARS - SDM
marsFit.sp <- function(commM, xData, train, iterNum, para){
    require(earth)
    taxonNames <- colnames(commM)
    marsList.sp <- foreach(i=1:iterNum, .packages=c("earth"), .verbose=F, .errorhandling="pass") %dopar% {
        iTrain <- train[,i]
        sub.commM <- subset(commM, iTrain)
        sub.xData <- subset(xData, iTrain)
        marsfit.sp <- lapply(1:19, FUN=function(m){
          if(para[m,7] == "NULL"){nPrune <- NULL}
          earth(x=sub.xData, y=sub.commM[,taxonNames[m]], nprune=nPrune, degree=para[m,8], trace=T, glm=list(family=binomial))
        })
        marsfit.sp
    }
    return(marsList.sp)
}


# MARS - CLM 
marsFit <- function(commM, xData, train, iterNum, para){
    require(earth)
    marsList <- foreach(i=1:iterNum, .packages=c("earth"), .verbose=F, .errorhandling="pass") %dopar% {
        iTrain <- train[,i]
        sub.commM <- subset(commM, iTrain)
        sub.xData <- subset(xData, iTrain)
        if(para[25,7] == "NULL"){nPrune <- NULL}
        marsfit <- earth(x=sub.xData, y=sub.commM, nprune=nPrune, degree=para[25,8], trace=T, glm=list(family=binomial))
        marsfit
    }
    return(marsList)
}


# CARTS - SDM
cartsFit.sp <- function(commM, xData, train, iterNum){
    require(rpart)
    cartsList.sp <- foreach(i=1:iterNum, .packages=c("rpart"), .verbose=F, .errorhandling="pass") %dopar% {
        iTrain <- train[,i]
        sub.commM <- subset(commM, iTrain)
        sub.xData <- subset(xData, iTrain)
        tNames <- colnames(commM)
        cartsfit.sp <- lapply(tNames,FUN=function(c){rpart(sub.commM[,c]~etr_year_ave + aet_year_ave + wdei_year_ave + tmax_high_quart + prcp_low_quart + prcp_high_quart,cbind(sub.commM, sub.xData))})
        cpValues <- lapply(cartsfit.sp, FUN=function(x){as.data.frame(printcp(x))})
        cpTarget <- lapply(cpValues, FUN=function(x){x[which.min(x[,"xerror"]), "CP"]})
        cartsfit.sp <- mapply(FUN=function(x, y){prune(x, cp=y)}, cartsfit.sp, cpTarget, SIMPLIFY=FALSE)
        cartsfit.sp
    }
    return(cartsList.sp)
}


# CARTS - CLM
cartsFit <- function(commM, xData, train, iterNum){
    require(mvpart)
    cartsList <- foreach(i=1:iterNum, .packages=c("mvpart"), .verbose=F, .errorhandling="pass") %dopar% {
        iTrain <- train[,i]
        sub.commM <- subset(commM, iTrain)
        sub.xData <- subset(xData, iTrain)
        form <- formula(sub.commM ~ etr_year_ave + aet_year_ave + wdei_year_ave + tmax_high_quart + prcp_low_quart + prcp_high_quart)
        cartsfit <- mvpart(form, cbind(sub.commM, sub.xData), plot.add=FALSE)
        cpValues <- as.data.frame(printcp(cartsfit))
        cpTarget <- cpValues[which.min(cpValues[,"xerror"]), "CP"]
        cartsfit <- prune(cartsfit, cp=cpTarget)
        cartsfit
    }
    return(cartsList)
}


# Nueral Network - SDM
snnFit <- function(commM, xData, train, iterNum, par){
    require(nnet)
    spNames <- colnames(commM)
    netList <- foreach(i=1:iterNum, .packages=c("nnet"), .verbose=F, .errorhandling="pass") %dopar% {
        iTrain <- train[,i]
        sub.commM <- subset(commM, iTrain)
        sub.xData <- subset(xData, iTrain)
        netfit <- lapply(spNames, FUN=function(m){nnet(x=sub.xData, y=sub.commM[,m], size=par[m,9], decay=par[m,10], linout=T)})
        netfit
    }
    return(netList)
}


# Nueral Network - CLM
mnnFit <- function(commM, xData, train, iterNum, para){
    require(nnet)
    netList <- foreach(i=1:iterNum, .packages=c("nnet"), .verbose=F, .errorhandling="pass") %dopar% {
        iTrain <- train[,i]
        sub.commM <- subset(commM, iTrain)
        sub.xData <- subset(xData, iTrain)
        netfit <- nnet(x=sub.xData, y=sub.commM, size=para[25,9], decay=para[25,10], linout=T)
        netfit
    }
    return(netList)
}


################################################################################
#
# FUNCTION TO EVALUATE SPECIES DISTRIBUTION IN THE TRAINING DATASET
#
# Function to evaluate models
# It takes a list of models (one for each iteration), the community matrix of
# observed values, the predictors matrix, a list with train/test splits, the
# type of response expected by the models, a model index for those cases in which
# models fitting could break and a logical value indicating whether the model is
# multiresponse or not. 
# The function give back list of matrices (one for each iteration) with AUC, TSS and# other metrics for each species
evalDistTrain <- function(modlist, commM, xData, train, respType, multiSp=T){
#modlist <- glmList
#commM <- pollenFit[[1]]
#xData <- climFit[[1]]
#train <- trainSetsFit[[1]]
#respType <- "response"
#multiSp <- F
#rm(modlist, commM, xData, train, respType, multiSp)
#rm(eDist, iTrain, modPred, trainDataset, testDataset, trainPredict, testPredict, edist)
    # Foreach loop to work with each model (on for each iteration)
    eDist <- foreach(i=1:length(modlist), .packages=c("paleoCLMs","VGAM", "dismo", "vegan", "betapart", "gam", "earth", "mvpart", "nnet"), .verbose=F, .errorhandling='pass') %dopar% {
        # modlist <- glmList[1]
        # commM <- pollenFit[[1]]
        # xData <- climFit[[1]]
        # train <- trainSetsFit[[1]]
        # respType <- "response"
        # multiSp <- FALSE
        
        # Get the current train/test splitting
        iTrain <- train[,i]
        
        # If statement to work with multiresponse models (CLMs)
        if(multiSp == TRUE){
            # Get species (communities) predictions to the predictor variables
            modPred <- predFun(xData, modlist[[i]], type=respType)
            # or with single species models (SDMs)
        }else{
            # Get species predictions to the predictor variables and assemblage them
            # in a community matrix
            modPred <- matrix(NA, ncol=ncol(commM), nrow=nrow(commM))
            for(j in 1:ncol(commM)){
                modPred[,j] <- predFun(xData, modlist[[i]][[j]], type=respType)
                colnames(modPred) <- colnames(commM)
            }
        }
        
        # Split the observed community matrix in training and testing
        trainDataset <- commM[iTrain,]
        testDataset <- commM[!iTrain,]
        
        # Split the predicted community matrix in training and testing
        trainPredict <- modPred[iTrain,]
        testPredict <- modPred[!iTrain,]
        
        # Evaluate the predictions with the evalDist function from the paleoCLMs package
        edist <- evalDist(obse=testDataset, pred=testPredict, obseTrain=trainDataset, predTrain=trainPredict, output="train")
        edist
    }
    
    # Return evaluation
    return(eDist)
}



################################################################################
#
# FUNCTION TO EVALUATE SPECIES DISTRIBUTION TO MULTIPLE PROJECTED TIME PERIODS
#
# Similar to the previous function but this time the models are evaluated projecting
# back on time to different time periods. The outputs is a list of lists (list for
# each iteration), in which each element is a matrix with AUC, TSS, etc for each
# species at each projected time period.
evalDistTest <- function(modlist, testCommM, trainCommM, xData, xDataTrain, train, trainTrain, respType, multiSp=T){
    #k is the time period the model was build in
    
    # Foreach loop to work with all the models (one for each iteration)
    eDist <- foreach(i=1:length(modlist), .packages=c("paleoCLMs","VGAM", "dismo", "vegan", "betapart", "gam", "earth", "mvpart", "nnet"), .verbose=F, .errorhandling='pass') %dopar% {
        
        # Getting the current train/test split
        iTrain <- lapply(train, function(x){x[,i]})
        iTrainTrain <- trainTrain[,i]
        
        # If statement to work with multiresponse models (CLMs) or ...
        if(multiSp == TRUE){
            # Get species (communities) predictions to the predictor variables
            modPred <- lapply(xData, FUN=predFun, modlist[[i]], type=respType)
            modPredTrain <- predFun(xDataTrain, modlist[[i]], type=respType)
            # with single species models (SDMs)
        }else{
            # Get species predictions to the predictor variables and assemblage them
            # in a community matrix. This step is repeated for each time period in which 
            # the models are projected to
            tempList <- list()
            for(l in 1:length(xData)){
                tempMat <- matrix(NA, ncol=ncol(testCommM[[l]]), nrow=nrow(xData[[l]]))
                for(j in 1:ncol(testCommM[[l]])){
                    tempMat[,j] <- predFun(xData[[l]], modlist[[i]][[j]], type=respType)
                }
                tempList[[l]] <- tempMat
            }
            modPred <- tempList
            modPredTrain <- matrix(NA, ncol=ncol(trainCommM), nrow=nrow(xDataTrain))
            for (j in 1:ncol(trainCommM)){
                modPredTrain[,j] <- predFun(xDataTrain, modlist[[i]][[j]], type=respType)
            }
        }
        
        # Splitting the observed community matrices in train/test
        trainDataset <- trainCommM[iTrainTrain,]
        testDataset <- mapply(function(x, y){x[!y,]}, testCommM, iTrain, SIMPLIFY=F)
        
        # The same for models predictions
        trainPredict <- modPredTrain[iTrainTrain,]
        testPredict <- mapply(function(x, y){x[!y,]}, modPred, iTrain, SIMPLIFY=F)
        
        # Evaluate the predictions with the evalDist function from the paleoCLMs package
        edist <- mapply(evalDist, obse=testDataset, pred=testPredict, MoreArgs=list(obseTrain=trainDataset, predTrain=trainPredict), SIMPLIFY=F, output="test")
        edist
    }
    
    # Return evaluation
    return(eDist)
}



################################################################################
#
# FUNCTION TO EVALUATE COMMUNITY ASSEMBLAGE TO THE TRAINING DATASET
#
# This function works similar to the first one but calculating metrics at the
# community level.
evalEnseTrain <- function(modList, commM, xData, train, respType, modDistEval, multiSp=T){
    # Foreach loop to work with all the models (one for each iteration)
    eEnse <- foreach(i=1:length(modList), .packages=c("paleoCLMs","VGAM", "dismo", "vegan", "betapart", "gam", "earth", "mvpart", "nnet"), .verbose=F, .errorhandling='pass') %dopar% {
        # modList <- cqoList[1]
        # commM <- pollenFit[[1]]
        # xData <- climFit[[1]]
        # train <- trainSetsFit[[1]]
        # respType <- "response"
        # modDistEval <- cqoDistTrain
        # multiSp <- TRUE
        # rm(modList, commM, xData, train, respType, modDistEval, multiSp)
        # rm(iTrain, modPred, eense1, tr, modPred, eense2, eense)
        
        # Getting the current train/test splitting
        iTrain <- train[,i]
        
        # If statement to work with multiresponse models (CLMs) or ...
        if(multiSp == TRUE){
            # Get species (communities) predictions to the predictor variables
            modPred <- predFun(xData, modList[[i]], type=respType)
            # with single species models (SDMs)
        }else{
            # Get species predictions to the predictor variables and assemblage them
            # in a community matrix.
            modPred <- matrix(NA, ncol=ncol(commM), nrow=nrow(commM))
            for(j in 1:ncol(commM)){
                modPred[,j] <- predFun(xData, modList[[i]][[j]], type=respType)
                colnames(modPred) <- colnames(commM)
            }
        }
        
        # If statement to fix an issue with the "rpart" package, that not assign species    # names to the prediction matrix
        if(class(modList[[1]]) == "rpart"){
            colnames(modPred) <- colnames(commM)
        }
        
        commM100 <- round(commM * 100)
        modPred100 <- round(modPred * 100)
        
        eense1 <- evalEnse(obse=commM100, pred=modPred100, tS=!iTrain, binary=FALSE)
        
        # Get the threshold that minimize specificity + sensitivity in the training
        # dataset for each species
        #nonErrorIter <- which(unlist(lapply(lapply(modDistEval, class), function(x){!any(x == "error")})))[[1]]
        #tr <- cbind(as.character(modDistEval[[nonErrorIter]][,"taxon"]), modDistEval[[nonErrorIter]][,"trmin"])
        tr <- modDistEval[[i]][,c("taxon","trmin")]
        # Apply the thresholds to the community matrix
        modPred <- applyThreshold(modPred, threshold=tr)
        # Evaluate the models with the evalEnse function from the "paleoCLMs" package
        eense2 <- evalEnse(obse=commM, pred=modPred, tS=!iTrain, binary=TRUE)
        
        eense <- c(eense1, eense2)
        names(eense) <- c("jacc.MEAN.prob","jacc.SD.prob","jacc.CAL.prob","jacc.DIS.prob","spRich.COR.prob","spRich.RSQ.prob","jacc.MEAN.thre","jacc.SD.thre","jacc.CAL.thre","jacc.DIS.thre","spRich.COR.thre","spRich.RSQ.thre")
        eense
    }
    
    # Return evaluation
    return(eEnse)
}




################################################################################
#
# FUNCTION TO EVALUATE COMMUNITY ASSEMBLAGE TO MULTIPLE PROJECTED TIME PERIODS
#
# This function works similar to the second one but calculating metrics at the
# community level.
evalEnseTest <- function(modList, testCommM, trainCommM, xData, train, respType, modDistEval=NULL, multiSp=T){    
    # Foreach loop to work with all the models (one for each iteration)
    eEnse <- foreach(i=1:length(modList), .packages=c("paleoCLMs","VGAM", "dismo", "vegan", "betapart", "gam", "earth", "mvpart", "nnet"), .verbose=F, .errorhandling='pass') %dopar% {
        # modList <- glmList
        # testCommM <- pollenPrj
        # trainCommM <- pollenFit[[1]]
        # xData <- climPrj
        # train <- trainSetsPrj
        # respType <- "response"
        # modDistEval <- glmDistTrain
        # multiSp <- F
        
        
        # Getting the current train/test split
        iTrain <- lapply(train, function(x){x[,i]})
        
        # If statement to work with multiresponse models (CLMs) or ...
        if(multiSp == TRUE){
            # Get species (communities) predictions to the predictor variables
            modPred <- lapply(xData, FUN=predFun, modList[[i]], type=respType)
            # with single species models (SDMs)
        }else{
            # Get species predictions to the predictor variables and assemblage them
            # in a community matrix.
            tempList <- list()
            for(l in 1:length(xData)){
                tempMat <- matrix(NA, ncol=ncol(testCommM[[l]]), nrow=nrow(xData[[l]]))
                for(j in 1:ncol(testCommM[[l]])){
                    tempMat[,j] <- predFun(xData[[l]], modList[[i]][[j]], type=respType)
                }
                colnames(tempMat) <- colnames(testCommM[[l]])
                tempList[[l]] <- tempMat
            }
            modPred <- tempList
        }  
        
        # If statement to fix an issue with the "rpart" package, that not assign species    # names to the prediction matrix
        if(class(modList[[1]]) == "rpart"){
            modPred <- lapply(modPred, FUN=function(x, y){colnames(x) <- colnames(y); return(x)}, trainCommM)
        }
        
        testCommM100 <- lapply(testCommM, FUN=function(x){round(x * 100)})
        modPred100 <- lapply(modPred, FUN=function(x){round(x * 100)})
        
        eense1 <- mapply(evalEnse, obse=testCommM100, pred=modPred100, tS=iTrain, binary=FALSE)
        
        # Get the threshold that minimize specificity + sensitivity in the training
        # dataset for each species
        tr <- modDistEval[[i]][,c("taxon","trmin")]
        # Apply the thresholds to the community matrix
        modPred <- lapply(modPred, FUN=applyThreshold, threshold=tr)
        # Evaluate the models with the evalEnse function from the "paleoCLMs" package
        eense2 <- mapply(evalEnse, obse=testCommM, pred=modPred, tS=iTrain, binary=TRUE)
        
        eense <- rbind(eense1, eense2)
        rownames(eense) <- c("jacc.MEAN.prob","jacc.SD.prob","jacc.CAL.prob","jacc.DIS.prob","spRich.COR.prob","spRich.RSQ.prob","jacc.MEAN.thre","jacc.SD.thre","jacc.CAL.thre","jacc.DIS.thre","spRich.COR.thre","spRich.RSQ.thre")
        eense
    }
    
    # Return evaluation
    return(eEnse)
}


################################################################################
#
# FUNCTION TO CALCULATE PROJECTED DISTRIBUTION RANGE TO DIFFERENT TIME PERIODS
#
# This function use the models to project to different timer periods, then using
# a species-model specific threshold would split the prediction in presence/absence
# estimation. Based in this estimation would compute two metrics of the distribution range (area and perimeter).
getDistRange <- function(modList, xStacks, modDistEval, respType, multiSp=T){
    
    modListShort <- modList
    
    # GET AVERAGE THRESHOLDS ACROSS ITERATIONS FOR EACH SPECIES
    # Foreach model to work with all the models. Here we get the specied-model specific thresholds
    eTr <- foreach(i=1:length(modListShort), .packages=c("paleoCLMs", "dismo", "vegan", "betapart", "SDMTools", "VGAM", "gam", "earth", "mvpart", "nnet"), .verbose=F, .errorhandling="pass") %dopar% {
        cbind(as.character(modDistEval[[i]][,"taxon"]), modDistEval[[i]][,"trmin"])
    }
    # Especify the column names for the thresholds matrices
    eTr <- lapply(eTr, FUN=function(x){colnames(x) <- c("taxon","trmin"); return(x)})
    # Combine them in a single dataframe from the previous list, and remove duplicated
    # columns ("taxon"), get numeric values instead of text, calculate the average and
    # get a final dataframe.
    eTr <- do.call(cbind, eTr)
    taxon <- eTr[,1]
    eTr <- eTr[,-c(1,3,5,7,9,11,13,15,17,19)]
    eTr <- matrix(as.numeric(eTr), ncol=10, nrow=17)
    eTr <- apply(eTr, MARGIN=1, FUN=mean)
    eTr <- data.frame(taxon, eTr)
    
    # Extract predictor values as dataframes from the raster stacks.
    xValues <- lapply(xStacks, getValues)
    xValues <- lapply(xValues, as.data.frame)
    
    # Get values only for land areas (remove seas and lakes)
    index <- lapply(xValues, function(x){which(!is.na(x[,1]))})
    xValuesNoNA <- mapply(function(x, y){x[y,]}, xValues, index, SIMPLIFY=F)
    
    # Foreach package to work with all the models
    eMod <- foreach(i=1:length(modListShort), .packages=c("paleoCLMs","VGAM", "dismo", "vegan", "betapart"), .verbose=T) %dopar% {
        
        # If statement to work with multiresponse models (CLMs) or ...
        if(multiSp == TRUE){
            # Get species (communities) predictions to the predictor variables
            modPred <- lapply(xValuesNoNA, FUN=predFun, modListShort[[i]], type=respType)
            # with single species models (SDMs)
        }else{
            # Get species predictions to the predictor variables and assemblage them
            # in a community matrix.
            tempList <- list()
            for(k in 1:length(xValuesNoNA)){
                tempMat <- matrix(NA, ncol=nrow(modDistEval[[i]]), nrow=nrow(xValuesNoNA[[k]]))
                for(j in 1:nrow(modDistEval[[i]])){
                    tempMat[,j] <- predFun(xValuesNoNA[[k]], modListShort[[i]][[j]], type=respType)
                }
                colnames(tempMat) <- modDistEval[[i]][,1]
                tempList[[k]] <- tempMat
            }
            modPred <- tempList
        }
        modPred
    }
    
    # For loop to combine predictions from different iterations of the same model
    # into a single average prediction.
    modPred <- list()
    for(j in 1:length(xStacks)){
        t1 <- as.data.frame(eMod[[1]][[j]])
        t2 <- as.data.frame(eMod[[2]][[j]])
        t3 <- as.data.frame(eMod[[3]][[j]])
        t4 <- as.data.frame(eMod[[4]][[j]])
        t5 <- as.data.frame(eMod[[5]][[j]])
        t6 <- as.data.frame(eMod[[6]][[j]])
        t7 <- as.data.frame(eMod[[7]][[j]])
        t8 <- as.data.frame(eMod[[8]][[j]])
        t9 <- as.data.frame(eMod[[9]][[j]])
        t10 <- as.data.frame(eMod[[10]][[j]])
        allData <- list(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10)
        modPred[[j]] <- Reduce(`+`, allData) / length(allData)
    }
    
    # Apply the average thresholds computed previously (eTr) to the average prediction
    modPred <- lapply(modPred, applyThreshold, eTr)
    
    # Project the predictions to rasters, by including NAs in the dataframe where 
    # there is water (seas or lakes) and rasterize these values.
    predValues <- matrix(NA, ncol=17, nrow=nrow(xValues[[1]]))
    predValues <- mapply(function(x,y,z){x[z,] <- y; return(x)}, y=modPred, z=index, MoreArgs=list(x=predValues), SIMPLIFY=F)
    predStacks <- mapply(setValues, xStacks, predValues)
    
    # Calculate distribution range statistics with the "PatchStat" function in the "SDMTools" package
    require(SDMTools)
    distRange <- lapply(predStacks, FUN=function(x){lapply(unstack(x), FUN=function(y){PatchStat(y, latlon=T)})})
    
    # Extract information regarding area and perimeter, and disregard all other statistics
    distRange <- lapply(distRange, FUN=function(x){lapply(x, FUN=function(y){y[2,c("area","perimeter")]})})
    
    return(distRange)
}


