################################################################################
##
## Wrapping functions to run models and save code afterwards
##
################################################################################

# GLM - SDM
glmFit <- function(commM, xData, train, iterNum){
  spNames <- colnames(commM)
  glmList <- foreach(i=1:iterNum, .verbose=T) %dopar% {
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
  cqoList <- foreach(i=1:numIter, .packages=c("VGAM"), .verbose=T, .errorhandling='pass') %dopar% {
    iTrain <- train[,i]
    sub.commM <- subset(commM, iTrain)
    sub.xData <- subset(xData, iTrain)
    form <- formula(sub.commM ~ etr_year_ave + aet_year_ave + wdei_year_ave + tmax_high_quart + prcp_low_quart + prcp_high_quart)
    cqofit <- cqo(form, family=binomialff(mv=TRUE), data=sub.xData, trace=TRUE, Rank=2, maxit=100, step=0.5)
    cqofit
  }
  return(cqoList)
}

# Nueral Network - SDM
snnFit <- function(commM, xData, train, iterNum){
  require(nnet)
  spNames <- colnames(commM)
  netList <- foreach(i=1:numIter, .packages=c("nnet"), .verbose=T) %dopar% {
    iTrain <- train[,i]
    sub.commM <- subset(commM, iTrain)
    sub.xData <- subset(xData, iTrain)
    netfit <- lapply(spNames, FUN=function(m){nnet(x=sub.xData, y=sub.commM[,m], size=15, linout=T)})
    netfit
  }
  return(netList)
}


# Nueral Network - CLM
mnnFit <- function(commM, xData, train, iterNum){
  require(nnet)
  netList <- foreach(i=1:numIter, .packages=c("nnet"), .verbose=T) %dopar% {
    iTrain <- train[,i]
    sub.commM <- subset(commM, iTrain)
    sub.xData <- subset(xData, iTrain)
    netfit <- nnet(x=sub.xData, y=sub.commM, size=15, linout=T)
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
#
#################################################################################
evalDistTrain <- function(modlist, commM, xData, train, respType, modelIndex=NULL, multiSp=T){

  # If modelIndex not specified then all the models are included
  if(is.null(modelIndex)){
    modelIndex <- 1:length(modlist)
  }

  # Foreach loop to work with each model (on for each iteration)
  eDist <- foreach(i=1:length(modlist), .packages=c("paleoCLMs","VGAM", "dismo", "vegan", "betapart", "gam", "earth", "mvpart", "nnet"), .verbose=T, .errorhandling='pass') %dopar% {

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

  # Return only evaluation for those models that are specified by modelIndex
  return(eDist[modelIndex])
}



################################################################################
#
# FUNCTION TO EVALUATE SPECIES DISTRIBUTION TO MULTIPLE PROJECTED TIME PERIODS
#
# Similar to the previous function but this time the models are evaluated projecting
# back on time to different time periods. The outputs is a list of lists (list for
# each iteration), in which each element is a matrix with AUC, TSS, etc for each
# species at each projected time period.
#
################################################################################

evalDistTest <- function(modlist, commM, xData, train, respType, modelIndex=NULL, multiSp=T){

  # If modelIndex not specified then all the models are included
  if(is.null(modelIndex)){
    modelIndex <- 1:length(modlist)
  }

  # Foreach loop to work with all the models (one for each iteration)
  eDist <- foreach(i=1:length(modlist), .packages=c("paleoCLMs","VGAM", "dismo", "vegan", "betapart", "gam", "earth", "mvpart", "nnet"), .verbose=T, .errorhandling='pass') %dopar% {

    # Getting the current train/test split
    iTrain <- lapply(train, function(x){x[,i]})
    
    # If statement to work with multiresponse models (CLMs) or ...
    if(multiSp == TRUE){
      # Get species (communities) predictions to the predictor variables
      modPred <- lapply(xData, FUN=predFun, modlist[[i]], type=respType)
    # with single species models (SDMs)
    }else{
      # Get species predictions to the predictor variables and assemblage them
      # in a community matrix. This step is repeated for each time period in which 
      # the models are projected to
      tempList <- list()
      for(k in 1:length(xData)){
        tempMat <- matrix(NA, ncol=ncol(commM[[k]]), nrow=nrow(xData[[k]]))
        for(j in 1:ncol(commM[[k]])){
          tempMat[,j] <- predFun(xData[[k]], modlist[[i]][[j]], type=respType)
        }
        tempList[[k]] <- tempMat
      }
    modPred <- tempList
    }

    # Splitting the observed community matrices in train/test. Here we use 0BP
    # as the training dataset, and other time periods as test
    trainDataset <- commM[[1]][iTrain[[1]],]
    testDataset <- commM

    # The same for models predictions
    trainPredict <- modPred[[1]][iTrain[[1]],]
    testPredict <- modPred
    
    # We change the 0BP matrix to correctly being splitted between training and testing, since the others timer periods are all records for testing.
    testDataset[[1]] <- commM[[1]][!iTrain[[1]],]
    testPredict[[1]] <- modPred[[1]][!iTrain[[1]],]

    # Evaluate the predictions with the evalDist function from the paleoCLMs package
    edist <- mapply(evalDist, obse=testDataset, pred=testPredict, MoreArgs=list(obseTrain=trainDataset, predTrain=trainPredict), SIMPLIFY=F, output="test")
    edist
  }

  # Return only evaluation for those models that are specified by modelIndex
  return(eDist[modelIndex])
}



################################################################################
#
# FUNCTION TO EVALUATE COMMUNITY ASSEMBLAGE TO THE TRAINING DATASET
#
# This function works similar to the first one but calculating metrics at the
# community level.
#
###############################################################################

evalEnseTrain <- function(modList, commM, xData, train, respType, modDistEval, modelIndex=NULL, multiSp=T){
  
  # If modelIndex not specified then all the models are included
  if(is.null(modelIndex)){
    modelIndex <- 1:length(modList)
  }

  # Foreach loop to work with all the models (one for each iteration)
  eEnse <- foreach(i=1:length(modList), .packages=c("paleoCLMs","VGAM", "dismo", "vegan", "betapart", "gam", "earth", "mvpart", "nnet"), .verbose=T, .errorhandling='pass') %dopar% {
  
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
    
    # Get the threshold that minimize specificity + sensitivity in the training
    # dataset for each species
    tr <- cbind(as.character(modDistEval[[1]][,"taxon"]), modDistEval[[1]][,"trmin"])
    
    # Apply the thresholds to the community matrix
    modPred <- applyThreshold(modPred, threshold=tr)
    
    # Evaluate the models with the evalEnse function from the "paleoCLMs" package
    eense <- evalEnse(obse=commM, pred=modPred, tS=!iTrain)
    eense
  }

  # Return only evaluation for those models that are specified by modelIndex
  return(eEnse[modelIndex])
}




################################################################################
#
# FUNCTION TO EVALUATE COMMUNITY ASSEMBLAGE TO MULTIPLE PROJECTED TIME PERIODS
#
# This function works similar to the second one but calculating metrics at the
# community level.
#
###############################################################################

evalEnseTest <- function(modList, commM, xData, train, respType, modDistEval, modelIndex=NULL, multiSp=T){

  # If modelIndex not specified then all the models are included
  if(is.null(modelIndex)){
    modelIndex <- 1:length(modList)
  }

  # Foreach loop to work with all the models (one for each iteration)
  eEnse <- foreach(i=1:length(modList), .packages=c("paleoCLMs","VGAM", "dismo", "vegan", "betapart", "gam", "earth", "mvpart", "nnet"), .verbose=T, .errorhandling='pass') %dopar% {
  
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
      for(k in 1:length(xData)){
        tempMat <- matrix(NA, ncol=ncol(commM[[k]]), nrow=nrow(xData[[k]]))
        for(j in 1:ncol(commM[[k]])){
          tempMat[,j] <- predFun(xData[[k]], modList[[i]][[j]], type=respType)
        }
        colnames(tempMat) <- colnames(commM[[k]])
        tempList[[k]] <- tempMat
      }
      modPred <- tempList
    }

    # If statement to fix an issue with the "rpart" package, that not assign species    # names to the prediction matrix
    if(class(modList[[1]]) == "rpart"){
      modPred <- lapply(modPred, FUN=function(x, y){colnames(x) <- colnames(y); return(x)}, commM[[1]])
    }

    # Get the threshold that minimize specificity + sensitivity in the training
    # dataset for each species
    tr <- cbind(as.character(modDistEval[[1]][,"taxon"]), modDistEval[[1]][,"trmin"])

    # Apply the thresholds to the community matrix
    modPred <- lapply(modPred, FUN=applyThreshold, threshold=tr)

    # Evaluate the models with the evalEnse function from the "paleoCLMs" package
    eense <- mapply(evalEnse, obse=commM, pred=modPred, tS=iTrain)
    eense
  }

  # Return only evaluation for those models that are specified by modelIndex
  return(eEnse[modelIndex])
}


################################################################################
#
# FUNCTION TO CALCULATE PROJECTED DISTRIBUTION RANGE TO DIFFERENT TIME PERIODS
#
# This function use the models to project to different timer periods, then using
# a species-model specific threshold would split the prediction in presence/absence
# estimation. Based in this estimation would compute two metrics of the distribution range (area and perimeter).
#
################################################################################

getDistRange <- function(modList, xStacks, modDistEval, respType, modelIndex=NULL, multiSp=T){
  
  # If modelIndex not specified then all the models are included
  if(is.null(modelIndex)){
    modelIndex <- 1:length(modList)
  }

  # Limit the work to those models included in modelIndex (success models)  
  modListShort <- modList[modelIndex]
  
  # GET AVERAGE THRESHOLDS ACROSS ITERATIONS FOR EACH SPECIES
  # Foreach model to work with all the models included in the modelIndex. Here we get the specied-model specific thresholds
  eTr <- foreach(i=1:length(modListShort), .packages=c("paleoCLMs", "dismo", "vegan", "betapart", "SDMTools", "VGAM", "gam", "earth", "mvpart", "nnet"), .verbose=T) %dopar% {
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

  # Foreach package to work with all the models especified in modelIndex
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


