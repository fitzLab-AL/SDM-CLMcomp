################################################################################
#
# LOAD DATA, PARAMETERS, AND FUNCTIONS
#
################################################################################

# Load the project library with functions to load and arrange data
library(paleoCLMs)

# Load wrapping functions
source("../R-package/00-WrappingFunctions.R")

# Define parameters
indir <- "../Data/" 

periods <- seq(0, 21000, by=500) # That generate a numerical vector with time slides. I use this afterwards to name files and objects
periods_red <- periods/1000
index <- 1:length(periods)

vars <- c("etr_year_ave","aet_year_ave","wdei_year_ave","tmax_high_quart","prcp_low_quart","prcp_high_quart") # This select the name of the variables that will be used to fit the models

numIter <- 10 # Here I specified how many iteractions will be runned. In this case I fitted 10 models with different random train-test datasets.

# Load pollen data and remove low quality data
pollenOriginal <- lapply(periods, loadPollen, indir)
pollenHQ <- lapply(pollenOriginal, removeLowQualitySamples, 0.75)

# Load a raster as template, it will be deleted a bit later.
rasTemp <- raster("../Data/Climate/Climate/CCSM/0BP/gdd5_year_ave.tif")

# Average values for those pollen sites in the same grid cell
pollenAbun <- lapply(pollenHQ, reduceDuplicated, rasTemp, weighted=F)
rm(rasTemp, pollenOriginal, pollenHQ)

# Load climate data
clim <- mapply(loadClim, period=periods, pollen=pollenAbun, MoreArgs=list(clim_model="CCSM", indir=paste(indir, "Climate/", sep=""), vars=vars))

clim <- lapply(clim, as.data.frame)

# Remove incomplete cases from the pollen and the climate datasets
compCasesIndex <- lapply(clim, FUN=function(x){which(complete.cases(x))})
pollenAbun <- mapply(FUN=function(x,i){x[i,]}, pollenAbun, compCasesIndex, SIMPLIFY=F)
clim <- lapply(clim, FUN=function(x){x[which(complete.cases(x)),]})
rm(compCasesIndex)

# Extract coordinates to a new object and remove them from the pollen data.frame
coord <- lapply(pollenAbun, FUN=function(x){x[,which(colnames(x) %in% c("x","y"))]})
pollenAbun <- lapply(pollenAbun, FUN=function(x){x[,-which(colnames(x) %in% c("x","y"))]})

# Convert pollen concentrations to presences/absences matrix given a specific pollen threshold(Community matrix)
# Calculate variable thresholds
var05 <- lapply(pollenAbun, pollenThreshold, 0.05)

# apply various thresholds
pollen <- mapply(applyThreshold, pollenAbun, threshold=var05, SIMPLIFY=F)
rm(var05)

# Get a taxon list for training the models and remove all species with less than 20 presences and/or absences. Then propagate the taxon list to all time periods
taxonNames <- colnames(removeAbsentSpecies(pollen[[1]], 19, applyTo="both"))

#remove taxa because too few datapoints
taxonNames <- taxonNames[-which(taxonNames %in% c("Rumex.Oxyria","Castanea","Tilia","Tsuga","Fagus"))]

# Propagate the taxon list to all pollen datasets
pollenAbun <- lapply(pollenAbun, FUN=function(x){x[,taxonNames]})
pollen <- lapply(pollen, FUN=function(x){x[,taxonNames]})

# Define time periods to fit and project the models
periods_fit <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "11.5", "12", "12.5", "13", "13.5", "14-15", "15.5-17.5", "15.5-21")
periods_prj <- c("0", "0.5", "3", "6", "7.5", "9.5", "12.5", "14.5", "15", "15.5-17.5", "15.5-21", "18-21")

index_fit <- as.list(which(periods_red %in% periods_fit))
index_fit <- append(index_fit, list(c(29:31), c(32:36), c(32:43))) 
index_prj <- as.list(which(periods_red %in% periods_prj))
index_prj <- append(index_prj, list(c(32:36), c(32:43), c(37:43)))


###############################################################################
##
## FIT THE MODELS
##
###############################################################################

seed <- c(1, 11, 21, 31, 41)

require(foreach)
require(doParallel)
registerDoParallel(cores=5)

## Some of the cqo and cao models crash at certain time periods. I run the cqo and cao models twice (Replicate1 and Replicate2) by hand. Then I replaced those that crashed in the first run with some from the second run. See the code below for merging the models

for(seed.i in 1:5){
  for(i in 1:length(periods_fit)){ # Only until 21 because the 0-21 models have to be run multiple times (see code below)
    j=index_fit[[i]]
    k=periods_fit[[i]]
    
    a <- Sys.time()
    print(paste("Seed:", seed[seed.i], sep=""))
    print(paste("Fitting period: ", k, "BP", sep=""))
    
    # Split training and testing datasets. THAT SHOULD BE BEFORE POOLING DATA
    trainSets <- lapply(pollen, trainSplit, numIter, 0, seed[seed.i])
    trainSets[j] <- lapply(pollen[j], trainSplit, numIter, 0.7, seed[seed.i])
    
    # Pooling data together
    # Binning 14-15
    binIndex <- which(periods %in% seq(14000, 15000,by=500))
    pbin1 <- poolingData(binIndex, pollen)
    cbin1 <- poolingData(binIndex, clim)
    tbin1 <- poolingData(binIndex, trainSets)
    
    # Binning 15.5 to 17.5
    binIndex <- which(periods %in% seq(15500, 17500, by=500))
    pbin2 <- poolingData(binIndex, pollen)
    cbin2 <- poolingData(binIndex, clim)
    tbin2 <- poolingData(binIndex, trainSets)
    
    # Binning 15.5 to 21
    binIndex <- which(periods %in% seq(15500, 21000, by=500))
    pbin3 <- poolingData(binIndex, pollen)
    cbin3 <- poolingData(binIndex, clim)
    tbin3 <- poolingData(binIndex, trainSets)
    
    # Binning 18 to 21
    binIndex <- which(periods %in% seq(18000, 21000, by=500))
    pbin4 <- poolingData(binIndex, pollen)
    cbin4 <- poolingData(binIndex, clim)
    tbin4 <- poolingData(binIndex, trainSets)
    
    # Appending the pooled dataframes to the data lists
    pollenPool <- append(pollen, list(pbin1, pbin2, pbin3, pbin4))
    climPool <- append(clim, list(cbin1, cbin2, cbin3, cbin4))
    trainSetsPool <- append(trainSets, list(tbin1, tbin2, tbin3, tbin4))
    
    pollenPool <- lapply(pollenPool, FUN=function(x){rownames(x) <- 1:nrow(x); return(x)})
    climPool <- lapply(climPool, FUN=function(x){rownames(x) <- 1:nrow(x); return(x)})
    trainSetsPool <- lapply(trainSetsPool, FUN=function(x){rownames(x) <- 1:nrow(x); return(x)})
  
    # Remove temporal objects
    rm(binIndex, cbin1, cbin2, cbin3, cbin4, pbin1, pbin2, pbin3, pbin4, tbin1, tbin2, tbin3, tbin4)
    
    # Change names
    names(pollenPool) <- names(climPool) <- names(trainSetsPool) <- c(periods_red, "14-15", "15.5-17.5", "15.5-21", "18-21")
    
    # Selecting and splitting data for fitting and projecting time periods
    pollenFit <- pollenPool[which(names(pollenPool) %in% periods_fit)]
    climFit <- climPool[which(names(climPool) %in% periods_fit)]
    trainSetsFit <- trainSetsPool[which(names(trainSetsPool) %in% periods_fit)]
    
    pollenPrj <- pollenPool[which(names(pollenPool) %in% periods_prj)]
    climPrj <- climPool[which(names(climPool) %in% periods_prj)]
    trainSetsPrj <- trainSetsPool[which(names(trainSetsPool) %in% periods_prj)]
    
    # Load tunned parameters for each taxon and model
    tunedParam <- read.csv("parameters.csv", header=TRUE, row.names=1)
    
    glmList <- glmFit(pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], numIter)
    cqoList <- cqoFit(pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], numIter)
    gamList <- gamFit(pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], numIter)
    caoList <- caoFit(pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], numIter)
    smarsList <- marsFit.sp(pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], numIter, tunedParam)
    marsList <- marsFit(pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], numIter, tunedParam)
    scartsList <- cartsFit.sp(pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], numIter)
    cartsList <- cartsFit(pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], numIter)
    snnList <- snnFit(data.frame(pollenFit[[i]]), climFit[[i]], trainSetsFit[[i]], numIter, tunedParam)
    mnnList <- mnnFit(pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], numIter, tunedParam)
  
    # Create directory to save models
    direct <- paste("Results/RObjects/Models/Replicate", seed.i, "/", sep="")
    dir.create(direct, recursive=T) 
  
    save(glmList, cqoList, gamList, caoList, scartsList, cartsList, smarsList, marsList, snnList, mnnList, file=paste(direct, periods_fit[[i]], "BP.RData", sep=""))
  
    b <- Sys.time()
    print(b-a)
  
    rm(j, k, a, trainSets, pollenPool, climPool, trainSetsPool, pollenFit, climFit, trainSetsFit, pollenPrj, climPrj, trainSetsPrj, tunedParam, direct, glmList, cqoList, gamList, caoList, smarsList, marsList, scartsList, cartsList, snnList, mnnList, b)
  }
}
rm(seed.i, i)


################################################################################
#
#  EVALUATE THE MODELS
#
################################################################################
require(foreach)
require(doParallel)
registerDoParallel(cores=5)

# Define which models should be evaluated in each replica
for(seed.i in 1:length(seed)){
  modelList <- 1:length(periods_fit)

  # Print some info for tracking the calculation
  print(paste("Seed:", seed[seed.i], sep=""))
  # FOR loop for each time period  
  for(i in modelList){
   
    # Define a couple of parameters
    j=index_fit[[i]]
    k=periods_fit[[i]]
    
    # Print some info for tracking the calculation
    a <- Sys.time()
    print(paste("  Fitting period: ", k, "BP", sep=""))
    
    # Split training and testing datasets. THAT SHOULD BE BEFORE POOLING DATA
    trainSets <- lapply(pollen, trainSplit, numIter, 0, seed[seed.i])
    trainSets[j] <- lapply(pollen[j], trainSplit, numIter, 0.7, seed[seed.i])
    
    # Pooling data together
    # Binning 14-15
    binIndex=which(periods %in% seq(14000, 15000,by=500))
    pbin1 <- poolingData(binIndex, pollen)
    cbin1 <- poolingData(binIndex, clim)
    tbin1 <- poolingData(binIndex, trainSets)
    
    # Binning 15.5 to 17.5
    binIndex=which(periods %in% seq(15500, 17500, by=500))
    pbin2 <- poolingData(binIndex, pollen)
    cbin2 <- poolingData(binIndex, clim)
    tbin2 <- poolingData(binIndex, trainSets)
    
    # Binning 15.5 to 21
    binIndex=which(periods %in% seq(15500, 21000, by=500))
    pbin3 <- poolingData(binIndex, pollen)
    cbin3 <- poolingData(binIndex, clim)
    tbin3 <- poolingData(binIndex, trainSets)
    
    # Binning 18 to 21
    binIndex=which(periods %in% seq(18000, 21000, by=500))
    pbin4 <- poolingData(binIndex, pollen)
    cbin4 <- poolingData(binIndex, clim)
    tbin4 <- poolingData(binIndex, trainSets)
    
    # Appending the pooled dataframes to the data lists
    pollenPool <- append(pollen, list(pbin1, pbin2, pbin3, pbin4))
    climPool <- append(clim, list(cbin1, cbin2, cbin3, cbin4))
    trainSetsPool <- append(trainSets, list(tbin1, tbin2, tbin3, tbin4))

    pollenPool <- lapply(pollenPool, FUN=function(x){rownames(x) <- 1:nrow(x); return(x)})
    climPool <- lapply(climPool, FUN=function(x){rownames(x) <- 1:nrow(x); return(x)})
    trainSetsPool <- lapply(trainSetsPool, FUN=function(x){rownames(x) <- 1:nrow(x); return(x)})

    # Remove temporal files
    rm(binIndex, cbin1, cbin2, cbin3, cbin4, pbin1, pbin2, pbin3, pbin4, tbin1, tbin2, tbin3, tbin4)
    
    # Names the element of the list
    names(pollenPool) <- names(climPool) <- names(trainSetsPool) <- c(periods_red, "14-15", "15.5-17.5", "15.5-21", "18-21")
    
    # Selecting and splitting data for fitting and projecting time periods
    pollenFit <- pollenPool[which(names(pollenPool) %in% periods_fit)]
    climFit <- climPool[which(names(climPool) %in% periods_fit)]
    trainSetsFit <- trainSetsPool[which(names(trainSetsPool) %in% periods_fit)]
    
    pollenPrj <- pollenPool[which(names(pollenPool) %in% periods_prj)]
    climPrj <- climPool[which(names(climPool) %in% periods_prj)]
    trainSetsPrj <- trainSetsPool[which(names(trainSetsPool) %in% periods_prj)]
  
    # Load the models
    direct <- paste("Results/RObjects/Models/Replicate", seed.i, "/", sep="")
    load(paste(direct, periods_fit[[i]], "BP.RData", sep=""), .GlobalEnv)
    
    ## EVALUATE THE MODELS
    # GLM - SDM
#    glmDistTrain <- evalDistTrain(glmList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="response", multiSp=F)
#    glmDistTest <- evalDistTest(glmList, pollenPrj, pollenFit[[i]], climPrj, climFit[[i]], trainSetsPrj, trainSetsFit[[i]], respType="response", multiSp=F)
    glmEnseTrain <- evalEnseTrain(glmList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="response", glmDistTrain, multiSp=F)
    glmEnseTest <- evalEnseTest(glmList, pollenPrj, pollenFit[[i]], climPrj, trainSetsPrj, respType="response", glmDistTrain, multiSp=F)
     
    # CQO - CLM
#    cqoDistTrain<- evalDistTrain(cqoList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="response")
#    cqoDistTest <- evalDistTest(cqoList, pollenPrj, pollenFit[[i]], climPrj, climFit[[i]], trainSetsPrj, trainSetsFit[[i]], respType="response")
    cqoEnseTrain <- evalEnseTrain(cqoList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="response", cqoDistTrain)
    cqoEnseTest <- evalEnseTest(cqoList, pollenPrj, pollenFit[[i]], climPrj, trainSetsPrj, respType="response", cqoDistTrain)
    
    # GAM - SDM
#    gamDistTrain <- evalDistTrain(gamList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="response", multiSp=F)
#    gamDistTest <- evalDistTest(gamList, pollenPrj, pollenFit[[i]], climPrj, climFit[[i]], trainSetsPrj, trainSetsFit[[i]], respType="response", multiSp=F)
    gamEnseTrain <- evalEnseTrain(gamList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="response", gamDistTrain, multiSp=F)
    gamEnseTest <- evalEnseTest(gamList, pollenPrj, pollenFit[[i]], climPrj, trainSetsPrj, respType="response", gamDistTrain, multiSp=F)
     
    # CAO - CLM
#    caoDistTrain<- evalDistTrain(caoList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="response")
#    caoDistTest <- evalDistTest(caoList, pollenPrj, pollenFit[[i]], climPrj, climFit[[i]], trainSetsPrj, trainSetsFit[[i]], respType="response")
    caoEnseTrain <- evalEnseTrain(caoList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="response", caoDistTrain)
    caoEnseTest <- evalEnseTest(caoList, pollenPrj, pollenFit[[i]], climPrj, trainSetsPrj, respType="response", caoDistTrain)

    # MARS - SDM
#    smarsDistTrain <- evalDistTrain(smarsList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="response", multiSp=F)
#    smarsDistTest <- evalDistTest(smarsList, pollenPrj, pollenFit[[i]], climPrj, climFit[[i]], trainSetsPrj, trainSetsFit[[i]], respType="response", multiSp=F)
    smarsEnseTrain <- evalEnseTrain(smarsList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="response", smarsDistTrain, multiSp=F)
    smarsEnseTest <- evalEnseTest(smarsList, pollenPrj, pollenFit[[i]], climPrj, trainSetsPrj, respType="response", smarsDistTrain, multiSp=F)
    
    # MARS - CLM
#    marsDistTrain <- evalDistTrain(marsList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="response")
#    marsDistTest <- evalDistTest(marsList, pollenPrj, pollenFit[[i]], climPrj, climFit[[i]], trainSetsPrj, trainSetsFit[[i]], respType="response")
    marsEnseTrain <- evalEnseTrain(marsList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="response", marsDistTrain)
    marsEnseTest <- evalEnseTest(marsList, pollenPrj, pollenFit[[i]], climPrj, trainSetsPrj, respType="response", marsDistTrain)
     
    # CARTs - SDM
#    scartsDistTrain <- evalDistTrain(scartsList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="matrix", multiSp=F)
#    scartsDistTest <- evalDistTest(scartsList, pollenPrj, pollenFit[[i]], climPrj, climFit[[i]], trainSetsPrj, trainSetsFit[[i]], respType="matrix", multiSp=F)
    scartsEnseTrain <- evalEnseTrain(scartsList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="matrix", scartsDistTrain, multiSp=F)
    scartsEnseTest <- evalEnseTest(scartsList, pollenPrj, pollenFit[[i]], climPrj, trainSetsPrj, respType="matrix", scartsDistTrain, multiSp=F)
    
    # CARTs - CLM
#    cartsDistTrain <- evalDistTrain(cartsList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="matrix")
#    cartsDistTest <- evalDistTest(cartsList, pollenPrj, pollenFit[[i]], climPrj, climFit[[i]], trainSetsPrj, trainSetsFit[[i]], respType="matrix")
    cartsEnseTrain <- evalEnseTrain(cartsList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="matrix", cartsDistTrain)
    cartsEnseTest <- evalEnseTest(cartsList, pollenPrj, pollenFit[[i]], climPrj, trainSetsPrj, respType="matrix", cartsDistTrain)
    
    # NEURAL NETWORKS - SDM 
    library (nnet)
#    snnDistTrain <- evalDistTrain(snnList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="raw", multiSp=F)
#    snnDistTest <- evalDistTest(snnList, pollenPrj, pollenFit[[i]], climPrj, climFit[[i]], trainSetsPrj, trainSetsFit[[i]], respType="raw", multiSp=F)
    snnEnseTrain <- evalEnseTrain(snnList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="raw", snnDistTrain, multiSp=F)
    snnEnseTest <- evalEnseTest(snnList, pollenPrj, pollenFit[[i]], climPrj, trainSetsPrj, respType="raw", snnDistTrain, multiSp=F)
     
    # NEURAL NETWORKS - CLM 
#    mnnDistTrain <- evalDistTrain(mnnList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="raw")
#    mnnDistTest <- evalDistTest(mnnList, pollenPrj, pollenFit[[i]], climPrj, climFit[[i]], trainSetsPrj, trainSetsFit[[i]], respType="raw")
    mnnEnseTrain <- evalEnseTrain(mnnList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="raw", mnnDistTrain)
    mnnEnseTest <- evalEnseTest(mnnList, pollenPrj, pollenFit[[i]], climPrj, trainSetsPrj, respType="raw", mnnDistTrain)

    # Save the RObjects in their respective folders
    # Distribution
#    direct <- paste("Results/RObjects/Evaluation/Distribution/Replicate", seed.i, "/", sep="")
#    dir.create(direct, recursive=T) 
#    save(glmDistTrain, glmDistTest, cqoDistTrain, cqoDistTest, gamDistTrain, gamDistTest, caoDistTrain, caoDistTest, scartsDistTrain, scartsDistTest, cartsDistTrain, cartsDistTest, smarsDistTrain, smarsDistTest, marsDistTrain, marsDistTest, snnDistTrain, snnDistTest, mnnDistTrain, mnnDistTest, file=paste(direct, periods_fit[[i]], "BP.RData", sep=""))
   
    # Assemblage
    direct <- paste("Results/RObjects/Evaluation/Assemblage/Replicate", seed.i, "/", sep="")
    dir.create(direct, recursive=T) 
    save(glmEnseTrain, glmEnseTest, cqoEnseTrain, cqoEnseTest, gamEnseTrain, gamEnseTest, caoEnseTrain, caoEnseTest, scartsEnseTrain, scartsEnseTest, cartsEnseTrain, cartsEnseTest, smarsEnseTrain, smarsEnseTest, marsEnseTrain, marsEnseTest, snnEnseTrain, snnEnseTest, mnnEnseTrain, mnnEnseTest, file=paste(direct, periods_fit[[i]], "BP.RData", sep=""))

    # Remove models to avoid problems
    rm(glmList, cqoList, gamList, caoList, scartsList, cartsList, smarsList, marsList, snnList, mnnList)
#    rm(glmDistTrain, cqoDistTrain, gamDistTrain, caoDistTrain, scartsDistTrain, cartsDistTrain, smarsDistTrain, marsDistTrain, snnDistTrain, mnnDistTrain)
#    rm(glmDistTest, cqoDistTest, gamDistTest, caoDistTest, scartsDistTest, cartsDistTest, smarsDistTest, marsDistTest, snnDistTest, mnnDistTest)
    rm(glmEnseTrain, cqoEnseTrain, gamEnseTrain, caoEnseTrain, scartsEnseTrain, cartsEnseTrain, smarsEnseTrain, marsEnseTrain, snnEnseTrain, mnnEnseTrain)
    rm(glmEnseTest, cqoEnseTest, gamEnseTest, caoEnseTest, scartsEnseTest, cartsEnseTest, smarsEnseTest, marsEnseTest, snnEnseTest, mnnEnseTest)

  }
  rm(j, k, a, trainSets, pollenPool, climPool, trainSetsPool, pollenFit, climFit, trainSetsFit, pollenPrj, climPrj, trainSetsPrj, direct)
}

rm(seed.i, i)


################################################################################
#
# MERGE RESULTS FROM THE DIFFERENTE REPLICATES
#
################################################################################

###############
## FUNCTIONS ##
###############

## Get the crashes index for the first random seed
getCrashes <- function(p.fit
                       , mod.names
                       , eval.dis.dir
                       , eval.ass.dir
                       )
  {
    #p.fit <- "9BP"
    #mod.names <- modNames
    #eval.dis.dir <- evalDisDir
    #eval.ass.dir <- evalAssDir
    #rm(p.fit, mod.names, eval.dis.dir, eval.ass.dir)
  
    load(paste(eval.dis.dir, "1/", p.fit, ".RData", sep=""))
    load(paste(eval.ass.dir, "1/", p.fit, ".RData", sep=""))
    crashesID <- c()    
    for(i in 1:length(mod.names))
    {
      nam <- mod.names[i]
    
      disTrainID <- which(sapply(get(paste(nam, "DistTrain", sep="")), FUN=function(x){any(class(x) == "error")}))
      disTestID <- which(sapply(get(paste(nam, "DistTest", sep="")), FUN=function(x){any(class(x) == "error")}))
      assTrainID <- which(sapply(get(paste(nam, "EnseTrain", sep="")), FUN=function(x){any(class(x) == "error")}))
      assTestID <- which(sapply(get(paste(nam, "EnseTest", sep="")), FUN=function(x){any(class(x) == "error")}))
    
      crashesID <- append(crashesID, c(disTrainID, disTestID, assTrainID, assTestID))
      crashesID <- unique(crashesID)
    }
    return(crashesID)
  }

## Get the crashes index for the first random seed
getSubstitutes <- function(p.fit
                       , mod.names
                       , eval.dis.dir
                       , eval.ass.dir
                       )
  {
    substitutesList <- list()
    for(j in 2:5){
      load(paste(eval.dis.dir, j, "/", p.fit,".RData", sep=""))
      load(paste(eval.ass.dir, j, "/", p.fit,".RData", sep=""))
      substID <- c()    
      for(i in 1:length(mod.names))
      {
        nam <- mod.names[i]
      
        disTrainID <- which(sapply(get(paste(nam, "DistTrain", sep="")), FUN=function(x){any(class(x) == "error")}))
        disTestID <- which(sapply(get(paste(nam, "DistTest", sep="")), FUN=function(x){any(class(x) == "error")}))
        assTrainID <- which(sapply(get(paste(nam, "EnseTrain", sep="")), FUN=function(x){any(class(x) == "error")}))
        assTestID <- which(sapply(get(paste(nam, "EnseTest", sep="")), FUN=function(x){any(class(x) == "error")}))
      
        substID <- append(substID, c(disTrainID, disTestID, assTrainID, assTestID))
      }
      substID <- unique(substID)
      substID <- setdiff(1:10, substID)
      substitutesList[[j - 1]] <- substID
    }
    return(substitutesList)
  }
  

getNeeds <- function(nCrash, nSubst)
  {  
    # Create an index with the substitutions to be used from each iteration
    nNeeds <- list(0, 0, 0, 0)
    if(nCrash <= nSubst[[1]])
    {
      nNeeds[[1]] <- nCrash
    }else{
      nNeeds[[1]] <- nSubst[[1]]
      nCrash2 <- nCrash - nSubst[[1]]
      if(nCrash2 <= nSubst[[2]])
      {
        nNeeds[[2]] <- nCrash2
      }else{
        nNeeds[[2]] <- nSubst[[2]]
        nCrash3 <- nCrash2 - nSubst[[2]]
        if(nCrash3 <= nSubst[[3]])
        {
          nNeeds[[3]] <- nCrash3
        }else{
          nNeeds[[3]] <- nSubst[[3]]
          nCrash4 <- nCrash3 - nSubst[[2]]
          if(nCrash4 <= nSubst[[4]])
          {
            nNeeds[[4]] <- nCrash4
          }else{
            nNeeds[[4]] <- nSubst[[4]]
          }
        }
      }
    }
    return(nNeeds)
  }

##########################
## RUN THE SUBSTITUTION ##
##########################

## Define some parameters: Directories and object names
sdmNames <- c("glm","gam","scarts","smars","snn")
clmNames <- c("cqo","cao","carts","mars","mnn")
modNames <- c(sdmNames, clmNames)

objDTrainNames <- c("glmDistTrain", "cqoDistTrain", "gamDistTrain", "caoDistTrain", "scartsDistTrain", "cartsDistTrain", "smarsDistTrain", "marsDistTrain", "snnDistTrain", "mnnDistTrain")
objDTestNames <- c("glmDistTest", "cqoDistTest", "gamDistTest", "caoDistTest", "scartsDistTest", "cartsDistTest", "smarsDistTest", "marsDistTest", "snnDistTest", "mnnDistTest")
objETrainNames <- c("glmEnseTrain", "cqoEnseTrain", "gamEnseTrain", "caoEnseTrain", "scartsEnseTrain", "cartsEnseTrain", "smarsEnseTrain", "marsEnseTrain", "snnEnseTrain", "mnnEnseTrain")
objETestNames <- c("glmEnseTest", "cqoEnseTest", "gamEnseTest", "caoEnseTest", "scartsEnseTest", "cartsEnseTest", "smarsEnseTest", "marsEnseTest", "snnEnseTest", "mnnEnseTest")

objNames <- c(objDTrainNames, objDTestNames, objETrainNames, objETestNames)

evalDisDir <- c("Results/RObjects/Evaluation/Distribution/Replicate")
evalAssDir <- c("Results/RObjects/Evaluation/Assemblage/Replicate")


### ORIGINAL MODELS
# FOR loops for each fitting period 
for(i in periods_fit)
{
  print(paste("Period fit:", i, "BP", sep=""))

  # Get an index of crashes in each model and potential substitutions
  crashID <- getCrashes(paste(i, "BP", sep=""), modNames, evalDisDir, evalAssDir)
  substID <- getSubstitutes(paste(i, "BP", sep=""), modNames, evalDisDir, evalAssDir)
             
  # Calculate the number of crashes and substitutions
  numCrash <- length(crashID)
  numSubst <- lapply(substID, length)
  
  # If there are 10 crashes and no substitution nothing is done
  if(numCrash == 10 & sum(unlist(numSubst)) == 0){next}

  numNeeds <- getNeeds(numCrash, numSubst)
  
  # Load the data from the first iteration and assign new name to the objects
  load(paste(evalDisDir, "1/", i, "BP.RData", sep=""))
  load(paste(evalAssDir, "1/", i, "BP.RData", sep=""))
  for(l in 1:length(objNames))
  {
    assign(paste(objNames[l], ".Merged", sep=""), get(objNames[l]))  
  }
  
  # If there are "0" crashes don't need to do anything. Go to the end and save the files
  # Otherwise, go through the whole section to replace the crashes with the substitutions
  if(numCrash != 0 & sum(unlist(numSubst)) != 0)
  {
    
    # Load second iteration, check if there are substitution from this iteration and 
    # extract their results
    if(numNeeds[[1]] != 0)
    {
      load(paste(evalDisDir, "2/", i, "BP.RData", sep=""))
      load(paste(evalAssDir, "2/", i, "BP.RData", sep=""))
      index <- substID[[1]][1:numNeeds[[1]]]
      for(l in 1:length(objNames))
      {
        assign(paste(objNames[l], ".Subst", sep=""), get(objNames[l])[index])
      }
    }

    # Load third iteration, check if there are substitution from this iteration, 
    # extract their results, and if necesary append them to those from the former iteration
    if(numNeeds[[2]] != 0)
    {
      load(paste(evalDisDir, "3/", i, "BP.RData", sep=""))
      load(paste(evalAssDir, "3/", i, "BP.RData", sep=""))
      index <- substID[[2]][1:numNeeds[[2]]]
      for(l in 1:length(objNames))
      {
        if(numNeeds[[1]] == 0)
        {
          assign(paste(objNames[l], ".Subst", sep=""), get(objNames[l])[index])
        }else{
          tmp <- append(get(paste(objNames[l], ".Subst", sep="")), get(objNames[l])[index])
          assign(paste(objNames[l], ".Subst", sep=""), tmp)
        }
      }
    }
      
    # Load fourth iteration, check if there are substitution from this iteration, 
    # extract there results, and if necesary append them to those from the former iterations
    if(numNeeds[[3]] != 0)
    {
      load(paste(evalDisDir, "4/", i, "BP.RData", sep=""))
      load(paste(evalAssDir, "4/", i, "BP.RData", sep=""))
      index <- substID[[3]][1:numNeeds[[3]]]
      for(l in 1:length(objNames))
      {
        if(numNeeds[[1]] == 0 & numNeeds[[2]] == 0)
        {
          assign(paste(objNames[l], ".Subst", sep=""), get(objNames[l])[index])
        }else{
          tmp <- append(get(paste(objNames[l], ".Subst", sep="")), get(objNames[l])[index])
          assign(paste(objNames[l], ".Subst", sep=""), tmp)
        }
      }
    }
    
    # Load fourth iteration, check if there are substitution from this iteration, 
    # extract there results, and if necesary append them to those from the former iterations
    if(numNeeds[[4]] != 0)
    {
      load(paste(evalDisDir, "5/", i, "BP.RData", sep=""))
      load(paste(evalAssDir, "5/", i, "BP.RData", sep=""))
      index <- substID[[4]][1:numNeeds[[4]]]
      for(l in 1:length(objNames))
      {
        if(numNeeds[[1]] == 0 & numNeeds[[2]] == 0 & numNeeds[[3]] == 0)
        {
          assign(paste(objNames[l], ".Subst", sep=""), get(objNames[l])[index])
        }else{
          tmp <- append(get(paste(objNames[l], ".Subst", sep="")), get(objNames[l])[index])
          assign(paste(objNames[l], ".Subst", sep=""), tmp)
        }
      }
    }

    # Get an index in case there are more crashes than substitutions. Then, replace crashes
    # with substitutions
    index <- crashID[1:length(glmDistTrain.Subst)] 
    for(l in 1:length(objNames))
    {
      assign("tmpObj.Subst", get(paste(objNames[l], ".Subst", sep="")))
      assign("tmpObj.Merged", get(paste(objNames[l], ".Merged", sep="")))
      
      tmpObj.Merged[index] <- tmpObj.Subst
      
      assign(paste(objNames[l], ".Merged", sep=""), tmpObj.Merged)
      
      rm(tmpObj.Merged, tmpObj.Subst)
    }
  }

  # Save assemblages results to the HDD 
  direct <- "Results/RObjects/Evaluation/Assemblage/Merged"
  dir.create(direct, recursive=T) 
  save(glmEnseTrain.Merged, glmEnseTest.Merged, cqoEnseTrain.Merged, cqoEnseTest.Merged, gamEnseTrain.Merged, gamEnseTest.Merged, caoEnseTrain.Merged, caoEnseTest.Merged, scartsEnseTrain.Merged, scartsEnseTest.Merged, cartsEnseTrain.Merged, cartsEnseTest.Merged, smarsEnseTrain.Merged, smarsEnseTest.Merged, marsEnseTrain.Merged, marsEnseTest.Merged, snnEnseTrain.Merged, snnEnseTest.Merged, mnnEnseTrain.Merged, mnnEnseTest.Merged, file=paste(direct, "/", i, "BP.RData", sep=""))

  # Save distribution results to the HDD
  direct <- "Results/RObjects/Evaluation/Distribution/Merged"
  dir.create(direct, recursive=T) 
  save(glmDistTrain.Merged, glmDistTest.Merged, cqoDistTrain.Merged, cqoDistTest.Merged, gamDistTrain.Merged, gamDistTest.Merged, caoDistTrain.Merged, caoDistTest.Merged, scartsDistTrain.Merged, scartsDistTest.Merged, cartsDistTrain.Merged, cartsDistTest.Merged, smarsDistTrain.Merged, smarsDistTest.Merged, marsDistTrain.Merged, marsDistTest.Merged, snnDistTrain.Merged, snnDistTest.Merged, mnnDistTrain.Merged, mnnDistTest.Merged, file=paste(direct, "/", i, "BP.RData", sep=""))
  
  # Remove temporary objects
  rm(direct, numNeeds, numCrash, numSubst, crashID, substID, l) 
  rm(list=ls(pattern="DistTest"))
  rm(list=ls(pattern="DistTrain"))
  rm(list=ls(pattern="EnseTest"))
  rm(list=ls(pattern="EnseTrain"))
}

rm(i)


################################################################################
#
# FORMAT RESULTS
#
################################################################################

# List to store results
resd.All <- list()
rese.All <- list()

for(i in c(1:length(periods_fit))){

  print(i)

  j=index_fit[[i]]
  k=periods_fit[[i]]

  
  load(paste("Results/RObjects/Evaluation/Assemblage/Merged/", periods_fit[[i]], "BP.RData", sep=""), .GlobalEnv)
  load(paste("Results/RObjects/Evaluation/Distribution/Merged/", periods_fit[[i]], "BP.RData", sep=""), .GlobalEnv)

  ##
  ## Change errors in the evaluation with -9999 values. That avoid an error when melting data frames with all NAs
  ##
  for(v in 1:length(objDTrainNames)){
    assign("tmpName", get(paste(objDTrainNames[v], ".Merged", sep=""))) 
    ind <- which(unlist(lapply(lapply(tmpName, class), function(x){any(x == "error")})))
    if(length(ind) > 0){
       for(w in ind){
         tmp <- as.data.frame(matrix(-9999, nrow=19, ncol=8))
         colnames(tmp) <- c("taxon","prev","auc","trmin","tss","sens","spec","brier")
         tmp$taxon <- c("Ambrosia.type","Artemisia","Alnus","Corylus","Ostrya.Carpinus","Betula","Quercus","Carya","Juglans","Fraxinus","Abies","Larix","Picea","Platanus","Populus","Salix","Acer","Ulmus","Pinus")
         tmpName[[w]] <- tmp
       }
    }
    assign(objDTrainNames[v], tmpName) 
  }
  for(v in 1:length(objDTestNames)){
    assign("tmpName", get(paste(objDTestNames[v], ".Merged", sep=""))) 
    ind <- which(unlist(lapply(lapply(tmpName, class), function(x){any(x == "error")})))
    if(length(ind) > 0){
       for(w in ind){
         tmp <- as.data.frame(matrix(-9999, nrow=19, ncol=8))
         colnames(tmp) <- c("taxon","prev","auc","trmin","tss","sens","spec","brier")
         tmp$taxon <- c("Ambrosia.type","Artemisia","Alnus","Corylus","Ostrya.Carpinus","Betula","Quercus","Carya","Juglans","Fraxinus","Abies","Larix","Picea","Platanus","Populus","Salix","Acer","Ulmus","Pinus")
         
         tmpName[[w]] <- rep(list(tmp), length(periods_prj))
         names(tmpName[[w]]) <- periods_prj
       }
    }
    assign(objDTestNames[v], tmpName) 
  }
  for(v in 1:length(objETrainNames)){
    assign("tmpName", get(paste(objETrainNames[v], ".Merged", sep=""))) 
    ind <- which(unlist(lapply(lapply(tmpName, class), function(x){any(x == "error")})))
    if(length(ind) > 0){
      for(w in ind){
        tmpName[[w]] <- rep(-9999, 12)
        names(tmpName[[w]]) <- c("jacc.MEAN.prob","jacc.SD.prob","jacc.CAL.prob","jacc.DIS.prob","spRich.RSQ.prob","spRich.COR.prob","jacc.MEAN.thre","jacc.SD.thre","jacc.CAL.thre","jacc.DIS.thre","spRich.RSQ.thre","spRich.COR.thre")
      }
    }
    assign(objETrainNames[v], tmpName) 
  }
  for(v in 1:length(objETestNames)){
    assign("tmpName", get(paste(objETestNames[v], ".Merged", sep=""))) 
    ind <- which(unlist(lapply(lapply(tmpName, class), function(x){any(x == "error")})))
    if(length(ind) > 0){
      for(w in ind){
        tmpName[[w]] <- matrix(-9999, nrow=12, ncol=12, dimnames=list(c("jacc.MEAN.prob","jacc.SD.prob","jacc.CAL.prob","jacc.DIS.prob","spRich.RSQ.prob","spRich.COR.prob","jacc.MEAN.thre","jacc.SD.thre","jacc.CAL.thre","jacc.DIS.thre","spRich.RSQ.thre","spRich.COR.thre"), periods_prj))
      }
    }
    assign(objETestNames[v], tmpName) 
  }
  
  rm(v, w, tmp, tmpName)
  
  library(reshape2)
  
  # Results for TRAIN in 0BP
  # Distribution
  glmdt <- melt(glmDistTrain)
  cqodt <- melt(cqoDistTrain)
  gamdt <- melt(gamDistTrain)
  caodt <- melt(caoDistTrain)
  smarsdt <-melt(smarsDistTrain)
  marsdt <- melt(marsDistTrain)
  scartsdt <-melt(scartsDistTrain)
  cartsdt <- melt(cartsDistTrain)
  snndt <-melt(snnDistTrain)
  mnndt <- melt(mnnDistTrain)
  colnames(glmdt) <- colnames(cqodt) <- colnames(gamdt) <- colnames(caodt) <- colnames(smarsdt) <- colnames(marsdt) <- colnames(scartsdt) <- colnames(cartsdt) <- colnames(snndt) <- colnames(mnndt) <- c("taxon","variable","value","iteration")
  
  resDt <- rbind(glmdt, cqodt, gamdt, caodt, smarsdt, marsdt, scartsdt, cartsdt, snndt, mnndt)
  resDt <- cbind(resDt, model=c(rep("glm", nrow(glmdt)), rep("cqo", nrow(cqodt)), rep("gam", nrow(gamdt)), rep("cao", nrow(caodt)), rep("smars", nrow(smarsdt)), rep("mars", nrow(marsdt)), rep("scarts", nrow(scartsdt)), rep("carts", nrow(cartsdt)), rep("snn", nrow(snndt)), rep("mnn", nrow(mnndt))))
  resDt$periodFit <- as.character(k)
  resDt$periodPrj <- as.character(k)
  resDt$test <- "train"
  rm(glmdt, gamdt, cqodt, caodt, marsdt, smarsdt, cartsdt, scartsdt, mnndt, snndt)
  
  # Assemblage
  glmet <- melt(glmEnseTrain)
  cqoet <- melt(cqoEnseTrain)
  gamet <- melt(gamEnseTrain)
  caoet <- melt(caoEnseTrain)
  smarset <- melt(smarsEnseTrain)
  marset <- melt(marsEnseTrain)
  scartset <-melt(scartsEnseTrain)
  cartset <- melt(cartsEnseTrain)
  snnet <- melt(snnEnseTrain)
  mnnet <- melt(mnnEnseTrain)
  colnames(glmet) <- colnames(cqoet) <- colnames(gamet) <- colnames(caoet) <- colnames(smarset) <- colnames(marset) <- colnames(scartset) <- colnames(cartset) <- colnames(snnet) <- colnames(mnnet) <- c("value","iteration")
  
  resEt <- rbind(glmet, cqoet, gamet, caoet, smarset, marset, scartset, cartset, snnet, mnnet)
  resEt <- cbind(variable=rep(names(glmEnseTrain[[1]]), nrow(resEt)/length(names(glmEnseTrain[[1]]))), resEt)
  resEt <- cbind(resEt, model=c(rep("glm", nrow(glmet)), rep("cqo", nrow(cqoet)), rep("gam", nrow(gamet)), rep("cao", nrow(caoet)), rep("smars", nrow(smarset)), rep("mars", nrow(marset)), rep("scarts", nrow(scartset)), rep("carts", nrow(cartset)), rep("snn", nrow(snnet)), rep("mnn", nrow(mnnet))))
  resEt$periodFit <- as.character(k)
  resEt$periodPrj <- as.character(k)
  resEt$test <- "train"
  rm(glmet, gamet, cqoet, caoet, marset, smarset, cartset, scartset, mnnet, snnet)
  
  # Results for TEST in all periods
  # Distribution
  glmd <- melt(glmDistTest)
  cqod <- melt(cqoDistTest)
  gamd <- melt(gamDistTest)
  caod <- melt(caoDistTest)
  smarsd <-melt(smarsDistTest)
  marsd <- melt(marsDistTest)
  scartsd <-melt(scartsDistTest)
  cartsd <- melt(cartsDistTest)
  snnd <-melt(snnDistTest)
  mnnd <- melt(mnnDistTest)
  colnames(glmd) <- colnames(cqod) <- colnames(gamd) <- colnames(caod) <- colnames(smarsd) <- colnames(marsd) <- colnames(scartsd) <- colnames(cartsd) <- colnames(mnnd) <- colnames(snnd) <- c("taxon","variable","value","periodPrj","iteration")
  
  resD <- rbind(glmd, cqod, gamd, caod, smarsd, marsd, scartsd, cartsd, snnd, mnnd)
  resD <- cbind(resD, model=c(rep("glm", nrow(glmd)), rep("cqo", nrow(cqod)), rep("gam", nrow(gamd)), rep("cao", nrow(caod)), rep("smars", nrow(smarsd)), rep("mars", nrow(marsd)), rep("scarts", nrow(scartsd)), rep("carts", nrow(cartsd)), rep("snn", nrow(snnd)), rep("mnn", nrow(mnnd))))
  resD$periodFit <- as.character(k)
  resD$test <- "test"
  rm(glmd, gamd, cqod, caod, marsd, smarsd, cartsd, scartsd, mnnd, snnd)
  
  # Assemblage
  glme <- melt(glmEnseTest)
  game <- melt(gamEnseTest)
  cqoe <- melt(cqoEnseTest)
  caoe <- melt(caoEnseTest)
  marse <- melt(marsEnseTest)
  smarse <- melt(smarsEnseTest)
  cartse <- melt(cartsEnseTest)
  scartse <- melt(scartsEnseTest)
  mnne <- melt(mnnEnseTest)
  snne <- melt(snnEnseTest)
  colnames(glme) <- colnames(cqoe) <- colnames(game) <- colnames(caoe) <- colnames(smarse) <- colnames(marse) <- colnames(scartse) <- colnames(cartse) <- colnames(snne) <- colnames(mnne) <- c("variable","periodPrj","value","iteration")
  
  resE <- rbind(glme, cqoe, game, caoe, smarse, marse, scartse, cartse, snne, mnne)
  resE <- cbind(resE, model=c(rep("glm", nrow(glme)), rep("cqo", nrow(cqoe)), rep("gam", nrow(game)), rep("cao", nrow(caoe)), rep("smars", nrow(smarse)), rep("mars", nrow(marse)), rep("scarts", nrow(scartse)), rep("carts", nrow(cartse)), rep("snn", nrow(snne)), rep("mnn", nrow(mnne))))
  resE$test <- "test"
  resE$periodFit <- as.character(k)
  rm(glme, game, cqoe, caoe, marse, smarse, cartse, scartse, mnne, snne)
  
  # Combining the TRAINING and TESTING datasets
  resd.All[[i]] <- rbind(resDt, resD)
  rese.All[[i]] <- rbind(resEt, resE)

  rm(list=ls(pattern=".Merged"))
  rm(list=ls(pattern="DistTest"))
  rm(list=ls(pattern="DistTrain"))
  rm(list=ls(pattern="EnseTest"))
  rm(list=ls(pattern="EnseTrain"))
}


rm(i, j, k, resDt, resEt, resD, resE)
rm(objNames, objDTestNames, objDTrainNames, objETestNames, objETrainNames)
rm(caoFit, cartsFit, cartsFit.sp, cqoFit, gamFit, glmFit, marsFit, marsFit.sp, mnnFit, snnFit, seed)
 
# Write data to file
save(resd.All, file="Results/RObjects/DisResults.RData")
save(rese.All, file="Results/RObjects/AssResults.RData")
