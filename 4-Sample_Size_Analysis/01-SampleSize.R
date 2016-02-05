################################################################################
#
# LOAD DATA, PARAMETERS, AND FUNCTIONS
#
################################################################################

# Load the project library with functions to load and arrange data
library(paleoCLMs)

# Load wrapping functions
source("00-WrappingFunctions.R")

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
periods_fit <- c("1", "9")
periods_prj <- c("0", "0.5", "3", "6", "7.5", "9.5", "12.5", "14.5", "15", "15.5-17.5", "15.5-21", "18-21")

index_fit <- as.list(which(periods_red %in% periods_fit))
index_prj <- as.list(which(periods_red %in% periods_prj))
index_prj <- append(index_prj, list(c(32:36), c(32:43), c(37:43)))


###############################################################################
##
## FIT THE MODELS
##
###############################################################################

require(foreach)
require(doParallel)
registerDoParallel(cores=5)

seed <- c(1, 11, 21, 31)

# Set the samples sizes to be tested
sampleSizes <- c(200, 150, 100, 50)
  

## Some of the cqo and cao models crash at certain time periods. I run the cqo and cao models twice (Replicate1 and Replicate2) by hand. Then I replaced those that crashed in the first run with some from the second run. See the code below for merging the models
for(i in 1:length(periods_fit)){

  for(seed.i in 1:4){

    # Get the index and the value itself for periods_fit
    j=index_fit[[i]]
    k=periods_fit[[i]]
    
    a <- Sys.time()
    print(paste("Seed:", seed[seed.i], sep=""))
    print(paste("  Fitting period: ", k, "BP", sep=""))
    
    # Split training and testing datasets. THAT SHOULD BE BEFORE POOLING DATA
    trainSets <- lapply(pollen, trainSplit, numIter, 0, seed[seed.i])
    trainSets[j] <- lapply(pollen[j], trainSplit, numIter, 0.7, seed[seed.i])
    
    # Pooling data together
    # Binning 15.5 to 17.5
    binIndex <- which(periods %in% seq(15500, 17500, by=500))
    pbin1 <- poolingData(binIndex, pollen)
    cbin1 <- poolingData(binIndex, clim)
    tbin1 <- poolingData(binIndex, trainSets)
    
    # Binning 15.5 to 21
    binIndex <- which(periods %in% seq(15500, 21000, by=500))
    pbin2 <- poolingData(binIndex, pollen)
    cbin2 <- poolingData(binIndex, clim)
    tbin2 <- poolingData(binIndex, trainSets)
    
    # Binning 18 to 21
    binIndex <- which(periods %in% seq(18000, 21000, by=500))
    pbin3 <- poolingData(binIndex, pollen)
    cbin3 <- poolingData(binIndex, clim)
    tbin3 <- poolingData(binIndex, trainSets)
    
    # Appending the pooled dataframes to the data lists
    pollenPool <- append(pollen, list(pbin1, pbin2, pbin3))
    climPool <- append(clim, list(cbin1, cbin2, cbin3))
    trainSetsPool <- append(trainSets, list(tbin1, tbin2, tbin3))
    
    # Remove temporal objects
    rm(binIndex, cbin1, cbin2, cbin3, pbin1, pbin2, pbin3, tbin1, tbin2, tbin3)
    
    # Change names
    names(pollenPool) <- names(climPool) <- names(trainSetsPool) <- c(periods_red, "15.5-17.5", "15.5-21", "18-21")
    
    # Selecting and splitting data for fitting and projecting time periods
    pollenFit <- pollenPool[which(names(pollenPool) %in% periods_fit)]
    climFit <- climPool[which(names(climPool) %in% periods_fit)]
    trainSetsFit <- trainSetsPool[which(names(trainSetsPool) %in% periods_fit)]
    
    pollenPrj <- pollenPool[which(names(pollenPool) %in% periods_prj)]
    climPrj <- climPool[which(names(climPool) %in% periods_prj)]
    trainSetsPrj <- trainSetsPool[which(names(trainSetsPool) %in% periods_prj)]
    
    # Load tunned parameters for each taxon and model
    tunedParam <- read.csv("parameters.csv", header=TRUE, row.names=1)
    
    # Fitting the original models
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
  
    # Create directory and save models
    directMod <- paste("Results/RObjects/Models/Replicate", seed.i, "/", sep="")
    dir.create(directMod, recursive=T) 
    save(glmList, cqoList, gamList, caoList, scartsList, cartsList, smarsList, marsList, snnList, mnnList, file=paste(directMod, periods_fit[[i]], "BP.RData", sep=""))
  
    ## EVALUATE THE MODELS
    # GLM - SDM
    glmDistTrain <- evalDistTrain(glmList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="response", multiSp=F)
    glmDistTest <- evalDistTest(glmList, pollenPrj, pollenFit[[i]], climPrj, climFit[[i]], trainSetsPrj, trainSetsFit[[i]], respType="response", multiSp=F)
    glmEnseTrain <- evalEnseTrain(glmList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="response", glmDistTrain, multiSp=F)
    glmEnseTest <- evalEnseTest(glmList, pollenPrj, pollenFit[[i]], climPrj, trainSetsPrj, respType="response", glmDistTrain, multiSp=F)
     
    # CQO - CLM
    cqoDistTrain<- evalDistTrain(cqoList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="response")
    cqoDistTest <- evalDistTest(cqoList, pollenPrj, pollenFit[[i]], climPrj, climFit[[i]], trainSetsPrj, trainSetsFit[[i]], respType="response")
    cqoEnseTrain <- evalEnseTrain(cqoList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="response", cqoDistTrain)
    cqoEnseTest <- evalEnseTest(cqoList, pollenPrj, pollenFit[[i]], climPrj, trainSetsPrj, respType="response", cqoDistTrain)
    
    # GAM - SDM
    gamDistTrain <- evalDistTrain(gamList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="response", multiSp=F)
    gamDistTest <- evalDistTest(gamList, pollenPrj, pollenFit[[i]], climPrj, climFit[[i]], trainSetsPrj, trainSetsFit[[i]], respType="response", multiSp=F)
    gamEnseTrain <- evalEnseTrain(gamList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="response", gamDistTrain, multiSp=F)
    gamEnseTest <- evalEnseTest(gamList, pollenPrj, pollenFit[[i]], climPrj, trainSetsPrj, respType="response", gamDistTrain, multiSp=F)
     
    # CAO - CLM
    caoDistTrain<- evalDistTrain(caoList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="response")
    caoDistTest <- evalDistTest(caoList, pollenPrj, pollenFit[[i]], climPrj, climFit[[i]], trainSetsPrj, trainSetsFit[[i]], respType="response")
    caoEnseTrain <- evalEnseTrain(caoList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="response", caoDistTrain)
    caoEnseTest <- evalEnseTest(caoList, pollenPrj, pollenFit[[i]], climPrj, trainSetsPrj, respType="response", caoDistTrain)
  
    # MARS - SDM
    smarsDistTrain <- evalDistTrain(smarsList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="response", multiSp=F)
    smarsDistTest <- evalDistTest(smarsList, pollenPrj, pollenFit[[i]], climPrj, climFit[[i]], trainSetsPrj, trainSetsFit[[i]], respType="response", multiSp=F)
    smarsEnseTrain <- evalEnseTrain(smarsList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="response", smarsDistTrain, multiSp=F)
    smarsEnseTest <- evalEnseTest(smarsList, pollenPrj, pollenFit[[i]], climPrj, trainSetsPrj, respType="response", smarsDistTrain, multiSp=F)
    
    # MARS - CLM
    marsDistTrain <- evalDistTrain(marsList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="response")
    marsDistTest <- evalDistTest(marsList, pollenPrj, pollenFit[[i]], climPrj, climFit[[i]], trainSetsPrj, trainSetsFit[[i]], respType="response")
    marsEnseTrain <- evalEnseTrain(marsList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="response", marsDistTrain)
    marsEnseTest <- evalEnseTest(marsList, pollenPrj, pollenFit[[i]], climPrj, trainSetsPrj, respType="response", marsDistTrain)
     
    # CARTs - SDM
    scartsDistTrain <- evalDistTrain(scartsList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="matrix", multiSp=F)
    scartsDistTest <- evalDistTest(scartsList, pollenPrj, pollenFit[[i]], climPrj, climFit[[i]], trainSetsPrj, trainSetsFit[[i]], respType="matrix", multiSp=F)
    scartsEnseTrain <- evalEnseTrain(scartsList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="matrix", scartsDistTrain, multiSp=F)
    scartsEnseTest <- evalEnseTest(scartsList, pollenPrj, pollenFit[[i]], climPrj, trainSetsPrj, respType="matrix", scartsDistTrain, multiSp=F)
    
    # CARTs - CLM
    cartsDistTrain <- evalDistTrain(cartsList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="matrix")
    cartsDistTest <- evalDistTest(cartsList, pollenPrj, pollenFit[[i]], climPrj, climFit[[i]], trainSetsPrj, trainSetsFit[[i]], respType="matrix")
    cartsEnseTrain <- evalEnseTrain(cartsList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="matrix", cartsDistTrain)
    cartsEnseTest <- evalEnseTest(cartsList, pollenPrj, pollenFit[[i]], climPrj, trainSetsPrj, respType="matrix", cartsDistTrain)
    
    # NEURAL NETWORKS - SDM 
    library (nnet)
    snnDistTrain <- evalDistTrain(snnList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="raw", multiSp=F)
    snnDistTest <- evalDistTest(snnList, pollenPrj, pollenFit[[i]], climPrj, climFit[[i]], trainSetsPrj, trainSetsFit[[i]], respType="raw", multiSp=F)
    snnEnseTrain <- evalEnseTrain(snnList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="raw", snnDistTrain, multiSp=F)
    snnEnseTest <- evalEnseTest(snnList, pollenPrj, pollenFit[[i]], climPrj, trainSetsPrj, respType="raw", snnDistTrain, multiSp=F)
     
    # NEURAL NETWORKS - CLM 
    mnnDistTrain <- evalDistTrain(mnnList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="raw")
    mnnDistTest <- evalDistTest(mnnList, pollenPrj, pollenFit[[i]], climPrj, climFit[[i]], trainSetsPrj, trainSetsFit[[i]], respType="raw")
    mnnEnseTrain <- evalEnseTrain(mnnList, pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], respType="raw", mnnDistTrain)
    mnnEnseTest <- evalEnseTest(mnnList, pollenPrj, pollenFit[[i]], climPrj, trainSetsPrj, respType="raw", mnnDistTrain)
  
    # Save the RObjects in their respective folders
    # Distribution
    directEvalDis <- paste("Results/RObjects/Evaluation/Distribution/Replicate", seed.i, "/", sep="")
    dir.create(directEvalDis, recursive=T) 
    save(glmDistTrain, glmDistTest, cqoDistTrain, cqoDistTest, gamDistTrain, gamDistTest, caoDistTrain, caoDistTest, scartsDistTrain, scartsDistTest, cartsDistTrain, cartsDistTest, smarsDistTrain, smarsDistTest, marsDistTrain, marsDistTest, snnDistTrain, snnDistTest, mnnDistTrain, mnnDistTest, file=paste(directEvalDis, periods_fit[[i]], "BP.RData", sep=""))
   
    # Assemblage
    directEvalAss <- paste("Results/RObjects/Evaluation/Assemblage/Replicate", seed.i, "/", sep="")
    dir.create(directEvalAss, recursive=T) 
    save(glmEnseTrain, glmEnseTest, cqoEnseTrain, cqoEnseTest, gamEnseTrain, gamEnseTest, caoEnseTrain, caoEnseTest, scartsEnseTrain, scartsEnseTest, cartsEnseTrain, cartsEnseTest, smarsEnseTrain, smarsEnseTest, marsEnseTrain, marsEnseTest, snnEnseTrain, snnEnseTest, mnnEnseTrain, mnnEnseTest, file=paste(directEvalAss, periods_fit[[i]], "BP.RData", sep=""))
  
    # Remove models to avoid problems
    rm(glmList, cqoList, gamList, caoList, scartsList, cartsList, smarsList, marsList, snnList, mnnList)
    rm(glmDistTrain, cqoDistTrain, gamDistTrain, caoDistTrain, scartsDistTrain, cartsDistTrain, smarsDistTrain, marsDistTrain, snnDistTrain, mnnDistTrain)
    rm(glmDistTest, cqoDistTest, gamDistTest, caoDistTest, scartsDistTest, cartsDistTest, smarsDistTest, marsDistTest, snnDistTest, mnnDistTest)
    rm(glmEnseTrain, cqoEnseTrain, gamEnseTrain, caoEnseTrain, scartsEnseTrain, cartsEnseTrain, smarsEnseTrain, marsEnseTrain, snnEnseTrain, mnnEnseTrain)
    rm(glmEnseTest, cqoEnseTest, gamEnseTest, caoEnseTest, scartsEnseTest, cartsEnseTest, smarsEnseTest, marsEnseTest, snnEnseTest, mnnEnseTest)
  
  
    #################################
    # FITTING THE SUBSETTING MODELS
    
    # FOR loop for each iteration  
    for(o in 1:length(sampleSizes)){
  
      # Assigning sample size
      if(sampleSizes[o] > nrow(pollenFit[[i]])){
        next
      }else{
        sampleSize <- sampleSizes[o]
      }
  
      print(paste("    Sample size: ", sampleSize, sep=""))
  
      # FOR loop for ten iterations each sample size
      for(m in 1:10){
      
        print(paste("      Random subset: ", m, sep=""))
    
        # Subsetting all the dataframes
        set.seed(m)
        randomSites <- sample(1:nrow(pollenFit[[i]]), sampleSize)
  
        subPollenFit <- pollenFit[[i]][randomSites,]
        subClimFit <- climFit[[i]][randomSites,]
        subTrainFit <- trainSetsFit[[i]][randomSites,]
      
        # Fitting the models
        glmList <- glmFit(subPollenFit, subClimFit, subTrainFit, numIter)
        cqoList <- cqoFit(subPollenFit, subClimFit, subTrainFit, numIter)
        gamList <- gamFit(subPollenFit, subClimFit, subTrainFit, numIter)
        caoList <- caoFit(subPollenFit, subClimFit, subTrainFit, numIter)
        smarsList <- marsFit.sp(subPollenFit, subClimFit, subTrainFit, numIter, tunedParam)
        marsList <- marsFit(subPollenFit, subClimFit, subTrainFit, numIter, tunedParam)
        scartsList <- cartsFit.sp(subPollenFit, subClimFit, subTrainFit, numIter)
        cartsList <- cartsFit(subPollenFit, subClimFit, subTrainFit, numIter)
        snnList <- snnFit(data.frame(subPollenFit), subClimFit, subTrainFit, numIter, tunedParam)
        mnnList <- mnnFit(subPollenFit, subClimFit, subTrainFit, numIter, tunedParam)
      
        save(glmList, cqoList, gamList, caoList, scartsList, cartsList, smarsList, marsList, snnList, mnnList, file=paste(directMod, periods_fit[[i]], "BP-",sampleSizes[o],"-",m,".RData", sep=""))
   
        
        ## EVALUATE THE MODELS
        # Only replica 1 have all the models, replicas 2, and 3 have only CQO and CAO
        # GLM - SDM
        glmDistTrain <- evalDistTrain(glmList, subPollenFit, subClimFit, subTrainFit, respType="response", multiSp=F)
        glmDistTest <- evalDistTest(glmList, pollenPrj, subPollenFit, climPrj, subClimFit, trainSetsPrj, subTrainFit, respType="response", multiSp=F)
        glmEnseTrain <- evalEnseTrain(glmList, subPollenFit, subClimFit, subTrainFit, respType="response", glmDistTrain, multiSp=F)
        glmEnseTest <- evalEnseTest(glmList, pollenPrj, subPollenFit, climPrj, trainSetsPrj, respType="response", glmDistTrain, multiSp=F)
         
        # CQO - CLM
        cqoDistTrain<- evalDistTrain(cqoList, subPollenFit, subClimFit, subTrainFit, respType="response")
        cqoDistTest <- evalDistTest(cqoList, pollenPrj, subPollenFit, climPrj, subClimFit, trainSetsPrj, subTrainFit, respType="response")
        cqoEnseTrain <- evalEnseTrain(cqoList, subPollenFit, subClimFit, subTrainFit, respType="response", cqoDistTrain)
        cqoEnseTest <- evalEnseTest(cqoList, pollenPrj, subPollenFit, climPrj, trainSetsPrj, respType="response", cqoDistTrain)
        
        # GAM - SDM
        gamDistTrain <- evalDistTrain(gamList, subPollenFit, subClimFit, subTrainFit, respType="response", multiSp=F)
        gamDistTest <- evalDistTest(gamList, pollenPrj, subPollenFit, climPrj, subClimFit, trainSetsPrj, subTrainFit, respType="response", multiSp=F)
        gamEnseTrain <- evalEnseTrain(gamList, subPollenFit, subClimFit, subTrainFit, respType="response", gamDistTrain, multiSp=F)
        gamEnseTest <- evalEnseTest(gamList, pollenPrj, subPollenFit, climPrj, trainSetsPrj, respType="response", gamDistTrain, multiSp=F)
  
        # CAO - CLM
        caoDistTrain<- evalDistTrain(caoList, subPollenFit, subClimFit, subTrainFit, respType="response")
        caoDistTest <- evalDistTest(caoList, pollenPrj, subPollenFit, climPrj, subClimFit, trainSetsPrj, subTrainFit, respType="response")
        caoEnseTrain <- evalEnseTrain(caoList, subPollenFit, subClimFit, subTrainFit, respType="response", caoDistTrain)
        caoEnseTest <- evalEnseTest(caoList, pollenPrj, subPollenFit, climPrj, trainSetsPrj, respType="response", caoDistTrain)
             
        # MARS - SDM
        smarsDistTrain <- evalDistTrain(smarsList, subPollenFit, subClimFit, subTrainFit, respType="response", multiSp=F)
        smarsDistTest <- evalDistTest(smarsList, pollenPrj, subPollenFit, climPrj, subClimFit, trainSetsPrj, subTrainFit, respType="response", multiSp=F)
        smarsEnseTrain <- evalEnseTrain(smarsList, subPollenFit, subClimFit, subTrainFit, respType="response", smarsDistTrain, multiSp=F)
        smarsEnseTest <- evalEnseTest(smarsList, pollenPrj, subPollenFit, climPrj, trainSetsPrj, respType="response", smarsDistTrain, multiSp=F)
        
        # MARS - CLM
        marsDistTrain <- evalDistTrain(marsList, subPollenFit, subClimFit, subTrainFit, respType="response")
        marsDistTest <- evalDistTest(marsList, pollenPrj, subPollenFit, climPrj, subClimFit, trainSetsPrj, subTrainFit, respType="response")
        marsEnseTrain <- evalEnseTrain(marsList, subPollenFit, subClimFit, subTrainFit, respType="response", marsDistTrain)
        marsEnseTest <- evalEnseTest(marsList, pollenPrj, subPollenFit, climPrj, trainSetsPrj, respType="response", marsDistTrain)
         
        # CARTs - SDM
        scartsDistTrain <- evalDistTrain(scartsList, subPollenFit, subClimFit, subTrainFit, respType="matrix", multiSp=F)
        scartsDistTest <- evalDistTest(scartsList, pollenPrj, subPollenFit, climPrj, subClimFit, trainSetsPrj, subTrainFit, respType="matrix", multiSp=F)
        scartsEnseTrain <- evalEnseTrain(scartsList, subPollenFit, subClimFit, subTrainFit, respType="matrix", scartsDistTrain, multiSp=F)
        scartsEnseTest <- evalEnseTest(scartsList, pollenPrj, subPollenFit, climPrj, trainSetsPrj, respType="matrix", scartsDistTrain, multiSp=F)
        
        # CARTs - CLM
        cartsDistTrain <- evalDistTrain(cartsList, subPollenFit, subClimFit, subTrainFit, respType="matrix")
        cartsDistTest <- evalDistTest(cartsList, pollenPrj, subPollenFit, climPrj, subClimFit, trainSetsPrj, subTrainFit, respType="matrix")
        cartsEnseTrain <- evalEnseTrain(cartsList, subPollenFit, subClimFit, subTrainFit, respType="matrix", cartsDistTrain)
        cartsEnseTest <- evalEnseTest(cartsList, pollenPrj, subPollenFit, climPrj, trainSetsPrj, respType="matrix", cartsDistTrain)
        
        # NEURAL NETWORKS - SDM 
        library (nnet)
        snnDistTrain <- evalDistTrain(snnList, subPollenFit, subClimFit, subTrainFit, respType="raw", multiSp=F)
        snnDistTest <- evalDistTest(snnList, pollenPrj, subPollenFit, climPrj, subClimFit, trainSetsPrj, subTrainFit, respType="raw", multiSp=F)
        snnEnseTrain <- evalEnseTrain(snnList, subPollenFit, subClimFit, subTrainFit, respType="raw", snnDistTrain, multiSp=F)
        snnEnseTest <- evalEnseTest(snnList, pollenPrj, subPollenFit, climPrj, trainSetsPrj, respType="raw", snnDistTrain, multiSp=F)
         
        # NEURAL NETWORKS - CLM 
        mnnDistTrain <- evalDistTrain(mnnList, subPollenFit, subClimFit, subTrainFit, respType="raw")
        mnnDistTest <- evalDistTest(mnnList, pollenPrj, subPollenFit, climPrj, subClimFit, trainSetsPrj, subTrainFit, respType="raw")
        mnnEnseTrain <- evalEnseTrain(mnnList, subPollenFit, subClimFit, subTrainFit, respType="raw", mnnDistTrain)
        mnnEnseTest <- evalEnseTest(mnnList, pollenPrj, subPollenFit, climPrj, trainSetsPrj, respType="raw", mnnDistTrain)
      
        # Save the RObjects in their respective folders
        # Distribution
        save(glmDistTrain, glmDistTest, cqoDistTrain, cqoDistTest, gamDistTrain, gamDistTest, caoDistTrain, caoDistTest, scartsDistTrain, scartsDistTest, cartsDistTrain, cartsDistTest, smarsDistTrain, smarsDistTest, marsDistTrain, marsDistTest, snnDistTrain, snnDistTest, mnnDistTrain, mnnDistTest, file=paste(directEvalDis, periods_fit[[i]], "BP-",sampleSizes[o],"-",m,".RData", sep=""))
  
        # Assemblage
        save(glmEnseTrain, glmEnseTest, cqoEnseTrain, cqoEnseTest, gamEnseTrain, gamEnseTest, caoEnseTrain, caoEnseTest, scartsEnseTrain, scartsEnseTest, cartsEnseTrain, cartsEnseTest, smarsEnseTrain, smarsEnseTest, marsEnseTrain, marsEnseTest, snnEnseTrain, snnEnseTest, mnnEnseTrain, mnnEnseTest, file=paste(directEvalAss, periods_fit[[i]], "BP-",sampleSizes[o],"-",m,".RData", sep=""))
  
  
        # Remove models to avoid problems
        rm(glmList, cqoList, gamList, caoList, scartsList, cartsList, smarsList, marsList, snnList, mnnList)
        rm(glmDistTrain, cqoDistTrain, gamDistTrain, caoDistTrain, scartsDistTrain, cartsDistTrain, smarsDistTrain, marsDistTrain, snnDistTrain, mnnDistTrain)
        rm(glmDistTest, cqoDistTest, gamDistTest, caoDistTest, scartsDistTest, cartsDistTest, smarsDistTest, marsDistTest, snnDistTest, mnnDistTest)
        rm(glmEnseTrain, cqoEnseTrain, gamEnseTrain, caoEnseTrain, scartsEnseTrain, cartsEnseTrain, smarsEnseTrain, marsEnseTrain, snnEnseTrain, mnnEnseTrain)
        rm(glmEnseTest, cqoEnseTest, gamEnseTest, caoEnseTest, scartsEnseTest, cartsEnseTest, smarsEnseTest, marsEnseTest, snnEnseTest, mnnEnseTest)
      }
    }
    
    b <- Sys.time()
    print(b-a)
  
    rm(j, k, a, o, m, trainSets, pollenPool, climPool, trainSetsPool, subPollenFit, subClimFit, subTrainFit, pollenFit, climFit, trainSetsFit, pollenPrj, climPrj, trainSetsPrj, tunedParam, b)
  }
}

rm(seed.i, i)             


rm(caoFit, cartsFit, cartsFit.sp, cqoFit, gamFit, glmFit, marsFit, marsFit.sp, mnnFit, snnFit)


################################################################################
#
# MERGE RESULTS FROM THE DIFFERENTE REPLICATES
#
################################################################################

## FUNCTIONS

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
    for(j in 2:4){
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



##########################
## RUN THE SUBSTITUTION

## Define some parameters: Directories and object names
sdmNames <- c("glm","gam","scarts","smars","snn")
clmNames <- c("cqo","cao","carts","mars","mnn")
modNames <- c(sdmNames, clmNames)

objNames <- c("glmDistTrain", "glmDistTest", "glmEnseTrain", "glmEnseTest"
            , "cqoDistTrain", "cqoDistTest", "cqoEnseTrain", "cqoEnseTest"
            , "gamDistTrain", "gamDistTest", "gamEnseTrain", "gamEnseTest"
            , "caoDistTrain", "caoDistTest", "caoEnseTrain", "caoEnseTest"
            , "scartsDistTrain", "scartsDistTest", "scartsEnseTrain", "scartsEnseTest"
            , "cartsDistTrain", "cartsDistTest", "cartsEnseTrain", "cartsEnseTest"
            , "smarsDistTrain", "smarsDistTest", "smarsEnseTrain", "smarsEnseTest"
            , "marsDistTrain", "marsDistTest", "marsEnseTrain", "marsEnseTest"
            , "snnDistTrain", "snnDistTest", "snnEnseTrain", "snnEnseTest"
            , "mnnDistTrain", "mnnDistTest", "mnnEnseTrain", "mnnEnseTest"
             )


evalDisDir <- c("Results/RObjects/Evaluation/Distribution/Replicate")
evalAssDir <- c("Results/RObjects/Evaluation/Assemblage/Replicate")


### ORIGINAL MODELS
# FOR loops for each fitting period, sample size, and iteration 
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

  # Create an index with the substitutions to be used from each iteration
  numNeeds <- list(0, 0, 0)
  if(numCrash <= numSubst[[1]])
  {
    numNeeds[[1]] <- numCrash
  }else{
    numNeeds[[1]] <- numSubst[[1]]
    numCrash2 <- numCrash - numSubst[[1]]
    if(numCrash2 <= numSubst[[2]])
    {
      numNeeds[[2]] <- numCrash2
    }else{
      numNeeds[[2]] <- numSubst[[2]]
      numCrash3 <- numCrash2 - numSubst[[2]]
      if(numCrash3 <= numSubst[[3]])
      {
        numNeeds[[3]] <- numCrash3
      }else{
        numNeeds[[3]] <- numSubst[[3]]
      }
      rm(numCrash3)
    }
    rm(numCrash2)  
  }
  
  # Load the data from the first iteration and assign new name to the objects
  load(paste(evalDisDir, "1/", i, "BP.RData", sep=""))
  load(paste(evalAssDir, "1/", i, "BP.RData", sep=""))
  for(l in 1:length(objNames))
  {
    assign(paste(objNames[l], ".Merged", sep=""), get(objNames[l]))  
  }
  
  # If there are "0" crashes don't need to do anything. Go to the end and save the files
  # Otherwise, go through the whole section to replace the crashes with the substitutions
  if(numCrash != 0)
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
    
    # Get an index in case there are more crashes than substitutions. Then, replace crashes
    # with substitutions
    index <- crashID[1:length(glmDistTrain.Subst)] 
    for(l in 1:length(objNames))
    {
      assign("tmpObj.Subst", get(paste(objNames[l], ".Subst", sep="")))
      assign("tmpObj.Merged", get(paste(objNames[l], ".Merged", sep="")))
      
      tmpObj.Merged[index] <- tmpObj.Subst
      
      assign(paste(objNames[l], ".Merged", sep=""), tmpObj.Merged)
      
      rm(tmpObj.Merged)
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
  rm(direct, tmpObj.Merged, tmpObj.Subst, index, numNeeds, numCrash, numSubst, crashID, substID, i, l) 
  rm(list=ls(pattern="DistTest"))
  rm(list=ls(pattern="DistTrain"))
  rm(list=ls(pattern="EnseTest"))
  rm(list=ls(pattern="EnseTrain"))
}

### SUBSET MODELS
# FOR loops for each fitting period, sample size, and iteration 
for(i in periods_fit)
{
  print(paste("Period fit:", i, "BP", sep=""))
  for(j in sampleSizes)
  {
    if(i == 9 & j == 200){next}
    print(paste("  Sample size:", j, sep=""))
    for(k in 1:10)
    {
      print(paste("    Random seed:", k, sep=""))

      # Get an index of crashes in each model and potential substitutions
      crashID <- getCrashes(paste(i, "BP-", j, "-", k, sep=""), modNames, evalDisDir, evalAssDir)
      substID <- getSubstitutes(paste(i, "BP-", j, "-", k, sep=""), modNames, evalDisDir, evalAssDir)
                 
      # Calculate the number of crashes and substitutions
      numCrash <- length(crashID)
      numSubst <- lapply(substID, length)
      
      # If there are 10 crashes and no substitution nothing is done
      if(numCrash == 10 & sum(unlist(numSubst)) == 0){next}

      # Create an index with the substitutions to be used from each iteration
      numNeeds <- list(0, 0, 0)
      if(numCrash <= numSubst[[1]])
      {
        numNeeds[[1]] <- numCrash
      }else{
        numNeeds[[1]] <- numSubst[[1]]
        numCrash2 <- numCrash - numSubst[[1]]
        if(numCrash2 <= numSubst[[2]])
        {
          numNeeds[[2]] <- numCrash2
        }else{
          numNeeds[[2]] <- numSubst[[2]]
          numCrash3 <- numCrash2 - numSubst[[2]]
          if(numCrash3 <= numSubst[[3]])
          {
            numNeeds[[3]] <- numCrash3
          }else{
            numNeeds[[3]] <- numSubst[[3]]
          }
          rm(numCrash3)
        }
        rm(numCrash2)  
      }
      
      # Load the data from the first iteration and assign new name to the objects
      load(paste(evalDisDir, "1/", i, "BP-", j, "-", k, ".RData", sep=""))
      load(paste(evalAssDir, "1/", i, "BP-", j, "-", k, ".RData", sep=""))
      for(l in 1:length(objNames))
      {
        assign(paste(objNames[l], ".Merged", sep=""), get(objNames[l]))  
      }
      
      # If there "0" crashes don't need to do anything. Go to the end and save the files
      # Otherwise, go through the whole section to replace the crashes with the substitutions
      if(numCrash != 0)
      {
        
        # Load second iteration, check if there are substitution from this iteration and 
        # extract there results
        if(numNeeds[[1]] != 0)
        {
          load(paste(evalDisDir, "2/", i, "BP-", j, "-", k, ".RData", sep=""))
          load(paste(evalAssDir, "2/", i, "BP-", j, "-", k, ".RData", sep=""))
          index <- substID[[1]][1:numNeeds[[1]]]
          for(l in 1:length(objNames))
          {
            assign(paste(objNames[l], ".Subst", sep=""), get(objNames[l])[index])
          }
        }
  
        # Load third iteration, check if there are substitution from this iteration, 
        # extract there results, and if necesary append them to those from the former iteration
        if(numNeeds[[2]] != 0)
        {
          load(paste(evalDisDir, "3/", i, "BP-", j, "-", k, ".RData", sep=""))
          load(paste(evalAssDir, "3/", i, "BP-", j, "-", k, ".RData", sep=""))
          index <- substID[[2]][1:numNeeds[[2]]]
          for(l in 1:length(objNames))
          {
            if(numNeeds[[1]] == 0)
            {
              assign(paste(objNames[l], ".Subst", sep=""), get(objNames[l])[index])
            }else{
              append(get(paste(objNames[l], ".Subst", sep="")), get(objNames[l])[index])
            }
          }
        }
          
        # Load fourth iteration, check if there are substitution from this iteration, 
        # extract there results, and if necesary append them to those from the former iterations
        if(numNeeds[[3]] != 0)
        {
          load(paste(evalDisDir, "4/", i, "BP-", j, "-", k, ".RData", sep=""))
          load(paste(evalAssDir, "4/", i, "BP-", j, "-", k, ".RData", sep=""))
          index <- substID[[3]][1:numNeeds[[3]]]
          for(l in 1:length(objNames))
          {
            if(numNeeds[[1]] == 0 & numNeeds[[2]] == 0)
            {
              assign(paste(objNames[l], ".Subst", sep=""), get(objNames[l])[index])
            }else{
              append(get(paste(objNames[l], ".Subst", sep="")), get(objNames[l])[index])
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
          
          rm(tmpObj.Merged)
        }
      }

      # Save assemblages results to the HDD 
      direct <- "Results/RObjects/Evaluation/Assemblage/Merged"
      dir.create(direct, recursive=T) 
      save(glmEnseTrain.Merged, glmEnseTest.Merged, cqoEnseTrain.Merged, cqoEnseTest.Merged, gamEnseTrain.Merged, gamEnseTest.Merged, caoEnseTrain.Merged, caoEnseTest.Merged, scartsEnseTrain.Merged, scartsEnseTest.Merged, cartsEnseTrain.Merged, cartsEnseTest.Merged, smarsEnseTrain.Merged, smarsEnseTest.Merged, marsEnseTrain.Merged, marsEnseTest.Merged, snnEnseTrain.Merged, snnEnseTest.Merged, mnnEnseTrain.Merged, mnnEnseTest.Merged, file=paste(direct, "/", i, "BP-", j, "-", k, ".RData", sep=""))
    
      # Save distribution results to the HDD
      direct <- "Results/RObjects/Evaluation/Distribution/Merged"
      dir.create(direct, recursive=T) 
      save(glmDistTrain.Merged, glmDistTest.Merged, cqoDistTrain.Merged, cqoDistTest.Merged, gamDistTrain.Merged, gamDistTest.Merged, caoDistTrain.Merged, caoDistTest.Merged, scartsDistTrain.Merged, scartsDistTest.Merged, cartsDistTrain.Merged, cartsDistTest.Merged, smarsDistTrain.Merged, smarsDistTest.Merged, marsDistTrain.Merged, marsDistTest.Merged, snnDistTrain.Merged, snnDistTest.Merged, mnnDistTrain.Merged, mnnDistTest.Merged, file=paste(direct, "/", i, "BP-", j, "-", k, ".RData", sep=""))
      
      # Remove temporary objects
      rm(direct, tmpObj.Merged, tmpObj.Subst, index, numNeeds, numCrash, numSubst, crashID, substID) 
      rm(list=ls(pattern="DistTest"))
      rm(list=ls(pattern="DistTrain"))
      rm(list=ls(pattern="EnseTest"))
      rm(list=ls(pattern="EnseTrain"))
    }
  }
}



################################################################################
#
# FORMAT RESULTS
#
################################################################################

objDTestNames <- c("glmDistTest", "cqoDistTest", "gamDistTest", "caoDistTest", "scartsDistTest", "cartsDistTest", "smarsDistTest", "marsDistTest", "snnDistTest", "mnnDistTest")
objDTrainNames <- c("glmDistTrain", "cqoDistTrain", "gamDistTrain", "caoDistTrain", "scartsDistTrain", "cartsDistTrain", "smarsDistTrain", "marsDistTrain", "snnDistTrain", "mnnDistTrain")
objETestNames <- c("glmEnseTest", "cqoEnseTest", "gamEnseTest", "caoEnseTest", "scartsEnseTest", "cartsEnseTest", "smarsEnseTest", "marsEnseTest", "snnEnseTest", "mnnEnseTest")
objETrainNames <- c("glmEnseTrain", "cqoEnseTrain", "gamEnseTrain", "caoEnseTrain", "scartsEnseTrain", "cartsEnseTrain", "smarsEnseTrain", "marsEnseTrain", "snnEnseTrain", "mnnEnseTrain")


# List to store results
resd.All <- list()
rese.All <- list()

sampleSizes <- c(300, 200, 150, 100, 50)

files <- strsplit(list.files("Results/RObjects/Evaluation/Distribution/Merged"), split='[.]')
files <- sapply(files, FUN=function(x){x <- x[[1]]; return(x)})

# FOR loop for each period fit, sample size and iteration
a <- 1
for(i in periods_fit)
{
  print(paste("Period fit:", i, "BP", sep=""))
  for(j in sampleSizes)
  {
    if(i == 9 & j == 200){next}
    print(paste("  Sample size:", j, sep=""))
    for(k in 1:10)
    {
      print(paste("    Random seed:", k, sep=""))

      # Skip
      if(j == 300)
      {
        if(k == 1)
        {
          # Load the results from the merged folders
          load(paste("Results/RObjects/Evaluation/Assemblage/Merged/", i, "BP.RData", sep=""), .GlobalEnv)
          load(paste("Results/RObjects/Evaluation/Distribution/Merged/", i, "BP.RData", sep=""), .GlobalEnv)
        }else{
          # For the original data set there is only one iteration, the others are skipped
          next
        }
      }else{
        if(paste(i, "BP-", j, "-", k, sep="") %in% files){
          # Load the results from the merged folders
          load(paste("Results/RObjects/Evaluation/Assemblage/Merged/", i, "BP-", j, "-", k, ".RData", sep=""), .GlobalEnv)
          load(paste("Results/RObjects/Evaluation/Distribution/Merged/", i, "BP-", j, "-", k, ".RData", sep=""), .GlobalEnv)
        }else{
          # There are a couple of combinations (periods fit, sample size and iteration) that have no data and wont be loaded.
          next
        }
      }
    
      # Reshape all the data in a dataframe. FOR loops are for all models (GLM, CQO, etc) for each type of results object (distribution vs assemblage and train vs test).
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
            names(tmpName[[w]]) <- c("jacc.MEAN.prob","jacc.SD.prob","jacc.CAL.prob","jacc.DIS.prob","spRich.COR.prob","spRich.RSQ.prob","jacc.MEAN.thre","jacc.SD.thre","jacc.CAL.thre","jacc.DIS.thre","spRich.COR.thre","spRich.RSQ.thre")
          }
        }
        assign(objETrainNames[v], tmpName) 
      }
      for(v in 1:length(objETestNames)){
        assign("tmpName", get(paste(objETestNames[v], ".Merged", sep=""))) 
        ind <- which(unlist(lapply(lapply(tmpName, class), function(x){any(x == "error")})))
        if(length(ind) > 0){
          for(w in ind){
            tmpName[[w]] <- matrix(-9999, nrow=12, ncol=12, dimnames=list(c("jacc.MEAN.prob","jacc.SD.prob","jacc.CAL.prob","jacc.DIS.prob","spRich.COR.prob","spRich.RSQ.prob","jacc.MEAN.thre","jacc.SD.thre","jacc.CAL.thre","jacc.DIS.thre","spRich.COR.thre","spRich.RSQ.thre"), periods_prj))
          }
        }
        assign(objETestNames[v], tmpName) 
      }
      
      library(reshape2)
      
      # Melting results for TRAIN
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
      resDt$periodFit <- as.character(i)
      resDt$periodPrj <- as.character(i)
      resDt$test <- "train"
      resDt$sampleSize <- j
      resDt$iteration2 <- k
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
      resEt$periodFit <- as.character(i)
      resEt$periodPrj <- as.character(i)
      resEt$test <- "train"
      resEt$sampleSize <- j
      resEt$iteration2 <- k
      rm(glmet, gamet, cqoet, caoet, marset, smarset, cartset, scartset, mnnet, snnet)
      
      # Melting results for TEST
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
      resD$periodFit <- as.character(i)
      resD$test <- "test"
      resD$sampleSize <- j
      resD$iteration2 <- k

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
      resE$periodFit <- as.character(i)
      resE$sampleSize <- j
      resE$iteration2 <- k

      rm(glme, game, cqoe, caoe, marse, smarse, cartse, scartse, mnne, snne)
      
      # Combining the TRAINING and TESTING datasets
      resd.All[[a]] <- rbind(resDt, resD)
      rese.All[[a]] <- rbind(resEt, resE)
    
      # Remove objects
      rm(list=ls(pattern="[.]Merged"))
      rm(list=ls(pattern="DistTest"))
      rm(list=ls(pattern="DistTrain"))
      rm(list=ls(pattern="EnseTest"))
      rm(list=ls(pattern="EnseTrain"))
      rm(resDt, resD, resEt, resE)
      
      a <- a+1
    }
  }
}

# Remove objects
rm(i, j, k, a, v, w, ind, tmp, tmpName)
rm(objNames, objDTestNames, objDTrainNames, objETestNames, objETrainNames)

# Write data to file
save(resd.All, file="Results/RObjects/DisResults.RData")
save(rese.All, file="Results/RObjects/AssResults.RData")


################################################################################
##
## ARRANGING DATA
##
################################################################################
##Load plotting libraries
library(reshape2)
library(ggplot2)

##Load data
load("Results/RObjects/DisResults.RData")
load("Results/RObjects/AssResults.RData")

##Remove the missing data and melt the data in a single dataframe
resd.all <- melt(resd.All, id.vars=c("taxon","variable","iteration","model","periodFit","periodPrj","test","sampleSize","iteration2"))
rese.all <- melt(rese.All, id.vars=c("variable","iteration","model","periodFit","periodPrj","test","sampleSize","iteration2"))

# Change -9999 values for NA
resd.all$value[which(resd.all$value == -9999)] <- NA
rese.all$value[which(rese.all$value == -9999)] <- NA

##Change column names in the new dataframe
resd.all <- resd.all[,-which(colnames(resd.all) == "L1")]
rese.all <- rese.all[,-which(colnames(rese.all) == "L1")]

# Create new column with modelType information (SDM vs CLM)
resd.all$modelType <- NA
rese.all$modelType <- NA
resd.all$modelType[which(resd.all$model %in% c("glm","gam","scarts","smars","snn"))] <- "SDM"
rese.all$modelType[which(rese.all$model %in% c("glm","gam","scarts","smars","snn"))] <- "SDM"
resd.all$modelType[which(resd.all$model %in% c("cqo","cao","carts","mars","mnn"))] <- "CLM"
rese.all$modelType[which(rese.all$model %in% c("cqo","cao","carts","mars","mnn"))] <- "CLM"

# Define several factors, that way they will be plotted in order afterwards
resd.all$model <- factor(resd.all$model, c("glm","gam","scarts","smars","snn","cqo","cao","carts","mars","mnn"))
rese.all$model <- factor(rese.all$model, c("glm","gam","scarts","smars","snn","cqo","cao","carts","mars","mnn"))             
resd.all$periodFit <- factor(resd.all$periodFit, periods_fit)
rese.all$periodFit <- factor(rese.all$periodFit, periods_fit)             
resd.all$periodPrj <- factor(resd.all$periodPrj, levels=periods_prj)
rese.all$periodPrj <- factor(rese.all$periodPrj, levels=periods_prj)

resd.all$modelType <- factor(resd.all$modelType, levels=c("SDM","CLM"))
rese.all$modelType <- factor(rese.all$modelType, levels=c("SDM","CLM"))

resd.all$variable <- factor(resd.all$variable, c("auc","tss","brier","sens","spec","prev","trmin"))
rese.all$variable <- factor(rese.all$variable, c("jacc.MEAN.prob","jacc.SD.prob","jacc.CAL.prob","jacc.DIS.prob","spRich.COR.prob","spRich.RSQ.prob","jacc.MEAN.thre","jacc.SD.thre","jacc.CAL.thre","jacc.DIS.thre","spRich.COR.thre","spRich.RSQ.thre"))

resd.all$test <- factor(resd.all$test, levels=c("train","test"))
rese.all$test <- factor(rese.all$test, levels=c("train","test"))

resd.all$sampleSize[which(resd.all$sampleSize == 300)] <- "Or."
rese.all$sampleSize[which(rese.all$sampleSize == 300)] <- "Or."
resd.all$sampleSize <- factor(resd.all$sampleSize, levels=c("50","100","150","200","Or."))
rese.all$sampleSize <- factor(rese.all$sampleSize, levels=c("50","100","150","200","Or."))

resd.all$iteration <- factor(resd.all$iteration, levels=1:10)
rese.all$iteration <- factor(rese.all$iteration, levels=1:10)
resd.all$iteration2 <- factor(resd.all$iteration2, levels=1:10)
rese.all$iteration2 <- factor(rese.all$iteration2, levels=1:10)

# Split the results in TRAIN and TESTING
resd.train <- resd.all[which(as.character(resd.all$periodFit) == as.character(resd.all$periodPrj)),]
rese.train <- rese.all[which(as.character(rese.all$periodFit) == as.character(rese.all$periodPrj)),]

resd.test <- resd.all[which(resd.all$test == "test"),]
rese.test <- rese.all[which(rese.all$test == "test"),]


#################################################################################
###
### SIGNIFICANCE TEST BETWEEN SDMs AND CLMs
###
#################################################################################
### function for calculating t-test for each metric
sdmStats <- function(dat
                    , p_prj
                    , p_fit
                    , measure
                    )
{
  tmp <- matrix(NA, length(p_prj), length(p_fit))
  rownames(tmp) <- p_prj
  colnames(tmp) <- p_fit
  stat.glm <- stat.gam <- stat.scarts <- stat.smars <- stat.snnet <- list(tmp, tmp, tmp, tmp, tmp)
  names(stat.glm) <- names(stat.gam) <- names(stat.scarts) <- names(stat.smars) <- names(stat.snnet) <- unique(dat$sampleSize)

  dat <- dat[which(dat$variable == measure),]
  for (l in 1:length(p_fit))
  {
    i <- p_fit[l]
    print(paste("Working on", i))
    for (j in 1:length(p_prj))
    {
      k <- p_prj[j]
      print(paste("  and", k))
      for (m in 1:length(unique(dat$sampleSize)))
      {
        n <- unique(dat$sampleSize)[m]
        print(paste("    Sample size:", m, sep="")) 

        dat1 <- dat[which(dat$periodPrj == k & dat$periodFit == i & dat$sampleSize == n),]
  
        if(length(unique(dat1$iteration)) < 3)
        {
          stat.glm[[m]][j,l] <- stat.gam[[m]][j,l] <- stat.scarts[[m]][j,l] <- stat.smars[[m]][j,l] <- stat.snnet[[m]][j,l] <- NA
        }else{
          sdm <- dat1[which(dat1$model == "glm"),]
          clm <- dat1[which(dat1$model == "cqo"),]
          results <- t.test(sdm$value, clm$value, paired=T)
          stat.glm[[m]][j,l] <- results$p.value
          
          sdm <- dat1[which(dat1$model == "gam"),]
          clm <- dat1[which(dat1$model == "cao"),]
          results <- t.test(sdm$value, clm$value, paired=T)
          stat.gam[[m]][j,l] <- results$p.value
          
          sdm <- dat1[which(dat1$model == "scarts"),]
          clm <- dat1[which(dat1$model == "carts"),]
          results <- t.test(sdm$value, clm$value, paired=T)
          stat.scarts[[m]][j,l] <- results$p.value
          
          sdm <- dat1[which(dat1$model == "smars"),]
          clm <- dat1[which(dat1$model == "mars"),]
          results <- t.test(sdm$value, clm$value, paired=T)
          stat.smars[[m]][j,l] <- results$p.value
          
          sdm <- dat1[which(dat1$model == "snn"),]
          clm <- dat1[which(dat1$model == "mnn"),]
          results <- t.test(sdm$value, clm$value, paired=T)
          stat.snnet[[m]][j,l] <- results$p.value
        }
      }
    }
  }        
 
  stat.glm <- melt(stat.glm)
  stat.glm$models <- "glm_cqo"
  stat.gam <- melt(stat.gam)
  stat.gam$models <- "gam_cao"
  stat.scarts <- melt(stat.scarts)
  stat.scarts$models <- "scarts_carts"
  stat.smars <- melt(stat.smars)
  stat.smars$models <- "smars_mars"
  stat.snnet <- melt(stat.snnet)
  stat.snnet$models <- "snn_mnn"
  
  stat.Data <- rbind(stat.glm, stat.gam, stat.scarts, stat.smars, stat.snnet)
  colnames(stat.Data) <- c("periodPrj","periodFit","p_value","sampleSize","models")
  stat.Data$models <- factor(stat.Data$models)
  stat.Data$sampleSize <- factor(stat.Data$sampleSize, levels=c("50","100","150","200","Or."))
  stat.Data$periodFit <- factor(stat.Data$periodFit)
  return(stat.Data)
}

## apply the function for each metric and save the results
auc.ttest <- sdmStats(resd.test, periods_prj, periods_fit, "auc")
save(auc.ttest, file="Results/RObjects/Ttest/auc.ttest")

brier.ttest <- sdmStats(resd.test, periods_prj, periods_fit, "brier")
save(brier.ttest, file="Results/RObjects/Ttest/brier.ttest")

tss.ttest <- sdmStats(resd.test, periods_prj, periods_fit, "tss")
save(tss.ttest, file="Results/RObjects/Ttest/tss.ttest")

sens.ttest <- sdmStats(resd.test, periods_prj, periods_fit, "sens")
save(sens.ttest, file="Results/RObjects/Ttest/sens.ttest")

spec.ttest <- sdmStats(resd.test, periods_prj, periods_fit, "spec")
save(spec.ttest, file="Results/RObjects/Ttest/spec.ttest")

jacc.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "jacc.MEAN.prob")
save(jacc.ttest, file="Results/RObjects/Ttest/jacc.ttest")

jacDis.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "jacc.DIS.prob")
save(jacDis.ttest, file="Results/RObjects/Ttest/jacDis.ttest")

jacCal.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "jacc.CAL.prob")
save(jacCal.ttest, file="Results/RObjects/Ttest/jacCal.ttest")

spRichCor.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "spRich.COR.prob")
save(spRichCor.ttest, file="Results/RObjects/Ttest/spRichCor.ttest")

spRichRsq.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "spRich.RSQ.prob")
save(spRichRsq.ttest, file="Results/RObjects/Ttest/spRichRsq.ttest")

jacct.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "jacc.MEAN.thre")
save(jacct.ttest, file="Results/RObjects/Ttest/jacct.ttest")

jacDist.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "jacc.DIS.thre")
save(jacDist.ttest, file="Results/RObjects/Ttest/jacDist.ttest")

jacCalt.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "jacc.CAL.thre")
save(jacCalt.ttest, file="Results/RObjects/Ttest/jacCalt.ttest")

spRichCort.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "spRich.COR.thre")
save(spRichCort.ttest, file="Results/RObjects/Ttest/spRichCort.ttest")

spRichRsqt.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "spRich.RSQ.thre")
save(spRichRsqt.ttest, file="Results/RObjects/Ttest/spRichRsqt.ttest")

### Function for testing the difference betweeen SDM and CLM using a t-test
sdmStats_ModelType <- function(dat
                              , p_prj
                              , p_fit
                              , measure
                              )
{
# dat <- resd.test
# p_fit <- periods_fit
# p_prj <- periods_prj
# measure <- "auc"
# rm(dat, p_fit, p_prj, measure)
# rm(tmp, stat.Data, l, i, j, k, m, n, dat1, sdm, clm, ind, results, stat.Data)
  tmp <- matrix(NA, length(p_prj), length(p_fit))
  rownames(tmp) <- p_prj
  colnames(tmp) <- p_fit

  stat.Data <- list(tmp, tmp, tmp, tmp, tmp)
  names(stat.Data) <- unique(dat$sampleSize)
  dat <- dat[which(dat$variable == measure),]

  for(l in 1:length(p_fit)){
    i <- p_fit[l]
    print(paste("Working on", i))
    for(j in 1:length(p_prj)){
      k <- p_prj[j]
      print(paste("  and", k))
      for(m in 1:length(unique(dat$sampleSize))){
        n <- unique(dat$sampleSize)[m]
        print(paste("    Sample size:", m, sep="")) 
  
        dat1 <- dat[which(dat$periodPrj == k & dat$periodFit == i & dat$sampleSize == n),]
  
        if(length(unique(dat1$iteration)) < 3)
        {
          stat.Data[[m]][j,l] <- NA
        }else{
          sdm <- dat1[which(dat1$modelType == "SDM"),]
          clm <- dat1[which(dat1$modelType == "CLM"),]
    
          ind <- intersect(unique(sdm$iteration), unique(clm$iteration))
          sdm <- sdm[which(sdm$iteration %in% ind),]
          clm <- clm[which(clm$iteration %in% ind),]
          results <- t.test(sdm$value, clm$value, paired=T)
          stat.Data[[m]][j,l] <- results$p.value
        }
      }
    }
  }
  stat.Data <- melt(stat.Data)
  colnames(stat.Data) <- c("periodPrj","periodFit","p_value","sampleSize")
  stat.Data$variable <- measure
  return(stat.Data)
}

## Apply the function for each metric
auc.modelType <- sdmStats_ModelType(resd.test, periods_prj, periods_fit, "auc")
brier.modelType <- sdmStats_ModelType(resd.test, periods_prj, periods_fit, "brier")
tss.modelType <- sdmStats_ModelType(resd.test, periods_prj, periods_fit, "tss")
sens.modelType <- sdmStats_ModelType(resd.test, periods_prj, periods_fit, "sens")
spec.modelType <- sdmStats_ModelType(resd.test, periods_prj, periods_fit, "spec")

jacc.modelType <- sdmStats_ModelType(rese.test, periods_prj, periods_fit, "jacc.MEAN.prob")
jacCal.modelType <- sdmStats_ModelType(rese.test, periods_prj, periods_fit, "jacc.CAL.prob")
jacDis.modelType <- sdmStats_ModelType(rese.test, periods_prj, periods_fit, "jacc.DIS.prob")
spRichCor.modelType <- sdmStats_ModelType(rese.test, periods_prj, periods_fit, "spRich.COR.prob")
spRichRsq.modelType <- sdmStats_ModelType(rese.test, periods_prj, periods_fit, "spRich.RSQ.prob")

jacct.modelType <- sdmStats_ModelType(rese.test, periods_prj, periods_fit, "jacc.MEAN.thre")
jacCalt.modelType <- sdmStats_ModelType(rese.test, periods_prj, periods_fit, "jacc.CAL.thre")
jacDist.modelType <- sdmStats_ModelType(rese.test, periods_prj, periods_fit, "jacc.DIS.thre")
spRichCort.modelType <- sdmStats_ModelType(rese.test, periods_prj, periods_fit, "spRich.COR.thre")
spRichRsqt.modelType <- sdmStats_ModelType(rese.test, periods_prj, periods_fit, "spRich.RSQ.thre")

### Combined results and save
stats.modelType <- rbind(auc.modelType, tss.modelType, sens.modelType, spec.modelType, brier.modelType, jacc.modelType, jacCal.modelType, jacDis.modelType, spRichCor.modelType, spRichRsq.modelType, jacct.modelType, jacCalt.modelType, jacDist.modelType, spRichCort.modelType, spRichRsqt.modelType)

save(stats.modelType, file="Results/RObjects/Ttest/stats.modelType.R")


load("Results/RObjects/Ttest/auc.ttest")
load("Results/RObjects/Ttest/tss.ttest")
load("Results/RObjects/Ttest/spec.ttest")
load("Results/RObjects/Ttest/sens.ttest")
load("Results/RObjects/Ttest/brier.ttest")
load("Results/RObjects/Ttest/jacc.ttest")
load("Results/RObjects/Ttest/jacCal.ttest")
load("Results/RObjects/Ttest/jacDis.ttest")
load("Results/RObjects/Ttest/spRichCor.ttest")
load("Results/RObjects/Ttest/spRichRsq.ttest")
load("Results/RObjects/Ttest/stats.modelType.R")
load("Results/RObjects/Ttest/jacct.ttest")
load("Results/RObjects/Ttest/jacCalt.ttest")
load("Results/RObjects/Ttest/jacDist.ttest")
load("Results/RObjects/Ttest/spRichCort.ttest")
load("Results/RObjects/Ttest/spRichRsqt.ttest")
load("Results/RObjects/Ttest/stats.modelType.R")

## Rearrange p-value data values for each model
## Describe the function to do so
orgData <- function(x, per, varName){
  tmpObj <- strsplit(as.character(x$models), "_")
  
  x$model.x <- factor(sapply(tmpObj, function(x){x[[1]]}))
  x$model.y <- factor(sapply(tmpObj, function(x){x[[2]]}))
  
  disLevels <- c("prev","auc","trmin","tss","sens","spec","brier")
  assLevels <- c("jacc.MEAN.prob","jacc.SD.prob","jacc.CAL.prob","jacc.DIS.prob","spRich.COR.prob","spRich.RSQ.prob","jacc.MEAN.thre","jacc.SD.thre","jacc.CAL.thre","jacc.DIS.thre","spRich.COR.thre","spRich.RSQ.thre")
  
  if(varName %in% disLevels){
    x$variable <- factor(varName, levels=disLevels)
  }
  if(varName %in% assLevels){
    x$variable <- factor(varName, levels=assLevels)
  }
  
  x$significant <- as.numeric(x$p_value < 0.001)
  x$significant2 <- x$significant
  x$significant2[!x$significant2] <- NA
  x$significant <- factor(x$significant, labels=c("n.s.","< 0.001"))
  x$significant2 <- factor(x$significant2, labels=c("< 0.001"))
  
  return(x)
}

## Apply the function to p-value files for each model from above
aucStats.Models <- orgData(auc.ttest, periods_prj, "auc")
tssStats.Models <- orgData(tss.ttest, periods_prj, "tss")
sensStats.Models <- orgData(sens.ttest, periods_prj, "sens")
specStats.Models <- orgData(spec.ttest, periods_prj, "spec")
brierStats.Models <- orgData(brier.ttest, periods_prj, "brier")

jaccStats.Models <- orgData(jacc.ttest, periods_prj, "jacc.MEAN.prob")
jacDisStats.Models <- orgData(jacDis.ttest, periods_prj, "jacc.DIS.prob")
jacCalStats.Models <- orgData(jacCal.ttest, periods_prj, "jacc.CAL.prob")
spRichCorStats.Models <- orgData(spRichCor.ttest, periods_prj, "spRich.COR.prob")
spRichRsqStats.Models <- orgData(spRichRsq.ttest, periods_prj, "spRich.RSQ.prob")

jacctStats.Models <- orgData(jacct.ttest, periods_prj, "jacc.MEAN.thre")
jacDistStats.Models <- orgData(jacDist.ttest, periods_prj, "jacc.DIS.thre")
jacCaltStats.Models <- orgData(jacCalt.ttest, periods_prj, "jacc.CAL.thre")
spRichCortStats.Models <- orgData(spRichCort.ttest, periods_prj, "spRich.COR.thre")
spRichRsqtStats.Models <- orgData(spRichRsqt.ttest, periods_prj, "spRich.RSQ.thre")

## Combine the data
disStats.Models <- rbind(aucStats.Models, tssStats.Models, sensStats.Models, specStats.Models, brierStats.Models)
assStats.Models <- rbind(jaccStats.Models, jacDisStats.Models, jacCalStats.Models, spRichCorStats.Models, spRichRsqStats.Models, jacctStats.Models, jacDistStats.Models, jacCaltStats.Models, spRichCortStats.Models, spRichRsqtStats.Models)

## Rearrange p-value data values for model type
## Describe the function to do so
orgData <- function(x, per){

  x$variable <- factor(x$variable)
  
  x$significant <- as.numeric(x$p_value < 0.001)
  x$significant2 <- x$significant
  x$significant2[!x$significant2] <- NA
  x$significant <- factor(x$significant, labels=c("n.s.","< 0.001"))
  x$significant2 <- factor(x$significant2, labels=c("< 0.001"))
  
  return(x)
}

## Apply the function to the p-value file for model type from above
stats.ModType <- orgData(stats.modelType, periods_prj)

rm(resd.All, rese.All, sdmStats, sdmStats_ModelType, orgData)
rm(auc.ttest, tss.ttest, sens.ttest, spec.ttest, brier.ttest, jacc.ttest, jacCal.ttest, jacDis.ttest, spRichCor.ttest, spRichRsq.ttest, jacct.ttest, jacCalt.ttest, jacDist.ttest, spRichCort.ttest, spRichRsqt.ttest)
rm(auc.modelType, tss.modelType, sens.modelType, spec.modelType, brier.modelType, jacc.modelType, jacCal.modelType, jacDis.modelType, spRichCor.modelType, spRichRsq.modelType, jacct.modelType, jacCalt.modelType, jacDist.modelType, spRichCort.modelType, spRichRsqt.modelType)
rm(aucStats.Models, tssStats.Models, sensStats.Models, specStats.Models, brierStats.Models, jaccStats.Models, jacCalStats.Models, jacDisStats.Models, spRichCorStats.Models, spRichRsqStats.Models, jacctStats.Models, jacCaltStats.Models, jacDistStats.Models, spRichCortStats.Models, spRichRsqtStats.Models)
rm(stats.modelType)


################################################################################
#
# TESTING ACROSS TIME DATA
#
################################################################################

## De-melt the data to calculate mean values. I calculated means across iterations, taxon, and model for ploting results in CLMvsSDM; and across iterations, taxon, and model_type for plotting results for each specific model (cqo, glm, ...). For assemblage assessments taxon was obviously disregarded.
index.dColumns <- which(colnames(resd.test) %in% c("iteration","taxon","test","modelType","iteration2"))
index.eColumns <- which(colnames(rese.test) %in% c("iteration","test","modelType","iteration2"))

resd.Models.Mean <- acast(resd.test[,-index.dColumns], periodPrj ~ periodFit ~ variable ~ model ~ sampleSize, mean, na.rm=T)
rese.Models.Mean <- acast(rese.test[,-index.eColumns], periodPrj ~ periodFit ~ variable ~ model ~ sampleSize, mean, na.rm=T)

index.dColumns <- which(colnames(resd.test) %in% c("iteration","taxon","test","model","iteration2"))
index.eColumns <- which(colnames(rese.test) %in% c("iteration","test","model","iteration2"))

resd.ModType.Mean <- acast(resd.test[,-index.dColumns], periodPrj ~ periodFit ~ variable ~ modelType ~ sampleSize, mean, na.rm=T)
rese.ModType.Mean <- acast(rese.test[,-index.eColumns], periodPrj ~ periodFit ~ variable ~ modelType ~ sampleSize, mean, na.rm=T)

## Melting again the new data (means and sd) to get the data in shape for the ggplot functions.
resd.Models.Mean <- melt(resd.Models.Mean)
rese.Models.Mean <- melt(rese.Models.Mean)

resd.ModType.Mean <- melt(resd.ModType.Mean)
rese.ModType.Mean <- melt(rese.ModType.Mean)

## Change colnames of the new dataframes
colnames(resd.Models.Mean) <- colnames(rese.Models.Mean) <- c("periodPrj","periodFit","variable","model","sampleSize","value")
colnames(resd.ModType.Mean) <- colnames(rese.ModType.Mean) <- c("periodPrj","periodFit","variable","modelType","sampleSize", "value")

### Define two new columns in the models dataframes. The idea is to have a CLM versus SDM column and a SDM equivalent column. 
newModelTypeColumn <- function(d){
  d$modelType <- NA
  d$modelType[which(d$model %in% c("cqo","cao","carts","mars","mnn"))] <- "CLM"
  d$modelType[which(d$model %in% c("glm","gam","scarts","smars","snn"))] <- "SDM"
  d$modelType <- factor(d$modelType)
  return(d)
}

newSdmEqColumn <- function(d){
  d$sdm_eq <- factor(NA, levels=(c("glm-cqo","gam-cao","scarts-carts","smars-mars","snn-mnn")))
  d$sdm_eq[which(d$model %in% c("glm", "cqo"))] <- "glm-cqo"
  d$sdm_eq[which(d$model %in% c("gam", "cao"))] <- "gam-cao"
  d$sdm_eq[which(d$model %in% c("smars", "mars"))] <- "smars-mars"
  d$sdm_eq[which(d$model %in% c("scarts", "carts"))] <- "scarts-carts"
  d$sdm_eq[which(d$model %in% c("snn", "mnn"))] <- "snn-mnn"
  return(d)
}

# Apply the functions
resd.Models.Mean <- newModelTypeColumn(resd.Models.Mean)
rese.Models.Mean <- newModelTypeColumn(rese.Models.Mean)

resd.Models.Mean <- newSdmEqColumn(resd.Models.Mean)
rese.Models.Mean <- newSdmEqColumn(rese.Models.Mean)

rm(index.dColumns, index.eColumns, newModelTypeColumn, newSdmEqColumn)

################################################################################
##
## DIFFERENCES - Calculate the difference between mean SDM and mean CLM values
##
################################################################################

## Rearrange the data
calDif <- function(d)
{
  if(ncol(d) == 8){
    p <- merge(d, d, by=c("periodPrj","periodFit","sdm_eq","variable","sampleSize"))
  }
  if(ncol(d) == 6){
    p <- merge(d, d, by=c("periodPrj","periodFit","variable","sampleSize"))
  }
  p <- p[which(p$modelType.x == "SDM" & p$modelType.y == "CLM"),]
  
  index1 <- which(p$variable %in% c("auc","tss","sens","spec","jacc.CAL.prob","jacc.CAL.thre","jacc.DIS.prob","jacc.DIS.thre","spRich.RSQ.prob","spRich.RSQ.thre"))
  index2 <- which(p$variable %in% c("spRich.COR.prob","spRich.COR.thre"))
  index3 <- which(p$variable %in% c("jacc.SD.prob","jacc.SD.thre"))
  index4 <- which(p$variable %in% c("brier","jacc.MEAN.prob","jacc.MEAN.thre"))
  
  p$dif <- NA
  
  p$dif[index1] <- p$value.x[index1] - p$value.y[index1]
  
  p$dif[index4] <- p$value.y[index4] - p$value.x[index4]
  
  tmp <- pmax(p$value.x[index2], p$value.y[index2])
  tmp[which(tmp == p$value.x[index2] & tmp != p$value.y[index2])] <- 1
  tmp[which(tmp != p$value.x[index2] & tmp == p$value.y[index2])] <- -1
  tmp[which(tmp == p$value.x[index2] & tmp == p$value.y[index2])] <- 0
  p$dif[index2] <- abs(p$value.x[index2] - p$value.y[index2]) * tmp
  
  tmp <- pmin(p$value.x[index3], p$value.y[index3])
  tmp[which(tmp == p$value.x[index3] & tmp != p$value.y[index3])] <- 1
  tmp[which(tmp != p$value.x[index3] & tmp == p$value.y[index3])] <- -1
  tmp[which(tmp == p$value.x[index3] & tmp == p$value.y[index3])] <- 0
  p$dif[index3] <- abs(p$value.x[index3] - p$value.y[index3]) * tmp
  
  return(p)
}

## Apply the functions
resd.Models.Dif <- calDif(resd.Models.Mean)
rese.Models.Dif <- calDif(rese.Models.Mean)

resd.ModType.Dif <- calDif(resd.ModType.Mean)
rese.ModType.Dif <- calDif(rese.ModType.Mean)

rm(calDif)

################################################################################
#
# Merge Difference and p-value datasets for plotting purposes
#
################################################################################
# Merge the data frames
resd.Models.Dif2 <- merge(resd.Models.Dif, disStats.Models, by=c("periodPrj", "periodFit", "sampleSize","model.x", "model.y", "variable"))
rese.Models.Dif2 <- merge(rese.Models.Dif, assStats.Models, by=c("periodPrj", "periodFit", "sampleSize","model.x", "model.y", "variable"))

resd.ModType.Dif2 <- merge(resd.ModType.Dif, stats.ModType, by=c("periodPrj", "periodFit", "sampleSize","variable"))
rese.ModType.Dif2 <- merge(rese.ModType.Dif, stats.ModType, by=c("periodPrj", "periodFit", "sampleSize","variable"))

# Save them to the HDD
save(resd.Models.Mean, file="Results/RObjects/resdModelsMean_sampleSize.RData")
save(rese.Models.Mean, file="Results/RObjects/reseModelsMean_sampleSize.RData")
save(resd.ModType.Mean, file="Results/RObjects/resdModTypeMean_sampleSize.RData")
save(rese.ModType.Mean, file="Results/RObjects/reseModTypeMean_sampleSize.RData")
save(resd.Models.Dif2, file="Results/RObjects/resdModelsDif2_sampleSize.RData")
save(rese.Models.Dif2, file="Results/RObjects/reseModelsDif2_sampleSize.RData")
save(resd.ModType.Dif2, file="Results/RObjects/resdModTypeDif2_sampleSize.RData")
save(rese.ModType.Dif2, file="Results/RObjects/reseModTypeDif2_sampleSize.RData")


################################################################################
##
## PLOTTING RESULTS
##
################################################################################
# Load objects with the results ready to be ploted
load("Results/RObjects/resdModelsDif2_sampleSize.RData", .GlobalEnv)
load("Results/RObjects/reseModelsDif2_sampleSize.RData", .GlobalEnv)
load("Results/RObjects/resdModTypeDif2_sampleSize.RData", .GlobalEnv)
load("Results/RObjects/reseModTypeDif2_sampleSize.RData", .GlobalEnv)
load("Results/RObjects/resdModelsMean_sampleSize.RData", .GlobalEnv)
load("Results/RObjects/reseModelsMean_sampleSize.RData", .GlobalEnv)
load("Results/RObjects/resdModTypeMean_sampleSize.RData", .GlobalEnv)
load("Results/RObjects/reseModTypeMean_sampleSize.RData", .GlobalEnv)

# Remove records for 9BP and n=200 of sample size
resd.Models.Mean <- resd.Models.Mean[!(resd.Models.Mean$sampleSize == 200 & resd.Models.Mean$periodFit == 9),]
rese.Models.Mean <- rese.Models.Mean[!(rese.Models.Mean$sampleSize == 200 & rese.Models.Mean$periodFit == 9),]
resd.ModType.Mean <- resd.ModType.Mean[!(resd.ModType.Mean$sampleSize == 200 & resd.ModType.Mean$periodFit == 9),]
rese.ModType.Mean <- rese.ModType.Mean[!(rese.ModType.Mean$sampleSize == 200 & rese.ModType.Mean$periodFit == 9),]
resd.Models.Dif2 <- resd.Models.Dif2[!(resd.Models.Dif2$sampleSize == 200 & resd.Models.Dif2$periodFit == 9),]
rese.Models.Dif2 <- rese.Models.Dif2[!(rese.Models.Dif2$sampleSize == 200 & rese.Models.Dif2$periodFit == 9),]
resd.ModType.Dif2 <- resd.ModType.Dif2[!(resd.ModType.Dif2$sampleSize == 200 & resd.ModType.Dif2$periodFit == 9),]
rese.ModType.Dif2 <- rese.ModType.Dif2[!(rese.ModType.Dif2$sampleSize == 200 & rese.ModType.Dif2$periodFit == 9),]

# Define some common details to all the plots
themeOpt <- theme(axis.text.x=element_text(angle=90, hjust=1)) + theme_bw()

axesLabels <- labs(x="Projecting period", y="Sample size")

aspRatio <- coord_fixed(ratio=1)


#################################################################################
## CLM versus SDM comparison
#################################################################################

########
# 1 BP

# AUC
pdf("Results/Graphs/AUC-1BP-mean.pdf", 8, 3)
guidesLegend <- guides(fill=guide_legend(title="Mean AUC"), colour=guide_legend(title="Mean AUC"))
ggplot(resd.ModType.Mean[which(resd.ModType.Mean$variable == "auc" & resd.ModType.Mean$periodFit == 1),]) +
      geom_tile(aes(x=periodPrj, y=sampleSize, fill=value, colour=value)) +
      scale_fill_gradient2(low="darkorange", high="darkgreen", mid="white", midpoint=0.75, limits=c(0.5, 1)) +
      scale_colour_gradient(low="gray90", high="black", limits=c(0.5, 1)) +
      facet_grid(. ~ modelType) +
      themeOpt + guidesLegend + axesLabels + aspRatio
dev.off()

# JACCARD
pdf("Results/Graphs/JACC-1BP-mean.pdf", 8, 3)
guidesLegend <- guides(fill=guide_legend(title="Mean \nJaccard index"), colour=guide_legend(title="Mean \nJaccard index"))
ggplot(rese.ModType.Mean[which(rese.ModType.Mean$variable %in% c("jacc.MEAN.prob") & rese.ModType.Mean$periodFit == 1),]) +
      geom_tile(aes(x=periodPrj, y=sampleSize, fill=value, colour=value)) +
      scale_fill_gradient2(low="darkgreen", high="darkorange", mid="white", midpoint=0.5, limits=c(0, 1)) +
      scale_colour_gradient2(low="black", mid="gray90", high="black", limits=c(0, 1)) +
      facet_grid(. ~ modelType) +
      themeOpt + guidesLegend + axesLabels + aspRatio
dev.off()

# SPECIES RICHNESS
pdf("Results/Graphs/SPRICH-1BP-mean.pdf", 8, 3)
guidesLegend <- guides(fill=guide_legend(title="Mean species \nrichness\ncorrelation"), colour=guide_legend(title="Mean species \nrichness\ncorrelation"))
ggplot(rese.ModType.Mean[which(rese.ModType.Mean$variable %in% c("spRich.COR.prob") & rese.ModType.Mean$periodFit == 1),]) +
      geom_tile(aes(x=periodPrj, y=sampleSize, fill=value, colour=value)) +
      scale_fill_gradient2(low="darkorange", high="darkgreen", mid="white", midpoint=0, limits=c(-1, 1)) +
      scale_colour_gradient2(low="black", mid="gray90", high="black", limits=c(-1, 1)) +
      facet_grid(. ~ modelType) +
      themeOpt + guidesLegend + axesLabels + aspRatio
dev.off()


########
# 9 BP

# AUC
pdf("Results/Graphs/AUC-9BP-mean.pdf", 8, 3)
guidesLegend <- guides(fill=guide_legend(title="Mean AUC"), colour=guide_legend(title="Mean AUC"))
ggplot(resd.ModType.Mean[which(resd.ModType.Mean$variable == "auc" & resd.ModType.Mean$periodFit == 9),]) +
      geom_tile(aes(x=periodPrj, y=sampleSize, fill=value, colour=value)) +
      scale_fill_gradient2(low="darkorange", high="darkgreen", mid="white", midpoint=0.75, limits=c(0.5, 1)) +
      scale_colour_gradient(low="gray90", high="black", limits=c(0.5, 1)) +
      facet_grid(. ~ modelType) +
      themeOpt + guidesLegend + axesLabels + aspRatio
dev.off()

# JACCARD
pdf("Results/Graphs/JACC-9BP-mean.pdf", 8, 3)
guidesLegend <- guides(fill=guide_legend(title="Mean \nJaccard index"), colour=guide_legend(title="Mean \nJaccard index"))
ggplot(rese.ModType.Mean[which(rese.ModType.Mean$variable %in% c("jacc.MEAN.prob") & rese.ModType.Mean$periodFit == 9),]) +
      geom_tile(aes(x=periodPrj, y=sampleSize, fill=value, colour=value)) +
      scale_fill_gradient2(low="darkgreen", high="darkorange", mid="white", midpoint=0.5, limits=c(0, 1)) +
      scale_colour_gradient2(low="black", mid="gray90", high="black", limits=c(0, 1)) +
      facet_grid(. ~ modelType) +
      themeOpt + guidesLegend + axesLabels + aspRatio
dev.off()

# SPECIES RICHNESS
pdf("Results/Graphs/SPRICH-9BP-mean.pdf", 8, 3)
guidesLegend <- guides(fill=guide_legend(title="Mean species \nrichness\ncorrelation"), colour=guide_legend(title="Mean species \nrichness\ncorrelation"))
ggplot(rese.ModType.Mean[which(rese.ModType.Mean$variable %in% c("spRich.COR.prob") & rese.ModType.Mean$periodFit == 9),]) +
      geom_tile(aes(x=periodPrj, y=sampleSize, fill=value, colour=value)) +
      scale_fill_gradient2(low="darkorange", high="darkgreen", mid="white", midpoint=0, limits=c(-1, 1)) +
      scale_colour_gradient2(low="black", mid="gray90", high="black", limits=c(-1, 1)) +
      facet_grid(. ~ modelType) +
      themeOpt + guidesLegend + axesLabels + aspRatio
dev.off()


################################################################################
##
## DIFFERENCES
##
################################################################################

guidesLegend <- guides(fill=guide_legend(title="Difference"), size=guide_legend(title="Significance"))

f_labeller <- function(var, value){
  value <- as.character(value)
  if (var=="variable") { 
    value <- "Diff."
  }
  return(value)
}

########
# 1 BP

# AUC
pdf("Results/Graphs/AUCDif-1BP-mean.pdf", 5.25, 3)
ggplot(resd.ModType.Dif2[which(resd.ModType.Dif2$variable %in% c("auc") & resd.ModType.Dif2$periodFit == 1),]) +
      geom_tile(aes(x=periodPrj, y=sampleSize, fill=dif), colour="gray90") +                             
      geom_tile(aes(x=periodPrj, y=sampleSize, size=significant2), colour="black", fill=NA) +
      scale_fill_gradient2(limits=c(-0.06, 0.06), low="#B2182B", high="#2166AC") +
      scale_size_discrete(range=c(0, 1)) +
      facet_grid(.~ variable, labeller=f_labeller) +
      themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

# JACCARD
pdf("Results/Graphs/jaccDif-1BP-mean.pdf", 5.25, 3)
ggplot(rese.ModType.Dif2[which(rese.ModType.Dif2$variable %in% c("jacc.MEAN.prob") & rese.ModType.Dif2$periodFit == 1),]) +
      geom_tile(aes(x=periodPrj, y=sampleSize, fill=dif), colour="gray90") +
      geom_tile(aes(x=periodPrj, y=sampleSize, size=significant2), colour="black", fill=NA) +
      scale_fill_gradient2(limits=c(-0.1, 0.1), low="#B2182B", high="#2166AC") +
      scale_size_discrete(range=c(0, 1)) +
      facet_grid(.~ variable, labeller=f_labeller) +
      themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

# SPECIES RICHNESS
pdf("Results/Graphs/spRichDif-1BP-mean.pdf", 5.25, 3)
ggplot(rese.ModType.Dif2[which(rese.ModType.Dif2$variable %in% c("spRich.COR.prob") & rese.ModType.Dif2$periodFit == 1),]) +
      geom_tile(aes(x=periodPrj, y=sampleSize, fill=dif), colour="gray90") +
      geom_tile(aes(x=periodPrj, y=sampleSize, size=significant2), colour="black", fill=NA) +
      scale_fill_gradient2(limits=c(-0.3, 0.3), low="#B2182B", high="#2166AC") +
      scale_size_discrete(range=c(0, 1)) +
      facet_grid(.~ variable, labeller=f_labeller) +
      themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

########
# 9 BP

# AUC
pdf("Results/Graphs/AUCDif-9BP-mean.pdf", 5.25, 3)
ggplot(resd.ModType.Dif2[which(resd.ModType.Dif2$variable %in% c("auc") & resd.ModType.Dif2$periodFit == 9),]) +
      geom_tile(aes(x=periodPrj, y=sampleSize, fill=dif), colour="gray90") +                             
      geom_tile(aes(x=periodPrj, y=sampleSize, size=significant2), colour="black", fill=NA) +
      scale_fill_gradient2(limits=c(-0.06, 0.06), low="#B2182B", high="#2166AC") +
      scale_size_discrete(range=c(0, 1)) +
      facet_grid(.~ variable, labeller=f_labeller) +
      themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

# JACCARD
pdf("Results/Graphs/jaccDif-9BP-mean.pdf", 5.25, 3)
ggplot(rese.ModType.Dif2[which(rese.ModType.Dif2$variable %in% c("jacc.MEAN.prob") & rese.ModType.Dif2$periodFit == 9),]) +
      geom_tile(aes(x=periodPrj, y=sampleSize, fill=dif), colour="gray90") +
      geom_tile(aes(x=periodPrj, y=sampleSize, size=significant2), colour="black", fill=NA) +
      scale_fill_gradient2(limits=c(-0.08, 0.08), low="#B2182B", high="#2166AC") +
      scale_size_discrete(range=c(0, 1)) +
      facet_grid(.~ variable, labeller=f_labeller) +
      themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

# SPECIES RICHNESS
pdf("Results/Graphs/spRichDif-9BP-mean.pdf", 5.25, 3)
ggplot(rese.ModType.Dif2[which(rese.ModType.Dif2$variable %in% c("spRich.COR.prob") & rese.ModType.Dif2$periodFit == 9),]) +
      geom_tile(aes(x=periodPrj, y=sampleSize, fill=dif), colour="gray90") +
      geom_tile(aes(x=periodPrj, y=sampleSize, size=significant2), colour="black", fill=NA) +
      scale_fill_gradient2(limits=c(-0.3, 0.3), low="#B2182B", high="#2166AC") +
      scale_size_discrete(range=c(0, 1)) +
      facet_grid(.~ variable, labeller=f_labeller) +
      themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()



################
##
## Line plots
##
################
f_labeller <- function(variable, value){
    value <- as.character(value)
    value <- paste(value, " BP", sep="")
    return(value)
}

##AUC
pdf("Results/Graphs/linePlots/auc-1BP.pdf", 17, 2)
ggplot(resd.test[which(resd.test$variable == "auc" & resd.test$periodFit == 1),], aes(y=value, x=sampleSize, colour=modelType, group=modelType)) + 
  stat_summary(fun.data="mean_cl_boot", geom="pointrange") +
  stat_summary(fun.data="mean_cl_boot", geom="line") +
  facet_grid(. ~ periodPrj, labeller=f_labeller) + labs(x="Sample size", y="AUC", colour="Legend") + 
    theme_bw()
dev.off()
pdf("Results/Graphs/linePlots/auc-9BP.pdf", 17, 2)
ggplot(resd.test[which(resd.test$variable == "auc" & resd.test$periodFit == 9),], aes(y=value, x=sampleSize, colour=modelType, group=modelType)) + 
  stat_summary(fun.data="mean_cl_boot", geom="pointrange") +
  stat_summary(fun.data="mean_cl_boot", geom="line") +
  facet_grid(. ~ periodPrj, labeller=f_labeller) + labs(x="Sample size", y="AUC", colour="Legend") + 
    theme_bw()    
dev.off()

##BRIER
pdf("Results/Graphs/linePlots/brier-1BP.pdf", 17, 2)
ggplot(resd.test[which(resd.test$variable == "brier" & resd.test$periodFit == 1),], aes(y=1-value, x=sampleSize, colour=modelType, group=modelType)) + 
  stat_summary(fun.data="mean_cl_boot", geom="pointrange") +
  stat_summary(fun.data="mean_cl_boot", geom="line") +
  facet_grid(. ~ periodPrj, labeller=f_labeller) + labs(x="Sample size", y="1 - Brier score", colour="Legend") + 
    theme_bw()    
dev.off()
pdf("Results/Graphs/linePlots/brier-9BP.pdf", 17, 2)
ggplot(resd.test[which(resd.test$variable == "brier" & resd.test$periodFit == 9),], aes(y=1-value, x=sampleSize, colour=modelType, group=modelType)) + 
  stat_summary(fun.data="mean_cl_boot", geom="pointrange") +
  stat_summary(fun.data="mean_cl_boot", geom="line") +
  facet_grid(. ~ periodPrj, labeller=f_labeller) + labs(x="Sample size", y="1 - Brier score", colour="Legend") + 
    theme_bw()    
dev.off()

pdf("Results/Graphs/linePlots/jacc_mean_prob-1BP.pdf", 17, 2)
ggplot(rese.test[which(rese.test$variable == "jacc.MEAN.prob" & rese.test$periodFit == 1),], aes(y=1-value, x=sampleSize, colour=modelType, group=modelType)) + 
  stat_summary(fun.data="mean_cl_boot", geom="pointrange") +
  stat_summary(fun.data="mean_cl_boot", geom="line") +
  facet_grid(. ~ periodPrj, labeller=f_labeller) + labs(x="Sample size", y="Bray Curtis\n simmilarity", colour="Legend") + 
    theme_bw()    
dev.off()
pdf("Results/Graphs/linePlots/jacc_mean_prob-9BP.pdf", 17, 2)
ggplot(rese.test[which(rese.test$variable == "jacc.MEAN.prob" & rese.test$periodFit == 9),], aes(y=1-value, x=sampleSize, colour=modelType, group=modelType)) + 
  stat_summary(fun.data="mean_cl_boot", geom="pointrange") +
  stat_summary(fun.data="mean_cl_boot", geom="line") +
  facet_grid(. ~ periodPrj, labeller=f_labeller) + labs(x="Sample size", y="Bray Curtis\n simmilarity", colour="Legend") + 
    theme_bw()    
dev.off()

pdf("Results/Graphs/linePlots/spRich_cor-1BP.pdf", 17, 2)
ggplot(rese.test[which(rese.test$variable == "spRich.COR.prob" & rese.test$periodFit == 1),], aes(y=value, x=sampleSize, colour=modelType, group=modelType)) + 
  stat_summary(fun.data="mean_cl_boot", geom="pointrange") +
  stat_summary(fun.data="mean_cl_boot", geom="line") +
  facet_grid(. ~ periodPrj, labeller=f_labeller) + labs(x="Sample size", y="Species Richness \n correlation", colour="Legend") + 
    theme_bw()    
dev.off()
pdf("Results/Graphs/linePlots/spRich_cor-9BP.pdf", 17, 2)
ggplot(rese.test[which(rese.test$variable == "spRich.COR.prob" & rese.test$periodFit == 9),], aes(y=value, x=sampleSize, colour=modelType, group=modelType)) + 
  stat_summary(fun.data="mean_cl_boot", geom="pointrange") +
  stat_summary(fun.data="mean_cl_boot", geom="line") +
  facet_grid(. ~ periodPrj, labeller=f_labeller) + labs(x="Sample size", y="Species Richness \n correlation", colour="Legend") + 
  theme_bw()
dev.off()

