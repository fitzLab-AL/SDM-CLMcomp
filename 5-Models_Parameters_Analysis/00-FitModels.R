################################################################################
#
# LOAD DATA, PARAMETERS, AND FUNCTIONS
#
################################################################################

# Load the project library with functions to load and arrange data
library(paleoCLMs)

# Load wrapping functions
source("SDMvsCLM/Working_Code/VarPar/00-Wrapping.R")
source("SDMvsCLM/Working_Code/VarPar/00-StepFunctions.R")

# Define parameters
indir <- "SDMvsCLM/PaleoCLM/Data/"

periods <- seq(0, 21000, by=500) # That generate a numerical vector with time slides. I use this afterwards to name files and objects
periods_red <- periods/1000
index <- 1:length(periods)

vars <- c("etr_year_ave","aet_year_ave","wdei_year_ave","tmax_high_quart","prcp_low_quart","prcp_high_quart") # This select the name of the variables that will be used to fit the models

numIter <- 10 # Here I specified how many iteractions will be runned. In this case I fitted 10 models with different random train-test datasets.

# Load pollen data and remove low quality data
pollenOriginal <- lapply(periods, loadPollen, indir)
pollenHQ <- lapply(pollenOriginal, removeLowQualitySamples, 0.75)

# Load a raster as template, it will be deleted a bit later.
rasTemp <- raster("SDMvsCLM/PaleoCLM/Data/Climate/CCSM/0BP/gdd5_year_ave.tif")

# Average values for those pollen sites in the same grid cell
pollenAbun <- lapply(pollenHQ, reduceDuplicated, rasTemp, weighted=F)
rm(rasTemp, pollenOriginal, pollenHQ)

# Load climate data
clim <- mapply(loadClim, period=periods, pollen=pollenAbun, MoreArgs=list(clim_model="CCSM", indir=paste(indir, "Climate-IceSheet/", sep=""), vars=vars))

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
require(foreach)
require(doParallel)
registerDoParallel(cores=4)

seed <- c(1, 11, 21, 31, 41)
seed.i=1

# Only testing this for 0BP so
i=1
j=index_fit[[i]]
k=periods_fit[[i]]

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

# Run step functions

glmList <- glmStep(pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], numIter)
cqoList <- cqoStep(pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], numIter)
gamList <- gamStep(pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], numIter)
caoList <- caoStep(pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], numIter)
smarsList <- smarsStep(pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], numIter)
marsList <- marsStep(pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], numIter)
scartsList <- scartsStep(pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], numIter)
cartsList <- cartsStep(pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], numIter)
snnList <- snnStep(data.frame(pollenFit[[i]]), climFit[[i]], trainSetsFit[[i]], numIter)
mnnList <- mnnStep(pollenFit[[i]], climFit[[i]], trainSetsFit[[i]], numIter)

save(glmList, cqoList, gamList, caoList, scartsList, cartsList, smarsList, marsList, snnList, mnnList, file="SDMvsCLM/R_objects/VaryingResults/Lists.RData")

################################################################################
#
#  EVALUATE THE MODELS
#
################################################################################
require(foreach)
require(doParallel)
registerDoParallel(cores=5)

seed <- c(1, 11, 21, 31, 41)
seed.i=1

# Only testing this for 0BP so
i=1
j=index_fit[[i]]
k=periods_fit[[i]]

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

# Evaluate the models
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

save(glmDistTest,glmDistTrain,glmEnseTest,glmEnseTrain,gamDistTest,gamDistTrain,gamEnseTest,gamEnseTrain,cqoDistTest,cqoDistTrain,cqoEnseTest,cqoEnseTrain,caoDistTest,caoDistTrain,caoEnseTrain,caoEnseTest,smarsDistTrain,smarsDistTest,smarsEnseTrain,smarsEnseTest,marsDistTrain,marsDistTest,marsEnseTrain,marsEnseTest,scartsDistTrain,scartsDistTest,scartsEnseTrain,scartsEnseTest,cartsDistTrain,cartsDistTest,cartsEnseTrain,cartsEnseTest,snnDistTrain,snnDistTest,snnEnseTrain,snnEnseTest,mnnDistTrain,mnnDistTest,mnnEnseTrain,mnnEnseTest,file="SDMvsCLM/R_objects/VaryingResults/DistEnseResults.RData")

################################################################################
#
# FORMAT RESULTS
#
################################################################################

objDTrainNames <- c("glmDistTrain", "cqoDistTrain", "gamDistTrain", "caoDistTrain", "scartsDistTrain", "cartsDistTrain", "smarsDistTrain", "marsDistTrain", "snnDistTrain", "mnnDistTrain")
objDTestNames <- c("glmDistTest", "cqoDistTest", "gamDistTest", "caoDistTest", "scartsDistTest", "cartsDistTest", "smarsDistTest", "marsDistTest", "snnDistTest", "mnnDistTest")
objETrainNames <- c("glmEnseTrain", "cqoEnseTrain", "gamEnseTrain", "caoEnseTrain", "scartsEnseTrain", "cartsEnseTrain", "smarsEnseTrain", "marsEnseTrain", "snnEnseTrain", "mnnEnseTrain")
objETestNames <- c("glmEnseTest", "cqoEnseTest", "gamEnseTest", "caoEnseTest", "scartsEnseTest", "cartsEnseTest", "smarsEnseTest", "marsEnseTest", "snnEnseTest", "mnnEnseTest")

objNames <- c(objDTrainNames, objDTestNames, objETrainNames, objETestNames)

  for(v in 1:length(objDTrainNames)){
    assign("tmpName", get(objDTrainNames[v])) 
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
    assign("tmpName", get(objDTestNames[v])) 
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
    assign("tmpName", get(objETrainNames[v])) 
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
    assign("tmpName", get(objETestNames[v])) 
    ind <- which(unlist(lapply(lapply(tmpName, class), function(x){any(x == "error")})))
    if(length(ind) > 0){
      for(w in ind){
        tmpName[[w]] <- matrix(-9999, nrow=12, ncol=12, dimnames=list(c("jacc.MEAN.prob","jacc.SD.prob","jacc.CAL.prob","jacc.DIS.prob","spRich.RSQ.prob","spRich.COR.prob","jacc.MEAN.thre","jacc.SD.thre","jacc.CAL.thre","jacc.DIS.thre","spRich.RSQ.thre","spRich.COR.thre"), periods_prj))
      }
    }
    assign(objETestNames[v], tmpName) 
  }
  
  
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
  resd.All <- rbind(resDt, resD)
  rese.All <- rbind(resEt, resE)
  
  rm(list=ls(pattern="DistTest"))
  rm(list=ls(pattern="DistTrain"))
  rm(list=ls(pattern="EnseTest"))
  rm(list=ls(pattern="EnseTrain"))


rm(i, j, k, resDt, resEt, resD, resE)
rm(objNames, objDTestNames, objDTrainNames, objETestNames, objETrainNames)

# Write data to file
save(resd.All, file="SDMvsCLM/R_objects/VaryingResults/DisResults.RData")
save(rese.All, file="SDMvsCLM/R_objects/VaryingResults/AssResults.RData")


################################################################################
##
## ARRANGING DATA
##
################################################################################

##Load plotting libraries
library(reshape2)
library(ggplot2)

##Load data
load("SDMvsCLM/R_objects/VaryingResults/DisResults.RData")
load("SDMvsCLM/R_objects/VaryingResults/AssResults.RData")

resd.all <- resd.All
rese.all <- rese.All

##Remove the missing data 
resd.all$value[which(resd.all$value == -9999)] <- NA
rese.all$value[which(rese.all$value == -9999)] <- NA

##Create new column for modely type
resd.all$modelType <- NA
rese.all$modelType <- NA
resd.all$modelType[which(resd.all$model %in% c("glm","gam","scarts","smars","snn"))] <- "SDM"
rese.all$modelType[which(rese.all$model %in% c("glm","gam","scarts","smars","snn"))] <- "SDM"
resd.all$modelType[which(resd.all$model %in% c("cqo","cao","carts","mars","mnn"))] <- "CLM"
rese.all$modelType[which(rese.all$model %in% c("cqo","cao","carts","mars","mnn"))] <- "CLM"

#Specify the fitting periods as a factor, that way they will be plotted in order afterwards
resd.all$model <- factor(resd.all$model, c("glm","gam","scarts","smars","snn","cqo","cao","carts","mars","mnn"))
rese.all$model <- factor(rese.all$model, c("glm","gam","scarts","smars","snn","cqo","cao","carts","mars","mnn"))             
resd.all$periodFit <- factor(resd.all$periodFit, rev(periods_fit))
rese.all$periodFit <- factor(rese.all$periodFit, rev(periods_fit))             
resd.all$periodPrj <- factor(resd.all$periodPrj)
rese.all$periodPrj <- factor(rese.all$periodPrj)

resd.all$modelType <- factor(resd.all$modelType, levels=c("SDM","CLM"))
rese.all$modelType <- factor(rese.all$modelType, levels=c("SDM","CLM"))

resd.all$variable <- factor(resd.all$variable, c("auc","brier","tss","sens","spec","prev","trmin"))
rese.all$variable <- factor(rese.all$variable, c("jacc.MEAN.prob","jacc.SD.prob","jacc.DIS.prob","jacc.CAL.prob","spRich.COR.prob","spRich.RSQ.prob","jacc.MEAN.thre","jacc.SD.thre","jacc.DIS.thre","jacc.CAL.thre","spRich.COR.thre","spRich.RSQ.thre"))

resd.all$test <- factor(resd.all$test, levels=c("train","test"))
rese.all$test <- factor(rese.all$test, levels=c("train","test"))

resd.train <- resd.all[which(as.character(resd.all$periodFit) == as.character(resd.all$periodPrj)),]
rese.train <- rese.all[which(as.character(rese.all$periodFit) == as.character(rese.all$periodPrj)),]

resd.test <- resd.all[which(resd.all$test == "test"),]
rese.test <- rese.all[which(rese.all$test == "test"),]

resd.test$periodPrj <- factor(resd.test$periodPrj, levels=rev(periods_prj))
rese.test$periodPrj <- factor(rese.test$periodPrj, levels=rev(periods_prj))

#################################################################################
###
### SIGNIFICANCE TEST BETWEEN SDM AND CLM
###
#################################################################################
### function for calculating t-test for each metric
sdmStats <- function(dat, p_prj, p_fit, measure){
  #dat <- rese.test
  #p_prj <- periods_prj
  #p_fit <- periods_fit
  #measure <- "jacc.DIS.thre"
  
  stat.glm <- stat.gam <- stat.scarts <- stat.smars <- stat.snnet <- matrix(NA, length(p_prj), length(p_fit))
  dat1 <- dat[which(dat$variable == measure),]
  for (l in 1:length(p_fit)){
    i <- p_fit[l]
    for (j in 1:length(p_prj)){
      k <- p_prj[j]
      
      print(paste("Working on", i))
      print(paste("  and", k))
      
      dat2 <- dat1[which(dat1$periodPrj == k & dat1$periodFit == i),]
      
      if(length(na.omit(unique(dat2$value[which(dat2$model == "gam")]))) < 3 | length(na.omit(unique(dat2$value[which(dat2$model == "cao")]))) < 3){
        stat.gam[j,l] <- NA
      }else{
        sdm <- dat2[which(dat2$model == "gam"),]
        clm <- dat2[which(dat2$model == "cao"),]
        results <- t.test(sdm$value, clm$value, paired=T)
        stat.gam[j,l] <- results$p.value
      }     
      
      if(length(na.omit(unique(dat2$value[which(dat2$model == "glm")]))) < 3 | length(na.omit(unique(dat2$value[which(dat2$model == "cqo")]))) < 3){
        stat.gam[j,l] <- NA
      }else{
        sdm <- dat2[which(dat2$model == "glm"),]
        clm <- dat2[which(dat2$model == "cqo"),]
        results <- t.test(sdm$value, clm$value, paired=T)
        stat.glm[j,l] <- results$p.value
      }
      
      sdm <- dat2[which(dat2$model == "scarts"),]
      clm <- dat2[which(dat2$model == "carts"),]
      results <- t.test(sdm$value, clm$value, paired=T)
      stat.scarts[j,l] <- results$p.value
      
      sdm <- dat2[which(dat2$model == "smars"),]
      clm <- dat2[which(dat2$model == "mars"),]
      results <- t.test(sdm$value, clm$value, paired=T)
      stat.smars[j,l] <- results$p.value
      
      sdm <- dat2[which(dat2$model == "snn"),]
      clm <- dat2[which(dat2$model == "mnn"),]
      results <- t.test(sdm$value, clm$value, paired=T)
      stat.snnet[j,l] <- results$p.value
    }
  }
  
  rownames(stat.glm) <- rownames(stat.gam) <- rownames(stat.scarts) <- rownames(stat.smars) <- rownames(stat.snnet) <- p_prj
  colnames(stat.glm) <- colnames(stat.gam) <- colnames(stat.scarts) <- colnames(stat.smars) <- colnames(stat.snnet) <- p_fit
  
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
  colnames(stat.Data) <- c("periodPrj","periodFit","p_value","models")
  stat.Data$models <- factor(stat.Data$models)
  return(stat.Data)
}



## apply the function for each metric and save the results
auc.ttest <- sdmStats(resd.test, periods_prj, periods_fit, "auc")
auc.ttest <- auc.ttest[complete.cases(auc.ttest),]
save(auc.ttest, file="SDMvsCLM/R_objects/VaryingResults/Ttest/auc.ttest")

brier.ttest <- sdmStats(resd.test, periods_prj, periods_fit, "brier")
brier.ttest <- brier.ttest[complete.cases(brier.ttest),]
save(brier.ttest, file="SDMvsCLM/R_objects/VaryingResults/Ttest/brier.ttest")

tss.ttest <- sdmStats(resd.test, periods_prj, periods_fit, "tss")
tss.ttest <- tss.ttest[complete.cases(tss.ttest),]
save(tss.ttest, file="SDMvsCLM/R_objects/VaryingResults/Ttest/tss.ttest")

sens.ttest <- sdmStats(resd.test, periods_prj, periods_fit, "sens")
sens.ttest <- sens.ttest[complete.cases(sens.ttest),]
save(sens.ttest, file="SDMvsCLM/R_objects/VaryingResults/Ttest/sens.ttest")

spec.ttest <- sdmStats(resd.test, periods_prj, periods_fit, "spec")
spec.ttest <- spec.ttest[complete.cases(spec.ttest),]
save(spec.ttest, file="SDMvsCLM/R_objects/VaryingResults/Ttest/spec.ttest")

jacc.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "jacc.MEAN.prob")
jacc.ttest <- jacc.ttest[complete.cases(jacc.ttest),]
save(jacc.ttest, file="SDMvsCLM/R_objects/VaryingResults/Ttest/jacc.ttest")

jacDis.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "jacc.DIS.prob")
jacDis.ttest <- jacDis.ttest[complete.cases(jacDis.ttest),]
save(jacDis.ttest, file="SDMvsCLM/R_objects/VaryingResults/Ttest/jacDis.ttest")

jacCal.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "jacc.CAL.prob")
jacCal.ttest <- jacCal.ttest[complete.cases(jacCal.ttest),]
save(jacCal.ttest, file="SDMvsCLM/R_objects/VaryingResults/Ttest/jacCal.ttest")

spRichCor.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "spRich.COR.prob")
spRichCor.ttest <- spRichCor.ttest[complete.cases(spRichCor.ttest),]
save(spRichCor.ttest, file="SDMvsCLM/R_objects/VaryingResults/Ttest/spRichCor.ttest")

spRichRsq.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "spRich.RSQ.prob")
spRichRsq.ttest <- spRichRsq.ttest[complete.cases(spRichRsq.ttest),]
save(spRichRsq.ttest, file="SDMvsCLM/R_objects/VaryingResults/Ttest/spRichRsq.ttest")

jacct.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "jacc.MEAN.thre")
jacct.ttest <- jacct.ttest[complete.cases(jacct.ttest),]
save(jacct.ttest, file="SDMvsCLM/R_objects/VaryingResults/Ttest/jacct.ttest")

jacDist.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "jacc.DIS.thre")
jacDist.ttest <- jacDist.ttest[complete.cases(jacDist.ttest),]
save(jacDist.ttest, file="SDMvsCLM/R_objects/VaryingResults/Ttest/jacDist.ttest")

jacCalt.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "jacc.CAL.thre")
jacCalt.ttest <- jacCalt.ttest[complete.cases(jacCalt.ttest),]
save(jacCalt.ttest, file="SDMvsCLM/R_objects/VaryingResults/Ttest/jacCalt.ttest")

spRichCort.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "spRich.COR.thre")
spRichCort.ttest <- spRichCort.ttest[complete.cases(spRichCort.ttest),]
save(spRichCort.ttest, file="SDMvsCLM/R_objects/VaryingResults/Ttest/spRichCort.ttest")

spRichRsqt.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "spRich.RSQ.thre")
spRichRsqt.ttest <- spRichRsqt.ttest[complete.cases(spRichRsqt.ttest),]
save(spRichRsqt.ttest, file="SDMvsCLM/R_objects/VaryingResults/Ttest/spRichRsqt.ttest")

### Function for testing the difference betweeen SDM and CLM using a t-test
sdmStats_ModelType <- function(dat, p_prj, p_fit, measure){
  stat.Data <- matrix(NA, length(p_prj), length(p_fit))
  dat1 <- dat[which(dat$variable == measure),]
  
  l=1
    for(j in 1:length(p_prj)){
      i <- p_fit[l]
      k <- p_prj[j]
      
      print(paste("Working on", i))
      print(paste("  and", k))
      
      dat2 <- dat1[which(dat1$periodPrj == k & dat1$periodFit == i),]
      
      sdm <- dat2[which(dat2$modelType == "SDM"),]
      clm <- dat2[which(dat2$modelType == "CLM"),]
      
      ind <- intersect(unique(sdm$iteration), unique(clm$iteration))
      sdm <- sdm[which(sdm$iteration %in% ind),]
      clm <- clm[which(clm$iteration %in% ind),]
      results <- t.test(sdm$value, clm$value, paired=T)
      stat.Data[j,l] <- results$p.value
    }
  
  rownames(stat.Data) <- p_prj
  colnames(stat.Data) <- p_fit
  stat.Data <- melt(stat.Data)
  colnames(stat.Data) <- c("periodPrj","periodFit","p_value")
  stat.Data$variable <- measure
  return(stat.Data)
}

## Apply the function for each metric
auc.modelType <- sdmStats_ModelType(resd.test, periods_prj, periods_fit, "auc")
auc.modelType <- auc.modelType[complete.cases(auc.modelType),]

brier.modelType <- sdmStats_ModelType(resd.test, periods_prj, periods_fit, "brier")
brier.modelType <- brier.modelType[complete.cases(brier.modelType),]

tss.modelType <- sdmStats_ModelType(resd.test, periods_prj, periods_fit, "tss")
tss.modelType <- tss.modelType[complete.cases(tss.modelType),]

sens.modelType <- sdmStats_ModelType(resd.test, periods_prj, periods_fit, "sens")
sens.modelType <- sens.modelType[complete.cases(sens.modelType),]

spec.modelType <- sdmStats_ModelType(resd.test, periods_prj, periods_fit, "spec")
spec.modelType <- spec.modelType[complete.cases(spec.modelType),]


jacc.modelType <- sdmStats_ModelType(rese.test, periods_prj, periods_fit, "jacc.MEAN.prob")
jacc.modelType <- jacc.modelType[complete.cases(jacc.modelType),]

jacCal.modelType <- sdmStats_ModelType(rese.test, periods_prj, periods_fit, "jacc.CAL.prob")
jacCal.modelType <- jacCal.modelType[complete.cases(jacCal.modelType),]

jacDis.modelType <- sdmStats_ModelType(rese.test, periods_prj, periods_fit, "jacc.DIS.prob")
jacDis.modelType <- jacDis.modelType[complete.cases(jacDis.modelType),]

spRichCor.modelType <- sdmStats_ModelType(rese.test, periods_prj, periods_fit, "spRich.COR.prob")
spRichCor.modelType <- spRichCor.modelType[complete.cases(spRichCor.modelType),]

spRichRsq.modelType <- sdmStats_ModelType(rese.test, periods_prj, periods_fit, "spRich.RSQ.prob")
spRichRsq.modelType <- spRichRsq.modelType[complete.cases(spRichRsq.modelType),]

jacct.modelType <- sdmStats_ModelType(rese.test, periods_prj, periods_fit, "jacc.MEAN.thre")
jacct.modelType <- jacct.modelType[complete.cases(jacct.modelType),]

jacCalt.modelType <- sdmStats_ModelType(rese.test, periods_prj, periods_fit, "jacc.CAL.thre")
jacCalt.modelType <- jacCalt.modelType[complete.cases(jacCalt.modelType),]

jacDist.modelType <- sdmStats_ModelType(rese.test, periods_prj, periods_fit, "jacc.DIS.thre")
jacDist.modelType <- jacDist.modelType[complete.cases(jacDist.modelType),]

spRichCort.modelType <- sdmStats_ModelType(rese.test, periods_prj, periods_fit, "spRich.COR.thre")
spRichCort.modelType <- spRichCort.modelType[complete.cases(spRichCort.modelType),]

spRichRsqt.modelType <- sdmStats_ModelType(rese.test, periods_prj, periods_fit, "spRich.RSQ.thre")
spRichRsqt.modelType <- spRichRsqt.modelType[complete.cases(spRichRsqt.modelType),]

### Combined results and save
stats.modelType <- rbind(auc.modelType, tss.modelType, sens.modelType, spec.modelType, brier.modelType, jacc.modelType, jacCal.modelType, jacDis.modelType, spRichCor.modelType, spRichRsq.modelType, jacct.modelType, jacCalt.modelType, jacDist.modelType, spRichCort.modelType, spRichRsqt.modelType)

save(stats.modelType, file="SDMvsCLM/R_objects/VaryingResults/Ttest/stats.modelType.R")


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
index.dColumns <- which(colnames(resd.test) %in% c("iteration","taxon","test","modelType"))
index.eColumns <- which(colnames(rese.test) %in% c("iteration","test","modelType"))

resd.Models.Mean <- acast(resd.test[,-index.dColumns], periodPrj ~ periodFit ~ variable ~ model , mean, na.rm=T)
rese.Models.Mean <- acast(rese.test[,-index.eColumns], periodPrj ~ periodFit ~ variable ~ model, mean, na.rm=T)

index.dColumns <- which(colnames(resd.test) %in% c("iteration","taxon","test","model"))
index.eColumns <- which(colnames(rese.test) %in% c("iteration","test","model"))

resd.ModType.Mean <- acast(resd.test[,-index.dColumns], periodPrj ~ periodFit ~ variable ~ modelType, mean, na.rm=T)
rese.ModType.Mean <- acast(rese.test[,-index.eColumns], periodPrj ~ periodFit ~ variable ~ modelType, mean, na.rm=T)

## Melting again the new data (means and sd) to get the data in shape for the ggplot functions.
resd.Models.Mean <- melt(resd.Models.Mean)
rese.Models.Mean <- melt(rese.Models.Mean)

resd.ModType.Mean <- melt(resd.ModType.Mean)
rese.ModType.Mean <- melt(rese.ModType.Mean)


## Change colnames of the new dataframes
colnames(resd.Models.Mean) <- colnames(rese.Models.Mean) <- c("periodPrj","periodFit","variable","model","value")

colnames(resd.ModType.Mean) <- colnames(rese.ModType.Mean) <- c("periodPrj","periodFit","variable","modelType","value")


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
calDif <- function(d){
  if(ncol(d) == 7){
    p <- merge(d, d, by=c("periodPrj","periodFit","sdm_eq","variable"))
  }
  if(ncol(d) == 5){
    p <- merge(d, d, by=c("periodPrj","periodFit","variable"))
  }
  p <- p[which(p$modelType.x == "SDM"),]
  p <- p[which(p$modelType.y == "CLM"),]
  
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
# Calculate and include sample size of the fitting period
#
################################################################################


# Pooling data together
# Binning 14-15
binIndex <- which(periods %in% seq(14000, 15000,by=500))
pbin1 <- poolingData(binIndex, pollen)

# Binning 15.5 to 17.5
binIndex <- which(periods %in% seq(15500, 17500, by=500))
pbin2 <- poolingData(binIndex, pollen)

# Binning 15.5 to 21
binIndex <- which(periods %in% seq(15500, 21000, by=500))
pbin3 <- poolingData(binIndex, pollen)

# Binning 18 to 21
binIndex <- which(periods %in% seq(18000, 21000, by=500))
pbin4 <- poolingData(binIndex, pollen)

# Appending the pooled dataframes to the data lists
pollenPool <- append(pollen, list(pbin1, pbin2, pbin3, pbin4))
pollenPool <- lapply(pollenPool, FUN=function(x){rownames(x) <- 1:nrow(x); return(x)})


# Change names
names(pollenPool) <- c(periods_red, "14-15", "15.5-17.5", "15.5-21", "18-21")

# Selecting and splitting data for fitting and projecting time periods
pollenFit <- pollenPool[which(names(pollenPool) %in% periods_fit)]

sample.size <- sapply(pollenFit, nrow)
sample.size <- data.frame(periodFit=periods_fit, sampleSize=sample.size)

# Remove temporal objects
rm(binIndex, pbin1, pbin2, pbin3, pbin4, pollenPool, pollenFit)

resd.Models.Mean <- merge(resd.Models.Mean, sample.size, by="periodFit")
rese.Models.Mean <- merge(rese.Models.Mean, sample.size, by="periodFit")
resd.Models.Dif <- merge(resd.Models.Dif, sample.size, by="periodFit")
rese.Models.Dif <- merge(rese.Models.Dif, sample.size, by="periodFit")

resd.ModType.Mean <- merge(resd.ModType.Mean, sample.size, by="periodFit")
rese.ModType.Mean <- merge(rese.ModType.Mean, sample.size, by="periodFit")
resd.ModType.Dif <- merge(resd.ModType.Dif, sample.size, by="periodFit")
rese.ModType.Dif <- merge(rese.ModType.Dif, sample.size, by="periodFit")

################################################################################
#
# Merge Difference and p-value datasets for plotting purposes
#
################################################################################

resd.Models.Dif2 <- merge(resd.Models.Dif, disStats.Models, by=c("periodPrj", "periodFit", "model.x", "model.y", "variable"))
rese.Models.Dif2 <- merge(rese.Models.Dif, assStats.Models, by=c("periodPrj", "periodFit", "model.x", "model.y", "variable"))

resd.ModType.Dif2 <- merge(resd.ModType.Dif, stats.ModType, by=c("periodPrj", "periodFit", "variable"))
rese.ModType.Dif2 <- merge(rese.ModType.Dif, stats.ModType, by=c("periodPrj", "periodFit", "variable"))


save(resd.Models.Dif2, file="SDMvsCLM/R_objects/VaryingResults/resd.Models.Dif2.RData")
save(rese.Models.Dif2, file="SDMvsCLM/R_objects/VaryingResults/rese.Models.Dif2.RData")
save(resd.ModType.Dif2, file="SDMvsCLM/R_objects/VaryingResults/resd.ModType.Dif2.RData")
save(rese.ModType.Dif2, file="SDMvsCLM/R_objects/VaryingResults/rese.ModType.Dif2.RData")


################################################################################
#
# Plotting Results
#
################################################################################

##Load plotting libraries
library(reshape2)
library(ggplot2)
library(scales)


# These plotting functions can be modified to reflect the fixed scales (min and max values) and use more attractive colors.

# fitLine <- geom_line(data=forecastLine, aes(x=x, y=y), colour="black", size=1.3, alpha=0.75)

themeOpt <- theme(axis.text.x=element_text(angle=90, hjust=1)) + theme_bw()

axesLabels <- labs(x="Projecting period", y="Fitting period")

aspRatio <- coord_fixed(ratio=1)



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

# DIFFERENCES IN MEAN DISTRIBUTION VALUES
pdf("SDMvsCLM/R_objects/VaryingResults/Graphs/acrossTime-Dif/auc-Dif.pdf", 4.5, 4.85)
ggplot(resd.ModType.Dif2[which(resd.ModType.Dif2$variable %in% c("auc")),]) +
  geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
  geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
  scale_fill_gradient2(limits=c(-0.08, 0.08), low="#B2182B", high="#2166AC") +
  scale_size_discrete(range=c(0, 1)) +
  facet_grid(.~ variable, labeller=f_labeller) +
  themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("SDMvsCLM/R_objects/VaryingResults/Graphs/acrossTime-Dif/sens-Dif.pdf", 4.5, 4.85)
ggplot(resd.ModType.Dif2[which(resd.ModType.Dif2$variable %in% c("sens")),]) +
  geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
  geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
  scale_fill_gradient2(limits=c(-0.11, 0.11), low="#B2182B", high="#2166AC") +
  scale_size_discrete(range=c(0, 1)) +
  facet_grid(.~ variable, labeller=f_labeller) +
  themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("SDMvsCLM/R_objects/VaryingResults/Graphs/acrossTime-Dif/spec-Dif.pdf", 4.5, 4.85)
ggplot(resd.ModType.Dif2[which(resd.ModType.Dif2$variable %in% c("spec")),]) +
  geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
  geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
  scale_fill_gradient2(limits=c(-0.12, 0.12), low="#B2182B", high="#2166AC") +
  scale_size_discrete(range=c(0, 1)) +
  facet_grid(.~ variable, labeller=f_labeller) +
  themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("SDMvsCLM/R_objects/VaryingResults/Graphs/acrossTime-Dif/brier-Dif.pdf", 4.5, 4.85)
ggplot(resd.ModType.Dif2[which(resd.ModType.Dif2$variable %in% c("brier")),]) +
  geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
  geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
  scale_fill_gradient2(limits=c(-0.1, 0.1), low="#B2182B", high="#2166AC") +
  scale_size_discrete(range=c(0, 1)) +
  facet_grid(.~ variable, labeller=f_labeller) +
  themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()


# DIFFERENCES IN MEAN ASSEMBLAGE VALUES
pdf("SDMvsCLM/R_objects/VaryingResults/Graphs/acrossTime-Dif/jac-thre-Dif.pdf", 4.5, 4.85)
ggplot(rese.ModType.Dif2[which(rese.ModType.Dif2$variable %in% c("jacc.MEAN.thre")),]) +
  geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
  geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
  scale_fill_gradient2(limits=c(-0.1, 0.1), low="#B2182B", high="#2166AC") +
  scale_size_discrete(range=c(0, 1)) +
  facet_grid(.~ variable, labeller=f_labeller) +
  themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("SDMvsCLM/R_objects/VaryingResults/Graphs/acrossTime-Dif/jac-prob-Dif.pdf", 4.5, 4.85)
ggplot(rese.ModType.Dif2[which(rese.ModType.Dif2$variable %in% c("jacc.MEAN.prob")),]) +
  geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
  geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
  scale_fill_gradient2(limits=c(-0.1, 0.1), low="#B2182B", high="#2166AC") +
  scale_size_discrete(range=c(0, 1)) +
  facet_grid(.~ variable, labeller=f_labeller) +
  themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("SDMvsCLM/R_objects/VaryingResults/Graphs/acrossTime-Dif/jac-cal-thre-Dif.pdf", 4.5, 4.85)
ggplot(rese.ModType.Dif2[which(rese.ModType.Dif2$variable %in% c("jacc.CAL.thre")),]) +
  geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
  geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
  scale_fill_gradient2(limits=c(-6, 6), low="#B2182B", high="#2166AC") +
  scale_size_discrete(range=c(0, 1)) +
  facet_grid(.~ variable, labeller=f_labeller) +
  themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("SDMvsCLM/R_objects/VaryingResults/Graphs/acrossTime-Dif/jac-dis-thre-Dif.pdf", 4.5, 4.85)
ggplot(rese.ModType.Dif2[which(rese.ModType.Dif2$variable %in% c("jacc.DIS.thre")),]) +
  geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
  geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
  scale_fill_gradient2(limits=c(-0.15, 0.15), low="#B2182B", high="#2166AC") +
  scale_size_discrete(range=c(0, 1)) +
  facet_grid(.~ variable, labeller=f_labeller) +
  themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("SDMvsCLM/R_objects/VaryingResults/Graphs/acrossTime-Dif/jac-cal-prob-Dif.pdf", 4.5, 4.85)
ggplot(rese.ModType.Dif2[which(rese.ModType.Dif2$variable %in% c("jacc.CAL.prob")),]) +
  geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
  geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
  scale_fill_gradient2(limits=c(-6, 6), low="#B2182B", high="#2166AC") +
  scale_size_discrete(range=c(0, 1)) +
  facet_grid(.~ variable, labeller=f_labeller) +
  themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("SDMvsCLM/R_objects/VaryingResults/Graphs/acrossTime-Dif/jac-dis-prob-Dif.pdf", 4.5, 4.85)
ggplot(rese.ModType.Dif2[which(rese.ModType.Dif2$variable %in% c("jacc.DIS.prob")),]) +
  geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
  geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
  scale_fill_gradient2(limits=c(-0.15, 0.15), low="#B2182B", high="#2166AC") +
  scale_size_discrete(range=c(0, 1)) +
  facet_grid(.~ variable, labeller=f_labeller) +
  themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("SDMvsCLM/R_objects/VaryingResults/Graphs/acrossTime-Dif/spRich-dis-Dif.pdf", 4.5, 4.85)
ggplot(rese.ModType.Dif2[which(rese.ModType.Dif2$variable %in% c("spRich.COR.thre")),]) +
  geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
  geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
  scale_fill_gradient2(limits=c(-0.3, 0.3), low="#B2182B", high="#2166AC") +
  scale_size_discrete(range=c(0, 1)) +
  facet_grid(.~ variable, labeller=f_labeller) +
  themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("SDMvsCLM/R_objects/VaryingResults/Graphs/acrossTime-Dif/spRich-cal-Dif.pdf", 4.5, 4.85)
ggplot(rese.ModType.Dif2[which(rese.ModType.Dif2$variable %in% c("spRich.RSQ.thre")),]) +
  geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
  geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
  scale_fill_gradient2(limits=c(-1.5, 1.5), low="#B2182B", high="#2166AC") +
  scale_size_discrete(range=c(0, 1)) +
  facet_grid(.~ variable, labeller=f_labeller) +
  themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()


# DIFFERENCES BY MODELS IN MEAN DISTRIBUTION VALUES
f_labeller <- function(var, value){
  value <- as.character(value)
  if (var=="variable") { 
    value <- paste("Diff:", value, sep="")
  }
  return(value)
}

pdf("SDMvsCLM/R_objects/VaryingResults/Graphs/acrossTime-Dif/auc-Dif-models.pdf", 12.8, 8.5)
ggplot(resd.Models.Dif2[which(resd.Models.Dif2$variable %in% c("auc")),]) +
  geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
  geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
  scale_fill_gradient2(limits=c(-0.26, 0.26), midpoint=0, low="#B2182B", high="#2166AC") +
  scale_size_discrete(range=c(0, 1)) +
  facet_grid(. ~ sdm_eq) +
  themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("SDMvsCLM/R_objects/VaryingResults/Graphs/acrossTime-Dif/sens-Dif-models.pdf", 12.8, 8.5)
ggplot(resd.Models.Dif2[which(resd.Models.Dif2$variable %in% c("sens")),]) +
  geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
  geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
  scale_fill_gradient2(limits=c(-0.3, 0.3), midpoint=0, low="#B2182B", high="#2166AC") +
  scale_size_discrete(range=c(0, 1)) +
  facet_grid(. ~ sdm_eq) +
  themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("SDMvsCLM/R_objects/VaryingResults/Graphs/acrossTime-Dif/spec-Dif-models.pdf", 12.8, 8.5)
ggplot(resd.Models.Dif2[which(resd.Models.Dif2$variable %in% c("spec")),]) +
  geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
  geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
  scale_fill_gradient2(limits=c(-0.31, 0.31), midpoint=0, low="#B2182B", high="#2166AC") +
  scale_size_discrete(range=c(0, 1)) +
  facet_grid(. ~ sdm_eq) +
  themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("SDMvsCLM/R_objects/VaryingResults/Graphs/acrossTime-Dif/brier-Dif-models.pdf", 12.8, 8.5)
ggplot(resd.Models.Dif2[which(resd.Models.Dif2$variable %in% c("brier")),]) +
  geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
  geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
  scale_fill_gradient2(limits=c(-0.2, 0.2), midpoint=0, low="#B2182B", high="#2166AC") +
  scale_size_discrete(range=c(0, 1)) +
  facet_grid(. ~ sdm_eq) +
  themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()


pdf("SDMvsCLM/R_objects/VaryingResults/Graphs/acrossTime-Dif/jac-thre-Dif-models.pdf", 12.8, 12.5)
ggplot(rese.Models.Dif2[which(rese.Models.Dif2$variable == "jacc.MEAN.thre"),]) +
  geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
  geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
  scale_fill_gradient2(limits=c(-0.2, 0.2), midpoint=0, low="#B2182B", high="#2166AC") +
  scale_size_discrete(range=c(0, 1)) + facet_grid(. ~ sdm_eq) +
  themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("SDMvsCLM/R_objects/VaryingResults/Graphs/acrossTime-Dif/jac-prob-Dif-models.pdf", 12.8, 12.5)
ggplot(rese.Models.Dif2[which(rese.Models.Dif2$variable == "jacc.MEAN.prob"),]) +
  geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
  geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
  scale_fill_gradient2(limits=c(-0.2, 0.2), midpoint=0, low="#B2182B", high="#2166AC") +
  scale_size_discrete(range=c(0, 1)) + facet_grid(. ~ sdm_eq) +
  themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("SDMvsCLM/R_objects/VaryingResults/Graphs/acrossTime-Dif/jac-cal-thre-Dif-models.pdf", 12.8, 12.5)
ggplot(rese.Models.Dif2[which(rese.Models.Dif2$variable == "jacc.CAL.thre"),]) +
  geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
  geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
  scale_fill_gradient2(limits=c(-20, 20), midpoint=0, low="#B2182B", high="#2166AC") +
  scale_size_discrete(range=c(0, 1)) + facet_grid(. ~ sdm_eq) +
  themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("SDMvsCLM/R_objects/VaryingResults/Graphs/acrossTime-Dif/jac-dis-thre-Dif-models.pdf", 12.8, 12.5)
ggplot(rese.Models.Dif2[which(rese.Models.Dif2$variable == "jacc.DIS.thre"),]) +
  geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
  geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
  scale_fill_gradient2(limits=c(-0.3, 0.3), midpoint=0, low="#B2182B", high="#2166AC") +
  scale_size_discrete(range=c(0, 1)) + facet_grid(. ~ sdm_eq) +
  themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("SDMvsCLM/R_objects/VaryingResults/Graphs/acrossTime-Dif/jac-cal-prob-Dif-models.pdf", 12.8, 12.5)
ggplot(rese.Models.Dif2[which(rese.Models.Dif2$variable == "jacc.CAL.prob"),]) +
  geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
  geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
  scale_fill_gradient2(limits=c(-20, 20), midpoint=0, low="#B2182B", high="#2166AC") +
  scale_size_discrete(range=c(0, 1)) + facet_grid(. ~ sdm_eq) +
  themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("SDMvsCLM/R_objects/VaryingResults/Graphs/acrossTime-Dif/jac-dis-prob-Dif-models.pdf", 12.8, 12.5)
ggplot(rese.Models.Dif2[which(rese.Models.Dif2$variable == "jacc.DIS.prob"),]) +
  geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
  geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
  scale_fill_gradient2(limits=c(-0.3, 0.3), midpoint=0, low="#B2182B", high="#2166AC") +
  scale_size_discrete(range=c(0, 1)) + facet_grid(. ~ sdm_eq) +
  themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("SDMvsCLM/R_objects/VaryingResults/Graphs/acrossTime-Dif/spRich-dis-Dif-models.pdf", 12.8, 12.5)
ggplot(rese.Models.Dif2[which(rese.Models.Dif2$variable == "spRich.COR.thre"),]) +
  geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
  geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
  scale_fill_gradient2(limits=c(-0.85, 0.85), midpoint=0, low="#B2182B", high="#2166AC") +
  scale_size_discrete(range=c(0, 1)) +
  facet_grid(. ~ sdm_eq) +
  themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("SDMvsCLM/R_objects/VaryingResults/Graphs/acrossTime-Dif/spRich-cal-Dif-models.pdf", 12.8, 12.5)
ggplot(rese.Models.Dif2[which(rese.Models.Dif2$variable == "spRich.RSQ.thre"),]) +
  geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
  geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
  scale_fill_gradient2(limits=c(-4, 4), midpoint=0, low="#B2182B", high="#2166AC") +
  scale_size_discrete(range=c(0, 1)) +
  facet_grid(. ~ sdm_eq) +
  themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

################################################################################
##
## MEANS
##
################################################################################

themeOpt <- theme(axis.text.x=element_text(angle=90, hjust=1)) + theme_bw()

axesLabels <- labs(x="Projecting period", y="Fitting period")

aspRatio <- coord_fixed(ratio=1)

resd.ModType.Mean$periodPrj <- factor(resd.ModType.Mean$periodPrj, rev(periods_prj))

# MEAN DISTRIBUTION VALUES
pdf("Results/Graphs/acrossTime-Means/auc.pdf", 6.5, 5)
guidesLegend <- guides(fill=guide_legend(title="Mean AUC"), colour=guide_legend(title="Mean AUC"))
ggplot(resd.ModType.Mean[which(resd.ModType.Mean$variable == "auc"),]) +
  geom_tile(aes(x=periodPrj, y=periodFit, fill=value, colour=value)) +
  scale_fill_gradient2(low="darkorange", high="darkgreen", mid="white", midpoint=0.75, limits=c(0.5, 1)) +
  scale_colour_gradient(low="gray90", high="black", limits=c(0.5, 1)) +
  facet_grid(. ~ modelType) +
  themeOpt + guidesLegend + axesLabels + aspRatio
dev.off()

# Specific Model Results

resd.Models.Mean$periodPrj <- factor(resd.Models.Mean$periodPrj, rev(periods_prj))

pdf("Results/Graphs/acrossTime-Means/auc-models.pdf", 12.8, 8.5)
guidesLegend <- guides(fill=guide_legend(title="Mean AUC"), colour=guide_legend(title="Mean AUC"))
ggplot(resd.Models.Mean[which(resd.Models.Mean$variable == "auc"),]) +
  geom_tile(aes(x=periodPrj, y=periodFit, fill=value, colour=value)) +
  scale_fill_gradient2(low="darkorange", high="darkgreen", mid="white", midpoint=0.75, limits=c(0.485, 1)) + 
  scale_colour_gradient(low="gray90", high="black", limits=c(0.485, 1)) +
  facet_grid(modelType ~ sdm_eq) +
  themeOpt + guidesLegend + axesLabels + aspRatio
dev.off()


# Line Plot
pdf("Results/Graphs/linePlots/auc-1BP.pdf", 17, 2)
ggplot(resd.Models.Mean[which(resd.Models.Mean$variable == "auc",], aes(y=value, x=periodPrj, colour=modelType, group=modelType)) + 
  stat_summary(fun.data="mean_cl_boot", geom="pointrange") +
  stat_summary(fun.data="mean_cl_boot", geom="line") +
  theme_bw()
dev.off()