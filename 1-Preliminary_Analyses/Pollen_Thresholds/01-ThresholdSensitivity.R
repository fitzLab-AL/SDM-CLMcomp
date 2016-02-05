################################################################################
#
# LOAD DATA, PARAMETERS, AND FUNCTIONS
#
################################################################################

# Load functions for the project. It will be changed by a R-package
library(paleoCLMs)

# Define parameters
indir <- "SDMvsCLM/PaleoCLM/Data/" # That is a general folder with a folder call "Pollen" and other call "Climate". Both folders have exactly the same data and structure as in Research Hub.
period <- seq(0, 21000, by=500) # That generate a numerical vector with time slides. I use this afterwards to name files and objects
vars <- c("etr_year_ave","aet_year_ave","wdei_year_ave","tmax_high_quart","prcp_low_quart","prcp_high_quart") # This select the name of the variables that will be used to fit the models
numIter <- 10 # Here I specified how many iteractions will be runned. In this case I fitted 10 models with different random train-test datasets.

# Load pollen data and remove low quality data
pollenOriginal <- lapply(period, loadPollen, indir)
pollenHQ <- lapply(pollenOriginal, removeLowQualitySamples, 0.75)

# Load a raster as template, it will be deleted a bit later.
library(raster)
rasTemp <- raster("SDMvsCLM/PaleoCLM/Data/Climate-Icesheet/CCSM/0BP/gdd5_year_ave.tif")

# Average values for those pollen sites in the same grid cell
pollen <- lapply(pollenHQ, reduceDuplicated, rasTemp, weighted=F)
rm(rasTemp)

# Load climate data
clim <- mapply(loadClim, period=period, pollen=pollen, MoreArgs=list(clim_model="CCSM", indir=indir, vars=vars))
clim <- lapply(clim, as.data.frame) # For some of the models I need climate data as data.frame


climLists <- lapply(period, function(x, y){paste(indir, "Climate-IceSheet/CCSM/", x, "BP/", y, ".tif", sep="")}, vars)
climStacks <- lapply(climLists, stack)

# Remove incomplete cases from the pollen and the climate datasets
compCasesIndex <- lapply(clim, FUN=function(x){which(complete.cases(x))})
pollen <- mapply(FUN=function(x,i){x[i,]}, pollen, compCasesIndex, SIMPLIFY=F)
clim <- lapply(clim, FUN=function(x){x[which(complete.cases(x)),]})

# Extract coordinates to a new object and remove them from the pollen data.frame
coord <- lapply(pollen, FUN=function(x){x[,which(colnames(x) %in% c("x","y"))]})
pollen <- lapply(pollen, FUN=function(x){x[,-which(colnames(x) %in% c("x","y"))]})

# Convert pollen concentrations to presences/absences matrix (Community matrix)
# This should be changed when we figure out a species-specific threshold
pollenThres <- read.csv("SDMvsCLM/PaleoCLM/Data/Pollen/PollenThresholds.csv")
pollenThres[,c(4:7)] <- pollenThres[,c(4:7)]/100

pollen1 <- lapply(pollen, applyThreshold, pollenThres[,c(1,4)]) #No threshold
pollen2 <- lapply(pollen, applyThreshold, pollenThres[,c(1,5)]) #Half percent threshold
pollen3 <- lapply(pollen, applyThreshold, pollenThres[,c(1,6)]) #Biomes threshold
pollen4 <- lapply(pollen, applyThreshold, pollenThres[,c(1,7)]) #Intercept threshold

# Calculate variable thresholds
var05 <- lapply(pollen, pollenThreshold, 0.05)
var01 <- lapply(pollen, pollenThreshold, 0.01)
pollen5 <- mapply(applyThreshold, pollen, threshold=var01, SIMPLIFY=F)
pollen6 <- mapply(applyThreshold, pollen, threshold=var05, SIMPLIFY=F)


# Get a taxon list for training the models and remove all species with less than 20 presences and/or absences. Then propagate the taxon list to all time periods
taxList1 <- colnames(removeAbsentSpecies(pollen1[[1]], 19, applyTo="both"))
taxList2 <- colnames(removeAbsentSpecies(pollen2[[1]], 19, applyTo="both"))
taxList3 <- colnames(removeAbsentSpecies(pollen3[[1]], 19, applyTo="both"))
taxList4 <- colnames(removeAbsentSpecies(pollen4[[1]], 19, applyTo="both"))
taxList5 <- colnames(removeAbsentSpecies(pollen5[[1]], 19, applyTo="both"))
taxList6 <- colnames(removeAbsentSpecies(pollen6[[1]], 19, applyTo="both"))

# Combine the taxon list to get only those pollen in all the lists
taxonNames <- Reduce(intersect, list(taxList1, taxList2, taxList3, taxList4,  taxList5, taxList6))

# Propagate the taxon list to all pollen datasets
pollen1 <- lapply(pollen1, FUN=function(x){x[,taxonNames]})
pollen2 <- lapply(pollen2, FUN=function(x){x[,taxonNames]})
pollen3 <- lapply(pollen3, FUN=function(x){x[,taxonNames]})
pollen4 <- lapply(pollen4, FUN=function(x){x[,taxonNames]})
pollen5 <- lapply(pollen5, FUN=function(x){x[,taxonNames]})
pollen6 <- lapply(pollen6, FUN=function(x){x[,taxonNames]})

# Get present data to fit the models
p0BP1 <- pollen1[[1]]
p0BP2 <- pollen2[[1]]
p0BP3 <- pollen3[[1]]
p0BP4 <- pollen4[[1]]
p0BP5 <- pollen5[[1]]
p0BP6 <- pollen6[[1]]

# Get a set of present climate conditions for training the models
c0BP <- clim[[1]]

# Set all the paleo data for testing but the current are splited in 70/30 for training/testing
trainSets <- lapply(pollen1, trainSplit, numIter, 0)
trainSets[[1]] <- trainSplit(p0BP6, numIter, 0.7)



################################################################################
#
# FIT AND EVALUATE MODELS
#
################################################################################

## Loading wrapping functions
source("SDMvsCLM/FINAL/Final_Code/Preliminary_Analyses/Pollen_Thresholds/00-ThreshWrappingFunctions.R")

## Selecting time periods to project the models
projIndex <- which(period %in% c(0, 6000, 12000, 14000, 15000))


################################################################################
## GLM - Different thresholds

require(foreach)
require(doParallel)
registerDoParallel(cores=4) # Check how many cores have your computer and set up a lower number here

# Fitting VGLM models
glmList.1 <- glmFit(p0BP1, c0BP, trainSets[[1]], numIter)
glmList.2 <- glmFit(p0BP2, c0BP, trainSets[[1]], numIter)
glmList.3 <- glmFit(p0BP3, c0BP, trainSets[[1]], numIter)
glmList.4 <- glmFit(p0BP4, c0BP, trainSets[[1]], numIter)
glmList.5 <- glmFit(p0BP5, c0BP, trainSets[[1]], numIter)
glmList.6 <- glmFit(p0BP6, c0BP, trainSets[[1]], numIter)

# Evaluating Distribution on train dataset
glmDistTrain.1 <- evalDistTrain(glmList.1, pollen1[[1]], clim[[1]], trainSets[[1]], respType="response", multiSp=F)
glmDistTrain.2 <- evalDistTrain(glmList.2, pollen2[[1]], clim[[1]], trainSets[[1]], respType="response", multiSp=F)
glmDistTrain.3 <- evalDistTrain(glmList.3, pollen3[[1]], clim[[1]], trainSets[[1]], respType="response", multiSp=F)
glmDistTrain.4 <- evalDistTrain(glmList.4, pollen4[[1]], clim[[1]], trainSets[[1]], respType="response", multiSp=F)
glmDistTrain.5 <- evalDistTrain(glmList.5, pollen5[[1]], clim[[1]], trainSets[[1]], respType="response", multiSp=F)
glmDistTrain.6 <- evalDistTrain(glmList.6, pollen6[[1]], clim[[1]], trainSets[[1]], respType="response", multiSp=F)

# Evaluating Distribution on test dataset
glmDistTest.1 <- evalDistTest(glmList.1, pollen1[projIndex], clim[projIndex], trainSets[projIndex], respType="response", multiSp=F)
glmDistTest.2 <- evalDistTest(glmList.2, pollen2[projIndex], clim[projIndex], trainSets[projIndex], respType="response", multiSp=F)
glmDistTest.3 <- evalDistTest(glmList.3, pollen3[projIndex], clim[projIndex], trainSets[projIndex], respType="response", multiSp=F)
glmDistTest.4 <- evalDistTest(glmList.4, pollen4[projIndex], clim[projIndex], trainSets[projIndex], respType="response", multiSp=F)
glmDistTest.5 <- evalDistTest(glmList.5, pollen5[projIndex], clim[projIndex], trainSets[projIndex], respType="response", multiSp=F)
glmDistTest.6 <- evalDistTest(glmList.6, pollen6[projIndex], clim[projIndex], trainSets[projIndex], respType="response", multiSp=F)
                              
# Evaluating Assamblage on train dataset
glmEnseTrain.1 <- evalEnseTrain(glmList.1, pollen1[[1]], clim[[1]], trainSets[[1]], respType="response", glmDistTrain.1, multiSp=F)
glmEnseTrain.2 <- evalEnseTrain(glmList.2, pollen2[[1]], clim[[1]], trainSets[[1]], respType="response", glmDistTrain.2, multiSp=F)
glmEnseTrain.3 <- evalEnseTrain(glmList.3, pollen3[[1]], clim[[1]], trainSets[[1]], respType="response", glmDistTrain.3, multiSp=F)
glmEnseTrain.4 <- evalEnseTrain(glmList.4, pollen4[[1]], clim[[1]], trainSets[[1]], respType="response", glmDistTrain.4, multiSp=F)
glmEnseTrain.5 <- evalEnseTrain(glmList.5, pollen5[[1]], clim[[1]], trainSets[[1]], respType="response", glmDistTrain.5, multiSp=F)
glmEnseTrain.6 <- evalEnseTrain(glmList.6, pollen6[[1]], clim[[1]], trainSets[[1]], respType="response", glmDistTrain.6, multiSp=F)

# Evaluating Assamblages on test dataset
glmEnseTest.1 <- evalEnseTest(glmList.1, pollen1[projIndex], clim[projIndex], trainSets[projIndex], respType="response", glmDistTrain.1, multiSp=F)
glmEnseTest.2 <- evalEnseTest(glmList.2, pollen2[projIndex], clim[projIndex], trainSets[projIndex], respType="response", glmDistTrain.2, multiSp=F)
glmEnseTest.3 <- evalEnseTest(glmList.3, pollen3[projIndex], clim[projIndex], trainSets[projIndex], respType="response", glmDistTrain.3, multiSp=F)
glmEnseTest.4 <- evalEnseTest(glmList.4, pollen4[projIndex], clim[projIndex], trainSets[projIndex], respType="response", glmDistTrain.4, multiSp=F)
glmEnseTest.5 <- evalEnseTest(glmList.5, pollen5[projIndex], clim[projIndex], trainSets[projIndex], respType="response", glmDistTrain.5, multiSp=F)
glmEnseTest.6 <- evalEnseTest(glmList.6, pollen6[projIndex], clim[projIndex], trainSets[projIndex], respType="response", glmDistTrain.6, multiSp=F)

# Evaluating Distributioni range
library(SDMTools)
glmDistRange.1 <- getDistRange(glmList.1, climStacks[projIndex], glmDistTrain.1, respType="response", multiSp=F) 
glmDistRange.2 <- getDistRange(glmList.2, climStacks[projIndex], glmDistTrain.2, respType="response", multiSp=F) 
glmDistRange.3 <- getDistRange(glmList.3, climStacks[projIndex], glmDistTrain.3, respType="response", multiSp=F) 
glmDistRange.4 <- getDistRange(glmList.4, climStacks[projIndex], glmDistTrain.4, respType="response", multiSp=F) 
glmDistRange.5 <- getDistRange(glmList.5, climStacks[projIndex], glmDistTrain.5, respType="response", multiSp=F) 
glmDistRange.6 <- getDistRange(glmList.6, climStacks[projIndex], glmDistTrain.6, respType="response", multiSp=F) 


################################################################################
## CQO - different threshold

require(foreach)
require(doParallel)
registerDoParallel(cores=4) # Check how many cores have your computer and set up a lower number here

cqoList.1 <- cqoFit(p0BP1, c0BP, trainSets[[1]], numIter)
cqoList.2 <- cqoFit(p0BP2, c0BP, trainSets[[1]], numIter)
cqoList.3 <- cqoFit(p0BP3, c0BP, trainSets[[1]], numIter)
cqoList.4 <- cqoFit(p0BP4, c0BP, trainSets[[1]], numIter)
cqoList.5 <- cqoFit(p0BP5, c0BP, trainSets[[1]], numIter)
cqoList.6 <- cqoFit(p0BP6, c0BP, trainSets[[1]], numIter)

cqoDistTrain.1 <- evalDistTrain(cqoList.1, pollen1[[1]], clim[[1]], trainSets[[1]], respType="response")
cqoDistTrain.2 <- evalDistTrain(cqoList.2, pollen2[[1]], clim[[1]], trainSets[[1]], respType="response")
cqoDistTrain.3 <- evalDistTrain(cqoList.3, pollen3[[1]], clim[[1]], trainSets[[1]], respType="response")
cqoDistTrain.4 <- evalDistTrain(cqoList.4, pollen4[[1]], clim[[1]], trainSets[[1]], respType="response")
cqoDistTrain.5 <- evalDistTrain(cqoList.5, pollen5[[1]], clim[[1]], trainSets[[1]], respType="response")
cqoDistTrain.6 <- evalDistTrain(cqoList.6, pollen6[[1]], clim[[1]], trainSets[[1]], respType="response")

cqoDistTest.1 <- evalDistTest(cqoList.1, pollen1[projIndex], clim[projIndex], trainSets[projIndex], respType="response")
cqoDistTest.2 <- evalDistTest(cqoList.2, pollen2[projIndex], clim[projIndex], trainSets[projIndex], respType="response")
cqoDistTest.3 <- evalDistTest(cqoList.3, pollen3[projIndex], clim[projIndex], trainSets[projIndex], respType="response")
cqoDistTest.4 <- evalDistTest(cqoList.4, pollen4[projIndex], clim[projIndex], trainSets[projIndex], respType="response")
cqoDistTest.5 <- evalDistTest(cqoList.5, pollen5[projIndex], clim[projIndex], trainSets[projIndex], respType="response")
cqoDistTest.6 <- evalDistTest(cqoList.6, pollen6[projIndex], clim[projIndex], trainSets[projIndex], respType="response")

cqoEnseTrain.1 <- evalEnseTrain(cqoList.1, pollen1[[1]], clim[[1]], trainSets[[1]], respType="response", cqoDistTrain.1)
cqoEnseTrain.2 <- evalEnseTrain(cqoList.2, pollen2[[1]], clim[[1]], trainSets[[1]], respType="response", cqoDistTrain.2)
cqoEnseTrain.3 <- evalEnseTrain(cqoList.3, pollen3[[1]], clim[[1]], trainSets[[1]], respType="response", cqoDistTrain.3)
cqoEnseTrain.4 <- evalEnseTrain(cqoList.4, pollen4[[1]], clim[[1]], trainSets[[1]], respType="response", cqoDistTrain.4)
cqoEnseTrain.5 <- evalEnseTrain(cqoList.5, pollen5[[1]], clim[[1]], trainSets[[1]], respType="response", cqoDistTrain.5)
cqoEnseTrain.6 <- evalEnseTrain(cqoList.6, pollen6[[1]], clim[[1]], trainSets[[1]], respType="response", cqoDistTrain.6)

cqoEnseTest.1 <- evalEnseTest(cqoList.1, pollen1[projIndex], clim[projIndex], trainSets[projIndex], respType="response", cqoDistTrain.1)
cqoEnseTest.2 <- evalEnseTest(cqoList.2, pollen2[projIndex], clim[projIndex], trainSets[projIndex], respType="response", cqoDistTrain.2)
cqoEnseTest.3 <- evalEnseTest(cqoList.3, pollen3[projIndex], clim[projIndex], trainSets[projIndex], respType="response", cqoDistTrain.3)
cqoEnseTest.4 <- evalEnseTest(cqoList.4, pollen4[projIndex], clim[projIndex], trainSets[projIndex], respType="response", cqoDistTrain.4)
cqoEnseTest.5 <- evalEnseTest(cqoList.5, pollen5[projIndex], clim[projIndex], trainSets[projIndex], respType="response", cqoDistTrain.5)
cqoEnseTest.6 <- evalEnseTest(cqoList.6, pollen6[projIndex], clim[projIndex], trainSets[projIndex], respType="response", cqoDistTrain.6)

cqoDistRange.1 <- getDistRange(cqoList.1, climStacks[projIndex], cqoDistTrain.1, respType="response") 
cqoDistRange.2 <- getDistRange(cqoList.2, climStacks[projIndex], cqoDistTrain.2, respType="response") 
cqoDistRange.3 <- getDistRange(cqoList.3, climStacks[projIndex], cqoDistTrain.3, respType="response") 
cqoDistRange.4 <- getDistRange(cqoList.4, climStacks[projIndex], cqoDistTrain.4, respType="response") 
cqoDistRange.5 <- getDistRange(cqoList.5, climStacks[projIndex], cqoDistTrain.5, respType="response") 
cqoDistRange.6 <- getDistRange(cqoList.6, climStacks[projIndex], cqoDistTrain.6, respType="response") 


################################################################################
## NNET - different threshold
require(foreach)
require(doParallel)
registerDoParallel(cores=4) # Check how many cores have your computer and set up a lower number here

mnnList.1 <- mnnFit(p0BP1, c0BP, trainSets[[1]], numIter)
mnnList.2 <- mnnFit(p0BP2, c0BP, trainSets[[1]], numIter)
mnnList.3 <- mnnFit(p0BP3, c0BP, trainSets[[1]], numIter)
mnnList.4 <- mnnFit(p0BP4, c0BP, trainSets[[1]], numIter)
mnnList.5 <- mnnFit(p0BP5, c0BP, trainSets[[1]], numIter)
mnnList.6 <- mnnFit(p0BP6, c0BP, trainSets[[1]], numIter)

mnnDistTrain.1 <- evalDistTrain(mnnList.1, pollen1[[1]], clim[[1]], trainSets[[1]], respType="raw")
mnnDistTrain.2 <- evalDistTrain(mnnList.2, pollen2[[1]], clim[[1]], trainSets[[1]], respType="raw")
mnnDistTrain.3 <- evalDistTrain(mnnList.3, pollen3[[1]], clim[[1]], trainSets[[1]], respType="raw")
mnnDistTrain.4 <- evalDistTrain(mnnList.4, pollen4[[1]], clim[[1]], trainSets[[1]], respType="raw")
mnnDistTrain.5 <- evalDistTrain(mnnList.5, pollen5[[1]], clim[[1]], trainSets[[1]], respType="raw")
mnnDistTrain.6 <- evalDistTrain(mnnList.6, pollen6[[1]], clim[[1]], trainSets[[1]], respType="raw")

mnnDistTest.1 <- evalDistTest(mnnList.1, pollen1[projIndex], clim[projIndex], trainSets[projIndex], respType="raw")
mnnDistTest.2 <- evalDistTest(mnnList.2, pollen2[projIndex], clim[projIndex], trainSets[projIndex], respType="raw")
mnnDistTest.3 <- evalDistTest(mnnList.3, pollen3[projIndex], clim[projIndex], trainSets[projIndex], respType="raw")
mnnDistTest.4 <- evalDistTest(mnnList.4, pollen4[projIndex], clim[projIndex], trainSets[projIndex], respType="raw")
mnnDistTest.5 <- evalDistTest(mnnList.5, pollen5[projIndex], clim[projIndex], trainSets[projIndex], respType="raw")
mnnDistTest.6 <- evalDistTest(mnnList.6, pollen6[projIndex], clim[projIndex], trainSets[projIndex], respType="raw")

mnnEnseTrain.1 <- evalEnseTrain(mnnList.1, pollen1[[1]], clim[[1]], trainSets[[1]], respType="raw", mnnDistTrain.1)
mnnEnseTrain.2 <- evalEnseTrain(mnnList.2, pollen2[[1]], clim[[1]], trainSets[[1]], respType="raw", mnnDistTrain.2)
mnnEnseTrain.3 <- evalEnseTrain(mnnList.3, pollen3[[1]], clim[[1]], trainSets[[1]], respType="raw", mnnDistTrain.3)
mnnEnseTrain.4 <- evalEnseTrain(mnnList.4, pollen4[[1]], clim[[1]], trainSets[[1]], respType="raw", mnnDistTrain.4)
mnnEnseTrain.5 <- evalEnseTrain(mnnList.5, pollen5[[1]], clim[[1]], trainSets[[1]], respType="raw", mnnDistTrain.5)
mnnEnseTrain.6 <- evalEnseTrain(mnnList.6, pollen6[[1]], clim[[1]], trainSets[[1]], respType="raw", mnnDistTrain.6)

mnnEnseTest.1 <- evalEnseTest(mnnList.1, pollen1[projIndex], clim[projIndex], trainSets[projIndex], respType="raw", mnnDistTrain.1)
mnnEnseTest.2 <- evalEnseTest(mnnList.2, pollen2[projIndex], clim[projIndex], trainSets[projIndex], respType="raw", mnnDistTrain.2)
mnnEnseTest.3 <- evalEnseTest(mnnList.3, pollen3[projIndex], clim[projIndex], trainSets[projIndex], respType="raw", mnnDistTrain.3)
mnnEnseTest.4 <- evalEnseTest(mnnList.4, pollen4[projIndex], clim[projIndex], trainSets[projIndex], respType="raw", mnnDistTrain.4)
mnnEnseTest.5 <- evalEnseTest(mnnList.5, pollen5[projIndex], clim[projIndex], trainSets[projIndex], respType="raw", mnnDistTrain.5)
mnnEnseTest.6 <- evalEnseTest(mnnList.6, pollen6[projIndex], clim[projIndex], trainSets[projIndex], respType="raw", mnnDistTrain.6)

mnnDistRange.1 <- getDistRange(mnnList.1, climStacks[projIndex], mnnDistTrain.1, respType="raw") 
mnnDistRange.2 <- getDistRange(mnnList.2, climStacks[projIndex], mnnDistTrain.2, respType="raw") 
mnnDistRange.3 <- getDistRange(mnnList.3, climStacks[projIndex], mnnDistTrain.3, respType="raw") 
mnnDistRange.4 <- getDistRange(mnnList.4, climStacks[projIndex], mnnDistTrain.4, respType="raw") 
mnnDistRange.5 <- getDistRange(mnnList.5, climStacks[projIndex], mnnDistTrain.5, respType="raw") 
mnnDistRange.6 <- getDistRange(mnnList.6, climStacks[projIndex], mnnDistTrain.6, respType="raw") 


################################################################################
## NNET - Single species - different threshold

snnList.1 <- snnFit(data.frame(p0BP1), c0BP, trainSets[[1]], numIter)
snnList.2 <- snnFit(data.frame(p0BP2), c0BP, trainSets[[1]], numIter)
snnList.3 <- snnFit(data.frame(p0BP3), c0BP, trainSets[[1]], numIter)
snnList.4 <- snnFit(data.frame(p0BP4), c0BP, trainSets[[1]], numIter)
snnList.5 <- snnFit(data.frame(p0BP5), c0BP, trainSets[[1]], numIter)
snnList.6 <- snnFit(data.frame(p0BP6), c0BP, trainSets[[1]], numIter)

snnDistTrain.1 <- evalDistTrain(snnList.1, pollen1[[1]], clim[[1]], trainSets[[1]], respType="raw", multiSp=F)
snnDistTrain.2 <- evalDistTrain(snnList.2, pollen2[[1]], clim[[1]], trainSets[[1]], respType="raw", multiSp=F)
snnDistTrain.3 <- evalDistTrain(snnList.3, pollen3[[1]], clim[[1]], trainSets[[1]], respType="raw", multiSp=F)
snnDistTrain.4 <- evalDistTrain(snnList.4, pollen4[[1]], clim[[1]], trainSets[[1]], respType="raw", multiSp=F)
snnDistTrain.5 <- evalDistTrain(snnList.5, pollen5[[1]], clim[[1]], trainSets[[1]], respType="raw", multiSp=F)
snnDistTrain.6 <- evalDistTrain(snnList.6, pollen6[[1]], clim[[1]], trainSets[[1]], respType="raw", multiSp=F)

snnDistTest.1 <- evalDistTest(snnList.1, pollen1[projIndex], clim[projIndex], trainSets[projIndex], respType="raw", multiSp=F)
snnDistTest.2 <- evalDistTest(snnList.2, pollen2[projIndex], clim[projIndex], trainSets[projIndex], respType="raw", multiSp=F)
snnDistTest.3 <- evalDistTest(snnList.3, pollen3[projIndex], clim[projIndex], trainSets[projIndex], respType="raw", multiSp=F)
snnDistTest.4 <- evalDistTest(snnList.4, pollen4[projIndex], clim[projIndex], trainSets[projIndex], respType="raw", multiSp=F)
snnDistTest.5 <- evalDistTest(snnList.5, pollen5[projIndex], clim[projIndex], trainSets[projIndex], respType="raw", multiSp=F)
snnDistTest.6 <- evalDistTest(snnList.6, pollen6[projIndex], clim[projIndex], trainSets[projIndex], respType="raw", multiSp=F)

snnEnseTrain.1 <- evalEnseTrain(snnList.1, pollen1[[1]], clim[[1]], trainSets[[1]], respType="raw", snnDistTrain.1, multiSp=F)
snnEnseTrain.2 <- evalEnseTrain(snnList.2, pollen2[[1]], clim[[1]], trainSets[[1]], respType="raw", snnDistTrain.2, multiSp=F)
snnEnseTrain.3 <- evalEnseTrain(snnList.3, pollen3[[1]], clim[[1]], trainSets[[1]], respType="raw", snnDistTrain.3, multiSp=F)
snnEnseTrain.4 <- evalEnseTrain(snnList.4, pollen4[[1]], clim[[1]], trainSets[[1]], respType="raw", snnDistTrain.4, multiSp=F)
snnEnseTrain.5 <- evalEnseTrain(snnList.5, pollen5[[1]], clim[[1]], trainSets[[1]], respType="raw", snnDistTrain.5, multiSp=F)
snnEnseTrain.6 <- evalEnseTrain(snnList.6, pollen6[[1]], clim[[1]], trainSets[[1]], respType="raw", snnDistTrain.6, multiSp=F)

snnEnseTest.1 <- evalEnseTest(snnList.1, pollen1[projIndex], clim[projIndex], trainSets[projIndex], respType="raw", snnDistTrain.1, multiSp=F)
snnEnseTest.2 <- evalEnseTest(snnList.2, pollen2[projIndex], clim[projIndex], trainSets[projIndex], respType="raw", snnDistTrain.2, multiSp=F)
snnEnseTest.3 <- evalEnseTest(snnList.3, pollen3[projIndex], clim[projIndex], trainSets[projIndex], respType="raw", snnDistTrain.3, multiSp=F)
snnEnseTest.4 <- evalEnseTest(snnList.4, pollen4[projIndex], clim[projIndex], trainSets[projIndex], respType="raw", snnDistTrain.4, multiSp=F)
snnEnseTest.5 <- evalEnseTest(snnList.5, pollen5[projIndex], clim[projIndex], trainSets[projIndex], respType="raw", snnDistTrain.5, multiSp=F)
snnEnseTest.6 <- evalEnseTest(snnList.6, pollen6[projIndex], clim[projIndex], trainSets[projIndex], respType="raw", snnDistTrain.6, multiSp=F)

snnDistRange.1 <- getDistRange(snnList.1, climStacks[projIndex], snnDistTrain.1, respType="raw", multiSp=F) 
snnDistRange.2 <- getDistRange(snnList.2, climStacks[projIndex], snnDistTrain.2, respType="raw", multiSp=F) 
snnDistRange.3 <- getDistRange(snnList.3, climStacks[projIndex], snnDistTrain.3, respType="raw", multiSp=F) 
snnDistRange.4 <- getDistRange(snnList.4, climStacks[projIndex], snnDistTrain.4, respType="raw", multiSp=F) 
snnDistRange.5 <- getDistRange(snnList.5, climStacks[projIndex], snnDistTrain.5, respType="raw", multiSp=F) 
snnDistRange.6 <- getDistRange(snnList.6, climStacks[projIndex], snnDistTrain.6, respType="raw", multiSp=F) 



################################################################################
#
# FUNCTIONS TO MELT DATA LISTS
#
################################################################################
 
############################
# DISTRIBUTION - TRAIN - OBP
meltDistTrain <- function(dataList, modelName){
  require(reshape2)

  # Melting data lists
  data.1 <- melt(dataList[[1]])
  data.2 <- melt(dataList[[2]])
  data.3 <- melt(dataList[[3]])
  data.4 <- melt(dataList[[4]])
  data.5 <- melt(dataList[[5]])
  data.6 <- melt(dataList[[6]])

  # Assign column names to the melted dataframe  
  colnames(data.1) <- colnames(data.2) <- colnames(data.3) <- colnames(data.4) <- colnames(data.5) <- colnames(data.6) <- c("taxon","variable","value","iteration")
  
  # Combining the different dataframes for each pollen threshold
  datadt <- rbind(data.1, data.2, data.3, data.4, data.5, data.6)
  
  # Add new column with information on the threshold used for each data
  datadt <- cbind(datadt, thres=c(rep("NO_THRES", nrow(data.1)), rep("HALF_P", nrow(data.2)), rep("BIOMES", nrow(data.3)), rep("INTERC", nrow(data.4)), rep("VAR_01", nrow(data.5)), rep("VAR_05", nrow(data.6))))
  
  # New columns for period, train, and the model
  datadt$period <- 0
  datadt$test <- "train"
  datadt$model <- modelName
  return(datadt)
}


##########################
# ASSEMBLAGE - TRAIN - OBP
meltEnseTrain <- function(dataList, modelName){

  require(reshape2)

  # Melting data lists
  data.1 <- melt(dataList[[1]])
  data.2 <- melt(dataList[[2]])
  data.3 <- melt(dataList[[3]])
  data.4 <- melt(dataList[[4]])
  data.5 <- melt(dataList[[5]])
  data.6 <- melt(dataList[[6]])
  
  # Assign column names to the melted dataframe  
  colnames(data.1) <- colnames(data.2) <- colnames(data.3) <- colnames(data.4) <- colnames(data.5) <- colnames(data.6) <- c("value","iteration")
  
  # Combining the different dataframes for each pollen threshold
  dataet <- rbind(data.1, data.2, data.3, data.4, data.5, data.6)

  # Add column with the name of the variable (AUC, TSS, etc)
  dataet <- cbind(variable=rep(names(dataList[[1]][[1]]), nrow(dataet)/length(names(dataList[[1]][[1]]))), dataet)
  
  #Add new column with information on the threshold used for each data
  dataet <- cbind(dataet, thres=c(rep("NO_THRES", nrow(data.1)), rep("HALF_P", nrow(data.2)), rep("BIOMES", nrow(data.3)), rep("INTERC", nrow(data.4)), rep("VAR_01", nrow(data.5)), rep("VAR_05", nrow(data.6))))
  
  # New columns for period, train, and the model
  dataet$period <- 0
  dataet$test <- "train"
  dataet$model <- modelName
  return(dataet)
}


####################################
# DISTRIBUTION RANGE - SEVERAL TIMES
meltDistRange <- function(dataList, modelName){

  require(reshape2)

  # Melting data lists
  data.1 <- melt(dataList[[1]])
  data.2 <- melt(dataList[[2]])
  data.3 <- melt(dataList[[3]])
  data.4 <- melt(dataList[[4]])
  data.5 <- melt(dataList[[5]])
  data.6 <- melt(dataList[[6]])
  
  # Assign column names to the melted dataframe  
  colnames(data.1) <- colnames(data.2) <- colnames(data.3) <- colnames(data.4) <- colnames(data.5) <- colnames(data.6) <- c("variable","value","species","period")
  
  # Combining the different dataframes for each pollen threshold
  datadr <- rbind(data.1, data.2, data.3, data.4, data.5, data.6)
  
  # Add new column with information on the threshold used for each data
  datadr <- cbind(datadr, thres=c(rep("NO_THRES", nrow(data.1)), rep("HALF_P", nrow(data.2)), rep("BIOMES", nrow(data.3)), rep("INTERC", nrow(data.4)), rep("VAR_01", nrow(data.5)), rep("VAR_05", nrow(data.6))))
  
  # New columns for period, train, and the model
  datadr$period <- as.factor(datadr$period)
  datadr$model <- modelName
  datadr$test <- "test"
  return(datadr)
}


#####################################
# DISTRIBUTION - TEST - SEVERAL TIMES
meltDistTest <- function(dataList, modelName){

  require(reshape2)

  # Melting data lists
  data.1 <- melt(dataList[[1]])
  data.2 <- melt(dataList[[2]])
  data.3 <- melt(dataList[[3]])
  data.4 <- melt(dataList[[4]])
  data.5 <- melt(dataList[[5]])
  data.6 <- melt(dataList[[6]])
  
  # Assign column names to the melted dataframe  
  colnames(data.1) <- colnames(data.2) <- colnames(data.3) <- colnames(data.4) <- colnames(data.5) <- colnames(data.6) <- c("taxon","variable","value","period","iteration")
  
  # Combining the different dataframes for each pollen threshold
  datadt <- rbind(data.1, data.2, data.3, data.4, data.5, data.6)
  
  # Add new column with information on the threshold used for each data
  datadt <- cbind(datadt, thres=c(rep("NO_THRES", nrow(data.1)), rep("HALF_P", nrow(data.2)), rep("BIOMES", nrow(data.3)), rep("INTERC", nrow(data.4)), rep("VAR_01", nrow(data.5)), rep("VAR_05", nrow(data.6))))

  # New columns for train, and the model
  datadt$test <- "test"
  datadt$model <- modelName
  return(datadt)
}


###################################
# ASSEMBLAGE - TEST - SEVERAL TIMES
meltEnseTest <- function(dataList, modelName){

  require(reshape2)

  # Melting data lists
  data.1 <- melt(dataList[[1]])
  data.2 <- melt(dataList[[2]])
  data.3 <- melt(dataList[[3]])
  data.4 <- melt(dataList[[4]])
  data.5 <- melt(dataList[[5]])
  data.6 <- melt(dataList[[6]])

  # Assign column names to the melted dataframe  
  colnames(data.1) <- colnames(data.2) <- colnames(data.3) <- colnames(data.4) <- colnames(data.5) <- colnames(data.6) <- c("variable","period","value","iteration")
  
  # Combining the different dataframes for each pollen threshold
  dataet <- rbind(data.1, data.2, data.3, data.4, data.5, data.6)

  # Add new column with information on the threshold used for each data
  colnames(dataet) <- c("variable","period","value","iteration")
  dataet <- cbind(dataet, thres=c(rep("NO_THRES", nrow(data.1)), rep("HALF_P", nrow(data.2)), rep("BIOMES", nrow(data.3)), rep("INTERC", nrow(data.4)), rep("VAR_01", nrow(data.5)), rep("VAR_05", nrow(data.6))))

  # New columns for train, and the model
  dataet$test <- "test"
  dataet$model <- modelName
  return(dataet)
}


##################################################
# REDUCED MAJOR AXIS (RMA) regression coefficients 
rmaFit <- function(dat, a, b){
  require(lmodel2)
  form <- as.formula(paste(b, "~", a, sep=""))
  mod <- lmodel2(form, data=dat, "interval", "interval", 1)
  reg <- mod$regression.results
  names(reg) <- c("method", "intercept", "slope", "angle", "p-value")
  reg$Rsquare <- mod$rsquare
  return(reg[4,])
}



################################################################################
#
# FORMAT RESULTS
#
################################################################################

#####################
# DISTRIBUTION TRAIN
glmdt <- meltDistTrain(list(glmDistTrain.1, glmDistTrain.2, glmDistTrain.3, glmDistTrain.4, glmDistTrain.5, glmDistTrain.6), "glm")
cqodt <- meltDistTrain(list(cqoDistTrain.1, cqoDistTrain.2, cqoDistTrain.3, cqoDistTrain.4, cqoDistTrain.5, cqoDistTrain.6), "cqo")
mnndt <- meltDistTrain(list(mnnDistTrain.1, mnnDistTrain.2, mnnDistTrain.3, mnnDistTrain.4, mnnDistTrain.5, mnnDistTrain.6), "mnn")
snndt <- meltDistTrain(list(snnDistTrain.1, snnDistTrain.2, snnDistTrain.3, snnDistTrain.4, snnDistTrain.5, snnDistTrain.6), "snn")

disTrain <- rbind(glmdt, cqodt, mnndt, snndt)
rm(glmdt, cqodt, mnndt, snndt)

disTrain2 <- merge(disTrain, disTrain, by=c("taxon","variable","iteration","period","test","model"))
disTrain2 <- disTrain2[which(disTrain2$thres.x != disTrain2$thres.y),]


###################
# ASSEMBLAGE TRAIN
glmet <- meltEnseTrain(list(glmEnseTrain.1, glmEnseTrain.2, glmEnseTrain.3, glmEnseTrain.4, glmEnseTrain.5, glmEnseTrain.6), "glm")
cqoet <- meltEnseTrain(list(cqoEnseTrain.1, cqoEnseTrain.2, cqoEnseTrain.3, cqoEnseTrain.4, cqoEnseTrain.5, cqoEnseTrain.6), "cqo")
mnnet <- meltEnseTrain(list(mnnEnseTrain.1, mnnEnseTrain.2, mnnEnseTrain.3, mnnEnseTrain.4, mnnEnseTrain.5, mnnEnseTrain.6), "mnn")
snnet <- meltEnseTrain(list(snnEnseTrain.1, snnEnseTrain.2, snnEnseTrain.3, snnEnseTrain.4, snnEnseTrain.5, snnEnseTrain.6), "snn")

enseTrain <- rbind(glmet, cqoet, mnnet, snnet)
rm(glmet, cqoet, mnnet, snnet)

enseTrain2 <- merge(enseTrain, enseTrain, by=c("variable", "iteration", "period","test","model"))
enseTrain2 <- enseTrain2[which(enseTrain2$thres.x != enseTrain2$thres.y),]


##########################
# DISTRIBUTION RANGE TEST
glmdr <- meltDistRange(list(glmDistRange.1, glmDistRange.2, glmDistRange.3, glmDistRange.4, glmDistRange.5, glmDistRange.6), "glm")
cqodr <- meltDistRange(list(cqoDistRange.1, cqoDistRange.2, cqoDistRange.3, cqoDistRange.4, cqoDistRange.5, cqoDistRange.6), "cqo")
mnndr <- meltDistRange(list(mnnDistRange.1, mnnDistRange.2, mnnDistRange.3, mnnDistRange.4, mnnDistRange.5, mnnDistRange.6), "mnn")
snndr <- meltDistRange(list(snnDistRange.1, snnDistRange.2, snnDistRange.3, snnDistRange.4, snnDistRange.5, snnDistRange.6), "snn")

distRange <- rbind(glmdr, cqodr, mnndr, snndr)
rm(glmdr, cqodr, mnndr, snndr)

distRange2 <- merge(distRange, distRange, by=c("variable", "species", "period","model"))

  # Calculate RMA coefficients
  pollenThres <- c("NO_THRES","HALF_P","BIOMES","INTERC","VAR_05","VAR_01")
  
  # Loop for to calculate one regression for each threshold/threshold combination
  rmaDR.model <- list()
  thres.x <- thres.y <- model <- NULL   
  n <- 1
  for(k in c("glm","cqo","mnn","snn")){
    for(i in 1:6){
      for(j in 1:6){
        print(k)
        print(i)
        print(j)
        dat <- distRange2[which(distRange2$thres.x == pollenThres[j] & distRange2$thres.y == pollenThres[i] & distRange2$model == k),]
        rmaDR.model[[n]] <- rmaFit(dat, "value.x", "value.y")
        thres.x[n] <- pollenThres[j]
        thres.y[n] <- pollenThres[i]
        model[n] <- k
        n <- n+1
       }
    }
  }
  # Combine the results
  rmaDR.model <- do.call(rbind, rmaDR.model)
  rmaDR.model <- cbind(rmaDR.model, thres.x, thres.y, model)
  
  # Loop for to calculate one regression for each threshold/threshold combination
  rmaDR.period <- list()
  thres.x <- thres.y <- model <- NULL   
  n <- 1
  for(k in 1:5){
    for(i in 1:6){
      for(j in 1:6){
        print(k)
        print(i)
        print(j)
        dat <- distRange2[which(distRange2$thres.x == pollenThres[j] & distRange2$thres.y == pollenThres[i] & distRange2$period == k),]
        rmaDR.period[[n]] <- rmaFit(dat, "value.x", "value.y")
        thres.x[n] <- pollenThres[j]
        thres.y[n] <- pollenThres[i]
        period[n] <- k
        n <- n+1
       }
    }
  }
  # Combine the results
  rmaDR.period <- do.call(rbind, rmaDR.period)
  rmaDR.period <- cbind(rmaDR.period, thres.x, thres.y, period)

distRange2 <- distRange2[which(distRange2$thres.x != distRange2$thres.y),]
rmaDR.model <- rmaDR.model[which(rmaDR.model$thres.x != rmaDR.model$thres.y),]
rmaDR.period <- rmaDR.period[which(rmaDR.period$thres.x != rmaDR.period$thres.y),]


####################
# DISTRIBUTION TEST
glmdt <- meltDistTest(list(glmDistTest.1, glmDistTest.2, glmDistTest.3, glmDistTest.4, glmDistTest.5, glmDistTest.6), "glm")
cqodt <- meltDistTest(list(cqoDistTest.1, cqoDistTest.2, cqoDistTest.3, cqoDistTest.4, cqoDistTest.5, cqoDistTest.6), "cqo")
mnndt <- meltDistTest(list(mnnDistTest.1, mnnDistTest.2, mnnDistTest.3, mnnDistTest.4, mnnDistTest.5, mnnDistTest.6), "mnn")
snndt <- meltDistTest(list(snnDistTest.1, snnDistTest.2, snnDistTest.3, snnDistTest.4, snnDistTest.5, snnDistTest.6), "snn")

distTest <- rbind(glmdt, cqodt, mnndt, snndt)
distTest2 <- merge(distTest, distTest, by=c("taxon","variable","period","iteration","test","model"))
distTest2 <- distTest2[which(distTest2$thres.x != distTest2$thres.y),]
rm(glmdt, cqodt, mnndt, snndt)


##################
# ASSEMBLAGE TEST
glmet <- meltEnseTest(list(glmEnseTest.1, glmEnseTest.2, glmEnseTest.3, glmEnseTest.4, glmEnseTest.5, glmEnseTest.6), "glm")
cqoet <- meltEnseTest(list(cqoEnseTest.1, cqoEnseTest.2, cqoEnseTest.3, cqoEnseTest.4, cqoEnseTest.5, cqoEnseTest.6), "cqo")
mnnet <- meltEnseTest(list(mnnEnseTest.1, mnnEnseTest.2, mnnEnseTest.3, mnnEnseTest.4, mnnEnseTest.5, mnnEnseTest.6), "mnn")
snnet <- meltEnseTest(list(snnEnseTest.1, snnEnseTest.2, snnEnseTest.3, snnEnseTest.4, snnEnseTest.5, snnEnseTest.6), "snn")

enseTest <- rbind(glmet, cqoet, mnnet, snnet)
rm(glmet, cqoet, mnnet, snnet)

enseTest2 <- merge(enseTest, enseTest, by=c("variable","period","iteration","test","model"))


  # Calculate RMA coefficients
  
  # Loop for to calculate one regression for each threshold/threshold combination
  rmaAT <- list()
  thres.x <- thres.y <- NULL   
  n <- 1
  for(i in 1:6){
    for(j in 1:6){
      print(i)
      print(j)
      dat <- enseTest2[which(enseTest2$thres.x == pollenThres[j] & enseTest2$thres.y == pollenThres[i] ),]
      rmaAT[[n]] <- rmaFit(dat, "value.x", "value.y")
      thres.x[n] <- pollenThres[j]
      thres.y[n] <- pollenThres[i]
      n <- n+1
     }
  }
  # Combine the results
  rmaAT <- do.call(rbind, rmaAT)
  rmaAT <- cbind(rmaAT, thres.x, thres.y)

enseTest2 <- enseTest2[which(enseTest2$thres.x != enseTest2$thres.y),]
rmaAT <- rmaAT[which(rmaAT$thres.x != rmaAT$thres.y),]


################################################################################
#
# PLOT RESULTS
#
################################################################################

library(ggplot2)

#####################
# DISTRIBUTION TRAIN
#
# AUC
#png("Results/Graphs/ThresholdSens/AUC-train.png", 800, 720)
pdf("SDMvsCLM/FINAL/Results/Graphs/ThresholdSens/AUC-train.pdf", 10, 9)
ggplot(disTrain2[which(disTrain2$variable == "auc"),], aes(x=value.x, y=value.y, colour=model)) + 
      geom_point(alpha=0.25) + 
      facet_grid(thres.x ~ thres.y) + 
      xlim(0,1) + 
      ylim(0,1) + 
      geom_abline(linetype="dashed", colour="red", size=0.75) + 
      xlab("AUC") + 
      ylab("AUC") + 
      guides(colour=guide_legend(override.aes=list(alpha = 1))) +
      theme_bw()
dev.off()

# TSS
pdf("SDMvsCLM/FINAL/Results/Graphs/ThresholdSens/TSS-train.pdf", 10, 9)
ggplot(disTrain2[which(disTrain2$variable == "tss"),], aes(x=value.x, y=value.y, colour=model)) + 
      geom_point(alpha=0.25) + 
      facet_grid(thres.x ~ thres.y) + 
      expand_limits(x=c(0.5,1), y=c(0.5,1)) + 
      geom_abline(linetype="dashed", colour="red", size=0.75) + 
      xlab("TSS") + 
      ylab("TSS") + 
      guides(colour=guide_legend(override.aes=list(alpha = 1))) +
      theme_bw()
dev.off()


###################
# ASSEMBLAGE TRAIN
# 
# Species Richness
pdf("SDMvsCLM/FINAL/Results/Graphs/ThresholdSens/spRich-train.pdf", 10, 9)
ggplot(enseTrain2[which(enseTrain2$variable == "spRich.COR"),], aes(x=value.x, y=value.y, colour=model)) + 
      geom_point() + 
      facet_grid(thres.y ~ thres.x) + 
      xlim(-1, 1) + 
      ylim(-1,1) + 
      geom_abline(linetype="dashed", colour="grey", size=0.75) + 
      xlab("Species richness correlation\n between observed and predicted communities") + 
      ylab("Species richness correlation\n between observed and predicted communities") + 
      theme_bw()
dev.off()

# Jaccard
pdf("SDMvsCLM/FINAL/Results/Graphs/ThresholdSens/jacc-train.pdf", 10, 9)
ggplot(enseTrain2[which(enseTrain2$variable == "jacc.MEAN"),], aes(x=value.x, y=value.y, colour=model)) + 
      geom_point() + 
      facet_grid(thres.y ~ thres.x) + 
      xlim(0, 1) + 
      ylim(0,1) + 
      geom_abline(linetype="dashed", colour="grey", size=0.75) + 
      xlab("Dissimilarity between predicted and observed communities\n(Jaccard)") + 
      ylab("Dissimilarity between predicted and observed communities\n(Jaccard)") + 
      theme_bw()
dev.off()

# Betadiversity
pdf("SDMvsCLM/FINAL/Results/Graphs/ThresholdSens/beta-train.pdf", 10, 9)
ggplot(enseTrain2[which(enseTrain2$variable == "beta.SOR.DIF"),], aes(x=value.x, y=value.y, colour=model)) + 
      geom_point() + 
      facet_grid(thres.y ~ thres.x) + 
      xlim(-1, 1) + 
      ylim(-1,1) + 
      geom_abline(linetype="dashed", colour="grey", size=0.75) + 
      xlab("Difference between predicted and observed betadiversity\n(Sorensen)") + 
      ylab("Difference between predicted and observed betadiversity\n(Sorensen)") +
      theme_bw()
dev.off()

# Turnover
pdf("SDMvsCLM/FINAL/Results/Graphs/ThresholdSens/turn-train.pdf", 10, 9)
ggplot(enseTrain2[which(enseTrain2$variable == "beta.SIM.DIF"),], aes(x=value.x, y=value.y, colour=model)) + 
      geom_point() + 
      facet_grid(thres.y ~ thres.x) + 
      xlim(-1, 1) + 
      ylim(-1,1) + 
      geom_abline(linetype="dashed", colour="grey", size=0.75) + 
      xlab("Difference between predicted and observed turnover\n(Sorensen)") + 
      ylab("Difference between predicted and observed turnover\n(Sorensen)") + 
      theme_bw()
dev.off()

# Nestedness
pdf("SDMvsCLM/FINAL/Results/Graphs/ThresholdSens/nest-train.pdf", 10, 9)
ggplot(enseTrain2[which(enseTrain2$variable == "beta.SNE.DIF"),], aes(x=value.x, y=value.y, colour=model)) + 
      geom_point() + 
      facet_grid(thres.y ~ thres.x) + 
      xlim(-1, 1) + 
      ylim(-1,1) + 
      geom_abline(linetype="dashed", colour="grey", size=0.75) + 
      xlab("Difference between predicted and observed nestedness\n(Sorensen)") +
      ylab("Difference between predicted and observed nestedness\n(Sorensen)") +
      theme_bw()
dev.off()


#################################
# DISTRIBUTION RANGE ACROSS TIME
#
# Distribution range by periods
pdf("SDMvsCLM/FINAL/Results/Graphs/ThresholdSens/DistRange-period.pdf", 10, 9)
ggplot(distRange2[which(distRange2$variable == "area"),], aes(x=value.x, y=value.y, color=period)) + 
      geom_point(size=1.5, alpha=0.25) + 
      facet_grid(thres.y ~ thres.x) + 
      xlab("Area") + 
      ylab("Area") + 
      geom_abline(aes(intercept=intercept, slope=slope, color=factor(period)), data=rmaDR.period, size=1) +
      geom_abline(linetype="dashed", colour="grey", size=0.75) + 
      guides(colour=guide_legend(override.aes=list(alpha = 1))) + 
      theme_bw()
dev.off()

# Distribution range by models
pdf("SDMvsCLM/FINAL/Results/Graphs/ThresholdSens/DistRange-model.pdf", 10, 9)
ggplot(distRange2[which(distRange2$variable == "area"),], aes(x=value.x, y=value.y, color=model)) +
      geom_point(size=1.5, alpha=0.25) + 
      facet_grid(thres.y ~ thres.x) + 
      xlab("Area") + 
      ylab("Area") + 
      geom_abline(aes(intercept=intercept, slope=slope, color=model), data=rmaDR.model, size=1) +
      geom_abline(linetype="dashed", colour="grey", size=0.75) + 
      guides(colour=guide_legend(override.aes=list(alpha = 1))) + 
      theme_bw()
dev.off()

# Distribution range BOXPLOT
pdf("SDMvsCLM/FINAL/Results/Graphs/ThresholdSens/DistRange-boxplot.pdf", 10, 9)
ggplot(distRange[which(distRange$variable == "area"),], aes(x=factor(species), y=value)) +
      geom_boxplot(outlier.shape=NA) + 
      facet_grid(period ~ .) + 
      geom_jitter(aes(colour=thres), alpha=0.25) + 
      ylab("Area") + 
      xlab("Species") + 
      guides(colour=guide_legend(override.aes=list(alpha = 1))) + 
      theme_bw() 
dev.off()


#####################
# DISTRIBUTION TEST
#
# AUC
pdf("SDMvsCLM/FINAL/Results/Graphs/ThresholdSens/AUC-test.pdf", 10, 9)
ggplot(distTest2[which(distTest2$variable == "auc"),], aes(x=value.x, y=value.y, colour=factor(period))) + 
      geom_point(alpha=0.2) + 
      facet_grid(thres.x ~ thres.y) + 
      geom_abline(linetype="dashed", colour="red", size=0.75) + 
      ylab("AUC") + 
      xlab("AUC") + 
      theme_bw() + 
      guides(colour=guide_legend(override.aes=list(alpha = 1))) 
dev.off()

# TSS
pdf("SDMvsCLM/FINAL/Results/Graphs/ThresholdSens/TSS-test.pdf", 10, 9)
ggplot(distTest2[which(distTest2$variable == "tss"),], aes(x=value.x, y=value.y, colour=factor(period))) + 
      geom_point(alpha=0.2) + 
      facet_grid(thres.x ~ thres.y) + 
      geom_abline(linetype="dashed", colour="red", size=0.75) + 
      ylab("TSS") + 
      xlab("TSS") + 
      guides(colour=guide_legend(override.aes=list(alpha = 1))) +
      theme_bw()
dev.off()


###################
# ASSEMBLAGE TEST
#
# Species richness
pdf("SDMvsCLM/FINAL/Results/Graphs/ThresholdSens/spRich-test.pdf", 10, 9)
ggplot(enseTest2[which(enseTest2$variable == "spRich.COR"),], aes(x=value.x, y=value.y, colour=factor(period))) + 
      geom_point(alpha=0.25) + 
      facet_grid(thres.y ~ thres.x) + 
      xlim(-1,1) + 
      ylim(-1,1) + 
      geom_abline(linetype="dashed", colour="grey", size=0.75) + 
      ylab("Species richness correlation") + 
      xlab("Species richness correlation") + 
      guides(colour=guide_legend(override.aes=list(alpha = 1))) + 
      theme_bw()
dev.off()

# Jaccard
pdf("SDMvsCLM/FINAL/Results/Graphs/ThresholdSens/jacc-test.pdf", 10, 9)
ggplot(enseTest2[which(enseTest2$variable == "jacc.MEAN"),], aes(x=value.x, y=value.y, colour=factor(period))) + 
      geom_point(alpha=0.25) + 
      facet_grid(thres.y ~ thres.x) + 
      xlim(0,1) + 
      ylim(0,1) + 
      geom_abline(linetype="dashed", colour="grey", size=0.75) + 
      ylab("Dissimilarity between predicted and observed communities\n(Jaccard)") +
      xlab("Dissimilarity between predicted and observed communities\n(Jaccard)") +
      guides(colour=guide_legend(override.aes=list(alpha = 1))) + 
      theme_bw()
dev.off()

# Betadiversity
pdf("SDMvsCLM/FINAL/Results/Graphs/ThresholdSens/beta-test.pdf", 10, 9)
ggplot(enseTest2[which(enseTest2$variable == "beta.SOR.DIF"),], aes(x=value.x, y=value.y, colour=factor(period))) + 
      geom_point(alpha=0.25) + 
      facet_grid(thres.y ~ thres.x) + 
      xlim(-1,1) + 
      ylim(-1,1) + 
      geom_abline(linetype="dashed", colour="grey", size=0.75) + 
      ylab("Difference between predicted and observed betadiversity\n(Sorensen)") +
      xlab("Difference between predicted and observed betadiversity\n(Sorensen)") +
      guides(colour=guide_legend(override.aes=list(alpha = 1))) + 
      theme_bw()
dev.off()

# Turnover
pdf("SDMvsCLM/FINAL/Results/Graphs/ThresholdSens/turn-test.pdf", 10, 9)
ggplot(enseTest2[which(enseTest2$variable == "beta.SIM.DIF"),], aes(x=value.x, y=value.y, colour=factor(period))) + 
      geom_point(alpha=0.25) + 
      facet_grid(thres.y ~ thres.x) + 
      xlim(-1,1) + 
      ylim(-1,1) + 
      geom_abline(linetype="dashed", colour="grey", size=0.75) + 
      ylab("Difference between predicted and observed turnover\n(Sorensen)") + 
      xlab("Difference between predicted and observed turnover\n(Sorensen)") +
      guides(colour=guide_legend(override.aes=list(alpha = 1))) + 
      theme_bw()
dev.off()

# Nestedness
pdf("SDMvsCLM/FINAL/Results/Graphs/ThresholdSens/nest-test.pdf", 10, 9)
ggplot(enseTest2[which(enseTest2$variable == "beta.SNE.DIF"),], aes(x=value.x, y=value.y, colour=factor(period))) + 
      geom_point(alpha=0.25) + 
      facet_grid(thres.y ~ thres.x) + 
      xlim(-1,1) + 
      ylim(-1,1) + 
      geom_abline(linetype="dashed", colour="grey", size=0.75) + 
      ylab("Difference between predicted and observed nestedness\n(Sorensen)") +
      xlab("Difference between predicted and observed nestedness\n(Sorensen)") +
      guides(colour=guide_legend(override.aes=list(alpha = 1))) + 
      theme_bw()
dev.off()















































