################################################################################
#
# LOAD DATA, PARAMETERS, AND FUNCTIONS
#
################################################################################

# Load functions for the project. It will be changed by a R-package
library(paleoCLMs) # This file doesn't have to be in the project folder. In that way I only keep one copy of the file and the changes will be automatically propagated to other scripts in the project

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
rasTemp <- raster("SDMvsCLM/PaleoCLM/Data/Climate/CCSM/0BP/gdd5_year_ave.tif")

# Average values for those pollen sites in the same grid cell
pollen <- lapply(pollenHQ, reduceDuplicated, rasTemp, weighted=F)
rm(rasTemp)

# Load climate data
clim <- mapply(loadClim, period=period, pollen=pollen, MoreArgs=list(clim_model="CCSM", indir=indir, vars=vars))
clim <- lapply(clim, as.data.frame) # For some of the models I need climate data as data.frame

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
pollenThres[,4:7] <- pollenThres[,4:7]/100

#pollen1 <- lapply(pollen, applyThreshold, pollenThres[,c(1,4)]) #No threshold
#pollen2 <- lapply(pollen, applyThreshold, pollenThres[,c(1,5)]) #Half percent threshold
#pollen3 <- lapply(pollen, applyThreshold, pollenThres[,c(1,6)]) #Biomes threshold
#pollen4 <- lapply(pollen, applyThreshold, pollenThres[,c(1,7)]) #Intercept threshold

# Calculate variable thresholds
var05 <- lapply(pollen, pollenThreshold, 0.05)
#var01 <- lapply(pollen, pollenThreshold, 0.01)
#pollen5 <- mapply(applyThreshold, pollen, threshold=var01, SIMPLIFY=F)
pollen6 <- mapply(applyThreshold, pollen, threshold=var05, SIMPLIFY=F)


# Get a taxon list for training the models and remove all species with less than 20 presences and/or absences. Then propagate the taxon list to all time periods
#taxList1 <- colnames(removeAbsentSpecies(pollen1[[1]], 19, applyTo="low"))
#taxList2 <- colnames(removeAbsentSpecies(pollen2[[1]], 19, applyTo="low"))
#taxList3 <- colnames(removeAbsentSpecies(pollen3[[1]], 19, applyTo="low"))
#taxList4 <- colnames(removeAbsentSpecies(pollen4[[1]], 19, applyTo="low"))
#taxList5 <- colnames(removeAbsentSpecies(pollen5[[1]], 19, applyTo="low"))
taxList6 <- colnames(removeAbsentSpecies(pollen6[[1]], 19, applyTo="low"))

# Combine the taxon list to get only those pollen in all the lists
#taxonNames <- Reduce(intersect, list(taxList1, taxList2, taxList3, taxList4,  taxList5, taxList6))
taxonNames <- taxList6

# Propagate the taxon list to all pollen datasets
#pollen1 <- lapply(pollen1, FUN=function(x){x[,taxonNames]})
#pollen2 <- lapply(pollen2, FUN=function(x){x[,taxonNames]})
#pollen3 <- lapply(pollen3, FUN=function(x){x[,taxonNames]})
#pollen4 <- lapply(pollen4, FUN=function(x){x[,taxonNames]})
#pollen5 <- lapply(pollen5, FUN=function(x){x[,taxonNames]})
pollen6 <- lapply(pollen6, FUN=function(x){x[,taxonNames]})

# Get present data to fit the models
#p0BP1 <- pollen1[[1]]
#p0BP2 <- pollen2[[1]]
#p0BP3 <- pollen3[[1]]
#p0BP4 <- pollen4[[1]]
#p0BP5 <- pollen5[[1]]
p0BP6 <- pollen6[[1]]

# Get a set of present climate conditions for training the models
c0BP <- clim[[1]]
pollen = pollen6
p0BP = p0BP6

# Set all the paleo data for testing but the current are splited in 70/30 for training/testing
trainSets <- lapply(pollen, trainSplit, numIter, 0)
trainSets[[1]] <- trainSplit(p0BP, numIter, 0.7)

################################################################################
#
# GAM TUNING for degrees of freedom
#
################################################################################
## Tune the GAM model for each species
library(mgcv)

commM=p0BP
xData=c0BP
iTrain <- trainSets[[1]][,1]
sub.commM <- subset(commM, iTrain)
sub.xData <- subset(xData, iTrain)
spNames=taxonNames

# species 1 (Ambrosia.type)
# k=5
form=formula(paste(spNames[[1]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=5) + s(wdei_year_ave, k=5) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp1_1=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp1_1)
rsd_sp1_1=residuals(gamfit_sp1_1)
gam(rsd_sp1_1~s((wdei_year_ave),k=10,bs="cs"),gamma=1.4,data=cbind(sub.commM,sub.xData))
# k=5 and 30
form2=formula(paste(spNames[[1]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=5) + s(wdei_year_ave, k=30) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp1_2=gam(form2,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp1_2)

# species 2 (Artemisia)
# k=5
form=formula(paste(spNames[[2]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=5) + s(wdei_year_ave, k=5) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp2_1=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp2_1)

# species 3 (Alnus)
# k=5
form=formula(paste(spNames[[3]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=5) + s(wdei_year_ave, k=5) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp3_1=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp3_1)

# species 4 (Corylus)
# k=5
form=formula(paste(spNames[[4]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=5) + s(wdei_year_ave, k=5) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp4_1=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp4_1)
rsd_sp4_1=residuals(gamfit_sp4_1)
gam(rsd_sp4_1~s((wdei_year_ave),k=10,bs="cs"),gamma=1.4,data=cbind(sub.commM,sub.xData))
# k=5 and 20
form=formula(paste(spNames[[4]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=5) + s(wdei_year_ave, k=20) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp4_2=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp4_2)

# species 5 (Ostrya)
# k=5
form=formula(paste(spNames[[5]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=5) + s(wdei_year_ave, k=5) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp5_1=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp5_1)

# species 6 (Betula)
# k=5
form=formula(paste(spNames[[6]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=5) + s(wdei_year_ave, k=5) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp6_1=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp6_1)

# species 7 (Castanea)
# k=5
form=formula(paste(spNames[[7]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=5) + s(wdei_year_ave, k=5) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp7_1=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp7_1)

# species 8 (Fagus)
# k=5
form=formula(paste(spNames[[8]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=5) + s(wdei_year_ave, k=5) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp8_1=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp8_1)

# species 9 (Quercus)
# k=5
form=formula(paste(spNames[[9]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=5) + s(wdei_year_ave, k=5) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp9_1=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp9_1)

# species 10 (Carya)
# k=5
form=formula(paste(spNames[[10]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=5) + s(wdei_year_ave, k=5) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp10_1=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp10_1)

# Species 11 (Juglans)
# k=5
form=formula(paste(spNames[[11]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=5) + s(wdei_year_ave, k=5) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp11_1=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp11_1)

# species 12 (Tilia)
# k=5
form=formula(paste(spNames[[12]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=5) + s(wdei_year_ave, k=5) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp12_1=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp12_1)
rsd_sp12_1=residuals(gamfit_sp12_1)
gam(rsd_sp12_1~s((aet_year_ave),k=10,bs="cs"),gamma=1.4,data=cbind(sub.commM,sub.xData))
# k=5 and 20
form=formula(paste(spNames[[12]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=20) + s(wdei_year_ave, k=5) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp12_2=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp12_2)

# species 13 (Fraxinus)
# k=5
form=formula(paste(spNames[[13]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=5) + s(wdei_year_ave, k=5) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp13_1=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp13_1)

# species 14 (Abies)
# k=5
form=formula(paste(spNames[[14]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=5) + s(wdei_year_ave, k=5) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp14_1=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp14_1)

# species 15 (Larix)
# k=5
form=formula(paste(spNames[[15]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=5) + s(wdei_year_ave, k=5) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp15_1=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp15_1)
rsd_sp15_1=residuals(gamfit_sp15_1)
gam(rsd_sp15_1~s((aet_year_ave),k=20,bs="cs"),gamma=1.4,data=cbind(sub.commM,sub.xData))
# k=5 and k=20
form=formula(paste(spNames[[15]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=20) + s(wdei_year_ave, k=5) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp15_2=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp15_2)

# species 16 (Picea)
# k=5
form=formula(paste(spNames[[16]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=5) + s(wdei_year_ave, k=5) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp16_1=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp16_1)

# species 17 (Tsuga)
# k=5
form=formula(paste(spNames[[17]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=5) + s(wdei_year_ave, k=5) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp17_1=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp17_1)

# species 18 (Platanus)
# k=5
form=formula(paste(spNames[[18]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=5) + s(wdei_year_ave, k=5) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp18_1=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp18_1)

# species 19 (Rumex)
# k=5
form=formula(paste(spNames[[19]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=5) + s(wdei_year_ave, k=5) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp19_1=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp19_1)
rsd_sp19_1=residuals(gamfit_sp19_1)
gam(rsd_sp19_1~s((aet_year_ave),k=10,bs="cs"),gamma=1.4,data=cbind(sub.commM,sub.xData))
gam(rsd_sp19_1~s((etr_year_ave),k=20,bs="cs"),gamma=1.4,data=cbind(sub.commM,sub.xData))
# k=5,10 and 20
form=formula(paste(spNames[[19]], "~ s(etr_year_ave, k=10) + s(aet_year_ave, k=30) + s(wdei_year_ave, k=5) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp19_2=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp19_2)

# species 20 (Populus)
# k=5
form=formula(paste(spNames[[20]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=5) + s(wdei_year_ave, k=5) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp20_1=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp20_1)

# species 21 (Salix)
# k=5
form=formula(paste(spNames[[21]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=5) + s(wdei_year_ave, k=5) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp21_1=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp21_1)

# species 22 (Acer)
# k=5
for2=formula(paste(spNames[[22]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=5) + s(wdei_year_ave, k=5) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp22_1=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp22_1)
rsd_sp22_1=residuals(gamfit_sp22_1)
gam(rsd_sp22_1~s((tmax_high_quart),k=10,bs="cs"),gamma=1.4,data=cbind(sub.commM,sub.xData))
# k=5 and 10
form=formula(paste(spNames[[22]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=5) + s(wdei_year_ave, k=5) + s(tmax_high_quart, k=10) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp22_2=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp22_2)

# species 23 (Ulmus)
# k=5
form=formula(paste(spNames[[23]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=5) + s(wdei_year_ave, k=5) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp23_1=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp23_1)

# species 24 (Pinus)
# k=5
form=formula(paste(spNames[[24]], "~ s(etr_year_ave, k=5) + s(aet_year_ave, k=5) + s(wdei_year_ave, k=5) + s(tmax_high_quart, k=5) + s(prcp_low_quart, k=5) + s(prcp_high_quart, k=5)"))
gamfit_sp24_1=gam(form,family=binomial,data=cbind(sub.commM,sub.xData))
gam.check(gamfit_sp24_1)

# RESULST FOR INDIVIDUAL SPECIES
# all species use k=5 except the following: 
# species 1: wdei_year)ave, k=30
# species 3: wdei_year_ave, k=20
# species 12: aet_year_ave k=20
# species 15: aet_year_ave k=20
# species 19: aet_year_ave k=30, etr_year_ave k=10
# species 22: tmax_high_quart k=10


################################################################################
#
# CAO TUNING for degrees of freedom
#
################################################################################
require(foreach)
require(doParallel)
registerDoParallel(cores=4) # Check how many cores have your computer and set up a lower number here

caoFit <- function(commM, xData, train, iterNum, dfree){
  require(foreach)
  require(doParallel)
  registerDoParallel(cores=4)
  require(VGAM)
  caoList <- foreach(i=1:iterNum, .packages=c("VGAM"), .verbose=T, .errorhandling='remove') %dopar% {
    iTrain <- train[,i]
    sub.commM <- subset(commM, iTrain)
    sub.xData <- subset(xData, iTrain)
    form <- formula(sub.commM ~ etr_year_ave + aet_year_ave + wdei_year_ave + tmax_high_quart + prcp_low_quart + prcp_high_quart)
    caofit <- cao(form, family=binomialff(mv=TRUE), data=sub.xData, trace=TRUE, df1.nl=dfree)
    caofit
  }
  return(caoList)
}

source("SDMvsCLM/FINAL/Final_Code/Preliminary_Analyses/Pollen_Thresholds/00-ThreshWrappingFunctions.R")
projIndex <-which(period %in% period)

cao1List <- caoFit(p0BP, c0BP, trainSets[[1]], numIter, 1)
cao1DistTrain<- evalDistTrain(cao1List, pollen[[1]], clim[[1]], trainSets[[1]], respType="response")
cao1DistTest <- evalDistTest(cao1List, pollen[projIndex], clim[projIndex], trainSets[projIndex], respType="response")
cao1EnseTrain <- evalEnseTrain(cao1List, pollen[[1]], clim[[1]], trainSets[[1]], respType="response", cao1DistTrain)
cao1EnseTest <- evalEnseTest(cao1List, pollen[projIndex], clim[projIndex], trainSets[projIndex], respType="response", cao1DistTrain)

cao2List <- caoFit(p0BP, c0BP, trainSets[[1]], numIter, 2)
cao2DistTrain<- evalDistTrain(cao2List, pollen[[1]], clim[[1]], trainSets[[1]], respType="response")
cao2DistTest <- evalDistTest(cao2List, pollen[projIndex], clim[projIndex], trainSets[projIndex], respType="response")
cao2EnseTrain <- evalEnseTrain(cao2List, pollen[[1]], clim[[1]], trainSets[[1]], respType="response", cao2DistTrain)
cao2EnseTest <- evalEnseTest(cao2List, pollen[projIndex], clim[projIndex], trainSets[projIndex], respType="response", cao2DistTrain)

cao3List <- caoFit(p0BP, c0BP, trainSets[[1]], numIter, 3)
cao3DistTrain<- evalDistTrain(cao3List, pollen[[1]], clim[[1]], trainSets[[1]], respType="response")
cao3DistTest <- evalDistTest(cao3List, pollen[projIndex], clim[projIndex], trainSets[projIndex], respType="response")
cao3EnseTrain <- evalEnseTrain(cao3List, pollen[[1]], clim[[1]], trainSets[[1]], respType="response", cao3DistTrain)
cao3EnseTest <- evalEnseTest(cao3List, pollen[projIndex], clim[projIndex], trainSets[projIndex], respType="response", cao3DistTrain)

cao4List <- caoFit(p0BP, c0BP, trainSets[[1]], numIter, 4)
cao4DistTrain<- evalDistTrain(cao4List, pollen[[1]], clim[[1]], trainSets[[1]], respType="response")
cao4DistTest <- evalDistTest(cao4List, pollen[projIndex], clim[projIndex], trainSets[projIndex], respType="response")
cao4EnseTrain <- evalEnseTrain(cao4List, pollen[[1]], clim[[1]], trainSets[[1]], respType="response", cao4DistTrain)
cao4EnseTest <- evalEnseTest(cao4List, pollen[projIndex], clim[projIndex], trainSets[projIndex], respType="response", cao4DistTrain)

cao5List <- caoFit(p0BP, c0BP, trainSets[[1]], numIter, 5)
cao5DistTrain<- evalDistTrain(cao5List, pollen[[1]], clim[[1]], trainSets[[1]], respType="response")
cao5DistTest <- evalDistTest(cao5List, pollen[projIndex], clim[projIndex], trainSets[projIndex], respType="response")
cao5EnseTrain <- evalEnseTrain(cao5List, pollen[[1]], clim[[1]], trainSets[[1]], respType="response", cao5DistTrain)
cao5EnseTest <- evalEnseTest(cao5List, pollen[projIndex], clim[projIndex], trainSets[projIndex], respType="response", cao5DistTrain)


################################################################################
#
# PLOTING CAO TUNING RESULTS
#
################################################################################

library(reshape2)

cao1d <- melt(cao1DistTest)
cao2d <- melt(cao2DistTest)
cao3d <- melt(cao3DistTest)
cao4d <- melt(cao4DistTest)
cao5d <- melt(cao5DistTest)
colnames(cao1d) <- colnames(cao2d) <- colnames(cao3d) <- colnames(cao4d) <- colnames(cao5d) <- c("taxon","variable","value","period","iteration")

resCAO <- rbind(cao1d, cao2d, cao3d, cao4d, cao5d)
colnames(resCAO) <- c("taxon","variable","value","period","iteration")
resCAO <- cbind(resCAO, model=c(rep("cao1", nrow(cao1d)), rep("cao2", nrow(cao2d)), rep("cao3", nrow(cao3d)), rep("cao4", nrow(cao4d)), rep("cao5", nrow(cao5d))))

resCAO <- resCAO[which(resCAO$period == 1),]

aucCAO <- resCAO[which(resCAO$variable == "auc"),]
tssCAO <- resCAO[which(resCAO$variable == "tss"),]
specCAO <- resCAO[which(resCAO$variable == "spec"),]
sensCAO <- resCAO[which(resCAO$variable == "sens"),]
prevCAO <- resCAO[which(resCAO$variable == "prev"),]
trminCAO <- resCAO[which(resCAO$variable == "trmin"),]


cao1e <- melt(cao1EnseTest)
cao2e <- melt(cao2EnseTest)
cao3e <- melt(cao3EnseTest)
cao4e <- melt(cao4EnseTest)
cao5e <- melt(cao5EnseTest)
colnames(cao1e) <- colnames(cao2e) <- colnames(cao3e) <- colnames(cao4e) <- colnames(cao5e) <- c("variable","period","value","iteration")

reseCAO <- rbind(cao1e, cao2e, cao3e, cao4e, cao5e)
reseCAO <- cbind(reseCAO, model=c(rep("cao1", nrow(cao1e)), rep("cao2", nrow(cao2e)), rep("cao3", nrow(cao3e)), rep("cao4", nrow(cao4e)), rep("cao5", nrow(cao5e))))

reseCAO <- reseCAO[which(reseCAO$period == 1),]

beta.sim.dif <- reseCAO[which(reseCAO$variable == "beta.SIM.DIF"),]
beta.sne.dif <- reseCAO[which(reseCAO$variable == "beta.SNE.DIF"),]
beta.sor.dif <- reseCAO[which(reseCAO$variable == "beta.SOR.DIF"),]
jacc.mean <- reseCAO[which(reseCAO$variable == "jacc.MEAN"),]
jacc.sd <- reseCAO[which(reseCAO$variable == "jacc.SD"),]
sp.rich.cor <- reseCAO[which(reseCAO$variable == "spRich.COR"),]

library(ggplot2)
p1 <- ggplot(data=aucCAO, aes(model, value)) + geom_boxplot(notch=T, outlier.shape=NA) + ylim(0, 1) + labs(title="AUC", x="", y="") + geom_abline(aes(intercept=0.5, slope=0, colour="red"))
p2 <- ggplot(data=tssCAO, aes(model, value)) + geom_boxplot(notch=T, outlier.shape=NA) + ylim(-0.5, 1) + labs(title="TSS", x="", y="") + geom_abline(aes(intercept=0, slope=0, colour="red"))
p3 <- ggplot(data=specCAO, aes(model, value)) + geom_boxplot(notch=T, outlier.shape=NA) + ylim(0, 1) + labs(title="Specificity", x="", y="")
p4 <- ggplot(data=sensCAO, aes(model, value)) + geom_boxplot(notch=T, outlier.shape=NA) + ylim(0, 1) + labs(title="Sensitivity", x="", y="")

pdf("SDMvsCLM/FINAL/Results/Graphs/Tuning/CAO-Test-Dist-Sensitivity_0BP.pdf", width=20/2.54, height=20/2.54)
multiPlot(p1, p4, p2, p3, cols=2)
dev.off()


p5 <- ggplot(data=jacc.mean, aes(model, value)) + geom_boxplot(outlier.shape=NA) + ylim(0, 1) + labs(title="Obs. vs Pred. communities", x="", y="Jaccard index")

p6 <- ggplot(data=sp.rich.cor, aes(model, value)) + geom_boxplot(outlier.shape=NA) + ylim(0, 1) + labs(title="Species richness", x="", y="Cor. Pearson")

p7 <- ggplot(data=beta.sor.dif, aes(model, value)) + geom_boxplot(outlier.shape=NA) + ylim(-0.025, 0.075) + labs(title="Reg. betadiversity", x="", y="Increase (Sorensen)") + geom_abline(aes(intercept=0, slope=0, colour="red"))
p8 <- ggplot(data=beta.sim.dif, aes(model, value)) + geom_boxplot(outlier.shape=NA) + ylim(-0.1, 0.5) + labs(title="Reg. turnover", x="", y="Increase (Sorensen)") + geom_abline(aes(intercept=0, slope=0, colour="red"))
p9 <- ggplot(data=beta.sne.dif, aes(model, value)) + geom_boxplot(outlier.shape=NA) + ylim(-0.025, 0.025) + labs(title="Reg. nestedness", x="", y="Increase (Sorensen)") + geom_abline(aes(intercept=0, slope=0, colour="red"))


pdf("SDMvsCLM/FINAL/Results/Graphs/Tuning/CAO-Test-Ensem-Sensitivity_0BP.pdf", width=20/2.54, height=30/2.54)
lay <- matrix(c(1,2,3,1,4,5), nrow=3, ncol=2)
multiPlot(p5, p6, p8, p7, p9, layout=lay)
dev.off()

################################################################################
#
# MARS TUNING FOR NPRUNE AND DEGREE - SDM and CLM
#
# nprune: number of terms in the model after pruning (including intercept)
# degree: degree of interactions, default is 1 (no interaction between predictors)
################################################################################
library(earth)
library(foreach)
library(caret)

## Tuning of SDM version of MARS
fitControl <- trainControl(method = "repeatedcv", number=10, repeats = 50)
marsGrid <- expand.grid(nprune=c(1, 2, 5, 7, 9, 10),degree=c(1, 2,3,4,5))
marsfits <- lapply(taxonNames,FUN=function(t){train(x=c0BP, y=p0BP[,t], method = "earth", metric="Rsquared", trControl=fitControl, tuneGrid=marsGrid,trace=T,glm=list(family=binomial))})

# plot the results
lapply(marsfits,FUN=function(x){plot(x)})
lapply(marsfits,FUN=function(x){plot(x, metric="Rsquared")}) 

# store the nprune and degree for the best model for each species
smarsSensMatrix1=matrix(NA, 24, 2)
for (i in 1:length(marsfits)){
  smarsSensMatrix1[i,1]=marsfits[[i]]$bestTune$nprune
  smarsSensMatrix1[i,2]=marsfits[[i]]$bestTune$degree
}
colnames(smarsSensMatrix1)=c("nprune","degree")
rownames(smarsSensMatrix1)=taxonNames

smarsSensMatrix1

# Because the train function produced unstable results, I ran it multiple times. I then used results from one of the runs in which repeats=50 and in cases where there is variance in results, chose the nprune/degree combo with the least number of parameters.

# Tuning of the CLM version
# nprune 5 - degree 1
library(foreach)
library(doParallel)
registerDoParallel(cores=4)

mars51List <- foreach(i=1:numIter, .packages=c("earth"), .verbose=T) %dopar% {
  iTrainSets <- lapply(trainSets, function(x){x[,i]})
  sub.p0BP <- subset(p0BP, iTrainSets[[1]])
  sub.c0BP <- subset(c0BP, iTrainSets[[1]])
  marsfit <- earth(x=sub.c0BP, y=sub.p0BP, nprune=5, degree=1, trace=T, glm=list(family=binomial))
  marsfit
}

mars51DistTest <- evalDistTest(mars51List, pollen[projIndex], clim[projIndex], trainSets[projIndex], respType="response")
mars51EnseTest <- evalEnseTest(mars51List, pollen[projIndex], clim[projIndex], trainSets[projIndex], respType="response", cao1DistTrain)


# nprune 5 - degree 2
mars52List <- foreach(i=1:numIter, .packages=c("earth"), .verbose=T) %dopar% {
  iTrainSets <- lapply(trainSets, function(x){x[,i]})
  sub.p0BP <- subset(p0BP, iTrainSets[[1]])
  sub.c0BP <- subset(c0BP, iTrainSets[[1]])
  marsfit <- earth(x=sub.c0BP, y=sub.p0BP, nprune=5, degree=2, trace=T, glm=list(family=binomial))
  marsfit
}

mars52DistTest <- evalDistTest(mars52List, pollen[projIndex], clim[projIndex], trainSets[projIndex], respType="response")
mars52EnseTest <- evalEnseTest(mars52List, pollen[projIndex], clim[projIndex], trainSets[projIndex], respType="response", cao1DistTrain)


# nprune 5 - degree 3
mars53List <- foreach(i=1:numIter, .packages=c("earth"), .verbose=T) %dopar% {
  iTrainSets <- lapply(trainSets, function(x){x[,i]})
  sub.p0BP <- subset(p0BP, iTrainSets[[1]])
  sub.c0BP <- subset(c0BP, iTrainSets[[1]])
  marsfit <- earth(x=sub.c0BP, y=sub.p0BP, nprune=5, degree=3, trace=T, glm=list(family=binomial))
  marsfit
}

mars53DistTest <- evalDistTest(mars53List, pollen[projIndex], clim[projIndex], trainSets[projIndex], respType="response")
mars53EnseTest <- evalEnseTest(mars53List, pollen[projIndex], clim[projIndex], trainSets[projIndex], respType="response", cao1DistTrain)


# nprune 9 - degree 1
mars91List <- foreach(i=1:numIter, .packages=c("earth"), .verbose=T) %dopar% {
  iTrainSets <- lapply(trainSets, function(x){x[,i]})
  sub.p0BP <- subset(p0BP, iTrainSets[[1]])
  sub.c0BP <- subset(c0BP, iTrainSets[[1]])
  marsfit <- earth(x=sub.c0BP, y=sub.p0BP, nprune=9, degree=1, trace=T, glm=list(family=binomial))
  marsfit
}

mars91DistTest <- evalDistTest(mars91List, pollen[projIndex], clim[projIndex], trainSets[projIndex], respType="response")
mars91EnseTest <- evalEnseTest(mars91List, pollen[projIndex], clim[projIndex], trainSets[projIndex], respType="response", cao1DistTrain)


# nprune 9 - degree 2
mars92List <- foreach(i=1:numIter, .packages=c("earth"), .verbose=T) %dopar% {
  iTrainSets <- lapply(trainSets, function(x){x[,i]})
  sub.p0BP <- subset(p0BP, iTrainSets[[1]])
  sub.c0BP <- subset(c0BP, iTrainSets[[1]])
  marsfit <- earth(x=sub.c0BP, y=sub.p0BP, nprune=9, degree=2, trace=T, glm=list(family=binomial))
  marsfit
}

mars92DistTest <- evalDistTest(mars92List, pollen[projIndex], clim[projIndex], trainSets[projIndex], respType="response")
mars92EnseTest <- evalEnseTest(mars92List, pollen[projIndex], clim[projIndex], trainSets[projIndex], respType="response", cao1DistTrain)


# nprune 9 - degree 3
mars93List <- foreach(i=1:numIter, .packages=c("earth"), .verbose=T) %dopar% {
  iTrainSets <- lapply(trainSets, function(x){x[,i]})
  sub.p0BP <- subset(p0BP, iTrainSets[[1]])
  sub.c0BP <- subset(c0BP, iTrainSets[[1]])
  marsfit <- earth(x=sub.c0BP, y=sub.p0BP, nprune=9, degree=3, trace=T, glm=list(family=binomial))
  marsfit
}

mars93DistTest <- evalDistTest(mars93List, pollen[projIndex], clim[projIndex], trainSets[projIndex], respType="response")
mars93EnseTest <- evalEnseTest(mars93List, pollen[projIndex], clim[projIndex], trainSets[projIndex], respType="response", cao1DistTrain)


# nprune 7 - degree 2
mars72List <- foreach(i=1:numIter, .packages=c("earth"), .verbose=T) %dopar% {
  iTrainSets <- lapply(trainSets, function(x){x[,i]})
  sub.p0BP <- subset(p0BP, iTrainSets[[1]])
  sub.c0BP <- subset(c0BP, iTrainSets[[1]])
  marsfit <- earth(x=sub.c0BP, y=sub.p0BP, nprune=7, degree=2, trace=T, glm=list(family=binomial))
  marsfit
}

mars72DistTest <- evalDistTest(mars72List, pollen[projIndex], clim[projIndex], trainSets[projIndex], respType="response")
mars72EnseTest <- evalEnseTest(mars72List, pollen[projIndex], clim[projIndex], trainSets[projIndex], respType="response", cao1DistTrain)


# PLOTING - CLM MARS


library(reshape2)

d51 <- melt(mars51DistTest)
d52 <- melt(mars52DistTest)
d53 <- melt(mars53DistTest)
d91 <- melt(mars91DistTest)
d92 <- melt(mars92DistTest)
d93 <- melt(mars93DistTest)
d72 <- melt(mars72DistTest)
colnames(d51) <- colnames(d52) <- colnames(d53) <- colnames(d91) <- colnames(d92) <- colnames(d93) <- colnames(d72) <- c("taxon","variable","value","iteration")

resMARS <- rbind(d51, d52, d53, d91, d92, d93, d72)
colnames(resMARS) <- c("taxon","variable","value","iteration")
resMARS <- cbind(resMARS, model=c(rep("mars51", nrow(d51)), rep("mars52", nrow(d52)), rep("mars53", nrow(d53)), rep("mars91", nrow(d91)), rep("mars92", nrow(d92)), rep("mars93", nrow(d93)), rep("mars72", nrow(d72))))

aucMARS <- resMARS[which(resMARS$variable == "auc"),]
tssMARS <- resMARS[which(resMARS$variable == "tss"),]
specMARS <- resMARS[which(resMARS$variable == "spec"),]
sensMARS <- resMARS[which(resMARS$variable == "sens"),]
prevMARS <- resMARS[which(resMARS$variable == "prev"),]
trminMARS <- resMARS[which(resMARS$variable == "trmin"),]


e51 <- melt(mars51EnseTest)
e52 <- melt(mars52EnseTest)
e53 <- melt(mars53EnseTest)
e91 <- melt(mars91EnseTest)
e92 <- melt(mars92EnseTest)
e93 <- melt(mars93EnseTest)
e72 <- melt(mars72EnseTest)
colnames(e51) <- colnames(e52) <- colnames(e53) <- colnames(e91) <- colnames(e92) <- colnames(e93) <- colnames(e72) <- c("variable","period","value","iteration")

reseMARS <- rbind(e51, e52, e53, e91, e92, e93, e72)
reseMARS <- cbind(reseMARS, model=c(rep("mars51", nrow(e51)), rep("mars52", nrow(e52)), rep("mars53", nrow(e53)), rep("mars91", nrow(e91)), rep("mars92", nrow(e92)), rep("mars93", nrow(e93)), rep("mars72", nrow(e72))))

beta.sim.dif.MARS <- reseMARS[which(reseMARS$variable == "beta.SIM.DIF"),]
beta.sne.dif.MARS <- reseMARS[which(reseMARS$variable == "beta.SNE.DIF"),]
beta.sor.dif.MARS <- reseMARS[which(reseMARS$variable == "beta.SOR.DIF"),]
jacc.mean.MARS <- reseMARS[which(reseMARS$variable == "jacc.MEAN"),]
jacc.sd.MARS <- reseMARS[which(reseMARS$variable == "jacc.SD"),]
sp.rich.cor.MARS <- reseMARS[which(reseMARS$variable == "spRich.COR"),]

library(ggplot2)
p1 <- ggplot(data=aucMARS, aes(factor(model), value)) + geom_boxplot(notch=T, outlier.shape=NA) + ylim(0, 1) + labs(title="AUC", x="", y="") + geom_abline(aes(intercept=0.5, slope=0, colour="red"))
p2 <- ggplot(data=tssMARS, aes(factor(model), value)) + geom_boxplot(notch=T, outlier.shape=NA) + ylim(-0.5, 1) + labs(title="TSS", x="", y="") + geom_abline(aes(intercept=0, slope=0, colour="red"))
p3 <- ggplot(data=specMARS, aes(factor(model), value)) + geom_boxplot(notch=T, outlier.shape=NA) + ylim(0, 1) + labs(title="Specificity", x="", y="")
p4 <- ggplot(data=sensMARS, aes(factor(model), value)) + geom_boxplot(notch=T, outlier.shape=NA) + ylim(0, 1) + labs(title="Sensitivity", x="", y="")

pdf("SDMvsCLM/FINAL/Results/Graphs/Tuning/MARS-Test-Dist-Sensitivity.pdf", width=20/2.54, height=20/2.54)
multiPlot(p1, p4, p2, p3, cols=2)
dev.off()


p5 <- ggplot(data=jacc.mean.MARS, aes(factor(model), value)) + geom_boxplot(outlier.shape=NA) + ylim(0, 1) + labs(title="Obs. vs Pred. communities", x="", y="Jaccard index")

p6 <- ggplot(data=sp.rich.cor.MARS, aes(factor(model), value)) + geom_boxplot(outlier.shape=NA) + ylim(0, 1) + labs(title="Species richness", x="", y="Cor. Pearson")

p7 <- ggplot(data=beta.sor.dif.MARS, aes(factor(model), value)) + geom_boxplot(outlier.shape=NA) + ylim(-0.025, 0.075) + labs(title="Reg. betadiversity", x="", y="Increase (Sorensen)") + geom_abline(aes(intercept=0, slope=0, colour="red"))
p8 <- ggplot(data=beta.sim.dif.MARS, aes(factor(model), value)) + geom_boxplot(outlier.shape=NA) + ylim(-0.1, 0.5) + labs(title="Reg. turnover", x="", y="Increase (Sorensen)") + geom_abline(aes(intercept=0, slope=0, colour="red"))
p9 <- ggplot(data=beta.sne.dif.MARS, aes(factor(model), value)) + geom_boxplot(outlier.shape=NA) + ylim(-0.025, 0.025) + labs(title="Reg. nestedness", x="", y="Increase (Sorensen)") + geom_abline(aes(intercept=0, slope=0, colour="red"))


pdf("SDMvsCLM/FINAL/Results/Graphs/Tuning/MARS-Test-Ensem-Sensitivity.pdf", width=20/2.54, height=30/2.54)
lay <- matrix(c(1,2,3,1,4,5), nrow=3, ncol=2)
multiPlot(p5, p6, p8, p7, p9, layout=lay)
dev.off()



################################################################################
#
# SDM CARTS TUNING FOR Complexity Parameter (CP)
#
################################################################################
library(caret)
library(rpart)
require(foreach)

## Tuning of SDM version of CARTS
fitControl <- trainControl(method = "repeatedcv", number=10, repeats = 50)
cartsGrid <- expand.grid(cp=c(0.01, 0.02, 0.05, 0.07, 0.09, .1,.5,.9))
cartsfits <- lapply(taxonNames,FUN=function(t){train(x=c0BP, y=p0BP[,t], method = "rpart", metric="RMSE", trControl=fitControl, tuneGrid=cartsGrid)})

# plot the results
lapply(marsfits,FUN=function(x){plot(x)})
lapply(marsfits,FUN=function(x){plot(x, metric="Rsquared")})

# store the results
cartsSensMatrix=matrix(NA, 24, 1)
for (i in 1:length(cartsfits)){
  cartsSensMatrix[i,1]=cartsfits[[i]]$bestTune$cp
}
colnames(cartsSensMatrix)=c("cp")
rownames(cartsSensMatrix)=taxonNames

cartsSensMatrix

################################################################################
#
# TUNING NEURAL NETWORKS FOR # OF HIDDEN LAYERS AND DECAY RATE 
#
################################################################################
library(nnet)
library(caret)
require(foreach)


## Tuning of SDM version of CARTS
fitControl <- trainControl(method = "repeatedcv", number=10, repeats = 50)
nnetGrid <- expand.grid(size=c(1, 5, 10, 15, 20, 25),decay=c(0, 1, 5, 10))
snnetfits <- lapply(taxonNames, FUN=function(t){train(x=c0BP, y=p0BP[,t], method = "nnet", metric="Rsquared", maxit=2000, skip=F, trControl=fitControl, tuneGrid=nnetGrid, verbose=FALSE)})

# Plot results
lapply(marsfits,FUN=function(x){plot(x)})
lapply(marsfits,FUN=function(x){plot(x, metric="RMSE")})

# Store results
snnetSensMatrix1=matrix(NA, 24, 2)
for (i in 1:length(snnetfits)){
  snnetSensMatrix1[i,1]=snnetfits[[i]]$bestTune$size
  snnetSensMatrix1[i,2]=snnetfits[[i]]$bestTune$decay
}
colnames(snnetSensMatrix1)=c("nprune","degree")
rownames(snnetSensMatrix1)=taxonNames

snnetSensMatrix1

# I ran the train function multiple times. Results: most models peak at size=5 and decay=1, there are only very suttle (<0.01) differences above size=5.
# All species should be set to decay=1.

## Tuning of CLM version of CARTS
# FINAL results: decay=1, size=20. According to Rsquared and RMSE.

fitControl <- trainControl(method = "repeatedcv", number=10, repeats = 50)
nnetGrid <- expand.grid(size=c(1, 5, 10, 15, 20, 25),decay=c(0, 1, 5, 10))
form=formula(Ambrosia.type+Artemisia+Alnus+Corylus+Ostrya.Carpinus+Betula+Castanea+Fagus+Quercus+Carya+Juglans+Tilia+Fraxinus+Abies+Larix+Picea+Tsuga+Platanus+Rumex.Oxyria+Populus+Salix+Acer+Ulmus+Pinus~etr_year_ave+aet_year_ave+wdei_year_ave+tmax_high_quart+prcp_low_quart+prcp_high_quart)
nnetfits <- train(form, data=cbind(p0BP,c0BP), method = "nnet", maxit=2000, skip=F, metric="Rsquared", trControl=fitControl, tuneGrid=nnetGrid, verbose=FALSE)

# Plot results
plot(nnetfits)
plot(nnetfits, metric="RMSE")



