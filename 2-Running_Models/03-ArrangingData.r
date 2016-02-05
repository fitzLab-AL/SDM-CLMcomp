##Load plotting libraries
library(reshape2)
library(ggplot2)

################################################################################
##
## ARRANGING DATA
##
################################################################################

##Load data
load("Results/RObjects/DisResults.RData")
load("Results/RObjects/AssResults.RData")

resd.all <- resd.All
rese.all <- rese.All

##Name the data
names(resd.all) <- names(rese.all) <- periods_fit

resd.all <- resd.all
rese.all <- rese.all

##Remove the missing data and melt the data in a single dataframe
resd.all <- melt(resd.all, id.vars=c("taxon","variable","iteration","model","periodFit","periodPrj","test"))
rese.all <- melt(rese.all, id.vars=c("variable","iteration","model","periodFit","periodPrj","test"))

resd.all$value[which(resd.all$value == -9999)] <- NA
rese.all$value[which(rese.all$value == -9999)] <- NA

##Change column names in the new dataframe
resd.all <- resd.all[,-which(colnames(resd.all) == "L1")]
rese.all <- rese.all[,-which(colnames(rese.all) == "L1")]

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
                #results <- t.test(sdm$value, clm$value, paired=T)
                results <- wilcox.test(sdm$value, clm$value, paired=T)
                stat.gam[j,l] <- results$p.value
            }
            
            if(length(na.omit(unique(dat2$value[which(dat2$model == "glm")]))) < 3 | length(na.omit(unique(dat2$value[which(dat2$model == "cqo")]))) < 3){
                stat.gam[j,l] <- NA
            }else{
                sdm <- dat2[which(dat2$model == "glm"),]
                clm <- dat2[which(dat2$model == "cqo"),]
                #results <- t.test(sdm$value, clm$value, paired=T)
                results <- wilcox.test(sdm$value, clm$value, paired=T)
                stat.glm[j,l] <- results$p.value
            }
            
            sdm <- dat2[which(dat2$model == "scarts"),]
            clm <- dat2[which(dat2$model == "carts"),]
            #results <- t.test(sdm$value, clm$value, paired=T)
            results <- wilcox.test(sdm$value, clm$value, paired=T)
            stat.scarts[j,l] <- results$p.value
            
            sdm <- dat2[which(dat2$model == "smars"),]
            clm <- dat2[which(dat2$model == "mars"),]
            #results <- t.test(sdm$value, clm$value, paired=T)
            results <- wilcox.test(sdm$value, clm$value, paired=T)
            stat.smars[j,l] <- results$p.value
            
            sdm <- dat2[which(dat2$model == "snn"),]
            clm <- dat2[which(dat2$model == "mnn"),]
            #results <- t.test(sdm$value, clm$value, paired=T)
            results <- wilcox.test(sdm$value, clm$value, paired=T)
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
# auc.ttest <- sdmStats(resd.test, periods_prj, periods_fit, "auc")
# save(auc.ttest, file="Results/RObjects/Ttest/auc.ttest")
# 
# brier.ttest <- sdmStats(resd.test, periods_prj, periods_fit, "brier")
# save(brier.ttest, file="Results/RObjects/Ttest/brier.ttest")
# 
# tss.ttest <- sdmStats(resd.test, periods_prj, periods_fit, "tss")
# save(tss.ttest, file="Results/RObjects/Ttest/tss.ttest")
# 
# sens.ttest <- sdmStats(resd.test, periods_prj, periods_fit, "sens")
# save(sens.ttest, file="Results/RObjects/Ttest/sens.ttest")
# 
# spec.ttest <- sdmStats(resd.test, periods_prj, periods_fit, "spec")
# save(spec.ttest, file="Results/RObjects/Ttest/spec.ttest")
# 
# jacc.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "jacc.MEAN.prob")
# save(jacc.ttest, file="Results/RObjects/Ttest/jacc.ttest")
# 
# jacDis.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "jacc.DIS.prob")
# save(jacDis.ttest, file="Results/RObjects/Ttest/jacDis.ttest")
# 
# jacCal.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "jacc.CAL.prob")
# save(jacCal.ttest, file="Results/RObjects/Ttest/jacCal.ttest")
# 
# spRichCor.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "spRich.COR.prob")
# save(spRichCor.ttest, file="Results/RObjects/Ttest/spRichCor.ttest")
# 
# spRichRsq.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "spRich.RSQ.prob")
# save(spRichRsq.ttest, file="Results/RObjects/Ttest/spRichRsq.ttest")
# 
# jacct.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "jacc.MEAN.thre")
# save(jacct.ttest, file="Results/RObjects/Ttest/jacct.ttest")
# 
# jacDist.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "jacc.DIS.thre")
# save(jacDist.ttest, file="Results/RObjects/Ttest/jacDist.ttest")
# 
# jacCalt.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "jacc.CAL.thre")
# save(jacCalt.ttest, file="Results/RObjects/Ttest/jacCalt.ttest")
# 
# spRichCort.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "spRich.COR.thre")
# save(spRichCort.ttest, file="Results/RObjects/Ttest/spRichCort.ttest")
# 
# spRichRsqt.ttest <- sdmStats(rese.test, periods_prj, periods_fit, "spRich.RSQ.thre")
# save(spRichRsqt.ttest, file="Results/RObjects/Ttest/spRichRsqt.ttest")

auc.wtest <- sdmStats(resd.test, periods_prj, periods_fit, "auc")
save(auc.wtest, file="Results/RObjects/wtest/auc.wtest")

brier.wtest <- sdmStats(resd.test, periods_prj, periods_fit, "brier")
save(brier.wtest, file="Results/RObjects/wtest/brier.wtest")

tss.wtest <- sdmStats(resd.test, periods_prj, periods_fit, "tss")
save(tss.wtest, file="Results/RObjects/wtest/tss.wtest")

sens.wtest <- sdmStats(resd.test, periods_prj, periods_fit, "sens")
save(sens.wtest, file="Results/RObjects/wtest/sens.wtest")

spec.wtest <- sdmStats(resd.test, periods_prj, periods_fit, "spec")
save(spec.wtest, file="Results/RObjects/wtest/spec.wtest")

jacc.wtest <- sdmStats(rese.test, periods_prj, periods_fit, "jacc.MEAN.prob")
save(jacc.wtest, file="Results/RObjects/wtest/jacc.wtest")

jacDis.wtest <- sdmStats(rese.test, periods_prj, periods_fit, "jacc.DIS.prob")
save(jacDis.wtest, file="Results/RObjects/wtest/jacDis.wtest")

jacCal.wtest <- sdmStats(rese.test, periods_prj, periods_fit, "jacc.CAL.prob")
save(jacCal.wtest, file="Results/RObjects/wtest/jacCal.wtest")

spRichCor.wtest <- sdmStats(rese.test, periods_prj, periods_fit, "spRich.COR.prob")
save(spRichCor.wtest, file="Results/RObjects/wtest/spRichCor.wtest")

spRichRsq.wtest <- sdmStats(rese.test, periods_prj, periods_fit, "spRich.RSQ.prob")
save(spRichRsq.wtest, file="Results/RObjects/wtest/spRichRsq.wtest")

jacct.wtest <- sdmStats(rese.test, periods_prj, periods_fit, "jacc.MEAN.thre")
save(jacct.wtest, file="Results/RObjects/wtest/jacct.wtest")

jacDist.wtest <- sdmStats(rese.test, periods_prj, periods_fit, "jacc.DIS.thre")
save(jacDist.wtest, file="Results/RObjects/wtest/jacDist.wtest")

jacCalt.wtest <- sdmStats(rese.test, periods_prj, periods_fit, "jacc.CAL.thre")
save(jacCalt.wtest, file="Results/RObjects/wtest/jacCalt.wtest")

spRichCort.wtest <- sdmStats(rese.test, periods_prj, periods_fit, "spRich.COR.thre")
save(spRichCort.wtest, file="Results/RObjects/wtest/spRichCort.wtest")

spRichRsqt.wtest <- sdmStats(rese.test, periods_prj, periods_fit, "spRich.RSQ.thre")
save(spRichRsqt.wtest, file="Results/RObjects/wtest/spRichRsqt.wtest")

### Function for testing the difference betweeen SDM and CLM using a t-test
sdmStats_ModelType <- function(dat, p_prj, p_fit, measure){
  stat.Data <- matrix(NA, length(p_prj), length(p_fit))
  dat1 <- dat[which(dat$variable == measure),]

  for(l in 1:length(p_fit)){
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
      #results <- t.test(sdm$value, clm$value, paired=T)
      results <- wilcox.test(sdm$value, clm$value, paired=T)
      stat.Data[j,l] <- results$p.value
    }
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

#save(stats.modelType, file="Results/RObjects/Ttest/stats.modelType.R")
save(stats.modelType, file="Results/RObjects/wtest/stats.modelType.R")


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
  x$significant <- factor(x$significant, levels=c(0,1), labels=c("n.s.", "< 0.001"))
  x$significant2 <- factor(x$significant2, labels=c("< 0.001"))
  
  return(x)
}

## Apply the function to p-value files for each model from above
# aucStats.Models <- orgData(auc.ttest, periods_prj, "auc")
# tssStats.Models <- orgData(tss.ttest, periods_prj, "tss")
# sensStats.Models <- orgData(sens.ttest, periods_prj, "sens")
# specStats.Models <- orgData(spec.ttest, periods_prj, "spec")
# brierStats.Models <- orgData(brier.ttest, periods_prj, "brier")
# 
# jaccStats.Models <- orgData(jacc.ttest, periods_prj, "jacc.MEAN.prob")
# jacDisStats.Models <- orgData(jacDis.ttest, periods_prj, "jacc.DIS.prob")
# jacCalStats.Models <- orgData(jacCal.ttest, periods_prj, "jacc.CAL.prob")
# spRichCorStats.Models <- orgData(spRichCor.ttest, periods_prj, "spRich.COR.prob")
# spRichRsqStats.Models <- orgData(spRichRsq.ttest, periods_prj, "spRich.RSQ.prob")
# 
# jacctStats.Models <- orgData(jacct.ttest, periods_prj, "jacc.MEAN.thre")
# jacDistStats.Models <- orgData(jacDist.ttest, periods_prj, "jacc.DIS.thre")
# jacCaltStats.Models <- orgData(jacCalt.ttest, periods_prj, "jacc.CAL.thre")
# spRichCortStats.Models <- orgData(spRichCort.ttest, periods_prj, "spRich.COR.thre")
# spRichRsqtStats.Models <- orgData(spRichRsqt.ttest, periods_prj, "spRich.RSQ.thre")

aucStats.Models <- orgData(auc.wtest, periods_prj, "auc")
tssStats.Models <- orgData(tss.wtest, periods_prj, "tss")
sensStats.Models <- orgData(sens.wtest, periods_prj, "sens")
specStats.Models <- orgData(spec.wtest, periods_prj, "spec")
brierStats.Models <- orgData(brier.wtest, periods_prj, "brier")

jaccStats.Models <- orgData(jacc.wtest, periods_prj, "jacc.MEAN.prob")
jacDisStats.Models <- orgData(jacDis.wtest, periods_prj, "jacc.DIS.prob")
jacCalStats.Models <- orgData(jacCal.wtest, periods_prj, "jacc.CAL.prob")
spRichCorStats.Models <- orgData(spRichCor.wtest, periods_prj, "spRich.COR.prob")
spRichRsqStats.Models <- orgData(spRichRsq.wtest, periods_prj, "spRich.RSQ.prob")

jacctStats.Models <- orgData(jacct.wtest, periods_prj, "jacc.MEAN.thre")
jacDistStats.Models <- orgData(jacDist.wtest, periods_prj, "jacc.DIS.thre")
jacCaltStats.Models <- orgData(jacCalt.wtest, periods_prj, "jacc.CAL.thre")
spRichCortStats.Models <- orgData(spRichCort.wtest, periods_prj, "spRich.COR.thre")
spRichRsqtStats.Models <- orgData(spRichRsqt.wtest, periods_prj, "spRich.RSQ.thre")

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
  x$significant <- factor(x$significant, levels=c(0,1), labels=c("n.s.","< 0.001"))
  x$significant2 <- factor(x$significant2, labels=c("< 0.001"))
  
  return(x)
}

## Apply the function to the p-value file for model type from above
stats.ModType <- orgData(stats.modelType, periods_prj)

rm(resd.All, rese.All, sdmStats, sdmStats_ModelType, orgData)
rm(auc.ttest, tss.ttest, sens.ttest, spec.ttest, brier.ttest, jacc.ttest, jacCal.ttest, jacDis.ttest, spRichCor.ttest, spRichRsq.ttest, jacct.ttest, jacCalt.ttest, jacDist.ttest, spRichCort.ttest, spRichRsqt.ttest)
rm(auc.wtest, tss.wtest, sens.wtest, spec.wtest, brier.wtest, jacc.wtest, jacCal.wtest, jacDis.wtest, spRichCor.wtest, spRichRsq.wtest, jacct.wtest, jacCalt.wtest, jacDist.wtest, spRichCort.wtest, spRichRsqt.wtest)
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

rm(index.dColumns, index.eColumns)


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


####################################################################################
#
# Merging metric data with climatic distance and biotic distance data
#
####################################################################################
min.bioDist <- read.csv(file="Results/MinBioDist.csv")
colnames(min.bioDist) <- c("value","variable","periodFit","periodPrj")

minEucDist.melt <- read.csv("Results/MinEucDist.csv")
colnames(minEucDist.melt) <- c("value","periodPrj","periodFit")

tmpD <- resd.ModType.Mean[which(resd.ModType.Mean$periodFit %in% periods_fit),]
tmpE <- rese.ModType.Mean[which(rese.ModType.Mean$periodFit %in% periods_fit),]

## Merge metric data with climatic distance data
resd.ModType.climDist <- merge(tmpD, minEucDist.melt, by=c("periodPrj","periodFit"))
rese.ModType.climDist <- merge(tmpE, minEucDist.melt, by=c("periodPrj","periodFit"))

## Merge new datasets together so that metric data, climatic distance data and biotic distance data are together.
resd.ModType.climbioDist <- merge(resd.ModType.climDist, min.bioDist, by=c("periodPrj","periodFit"))
rese.ModType.climbioDist <- merge(rese.ModType.climDist, min.bioDist, by=c("periodPrj","periodFit"))

colnames(resd.ModType.climbioDist) <- c("periodPrj","periodFit","variable","modelType","value","sampleSize","Clim.Dist","Bio.Dist","Bio.variable")
colnames(rese.ModType.climbioDist) <- c("periodPrj","periodFit","variable","modelType","value","sampleSize","Clim.Dist","Bio.Dist","Bio.variable")

rm(resd.all, rese.all, min.bioDist, minEucDist.melt, tmpD, tmpE, resd.ModType.climDist, rese.ModType.climDist)
rm(disStats.Models, assStats.Models, stats.ModType)

save(resd.ModType.climbioDist, file="Results/RObjects/resd.ModType.climbioDist.RData")
save(rese.ModType.climbioDist, file="Results/RObjects/rese.ModType.climbioDist.RData")
save(resd.Models.Dif2, file="Results/RObjects/resd.Models.Dif2.RData")
save(rese.Models.Dif2, file="Results/RObjects/rese.Models.Dif2.RData")
save(resd.ModType.Dif2, file="Results/RObjects/resd.ModType.Dif2.RData")
save(rese.ModType.Dif2, file="Results/RObjects/rese.ModType.Dif2.RData")

