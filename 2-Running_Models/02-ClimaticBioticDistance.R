library(pdist)
library(raster)
library(reshape2)
library(paleoCLMs)

####################################################################################
#
# PCA Analysis on Occupied Grid Cells
#
####################################################################################
#combine climate data into one matrix
climAll <- do.call(rbind,clim)

#principal components analysis in R mode
pca <- prcomp(climAll,center=TRUE,scale=TRUE)

####################################################################################
#
# Distances between occupied grid cells function
# - this measures the distance between an occupied grid cell in one time period to all of the occupied cells in a second time period
# - distance can be measured 4 different ways
#     1) PCA analysis followed by minimum euclidean distance
#     2) PCA analysis followed by mean euclidean distance
#     3) Mahalanobis Distance
#     4) Standard Euclidean Distance
#
####################################################################################

## Define a function that make the calculations
climDist <- function(from, to, pcaObject, distType, distFun, returnAverage=T){

  if(distFun == "Euclidean"){
    #predict pca relationship for each time period
    from.new <- predict(pcaObject, from)
    to.new <- predict(pcaObject, to)
    
    # Calculate the distance
    d <- pdist(from.new[,1:2], to.new[,1:2])
    
    # Get the min or mean value for each cell in "from" to all the cells in "to"
    distance <- apply(as.matrix(d), MARGIN=1, FUN=distType)
  }

  if(distFun == "Mahalanobis"){
    m <- mahal(to)
    d <- predict(m, from)
    distance <- 1 - d
  }
  if(distFun=="StandEuclidean"){
    to.mean <- apply(to, 2, mean)
    to.sd <- apply(to, 2, sd)
  
    from.new <- t((t(from)-to.mean)/to.sd)
    to.new <- t((t(to)-to.mean)/to.sd)
    
    d <- pdist(from.new, to.new)
  
    distance <- apply(as.matrix(d), MARGIN=1, FUN=distType)
  }

  # Return the mean value of the previous distances
  if(returnAverage == T){
    return(mean(distance))
  }else{
    return(distance)
  }
}

## Defining periods
periodFitNum <- (c(0:10, seq(10.5, 13.5, by=0.5))*2)+1
periodPrjNum <- c(1,2,7,13,16,20,26,30,31)
 
## Selecting and Pooling data
climFit <- clim[periodFitNum]
climFit[[19]] <- poolingData(c(29:31), clim)
climFit[[20]] <- poolingData(c(32:36), clim)
climFit[[21]] <- poolingData(c(32:43), clim)

climPrj <- clim[periodPrjNum]
climPrj[[10]] <- poolingData(c(32:36), clim)
climPrj[[11]] <- poolingData(c(32:43), clim)
climPrj[[12]] <- poolingData(c(37:43), clim)

## Apply the function climDist to make the calculation across all time periods.
minEucDist <- lapply(climFit, FUN=function(x, y, z, w, j){lapply(y, FUN=climDist, x, z, w, j)}, climPrj, pca, min, "Euclidean")

## Melting and ploting these results in a HeatMap should give us the results that are comparable to the AUC, Jaccard, etc from our models (CLMs and SDMs)
minEucDist.melt <- melt(minEucDist)

## function to correct to factors
correctFactors <- function(x, pF, pP){
  x$L1 <- factor(x$L1, labels=pF)
  x$L2 <- factor(x$L2, labels=pP)
  x$L1 <- factor(x$L1, levels=rev(pF))
  x$L2 <- factor(x$L2, levels=rev(pP))
  return(x)
}

## apply the factor correction function
minEucDist.melt <- correctFactors(minEucDist.melt, periods_fit, periods_prj)
colnames(minEucDist.melt) <- c("values", "period", "period_fit")

## save distances
write.csv(minEucDist.melt, "Results/MinEucDist.csv",row.names=FALSE)

## Plotting results to double check
library(ggplot2)
forecastLine <- data.frame(x=c(0,3,5,5,6,6,7,7,8,8,9,9,10,10,11,12)+0.5, y=c(0,2,3,5,6,11,11,13,13,14,15,17,18,20,20,21)+0.5)

fitLine <- geom_line(data=forecastLine, aes(x=x, y=y), colour="black", size=1.25, alpha=0.75)
themeOpt <- theme(axis.text.x=element_text(angle=90, hjust=1)) + theme_bw()
axesLabels <- labs(x="Projecting period", y="Fitting period")
guidesLegend <- guides(fill=guide_legend(title="Legend"), colour=guide_legend(title="Legend"))
aspRatio <- coord_fixed(ratio=1)

ggplot(minEucDist.melt) + geom_tile(aes(x=period, y=period_fit, fill=values, colour=values)) + scale_fill_gradient(low="white", high="darkblue", limits=c(0, 2)) + scale_colour_gradient(low="gray90", high="black", limits=c(0, 2)) + fitLine + themeOpt + guidesLegend + axesLabels + aspRatio


####################################################################################
#
# Biotic Distance
#
####################################################################################


## Selecting and Pooling data
pollenFit <- pollen[periodFitNum]
pollenFit[[19]] <- poolingData(c(29:31), pollen)
pollenFit[[20]] <- poolingData(c(32:36), pollen)
pollenFit[[21]] <- poolingData(c(32:43), pollen)

pollenPrj <- pollen[periodPrjNum]
pollenPrj[[10]] <- poolingData(c(32:36), pollen)
pollenPrj[[11]] <- poolingData(c(32:43), pollen)
pollenPrj[[12]] <- poolingData(c(37:43), pollen)


## Calculate distance (this is done for minimum and mean)
## Define the function
bioticDistance <- function(dat1, dat2){
  require(vegan)
  require(raster)
  require(reshape2)
  
  minJaccard <- list()
  meanJaccard <- list()
  for(i in 1:nrow(dat1)){
    tmp <- vegdist(rbind(dat1[i,], dat2), method="jaccard", binary=F)
    minJaccard[[i]] <- min(as.matrix(tmp)[1,-1])
    meanJaccard[[i]] <- mean(as.matrix(tmp)[1,-1])
  }
  
  minJaccard <- mean(melt(minJaccard)$value)
  meanJaccard <- mean(melt(meanJaccard)$value)
  
  return(list(min.dissimilarity=minJaccard, mean.dissimilarity=meanJaccard))
}

## Apply the function for all the projecting periods and all the fitting periods
results.list <- lapply(pollenPrj, function(x, y){lapply(y, bioticDistance, x)}, pollenFit)

## Melt the results in a data frame
results <- melt(results.list)

## Define the factors, its levels
results$L1 <- factor(results$L1, labels=periods_prj)
results$L1 <- factor(results$L1, levels=rev(periods_prj))
results$L2 <- factor(results$L2, labels=periods_fit)
results$L2 <- factor(results$L2, levels=rev(periods_fit))
results$L3 <- factor(results$L3)


colnames(results) <- c("value","variable","period_fit","period_prj")


min.bioDist <- results[which(results$variable == "min.dissimilarity"),]
mean.bioDist <- results[which(results$variable == "mean.dissimilarity"),]

## Plot results to check
ggplot(mean.bioDist) + geom_tile(aes(x=period_prj, y=period_fit, fill=value, colour=value)) + scale_fill_gradient(low="white", high="darkblue") + scale_colour_gradient(low="gray90", high="black") + facet_grid(.~variable) + fitLine + themeOpt + guidesLegend + axesLabels + aspRatio

ggplot(min.bioDist) + geom_tile(aes(x=period_prj, y=period_fit, fill=value, colour=value)) + scale_fill_gradient(low="white", high="darkblue") + scale_colour_gradient(low="gray90", high="black") + facet_grid(.~variable) + fitLine + themeOpt + guidesLegend + axesLabels + aspRatio

write.csv(mean.bioDist, file="Results/MeanBioDist.csv",row.names=FALSE)
write.csv(min.bioDist, file="Results/MinBioDist.csv",row.names=FALSE)


plot(min.bioDist$period_prj[min.bioDist$period_fit == 0], min.bioDist$value[min.bioDist$period_fit == 0], border="white")
lines(min.bioDist$period_prj[min.bioDist$period_fit == 0], min.bioDist$value[min.bioDist$period_fit == 0])


rm(periodFitNum, periodPrjNum, climAll, pca, climDist, climFit, climPrj, pollenFit, pollenPrj, minEucDist, correctFactors, minEucDist.melt, results, results.list, min.bioDist, mean.bioDist, bioticDistance, forecastLine, fitLine, themeOpt, axesLabels, guidesLegend, aspRatio)
