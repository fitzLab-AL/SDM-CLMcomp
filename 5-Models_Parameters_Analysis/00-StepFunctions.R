################################################################################
##
## Wrapping step functions to fit models with varying variables and parameters
##
################################################################################


# GLM - SDM
glmStep <- function(commM, xData, train, iterNum){
  spNames <- colnames(commM)
  glmList <- foreach(i=1:iterNum, .verbose=T) %dopar% {
    iTrain <- train[,i]
    sub.commM <- subset(commM, iTrain)
    sub.xData <- subset(xData, iTrain)
    d=cbind(sub.commM,sub.xData)
    forms <- lapply(spNames, FUN=function(x){formula(paste(x, "~ etr_year_ave + aet_year_ave + wdei_year_ave + tmax_high_quart + prcp_low_quart + prcp_high_quart + I(etr_year_ave^2) + I(aet_year_ave^2) + I(wdei_year_ave^2) + I(tmax_high_quart^2) + I(prcp_low_quart^2) + I(prcp_high_quart^2)"))})
    glmfit <- lapply(forms, FUN=function(m){glm(m, family=binomial, data=d)})
    stepR=lapply(glmfit, FUN=function(x){step(x,direction="both")})
    forms2=lapply(stepR, FUN=function(x){x$formula})
    glmfit2=lapply(forms2, FUN=function(m){glm(m, family=binomial, data=d)})
    glmfit2
  }
  return(glmList)
}


# GAM - SDM
# There is something wrong with this function so that it doesn't work when I use %dopar% but does work when I use %do%. With %dopar% it doesn't keep track of 'd'.
gamStep <- function(commM, xData, train, iterNum){
  require(gam)
  spNames <- colnames(commM)
  gamList <- foreach(i=1:iterNum, .packages=c("gam"), .verbose=T) %dopar% {
    iTrain <- train[,i]
    sub.commM <- subset(commM, iTrain)
    sub.xData <- subset(xData, iTrain)
    d <- cbind(sub.commM, sub.xData)
    d2 <- lapply(spNames, FUN=function(x, y, z){cbind(y[,x], z)}, sub.commM, sub.xData)
    forms <- lapply(spNames, FUN=function(x){formula(paste(x,"~ etr_year_ave + aet_year_ave + wdei_year_ave + tmax_high_quart + prcp_low_quart + prcp_high_quart"))})
    gamfit <- lapply(forms, FUN=function(m){gam(m, family=binomial, data=cbind(sub.commM,sub.xData))})
    gamscope <- lapply(d2, FUN=function(x){gam.scope(x, response=1, smoother="s",arg=c("df=1","df=2","df=3","df=4","df=5"), form=F)})
    stepR <- mapply(FUN=function(j, l){step.gam(j, scope=l, direction="both")}, j=gamfit, l=gamscope, SIMPLIFY=F)
    forms2 <- lapply(stepR, FUN=function(x){x$formula})
    gamfit2 <- lapply(forms2, FUN=function(m){gam(m, family=binomial, data=d)})
    gamfit2
  }
  return(gamList)
}


# CQO - CLM
cqoStep <- function(commM, xData, train, iterNum){
  require(VGAM)
  require(gtools)
  cqoList <- foreach(i=1:numIter, .packages=c("VGAM"), .verbose=T, .errorhandling='pass') %dopar% {
    iTrain <- train[,i]
    sub.commM <- subset(commM, iTrain)
    sub.xData <- subset(xData, iTrain)
    
    vars.all <- c("etr_year_ave","aet_year_ave","wdei_year_ave","tmax_high_quart","prcp_low_quart","prcp_high_quart")
    form <- formula(sub.commM ~ etr_year_ave + aet_year_ave + wdei_year_ave + tmax_high_quart + prcp_low_quart + prcp_high_quart)
    cqo.obj <- cqo(form, family=binomialff(multiple.responses=T), data=sub.xData, trace=TRUE)
    
    cqo.aic <- AIC(cqo.obj)
    
    mods.rec <- list()
    aic.rec <- list()
    delta.rec <- list()
    min.delta <- -3
    i <- 0
    
    while(min.delta < -2){
      # increase index at each iteration
      i <- i + 1
      
      # get variables used in the actual model
      cqo.terms <- terms(cqo.obj)
      vars <- attr(cqo.terms, "term.labels")
      start.terms <- vars.all %in% vars
      
      # create the candidate formulas, changing one variable at a time (adding it if it was absent and substracting if it was present)
      candidates <- abs(t(diag(length(vars.all)) - start.terms))
      candidates <- apply(candidates, 2, as.logical)
      candidates.vars <- lapply(1:nrow(candidates), FUN=function(x, y, z){y[z[x,]]}, vars.all, candidates)
      candidates.forms <- lapply(candidates.vars, FUN=function(x){formula(paste("sub.commM ~", paste(x, collapse="+")))})
      
      # fit the candidate models
      candidates.mods <- lapply(candidates.forms, FUN=function(x, y){cqo(x, family=binomialff(multiple.responses=TRUE), data=y, trace=FALSE)}, sub.xData)
      
      # calculate candidates aic and deltas
      candidates.aic <- sapply(candidates.mods, FUN=function(x){AIC(x)})
      candidates.delta <- candidates.aic - cqo.aic
      min.aic.i <- which.min(candidates.delta)
      min.delta <- candidates.delta[min.aic.i]
      
      # substitute actual model by the winning candidate
      if(candidates.delta[min.aic.i] < -2){
        cqo.obj <- candidates.mods[[min.aic.i]]
        cqo.aic <- AIC(cqo.obj)
      }
      
      # store information in storing variables
      mods.rec[[i]] <- candidates.mods
      aic.rec[[i]] <- candidates.aic
      delta.rec[[i]] <- candidates.delta
      
      # print iteration info for tracking purpouses
      print(paste("ITERATION: ", i, sep=""))
      print(paste("AIC: ", round(cqo.aic, digit=2)))
      print(paste("AIC cand: ", paste(round(candidates.aic, digit=2), collapse=", ")))
      print(paste("AIC deltas: ", paste(round(candidates.delta, digit=2), collapse=", ")))
    }
    
    cqo.obj
  }
  return(cqoList)
}


# CAO - CLM
caoStep <- function(commM, xData, train, iterNum){
  #commM <- pollenList[[a]]
  #xData <- climList[[a]]
  #train <- trainSets[[k]]
  #iterNum <- 10
  require(VGAM)
  
  caoList <- foreach(ii=1:numIter, .packages=c("VGAM"), .verbose=T, .errorhandling='pass') %dopar% {
    iTrain <- train[,ii]
    sub.commM <- subset(commM, iTrain)
    sub.xData <- subset(xData, iTrain)
    
    vars.all <- colnames(sub.xData)
    df.candidates <- c(0.5, 1.5, 2.5, 3.5, 4.5)
    species <- colnames(sub.commM)
    
    form <- formula(paste("sub.commM", paste(vars.all, collapse="+"), sep="~"))
    cao.obj <- cao(form, family=binomialff(multiple.responses=T), data=sub.xData, trace=FALSE)
    cao.aic <- AIC(cao.obj)
    
    min.delta <- -3
    i <- 0
    while(min.delta < -2){
      # increase index at each iteration
      i <- i + 1
      
      # print iteration info for tracking purpouses
      print(paste("ITERATION: ", i, sep=""))
      print(paste("AIC: ", round(cao.aic, digit=2)))
      
      # get variables used in the actual model
      cao.terms <- terms(cao.obj)
      vars <- attr(cao.terms, "term.labels")
      start.terms <- vars.all %in% vars
      
      # create the candidate formulas, changing one variable at a time (adding it if it was absent and substracting if it was present)
      candidates <- abs(t(diag(length(vars.all)) - start.terms))
      candidates <- apply(candidates, 2, as.logical)
      candidates.vars <- lapply(1:nrow(candidates), FUN=function(x, y, z){y[z[x,]]}, vars.all, candidates)
      candidates.forms <- lapply(candidates.vars, FUN=function(x){formula(paste("sub.commM ~", paste(x, collapse="+")))})
      
      # fit the candidate models
      candidates.mods <- lapply(candidates.forms, FUN=function(x, y){cao(x, family=binomialff(multiple.responses=TRUE), data=y, trace=FALSE)}, sub.xData)
      
      # calculate candidates aic and deltas
      candidates.aic <- sapply(candidates.mods, AIC)
      candidates.delta <- candidates.aic - cao.aic
      
      
      if(is.na(cao.aic)){
        if(all(is.na(candidates.aic))){
          min.aic.i <- sample(1:length(candidates.aic), 1)
        }else{
          min.aic.i <- which.min(candidates.aic)
        }
      }else{
        if(all(is.na(candidates.aic))){
          min.delta <- 0
        }else{
          min.aic.i <- which.min(candidates.delta)
          min.delta <- candidates.delta[min.aic.i]
        }
      }
      
      # print output
      print(paste("AIC cand: ", paste(round(candidates.aic, digit=2), collapse=", ")))
      print(paste("AIC deltas: ", paste(round(candidates.delta, digit=2), collapse=", ")))
      
      # substitute actual model by the winning candidate
      if(min.delta < -2){
        cao.obj <- candidates.mods[[min.aic.i]]
        cao.aic <- AIC(cao.obj)
      }
    }
    
    # get variables used in the actual model
    cao.terms <- terms(cao.obj)
    vars <- attr(cao.terms, "term.labels")
    form <- formula(paste("sub.commM", paste(vars, collapse="+"), sep="~"))
    
    species.df <- c()
    for(j in 1:length(species)){
      # remove previous temporary objects
      rm(candidate.mod, candidate.aic)
      
      # print info
      cat("Species: ", species[j], "\n", sep="")
      cat("  Iteration: 0 (DF candidate: 0.5)\n")
      
      # fit initial model with 0.5 degrees of freedom
      cao.obj <- cao(form, family=binomialff(multiple.responses=TRUE), data=sub.xData, df1.nl=eval(parse(text=paste("c(", species[j], "=", 0.5, ", 2.5)", sep=""))), trace=FALSE)
      
      # calculate initial AIC and print info
      cao.aic <- AIC(cao.obj)
      cat("    AIC: ", round(cao.aic, digit=2), "\n", sep="")
      
      # set up initial values for the search of optimal degrees of freedom
      species.df[j] <- 0.5
      candidate.delta <- -3
      i <- 0
      
      # while loop to iterate and search optimal number of degrees of freedom
      while(candidate.delta < -2){
        
        # increase index at each iteration
        i <- i + 1
        
        # stop the loop in the last iteration
        if(i == 5){break}
        
        # pick up the degree of freedom to be tested
        df.test <- df.candidates[i + 1]
        
        # print iteration info for tracking purpouses
        cat("  Iteration: ", i, " (DF candidate: ", df.test, ")\n", sep="")
        cat("    AIC: ", round(cao.aic, digit=2), "\n", sep="")
        
        # fit the candidate model
        df.cand <- paste("c(", species[j], "=", df.test, ", 2.5)", sep="")
        try(candidate.mod <- cao(form, family=binomialff(multiple.responses=TRUE), data=sub.xData, df1.nl=eval(parse(text=df.cand)), trace=FALSE), silent=TRUE)
        
        # if the model didn't work skip to the next degree of freedom
        if(is.null(candidate.mod)){next}
        
        # calculate candidates aic and deltas
        # if try(cao) didn't work, the candidate.mod will be the same as in previous iteration
        # hence, the aic will be same and delta will be 0. Taking the last valid model as the right one.
        candidate.aic <- AIC(candidate.mod)
        
        # print output
        cat("    AIC cand: ", round(candidate.aic, digit=2), sep="")
        
        # if the aic of the candidate is NA skip to the next degree of freedom
        if(is.na(candidate.aic)){
          cat("\n")
          next
          # if aic is not NA calculate delta and update model and aic if necessary
        }else{
          candidate.delta <- candidate.aic - cao.aic
          cat(" (delta=", round(candidate.delta, digit=2), ")\n", sep="")
          
          # substitute actual model by the winning candidate
          if(candidate.delta < -2){
            cao.obj <- candidate.mod
            cao.aic <- candidate.aic
            species.df[j] <- df.test
          }
        }
      }
    }
    
    # create object with the degrees of freedom for each species in the cao format
    df.mod <- paste(species, species.df, sep="=")
    df.mod <- paste(df.mod, collapse=",")
    df.mod <- paste("c(", df.mod, ")", sep="")
    
    # fit final model with optimal degrees of freedom for each species
    cao.obj <- cao(form, family=binomialff(multiple.responses=TRUE), data=sub.xData, df1.nl=eval(parse(text=df.mod)), trace=F)
    
    # return the final model
    cao.obj
  }
  return(caoList)
}


smarsStep <- function(commM, xData, train, iterNum, par){
  require(earth)
  require(caret)
  taxonNames <- colnames(commM)
  marsList.sp <- foreach(i=1:numIter, .packages=c("earth"), .verbose=T) %do% {
    iTrain <- train[,i]
    sub.commM <- subset(commM, iTrain)
    sub.xData <- subset(xData, iTrain)
    fitControl <- trainControl(method = "repeatedcv", number=10, repeats = 50)
    marsGrid <- expand.grid(nprune=c(1, 2, 5, 7, 9, 10),degree=c(1,2,3,4,5))
    marsfits <- lapply(taxonNames,FUN=function(t){train(x=sub.xData, y=sub.commM[,t], method = "earth", metric="Rsquared", trControl=fitControl, tuneGrid=marsGrid,trace=T,glm=list(family=binomial))})
    nprune <- lapply(marsfits, FUN=function(x){x$bestTune$nprune})
    degree <- lapply(marsfits, FUN=function(x){x$bestTune$degree})
    marsfit.sp <- mapply(taxonNames, FUN=function(m, n){earth(x=sub.xData, y=sub.commM[,m], nprune=n, degree=1, trace=T, glm=list(family=binomial))}, nprune)
    marsfit.sp=list()
    for (j in 1:length(taxonNames)){
      marsfit.sp[[j]]=earth(x=sub.xData, y=sub.commM[,j], nprune=nprune[[j]], degree=degree[[j]], trace=T, glm=list(family=binomial))
    }
    marsfit.sp
  }
  return(marsList.sp)
}



# MARS - CLM 
marsStep <- function(commM, xData, train, iterNum){
  require(earth)
  marsList <- foreach(i=1:numIter, .packages=c("earth"), .verbose=T) %dopar% {
    iTrain <- train[,i]
    sub.commM <- subset(commM, iTrain)
    sub.xData <- subset(xData, iTrain)
    marsGrid <- expand.grid(nprune=c(1, 2, 5, 7, 9, 10),degree=c(1,2,3,4,5))
    marsfits=list()
    for (j in 1:nrow(marsGrid)){
      marsfits[[j]] <- earth(x=sub.xData, y=sub.commM, nprune=marsGrid[j,1], degree=marsGrid[j,2], trace=T, glm=list(family=binomial))
    }
    rsquared=lapply(marsfits,FUN=function(x){x$rsq})
    x=which.max(rsquared)
    marsfits[[x]]
  }
  return(marsList)
}



# Nueral Network - SDM
snnStep <- function(commM, xData, train, iterNum, par){
  require(nnet)
  spNames <- colnames(commM)
  netList <- foreach(i=1:numIter, .packages=c("nnet"), .verbose=T) %do% {
    iTrain <- train[,i]
    sub.commM <- subset(commM, iTrain)
    sub.xData <- subset(xData, iTrain)
    fitControl <- trainControl(method = "repeatedcv", number=10, repeats = 50)
    decayseq <- 10^seq(-3, 0, length = 10)
    nnetGrid <- expand.grid(size=20,decay=decayseq)
    snnetfits <-lapply(taxonNames, FUN=function(t){train(x=sub.xData, y=sub.commM[,t], method = "nnet", metric="Rsquared", maxit=2000, skip=F, trControl=fitControl, tuneGrid=nnetGrid, verbose=FALSE)})
    snnetfit=list()
    for (j in 1:length(snnetfits)){
      taxon=taxonNames[[j]]
      tempsnnetfit=nnet(x=sub.xData,y=sub.commM[,taxon], size=snnetfits[[j]]$bestTune$size, decay=snnetfits[[j]]$bestTune$decay,lineout=T)
      snnetfit[[j]]=tempsnnetfit
    }
    snnetfit
  }    
  return(netList)
}



# Nueral Network - CLM
mnnStep <- function(commM, xData, train, iterNum){
  require(nnet)
  netList <- foreach(i=1:numIter, .packages=c("nnet"), .verbose=T) %dopar% {
    iTrain <- train[,i]
    sub.commM <- subset(commM, iTrain)
    sub.xData <- subset(xData, iTrain)
    fitControl <- trainControl(method = "repeatedcv", number=10, repeats = 50)
    decayseq <- 10^seq(-3, 0, length = 10)
    nnetGrid <- expand.grid(size=20,decay=decayseq)
    form=formula(Ambrosia.type+Artemisia+Alnus+Corylus+Ostrya.Carpinus+Betula+Quercus+Carya+Juglans+Fraxinus+Abies+Larix+Picea+Platanus+Populus+Salix+Acer+Ulmus+Pinus~etr_year_ave+aet_year_ave+wdei_year_ave+tmax_high_quart+prcp_low_quart+prcp_high_quart)
    nnetfits <- train(form, data=cbind(sub.commM,sub.xData), method = "nnet", maxit=2000, skip=F, metric="Rsquared", trControl=fitControl, tuneGrid=nnetGrid, verbose=FALSE)
    netfit <- nnet(x=sub.xData, y=sub.commM, size=nnetfits$bestTune$size, decay=nnetfits$bestTune$decay, linout=T)
    netfit
  }
  return(netList)
}

#For CARTS we do not tune anything and the model itself picks the variables

# CARTS - SDM
scartsStep <- function(commM, xData, train, iterNum){
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
cartsStep <- function(commM, xData, train, iterNum){
  require(mvpart)
  cartsList <- foreach(i=1:numIter, .packages=c("mvpart"), .verbose=T) %dopar% {
    iTrain <- train[,i]
    sub.commM <- subset(commM, iTrain)
    sub.xData <- subset(xData, iTrain)
    form <- formula(sub.commM ~ etr_year_ave + aet_year_ave + wdei_year_ave + tmax_high_quart + prcp_low_quart + prcp_high_quart)
    cartsfit <- mvpart(form, cbind(sub.commM, sub.xData),plot.add=FALSE)
    cartsfit
  }
  return(cartsList)
}

