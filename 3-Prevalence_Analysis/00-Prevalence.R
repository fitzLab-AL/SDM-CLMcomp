library(reshape2)
library(ggplot2)

load("SDMvsCLM/FINAL/Objects/DisResults.RData")
load("FINAL/Objects/AssResults.RData")

periodFit = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "11.5", "12", "12.5", "13", "13.5", "14-15", "15.5-17.5", "15.5-21")
periodPrj=c("0", "0.5", "3", "6", "7.5", "9.5", "12.5", "14.5", "15", "15.5-17.5", "18-21", "15.5-21")

resd.all <- melt(resd.All, id.vars=c("taxon","variable","value", "iteration","model","periodFit","periodPrj","test"))
rese.all <- melt(rese.All, id.vars=c("variable","value","iteration","model","periodFit","periodPrj","test"))



resd.all=resd.all[,-9]
rese.all=rese.all[,-8]
resd.all$value[which(resd.all$value == -9999)] <- NA
rese.all$value[which(rese.all$value == -9999)] <- NA
#colnames(resd.all) <- c("taxon","variable","iteration","model","period_fit","period","test","value")
#colnames(rese.all) <- c("variable","iteration","model","period_fit","period","test","value")


resd.all$periodFit <- factor(resd.all$periodFit, rev(periodFit))
resd.all$periodPrj <- factor(resd.all$periodPrj)
rese.all$periodFit <- factor(rese.all$periodFit, rev(periodFit))
rese.all$periodPrj <- factor(rese.all$periodPrj)

newModelTypeColumn <- function(d){
  d$model_type <- NA
  d$model_type[which(d$model %in% c("cqo","cao","carts","mars","mnn"))] <- "CLM"
  d$model_type[which(d$model %in% c("glm","gam","scarts","smars","snn"))] <- "SDM"
  d$model_type <- factor(d$model_type)
  return(d)
}

resd.all=newModelTypeColumn(resd.all)
rese.all=newModelTypeColumn(rese.all)

## Calculate average prevalence for each species at the time period the model was fitted
prev.train <- resd.all[which(resd.all$test == "train" & resd.all$variable == "prev"),]
prev.train <- acast(prev.train[,-which(colnames(prev.train) %in% c("test","model_type"))], taxon ~ periodFit, mean, na.rm=T)
prev.train <- melt(prev.train)
colnames(prev.train) <- c("taxon","periodFit","prevFit")

## Extract distribution assessment and reshape them, similarly to how we did before
resd.Prev <- acast(resd.all[which(resd.all$test == "test"),], taxon ~ variable ~ model ~ model_type ~ periodPrj ~ periodFit, mean, na.rm=T)
resd.Prev <- melt(resd.Prev)
resd.Prev <- na.omit(resd.Prev)
colnames(resd.Prev) <- c("taxon","variable","model","model_type","periodPrj","periodFit","value")
resd.Prev <- resd.Prev[-which(resd.Prev$variable %in% c("prev","trmin")),]

p <- merge(resd.Prev, prev.train, by=c("taxon","periodFit"))

# AUC
auc=p[which(p$variable=="auc"),]

prev.res=list()
taxa=unique(resd.all$taxon)
for (i in 1:length(taxa)){
  taxon=taxa[i]
  sp=auc[which(auc$taxon==taxon),]
  sp1=matrix(NA,20,3)
  for (j in 1:length(periodFit)){
    age=sp[which(sp$periodFit==periodFit[j]),]
    SDMage=age[which(age$model_type=="SDM"),]
    CLMage=age[which(age$model_type=="CLM"),]
    x=mean(SDMage$value)-mean(CLMage$value)
    y=mean(age$prevFit)
    sp1[j,1]=x
    sp1[j,2]=y
    sp1[j,3]=as.character(unique(age$periodFit))
  }
  prev.res[[i]]=sp1
}


for (i in 1:length(prev.res)){
  colnames(prev.res[[i]])=c("Diff","Prev","PeriodFit")  
}
names(prev.res)=taxa
for (i in 1:length(prev.res)){
  prev.res[[i]]=data.frame(prev.res[[i]])
}

prev.res.melt=melt(prev.res,id.vars=c("Diff","Prev","PeriodFit"))
colnames(prev.res.melt)=c("Diff","Prev","PeriodFit","Taxa")

prev.res.melt$Prev=as.numeric(levels(prev.res.melt$Prev))[prev.res.melt$Prev]
prev.res.melt$Diff=as.numeric(levels(prev.res.melt$Diff))[prev.res.melt$Diff]
prev.res.melt$Period_Fit=factor(prev.res.melt$PeriodFit,periodFit)

forecastLine <- data.frame(x=c(0,0.25,0.5,0.75,1), y=c(0,0,0,0,0))
fitLine <- geom_line(data=forecastLine, aes(x=x, y=y), colour="black", size=1.3, alpha=0.75)


ggplot(prev.res.melt, aes(x=Prev,y=Diff)) + geom_point() + facet_wrap(~Taxa,ncol=5) + fitLine + stat_smooth(method="lm", size=1.25, fullrange = TRUE, color="red")


ggplot(prev.res.melt, aes(x=Prev,y=Diff,col=Taxa)) + geom_point() + facet_wrap(~Period_Fit,ncol=5) + fitLine + stat_smooth(method="lm", size=1.25, fullrange = TRUE, color="red")
*********************************************
# Brier Score
brier=p[which(p$variable=="brier"),]

prev.res=list()
taxa=unique(resd.all$taxon)
for (i in 1:length(taxa)){
  taxon=taxa[i]
  sp=brier[which(brier$taxon==taxon),]
  sp1=matrix(NA,20,3)
  for (j in 1:length(periodFit)){
    age=sp[which(sp$periodFit==periodFit[j]),]
    SDMage=age[which(age$model_type=="SDM"),]
    CLMage=age[which(age$model_type=="CLM"),]
    x=mean(SDMage$value)-mean(CLMage$value)
    y=mean(age$prevFit)
    sp1[j,1]=x
    sp1[j,2]=y
    sp1[j,3]=as.character(unique(age$periodFit))
  }
  prev.res[[i]]=sp1
}


for (i in 1:length(prev.res)){
  colnames(prev.res[[i]])=c("Diff","Prev","PeriodFit")  
}
names(prev.res)=taxa
for (i in 1:length(prev.res)){
  prev.res[[i]]=data.frame(prev.res[[i]])
}

prev.res.melt=melt(prev.res,id.vars=c("Diff","Prev","PeriodFit"))
colnames(prev.res.melt)=c("Diff","Prev","PeriodFit","Taxa")

prev.res.melt$Prev=as.numeric(levels(prev.res.melt$Prev))[prev.res.melt$Prev]
prev.res.melt$Diff=as.numeric(levels(prev.res.melt$Diff))[prev.res.melt$Diff]
prev.res.melt$Period_Fit=factor(prev.res.melt$PeriodFit,periodFit)

forecastLine <- data.frame(x=c(0,0.25,0.5,0.75,1), y=c(0,0,0,0,0))
fitLine <- geom_line(data=forecastLine, aes(x=x, y=y), colour="black", size=1.3, alpha=0.75)


ggplot(prev.res.melt, aes(x=Prev,y=Diff)) + geom_point() + facet_wrap(~Taxa,ncol=5) + fitLine + stat_smooth(method="lm", size=1.25, fullrange = TRUE, color="red")


ggplot(prev.res.melt, aes(x=Prev,y=Diff,col=Taxa)) + geom_point() + facet_wrap(~Period_Fit,ncol=5) + fitLine + stat_smooth(method="lm", size=1.25, fullrange = TRUE, color="red")
*********************************************
