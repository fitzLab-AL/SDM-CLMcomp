##Load plotting libraries
library(reshape2)
library(ggplot2)
library(scales)

################################################################################
#
# TRAINING AND TESTING DATA FOR THE FITTING TIME PERIOD
#
################################################################################


## AUC
d.auc <- resd.train[which(resd.train$variable == "auc"),]

pdf("Results/Graphs/Train/auc-train.pdf", 4, 4)
    ggplot(data=d.auc, aes(modelType, value)) +
    geom_boxplot(aes(fill=test), width=1, notch=T, outlier.shape=NA) +
    scale_fill_manual(name="Legend", values=c("grey70", "grey30")) +
    coord_cartesian(ylim=c(0.5, 1)) +
    theme_bw() +
    theme(legend.justification=c(0, 0), legend.position=c(0, 0)) +
    labs(title="", x="", y="AUC")
dev.off()

pdf("Results/Graphs/Train/auc-train-models.pdf", 6, 4)
    ggplot(data=d.auc, aes(model, value)) +
    geom_boxplot(aes(fill=test), width=1, notch=T, outlier.shape=NA) +
    scale_fill_manual(name="Legend", values=c("grey70", "grey30")) +
    coord_cartesian(ylim=c(0.5, 1)) +
    theme_bw() +
    theme(legend.justification=c(0, 0), legend.position=c(0, 0)) +
    labs(title="", x="", y="AUC")
dev.off()


## BRIER
d.brier <- resd.train[which(resd.train$variable == "brier"),]

pdf("Results/Graphs/Train/brier-train.pdf", 4, 4)
    ggplot(data=d.brier, aes(modelType, 1-value)) +
    geom_boxplot(aes(fill=test), width=1, notch=T, outlier.shape=NA) +
    scale_fill_manual(name="Legend", values=c("grey70", "grey30")) +
    coord_cartesian(ylim=c(0.6, 1)) +
    theme_bw() +
    theme(legend.justification=c(0, 0), legend.position=c(0, 0)) +
    labs(title="", x="", y="1 - Brier score")
dev.off()

pdf("Results/Graphs/Train/brier-train-models.pdf", 6, 4)
    ggplot(data=d.brier, aes(model, 1-value)) +
    geom_boxplot(aes(fill=test), width=1, notch=T, outlier.shape=NA) +
    scale_fill_manual(name="Legend", values=c("grey70", "grey30")) +
    coord_cartesian(ylim=c(0.6, 1)) +
    theme_bw() +
    theme(legend.justification=c(0, 0), legend.position=c(0, 0)) +
    labs(title="", x="", y="1 - Brier score")
dev.off()


## Jaccard
e.jacc <- rese.train[which(rese.train$variable == "jacc.MEAN.thre"),]

pdf("Results/Graphs/Train/jac-thre-train.pdf", 4, 4)
    ggplot(data=e.jacc, aes(modelType, 1-value)) +
    geom_boxplot(aes(fill=test), width=1, notch=T, outlier.shape=NA) +
    labs(title="", x="", y="Jaccard similarity index") +
    scale_fill_manual(name="Legend", values=c("grey70", "grey30")) +
    coord_cartesian(ylim=c(0, 1)) +
    theme_bw() +
    theme(legend.justification=c(1, 0), legend.position=c(1, 0))
dev.off()

pdf("Results/Graphs/Train/jac-thre-train-models.pdf", 6, 4)
    ggplot(data=e.jacc, aes(model, 1-value)) +
    geom_boxplot(aes(fill=test), width=1, notch=T, outlier.shape=NA) +
    labs(title="", x="", y="Jaccard similarity index") +
    scale_fill_manual(name="Legend", values=c("grey70", "grey30")) +
    coord_cartesian(ylim=c(0, 1)) +
    theme_bw() +
    theme(legend.justification=c(1, 0), legend.position=c(1, 0))
dev.off()

## Jaccard
e.jacc <- rese.train[which(rese.train$variable == "jacc.MEAN.prob"),]

pdf("Results/Graphs/Train/jac-prob-train.pdf", 4, 4)
    ggplot(data=e.jacc, aes(modelType, 1-value)) +
    geom_boxplot(aes(fill=test), width=1, notch=T, outlier.shape=NA) +
    labs(title="", x="", y="Jaccard similarity index") +
    scale_fill_manual(name="Legend", values=c("grey70", "grey30")) +
    coord_cartesian(ylim=c(0, 1)) +
    theme_bw() +
    theme(legend.justification=c(1, 0), legend.position=c(1, 0))
dev.off()

pdf("Results/Graphs/Train/jac-prob-train-models.pdf", 6, 4)
    ggplot(data=e.jacc, aes(model, 1-value)) +
    geom_boxplot(aes(fill=test), width=1, notch=T, outlier.shape=NA) +
    labs(title="", x="", y="Jaccard similarity index") +
    scale_fill_manual(name="Legend", values=c("grey70", "grey30")) +
    coord_cartesian(ylim=c(0, 1)) +
    theme_bw() +
    theme(legend.justification=c(1, 0), legend.position=c(1, 0))
dev.off()

## Jaccard Calibration
e.jacC <- rese.train[which(rese.train$variable == "jacc.CAL.prob"),]

pdf("Results/Graphs/Train/jac-cal-prob-train.pdf", 4, 4)
    ggplot(data=e.jacC, aes(modelType, value)) +
    geom_boxplot(aes(fill=test), width=1, notch=T, outlier.shape=NA) +
    labs(title="", x="", y="Jaccard calibration\n(model efficiency)") +
    scale_fill_manual(name="Legend", values=c("grey70", "grey30")) +
    coord_cartesian(ylim=c(-10, 1)) +
    theme_bw() +
    theme(legend.justification=c(0,0), legend.position=c(0,0))
dev.off()

pdf("Results/Graphs/Train/jac-cal-prob-train-models.pdf", 6, 4)
    ggplot(data=e.jacC, aes(model, value)) +
    geom_boxplot(aes(fill=test), width=1, notch=T, outlier.shape=NA) +
    labs(title="", x="", y="Jaccard calibration\n(model efficiency)") +
    scale_fill_manual(name="Legend", values=c("grey70", "grey30")) +
    coord_cartesian(ylim=c(-10, 1)) +
    theme_bw() +
    theme(legend.justification=c(0,0), legend.position=c(0,0))
dev.off()


## Jaccard Discrimination
e.jacD <- rese.train[which(rese.train$variable == "jacc.DIS.prob"),]

pdf("Results/Graphs/Train/jac-dis-prob-train.pdf", 4, 4)
    ggplot(data=e.jacD, aes(modelType, value)) +
    geom_boxplot(aes(fill=test), width=1, notch=T, outlier.shape=NA) +
    labs(title="", x="", y="Jaccard discrimination\n(correlation)") +
    scale_fill_manual(name="Legend", values=c("grey70", "grey30")) +
    coord_cartesian(ylim=c(0, 1)) +
    theme_bw() +
    theme(legend.justification=c(0, 0), legend.position=c(0, 0))
dev.off()

pdf("Results/Graphs/Train/jac-dis-prob-train-models.pdf", 6, 4)
    ggplot(data=e.jacD, aes(model, value)) +
    geom_boxplot(aes(fill=test), width=1, notch=T, outlier.shape=NA) +
    labs(title="", x="", y="Jaccard discrimination\n(correlation)") +
    scale_fill_manual(name="Legend", values=c("grey70", "grey30")) +
    coord_cartesian(ylim=c(0, 1)) +
    theme_bw() +
    theme(legend.justification=c(0, 0), legend.position=c(0, 0))
dev.off()


## Species richness correlation
e.spRCor <- rese.train[which(rese.train$variable == "spRich.COR.prob"),]

pdf("Results/Graphs/Train/spr-dis-train.pdf", 4, 4)
    ggplot(data=e.spRCor, aes(modelType, value)) +
    geom_boxplot(aes(fill=test), width=1, notch=T, outlier.shape=NA) +
    labs(title="", x="", y="Species richness discrimination\n(correlation)") +
    scale_fill_manual(name="Legend", values=c("grey70", "grey30")) +
    coord_cartesian(ylim=c(0, 1)) +
    theme_bw() +
    theme(legend.justification=c(0,0), legend.position=c(0,0))
dev.off()

pdf("Results/Graphs/Train/spr-dis-train-models.pdf", 6, 4)
    ggplot(data=e.spRCor, aes(model, value)) +
    geom_boxplot(aes(fill=test), width=1, notch=T, outlier.shape=NA) +
    labs(title="", x="", y="Species richness discrimination\n(correlation)") +
    scale_fill_manual(name="Legend", values=c("grey70", "grey30")) +
    coord_cartesian(ylim=c(0, 1)) +
    theme_bw() +
    theme(legend.justification=c(0,0), legend.position=c(0,0))
dev.off()


## Species richness Rsquared
e.spRRsq <- rese.train[which(rese.train$variable == "spRich.RSQ.prob"),]

pdf("Results/Graphs/Train/spr-cal-train.pdf", 4, 4)
    ggplot(data=e.spRRsq, aes(modelType, value)) +
    geom_boxplot(aes(fill=test), width=1, notch=T, outlier.shape=NA) +
    labs(title="", x="", y="Species richness calibration\n(model efficiency)") +
    scale_fill_manual(name="Legend", values=c("grey70", "grey30")) +
    coord_cartesian(ylim=c(-0.5, 1)) +
    theme_bw() +
    theme(legend.justification=c(0,0), legend.position=c(0,0))
dev.off()

pdf("Results/Graphs/Train/spr-cal-train-models.pdf", 6, 4)
    ggplot(data=e.spRRsq, aes(model, value)) +
    geom_boxplot(aes(fill=test), width=1, notch=T, outlier.shape=NA) +
    labs(title="", x="", y="Species richness calibration\n(model efficiency)") +
    scale_fill_manual(name="Legend", values=c("grey70", "grey30")) +
    coord_cartesian(ylim=c(-0.5, 1)) +
    theme_bw() +
    theme(legend.justification=c(0,0), legend.position=c(0,0))
dev.off()


rm(d.auc, d.tss, d.brier, e.jacc, e.jacC, e.jacD, e.spRCor, e.spRRsq)


################################################################################
##
## PLOTTING RESULTS
##
################################################################################

forecastLine <- data.frame(x=c(0,3,5,5,6,6,7,7,8,8,9,9,10,10,11,12)+0.5, y=c(0,2,3,5,6,10,10,12,12,13,14,16,17,19,19,20)+0.5)

# These plotting functions can be modified to reflect the fixed scales (min and max values) and use more attractive colors.

fitLine <- geom_line(data=forecastLine, aes(x=x, y=y), colour="black", size=1.3, alpha=0.75)

themeOpt <- theme(axis.text.x=element_text(angle=90, hjust=1)) + theme_bw()

axesLabels <- labs(x="Projecting period", y="Fitting period")

aspRatio <- coord_fixed(ratio=1)

#################################################################################
## CLM versus SDM comparison
#################################################################################

# MEAN DISTRIBUTION VALUES
pdf("Results/Graphs/acrossTime-Means/auc.pdf", 6.5, 5)
guidesLegend <- guides(fill=guide_legend(title="Mean AUC", reverse=T), colour=guide_legend(title="Mean AUC", reverse=T))
ggplot(resd.ModType.Mean[which(resd.ModType.Mean$variable == "auc"),]) +
      geom_tile(aes(x=periodPrj, y=periodFit, fill=value, colour=value)) +
      scale_fill_gradient2(low="darkorange", high="darkgreen", mid="white", midpoint=0.75, limits=c(0.5, 1)) +
      scale_colour_gradient(low="gray90", high="black", limits=c(0.5, 1)) +
      facet_grid(. ~ modelType) +
      themeOpt + guidesLegend + axesLabels + aspRatio + fitLine 
dev.off()

pdf("Results/Graphs/acrossTime-Means/sens.pdf", 6.5, 5)
guidesLegend <- guides(fill=guide_legend(title="Mean Sens.", reverse=T), colour=guide_legend(title="Mean Sens.", reverse=T))
ggplot(resd.ModType.Mean[which(resd.ModType.Mean$variable == "sens"),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=value, colour=value)) +
    scale_fill_gradient2(low="darkorange", high="darkgreen", mid="white", midpoint=0.5, limits=c(0, 1)) +
    scale_colour_gradient(low="gray90", high="black", limits=c(0, 1)) +
    facet_grid(. ~ modelType) +
    themeOpt + guidesLegend + axesLabels + aspRatio + fitLine 
dev.off()

pdf("Results/Graphs/acrossTime-Means/spec.pdf", 6.5, 5)
guidesLegend <- guides(fill=guide_legend(title="Mean Spec.", reverse=T), colour=guide_legend(title="Mean Spec.", reverse=T))
ggplot(resd.ModType.Mean[which(resd.ModType.Mean$variable == "spec"),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=value, colour=value)) +
    scale_fill_gradient2(low="darkorange", high="darkgreen", mid="white", midpoint=0.5, limits=c(0, 1)) +
    scale_colour_gradient(low="gray90", high="black", limits=c(0, 1)) +
    facet_grid(. ~ modelType) +
    themeOpt + guidesLegend + axesLabels + aspRatio + fitLine 
dev.off()

pdf("Results/Graphs/acrossTime-Means/brier.pdf", 6.5, 5)
guidesLegend <- guides(fill=guide_legend(title="Mean 1-Brier", reverse=T), colour=guide_legend(title="Mean 1-Brier", reverse=T))
ggplot(resd.ModType.Mean[which(resd.ModType.Mean$variable == "brier"),]) +
      geom_tile(aes(x=periodPrj, y=periodFit, fill=1-value, colour=1-value)) +
      scale_fill_gradient2(low="darkorange", high="darkgreen", mid="white", midpoint=0.5, limits=c(0, 1)) +
      scale_colour_gradient(low="gray90", high="black", limits=c(0, 1)) +
      facet_grid(. ~ modelType) +
      themeOpt + guidesLegend + axesLabels + aspRatio + fitLine 
dev.off()


# MEAN ASSEMBLAGE VALUES
pdf("Results/Graphs/acrossTime-Means/jac-thre.pdf", 6.5, 5)
guidesLegend <- guides(fill=guide_legend(title="Mean Jaccard\nsimilarity index", reverse=T), colour=guide_legend(title="Mean Jaccard\nsimilarity index", reverse=T))
ggplot(rese.ModType.Mean[which(rese.ModType.Mean$variable %in% c("jacc.MEAN.thre")),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=1-value, colour=1-value)) +
    scale_fill_gradient2(low="darkorange", high="darkgreen", mid="white", midpoint=0.5, limits=c(0, 1)) +
    scale_colour_gradient2(low="black", mid="gray90", high="black", limits=c(0, 1)) +
    facet_grid(. ~ modelType) +
    themeOpt + guidesLegend + axesLabels + aspRatio + fitLine
dev.off()

pdf("Results/Graphs/acrossTime-Means/jac-prob.pdf", 6.5, 5)
guidesLegend <- guides(fill=guide_legend(title="Mean Jaccard\nsimilarity index", reverse=T), colour=guide_legend(title="Mean Jaccard\nsimilarity index", reverse=T))
ggplot(rese.ModType.Mean[which(rese.ModType.Mean$variable %in% c("jacc.MEAN.prob")),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=1-value, colour=1-value)) +
    scale_fill_gradient2(low="darkorange", high="darkgreen", mid="white", midpoint=0.5, limits=c(0, 1)) +
    scale_colour_gradient2(low="black", mid="gray90", high="black", limits=c(0, 1)) +
    facet_grid(. ~ modelType) +
    themeOpt + guidesLegend + axesLabels + aspRatio + fitLine
dev.off()

pdf("Results/Graphs/acrossTime-Means/jac-cal-thre.pdf", 6.5, 5)
guidesLegend <- guides(fill=guide_legend(title="Mean Jaccard\ncalibration", reverse=T), colour=guide_legend(title="Mean Jaccard\ncalibration", reverse=T))
ggplot(rese.ModType.Mean[which(rese.ModType.Mean$variable %in% c("jacc.CAL.thre")),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=value, colour=value)) +
    scale_fill_gradientn(colours=c("darkorange","white","darkgreen"), limits=c(-5.5, 1), values=rescale(c(-5.5, 0, 1))) +
    scale_colour_gradientn(colours=c("black","gray90","black"), limits=c(-5.5, 1), values=rescale(c(-5.5, 0, 1))) +
    facet_grid(. ~ modelType) +
    themeOpt + guidesLegend + axesLabels + aspRatio + fitLine
dev.off()

pdf("Results/Graphs/acrossTime-Means/jac-dis-thre.pdf", 6.5, 5)
guidesLegend <- guides(fill=guide_legend(title="Mean Jaccard\ndiscrimination", reverse=T), colour=guide_legend(title="Mean Jaccard\ndiscrimination", reverse=T))
ggplot(rese.ModType.Mean[which(rese.ModType.Mean$variable %in% c("jacc.DIS.thre")),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=value, colour=value)) +
    scale_fill_gradient2(low="darkorange", high="darkgreen", mid="white", midpoint=0.5, limits=c(0, 1)) +
    scale_colour_gradient2(low="black", mid="gray90", high="black", limits=c(0, 1)) +
    facet_grid(. ~ modelType) +
    themeOpt + guidesLegend + axesLabels + aspRatio + fitLine
dev.off()

pdf("Results/Graphs/acrossTime-Means/jac-cal-prob.pdf", 6.5, 5)
guidesLegend <- guides(fill=guide_legend(title="Mean Jaccard\ncalibration", reverse=T), colour=guide_legend(title="Mean Jaccard\ncalibration", reverse=T))
ggplot(rese.ModType.Mean[which(rese.ModType.Mean$variable %in% c("jacc.CAL.prob")),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=value, colour=value)) +
    scale_fill_gradientn(colours=c("darkorange","white","darkgreen"), limits=c(-5.5,1),values=rescale(c(-5.5,0,1))) +
    scale_colour_gradientn(colours=c("black","gray90","black"), limits=c(-5.5, 1), values=rescale(c(-5.5, 0, 1))) +
    facet_grid(. ~ modelType) +
    themeOpt + guidesLegend + axesLabels + aspRatio + fitLine
dev.off()

pdf("Results/Graphs/acrossTime-Means/jac-dis-prob.pdf", 6.5, 5)
guidesLegend <- guides(fill=guide_legend(title="Mean Jaccard\ndiscrimination", reverse=T), colour=guide_legend(title="Mean Jaccard\ndiscrimination", reverse=T))
ggplot(rese.ModType.Mean[which(rese.ModType.Mean$variable %in% c("jacc.DIS.prob")),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=value, colour=value)) +
    scale_fill_gradient2(low="darkorange", high="darkgreen", mid="white", midpoint=0.5, limits=c(0, 1)) +
    scale_colour_gradient2(low="black", mid="gray90", high="black", limits=c(0, 1)) +
    facet_grid(. ~ modelType) +
    themeOpt + guidesLegend + axesLabels + aspRatio + fitLine
dev.off()


pdf("Results/Graphs/acrossTime-Means/spr-dis.pdf", 6.5, 5)
guidesLegend <- guides(fill=guide_legend(title="Mean species richness\ndiscrimination", reverse=T), colour=guide_legend(title="Mean species richness\ndiscrimination", reverse=T))
ggplot(rese.ModType.Mean[which(rese.ModType.Mean$variable %in% c("spRich.COR.thre")),]) +
      geom_tile(aes(x=periodPrj, y=periodFit, fill=value, colour=value)) +
      scale_fill_gradient2(low="darkorange", high="darkgreen", mid="white", midpoint=0, limits=c(-1, 1)) +
      scale_colour_gradient2(low="black", mid="gray90", high="black", limits=c(-1, 1)) +
      facet_grid(. ~ modelType) +
      themeOpt + guidesLegend + axesLabels + aspRatio + fitLine
dev.off()

pdf("Results/Graphs/acrossTime-Means/spr-cal.pdf", 6.5, 5)
guidesLegend <- guides(fill=guide_legend(title="Mean species richness\ncalibration", reverse=T), colour=guide_legend(title="Mean species richness\ncalibration", reverse=T))
ggplot(rese.ModType.Mean[which(rese.ModType.Mean$variable %in% c("spRich.RSQ.thre")),]) +
      geom_tile(aes(x=periodPrj, y=periodFit, fill=value, colour=value)) +
      scale_fill_gradientn(colours=c("darkorange","white","darkgreen"), limits=c(-4,1),values=rescale(c(-4,0,1))) +
      scale_colour_gradientn(colours=c("black","gray90","black"), limits=c(-4, 1), values=rescale(c(-4, 0, 1))) +
      facet_grid(. ~ modelType) +
      themeOpt + guidesLegend + axesLabels + aspRatio + fitLine
dev.off()


##############################################################################
## Specific models results 
################################################################################

# AUC
pdf("Results/Graphs/acrossTime-Means/auc-models.pdf", 12.8, 8.5)
guidesLegend <- guides(fill=guide_legend(title="Mean AUC", reverse=T), colour=guide_legend(title="Mean AUC", reverse=T))
ggplot(resd.Models.Mean[which(resd.Models.Mean$variable == "auc"),]) +
      geom_tile(aes(x=periodPrj, y=periodFit, fill=value, colour=value)) +
      scale_fill_gradient2(low="darkorange", high="darkgreen", mid="white", midpoint=0.75, limits=c(0.485, 1)) + 
      scale_colour_gradient(low="gray90", high="black", limits=c(0.485, 1)) +
      facet_grid(modelType ~ sdm_eq) +
      themeOpt + guidesLegend + axesLabels + aspRatio + fitLine
dev.off()

# Sens.
pdf("Results/Graphs/acrossTime-Means/sens-models.pdf", 12.8, 8.5)
guidesLegend <- guides(fill=guide_legend(title="Mean Sens.", reverse=T), colour=guide_legend(title="Mean Sens.", reverse=T))
ggplot(resd.Models.Mean[which(resd.Models.Mean$variable == "sens"),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=value, colour=value)) +
    scale_fill_gradient2(low="darkorange", high="darkgreen", mid="white", midpoint=0.5, limits=c(0, 1)) + 
    scale_colour_gradient(low="gray90", high="black", limits=c(0, 1)) +
    facet_grid(modelType ~ sdm_eq) +
    themeOpt + guidesLegend + axesLabels + aspRatio + fitLine
dev.off()

# Spec.
pdf("Results/Graphs/acrossTime-Means/spec-models.pdf", 12.8, 8.5)
guidesLegend <- guides(fill=guide_legend(title="Mean Spec.", reverse=T), colour=guide_legend(title="Mean Spec.", reverse=T))
ggplot(resd.Models.Mean[which(resd.Models.Mean$variable == "spec"),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=value, colour=value)) +
    scale_fill_gradient2(low="darkorange", high="darkgreen", mid="white", midpoint=0.5, limits=c(0, 1)) + 
    scale_colour_gradient(low="gray90", high="black", limits=c(0, 1)) +
    facet_grid(modelType ~ sdm_eq) +
    themeOpt + guidesLegend + axesLabels + aspRatio + fitLine
dev.off()

# BRIER
pdf("Results/Graphs/acrossTime-Means/brier-models.pdf", 12.8, 8.5)
guidesLegend <- guides(fill=guide_legend(title="Mean 1-Brier", reverse=T), colour=guide_legend(title="Mean 1-Brier", reverse=T))
ggplot(resd.Models.Mean[which(resd.Models.Mean$variable == "brier"),]) +
      geom_tile(aes(x=periodPrj, y=periodFit, fill=1-value, colour=1-value)) +
      scale_fill_gradient2(low="darkorange", high="darkgreen", mid="white", midpoint=0.5, limits=c(0, 1)) +
      scale_colour_gradient(low="gray90", high="black", limits=c(0, 1)) +
      facet_grid(modelType ~ sdm_eq) +
      themeOpt + guidesLegend + axesLabels + aspRatio + fitLine
dev.off()


# MEAN JACCARD INDEX
pdf("Results/Graphs/acrossTime-Means/jac-thre-models.pdf", 12.8, 8.5)
guidesLegend <- guides(fill=guide_legend(title="Mean Jaccard\nsimilarity index", reverse=T), colour=guide_legend(title="Mean Jaccard\nsimilarity index", reverse=T))
ggplot(rese.Models.Mean[which(rese.Models.Mean$variable == "jacc.MEAN.thre"),]) +
      geom_tile(aes(x=periodPrj, y=periodFit, fill=1-value, colour=1-value)) +
      scale_fill_gradient2(low="darkorange", high="darkgreen", mid="white", midpoint=0.5, limits=c(0, 1)) +
      scale_colour_gradient2(low="black", mid="gray90", high="black", limits=c(0, 1)) +
      facet_grid(modelType ~ sdm_eq) +
      themeOpt + guidesLegend + axesLabels + aspRatio + fitLine
dev.off()

pdf("Results/Graphs/acrossTime-Means/jac-prob-models.pdf", 12.8, 8.5)
guidesLegend <- guides(fill=guide_legend(title="Mean Jaccard\nsimilarity index", reverse=T), colour=guide_legend(title="Mean Jaccard\nsimilarity index", reverse=T))
ggplot(rese.Models.Mean[which(rese.Models.Mean$variable == "jacc.MEAN.prob"),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=1-value, colour=1-value)) +
    scale_fill_gradient2(low="darkorange", high="darkgreen", mid="white", midpoint=0.5, limits=c(0, 1)) +
    scale_colour_gradient2(low="black", mid="gray90", high="black", limits=c(0, 1)) +
    facet_grid(modelType ~ sdm_eq) +
    themeOpt + guidesLegend + axesLabels + aspRatio + fitLine
dev.off()

# JACCARD CALIBRATION
pdf("Results/Graphs/acrossTime-Means/jac-cal-thre-models.pdf", 12.8, 8.5)
guidesLegend <- guides(fill=guide_legend(title="Mean Jaccard\ncalibration", reverse=T), colour=guide_legend(title="Mean Jaccard\ncalibration", reverse=T))
ggplot(rese.Models.Mean[which(rese.Models.Mean$variable %in% c("jacc.CAL.thre")),]) +
      geom_tile(aes(x=periodPrj, y=periodFit, fill=value, colour=value)) +
      scale_fill_gradientn(colours=c("darkorange","white","darkgreen"), limits=c(-21,1),values=rescale(c(-21,0,1))) +
      scale_colour_gradientn(colours=c("black","gray90","black"), limits=c(-21, 1), values=rescale(c(-21, 0, 1))) +
      facet_grid(modelType ~ sdm_eq) +
      themeOpt + guidesLegend + axesLabels + aspRatio + fitLine
dev.off()

pdf("Results/Graphs/acrossTime-Means/jac-dis-thre-models.pdf", 12.8, 8.5)
guidesLegend <- guides(fill=guide_legend(title="Mean Jaccard\ndiscrimination", reverse=T), colour=guide_legend(title="Mean Jaccard\ndiscrimination", reverse=T))
ggplot(rese.Models.Mean[which(rese.Models.Mean$variable %in% c("jacc.DIS.thre")),]) +
      geom_tile(aes(x=periodPrj, y=periodFit, fill=value, colour=value)) +
      scale_fill_gradient2(low="darkorange", high="darkgreen", mid="white", midpoint=0.5, limits=c(0, 1)) +
      scale_colour_gradient2(low="black", mid="gray90", high="black", limits=c(0, 1)) +
      facet_grid(modelType ~ sdm_eq) +
      themeOpt + guidesLegend + axesLabels + aspRatio + fitLine
dev.off()

# JACCARD CALIBRATION - PROB
pdf("Results/Graphs/acrossTime-Means/jac-cal-prob-models.pdf", 12.8, 8.5)
guidesLegend <- guides(fill=guide_legend(title="Mean Jaccard\ncalibration", reverse=T), colour=guide_legend(title="Mean Jaccard\ncalibration", reverse=T))
ggplot(rese.Models.Mean[which(rese.Models.Mean$variable %in% c("jacc.CAL.prob")),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=value, colour=value)) +
    scale_fill_gradientn(colours=c("darkorange","white","darkgreen"), limits=c(-21,1),values=rescale(c(-21,0,1))) +
    scale_colour_gradientn(colours=c("black","gray90","black"), limits=c(-21, 1), values=rescale(c(-21, 0, 1))) +
    facet_grid(modelType ~ sdm_eq) +
    themeOpt + guidesLegend + axesLabels + aspRatio + fitLine
dev.off()

pdf("Results/Graphs/acrossTime-Means/jac-dis-prob-models.pdf", 12.8, 8.5)
guidesLegend <- guides(fill=guide_legend(title="Mean Jaccard\ndiscrimination", reverse=T), colour=guide_legend(title="Mean Jaccard\ndiscrimination", reverse=T))
ggplot(rese.Models.Mean[which(rese.Models.Mean$variable %in% c("jacc.DIS.prob")),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=value, colour=value)) +
    scale_fill_gradient2(low="darkorange", high="darkgreen", mid="white", midpoint=0.5, limits=c(0, 1)) +
    scale_colour_gradient2(low="black", mid="gray90", high="black", limits=c(0, 1)) +
    facet_grid(modelType ~ sdm_eq) +
    themeOpt + guidesLegend + axesLabels + aspRatio + fitLine
dev.off()


pdf("Results/Graphs/acrossTime-Means/spr-dis-models.pdf", 12.8, 8.5)
guidesLegend <- guides(fill=guide_legend(title="Mean species richness\ndiscrimination", reverse=T), colour=guide_legend(title="Mean species richness\ndiscrimination", reverse=T))
ggplot(rese.Models.Mean[which(rese.Models.Mean$variable %in% c("spRich.COR.prob")),]) +
      geom_tile(aes(x=periodPrj, y=periodFit, fill=value, colour=value)) +
      scale_fill_gradient2(low="darkorange", high="darkgreen", mid="white", midpoint=0, limits=c(-1, 1)) +
      scale_colour_gradient2(low="black", mid="gray90", high="black", limits=c(-1, 1)) +
      facet_grid(modelType ~ sdm_eq) +
      themeOpt + guidesLegend + axesLabels + aspRatio + fitLine
dev.off()

pdf("Results/Graphs/acrossTime-Means/spr-cal-models.pdf", 12.8, 8.5)
guidesLegend <- guides(fill=guide_legend(title="Mean species richness\ncalibration", reverse=T), colour=guide_legend(title="Mean species richness\ncalibration", reverse=T))
ggplot(rese.Models.Mean[which(rese.Models.Mean$variable %in% c("spRich.RSQ.prob")),]) +
      geom_tile(aes(x=periodPrj, y=periodFit, fill=value, colour=value)) +
      scale_fill_gradientn(colours=c("darkorange","white","darkgreen"), limits=c(-6,1),values=rescale(c(-6,0,1))) +
      scale_colour_gradientn(colours=c("black","gray90","black"), limits=c(-6, 1), values=rescale(c(-6, 0, 1))) +
      facet_grid(modelType ~ sdm_eq) +
      themeOpt + guidesLegend + axesLabels + aspRatio + fitLine
dev.off()



################################################################################
##
## DIFFERENCES
##
################################################################################

guidesLegend <- guides(fill=guide_legend(title="Difference", reverse=T, order=1), size=guide_legend(title="Significance", reverse=T, order=2))

f_labeller <- function(var, value){
  value <- as.character(value)
  if (var=="variable") { 
    value <- "Diff."
  }
  return(value)
}

# DIFFERENCES IN MEAN DISTRIBUTION VALUES
pdf("Results/Graphs/acrossTime-Dif/auc-Dif.pdf", 4.5, 4.85)
ggplot(resd.ModType.Dif2[which(resd.ModType.Dif2$variable %in% c("auc")),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
    geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
    scale_fill_gradient2(limits=c(-0.08, 0.08), low="#B2182B", high="#2166AC") +
    scale_size_discrete(range=c(0, 1)) +
    facet_grid(.~ variable, labeller=f_labeller) +
    fitLine + themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("Results/Graphs/acrossTime-Dif/sens-Dif.pdf", 4.5, 4.85)
ggplot(resd.ModType.Dif2[which(resd.ModType.Dif2$variable %in% c("sens")),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
    geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
    scale_fill_gradient2(limits=c(-0.11, 0.11), low="#B2182B", high="#2166AC") +
    scale_size_discrete(range=c(0, 1)) +
    facet_grid(.~ variable, labeller=f_labeller) +
    fitLine + themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("Results/Graphs/acrossTime-Dif/spec-Dif.pdf", 4.5, 4.85)
ggplot(resd.ModType.Dif2[which(resd.ModType.Dif2$variable %in% c("spec")),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
    geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
    scale_fill_gradient2(limits=c(-0.12, 0.12), low="#B2182B", high="#2166AC") +
    scale_size_discrete(range=c(0, 1)) +
    facet_grid(.~ variable, labeller=f_labeller) +
    fitLine + themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("Results/Graphs/acrossTime-Dif/brier-Dif.pdf", 4.5, 4.85)
ggplot(resd.ModType.Dif2[which(resd.ModType.Dif2$variable %in% c("brier")),]) +
      geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
      geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
      scale_fill_gradient2(limits=c(-0.1, 0.1), low="#B2182B", high="#2166AC") +
      scale_size_discrete(range=c(0, 1)) +
      facet_grid(.~ variable, labeller=f_labeller) +
      fitLine + themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()


# DIFFERENCES IN MEAN ASSEMBLAGE VALUES
pdf("Results/Graphs/acrossTime-Dif/jac-thre-Dif.pdf", 4.5, 4.85)
ggplot(rese.ModType.Dif2[which(rese.ModType.Dif2$variable %in% c("jacc.MEAN.thre")),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
    geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
    scale_fill_gradient2(limits=c(-0.12, 0.12), low="#B2182B", high="#2166AC") +
    scale_size_discrete(range=c(0, 1)) +
    facet_grid(.~ variable, labeller=f_labeller) +
    fitLine + themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("Results/Graphs/acrossTime-Dif/jac-prob-Dif.pdf", 4.5, 4.85)
ggplot(rese.ModType.Dif2[which(rese.ModType.Dif2$variable %in% c("jacc.MEAN.prob")),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
    geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
    scale_fill_gradient2(limits=c(-0.1, 0.1), low="#B2182B", high="#2166AC") +
    scale_size_discrete(range=c(0, 1)) +
    facet_grid(.~ variable, labeller=f_labeller) +
    fitLine + themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("Results/Graphs/acrossTime-Dif/jac-cal-thre-Dif.pdf", 4.5, 4.85)
ggplot(rese.ModType.Dif2[which(rese.ModType.Dif2$variable %in% c("jacc.CAL.thre")),]) +
      geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
      geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
      scale_fill_gradient2(limits=c(-4, 4), low="#B2182B", high="#2166AC") +
      scale_size_discrete(range=c(0, 1)) +
      facet_grid(.~ variable, labeller=f_labeller) +
      fitLine + themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()
 
pdf("Results/Graphs/acrossTime-Dif/jac-dis-thre-Dif.pdf", 4.5, 4.85)
ggplot(rese.ModType.Dif2[which(rese.ModType.Dif2$variable %in% c("jacc.DIS.thre")),]) +
      geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
      geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
      scale_fill_gradient2(limits=c(-0.15, 0.15), low="#B2182B", high="#2166AC") +
      scale_size_discrete(range=c(0, 1)) +
      facet_grid(.~ variable, labeller=f_labeller) +
      fitLine + themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()
 
pdf("Results/Graphs/acrossTime-Dif/jac-cal-prob-Dif.pdf", 4.5, 4.85)
ggplot(rese.ModType.Dif2[which(rese.ModType.Dif2$variable %in% c("jacc.CAL.prob")),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
    geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
    scale_fill_gradient2(limits=c(-6, 6), low="#B2182B", high="#2166AC") +
    scale_size_discrete(range=c(0, 1)) +
    facet_grid(.~ variable, labeller=f_labeller) +
    fitLine + themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("Results/Graphs/acrossTime-Dif/jac-dis-prob-Dif.pdf", 4.5, 4.85)
ggplot(rese.ModType.Dif2[which(rese.ModType.Dif2$variable %in% c("jacc.DIS.prob")),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
    geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
    scale_fill_gradient2(limits=c(-0.15, 0.15), low="#B2182B", high="#2166AC") +
    scale_size_discrete(range=c(0, 1)) +
    facet_grid(.~ variable, labeller=f_labeller) +
    fitLine + themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("Results/Graphs/acrossTime-Dif/spRich-dis-Dif.pdf", 4.5, 4.85)
ggplot(rese.ModType.Dif2[which(rese.ModType.Dif2$variable %in% c("spRich.COR.thre")),]) +
      geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
      geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
      scale_fill_gradient2(limits=c(-0.3, 0.3), low="#B2182B", high="#2166AC") +
      scale_size_discrete(range=c(0, 1)) +
      facet_grid(.~ variable, labeller=f_labeller) +
      fitLine + themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("Results/Graphs/acrossTime-Dif/spRich-cal-Dif.pdf", 4.5, 4.85)
ggplot(rese.ModType.Dif2[which(rese.ModType.Dif2$variable %in% c("spRich.RSQ.thre")),]) +
      geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
      geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
      scale_fill_gradient2(limits=c(-2, 2), low="#B2182B", high="#2166AC") +
      scale_size_discrete(range=c(0, 1)) +
      facet_grid(.~ variable, labeller=f_labeller) +
      fitLine + themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()


# DIFFERENCES BY MODELS IN MEAN DISTRIBUTION VALUES
f_labeller <- function(var, value){
  value <- as.character(value)
  if (var=="variable") { 
    value <- paste("Diff:", value, sep="")
  }
  return(value)
}

pdf("Results/Graphs/acrossTime-Dif/auc-Dif-models.pdf", 12.8, 8.5)
ggplot(resd.Models.Dif2[which(resd.Models.Dif2$variable %in% c("auc")),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
    geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
    scale_fill_gradient2(limits=c(-0.26, 0.26), midpoint=0, low="#B2182B", high="#2166AC") +
    scale_size_discrete(range=c(0, 1)) +
    facet_grid(. ~ sdm_eq) +
    fitLine + themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("Results/Graphs/acrossTime-Dif/sens-Dif-models.pdf", 12.8, 8.5)
ggplot(resd.Models.Dif2[which(resd.Models.Dif2$variable %in% c("sens")),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
    geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
    scale_fill_gradient2(limits=c(-0.3, 0.3), midpoint=0, low="#B2182B", high="#2166AC") +
    scale_size_discrete(range=c(0, 1)) +
    facet_grid(. ~ sdm_eq) +
    fitLine + themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("Results/Graphs/acrossTime-Dif/spec-Dif-models.pdf", 12.8, 8.5)
ggplot(resd.Models.Dif2[which(resd.Models.Dif2$variable %in% c("spec")),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
    geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
    scale_fill_gradient2(limits=c(-0.31, 0.31), midpoint=0, low="#B2182B", high="#2166AC") +
    scale_size_discrete(range=c(0, 1)) +
    facet_grid(. ~ sdm_eq) +
    fitLine + themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("Results/Graphs/acrossTime-Dif/brier-Dif-models.pdf", 12.8, 8.5)
ggplot(resd.Models.Dif2[which(resd.Models.Dif2$variable %in% c("brier")),]) +
      geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
      geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
      scale_fill_gradient2(limits=c(-0.2, 0.2), midpoint=0, low="#B2182B", high="#2166AC") +
      scale_size_discrete(range=c(0, 1)) +
      facet_grid(. ~ sdm_eq) +
      fitLine + themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()


pdf("Results/Graphs/acrossTime-Dif/jac-thre-Dif-models.pdf", 12.8, 12.5)
ggplot(rese.Models.Dif2[which(rese.Models.Dif2$variable == "jacc.MEAN.thre"),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
    geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
    scale_fill_gradient2(limits=c(-0.25, 0.25), midpoint=0, low="#B2182B", high="#2166AC") +
    scale_size_discrete(range=c(0, 1)) + facet_grid(. ~ sdm_eq) +
    fitLine + themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("Results/Graphs/acrossTime-Dif/jac-prob-Dif-models.pdf", 12.8, 12.5)
ggplot(rese.Models.Dif2[which(rese.Models.Dif2$variable == "jacc.MEAN.prob"),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
    geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
    scale_fill_gradient2(limits=c(-0.2, 0.2), midpoint=0, low="#B2182B", high="#2166AC") +
    scale_size_discrete(range=c(0, 1)) + facet_grid(. ~ sdm_eq) +
    fitLine + themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("Results/Graphs/acrossTime-Dif/jac-cal-thre-Dif-models.pdf", 12.8, 12.5)
ggplot(rese.Models.Dif2[which(rese.Models.Dif2$variable == "jacc.CAL.thre"),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
    geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
    scale_fill_gradient2(limits=c(-8, 8), midpoint=0, low="#B2182B", high="#2166AC") +
    scale_size_discrete(range=c(0, 1)) + facet_grid(. ~ sdm_eq) +
    fitLine + themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("Results/Graphs/acrossTime-Dif/jac-dis-thre-Dif-models.pdf", 12.8, 12.5)
ggplot(rese.Models.Dif2[which(rese.Models.Dif2$variable == "jacc.DIS.thre"),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
    geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
    scale_fill_gradient2(limits=c(-0.3, 0.3), midpoint=0, low="#B2182B", high="#2166AC") +
    scale_size_discrete(range=c(0, 1)) + facet_grid(. ~ sdm_eq) +
    fitLine + themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("Results/Graphs/acrossTime-Dif/jac-cal-prob-Dif-models.pdf", 12.8, 12.5)
ggplot(rese.Models.Dif2[which(rese.Models.Dif2$variable == "jacc.CAL.prob"),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
    geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
    scale_fill_gradient2(limits=c(-20, 20), midpoint=0, low="#B2182B", high="#2166AC") +
    scale_size_discrete(range=c(0, 1)) + facet_grid(. ~ sdm_eq) +
    fitLine + themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("Results/Graphs/acrossTime-Dif/jac-dis-prob-Dif-models.pdf", 12.8, 12.5)
ggplot(rese.Models.Dif2[which(rese.Models.Dif2$variable == "jacc.DIS.prob"),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
    geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
    scale_fill_gradient2(limits=c(-0.3, 0.3), midpoint=0, low="#B2182B", high="#2166AC") +
    scale_size_discrete(range=c(0, 1)) + facet_grid(. ~ sdm_eq) +
    fitLine + themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("Results/Graphs/acrossTime-Dif/spRich-dis-Dif-models.pdf", 12.8, 12.5)
ggplot(rese.Models.Dif2[which(rese.Models.Dif2$variable == "spRich.COR.thre"),]) +
      geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
      geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
      scale_fill_gradient2(limits=c(-0.85, 0.85), midpoint=0, low="#B2182B", high="#2166AC") +
      scale_size_discrete(range=c(0, 1)) +
      facet_grid(. ~ sdm_eq) +
      fitLine + themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()

pdf("Results/Graphs/acrossTime-Dif/spRich-cal-Dif-models.pdf", 12.8, 12.5)
ggplot(rese.Models.Dif2[which(rese.Models.Dif2$variable == "spRich.RSQ.thre"),]) +
      geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
      geom_tile(aes(x=periodPrj, y=periodFit, size=significant2), colour="black", fill=NA) +
      scale_fill_gradient2(limits=c(-4, 4), midpoint=0, low="#B2182B", high="#2166AC") +
      scale_size_discrete(range=c(0, 1)) +
      facet_grid(. ~ sdm_eq) +
      fitLine + themeOpt + axesLabels + aspRatio + guidesLegend
dev.off()




####################################################################################
##
##  CLIMATIC DISTANCE PLOTS
##
####################################################################################

pdf("Results/Graphs/climaticDistance/auc-minDist.pdf", 4, 4)
ggplot(resd.ModType.climbioDist[which(resd.ModType.climbioDist$variable %in% c("auc")),], aes(x=Clim.Dist, y=value, colour=modelType)) +
      geom_point(shape=21) +
      geom_smooth(aes(linetype=modelType), method=loess, colour="black", size=1) +
      scale_colour_manual(values=c("gray20", "gray80")) +
      labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="AUC") +
      theme_bw() +
      guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type")) +
      theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/climaticDistance/brier-minDist.pdf", 4, 4)
ggplot(resd.ModType.climbioDist[which(resd.ModType.climbioDist$variable %in% c("brier")),], aes(x=Clim.Dist, y=1-value, colour=modelType)) +
      geom_point(shape=21) +
      geom_smooth(aes(linetype=modelType), method=loess, colour="black", size=1) +
      scale_colour_manual(values=c("gray20", "gray80")) +
      labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="1 - Brier score") +
      theme_bw() +
      guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type")) +
      theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/climaticDistance/jac-thre-minDist.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.MEAN.thre")),], aes(x=Clim.Dist, y=1-value, colour=modelType)) +
    geom_point(shape=21) +
    geom_smooth(aes(linetype=modelType), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="Jaccard similarity index\n(obs. vs pred.)") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type")) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/climaticDistance/jac-prob-minDist.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.MEAN.prob")),], aes(x=Clim.Dist, y=1-value, colour=modelType)) +
    geom_point(shape=21) +
    geom_smooth(aes(linetype=modelType), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="Jaccard similarity index\n(obs. vs pred.)") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type")) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/climaticDistance/jac-dis-thre-minDist.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.DIS.thre")),], aes(x=Clim.Dist, y=value, colour=modelType)) +
    geom_point(shape=21) +
    geom_smooth(aes(linetype=modelType), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="Jaccard discrimination") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type")) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/climaticDistance/jac-cal-thre-minDist.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.CAL.thre")),], aes(x=Clim.Dist, y=value, colour=modelType)) +
    geom_point(shape=21) +
    geom_smooth(aes(linetype=modelType), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="Jaccard calibration") +
    coord_cartesian(ylim=c(-6,1)) +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type")) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/climaticDistance/jac-dis-prob-minDist.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.DIS.prob")),], aes(x=Clim.Dist, y=value, colour=modelType)) +
      geom_point(shape=21) +
      geom_smooth(aes(linetype=modelType), method=loess, colour="black", size=1) +
      scale_colour_manual(values=c("gray20", "gray80")) +
      labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="Jaccard discrimination") +
      theme_bw() +
      guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type")) +
      theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/climaticDistance/jac-cal-prob-minDist.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.CAL.prob")),], aes(x=Clim.Dist, y=value, colour=modelType)) +
      geom_point(shape=21) +
      geom_smooth(aes(linetype=modelType), method=loess, colour="black", size=1) +
      scale_colour_manual(values=c("gray20", "gray80")) +
      labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="Jaccard calibration") +
      scale_y_continuous(limits=c(-6,1)) +
      theme_bw() +
      guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type")) +
      theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/climaticDistance/spr-dis-minDist.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("spRich.COR.prob")),], aes(x=Clim.Dist, y=value, colour=modelType)) +
      geom_point(shape=21) +
      geom_smooth(aes(linetype=modelType), method=loess, colour="black", size=1) +
      scale_colour_manual(values=c("gray20", "gray80")) +
      labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="Species richness discrimination") +
      theme_bw() +
      guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type")) +
      theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/climaticDistance/spr-cal-minDist.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("spRich.RSQ.prob")),], aes(x=Clim.Dist, y=value, colour=modelType)) +
      geom_point(shape=21) +
      geom_smooth(aes(linetype=modelType), method=loess, colour="black", size=1) +
      scale_colour_manual(values=c("gray20", "gray80")) +
      labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="Species richness calibration") +
      theme_bw() +
      guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type")) +
      theme(legend.justification=c(0,0), legend.position=c(0,0))
dev.off()

####################################################################################
#
# Climatic Distance Plots with Weights
#
####################################################################################
pdf("Results/Graphs/climaticDistance/auc-minDist-weight.pdf", 4, 4)
ggplot(resd.ModType.climbioDist[which(resd.ModType.climbioDist$variable %in% c("auc")),], aes(x=Clim.Dist, y=value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="AUC") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/climaticDistance/brier-minDist-weight.pdf", 4, 4)
ggplot(resd.ModType.climbioDist[which(resd.ModType.climbioDist$variable %in% c("brier")),], aes(x=Clim.Dist, y=1-value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="1 - Brier score") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/climaticDistance/jac-thre-minDist-weight.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.MEAN.thre")),], aes(x=Clim.Dist, y=1-value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="Jaccard similarity index\n(obs. vs pred.)") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/climaticDistance/jac-prob-minDist-weight.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.MEAN.prob")),], aes(x=Clim.Dist, y=1-value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="Jaccard similarity index\n(obs. vs pred.)") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/climaticDistance/jac-dis-thre-minDist-weight.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.DIS.thre")),], aes(x=Clim.Dist, y=value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="Jaccard discrimination") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/climaticDistance/jac-cal-thre-minDist-weight.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.CAL.thre")),], aes(x=Clim.Dist, y=value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="Jaccard calibration") +
    coord_cartesian(ylim=c(-6,1)) +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/climaticDistance/jac-dis-prob-minDist-weight.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.DIS.prob")),], aes(x=Clim.Dist, y=value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="Jaccard discrimination") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/climaticDistance/jac-cal-prob-minDist-weight.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.CAL.prob")),], aes(x=Clim.Dist, y=value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="Jaccard calibration") +
    scale_y_continuous(limits=c(-6,1)) +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/climaticDistance/spr-dis-minDist-weight.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("spRich.COR.prob")),], aes(x=Clim.Dist, y=value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="Species richness discrimination") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/climaticDistance/spr-cal-minDist-weight.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("spRich.RSQ.prob")),], aes(x=Clim.Dist, y=value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="Species richness calibration") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.justification=c(0,0), legend.position=c(0,0))
dev.off()




####################################################################################
#
# Climatic Distance Plots with Weights without legend
#
####################################################################################
pdf("Results/Graphs/climaticDistance/auc-minDist-weight-wo_legend.pdf", 4, 4)
ggplot(resd.ModType.climbioDist[which(resd.ModType.climbioDist$variable %in% c("auc")),], aes(x=Clim.Dist, y=value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="AUC") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.position="none")
dev.off()

pdf("Results/Graphs/climaticDistance/brier-minDist-weight-wo_legend.pdf", 4, 4)
ggplot(resd.ModType.climbioDist[which(resd.ModType.climbioDist$variable %in% c("brier")),], aes(x=Clim.Dist, y=1-value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="1 - Brier score") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.position="none")
dev.off()

pdf("Results/Graphs/climaticDistance/jac-thre-minDist-weight-wo_legend.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.MEAN.thre")),], aes(x=Clim.Dist, y=1-value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="Jaccard similarity index\n(obs. vs pred.)") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.position="none")
dev.off()

pdf("Results/Graphs/climaticDistance/jac-prob-minDist-weight-wo_legend.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.MEAN.prob")),], aes(x=Clim.Dist, y=1-value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="Jaccard similarity index\n(obs. vs pred.)") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.position="none")
dev.off()

pdf("Results/Graphs/climaticDistance/jac-dis-thre-minDist-weight-wo_legend.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.DIS.thre")),], aes(x=Clim.Dist, y=value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="Jaccard discrimination") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.position="none")
dev.off()

pdf("Results/Graphs/climaticDistance/jac-cal-thre-minDist-weight-wo_legend.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.CAL.thre")),], aes(x=Clim.Dist, y=value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="Jaccard calibration") +
    coord_cartesian(ylim=c(-6,1)) +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.position="none")
dev.off()

pdf("Results/Graphs/climaticDistance/jac-dis-prob-minDist-weight-wo_legend.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.DIS.prob")),], aes(x=Clim.Dist, y=value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="Jaccard discrimination") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.position="none")
dev.off()

pdf("Results/Graphs/climaticDistance/jac-cal-prob-minDist-weight-wo_legend.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.CAL.prob")),], aes(x=Clim.Dist, y=value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="Jaccard calibration") +
    scale_y_continuous(limits=c(-6,1)) +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.position="none")
dev.off()

pdf("Results/Graphs/climaticDistance/spr-dis-minDist-weight-wo_legend.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("spRich.COR.prob")),], aes(x=Clim.Dist, y=value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="Species richness discrimination") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.position="none")
dev.off()

pdf("Results/Graphs/climaticDistance/spr-cal-minDist-weight-wo_legend.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("spRich.RSQ.prob")),], aes(x=Clim.Dist, y=value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Climatic distance \n(Euclidean distance in PCA space)", y="Species richness calibration") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.position="none")
dev.off()




####################################################################################
##
## BIOTIC DISTANCE PLOTS
##
####################################################################################

pdf("Results/Graphs/bioticDistance/auc-minDist.pdf", 4, 4)
ggplot(resd.ModType.climbioDist[which(resd.ModType.climbioDist$variable %in% c("auc")),], aes(x=Bio.Dist, y=value, colour=modelType)) +
      geom_point(shape=21) +
      geom_smooth(aes(linetype=modelType), method=loess, colour="black", size=1) +
      scale_colour_manual(values=c("gray20", "gray80")) +
      labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="AUC") +
      theme_bw() +
      guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type")) +
      theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/bioticDistance/brier-minDist.pdf", 4, 4)
ggplot(resd.ModType.climbioDist[which(resd.ModType.climbioDist$variable %in% c("brier")),], aes(x=Bio.Dist, y=1-value, colour=modelType)) +
      geom_point(shape=21) +
      geom_smooth(aes(linetype=modelType), method=loess, colour="black", size=1) +
      scale_colour_manual(values=c("gray20", "gray80")) +
      labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="1 - Brier score") +
      theme_bw() +
      guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type")) +
      theme(legend.justification=c(0,0), legend.position=c(0,0))
dev.off()

pdf("Results/Graphs/bioticDistance/jac-thre-minDist.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.MEAN.thre")),], aes(x=Bio.Dist, y=1-value, colour=modelType)) +
      geom_point(shape=21) +
      geom_smooth(aes(linetype=modelType), method=loess, colour="black", size=1) +
      scale_colour_manual(values=c("gray20", "gray80")) +
      labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="Jaccard similarity index\n(obs. vs pred.)") +
      theme_bw() +
      guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type")) +
      theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/bioticDistance/jac-prob-minDist.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.MEAN.prob")),], aes(x=Bio.Dist, y=1-value, colour=modelType)) +
    geom_point(shape=21) +
    geom_smooth(aes(linetype=modelType), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="Jaccard similarity index\n(obs. vs pred.)") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type")) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/bioticDistance/jac-dis-thre-minDist.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.DIS.thre")),], aes(x=Bio.Dist, y=value, colour=modelType)) +
    geom_point(shape=21) +
    geom_smooth(aes(linetype=modelType), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="Jaccard discrimination") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type")) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/bioticDistance/jac-cal-thre-minDist.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.CAL.thre")),], aes(x=Bio.Dist, y=value, colour=modelType)) +
    geom_point(shape=21) +
    geom_smooth(aes(linetype=modelType), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    coord_cartesian(ylim=c(-6,1)) +
    labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="Jaccard calibration") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type")) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/bioticDistance/jac-dis-prob-minDist.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.DIS.prob")),], aes(x=Bio.Dist, y=value, colour=modelType)) +
      geom_point(shape=21) +
      geom_smooth(aes(linetype=modelType), method=loess, colour="black", size=1) +
      scale_colour_manual(values=c("gray20", "gray80")) +
      labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="Jaccard discrimination") +
      theme_bw() +
      guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type")) +
      theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/bioticDistance/jac-cal-prob-minDist.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.CAL.prob")),], aes(x=Bio.Dist, y=value, colour=modelType)) +
      geom_point(shape=21) +
      geom_smooth(aes(linetype=modelType), method=loess, colour="black", size=1) +
      scale_colour_manual(values=c("gray20", "gray80")) +
      scale_y_continuous(limits=c(-6,1)) +
      labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="Jaccard calibration") +
      theme_bw() +
      guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type")) +
      theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()
             
pdf("Results/Graphs/bioticDistance/spr-dis-minDist.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("spRich.COR.prob")),], aes(x=Bio.Dist, y=value, colour=modelType)) +
      geom_point(shape=21) +
      geom_smooth(aes(linetype=modelType), method=loess, colour="black", size=1) +
      scale_colour_manual(values=c("gray20", "gray80")) +
      labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="Species richness discrimination") +
      theme_bw() +
      guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type")) +
      theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/bioticDistance/spr-cal-minDist.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("spRich.RSQ.prob")),], aes(x=Bio.Dist, y=value, colour=modelType)) +
      geom_point(shape=21) +
      geom_smooth(aes(linetype=modelType), method=loess, colour="black", size=1) +
      scale_colour_manual(values=c("gray20", "gray80")) +
      labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="Species richness calibration") +
      theme_bw() +
      guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type")) +
      theme(legend.justification=c(0,0), legend.position=c(0,0))
dev.off()




####################################################################################
#
# Biotic Distance Plots with weights (sample size)
#
####################################################################################

pdf("Results/Graphs/bioticDistance/auc-minDist-weight.pdf", 4, 4)
ggplot(resd.ModType.climbioDist[which(resd.ModType.climbioDist$variable %in% c("auc")),], aes(x=Bio.Dist, y=value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="AUC") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/bioticDistance/brier-minDist-weight.pdf", 4, 4)
ggplot(resd.ModType.climbioDist[which(resd.ModType.climbioDist$variable %in% c("brier")),], aes(x=Bio.Dist, y=1-value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="1 - Brier score") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.justification=c(0,0), legend.position=c(0,0))
dev.off()

pdf("Results/Graphs/bioticDistance/jac-thre-minDist-weight.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.MEAN.thre")),], aes(x=Bio.Dist, y=1-value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="Jaccard similarity index\n(obs. vs pred.)") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/bioticDistance/jac-prob-minDist-weight.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.MEAN.prob")),], aes(x=Bio.Dist, y=1-value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="Jaccard similarity index\n(obs. vs pred.)") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/bioticDistance/jac-dis-thre-minDist-weight.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.DIS.thre")),], aes(x=Bio.Dist, y=value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="Jaccard discrimination") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/bioticDistance/jac-cal-thre-minDist-weight.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.CAL.thre")),], aes(x=Bio.Dist, y=value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    scale_y_continuous(limits=c(-6,1)) +
    labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="Jaccard calibration") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/bioticDistance/jac-dis-prob-minDist-weight.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.DIS.prob")),], aes(x=Bio.Dist, y=value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="Jaccard discrimination") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/bioticDistance/jac-cal-prob-minDist-weight.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.CAL.prob")),], aes(x=Bio.Dist, y=value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    scale_y_continuous(limits=c(-6,1)) +
    labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="Jaccard calibration") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/bioticDistance/spr-dis-minDist-weight.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("spRich.COR.prob")),], aes(x=Bio.Dist, y=value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="Species richness discrimination") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
dev.off()

pdf("Results/Graphs/bioticDistance/spr-cal-minDist-weight.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("spRich.RSQ.prob")),], aes(x=Bio.Dist, y=value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="Species richness calibration") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.justification=c(0,0), legend.position=c(0,0))
dev.off()



####################################################################################
#
# Biotic Distance Plots with weights (sample size) without legend
#
####################################################################################

pdf("Results/Graphs/bioticDistance/auc-minDist-weight-wo_legend.pdf", 4, 4)
ggplot(resd.ModType.climbioDist[which(resd.ModType.climbioDist$variable %in% c("auc")),], aes(x=Bio.Dist, y=value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="AUC") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.position="none")
dev.off()

pdf("Results/Graphs/bioticDistance/brier-minDist-weight-wo_legend.pdf", 4, 4)
ggplot(resd.ModType.climbioDist[which(resd.ModType.climbioDist$variable %in% c("brier")),], aes(x=Bio.Dist, y=1-value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="1 - Brier score") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.position="none")
dev.off()

pdf("Results/Graphs/bioticDistance/jac-thre-minDist-weight-wo_legend.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.MEAN.thre")),], aes(x=Bio.Dist, y=1-value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="Jaccard similarity index\n(obs. vs pred.)") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.position="none")
dev.off()

pdf("Results/Graphs/bioticDistance/jac-prob-minDist-weight-wo_legend.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.MEAN.prob")),], aes(x=Bio.Dist, y=1-value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="Jaccard similarity index\n(obs. vs pred.)") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.position="none")
dev.off()

pdf("Results/Graphs/bioticDistance/jac-dis-thre-minDist-weight-wo_legend.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.DIS.thre")),], aes(x=Bio.Dist, y=value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="Jaccard discrimination") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.position="none")
dev.off()

pdf("Results/Graphs/bioticDistance/jac-cal-thre-minDist-weight-wo_legend.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.CAL.thre")),], aes(x=Bio.Dist, y=value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    coord_cartesian(ylim=c(-6,1)) +
    labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="Jaccard calibration") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.position="none")
dev.off()

pdf("Results/Graphs/bioticDistance/jac-dis-prob-minDist-weight-wo_legend.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.DIS.prob")),], aes(x=Bio.Dist, y=value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="Jaccard discrimination") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.position="none")
dev.off()

pdf("Results/Graphs/bioticDistance/jac-cal-prob-minDist-weight-wo_legend.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("jacc.CAL.prob")),], aes(x=Bio.Dist, y=value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    scale_y_continuous(limits=c(-6,1)) +
    labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="Jaccard calibration") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.position="none")
dev.off()

pdf("Results/Graphs/bioticDistance/spr-dis-minDist-weight-wo_legend.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("spRich.COR.prob")),], aes(x=Bio.Dist, y=value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="Species richness discrimination") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.position="none")
dev.off()

pdf("Results/Graphs/bioticDistance/spr-cal-minDist-weight-wo_legend.pdf", 4, 4)
ggplot(rese.ModType.climbioDist[which(rese.ModType.climbioDist$variable %in% c("spRich.RSQ.prob")),], aes(x=Bio.Dist, y=value, colour=modelType)) +
    geom_point(shape=21, aes(size=sampleSize)) +
    geom_smooth(aes(linetype=modelType, weight=1/sampleSize), method=loess, colour="black", size=1) +
    scale_colour_manual(values=c("gray20", "gray80")) +
    labs(x="Ecological distance\n(Jaccard dissimilarity index)", y="Species richness calibration") +
    theme_bw() +
    guides(linetype=guide_legend(title="Model type"), colour=guide_legend(title="Model type"), size=guide_legend(title="Sample size")) +
    theme(legend.position="none")
dev.off()




################################################################################
##
## PREVALENCE PARTITIONING
##
################################################################################

forecastLine <- data.frame(x=c(0,3,5,5,6,6,7,7,8,8,9,9,10,10,11,12)+0.5, y=c(0,2,3,5,6,10,10,12,12,13,14,16,17,19,19,20)+0.5)

# These plotting functions can be modified to reflect the fixed scales (min and max values) and use more attractive colors.

fitLine <- geom_line(data=forecastLine, aes(x=x, y=y), colour="black", size=1.3, alpha=0.75)

themeOpt <- theme(axis.text.x=element_text(angle=90, hjust=1)) + theme_bw()

axesLabels <- labs(x="Projecting period", y="Fitting period")

aspRatio <- coord_fixed(ratio=1)

#################################################################################
## CLM versus SDM comparison
#################################################################################

guidesLegend <- guides(fill=guide_legend(title="Difference", reverse=T))

f_labeller <- function(var, value){
    value <- as.character(value)
    if (var=="variable") { 
        value <- "Diff."
    }
    return(value)
}

# DIFFERENCES IN MEAN DISTRIBUTION VALUES
pdf("Results/Graphs/Prevalence/auc-lowPrevalence.pdf", 4.5, 4.85)
ggplot(resd.ModType.Dif.lowPrev[which(resd.ModType.Dif.lowPrev$variable %in% c("auc")),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
    scale_fill_gradient2(limits=c(-0.1, 0.1), low="#B2182B", high="#2166AC") +
    facet_grid(.~ variable, labeller=f_labeller) +
    fitLine + themeOpt + axesLabels + aspRatio + guidesLegend + ggtitle("low prevalence (< 0.5)")
dev.off()
pdf("Results/Graphs/Prevalence/auc-highPrevalence.pdf", 4.5, 4.85)
ggplot(resd.ModType.Dif.highPrev[which(resd.ModType.Dif.highPrev$variable %in% c("auc")),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
    scale_fill_gradient2(limits=c(-0.1, 0.1), low="#B2182B", high="#2166AC") +
    facet_grid(.~ variable, labeller=f_labeller) +
    fitLine + themeOpt + axesLabels + aspRatio + guidesLegend + ggtitle("high prevalence (> 0.5)")
dev.off()

# DIFFERENCES BY MODELS IN MEAN DISTRIBUTION VALUES
f_labeller <- function(var, value){
    value <- as.character(value)
    if (var=="variable") { 
        value <- paste("Diff:", value, sep="")
    }
    return(value)
}
pdf("Results/Graphs/Prevalence/auc-models-lowPrevalence.pdf", 12.8, 8.5)
ggplot(resd.Models.Dif.lowPrev[which(resd.Models.Dif.lowPrev$variable %in% c("auc")),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
    scale_fill_gradient2(limits=c(-0.26, 0.26), midpoint=0, low="#B2182B", high="#2166AC") +
    facet_grid(. ~ sdm_eq) +
    fitLine + themeOpt + axesLabels + aspRatio + guidesLegend + ggtitle("low prevalence (< 0.5)")
dev.off()
pdf("Results/Graphs/Prevalence/auc-models-highPrevalence.pdf", 12.8, 8.5)
ggplot(resd.Models.Dif.highPrev[which(resd.Models.Dif.highPrev$variable %in% c("auc")),]) +
    geom_tile(aes(x=periodPrj, y=periodFit, fill=dif), colour="gray90") +
    scale_fill_gradient2(limits=c(-0.26, 0.26), midpoint=0, low="#B2182B", high="#2166AC") +
    facet_grid(. ~ sdm_eq) +
    fitLine + themeOpt + axesLabels + aspRatio + guidesLegend + ggtitle("high prevalence (> 0.5)")
dev.off()

