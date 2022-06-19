#!/usr/bin/R

### Packages -------------------------------------------------------------------
library(edgeR)
library(ggplot2)
library(singscore)
library(pheatmap)
library(wesanderson)
library(survival)
library(ggpubr)
library(cowplot)
library(webr)
library(dplyr)

### Functions ------------------------------------------------------------------
getConsensusClass <- function(x, minCor = .2, gene_id = c("entrezgene", "ensembl_gene_id", "hgnc_symbol")[1]){
  
  centroids <- Centroids
  lev.cs <- c("Classical","Basal","Exo")
  
  if(is.vector(x)) {
    if(is.null(names(x))) stop("Input vector of gene expression is missing names.\n The names must be the type of gene identifiers specified by the gene_id argument.")
    x <- data.frame(ss = x, row.names = names(x))
  }
  
  gkeep <- intersect(centroids[, gene_id], rownames(x))
  if (length(gkeep) == 0) stop("Empty intersection between profiled genes and the genes used for consensus classification.\n Make sure that gene names correspond to the type of identifiers specified by the gene_id argument")
  if (length(gkeep) < 0.5 * nrow(centroids)) warning("Input gene expression profile(s) include(s) less than half of the genes used for consensus classification. Results may not be relevant") 
  cor.dat <- as.data.frame(cor(x[gkeep, ], centroids[match(gkeep, centroids[, gene_id]), lev.cs], use = "complete.obs"), row.names = colnames(x))
  
  # Best correlated centroid
  cor.dat$nearestCentroid <- apply(cor.dat, 1, function(y){lev.cs[which.max(y)]})
  cor.dat$corToNearest <- apply(cor.dat[, lev.cs], 1, max)
  cor.dat$cor_pval <- sapply(colnames(x), function(smp){
    cor.test(x[gkeep, smp], centroids[match(gkeep, centroids[, gene_id]), cor.dat[smp, "nearestCentroid"]])$p.value
  })
  
  # Separation level metrics
  cor.dat$deltaSecondNearest <- apply(cor.dat$corToNearest - cor.dat[, lev.cs], 1, function(x){sort(x)[2]})
  cor.dat$deltaMed <- apply(cor.dat$corToNearest - cor.dat[, lev.cs], 1, median)
  cor.dat$separationLevel <- cor.dat$deltaSecondNearest/cor.dat$deltaMed
  
  cor.dat$consensusClass <- cor.dat$nearestCentroid
  
  # Set to NA if best correlation < minCor
  try(cor.dat[which(cor.dat$corToNearest < minCor), "consensusClass"] <-  NA)
  try(cor.dat[which(cor.dat$corToNearest < minCor), "separationLevel"] <-  NA)
  
  cor.dat <- cor.dat[, c("consensusClass" , "cor_pval", "separationLevel", lev.cs)]
  return(cor.dat)
}

cluster_classify <- function(data, centroid, method = "pearson") {
  R <- stats::cor(data, centroid, method = method)
  scores <- apply(R, 1, which.max)
}

### Data Loading ---------------------------------------------------------------
## Base directory
basedir <-  'C:/Users/maqju/Desktop/Projects/Panchain/'
## Metadata
Metadata <- read.table(paste0(basedir,'Input/PGyMB_DataForAnalysis_23032022.csv'),
                       header = TRUE, sep = ',', row.names = 'sampleid')
## Centroids
Centroids <- read.table(paste0(basedir,'Centroids/Centroids_4C_3S_0.9AUC.tsv'),sep = '\t')
# Centroids <- Centroids[,c(1,2,ncol(Centroids))]
## Data and Symbols
HTSeq_Data <- read.table(paste0(basedir,'Input/HTSeq66_ENS.tsv'), header = TRUE)
TCGA_Data <- read.table(paste0(basedir,'Input/TCGA_Counts_PAAD.tsv'), header = TRUE)
Biolinks <-  readRDS(paste0(basedir,'Input/PAAD_TCGA_Biolinks.rds'))
# HTSeq_Data <- read.table(paste0(basedir,'Input/TCGA_Counts_PAAD.tsv'), header = TRUE)
Symbols_DB <- read.table(paste0(basedir,'Symbols/ENSG2_SYMBOL.tsv'), header = TRUE)

### Genesets loading -----------------------------------------------------------
## Collison Signatures 
Collisson_Classical <- read.table(paste0(basedir,'Genesets/Collison_Classical.txt'), row.names = 1)
Collisson_Exocrine <- read.table(paste0(basedir,'Genesets/Collison_Exocrine.txt'), row.names = 1)
Collisson_Quasimesenchymal <- read.table(paste0(basedir,'Genesets/Collison_Quasimesenchymal.txt'), row.names = 1)
## Bailey Signatures
Bailey_Progenitor_UP <- read.table(paste0(basedir,'Genesets/Bailey_Progenitor_Up.txt'), row.names = 1)
Bailey_Progenitor_DOWN <- read.table(paste0(basedir,'Genesets/Bailey_Progenitor_Down.txt'), row.names = 1)
Bailey_Squamous_UP <- read.table(paste0(basedir,'Genesets/Bailey_Squamous_Up.txt'), row.names = 1)
Bailey_Squamous_DOWN <- read.table(paste0(basedir,'Genesets/Bailey_Squamous_Down.txt'), row.names = 1)
Bailey_ADEX_UP <- read.table(paste0(basedir,'Genesets/Bailey_ADEX_Up.txt'), row.names = 1)
Bailey_ADEX_DOWN <- read.table(paste0(basedir,'Genesets/Bailey_ADEX_Down.txt'), row.names = 1)
Bailey_Immunogenic_UP <- read.table(paste0(basedir,'Genesets/Bailey_Immunogenic_Up.txt'), row.names = 1)
Bailey_Immunogenic_DOWN <- read.table(paste0(basedir,'Genesets/Bailey_Immunogenic_Down.txt'), row.names = 1)
## Moffitt Signatures
Moffitt_Basal <- read.table(paste0(basedir,'Genesets/Validation_Factors/Basal_Factor.txt'), row.names = 1)
Moffitt_Classical <- read.table(paste0(basedir,'Genesets/Validation_Factors/Classical_Factor.txt'), row.names = 1)
Moffitt_Exocrine <- read.table(paste0(basedir,'Genesets/Validation_Factors/Exocrine_Factor.txt'), row.names = 1)
Moffitt_Immunogenic <- read.table(paste0(basedir,'Genesets/Validation_Factors/Immune_Factor.txt'), row.names = 1)
Moffitt_Stroma_A <- read.table(paste0(basedir,'Genesets/Validation_Factors/Stroma_A_Factor.txt'), row.names = 1)

### Data Annotation and Deduplication for PanGenEU -----------------------------
## Annotation for PanGenEU
Annotated_Data <- merge(Symbols_DB,HTSeq_Data,by = 1,1)
## Deduplication for PanGenEU
Annotated_Data <- Annotated_Data[which(!duplicated
                                       (Annotated_Data$Gene_symbol)=='TRUE'),]
## Parsing for PanGenEU
rownames(Annotated_Data) <- Annotated_Data$Gene_symbol
Annotated_Data[,c(1,2)] <- NULL

### Data Annotation and Deduplication for TCGA ---------------------------------
## Annotation for TCGA Data
Annotated_Data_TCGA <- merge(Symbols_DB,TCGA_Data,by = 0,1)
## Deduplication for TCGA Data
Annotated_Data_TCGA <- Annotated_Data_TCGA[which(!duplicated
                                       (Annotated_Data_TCGA$Gene_symbol)=='TRUE'),]
## Parsing for TCGA DAta
rownames(Annotated_Data_TCGA) <- Annotated_Data_TCGA$Gene_symbol
Annotated_Data_TCGA[,c(1,2)] <- NULL

### Data Normalization and Transformation for PanGenEU -------------------------
## TMM Normalization
Norm_Data <- as.data.frame(cpm
                           (calcNormFactors
                             (DGEList(counts=as.matrix(Annotated_Data)))))
## Log Tranformation of Data
Log_Data <- log2(Norm_Data + 1)

### Data Normalization and Transformation for TCGA -----------------------------
## TMM Normalization
Norm_Data_TCGA <- as.data.frame(cpm
                           (calcNormFactors
                             (DGEList(counts=as.matrix(Annotated_Data_TCGA)))))
## Log Tranformation of Data
Log_Data_TCGA <- log2(Norm_Data_TCGA + 1)


### Singscore Only PanGenEu ----------------------------------------------------
## Ranking of Data 
Rank_Data <- rankGenes(Norm_Data)
## Moffitt Scores
Moffitt_Basal_Score <- simpleScore(Rank_Data, upSet = Moffitt_Basal)
Moffitt_Classical_Score <- simpleScore(Rank_Data, upSet = Moffitt_Classical)
Moffitt_Immunogenic_Score <- simpleScore(Rank_Data, upSet = Moffitt_Immunogenic)
Moffitt_Exocrine_Score <- simpleScore(Rank_Data,upSet = Moffitt_Exocrine)
Moffitt_Stroma_A_Score <- simpleScore(Rank_Data,upSet = Moffitt_Stroma_A)
## Collissson Scores
Collisson_Classical_Score <- simpleScore(Rank_Data, upSet = Collisson_Classical)
Collisson_Quasimesenchymal_Score <- simpleScore(Rank_Data, upSet = Collisson_Quasimesenchymal)
Collisson_Exocrine_Score <- simpleScore(Rank_Data, upSet = Collisson_Exocrine)
## Bailey
Bailey_Progenitor_Score <- simpleScore(Rank_Data, upSet = Bailey_Progenitor_UP, downSet = Bailey_Progenitor_DOWN)
Bailey_Squamous_Score <- simpleScore(Rank_Data, upSet = Bailey_Squamous_UP, downSet = Bailey_Squamous_DOWN)
Bailey_ADEX_Score <- simpleScore(Rank_Data, upSet = Bailey_ADEX_UP, downSet = Bailey_ADEX_DOWN)
Bailey_Immunogenic_Score <- simpleScore(Rank_Data, upSet = Bailey_Immunogenic_UP, downSet = Bailey_Immunogenic_DOWN)
## Create a Table of all Scores
Scores <- data.frame(row.names = rownames(Moffitt_Basal_Score),
                     Collisson_Classical_Score = Collisson_Classical_Score$TotalScore,
                     Collisson_Quasimesenchymal_Score = Collisson_Quasimesenchymal_Score$TotalScore,
                     Collisson_Exocrine_Score = Collisson_Exocrine_Score$TotalScore,
                     Moffitt_Basal_Score = Moffitt_Basal_Score$TotalScore,
                     Moffitt_Classical_Score = Moffitt_Classical_Score$TotalScore,
                     Moffitt_Exocrine_Score = Moffitt_Exocrine_Score$TotalScore,
                     Moffitt_Immunogenic_Score = Moffitt_Immunogenic_Score$TotalScore,
                     Moffitt_Stroma_A_Score = Moffitt_Stroma_A_Score$TotalScore,
                     Bailey_Progenitor_Score = Bailey_Progenitor_Score$TotalScore,
                     Bailey_Squamous_Score = Bailey_Squamous_Score$TotalScore,
                     Bailey_ADEX_Score = Bailey_ADEX_Score$TotalScore,
                     Bailey_Immunogenic_Score = Bailey_Immunogenic_Score$TotalScore)


### Prediction (Both) ----------------------------------------------------------
Consensus_Classes <- getConsensusClass(Log_Data,gene_id = "hgnc_symbol",minCor = 0.2)
Consensus_Classes_TCGA <- getConsensusClass(Log_Data_TCGA,gene_id = "hgnc_symbol",minCor = 0.2)

### All Classification -------------------------------------------------------
Puleo_Centroids <- read.table(paste0
                              (basedir,'Centroids/Puleo_Centroids_HC.tsv'),
                              header = TRUE, sep = '\t', row.names = 'GeneSymbol')
Puleo_TCGA_Sorted <-  Norm_Data_TCGA[rownames(Puleo_Centroids),]
Puleo <- cluster_classify(Puleo_TCGA_Sorted,
                          Puleo_Centroids,
                          method = 'spearman')
Puleo <- as.data.frame(Puleo)
Puleo[Puleo == 1] <- 'Basal'
Puleo[Puleo == 2] <- 'Classical'
# Moffit et al. (2015)
Moffit <- Biolinks$`paper_mRNA Moffitt clusters (All 150 Samples) 1basal  2classical`
Moffit[Moffit == 1] <- 'Basal'
Moffit[Moffit == 2] <- 'Classical'
#Collison et al. (2011)
Collison <- Biolinks$`paper_mRNA Collisson clusters (All 150 Samples) 1classical 2exocrine 3QM`
Collison[Collison == 1] <- 'Classical'
Collison[Collison == 2] <- 'Exocrine'
Collison[Collison == 3] <- 'Quasimesenchymal'
#Bailey et al. (2016)
Bailey <- Biolinks$`paper_mRNA Bailey Clusters (All 150 Samples) 1squamous 2immunogenic 3progenitor 4ADEX`
Bailey[Bailey == 1] <- 'Squamous'
Bailey[Bailey == 2] <- 'Immunogenic'
Bailey[Bailey == 3] <- 'Progenitor'
Bailey[Bailey == 4] <- 'ADEX'

### Annotated Pheatmap ---------------------------------------------------------
# classes_annotation <- data.frame(row.names = colnames(TCGA_Data), Consensus = Consensus_Classes_TCGA$consensusClass, Puleo = Puleo, Bailey = Bailey, Collison = Collison, Moffit = Moffit)
# png('PHeatmap_All_182.png',width = 200, height = 400, units = 'mm', res = 300)
# pheatmap(Consensus_Classes_TCGA[,c(4,5,6)], show_rownames = F, annotation_row = classes_annotation, cutree_rows = 4)
# dev.off()

### Classification -------------------------------------------------------------
Classification <- data.frame(row.names = colnames(TCGA_Data), Consensus = Consensus_Classes_TCGA$consensusClass, Puleo = Puleo, Bailey = Bailey, Collison = Collison, Moffit = Moffit)
## Puleo
# PD <- Classification %>% group_by(Consensus, Puleo) %>% summarise(n = n())
# PD <- PD[complete.cases(PD),]
# png('Donut_Puleo.png',width = 100, height = 100, units = 'mm', res = 300)
# PieDonut(PD, aes(Consensus, Puleo, count=n), maxx=1.7)
# dev.off()
## Bailey
# PD <- Classification %>% group_by(Consensus, Bailey) %>% summarise(n = n())
# PD <- PD[complete.cases(PD),]
# png('Donut_Bailey.png',width = 100, height = 100, units = 'mm', res = 300)
# PieDonut(PD, aes(Consensus, Bailey, count=n),labelpositionThreshold = 0, maxx=1.7)
# dev.off()
## Collisson
# PD <- Classification %>% group_by(Consensus, Collison) %>% summarise(n = n())
# PD <- PD[complete.cases(PD),]
# png('Donut_Collisson.png',width = 100, height = 100, units = 'mm', res = 300)
# PieDonut(PD, aes(Consensus, Collison, count=n),labelpositionThreshold = 0, maxx=1.7)
# dev.off()
## Moffitt
# PD <- Classification %>% group_by(Consensus, Moffit) %>% summarise(n = n())
# PD <- PD[complete.cases(PD),]
# png('Donut_Moffitt.png',width = 100, height = 100, units = 'mm', res = 300)
# PieDonut(PD, aes(Consensus, Moffit, count=n),labelpositionThreshold = 0, maxx=1.7)
# dev.off()


### Metadata Preparation -------------------------------------------------------
## Merging and Parsing
Full_Data <- cbind(Consensus_Classes,Scores)
Full_Data <- Full_Data[complete.cases(Full_Data),]
Full_Data <- merge(Full_Data,Metadata,by=0,0)
rownames(Full_Data) <- Full_Data[,1]
Full_Data[,1] <- NULL
## Clustering by Tissue
Tissue <- Full_Data$consensusClass
Tissue[Tissue == 'Classical'] <- 'Tumor'
Tissue[Tissue == 'Basal'] <- 'Tumor'
Tissue[Tissue == 'Exo'] <- 'Pancreatic'
Full_Data <- cbind(Full_Data,Tissue)
## Setting a comparision list
Comparisons <- list(c('Basal','Classical'), c('Classical','Exo'), c('Basal', 'Exo'))
Comparisons_Classical <- list(c('Basal','Classical'), c('Classical','Exo'))
Comparisons_Basal <- list(c('Basal','Classical'), c('Basal','Exo'))
Comparisons_Exocrine <- list(c('Classical','Exo'), c('Basal', 'Exo'))

### Correlation plots ----------------------------------------------------------
## Classical and Basal
ggplot(Full_Data,aes(x=Full_Data$Classical,y=Full_Data$Basal)) +
  geom_point(aes(col = Full_Data$consensusClass)) +
  labs(col = 'Consensus Class') +
  scale_color_manual(values = wes_palette('Darjeeling1')) +
  ylab('Basal Correlation') +
  xlab('Classical Correlation')
## Basal and Exocrine
exo_basal <- ggplot(Full_Data,aes(x=Full_Data$Exo,y=Full_Data$Basal)) +
  geom_point(aes(col = Full_Data$consensusClass)) +
  labs(col = 'Consensus Class') +
  scale_color_manual(values = 'jco') +
  ylab('Basal Correlation') +
  xlab('Exocrine Correlation')
## Classical and Exocrine
exo_classical <- ggplot(Full_Data,aes(x=Full_Data$Exo,y=Full_Data$Classical)) +
  geom_point(aes(col = Full_Data$consensusClass)) +
  labs(col = 'Consensus Class') +
  scale_color_manual(values = 'jco') +
  ylab('Classical Correlation') +
  xlab('Exocrine Correlation')

### Singscore Plots ------------------------------------------------------------
## Moffitt Classical
Moffitt_Classical_Singscore_Plot <- ggboxplot(Full_Data, x = 'consensusClass', y = 'Moffitt_Classical_Score', color = 'consensusClass',
                                             pallette = 'jco', add = 'jitter') + 
  stat_compare_means(comparisons = Comparisons_Classical, label.y = c(0.35, 0.3675)) + 
  labs(color = 'Consensus Class') +
  ylab('Moffit Classical') +
  xlab('Consensus Class') +
  theme(axis.title.x = element_blank()) +
  theme(legend.position="none")
## Moffitt Basal
Moffitt_Basal_Singscore_Plot <- ggboxplot(Full_Data, x = 'consensusClass', y = 'Moffitt_Basal_Score', color = 'consensusClass',
          pallette = 'jco', add = 'jitter') +
  stat_compare_means(comparisons = Comparisons_Basal, label.y = c(0.355, 0.38)) + 
  labs(color = 'Consensus Class') +
  ylab('Moffit Basal') +
  xlab('Consensus Class') +
  theme(axis.title.x = element_blank()) +
  theme(legend.position="none")
## Moffitt Exocrine
Moffitt_Exocrine_Singscore_Plot <- ggboxplot(Full_Data, x = 'consensusClass', y = 'Moffitt_Exocrine_Score', color = 'consensusClass',
          pallette = 'jco', add = 'jitter') +
  stat_compare_means(comparisons = Comparisons_Exocrine, label.y = c(0.32, 0.35)) + 
  labs(color = 'Consensus Class') +
  ylab('Moffit Exocrine') +
  xlab('Consensus Class') +
  theme(legend.position="none")
## Collison Exocrine
Collisson_Exocrine_Singscore_Plot <- ggboxplot(Full_Data, x = 'consensusClass', y = 'Collisson_Exocrine_Score', color = 'consensusClass',
          pallette = 'jco', add = 'jitter') +
  stat_compare_means(comparisons = Comparisons_Exocrine, label.y = c(0.5, 0.56)) + 
  labs(color = 'Consensus Class') +
  ylab('Collisson Exocrine') +
  xlab('Consensus Class') +
  theme(legend.position="none")
##Collisson Classical
Collisson_Classical_Singscore_Plot <- ggboxplot(Full_Data, x = 'consensusClass', y = 'Collisson_Classical_Score', color = 'consensusClass',
          pallette = 'jco', add = 'jitter') +
  stat_compare_means(comparisons = Comparisons_Classical, label.y = c(0.5, 0.527)) + 
  labs(color = 'Consensus Class') +
  ylab('Collisson Classical') +
  xlab('Consensus Class') +
  theme(axis.title.x = element_blank()) +
  theme(legend.position="none")
## Collisson Quasimensenchymal
Collisson_Quasimesenchymal_Singscore_Plot <- ggboxplot(Full_Data, x = 'consensusClass', y = 'Collisson_Quasimesenchymal_Score', color = 'consensusClass',
          pallette = 'jco', add = 'jitter') +
  stat_compare_means(comparisons = Comparisons_Basal, label.y = c(0.45, 0.48)) + 
  labs(color = 'Consensus Class') +
  ylab('Collisson QM-PDA') +
  xlab('Consensus Class') +
  theme(axis.title.x = element_blank()) +
  theme(legend.position="none")
## Bailey Progenitor
Bailey_Progenitor_Singscore_Plot <- ggboxplot(Full_Data, x = 'consensusClass', y = 'Bailey_Progenitor_Score', color = 'consensusClass',
          pallette = 'jco', add = 'jitter') +
  stat_compare_means(comparisons = Comparisons_Classical, label.y = c(0.1, 0.14)) + 
  labs(color = 'Consensus Class') +
  ylab('Bailey Progenitor') +
  xlab('Consensus Class') +
  theme(axis.title.x = element_blank()) +
  theme(legend.position="none")
## Bailey Squamous
Bailey_Squamous_Singscore_Plot <- ggboxplot(Full_Data, x = 'consensusClass', y = 'Bailey_Squamous_Score', color = 'consensusClass',
          pallette = 'jco', add = 'jitter') +
  stat_compare_means(comparisons = Comparisons_Basal, label.y = c(0.2, 0.232)) + 
  labs(color = 'Consensus Class') +
  ylab('Bailey Squamous') +
  xlab('Consensus Class') +
  theme(axis.title.x = element_blank()) +
  theme(legend.position="none")
## Bailey ADEX
Bailey_ADEX_Singscore_Plot <- ggboxplot(Full_Data, x = 'consensusClass', y = 'Bailey_ADEX_Score', color = 'consensusClass',
          pallette = 'jco', add = 'jitter') +
  stat_compare_means(comparisons = Comparisons_Exocrine, label.y = c(-0.07, -0.025)) + 
  labs(color = 'Consensus Class') +
  ylab('Bailey ADEX') +
  xlab('Consensus Class') +
  theme(legend.position="none")
## Grid Plot
# png("rplot.png",width = 1000,height = 1000)
# plot_grid(Collisson_Classical_Singscore_Plot,Moffitt_Classical_Singscore_Plot,
#           Bailey_Progenitor_Singscore_Plot,Collisson_Quasimesenchymal_Singscore_Plot,
#           Moffitt_Basal_Singscore_Plot,Bailey_Squamous_Singscore_Plot,
#           Collisson_Exocrine_Singscore_Plot,Moffitt_Exocrine_Singscore_Plot,
#           Bailey_ADEX_Singscore_Plot)
# dev.off()

### Other Plots ----------------------------------------------------------------
## Tumor size
ggboxplot(Full_Data, x = 'consensusClass', y = 'tumour_size_2', color = 'consensusClass',
          pallette = 'jco', add = 'jitter') +
  stat_compare_means(comparisons = Comparisons, label.y = c(7, 8, 9)) + 
  stat_compare_means(label.y = 10) +
  labs(color = 'Consensus Class') +
  ylab('Tumour Size') +
  xlab('Consensus Class')
  
# ggplot(Full_Data, aes(x=as.factor(Full_Data$consensusClass),y=Full_Data$tumour_size_2)) +
#   geom_boxplot(aes(col = Full_Data$consensusClass), outlier.shape = NA) +
#   geom_jitter(aes(col = Full_Data$consensusClass)) +
#   ylab('Tumor size') +
#   xlab('Consensus Class') +
#   labs(col = 'Consensus Class') +
#   scale_color_manual(values = wes_palette('Darjeeling1')) +
#   stat_compare_means(label = 'p.format', label.y = 7.25) +
#   stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.",label.y = 6.9)

## Blocktype Plot
# ggplot(Full_Data,aes(x=Full_Data$Basal,y=Full_Data$Classical)) +
#   geom_point(aes(col = as.factor(Full_Data$blocktype))) +
#   labs(col = 'Blocktype')

### Heatmaps -------------------------------------------------------------------
# annotation_row <- data.frame(row.names = rownames(Full_Data), Consensus = Full_Data$consensusClass)
# pheatmap(Log_Data[rownames(Centroids),],show_rownames = T,show_colnames = F,cluster_rows = T)
# pheatmap(Full_Data[,c(4,5,6)],show_rownames = FALSE, annotation_row = annotation_row)

### Kaplan Curve ---------------------------------------------------------------
## Cox test
coxSurv <- coxph(Surv(Full_Data$triskMonths, Full_Data$vitalstatusnuevo) ~ Full_Data$consensusClass, data=Full_Data)
pvalCoxSurv = paste0("Log-rank p value < ", round(summary(coxSurv)$logtest[[3]], digits=4))
textToPlot = paste(pvalCoxSurv)
## Surv fitting
objectSurvFit <- survfit(Surv(Full_Data$triskMonths, Full_Data$vitalstatusnuevo) ~ Full_Data$consensusClass, data=Full_Data, conf.type="log-log")
## Plotting
Kaplan_PanGenEU <- plot(objectSurvFit,col=c("blue","red",'orange'), lwd=2, cex=1.5, yscale=1, fun="log", xlab="Months", xscale=365.25, ylab="Proportion surviving", mark.time=FALSE, xaxs="S", main="Survival time by consensus classes", cex.main=1)
legend("topright", c("Basal", "Classical", "Exocrine"), text.col=c("blue","red","orange"), fill=c("blue","red","orange"), border=c("blue","red","orange"), bg="transparent", box.col="transparent")
legend("top", textToPlot, text.col="forestgreen", bg="transparent", box.col="transparent")

### Kaplan Curve 2 -------------------------------------------------------------
# TCGA_Data <- read.table('TCGA_Data_Final.tsv')
## Cox test
# coxSurv <- coxph(Surv(TCGA_Data$Days_Death, TCGA_Data$vital) ~ TCGA_Data$Consensus, data=TCGA_Data)
# pvalCoxSurv = paste0("Log-rank p value < ", round(summary(coxSurv)$logtest[[3]], digits=4))
# textToPlot = paste(pvalCoxSurv)
## Surv fitting
# objectSurvFit_TCGA <- survfit(Surv(TCGA_Data$Days_Death, TCGA_Data$vital) ~ TCGA_Data$Consensus, data=TCGA_Data, conf.type="log-log")
## Plotting
# Kaplan_TCGA <- plot(objectSurvFit_TCGA,col=c("blue","red",'orange'), lwd=2, cex=1.5, yscale=1, fun="log", xlab="Months", xscale=365.25, ylab="Proportion surviving", mark.time=FALSE, xaxs="S", main="Survival time by consensus classes", cex.main=1)
# legend("topright", c("Basal", "Classical", "Exocrine"), text.col=c("blue","red","orange"), fill=c("blue","red","orange"), border=c("blue","red","orange"), bg="transparent", box.col="transparent")
# legend("top", textToPlot, text.col="forestgreen", bg="transparent", box.col="transparent")

