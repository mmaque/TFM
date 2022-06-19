#!/usr/bin/R

### Packages ###################################################################
library(edgeR)
library(BuildConsensus)
library(MKmisc)
library(pROC)

### Functions ##################################################################
cluster_classify <- function(data, centroid, method = "pearson") {
  R <- stats::cor(data, centroid, method = method)
  scores <- apply(R, 1, which.max)
}

### Base Directory #############################################################
basedir <-  '/home/miguel/Escritorio/Documents/Panchain/'

### Data Loading ###############################################################
TCGA_Data <- read.table(paste0(basedir,'Input/TCGA_Counts_PAAD.tsv'), 
                        header = TRUE)
Biolinks <-  readRDS(paste0(basedir,'Input/PAAD_TCGA_Biolinks.rds'))
Symbols_DB <- read.table(paste0(basedir,'Symbols/ENSG2_SYMBOL.tsv'), 
                         header = TRUE)

### Data Annotation and Deduplication ##########################################
## Annotation
Annotated_Data <- merge(Symbols_DB,TCGA_Data,by = 0,1)
## Deduplication
Annotated_Data <- Annotated_Data[which(!duplicated
                                       (Annotated_Data$Gene_symbol)=='TRUE'),]
## Parsing
rownames(Annotated_Data) <- Annotated_Data$Gene_symbol
Annotated_Data[,c(1,2)] <- NULL

### Data Normalization #########################################################
Norm_Data <- as.data.frame(cpm
                           (calcNormFactors
                             (DGEList(counts=as.matrix(Annotated_Data)))))

### Loading Classifications and Parsing from Biolinks ##########################
## Sample Names
identical(colnames(TCGA_Data),gsub('-','.',colnames(Biolinks)))
Sample_Names <- colnames(TCGA_Data)
## Moffit et al. (2015)
Moffit <- data.frame(row.names = Sample_Names,Biolinks$`paper_mRNA Moffitt clusters (All 150 Samples) 1basal  2classical`)
colnames(Moffit) <- 'Moffit'
Moffit[Moffit == 1] <- 'Basal'
Moffit[Moffit == 2] <- 'Classical'

##Collison et al. (2011)
Collison <- data.frame(row.names = Sample_Names,Biolinks$`paper_mRNA Collisson clusters (All 150 Samples) 1classical 2exocrine 3QM`)
colnames(Collison) <- 'Collison'
Collison[Collison == 1] <- 'Classical'
Collison[Collison == 2] <- 'Exocrine'
Collison[Collison == 3] <- 'Quasimesenchymal'
##Bailey et al. (2016)
Bailey <- data.frame(row.names = Sample_Names,Biolinks$`paper_mRNA Bailey Clusters (All 150 Samples) 1squamous 2immunogenic 3progenitor 4ADEX`)
colnames(Bailey) <- 'Bailey'
Bailey[Bailey == 1] <- 'Squamous'
Bailey[Bailey == 2] <- 'Immunogenic'
Bailey[Bailey == 3] <- 'Progenitor'
Bailey[Bailey == 4] <- 'ADEX'

### Subclassification from Puleo et al. (2019) #################################
Puleo_Centroids <- read.table(paste0
                              (basedir,'Centroids/Puleo_Centroids_HC.tsv'),
                              header = TRUE, sep = '\t', row.names = 'GeneSymbol')


Puleo_TCGA_Sorted <-  Norm_Data[rownames(Puleo_Centroids),]
                              
Puleo <- cluster_classify(Puleo_TCGA_Sorted,
                          Puleo_Centroids,
                          method = 'spearman')

Puleo <- as.data.frame(Puleo)
Puleo[Puleo == 1] <- 'Basal'
Puleo[Puleo == 2] <- 'Classical'

### Subclassifications Dataframe ###############################################
## Creation

Classification <- data.frame(row.names = Sample_Names,
                             Collison = Collison$Collison,
                             Moffit = Moffit$Moffit,
                             Bailey = Bailey$Bailey,
                             Puleo = Puleo$Puleo)

# Classification_All_Cases[is.na(Classification_All_Cases)] <- 'Not_Available'
# Discordant_Biolinks <- Biolinks[,gsub('-','.',colnames(Biolinks)) %in% rownames(Classification_All_Cases[Classification_All_Cases$Moffit!=Classification_All_Cases$Puleo,])]
# Concordant_Biolinks <- Biolinks[,!(gsub('-','.',colnames(Biolinks)) %in% rownames(Classification_All_Cases[Classification_All_Cases$Moffit!=Classification_All_Cases$Puleo,]))]

## Removal of Samples with NA
Classification <- Classification[complete.cases(Classification),]

### MCL Consensus Network ######################################################
MCL_Res <- consensus.MCL(Classification,
                         I.values = 3,
                         outdir = "HC_Consensus_Solution",
                         n.iter = 500,
                         sim.method = "CohenKappa")

### Hypergeometric Testing for Selection of Consensus Samples ##################
## Count Number of Successes for each Consensus Group
Consensus_Classical_Counts <- as.data.frame(rowSums(Classification == 'Classical') + 
                                              rowSums(Classification == 'Progenitor'))

## Count Immunogenic as Basal (Comment if not)
# Consensus_Basal_Counts <- as.data.frame(rowSums(Classification == 'Basal') + 
#                                           rowSums(Classification == 'Quasimesenchymal') +
#                                           rowSums(Classification == 'Squamous') +
#                                           rowSums(Classification == 'Immunogenic'))

Consensus_Basal_Counts <- as.data.frame(rowSums(Classification == 'Basal') + 
                                          rowSums(Classification == 'Quasimesenchymal') +
                                          rowSums(Classification == 'Squamous'))

Consensus_Exo_Counts <- as.data.frame(rowSums(Classification == 'Exocrine') + 
                                        rowSums(Classification == 'ADEX'))
## Hypergeometric Test for each Consensus Group
Consensus_Classical_Samples_Names <- rownames(Consensus_Classical_Counts)[phyper(as.vector(Consensus_Classical_Counts[,1]),4,4,4,lower.tail = FALSE) <= 0 ]
Consensus_Basal_Samples_Names <- rownames(Consensus_Basal_Counts)[phyper(as.vector(Consensus_Basal_Counts[,1]),4,4,4,lower.tail = FALSE) <= 0 ]
Consensus_Exo_Samples_Names <- rownames(Consensus_Exo_Counts)[phyper(as.vector(Consensus_Exo_Counts[,1]),4,4,4,lower.tail = FALSE) <= 0.25 ]

### Core Sample RNA-Seq Data ###################################################
Core_Classical <- Norm_Data[colnames(Norm_Data) %in% Consensus_Classical_Samples_Names]
Core_Basal <- Norm_Data[colnames(Norm_Data) %in% Consensus_Basal_Samples_Names]
Core_Exo <- Norm_Data[colnames(Norm_Data) %in% Consensus_Exo_Samples_Names]

Core_Samples <- cbind(Core_Classical,Core_Basal,Core_Exo)
Core_Samples <- as.matrix(Core_Samples)

### T Testing ##################################################################
Factor_Classical <- factor(c(rep(1,ncol(Core_Classical)),rep(0,ncol(Core_Basal) + ncol(Core_Exo))))
Factor_Basal <- factor(c(rep(0,ncol(Core_Classical)),rep(1,ncol(Core_Basal)),rep(0,ncol(Core_Exo))))
Factor_Exo <- factor(c(rep(0,ncol(Core_Classical) + ncol(Core_Basal)),rep(1,ncol(Core_Exo))))

T_Test_Classical <- mod.t.test(Core_Samples,group = Factor_Classical)
T_Test_Basal <- mod.t.test(Core_Samples,group = Factor_Basal)
T_Test_Exo <- mod.t.test(Core_Samples,group = Factor_Exo)

Core_Samples <- as.data.frame(t(Core_Samples))

PVal_Classical_Filtered_Genes <- rownames(T_Test_Classical[T_Test_Classical$p.value <= 0.001,])
PVal_Classical_Filtered_Core_Samples <- Core_Samples[colnames(Core_Samples) %in% PVal_Classical_Filtered_Genes]

PVal_Basal_Filtered_Genes <- rownames(T_Test_Basal[T_Test_Basal$p.value < 0.001,])
PVal_Basal_Filtered_Core_Samples <- Core_Samples[colnames(Core_Samples) %in% PVal_Basal_Filtered_Genes]

PVal_Exo_Filtered_Genes <- rownames(T_Test_Exo[T_Test_Exo$p.value < 0.001,])
PVal_Exo_Filtered_Core_Samples <- Core_Samples[colnames(Core_Samples) %in% PVal_Exo_Filtered_Genes]

### AUC threshold ##############################################################
AUC_Threshold <- 0.9

###Classical AUC Selection######################################################
Classical_Fitted_Values <- list()
for (i in c(1:ncol(PVal_Classical_Filtered_Core_Samples))) {
  vector <- fitted.values(glm(Factor_Classical~PVal_Classical_Filtered_Core_Samples[,i],family = binomial))
  Classical_Fitted_Values[[i]] <- vector
}

Classical_Aucs <- vector()
for (i in c(1:length(Classical_Fitted_Values))) {
  result <- auc(roc(Factor_Classical,Classical_Fitted_Values[[i]]))
  Classical_Aucs[[i]] <- result
}

AUC_Classical_Genes <- data.frame(colnames(PVal_Classical_Filtered_Core_Samples),Classical_Aucs)
rownames(AUC_Classical_Genes) <- AUC_Classical_Genes$colnames.PVal_Classical_Filtered_Core_Samples.
Final_Classical_Genes <- rownames(AUC_Classical_Genes[AUC_Classical_Genes$Classical_Aucs >= AUC_Threshold,])

###Basal AUC Selection######################################################
Basal_Fitted_Values <- list()
for (i in c(1:ncol(PVal_Basal_Filtered_Core_Samples))) {
  vector <- fitted.values(glm(Factor_Basal~PVal_Basal_Filtered_Core_Samples[,i],family = binomial))
  Basal_Fitted_Values[[i]] <- vector
}

Basal_Aucs <- vector()
for (i in c(1:length(Basal_Fitted_Values))) {
  result <- auc(roc(Factor_Basal,Basal_Fitted_Values[[i]]))
  Basal_Aucs[[i]] <- result
}

AUC_Basal_Genes <- data.frame(colnames(PVal_Basal_Filtered_Core_Samples),Basal_Aucs)
rownames(AUC_Basal_Genes) <- AUC_Basal_Genes$colnames.PVal_Basal_Filtered_Core_Samples.
Final_Basal_Genes <- rownames(AUC_Basal_Genes[AUC_Basal_Genes$Basal_Aucs >= AUC_Threshold,])

###Exo AUC Selection########################################################
Exo_Fitted_Values <- list()
for (i in c(1:ncol(PVal_Exo_Filtered_Core_Samples))) {
  vector <- fitted.values(glm(Factor_Exo~PVal_Exo_Filtered_Core_Samples[,i],family = binomial))
  Exo_Fitted_Values[[i]] <- vector
}

Exo_Aucs <- vector()
for (i in c(1:length(Exo_Fitted_Values))) {
  result <- auc(roc(Factor_Exo,Exo_Fitted_Values[[i]]))
  Exo_Aucs[[i]] <- result
}

AUC_Exo_Genes <- data.frame(colnames(PVal_Exo_Filtered_Core_Samples),Exo_Aucs)
rownames(AUC_Exo_Genes) <- AUC_Exo_Genes$colnames.PVal_Exo_Filtered_Core_Samples.
Final_Exo_Genes <- rownames(AUC_Exo_Genes[AUC_Exo_Genes$Exo_Aucs >= AUC_Threshold,])

###Final total genes############################################################
Final_Total_Genes <- c(Final_Classical_Genes,Final_Basal_Genes,Final_Exo_Genes)
Final_Total_Genes <- unique(Final_Total_Genes)
Core_Samples_Filtered <- as.data.frame(t(Core_Samples[,colnames(Core_Samples) %in% Final_Total_Genes]))

###Obtention of logfc for classical samples#####################################
Classical_DGE <- DGEList(counts=Core_Samples_Filtered, group=Factor_Classical)
Classical_DGE <- estimateDisp(Classical_DGE)
Classical_LogFC <- exactTest(Classical_DGE,dispersion = 'common')

###Obtention of logfc for basal samples#########################################
Basal_DGE <- DGEList(counts=Core_Samples_Filtered, group=Factor_Basal)
Basal_DGE <- estimateDisp(Basal_DGE)
Basal_LogFC <- exactTest(Basal_DGE,dispersion = 'common')

###Obtention of logfc for basal samples#########################################
Exo_DGE <- DGEList(counts=Core_Samples_Filtered, group=Factor_Exo)
Exo_DGE <- estimateDisp(Exo_DGE)
Exo_LogFC <- exactTest(Exo_DGE,dispersion = 'common')

###Creation of centroids########################################################
# identical(rownames(Exo_LogFC$table),rownames(Basal_LogFC$table))
Centroids <- data.frame('Classical' = Classical_LogFC$table$logFC,
                        'Basal' = Basal_LogFC$table$logFC,
                        'Exo' = Exo_LogFC$table$logFC,
                        'hgnc_symbol' = rownames(Exo_LogFC$table),
                        row.names = rownames(Exo_LogFC$table))



# Best_Genes <- read.table('Best_Genes.txt',header =TRUE)
# Filtered_Centroid <- Centroids[rownames(Centroids) %in% Best_Genes$Best_Genes,]
# write.table(Centroids,'Centroids_4C_3S_0.9AUC.tsv',sep = '\t')
# write.table(Filtered_Centroid,'Filtered_Centroids_No_Imm_HC_0.9_3.tsv',sep = '\t')
