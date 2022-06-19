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
## Loading
TCGA_Data <- read.table(paste0(basedir,'Input/TCGA_Counts_PAAD.tsv'), 
                        header = TRUE)
Biolinks <-  readRDS(paste0(basedir,'Input/PAAD_TCGA_Biolinks.rds'))
Symbols_DB <- read.table(paste0(basedir,'Symbols/ENSG2_SYMBOL.tsv'), 
                         header = TRUE)
## Checking
identical(gsub('-','.',colnames(Biolinks)), colnames(TCGA_Data))

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
Sample_Names <- colnames(TCGA_Data)
## Moffit et al. (2015)
Moffit <- Biolinks$`paper_mRNA Moffitt clusters (All 150 Samples) 1basal  2classical`
Moffit[Moffit == 1] <- 'Basal'
Moffit[Moffit == 2] <- 'Classical'
##Collison et al. (2011)
Collison <- Biolinks$`paper_mRNA Collisson clusters (All 150 Samples) 1classical 2exocrine 3QM`
Collison[Collison == 1] <- 'Classical'
Collison[Collison == 2] <- 'Exocrine'
Collison[Collison == 3] <- 'Quasimesenchymal'
##Bailey et al. (2016)
Bailey <- Biolinks$`paper_mRNA Bailey Clusters (All 150 Samples) 1squamous 2immunogenic 3progenitor 4ADEX`
Bailey[Bailey == 1] <- 'Squamous'
Bailey[Bailey == 2] <- 'Immunogenic'
Bailey[Bailey == 3] <- 'Progenitor'
Bailey[Bailey == 4] <- 'ADEX'

### Subclassification from Puleo et al. (2019) #################################
Puleo_Centroids <- read.table(paste0
                              (basedir,'Centroids/Puleo_Centroids.tsv'),
                              header = TRUE, sep = '\t', row.names = 'GeneSymbol')

Puleo_TCGA_Sorted <-  Norm_Data[rownames(Puleo_Centroids),]

Puleo <- cluster_classify(Puleo_TCGA_Sorted,
                          Puleo_Centroids,
                          method = 'spearman')

Puleo[Puleo == 1] <- 'Desmoplastic'
Puleo[Puleo == 2] <- 'Stroma'
Puleo[Puleo == 3] <- 'Immune_Classical'
Puleo[Puleo == 4] <- 'Classical'
Puleo[Puleo == 5] <- 'Basal'

### Subclassifications Dataframe ###############################################
## Creation
Classification <- data.frame(row.names = Sample_Names, 
                             Collison = Collison, 
                             Moffit = Moffit, 
                             Bailey = Bailey,
                             Puleo = Puleo)

## Removal of Samples with NA
Classification <- Classification[complete.cases(Classification),]

### MCL Consensus Network ######################################################
# MCL_Res <- consensus.MCL(Classification,
#                          I.values = seq(1,15,1),
#                          outdir = "Test_0.01_4C",
#                          n.iter = 500,
#                          sim.method = "CohenKappa",
#                          pval.cut = 0.01)

# Silhouette <- data.frame()
# 
# for (i in 1:length(MCL_Res)) {
#    Silhouette[i,"minweighted"] = MCL_Res[[i]][["wsil.mean"]]
# }
# 
# optimalInflation= which.max (Silhouette$minweighted)
# MCL_Res[[optimalInflation]][['cl']]
# Consensus_Classes_Types <- data.frame(class = as.numeric(MCL_Res[[optimalInflation]][['cl']]), subtypes_names = names(MCL_Res[[optimalInflation]][['cl']]))


### Hypergeometric Testing for Selection of Consensus Samples ##################
## Count Number of Successes for each Consensus Group
Consensus_Classical_Counts <- as.data.frame(rowSums(Classification == 'Classical') + 
                                              rowSums(Classification == 'Progenitor'))

Consensus_Basal_Counts <- as.data.frame(rowSums(Classification == 'Basal') + 
                                          rowSums(Classification == 'Squamous') +
                                          rowSums(Classification == 'Stroma') +
                                          rowSums(Classification == 'Quasimesenchymal'))

Consensus_Exo_Counts <- as.data.frame(rowSums(Classification == 'Exocrine') + 
                                        rowSums(Classification == 'ADEX') +
                                        rowSums(Classification == 'Immune_Classical'))

Consensus_Imm_Counts <- as.data.frame(rowSums(Classification == 'Desmoplastic') + 
                                          rowSums(Classification == 'Immunogenic'))

## Hypergeometric Test for each Consensus Group
Consensus_Classical_Samples_Names <- rownames(Consensus_Classical_Counts)[phyper(as.vector(Consensus_Classical_Counts[,1]),4,4,4,lower.tail = FALSE) <= 0.015 ]
Consensus_Basal_Samples_Names <- rownames(Consensus_Basal_Counts)[phyper(as.vector(Consensus_Basal_Counts[,1]),4,4,4,lower.tail = FALSE) <= 0.015 ]
Consensus_Exo_Samples_Names <- rownames(Consensus_Exo_Counts)[phyper(as.vector(Consensus_Exo_Counts[,1]),4,4,4,lower.tail = FALSE) <= 0.15 ]
Consensus_Imm_Samples_Names <- rownames(Consensus_Imm_Counts)[phyper(as.vector(Consensus_Imm_Counts[,1]),4,4,4,lower.tail = FALSE) <= 0.25 ]


####Lolus Test####

# Norm_Data[,Consensus_Classical_Samples_Names] 
# Norm_Data[,c(Consensus_Basal_Samples_Names, Consensus_Exo_Samples_Names, Consensus_Imm_Samples_Names)]
# 
# 
# mod.t.test(cbind(Norm_Data[,c(Norm_Data[,Consensus_Classical_Samples_Names], Consensus_Basal_Samples_Names, Consensus_Exo_Samples_Names, Consensus_Imm_Samples_Names)]))
# 
# 
# 
# group1=length(Consensus_Classical_Samples_Names)
# group2=length(c(Consensus_Basal_Samples_Names, Consensus_Exo_Samples_Names, Consensus_Imm_Samples_Names))


### Core Sample RNA-Seq Data ###################################################
Core_Classical <- Norm_Data[colnames(Norm_Data) %in% Consensus_Classical_Samples_Names]
Core_Basal <- Norm_Data[colnames(Norm_Data) %in% Consensus_Basal_Samples_Names]
Core_Exo <- Norm_Data[colnames(Norm_Data) %in% Consensus_Exo_Samples_Names]
Core_Imm <- Norm_Data[colnames(Norm_Data) %in% Consensus_Imm_Samples_Names]

# Order is important for gene filtering
Core_Samples <- cbind(Core_Classical,Core_Basal,Core_Exo,Core_Imm)
Core_Samples <- as.matrix(Core_Samples)

### T Testing ##################################################################
Factor_Classical <- factor(c(rep(1,ncol(Core_Classical)),
                             rep(0,ncol(Core_Basal) + ncol(Core_Exo) + ncol(Core_Imm))))

Factor_Basal <- factor(c(rep(0,ncol(Core_Classical)),
                         rep(1,ncol(Core_Basal)),
                         rep(0,ncol(Core_Exo) + ncol(Core_Imm))))

Factor_Exo <- factor(c(rep(0,ncol(Core_Classical) + ncol(Core_Basal)),
                       rep(1,ncol(Core_Exo)),
                       rep(0,ncol(Core_Imm))))

Factor_Imm <-  factor(c(rep(0,ncol(Core_Classical) + ncol(Core_Basal) + ncol(Core_Exo)),
                        rep(1,ncol(Core_Imm))))

T_Test_Classical <- mod.t.test(Core_Samples,group = Factor_Classical)
T_Test_Basal <- mod.t.test(Core_Samples,group = Factor_Basal)
T_Test_Exo <- mod.t.test(Core_Samples,group = Factor_Exo)
T_Test_Imm <- mod.t.test(Core_Samples,group = Factor_Imm)

Core_Samples <- as.data.frame(t(Core_Samples))

PVal_Classical_Filtered_Genes <- rownames(T_Test_Classical[T_Test_Classical$p.value < 0.00001,])
PVal_Classical_Filtered_Core_Samples <- Core_Samples[colnames(Core_Samples) %in% PVal_Classical_Filtered_Genes]

PVal_Basal_Filtered_Genes <- rownames(T_Test_Basal[T_Test_Basal$p.value < 0.001,])
PVal_Basal_Filtered_Core_Samples <- Core_Samples[colnames(Core_Samples) %in% PVal_Basal_Filtered_Genes]

PVal_Exo_Filtered_Genes <- rownames(T_Test_Exo[T_Test_Exo$p.value < 0.001,])
PVal_Exo_Filtered_Core_Samples <- Core_Samples[colnames(Core_Samples) %in% PVal_Exo_Filtered_Genes]

PVal_Imm_Filtered_Genes <- rownames(T_Test_Imm[T_Test_Imm$p.value < 0.001,])
PVal_Imm_Filtered_Core_Samples <- Core_Samples[colnames(Core_Samples) %in% PVal_Imm_Filtered_Genes]

### AUC threshold ##############################################################
AUC_Threshold <- 0.9

### Classical AUC Selection ####################################################
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

### Basal AUC Selection ########################################################
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

### Exo AUC Selection ##########################################################
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

### Imm AUC Selection ##########################################################
Imm_Fitted_Values <- list()
for (i in c(1:ncol(PVal_Imm_Filtered_Core_Samples))) {
  vector <- fitted.values(glm(Factor_Imm~PVal_Imm_Filtered_Core_Samples[,i],family = binomial))
  Imm_Fitted_Values[[i]] <- vector
}

Imm_Aucs <- vector()
for (i in c(1:length(Imm_Fitted_Values))) {
  result <- auc(roc(Factor_Imm,Imm_Fitted_Values[[i]]))
  Imm_Aucs[[i]] <- result
}

AUC_Imm_Genes <- data.frame(colnames(PVal_Imm_Filtered_Core_Samples),Imm_Aucs)
rownames(AUC_Imm_Genes) <- AUC_Imm_Genes$colnames.PVal_Imm_Filtered_Core_Samples.
Final_Imm_Genes <- rownames(AUC_Imm_Genes[AUC_Imm_Genes$Imm_Aucs >= AUC_Threshold,])

###Final total genes############################################################
Final_Total_Genes <- c(Final_Classical_Genes,Final_Basal_Genes,Final_Exo_Genes,Final_Imm_Genes)
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

###Obtention of logfc for exocrine samples######################################
Exo_DGE <- DGEList(counts=Core_Samples_Filtered, group=Factor_Exo)
Exo_DGE <- estimateDisp(Exo_DGE)
Exo_LogFC <- exactTest(Exo_DGE,dispersion = 'common')

###Obtention of logfc for immunogenic samples###################################
Imm_DGE <- DGEList(counts=Core_Samples_Filtered, group=Factor_Imm)
Imm_DGE <- estimateDisp(Imm_DGE)
Imm_LogFC <- exactTest(Imm_DGE,dispersion = 'common')

###Creation of centroids########################################################
# identical(rownames(Exo_LogFC$table),rownames(Basal_LogFC$table))
Centroids <- data.frame('Classical' = Classical_LogFC$table$logFC,
                        'Basal' = Basal_LogFC$table$logFC,
                        'Exo' = Exo_LogFC$table$logFC,
                        'Imm' = Imm_LogFC$table$logFC,
                        'hgnc_symbol' = rownames(Exo_LogFC$table),
                        row.names = rownames(Exo_LogFC$table))



# write.table(Centroids,'Centroids_4S_4C_0.9AUC',sep = '\t')
