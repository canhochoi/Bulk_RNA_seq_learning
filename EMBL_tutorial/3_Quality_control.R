#follow workshop https://www.ebi.ac.uk/training/materials/bioinformatics-immunologists-materials/rna-sequencing/bulk-rna-seq-analysis/
#follow workshop at https://mperalc.gitlab.io/bulk_RNA-seq_workshop_2021/qc.html
#https://www.ebi.ac.uk/training/materials/bioinformatics-immunologists-materials/rna-sequencing/bulk-rna-seq-analysis/


library(airway)
library(reshape2)
library(DESeq2)
library(edgeR)
library(vsn) # for meanSdPlot
library(tidyr)

data(airway)

head(airway)

metadata(airway)

#information about how the gene counts were generated
metadata(rowRanges(airway))

#raw count matrix
#row is gene names (represented as ENSEMBL gene IDs)
head(assay(airway))

#sample information such as cell line , treatment with dex (dexamethasone)
head(as.data.frame(colData(airway)))

#summary of all columns by 
#1st summing all rows in each column
#2nd perform statistic across all columns
summary(colSums(assay(airway)))


#look at each sample
summary(assay(airway))

#plot histogram

library(ggplot2)

df <- as.data.frame(assay(airway))
df$gene_ID <- rownames(df)

df_long <- df %>% pivot_longer(cols = 1:8, names_to = "Sample", values_to = "Count")
#altertnatively, use melt in dshape2
#df_long <- melt(df, id = ncol(df))
#colnames(df_long)[2:3] = c("Sample", "Count")
summary(df_long)


#raw count follows negative binomial or Poisson distribution
ggplot(df_long, aes(x = Count)) + 
  geom_histogram() +
  labs(x = "gene count", y = "Occurence", title = "Histogram") +
#  scale_x_log10() +  
  theme_minimal()


#builiding DESeqDataSet object

#this way
#dds <- DESeqDataSet(airway, design = ~ cell + dex)

#or this way
dds <- DESeqDataSetFromMatrix(countData = assay(airway),
                              colData = colData(airway),
                              design = ~ cell + dex)

dds

#filter lowly expressed gene which has less than 1 count per million in at least 2 samples

dds <- dds[rowSums(cpm(counts(dds)) > 1) >= 2, ]

# normalize
# use variance-stabilizing transformation, log-transformed 
# blind = FALSE --> transformations are affected by experimental design 
# --> important for downstream analysis

dds_rlog <- rlog(dds,blind = FALSE)
dds_vst <- vst(dds,blind = FALSE)
head(assay(dds_rlog), 3)
head(assay(dds_vst), 3)

dds <-  estimateSizeFactors(dds)

#  log2(n+1) of raw counts
meanSdPlot(assay(normTransform(dds))) 

meanSdPlot(assay(dds_rlog))

meanSdPlot(assay(dds_vst))

# do PCA

getPCs <-  function(dds){
  #function to do PCA
  df <- assay(dds)
  #transpose so that row is number of samples (i.e., trials)
  #column is number of measurement (i.e., each gene experssion)
  pca <- prcomp(t(df), retx = TRUE)
  percentvariance <- round(100*(pca$sdev^2)/sum(pca$sdev^2))
  pcs <- as.data.frame(pca$x)
  pcs <- cbind(pcs,colData(dds))
  pcs <- as.data.frame(pcs)
  pcs <- list(pcs, percentvariance)
  names(pcs) <- c("pcs","percentVar")
  return(pcs)
}

pca_rlog <- getPCs(dds_rlog)

ggplot(pca_rlog$pcs, aes(PC1, PC2, color = dex, shape = cell)) + 
  geom_point() + 
  xlab(paste("PC1 of", pca_rlog$percentVar[1], "% variance")) +
  ylab(paste("PC2 of", pca_rlog$percentVar[2], "% variance")) 

#without any normalization
pca_raw <- getPCs(dds)

ggplot(pca_raw$pcs, aes(PC1, PC2, color = dex, shape = cell)) + 
  geom_point() + 
  xlab(paste("PC1 of", pca_raw$percentVar[1], "% variance")) +
  ylab(paste("PC2 of", pca_raw$percentVar[2], "% variance")) 


#get dominant genes influencing PCAs
getLoadings = function(dds){
  #return eigen vectors 
  df <- assay(dds)
  pca <- prcomp(t(df), retx = TRUE)
  df_rot <- as.data.frame(pca$rotation)
  return(df_rot)
}

loadings_rlog <- getLoadings(dds_rlog) 
# Annotating gene names
library(AnnotationDbi) # gene annotation
library(org.Hs.eg.db) # gene annotation

loadings_rlog$symbol <- mapIds(org.Hs.eg.db,
                              keys=row.names(loadings_rlog),
                              column="SYMBOL",
                              keytype="ENSEMBL",
                              multiVals="first")

loadings_PC1 <- pivot_longer(loadings_rlog[,c("PC1","symbol")], cols = "PC1", names_to = "PC1", values_to = "loadings")
loadings_PC1 <- loadings_PC1[order(loadings_PC1$loadings, decreasing = TRUE),]

#select top 10 genes
loadings_PC1$symbol[1:10]