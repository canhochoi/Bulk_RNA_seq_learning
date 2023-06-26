library(airway)
library(DESeq2)


#https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

geo_mean <- function(mat){
  #geometric mean
  geo_m <- numeric(nrow(mat))
  for (i in 1:nrow(mat)){
    arr <- mat[i,]
    geo_m[i] <- exp(mean(log(arr[arr>0])))
  }
  return(geo_m)
}

data(airway)

df <- as.data.frame(assay(airway))

#1 Calculate geometric mean
geo_mean_array <- geo_mean(df)

#2 for each sample (i.e., column), calculate the normalization factor
# by dividing with its respective row in geometric mean

#create geometric mean matrix
repeated_geo_mean <- rep(geo_mean_array, times = ncol(df))

df_geo_mean <- as.data.frame(matrix(repeated_geo_mean, nrow = nrow(df), ncol = ncol(df)))

#divide to get normalization factor

df_norm_factor <- df/df_geo_mean

#3 get size factor for each sample

median_sample <- apply(df_norm_factor, 2, function(x) median(x, na.rm = TRUE))


#compare with DESeq package
#or this way
dds <- DESeqDataSetFromMatrix(countData = assay(airway),
                              colData = colData(airway),
                              design = ~ cell + dex)

dds <- DESeq(dds)

#not the same 
#need to check
#source code for DESseq
#https://github.com/mikelove/DESeq2/blob/devel/R/core.R



dds$sizeFactor/median_sample

s <- estimateSizeFactors(dds)

#plot dispersion
plotDispEsts(dds)

#log2FoldChange and p-values 
res_condition <- results(dds)
res_condition

#meaning of results column
mcols(res_condition)

#specify the untreat as baseline
res_condition = results(dds, contrast=c("dex","trt","untrt"))
res_condition

summary(res_condition)

#False Discovery Rate of 10% --> adjusted p-value < 0.1

## Of those with abs(log2FC)>1, how many genes are significant at 5% FDR level?

sum(res_condition$padj<0.05 & abs(res_condition$log2FoldChange) > 1, na.rm = TRUE)

## Summary of results setting FDR to 5% and abs(log2FC) to 1

res_condition = results(dds, contrast=c("dex","trt","untrt"), alpha=0.05, lfcThreshold=1)
res_condition

## After setting lfcThreshold and alpha: of those with abs(log2FC)>1, how many genes are significant at 5% FDR level?

sum(res_condition$padj<0.05 & abs(res_condition$log2FoldChange) > 1, na.rm = TRUE)


sum(res_condition$padj<0.01 & abs(res_condition$log2FoldChange) > 1, na.rm = TRUE)

sum(res_condition$padj<0.01 , na.rm = TRUE)


#log2 fold change is not accuarate for low gene expression and high variance within a condition

res_condition_shrink <- lfcShrink(dds, contrast = c("dex", "trt", "untrt"), type = "normal", 
                                  res = res_condition[order(rownames(res_condition)),],
                                  lfcThreshold = 1)

sum(res_condition_shrink$padj<0.05 & abs(res_condition_shrink$log2FoldChange) > 1, na.rm = TRUE)

res_condition_shrink <- res_condition_shrink[order(res_condition_shrink$padj),]
res_condition_shrink


#compare and contrast two cell types with contrast option
res_cell <- results(dds, contrast=c("cell", "N061011", "N61311"),lfcThreshold = 1, alpha = 0.05)
res_cell_shrink <- lfcShrink(
  dds,
  contrast = c("cell", "N061011", "N61311"),
  type = "normal",
  res = res_condition[order(rownames(res_condition)), ],
  lfcThreshold = 1
)

res_cell_shrink[order(res_cell_shrink$padj),]

res_condition_shrink$symbol = mapIds(org.Hs.eg.db,
                                     keys=row.names(res_condition_shrink),
                                     column="SYMBOL",
                                     keytype="ENSEMBL",
                                     multiVals="first")

res_condition_shrink$Ensembl_ID = row.names(res_condition_shrink)

head(res_condition_shrink)

#export txt result
# write.table(res_condition_shrink, file = "DiffExpr_all_genes.txt", quote = F, sep = "\t",row.names = F, col.names = T)
# res_condition_shrink %>%  
#   as.data.frame() %>%
#   filter(padj<0.05 & log2FoldChange > 1) %>%
#   write.table(file = "DiffExpr_padj05_log2FC1.txt", quote = F, sep = "\t",row.names = F, col.names = T)


# Getting the count table in ggplot-friendly format
library(ggplot2)
library(ggrepel) # repel labels

geneCounts <- plotCounts(dds, gene=res_condition_shrink$Ensembl_ID[1], intgroup=c("dex","cell"), returnData=TRUE)

ggplot(geneCounts, aes(x = dex, y = count, color = cell)) + 
  scale_y_log10() + 
  geom_point(position=position_jitter(width=.1,height=0), size=3) + 
  geom_text_repel(aes(label =  cell)) +
  theme_bw() + 
  ggtitle( res_condition_shrink$symbol[1]) + 
  ylab("log10(counts)")



geneCounts = plotCounts(dds, gene=res_condition_shrink$Ensembl_ID[10800], intgroup=c("dex","cell"), returnData=TRUE)

ggplot(geneCounts, aes(x=dex, y=count, color=cell)) +
  scale_y_log10() +
  geom_point(position=position_jitter(width=.1,height=0), size=3) + 
  geom_text_repel(aes(label = cell)) +
  theme_bw() + 
  ggtitle( res_condition_shrink$symbol[10800]) +
  ylab("log10(counts)")


# ordering and mapping gene names
res_cell = results(dds, contrast=c("cell", "N080611", "N61311"))
res_cell = lfcShrink(dds, contrast=c("cell", "N080611", "N61311"), 
                     type = "normal",
                     res=res_cell[order(rownames(res_cell)),])

res_cell = res_cell[order(res_cell$padj),]
res_cell$symbol = mapIds(org.Hs.eg.db,
                         keys=row.names(res_cell),
                         column="SYMBOL",
                         keytype="ENSEMBL",
                         multiVals="first")

res_cell$Ensembl_ID = row.names(res_cell)
# getting counts
geneCounts = plotCounts(dds, gene=res_cell$Ensembl_ID[1], intgroup=c("dex","cell"), returnData=TRUE)
# plotting
ggplot(geneCounts, aes(x=cell, y=count, color=dex)) +
  scale_y_log10() +
  geom_point(position=position_jitter(width=.1,height=0), size=3) + 
  geom_text_repel(aes(label = dex)) +
  theme_bw() + 
  ggtitle( res_cell$symbol[1]) +
  ylab("log10(counts)")


#to get several gene
norm_counts = as.data.frame(counts(dds, normalized=TRUE)[res_condition_shrink$Ensembl_ID[1:10], ])
norm_counts$Ensembl_ID = rownames(norm_counts)

# Getting to long format
# norm_counts = norm_counts %>%
#   as.tibble() %>%
#   pivot_longer(cols = starts_with("SRR"),names_to = "Run",values_to = "Normalised_counts")

norm_counts = norm_counts %>%
  pivot_longer(cols = starts_with("SRR"),names_to = "Run",values_to = "Normalised_counts")

# Mapping gene names
norm_counts$symbol = mapIds(org.Hs.eg.db,
                            keys=norm_counts$Ensembl_ID,
                            column="SYMBOL",
                            keytype="ENSEMBL",
                            multiVals="first")  

df = merge(norm_counts,colData(dds),by="Run")

df %>% 
  as.data.frame() %>%
  ggplot( aes(x=symbol, y=Normalised_counts, color=dex)) +
  scale_y_log10() +
  geom_point(position=position_jitter(width=.1,height=0), size=3) + 
  theme_bw() + 
  ylab("log10(counts)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))


#MA plot
DESeq2::plotMA(res_condition_shrink, ylim=c(-5,5))

#volcano plot
# adding another section for coloring based on p-value and log2FC threshold with ggplot
res_color =  res_condition_shrink %>% 
  as.tibble() %>%
  mutate(color = padj < 0.01 & abs(log2FoldChange) >= 1)
res_color[which(is.na(res_color$color)),"color"] = FALSE # setting the NAs to FALSE

# adding labels to plot
res_color = res_color %>% 
  arrange(padj) %>% 
  mutate(labels = "")

res_color$labels[1:10] = res_color$symbol[1:10]

ggplot(res_color) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = color, label = labels)) +
  geom_text_repel(aes(label = labels, x = log2FoldChange, y = -log10(padj))) +
  ggtitle("Mov10 overexpression") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  +
  theme_bw()


#heatmap
library(pheatmap) # heatmaps
mat = as.data.frame(counts(dds, normalized=TRUE)[res_condition_shrink$Ensembl_ID[1:10], ])
mat = log10(mat)
mat = t(scale(t(mat))) # Need to scale per gene, hence the transposition
rownames(mat) = mapIds(org.Hs.eg.db,
                       keys=row.names(mat),
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")

dds_rlog <- rlog(dds,blind = FALSE)
df = as.data.frame(colData(dds_rlog)[,c("cell","dex")])
pheatmap(mat, annotation_col=df)


#batch effect
library(sva) # remove unwanted variation
mat <- counts(dds, normalized=TRUE)

threshold <- mat > 3
keep <- rowSums(threshold) > 3

mat_filtered <- mat[keep,]

mod <- model.matrix(~ dex, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))

#if not filtered, there will rows with all 0 in all columns
# --> cause error in svaseq package
svseq <- svaseq(mat_filtered, mod, mod0)
#number of surrogate variable is also called latent variable
str(svseq)

par(mfrow=c(3,1),mar=c(2,4,2,1))
stripchart(svseq$sv[,1] ~ dds$cell,vertical=TRUE,main="SV1")
abline(h=0)
stripchart(svseq$sv[,2] ~ dds$cell,vertical=TRUE,main="SV2")
abline(h=0)
stripchart(svseq$sv[,3] ~ dds$cell,vertical=TRUE,main="SV3")
abline(h=0)


ddssva <- dds # making a copy
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
ddssva$SV3 <- svseq$sv[,3]
design(ddssva) <- ~ SV1 + SV2 + SV3 + dex

ddssva <- DESeq(ddssva)

res_sva <- results(ddssva, contrast=c("dex","trt","untrt"),
                  lfcThreshold = 1, alpha = 0.05)

res_sva <- lfcShrink(ddssva, 
                    contrast=c("dex","trt","untrt"), 
                    type = "normal",
                    res=res_sva[order(rownames(res_sva)),],
                    lfcThreshold = 1)

res_sva <- res_sva[order(res_sva$padj),]

res_sva$symbol <- mapIds(org.Hs.eg.db,
                        keys=row.names(res_sva),
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")

res_sva$Ensembl_ID = row.names(res_sva)

#only the treatment

dds_reduced <- dds

design(dds_reduced) <- ~  dex

dds_reduced <- DESeq(dds_reduced)


res_reduced$Ensembl_ID <- row.names(res_reduced)

res_reduced <- results(dds_reduced, contrast=c("dex","trt","untrt"),
                      lfcThreshold = 1, alpha = 0.05)

res_reduced <- lfcShrink(dds_reduced, 
                        contrast=c("dex","trt","untrt"), 
                        type = "normal",
                        res=res_reduced[order(rownames(res_reduced)),],
                        lfcThreshold = 1)

res_reduced <- res_reduced[order(res_reduced$padj),]

res_reduced$symbol <- mapIds(org.Hs.eg.db,
                            keys=row.names(res_reduced),
                            column="SYMBOL",
                            keytype="ENSEMBL",
                            multiVals="first")

res_reduced$Ensembl_ID <- row.names(res_reduced)

# Now compare all three results

res_condition_shrink

res_sva

res_reduced

cell <- colData(ddssva)$cell
# We use the normalised data, as previously with PCA
adjusted <- ComBat(assay(ddssva), batch=cell)


df <- adjusted
pca <- prcomp(t(df), retx = TRUE)

percentVar <- (pca$sdev)^2 / sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pcs <- as.data.frame(pca$x)
pcs <- cbind(pcs,colData(dds_rlog))
pcs <- as.data.frame(pcs)
pcs <- list(pcs, percentVar)

names(pcs) <- c("pcs","percentVar")

ggplot(pcs$pcs, aes(PC1, PC2, color=dex, shape=cell)) +   geom_point(size=3) +
  xlab(paste0("PC1: ",pcs$percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",pcs$percentVar[2],"% variance")) + theme_bw()