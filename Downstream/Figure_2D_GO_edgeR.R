library(ggplot2)
library(edgeR)
library(readxl)
library(ggpubr)
library(sva)
library(biomaRt)
library(msigdbr)
library(ape)
library(GSVA)
library(dplyr)

#has patient information
#https://gtpb.github.io/ADER18F/pages/tutorial2.html



setwd("D:/OneDrive - Nanyang Technological University/Postdoc/Job_application/Australia/David Lynn/to_work/lynnlab-covirs")

######## Read in files ######## 
raw_counts.df = readRDS("raw_RNASeq_counts.RDS")
meta_data = readRDS("Sample_meta_data_RNASeq.RDS")

all(colnames(raw_counts.df)==rownames(meta_data))


#keep data for V0 and V1
meta_data = meta_data[meta_data$Timepoint %in% c("V0","V1"),]
#meta_data = meta_data[meta_data$Run %in% c("Run_1"),]

#keep participants that have results for both V0 and V1
tab_list = table(meta_data$Participant.ID)
tab_list = tab_list[tab_list == 2]
rownames(tab_list)
meta_data = meta_data[meta_data$Participant.ID %in% rownames(tab_list), ]

#keep participants that use ChAdOx1-S
meta_data = meta_data[meta_data$first.second.vaccine == "ChAdOx1", ]

# Function to swap rows i and j in a data frame and their row names
swap_rows <- function(df, row_names, i, j) {
  temp_row <- df[i, ]
  temp_name <- row_names[i]
  df[i, ] <- df[j, ]
  row_names[i] <- row_names[j]
  df[j, ] <- temp_row
  row_names[j] <- temp_name
  return(list(df = df, row_names = row_names))
}

# Loop to perform swapping
for (i in seq(1, nrow(meta_data), by = 2)) {
  if (i + 1 <= nrow(meta_data) && !(i %in% c(9,11,21,23))) {
    result <- swap_rows(meta_data, row.names(meta_data), i, i + 1)
    meta_data <- result$df
    row.names(meta_data) <- result$row_names
  }
}

counts_select = raw_counts.df[,rownames(meta_data)]

#build a DGE class
y = DGEList(counts= counts_select,group=meta_data$Timepoint)

#obtain cpm = each count/total count *1e6
y_cpm = cpm(y)

# choose values in cpm are greater than 3
thresh =  y_cpm > 3

#keep genes that have more than 15 TRUES in each row of thresh
keep = rowSums(thresh) > 5

#subset the rows of count data to keep the more highly expressed genes
y_keep = y[keep, ,keep.lib.sizes=FALSE]

#normalize
y_keep = calcNormFactors(y_keep,method = "TMM") 

#study the effect of Timepoint
meta_data$Participant.ID = factor(meta_data$Participant.ID)
meta_data$Timepoint = factor(meta_data$Timepoint)

#design = model.matrix(~ Participant.ID + Timepoint, meta_data)
#Participant.ID = as.factor(meta_data$Participant.ID)

#this works ok
#design = model.matrix(~ Timepoint, meta_data)

design = model.matrix(~Participant.ID  + Timepoint, meta_data)

design

#colnames(design)[2:15] = levels(meta_data$Participant.ID)[2:15]
colnames(design)[1:15] = levels(meta_data$Participant.ID)[1:15]

all(rownames(design) == colnames(y_keep))

design

y_keep = estimateDisp(y_keep, design, robust=TRUE)
#Type ?estimateDisp() to see more information

#This shows the curve fitting to reestimate the dispersion
plotBCV(y_keep)

fit = glmQLFit(y_keep, design, robust = TRUE)
qlf = glmQLFTest(fit, coef = "TimepointV1")
topTags(qlf, nrow(qlf))

idx = order(qlf$table$PValue)
cpm(y_keep)[idx[1:10],]
#fit = glmFit(y_keep, design, robust = TRUE)
#qlf = glmLRT(fit, coef = "TimepointV1")

is.de <- decideTests(qlf, p.value=0.05)
summary(is.de)

plotMD(qlf)
abline(h=c(-1,1), col="blue")


qlf.df = data.frame(qlf$table)
qlf.df$padj <- p.adjust(qlf.df$PValue, method="BH")

#extracting significant gene only
qlf.df_significant = qlf.df[(qlf.df$padj < 0.05),]
qlf.df_significant = qlf.df_significant[order(qlf.df_significant$logFC, decreasing = TRUE), ]



library(org.Hs.eg.db)
library(clusterProfiler)
library(AnnotationHub)
library(ensembldb)
library(biomaRt)

#get the background genes

ensembl = useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

ensids_background = rownames(y)
ens2gene_background = getBM(attributes=c('external_gene_name','ensembl_gene_id','gene_biotype','entrezgene_id'), 
                            filters = 'ensembl_gene_id', 
                            values = ensids_background, 
                            mart = ensembl) 
ens2gene_background = ens2gene_background[!duplicated(ens2gene_background$ensembl_gene_id),]
#remove NA in gene annotation
#ens2gene_background = ens2gene_background[!is.na(ens2gene_background$entrezgene_id), ]

rownames(ens2gene_background) = ens2gene_background$ensembl_gene_id
head(ens2gene_background)


#get the significant genes ID from ensemble
library(biomaRt)
ensembl = useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes = listAttributes(ensembl)
ensids = rownames(qlf.df_significant)
ens2gene = getBM(attributes=c('external_gene_name','ensembl_gene_id','gene_biotype','entrezgene_id'), 
                 filters = 'ensembl_gene_id', 
                 values = ensids, 
                 mart = ensembl) 
ens2gene = ens2gene[!duplicated(ens2gene$ensembl_gene_id),]
rownames(ens2gene) = ens2gene$ensembl_gene_id
head(ens2gene)

#this is important to match gene names and ID with qlf.df_significant
ens2gene_reorder = ens2gene[rownames(qlf.df_significant),]

qlf.df_significant <- cbind(qlf.df_significant, ens2gene_reorder)

all(rownames(qlf.df_significant) == qlf.df_significant$ensembl_gene_id)
  
##Look at FDR spread
ggplot(qlf.df_significant, aes(x=padj)) +
  geom_histogram(bins=100) +
  theme_classic()

## Run GO enrichment analysis 
# https://hbctraining.github.io/DGE_workshop/lessons/09_functional_analysis.html

#### subset into up and down-regulated genes ####

#separate into up-regulated only
gene_sig_up <- qlf.df_significant %>% 
  dplyr::filter(padj <= 5E-2 & logFC > 0)

#separate into down-regulated only
gene_sig_down <- qlf.df_significant %>% 
  dplyr::filter(padj <= 5E-2 & logFC <= 0)

#run enrichment for up-regulated gene
ego_up <- enrichGO(gene = gene_sig_up$ensembl_gene_id, 
                   universe = ens2gene_background$ensembl_gene_id,
                   keyType = "ENSEMBL",
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)

## Output results from GO analysis to a table
cluster_summary_up <- data.frame(ego_up)

#run enrichment for down-regulated gene
ego_down <- enrichGO(gene = gene_sig_down$ensembl_gene_id, 
                     universe = ens2gene_background$ensembl_gene_id,
                     keyType = "ENSEMBL",
                     OrgDb = org.Hs.eg.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05, 
                     readable = TRUE)

## Output results from GO analysis to a table
cluster_summary_down <- data.frame(ego_down)


#plot bar plot
# https://bioc.ism.ac.jp/packages/3.7/bioc/vignettes/enrichplot/inst/doc/enrichplot.html

cluster_summary_up %>% 
  dplyr::filter(p.adjust <= 0.05) %>% 
  
  ggplot(aes(x=reorder(Description, -log10(qvalue)), 
             y=-log10(qvalue))) +
  geom_col() +
  theme_classic() +
  #Some more customization to pretty it up
  #Flip x and y so long labels can be read
  coord_flip() +
  #fix labels
  labs(y="-log(qvalue)",
       x="",
       title = "Significant pathway enriched")

cluster_summary_down %>% 
  dplyr::filter(p.adjust <= 0.05) %>% 
  
  ggplot(aes(x=reorder(Description, -log10(qvalue)), 
             y=-log10(qvalue))) +
  geom_col() +
  theme_classic() +
  #Some more customization to pretty it up
  #Flip x and y so long labels can be read
  coord_flip() +
  #fix labels
  labs(y="-log(qvalue)",
       x="",
       title = "Significant pathway reduced")


#### end ####


#extract the first two rows of up-regulated gene

#cluster_combined = rbind(cluster_summary_up[1:3,], cluster_summary_down[c(1,24,26,31,36),])

cluster_combined = rbind(cluster_summary_up[1:3,], cluster_summary_down[c(1,22,23,34,40),])
cluster_combined$regulation = c(rep("up",3), rep("down",5))

ggplot(cluster_combined, aes(x=reorder(Description, -log10(qvalue)),
                             y=-log10(qvalue), fill = regulation)) +
  geom_col() +
  theme_classic() +
  scale_fill_manual(values = c("up" = "brown", "down" = "lightblue")) +
  coord_flip() + 
  #fix labels
  labs(y="-log(qvalue)",
       x="",
       title = "Pathway analysis")

