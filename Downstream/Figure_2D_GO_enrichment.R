library(DESeq2)
library(edgeR)
library(AnnotationDbi) # gene annotation
library(org.Hs.eg.db) # gene annotation

library(ggplot2)
library(edgeR)
library(readxl)
library(ggpubr)
#library(sva)
library(biomaRt)
library(msigdbr)
library(ape)
library(GSVA)
#library(PCAtools)
library(dplyr)

library(clusterProfiler)

#http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#de
#https://mperalc.gitlab.io/bulk_RNA-seq_workshop_2021/diffExpr.html#inspecting-the-differential-expression-results


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


meta_data$Participant.ID = as.factor(meta_data$Participant.ID)
meta_data$Timepoint = as.factor(meta_data$Timepoint)

data = DESeqDataSetFromMatrix(countData = raw_counts.df[,rownames(meta_data)],
                       colData = meta_data,
                       design = ~ Participant.ID + Timepoint)

head(counts(data))

data = data[ rowSums(edgeR::cpm(counts(data)) > 3) > 5, ]
nrow(data)


data <- DESeq(data)

plotDispEsts(data)

res_condition = results(data)
res_condition

levels(data$Timepoint)

summary(res_condition)

res_condition$symbol = mapIds(org.Hs.eg.db,
                              keys=row.names(res_condition),
                              column="SYMBOL",
                              keytype="ENSEMBL",
                              multiVals="first")

res_condition$Ensembl_ID = row.names(res_condition)

head(res_condition)

#remove NA in gene annotation

res_condition = res_condition[!is.na(res_condition$symbol), ]

#separate into up-regulated only
gene_sig_up <- res_condition[res_condition$padj <= 5E-2 & res_condition$log2FoldChange > 0, ]


#separate into down-regulated only
gene_sig_down <- res_condition[res_condition$padj <= 5E-2 & res_condition$log2FoldChange <= 0, ]


#get the background genes
ensembl = useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
ensids_background = rownames(res_condition)
ens2gene_background = getBM(attributes=c('external_gene_name','ensembl_gene_id','gene_biotype','entrezgene_id'), 
                            filters = 'ensembl_gene_id', 
                            values = ensids_background, 
                            mart = ensembl) 
ens2gene_background = ens2gene_background[!duplicated(ens2gene_background$ensembl_gene_id),]
rownames(ens2gene_background) = ens2gene_background$ensembl_gene_id
head(ens2gene_background)


ego_up <- enrichGO(gene = gene_sig_up$Ensembl_ID, 
                   universe = ens2gene_background$ensembl_gene_id,
                   keyType = "ENSEMBL",
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)

## Output results from GO analysis to a table
cluster_summary_up <- data.frame(ego_up)

ego_down <- enrichGO(gene = gene_sig_down$Ensembl_ID, 
                   universe = ens2gene_background$ensembl_gene_id,
                   keyType = "ENSEMBL",
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)
## Output results from GO analysis to a table
cluster_summary_down <- data.frame(ego_down)


cluster_summary_up %>% 
  dplyr::filter(p.adjust <= 0.05 & qvalue < 1e-15) %>% 
  
  ggplot(aes(x=reorder(Description, -log10(qvalue)), #Reorder gene sets by k/K values
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
  dplyr::filter(p.adjust <= 0.05 & qvalue < 1e-1) %>% 
  
  ggplot(aes(x=reorder(Description, -log10(qvalue)), #Reorder gene sets by k/K values
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

#extract the first two rows of up-regulated gene

cluster_combined = rbind(cluster_summary_up[c(1,2,23),], cluster_summary_down[c(1,20,26,50,53),])
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
       