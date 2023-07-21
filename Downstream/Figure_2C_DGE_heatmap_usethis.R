library(ggplot2)
library(edgeR)
library(readxl)
library(ggpubr)
library(sva)
library(biomaRt)
library(msigdbr)
library(ape)
library(GSVA)
library(PCAtools)
library(dplyr)

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


qlf.df = data.frame(qlf$table)
qlf.df$FDR <- p.adjust(qlf.df$PValue, method="BH")
qlf.df_significant = qlf.df[(qlf.df$FDR < 0.05),]
qlf.df_significant = qlf.df_significant[order(qlf.df_significant$logFC, decreasing = FALSE), ]


# Get log2 counts per million for filtered gene from qlf.df_significant 
logcounts_y_keep = cpm(y_keep, log=TRUE)
logcounts_y_keep = logcounts_y_keep[rownames(logcounts_y_keep) %in% rownames(qlf.df_significant), ]

logcounts_y_keep = logcounts_y_keep[rownames(qlf.df_significant), ]

#change name of column to just V1 or V0
c_name = c()
for (name in colnames(logcounts_y_keep)){
  if (grepl("V1", name)){
    c_name = append(c_name, "V1")
  }
  else {
    c_name = append(c_name, "V0")
  }
}
colnames(logcounts_y_keep) = c_name


#center and scale each column to get Z-score then transpose
logcounts_scaled = t(apply(logcounts_y_keep, 1, scale))
colnames(logcounts_scaled) = c_name

#same as

#z_score = t(scale(t(logcounts_y_keep)))


library(ComplexHeatmap)

h <- Heatmap(logcounts_scaled, cluster_rows = T, 
             column_labels = colnames(logcounts_scaled), name="Z-score",
             show_row_names = FALSE,
             cluster_columns = T)
h

ht = draw(h)
rorder = row_order(ht)

#to see which row has change in Z score
c_array = array(NA, dim = c(nrow(logcounts_scaled)))
for (i in 1:nrow(logcounts_scaled)){
  c_array[i] = sum(logcounts_scaled[rorder[i], seq(1, ncol(logcounts_scaled), by = 2)] < 0)
  
}

#reorganzie the row order
row_order = c(rorder[282:438], rorder[1:281])

h_rorder = Heatmap(logcounts_scaled, row_order = row_order, 
                   show_row_names = FALSE, cluster_columns = T, 
                   name="Z-score")
h_rorder