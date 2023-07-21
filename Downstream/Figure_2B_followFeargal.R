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

setwd("D:/OneDrive - Nanyang Technological University/Postdoc/Job_application/Australia/David Lynn/to_work/lynnlab-covirs")

######## Read in files ######## 
raw_counts.df = readRDS("raw_RNASeq_counts.RDS")
meta_data = readRDS("Sample_meta_data_RNASeq.RDS")

all(colnames(raw_counts.df)==rownames(meta_data))

table(meta_data$Timepoint,meta_data$first.second.vaccine)

meta_data = meta_data[meta_data$Timepoint %in% c("V0","V1"),]
meta_data = meta_data[meta_data$first.second.vaccine %in% c("ChAdOx1"),]


col_id <- c()
colname_remove <- c('COVIRS_4_V0')
row_id <- c()

for (colname in colname_remove){
  #col_id <- append(col_id, which(colnames(raw_counts_keep.df) == colname))
  row_id <- append(row_id, which(rownames(meta_data) == colname))
}

#subset meta_data too 
meta_data <- meta_data[-row_id,]

#only choose columns correspond to new meta_data
raw_counts_keep.df <- raw_counts.df[, rownames(meta_data)]

all(colnames(raw_counts_keep.df) == rownames(meta_data))


#convert to DGElist object
dgeObj <- DGEList(raw_counts_keep.df)

keep <- rowSums(cpm(dgeObj)>3) > 10
table(keep)
dgeObj <- dgeObj[keep, , keep.lib.sizes=FALSE]

#normalize
dgeObj <- calcNormFactors(dgeObj,method = "TMM") 

# Get log2 counts per million
#logcounts <- cpm(dgeObj, log=TRUE, normalized.lib.sizes = TRUE)

logcounts <- cpm(dgeObj, log=TRUE)

mod = model.matrix(~ Timepoint,data=meta_data)
mod0 = model.matrix(~1,data=meta_data)

#number of surrogate variable
n.sv = num.sv(logcounts,mod,method = "be")

svobj = sva(logcounts,mod,mod0,n.sv=n.sv)
sv_matrix = svobj$sv
rownames(sv_matrix) = colnames(logcounts)

#do adjustment on logcounts 
X = cbind(mod, sv_matrix)
#solve(t(X) %*% X) is inverse of (transpose of X times dot product of X)
Hat = solve(t(X) %*% X) %*% t(X)
beta = (Hat %*% t(logcounts))
rm(Hat)
P = ncol(mod)
adj_logcounts = logcounts - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),])

mdsplot_adjusted = plotMDS(adj_logcounts, main="Vaccine day_adjusted") 

pcoa_out_corrected = pcoa(as.dist(mdsplot_adjusted$distance.matrix))

pf_data.df = data.frame(meta_data,"PC1" = mdsplot_adjusted$x,"PC2" = mdsplot_adjusted$y)

#Plots MDS plot for adjusted data
ggplot(pf_data.df,aes(x = PC1 , y = PC2,color=Timepoint,shape=first.second.vaccine))+geom_point()+
  xlab(paste("PC1:",round(pcoa_out_corrected$values$Relative_eig[1]*100,2),"% variation"))+
  ylab(paste("PC2:",round(pcoa_out_corrected$values$Relative_eig[2]*100,2),"% variation"))+
  scale_shape_manual(values=c(16,23,24))+
  scale_alpha_manual(values=c(0.75,0.9))+
  scale_size_manual(values=c(2,2))+theme_classic()