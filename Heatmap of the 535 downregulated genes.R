#################0.init environment#################
rm(list = ls())
suppressMessages(library(pheatmap))
suppressMessages(library(ComplexHeatmap))

#################1.plot heatmap#################
rawdata <- read.csv('sample_cpm.csv',header = T)
rawdata <- as.data.frame(rawdata)
rownames(rawdata) <- rawdata[,1]
rawdata <- rawdata[,-1]
rawdata <- as.matrix(rawdata)
pheatmap(rawdata, scale = "row",show_rownames = T,cluster_cols = F,border_color = F)
heat <- pheatmap(rawdata,scale = "row",show_rownames = F,cluster_cols = F,border_color = F)

