#! /usr/bin/Rscript

####
# Working version by E. Corel
# Further adapted by J. Teulière
# Last modified on 03/11/2021
####

## Parameters:
setwd("/media/jerome/Seagate Portable Drive/AIRE_backup/GCN_papers")

file_mat = "./HUMAN/GTEx_dataset/large_read_matrices/by_organ/13756__Whole_Blood_reads_master.csv"
file_annot = "./HUMAN/GTEx_dataset/annotations/13756__Whole_annotation.csv"
min = 20

dataname = "AGE"

#######################################################################
# 0) Load libraries 
#######################################################################

library("DESeq2")
library("ggplot2")
library("factoextra")
library("NbClust")
  
##### LOAD DATA #######################################################
batCts <- as.matrix(read.csv(file=file_mat,
                               header=TRUE, sep=",", row.names = "gene_id"))          ## To do: parameterize row.names -- h-coded to "gene_id" ***
batcoldata <- read.csv(file=file_annot,
                         header = TRUE, sep =",", row.names = "samples")              ## To do: parameterize row.names -- h-coded to "samples" ***

head(batcoldata)
#colnames(batCts) <- sub("_raw","", colnames(batCts))				      # alt 1 : unused

colnames(batCts) <- sub("\\.","-", colnames(batCts))
 
colnames(batCts) <- lapply(colnames(batCts) , function(x){gsub("\\.", "-", x)})      # alt2 : unused
#head(batCts)
ncol(batCts)
#all(rownames(batcoldata) %in% colnames(batCts))
  
#all(rownames(batcoldata) == colnames(batCts))
  
#batCts <- batCts[, rownames(batcoldata)]
all(rownames(batcoldata) == colnames(batCts))
dds <- DESeqDataSetFromMatrix(countData = batCts,
                                colData = batcoldata,
                                design = ~ 1)

#stop("Ok: arguments checked")

dds <- DESeq(dds)
  
keep <- rowSums(counts(dds)) >= min							# Parameter: keep only counts >= 10 ***
dds <- dds[keep,]
##### END LOAD DATA #####################################################
  
##### COMPUTE NORMALISATION #############################################
vst <- vst(dds, blind=FALSE)

df2=assay(vst)
df=t(scale(df2))

nc <- NbClust(df, distance="euclidean", 
              min.nc=2, max.nc=4, method="kmeans")


z <- plotPCA(vst, intgroup = dataname)
z + geom_text(aes(label = name), size = 2, colour = "black")

##### END COMPUTE NORMALISATION #########################################
  
