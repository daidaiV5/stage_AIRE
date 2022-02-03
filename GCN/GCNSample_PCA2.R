#! /usr/bin/Rscript

####
# Working version by E. Corel
# Further adapted by J. Teulière
# Last modified on 03/11/2021
####

## Parameters:
setwd("/Users/daiyuping/Desktop/stageM2/echantionage")


file_mat = "/Users/daiyuping/Desktop/stageM2/echantionage/dataset_test/male_20_sub/matrices/wholeblood_male_20-29_sub_20.csv"
file_annot = "/Users/daiyuping/Desktop/stageM2/echantionage/dataset_test/male_20_sub/annot/wholeblood_male_20-29_annotation.csv"
min = 20

dataname = "AGE"

#######################################################################
# 0) Load libraries 
#######################################################################

library("DESeq2")
library("ggplot2")
library("sp")

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

batcoldata['AGE']='A'
dds <- DESeqDataSetFromMatrix(countData = batCts,
                              colData = batcoldata,
                              design = ~1)

#stop("Ok: arguments checked")

dds <- DESeq(dds)

keep <- rowSums(counts(dds)) >= min							# Parameter: keep only counts >= 10 ***
dds <- dds[keep,]
##### END LOAD DATA #####################################################

##### COMPUTE NORMALISATION #############################################
vst <- vst(dds, blind=FALSE)
vst

data=plotPCA(vst, intgroup = dataname, returnData=TRUE)
p=ggplot(data, aes(PC1, PC2,  shape=dataname))+   geom_point(size=3) +
  stat_ellipse(type = "t",level = 0.95)

ggplot(data, aes(PC1, PC2,  shape=dataname))+   geom_point(size=3) +
  stat_ellipse(type = "norm",level = 0.95)

ggplot(data, aes(PC1, PC2,  shape=dataname))+   geom_point(size=3) +
  stat_ellipse(type = "t",level = 0.95)



# Extract components
build <- ggplot_build(p)$data
points <- build[[1]]
ell <- build[[2]]

# Find which points are inside the ellipse, and add this to the data
dat <- data.frame(
  in.ell = as.logical(point.in.polygon(points$x, points$y, ell$x, ell$y))
)
data['in.ell']=dat

# Plot the result
p=ggplot(data, aes(PC1, PC2,  shape=dataname)) +
  geom_point(aes(col = in.ell)) +
  stat_ellipse()+geom_text(aes(label=name),size=2)


ggsave(p,filename = paste(name_test,".pdf",sep=""),width = 12,height = 9)

#### END ACP analysis ########################################################

data_list_supp=c(row.names(data)[which(data$in.ell==FALSE)])
data_list_neto=c(row.names(data)[which(data$in.ell==TRUE)])

batCts_clear=batCts[,data_list_neto]
batcoldata_clear=batcoldata[data_list_neto,]



##### END load data clear #########################################
  
