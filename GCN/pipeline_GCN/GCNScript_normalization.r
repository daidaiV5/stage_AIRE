#! /usr/bin/Rscript

####
# Working version by yuping.DAI
# R version 3.6.3
# Last modified by yuping on 29/06/2021
####

Args <- commandArgs(TRUE)

############## Ajout vérification paramètres
## Recall the arguments for checking ##############################
for (i in 1:length(Args)){
  print(Args[i])
}

## Default setting when no arguments passed
if(length(Args) < 1) {
  Args <- c("--help")
}
 
## Help section
if("--help" %in% Args) {
  cat("The R Script gcnScript.r
 
      Arguments:
      --in_mat=mat                     - file
      --out=radical                    - out file radical
      --dir=directory		               - out directory
      --number_sample               - the number of samples
      --norm=normal                    - bool if we normalize or not
      --help                           - print this text
      
 
      Example:
      GCNScript_nomalization.r --in_mat=matrice.csv --out=output --dir=OUTDIR \n\n")
 
  q(save="no")
}

##############
## Parse arguments (we expect the form --arg=value)
#
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(Args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
print(Args)
print("###") 
print(argsL)
print("###") 
print(argsDF)
print("###") 
parseMult <- function(x) strsplit(x, ",")
#
############

#############
# Data default MATRIX
#
if(is.null(argsL$in_mat)) {
                        print("no input matrix")
                        print(argsL$in_mat)
                        stop("--in-mat required")
} else {file_mat <- argsL$in_mat}

######
# Output default
#
if(is.null(argsL$out)) {
                        print("no radical given: using input name")
                        radical <- strsplit(file_mat,"_")[[1]][1]       ## RADICAL of output files
} else {radical <- argsL$out}
#############
######
# Output folder default
#
if(is.null(argsL$dir)) {
                        print("no output directory: using current dir")
                        folder = "./"                                   ## Current folder
} else {folder <- paste(argsL$dir, "/", sep="")
       dir.create(folder)
}

# number of sampling

if(is.null(argsL$number_sample)) {
  print("number of samping is null")

} else {number_sample <- argsL$number_sample}


if(is.null(argsL$norm)) {
  print("normalize")
  
} else {norm <- argsL$norm}


# make sure it's an integer
# samping <- round(as.double(samping))

############# END INPUT STEP



#######################################################################
# 0) Load libraries 
#######################################################################
# if (!requireNamespace("BiocManager", quietly = TRUE) || !require("DESeq2")){
#     install.packages("BiocManager", dependencies = TRUE)
#     BiocManager::install("DESeq2", dependencies = TRUE)
#    }
# if (!require("ggplot2")) {
#    install.packages("ggplot2", dependencies = TRUE)
#    library("ggplot2")
#    }
# if (!require("WGCNA")){
#     BiocManager::install("WGCNA", dependencies = TRUE) 
# }

library("DESeq2")
#library("ggplot2")

  
##### LOAD DATA #######################################################
cts <- as.matrix(read.csv(file_mat,
                               header=TRUE, sep=",", row.names = 1)) ## To do: parameterize row.names -- h-coded to "gene_id" ***

col.name=colnames(cts)
cts=cts[,sample(col.name,number_sample,replace = FALSE)]

coldata=data.frame(colnames(cts),nul_col='NULL',row.names = 1)

all(rownames(coldata) == colnames(cts))

##### Normale the DATA #######################################################

if(norm=='True'){
print('begin to normalize')
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ 1)

dds <- DESeq(dds)
resMF <- results(dds)

keep <- rowSums(counts(dds)) >= 20							# Parameter: keep only counts >= 20, if keep only counts >=10,removeBatchEffect will return value negetif
dds <- dds[keep,]
vst <- vst(dds, blind=FALSE)
matrix_final=assay(vst)
}else{
matrix_final=cts
}

file = paste(radical,"_matrice_nor",".csv", sep = "")	
write.table(matrix_final, paste(folder,file, sep = ""), sep = ",")

