#! /usr/bin/Rscript

####
# Working version by E. Corel
# Further adapted by J. Teulière
# Last modified by yuping on 97/05/2021
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
      --thr=threshold		               - 1-Pval Student threshold for significqnt pearson correlations (default 0.95)
      --pthr=pthreshold                - Pearson correlation threshold above which to keep edges (default 0.5)
      --min=minimum_reads              - Minimum number of reads to determine correlations (default 20)
      --random_result                  -   
      --repeats
      --help                           - print this text
 
      Example:
      gcnScript.r --in_mat=EXAMPLE_raw_count_sub.csv --annot=EXAMPLE_annot_sub.csv --out=output --dir=OUTDIR \n\n")
 
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
                        stop("--in-mat required")
} else {file_mat <- argsL$in_mat}
#############

if(is.null(argsL$repeats)) {
  print("poblem")
  stop("--in-mat required")
} else {repeat_time<- argsL$repeats}


if(is.null(argsL$random_result )) {
print("no input ")
stop("--random_result required")
  } else {random_result <- argsL$random_result}


# Data default ANNOTATION
#
# if(is.null(argsL$annot)) {
#                         print("no input annot")
#                         stop("--annot required")
# } else {file_annot <- argsL$annot}

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



######
# Threshold for Student test of Pearson correlations
#
if(is.null(argsL$thr)) {
  print("Student threshold default value: 0.95")
  thr <- 0.95
} else {thr <- argsL$thr}
# make sure it's a float
thr <- as.double(thr)

#####for PCA : type of ellipse 
# if(is.null(argsL$type_ellipse)) {
#   print("type_ellipse:t")
#   type_ellipse <- 't'
# } else {
#   type_ellipse <- argsL$type_ellipse
#   print(type_ellipse)
#   }


# if(is.null(argsL$thr_ellipse)) {
#   print("threshold default value : 0.95")
#   thr_ellipse <- 0.95
# } else {thr_ellipse <- argsL$thr_ellipse}
# thr_ellipse <- as.double(thr_ellipse)


######
# Threshold for Pearson correlations
#
if(is.null(argsL$pthr)) {
  print("Pearson threshold default value: 0.5") 
  pthr <- 0.50
} else {pthr <- argsL$pthr}
# make sure it's a float
pthr <- as.double(pthr)

######
# Threshold for reads
#
if(is.null(argsL$min)) {
  print("Minimum reads threshold default value: 20") 
  min <- 20
} else {min <- argsL$min}
# make sure it's an integer
min <- round(as.double(min))

# if(is.null(argsL$samping)) {
#   print("number of samping: 15") 
#   samping <- 5
# } else {samping <- argsL$samping}
# make sure it's an integer
# samping <- round(as.double(samping))

############# END INPUT STEP

initFolder = "./"

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
if (!require("WGCNA")){
    BiocManager::install("WGCNA", dependencies = TRUE) 
}
# if (!require("sp")){
#   BiocManager::install("sp", dependencies = TRUE) 
# }
#if (!require("vegan")){
#  BiocManager::install("vegan", dependencies = TRUE) 
#}


#library("DESeq2")
#library("ggplot2")

  
##### LOAD DATA #######################################################
vst2 <- as.matrix(read.csv(file=paste(initFolder, file_mat, sep = ""),
                               header=TRUE, sep=",", row.names = 1)) ## To do: parameterize row.names -- h-coded to "gene_id" ***


# batcoldata <- read.csv(file=paste(initFolder, file_annot, sep = ""),
#                          header = TRUE, sep =",", row.names = 1)              ## To do: parameterize row.names -- h-coded to "samples" ***
# 
# 
#   
# #colnames(batCts) <- sub("_raw","", colnames(batCts))				      # alt 1 : unused
# 
# colnames(batCts) <- sub("\\.","-", colnames(batCts))
# 
# colnames(batCts) <- lapply(colnames(batCts) , function(x){gsub("\\.", "-", x)})      # alt2 : unused
#   
# #all(rownames(batcoldata) %in% colnames(batCts))                                      # tests for same order in col names in the matrix, 
#                                                                                       # and row names in the annotation file
# #all(rownames(batcoldata) == colnames(batCts))
#   
# #batCts <- batCts[, rownames(batcoldata)]
# 
# #head(batCts)
# all(rownames(batcoldata) == colnames(batCts))
# dds <- DESeqDataSetFromMatrix(countData = batCts,
#                                 colData = batcoldata,
#                                 design = ~ 1)
# 
# #stop("Ok: arguments checked")
# 
# dds <- DESeq(dds)
#   
# keep <- rowSums(counts(dds)) >= min							# Parameter: keep only rows with total counts >= min ***
# dds <- dds[keep,]
# ##### END LOAD DATA #####################################################
# 
# ##### COMPUTE NORMALISATION #############################################
# data_list=c(row.names(batcoldata))
# a=sample(data_list,samping)
# dds=dds[,a]
# 
# vst <- vst(dds, blind=FALSE)

#rld <- rlog(dds, blind=FALSE)								# alt

# load WGCNA library #########



library("WGCNA")
options(stringsAsFactors = FALSE);
allowWGCNAThreads(nThreads=20)
## END load WGCNA library ####

# vst <- assay(vst)

repeat_time=as.integer(repeat_time)

df=read.csv(random_result,sep = ',',header = FALSE)

for (i in 1:repeat_time) {
  x=c(unlist(df[i,]))
  vst=vst2[,x]
  
  i=as.character(i)
  vst <- t(vst)
  dir.create(paste(folder,"/",i,"/",sep = ""))
  file = paste(radical,"_vst.txt", sep = "")	
  write.table(vst, paste(folder,"/",i,"/",file, sep = ""))
  
  cor_pval <- corAndPvalue(vst,alternative = "two.sided", method = "pearson", nThreads = 24)
  
  file = paste(radical,"_gcn_edges_cor_",pthr,".csv", sep = "")					# To do :  parameterize output name ***
  exportNetworkToCytoscape(cor_pval$cor, 
                           edgeFile = paste(folder, file, sep = ""),
                           threshold = pthr)
  cor_pval <- NULL
  #pvalue_mat <- NULL
  gc()
  
}


##### END COMPUTE NORMALISATION #########################################
  
##### PEARSON CORRELATIONS NETWORK AND PVALUES ##########################
# STEP 1 : compute Pearson correlation and pvalues ######################

  
#threshold = 0.8									# alt Section : unused
#cor_t <- cor_pval$cor[(abs(cor_pval$cor[]) > threshold)]				#
#qplot(cor_t, xlim = c(0.8,1))		      						#
#file = "young_gcn_final_cor_0.9.txt.cc_distrib"					#
#cc_distri <- read.csv(paste(folder, file, sep = ""), sep ="\t")			# end alt section
  
# STEP 2 : export Pearson's correlation network #########################

# Substep a : Export edges # threshold modified #########################
										# H-Coded parameter: thr ***
# file = paste(radical,"_gcn_edges_cor_",pthr,".csv", sep = "")					# To do :  parameterize output name ***
# exportNetworkToCytoscape(cor_pval$cor, 
#                            edgeFile = paste(folder, file, sep = ""),
#                            threshold = pthr)

# Substep b : Export nodes #############################################
# file = paste(radical,"_gcn_nodes_cor_",pthr,".txt", sep = "")					# To do :  parameterize output name ***
# exportNetworkToCytoscape(cor_pval$cor, 
#                            nodeFile = paste(folder, file, sep = ""),
#                            threshold = pthr)

# Substep c : export pvalues matrix ####################################
#pvalue_mat <- cor_pval$p
#dimnames(pvalue_mat) <- dimnames(cor_pval$cor)

#### !! P-values are expressed as 1 - Pval, all values should be > 0.95:
#pvalue_mat <- 1- pvalue_mat 

#file = paste(radical,"_gcn_edges_pval_",thr,".txt", sep = "")			# To do :  parameterize output name ***
#exportNetworkToCytoscape(pvalue_mat, 
#                           edgeFile = paste(folder, file, sep = ""), 
#                           threshold = thr)						
# cor_pval <- NULL
# pvalue_mat <- NULL
# gc()

##### END PEARSON CORRELATIONS NETWORK AND PVALUES #####################
  
# ##### SOFT THRESHOLD ###################################################
# # STEP 1 : define power range ##########################################
# powers = c(c(1:10), seq(from = 12, to=20, by=2))					# h-coded power range 1-10 and 12-20 by 2 ***
# # STEP 2 : find soft threshold #########################################
# soft_t <- pickSoftThreshold(vst, networkType = "unsigned", RsquaredCut = 0.9,		# h-coded Rsq = 0.9 ***
#                                powerVector = powers)
# s_t <- soft_t$powerEstimate
#    
# # if (is.na(s_t)){									# alt section : unused
# #   soft_t <- pickSoftThreshold(vst, networkType = "unsigned", RsquaredCut = 0.8,	#
# # #                               powerVector = powers)	       		     		#
# # #   s_t <- soft_t$powerEstimate							#
# # # } 	     										#
# # # if (is.na(s_t)){									#
# # #   soft_t <- pickSoftThreshold(vst, networkType = "unsigned", RsquaredCut = 0.70,	#
# # #                               powerVector = powers)		 	       		#
# # #   s_t <- soft_t$powerEstimate 	      						#
# # # } 	     										#
# #   											# end alt section
# 
# # ENDING : Plot the results ############################################
# sizeGrWindow(9, 5)
# par(mfrow = c(1,2));
# cex1 = 0.9;
# # Scale-free topology fit index as a function of the soft-thresholding power
# plot(soft_t$fitIndices[,1], -sign(soft_t$fitIndices[,3])*soft_t$fitIndices[,2],
#         xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,unsigned R^2",
#         type="n",main = paste("Scale independence"));text(soft_t$fitIndices[,1],
#                                                           -sign(soft_t$fitIndices[,3])*soft_t$fitIndices[,2],labels=powers,cex=cex1,col="red");
# # this line corresponds to using an R^2 cut-off of h = 0.90
# abline(h=0.90,col="red")
# # Mean connectivity as a function of the soft-thresholding power
# plot(soft_t$fitIndices[,1], soft_t$fitIndices[,5],xlab="Soft Threshold (power)",
#         ylab="Mean Connectivity", type="n",
#         main = paste("Mean connectivity"))
# text(soft_t$fitIndices[,1], soft_t$fitIndices[,5],labels=powers, cex=cex1,col="red")
# ######
# # 
# s_t = 14 #make sure its the same for all categories					# h-coded value of soft threshold = 14 ***
# # 
# # ##### ADJACENCY NETWORK
#  #compute adjacency network
#  adj <- adjacency(vst, selectCols = NULL, type = "unsigned",
#                   power = s_t, corFnc = "cor",
#                   corOptions = list(use = 'p'))
#  #####
#  
#  ##### TOM NETWORK AND MODULES
#  #Compute TOM and dissTOM
#  tom <- TOMsimilarity(adj, TOMType = "unsigned")
#  
#  # Call the hierarchical clustering function
#  geneTree = hclust(as.dist(1-tom), method = "average");
#  # Plot the resulting clustering tree (dendrogram)
#  sizeGrWindow(12,9)
#  plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
#       labels = FALSE, hang = 0.04);
#  
#  minModuleSize = 20;
#  # Module identification using dynamic tree cut:
#  dynamicMods = cutreeDynamic(dendro = geneTree, distM = 1-tom,
#                              deepSplit = 2, pamRespectsDendro = FALSE,
#                              minClusterSize = minModuleSize);
#  table(dynamicMods)
#  
#  # Convert numeric lables into colors
#  dynamicColors = labels2colors(dynamicMods)
#  table(dynamicColors)
#  # Plot the dendrogram and colors underneath
#  sizeGrWindow(8,6)
#  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
#                      dendroLabels = FALSE, hang = 0.03,
#                      addGuide = TRUE, guideHang = 0.05,
#                      main = "Gene dendrogram and module colors")
#  
#  #Select all modules in TOM to export 
#  probes = rownames(adj)
#  inModule = is.finite(match(dynamicColors, dynamicColors))
#  modProbes = probes[inModule]
#  dimnames(tom) <- dimnames(adj)
#  
#  #export TOM network
#  percentile = .999									# h-coded percentile = 0.999 ***
#  h_t <- quantile(as.vector(tom), percentile)
#  h_t <- unname(h_t)
#  
#  file <- paste(radical,"_gcn_tom_edges_s%d_h%.3f.txt", sep = "" )			# To do : parameterize output name ***
#  edges <- sprintf(file, s_t, percentile)
#  
#  exportNetworkToCytoscape(tom,
#                                 edgeFile = paste(folder, edges, sep=""),
#                                 weighted = TRUE,
#                                 threshold = h_t);
#  
#  file <- paste(radical,"_gcn_tom_nodes_s%d.txt", sep = "")				# To do : parameterize output name ***
#  nodes <- sprintf(file, s_t, percentile)
#  exportNetworkToCytoscape(tom,
#                           nodeFile = paste(folder, nodes, sep=""),
#                           weighted = TRUE,
#                           threshold = 0,
#                           nodeNames = modProbes,
#                           nodeAttr = dynamicColors[inModule]);
# 
# ##### END of MAIN Procedure
#  
#  adj <- NULL
#  tom <- NULL
#  gc()
  

