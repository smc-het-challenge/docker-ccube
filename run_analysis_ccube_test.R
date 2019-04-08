#!/usr/bin/env Rscript

rm(list = ls())

library(dplyr)
library(ccube)
library(doParallel)
library(foreach)

registerDoParallel(cores = detectCores())



args <- commandArgs(trailingOnly = TRUE)
vcfFile <- as.character(args[1])
batternbergFile <- as.character(args[2])

vcfFile <- "~/Downloads/P1-noXY/P1-noXY.mutect.vcf"
batternbergFile <- "~/Downloads/P1-noXY/P1-noXY.battenberg.txt"

ssm_file <- "ssm_data.txt"

# Parse vcf file
vcfParserPath <- dir(path = getwd(), pattern = "create_ccfclust_inputs.py", full.names = T)
shellCommandMutectSmcHet <- paste(
  vcfParserPath,
  " -v mutect_smchet",
  " -c ", 1,
  " --output-variants ", ssm_file,
  " ", vcfFile, sep = ""
)
system(shellCommandMutectSmcHet, intern = TRUE)
ssm <- read.delim(ssm_file, stringsAsFactors = F)

# Parse Battenberg CNA data
cna <- read.delim(batternbergFile, stringsAsFactors = F)
ssm <- ParseSnvCnaBattenberg(ssm, cna) 


# Estimate purity and write 1A.txt
cellularity <- GetPurity(ssm)
ssm$purity <- cellularity
write.table(cellularity, file = "1A.txt", sep = "\t", row.names = F, col.names = F, quote = F)


# Run Ccube 
numOfClusterPool = 1:10
numOfRepeat = 1
ccubeRes <- RunCcubePipeline(ssm = ssm, 
                             numOfClusterPool = numOfClusterPool, numOfRepeat = numOfRepeat,
                             runAnalysis = T, runQC = T, returnAll = T)

# write 1B.txt
uniqLabels <- unique(ccubeRes$res$label)
write.table(length(uniqLabels), file = "1B.txt", sep = "\t", row.names = F, col.names=F, quote = F)

if (!is.matrix(ccubeRes$res$full.model$responsibility)) {
  mutR <- data.frame(ccubeRes$res$full.model$responsibility)
  colnames(mutR) <- "cluster_1"
} else {
  mutR <- data.frame(ccubeRes$res$full.model$responsibility[, sort(uniqLabels)]) 
  colnames(mutR) <- paste0("cluster_", seq_along(uniqLabels))
}

ssmCcube <- ccubeRes$ssm
ssmCcube <- cbind(ssmCcube, mutR)
ssmCcube$label <- apply(mutR, 1, which.max)

tt <- ssmCcube
tt <- tt[order(as.numeric(gsub("[^\\d]+", "", tt$id, perl=TRUE))), ]

clusterCertainty <- as.data.frame(table(tt$label), stringsAsFactors = F)
clusterCertainty <- rename(clusterCertainty, cluster = Var1, n_ssms = Freq)
clusterCertainty$proportion <- ccubeRes$res$full.model$ccfMean[sort(uniqLabels)][as.integer(clusterCertainty$cluster)] * cellularity
clusterCertainty$cluster <- seq_along(uniqLabels) 
write.table(clusterCertainty, file = "1C.txt", sep = "\t", row.names = F, col.names=F, quote = F)

write.table(tt$label, file = "2A.txt", sep = "\t", row.names = F, col.names=F, quote = F)

RR = as.matrix(tt[, paste0("cluster_", seq_along(sort(uniqLabels)))])
coAssign <- Matrix::tcrossprod(RR)
diag(coAssign) <- 1
write.table(coAssign, file = "2B.txt", sep = "\t", row.names = F, col.names=F, quote = F)

