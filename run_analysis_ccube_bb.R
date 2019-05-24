#!/usr/bin/env Rscript

rm(list = ls())

library(dplyr)
library(ccube)
library(doParallel)
library(foreach)
options(stringsAsFactors = F)

registerDoParallel(cores = detectCores())

args <- commandArgs(trailingOnly = TRUE)
vcfFile <- as.character(args[1])
batternbergFile <- as.character(args[2])
purityFile <- as.character(args[3])


ssm_file <- "ssm_data.txt"

# Parse vcf file
vcfParserPath <- "create_ccfclust_inputs.py"
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

# Get Battenberg purity
cellularity <-read.delim(purityFile, stringsAsFactors=FALSE)$cellularity
ssm$purity <- cellularity
write.table(cellularity, file = "1A.txt", sep = "\t", row.names = F, col.names = F, quote = F)

ssm$vaf = ssm$var_counts/ssm$total_counts
allSsm <- ccube:::CheckAndPrepareCcubeInupts(ssm)

# 1st filter: major_cn == 0
ProblemSsmIDs <- filter(allSsm, major_cn == 0 )$id
ssm <- filter( allSsm, major_cn > 0 )


# 2nd filter: fp_qval, rough_ccf1, rough_ccf0
ssm<- mutate(rowwise(ssm), 
             fp_pval = binom.test(var_counts, total_counts, 5e-2, alternative = "great")$p.value, 
             fp_qval = p.adjust(fp_pval, method = "hochberg")
             )

ssm<- mutate(rowwise(ssm), 
             rough_ccf1 =  
               MapVaf2CcfLinear( vaf, purity, normal_cn, total_cn, 1) 
             )

ssm<- mutate(rowwise(ssm), 
             rough_ccf0 =  
               MapVaf2CcfLinear( vaf, purity, normal_cn, total_cn, 0) 
             )


ProblemSsmIDs <- c(ProblemSsmIDs, filter(ssm, fp_qval > 0.05  &
                                              rough_ccf1 < 0.2 &
                                              rough_ccf0 > -90
                                         )$id 
                   )

ssm = filter( ssm, !id %in% ProblemSsmIDs)

# Run Ccube 
numOfClusterPool = 1:7
numOfRepeat = 5
ccubeRes <- RunCcubePipeline(ssm = ssm, 
                             numOfClusterPool = numOfClusterPool, numOfRepeat = numOfRepeat,
                             runAnalysis = T, runQC = T, multiCore = T)

MakeCcubeStdPlot(ccubeRes$ssm, ccubeRes$res)


# write 1B.txt
uniqLabels <- unique(ccubeRes$res$label)
write.table(length(uniqLabels), file = "1B.txt", sep = "\t", row.names = F, col.names=F, quote = F)

if (!is.matrix(ccubeRes$res$full.model$responsibility)) {
  mutR <- data.frame(ccubeRes$res$full.model$responsibility)
  colnames(mutR) <- "cluster_1"
} else {
  mutR <- data.frame(ccubeRes$res$full.model$responsibility[, sort(uniqLabels)]) 
  colnames(mutR) <- paste0("cluster_", seq_along(uniqLabels) )
}

ssmCcube <- ccubeRes$ssm[, "id"]
ssmCcube <- cbind(ssmCcube, mutR)
ssmCcube <- cbind(ssmCcube, matrix(0, nrow = nrow(ssmCcube), 1 ) )
colnames(ssmCcube)[ncol(ssmCcube)] <- paste0("cluster_", length(uniqLabels) + 1 )
ssmCcube$label <- apply(mutR, 1, which.max)

ProblemSsm <- data.frame( id =  ProblemSsmIDs )
ProblemSsmMutR <-  data.frame ( matrix(0, nrow = nrow(ProblemSsm), ncol = length(uniqLabels)  ) )
colnames(ProblemSsmMutR) <- paste0("cluster_", seq_along(uniqLabels))
ProblemSsm <-  cbind(ProblemSsm, ProblemSsmMutR)
ProblemSsm <- cbind(ProblemSsm, matrix(1, nrow = nrow(ProblemSsm), 1 ) )
colnames(ProblemSsm)[ncol(ProblemSsm)] <- paste0("cluster_", length(uniqLabels) + 1 )
ProblemSsm$label <- length(uniqLabels) + 1

tt <- rbind( ssmCcube, ProblemSsm)
tt <- tt[order(as.numeric(gsub("[^\\d]+", "", tt$id, perl=TRUE))), ]
rownames(tt)<-NULL

clusterCertainty <- as.data.frame(table(tt$label), stringsAsFactors = F)
clusterCertainty <- rename(clusterCertainty, cluster = Var1, n_ssms = Freq)
clusterCertainty$proportion <- c ( ccubeRes$res$full.model$ccfMean[as.integer(clusterCertainty$cluster[-nrow(clusterCertainty)] )] * cellularity, 
                                   0)
write.table(clusterCertainty, file = "1C.txt", sep = "\t", row.names = F, col.names=F, quote = F)
write.table(tt$label, file = "2A.txt", sep = "\t", row.names = F, col.names=F, quote = F)

tt$id <- NULL
tt$label <- NULL
RR = as.matrix(tt)
coAssign <- Matrix::tcrossprod(RR)
diag(coAssign) <- 1
write.table(coAssign, file = "2B.txt", sep = "\t", row.names = F, col.names=F, quote = F)

