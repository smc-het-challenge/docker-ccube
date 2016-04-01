#!/usr/bin/env Rscript

library(dplyr)
library(ccube)

library(colorspace)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
myColors <- gg_color_hue(10)

# args <- commandArgs(trailingOnly = TRUE)
# vcfFile <- as.character(args[1])
# batternbergFile <- as.character(args[2])

vcfFile <- "Tumour2/Tumour2.mutect.vcf"
batternbergFile <- "Tumour2/Tumour2.battenberg.txt"


ssm_file <- "ssm_data.txt"
cnv_file <- "cnv_data.txt"

shellCommandMutectSmcHet <- paste(
  "create_ccfclust_inputs.py -v mutect_smchet",
  " -b ", batternbergFile,
  " -c ", 1,
  " --output-cnvs ", cnv_file,
  " --output-variants ", ssm_file,
  " ", vcfFile, sep = ""
)
system(shellCommandMutectSmcHet, intern = TRUE)
ssm <- read.delim(ssm_file,
                  stringsAsFactors = F)
cnv <- read.delim(cnv_file,
                  stringsAsFactors = F)

ssm$major_cn = 1
ssm$minor_cn = 1
ssm$cn_frac = 1
ssm$mu_r <- NULL
ssm$mu_v <- NULL
ssm$cn_ref <- NULL
ssm$cn_tot <- NULL

cnv <- na.omit(cnv)
cnv <- dplyr::filter(cnv, ssms!="" )

if (nrow(cnv)>0) {
  cnvtmp1 <- strsplit(as.character(cnv$ssms), ";")
  for (j in seq_len(nrow(cnv))) {
    if (length(cnvtmp1[[j]])==0) { next }
    cnvtmp1[[j]] = paste(cnvtmp1[[j]], cnv[j,]$frac, sep="," )
  }
  cnvtmp1 <- unlist(cnvtmp1)
  cnvtmp2 <- Reduce(
    rbind, strsplit(cnvtmp1, ",")
  )

  if (is.null(dim(cnvtmp2) )) {
    cnvtmp2 = as.data.frame(t(cnvtmp2), stringsAsFactors=F)
  } else {
    cnvtmp2 = as.data.frame(cnvtmp2, stringsAsFactors=F)
  }

  for (j in 2:ncol(cnvtmp2)) {
    cnvtmp2[,j] = as.numeric(cnvtmp2[,j])
  }

  ssm <- dplyr::left_join(ssm, cnvtmp2, by=c("id"="V1"))
  ssm$major_cn <- ssm$V3
  ssm$minor_cn <- ssm$V2
  ssm$cn_frac <- ssm$V4

  ssm$V2 <- NULL
  ssm$V3 <- NULL
  ssm$V4 <- NULL

  ssm[is.na(ssm[,5]), 5] = 1
  ssm[is.na(ssm[,6]), 6] = 1
  ssm[is.na(ssm[,7]), 7] = 1
}

clonalCnFrac <- sum(ssm$cn_frac==1)/nrow(ssm)
ssm <- dplyr::filter(ssm, cn_frac==1)


# maxSnv <- 30000
# if (nrow(ssm) > maxSnv) {
#   ssm <- dplyr::sample_n(ssm, maxSnv)
# }

ssm$normal_cn = 2
ssm <- dplyr::rename(ssm, ref_counts=a, total_counts=d)
ssm <- dplyr::mutate(ssm, var_counts=total_counts-ref_counts, mutation_id = gene)
ssm$purity <- GetPurity(ssm)
cellularity <- unique(ssm$purity)
write.table(cellularity, file = "1A.txt", sep = "\t", row.names = F, 
  col.names = F, quote = F)