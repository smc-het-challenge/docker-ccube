#!/usr/bin/env Rscript

library(dplyr)
library(ccube)

library(colorspace)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
myColors <- gg_color_hue(10)

args <- commandArgs(trailingOnly = TRUE)
vcfFile <- as.character(args[1])
batternbergFile <- as.character(args[2])
purityFile <- as.character(args[3])



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

allSsm <- dplyr::mutate(rowwise(ssm), chr =  strsplit(gene, split = "_")[[1]][1],
                        pos = strsplit(gene, split = "_")[[1]][2]) 

clonalCnFrac <- sum(allSsm$cn_frac==1)/nrow(allSsm)
ssm <- dplyr::filter(allSsm, cn_frac==1 & !chr %in% c("x", "y") )
problemSsm <- dplyr::filter(allSsm, cn_frac!=1 | chr %in% c("x", "y") )

# maxSnv <- 30000
# if (nrow(ssm) > maxSnv) {
#   ssm <- dplyr::sample_n(ssm, maxSnv)
# }

ssm$normal_cn = 2
ssm <- dplyr::rename(ssm, ref_counts=a, total_counts=d)
ssm <- dplyr::mutate(ssm, var_counts=total_counts-ref_counts, mutation_id = gene)
#ssm$purity <- GetPurity(ssm)
#cellularity <- unique(ssm$purity)
cellularity <-read.delim(purityFile, stringsAsFactors=FALSE)$cellularity
ssm$purity <- cellularity
write.table(cellularity, file = "1A.txt", sep = "\t", row.names = F, 
  col.names = F, quote = F)


library(doParallel)
library(foreach)

registerDoParallel(cores = detectCores())

kk <- 6
if (kk > nrow(ssm)){
  kk <- nrow(ssm) - 1
}
rr <- 5
iterSetting <- data.frame( sort(rep(seq(1, kk, length.out = kk), rr)) )

results <- foreach(n = 1:nrow(iterSetting), 
  .combine = c, .packages = "ccube") %dopar% {
  library(ccube)
  k <- iterSetting[n, ]
  list(ccube_m6(ssm, epi=1e-3,
         init=k, tol = 1e-10, maxiter = 1e3,
         fit_mult = T, fit_hyper = F, use = "use_base", verbose = F))
}

maxLbIndex <- which.max(Map( function(x) max(x$L), results))
lb <- unlist(Map( function(x) max(x$L), results))
res <- results[[maxLbIndex]]

ssm$ccube_ccf_mean <- res$mu[res$label]
ssm$ccube_mult <- res$full.model$bv  
ssm <- mutate(rowwise(ssm),
              vaf = var_counts/(var_counts+ref_counts),
              ccube_ccf = MapVaf2CcfPyClone(vaf,
                          purity,
                          normal_cn,
                          major_cn + minor_cn,
                          major_cn + minor_cn,
                          ccube_mult,
                          constraint=F))

uniqLabels <- unique(res$label)

if (length(uniqLabels) > 1 & min(ssm$ccube_ccf_mean) < 0.3) {
  write.table(length(uniqLabels)-1, file = "1B.txt", sep = "\t", row.names = F, col.names=F, quote = F)

} else {
  write.table(length(uniqLabels), file = "1B.txt", sep = "\t", row.names = F, col.names=F, quote = F)
}



if (!is.matrix(res$R)) {
  mutR <- data.frame(res$R)
  colnames(mutR) <- "cluster_1"
} else {
  mutR <- data.frame(res$R[, sort(uniqLabels)]) 
  colnames(mutR) <- paste0("cluster_", seq_along(uniqLabels))
}

ssm <- cbind(ssm, mutR)
ssm$label <- apply(mutR, 1, which.max)


## Post assign problem SSMs-- temporal solution

MapVaf2CcfTest <- function(x, t, cv, bv, frac,
                              epi = 1e-3, constraint = T,
                              lower = -Inf, upper = Inf) {
  
  if(bv==0) {
    return(0)
  }
  
  cn2 = 2
  zz = (1-t)*2 + t*(frac*cv + (1-frac)*2 )
  
  
  ccf <- ( x * zz - t*(1-frac)*2*epi - (1-t) *2*epi - t*frac*cv*epi  ) /
    ( t*frac*( bv *(1-epi) - cv*epi ) )
  
  if (constraint) {
    if (is.na(ccf)) {
      return(as.numeric(NA))
    } else if (ccf < 0.9 && bv > 1) {
      return(as.numeric(NA))
    } else if (ccf < upper && ccf > lower) {
      return(ccf)
    } else {
      return(as.numeric(NA))
    }
  } else {
    return(ccf)
  }
}

Assign <- function(x, centers, s) {
  
  n <- length(x)
  
  k <- length(centers)
  
  logRho <- array(0, dim= c(n ,k))
  
  for (ii in 1:k) {
    logRho[,ii] = ccube::bsxfun.se("-", -(x-centers[ii])^2/(2*s[ii]), log(s[ii]))
  }
  
  if (n==k) {
    logR <- ccube::bsxfun.se("-", logRho, ccube:::logsumexp(logRho, 1), expandByRow = F)  # 10.49
  } else {
    logR <- ccube::bsxfun.se("-", logRho, ccube:::logsumexp(logRho, 1)) # 10.49
  }
  
  R <- exp(logR)
  
} 

if (nrow(problemSsm) > 0) {
  
  problemSsm$normal_cn <- 2
  problemSsm$purity <- cellularity
  problemSsm$mutation_id <- problemSsm$gene
  problemSsm <- rename(problemSsm, ref_counts = a, total_counts = d)
  problemSsm <- mutate(rowwise(problemSsm), 
                       var_counts = total_counts - ref_counts, 
                       vaf = var_counts/total_counts,
                       ccf1 = MapVaf2CcfTest(var_counts/total_counts, 
                                             purity, 
                                             major_cn+minor_cn, 
                                             major_cn, cn_frac, constraint = F), 
                       ccf2 = MapVaf2CcfTest(var_counts/total_counts, 
                                             purity, 
                                             major_cn+minor_cn, 
                                             minor_cn, cn_frac, constraint = F),
                       ccf3 = MapVaf2CcfTest(var_counts/total_counts, 
                                             purity, 
                                             major_cn+minor_cn, 
                                             1, cn_frac, constraint = F), 
                       ccube_ccf = mean( unique( c(ccf1, ccf2, ccf3)) ) )



  postR <- Assign(problemSsm$ccube_ccf, res$mu[sort(unique(res$label))], 
                       rep(res$full.model$invWhishartScale, length(unique(res$label))))
  postLabel <- apply(postR, 1, which.max)
  problemSsm$ccube_ccf_mean <- res$mu[sort(uniqLabels)][postLabel]
  problemSsm <- mutate(problemSsm, ccube_mult = mean( unique( c(major_cn, minor_cn, 1) ) ) )
  problemSsm$ccf1 <- NULL
  problemSsm$ccf2 <- NULL
  problemSsm$ccf3 <- NULL
  postR <- data.frame(postR) 
  colnames(postR) <- paste0("cluster_", seq_along(sort(uniqLabels)))
  problemSsm <- cbind(problemSsm, postR)

  problemSsm$label <- postLabel

  tt <- rbind(ssm, problemSsm)

} else {

  tt <- ssm

}

tt <- tt[order(as.numeric(gsub("[^\\d]+", "", tt$id, perl=TRUE))), ]

clusterCertainty <- as.data.frame(table(tt$label), stringsAsFactors = F)
clusterCertainty <- rename(clusterCertainty, cluster = Var1, n_ssms = Freq)
clusterCertainty$proportion <- res$mu[sort(uniqLabels)][as.integer(clusterCertainty$cluster)] * cellularity
clusterCertainty$cluster <- seq_along(uniqLabels) 

if (length(uniqLabels) > 1 & min(ssm$ccube_ccf_mean) < 0.3) {
  clusterCertainty[1, ]$proportion <- 0
} 

write.table(clusterCertainty, file = "1C.txt", sep = "\t", row.names = F, col.names=F, quote = F)

write.table(tt$label, file = "2A.txt", sep = "\t", row.names = F, col.names=F, quote = F)

RR = as.matrix(tt[, paste0("cluster_", seq_along(sort(uniqLabels)))])
coAssign <- Matrix::tcrossprod(RR)
diag(coAssign) <- 1
write.table(coAssign, file = "2B.txt", sep = "\t", row.names = F, col.names=F, quote = F)

# summary graph
fn = "clonal_results_summary.pdf"
ppi <- 500
pdf(fn, width=8, height=8)
par(mfrow=c(2,2))

plot(ssm$ccube_ccf, ssm$vaf, col = myColors[res$label], 
     xlab = "cancer cell fraction", ylab = "variant allele frequecy", 
     main = "ccf vs vaf (colored by cluster memebership)")
ssm$total_cn =ssm$major_cn+ssm$minor_cn
uniqueTotCn = unique(ssm$total_cn)
xx = seq(0,2, length.out = 100)
for (cn in uniqueTotCn) {
  for (i in 1:cn) {
    points(MapVaf2CcfPyClone(xx, cellularity, 2, cn, cn, i, constraint = F), xx, type = 'l')
  }
}

Emu <- res$full.model$ccfMean
Esigma <- res$full.model$ccfCov
Epi <- res$full.model$Epi

params <- data.frame(Emu, Esigma, Epi)
xx <- seq(0,2,  length.out = 1000)
ll <- 0
ll1 <- 0 

for (j in seq_len(nrow(params))) {
  ll <- ll + params[j,]$Epi * dnorm(xx, mean = params[j,]$Emu, sd = sqrt(params[j,]$Esigma))
}

hist(ssm$ccube_ccf, density=20, breaks=20, prob=TRUE, 
     main = "ccf histogram +
       cluster uncertainties",
     xlab = "cancer cell fraction")
lines(xx,ll, lwd=2, col = "darkred")

numSnv <- table(res$label)
names(numSnv) <- as.character(format(round(Emu[sort(uniqLabels)], 2), nsmall = 2))
barplot(numSnv, las = 2, col = myColors[sort(uniqLabels)], 
        xlab = "cluster mean", ylab="number of variants", 
        main = "cluster prevalence")

dev.off()


