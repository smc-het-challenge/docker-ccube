#!/Library/Frameworks/R.framework/Resources/Rscript

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

#dir.create(sampleName, recursive = T)

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


library(doParallel)
library(foreach)

registerDoParallel(cores = detectCores())

kk <- 6
if (kk > nrow(ssm)){
  kk <- nrow(ssm) - 1
}
rr <- 1
iterSetting <- data.frame( sort(rep(seq(1, kk, length.out = kk), rr)) )

results <- foreach(n = 1:nrow(iterSetting), 
  .combine = c, .packages = "ccube") %dopar% {
  library(ccube)
  k <- iterSetting[n, ]
  list(ccube_m6(ssm, epi=1e-3,
         init=k, tol = 1e-10, maxiter = 1e3,
         fit_mult = T, fit_hyper = T, use = "use_base", verbose = F))
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
                          constraint=F) )

uniqLabels <- unique(res$label)
write.table(length(uniqLabels), file = "1B.txt", sep = "\t", row.names = F, col.names=F, quote = F)


clusterCertainty <- as.data.frame(table(res$label), stringsAsFactors = F)
clusterCertainty <- rename(clusterCertainty, cluster = Var1, n_ssms = Freq)
clusterCertainty$proportion <- res$mu[as.integer(clusterCertainty$cluster)] * cellularity
clusterCertainty$cluster <- seq_along(uniqLabels) 

write.table(clusterCertainty, file = "1C.txt", sep = "\t", row.names = F, col.names=F, quote = F)



if (length(uniqLabels) == 1) {
  mutR = data.frame(res$R)
  colnames(mutR) <- 1
} else {
  mutR <- data.frame(res$R[, sort(uniqLabels)]) 
  colnames(mutR) <- seq_along(uniqLabels)
}

label <- apply(mutR, 1, which.max)
write.table(label, file = "2A.txt", sep = "\t", row.names = F, col.names=F, quote = F)


coAssign <- Matrix::tcrossprod(res$R[, sort(uniqLabels)])
diag(coAssign) <- 1
write.table(coAssign, file = "2B.txt", sep = "\t", row.names = F, col.names=F, quote = F)

# summary graph
fn = "results_summary.pdf"
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


