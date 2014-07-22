library(HPQTL)

load("../svenson_hap_probs.Rdata")

# load phenotype, exclude mice with missing genotype
pheno <- read.csv("data_do/svenson.bodyweight.csv", as.is=TRUE)
idx <- match(dimnames(probs)[[1]], pheno$Sample)
pheno <- pheno[idx,]
stopifnot(pheno$Sample == dimnames(probs)[[1]])

# snpmap <- read.csv("H:\\SmallProjects\\1309 svensonDO\\kindship pred\\SNP_Map.txt", sep="\t", as.is=TRUE)
# idx <- match(dimnames(probs)[[3]], snpmap$Name)
# snpmap <- snpmap[idx,]
# stopifnot(snpmap$Name == dimnames(probs)[[3]])

# markers <- data.frame(marker = snpmap$Name, chr = snpmap$Chromosome, pos = snpmap$Position)
# saveRDS(markers, "data_do/snpmap.rds")
markers <- readRDS("data_do/snpmap.rds")
attr(probs, "markers") <- markers

# making 'genotype.probs' object
geno <- extract.geno(probs)
# G <- gensim.matrix(geno, method='allele-multif-cosine')
# Glist <- gensim.matrix(geno, method='allele-multif-cosine', procedure="LMM-L1O")
# save(G, Glist, file="../kinship.rdata")
load("../kinship.rdata")

# heritability calculations
obs <- which(!is.na(pheno$Weight1) & !is.na(pheno$Weight2))
h2 <- function(x) x$sigma[1] / sum(x$sigma)

fit1 <- regress(Weight1 ~ Sex*Diet + Gen, data=pheno[obs,],~G[obs,obs])
fit2 <- regress(Weight2 ~ Sex*Diet + Gen, data=pheno[obs,],~G[obs,obs])

n <- which(names(pheno) == "Weight2")
covar <- model.matrix(~ Sex*Diet + Gen, data=pheno)[,-1]

# scan1
system.time(fit.lm <- scan1(geno, pheno, covar=covar, pheno.cols=n, procedure="LM"))
system.time(fit.lmm <- scan1(geno, pheno, covar=covar, pheno.cols=n, G=G, procedure="LMM"))
system.time(fit.l1o <- scan1(geno, pheno, covar=covar, pheno.cols=n, G=Glist, procedure="LMM-L1O"))

# save results
saveRDS(fit.lm, file="results_do/do_Weight2_lm.rds")
saveRDS(fit.lmm, file="results_do/do_Weight2_lmm.rds")
saveRDS(fit.l1o, file="results_do/do_Weight2_l1o.rds")

# system.time(trhold.lm <- scan1.threshold(geno, pheno, covar=covar, pheno.cols=n, procedure="LM", n.perm=3, keep.lods=TRUE))
load("results_do/do_Weight2_thresholds.rdata")

plot(fit.lm, main="LM", incl.markers=FALSE)
abline(h=trhold.lm, col="red", lty=2)
plot(fit.lmm, main="LMM", incl.markers=FALSE)
abline(h=trhold.lmm, col="red", lty=2)
plot(fit.l1o, main="LMM-L1O", incl.markers=FALSE, bandcol="gray70")
abline(h=trhold.l1o, col="red", lty=2)

# regressing out
most.significant.snp <- probs[,,3951]
dim(most.significant.snp)
obs <- which(!is.na(pheno$Weight2))
regress(Weight2 ~ Sex*Diet + Gen + most.significant.snp[obs,-1], data=pheno[obs,],~G[obs,obs])
covar2 <- model.matrix(~ Sex*Diet + Gen + most.significant.snp[,-1], data=pheno)[,-1]

system.time(fit.lm2 <- scan1(geno, pheno, covar=covar2, pheno.cols=n, procedure="LM"))
system.time(fit.l1o2 <- scan1(geno, pheno, covar=covar2, pheno.cols=n, G=Glist, procedure="LMM-L1O"))
saveRDS(fit.lm2, file="results_do/do_Weight2_lm2.rds")
saveRDS(fit.l1o2, file="results_do/do_Weight2_l1o2.rds")