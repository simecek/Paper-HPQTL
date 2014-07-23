# read LXS MDA haplotype reconstruction, k

lxsmap <- read.csv("C:\\Users\\petrs\\Dropbox\\High Resolution Mapping of Reference Populations\\High Resolution Mapping of Reference Populations\\web\\consensus genetic maps - CSVr files\\geno.lxs.csv", as.is=TRUE)
signature <- function(x) paste0(x[-(1:4)],collapse="")
lxsmap$sign <- apply(lxsmap,1,signature)
breakpoints <- which(lxsmap$sign[-nrow(lxsmap)] != lxsmap$sign[-1])
lxsmap <- lxsmap[breakpoints,]
lxsmap$positionCM <- as.numeric(lxsmap$positionCM)

#########
# Bennett
#########

pheno <- read.csv("data_ris/Bennett1_table.csv", skip=5, as.is=TRUE)
pheno$strain <- sub("/TejJ", "", pheno$strain)
pheno$strain <- sub("ILSXISS", "LXS", pheno$strain)
pheno <- pheno[pheno$strain %in% names(lxsmap),]
idx <- match(pheno$strain, names(lxsmap))
lxs <- as.data.frame(t(lxsmap[,c(2,4,idx)]))
names(lxs) <- lxsmap[,1]
output <- cbind(rbind(NA,NA,pheno), lxs)

write.csv(output, "data_ris/Bennett1.csv", row.names = FALSE, na="")


# real processing

cross <- read.cross(file="data_ris/Bennett1.csv", format="csv", genotypes = c("S", "H", "L"), na.strings=c("N",""), 
           crosstype="risib")
cross <- calc.genoprob(cross, step=0)
for (i in 1:ncol(cross$pheno)) cross$pheno[,i] <- as.numeric(cross$pheno[,i])
class(cross) <- c("f2", "cross")

geno <- extract.geno(cross)
G <- gensim.matrix(geno, method="allele-2f-additive")

n=9 # maxi
Glist <- gensim.matrix(geno, method="allele-2f-additive", procedure="LMM-L1O")

fit.lm <- scan1(geno, cross$pheno, pheno.cols=n, procedure="LM")
fit.lmm <- scan1(geno, cross$pheno, pheno.cols=n, procedure="LMM", G=G)
fit.l1o <- scan1(geno, cross$pheno, pheno.cols=n, procedure="LMM-L1O", G=Glist)
saveRDS(fit.lm, file="results_ris/Bennett1_lm.rds")
saveRDS(fit.lmm, file="results_ris/Bennett1_lmm.rds")
saveRDS(fit.l1o, file="results_ris/Bennett1_l1o.rds")

system.time(trhold.lm <- scan1.threshold(geno, cross$pheno, pheno.cols=n, procedure="LM", n.perm=1000))
system.time(trhold.lmm <- scan1.threshold(geno, cross$pheno, pheno.cols=n, procedure="LMM", G=G, n.perm=100))
system.time(trhold.l1o <- scan1.threshold(geno, cross$pheno, pheno.cols=n, procedure="LM", G=Glist,n.perm=100))

saveRDS(trhold.lm, file="results_ris/Bennett1_lm_trhold.rds")
saveRDS(trhold.lmm, file="results_ris/Bennett1_lmm_trhold.rds")
saveRDS(trhold.l1o, file="results_ris/Bennett1_l1o_trhold.rds")

### Nelson1_table.csv

pheno <- read.csv("data_ris/Nelson1_table.csv", skip=5, as.is=TRUE)
pheno$strain <- sub("/TejJ", "", pheno$strain)
pheno$strain <- sub("ILSXISS", "LXS", pheno$strain)
pheno <- pheno[pheno$strain %in% names(lxsmap),]
idx <- match(pheno$strain, names(lxsmap))
lxs <- as.data.frame(t(lxsmap[,c(2,4,idx)]))
names(lxs) <- lxsmap[,1]
output <- cbind(rbind(NA,NA,pheno), lxs)

write.csv(output, "data_ris/Nelson1.csv", row.names = FALSE, na="")


# real processing

cross <- read.cross(file="data_ris/Nelson1.csv", format="csv", genotypes = c("S", "H", "L"), na.strings=c("N",""), 
                    crosstype="f2")
cross <- calc.genoprob(cross, step=0)
for (i in 1:ncol(cross$pheno)) cross$pheno[,i] <- as.numeric(cross$pheno[,i])
class(cross) <- c("f2", "cross")

cross[,!is.na(cross$pheno[,5])]

geno <- extract.geno(cross)
G <- gensim.matrix(geno, method="allele-2f-additive")

# obs <- which(!is.na(cross$pheno$bw_AL))
# h2(regress(bw_DR~sex*diet, ~G[obs,obs], data=cross$pheno[obs,]))

n=5 # bw_AL
Glist <- gensim.matrix(geno, method="allele-2f-additive", procedure="LMM-L1O")
covar <- model.matrix(~sex, data=cross$pheno)[,-1]

fit.lm <- scan1(geno, cross$pheno, pheno.cols=n, covar=covar, procedure="LM")
fit.lmm <- scan1(geno, cross$pheno, pheno.cols=n, covar=covar, procedure="LMM", G=G)
fit.l1o <- scan1(geno, cross$pheno, pheno.cols=n, covar=covar, procedure="LMM-L1O", G=Glist)
saveRDS(fit.lm, file="results_ris/Nelson1_lm.rds")
saveRDS(fit.lmm, file="results_ris/Nelson1_lmm.rds")
saveRDS(fit.l1o, file="results_ris/Nelson1_l1o.rds")

system.time(trhold.lm <- scan1.threshold(geno, cross$pheno, pheno.cols=n, procedure="LM", n.perm=1000))
system.time(trhold.lmm <- scan1.threshold(geno, cross$pheno, pheno.cols=n, procedure="LMM", G=G, n.perm=100))
system.time(trhold.l1o <- scan1.threshold(geno, cross$pheno, pheno.cols=n, procedure="LM", G=Glist,n.perm=100))

saveRDS(trhold.lm, file="results_ris/Bennett1_lm_trhold.rds")
saveRDS(trhold.lmm, file="results_ris/Bennett1_lmm_trhold.rds")
saveRDS(trhold.l1o, file="results_ris/Bennett1_l1o_trhold.rds")

# read BXD MDA haplotype reconstruction

bxdmap <- read.csv("C:\\Users\\petrs\\Dropbox\\High Resolution Mapping of Reference Populations\\High Resolution Mapping of Reference Populations\\web\\consensus genetic maps - CSVr files\\geno.bxd.csv", as.is=TRUE)
bxdmap <- bxdmap[!is.na(bxdmap[,4]),]
signature <- function(x) paste0(x[-(1:4)],collapse="")
bxdmap$sign <- apply(bxdmap,1,signature)
breakpoints <- which(bxdmap$sign[-nrow(bxdmap)] != bxdmap$sign[-1])
bxdmap <- bxdmap[breakpoints,]
bxdmap$positionCM <- as.numeric(bxdmap$positionCM)

pheno <- read.csv("data_ris/Schughart2_table.csv", skip=5, as.is=TRUE)
pheno$strain <- gsub("(BXD[0-9]*)/.*", "\\1", pheno$strain)
names(bxdmap) <- gsub("(BXD[0-9]*)[^0-9]*.*", "\\1", names(bxdmap))
names(bxdmap) <- gsub("[.]", "/", names(bxdmap))
pheno <- pheno[pheno$strain %in% names(bxdmap),]
idx <- match(pheno$strain, names(bxdmap))
bxd <- as.data.frame(t(bxdmap[,c(2,4,idx)]))
names(bxd) <- bxdmap[,1]
output <- cbind(rbind(NA,NA,pheno), bxd)

write.csv(output, "data_ris/Schughart2.csv", row.names = FALSE, na="")


# real processing

cross <- read.cross(file="data_ris/Schughart2.csv", format="csv", genotypes = c("B", "H", "D"), na.strings=c("N",""), 
                    crosstype="risib")
cross <- calc.genoprob(cross, step=0)
for (i in 1:ncol(cross$pheno)) cross$pheno[,i] <- as.numeric(cross$pheno[,i])
class(cross) <- c("f2", "cross")

geno <- extract.geno(cross)
G <- gensim.matrix(geno, method="allele-2f-additive")

n=5 # maxi
Glist <- gensim.matrix(geno, method="allele-2f-additive", procedure="LMM-L1O")
covar <- model.matrix(~sex, data=cross$pheno)[,-1]


fit.lm <- scan1(geno, cross$pheno, pheno.cols=n, covar=covar, procedure="LM")
fit.lmm <- scan1(geno, cross$pheno, pheno.cols=n, covar=covar,procedure="LMM", G=G)
fit.l1o <- scan1(geno, cross$pheno, pheno.cols=n, covar=covar,procedure="LMM-L1O", G=Glist)
saveRDS(fit.lm, file="results_ris/Schughart2_lm.rds")
saveRDS(fit.lmm, file="results_ris/Schughart2_lmm.rds")
saveRDS(fit.l1o, file="results_ris/Schughart2_l1o.rds")

system.time(trhold.lm <- scan1.threshold(geno, cross$pheno, pheno.cols=n, covar=covar,procedure="LM", n.perm=1000))
system.time(trhold.lmm <- scan1.threshold(geno, cross$pheno, pheno.cols=n, covar=covar,procedure="LMM", G=G, n.perm=100))
system.time(trhold.l1o <- scan1.threshold(geno, cross$pheno, pheno.cols=n, covar=covar,procedure="LM", G=Glist,n.perm=100))

saveRDS(trhold.lm, file="results_ris/Schughart2_lm_trhold.rds")
saveRDS(trhold.lmm, file="results_ris/Schughart2_lmm_trhold.rds")
saveRDS(trhold.l1o, file="results_ris/Schughart2_l1o_trhold.rds")
