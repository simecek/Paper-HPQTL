library(HPQTL)

set.seed(123)

master.table <- read.csv("data_qtlarchive/_mastersheet_qtlarchive_heritability.csv", as.is=TRUE)

### 'scanone' objects for all 3 methods

for (i in 1:nrow(master.table)) {
  print(i)
  csv.name <- paste0("data_qtlarchive/",master.table$dataset[i],".csv")
  
  # guess calls
  tmp <- read.csv(csv.name, na.strings = c("NA", "-", "H", ""))
  calls <- names(table(as.character(unlist(tmp[-c(1,2,3),which(tmp[1,]!="")]))))
  
  cross <- read.cross(file = csv.name, format = "csv", genotypes = c(calls[1],"H",calls[2]))
  cross <- calc.genoprob(cross, step=0)
  
  n = which(make.names(names(cross$pheno), allow_ = FALSE) == make.names(master.table$trait[i], allow_ = FALSE))
  stopifnot(length(n)==1)
  sex.column <- grep("sex", tolower(names(cross$pheno)))
  stopifnot(length(sex.column)<=1)
  
  geno <- extract.geno(cross)
  G <- gensim.matrix(geno, method="allele-2f-additive")
  Glist <- gensim.matrix(geno, method="allele-2f-additive", procedure="LMM-L1O")
  
  if (length(sex.column) == 0 || length(table(cross$pheno[, sex.column]))==1) {
    fit.lm <- scan1(geno, cross$pheno, pheno.cols=n, procedure="LM")
    fit.lmm <- scan1(geno, cross$pheno, pheno.cols=n, procedure="LMM", G=G)
    fit.l1o <- scan1(geno, cross$pheno, pheno.cols=n, procedure="LMM-L1O", G=Glist)
  } else {
    sex <- as.matrix(cross$pheno[, sex.column, drop=FALSE])
    fit.lm <- scan1(geno, cross$pheno, covar=sex, pheno.cols=n, procedure="LM")
    fit.lmm <- scan1(geno, cross$pheno, covar=sex, pheno.cols=n, procedure="LMM", G=G)
    fit.l1o <- scan1(geno, cross$pheno, covar=sex, pheno.cols=n, procedure="LMM-L1O", G=Glist)
  }
  
  saveRDS(fit.lm, file=paste0("results_qtlarchive/",master.table$dataset[i],"_lm.rds"))
  saveRDS(fit.lmm, file=paste0("results_qtlarchive/",master.table$dataset[i],"_lmm.rds"))
  saveRDS(fit.l1o, file=paste0("results_qtlarchive/",master.table$dataset[i],"_l1o.rds"))
}

### LM thresholds calculation (alpha=5%)

for (i in 1:nrow(master.table)) {
  print(i)
  csv.name <- paste0("data_qtlarchive/",master.table$dataset[i],".csv")
  
  # guess calls
  tmp <- read.csv(csv.name, na.strings = c("NA", "-", "H", ""))
  calls <- names(table(as.character(unlist(tmp[-c(1,2,3),which(tmp[1,]!="")]))))
  
  cross <- read.cross(file = csv.name, format = "csv", genotypes = c(calls[1],"H",calls[2]))
  cross <- calc.genoprob(cross, step=0)
  
  n = which(make.names(names(cross$pheno), allow_ = FALSE) == make.names(master.table$trait[i], allow_ = FALSE))
  stopifnot(length(n)==1)
  sex.column <- grep("sex", tolower(names(cross$pheno)))
  stopifnot(length(sex.column)<=1)
  
  geno <- extract.geno(cross)
  
  if (length(sex.column) == 0 || length(table(cross$pheno[, sex.column]))==1) {
    trhold.lm <- scan1.threshold(geno, cross$pheno, pheno.cols=n, procedure="LM", n.perm=1000)
  } else {
    sex <- as.matrix(cross$pheno[, sex.column, drop=FALSE])
    trhold.lm <- scan1.threshold(geno, cross$pheno, covar=sex, pheno.cols=n, procedure="LM", n.perm=1000)
  }
  
  saveRDS(trhold.lm, file=paste0("results_qtlarchive/",master.table$dataset[i],"_lm_trhold.rds"))
}

### LMM thresholds calculation (alpha=5%)

for (i in 1:nrow(master.table)) {
  print(i)
  csv.name <- paste0("data_qtlarchive/",master.table$dataset[i],".csv")
  
  # guess calls
  tmp <- read.csv(csv.name, na.strings = c("NA", "-", "H", ""))
  calls <- names(table(as.character(unlist(tmp[-c(1,2,3),which(tmp[1,]!="")]))))
  
  cross <- read.cross(file = csv.name, format = "csv", genotypes = c(calls[1],"H",calls[2]))
  cross <- calc.genoprob(cross, step=0)
  
  n = which(make.names(names(cross$pheno), allow_ = FALSE) == make.names(master.table$trait[i], allow_ = FALSE))
  stopifnot(length(n)==1)
  sex.column <- grep("sex", tolower(names(cross$pheno)))
  stopifnot(length(sex.column)<=1)
  
  geno <- extract.geno(cross)
  G <- gensim.matrix(geno, method="allele-2f-additive")
  
  if (length(sex.column) == 0 || length(table(cross$pheno[, sex.column]))==1) {
    trhold.lmm <- try(scan1.threshold(geno, cross$pheno, pheno.cols=n, procedure="LMM", G=G, n.perm=1000), TRUE)
  } else {
    sex <- as.matrix(cross$pheno[, sex.column, drop=FALSE])
    trhold.lmm <- try(scan1.threshold(geno, cross$pheno, covar=sex, pheno.cols=n, procedure="LMM", G=G, n.perm=1000), TRUE)
  }
  
  if (class(trhold.lmm)!="try-error") saveRDS(trhold.lmm, file=paste0("results_qtlarchive/",master.table$dataset[i],"_lmm_trhold.rds"))
}

### LMM-L1O thresholds calculation (alpha=5%)

for (i in 1:nrow(master.table)) {
  set.seed(123)
  print(i)
  csv.name <- paste0("data_qtlarchive/",master.table$dataset[i],".csv")
  
  # guess calls
  tmp <- read.csv(csv.name, na.strings = c("NA", "-", "H", ""))
  calls <- names(table(as.character(unlist(tmp[-c(1,2,3),which(tmp[1,]!="")]))))
  
  cross <- read.cross(file = csv.name, format = "csv", genotypes = c(calls[1],"H",calls[2]))
  cross <- calc.genoprob(cross, step=0)
  
  n = which(make.names(names(cross$pheno), allow_ = FALSE) == make.names(master.table$trait[i], allow_ = FALSE))
  stopifnot(length(n)==1)
  sex.column <- grep("sex", tolower(names(cross$pheno)))
  stopifnot(length(sex.column)<=1)
  
  geno <- extract.geno(cross)
  Glist <- gensim.matrix(geno, method="allele-2f-additive", procedure="LMM-L1O")
  
  if (length(sex.column) == 0 || length(table(cross$pheno[, sex.column]))==1) {
    trhold.l1o <- try(scan1.threshold(geno, cross$pheno, pheno.cols=n, procedure="LMM-L1O", G=Glist, n.perm=1000), TRUE)
  } else {
    sex <- as.matrix(cross$pheno[, sex.column, drop=FALSE])
    trhold.l1o <- try(scan1.threshold(geno, cross$pheno, covar=sex, pheno.cols=n, procedure="LMM-L1O", G=Glist, n.perm=1000), TRUE)
  }
  
  if (class(trhold.l1o)!="try-error") saveRDS(trhold.l1o, file=paste0("results_qtlarchive/",master.table$dataset[i],"_l1o_trhold.rds"))
}

### SUMMARY