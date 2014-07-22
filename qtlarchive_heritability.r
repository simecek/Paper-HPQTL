library(HPQTL)

h2 <- function(x) x$sigma[1] / sum(x$sigma)

master.table <- read.csv("data_qtlarchive/_mastersheet_qtlarchive.csv", as.is=TRUE)
master.table$h2 <- NA
master.table$totmar <- NA
master.table$nind <- NA
master.table$observed <- NA


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
  obs <- which(!is.na(cross$pheno[,n]))
  y <- cross$pheno[obs,n]
  
  geno <- extract.geno(cross)
  G <- gensim.matrix(geno, method="allele-2f-additive")[obs,obs]
  
  sex.column <- grep("sex", tolower(names(cross$pheno)))
  
  if (length(sex.column) == 0 || length(table(cross$pheno[, sex.column]))==1) {
    master.table$h2[i] <- h2(regress(y~1, ~G))
  } else {
    selsex <- cross$pheno[obs, sex.column]
    master.table$h2[i] <- h2(regress(y~selsex, ~G))
  }  
  
  master.table$totmar[i] <- totmar(cross)
  master.table$nind[i] <- nind(cross)
  master.table$observed[i] <- length(obs)/nind(cross)
}

# maximum 20% of missing observations
master.table <- subset(master.table, observed >= 0.8)

# for each dataset, take trait with maximum heritability
master.table <- master.table[order(master.table$dataset, -master.table$h2),]
master.table <- subset(master.table, !duplicated(dataset))
write.csv(master.table, "data_qtlarchive/_mastersheet_qtlarchive_heritability.csv", row.names=FALSE)
