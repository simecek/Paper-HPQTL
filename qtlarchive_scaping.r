library(XML)

#############################
# QTL ARCHIVE SCRAPING
#############################

download.category <- function(url) {
  tables <- readHTMLTable(url)
  col <- as.character(tables[[1]][,1]) # first column
  headers <- grep("QTL data set:", col) # rows with dataset headers
  headers <- c(headers, length(col)+1)
  
  output <- NULL
  for (i in 1:(length(headers)-1)) {
    cross.name <- tolower(sub("QTL data set: ([^ ]*) .*", "\\1",  col[headers[i]]))
    cross.type <- sub("^\\s+","",sub("\\s+$","",strsplit(col[headers[i]], "Â")[[1]][6]))
    founder1 <- sub("^\\s+","",sub("\\s+$","",strsplit(col[headers[i]], "Â")[[1]][3]))
    founder2 <- sub("^\\s+","",sub("\\s+$","",strsplit(col[headers[i]], "Â")[[1]][5]))
    
    for (j in (headers[i]+1):(headers[i+1]-2)) {
      tmp <- sub("^\\s+","",sub("\\s+$","",strsplit(col[j], "Â")[[1]]))
      output <- rbind(output, data.frame(database = "QTL Archive", dataset = cross.name, type = cross.type,
                                         founders = paste(founder1, founder2), trait = tmp[1],
                                         description = tmp[2]))
    }
    
    cross.url <- paste0("http://qtlarchive.org/grpdoc/", cross.name, "/")
    cross.files <- readLines(cross.url)
    csv.line <- grep("csv", cross.files)
    stopifnot(length(csv.line) == 1)
    cross.files[csv.line]
    file.name <- gsub('.*href=\"(.*)\">.*',"\\1", cross.files[csv.line])
    
    download.file(paste0(cross.url, file.name), paste0("data_qtlarchive/",cross.name,".csv"))
    
  }

  return(output) 
}

master.table <- download.category("http://qtlarchive.org/db/q?pg=phenolist&reqcat1=body%20weight%20size%20and%20growth&reqcat2=body%20weight")
master.table <- rbind(master.table,
                      download.category("http://qtlarchive.org/db/q?pg=phenolist&reqcat1=body%20composition&reqcat2=total"))
master.table <- rbind(master.table,
                      download.category("http://qtlarchive.org/db/q?pg=phenolist&reqcat1=body%20weight%20size%20and%20growth&reqcat2=body%20weight%20growth%20curve"))

ignore.list <- read.csv("data_mastersheets/qtlarchive_ignore.csv", as.is=TRUE)
ignore.list$hash <- paste(ignore.list$dataset, ignore.list$trait)

master.table <- subset(master.table, !(paste(dataset, trait) %in% ignore.list$hash))

write.csv(master.table, "data_mastersheets/qtlarchive.csv", row.names=FALSE)

#####################################
# MANUAL CURATION
#####################################

# leitner_2009 - ChrX deleted (all genotypes are N)
# zhang_2012 - ^M deleted
# beamer should be "body_wt", not "bw"


#####################################
# TRY TO LOAD THEM
#####################################

library(qtl)
datalist <- unique(as.character(master.table$dataset))

for (i in 1:length(datalist)) {
  csv.name <- paste0("data_qtlarchive/",datalist[i],".csv")
  tmp <- read.csv(csv.name, na.strings = c("NA", "-", "H", ""))
  (calls <- names(table(as.character(unlist(tmp[-c(1,2,3),which(tmp[1,]!="")])))))

  read.cross(file = csv.name, format = "csv", genotypes = c(calls[1],"H",calls[2]))
}

####################################
# TRY TO PLOT LODSCORE
####################################

master.table <- read.csv("data_mastersheets/qtlarchive.csv", as.is=TRUE)

for (i in 8:nrow(master.table)) {
  csv.name <- paste0("data_qtlarchive/",master.table$dataset[i],".csv")
  
  # guess calls
  tmp <- read.csv(csv.name, na.strings = c("NA", "-", "H", ""))
  calls <- names(table(as.character(unlist(tmp[-c(1,2,3),which(tmp[1,]!="")]))))
  
  cross <- read.cross(file = csv.name, format = "csv", genotypes = c(calls[1],"H",calls[2]))
  cross <- calc.genoprob(cross, step=0)
  n = which(make.names(names(cross$pheno), allow_ = FALSE) == make.names(master.table$trait[i], allow_ = FALSE))
  stopifnot(length(n)==1)
  plot(scanone(cross, pheno.col=n, method="hk"), main= paste(master.table$dataset[i],":",master.table$trait[i]))
}
<<<<<<< HEAD

######################################
# MARKDOWN TABLE
######################################

library(knitr)
master.table$description<- sub("\n", " ", master.table$description)
master.table$dataset <- paste0("[",master.table$dataset,"](http://qtlarchive.org/db/q?pg=projdetails&proj=",master.table$dataset,")")
kable(master.table[,-1], format = "markdown")
=======
>>>>>>> e0546bdeccad5e14b3c50adf2498df9805ece8dc
