library(XML)

#############################
# QTL ARCHIVE SCRAPING
#############################

# download the table
theurl <- "http://qtlarchive.org/db/q?pg=phenolist&reqcat1=body%20weight%20size%20and%20growth&reqcat2=body%20weight"
tables <- readHTMLTable(theurl)
n.rows <- unlist(lapply(tables, function(t) dim(t)[1]))
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
  
write.csv(output, "data_mastersheets/qtlarchive.csv", row.names=FALSE)
