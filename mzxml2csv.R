######################################################################
# AUTH  Jeff Jones
# DATE  2018.07.23
# DESC  parse mzXML post TPP: msconvert
#
######################################################################

rm(list=ls())
require(XML)

path_xml <- NULL
path_csv <- NULL

text_cmd <- "R --vanilla mzxml2csv.R <path_to.mzXML> <path_to.csv>"

args = commandArgs(trailingOnly=TRUE)
#
# test if there are 2 arguments
#
if (length(args) !=2 ) {

  stop(text_cmd, call.=FALSE)

} else {

  path_xml <- args[1]
  path_csv <- args[2]

  if(!grepl("mzXML$", path_xml, ignore.case = T) | !grepl("csv$", path_csv, ignore.case = T))

    stop(text_cmd, call.=FALSE)

}


#
# Read in the data
#
data <- xmlParse(path_xml)
file_xml <- sub('.*/', '', path_xml)

#
# convert to a master list
#
xml_data <- xmlToList(data)

scans <- list()
scan_i <- 0
#
# traverse the master list
#
for ( scan_n in which(names(xml_data$msRun) == 'scan') ){
  
  scan <- xml_data$msRun[scan_n]
  
  #
  # grab the spectrum attributes
  #
  scan_att <- append(unlist(scan$scan$.attrs), file_xml)
  names(scan_att)[length(scan_att)] <- 'fileName'
  
  scan_pre_att <- unlist(scan$scan$precursorMz$.attrs)
  scan_pre_att <- append(scan_pre_att, scan$scan$precursorMz$text)
  names(scan_pre_att)[length(scan_pre_att)] <- 'precursorMz'
  
  scan_i <- scan_i + 1
  
  scans[[scan_i]] <- list()
  scans[[scan_i]] <- append(scans[[scan_i]], as.list(scan_att))
  scans[[scan_i]] <- append(scans[[scan_i]], as.list(scan_pre_att))
  
}

#
# use lapply to converge the targets list into a matrix with named columns
#
data <- do.call(rbind, lapply(lapply(scans, unlist), "[",
                              unique(unlist(c(sapply(scans,names))))))

#
# convert the matrix to a data frame
#
data <- as.data.frame(data)
colnames(data) <- unique(unlist(c(sapply(scans,names))))

w_num <- which(colnames(data) == 'num')
colnames(data)[w_num] <- 'scan'
write.csv(data, path_csv, row.names = FALSE)
