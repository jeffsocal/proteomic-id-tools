######################################################################
# AUTH  Jeff Jones
# DATE  2018.07.23
# DESC  parse mzXML post TPP: msconvert
#
######################################################################

rm(list=ls())
require(XML)
library(progress)

options(warn=-1)

help_text <- "
 NAME
    mzxml2csv.R

 SYNOPSIS
    mzxml2csv.R --xml=<path_mzXML> --csv=<path_csv>

 DESCRIPTION
    extract MSn meta data from mzXML files for downstream analysis

 COMMAND LINE

    --xml <path_to.mzXML>

    --csv <path_to.csv> (optional)

 EXAMPLE

    Rscript mzxml2csv.R --xml=<path_mzXML>

"

###############################################################################
# USER INPUT
path_xml                      <- NULL
path_csv                      <- NULL

for (arg in commandArgs()){
    arg_value <- as.character(sub("-+[a-z]*\\=", "", arg))
    if( grepl("--xml", arg) ) path_xml <- arg_value
    if( grepl("--csv", arg) ) path_csv <- arg_value
    if( grepl("--help", arg) ) stop(help_text)
}

###############################################################################
# INPUT VALIDATION
message <- NULL
if(is.null(path_xml)) message <- stop("ERROR\n", "  no mzXML file declared\n")
if(!grepl("mzXML$", path_xml)) message <- paste0(message, "  mz file (--xml) not a supported format\n")
if(is.null(path_csv))
    path_csv = paste0(path_xml, ".csv")
if(!grepl(".csv$", path_csv)) message <- paste0(message, "  csv file (--csv) not a supported format\n")

if(!is.null(message)) stop("ERROR\n", message)


cat("mzXML to CSV started\n")
cat(" xml file:                        ", path_xml, "\n")
cat(" csv file:                        ", path_csv, "\n")


#
# Read in the data
#
cat(" reading xml file ...")
data <- xmlParse(path_xml)
file_xml <- sub('.*/', '', path_xml)

#
# convert to a master list
#
xml_data <- xmlToList(data)

w_scans <- which(names(xml_data$msRun) == 'scan')
cat("\n")

scans <- list()
scan_i <- 0

pb <- progress_bar$new(
    format = " extracting data [:bar] :percent eta: :eta",
    total = length(w_scans), clear = FALSE, width= 60)
#
# traverse the master list
#
for ( scan_n in w_scans ){
    
    pb$tick()
    
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

cat("\nFINISHED:\n")
cat(" data written to                  ", path_csv, "\n")
