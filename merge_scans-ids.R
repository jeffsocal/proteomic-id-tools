######################################################################
# AUTH  Jeff Jones
# DATE  2018.07.23
# DESC  parse mzXML post TPP: msconvert
#
######################################################################

rm(list=ls())
require(XML)

path_ms2 <- NULL
path_pep <- NULL

text_cmd <- "R --vanilla merge_scans-ids.R <path_to.mzXML.csv> <path_to.pepXML.csv>"

args = commandArgs(trailingOnly=TRUE)
#
# test if there are 2 arguments
#
if (length(args) !=2 ) {

  stop(text_cmd, call.=FALSE)

} else {

  path_ms2 <- args[1]
  path_pep <- args[2]

  if(!grepl("mzXML.csv$", path_ms2, ignore.case = T) | !grepl("pepXML.csv$", path_pep, ignore.case = T))

    stop(text_cmd, call.=FALSE)

}


#
# Read in the data
#
data_ms2 <- read.csv(path_ms2)
data_pep <- read.csv(path_pep)

#
# use merge the two datasets, assuming no scans have been co-added
# keep all scans
#
data <- merge(data_ms2, data_pep, 
              by.x='scan', 
              by.y='start_scan',
              all.x = TRUE,
              suffixes=c("_scan", "_id"))

path_merged <- sub('.csv', '.pepXML.csv', path_ms2)

write.csv(data, path_merged, row.names = FALSE)
