######################################################################
# AUTH  Jeff Jones
# DATE  2018.07.23
# DESC  parse pepXML post TPP: xtandem > Tandem2XML > PeptideProphet
#
######################################################################

rm(list=ls())
require(XML)

path_xml <- NULL
path_csv <- NULL

text_cmd <- "R --vanilla pepxml2csv.R <path_to.pepXML> <path_to.csv>"

args = commandArgs(trailingOnly=TRUE)
#
# test if there are 2 arguments
#
if (length(args) !=2 ) {
  
  stop(text_cmd, call.=FALSE)
  
} else {
  
  path_xml <- args[1]
  path_csv <- args[2]
  
  if(!grepl("pepXML$", path_xml, ignore.case = T) | !grepl("csv$", path_csv, ignore.case = T))
    
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

hits <- list()

#
# traverse the master list
#
for ( spec_i in which(names(xml_data$msms_run_summary) == 'spectrum_query') ){
  
  hits[[spec_i]] <- list()
  
  spec <- xml_data$msms_run_summary[spec_i]
  
  #
  # grab the spectrum attributes
  #
  spec_att <- append(unlist(spec$spectrum_query$.attrs), file_xml)
  names(spec_att)[length(spec_att)] <- 'fileName'
  
  for(result in spec$spectrum_query$search_result){
    
    #
    # grab the result attributes
    # 
    result_att <- unlist(result$.attrs)
    
    #
    # grab the X!T scoring values 
    #    
    w_search_scores <- which(names(result) == 'search_score')
    score_name <- as.character(as.data.frame(t(as.data.frame(result[w_search_scores])))$name)
    score_vals <- as.character(as.data.frame(t(as.data.frame(result[w_search_scores])))$value)
    names(score_vals) <- score_name  
    
    #
    # grab the Prophet attributes
    #
    prophet_att <- unlist(result$analysis_result$peptideprophet_result$.attrs)
    
    #
    # grab the Prophet scores
    #    
    prophet_scores <- result$analysis_result$peptideprophet_result$search_score_summary
    w_prophet_param <- which(names(result) == 'parameter')
    prophet_name <- as.character(as.data.frame(t(as.data.frame(prophet_scores[w_prophet_param])))$name)
    prophet_vals <- as.character(as.data.frame(t(as.data.frame(prophet_scores[w_prophet_param])))$value)
    names(prophet_vals) <- prophet_name  
    
    #
    # append all the values to our new list
    #
    hits[[spec_i]] <- append(hits[[spec_i]], as.list(spec_att))
    hits[[spec_i]] <- append(hits[[spec_i]], as.list(result_att))
    hits[[spec_i]] <- append(hits[[spec_i]], as.list(score_vals))
    hits[[spec_i]] <- append(hits[[spec_i]], as.list(prophet_att))
    hits[[spec_i]] <- append(hits[[spec_i]], as.list(prophet_vals))
  }
}

#
# use lapply to converge the hits list into a matrix with named columns
#
data <- do.call(rbind, lapply(lapply(hits, unlist), "[",
                              unique(unlist(c(sapply(hits,names))))))

#
# convert the matrix to a data frame
#
data <- as.data.frame(data)
colnames(data) <- unique(unlist(c(sapply(hits,names))))

write.csv(data, path_csv, row.names = FALSE)
