######################################################################
# AUTH  Jeff Jones
# DATE  2018.07.23
# DESC  parse protXML post TPP: ProteinProphet
#
######################################################################

rm(list=ls())
require(XML)

path_xml <- NULL
path_csv <- NULL

text_cmd <- "R --vanilla protxml2csv.R <path_to.protXML> <path_to.csv>"

args = commandArgs(trailingOnly=TRUE)
#
# test if there are 2 arguments
#
if (length(args) !=2 ) {
  
  stop(text_cmd, call.=FALSE)
  
} else {
  
  path_xml <- args[1]
  path_csv <- args[2]
  
  if(!grepl("protXML$", path_xml, ignore.case = T) | !grepl("csv$", path_csv, ignore.case = T))
    
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
hit_i <- 0
#
# traverse the master list
#
for ( prot_i in which(names(xml_data) == 'protein_group') ){
  
  prot <- xml_data[prot_i]
  
  #
  # grab the protein attributes
  #
  prot_att <- append(unlist(prot$protein_group$.attrs), file_xml)
  names(prot_att)[length(prot_att)] <- 'fileName'
  
  w_proteins <- which(names(prot$protein_group) == 'protein')
  
  for(result in prot$protein_group[w_proteins]){
    
    #
    # grab the result attributes
    # 
    result_att <- unlist(result$.attrs)
    
    w_peptides <- which(names(result) == 'peptide')
    
    for(peptide in result[w_peptides]){
      
      peptide_vals <- unlist(peptide$.attrs)
      
      #
      # append all the values to our new list
      #
      
      hit_i <- hit_i + 1
      
      hits[[hit_i]] <- list()
      hits[[hit_i]] <- append(hits[[hit_i]], as.list(prot_att))
      hits[[hit_i]] <- append(hits[[hit_i]], as.list(result_att))
      hits[[hit_i]] <- append(hits[[hit_i]], as.list(peptide_vals))
    }
    
    
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
