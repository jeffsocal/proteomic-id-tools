######################################################################
# AUTH  Jeff Jones
# DATE  2018.07.23
# DESC  parse protXML post TPP: ProteinProphet
#
######################################################################

rm(list=ls())
require(XML)
library(progress)

path_xml <- NULL
path_csv <- NULL

text_cmd <- "
R --vanilla --slave < protxml2csv.R --args <path_to.protXML> <path_to.csv> <FDR_cutoff> <min_n_peptides>
    FDR_cutoff:        value between 0-1 (default: 0.1 ~ 10% FDR) 
    min_n_peptides:    value >= 1 (default: 2) assumes user avoids 1-hit proteins
"

args = commandArgs(trailingOnly=TRUE)
#
# test if there are 2 arguments
#
if (length(args) !=4 ) {
    
    stop(text_cmd, call.=FALSE)
    
} else {
    
    path_xml <- args[1]
    path_csv <- args[2]
    fdr_cutf <- args[3]
    min_npep <- args[4]
    
    if(!grepl("protXML$", path_xml, ignore.case = T) | !grepl("csv$", path_csv, ignore.case = T))
        
        stop(text_cmd, call.=FALSE)
    
}


cat("pepXML to CSV started\n")
cat(" xml file:                        ", path_xml, "\n")
cat(" csv file:                        ", path_csv, "\n")
cat(" fdr cut off:                     ", fdr_cutf, "\n")
cat(" min n peptides:                  ", min_npep, "\n")

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

w_protg <- which(names(xml_data) == 'protein_group')
cat("\n")

pb <- progress_bar$new(
    format = " extracting data [:bar] :percent eta: :eta",
    total = length(w_protg), clear = FALSE, width= 60)
hits <- list()
hit_i <- 0
#
# traverse the master list
#
for ( prot_i in w_protg ){
    
    pb$tick()
    
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
data$probability <- as.numeric(as.character(data$probability))
data$total_number_distinct_peptides <- as.numeric(as.character(data$total_number_distinct_peptides))
data$percent_coverage <- as.numeric(as.character(data$percent_coverage))

#
# probe for model performance values
#
roc_dat <- xml_data$protein_summary_header$program_details$proteinprophet_details
w_dp <- which(names(roc_dat) == 'protein_summary_data_filter')
w_ep <- which(names(roc_dat) == 'error_point')

roc_df <- do.call(rbind, lapply(lapply(roc_dat[w_dp], unlist), "[",
                                unique(unlist(c(sapply(roc_dat[w_dp],names))))))
row.names(roc_df) <- 1:dim(roc_df)[1]
roc_df <- as.data.frame(roc_df)
roc_df$analysts <- 'roc_data_point'
colnames(roc_df)[c(1,3,4,5)] <- c('min_prob','error','num_corr','num_incorr')

err_df <- do.call(rbind, lapply(lapply(roc_dat[w_ep], unlist), "[",
                                unique(unlist(c(sapply(roc_dat[w_ep],names))))))
row.names(err_df) <- 1:dim(err_df)[1]
err_df <- as.data.frame(err_df)
err_df$analysts <- 'error_point'
err_df$sensitivity <- NA

pdf <- rbind(roc_df, err_df)
row.names(pdf) <- 1:dim(pdf)[1]

pdf$min_prob <- as.numeric(as.character(pdf$min_prob))
pdf$sensitivity <- as.numeric(as.character(pdf$sensitivity))
pdf$error <- as.numeric(as.character(pdf$error))
pdf$num_corr <- as.numeric(as.character(pdf$num_corr))
pdf$num_incorr <- as.numeric(as.character(pdf$num_incorr))

w_pcf <- which(pdf$error == max(pdf[pdf$error <= fdr_cutf &
                                        pdf$analysts == 'error_point',]$error) &
                   pdf$analysts == 'error_point')
min_prob <- pdf[w_pcf,]$min_prob
cat("\n")
print(pdf[w_pcf,], row.names = FALSE)
cat("\n")

#
# filter data by cutoff min_prob and rank
#
data <- data[data$probability >= min_prob & data$total_number_distinct_peptides >= min_npep & !is.na(data$group_number),]

pdata <- unique(data[,c('group_number','probability','protein_name',
                 'n_indistinguishable_proteins','percent_coverage',
                 'total_number_peptides','total_number_distinct_peptides')])

pdata <- pdata[order(-pdata$percent_coverage),]

print(head(pdata[,c(3,5,7)]), row.names = FALSE)

write.csv(data, path_csv, row.names = FALSE)
cat("\nFINISHED:\n")
cat(" data written to                  ", path_csv, "\n")

path_csv <- sub(".csv", "_summary.csv", path_csv)
write.csv(pdata, path_csv, row.names = FALSE)
cat("\nFINISHED:\n")
cat(" protein summary written to       ", path_csv, "\n")

path_csv <- sub(".csv", "_performance.csv", path_csv)
write.csv(pdf, path_csv, row.names = FALSE)
cat(" performance roc written to       ", path_csv, "\n")