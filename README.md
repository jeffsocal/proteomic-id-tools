# TPP Peptide-Protein ID Pipeline - R Helper Tools

## Introduction
These R scripts are an attempt to do a pure XML to CSV conversion for TPP files originating out of XTandem, ProteinProphet, and ProteinProphet, relying only on a single library R:XML. 

## Set-up

1. Install and setup R

          sudo apt-get install libxml2-dev  
          R
          > install.packages('XML')
          > q()

2. Checkout the R helper tools repository.

        git clone https://github.com/jeffsocal/tpp-post-id-tools.git
        
## Running the Scripts
You can run the scripts detached from the R terminal, as a command line execution or even bash scripted. The inputs will be looking for the specific _pepXML_, _protXML_, and _mzXML_ as a super-simple method of validating the input.

These scripts attempt to dump the XML out as a .csv file that can be universally utilized in down stream analyses. 

          Rscript --slave pepxml2csv.R 		<path_to.pepXML> 	<path_to.pepXML.csv>
          Rscript --slave proxml2csv.R 		<path_to.protXML> 	<path_to.protXML.csv>
          Rscript --slave mzxml2csv.R 		<path_to.mzXML> 	<path_to.mzXML.csv>
          
This script is intending to match up the ids with the originating scans so that accounting of search utility (how many spectra were id'd) can be assessed. Along with providing an ability to tune CID parameters, etc. 

          Rscript --slave merge_scans-ids.R 	<path_to.mzXML.csv> <path_to.pepXML.csv>
