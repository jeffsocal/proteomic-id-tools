# Peptide-Protein ID - R Helper Tools

## Introduction
These R scripts are an attempt to do a pure XML to CSV conversion for TPP files originating out of ProteoWizard's msconvert, XTandem, PeptideProphet, and ProteinProphet, relying only on a single library R:XML. 

## Set-up

1. Install and setup R

          sudo apt-get install libxml2-dev  
          R
          > install.packages(c('XML','progress'))
          > q()

2. Checkout the R helper tools repository.

        git clone https://github.com/jeffsocal/tpp-post-id-tools.git
        
## Running the Scripts
You can run the scripts detached from the R terminal, as a command line execution or even bash scripted. The inputs will be looking for the specific _pepXML_, _protXML_, and _mzXML_ as a super-simple method of validating the input.

These scripts attempt to dump the XML out as a .csv file that can be universally utilized in down stream analyses. 

          Rscript --slave pepxml2csv.R <path_to.pepXML> <path_to.pepXML.csv>
          Rscript --slave proxml2csv.R <path_to.protXML> <path_to.protXML.csv>
          Rscript --slave mzxml2csv.R <path_to.mzXML> <path_to.mzXML.csv>
          

### MERGE ms2 scans with sequence ids
This script is intending to match up the ids with the originating scans so that accounting of search utility (how many spectra were id'd) can be assessed. Along with providing an ability to tune CID parameters, etc. 

          Rscript --slave < merge_scans-ids.R --args <path_to.mzXML.csv> <path_to.pepXML.csv>

### CLUSTER ms2 scans with ms1 features
This script is intending to cluster ids with ms1 features so that accounting of search utility (how many features were id'd) can be assessed, etc. 

          R --slave < cluster_ms1-ms2.R --args \
          --out=<file_out.csv> \
          --xf=<path_to.ms1features.csv> \
          --yf=<path_to.protXML.csv> \
          --xmz=<x_col_mz> \                      # x's table column of m/z values
          --xrt=<x_col_rt> \                      # x's table column of rt values
          --xid=<x_col_id> \                      # x's table column of feature ids
          --ymz=<y_col_mz> \                      # y's table column of m/z values
          --yrt=<y_col_rt> \                      # y's table column of rt values
          --yid=<y_col_id> \                      # y's table column of scans ids
          --tmz=<tolerance_mz_daltons> \          # m/z clustering tolerance (default: 0.1 daltons)
          --trt=<tolerance_rt_sec> \              # rt clustering tolerance (default: 10 sec)
          --seg=<data_segments>                   # data segmenting for large files (default: 1000)
          
Typical stdout 

          cluster ms1-ms2 started
          reference file:       data/ms1_values.csv
          cluster file:         data/ms2_values.csv
          output file:          data/cluster_out.csv
          clustering [============>-------------------]  40% eta:  3s
          ...
          clustering [================================] 100% eta:  0s
           
          FINISHED: data written to  data/test.csv
          
          
