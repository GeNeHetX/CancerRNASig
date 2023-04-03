# CancerRNASig
RNA signature functions

## Installation

Install CancerRNASig_Rpackage by downloading the zip code : 

```bash
 Rscript -e "devtools::install_local('nameZipFile.zip')"
```

## Usage

 * **Purist** : Function compute the PurIST score (Purity Independent Subtyping of Tumors)
 
    → return a data frame with row names as colnames of newexp, first column is purist subtype, second is purist score
<br></br>

 * **McpCount** : 

     → return a data frame with samples in row and MCPcounter population quantification in columns
 <br></br>    

 * **ProjICA** : Function compute simple ICA projection (default using Puleo et al components)

     → return projected components
 <br></br>    

 * **Estimate** : 
 
     → return score data
