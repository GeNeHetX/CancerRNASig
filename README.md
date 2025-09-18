# CancerRNASig
CancerRNASig is an R package designed to facilitate the application of RNA-seq signatures published in the literature to your own datasets. It also provides an organized collection of published gene signatures related to cancer (mainly PDAC), stroma, and immunity.

## Prerequisites
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

install.packages("devtools")

install.packages("utils", repos="http://r-forge.r-project.org", dependencies=TRUE) 

install.packages("estimate", repos="http://r-forge.r-project.org", dependencies=TRUE)

BiocManager::install("DESeq2")

devtools::install_github("GeNeHetX/qutils")
```
## Installation & Usage
The CancerRNASig R package can be installed directly from GitHub in a terminal using Rscript:
```bash
Rscript -e "if (!require('devtools')) install.packages('devtools', repos='https://cloud.r-project.org'); devtools::install_github('GeNeHetX/CancerRNASig')"

```
or within an R environment:
```r
# Install devtools if needed
install.packages("devtools")

# Install CancerRNASig from GitHub
devtools::install_github("GeNeHetX/CancerRNASig")
```

### Input formats

**counts_matrix** : data frame or matrix containing gene expression values, with samples in columns and genes in rows. 
```r
### example of counts_matrix
counts_matrix <- read.csv("my_dataset_expression_matrix.tsv.gz", sep="\t",header=TRUE)

head(counts_matrix)
                SAMPLE_001 SAMPLE_002 SAMPLE_003 SAMPLE_004
ENSG00000160072      4          8         21          4
ENSG00000234396      0          0          1          0
ENSG00000225972    515        507        690        148
ENSG00000224315      0          0          0          0
ENSG00000198744     57        496        627         66
ENSG00000279928    125        108       1499        681
```

**geneSymbols** : vector of gene symbols corresponding to your reference genome. This vector must have the **same length as the number of rows in counts_matrix**.
```r
### example of geneSymbols
geneannot <- read.delim('geneAnnot.tsv', sep="\t")
head(geneannot)
                seqname   start     end strand          GeneID        GeneName
ENSG00000160072       1 1471765 1497848      + ENSG00000160072          ATAD3B
ENSG00000234396       1 2212523 2220738      + ENSG00000234396 ENSG00000234396
ENSG00000225972       1  629062  629433      + ENSG00000225972        MTND1P23
ENSG00000224315       1 8786211 8786913      - ENSG00000224315          RPL7P7
ENSG00000198744       1  634376  634922      + ENSG00000198744        MTCO3P12
ENSG00000279928       1  182696  184174      + ENSG00000279928        DDX11L17


geneSymbols <- geneannot$GeneName
head(geneSymbols)
> "ATAD3B" "ENSG00000234396" "MTND1P23" "RPL7P7" "MTCO3P12" "DDX11L17"

length(geneSymbols)
> 61809
dim(counts_matrix)[1]
> 61809
```

**ICAgw** : data frame with your components in columns and your genes (ENSG IDs or gene symbols) in rows.
```r
### example of ICAgw
head(CancerRNASig:::puleoICAgw[,1:3])
>                 Exocrine  Endocrine    Classic
CPB1              23.80349 -0.3141979  0.5453618
CELA3A /// CELA3B 22.51492 -1.3706352  0.5733326
PLA2G1B           21.70120 -0.6311313 -0.4521336
CEL               21.36494  0.2416872 -0.8614301
REG1B             21.04806  2.7321668 -1.1594257
PNLIP             21.02681 -2.3157584  0.1250148

```

### Usage
Load the package and use the main function **callSignature** :
```r
library(CancerRNASig)

res <- callSignature(
  matrix       = counts_matrix,       # expression matrix (ENSG in rows, samples in columns)
  geneSymbols  = geneSymbols,        # vector of gene symbols
  signature    = "Gempred",        # signature to apply
  normType       = "uq",             # normalization (raw, uq, vst)
  scaleType      = "sc"             # scaling (raw, sc, gc, gsc, ssc)
)
```


To project your dataset onto custom ICA components, simply use the **.qProjICA** function.
```r
res <- .qProjICA(
  newexp        = raw_counts,        # expression matrix (genes in rows, samples in columns)
  ICAgw         = yourICAgw,           # ICA gene weights (default: Puleo et al.)
  geneNormType  = "sc",                 # normalization of gene expression (sc: sample scale)
  projNormType  = "raw",                # normalization of components (raw by default)
  ming          = 500                   # minimum number of overlapping genes
)
```

## Available Signatures
The main function **callSignature()** applies several transcriptomic signatures.
Each subsection below describes the method, the expected output, and provides the publication reference.

 * **Purist** : 

     → **Description**: PurIST (Purity Independent Subtyping of Tumors) is a transcriptomic classifier for pancreatic ductal adenocarcinoma that separates tumors into basal-like and classical subtypes independently of tumor purity. These subtypes have distinct prognostic and therapeutic implications.<br>
    → **Output**: a data frame with samples as rows and two columns: the **purist score (NumPurist)** and the **purist class** (classic/basal) <br>
    → [DOI: 10.1158/1078-0432.CCR-19-1467](https://pmc.ncbi.nlm.nih.gov/articles/PMC6942634/)
<br></br>

 * **Mcpcount** : 

     → **Description**: MCP-counter is a transcriptomic method that estimates the abundance of key immune and stromal cell populations from bulk RNA-seq data. It provides robust, sample-independent scores that enable cross-sample comparison and tumor microenvironment characterization.<br>
     → **Output**: a data frame with samples as rows and the relative proportion of several immune and stroma cell populations (Tcells, CD8Tcells, Cytotox.lymph, NK, B.lineage, Mono.lineage, Myeloid.dendritic, Neutrophils, Endothelial, Fibroblasts) in columns.<br>
     → [DOI: 10.1186/s13059-016-1070-5](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1070-5)
 <br></br>    

 * **Puleo** : 

    → **Description**: A 403-gene transcriptomic signature derived from pancreatic ductal adenocarcinoma (PDAC) samples, designed to stratify tumors into five molecular subtypes: **Pure Classical, Immune Classical, Desmoplastic, Stroma-Activated, and Pure Basal-Like**. This classification integrates both tumor cell–intrinsic programs and tumor microenvironmental (stromal and immune) components. <br>
     → **Output**: a data frame with samples as rows and the relative proportion of several immune and stroma cell populations (Exocrine, Endocrine, Classic, StromaActiv, Basal, StromaActivInflam, Immune, StromaInactive, ICA9) in columns. <br>
     → [DOI: 10.1053/j.gastro.2018.08.033](https://pubmed.ncbi.nlm.nih.gov/30165049/)
 <br></br>

 * **Gempred** : 

     → **Description**: GemPred is an RNA-based transcriptomic signature developed for pancreatic ductal adenocarcinoma that predicts sensitivity to adjuvant gemcitabine in PDAC patients. Patients with GemPred-positive tumors benefit significantly more from gemcitabine than GemPred-negative cases.<br>
     → **Output**: a data frame with samples as rows and a single column containing the GemPred score.<br>
     → [DOI: 10.1016/j.annonc.2020.10.601](https://www.sciencedirect.com/science/article/pii/S092375342043131X)

* **tGempred** : 

     → **Description**: tGemPred is an improved GemPred signature optimized for biopsy samples, using ICA deconvolution to remove contamination signals (e.g., blood) and enhance predictive performance.<br>
     → **Output**: a data frame with samples as rows and a single column containing the tGemPred score.<br>
     → DOI: coming soon

### Note for scRNA-seq data
These signatures are designed for **bulk RNA-seq** but can also be applied to **pseudo-bulk scRNA-seq data**.
To use scRNA-seq data, first aggregate counts per sample or cluster (pseudo-bulk), normalize appropriately, and then apply callSignature

## Available Gene Sets
All gene sets and their annotations are available in the **signatures** object included in the package.

<table style="border-bottom:0; width: auto !important; margin-left: auto; margin-right: auto;border-bottom: 0;" class="table">
 <thead>
  <tr>
   <th style="text-align:left;"> Annotation </th>
   <th style="text-align:left;"> Source </th>
   <th style="text-align:left;"> Type </th>
   <th style="text-align:right;"> NumberOfSignatures </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> BUSSLINGER.HUMAN </td>
   <td style="text-align:left;"> G.Busslinger.etal;PMID: 33691112 </td>
   <td style="text-align:left;"> Normal Digestive scRNA-seq </td>
   <td style="text-align:right;"> 20 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BockerstettGut </td>
   <td style="text-align:left;"> K.Bockerstett;PMID:31481545 </td>
   <td style="text-align:left;"> Normal Digestive scRNA-seq </td>
   <td style="text-align:right;"> 20 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAF_FMG20 </td>
   <td style="text-align:left;"> Kieffer.etal;PMID.32434947 </td>
   <td style="text-align:left;"> Cancer-Associated Fibroblasts (CAF) </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAF_Neuzillet22 </td>
   <td style="text-align:left;"> Neuzillet.etal;PMID.36102377 </td>
   <td style="text-align:left;"> Cancer-Associated Fibroblasts (CAF) </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAF_Turley20 </td>
   <td style="text-align:left;"> Dominguez.etal;PMID.31699795 </td>
   <td style="text-align:left;"> Cancer-Associated Fibroblasts (CAF) </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAF_Tuveson19 </td>
   <td style="text-align:left;"> Elyada.etal;PMID.31197017 </td>
   <td style="text-align:left;"> Cancer-Associated Fibroblasts (CAF) </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CCCA_MetaProg </td>
   <td style="text-align:left;"> Gavish.etal;PMID.37258682 </td>
   <td style="text-align:left;"> Cancer </td>
   <td style="text-align:right;"> 41 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CCK_STIM </td>
   <td style="text-align:left;"> MartinSerrano.etal;PMID.35584893 </td>
   <td style="text-align:left;"> Cholangiocarcinoma (CCK) </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CCK_Sia13 </td>
   <td style="text-align:left;"> Sia.etal;PMID.23295441 </td>
   <td style="text-align:left;"> Cholangiocarcinoma (CCK) </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ECM_Helms22 </td>
   <td style="text-align:left;"> Helms.etal;PMID.34548310 </td>
   <td style="text-align:left;"> Extracellular Matrix (ECM) </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FibroAtlasGao </td>
   <td style="text-align:left;"> Yang-Gao.etal;PMID.39303725 </td>
   <td style="text-align:left;"> Fibroblast </td>
   <td style="text-align:right;"> 20 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IMMU_GenJCI121924 </td>
   <td style="text-align:left;"> Rodrigues.etal;PMID.30179225 </td>
   <td style="text-align:left;"> Immune Cells </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IMMU_MCPcounter </td>
   <td style="text-align:left;"> Becht.etal;PMID.27765066 </td>
   <td style="text-align:left;"> Immune Cells </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IMMU_Neutroatlas </td>
   <td style="text-align:left;"> Wu.etal;PMID.38447573 </td>
   <td style="text-align:left;"> Immune Cells </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IMMU_Tcellatlas </td>
   <td style="text-align:left;"> Chu.etal;PMID.37248301 </td>
   <td style="text-align:left;"> Immune Cells </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KIM.scRNAGatricCarcniogegenisis.cell </td>
   <td style="text-align:left;"> J.Kim;PMID:35087207 </td>
   <td style="text-align:left;"> Normal Digestive scRNA-seq </td>
   <td style="text-align:right;"> 20 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MA.MOUSE.STOMACH </td>
   <td style="text-align:left;"> Z.Ma.etal;PMID: 34695382 </td>
   <td style="text-align:left;"> Normal Digestive scRNA-seq </td>
   <td style="text-align:right;"> 20 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> OrganoidAtlas </td>
   <td style="text-align:left;"> Xu.etal;PMID: 40355592 </td>
   <td style="text-align:left;"> Organoid </td>
   <td style="text-align:right;"> 48 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PDAC_Bailey16 </td>
   <td style="text-align:left;"> Bailey.etal;PMID.26909576 </td>
   <td style="text-align:left;"> Pancreatic Ductal Adenocarcinoma (PDAC) </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PDAC_CSY20 </td>
   <td style="text-align:left;"> Chan-Seng-Yue.etal;PMID.31932696 </td>
   <td style="text-align:left;"> Pancreatic Ductal Adenocarcinoma (PDAC) </td>
   <td style="text-align:right;"> 12 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PDAC_Hwang22 </td>
   <td style="text-align:left;"> Hwang.etal;PMID.35902743 </td>
   <td style="text-align:left;"> Pancreatic Ductal Adenocarcinoma (PDAC) </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PDAC_Moffitt15 </td>
   <td style="text-align:left;"> Moffitt.etal;PMID.26343385 </td>
   <td style="text-align:left;"> Pancreatic Ductal Adenocarcinoma (PDAC) </td>
   <td style="text-align:right;"> 14 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PDAC_PDAssigner </td>
   <td style="text-align:left;"> Collisson.etal;PMID.21460848 </td>
   <td style="text-align:left;"> Pancreatic Ductal Adenocarcinoma (PDAC) </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PDAC_PDXph1 </td>
   <td style="text-align:left;"> Nicolle.etal;PMID.29186684 </td>
   <td style="text-align:left;"> Pancreatic Ductal Adenocarcinoma (PDAC) </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PDAC_Puleo </td>
   <td style="text-align:left;"> Puleo.etal;PMID.30165049 </td>
   <td style="text-align:left;"> Pancreatic Ductal Adenocarcinoma (PDAC) </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SCHLESINGER.MOUSE </td>
   <td style="text-align:left;"> Y.Schlesinger;PMID:32908137 </td>
   <td style="text-align:left;"> Normal Digestive scRNA-seq </td>
   <td style="text-align:right;"> 20 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Sathe.scrna </td>
   <td style="text-align:left;"> A.Sathe;PMID:32060101 </td>
   <td style="text-align:left;"> Normal Digestive scRNA-seq </td>
   <td style="text-align:right;"> 20 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ductal_pancreas_mouse_atlas </td>
   <td style="text-align:left;"> Y.Schlesinger;PMID:38908487 </td>
   <td style="text-align:left;"> Normal Digestive scRNA-seq </td>
   <td style="text-align:right;"> 20 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> scIBD </td>
   <td style="text-align:left;"> H.Nie;PMID:38177426 </td>
   <td style="text-align:left;"> Normal Digestive scRNA-seq </td>
   <td style="text-align:right;"> 20 </td>
  </tr>
</tbody>
<tfoot>
<tr><td style="padding: 0; " colspan="100%"><span style="font-style: italic;">Note: </span></td></tr>
<tr><td style="padding: 0; " colspan="100%">
<sup></sup> last update: 10/09/2025</td></tr>
</tfoot>
</table>

### For non-R users
A [json file](https://github.com/GeNeHetX/CancerRNASig/blob/main/data-raw/geneSetSignatures.json) is available with all the signatures. It contains the same data as the signatures.rda file and can be openned with other programming langages such as **Python**.

