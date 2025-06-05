# CancerRNASig
RNA signature functions

## Installation

Install CancerRNASig_Rpackage by downloading the zip code : 

```bash
 Rscript -e "devtools::install_github('GeNeHetX/CancerRNASig')"
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

## For none-R users
A [json file](https://github.com/GeNeHetX/CancerRNASig/blob/main/data-raw/geneSetSignatures.json) is available with all the signatures. It contains the same data as the signatures.rda file and can be openned with other programming langages (Python).

## 🧬 Table des Signatures

| Nom de la Signature | Source                          | Type         | Nombre de Signatures |
|:--------------------|:--------------------------------|:-------------|----------------------:|
| CAF_FMG20           | Kieffer *et al.*; PMID: 32434947 | CAF          |                     8 |
| CAF_Neuzillet22     | Neuzillet *et al.*; PMID: 36102377 | CAF        |                     4 |
| CAF_Turley20        | Dominguez *et al.*; PMID: 31699795 | CAF        |                     3 |
| CAF_Tuveson19       | Elyada *et al.*; PMID: 31197017   | CAF        |                     2 |
| CCCA_MetaProg       | Gavish *et al.*; PMID: 37258682   | Cancer     |                    41 |
| CCK_STIM            | Martin-Serrano *et al.*; PMID: 35584893 | CCK    |                     5 |
| CCK_Sia13           | Sia *et al.*; PMID: 23295441      | CCK        |                     2 |
| DrugBank            | DrugBank Dec 2022                 | Drug       |                     5 |
| ECM_Helms22         | Helms *et al.*; PMID: 34548310    | ECM        |                     1 |
| FibroAtlasGao       | Yang-Gao *et al.*; PMID: 39303725 | Fibroblast |                    20 |
| IMMU_GenJCI121924   | Rodrigues *et al.*; PMID: 30179225 | Immu       |                     2 |
| IMMU_MCPcounter     | Becht *et al.*; PMID: 27765066    | Immu       |                    10 |
| IMMU_Neutroatlas    | Wu *et al.*; PMID: 38447573       | Immu       |                    10 |
| IMMU_Tcellatlas     | Chu *et al.*; PMID: 37248301      | Immu       |                     8 |
| PDAC_Bailey16       | Bailey *et al.*; PMID: 26909576   | PDAC       |                     4 |
| PDAC_CSY20          | Chan-Seng-Yue *et al.*; PMID: 31932696 | PDAC   |                    12 |
| PDAC_Hwang22        | Hwang *et al.*; PMID: 35902743    | PDAC       |                    18 |
| PDAC_Moffitt15      | Moffitt *et al.*; PMID: 26343385  | PDAC       |                    14 |
| PDAC_PDAssigner     | Collisson *et al.*; PMID: 21460848 | PDAC      |                     3 |
| PDAC_PDXph1         | Nicolle *et al.*; PMID: 29186684  | PDAC       |                     2 |
| PDAC_Puleo          | Puleo *et al.*; PMID: 30165049    | PDAC       |                    10 |

last update: 24/01/2025

## Table des Signatures
|Annotation                           |Source                           |Type                 | NumberOfSignatures| |:------------------------------------|:--------------------------------|:--------------------|------------------:| |BUSSLINGER.HUMAN                     |G.Busslinger.etal;PMID: 33691112 |scrnaNormalDigestive |                 20| |BockerstettGut                       |K.Bockerstett;PMID:31481545      |scrnaNormalDigestive |                 20| |CAF_FMG20                            |Kieffer.etal;PMID.32434947       |CAF                  |                  8| |CAF_Neuzillet22                      |Neuzillet.etal;PMID.36102377     |CAF                  |                  4| |CAF_Turley20                         |Dominguez.etal;PMID.31699795     |CAF                  |                  3| |CAF_Tuveson19                        |Elyada.etal;PMID.31197017        |CAF                  |                  2| |CCCA_MetaProg                        |Gavish.etal;PMID.37258682        |Cancer               |                 41| |CCK_STIM                             |MartinSerrano.etal;PMID.35584893 |CCK                  |                  5| |CCK_Sia13                            |Sia.etal;PMID.23295441           |CCK                  |                  2| |DrugBank                             |DrugBankDec2022                  |Drug                 |                  5| |ECM_Helms22                          |Helms.etal;PMID.34548310         |ECM                  |                  1| |FibroAtlasGao                        |Yang-Gao.etal;PMID.39303725      |Fibroblast           |                 20| |IMMU_GenJCI121924                    |Rodrigues.etal;PMID.30179225     |Immu                 |                  2| |IMMU_MCPcounter                      |Becht.etal;PMID.27765066         |Immu                 |                 10| |IMMU_Neutroatlas                     |Wu.etal;PMID.38447573            |Immu                 |                 10| |IMMU_Tcellatlas                      |Chu.etal;PMID.37248301           |Immu                 |                  8| |KIM.scRNAGatricCarcniogegenisis.cell |J.Kim;PMID:35087207              |scrnaNormalDigestive |                 20| |MA.MOUSE.STOMACH                     |Z.Ma.etal;PMID: 34695382         |scrnaNormalDigestive |                 20| |OrganoidAtlas                        |Xu.etal;PMID: 40355592           |Organoid             |                 48| |PDAC_Bailey16                        |Bailey.etal;PMID.26909576        |PDAC                 |                  4| |PDAC_CSY20                           |Chan-Seng-Yue.etal;PMID.31932696 |PDAC                 |                 12| |PDAC_Hwang22                         |Hwang.etal;PMID.35902743         |PDAC                 |                 18| |PDAC_Moffitt15                       |Moffitt.etal;PMID.26343385       |PDAC                 |                 14| |PDAC_PDAssigner                      |Collisson.etal;PMID.21460848     |PDAC                 |                  3| |PDAC_PDXph1                          |Nicolle.etal;PMID.29186684       |PDAC                 |                  2| |PDAC_Puleo                           |Puleo.etal;PMID.30165049         |PDAC                 |                 10| |SCHLESINGER.MOUSE                    |Y.Schlesinger;PMID:32908137      |scrnaNormalDigestive |                 20| |Sathe.scrna                          |A.Sathe;PMID:32060101            |scrnaNormalDigestive |                 20| |mmc4                                 |Y.Schlesinger;PMID:38908487      |scrnaNormalDigestive |                 20| |scIBD                                |H.Nie;PMID:38177426              |scrnaNormalDigestive |                 20|

last update: 05/06/2025
