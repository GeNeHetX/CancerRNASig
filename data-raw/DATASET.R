# ---
# source("data-raw/DATASET.R")

# --- --- --- --- --- --- --- --- --- --- --- ---
#load and func
{
  library(tidyr)
  library(jsonlite)
  library(qutils)
  library(openxlsx)

  addgs=function(gsobj,geneset,type,src,id,addIDtoNameList=T){
    if(length(type)==1){type=rep(type,length(geneset))}
    if(length(src)==1){src=rep(src,length(geneset))}
    if(length(id)==1){id=rep(id,length(geneset))}

    gsid=names(geneset) 
    if(addIDtoNameList){
      gsid=paste0(id[1],"_",names(geneset)) 
    }
    gsobj$geneset=c(setNames(geneset,gsid),gsobj$geneset)
    gsobj$type=c(setNames(type,gsid),gsobj$type)
    gsobj$src=c(setNames(src,gsid),gsobj$src)
    gsobj$id=c(setNames(id,gsid),gsobj$id)
    gsobj
  }

  cleangid=function(x){unique(x[which(!is.na(x)& !x %in% c(""," "))])}
}

.refpath = 'data-raw/refData'

# --- --- --- --- --- --- --- --- --- --- --- ---
# PDAC
{ 
  # print(load(file.path(.refpath,"ChanSengYueSigs.RData")))
  # write_json(toJSON(ChanSengYueSigs),file.path(.refpath, "ChanSengYueSigs.json"), pretty = T) 
  ChanSengYueSigs=fromJSON(read_json(file.path(.refpath, "ChanSengYueSigs.json"),simplifyVector=T))


  rnami=c(1,2,6,10)
  rnam=c("ClassicalA","BasallikeA","ClassicalB","BasallikeB")
  names(ChanSengYueSigs)[rnami]=rnam

  # print(load(file.path(.refpath,"pdxGeneL.RData")))
  # write_json(toJSON(pdxGeneL),file.path(.refpath, "pdxGeneL.json"), pretty = T) 
  pdxGeneL=fromJSON(read_json(file.path(.refpath, "pdxGeneL.json"),simplifyVector=T))

  # print(load(file.path(.refpath,"puleoG.RData")))
  # write_json(toJSON(puleoG),file.path(.refpath, "puleoG.json"), pretty = T) 
  puleoG=fromJSON(read_json(file.path(.refpath, "puleoG.json"),simplifyVector=T))





  puleoG=lapply(puleoG,\(x){ unique(unlist(strsplit(x," /// ")))})

  a=read.delim(file.path(.refpath,"PDAssigner_human_annot.txt"),as.is=T,header=F,sep="\t")[,2:1]
  pdassigner=split(a[,2],a[,1])

  baileyMarkG=lapply(list(
    ADEX=file.path(.refpath,"PDAC_ADEXvRest_suptab2.txt"),
    Immunogenic=file.path(.refpath,"PDAC_ImmunovRest_suptab2.txt"),
    Progenitor=file.path(.refpath,"PDAC_ProgenitorvRest_suptab2.txt"),
    Squamous=file.path(.refpath,"PDAC_SquamousvRest_suptab2.txt")),function(f){
    x=read.delim(f,header=T,sep="\t",as.is=T)
    x=x[which(x$logFC>0&x$adj.P.Val<0.01),]
    rownames(x)=x$StableId
    cleangid(x$Symbol)
  }
  )

  allmoffsigtab=read.delim(file.path(.refpath,"moffittMetagene.txt"),sep="\t",as.is=T)
  allmoffsig=setNames(split(allmoffsigtab$symbol,allmoffsigtab$factor),colnames(allmoffsigtab)[3:16])
  names(allmoffsig) = gsub("^F.+_","",names(allmoffsig))

  #Hwang
  {
    a=read.xlsx(file.path(.refpath,"41588_2022_1134_MOESM4_ESM.xlsx"), sheet=2,startRow = 3)[,-1]

    HwangMarkG=as.list(a[-1,])


    b=sub("^X.+","",colnames(a));for(i in 1:length(b)) if(b[i]=="")b[i]=b[i-1];


    names(HwangMarkG)=sub("\\/","", sub(" ","",sub(")","",sub("-|\\/|\\(","",
      paste(sub("Fibroblast","Fibro",sub("Malignant.","Malign",sub("\\.programs","",b))),
        as.character(a[1,]),sep=".")))))

    

  }

}
# --- --- --- --- --- --- --- --- --- --- --- ---
# CAF
{
  turleyCaf=openxlsx::read.xlsx(file.path(.refpath,"TurleySupTab5.xlsx"))
  turleyCaf=split(turleyCaf[,2],turleyCaf[,1])
  names(turleyCaf)=c("hCAF0TGFB","hCAF1early","hCAF2il6lif")

  FMGcd2020=openxlsx::read.xlsx(file.path(.refpath,"6_DataS2_22March2020.xlsx"),sheet=2)
  colnames(FMGcd2020)=gsub(" |-",".", FMGcd2020[1,]);FMGcd2020=FMGcd2020[-1,]
  FMGcd2020=lapply(as.list(FMGcd2020),narm)

  TuvesonCAFsc=list(iCAF=read.xlsx(file.path(.refpath,"215864_2_supp_5552227_ps91bp.xlsx"),sheet=1)[,2],
    myo=read.xlsx(file.path(.refpath,"215864_2_supp_5552227_ps91bp.xlsx"),sheet=2)[,2])

  CAFgenesigsymL=readRDS(file.path(.refpath,"CAFgenesigsymL.rds"))[c("POSTN","MYH11","PDPN","subPOSTN")]

}
# --- --- --- --- --- --- --- --- --- --- --- ---
# IMMUNE
{


  # print(load(file.path(.refpath,"mcpgenes.RData")))
  # write_json(toJSON(mcpgenes),file.path(.refpath, "mcpgenes.json"), pretty = T) 
  mcpgenes=fromJSON(read_json(file.path(.refpath, "mcpgenes.json"),simplifyVector=T))


  mcpgenes=setNames(mcpgenes[,c(2,1,3)],c("class","gene","ENTREZID"))


  ImmuneL=list(ICKrelated=toupper(scan(file.path(.refpath, "YLickrelated.txt"),what="character",sep="\n")),
    GeneralImmu=toupper(read.xlsx(file.path(.refpath,"JCI121924.sdt1-6.xlsx"),sheet=6)[,1]))



  # Pan-cancer T cell atlas Y.CHU 

  table=read.csv(file.path(.refpath,'41591_2023_2371_MOESM3_ESM.csv'), sep=';')
  # table=read.csv('41591_2023_2371_MOESM3_ESM.csv', sep=';')
  table=table[1:25,2:9]
  colnames(table)=table[1,]
  table=table[-1,]
  PanCK_Tcellatlas=as.list(table)
  PanCK_Tcellatlas =lapply(PanCK_Tcellatlas, function(x) x[x!=""])


  # Pan-cancer Neutrophils, Wu PMID.38447573
  # table=read.xlsx(file.path(.refpath,'1-s2.0-S0092867424001260-mmc2.xlsx'),sheet=2)
  tab=read.delim(file.path(.refpath,'1-s2.0-S0092867424001260-mmc2_sheet2.txt'), sep='\t')
  # tab$ok=tab$Adjusted.P.value< 10^-20 #& tab$Log2.Fold.Change> log2(1.5)
  # tapply(  tab$ok, tab$Cluster, sum)
  PanCanNeutroWu=lapply(split(tab,tab$Cluster),\(x)x[which(rank(x$Adjusted.P.value,ties.method="first")<100),"Gene"] )


}
# --- --- --- --- --- --- --- --- --- --- --- ---
# CCK: Cjholangiocarcnoma
{
  stimclassg=read.delim(file.path(.refpath,"STIMClassifier.txt"),sep="\t",as.is=T,header=T)
  STIMgenes=split(stimclassg[,1],stimclassg[,2])
  names(STIMgenes)=gsub(" ","_",names(STIMgenes))



  Biclassgenes=list(PROLIFGenes=scan(file.path(.refpath,"proliferationGenes.txt"),what="character"),
    INFLAMGenes=scan(file.path(.refpath,"inflammationGenes.txt"),what="character"))


}

# --- --- --- --- --- --- --- --- --- --- --- ---
# DrugBank 
{
  #setwd(.refpath)
  fils=list.files(file.path(.refpath,"drugBank"))
  drugbank=setNames(lapply(fils,\(f){scan(file.path(.refpath,"drugBank",f),what="character",sep="\n")}),sub(".txt$","",fils))
}


# --- --- --- --- --- --- --- --- --- --- --- ---
# cancer atlas
{
  canceratlasMP=as.list(openxlsx::read.xlsx(file.path(.refpath,"41586_2023_6130_MOESM6_ESM.xlsx"),sheet=1))
  names(canceratlasMP)=gsub("\\.{2,3}",".",  gsub("-|\\/|\\(|\\)",".",names(canceratlasMP)))

  canceratlasRobNMF=as.list(openxlsx::read.xlsx(file.path(.refpath,"41586_2023_6130_MOESM6_ESM.xlsx"),sheet=2))

  }

{
  #Cell Atlas Organoid :
  atlas_organoid_df = read.xlsx(file.path(.refpath, "Cell_Atlas_Organoid.xlsx"),sheet=5)
  atlas_organoidVec = lapply(atlas_organoid_df, function(x) as.character(x))
  names(atlas_organoidVec) = paste0("Atlas_OrganoidVec_", colnames(atlas_organoid_df))
}
}

# --- --- --- --- --- --- --- --- --- --- --- ---
# ECM
{
  ECMHELMS=fromJSON(read_json(file.path(.refpath,"ECM_Helms_genesets.json"),simplifyVector=T))
}

# --- --- --- --- --- --- --- --- --- --- --- ---
# Fibroblast Atlas

tab_fibro=openxlsx::read.xlsx(file.path(.refpath,"mmc3.xlsx"),sheet=2)

### Nomage des clusters à leur signature biologique
fibroSigs <- c("c01", "c02", "c03", "c04", "c05", "c06", "c07", "c08", "c09", "c10",
          "c11", "c12", "c13", "c14", "c15", "c16", "c17", "c18", "c19", "c20")

#SMC = smooth muscle cells
names(fibroSigs) <- c("SMC_MYH11", "LaminaPropriaFibro_ADAMDEC",
                      "Progenitor-like-fibro_COL15A1", "MyoFibro_LRRC15",
                      "Progenitor-like-fibroblast_PI16", "AlveolarFibro_ADH1B+",
                      "InflaFibro_IL6", "Pericyte", "Fibro_CTNNB1",
                      "SynovialLiningFibro_PRG4", "MyoFibro_HOPX", "Mesothelial",
                      "InflaFibro_HGF", "InflaFibro_HSPA6", "EpithelialCryptFibro_SOX6",
                      "MyoFibro_SFRP2", "ProlifevrativeFibro_STMN1",
                      "SMC_HHIP", "Myofibro_MMP1", "ApFibro_CD74"
                      )

tab_fibro$clusterName <- names(fibroSigs[as.numeric(sub("c", "", tab_fibro$cluster))])


FibroAtlasSigs <- split(tab_fibro[,"gene"],tab_fibro$clusterName)

# --- --- --- --- --- --- --- --- --- --- --- ---
# Aggreg sigs

gsignatures=list(
  geneset=list(),
  type=c(),
  src=c(),
  id=c()
  )

gsignatures=gsignatures%>%
addgs(geneset=ChanSengYueSigs,type="PDAC",src="Chan-Seng-Yue.etal;PMID.31932696",id="PDAC_CSY20")%>%
addgs(geneset=pdxGeneL,type="PDAC",src="Nicolle.etal;PMID.29186684",id="PDAC_PDXph1")%>%
addgs(geneset=puleoG,type="PDAC",src="Puleo.etal;PMID.30165049",id="PDAC_Puleo")%>%
addgs(geneset=pdassigner,type="PDAC",src="Collisson.etal;PMID.21460848",id="PDAC_PDAssigner")%>%
addgs(geneset=allmoffsig,type="PDAC",src="Moffitt.etal;PMID.26343385",id="PDAC_Moffitt15")%>%
addgs(geneset=baileyMarkG,type="PDAC",src="Bailey.etal;PMID.26909576",id="PDAC_Bailey16")%>%
addgs(geneset=HwangMarkG,type="PDAC",src="Hwang.etal;PMID.35902743",id="PDAC_Hwang22")%>%

addgs(geneset=turleyCaf,type="CAF",src="Dominguez.etal;PMID.31699795",id="CAF_Turley20")%>%
addgs(geneset=FMGcd2020,type="CAF",src="Kieffer.etal;PMID.32434947",id="CAF_FMG20")%>%
addgs(geneset=TuvesonCAFsc,type="CAF",src="Elyada.etal;PMID.31197017",id="CAF_Tuveson19")%>%
addgs(geneset=CAFgenesigsymL,type="CAF",src="Neuzillet.etal;PMID.36102377",id="CAF_Neuzillet22")%>%

addgs(geneset=split(mcpgenes[,2],mcpgenes[,1]),type="Immu",src="Becht.etal;PMID.27765066",id="IMMU_MCPcounter")%>%
addgs(geneset=ImmuneL,type="Immu",src="Rodrigues.etal;PMID.30179225",id="IMMU_GenJCI121924")%>%


addgs(geneset=STIMgenes,type="CCK",src="MartinSerrano.etal;PMID.35584893",id="CCK_STIM")%>%
addgs(geneset=Biclassgenes,type="CCK",src="Sia.etal;PMID.23295441",id="CCK_Sia13")%>%

addgs(geneset=drugbank,type="Drug",src="DrugBankDec2022",id="DrugBank")%>%


addgs(geneset=list(PSCcaf=scan(file.path(.refpath,"ECMsignature_PMID34548310.txt"),what="character",sep="\n")),type="ECM",src="Helms.etal;PMID.34548310",id="ECM_Helms22")%>%

addgs(geneset=canceratlasMP,type="Cancer",src="Gavish.etal;PMID.37258682",id="CCCA_MetaProg")%>%

addgs(geneset=PanCK_Tcellatlas,type='Immu', src='Chu.etal;PMID.37248301',id="IMMU_Tcellatlas")%>%

addgs(geneset=PanCanNeutroWu,type='Immu', src='Wu.etal;PMID.38447573',id="IMMU_Neutroatlas")%>%

addgs(geneset=FibroAtlasSigs, type="Fibroblast", src="Yang-Gao.etal;PMID.39303725", id="FibroAtlasGao")




signatures=list(geneset=gsignatures$geneset,
  annotation=as.data.frame(gsignatures[2:4]),
  CancerAtlasAll=canceratlasRobNMF
  )

usethis::use_data(signatures,internal=FALSE,overwrite=T)
# write_json(toJSON(listSig), "geneSetSignatures.json", pretty = T) 
# .geneSetSignatures=fromJSON(read_json("geneSetSignatures.json",simplifyVector=T))
# usethis::use_data(.geneSetSignatures,internal=T,overwrite=T)
# write.table(estimate::SI_geneset,file="data-raw/refData/SI_geneset.tsv",sep="\t")
estimategenes=read.delim(file.path(.refpath,'SI_geneset.tsv'),sep="\t",header=T,colClasses="character")

mcpgenes=read.delim(file.path(.refpath,'mcpgenes.tsv'),sep="\t",header=T,colClasses="character")

puleoICAgw=read.delim(file.path(.refpath,"puleogw.tsv"),sep="\t",header=T)


#usethis::use_data(signatures,.puleoICAgw,.estimategenes,.mcpgenes,internal=FALSE,overwrite=T)

usethis::use_data(puleoICAgw,estimategenes,mcpgenes,internal=TRUE,overwrite=T)

# --- --- --- --- --- --- --- --- --- --- --- ---
# UPDATE README AND DATA FILES
### Write signature in json file for python or other langages
write_json(toJSON(signatures), "./data-raw/geneSetSignatures.json", pretty = T) 

update_RMD <- function(tab_summary){
  # Read README + keep all texte before "Table des Signatures"
  readme_file <- "./README.md"
  readme_content <- readLines(readme_file)
  title_index <- grep("## Table des Signatures", readme_content)
  if (length(title_index) > 0) {
    readme_content <- readme_content[1:(title_index[1]-2)]
  }
  
  #  Add Signatures informations to README
  md_table <- knitr::kable(tab_summary, format = "markdown", col.names = c("Annotation", "Source", "Type", "NumberOfSignatures"))
  cat(readme_content, file = readme_file, sep = "\n")
  cat("\n## Table des Signatures\n", file = readme_file, append = TRUE)
  cat(md_table, file = readme_file, append = TRUE)
  cat(paste0("\n\nlast update: ",format(Sys.Date(), "%d/%m/%Y"),"\n"), file = readme_file, append = TRUE)
}

# Summary table with signatures informations
tab_summary <- signatures$annotation %>%
  group_by(id, src, type) %>%
  summarise(
    type = unique(type),
    nb_signatures = n()
  ) %>%
  ungroup()
colnames(tab_summary) <- c("Annotation","Source","Type","NumberOfSignatures")

# Update files with sigs infos
write.csv(tab_summary, "./data-raw/sigs_summary.csv", row.names = FALSE)
update_RMD(tab_summary)

#library(devtools)
#source("data-raw/DATASET.R")
#build()
#install()
#reload(pkg = ".", quiet = FALSE)
