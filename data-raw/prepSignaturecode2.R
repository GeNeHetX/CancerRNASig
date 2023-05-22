# ---


# --- --- --- --- --- --- --- --- --- --- --- ---
#load and func
{
  library(tidyr)
  library(jsonlite)
  library(qutils)
  library(openxlsx)

  addgs=function(gsobj,gset,type,src,id,addIDtoNameList=T){
    if(length(type)==1){type=rep(type,length(gset))}
    if(length(src)==1){src=rep(src,length(gset))}
    if(length(id)==1){id=rep(id,length(gset))}

    gsid=names(gset) 
    if(addIDtoNameList){
      gsid=paste0(id[1],"_",names(gset)) 
    }
    gsobj$gset=c(setNames(gset,gsid),gsobj$gset)
    gsobj$type=c(setNames(type,gsid),gsobj$type)
    gsobj$src=c(setNames(src,gsid),gsobj$src)
    gsobj$id=c(setNames(id,gsid),gsobj$id)
    gsobj
  }

  cleangid=function(x){unique(x[which(!is.na(x)& !x %in% c(""," "))])}
}



# --- --- --- --- --- --- --- --- --- --- --- ---
# PDAC
{ setwd('C:/Users/cpignolet/Documents/VisualStudioCode/PackageR/CancerRNASig/data-raw')
  print(load("ChanSengYueSigs.RData"))
  rnami=c(1,2,6,10)
  rnam=c("ClassicalA","BasallikeA","ClassicalB","BasallikeB")
  names(ChanSengYueSigs)[rnami]=rnam

  print(load("pdxGeneL.RData"))
  print(load("puleoG.RData"))
  puleoG=lapply(puleoG,\(x){ unique(unlist(strsplit(x," /// ")))})

  a=read.delim("PDAssigner_human_annot.txt",as.is=T,header=F,sep="\t")[,2:1]
  pdassigner=split(a[,2],a[,1])

  baileyMarkG=lapply(list(
    ADEX="PDAC_ADEXvRest_suptab2.txt",
    Immunogenic="PDAC_ImmunovRest_suptab2.txt",
    Progenitor="PDAC_ProgenitorvRest_suptab2.txt",
    Squamous="PDAC_SquamousvRest_suptab2.txt"),function(f){
    x=read.delim(f,header=T,sep="\t",as.is=T)
    x=x[which(x$logFC>0&x$adj.P.Val<0.01),]
    rownames(x)=x$StableId
    cleangid(x$Symbol)
  }
  )

  allmoffsigtab=read.delim("moffittMetagene.txt",sep="\t",as.is=T)
  allmoffsig=setNames(split(allmoffsigtab$symbol,allmoffsigtab$factor),colnames(allmoffsigtab)[3:16])
  names(allmoffsig) = gsub("^F.+_","",names(allmoffsig))

}
# --- --- --- --- --- --- --- --- --- --- --- ---
# CAF
{
  turleyCaf=openxlsx::read.xlsx("TurleySupTab5.xlsx")
  turleyCaf=split(turleyCaf[,2],turleyCaf[,1])
  names(turleyCaf)=c("hCAF0TGFB","hCAF1early","hCAF2il6lif")

  FMGcd2020=openxlsx::read.xlsx("6_DataS2_22March2020.xlsx",sheet=2)
  colnames(FMGcd2020)=gsub(" |-",".", FMGcd2020[1,]);FMGcd2020=FMGcd2020[-1,]
  FMGcd2020=lapply(as.list(FMGcd2020),narm)

  TuvesonCAFsc=list(iCAF=read.xlsx("215864_2_supp_5552227_ps91bp.xlsx",sheet=1)[,2],
    myo=read.xlsx("215864_2_supp_5552227_ps91bp.xlsx",sheet=2)[,2])

  CAFgenesigsymL=readRDS("CAFgenesigsymL.rds")[c("POSTN","MYH11","PDPN","subPOSTN")]

}
# --- --- --- --- --- --- --- --- --- --- --- ---
# IMMUNE
{
  print(load("mcpgenes.RData"))
  mcpgenes=setNames(mcpgenes[,c(2,1,3)],c("class","gene","ENTREZID"))


  ImmuneL=list(ICKrelated=toupper(scan( "YLickrelated.txt",what="character",sep="\n")),
    GeneralImmu=toupper(read.xlsx( "JCI121924.sdt1-6.xlsx",sheet=6)[,1]))


}
# --- --- --- --- --- --- --- --- --- --- --- ---
# CCK: Cjholangiocarcnoma
{
  stimclassg=read.delim("STIMClassifier.txt",sep="\t",as.is=T,header=T)
  STIMgenes=split(stimclassg[,1],stimclassg[,2])
  names(STIMgenes)=gsub(" ","_",names(STIMgenes))



  Biclassgenes=list(PROLIFGenes=scan("proliferationGenes.txt",what="character"),
    INFLAMGenes=scan("inflammationGenes.txt",what="character"))


}

# --- --- --- --- --- --- --- --- --- --- --- ---
# DrugBank 
{
  fils=list.files("drugBank")
  drugbank=setNames(lapply(fils,\(f){scan(paste0("drugBank/",f),what="character",sep="\n")}),sub(".txt$","",fils))
}





gsignatures=list(
  gset=list(),
  type=c(),
  src=c(),
  id=c()
  )
gsignatures=gsignatures%>%
addgs(gset=ChanSengYueSigs,type="PDAC",src="Chan-Seng-Yue.etal;PMID.31932696",id="PDAC_CSY20")%>%
addgs(gset=pdxGeneL,type="PDAC",src="Nicolle.etal;PMID.29186684",id="PDAC_PDXph1")%>%
addgs(gset=puleoG,type="PDAC",src="Puleo.etal;PMID.30165049",id="PDAC_Puleo")%>%
addgs(gset=pdassigner,type="PDAC",src="Collisson.etal;PMID.21460848",id="PDAC_PDAssigner")%>%
addgs(gset=allmoffsig,type="PDAC",src="Moffitt.etal;PMID.26343385",id="PDAC_Moffitt15")%>%

addgs(gset=turleyCaf,type="CAF",src="Dominguez.etal;PMID.31699795",id="CAF_Turley20")%>%
addgs(gset=FMGcd2020,type="CAF",src="Kieffer.etal;PMID.32434947",id="CAF_FMG20")%>%
addgs(gset=TuvesonCAFsc,type="CAF",src="Elyada.etal;PMID.31197017",id="CAF_Tuveson19")%>%
addgs(gset=CAFgenesigsymL,type="CAF",src="Neuzillet.etal;PMID.36102377",id="CAF_Neuzillet22")%>%

addgs(gset=split(mcpgenes[,2],mcpgenes[,1]),type="Immu",src="Becht.etal;PMID.27765066",id="IMMU_MCPcounter")%>%
addgs(gset=ImmuneL,type="Immu",src="Rodrigues.etal;PMID.30179225",id="IMMU_GenJCI121924")%>%


addgs(gset=STIMgenes,type="CCK",src="MartinSerrano.etal;PMID.35584893",id="CCK_STIM")%>%
addgs(gset=Biclassgenes,type="CCK",src="Sia.etal;PMID.23295441",id="CCK_Sia13")%>%

addgs(gset=drugbank,type="Drug",src="DrugBankDec2022",id="DrugBank")%>%


addgs(gset=list(PSCcaf=scan("ECMsignature_PMID34548310.txt",what="character",sep="\n")),type="ECM",src="Helms.etal;PMID.34548310",id="ECM_Helms22")






A=list(gset=gsignatures$gset,annot=as.data.frame(gsignatures[2:4]))

write_json(toJSON(A), "geneSetSignatures.json", pretty = T) 

