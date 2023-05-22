# ---


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
{ print(load(file.path(.refpath,"ChanSengYueSigs.RData")))
  rnami=c(1,2,6,10)
  rnam=c("ClassicalA","BasallikeA","ClassicalB","BasallikeB")
  names(ChanSengYueSigs)[rnami]=rnam

  print(load(file.path(.refpath,"pdxGeneL.RData")))
  print(load(file.path(.refpath,"puleoG.RData")))
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
  print(load(file.path(.refpath,"mcpgenes.RData")))
  mcpgenes=setNames(mcpgenes[,c(2,1,3)],c("class","gene","ENTREZID"))


  ImmuneL=list(ICKrelated=toupper(scan(file.path(.refpath, "YLickrelated.txt"),what="character",sep="\n")),
    GeneralImmu=toupper(read.xlsx(file.path(.refpath,"JCI121924.sdt1-6.xlsx"),sheet=6)[,1]))


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
  setwd(.refpath)
  fils=list.files("drugBank")
  drugbank=setNames(lapply(fils,\(f){scan(paste0("drugBank/",f),what="character",sep="\n")}),sub(".txt$","",fils))
}





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

addgs(geneset=turleyCaf,type="CAF",src="Dominguez.etal;PMID.31699795",id="CAF_Turley20")%>%
addgs(geneset=FMGcd2020,type="CAF",src="Kieffer.etal;PMID.32434947",id="CAF_FMG20")%>%
addgs(geneset=TuvesonCAFsc,type="CAF",src="Elyada.etal;PMID.31197017",id="CAF_Tuveson19")%>%
addgs(geneset=CAFgenesigsymL,type="CAF",src="Neuzillet.etal;PMID.36102377",id="CAF_Neuzillet22")%>%

addgs(geneset=split(mcpgenes[,2],mcpgenes[,1]),type="Immu",src="Becht.etal;PMID.27765066",id="IMMU_MCPcounter")%>%
addgs(geneset=ImmuneL,type="Immu",src="Rodrigues.etal;PMID.30179225",id="IMMU_GenJCI121924")%>%


addgs(geneset=STIMgenes,type="CCK",src="MartinSerrano.etal;PMID.35584893",id="CCK_STIM")%>%
addgs(geneset=Biclassgenes,type="CCK",src="Sia.etal;PMID.23295441",id="CCK_Sia13")%>%

addgs(geneset=drugbank,type="Drug",src="DrugBankDec2022",id="DrugBank")%>%


addgs(geneset=list(PSCcaf=scan("ECMsignature_PMID34548310.txt",what="character",sep="\n")),type="ECM",src="Helms.etal;PMID.34548310",id="ECM_Helms22")


signatures=list(geneset=gsignatures$geneset,annotation=as.data.frame(gsignatures[2:4]))

# write_json(toJSON(listSig), "geneSetSignatures.json", pretty = T) 
# .geneSetSignatures=fromJSON(read_json("geneSetSignatures.json",simplifyVector=T))
# usethis::use_data(.geneSetSignatures,internal=T,overwrite=T)

usethis::use_data(signatures,internal=FALSE,overwrite=T)