
#' Title Simple ICA projection, default using Puleo et al components
#'
#' @param newexp A new expression dataset, should be gene symbol in rows
#' @param ICAgw ICA gene weights (S matrix)
#' @param geneNormType Normalization of gene expression matrix (sc: sample scale is default)
#' @param projNormType Normalization of components (raw:  is default)
#' @param ming Minimum number of overlapping genes
#'
#' @return projected components
#' @keywords internal
#'
#' @examples
.qProjICA=function(newexp,ICAgw=CancerRNASig:::puleoICAgw,geneNormType="sc",projNormType="raw",ming=500){

  comg = intersect(rownames(newexp), rownames(ICAgw))

  if(length(comg)<ming){
    stop("Too few genes in common in the expression dataset and the ICA weights")
  }

  scexp = CancerRNASig:::.qNormalize(newexp[comg, ],type=geneNormType)

  invs = MASS::ginv(as.matrix(ICAgw[comg,]))

  proj=t(CancerRNASig:::.qNormalize(invs %*%  scexp,projNormType))

  colnames(proj)=colnames(ICAgw)
  proj

}


.qNumerote=function(n,l=max(nchar(1:n))) {
  str_pad(1:n,width=l,side="left",pad=0)
}

.qNormalize=function(x,type){
  switch(type,
         raw={x},
         gsc={t(scale(t(x)))},
         gc={t(scale(t(x),scale=F))},
         sc={(scale((x),scale=F))},
         ssc={(scale((x)))}
         ,x
  )
}



.qICA=function(X,k,maxiter = 10^6,eps = 10^-6){
  resICA=NULL
  try({
    resICA <- JADE::JADE(X, n.comp = k, maxiter = maxiter, eps = eps)
    rownames(resICA$A) <- colnames(X)
    rownames(resICA$S) <- rownames(X)
    colnames(resICA$S)=paste("ICA",CancerRNASig:::.qNumerote(ncol(resICA$S)),sep="")
    colnames(resICA$A)=paste("ICA",CancerRNASig:::.qNumerote(ncol(resICA$A)),sep="")
  })
  resICA
}






.qProjICA.ds=function(icarez,dataset,geneNormType="sc",projNormType="raw"){

  # expg=getUniqueGeneMat(dataset$exp,dataset$probeannot[rownames(dataset$exp),dataset$genecol],rowSds(as.matrix(dataset$exp)))
   comg = intersect(rownames(icarez$S), rownames(dataset$exp))
  scexp = CancerRNASig:::.qNormalize(dataset$exp[comg, ],type=geneNormType)


  invs = MASS::ginv(as.matrix(icarez$S[comg,]))

  proj=t(CancerRNASig:::.qNormalize(invs %*%  scexp,projNormType))

  # proj = scale(t(t((t(scexp) %*% t(invs)))))

  colnames(proj)=colnames(icarez$S)
  proj

}



# rotation test
if(F){


  icas=readRDS("/Users/remy.nicolle/Downloads/ica_components.rds")

  edat=readRDS("/Users/remy.nicolle/Downloads/expression_matrix.rds")

  pc <- prcomp(t(edat),center = F,scale. = F)
  sapply(1:10,\(i)summary(lm(icas[,i]~pc$rotation))$r.squared)

  pc0 <- prcomp(t(edat),center = F,scale. = F)
  pc1 <- prcomp(t(edat[sample(1:nrow(edat),100),]),center = F,scale. = F)
  sapply(1:10,\(i)summary(lm(pc1$rotation[,i]~pc0$rotation[rownames(pc1$rotation),]))$r.squared)

}

# optimize txBenefit
if(F){
  devtools::install_github("msadatsafavi/txBenefit")
  library(txBenefit)
  library(survival)
  library(survminer)

  #simulate
  H0=T
  n=50 ; pevent=0.6 ; reff=2
  sig0=runif(n*2)
  time=c( (sig0[1:n]*reff +runif(n)),runif(n)*reff)
  event=1*(runif(n*2)> pevent)
  
  tx=rep(c(1,0),each=n)
  if(H0){tx=runif(n*2)>0.5}
  
  datcox=data.frame(time,event,tx,sig0)
  dat0cox=data.frame(time,event,tx=0,sig0)
  dat1cox=data.frame(time,event,tx=1,sig0)




  fit=coxph(Surv(time=time,event=event) ~ tx +sig0+ tx:sig0  , data=datcox, model=TRUE)
  summary(fit)$coefficients[3,]
  pt=predict(fit,type="expected")
  p0=predict(fit,newdata=dat0cox,type="expected")
  p1=predict(fit,newdata=dat1cox,type="expected")
  pdiff=0-p1
  boxplot(pdiff~datcox$tx)
  wilcox.test(p0-p1~datcox$tx)

  # ggsurvplot(survfit(Surv(time, event) ~ tx, data = data.frame(datcox,pgp=sig0>median(sig0))))  
  ggsurvplot(survfit(Surv(time, event) ~ tx+pgp, data = data.frame(datcox,pgp=sig0>median(sig0))))

  Cb.cox(coxph(Surv(time=time,event=event) ~ tx + tx:sig0 +sig0 , data=datcox, model=TRUE),tx_var = "tx")
  Cb.cox(coxph(Surv(time=time,event=event) ~ tx + tx:sig0  , data=datcox, model=TRUE),tx_var = "tx")
  Cb.cox(coxph(Surv(time=time,event=event) ~ tx:sig0  , data=datcox, model=TRUE),tx_var = "tx")




cor(pdiff[datcox$tx==0],datcox$sig0[datcox$tx==0],method = "spearman")
cor(pdiff[datcox$tx==1],datcox$sig0[datcox$tx==1],method = "spearman")

plot(p1,p0,col=datcox$tx+1)
plot(pt,p0,col=datcox$tx+1)
plot(pt,p1,col=datcox$tx+1)
plot(pt,p0-p1,col=datcox$tx+1)


predict(coxph(Surv(time=time,event=event) ~ tx + tx:sig0 +sig0 , data=datcox, model=TRUE))


fittxbenef=function(survdf,vardf,init=rnorm(ncol(vardf))){
  if(!is.data.frame(survdf)) stop("survdf should be a data frame")
  if(!all(c("time","event","trt") %in% colnames(survdf))) stop("survdf should have columns time , event and trt")
  if(length(init)!=ncol(vardf)) stop("init should have the same length as the number of columns of vardf")
  if(nrow(vardf)!=nrow(survdf)) stop("vardf and survdf should have the same number of rows")
  evaltxbenefit=function(w,datvoi,datcox){
    datcox$newvar=as.matrix(datvoi)%*%t(t(w))
    # -Cb.cox(coxph(Surv(time=tte,event=event) ~ tx + tx:newvar +newvar , data=datcox, model=TRUE),tx_var = "tx")$Cb
    summary(coxph(Surv(time=time,event=event) ~ trt + trt:newvar +newvar , data=datcox))$coefficients[3,4]
  }
  res=optim(init,fn=evaltxbenefit,datvoi=vardf,datcox=survdf)
  list(w=res$par,
  z=res$value ,
  newsig=tcrossprod(res$par,vardf)[1,]
  )
}


  H0=F
  n=100 ; pevent=0.6 ;reff=1
  k=4
  X=matrix(runif(n*k*2),ncol=k)
  w=runif(k)
  sig0=tcrossprod(w,X)[1,]
  timetx1=(sig0[1:n]/reff) +runif(n)
  
  time=c(timetx1, runif(n,min=min(timetx1),max=max(timetx1)))
  event=1*(runif(n*2)> pevent)
  
  tx=rep(c(1,0),each=n)
  if(H0){tx=runif(n*2)>0.5}
  survdat=data.frame(time,event,trt=tx)


    # ggsurvplot(survfit(Surv(time, event) ~ tx, data =survdat))

  customfit=fittxbenef(survdat,X)  
  
  fitsig=tcrossprod(customfit$w,X)[1,]
  cor(fitsig,sig0)
  customfit$z

  tesdat=data.frame(survdat,sig=sig0,sugcut=sig0>median(sig0),trainsig=fitsig,tsigcut=fitsig>median(fitsig))
  ggsurvplot(survfit(Surv(time, event) ~ tx+sugcut, data =tesdat))

  (coxph(Surv(time=time,event=event) ~ tx + tx:sig +sig , data=tesdat))
  (coxph(Surv(time=time,event=event) ~ tx + tx:trainsig +trainsig , data=tesdat))

  ggsurvplot(survfit(Surv(time, event) ~ tx+tsigcut, data =tesdat))

}