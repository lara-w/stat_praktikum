library(readxl)
library(ggplot2)
library(glmnet)
library(survival)
library(survminer)

###Daten einlesen
FIRE_3_EMT_genes_final <- read_excel("FIRE-3_EMT_genes_final.xlsx")
FIRE_3_clin_data <- read_excel("FIRE-3_clin_data.xlsx")

### Umbenennen von CRF in pat_nr fuer merge
FIRE_3_EMT_genes_final$pat_nr <- FIRE_3_EMT_genes_final$CRF   

### merge by pat_nr, NA aussortiert
fire3_complete <- merge(FIRE_3_clin_data, FIRE_3_EMT_genes_final, by = "pat_nr")

###fire3mit4Genes 
fire3_Genes <- read_excel("fire3_Genes.xlsx")
cols_remain <- c("pat_nr", "pfs_time","pfs_status","CDKN1A","EIF3B","FBLN1","MEF2D")
rt <- fire3_complete[,colnames(fire3_complete)%in% cols_remain]

###multivariate Cox
multiCox=coxph(Surv(pfs_time, pfs_status) ~ ., data = rt)
multiCox=step(multiCox, direction="both")
multiCoxSum=summary(multiCox)

###Coeffizienten,HR und p-value berechnen
outMultiTab=data.frame()
outMultiTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outMultiTab=cbind(id=row.names(outMultiTab),outMultiTab)
write.table(outMultiTab, file="multiCox.txt", sep="\t", row.names=F, quote=F)

#risk score
score=predict(multiCox, type="risk", newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`", "", coxGene)
outCol=c("pfs_time", "pfs_status", coxGene)
risk=as.vector(ifelse(score>median(score), "high", "low"))
outTab=cbind(rt[,outCol], riskScore=as.vector(score), risk)
write.table(cbind(id=rownames(outTab),outTab), file="risk.txt", sep="\t", quote=F, row.names=F)


bioForest=function(coxFile=null, forestFile=null, forestCol=null){
  ###data Einlesen
  rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  #plot 
  pdf(file=forestFile, width=9, height=7)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #forestplot linke Seite
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.05-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
  
  ##forestplot rechte Seite
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol[1], forestCol[2])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.6)
  axis(1)
  dev.off()
}
#forest plot multicox
bioForest(coxFile="multiCox.txt", forestFile="model.multiForest.pdf", forestCol=c("red","green"))
