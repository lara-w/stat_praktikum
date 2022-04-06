library(readxl)
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
### Daten einlesen
fire3 <- readRDS("fire3.Rds")
### EMTGene
EMT_genes<- fire3[,c(12:204)]



multiCox=coxph(Surv(os_time, os_status) ~ ., data = EMT_genes)
multiCoxSum=summary(multiCox)
outMultiTab=data.frame()
outMultiTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outMultiTab=cbind(id=row.names(outMultiTab),outMultiTab)
outMultiTab <- as.data.frame(outMultiTab)
write.table(outMultiTab, file="multiCox.txt", sep="\t", row.names=F, quote=F)
genes <- outMultiTab$id[outMultiTab$pvalue<0.05]
genes
###"ANXA3"    "BHLHE40"  "CA9"      "CEACAM1"  "CTHRC1"   "CXCL8"    "DCN"      "EGLN3"   
###"EMP1"     "EPHB2"    "EPHX2"    "ETV4"     "EZH2"     "GTF3A"    "HAND1"    "KAT2B"   
###"KDM1A"    "KLF4"     "KLF6"     "KRT17"    "LGALS3"   "MARVELD3" "MCM7"     "MMP11"   
####"MSH6"     "MYL9"     "NCS1"     "PDCD4"   

##new_multicox<- outMultiTab[,colnames(outMultiTab)%in%genes]
##write.table(new_multicox, file="newmultiCox.txt", sep="\t", row.names=F, quote=F)
cols_remain<- c("os_status","os_time",genes)


score=predict(multiCox, type="risk", newdata=EMT_genes)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`", "", coxGene)
outCol=c("os_time", "os_status", coxGene)
risk=as.vector(ifelse(score>median(score), "high", "low"))
outTab=cbind(EMT_genes[,outCol], riskScore=as.vector(score), risk)
newcol<- c(cols_remain,"riskScore","risk")
outTab <- outTab[,colnames(outTab)%in%newcol]
write.table(cbind(id=rownames(outTab),outTab), file="risk.txt", sep="\t", quote=F, row.names=F)

###Kaplan-Meier-Schätzer
fit <- survfit(Surv(os_time, os_status) ~ risk, data=outTab)
ggsurvplot(fit, data = outTab,
           pval = TRUE,
           pval.coord = c(0, 0.03), 
           surv.median.line = "hv", 
           palette=c("red", "blue"),  
           legend.labs=c("high-risk","low-risk"), 
           legend.title="groups", 
           ylab="Survival probability",xlab = " Time (Days)",
           censor.shape = 124,censor.size = 2,conf.int = FALSE, 
           break.x.by = 500,
           risk.table = TRUE,tables.height = 0.2,
           tables.theme = theme_cleantable(),
           ggtheme = theme_bw())

