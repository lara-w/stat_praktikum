library(readxl)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(glmnet)
library(carnet)
library(survival)
library(forestplot)

###Daten einlesen
FIRE_3_EMT_genes_final <- read_excel("FIRE-3_EMT_genes_final.xlsx")
FIRE_3_clin_data <- read_excel("FIRE-3_clin_data.xlsx")

###Umbenennen von CRF in pat_nr fuer merge
FIRE_3_EMT_genes_final$pat_nr <- FIRE_3_EMT_genes_final$CRF   

###merge by pat_nr, NA aussortiert
fire3_complete <- merge(FIRE_3_clin_data, FIRE_3_EMT_genes_final, by = "pat_nr")

#fire3_complete <-na.omit(fire3_complete)
str(fire3_complete)

#univariate cox regression
options(forestplot_new_page = FALSE)
clrs <- fpColors(box="green",line="darkblue", summary="royalblue")            
fire3_complete$pfs_time <- fire3_complete$pfs_time/30
fire3_complete <- fire3_complete %>% relocate(pfs_time,pfs_status, .after = pat_nr)
fire3_complete <- fire3_complete[,-(4:15)]
head(fire3_complete)


###HR,Pvalue berechnen 
outTab=data.frame()
for(i in colnames(fire3_complete[,4:ncol(fire3_complete)])){
  cox <- coxph(Surv(pfs_time, pfs_status) ~ fire3_complete[,i], data = fire3_complete)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
outTab
head(outTab)

###Genen mit p < 0,05 auswaehlen 
outTab$id[outTab$pvalue<0.05]
##[1] "CDKN1A" "EIF3B"  "FBLN1"  "MEF2D" 
rows_remain <-c("CDKN1A","EIF3B","FBLN1","MEF2D")
outTab4gene <-subset(outTab, outTab$id %in% c("CDKN1A","EIF3B","FBLN1","MEF2D")) 
write.table(outTab4gene,file="uniCox.xls",sep="\t",row.names=F,quote=F)


###uniCox mit forestplot
rt=read.table("uniCox.xls",header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.05, "<0.05", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          
pdf(file="forest.pdf",
    width = 6,             
    height = 7,            
)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           boxsize=0.4,
           xlab="Hazard ratio"
)
dev.off()


