##data einfuegen 
fire3 <- readRDS("fire3.Rds")
### Clin_Daten
clin_data <- fire3[1:13]
str(clin_data)

cols<- colnames(clin_data)
cols
new_cols<- c(cols[1],cols[12:13],cols[5:(length(cols)-2)])
clin_data<- clin_data[,new_cols]


####multiCox
str(clin_data)
res.cox <- coxph(Surv(os_time, os_status)~.,data =  clin_data)
res.cox

###Kaplan-Meier-Schätzer
km1<-survfit(Surv(os_time,os_status)~Liver_limited_disease,data= clin_data)

ggsurvplot(km1, data = clin_data,
           pval = TRUE,
           pval.coord = c(0, 0.03), 
           surv.median.line = "hv", 
           palette=c("red", "blue"),  
           legend.labs=c("no","yes"), 
           legend.title="Liver_limited_disease", 
           title="overall survival", 
           ylab="Survival probability",xlab = " Time (Days)",
           censor.shape = 124,censor.size = 2,conf.int = FALSE, 
           break.x.by = 500,
           risk.table = TRUE,tables.height = 0.2,
           tables.theme = theme_cleantable(),
           ggtheme = theme_bw())


km2<-survfit(Surv(os_time,os_status)~sex,data = clin_data)
ggsurvplot(km2, data = clin_data,
           pval = TRUE,
           pval.coord = c(0, 0.03), 
           surv.median.line = "hv", 
           palette=c("red", "blue"),  
           legend.labs=c("male","female"), 
           legend.title="Sex", 
           title="overall survival", 
           ylab="Survival probability",xlab = " Time (Days)",
           censor.shape = 124,censor.size = 2,conf.int = FALSE, 
           break.x.by = 500,
           risk.table = TRUE,tables.height = 0.2,
           tables.theme = theme_cleantable(),
           ggtheme = theme_bw())

km3<-survfit(Surv(os_time,os_status)~primeloc_side,data = clin_data)
ggsurvplot(km3, data = clin_data,
           pval = TRUE,
           pval.coord = c(0, 0.03), 
           surv.median.line = "hv", 
           palette=c("red", "blue","purple"),  
           legend.labs=c("left","right","both"), 
           legend.title="Primary_tumor_sidedness", 
           title="overall survival", 
           ylab="Survival probability",xlab = " Time (Days)",
           censor.shape = 124,censor.size = 2,conf.int = FALSE, 
           break.x.by = 500,
           risk.table = TRUE,tables.height = 0.2,
           tables.theme = theme_cleantable(),
           ggtheme = theme_bw())

km4<-survfit(Surv(os_time,os_status)~braf_wild,data = clin_data)
ggsurvplot(km4, data = clin_data,
           pval = TRUE,
           pval.coord = c(0, 0.03), 
           surv.median.line = "hv", 
           palette=c("red", "blue"),  
           legend.labs=c("Wildtype","nicht getestet/mutiert "), 
           legend.title="Braf", 
           title="overall survival", 
           ylab="Survival probability",xlab = " Time (Days)",
           censor.shape = 124,censor.size = 2,conf.int = FALSE, 
           break.x.by = 500,
           risk.table = TRUE,tables.height = 0.2,
           tables.theme = theme_cleantable(),
           ggtheme = theme_bw())


km5<-survfit(Surv(os_time,os_status)~braf_mut,data = clin_data)
ggsurvplot(km5, data = clin_data,
           pval = TRUE,
           pval.coord = c(0, 0.03), 
           surv.median.line = "hv", 
           palette=c("red", "blue"),  
           legend.labs=c("nicht getestet/Wildtype","mutiert "), 
           legend.title="Braf", 
           title="overall survival", 
           ylab="Survival probability",xlab = " Time (Days)",
           censor.shape = 124,censor.size = 2,conf.int = FALSE, 
           break.x.by = 500,
           risk.table = TRUE,tables.height = 0.2,
           tables.theme = theme_cleantable(),
           ggtheme = theme_bw())
