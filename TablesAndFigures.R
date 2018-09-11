rm(list=ls(all=TRUE))
library(xtable)
#1. 

setwd('C:/Users/razuero/Dropbox/OptmalTaxationShared/Informality/TablesandFigures')


r1<-c("$\\beta$",0.2135)
r2<-c("$\\chi$",2.0192)
r3<-c("$\\gamma$",0.7341)
r4<-c("$\\delta$",0.12873)
r5<-c("$\\kappa$",0.1021)
r6<-c("$\\psi$",0.4528)
r7<-c("$\\rho$",0.0912)
r8<-c("$\\sigma$",0.1827)
r9<-c("$\\mu_e$",1.2528)
r10<-c("$\\mu_w$",1.7626)
r11<-c("$\\sigma_w^2$",1.0921)
r12<-c("$\\sigma_e^2$",1.1675)
r13<-c("$\\sigma_{w,e}$",0.2782)


t<-rbind(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13)

#Anualizar transferencias del gobierno

colnames(t)<-c("Parameter", "Estimate")
rownames(t)<-NULL
t<-as.data.frame(t)
rownames(t)<-NULL
write(print(xtable(t,align="ccc",caption="Calibration results",table.placement="",label="tab:CalibrationResults"),table.placement="H",caption.placement="top", sanitize.text.function=function(x){x},include.rownames = FALSE),file="CalibrationResults.tex")



write(print(xtable(t,align="cccccc",caption="Calibration results",table.placement="",label="tab:CalibrationResults"),table.placement="H",caption.placement="top", sanitize.text.function=function(x){x}),file="CalibrationResults.tex")

write(print(xtable(TT,align="cccccc",caption="Calibration results",table.placement="",label="tab:CalibrationResults"),table.placement="H",caption.placement="top", sanitize.text.function=function(x){x}),file="CalibrationResults.tex")



