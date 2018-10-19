rm(list=ls(all=TRUE))
# IMPORTAR DATOS ----------------------------------------------------------
library(dplyr)
library(StatMeasures)
library(xtable)
library(ggplot2)
library(scales)
library(gridExtra)
setwd('/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/OptimalTaxation/EmpiricalMoments')
ENAHO<-read.csv("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/DataAnalysis/All/ENAHO/ENAHOARM/OptimaltaxationSubSampleENAHO.csv", header = T, sep=",")
CENSO<-read.csv("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/DataAnalysis/All/Census/Modified/OptimaltaxationSubSampleCenso.csv", header = T, sep=",")
graphdirectory<-"/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/OptimalTaxation/EmpiricalMoments/Graphs"
# TRIMMING ----------------------------------------------------------------
#A)Variable de interes (CENSO)
#Ventas (Ventas netas + prestacion de servicios)
CENSO$VENTAST<-CENSO$CAP5MONTO1+CENSO$CAP5MONTO5

#B)Sub-muestra
CENSO<-CENSO[which(CENSO$VENTAST>0 & !is.na(CENSO$VENTAST) & !is.na(CENSO$CITAX)),]

#Se selecciona la variable con la que se realiza el trimming.

#CITAX O VENTAS
CENSO<-CENSO[which(!(CENSO$VENTAST>quantile(CENSO$VENTAST,0.99)|CENSO$CITAX>quantile(CENSO$CITAX,0.99))),]

#CITAX
#CENSO<-CENSO[which(CENSO$CITAX<quantile(CENSO$CITAX,0.99)),]

#VENTAS
#CENSO<-CENSO[which(!(CENSO$VENTAST>quantile(CENSO$VENTAST,0.99))),]

#BENEFICIOS
#CENSO<-CENSO[which(CENSO$CAP5MONTO34<quantile(CENSO$CAP5MONTO34,0.99)),]

#C)Quantiles 
#Se define el numero de grupos
nq<-10

# MOMENTO 1 ---------------------------------------------------------------
# Informalidad y numero de trabajadores por tamaño de la empresa (ENAHO).

#Se agrupan las empresas con mas de 50 trabajadores en una misma categoria,
#si se quiere la informacion sin restringir el tamaÃ±o de la empresa se debe 
#usar la variable "numero_trabajadores" en el group by.

ENAHO$informal_empleado<-1-ENAHO$formal_empleado
ENAHO$nn<-1
ENAHO$NUMTRAB<-ENAHO$numero_trabajadores
ENAHO$NUMTRAB[ENAHO$NUMTRAB>50]<-50
informtotal<-round(mean(ENAHO$informal_empleado, na.rm = T)*100, 4)

#Restringiendo la muestra a quienes son empleados
ENAHOT<-ENAHO[which(ENAHO$categopri_ci=="Empleado"),]

MOMENTO1<-ENAHOT%>%
  group_by(NUMTRAB)%>%
  summarise(informalidad=mean(informal_empleado,na.rm = T), trabajadores=sum(nn))

MOMENTO1$informalidad<-MOMENTO1$informalidad*100
names(MOMENTO1)<-c("Tamano de la empresa","Informalidad","Numero de trabajadores")
M1<-as.matrix.data.frame(MOMENTO1,rownames.force = F)
print(xtable(M1,digits = c(0,0,3,0)), include.rownames=F)


# MOMENTO 2  --------------------------------------------------------------
#Numero de trabajadores por decil de ventas (CENSO)

CENSO$TOTALTRAB<-CENSO$NTRABAJADORES8+CENSO$NTRABAJADORES9+CENSO$NTRABAJADORES10
sub<-subset(CENSO, CAP5MONTO1>0)
sub$deciles<-decile(sub$CAP5MONTO1)
MOMENTO2<-sub%>%
  group_by(deciles)%>%
  summarise(trabajadores=sum(TOTALTRAB))
names(MOMENTO2)<-c("Decil de ventas","Numero de trabajadores")
M2<-as.matrix.data.frame(MOMENTO2,rownames.force = F)
print(xtable(M2), include.rownames=F)


# MOMENTO 3 ---------------------------------------------------------------
#Informalidad, numero de firmas y de trabajadores por tamaño de la empresa
#(CENSO+ENAHO)

#Se agrupan las empresas con mas de 50 trabajadores en una misma categoria,
#si se quiere la informacion sin restringir el tamaÃ±o de la empresa se debe 
#usar la variable "TOTALTRAB" en el group by.

CENSO$nn<-1
CENSO$NUMTRAB<-CENSO$TOTALTRAB
CENSO$NUMTRAB[CENSO$NUMTRAB>50]<-50
MM3<-CENSO%>%
  group_by(NUMTRAB)%>%
  summarise(firmas=sum(nn))

totalfirmas<-sum(MM3$firmas)
MM3$proporcion<-MM3$firmas/totalfirmas*100
MM3$propcum<-cumsum(MM3$proporcion)
names(MM3)<-c("Tamano de la empresa","Numero de firmas","% del total de firmas", "% acumulado del total de firmas")

MOMENTO3<-left_join(MM3,MOMENTO1,by=c("Tamano de la empresa","Tamano de la empresa"))
T15<-MOMENTO3

porctotal<-sum(MOMENTO3$`% del total de firmas`)

MOMENTO3$Informalidad<-round(MOMENTO3$Informalidad,3)
MOMENTO3$`% del total de firmas`<-round(MOMENTO3$`% del total de firmas`,3)
MOMENTO3$`% acumulado del total de firmas`<-round(MOMENTO3$`% acumulado del total de firmas`,3)
trabajadorestotal<-sum(MOMENTO3$`Numero de trabajadores`, na.rm = T)

row<-c("Total", totalfirmas,porctotal,"",informtotal,trabajadorestotal)
M3<-rbind(as.matrix.data.frame(MOMENTO3),t(row))

print(xtable(M3, digits =c(0,0,0,3,3,3,0)), include.rownames=F)


# MOMENTO 4 ---------------------------------------------------------------
# Proporcion de trabajadores empleados e independientes en la ENAHO

ENAHO$clasificacion<-NA   
ENAHO$clasificacion[ENAHO$categopri_ci=="Empleado"]<-1
ENAHO$clasificacion[ENAHO$categopri_ci=="Cuenta propia"|ENAHO$categopri_ci=="Patron"]<-2
ENAHO$clasificacion[ENAHO$categopri_ci=="Otro"|ENAHO$categopri_ci=="No_remunerado"]<-3
obsENAHO<-sum(ENAHO$clasificacion==1|ENAHO$clasificacion==2|ENAHO$clasificacion==3, na.rm = T)
ENAHO$clasificacion<-factor(ENAHO$clasificacion, levels = c(1,2,3), labels = c("Empleados","Independientes","Otros + No remunerados"))

M4<-cbind(t(table(ENAHO$clasificacion)),obsENAHO)
`%`<-c(M4[1,1]/M4[1,4]*100,M4[1,2]/M4[1,4]*100,M4[1,3]/M4[1,4]*100,M4[1,4]/M4[1,4]*100)
M4<-rbind(M4,`%`)

print(xtable(M4, digits = matrix(c(rep(0,6),rep(2,3),0),ncol=5, byrow = T)))

# MOMENTO 5 ---------------------------------------------------------------
# Distribucion de los salarios (ENAHO)

#Se define una funcion que arroja el promedio, desviacion estandar y deciles de cualquier
#variable
distribuciones<-function(x){
  Promedio<-round(mean(x,na.rm = T),3)
  SD<-round(sd(x,na.rm = T),3)
  deciles<-round(quantile(x, c(seq(from=0, to=1, by=0.1)), na.rm = T),3)
  dsum<-cbind(c("Promedio","Desviacion estandar","Min (0%)", "10%","20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%","Max (100%)"),as.numeric(c(Promedio,SD,deciles)))
  colnames(dsum)<-c(".","x")
  return(dsum)
}

ENAHO$salarioUSD<-(ENAHO$ylmpri_ci+ ENAHO$ylnmpri_ci + ENAHO$ylmsec_ci + ENAHO$ylnmsec_ci)*0.315
ENAHOT$salarioUSD<-(ENAHOT$ylmpri_ci+ ENAHOT$ylnmpri_ci + ENAHOT$ylmsec_ci + ENAHOT$ylnmsec_ci)*0.315

#Muestra total de la ENAHO
M5.1<-distribuciones(ENAHO$salarioUSD)  
print(xtable(M5.1, digits = matrix(rep(3,39),ncol=3)),include.rownames=F)

M5.2<-distribuciones(ENAHO$salarioUSD[which(ENAHO$salarioUSD>0)])
print(xtable(M5.2, digits = matrix(rep(3,39),ncol=3)),include.rownames=F)

#Sub- muestra: Trabajadores empleados
M5.3A<-distribuciones(ENAHOT$salarioUSD)  

#Trabajadores informales 
M5.3B<-distribuciones(ENAHOT$salarioUSD[ENAHOT$informal_empleado==1])

#Trabajadores formales 
M5.3C<-distribuciones(ENAHOT$salarioUSD[ENAHOT$informal_empleado==0])

M5.3<-rbind(c(".","Todos","Informales","Formales"),cbind(M5.3A,M5.3B[,2],M5.3C[,2])) 
print(xtable(M5.3),include.rownames=F, include.colnames = F)


# MOMENTO 6 ---------------------------------------------------------------
#Distribución de la produccion de las firmas (CENSO)

CENSO$PRODUCCIONUSD<-CENSO$CAP5MONTO9*0.315

M6<-distribuciones(CENSO$PRODUCCIONUSD)
print(xtable(M6),include.rownames=F)


# MOMENTO 7 ---------------------------------------------------------------
#Pago de impuestos y nivel de producción de la firma (CENSO)

CENSO$PROFITSUSD<-CENSO$CAP5MONTO34*0.315

#A)Impuesto a la renta (todos los regimenes)-CITAX

CENSO$TAX_PRODA<-CENSO$CITAXUSD/CENSO$PRODUCCIONUSD*100
CENSO$TAX_PROFA<-CENSO$CITAXUSD/CENSO$PROFITSUSD*100
CENSO$TAX_PROFA[CENSO$PROFITSUSD==0]<-NA

M7A.1<-distribuciones(CENSO$CITAXUSD)
M7A.2<-distribuciones(CENSO$TAX_PRODA)
M7A.3<-distribuciones(CENSO$TAX_PROFA)

M7A<-cbind(M7A.1,M7A.2[,2],M7A.3[,2])
colnames(M7A)<-c(".", "Impuestos","(%)Impuestos/Producción", "(%)Impuestos/Beneficios" )
print(xtable(M7A),include.rownames=F)

CENSO$decilesprod<-decile(CENSO$PRODUCCIONUSD)
MOMENTO7A<-CENSO%>%
  group_by(decilesprod)%>%
  summarise(produccion=sum(PRODUCCIONUSD), beneficios=sum(PROFITSUSD), impuestos=sum(CITAXUSD))

MOMENTO7A$propprod<-MOMENTO7A$impuestos/MOMENTO7A$produccion*100
MOMENTO7A$propprof<-MOMENTO7A$impuestos/MOMENTO7A$beneficios*100
MOMENTO7A<-cbind(MOMENTO7A[,1], MOMENTO7A[,2:4]/1000, MOMENTO7A[,5:6])

names(MOMENTO7A)<-c("Decil de produccion","Produccion","Beneficios","Impuestos","Impuestos/Produccion (%)","Impuestos/Beneficios (%)")
print(xtable(MOMENTO7A),include.rownames=F)


#B)Impuestos totales (CITAX + Transferencias de utilidades a trabajadores)
CENSO$TAXTOTALUSD<-(CENSO$CITAX+CENSO$CAP5MONTO35)*0.315

CENSO$TAX_PRODB<-CENSO$TAXTOTALUSD/CENSO$PRODUCCIONUSD*100
CENSO$TAX_PROFB<-CENSO$TAXTOTALUSD/CENSO$PROFITSUSD*100
CENSO$TAX_PROFB[CENSO$PROFITSUSD==0]<-NA

M7B.1<-distribuciones(CENSO$TAXTOTALUSD)
M7B.2<-distribuciones(CENSO$TAX_PRODB)
M7B.3<-distribuciones(CENSO$TAX_PROFB)

M7B<-cbind(M7B.1,M7B.2[,2],M7B.3[,2])
colnames(M7B)<-c(".", "Impuestos","(%)Impuestos/Producción", "(%)Impuestos/Beneficios" )
print(xtable(M7B),include.rownames=F)


MOMENTO7B<-CENSO%>%
  group_by(decilesprod)%>%
  summarise(produccion=sum(PRODUCCIONUSD), beneficios=sum(PROFITSUSD), impuestos=sum(TAXTOTALUSD))

MOMENTO7B$propprod<-MOMENTO7B$impuestos/MOMENTO7B$produccion*100
MOMENTO7B$propprof<-MOMENTO7B$impuestos/MOMENTO7B$beneficios*100
MOMENTO7B<-cbind(MOMENTO7B[,1], MOMENTO7B[,2:4]/1000, MOMENTO7B[,5:6])

names(MOMENTO7B)<-c("Decil de produccion","Produccion","Beneficios","Impuestos","Impuestos/Produccion (%)","Impuestos/Beneficios (%)")
print(xtable(MOMENTO7B),include.rownames=F)


# MOMENTO 8 ---------------------------------------------------------------
#Nivel de ingresos e informalidad para trabajadores de la ENAHO

ENAHOT$decilesing<-decile(ENAHOT$salarioUSD)

MOMENTO8<-ENAHOT%>%
  group_by(decilesing)%>%
  summarise(informalidad=mean(informal_empleado,na.rm = T)*100, trabajadores=sum(nn))

names(MOMENTO8)<-c("Decil de ingresos" ,"Informalidad (%)", "Observaciones")

print(xtable(MOMENTO8,digits=c(0,0,2,0)),include.rownames=F)


# INFORMACIÓN POR DECILES DE VENTAS-----------------------------------------

CENSO$TOTALTRAB<-CENSO$NTRABAJADORES8+CENSO$NTRABAJADORES9+CENSO$NTRABAJADORES10
CENSO$NUMTRAB<-CENSO$TOTALTRAB
CENSO$NUMTRAB[CENSO$NUMTRAB>50]<-50
CENSO$nn<-1


CENSO$VENTASTUSD<-CENSO$VENTAST*0.315
CENSO$PROFITSBEFOREUSD<-CENSO$CAP5MONTO34*0.315
CENSO$PROFITSAFTERUSD<-CENSO$CAP5MONTO37*0.315
CENSO$TAXTOTALUSD<-(CENSO$CITAX+CENSO$CAP5MONTO35)*0.315

CENSO$decilesventast<-ntile(CENSO$VENTASTUSD,nq)

DECILES_VENTAS<-CENSO%>%
  group_by(decilesventast)%>%
  summarise(Ventas=mean(VENTASTUSD, na.rm = T), BeneficiosAntes=mean(PROFITSBEFOREUSD, na.rm = T), BeneficiosDespues=mean(PROFITSAFTERUSD, na.rm = T), Impuestos1=mean(CITAXUSD, na.rm = T), Impuestos2=mean(TAXTOTALUSD, na.rm = T))

DECILES_VENTAS$TAX1_VENTAS<-DECILES_VENTAS$Impuestos1/DECILES_VENTAS$Ventas
DECILES_VENTAS$TAX2_VENTAS<-DECILES_VENTAS$Impuestos2/DECILES_VENTAS$Ventas

DECILES_VENTAS$TAX1_BENEFICIOS<-DECILES_VENTAS$Impuestos1/DECILES_VENTAS$BeneficiosDespues
DECILES_VENTAS$TAX2_BENEFICIOS<-DECILES_VENTAS$Impuestos2/DECILES_VENTAS$BeneficiosDespues

DECILES_VENTAS$Ventas_decil<-quantile(CENSO$VENTASTUSD, c(seq(1/nq,1,1/nq)))

# INFORMACION POR DECILES DE BENEFICIOS -----------------------------------

CENSO$decilesbeneficios<-ntile(CENSO$PROFITSBEFOREUSD,nq)

DECILES_BENEFICIOS<-CENSO%>%
  group_by(decilesbeneficios)%>%
  summarise(Ventas=mean(VENTASTUSD, na.rm = T), BeneficiosAntes=mean(PROFITSBEFOREUSD, na.rm = T), BeneficiosDespues=mean(PROFITSAFTERUSD, na.rm = T), Impuestos1=mean(CITAXUSD, na.rm = T), Impuestos2=mean(TAXTOTALUSD, na.rm = T))

DECILES_BENEFICIOS$TAX1_VENTAS<-DECILES_BENEFICIOS$Impuestos1/DECILES_BENEFICIOS$Ventas
DECILES_BENEFICIOS$TAX2_VENTAS<-DECILES_BENEFICIOS$Impuestos2/DECILES_BENEFICIOS$Ventas

DECILES_BENEFICIOS$TAX1_BENEFICIOS<-DECILES_BENEFICIOS$Impuestos1/DECILES_BENEFICIOS$BeneficiosDespues
DECILES_BENEFICIOS$TAX2_BENEFICIOS<-DECILES_BENEFICIOS$Impuestos2/DECILES_BENEFICIOS$BeneficiosDespues

DECILES_BENEFICIOS$Beneficios_decil<-quantile(CENSO$PROFITSBEFOREUSD, c(seq(1/nq,1,1/nq)))

# GRAFICA 1 ---------------------------------------------------------------
#Proporción de trabajadores informales por deciles de ventas de las firmas
#“InformalProportionDemandPercentile”

TS<-CENSO%>%
  group_by(NUMTRAB)%>%
  summarise(firmas=sum(nn), ventas=mean(VENTASTUSD, na.rm = T))
TS$decil<- ntile(TS$ventas,nq)

TS<-left_join(TS,MOMENTO1,by=c("NUMTRAB"="Tamano de la empresa"))
TS$`Numero de trabajadores`<- TS$firmas*TS$NUMTRAB
TS$Informales<-TS$firmas*TS$NUMTRAB*TS$Informalidad/100
TS$Formales<-TS$`Numero de trabajadores`-TS$Informales

INF<-TS%>%
  group_by(decil)%>%
  summarise(Total=sum(`Numero de trabajadores`), Informales=sum(Informales,na.rm = T), Formales=sum(Formales, na.rm = T))

INF$Informalidad<-INF$Informales/INF$Total
INF$Formalidad<-1-INF$Informalidad
INF$Percentil<-c(seq(1/nq,1,1/nq))


G1<-ggplot(data=INF,aes(x=Percentil,y=Informalidad))+
  geom_line()+
  labs(x=expression(theta~e),y="Proportion Informal Demand")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))

setwd(graphdirectory)
dev.set()
png(file="InformalProportionDemandPercentileE.png",width=1600,height=850)
G1
dev.off()


# GRAFICA 2 ---------------------------------------------------------------
#Proporción de trabajadores formales por deciles de ventas de las firmas
#“FormalProportionDemandPercentile”

G2<-ggplot(data=INF,aes(x=Percentil,y=Formalidad))+
  geom_line()+
  labs(x=expression(theta~e),y="Proportion Formal Demand")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))

dev.set()
png(file="FormalProportionDemandPercentileE.png",width=1600,height=850)
G2
dev.off()


# GRAFICA 3 ---------------------------------------------------------------
#Número de trabajadores totales (formales e informales) por deciles de ventas de las firmas
#"TotalLaborDemandPercentile"

T3<-CENSO%>%
  group_by(decilesventast)%>%
  summarise(trabajadores=mean(TOTALTRAB))

T3$Percentil<-c(seq(1/nq,1,1/nq))

G3<-ggplot(data=T3,aes(x=Percentil,y=trabajadores))+
  geom_line()+
  labs(x=expression(theta~e),y="Total Labor Demand (Average number of workers)")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))

dev.set()
png(file="TotalLaborDemandE.png",width=1600,height=850)
G3
dev.off()



# GRAFICA 4 ---------------------------------------------------------------
#Niveles de beneficios antes de impuestos por deciles de ventas.
#“PretaxProfitPercentile”.

DECILES_VENTAS$Percentil<-c(seq(1/nq,1,1/nq))

G4<-ggplot(data=DECILES_VENTAS,aes(x=Percentil,y=BeneficiosAntes))+
  geom_line()+
  labs(x=expression(theta~e),y="Profits Distribution (USD dollars)")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))

dev.set()
png(file="PretaxProfitPercentileE.png",width=1600,height=850)
G4
dev.off()


# GRAFICA 5 ---------------------------------------------------------------
#Niveles de beneficios después de impuestos por deciles de ventas.
#“AfterTaxProfitPercentile”

G5<-ggplot(data=DECILES_VENTAS,aes(x=Percentil,y=BeneficiosDespues))+
  geom_line()+
  labs(x=expression(theta~e),y="After Tax Profits Distribution (USD dollars)")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))


dev.set()
png(file="AfterTaxProfitPercentileE.png",width=1600,height=850)
G5
dev.off()

# GRAFICA 6 ---------------------------------------------------------------
#Niveles de venta (valor) en cada decil de ventas.
#“ProductionPercentile”

DECILES_VENTAS$VM<-DECILES_VENTAS$Ventas_decil/1000

G6<-ggplot(data=DECILES_VENTAS,aes(x=Percentil,y=VM))+
  geom_line()+
  labs(x=expression(theta~e),y="Production Distribution (Thousands of USD dollars)")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))

dev.set()
png(file="ProductionPercentileE.png",width=1600,height=850)
G6
dev.off()


# GRAFICA 7 ---------------------------------------------------------------
#Valor de pago de impuesto en cada decil de ventas (CITAX)
#Valor de pago de impuesto en cada decil de ventas (CITAX)
G7.1<-ggplot(data=DECILES_VENTAS,aes(x=Percentil,y=Impuestos1))+
  geom_line()+
  labs(x=expression(theta~e),y="Taxes Payed (CITAX - USD dollars)")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))

dev.set()
png(file="Taxes1E.png",width=1600,height=850)
G7.1
dev.off()

#Valor de pago de impuesto en cada decil de ventas (CITAX+Transferencias)
G7.2<-ggplot(data=DECILES_VENTAS,aes(x=Percentil,y=Impuestos2))+
  geom_line()+
  labs(x=expression(theta~e),y="Taxes Payed (CITAX+Transferencias - USD dollars)")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))

dev.set()
png(file="Taxes2E.png",width=1600,height=850)
G7.2
dev.off()


# GRAFICA 8 ---------------------------------------------------------------
#Valor de pago de impuesto como proporción de ventas, en cada decil de ventas (CITAX)
G8.1<-ggplot(data=DECILES_VENTAS,aes(x=Percentil,y=TAX1_VENTAS))+
  geom_line()+
  labs(x=expression(theta~e),y="Taxes Payed (CITAX) / Production")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))
dev.set()
G8.1
dev.off()

#Valor de pago de impuesto como proporciÃÂ³n de ventas, en cada decil de ventas (CITAX+TRANSFERENCIAS)
G8.2<-ggplot(data=DECILES_VENTAS,aes(x=Percentil,y=TAX2_VENTAS))+
  geom_line()+
  labs(x=expression(theta~e),y="Taxes Payed (CITAX+Transferencias) / Production")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))

dev.set()
png(file="TaxesSales2E.png",width=1600,height=850)
G8.2
dev.off()


# GRAFICA 9 ---------------------------------------------------------------
#Valor de pago de impuesto como proporción de los beneficios, en cada decil de ventas
G9.1<-ggplot(data=DECILES_VENTAS,aes(x=Percentil,y=TAX1_BENEFICIOS))+
  geom_line()+
  labs(x=expression(theta~e),y="Taxes Payed (CITAX)  / Profits")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))

dev.set()
png(file="TaxesProfits1E.png",width=1600,height=850)
G9.1
dev.off()

#Valor de pago de impuesto en cada decil de ventas (CITAX+Transferencias)
G9.2<-ggplot(data=DECILES_VENTAS,aes(x=Percentil,y=TAX2_BENEFICIOS))+
  geom_line()+
  labs(x=expression(theta~e),y="Taxes Payed  (CITAX+Transferencias)  / Profits")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))
dev.set()
png(file="TaxesProfits2E.png",width=1600,height=850)
G9.2
dev.off()


# GRAFICA 10 --------------------------------------------------------------
#Distribución de ingresos de los trabajadores formales

mins<-min(ENAHOT$salarioUSD)
maxs<-max(ENAHOT$salarioUSD)
ENAHOT$salarionorm<-(ENAHOT$salarioUSD-mins)/(maxs-mins)


G10<-ggplot(ENAHOT[which(ENAHOT$informal_empleado==0),],aes(x=salarionorm))+
  geom_histogram(aes(y=..ncount..))+
  geom_density(aes(y=..scaled..))+
  labs(x="Wage (Formal workers)", y="Density")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))

dev.set()
png(file="FormalSupplyPercentiles.png",width=1600,height=850)
G10
dev.off()


# GRAFICA 11 --------------------------------------------------------------
#Distribución de ingresos de los trabajadores informales

G11<-ggplot(ENAHOT[which(ENAHOT$informal_empleado==1),],aes(x=salarionorm))+
  geom_histogram(aes(y=..ncount..))+
  geom_density(aes(y=..scaled..))+
  labs(x="Wage (Informal workers)", y="Density")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))
dev.set()
png(file="InformalSupplyPercentiles.png",width=1600,height=850)
G11
dev.off()

# GRAFICA 12 --------------------------------------------------------------
#Distribución de ingresos de los trabajadores (formales e informales)

G12<-ggplot(ENAHOT,aes(x=salarionorm))+
  geom_histogram(aes(y=..ncount..))+
  geom_density(aes(y=..scaled..))+
  labs(x="Wage (Total labor supply)", y="Density")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))

dev.set()
png(file="TotalLaborSupplyPercentiles.png",width=1600,height=850)
G12
dev.off()

# GRAFICA 13 --------------------------------------------------------------
#Distribución de trabajadores informales como proporción de todos los trabajadores organizados por deciles de ingresos.

ENAHOT$percentiling<-ntile(ENAHOT$salarioUSD,nq)

T13<-ENAHOT%>%
  group_by(percentiling)%>%
  summarise(informalidad=mean(informal_empleado,na.rm = T)*100, trabajadores=sum(nn))
T13$Percentil<-c(seq(1/nq,1,1/nq))

G13<-ggplot(data=T13,aes(x=Percentil,y=informalidad/100))+
  geom_line()+
  labs(x=expression(theta~w),y="Proportion Informal Supply")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))


dev.set()
png(file="InformalProportionSupplyPercentiles.png",width=1600,height=850)
G13
dev.off()

# GRAFICA 14 --------------------------------------------------------------
#Proporción de informales por tamaño de firma

G14<-ggplot(data=MOMENTO1,aes(x=`Tamano de la empresa`,y=Informalidad/100))+
  geom_line()+
  labs(x="Firm Size (Number of workers)",y="Proportion Informal Supply")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))

dev.set()
png(file="InformalFirmSize.png",width=1600,height=850)
G14
dev.off()

# GRAFICA 15 --------------------------------------------------------------
#Informalidad por percentiles de tamaño de la firma

G15<-ggplot(data=T15,aes(x=`% acumulado del total de firmas`/100,y=Informalidad/100))+
  geom_line()+
  labs(x="Firm Size Percentile",y="Proportion Informal Supply")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))

dev.set()
png(file="InformalFirmSizePercentiles.png",width=1600,height=850)
G15
dev.off()

# GRAFICA 16 --------------------------------------------------------------
#Niveles de beneficios (valor) en cada decil de beneficios.


DECILES_BENEFICIOS$Percentil<-c(seq(1/nq,1,1/nq))
DECILES_BENEFICIOS$BM<-DECILES_BENEFICIOS$Beneficios_decil/1000

G16<-ggplot(data=DECILES_BENEFICIOS,aes(x=Percentil,y=BM))+
  geom_line()+
  labs(x=expression("Profits Percentile"),y="Profits Distribution (Thousands of USD dollars)")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))

dev.set()
png(file="ProfitsPercentileE.png",width=1600,height=850)
G16
dev.off()


# GRAFICA 17 --------------------------------------------------------------
#Valor de pago de impuesto como proporción de ventas, en cada decil de beneficios (CITAX)
G17.1<-ggplot(data=DECILES_BENEFICIOS,aes(x=Percentil,y=TAX1_VENTAS))+
  geom_line()+
  labs(x="Profits Percentile",y="Taxes Payed (CITAX) / Production")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))
dev.set()
png(file="PTaxesSales1E.png",width=1600,height=850)
G17.1
dev.off()

#Valor de pago de impuesto como proporción de ventas, en cada decil de beneficios (CITAX+TRANSFERENCIAS)
G17.2<-ggplot(data=DECILES_BENEFICIOS,aes(x=Percentil,y=TAX2_VENTAS))+
  geom_line()+
  labs(x="Profits Percentile",y="Taxes Payed (CITAX+Transferencias) / Production")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))

dev.set()
png(file="PTaxesSales2E.png",width=1600,height=850)
G17.2
dev.off()

# GRAFICA 18 ---------------------------------------------------------------
#Valor de pago de impuesto como proporción de los beneficios, en cada decil de beneficios (CITAX)
G18.1<-ggplot(data=DECILES_BENEFICIOS,aes(x=Percentil,y=TAX1_BENEFICIOS))+
  geom_line()+
  labs(x="Profits Percentile",y="Taxes Payed (CITAX)  / Profits")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))

dev.set()
png(file="PTaxesProfits1E.png",width=1600,height=850)
G18.1
dev.off()

#Valor de pago de impuesto como proporción de los beneficios, en cada decil de beneficios (CITAX+Transferencias)
G18.2<-ggplot(data=DECILES_BENEFICIOS,aes(x=Percentil,y=TAX2_BENEFICIOS))+
  geom_line()+
  labs(x="Profits Percentile",y="Taxes Payed  (CITAX+Transferencias)  / Profits")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))
dev.set()
png(file="PTaxesProfits2E.png",width=1600,height=850)
G18.2
dev.off()

# MOMENTO 19 --------------------------------------------------------------
#Horas promedio (semanales) de participación en el mercado laboral por decil de ingresos
#Version 1 (Trabajadores)
T19<-ENAHOT[which(ENAHOT$horastot_ci>0),]%>%
  group_by(percentiling)%>%
  summarise(horas=mean(horastot_ci,na.rm = T))

T19$Percentil<-c(seq(1/nq,1,1/nq))

G19<-ggplot(data=T19,aes(x=Percentil,y=horas))+
  geom_line()+
  labs(x=expression(theta~w),y="Number of hours worked per week")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))

dev.set()
png(file="HoursWorked.png",width=1600,height=850)
G19
dev.off()

#Versión 2 (PEA (incluye a quienes no trabajan) // Nivel de ingresos del hogar)
ENAHO$ing_ch<- ENAHO$ylm_ch + ENAHO$ylnm_ch
ENAHO$ing_ch_pc<-ENAHO$ing_ch/ENAHO$nmiembros_ch

ENAHO$horas_trab<- ENAHO$horastot_ci
ENAHO$horas_trab[is.na(ENAHO$horas_trab)]<-0

ENAHOPEA<-ENAHO[which(ENAHO$condocup_ci=="Ocupado"|ENAHO$condocup_ci=="Desocupado"),]
ENAHOPEA$percentil_hogar<-ntile(ENAHOPEA$ing_ch_pc,nq)
T19.2<-ENAHOPEA[!is.na(ENAHOPEA$ing_ch_pc),]%>%
  group_by(percentil_hogar)%>%
  summarise(horas=mean(horas_trab,na.rm = T))

T19.2$Percentil<-c(seq(1/nq,1,1/nq))

G19.2<-ggplot(data=T19.2,aes(x=Percentil,y=horas))+
  geom_line()+
  labs(x=expression(theta~w),y="Number of hours worked per week")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))

dev.set()
png(file="HoursWorkedPEA.png",width=1600,height=850)
G19.2
dev.off()

# MOMENTO 20 --------------------------------------------------------------
#Participación en el mercado laboral por decil de ingresos

ENAHO$ingreso_totUSD<-ENAHO$salarioUSD+(ENAHO$ynlm_ci+ENAHO$ynlnm_ci)*0.315
ENAHO$perc_ingt<-ntile(ENAHO$ingreso_totUSD,nq) 

T20<-ENAHO%>%
  group_by(perc_ingt)%>%
  summarise(ocupados=sum(emp_ci), desocupados=sum(desemp_ci))

T20$Percentil<-c(seq(1/nq,1,1/nq))
T20$Participacion<-T20$ocupados/(T20$ocupados+T20$desocupados)

G20<-ggplot(data=T20,aes(x=Percentil,y=Participacion))+
  geom_line()+
  labs(x=expression(theta~w),y="Participation in labor market")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))

dev.set()
png(file="ParticipationLaborMarket.png",width=1600,height=850)
G20
dev.off()
