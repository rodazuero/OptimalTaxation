rm(list=ls(all=TRUE))
# IMPORTAR DATOS ----------------------------------------------------------
library(dplyr)
library(StatMeasures)
library(xtable)
library(ggplot2)
library(scales)
library(gridExtra)
setwd('/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/EmpiricalMoments')
setwd("C:/Users/mr.porras10/OneDrive - Universidad de Los Andes/W/New/Graficas")
ENAHO<-read.csv("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/DataAnalysis/All/ENAHO/ENAHOARM/OptimaltaxationSubSampleENAHO.csv", header = T, sep=",")
CENSO<-read.csv("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/DataAnalysis/All/Census/Modified/OptimaltaxationSubSampleCenso.csv", header = T, sep=",")
CENSO<-CENSO[which(CENSO$CAP5MONTO1>0),]
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

porctotal<-sum(MOMENTO3$`% del total de firmas`)

MOMENTO3$Informalidad<-round(MOMENTO3$Informalidad,3)
MOMENTO3$`% del total de firmas`<-round(MOMENTO3$`% del total de firmas`,3)
MOMENTO3$`% acumulado del total de firmas`<-round(MOMENTO3$`% acumulado del total de firmas`,3)
trabajadorestotal<-sum(MOMENTO3$`Numero de trabajadores`, na.rm = T)

row<-c("Total", totalfirmas,porctotal,"",informtotal,trabajadorestotal)
M3<-rbind(as.matrix.data.frame(MOMENTO3),t(row))

print(xtable(M3, digits =c(0,0,0,3,3,3,0)), include.rownames=F)


#ggplot(MOMENTO3, aes(x=`Tamano de la empresa`))+geom_bar(stat = "identity",aes(y=`% del total de firmas`))+
#geom_line(aes(y= Informalidad))+  scale_y_continuous(sec.axis = sec_axis(~., name = "Informalidad"))


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

CENSO$VENTASUSD<-CENSO$CAP5MONTO1*0.315
CENSO$PROFITSBEFOREUSD<-CENSO$CAP5MONTO34*0.315
CENSO$PROFITSAFTERUSD<-CENSO$CAP5MONTO37*0.315

CENSO$decilesventas<-decile(CENSO$CAP5MONTO1)

DECILES_VENTAS<-CENSO%>%
  group_by(decilesventas)%>%
  summarise(Ventas=sum(VENTASUSD)/1000000, BeneficiosAntes=sum(PROFITSBEFOREUSD)/1000000, BeneficiosDespues=sum(PROFITSAFTERUSD)/1000000, Impuestos1=sum(CITAXUSD)/1000000, Impuestos2=sum(TAXTOTALUSD)/1000000)

DECILES_VENTAS$TAX1_VENTAS<-DECILES_VENTAS$Impuestos1/DECILES_VENTAS$Ventas
DECILES_VENTAS$TAX2_VENTAS<-DECILES_VENTAS$Impuestos2/DECILES_VENTAS$Ventas

DECILES_VENTAS$TAX1_BENEFICIOS<-DECILES_VENTAS$Impuestos1/DECILES_VENTAS$BeneficiosDespues
DECILES_VENTAS$TAX2_BENEFICIOS<-DECILES_VENTAS$Impuestos2/DECILES_VENTAS$BeneficiosDespues


# GRAFICA 1 ---------------------------------------------------------------
#Proporción de trabajadores informales por deciles de ventas de las firmas
#“InformalProportionDemandPercentile”

TS<-CENSO%>%
  group_by(NUMTRAB)%>%
  summarise(firmas=sum(nn), ventas=mean(CAP5MONTO1, na.rm = T))
TS$decil<- decile(TS$ventas)
TS$NUMTRAB
TS<-left_join(TS,MOMENTO1,by=c("NUMTRAB"="Tamano de la empresa"))
TS$`Numero de trabajadores`<- TS$firmas*TS$NUMTRAB
TS$Informales<-TS$firmas*TS$NUMTRAB*TS$Informalidad/100
TS$Formales<-TS$`Numero de trabajadores`-TS$Informales

INF<-TS%>%
  group_by(decil)%>%
  summarise(Total=sum(`Numero de trabajadores`), Informales=sum(Informales,na.rm = T), Formales=sum(Formales, na.rm = T))

INF$Informalidad<-INF$Informales/INF$Total
INF$Formalidad<-1-INF$Informalidad
INF$Percentil<-c(seq(0.1,1,0.1))

G1<-ggplot(data=INF,aes(x=Percentil,y=Informalidad))+
  geom_line()+
  labs(x=expression(theta~e),y="Proportion Informal Demand")
setwd("C:/Users/mr.porras10/OneDrive - Universidad de Los Andes/W/New/Graficas")
dev.set()
png(file="InformalProportionDemandPercentileE.png",width=1600,height=850)
G1
dev.off()


# GRAFICA 2 ---------------------------------------------------------------
#Proporción de trabajadores formales por deciles de ventas de las firmas
#“FormalProportionDemandPercentile”

G2<-ggplot(data=INF,aes(x=Percentil,y=Formalidad))+
  geom_line()+
  labs(x=expression(theta~e),y="Proportion Formal Demand")
setwd("C:/Users/mr.porras10/OneDrive - Universidad de Los Andes/W/New/Graficas")
dev.set()
png(file="FormalProportionDemandPercentileE.png",width=1600,height=850)
G2
dev.off()


# GRAFICA 3 ---------------------------------------------------------------
#Número de trabajadores totales (formales e informales) por deciles de ventas de las firmas
#"TotalLaborDemandPercentile"

MOMENTO2$Percentil<-c(seq(0.1,1,0.1))

G3<-ggplot(data=MOMENTO2,aes(x=Percentil,y=`Numero de trabajadores`))+
  geom_line()+
  labs(x=expression(theta~e),y="Total Labor Demand")

setwd("C:/Users/mr.porras10/OneDrive - Universidad de Los Andes/W/New/Graficas")
dev.set()
png(file="TotalLaborDemandE.png",width=1600,height=850)
G3
dev.off()



# GRAFICA 4 ---------------------------------------------------------------
#Niveles de beneficios antes de impuestos por deciles de ventas.
#“PretaxProfitPercentile”.

DECILES_VENTAS$Percentil<-c(seq(0.1,1,0.1))

G4<-ggplot(data=DECILES_VENTAS,aes(x=Percentil,y=BeneficiosAntes))+
  geom_line()+
  labs(x=expression(theta~e),y="Profits Distribution (Millions of dollars)")

setwd("C:/Users/mr.porras10/OneDrive - Universidad de Los Andes/W/New/Graficas")
dev.set()
png(file="PretaxProfitPercentileE.png",width=1600,height=850)
G4
dev.off()


# GRAFICA 5 ---------------------------------------------------------------
#Niveles de beneficios después de impuestos por deciles de ventas.
#“AfterTaxProfitPercentile”

G5<-ggplot(data=DECILES_VENTAS,aes(x=Percentil,y=BeneficiosDespues))+
  geom_line()+
  labs(x=expression(theta~e),y="After Tax Profits Distribution (Millions of dollars)")

setwd("C:/Users/mr.porras10/OneDrive - Universidad de Los Andes/W/New/Graficas")
dev.set()
png(file="AfterTaxProfitPercentileE.png",width=1600,height=850)
G5
dev.off()



# GRAFICA 6 ---------------------------------------------------------------
#Niveles de venta (valor) en cada decil de ventas.
#“ProductionPercentile”

G6<-ggplot(data=DECILES_VENTAS,aes(x=Percentil,y=Ventas))+
  geom_line()+
  labs(x=expression(theta~e),y="Production Distribution")

setwd("C:/Users/mr.porras10/OneDrive - Universidad de Los Andes/W/New/Graficas")
dev.set()
png(file="ProductionPercentileE.png",width=1600,height=850)
G6
dev.off()


# GRAFICA 7 ---------------------------------------------------------------
#Valor de pago de impuesto en cada decil de ventas (CITAX)
G7.1<-ggplot(data=DECILES_VENTAS,aes(x=Percentil,y=Impuestos1))+
  geom_line()+
  labs(x=expression(theta~e),y="Taxes Payed (CITAX - Millions of dollars)")

setwd("C:/Users/mr.porras10/OneDrive - Universidad de Los Andes/W/New/Graficas")
dev.set()
png(file="Taxes1E.png",width=1600,height=850)
G7.1
dev.off()

#Valor de pago de impuesto en cada decil de ventas (CITAX+Transferencias)
G7.2<-ggplot(data=DECILES_VENTAS,aes(x=Percentil,y=Impuestos2))+
  geom_line()+
  labs(x=expression(theta~e),y="Taxes Payed (CITAX+Transferencias - Millions of dollars)")

setwd("C:/Users/mr.porras10/OneDrive - Universidad de Los Andes/W/New/Graficas")
dev.set()
png(file="Taxes2E.png",width=1600,height=850)
G7.2
dev.off()


# GRAFICA 8 ---------------------------------------------------------------
#Valor de pago de impuesto como proporción de ventas, en cada decil de ventas (CITAX)
G8.1<-ggplot(data=DECILES_VENTAS,aes(x=Percentil,y=TAX1_VENTAS))+
  geom_line()+
  labs(x=expression(theta~e),y="Taxes Payed (CITAX) / Production")

setwd("C:/Users/mr.porras10/OneDrive - Universidad de Los Andes/W/New/Graficas")
dev.set()
png(file="Taxes_Sales1E.png",width=1600,height=850)
G8.1
dev.off()

#Valor de pago de impuesto como proporción de ventas, en cada decil de ventas (CITAX+TRANSFERENCIAS)
G8.2<-ggplot(data=DECILES_VENTAS,aes(x=Percentil,y=TAX2_VENTAS))+
  geom_line()+
  labs(x=expression(theta~e),y="Taxes Payed (CITAX+Transferencias) / Production")

setwd("C:/Users/mr.porras10/OneDrive - Universidad de Los Andes/W/New/Graficas")
dev.set()
png(file="Taxes_Sales2E.png",width=1600,height=850)
G8.2
dev.off()

# GRAFICA 9 ---------------------------------------------------------------
#Valor de pago de impuesto como proporción de los beneficios, en cada decil de ventas
G9.1<-ggplot(data=DECILES_VENTAS,aes(x=Percentil,y=TAX1_BENEFICIOS))+
  geom_line()+
  labs(x=expression(theta~e),y="Taxes Payed (CITAX)  / Profits")

setwd("C:/Users/mr.porras10/OneDrive - Universidad de Los Andes/W/New/Graficas")
dev.set()
png(file="Taxes_Profits1E.png",width=1600,height=850)
G9.1
dev.off()

#Valor de pago de impuesto en cada decil de ventas (CITAX+Transferencias)
G9.2<-ggplot(data=DECILES_VENTAS,aes(x=Percentil,y=TAX2_BENEFICIOS))+
  geom_line()+
  labs(x=expression(theta~e),y="Taxes Payed  (CITAX+Transferencias)  / Profits")

setwd("C:/Users/mr.porras10/OneDrive - Universidad de Los Andes/W/New/Graficas")
dev.set()
png(file="Taxes_Profits2E.png",width=1600,height=850)
G9.2
dev.off()


# GRAFICA 10 --------------------------------------------------------------
#Distribución de ingresos de los trabajadores formales

mins<-min(ENAHOT$salarioUSD)
maxs<-max(ENAHOT$salarioUSD)
ENAHOT$salarionorm<-(ENAHOT$salarioUSD-mins)/(maxs-mins)


G10<-ggplot(ENAHOT[which(ENAHOT$informal_empleado==0),],aes(x=salarionorm))+
  geom_histogram(aes(y=..density..))+
  stat_density(geom = "line")+
  labs(x="w (Formal workers)", y="")

setwd("C:/Users/mr.porras10/OneDrive - Universidad de Los Andes/W/New/Graficas")
dev.set()
png(file="FormalSupplyPercentiles.png",width=1600,height=850)
G10
dev.off()


# GRAFICA 11 --------------------------------------------------------------
#Distribución de ingresos de los trabajadores informales

G11<-ggplot(ENAHOT[which(ENAHOT$informal_empleado==1),],aes(x=salarionorm))+
  geom_histogram(aes(y=..density..))+
  stat_density(geom = "line")+
  labs(x="w (Informal workers)", y="")

setwd("C:/Users/mr.porras10/OneDrive - Universidad de Los Andes/W/New/Graficas")
dev.set()
png(file="InformalSupplyPercentiles.png",width=1600,height=850)
G11
dev.off()

# GRAFICA 12 --------------------------------------------------------------
#Distribución de ingresos de los trabajadores (formales e informales)

G12<-ggplot(ENAHOT,aes(x=salarionorm))+
  geom_histogram(aes(y=..density..))+
  stat_density(geom = "line")+
  labs(x="w (Total labor supply)", y="")

setwd("C:/Users/mr.porras10/OneDrive - Universidad de Los Andes/W/New/Graficas")
dev.set()
png(file="TotalLaborSupplyPercentiles.png",width=1600,height=850)
G12
dev.off()


# GRAFICA 13 --------------------------------------------------------------
#Distribución de trabajadores informales como proporción de todos los trabajadores organizados por deciles de ingresos.

MOMENTO8$Informales<-MOMENTO8$Observaciones*MOMENTO8$`Informalidad (%)`/100
MOMENTO8$Propinformales<-MOMENTO8$Informales/sum(MOMENTO8$Observaciones)
MOMENTO8$Percentil<-c(seq(0.1,1,0.1))

G13<-ggplot(data=MOMENTO8,aes(x=Percentil,y=Propinformales))+
  geom_line()+
  labs(x=expression(theta~w),y="Proportion Informal Supply")

setwd("C:/Users/mr.porras10/OneDrive - Universidad de Los Andes/W/New/Graficas")
dev.set()
png(file="InformalProportionSupplyPercentiles.png",width=1600,height=850)
G13
dev.off()
