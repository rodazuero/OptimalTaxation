
# IMPORTAR DATOS ----------------------------------------------------------
library(dplyr)
library(StatMeasures)
library(xtable)
library(ggplot2)
setwd(choose.dir())
ENAHO<-read.csv("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/DataAnalysis/All/Census/Modified/OptimaltaxationSubSampleENAHO.csv", header = T, sep=",")
CENSO<-read.csv("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/DataAnalysis/All/Census/Modified/OptimaltaxationSubSampleCenso.csv", header = T, sep=",")

# MOMENTO 1 ---------------------------------------------------------------
# Informalidad y numero de trabajadores por tama?o de la empresa (ENAHO)

#Se agrupan las empresas con mas de 50 trabajadores en una misma categoria,
#si se quiere la informacion sin restringir el tamaño de la empresa se debe 
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
#Informalidad, numero de firmas y de trabajadores por tama?o de la empresa
#(CENSO+ENAHO)

#Se agrupan las empresas con mas de 50 trabajadores en una misma categoria,
#si se quiere la informacion sin restringir el tamaño de la empresa se debe 
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
#Distribuci?n de la produccion de las firmas (CENSO)

CENSO$PRODUCCIONUSD<-CENSO$CAP5MONTO9*0.315

M6<-distribuciones(CENSO$PRODUCCIONUSD)
print(xtable(M6),include.rownames=F)


# MOMENTO 7 ---------------------------------------------------------------
#Pago de impuestos y nivel de producci?n de la firma (CENSO)

CENSO$PROFITSUSD<-CENSO$CAP5MONTO34*0.315

#A)Impuesto a la renta (todos los regimenes)-CITAX

CENSO$TAX_PRODA<-CENSO$CITAXUSD/CENSO$PRODUCCIONUSD*100
CENSO$TAX_PROFA<-CENSO$CITAXUSD/CENSO$PROFITSUSD*100
CENSO$TAX_PROFA[CENSO$PROFITSUSD==0]<-NA

M7A.1<-distribuciones(CENSO$CITAXUSD)
M7A.2<-distribuciones(CENSO$TAX_PRODA)
M7A.3<-distribuciones(CENSO$TAX_PROFA)

M7A<-cbind(M7A.1,M7A.2[,2],M7A.3[,2])
colnames(M7A)<-c(".", "Impuestos","(%)Impuestos/Producci?n", "(%)Impuestos/Beneficios" )
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
colnames(M7B)<-c(".", "Impuestos","(%)Impuestos/Producci?n", "(%)Impuestos/Beneficios" )
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
