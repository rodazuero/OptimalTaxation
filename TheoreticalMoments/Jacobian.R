
#Author: Rodrigo Azuero

#Created: Nov, 2018. 

#Description:
#This files analyzes the relationship between parameters and the moments, including how likely are 
#combination of parameters to generate an equilibrium. 

#Once a combination of parameters is generated, together with the predicted moments and the
#corresponding distances, it can be analyzed. 

#------------#
#Housekeeping#
#------------#


#install.packages("ggplot2", dependencies=TRUE, INSTALL_opts = c('--no-lock'))
#install.packages('tidyverse',dependencies=TRUE)
#install.packages('rlang',dependencies=TRUE)
#devtools::install_github("r-lib/rlang", build_vignettes = TRUE)
library(nleqslv)
library(lattice)
library("gridExtra")
library("cowplot")
library(ggpubr)
library(ggplot2)
library(reshape)
library(grid)

library(data.table)
library(pastecs)


rm(list=ls(all=TRUE))



#Obtaining the empirical moments#
source('/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/LocalCopy/OptimalTaxation/EmpiricalMoments/Momentos.R')



#In this section we modify the directory to specify where are we going to extract the parameters and theoretical
#moments that will be compared with the empircal ones. 

CD<-'/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/LocalCopy/OptimalTaxation/AWS/InAws/Output/V1/'
setwd(CD)




#----------------------#
#Loading the files     #
#----------------------#


#Loading the parameters
Parameters<-read.csv("ParametersCSV.csv", header = F, sep=",")
names<-c("aalpha","ggamma","ddelta","bbeta","ssigma","kkappa","rrho","psi",
         "chi","mmu1","mmu2","ssigma1","ssigma2","rho12","c","NA")
colnames(Parameters)<-names
#Loading distance from theoretical to empirical moments
DistanceMom<-read.csv("DistanceMoments.csv", header = F, sep=",")
#stat.desc(subset(DistanceMom[,1],DistanceMom[,1]!=5.71312e+01)) 

colnames(DistanceMom)<-"Distance"

#Loading equilibrium //Sometimes this file might not be present.
EqValues<-read.csv("EquilibriumValue.csv", header = F, sep=",")
EqValues<-EqValues[,1]
#Remove the first equilibrium value
#EqValues<-EqValues[-1]


#Loading the moments
ProductionMoments<-read.csv("ThMoments0CSV.csv", header = F, sep=",")
colnames(ProductionMoments)<-paste0(rep("Production",10),colnames(ProductionMoments))


Taxesproportionally<-read.csv("ThMoments1CSV.csv", header = F, sep=",")
colnames(Taxesproportionally)<-paste0(rep("Taxesproportionally",10),colnames(Taxesproportionally))


TotalWorkersDemanded<-read.csv("ThMoments2CSV.csv", header = F, sep=",")
colnames(TotalWorkersDemanded)<-paste0(rep("TotalWorkersDemanded",10),colnames(TotalWorkersDemanded))

InformalDemandProportion<-read.csv("ThMoments3CSV.csv", header = F, sep=",")
colnames(InformalDemandProportion)<-paste0(rep("InformalDemandProportion",10),colnames(InformalDemandProportion))

IncomeDistribution<-read.csv("ThMoments4CSV.csv", header = F, sep=",")
colnames(IncomeDistribution)<-paste0(rep("IncomeDistribution",10),colnames(IncomeDistribution))

InformalLaborSupplyProp<-read.csv("ThMoments5CSV.csv", header = F, sep=",")
colnames(InformalLaborSupplyProp)<-paste0(rep("InformalLaborSupplyProp",10),colnames(InformalLaborSupplyProp))

TotalLaborSupply<-read.csv("ThMoments6CSV.csv", header = F, sep=",")
colnames(TotalLaborSupply)<-paste0(rep("TotalLaborSupply",10),colnames(TotalLaborSupply))



#We only need first column of prop entrepreneurs
PropEntrepreneurs<-read.csv("ThMoments7CSV.csv", header = F, sep=",")
PropEntrepreneurs<-PropEntrepreneurs[,1]


#Moment 8: proportion of income going to wages
PropAlpha<-read.csv("ThMoments8CSV.csv", header = F, sep=",")
PropAlpha<-PropAlpha[,1]

#EquilibriumValues has 1046. Take one less\
#EqValues<-EqValues[1:10465]

AllMmoments<-data.table(ProductionMoments,Taxesproportionally,
                        TotalWorkersDemanded,InformalDemandProportion,
                        IncomeDistribution,InformalLaborSupplyProp,
                        TotalLaborSupply,PropEntrepreneurs,PropAlpha)





#We need a column indicating if the equilibrium is successful or not. We will call it Indic. 
#This will be if: 
#1. TotalLaborSupply!=0.
AllMmoments[,Indic1 :=1-(TotalLaborSupplyV1+TotalLaborSupplyV2+TotalLaborSupplyV3
                        +TotalLaborSupplyV4+ TotalLaborSupplyV5 + TotalLaborSupplyV6+
                          TotalLaborSupplyV7+TotalLaborSupplyV8+TotalLaborSupplyV9==0)]

#2. If Total Labor supply ==1000
AllMmoments[,Indic2 :=1-(TotalLaborSupplyV1==1000)]

#3. If total Production ==0
AllMmoments[,Indic3 :=1-(ProductionV1+ProductionV2+ProductionV3
                        +ProductionV4+ ProductionV5 + ProductionV6+
                          ProductionV7+ProductionV8+ProductionV9==0)]

#4. If total taxes are zero:
AllMmoments[,Indic4:=1-(TaxesproportionallyV1+TaxesproportionallyV2+TaxesproportionallyV3
                        +TaxesproportionallyV4+TaxesproportionallyV5+TaxesproportionallyV6
                        +TaxesproportionallyV7+TaxesproportionallyV8+TaxesproportionallyV9==0)]

#5. If total labor demand is zero:
AllMmoments[,Indic5:=1-(TotalWorkersDemandedV1+TotalWorkersDemandedV2+TotalWorkersDemandedV3
                        +TotalWorkersDemandedV4+TotalWorkersDemandedV5+TotalWorkersDemandedV6
                        +TotalWorkersDemandedV7+TotalWorkersDemandedV8+TotalWorkersDemandedV9==0)]

#6. Replace informal proportion demand if -nan to zero. 
AllMmoments[is.na(InformalDemandProportionV1), InformalDemandProportionV1 := 0]
AllMmoments[is.na(InformalDemandProportionV2), InformalDemandProportionV2 := 0]
AllMmoments[is.na(InformalDemandProportionV3), InformalDemandProportionV3 := 0]
AllMmoments[is.na(InformalDemandProportionV4), InformalDemandProportionV4 := 0]
AllMmoments[is.na(InformalDemandProportionV5), InformalDemandProportionV5 := 0]
AllMmoments[is.na(InformalDemandProportionV6), InformalDemandProportionV6 := 0]
AllMmoments[is.na(InformalDemandProportionV7), InformalDemandProportionV7 := 0]
AllMmoments[is.na(InformalDemandProportionV8), InformalDemandProportionV8 := 0]
AllMmoments[is.na(InformalDemandProportionV9), InformalDemandProportionV9 := 0]


#4. Total
AllMmoments[,Indic :=1-1*(Indic1==0|Indic2==0|Indic3==0|Indic4==0|Indic5==0)]

#5. And remove Indics
AllMmoments[,Indic1:=NULL]
AllMmoments[,Indic2:=NULL]
AllMmoments[,Indic3:=NULL]
AllMmoments[,Indic4:=NULL]
AllMmoments[,Indic5:=NULL]


##Names modified to have C, the constant in the production function as additional parameter
names<-c("aalpha","ggamma","ddelta","bbeta","ssigma","kkappa","psi","chi",
         "rho","mmu1","mmu2","ssigma1","ssigma2","rho12","c","NA")
colnames(Parameters)<-names

Everything<-data.table(Parameters,DistanceMom,AllMmoments)

#----------------------------------------------------------------------#
#We only take those where positive number of entrepreneurs where found #
#----------------------------------------------------------------------#

#The combinations with problems throw 666 into  InformalDemandProportion[,1]. 
#We take only the moments where InformalDemandProprtion[,1]!=666


#-----------------------------------#
#Analyzing Parameters and equilibria#
#-----------------------------------#


Parameters<-data.table(Parameters,AllMmoments$Indic)
setnames(Parameters,17,"Indic")
#For those where there is equilibrium, set Indic =1. 

#Subsetting the two 
ParametersEquilibria=subset(Parameters,Indic!=0)
ParametersNoEquilibria=subset(Parameters,Indic==0)


#----------------------#
#Analyzing the distance#
#----------------------#


#Equilibrium
EverythingEqOriginal<-Everything
EverythingEq<-subset(Everything,Indic!=0)
EverythingEquilibrium<-EverythingEq





#Specify if you want to do a subset of the moments. For example:
#limiting equilibria where entrepreneurs >0.2, etc. 

SUBSET=TRUE

if(SUBSET==TRUE){
  #Exclude predictions that do not have positive labor supply informal
  #EverythingEq<-subset(EverythingEq,InformalLaborSupplyPropV1>0.2)
  
  #Exclude predictions that do not have positive informal labor demand
  #EverythingEq<-subset(EverythingEq,InformalDemandProportionV4>0.5)
  
  #Leave those with high entrepreneurs
  dim(subset(EverythingEq,PropEntrepreneurs>0.20))
  EverythingEq<-subset(EverythingEq,PropEntrepreneurs>0.20)
  
  
  #Alpha more than 0.45
  #EverythingEq<-subset(EverythingEq,aalpha>0.46)
  
  #Exclude those with very few informal demand
  EverythingEq<-subset(EverythingEq,InformalLaborSupplyPropV1>0.5)
  
  #Doing the restrictions necessary in the data
  EverythingEq<-subset(EverythingEq,PropEntrepreneurs>0.25)
}


#Sort by the distance
EverythingEqDistance<-EverythingEq[order(Distance)]







#------------------------------------------------#
# Comparing the theoretical and empirical moments#
#------------------------------------------------#


#Generating the data.table to create graphs
Comparing<-data.table(as.factor(c(rep(1,9),rep(0,9))))
colnames(Comparing)<-c("Sample")
Comparing$Decile<-c(seq(1,9,1),seq(1,9,1))



for(i in 1:46){
  #Deciding which observation to be analyzed
  #i=1

  
  Folder=paste('Model',i, sep = "")
  cd=paste(CD,"Modelfit/",sep="")
  setwd(cd)
  newcd=paste(cd,Folder,sep="")
  dir.create(Folder)
  setwd(newcd)
  
  #Including the theoretical and empirical moments
  
  #I will do a subsample of sobols for which entrepreneurs where relevant:
  #And store the original EverythingEqDistance in other stuff
  
  #Only run it once to keep it!:
  
  #EverythingEqDistanceOriginal<-EverythingEqDistance
  #If you want to subset the equilibria analyzed this is the place to do it. 
  #Otherwise, just keep going. 
  
  #And the rest:
  #EverythingEqDistanceEntrep<-subset(EverythingEqDistanceOriginal,PropEntrepreneurs>0.19)
  #max(EverythingEqDistanceEntrep$InformalDemandProportionV4)
  #EverythingEqDistanceEntrep$InformalDemandProportionV5[i]
  #EverythingEqDistance<-EverythingEqDistanceEntrep
  
  
  
  #0. Printing the parameters
  write.csv(Parameters[i,],file='Param.csv',row.names = FALSE)
  
  
  write.csv(EverythingEqDistance[i,],file='Everything.csv')
  
  
  #1. Production 
  
  
  
  ProductionTheoretical<-c(EverythingEqDistance$ProductionV1[i],
                EverythingEqDistance$ProductionV2[i],
                EverythingEqDistance$ProductionV3[i],
                EverythingEqDistance$ProductionV4[i],
                EverythingEqDistance$ProductionV5[i],
                EverythingEqDistance$ProductionV6[i],
                EverythingEqDistance$ProductionV7[i],
                EverythingEqDistance$ProductionV8[i],
                EverythingEqDistance$ProductionV9[i])
  
  
  Comparing$Production <-  c(ProductionTheoretical/EverythingEqDistance$ProductionV5[i], MOMENTO7A$Produccion[1:9]/MOMENTO7A$Produccion[5])
  
  
  Comparing$ProductionReal <-  c(ProductionTheoretical, MOMENTO7A$Produccion[1:9])
  
  #2. Taxes Proportionally cumulative
  
  
  
  Taxesproportionally<-c(EverythingEqDistance$TaxesproportionallyV1[i],
                           EverythingEqDistance$TaxesproportionallyV2[i],
                           EverythingEqDistance$TaxesproportionallyV3[i],
                           EverythingEqDistance$TaxesproportionallyV4[i],
                           EverythingEqDistance$TaxesproportionallyV5[i],
                           EverythingEqDistance$TaxesproportionallyV6[i],
                           EverythingEqDistance$TaxesproportionallyV7[i],
                           EverythingEqDistance$TaxesproportionallyV8[i],
                           EverythingEqDistance$TaxesproportionallyV9[i])
  
  
  Comparing$TaxesPayedProportionally <-  c(Taxesproportionally, MOMENTO7A$`Proporción de pago de impuestos acumulada`[1:9]/sum(MOMENTO7A$`Proporción de pago de impuestos acumulada`[1:9]))
  
  
  #3. Total workers demanded
  
  DemandTotalWorkers<-c(EverythingEqDistance$TotalWorkersDemandedV1[i],
                         EverythingEqDistance$TotalWorkersDemandedV2[i],
                         EverythingEqDistance$TotalWorkersDemandedV3[i],
                         EverythingEqDistance$TotalWorkersDemandedV4[i],
                         EverythingEqDistance$TotalWorkersDemandedV5[i],
                         EverythingEqDistance$TotalWorkersDemandedV6[i],
                         EverythingEqDistance$TotalWorkersDemandedV7[i],
                         EverythingEqDistance$TotalWorkersDemandedV8[i],
                         EverythingEqDistance$TotalWorkersDemandedV9[i])
  
  
  Comparing$TotalDemandWorkers <-  c(DemandTotalWorkers/DemandTotalWorkers[5],MOMENTO2$`Numero de trabajadores`[1:9]/MOMENTO2$`Numero de trabajadores`[5] )
  
  Comparing$TotalDemandWorkersReal <-  c(DemandTotalWorkers,MOMENTO2$`Numero de trabajadores`[1:9] )
  
  
  #4. Informal Demand proportion
  InformalDemandProportion<-c(EverythingEqDistance$InformalDemandProportionV4[i],
                              EverythingEqDistance$InformalDemandProportionV7[i],
                              EverythingEqDistance$InformalDemandProportionV8[i])
  
  InformalPROPDemand <- data.table(c(MOMENTO3$Informalidad[1]/100,MOMENTO3$Informalidad[2]/100,MOMENTO3$Informalidad[3]/100,
                                     InformalDemandProportion))
  
  InformalPROPDemand$Decile<-c(4,7,8,4,7,8)
  InformalPROPDemand$Sample<-as.factor(c(0,0,0,1,1,1))
  
  
  
  #5. Income distribution
  
  
  Income<-c(EverythingEqDistance$IncomeDistributionV1[i],
            EverythingEqDistance$IncomeDistributionV2[i],
            EverythingEqDistance$IncomeDistributionV3[i],
            EverythingEqDistance$IncomeDistributionV4[i],
            EverythingEqDistance$IncomeDistributionV5[i],
            EverythingEqDistance$IncomeDistributionV6[i],
            EverythingEqDistance$IncomeDistributionV7[i],
            EverythingEqDistance$IncomeDistributionV8[i],
            EverythingEqDistance$IncomeDistributionV9[i])
  
  
  
  Comparing$IncomeDistribution <-  c(Income/EverythingEqDistance$IncomeDistributionV5[i],as.numeric(M5.3A[4:12,2])/as.numeric(M5.3A[8,2]) )
  
  
  Comparing$IncomeDistributionReal <-  c(Income,as.numeric(M5.3A[4:12,2]) )
  
  
  #6. Informal labor supply
  
  
  InformalLaborSupply<-c(EverythingEqDistance$InformalLaborSupplyPropV1[i],
            EverythingEqDistance$InformalLaborSupplyPropV2[i],
            EverythingEqDistance$InformalLaborSupplyPropV3[i],
            EverythingEqDistance$InformalLaborSupplyPropV4[i],
            EverythingEqDistance$InformalLaborSupplyPropV5[i],
            EverythingEqDistance$InformalLaborSupplyPropV6[i],
            EverythingEqDistance$InformalLaborSupplyPropV7[i],
            EverythingEqDistance$InformalLaborSupplyPropV8[i],
            EverythingEqDistance$InformalLaborSupplyPropV9[i])
  
  
  
  Comparing$InformalLaborSupply <-  c(InformalLaborSupply,MOMENTO8$`Informalidad (%)`[1:9]/100 )
  
  
  
  
  
  #7. Total labor supply
  
  
  TotalLaborSupply<-c(EverythingEqDistance$TotalLaborSupplyV1[i],
                         EverythingEqDistance$TotalLaborSupplyV2[i],
                         EverythingEqDistance$TotalLaborSupplyV3[i],
                         EverythingEqDistance$TotalLaborSupplyV4[i],
                         EverythingEqDistance$TotalLaborSupplyV5[i],
                         EverythingEqDistance$TotalLaborSupplyV6[i],
                         EverythingEqDistance$TotalLaborSupplyV7[i],
                         EverythingEqDistance$TotalLaborSupplyV8[i],
                         EverythingEqDistance$TotalLaborSupplyV9[i])
  
  
  
  Comparing$TotalLaborSupply <-  c(TotalLaborSupply,MOMENTO19$`Version 3A`[1:9] )
  
  Comparing$TotalLaborSupply2<-c(TotalLaborSupply,MOMENTO19$`Version 3A`[1:9])
  
  #8. Proportion of entrepreneurs
  PropEntrep<-data.table(c(EverythingEqDistance$PropEntrepreneurs[i],(M4[2,2]/(M4[2,1]+M4[2,2]))))
  PropEntrep$Sample<-c(as.factor(c(0,1)))
  print(PropEntrep)
  InformalPROPDemand
  
  
  #Subsetting the ones where entrepreneurs are more than 25%:
  
  
  
  #-----#
  #Plots#
  #-----#
  
  
  
  
  #Set directory
  #setwd('/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolsGenerated/Dim15Secondary/ModelFit')
  print("hahahahaha")
  #Size of line
  sizeline=10
  sizerel=7
  
  
  #Production relative to median
  Production<-ggplot(data=Comparing,aes(x=Decile,y=Production,colour=Sample))+geom_line(size=sizeline)+geom_point()
  Production<-Production+scale_colour_discrete(labels=c("Data","Model")  )
  Production<-Production+ theme_bw()
  Production<-Production+scale_x_continuous(breaks = seq(1,9,1))
  Production<-Production + theme(
    plot.title = element_text(hjust=0.5,size = rel(sizerel)),
    axis.title=element_text(size = rel(sizerel)),
    axis.text.x=element_text(size = rel(sizerel)),
    axis.text.y=element_text(size = rel(sizerel)),
    legend.text=element_text(size = rel(sizerel)),
    legend.title=element_text(size=rel(0)))
  
  Production<-Production+ggtitle(" ") +ylab("Relative to median")
  Production
  
  dev.set()
  png(file="Production.png",width=1600,height=850)
  print(Production)
  dev.off()


  #Production real
  Production<-ggplot(data=Comparing,aes(x=Decile,y=ProductionReal,colour=Sample))+geom_line(size=sizeline)+geom_point()
  Production<-Production+scale_colour_discrete(labels=c("Data","Model")  )
  Production<-Production+ theme_bw()
  Production<-Production+scale_x_continuous(breaks = seq(1,9,1))
  Production<-Production + theme(
    plot.title = element_text(hjust=0.5,size = rel(sizerel)),
    axis.title=element_text(size = rel(sizerel)),
    axis.text.x=element_text(size = rel(sizerel)),
    axis.text.y=element_text(size = rel(sizerel)),
    legend.text=element_text(size = rel(sizerel)),
    legend.title=element_text(size=rel(0)))
  
  Production<-Production+ggtitle(" ") +ylab("Thousands of S/.")
  Production
  
  dev.set()
  png(file="ProductionReal.png",width=1600,height=850)
  print(Production)
  dev.off()
  
  
  
  
  
  #Taxes
  Taxes<-ggplot(data=Comparing,aes(x=Decile,y=TaxesPayedProportionally,colour=Sample))+geom_line(size=sizeline)+geom_point()
  Taxes<-Taxes+scale_colour_discrete(labels=c("Data","Model")  )
  Taxes<-Taxes+ theme_bw()
  Taxes<-Taxes+scale_x_continuous(breaks = seq(1,9,1))
  Taxes<-Taxes + theme(
    plot.title = element_text(hjust=0.5,size = rel(sizerel)),
    axis.title=element_text(size = rel(sizerel)),
    axis.text.x=element_text(size = rel(sizerel)),
    axis.text.y=element_text(size = rel(sizerel)),
    legend.text=element_text(size = rel(sizerel)),
    legend.title=element_text(size=rel(0)))
  
  Taxes<-Taxes+ggtitle(" ") +ylab("%")
  Taxes
  
  dev.set()
  png(file="Taxes.png",width=1600,height=850)
  print(Taxes)
  dev.off()
  
  
  #Total demand workers
  
  #DemandWorkers median
  DemandWorkers<-ggplot(data=Comparing,aes(x=Decile,y=TotalDemandWorkers,colour=Sample))+geom_line(size=sizeline)+geom_point()
  DemandWorkers<-DemandWorkers+scale_colour_discrete(labels=c("Data","Model")  )
  DemandWorkers<-DemandWorkers+ theme_bw()
  DemandWorkers<-DemandWorkers+scale_x_continuous(breaks = seq(1,9,1))
  DemandWorkers<-DemandWorkers + theme(
    plot.title = element_text(hjust=0.5,size = rel(sizerel)),
    axis.title=element_text(size = rel(sizerel)),
    axis.text.x=element_text(size = rel(sizerel)),
    axis.text.y=element_text(size = rel(sizerel)),
    legend.text=element_text(size = rel(sizerel)),
    legend.title=element_text(size=rel(0)))
  
  DemandWorkers<-DemandWorkers+ggtitle(" ") +ylab("relative to median")
  DemandWorkers
  
  dev.set()
  png(file="DemandWorkers.png",width=1600,height=850)
  print(DemandWorkers)
  dev.off()
  
  
  DemandWorkers<-ggplot(data=Comparing,aes(x=Decile,y=TotalDemandWorkersReal,colour=Sample))+geom_line(size=sizeline)+geom_point()
  DemandWorkers<-DemandWorkers+scale_colour_discrete(labels=c("Data","Model")  )
  DemandWorkers<-DemandWorkers+ theme_bw()
  DemandWorkers<-DemandWorkers+scale_x_continuous(breaks = seq(1,9,1))
  DemandWorkers<-DemandWorkers + theme(
    plot.title = element_text(hjust=0.5,size = rel(sizerel)),
    axis.title=element_text(size = rel(sizerel)),
    axis.text.x=element_text(size = rel(sizerel)),
    axis.text.y=element_text(size = rel(sizerel)),
    legend.text=element_text(size = rel(sizerel)),
    legend.title=element_text(size=rel(0)))
  
  DemandWorkers<-DemandWorkers+ggtitle(" ") +ylab("relative to median")
  DemandWorkers
  
  dev.set()
  png(file="DemandWorkersReal.png",width=1600,height=850)
  print(DemandWorkers)
  dev.off()
  
  
  #IncomeDistribution relative to median
  IncomeDistribution<-ggplot(data=Comparing,aes(x=Decile,y=IncomeDistribution,colour=Sample))+geom_line(size=sizeline)+geom_point()
  IncomeDistribution<-IncomeDistribution+scale_colour_discrete(labels=c("Data","Model")  )
  IncomeDistribution<-IncomeDistribution+ theme_bw()
  IncomeDistribution<-IncomeDistribution+scale_x_continuous(breaks = seq(1,9,1))
  IncomeDistribution<-IncomeDistribution + theme(
    plot.title = element_text(hjust=0.5,size = rel(sizerel)),
    axis.title=element_text(size = rel(sizerel)),
    axis.text.x=element_text(size = rel(sizerel)),
    axis.text.y=element_text(size = rel(sizerel)),
    legend.text=element_text(size = rel(sizerel)),
    legend.title=element_text(size=rel(0)))
  
  IncomeDistribution<-IncomeDistribution+ggtitle(" ") +ylab("Relative to median")
  IncomeDistribution
  
  dev.set()
  png(file="IncomeDistribution.png",width=1600,height=850)
  print(IncomeDistribution)
  dev.off()
  
  #IncomeDistribution real
  IncomeDistribution<-ggplot(data=Comparing,aes(x=Decile,y=IncomeDistributionReal,colour=Sample))+geom_line(size=sizeline)+geom_point()
  IncomeDistribution<-IncomeDistribution+scale_colour_discrete(labels=c("Data","Model")  )
  IncomeDistribution<-IncomeDistribution+ theme_bw()
  IncomeDistribution<-IncomeDistribution+scale_x_continuous(breaks = seq(1,9,1))
  IncomeDistribution<-IncomeDistribution + theme(
    plot.title = element_text(hjust=0.5,size = rel(sizerel)),
    axis.title=element_text(size = rel(sizerel)),
    axis.text.x=element_text(size = rel(sizerel)),
    axis.text.y=element_text(size = rel(sizerel)),
    legend.text=element_text(size = rel(sizerel)),
    legend.title=element_text(size=rel(0)))
  
  IncomeDistribution<-IncomeDistribution+ggtitle(" ") +ylab("S/.")
  IncomeDistribution
  
  dev.set()
  png(file="IncomeDistributionReal.png",width=1600,height=850)
  print(IncomeDistribution)
  dev.off()
  
  #InformalLaborSupply
  InformalLaborSupply<-ggplot(data=Comparing,aes(x=Decile,y=InformalLaborSupply,colour=Sample))+geom_line(size=sizeline)+geom_point()
  InformalLaborSupply<-InformalLaborSupply+scale_colour_discrete(labels=c("Data","Model")  )
  InformalLaborSupply<-InformalLaborSupply+ theme_bw()
  InformalLaborSupply<-InformalLaborSupply+scale_x_continuous(breaks = seq(1,9,1))
  InformalLaborSupply<-InformalLaborSupply + theme(
    plot.title = element_text(hjust=0.5,size = rel(sizerel)),
    axis.title=element_text(size = rel(sizerel)),
    axis.text.x=element_text(size = rel(sizerel)),
    axis.text.y=element_text(size = rel(sizerel)),
    legend.text=element_text(size = rel(sizerel)),
    legend.title=element_text(size=rel(0)))
  
  InformalLaborSupply<-InformalLaborSupply+ggtitle(" ") +ylab("Proportion")
  InformalLaborSupply
  
  dev.set()
  png(file="InformalLaborSupply.png",width=1600,height=850)
  print(InformalLaborSupply)
  dev.off()
  
  
  #TotalLaborSupply 1 
  TotalLaborSupply<-ggplot(data=Comparing,aes(x=Decile,y=TotalLaborSupply,colour=Sample))+geom_line(size=sizeline)+geom_point()
  TotalLaborSupply<-TotalLaborSupply+scale_colour_discrete(labels=c("Data","Model")  )
  TotalLaborSupply<-TotalLaborSupply+ theme_bw()
  TotalLaborSupply<-TotalLaborSupply+scale_x_continuous(breaks = seq(1,9,1))
  TotalLaborSupply<-TotalLaborSupply + theme(
    plot.title = element_text(hjust=0.5,size = rel(sizerel)),
    axis.title=element_text(size = rel(sizerel)),
    axis.text.x=element_text(size = rel(sizerel)),
    axis.text.y=element_text(size = rel(sizerel)),
    legend.text=element_text(size = rel(sizerel)),
    legend.title=element_text(size=rel(0)))
  
  TotalLaborSupply<-TotalLaborSupply+ggtitle(" ") +ylab("#")
  
  
  dev.set()
  png(file="TotalLaborSupply.png",width=1600,height=850)
  print(TotalLaborSupply)
  dev.off()
  
  
  #TotalLaborSupply 2 
  TotalLaborSupply<-ggplot(data=Comparing,aes(x=Decile,y=TotalLaborSupply2,colour=Sample))+geom_line(size=sizeline)+geom_point()
  TotalLaborSupply<-TotalLaborSupply+scale_colour_discrete(labels=c("Data","Model")  )
  TotalLaborSupply<-TotalLaborSupply+ theme_bw()
  TotalLaborSupply<-TotalLaborSupply+scale_x_continuous(breaks = seq(1,9,1))
  TotalLaborSupply<-TotalLaborSupply + theme(
    plot.title = element_text(hjust=0.5,size = rel(sizerel)),
    axis.title=element_text(size = rel(sizerel)),
    axis.text.x=element_text(size = rel(sizerel)),
    axis.text.y=element_text(size = rel(sizerel)),
    legend.text=element_text(size = rel(sizerel)),
    legend.title=element_text(size=rel(0)))
  
  TotalLaborSupply<-TotalLaborSupply+ggtitle(" ") +ylab("#")
  
  
  dev.set()
  png(file="TotalLaborSupply2.png",width=1600,height=850)
  print(TotalLaborSupply)
  dev.off()
  
  
  
  #Informal proportion demand
  InformalDemand<-ggplot(data=InformalPROPDemand,aes(x=Decile,y=V1,colour=Sample))+geom_line(size=sizeline)+geom_point()
  InformalDemand<-InformalDemand+scale_colour_discrete(labels=c("Data","Model")  )
  InformalDemand<-InformalDemand+ theme_bw()
  InformalDemand<-InformalDemand+scale_x_continuous(breaks = seq(1,9,1))
  InformalDemand<-InformalDemand + theme(
    plot.title = element_text(hjust=0.5,size = rel(sizerel)),
    axis.title=element_text(size = rel(sizerel)),
    axis.text.x=element_text(size = rel(sizerel)),
    axis.text.y=element_text(size = rel(sizerel)),
    legend.text=element_text(size = rel(sizerel)),
    legend.title=element_text(size=rel(0)))
  
  InformalDemand<-InformalDemand+ggtitle(" ") +ylab("#")
  InformalDemand
  
  dev.set()
  png(file="InformalLaborDemand.png",width=1600,height=850)
  print(InformalDemand)
  dev.off()
  
  
  
  #Informal proportion demand ALL IN THEORETICAL PREDICTION
  InformalAll<-c(EverythingEqDistance$InformalDemandProportionV1[i],
                 EverythingEqDistance$InformalDemandProportionV2[i],
                 EverythingEqDistance$InformalDemandProportionV3[i],
                 EverythingEqDistance$InformalDemandProportionV4[i],
                 EverythingEqDistance$InformalDemandProportionV5[i],
                 EverythingEqDistance$InformalDemandProportionV6[i],
                 EverythingEqDistance$InformalDemandProportionV7[i],
                 EverythingEqDistance$InformalDemandProportionV8[i],
                 EverythingEqDistance$InformalDemandProportionV9[i])
  InformalAll<-cbind(InformalAll,seq(1,9,1))
  colnames(InformalAll)<-c("Informalidad","decile")
  InformalAll<-as.data.frame(InformalAll)
  PlotInformal<-ggplot(data=InformalAll,aes(y=Informalidad,x=decile))+geom_line()
  
  
  
  #Proportion of entrepreneurs
  PropEntrep$Sample<-c("Model","Data")
  PropEntrepGraph<-ggplot(data=PropEntrep,aes(y=V1,x=Sample))+geom_bar(stat="identity")
  PropEntrepGraph<-PropEntrepGraph+ggtitle("Proportion entrepreneurs") +ylab("%")
  dev.set()
  png(file="PropEntrepGraph.png",width=1600,height=850)
  print(PropEntrepGraph)
  dev.off()
  
  
  
  #Informal proportion demand total

  InformalDemand<-ggplot(data=InformalPROPDemand,aes(x=Decile,y=V1,colour=Sample))+geom_line(size=sizeline)+geom_point()
  InformalDemand<-InformalDemand+scale_colour_discrete(labels=c("Data","Model")  )
  InformalDemand<-InformalDemand+ theme_bw()
  InformalDemand<-InformalDemand+scale_x_continuous(breaks = seq(1,9,1))
  InformalDemand<-InformalDemand + theme(
    plot.title = element_text(hjust=0.5,size = rel(sizerel)),
    axis.title=element_text(size = rel(sizerel)),
    axis.text.x=element_text(size = rel(sizerel)),
    axis.text.y=element_text(size = rel(sizerel)),
    legend.text=element_text(size = rel(sizerel)),
    legend.title=element_text(size=rel(0)))
  
  InformalDemand<-InformalDemand+ggtitle(" ") +ylab("#")
  InformalDemand
  
  dev.set()
  png(file="InformalLaborDemand.png",width=1600,height=850)
  print(InformalDemand)
  dev.off()
  
  #Proportion of entrepreneurs
  PropEntrep$Sample<-c("Model","Data")
  PropEntrepGraph<-ggplot(data=PropEntrep,aes(y=V1,x=Sample))+geom_bar(stat="identity")
  PropEntrepGraph<-PropEntrepGraph+ggtitle("Proportion entrepreneurs") +ylab("%")
  dev.set()
  png(file="PropEntrepGraph.png",width=1600,height=850)
  print(PropEntrepGraph)
  dev.off()

}



#----------------------------------------------------------#
#Conditional distributions of each parameter on equilibrium#
#----------------------------------------------------------#
RUNTEST=FALSE
if (RUNTEST==TRUE){



#1. Aalpha
#1.1 Equilibrium


#Extra text
Desc<-stat.desc(ParametersEquilibria$aalpha)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))


#Plot
m <- ggplot(ParametersEquilibria, aes(x = aalpha))
m<-m + geom_histogram(binwidth = 0.02) 
m<-m+labs(title="Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=1,ymin=-6,ymax=-3)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gt <- ggplot_gtable(ggplot_build(m))
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gt)


#1.2. No equilibrium


#Extra text
Desc<-stat.desc(ParametersNoEquilibria$aalpha)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersNoEquilibria, aes(x = aalpha))
m<-m + geom_histogram(binwidth = 0.02) 
m<-m+labs(title="No Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=1,ymin=-100,ymax=-90)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gtN <- ggplot_gtable(ggplot_build(m))
gtN$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gtN)



dev.set()
png(file="Figures/aalpha.png",width=398,height=389)
ggarrange(gtN, gt)
dev.off()


#2 ggamma
#1.1 Equilibrium


#Extra text
Desc<-stat.desc(ParametersEquilibria$ggamma)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersEquilibria, aes(x = ggamma))
m<-m + geom_histogram(binwidth = 0.02) 
m<-m+labs(title="Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=10,ymin=-1.5,ymax=-0.5)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gt <- ggplot_gtable(ggplot_build(m))
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gt)


#1.2. No equilibrium


#Extra text
Desc<-stat.desc(ParametersNoEquilibria$ggamma)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersNoEquilibria, aes(x = ggamma))
m<-m + geom_histogram(binwidth = 0.02) 
m<-m+labs(title="No Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=10,ymin=-15,ymax=-10)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gtN <- ggplot_gtable(ggplot_build(m))
gtN$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gtN)



dev.set()
png(file="Figures/ggamma.png",width=398,height=389)
ggarrange(gtN, gt)
dev.off()


#3 ddelta
#1.1 Equilibrium


#Extra text
Desc<-stat.desc(ParametersEquilibria$ddelta)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersEquilibria, aes(x = ddelta))
m<-m + geom_histogram(binwidth = 10) 
m<-m+labs(title="Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=10,ymin=-1.5,ymax=-0.5)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gt <- ggplot_gtable(ggplot_build(m))
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gt)


#1.2. No equilibrium


#Extra text
Desc<-stat.desc(ParametersNoEquilibria$ddelta)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersNoEquilibria, aes(x = ddelta))
m<-m + geom_histogram(binwidth = 10) 
m<-m+labs(title="No Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=10,ymin=-15,ymax=-10)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gtN <- ggplot_gtable(ggplot_build(m))
gtN$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gtN)



dev.set()
png(file="Figures/ddelta.png",width=398,height=389)
ggarrange(gtN, gt)
dev.off()


#4 bbeta
#1.1 Equilibrium


#Extra text
Desc<-stat.desc(ParametersEquilibria$bbeta)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersEquilibria, aes(x = bbeta))
m<-m + geom_histogram(binwidth = 10) 
m<-m+labs(title="Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=10,ymin=-1.5,ymax=-0.5)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gt <- ggplot_gtable(ggplot_build(m))
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gt)


#1.2. No equilibrium


#Extra text
Desc<-stat.desc(ParametersNoEquilibria$bbeta)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersNoEquilibria, aes(x = bbeta))
m<-m + geom_histogram(binwidth = 10) 
m<-m+labs(title="No Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=10,ymin=-15,ymax=-10)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gtN <- ggplot_gtable(ggplot_build(m))
gtN$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gtN)



dev.set()
png(file="Figures/bbeta.png",width=398,height=389)
ggarrange(gtN, gt)
dev.off()


#5. ssigma
#1.1 Equilibrium


#Extra text
Desc<-stat.desc(ParametersEquilibria$ssigma)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersEquilibria, aes(x = ssigma))
m<-m + geom_histogram(binwidth = 0.1) 
m<-m+labs(title="Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=10,ymin=-1.5,ymax=-0.5)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gt <- ggplot_gtable(ggplot_build(m))
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gt)


#1.2. No equilibrium


#Extra text
Desc<-stat.desc(ParametersNoEquilibria$ssigma)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersNoEquilibria, aes(x = ssigma))
m<-m + geom_histogram(binwidth = 0.1) 
m<-m+labs(title="No Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=10,ymin=-15,ymax=-10)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gtN <- ggplot_gtable(ggplot_build(m))
gtN$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gtN)



dev.set()
png(file="Figures/ssigma.png",width=398,height=389)
ggarrange(gtN, gt)
dev.off()



#6. kkappa
#1.1 Equilibrium


#Extra text
Desc<-stat.desc(ParametersEquilibria$kkappa)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersEquilibria, aes(x = kkappa))
m<-m + geom_histogram(binwidth = 10) 
m<-m+labs(title="Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=10,ymin=-1.5,ymax=-0.5)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gt <- ggplot_gtable(ggplot_build(m))
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gt)


#1.2. No equilibrium


#Extra text
Desc<-stat.desc(ParametersNoEquilibria$kkappa)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersNoEquilibria, aes(x = kkappa))
m<-m + geom_histogram(binwidth = 10) 
m<-m+labs(title="No Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=10,ymin=-15,ymax=-10)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gtN <- ggplot_gtable(ggplot_build(m))
gtN$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gtN)



dev.set()
png(file="Figures/kkappa.png",width=398,height=389)
ggarrange(gtN, gt)
dev.off()




#7. psi
#1.1 Equilibrium


#Extra text
Desc<-stat.desc(ParametersEquilibria$psi)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersEquilibria, aes(x = psi))
m<-m + geom_histogram(binwidth = 0.02) 
m<-m+labs(title="Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=10,ymin=-1.5,ymax=-0.5)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gt <- ggplot_gtable(ggplot_build(m))
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gt)


#1.2. No equilibrium


#Extra text
Desc<-stat.desc(ParametersNoEquilibria$psi)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersNoEquilibria, aes(x = psi))
m<-m + geom_histogram(binwidth = 0.02) 
m<-m+labs(title="No Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=10,ymin=-15,ymax=-10)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gtN <- ggplot_gtable(ggplot_build(m))
gtN$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gtN)



dev.set()
png(file="Figures/psi.png",width=398,height=389)
ggarrange(gtN, gt)
dev.off()




#8. chi
#1.1 Equilibrium


#Extra text
Desc<-stat.desc(ParametersEquilibria$chi)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersEquilibria, aes(x = chi))
m<-m + geom_histogram(binwidth = 10) 
m<-m+labs(title="Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=10,ymin=-1.5,ymax=-0.5)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gt <- ggplot_gtable(ggplot_build(m))
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gt)


#1.2. No equilibrium


#Extra text
Desc<-stat.desc(ParametersNoEquilibria$chi)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersNoEquilibria, aes(x = chi))
m<-m + geom_histogram(binwidth = 10) 
m<-m+labs(title="No Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=10,ymin=-15,ymax=-10)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gtN <- ggplot_gtable(ggplot_build(m))
gtN$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gtN)



dev.set()
png(file="Figures/chi.png",width=398,height=389)
ggarrange(gtN, gt)
dev.off()




#9. rho
#1.1 Equilibrium


#Extra text
Desc<-stat.desc(ParametersEquilibria$rho)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersEquilibria, aes(x = rho))
m<-m + geom_histogram(binwidth = 0.02) 
m<-m+labs(title="Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=10,ymin=-1.5,ymax=-0.5)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gt <- ggplot_gtable(ggplot_build(m))
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gt)


#1.2. No equilibrium


#Extra text
Desc<-stat.desc(ParametersNoEquilibria$rho)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersNoEquilibria, aes(x = rho))
m<-m + geom_histogram(binwidth = 0.02) 
m<-m+labs(title="No Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=10,ymin=-15,ymax=-10)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gtN <- ggplot_gtable(ggplot_build(m))
gtN$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gtN)



dev.set()
png(file="Figures/rho.png",width=398,height=389)
ggarrange(gtN, gt)
dev.off()


#10. mmu1
#1.1 Equilibrium


#Extra text
Desc<-stat.desc(ParametersEquilibria$mmu1)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersEquilibria, aes(x = mmu1))
m<-m + geom_histogram(binwidth = 0.02) 
m<-m+labs(title="Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=10,ymin=-1.5,ymax=-0.5)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gt <- ggplot_gtable(ggplot_build(m))
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gt)


#1.2. No equilibrium


#Extra text
Desc<-stat.desc(ParametersNoEquilibria$mmu1)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersNoEquilibria, aes(x = mmu1))
m<-m + geom_histogram(binwidth = 0.02) 
m<-m+labs(title="No Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=10,ymin=-15,ymax=-10)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gtN <- ggplot_gtable(ggplot_build(m))
gtN$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gtN)



dev.set()
png(file="Figures/mmu1.png",width=398,height=389)
ggarrange(gtN, gt)
dev.off()


#11. mmu2
#1.1 Equilibrium


#Extra text
Desc<-stat.desc(ParametersEquilibria$mmu2)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersEquilibria, aes(x = mmu2))
m<-m + geom_histogram(binwidth = 0.02) 
m<-m+labs(title="Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=10,ymin=-1.5,ymax=-0.5)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gt <- ggplot_gtable(ggplot_build(m))
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gt)


#1.2. No equilibrium


#Extra text
Desc<-stat.desc(ParametersNoEquilibria$mmu2)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersNoEquilibria, aes(x = mmu2))
m<-m + geom_histogram(binwidth = 0.02) 
m<-m+labs(title="No Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=10,ymin=-15,ymax=-10)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gtN <- ggplot_gtable(ggplot_build(m))
gtN$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gtN)



dev.set()
png(file="Figures/mmu2.png",width=398,height=389)
ggarrange(gtN, gt)
dev.off()




#12. ssigma1
#1.1 Equilibrium


#Extra text
Desc<-stat.desc(ParametersEquilibria$ssigma1)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersEquilibria, aes(x = ssigma1))
m<-m + geom_histogram(binwidth = 0.02) 
m<-m+labs(title="Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=10,ymin=-1.5,ymax=-0.5)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gt <- ggplot_gtable(ggplot_build(m))
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gt)


#1.2. No equilibrium


#Extra text
Desc<-stat.desc(ParametersNoEquilibria$ssigma1)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersNoEquilibria, aes(x = ssigma1))
m<-m + geom_histogram(binwidth = 0.02) 
m<-m+labs(title="No Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=10,ymin=-15,ymax=-10)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gtN <- ggplot_gtable(ggplot_build(m))
gtN$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gtN)



dev.set()
png(file="Figures/ssigma1.png",width=398,height=389)
ggarrange(gtN, gt)
dev.off()




#13. ssigma2
#1.1 Equilibrium


#Extra text
Desc<-stat.desc(ParametersEquilibria$ssigma2)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersEquilibria, aes(x = ssigma2))
m<-m + geom_histogram(binwidth = 0.02) 
m<-m+labs(title="Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=10,ymin=-1.5,ymax=-0.5)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gt <- ggplot_gtable(ggplot_build(m))
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gt)


#1.2. No equilibrium


#Extra text
Desc<-stat.desc(ParametersNoEquilibria$ssigma2)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersNoEquilibria, aes(x = ssigma2))
m<-m + geom_histogram(binwidth = 0.02) 
m<-m+labs(title="No Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=10,ymin=-15,ymax=-10)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gtN <- ggplot_gtable(ggplot_build(m))
gtN$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gtN)



dev.set()
png(file="Figures/ssigma2.png",width=398,height=389)
ggarrange(gtN, gt)
dev.off()




#14. rho12
#1.1 Equilibrium


#Extra text
Desc<-stat.desc(ParametersEquilibria$rho12)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersEquilibria, aes(x = rho12))
m<-m + geom_histogram(binwidth = 0.02) 
m<-m+labs(title="Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=10,ymin=-1.5,ymax=-0.5)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gt <- ggplot_gtable(ggplot_build(m))
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gt)


#1.2. No equilibrium


#Extra text
Desc<-stat.desc(ParametersNoEquilibria$rho12)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersNoEquilibria, aes(x = rho12))
m<-m + geom_histogram(binwidth = 0.02) 
m<-m+labs(title="No Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=10,ymin=-15,ymax=-10)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gtN <- ggplot_gtable(ggplot_build(m))
gtN$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gtN)



dev.set()
png(file="Figures/rho12.png",width=398,height=389)
ggarrange(gtN, gt)
dev.off()






#14.2 c
#1.1 Equilibrium


#Extra text
Desc<-stat.desc(ParametersEquilibria$c)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersEquilibria, aes(x = c))
m<-m + geom_histogram(binwidth = 0.02) 
m<-m+labs(title="Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=10,ymin=-1.5,ymax=-0.5)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gt <- ggplot_gtable(ggplot_build(m))
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gt)


#1.2. No equilibrium


#Extra text
Desc<-stat.desc(ParametersNoEquilibria$c)
Desc<-data.frame(as.list(Desc))
Extratext<-paste("Min:",round(Desc$min,3), "Max:",round(Desc$ma,3),"Mean:",round(Desc$mean,3))



#Plot
m <- ggplot(ParametersNoEquilibria, aes(x = c))
m<-m + geom_histogram(binwidth = 0.02) 
m<-m+labs(title="No Equilibria")
m<-m+annotation_custom(grob=textGrob(Extratext,gp = gpar(fontsize = 8)),xmin=0,xmax=10,ymin=-15,ymax=-10)
m<-m+theme(plot.margin = unit(c(1,1,2,1), "lines"))
gtN <- ggplot_gtable(ggplot_build(m))
gtN$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gtN)



dev.set()
png(file="Figures/c.png",width=398,height=389)
ggarrange(gtN, gt)
dev.off()





#Probit model to analyze the relationship between parameters and equilibria


mod<-lm( Indic~ aalpha+ggamma+ddelta+bbeta+ssigma+kkappa+psi+chi+rho+mmu1+mmu2+ssigma1+ssigma2+rho12+c, Parameters)



print(coef(summary(mod)))
write(print(coef(summary(mod))),file="ProbabilityParameters.txt")
lapply(coef(summary(mod)), write, "ProbabilityParameters.txt", append=TRUE)
coef(summary(mod))[, "Std. Error"]
sink("ProbabilityParameters.txt")
print(coef(summary(mod)))
sink()

#------------------------------#
#Conclusion after the analysis:
#------------------------------#

#1. 


#------------------------------------------------#
#Plot relationship between moments and parameters#
#------------------------------------------------#

setwd('/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolsGenerated/Dim15Secondary/RelationshipParameters')

#1.Parameters and entrepreneurs
ParamPropentrep<-data.table(Parameters,PropEntrepreneurs)
setnames(ParamPropentrep,15,"Entrepreneurs")

#Graphs with regressions

#1. Aalpha
p<-ggplot(data=EverythingEquilibrium,aes(x=aalpha,y=PropEntrepreneurs))+geom_point(color='blue')
p<-p+geom_smooth(method="lm")
dev.set()
png(file="Aalpha.png",width=1600,height=850)
p
dev.off()


#2. Ggamma
p<-ggplot(data=EverythingEquilibrium,aes(x=ggamma,y=PropEntrepreneurs))+geom_point(color='blue')
p<-p+geom_smooth(method="lm")
dev.set()
png(file="Ggama.png",width=1600,height=850)
p
dev.off()

#3. Ddelta
p<-ggplot(data=EverythingEquilibrium,aes(x=ddelta,y=PropEntrepreneurs))+geom_point(color='blue')
p<-p+geom_smooth(method="lm")
dev.set()
png(file="Ddelta.png",width=1600,height=850)
p
dev.off()


#4. Bbeta
p<-ggplot(data=EverythingEquilibrium,aes(x=bbeta,y=PropEntrepreneurs))+geom_point(color='blue')
p<-p+geom_smooth(method="lm")
dev.set()
png(file="bbeta.png",width=1600,height=850)
p
dev.off()



#5. ssigma
p<-ggplot(data=EverythingEquilibrium,aes(x=ssigma,y=PropEntrepreneurs))+geom_point(color='blue')
p<-p+geom_smooth(method="lm")
dev.set()
png(file="ssigma.png",width=1600,height=850)
p
dev.off()


#6. kkappa
p<-ggplot(data=EverythingEquilibrium,aes(x=kkappa,y=PropEntrepreneurs))+geom_point(color='blue')
p<-p+geom_smooth(method="lm")
dev.set()
png(file="kkappa.png",width=1600,height=850)
p
dev.off()

#7. psi
p<-ggplot(data=EverythingEquilibrium,aes(x=psi,y=PropEntrepreneurs))+geom_point(color='blue')
p<-p+geom_smooth(method="lm")
dev.set()
png(file="psi.png",width=1600,height=850)
p
dev.off()


#8. chi
p<-ggplot(data=EverythingEquilibrium,aes(x=chi,y=PropEntrepreneurs))+geom_point(color='blue')
p<-p+geom_smooth(method="lm")
dev.set()
png(file="chi.png",width=1600,height=850)
p
dev.off()


#9. rho
p<-ggplot(data=EverythingEquilibrium,aes(x=rho,y=PropEntrepreneurs))+geom_point(color='blue')
p<-p+geom_smooth(method="lm")
dev.set()
png(file="rho.png",width=1600,height=850)
p
dev.off()


#10. mmu1
p<-ggplot(data=EverythingEquilibrium,aes(x=mmu1,y=PropEntrepreneurs))+geom_point(color='blue')
p<-p+geom_smooth(method="lm")
dev.set()
png(file="mmu1.png",width=1600,height=850)
p
dev.off()


#11. mmu2
p<-ggplot(data=EverythingEquilibrium,aes(x=mmu2,y=PropEntrepreneurs))+geom_point(color='blue')
p<-p+geom_smooth(method="lm")
dev.set()
png(file="mmu2.png",width=1600,height=850)
p
dev.off()


#12. ssigma1
p<-ggplot(data=EverythingEquilibrium,aes(x=ssigma1,y=PropEntrepreneurs))+geom_point(color='blue')
p<-p+geom_smooth(method="lm")
dev.set()
png(file="ssigma1.png",width=1600,height=850)
p
dev.off()


#13. ssigma2
p<-ggplot(data=EverythingEquilibrium,aes(x=ssigma2,y=PropEntrepreneurs))+geom_point(color='blue')
p<-p+geom_smooth(method="lm")
dev.set()
png(file="ssigma2.png",width=1600,height=850)
p
dev.off()


#14. rho12
p<-ggplot(data=EverythingEquilibrium,aes(x=rho12,y=PropEntrepreneurs))+geom_point(color='blue')
p<-p+geom_smooth(method="lm")
dev.set()
png(file="rho12.png",width=1600,height=850)
p
dev.off()





lm( PropEntrepreneurs~ aalpha+bbeta+ddelta+ggamma+kkappa+ssigma+rho+ssigma1+ssigma2+rho12, EverythingEquilibrium)


mod <- lm(PropEntrepreneurs ~ aalpha + 
            ggamma +
            ddelta+
            bbeta+
            ssigma+
            kkappa+
            psi+
            chi+
            rho+
            mmu1+
            mmu2+
            ssigma1+
            ssigma2+
            rho12
            ,EverythingEquilibrium)

write(print(coef(summary(mod))),file="Regression.txt")



#2. Moments of informality. 


setwd('/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolsGenerated/Dim15Secondary/Informality')


#Graphs with regressions

#1. Aalpha
p<-ggplot(data=EverythingEquilibrium,aes(x=aalpha,y=InformalDemandProportionV4))+geom_point(color='blue')
p<-p+geom_smooth(method="lm")
dev.set()
png(file="Aalpha.png",width=1600,height=850)
p
dev.off()


#2. Ggamma
p<-ggplot(data=EverythingEquilibrium,aes(x=ggamma,y=InformalDemandProportionV4))+geom_point(color='blue')
p<-p+geom_smooth(method="lm")
dev.set()
png(file="Ggama.png",width=1600,height=850)
p
dev.off()

#3. Ddelta
p<-ggplot(data=EverythingEquilibrium,aes(x=ddelta,y=InformalDemandProportionV4))+geom_point(color='blue')
p<-p+geom_smooth(method="lm")
dev.set()
png(file="Ddelta.png",width=1600,height=850)
p
dev.off()


#4. Bbeta
p<-ggplot(data=EverythingEquilibrium,aes(x=bbeta,y=InformalDemandProportionV4))+geom_point(color='blue')
p<-p+geom_smooth(method="lm")
dev.set()
png(file="bbeta.png",width=1600,height=850)
p
dev.off()



#5. ssigma
p<-ggplot(data=EverythingEquilibrium,aes(x=ssigma,y=InformalDemandProportionV4))+geom_point(color='blue')
p<-p+geom_smooth(method="lm")
dev.set()
png(file="ssigma.png",width=1600,height=850)
p
dev.off()


#6. kkappa
p<-ggplot(data=EverythingEquilibrium,aes(x=kkappa,y=InformalDemandProportionV4))+geom_point(color='blue')
p<-p+geom_smooth(method="lm")
dev.set()
png(file="kkappa.png",width=1600,height=850)
p
dev.off()

#7. psi
p<-ggplot(data=EverythingEquilibrium,aes(x=psi,y=InformalDemandProportionV4))+geom_point(color='blue')
p<-p+geom_smooth(method="lm")
dev.set()
png(file="psi.png",width=1600,height=850)
p
dev.off()


#8. chi
p<-ggplot(data=EverythingEquilibrium,aes(x=chi,y=InformalDemandProportionV4))+geom_point(color='blue')
p<-p+geom_smooth(method="lm")
dev.set()
png(file="chi.png",width=1600,height=850)
p
dev.off()


#9. rho
p<-ggplot(data=EverythingEquilibrium,aes(x=rho,y=InformalDemandProportionV4))+geom_point(color='blue')
p<-p+geom_smooth(method="lm")
dev.set()
png(file="rho.png",width=1600,height=850)
p
dev.off()


#10. mmu1
p<-ggplot(data=EverythingEquilibrium,aes(x=mmu1,y=InformalDemandProportionV4))+geom_point(color='blue')
p<-p+geom_smooth(method="lm")
dev.set()
png(file="mmu1.png",width=1600,height=850)
p
dev.off()


#11. mmu2
p<-ggplot(data=EverythingEquilibrium,aes(x=mmu2,y=InformalDemandProportionV4))+geom_point(color='blue')
p<-p+geom_smooth(method="lm")
dev.set()
png(file="mmu2.png",width=1600,height=850)
p
dev.off()


#12. ssigma1
p<-ggplot(data=EverythingEquilibrium,aes(x=ssigma1,y=InformalDemandProportionV4))+geom_point(color='blue')
p<-p+geom_smooth(method="lm")
dev.set()
png(file="ssigma1.png",width=1600,height=850)
p
dev.off()


#13. ssigma2
p<-ggplot(data=EverythingEquilibrium,aes(x=ssigma2,y=InformalDemandProportionV4))+geom_point(color='blue')
p<-p+geom_smooth(method="lm")
dev.set()
png(file="ssigma2.png",width=1600,height=850)
p
dev.off()


#14. rho12
p<-ggplot(data=EverythingEquilibrium,aes(x=rho12,y=InformalDemandProportionV4))+geom_point(color='blue')
p<-p+geom_smooth(method="lm")
dev.set()
png(file="rho12.png",width=1600,height=850)
p
dev.off()





lm( InformalDemandProportionV4~ aalpha+bbeta+ddelta+ggamma+kkappa+ssigma+rho+ssigma1+ssigma2+rho12, EverythingEquilibrium)


mod <- lm(InformalDemandProportionV4 ~ aalpha + 
            ggamma +
            ddelta+
            bbeta+
            ssigma+
            kkappa+
            psi+
            chi+
            rho+
            mmu1+
            mmu2+
            ssigma1+
            ssigma2+
            rho12
          ,EverythingEquilibrium)

write(print(coef(summary(mod))),file="Regression.txt")

}
















