rm(list=ls(all=TRUE))
#install.packages("randtoolbox")
library(randtoolbox)
#1. Run empirical moments:
source('/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/OptimalTaxation/EmpiricalMoments/Momentos.R')

#2. Run functions needed to generate theoretical moments 
source('/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/OptimalTaxation/TheoreticalMoments/Equilibrium.R')


#3. Generating the theoretical moments. Function
TheoMoments<-function(ParamsDecisionExcessDemand,WagesInitialGuess){
  
  #1. Finding wages of equilibrium 
  WEq=EqWagesNumericVector(ParamsDecisionExcessDemand,WagesInitialGuess)
  wi=WEq[1]
  wf=WEq[2]
  ExcessDemand=WEq[3]
  
  #2. Loading parameters correspondingly
  aalpha<-ParamsDecisionExcessDemand[1]
  ddelta<-ParamsDecisionExcessDemand[2]
  ggamma<-ParamsDecisionExcessDemand[3]
  bbeta<-ParamsDecisionExcessDemand[4]
  ssigma<-ParamsDecisionExcessDemand[5]
  kkappa<-ParamsDecisionExcessDemand[6]
  rrho<-ParamsDecisionExcessDemand[7]
  psi<-ParamsDecisionExcessDemand[8]
  chi<-ParamsDecisionExcessDemand[9]
  mmu1<-ParamsDecisionExcessDemand[10]
  mmu2<-ParamsDecisionExcessDemand[11]
  ssigma1<-ParamsDecisionExcessDemand[12]
  ssigma2<-ParamsDecisionExcessDemand[13]
  rho12<-ParamsDecisionExcessDemand[14]
  li<-ParamsDecisionExcessDemand[15]
  lf<-ParamsDecisionExcessDemand[16]
  ni<-ParamsDecisionExcessDemand[17]
  nf<-ParamsDecisionExcessDemand[18]
  nf<-ParamsDecisionExcessDemand[19]
  
  
  #3. Calculating number of entrepreneurs
  Sigma <- matrix(c(ssigma1,rho12,rho12,ssigma2),2,2)
  mu=c(mmu1,mmu2)
  set.seed(257)
  
  logtthetavec<-rmvnorm(100000, mean = mu, Sigma)
  tthetaw<-exp(logtthetavec[,1])
  tthetae<-exp(logtthetavec[,2])
  mean(tthetaw)
  mean(tthetae)
  var(tthetaw)
  var(tthetae)
  minTtw<-min(tthetaw)
  minTte<-min(tthetae)
  max(tthetaw)
  max(tthetae)
  
  #First, doing the analysis of who works and who doesn't. 
  #Define number of bins you want in each case. 

  
  
  
  
  params<-c(wiEq,wfEq,aalpha,ddelta,ggamma,bbeta,ssigma,kkappa,rrho,psi,chi)
  
  
  tthetavec<-c(2,1)
  
  
  
  
  
  #4. Calculating moments of workers and entrepreneurs
  set.seed(2581633)
  logtthetavec<-rmvnorm(100, mean = mu, Sigma)
  tthetae_Sample<-sort(exp(logtthetavec[,2]))
  tthetaw_Sample<-sort(exp(logtthetavec[,1]))
  le<-length(tthetae_Sample)
  lw<-length(tthetaw_Sample)
  
  
  DecisionMatrix<-matrix(0,lw,le)
  x.m <- melt(t(DecisionMatrix))
  DecisionMatrixVer=matrix(0,le*lw,3)
  
  it=1
  for (tte in 1:lw){
    for (ttw in 1:le){
      #Generating the ttheta vector
      tthetavec<-c(tthetaw_Sample[ttw],tthetae_Sample[tte])
      DecisionMatrix[ttw,tte]<-iDecision(tthetavec,params,InitLWorkers,InitProf)$Decission
      DecisionMatrixVer[it,1]=tthetae_Sample[tte]
      DecisionMatrixVer[it,2]=tthetaw_Sample[ttw]
      DecisionMatrixVer[it,3]=iDecision(tthetavec,params,InitLWorkers,InitProf)$Decission
      
      #Identifying ranges of entrepreneurs and workers
      #which(DecisionMatrix==1)
      #which(DecisionMatrix==0)
      #print(it)
      it=it+1
      
    }
  }
  
  #Number of entrepreneurs:
  
  logtthetavec<-rmvnorm(100000, mean = mu, Sigma)
  EntrepV<-numeric(100000)
  for(i in 1:10000){
    EntrepV[i]=iDecision(exp(logtthetavec[i,]),params,InitLWorkers,InitProf)$Decission
  }
  PropEntrep<-sum(EntrepV)/length(EntrepV)
  
  EntreP<-100*PropEntrep
  
  
  #Identifying the range of workers and entrepreneurs
  
  logtthetavec<-rmvnorm(10000, mean = mu, Sigma)
  tthetae_Sample<-sort(exp(logtthetavec[,2]))
  tthetaw_Sample<-sort(exp(logtthetavec[,1]))
  le<-length(tthetae_Sample)
  lw<-length(tthetaw_Sample)
  
  
  Dec=0
  tte=0
  while(Dec==0){
    tte<-tte+1
    tthetavec<-c(minTtw,tthetae_Sample[tte])
    Dec<-iDecision(tthetavec,params,InitLWorkers,InitProf)$Decission
  }
  minEnt<-tthetae_Sample[tte]
  maxEnt<-max(tthetae_Sample)
  
  #Identifying the range of workers and entrepreneurs
  Dec=1
  ttw=0
  while(Dec==1){
    ttw<-ttw+1
    tthetavec<-c(tthetaw_Sample[ttw],minTte)
    Dec<-iDecision(tthetavec,params,InitLWorkers,InitProf)$Decission
  }
  minWork<-tthetaw_Sample[ttw]
  maxWork<-max(tthetaw_Sample)
  #Their corresponding pdf
  tthetae_Dist<-pnorm(log(tthetae_Sample),mean=mu[2],sd=sqrt(Sigma[2,2]))
  
  
  #Normalization of the truncated distribution
  Zentrep=pnorm(log(maxEnt),mean=mu[2],sd=sqrt(Sigma[2,2]))-
                  pnorm(log(minEnt),mean=mu[2],sd=sqrt(Sigma[2,2]))
  
  PPHI_MINENTREP=pnorm(log(minEnt),mean=mu[2],sd=sqrt(Sigma[2,2]))
  
  Zoptimal<-numeric(le)
  InformalDemand<-numeric(le)
  FormalDemand<-numeric(le)
  PretaxProfit<-numeric(le)
  Production<-numeric(le)
  Zproportion<-numeric(le)
  Zproportion2<-numeric(le)
  TotalLaborForce<-numeric(le)
  InformalProportion<-numeric(le)
  AfterTaxProfit<-numeric(le)
  Taxpayed<-numeric(le)
  TaxSales<-numeric(le)
  TaxProfits<-numeric(le)
  tthetae_Trunc<-numeric(le)
  
  
  
  
  for (tte in 1:le){
    tthetae=tthetae_Sample[tte]
    tthetavec<-c(0.0,tthetae)
    
    #Decision =1 if entrepreneur
    Decision=iDecision(tthetavec,params,InitLWorkers,InitProf)$Decission
    
    #Optimal evasion levels
    zevasion=iDecision(tthetavec,params,InitLWorkers,InitProf)$Zoptimal
    Zoptimal[tte]=zevasion
    
    #Demand of informal labor
    ni=iDecision(tthetavec,params,InitLWorkers,InitProf)$InformalDemand
    InformalDemand[tte]=ni
    
    #Demand of formal labor
    nf=iDecision(tthetavec,params,InitLWorkers,InitProf)$FormalDemand
    FormalDemand[tte]=nf
    
    #Total Labor force
    TotalLaborForce[tte]=ni+nf
    
    #Informal proportion
    InformalProportion[tte]=ni/(ni+nf)
    
    #Pre tax profits
    prod=profm(ni,nf,aalpha,tthetae,wi,wf,z)
    PretaxProfit[tte]=profm(ni,nf,aalpha,tthetae,wi,wf,z)
    
    #Production level
    Production[tte]<-tthetae*((ni+nf)^(aalpha))
    
    #After Tax Profit
    AfterTaxProfit[tte]=PretaxProfit[tte]-TcActual(zevasion,ni,nf,aalpha,tthetae,wi,wf)
    
    #Evasion as proportion of profits
    Zproportion[tte]=zevasion/Production[tte]
    Zproportion2[tte]=zevasion/PretaxProfit[tte]
    
    FinProfits(c(ni,nf,zevasion),c(wi,wf,aalpha,ddelta,ggamma,tthetae,bbeta,ssigma))
    
    
    
    
    #Taxes payed
    Taxpayed[tte]=TcActual(zevasion,ni,nf,aalpha,tthetae,wi,wf)
    
    
    #Taxes payed as a proportion of production
    TaxSales[tte]=Taxpayed[tte]/Production[tte]
    
    #Taxes payed as proportion of benefits
    TaxProfits[tte]=Taxpayed[tte]/PretaxProfit[tte]
    
    #PDf of the truncated distribution
    tthetae_Trunc[tte]<-(tthetae_Dist[tte]-PPHI_MINENTREP)/(Zentrep)
  }
  
  
  #Total revenue
  TaxpayedProp<-Taxpayed/sum(Taxpayed)
  
  
  #Ploting the corresponding relationships
  zevasion1<-as.data.frame(cbind(tthetae_Sample,tthetae_Dist,Zoptimal,InformalDemand,
                                 FormalDemand,TotalLaborForce,InformalProportion,PretaxProfit,Zproportion,
                                 Production,AfterTaxProfit,Taxpayed,TaxSales,TaxProfits,Zproportion2,
                                 tthetae_Trunc))
  
  
  
  #Keeping only active entrepreneurs
  zevasion1<-subset(zevasion1,tthetae_Sample>=minEnt)
  zevasion1<-subset(zevasion1,tthetae_Sample<=maxEnt)
  
  
  #Obtaining percentiles
  perc<-seq(0.1,0.9,0.1)
  length_perc<-length(perc)
  
  A1<-subset(zevasion1,tthetae_Trunc==quantile(tthetae_Trunc,c(perc[1]),type=3))
  
  for(p in 2:length_perc){
    At<-subset(zevasion1,tthetae_Trunc==quantile(tthetae_Trunc,c(perc[p]),type=3))
    A1<-rbind(A1,At)
  }
  
  #Keeping the original
  zevasionALL<-zevasion1
  
  #If want to do moments based on percentiles, not the whole data, run the following line:
  zevasion1<-A1

  
  #Generating the tax burden for each percentile
  TaxTotalPayment<-sum(zevasion1$Taxpayed)
  TaxTotalProportion<-zevasion1$Taxpayed/TaxTotalPayment
  zevasion1<-cbind(zevasion1,TaxTotalProportion)
  #-------------------------#
  #Analysis   of   workers  #
  #-------------------------#
  
  #Draws of the worker skill
  logtthetavec<-rmvnorm(1000, mean = mu, Sigma)
  tthetaw_Sample<-sort(exp(logtthetavec[,1]))
  lw<-length(tthetaw_Sample)
  
  
  
  #Pdf of the distribution
  tthetaw_Dist<-pnorm(log(tthetaw_Sample),mean=mu[1],sd=sqrt(Sigma[1,1]))
  
  
  #Normalization for the truncated distribution
  Zworker=pnorm(log(maxWork),mean=mu[1],sd=sqrt(Sigma[1,1]))-pnorm(log(minWork),mean=mu[1],sd=sqrt(Sigma[1,1]))
  
  #CDF of the min:
  PPHI_MIN=pnorm(log(minWork),mean=mu[1],sd=sqrt(Sigma[1,1]))
  
  
  InformalSupply<-numeric(lw)
  FormalSupply<-numeric(lw)
  ValueWorker<-numeric(lw)
  TotalLaborSupply<-numeric(lw)
  InformalProportion<-numeric(lw)
  InformalIncome<-numeric(lw)
  FormalIncome<-numeric(lw)
  TotalIncome<-numeric(lw)
  tthetaw_Trunc<-numeric(lw)
  
  
  
  for (ttw in 1:lw){
    tthetaw=tthetaw_Sample[ttw]
    tthetavec<-c(tthetaw,1)
    #Decision
    Decision=iDecision(tthetavec,params,InitLWorkers,InitProf)$Decission
    
    #Informal Supply
    InformalSupply[ttw]=iDecision(tthetavec,params,InitLWorkers,InitProf)$InformalSupply
    
    #Formal Supply
    FormalSupply[ttw]=iDecision(tthetavec,params,InitLWorkers,InitProf)$FormalSupply
    
    #Value of worker
    ValueWorker[ttw]=iDecision(tthetavec,params,InitLWorkers,InitProf)$ValueWorker
    
    #Total Labor force
    TotalLaborSupply[ttw]=FormalSupply[ttw]+InformalSupply[ttw]
    
    #Informal proportion
    InformalProportion[ttw]=InformalSupply[ttw]/TotalLaborSupply[ttw]
    
    #Truncated pdf
    tthetaw_Trunc[ttw]=(tthetaw_Dist[ttw]-PPHI_MIN)/Zworker
    
    #Labor income from informal activities
    InformalIncome[ttw]=wi*tthetaw*InformalSupply[ttw]
    
    #Labor income from formal activities
    FormalIncome[ttw]=wf*tthetaw*FormalSupply[ttw]
    
    #Total Income
    TotalIncome[ttw]= FormalIncome[ttw]+ InformalIncome[ttw]
    
    
  }
  
  
  
  
  #Ploting the corresponding relationships
  Worker<-as.data.frame(cbind(tthetaw_Sample,InformalSupply,FormalSupply,ValueWorker,
                              TotalLaborSupply,InformalProportion,tthetaw_Dist,tthetaw_Trunc,
                              InformalIncome,FormalIncome,TotalIncome))
  
  
  #Keeping only the relevant workers
  Worker<-subset(Worker,tthetaw_Sample>=minWork)
  Worker<-subset(Worker,tthetaw_Sample<=maxWork)
  
  #Obtaining percentiles
  perc<-seq(0.1,0.9,0.1)
  length_perc<-length(perc)
  
  A1<-subset(Worker,tthetaw_Trunc==quantile(tthetaw_Trunc,c(perc[1]),type=3))
  
  for(p in 2:length_perc){
    At<-subset(Worker,tthetaw_Trunc==quantile(tthetaw_Trunc,c(perc[p]),type=3))
    A1<-rbind(A1,At)
  }
  
  
  #Keeping the original
  WorkerALL<-Worker
  #Keeping the corresponding deciles:
  Worker<-A1
  
  
  #5. Finally, generate the output of the matrix
  TheoMoments<-cbind(zevasion1,Worker,EntreP,ExcessDemand)
  return(TheoMoments)
}



#Generating function of distance between empirical and theoretical moments
DistanceMoments<-function(ParamsDecisionExcessDemand){
  
  #Obtain the theoretical moments
  ThMoments<-TheoMoments(ParamsDecisionExcessDemand,WagesInitialGuess)



  
  #Empirical Moments
  
  #Number of workers by sales deciles
  M2
  
  #Proportion of workers and entrepreneurs
  d1<-((M4[2,2]/(M4[2,1]+M4[2,2]))-ThMoments$EntreP[1])^2
  
  
  #Income from work
  #Sub- muestra: Trabajadores empleados
  PropIncomeTheoretical<-as.numeric(M5.3A[4:12,2])
  PropIncomeEmpirical<-ThMoments$TotalIncome
  d2<-sum((PropIncomeTheoretical-PropIncomeEmpirical)^2)
  
  #Trabajadores informales. Ingreso
  M5.3B
  #Trabajadores formales. Ingreso
  M5.3C
  
  #Distribucion ventas firma
  #PropSalesTheoretical<-as.numeric(M6[4:12,2])/sum(as.numeric(M6[4:12,2]))
  #PropSalesEmpirical<-ThMoments$Production/sum(ThMoments$Production)
  
  
  PropSalesEmpirical<-as.numeric(M6[4:12,2])
  PropSalesTheoretical<-ThMoments$Production*100
  
  d3<-sum((PropSalesTheoretical-PropSalesEmpirical)^2)
  
  #plot(x,PropSalesTheoretical,type="l",col="red")
  #lines(x,PropSalesEmpirical,col="blue")
  
  #Impuestos
  MOMENTO7A$PropPagoImpuestos=MOMENTO7A$Impuestos/sum(MOMENTO7A$Impuestos)
  d4<-sum((MOMENTO7A$PropPagoImpuestos[1:9]-ThMoments$TaxTotalProportion)^2)
  
  
  #Informalidad por nivel de ingresos de los trabajadores
  d5<-sum((MOMENTO8$`Informalidad (%)`[1:9]/100-ThMoments$InformalProportion)^2)
  
  #Informalidad y tamaño de la firma
  EmpiricInformSize<-T15$`% acumulado del total de firmas`[1:3]/100
  TheoInformSize<-cbind(ThMoments$InformalProportion[4],ThMoments$InformalProportion[7],ThMoments$InformalProportion[8])
  d6<-sum((EmpiricInformSize-TheoInformSize)^2)
  
  
  #Labor supply. Decil como propoción de la media. 
  TheoLaborSupply<-log(ThMoments$TotalLaborSupply )/log(ThMoments$TotalLaborSupply )[5]
  EmpLaborSupply<-T19$horas/T19$horas[5]
  EmpLaborSupply<-EmpLaborSupply[1:9]
  d7<-sum((TheoLaborSupply-EmpLaborSupply)^2)
  
  #Weights
  distance<-d1+d2+d3+d4+d5+d6+d7
  
  DistanceV<-rep(distance,9)
  ExcessDemand<-rep(ThMoments$ExcessDemand)
  Ans<-cbind(ThMoments,DistanceV,ExcessDemand)
  return(Ans)
}


#Sobol calibration
aalpha=0.8
wi=8
wf=8
ni=2.3
nf=0.53*70
ggamma=0.28
ddelta=0.12
bbeta=0.15
ssigma=0.2
kkappa=0.2
psi=0.4
chi=2.5
rrho=0.9
z=24
li=2.1
lf=2.1
mmu1<-0.2
mmu2<-1.3
ssigma1<-0.1
ssigma2<-0.9
rho12<-0.08




#13 sobol between 0 and 1. Need to transform correctly:
Rand<-sobol(1000,dim=13,seed=2581633)
Rand[,1]<-Rand[,1]*5+0.1 #ggamma between 0.1 and 2
Rand[,2]<-Rand[,2]*5+0.1 #ddelta between 0.1 and 5
Rand[,3]<-Rand[,3]*5+0.1 #bbeta  between 0.1 and 2
Rand[,4]<-Rand[,4]*5+0.1 #ssigma  between 0.1 and 2
Rand[,5]<-Rand[,5]*5+0.1 #Kappa between 0.1 and 2
Rand[,6]<-Rand[,6]*5+0.1 #Psi between 0.1 and 2
Rand[,7]<-Rand[,7]*5+0.1 #chi between 0.1 and 4
Rand[,8]<-Rand[,8]*5+0.1 #rrho between 0.1 and 2
Rand[,9]<-Rand[,9]*5+0.5 #Mmu1 between 0.5 and 3
Rand[,10]<-Rand[,10]*5 +0.5#mmu3 between 0.5 and 3. 
Rand[,11]<-Rand[,11]*5+0.5 #ssigma1 between 0.5 and 3
Rand[,12]<-Rand[,12]*5+0.5 #ssigma2 between 0.5 and 3
Rand[,13]<-Rand[,13]*2+0.5   #rho12 between 0 and 1 

#Doing the distance estimator sobol
TotalDistance<-data.matrix(read.csv(file = "/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/OptimalTaxation/TheoreticalMoments/MOM2.csv"))

#Problem with 61-63 why?



for(i in 1:500){
  Par<-seq(1,19)
  Par[1]<-0.8
  Par[2]<-Rand[i,1]
  Par[3]<-Rand[i,2]
  Par[4]<-Rand[i,3]
  Par[5]<-Rand[i,4]
  Par[6]<-Rand[i,5]
  Par[7]<-Rand[i,6]
  Par[8]<-Rand[i,7]
  Par[9]<-Rand[i,8]
  Par[10]<-Rand[i,9]
  Par[11]<-Rand[i,10]
  Par[12]<-Rand[i,11]
  Par[13]<-Rand[i,12]
  Par[14]<-Rand[i,13]
  Par[15]<-li
  Par[16]<-lf
  Par[17]<-ni
  Par[18]<-nf
  Par[19]<-z
  WagesInit<-WagesInitialGuess
  InitLWorkers
  DistanceEstimator(Par,WagesInit,InitLWorkers,InitProf)
}



testc<-seq(1,19)
testc[1]=0.8
testc[2]=0.396875
testc[3]=0.177344
testc[4]=0.157031
testc[5]=0.459375
testc[6]=1.27188
testc[7]=2.08438
testc[8]=0.170313
testc[9]=1.74063
testc[10]=3.10156
testc[11]=3.47656
testc[12]=0.664062
testc[13]=1.83594


testc[15]<-li
testc[16]<-lf
testc[17]<-ni
testc[18]<-nf
testc[19]<-z
DistanceEstimator(testc,WagesInit,InitLWorkers,InitProf)


Total<-data.matrix(read.csv(file = "/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/OptimalTaxation/TheoreticalMoments/MOM2.csv"))

for(i in 360:400){
  Par<-seq(1,19)
  Par[1]<-aalpha
  Par[2]<-Rand[i,1]
  Par[3]<-Rand[i,2]
  Par[4]<-Rand[i,3]
  Par[5]<-Rand[i,4]
  Par[6]<-Rand[i,5]
  Par[7]<-Rand[i,6]
  Par[8]<-Rand[i,7]
  Par[9]<-Rand[i,8]
  Par[10]<-Rand[i,9]
  Par[11]<-Rand[i,10]
  Par[12]<-Rand[i,11]
  Par[13]<-Rand[i,12]
  Par[14]<-Rand[i,13]
  Par[15]<-li
  Par[16]<-lf
  Par[17]<-ni
  Par[18]<-nf
  Par[19]<-z
  Calibration<-DistanceMoments(Par)
  #Mom[1:9,1:30,i]<-data.matrix(Calibration)
  MOM2=data.matrix(Calibration)
  if(i==1){
    Total=MOM2
  }else{
    Total<-data.matrix(read.csv(file = "/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/OptimalTaxation/TheoreticalMoments/MOM2.csv"))
    Total=cbind(Total,MOM2)
    
  }
  write.csv(Total, file = "/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/OptimalTaxation/TheoreticalMoments/MOM2.csv",row.names=FALSE)
  rm(Total)
}


                 
#Proportion of entrepreneurs is in:   
dim(Total)[2]/32



Total[1:9,1:32]
Total[1:9,33:64]
Total[1:9,29]
Total[1:9,61]
Total[1:9,92:93]

Total$EntreP.1
TTT<-as.data.frame.table(Total)
TTT$EntreP1

                 
                 
                 