#install.packages('nleqslv')
rm(list=ls(all=TRUE))
library(nleqslv)
library(lattice)
library(ggplot2)
library(reshape)
#Compile Rcpp Functions
library(Rcpp)
library(RcppArmadillo)
library(mvtnorm)
Rcpp::sourceCpp('/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/main.cpp')
GRAPHS="/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/OptimalTaxation/TheoreticalMoments"



#------------#
#Housekeeping#
#------------#

#setwd("C:/Users/razuero/Documents/Informalidad/OptimalTaxation/Data")
  
  
  #Parameter definition
  aalpha=0.8
  wi=7.2
  wf=7.59
  ni=2.3
  nf=0.53*70
  ggamma=0.28
  ddelta=0.12
  bbeta=0.15
  ssigma=0.2
  kkappa=0.1
  psi=0.4
  chi=1.5
  rrho=0.9
  z=24
  li=2.1
  lf=2.1
  mmu1<-1.1
  mmu2<-2
  ssigma1<-0.5
  ssigma2<-1.1
  rho12<-0.3
  #Obtaining the mean ttheta:
  Sigma <- matrix(c(ssigma1,rho12,rho12,ssigma2),2,2)
  mu=c(mmu1,mmu2)
  set.seed(257)
  logtthetavec<-rmvnorm(100000, mean = mu, Sigma)
  tthetaw<-exp(logtthetavec[,1])
  tthetae<-exp(logtthetavec[,2])
  tthetae<-mean(tthetae)
  tthetae
  tthetaw<-mean(tthetaw)
  tthetaw
  min(exp(logtthetavec[,1]))
  
  
  
  
  #Auxiliar function definitions
  
  #Payroll taxes marginal
  Tn<-function(nf){
    return(0.09)
  }
  
  #Payroll taxes actual
  TnActual<-function(nf){
    term1=-1/(1+nf^2)
    term2=term1+1
    term3=0.45*term2
    ans=term3+0.05
    return(ans)
  }
  
  #Payroll taxes actual
  TnActual<-function(nomina){
    term1=0.09*nomina
    return(term1)
  }
  
  #PreTaxProfits
  profm<-function(ni,nf,aalpha,tthetae,wi,wf,z){
    ppi<-tthetae*(ni+nf)^aalpha-wi*ni-wf*nf-TnActual(nf*wf)
    return(ppi)
  }
  
  #Corporate taxes marginal
  Tc<-function(z,ni,nf,aalpha,tthetae,wi,wf){
    profm<-tthetae*(ni+nf)^aalpha-wi*ni-wf*nf-TnActual(nf*wf)
    arg<-profm-z
    firsterm=-exp(-arg)
    return(max(firsterm+1,0))
  }
  
  
  #Tc Actual corporate taxes
  #Tc Actual corporate taxes
  TcActual<-function(z,ni,nf,aalpha,tthetae,wi,wf){
    profm<-tthetae*(ni+nf)^aalpha-wi*ni-wf*nf-TnActual(nf*wf)
    prod<-tthetae*(ni+nf)^aalpha
    if(prod-z<=0){
      ans=0
    }else if (prod-z<189){
      ans=0.756
    } else if (prod-z<302.4){
      ans=1.89
    } else if (prod-z<491.4){
      ans=7.56
    } else if (prod-z<756){
      ans=15.12
    } else if (prod-z<1134){
      ans=22.68
    }else if (prod-z>=1134){
      arg<-profm-z
      tax<-0.3
      ans=arg*tax
      ans=max(ans,0)
    }
    return(ans)
  }
  #Personal Income Taxes marginal
  PIT<-function(tthetaw,wf,lf){
    taxburden<-tthetaw*wf*lf
    taxrate<-0.3
    return(taxrate*taxburden)
  }
  
  #Personal Income Taxes marginal
  PIT<-function(tthetaw,wf,lf){
    taxburden<-tthetaw*wf*lf
    taxrate<-taxburden^0.015-1
    return(taxrate*taxburden)
  }
  
  
  PIT<-function(tthetaw,wf,lf){
    x<-tthetaw*wf*lf
    a=1000
    if(x<a){
      ans=(100/a)*x-100
    }
    else if(x<36000){
      ans=0
    }
    else if(x<66000){
      ans=(x^2)/100000-9*x/25
    }
    else if(x>=66000){
      ans=0.3*x
    }
    return(ans)
  }
  
  #Total profits
  FinProfits<-function(Args,params){
    #0. Loading arguments and parameters
    
    
    wi=params[1]
    wf=params[2]
    aalpha=params[3]
    ddelta=params[4]
    ggamma=params[5]
    tthetae=params[6]
    bbeta=params[7]
    ssigma=params[8]
    
    ni<-Args[1] #If need to set positive 
    nf<-Args[2]
    z<-Args[3]
    
    #First, operational profits
    term1<-profm(ni,nf,aalpha,tthetae,wi,wf,z)
    
    #Second, corp. taxes
    taxes<-TcActual(z,ni,nf,aalpha,tthetae,wi,wf)
    
    #Third, costs of informality
    infcost<-(ddelta/(1+ggamma))*(ni^(1+ggamma))
    
    #Fourth, evasion cost
    evcost<-(bbeta/(1+ssigma))*(z^(1+ssigma))
    
    #Compute the answer
    ans<-term1-taxes-infcost-evcost
    return(ans)
  }
  
  params=c(wi,wf,aalpha,ddelta,ggamma,tthetae,bbeta,ssigma)
  Args=c(ni,nf,z)
  FinProfits(Args,params)
  
  
  
  
  #Given parameters, find profits of individual through maximizing the objective function
  profitsFinMaxim<-function(Init,params){
    
    #Define profit function in terms of ni, nf, z only:
    ProfitoptFunc<-function(X){
      return(-FinProfits(exp(X),params))
    }
    
    #Maximizing function
    F<-optim(log(Init),ProfitoptFunc)
    ans<-c(0,0,0,0)
    ans[1]<-exp(F$par[1])
    ans[2]<-exp(F$par[2])
    ans[3]<-exp(F$par[3])
    ans[4]<--F$value
    return(ans)
  }
  
  
  #Testing this final function through FOC
  
  params=c(wi,wf,aalpha,ddelta,ggamma,tthetae,bbeta,ssigma)
  Init<-c(ni, nf, 0.1)
  FinProfits(c(0.004,0.009,0),params)
  TestProfOpt<-profitsFinMaxim(Init,params)
  TestProfOpt
  
  
  #Now obtaining the value of workers
  ValueWorkers<-function(Args,ParamsWorkers){
    wf<-ParamsWorkers[1]
    wi<-ParamsWorkers[2]
    kkappa<-ParamsWorkers[3]
    rrho<-ParamsWorkers[4]
    psi<-ParamsWorkers[5]
    chi<-ParamsWorkers[6]
    tthetaw<-ParamsWorkers[7]
    #li=0.1929058 //Found using max
    #li=0.4705325 //Found using derivative
    li<-Args[1]
    lf<-Args[2]
    
    
    #Term 1 is the total income derived from both informal and formal activities
    term1<-(wf*lf+wi*li)*tthetaw
    
    #Term 2 is the total disutility from providing labor supply of both types
    term2<-(chi/(1+psi))*((li+lf)^(1+psi))
    
    #Term 3 is the consumption penalty from the informal labor income
    term3<-(kkappa/(1+rrho))*((tthetaw*li)^(1+rrho))
    
    #Term 4 is the personal income tax
    Taxes<-PIT(tthetaw,wf,lf)
    
    ans<-term1-term2-term3-Taxes
    return(ans)
  }
  
  
  #Testing the function
  ParamsWorkers=c(wf,wi,kkappa,rrho,psi,chi,tthetaw)
  Args=c(li,lf)
  ValueWorkers(Args,ParamsWorkers)
  
  #Setting to maximize the function
  ValueMax<-function(InitL,ParamsWorkers){
    
    #Defining the function to maximize
    MaxValueWorker<-function(X){
      return(-ValueWorkers(exp(X),ParamsWorkers))
    }
    
    #Maximizing the function
    F<-optim(log(InitL),MaxValueWorker,method="BFGS")
    ans<-c(0,0,0)
    ans[1]<-exp((F$par[1]))
    ans[2]<-exp(F$par[2])
    ans[3]<--F$value
    return(ans)
  }
  
  
  
  ParamsWorkers=c(wf,wi,kkappa,rrho,psi,chi,tthetaw)
  InitL=c(li,lf)
  ValueMax(InitL,ParamsWorkers)
  F<-ValueMax(InitL,ParamsWorkers)
  F
  
  
  
  
  
  InitLWorkers<-c(li,lf)
  InitProf<-c(ni,nf,z)
  Ttheta<-c(tthetaw,tthetae)
  params<-c(wi,wf,aalpha,ddelta,ggamma,bbeta,ssigma,kkappa,rrho,psi,chi) 
  
  
  #Function defining if an individual is worker or an entrepreneur
  iDecision<-function(Ttheta,params,InitLWorkers,InitProf){
    
    #Loading parameters
    wi=params[1]
    wf=params[2]
    aalpha=params[3]
    ddelta=params[4]
    ggamma=params[5]
    bbeta=params[6]
    ssigma=params[7]
    kkappa=params[8]
    rrho=params[9]
    psi=params[10]
    chi=params[11]
    
    
    #Loading ttheta
    tthetaw<-Ttheta[1]
    tthetae<-Ttheta[2]
    
    #Matrix Dec will be the output. Entry [1,1] will be decission
    #Second column devoted to problem of worker (Informal supply, formal supply, value)
    #Third column is problem of entrepreneur (Formal demmand, informal demmand)
    Dec<-list()
    
    #Computing the value of workers
    #First, define the parameters to be used as input in the value of the workers
    ParamsWorkers<-c(wf,wi,kkappa,rrho,psi,chi,tthetaw)
    
    
    Worker<-ValueMax(InitLWorkers,ParamsWorkers)
    Dec$InformalSupply<-Worker[1]
    Dec$FormalSupply<-Worker[2]
    Dec$ValueWorker<-Worker[3]
    ValueWorker<-Worker[3]
    
    
    #Computing maximum possible profits
    #Defining the parameters of the entrepreneur
    ParamsEntrep=c(wi,wf,aalpha,ddelta,ggamma,tthetae,bbeta,ssigma)
    
    
    
    
    ni<-InitProf[1] #If need to set positive 
    nf<-InitProf[2]
    z<-InitProf[3]
    Entrepren<-profitsFinMaxim(InitProf,ParamsEntrep)
    Dec$InformalDemand<-Entrepren[1]
    Dec$FormalDemand<-Entrepren[2]
    Dec$Zoptimal<-Entrepren[3]
    ValueProfit<-Entrepren[4]
    Dec$ValueEntrepren<-ValueProfit
    
    
    
    #Entry (1,1) will say the decision
    
    Dec$Decission<-(ValueProfit>ValueWorker)*1
    return(Dec)
  }
  
  
  #Testing the function
  
  params=c(wi,wf,aalpha,ddelta,ggamma,bbeta,ssigma,kkappa,rrho,psi,chi)
  
  
  
  #Loading ttheta
  Ttheta<-c(tthetaw,tthetae)
  tthetaw<-Ttheta[1]
  tthetae<-Ttheta[2]
  InitProf=c(ni,nf,z)
  InitLWorkers=c(li,lf)
  
  iDecision(Ttheta,params,InitLWorkers,InitProf)
  
  #Parameters for excess demand
  Params=c(aalpha,ddelta,ggamma,bbeta,ssigma,kkappa,rrho,psi,chi,mmu1,mmu2,ssigma1,ssigma2,rho12)
  
  #Total demand of formal workers according to equation 2b. M denotes the number of monteCarlo 
  #draws
  ExcessDemandFunctions<-function(Wages,Params,M,InitLWorkers,InitProf){
    #First step will be to load the parameters
    wi<-Wages[1]
    wf<-Wages[2]
    
    aalpha<-Params[1]
    ddelta<-Params[2]
    ggamma<-Params[3]
    bbeta<-Params[4]
    ssigma<-Params[5]
    kkappa<-Params[6]
    rrho<-Params[7]
    psi<-Params[8]
    chi<-Params[9]
    mmu1<-Params[10]
    mmu2<-Params[11]
    ssigma1<-Params[12]
    ssigma2<-Params[13]
    rho12<-Params[14]
    
    #Loading correspondingly the inputs of iDecission
    params=c(wi,wf,aalpha,ddelta,ggamma,bbeta,ssigma,kkappa,rrho,psi,chi)
    
    #Although doing it in loop is slow, I will start in this to 
    #make sure in case need to readapt
    set.seed(257)
    #Creating empty vectors for each
    InformalDemand<-seq(1,M,1)
    InformalSupply<-seq(1,M,1)
    FormalDemand<-seq(1,M,1)
    FormalSupply<-seq(1,M,1)
    Entrepreneurs<-seq(1,M,1)
    
    
    #Define the variance covariance matrix
    
    #Old attempt
    #Sigma <- matrix(c(0.01,0.003,0.003,0.01),2,2)
    #mu=c(0,0)
    
    #New attempt
    Sigma <- matrix(c(ssigma1,rho12,rho12,ssigma2),2,2)
    mu=c(mmu1,mmu2)
    logtthetavec<-rmvnorm(M, mean = mu, Sigma)
    for(i in 1:M){
      #Obtaining draw of distribution of ttheta
      #First, define the variance
      
      tthetaw<-exp(logtthetavec[i,1])
      tthetae<-exp(logtthetavec[i,2])
      #tthetaw<-runif(1)+1
      #tthetae<-runif(1)+1
      #tthetaw<-1.2
      #tthetae<-5.5
      Ttheta<-c(tthetaw,tthetae)
      #Computing decisions
      Dec<-iDecision(Ttheta,params,InitLWorkers,InitProf)
      #1. Informal demmand
      InformalDemand[i]<-(Dec$InformalDemand)*(Dec$Decission)
      InformalSupply[i]<-(Dec$InformalSupply)*(1-Dec$Decission)*tthetaw
      FormalDemand[i]<-(Dec$FormalDemand)*(Dec$Decission)
      FormalSupply[i]<-(Dec$FormalSupply)*(1-Dec$Decission)*tthetaw
      Entrepreneurs[i]<-Dec$Decission
      
      
    }
    InformalExcessDemand<-sum(InformalDemand)-sum(InformalSupply)
    FormalExcessDemand<-sum(FormalDemand)-sum(FormalSupply)
    sum(InformalDemand)/M
    sum(InformalSupply)/M
    sum(FormalDemand)/M
    sum(FormalSupply)/M
    FormalExcessDemand
    InformalExcessDemand
    InformalDemand0=InformalDemand
    FormalDemand0=FormalDemand
    InformalSupply0=InformalSupply
    FormalSupply0=FormalSupply
    ans<-InformalExcessDemand^2+FormalExcessDemand^2
    ans<-ans/M
    return(ans)
    
  }
  
  M=100
  #ExcessDemandFunctions(c(wi,wf),Params,M,InitLWorkers,InitProf)
  
  
  
  Excess<-function(W){
    ans<-ExcessDemandFunctions(W,Params,M,InitLWorkers,InitProf)
    print(ans)
    return(ans)
  }
  
  
  #Trying the equilibrium function in C++
  ParamsDecisionExcessDemand=seq(1,1,19)
  ParamsDecisionExcessDemand[1]=aalpha
  ParamsDecisionExcessDemand[2]=ddelta
  ParamsDecisionExcessDemand[3]=ggamma
  ParamsDecisionExcessDemand[4]=bbeta
  ParamsDecisionExcessDemand[5]=ssigma
  ParamsDecisionExcessDemand[6]=kkappa
  ParamsDecisionExcessDemand[7]=rrho
  ParamsDecisionExcessDemand[8]=psi
  ParamsDecisionExcessDemand[9]=chi
  ParamsDecisionExcessDemand[10]=mmu1
  ParamsDecisionExcessDemand[11]=mmu2
  ParamsDecisionExcessDemand[12]=ssigma1
  ParamsDecisionExcessDemand[13]=ssigma2
  ParamsDecisionExcessDemand[14]=rho12
  ParamsDecisionExcessDemand[15]=li
  ParamsDecisionExcessDemand[16]=lf
  ParamsDecisionExcessDemand[17]=ni
  ParamsDecisionExcessDemand[18]=nf
  ParamsDecisionExcessDemand[19]=z
  
  
  #Wages initial guess
  WagesInitialGuess=c(wi,wf)
  WagesInitialGuess=c(wi,wf)
  
  #Trying the equilibrium function
  #WEq=EqWagesNumericVector(ParamsDecisionExcessDemand,WagesInitialGuess)
  wiEq=WEq[1]
  wfEq=WEq[2]
  wiEq=wi
  wfEq=wf
  
  
  #If you want to find the equilibrium in R:
  
  if(1==2){
    
    W=c(5.63,6.36)
    W=c(wiEq,wfEq)
    #W=c(wi,wf)
    #Excess(c(15.77,18.0332))
    Res<-optim(W,Excess,control=list(abstol=1.0e-2))
    wiEq<-Res$par[1]
    wfEq<-Res$par[2]
    WEQ=c(wiEq,wfEq)
    
    
    
  }
  
  
  #Identifying range of ttheta for the plot
  #Make sure this matrix coincides with the one defined in excess demand
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
  numbins=100
  stepttw=(max(tthetaw)-min(tthetaw))/numbins
  steptte=(max(tthetae)-min(tthetae))/numbins
  
  tthetaw_Sample<-seq(min(tthetaw),max(tthetaw),stepttw)
  tthetae_Sample<-seq(min(tthetae),max(tthetae),steptte)
  #Vector storing the decisions of individuals
  lw<-length(tthetaw_Sample)
  le<-length(tthetae_Sample)
  DecisionMatrix<-matrix(0,lw,le)
  x.m <- melt(t(DecisionMatrix))
  DecisionMatrixVer=matrix(0,le*lw,3)
  
  
  wi=wiEq
  wf=wfEq
  
  params<-c(wiEq,wfEq,aalpha,ddelta,ggamma,bbeta,ssigma,kkappa,rrho,psi,chi)
  
  
  tthetavec<-c(2,1)
  iDecision(tthetavec,params,InitLWorkers,InitProf)$Decission
  
  #Graph of entrepreneur-worker takes a while. Put it in conditional while testing parameters
  if(1==1){
  
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
      
      it=it+1
    }
  }
  
  
  
  
  #Identifying range of entrepreneurs and workers
  
  
  
  
  EntSubset<-subset(DecisionMatrixVer,DecisionMatrixVer[,3]==1)
  minEnt<-min(EntSubset[,1])
  maxEnt<-max(EntSubset[,1])
  
  WorkSubset<-subset(DecisionMatrixVer,DecisionMatrixVer[,3]==0)
  minWork<-min(WorkSubset[,2])
  maxWork<-max(WorkSubset[,2])
  
  #Convert it to a data.table
  x.m <- melt(t(DecisionMatrix))
  colnames(x.m)<-c("ttheta_e","ttheta_w","Decision")
  DecisionMatrixVer2=as.data.frame(DecisionMatrixVer)
  colnames(DecisionMatrixVer2)<-c("ttheta_e","ttheta_w","Decision")
  
  p<-ggplot(DecisionMatrixVer2, aes(ttheta_w, ttheta_e, fill = factor(Decision)))+geom_tile()
  p<-p+labs(x=expression(theta~w),y=expression(theta~e))
  p<-p+ scale_fill_discrete(name="Decision",labels=c("Workers","Entrepreneurs"))
  #p<-p+scale_x_continuous(limits = c(min(tthetaw), max(tthetaw)))+scale_y_continuous(limits=c(min(tthetae),max(tthetae)))
  p<-p+theme_bw()
  p<-p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  
  setwd(GRAPHS)
  dev.set()
  png(file="EntreprenDecision.png",width=1600,height=850)
  p
  dev.off()
  }
  #-------------------------------------#
  #Now doing analysis in equilibrium    # 
  #-------------------------------------#
  
  
  #-------------------------#
  #Analysis of Entrepreneurs#
  #-------------------------#
  
  
  
  #We will do a different sample of tthetas for this purpose. 
  #In the analysis of entrepreneur-worker it is not important to have
  #the actual distribution. However, we do want the distribution. I will sort
  #the tthetas based on their distribution. 
  
  logtthetavec<-rmvnorm(100000, mean = mu, Sigma)
  tthetae_Sample<-sort(exp(logtthetavec[,2]))
  maxEnt<-max(tthetae_Sample)
  maxWork<-max(tthetaw_Sample)
  
  
  logtthetavec<-rmvnorm(1000, mean = mu, Sigma)
  tthetae_Sample<-sort(exp(logtthetavec[,2]))
  
  le<-length(tthetae_Sample)
  
  
  #Identifying the range of workers and entrepreneurs
  Dec=0
  tte=0
  while(Dec==0){
    tte<-tte+1
    tthetavec<-c(minTtw,tthetae_Sample[tte])
    tthetavec<-c(0.0001,tthetae_Sample[tte])
    Dec<-iDecision(tthetavec,params,InitLWorkers,InitProf)$Decission
  }
  minEnt<-tthetae_Sample[tte]
  
  
  #Identifying the range of workers and entrepreneurs
  Dec=1
  ttw=0
  while(Dec==1){
    ttw<-ttw+1
    tthetavec<-c(tthetaw_Sample[ttw],minTte)
    Dec<-iDecision(tthetavec,params,InitLWorkers,InitProf)$Decission
  }
  minWork<-tthetaw_Sample[ttw]
  
  
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
                                 Production,AfterTaxProfit,Taxpayed,TaxSales,TaxProfits,tthetae_Trunc,Zproportion2,
                                 TaxpayedProp))
  
  
  
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
  
  #If want to do moments based on percentiles, not the whole data, run the following line:
  #zevasion1<-A1
  
  
  
  
  #Informal Demand
  p<-ggplot(data=zevasion1,aes(x=tthetae_Sample,y=InformalDemand))+geom_line()
  p<-p+labs(x=expression(theta~e),y="Demand Informal Workers")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="InformalDemand.png",width=1600,height=850)
  p
  dev.off()
  
  
  
  #Informal Demand percentile
  p<-ggplot(data=zevasion1,aes(x=tthetae_Trunc,y=InformalDemand))+geom_line()
  p<-p+labs(x=expression(theta~e),y="Demand Informal Workers")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  
  dev.set()
  png(file="InformalDemandPercentile.png",width=1600,height=850)
  p
  dev.off()
  
  #Formal Demand
  p<-ggplot(data=zevasion1,aes(x=tthetae_Sample,y=FormalDemand))+geom_line()
  p<-p+labs(x=expression(theta~e),y="Demand Formal Workers")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="FormalDemand.png",width=1600,height=850)
  p
  dev.off()
  
  
  #Formal Demand Percentiles
  p<-ggplot(data=zevasion1,aes(x=tthetae_Trunc,y=FormalDemand))+geom_line()
  p<-p+labs(x=expression(theta~e),y="Demand Formal Workers")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="FormalDemandPercentile.png",width=1600,height=850)
  p
  dev.off()
  
  #Total Labor Force
  p<-ggplot(data=zevasion1,aes(x=tthetae_Sample,y=TotalLaborForce))+geom_line()
  p<-p+labs(x=expression(theta~e),y="Total Labor Demand")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="TotalLaborDemand.png",width=1600,height=850)
  p
  dev.off()
  
  
  #Total Labor Force Percentiles
  p<-ggplot(data=zevasion1,aes(x=tthetae_Trunc,y=TotalLaborForce))+geom_line()
  p<-p+labs(x=expression(theta~e),y="Total Labor Demand")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="TotalLaborDemandPercentile.png",width=1600,height=850)
  p
  dev.off()
  
  
  #InformalProportion
  p<-ggplot(data=zevasion1,aes(x=tthetae_Sample,y=InformalProportion))+geom_line()
  p<-p+labs(x=expression(theta~e),y="Proportion Informal Demand")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="InformalProportionDemand.png",width=1600,height=850)
  p
  dev.off()
  
  #InformalProportion percentiles
  p<-ggplot(data=zevasion1,aes(x=tthetae_Trunc,y=InformalProportion))+geom_line()
  p<-p+labs(x=expression(theta~e),y="Proportion Informal Demand")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="InformalProportionDemandPercentile.png",width=1600,height=850)
  p
  dev.off()
  
  #Pre tax profit
  p<-ggplot(data=zevasion1,aes(x=tthetae_Sample,y=PretaxProfit))+geom_line()
  p<-p+labs(x=expression(theta~e),y="Pre tax profit")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="PretaxProfit.png",width=1600,height=850)
  p
  dev.off()
  
  
  #Pre tax profit percentile
  p<-ggplot(data=zevasion1,aes(x=tthetae_Trunc,y=PretaxProfit))+geom_line()
  p<-p+labs(x=expression(theta~e),y="Pre tax profit percentile")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="PretaxProfitPercentile.png",width=1600,height=850)
  p
  dev.off()
  
  
  #After tax profit
  p<-ggplot(data=zevasion1,aes(x=tthetae_Sample,y=AfterTaxProfit))+geom_line()
  p<-p+labs(x=expression(theta~e),y="After Tax Profit")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="AfterTaxProfit.png",width=1600,height=850)
  p
  dev.off()
  
  
  #After tax profit
  p<-ggplot(data=zevasion1,aes(x=tthetae_Trunc,y=AfterTaxProfit))+geom_line()
  p<-p+labs(x=expression(theta~e),y="After Tax Profit percentile")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="AfterTaxProfitPercentile.png",width=1600,height=850)
  p
  dev.off()
  
  #Production
  p<-ggplot(data=zevasion1,aes(x=tthetae_Sample,y=Production))+geom_line()
  p<-p+labs(x=expression(theta~e),y="Sales")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="Production.png",width=1600,height=850)
  p
  dev.off()
  
  
  #Production percentiles
  p<-ggplot(data=zevasion1,aes(x=tthetae_Trunc,y=Production))+geom_line()
  p<-p+labs(x=expression(theta~e),y="Sales percentile")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="ProductionPercentile.png",width=1600,height=850)
  p
  dev.off()
  
  
  
  #Evasion proportion 
  p<-ggplot(data=zevasion1,aes(x=tthetae_Sample,y=Zoptimal))+geom_line()
  p<-p+labs(x=expression(theta~e),y="Evasion")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="EvasionProportion.png",width=1600,height=850)
  p
  dev.off()
  
  #Evasion proportion percentiles 
  p<-ggplot(data=zevasion1,aes(x=tthetae_Trunc,y=Zproportion))+geom_line()
  p<-p+labs(x=expression(theta~e),y="Evasion proportion percentiles")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="EvasionProportionPercentiles.png",width=1600,height=850)
  p
  dev.off()
  
  
  
  #Evasion values percentiles 
  p<-ggplot(data=zevasion1,aes(x=tthetae_Trunc,y=Zoptimal))+geom_line()
  p<-p+labs(x=expression(theta~e),y="Evasion values percentiles")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="EvasionValuesPercentiles.png",width=1600,height=850)
  p
  dev.off()
  
  
  #Taxes payed as proportino of sales
  p<-ggplot(data=zevasion1,aes(x=tthetae_Trunc,y=TaxSales))+geom_line()
  p<-p+labs(x=expression(theta~e),y="Evasion proportion percentiles")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="TaxesPayed.png",width=1600,height=850)
  p
  dev.off()
  
  
  #Taxes payed
  p<-ggplot(data=zevasion1,aes(x=tthetae_Trunc,y=TaxSales))+geom_line()
  p<-p+labs(x=expression(theta~e),y="Taxes payed percentiles")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="TaxesPayed.png",width=1600,height=850)
  p
  dev.off()
  
  
  #Taxes payed proportion
  p<-ggplot(data=zevasion1,aes(x=tthetae_Trunc,y=TaxpayedProp))+geom_line()
  p<-p+labs(x=expression(theta~e),y="Taxes payed proportion percentiles")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="TaxesPayedProportion.png",width=1600,height=850)
  p
  dev.off()
  
  
  
  
  
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
  #Worker<-A1
  
  
  
  #Formal supply
  p<-ggplot(data=Worker,aes(x=tthetaw_Sample,y=FormalSupply))+geom_line()
  p<-p+labs(x=expression(theta~w),y="Formal Labor Supply")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="FormalSupply.png",width=1600,height=850)
  p
  dev.off()
  
  #Formal supply percentiles
  p<-ggplot(data=Worker,aes(x=tthetaw_Trunc,y=FormalSupply))+geom_line()
  p<-p+labs(x=expression(theta~w),y="Formal Labor Supply")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  
  dev.set()
  png(file="FormalSupplyPercentiles.png",width=1600,height=850)
  p
  dev.off()
  
  #Informal Supply
  p<-ggplot(data=Worker,aes(x=tthetaw_Sample,y=InformalSupply))+geom_line()
  p<-p+labs(x=expression(theta~w),y="Informal Labor Supply")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="InformalSupply.png",width=1600,height=850)
  p
  dev.off()
  
  
  #Informal Supply Percentiles
  p<-ggplot(data=Worker,aes(x=tthetaw_Trunc,y=InformalSupply))+geom_line()
  p<-p+labs(x=expression(theta~w),y="Informal Labor Supply")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="InformalSupplyPercentiles.png",width=1600,height=850)
  p
  dev.off()
  
  #Value Worker
  p<-ggplot(data=Worker,aes(x=tthetaw_Sample,y=ValueWorker))+geom_line()
  p<-p+labs(x=expression(theta~w),y="Value Worker")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="ValueWorker.png",width=1600,height=850)
  p
  dev.off()
  
  
  #Value Worker Percentiles
  p<-ggplot(data=Worker,aes(x=tthetaw_Trunc,y=ValueWorker))+geom_line()
  p<-p+labs(x=expression(theta~w),y="Value Worker")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="ValueWorkerPercentiles.png",width=1600,height=850)
  p
  dev.off()
  
  
  
  #Total Labor force
  p<-ggplot(data=Worker,aes(x=tthetaw_Sample,y=TotalLaborSupply))+geom_line()
  p<-p+labs(x=expression(theta~w),y="Total Labor Supply")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="TotalLaborSupply.png",width=1600,height=850)
  p
  dev.off()
  
  
  #Total Labor force Percentiles
  p<-ggplot(data=Worker,aes(x=tthetaw_Trunc,y=TotalLaborSupply))+geom_line()
  p<-p+labs(x=expression(theta~w),y="Total Labor Supply")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="TotalLaborSupplyPercentiles.png",width=1600,height=850)
  p
  dev.off()
  
  #Informal Proportion
  p<-ggplot(data=Worker,aes(x=tthetaw_Sample,y=InformalProportion))+geom_line()
  p<-p+labs(x=expression(theta~w),y="Informal Labor Supply proportion")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="InformalProportionSupply.png",width=1600,height=850)
  p
  dev.off()
  
  
  #Informal Labor Supply Proportion Percentiles
  p<-ggplot(data=Worker,aes(x=tthetaw_Trunc,y=InformalProportion))+geom_line()
  p<-p+labs(x=expression(theta~w),y="Informal Labor Supply proportion")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="InformalProportionSupplyPercentiles.png",width=1600,height=850)
  p
  dev.off()
  
  
  #Formal income
  p<-ggplot(data=Worker,aes(x=tthetaw_Trunc,y=FormalIncome))+geom_line()
  p<-p+labs(x=expression(theta~w),y="Formal Income")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="FormalIncomePercentiles.png",width=1600,height=850)
  p
  dev.off()
  
  
  #InFormal income
  p<-ggplot(data=Worker,aes(x=tthetaw_Trunc,y=InformalIncome))+geom_line()
  p<-p+labs(x=expression(theta~w),y="Informal Income")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="InformalIncomePercentiles.png",width=1600,height=850)
  p
  dev.off()
  
  
  
  
  #Total income
  p<-ggplot(data=Worker,aes(x=tthetaw_Trunc,y=TotalIncome))+geom_line()
  p<-p+labs(x=expression(theta~w),y="Total Income")
  p<-p+theme(axis.text=element_text(size=48),
             axis.title=element_text(size=48,face="bold"),
             legend.title=element_text(size=48),
             legend.text=element_text(size=48))
  p
  dev.set()
  png(file="TotalIncomePercentiles.png",width=1600,height=850)
  p
  dev.off()


