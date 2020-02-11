//
//  main.cpp
//  Equilibrium
//
//  Created by Rodrigo Azuero Melo on 5/11/18.
//  Copyright © 2018 Rodrigo Azuero Melo. All rights reserved.
//

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]

#include <iostream>
#include <boost/math/distributions/normal.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/random.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/math/distributions/uniform.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <ctime>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <sstream>
#include <string>
#include <stdio.h>
#include <vector>
#include <unistd.h>
#include <stdio.h>
#include <errno.h>
#include <random>
#include <unistd.h>
#include <nlopt.hpp>
#include <utility>
//#include <Rcpp.h> -> Not necessary if rcpparmadillo included
#include <RcppArmadillo.h>
//#include <bits/stdc++.h>

using std::vector;
using namespace std;
using namespace Rcpp;

ofstream ThMoments0CSV("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolsGenerated/ThMoments0CSV.csv");
ofstream ThMoments1CSV("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolsGenerated/ThMoments1CSV.csv");
ofstream ThMoments2CSV("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolsGenerated/ThMoments2CSV.csv");
ofstream ThMoments3CSV("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolsGenerated/ThMoments3CSV.csv");
ofstream ThMoments4CSV("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolsGenerated/ThMoments4CSV.csv");
ofstream ThMoments5CSV("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolsGenerated/ThMoments5CSV.csv");
ofstream ThMoments6CSV("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolsGenerated/ThMoments6CSV.csv");
ofstream ThMoments7CSV("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolsGenerated/ThMoments7CSV.csv");

ofstream ThMoments8CSV("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolsGenerated/ThMoments8CSV.csv");

ofstream ParametersCSV("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolsGenerated/ParametersCSV.csv");

ofstream DistanceMoments("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolsGenerated/DistanceMoments.csv");

ofstream EquilibriumValue("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolsGenerated/EquilibriumValue.csv");


typedef std::pair<double,double> mypair;
bool comparator ( const mypair& l, const mypair& r)
{ return l.first < r.first; }

//Function for cholesky decomposition
double *cholesky(double *A, int n) {
    double *L = (double*)calloc(n * n, sizeof(double));
    if (L == NULL)
        exit(EXIT_FAILURE);
    
    for (int i = 0; i < n; i++)
        for (int j = 0; j < (i+1); j++) {
            double s = 0;
            for (int k = 0; k < j; k++)
                s += L[i * n + k] * L[j * n + k];
            L[i * n + j] = (i == j) ?
            sqrt(A[i * n + i] - s) :
            (1.0 / L[j * n + j] * (A[i * n + j] - s));
        }
    
    return L;
}

//Function to show matrix decomposition cholesky
void show_matrix(double *A, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            printf("%2.5f ", A[i * n + j]);
        printf("\n");
    }
}


//Templates def
//Define template for random number generator
template<class T>
double gen_normal_3(T &generator)
{
    return generator();
}


//0. Payroll taxes definition (marginal rate)

long double Tn(double nf){
    return(0.09);
}

//1. Payroll taxes definition (actual payroll)
//long double TnActual(double nf){
//    return(0.09*nf);
//}


//1. Payroll taxes definition (actual payroll)
long double TnActual(double nomina){
    
    
    //double term1=pow(nf,2);
    //double term2=1+term1;
    //double term3=-1/term2;
    //double term4=term3+1;
    //double term5=0.45*term4;
    //double rate=term5+0.05;
    
    
    return(0.09*nomina);
}

//1.9. Production. C*0 to eliminate c
long double production(double ni, double nf, double aalpha, double tthetae,double c){
    double prod=tthetae*pow((c*1+ni+nf),aalpha);
    return(prod);
}

//2. Pre-tax profits.Checked. C*0 to eliminate c
long double profm(double ni, double nf, double aalpha, double tthetae, double wi, double wf, double c){
    double pi1=tthetae*pow((c*1+ni+nf),aalpha)-wi*ni-wf*nf-TnActual(nf*wf);
    return(pi1);
    
}

//3. Corporate tax profits, marginal. C*0 to eliminate c
long double Tc(double z, double ni, double nf, double aalpha,double tthetae, double wi, double wf,double c){
    double profm=tthetae*pow((c*1+ni+nf),aalpha)-wi*ni-wf*nf-TnActual(nf*wf);
    double arg=profm-z;
    double firsterm=-exp(-arg);
    //If negative profits, set zero marginal rate
    if (profm<0){
        firsterm=0;
    }
    return(firsterm);
}

//4. Tc Actual corporate taxes
//4. Tc Actual corporate taxes. In hundreds of dollars.
// For instance. 189=5000*12*0.315/100 c*0 to eliminate c
long double TcActual(double z, double ni, double nf, double aalpha, double tthetae, double wi, double wf,double c){
    
    double profm=tthetae*pow((c*1+ni+nf),aalpha)-wi*ni-wf*nf-TnActual(nf*wf);
    double arg=profm-z;
    double ans=0;
    double tax=0;
    double prod=tthetae*pow((c*1+ni+nf),aalpha);
    if(prod-z<=0){
        ans=0;
    }
    else if (prod-z<189){
        ans=0.756;
    } else if (prod-z<302.4){
        ans=1.89;
    } else if (prod-z<491.4){
        ans=7.56;
    } else if (prod-z<756){
        ans=15.12;
    } else if (prod-z<1134){
        ans=22.68;
    }else if (prod-z>=1134){
        arg=profm-z;
        tax=0.3;
        ans=arg*tax;
        if(ans<0){
            ans=0;
        }
    }
    return(ans);
    //double firsterm=-exp(-arg);
    //double tax=firsterm+1;
    //double tax=0.5*(exp(arg)-1);
    
    //if(tax<0){
    //    tax=0;
    //}
    
    
}

//5. Personal income tax
long double PIT(double tthetaw, double wf, double lf){
    double x=tthetaw*wf*lf;
    double a=1000;
    double ans=0;
    if(x<a){
        ans=(100/a)*x-100;
    }else if(x<36000){
        ans=0.0;
    }else if(x<66000){
        ans=pow(x,2)/100000-9*x/25;
    }else if(x>=66000){
        ans=0.3*x;
    }
    return(ans);
}

//6. Personal income tax marginal
long double PITM(double tthetaw, double wf, double lf){
    return(0.5);
}

//7. Final profits
long double FinProfits(const double *Args, double paramvec[9]){
    //vector<double> paramvec=*(vector<double>*)params;
    //Before, it was (const double *Args, void *params){
    
    double wi=paramvec[0];
    double wf=paramvec[1];
    double aalpha=paramvec[2];
    double ddelta=paramvec[3];
    double ggamma=paramvec[4];
    double bbeta=paramvec[5];
    double ssigma=paramvec[6];
    double tthetae=paramvec[7];
    double c=paramvec[8];//added
    
    
   // cout << ggamma << " ggamma "<< endl;
   // cout << ddelta << " ddelta "<< endl;
    
    int test=2;
    if(test==1){
        cout << aalpha << " aalpha"<< endl;
        cout << ddelta << " ddelta"<< endl;
        cout << ggamma << " ggamma"<< endl;
        cout << bbeta << " bbeta"<< endl;
        cout << ssigma << " ssigma"<< endl;
        cout << tthetae << " tthetae"<< endl;
        cout << c << " c"<< endl;
    }
 
    
    
    
    double ni=Args[0];
    double nf=Args[1];
    double z=Args[2];
    
    //Operational
    double term1=profm(ni, nf, aalpha, tthetae, wi, wf,c);
    
    //Corporate taxes
    double taxes=TcActual(z,ni,nf,aalpha,tthetae,wi,wf,c);
    
    //Cost of evacion
    double evcost=pow(z,1+ssigma)*(bbeta/(1+ssigma));
    
    //Evasion costs
    double infcost=pow(ni,1+ggamma)*(ddelta/(1+ggamma));
    //cout << ggamma << " ggamma "<< endl;
    //cout << ggamma << " ggamma "<< endl;
    double ans=term1-taxes-infcost-evcost;
    
    return(ans);
    
    
}

//8. Maximize profits
//Defining the profit function to be maximized by object nlopt
int iterat=0;
double profitsFinMaxim(unsigned n, const double *x, double *grad, void *profitsFinMaxim_data){
    //Void to double
    //cout << profitsFinMaxim_data<< " this stuff"<< endl;
    double *Params = (double *)profitsFinMaxim_data;
    
    
    
    
    
    
    double result=FinProfits(x, Params);
    //cout << result << " evalF"<< endl;
    iterat=iterat+1;
    //printf("Iteration=(%d); Feval=%0.10g\n", iterat, result);
    return(result);
}

//9. Defininf a function such that for given parameters, returns the optimal inputs and the maximum profits
vector<double> profitsFinMaxim(double InitialCond[3], double wi,
                               double wf,
                               double aalpha,
                               double ddelta,
                               double ggamma,
                               double bbeta,
                               double ssigma,
                               double tthetae,
                               double c){
    
    double Params[9]={};
    Params[0]=wi;
    Params[1]=wf;
    Params[2]=aalpha;
    Params[3]=ddelta;
    Params[4]=ggamma;
    Params[5]=bbeta;
    Params[6]=ssigma;
    Params[7]=tthetae;
    Params[8]=c;
    
    //Defining object to be optimized
    //First. lower bounds for ni, nf, z.
    double lb[3];
    lb[0]=0;
    lb[1]=0;
    lb[2]=0;
    
    //Setting up optimization
    nlopt_opt optProfits;
    optProfits = nlopt_create(NLOPT_LN_NELDERMEAD, 3);
    nlopt_set_lower_bounds(optProfits, lb);
    nlopt_set_max_objective(optProfits, profitsFinMaxim,(void *)&Params);
    
    //Tolerance, abs or rel.
    nlopt_set_xtol_rel(optProfits, 1.0e-8);
    nlopt_set_ftol_abs(optProfits,1.0e-8);
    double minf3; /* the minimum objective value, upon return */
    
    //Setting up initial conditions
    double x3[3];
    x3[0]=InitialCond[0];
    x3[1]=InitialCond[1];
    x3[2]=InitialCond[2];
    
    //Finding optimal
    nlopt_optimize(optProfits, x3, &minf3);
    
    //if (nlopt_optimize(optProfits, x3, &minf3) < 0) {
    //    printf("nlopt failed!\n");
    //}
    //else {
    //   printf("found minimum at f(%g) = %0.10g\n", x3[0], minf3);
    //}
    
    
    double ni=x3[0];
    double nf=x3[1];
    double z=x3[2];
    double Args[3];
    
    
    Args[0]=ni;
    Args[1]=nf;
    Args[2]=z;
    
    
    //Args[0]=ni;
    //Args[1]=nf;
    //Args[2]=z;
    
    //Findin the profits
    double Prof=FinProfits(Args, Params);
    
    vector<double> Ans;
    Ans.resize(4);
    Ans[0]=ni;
    Ans[1]=nf;
    Ans[2]=z;
    Ans[3]=Prof;
    nlopt_destroy(optProfits);
    
    
    return(Ans);
}





//10. Value of workers
double ValueWorkers(const double *Args, double ParamWorkers[7]){
    double wf=ParamWorkers[0];
    double wi=ParamWorkers[1];
    double kkappa=ParamWorkers[2];
    double rrho=ParamWorkers[3];
    double psi=ParamWorkers[4];
    double chi=ParamWorkers[5];
    double tthetaw=ParamWorkers[6];
    
    double li=Args[0];
    const double lf=Args[1];
    

    
    //Term1 is total income
    double term1=(wf*lf+wi*li)*tthetaw;
    
    //Term2 total disutility from working
    double term2=pow(li+lf,1+psi)*(chi/(1+psi));
    
    //Term 3 is consumption penalty from informality
    double term3=kkappa*pow(tthetaw*li,1+rrho)/(1+rrho);
    
    //Term4 is net transfers to.from gogernment
    double Taxes=PIT(tthetaw,wf,lf);
    
    
    double ans=term1-term2-term3-Taxes;
    return(ans);
}



//11. Maximize value of workers
int iteratvalworkers=0;
double valueWorkerMax(unsigned n, const double *x, double *grad, void *valueWorkerMax_data){
    //Void to double
    //cout << profitsFinMaxim_data<< " this stuff"<< endl;
    double *Params = (double *)valueWorkerMax_data;
    
    double result=ValueWorkers(x, Params);
    //cout << result << " evalF"<< endl;
    iteratvalworkers=iteratvalworkers+1;
    //printf("Iteration=(%d); Feval=%0.10g\n", iteratvalworkers, result);
    return(result);
}



vector<double> valueWorkerFinMaxim(double InitialCond[2],
                                   double wf,
                                   double wi,
                                   double kkappa,
                                   double rrho,
                                   double psi,
                                   double chi,
                                   double tthetaw){
    double Params[7]={};
    Params[0]=wf;
    Params[1]=wi;
    Params[2]=kkappa;
    Params[3]=rrho;
    Params[4]=psi;
    Params[5]=chi;
    Params[6]=tthetaw;
    
    
    
    
    //Defining the object to be maximized
    double lbWorker[2];
    lbWorker[0]=0;
    lbWorker[1]=0;
    
    
    nlopt_opt optWorkers;
    optWorkers = nlopt_create(NLOPT_LN_NELDERMEAD, 2);
    
    nlopt_set_lower_bounds(optWorkers, lbWorker);
    nlopt_set_max_objective(optWorkers, valueWorkerMax,(void *)&Params);
    
    
    nlopt_set_xtol_rel(optWorkers, 1.0e-8);
    nlopt_set_ftol_abs(optWorkers,1.0e-8);
    double MaxfWorker; /* the minimum objective value, upon return */
    double x3Worker[2];
    x3Worker[0]=InitialCond[0];
    x3Worker[1]=InitialCond[1];
    
    
    //memcpy(Params, (x3**)pointer, sizeof Params);
    
    //cout <<FinProfits(x3,Params)<< " test inicial"<< endl;
    //cout << x3[2]<< " x[2]"<< endl;
    nlopt_optimize(optWorkers, x3Worker, &MaxfWorker);
    nlopt_destroy(optWorkers);
    //if (nlopt_optimize(optWorkers, x3Worker, &MaxfWorker) < 0) {
    //    printf("nlopt failed!\n");
    //}
    //else {
    //    printf("found minimum at f(%g) = %0.10g\n",x3Worker[0], MaxfWorker);
    //}
    
    //Actual value of workers
    
    vector<double> ans;
    ans.resize(3);
    ans[0]=x3Worker[0];
    ans[1]=x3Worker[1];
    
    //cout << ans[0] << " ans[0]"<< endl;
    //wcout << ans[1] << " ans[1]"<< endl;
    ans[2]=MaxfWorker;
    return(ans);
    
    
}

//Vector de decisión
vector<vector<double> > iDecision(vector<double> Ttheta,
                                  vector<double> Params,
                                  vector<double> InitLWorkers,
                                  vector<double> InitProf){
    
    //Loading paradmeters
    double wi=Params[0];
    double wf=Params[1];
    double aalpha=Params[2];
    double ggamma=Params[3];
    double ddelta=Params[4];
    double bbeta=Params[5];
    double ssigma=Params[6];
    double kkappa=Params[7];
    double rrho=Params[8];
    double psi=Params[9];
    double chi=Params[10];
    double c=Params[11];
    
    
    
    
    
    //Loading ttheta
    double tthetaw=Ttheta[0];
    double tthetae=Ttheta[1];
    
    //Computing the value of the workers
    
    vector<double> ansValWorker;
    ansValWorker.resize(3);
    double InitialVworkers[2];
    InitialVworkers[0]=InitLWorkers[0];
    InitialVworkers[1]=InitLWorkers[1];
    
    
    
    
    if(1==2){
        cout << " inside of the idecision stuff "<< endl;
        cout << wi<< " wi "<< endl;
        cout <<  wf<< " wf "<< endl;
        cout <<  aalpha<< " aalpha "<< endl;
        cout <<  ddelta<< " ddelta "<< endl;
        cout <<  ggamma<< " ggamma "<< endl;
        cout <<  bbeta<< " bbeta "<< endl;
        cout <<  ssigma<< " ssigma "<< endl;
        cout <<  kkappa<< " kkappa "<< endl;
        cout <<  rrho<< " rrho "<< endl;
        cout <<  psi<< " psi "<< endl;
        cout <<  chi<< " chi "<< endl;
        cout << c << " c " << endl;
        cout << tthetaw << " tthetaw "<< endl;
        cout << tthetae << " tthetae "<< endl;
        cout << InitialVworkers[0] << " InitialVworkers[0]"<< endl;
        cout <<InitialVworkers[1] << " InitialVworkers[1]"<< endl;
    }
    
    
    
    
    ansValWorker=valueWorkerFinMaxim(InitialVworkers,wf,wi, kkappa,rrho,psi, chi, tthetaw);
    
    //Loading answers of value of workers
    double li=ansValWorker[0];
    
    double lf=ansValWorker[1];
    double ValueofWorker=ansValWorker[2];
    
    
    //Computing the Profits answer
    vector<double> ansProfits;
    ansProfits.resize(4);
    double InitialVProfits[3];
    InitialVProfits[0]=InitProf[0];
    InitialVProfits[1]=InitProf[1];
    InitialVProfits[2]=InitProf[2];
    ansProfits=profitsFinMaxim(InitialVProfits,  wi,wf, aalpha,ddelta,ggamma,bbeta,ssigma,tthetae,c);
    //Loading the answer of vlaue of profits
    double ni=ansProfits[0];
    double nf=ansProfits[1];
    double z=ansProfits[2];
    double prof=ansProfits[3];
    
    vector<vector<double> > ans;
    ans.resize(4);
    for (int it=0; it<4; it++){
        ans[it].resize(4);
    }
    //Loading to the matrix ans the solutions
    ans[0][0]=0;
    if(prof>ValueofWorker){
        ans[0][0]=1;
    }
    
    ans[0][1]=ni;
    ans[1][1]=nf;
    ans[2][1]=z;
    ans[3][1]=prof;
    
    ans[0][2]=li;
    ans[1][2]=lf;
    ans[2][2]=ValueofWorker;
    
    
    return(ans);
}






double ExcessDemandFunctions(vector<double>Wages,vector<double>Params,vector<double> InitLWorkers,vector<double> InitProf){
    
    
    //Checking inputs
    
    
    
    
    
    
    int M=1000;
    double InformalDemand[M];
    double InformalSupply[M];
    double FormalDemand[M];
    double FormalSupply[M];
    double Entrepreneurs[M];
    
    for(int it=0; it<M; it++){
        InformalDemand[it]=0;
        InformalSupply[it]=0;
        FormalDemand[it]=0;
        FormalSupply[it]=0;
    }
    
    
    
    
    //double x=unidistrib(rng);
    //cout << x << " unidit"<< endl;
    
    //Initialize the variables to iterate over
    vector<double> Ttheta;
    Ttheta.resize(2);
    
    int Decision=0;
    
    //Params for decission has the following structure:
    vector<double> ParamsDecision;
    ParamsDecision.resize(12);
    ParamsDecision[0]=Wages[0];
    ParamsDecision[1]=Wages[1];
    ParamsDecision[2]=Params[0];
    ParamsDecision[3]=Params[1];
    ParamsDecision[4]=Params[2];
    ParamsDecision[5]=Params[3];
    ParamsDecision[6]=Params[4];
    ParamsDecision[7]=Params[5];
    ParamsDecision[8]=Params[6];
    ParamsDecision[9]=Params[7];
    ParamsDecision[10]=Params[8];
    ParamsDecision[11]=Params[14];
    
    
    vector<vector<double> > DecVector;
    DecVector.resize(4);
    for (int it=0; it<4; it++){
        DecVector[it].resize(4);
    }
    
    //Initializing random number generator for the distribution
    int SEED=2581633;
    typedef boost::mt19937 RNGType;
    RNGType rng(SEED);
    boost::random::uniform_real_distribution<> unidistrib(1,2);
    boost::variate_generator<RNGType, boost::normal_distribution<> >
    generator(rng,
              boost::normal_distribution<>());
    
    double S0rr=gen_normal_3(generator);
    S0rr=gen_normal_3(generator);
    
    //We need to obtain the cholesky decomposition of the variance-covariance matrix
    int n=2;
    
    
    //double m2[] = {0.01, 0.003,  0.003,  0.01};
    double mmu1=Params[9];
    double mmu2=Params[10];
    double ssigma1=Params[11];
    double ssigma2=Params[12];
    double rho12=Params[13];
    
    
    
    
    double m2[] = {ssigma1, rho12,  rho12,  ssigma2};
    double *c2 = cholesky(m2, n);
    
    //Obtaining draws from multivariate pareto distribution
    
    //Draws from the normal distribution
    double z1=0;
    double z2=0;
    //Mean of log-normal thing
    
    
    for(int it=0; it<M; it++){
        
        //Obtaining the draws from a standard normal distribution
        z1=gen_normal_3(generator);
        z2=gen_normal_3(generator);
        
        //Doing the corresponding transformation to obtain them from bivariate distribution
        Ttheta[0]=exp(c2[0]*z1+c2[1]*z2+mmu1);
        Ttheta[1]=exp(c2[2]*z1+c2[3]*z2+mmu2);
        
        //cout << Ttheta[0] << " ttheta0"<< endl;
        //cout << Ttheta[1]<< " ttheta1"<< endl;
        
        //Ttheta[0]=1.2;
        //Ttheta[1]=5.5;
        //Ttheta[0]=unidistrib(rng);
        //Ttheta[1]=unidistrib(rng);
        
        //cout << Ttheta[0] << " Ttheta[0]"<< endl;
        //cout << Ttheta[1] << " Ttheta[1]"<< endl;
        //Computing the decissions
        DecVector=iDecision(Ttheta, ParamsDecision, InitLWorkers, InitProf);
        Entrepreneurs[it]=DecVector[0][0];
        Decision=DecVector[0][0];
        
        
        
        
        InformalDemand[it]=DecVector[0][1]*Decision;
        InformalSupply[it]=DecVector[0][2]*(1-Decision)*Ttheta[0];
        FormalDemand[it]=DecVector[1][1]*Decision;
        FormalSupply[it]=DecVector[1][2]*(1-Decision)*Ttheta[0];
        if(1==2){
            cout << "---- "<< endl;
            cout << it << " it "<< endl;
            cout << Ttheta[0] << " Ttheta0"<< endl;
            cout << Ttheta[1] << " Ttheta[1]" << endl;
            cout << Wages[0]<< " wages[0]"<< endl;
            cout << Wages[1]<< " wages[1]"<< endl;
            cout << Decision << " Decision "<< endl;
            cout << InformalDemand[it] << " InformalDemand[it] "<< endl;
            cout << InformalSupply[it] << " InformalSupply[it] "<< endl;
            cout << FormalDemand[it] << " FormalDemand[it] "<< endl;
            cout << FormalSupply[it] << " FormalSupply[it] "<< endl;
        }
        if(1==2){
            
            cout << " inside of the loop the decision -----"<<endl;
            cout << InitLWorkers[0] << " initlworkers0"<< endl;
            cout << InitLWorkers[1] << " initlworkers1"<< endl;
            cout << InitProf[0] << " InitProf0"<< endl;
            cout << InitProf[1] << " InitProf1"<< endl;
            cout << InitProf[2] << " InitProf2"<< endl;
            cout << DecVector[0][0]<< " Decision "<< endl;
            cout << DecVector[0][1]<< " ni "<< endl;
            cout << DecVector[1][1]<< " nf "<< endl;
            cout << DecVector[2][1]<< " z "<< endl;
            cout << DecVector[3][1]<< " prof "<< endl;
            cout << DecVector[0][2]<< " li "<< endl;
            cout << DecVector[1][2]<< " lf "<< endl;
            cout << DecVector[2][2]<< " VWorker "<< endl;
            
            
        }
        
        
        
    }
    
    //Obtaining the sum of each
    double InformalExcessDemand=0;
    double FormalExcessDemand=0;
    
    
    
    for(int it=0; it<M; it++){
        InformalExcessDemand=InformalExcessDemand+InformalDemand[it]-InformalSupply[it];
        FormalExcessDemand=FormalExcessDemand+FormalDemand[it]-FormalSupply[it];
    }
    
    double ExcessTotal=0;
    ExcessTotal=pow(InformalExcessDemand,2)+pow(FormalExcessDemand,2);
    ExcessTotal=ExcessTotal/M;
    //cout << InformalExcessDemand << " InformalExcessDemand "<< endl;
    //cout << FormalExcessDemand << " FormalExcessDemand "<< endl;
    //cout << ExcessTotal << " ExcessTotal"<< endl;
    return(ExcessTotal);
}





// Standardizing the excess demand to be minimized

double StandardizedExcessDemands(const double *Wages, double Others[21]){
    //Loaging everything from others to parameters
    vector<double>WagesVector;
    WagesVector.resize(2);
    WagesVector[0]=Wages[0];
    WagesVector[1]=Wages[1];
    
    vector<double>Params;
    Params.resize(16);
    
    for(int it=0;it<16;it++){
        Params[it]=Others[it];
    }
    
    
    vector<double>InitLWorkers;
    InitLWorkers.resize(2);
    InitLWorkers[0]=Others[15];
    InitLWorkers[1]=Others[16];
    
    
    vector<double> InitProf;
    InitProf.resize(3);
    InitProf[0]=Others[17];
    InitProf[1]=Others[18];
    InitProf[2]=Others[19];
    
    
    //Return the function
    
    
    
    return(pow(ExcessDemandFunctions(WagesVector,Params,InitLWorkers,InitProf),0.5));
}


//11. MinimizeExcessDemands
int iteratExcessDemands=0;

double ExcessDemandsTotal(unsigned n, const double *x, double *grad, void *ExcessDemandsTotal_data){
    //Void to double
    //cout << profitsFinMaxim_data<< " this stuff"<< endl;
    double *Params = (double *)ExcessDemandsTotal_data;
    
    double result=StandardizedExcessDemands(x, Params);
    //cout << result << " evalF"<< endl;
    iteratExcessDemands=iteratExcessDemands+1;
    printf("Iteration=(%d); Feval=%0.10g\n", iteratExcessDemands, result);
    return(result);
}

//12. Creating function to find equilibrium wages conditional on parameters
// [[Rcpp::export]]
vector<double> EqWages(vector <double> Others, vector<double> WagesInit){
    
    //Others: Vector of parameters used to find equilibrium.
    //        Need to reconvert to double.
    
    double DOthers[21];
    for(int it=0; it<20;it++){
        DOthers[it]=Others[it];
    }
    
    
    
    
    //WagesInit: Vector of initial wages to start looking for equilibrium.
    //Need to translate from vector to double
    double Winitial[2];
    Winitial[0]=WagesInit[0];
    Winitial[1]=WagesInit[1];
    
    
    double lbExcessDemmand[2];
    lbExcessDemmand[0]=0;
    lbExcessDemmand[1]=0;
    
    
    nlopt_opt ExcessDemand;
    ExcessDemand = nlopt_create(NLOPT_LN_NELDERMEAD, 2);
    
    nlopt_set_lower_bounds(ExcessDemand, lbExcessDemmand);
    nlopt_set_min_objective(ExcessDemand, ExcessDemandsTotal,(void *)&DOthers);
    
    
    nlopt_set_xtol_rel(ExcessDemand, 1.0e-8);
    nlopt_set_ftol_abs(ExcessDemand,1.0e-8);
    double ExcessDemandValueFinal; /* the minimum objective value, upon return */
    
    
    
    
    
    //A
    
    
    
    
    nlopt_optimize(ExcessDemand, Winitial, &ExcessDemandValueFinal);
    if (nlopt_optimize(ExcessDemand, Winitial, &ExcessDemandValueFinal) < 0) {
        printf("nlopt failed!\n");
    }
    else {
        printf("found minimum at f(%g) = %0.10g\n",Winitial[0], ExcessDemandValueFinal);
    }
    nlopt_destroy(ExcessDemand);
    vector<double> WagesVector;
    WagesVector.resize(3);
    WagesVector[0]=Winitial[0];
    WagesVector[1]=Winitial[1];
    WagesVector[2]=ExcessDemandValueFinal;
    
    return(WagesVector);
    
    
    
}


// [[Rcpp::export]]
arma::vec EqWagesNumericVector(arma::vec Others, arma::vec WagesInit){
    
    
    //Others: Vector of parameters used to find equilibrium.
    //        Need to reconvert to double.
    cout << Others[15] << " Others[15] in "<< endl;
    double DOthers[21];
    for(int it=0; it<20;it++){
        DOthers[it]=Others[it];
        cout << it << " it "<< endl;
        cout << DOthers[it]<< " DOthers[it] in eqwagesnumeric" << endl;
    }
    
    
    
    
    //WagesInit: Vector of initial wages to start looking for equilibrium.
    //Need to translate from vector to double
    double Winitial[2];
    Winitial[0]=WagesInit[0];
    Winitial[1]=WagesInit[1];
    
    
    double lbExcessDemmand[2];
    lbExcessDemmand[0]=0;
    lbExcessDemmand[1]=0;
    
    
    nlopt_opt ExcessDemand;
    ExcessDemand = nlopt_create(NLOPT_LN_NELDERMEAD, 2);
    
    nlopt_set_lower_bounds(ExcessDemand, lbExcessDemmand);
    nlopt_set_min_objective(ExcessDemand, ExcessDemandsTotal,(void *)&DOthers);
    
    //Originally:
    //nlopt_set_xtol_rel(ExcessDemand, 1.0e-8);
    nlopt_set_ftol_abs(ExcessDemand,1.0e-8);
    
    //nlopt_set_stopval(ExcessDemand, 0.1);
    double ExcessDemandValueFinal; /* the minimum objective value, upon return */
    
    
    
    
    
    
    
    
    
    
    nlopt_optimize(ExcessDemand, Winitial, &ExcessDemandValueFinal);
    if (nlopt_optimize(ExcessDemand, Winitial, &ExcessDemandValueFinal) < 0) {
        printf("nlopt failed!\n");
    }
    else {
        printf("found minimum at f(%g) = %0.10g\n",Winitial[0], ExcessDemandValueFinal);
    }
    
    EquilibriumValue.open("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2SobolsGenerated/EquilibriumValue.csv", ios::out | ios::app);
    
    EquilibriumValue<<ExcessDemandValueFinal << " , " ;
    EquilibriumValue << endl;
    EquilibriumValue.close();
    
    nlopt_destroy(ExcessDemand);
    
    arma::vec WagesVector=arma::zeros(3);
    WagesVector[0]=Winitial[0];
    WagesVector[1]=Winitial[1];
    WagesVector[2]=ExcessDemandValueFinal;
    
    return(WagesVector);
    
}




//14. Obtaining the theoretical moments
vector<vector<double> > TheoMoments(arma::vec Others, arma::vec WagesEquilibrium,
                                    arma::vec armaInitLWorkers,arma::vec armaInitProf){
    //Others: Vector of parameters used to determine equilibrium stuff.
    //It is used as arma::vec, translate it to a double vector;
    double DOthers[21];
    for(int it=0; it<20;it++){
        DOthers[it]=Others[it];
    }
    
    
    //WagesInit: Vector of equilibrium wages to
    //Need to translate from vector to double
    
    vector<double> WagesVector;
    WagesVector.resize(2);
    WagesVector[0]=WagesEquilibrium[0];
    WagesVector[1]=WagesEquilibrium[1];
    
    //Initial values to run the optimizers need to be translated from arma vec to double:
    vector<double> InitLWorkers;
    InitLWorkers.resize(2);
    InitLWorkers[0]=armaInitLWorkers[0];
    InitLWorkers[1]=armaInitLWorkers[1];
    
    
    double InitLWorkersDouble[2];
    InitLWorkersDouble[0]=InitLWorkers[0];
    InitLWorkersDouble[1]=InitLWorkers[1];
    
    
    vector<double> InitProf;
    InitProf.resize(3);
    InitProf[0]=armaInitProf[0];
    InitProf[1]=armaInitProf[1];
    InitProf[2]=armaInitProf[2];
    
    double InitProfdouble[3];
    InitProfdouble[0]=InitProf[0];
    InitProfdouble[1]=InitProf[1];
    InitProfdouble[2]=InitProf[2];
    
    //Generate how many draws to obtain the moments
    int M=1000;
    
    
    //Initialize vector of moments
    double InformalDemand[M];
    double InformalSupply[M];
    double FormalDemand[M];
    double FormalSupply[M];
    double Entrepreneurs[M];
    
    for(int it=0; it<M; it++){
        InformalDemand[it]=0;
        InformalSupply[it]=0;
        FormalDemand[it]=0;
        FormalSupply[it]=0;
    }
    
    
    //Initialize the variables to iterate over
    vector<double> Ttheta;
    Ttheta.resize(2);
    
    int Decision=0;
    
    //Params for decission has the following structure:
    vector<double> ParamsDecision;
    ParamsDecision.resize(12);
    ParamsDecision[0]=WagesVector[0];
    ParamsDecision[1]=WagesVector[1];
    ParamsDecision[2]=DOthers[0];
    ParamsDecision[3]=DOthers[1];
    ParamsDecision[4]=DOthers[2];
    ParamsDecision[5]=DOthers[3];
    ParamsDecision[6]=DOthers[4];
    ParamsDecision[7]=DOthers[5];
    ParamsDecision[8]=DOthers[6];
    ParamsDecision[9]=DOthers[7];
    ParamsDecision[10]=DOthers[8];
    ParamsDecision[11]=DOthers[14];
    
    
    
    
    //double m2[] = {0.01, 0.003,  0.003,  0.01};
    double aalpha=DOthers[0];
    double ggamma=DOthers[1];
    double ddelta=DOthers[2];
    double bbeta=DOthers[3];
    double ssigma=DOthers[4];
    double kkappa=DOthers[5];
    double rrho=DOthers[6];
    double psi=DOthers[7];
    double chi=DOthers[8];
    double mmu1=DOthers[9];
    double mmu2=DOthers[10];
    double ssigma1=DOthers[11];
    double ssigma2=DOthers[12];
    double rho12=DOthers[13];
    double c=DOthers[14];
    
    
    cout << aalpha << " aalpha in theomoments!!!"<< endl;
    cout << ddelta << " ddelta in theomoments!!!"<< endl;
    cout << ggamma << " ggamma in theomoments!!!"<< endl;
    cout << bbeta << " bbeta in theomoments!!!"<< endl;
    cout << ssigma << " ssigma in theomoments!!!"<< endl;
    cout << kkappa << " kkappa in theomoments!!!"<< endl;
    cout << rrho << " rrho in theomoments!!!"<< endl;
    cout << psi << " psi in theomoments!!!"<< endl;
    cout << chi << " chi in theomoments!!!"<< endl;
    cout << mmu1 << " mmu1 in theomoments!!!"<< endl;
    cout << mmu2 << " mmu2 in theomoments!!!"<< endl;
    cout << ssigma1 << " ssigma1 in theomoments!!!"<< endl;
    cout << ssigma2 << " ssigma2 in theomoments!!!"<< endl;
    cout << rho12 << " rho12 in theomoments!!!"<< endl;
    cout << c << " c in theomoments!!!"<< endl;
    
    
    
    //Initialize vector of moments and decisions
    vector<vector<double> > DecVector;
    DecVector.resize(4);
    for (int it=0; it<4; it++){
        DecVector[it].resize(4);
    }
    
    //Initializing random number generator for the distribution
    int SEED=2581633;
    typedef boost::mt19937 RNGType;
    RNGType rng(SEED);
    boost::random::uniform_real_distribution<> unidistrib(1,2);
    boost::variate_generator<RNGType, boost::normal_distribution<> >
    generator(rng,
              boost::normal_distribution<>());
    
    double S0rr=gen_normal_3(generator);
    S0rr=gen_normal_3(generator);
    
    //We need to obtain the cholesky decomposition of the variance-covariance matrix
    int nnn=2;
    double m2[] = {ssigma1, rho12,  rho12,  ssigma2};
    double *c2 = cholesky(m2, nnn);
    
    //Obtaining draws from multivariate pareto distribution
    
    //Draws from the normal distribution
    vector<double> z1;
    z1.resize(M);
    vector<double> z2;
    z2.resize(M);
    
    
    //We will generate a vector of indeces for workers and entrepreneurs with the corresponding tthetas
    vector<int> EntreprenIndex;
    EntreprenIndex.resize(M);
    vector<int> WorkerIndex;
    WorkerIndex.resize(M);
    
    vector<double> TthetasEntreprenIndex;
    TthetasEntreprenIndex.resize(M);
    vector<double> TthetaWorkerIndex;
    TthetaWorkerIndex.resize(M);
    
    //And the number of workers and entrepreneurs
    int NumberWorkers=0;
    int NumberEntrep=0;
    
    
    //Moments to obtain returns to wages and production total
    double WageTotal=0;
    double TotalProduction=0;
    
    
    //First we will identify the individuals who are workers or entrepreneurs
    for(int it=0; it<M; it++){
        
        //Obtaining the draws from a standard normal distribution
        z1[it]=gen_normal_3(generator);
        z2[it]=gen_normal_3(generator);
        
        //Doing the corresponding transformation to obtain them from bivariate distribution
        Ttheta[0]=exp(c2[0]*z1[it]+c2[1]*z2[it]+mmu1);
        Ttheta[1]=exp(c2[2]*z1[it]+c2[3]*z2[it]+mmu2);
        
        //Computing the decissions
        DecVector=iDecision(Ttheta, ParamsDecision, InitLWorkers, InitProf);
        Entrepreneurs[it]=DecVector[0][0];
        Decision=DecVector[0][0];
        
        
        
        if(DecVector[0][0]==0){
            TthetaWorkerIndex[NumberWorkers]=Ttheta[0];
            WorkerIndex[NumberWorkers]=it;
            WageTotal+=Ttheta[0]*(ParamsDecision[0]*DecVector[0][2]+ParamsDecision[1]*DecVector[1][2]);
            
            
            NumberWorkers=NumberWorkers+1;
        }
        
        
        if(DecVector[0][0]==1){
            TthetasEntreprenIndex[NumberEntrep]=Ttheta[1];
            EntreprenIndex[NumberEntrep]=it;
            NumberEntrep=NumberEntrep+1;
            TotalProduction+=production(DecVector[0][1],  DecVector[1][1],  aalpha,  Ttheta[1], c);
        }
        
    }

    double propWages=0;
    propWages=WageTotal/TotalProduction;
    cout <<propWages << " propWages "<< endl;
    
    
    //We have now the indeces of those who are workers and those who are entrepreneurs with their corresponding tthetas. Now we need to: 1. Sort them, 2. Generate the corresponding moments.
    //Need to generate pairs between entrepr ad their skils. Sort them by ttheta_e(TthetasEntreprenIndex
    //and then obtain the corresponding indeces.
    //int N=10;
    double NumberOriginalEntrep=NumberEntrep;
    if(NumberEntrep==0){
        NumberEntrep=10;
    }
    vector<pair<double, double> > EntreprenVecPair (NumberEntrep, std::make_pair(0, 0));
    
    for(int it=0; it<NumberEntrep;it++){
        EntreprenVecPair[it].first=TthetasEntreprenIndex[it];
        EntreprenVecPair[it].second=it;
        
    }
    
    vector<pair<double, double> > WorkersVecPair (NumberWorkers, std::make_pair(0, 0));
    
    for(int it=0; it<NumberWorkers;it++){
        WorkersVecPair[it].first=TthetaWorkerIndex[it];
        WorkersVecPair[it].second=it;
        
    }
    
    //It will be sorted by the second element
    
    
    
    
    std::sort(EntreprenVecPair.begin(), EntreprenVecPair.end());
    std::sort(WorkersVecPair.begin(), WorkersVecPair.end());
    
    
    //Sorting the entrepreneurs:
    vector<int> EntrepIndexFinal(NumberEntrep);
    
    
    
    
    
    
    
    
    //Vectors already sorted out.
    // Generating the corresponding deciles:
    vector<int> decileWorkers(9);
    vector<int> decileEntrep(9);
    vector<double>Tthetaw(9);
    vector<double>Tthetae(9);
    
    //Loading the parameters from the stuff
    double wi=WagesEquilibrium[0];
    double wf=WagesEquilibrium[1];
    
    
    //Obtaining the indices of the corresponding deciles:
    for(int it=0;it<9;it++){
        decileWorkers[it]=NumberWorkers*(it+1)/10;
        decileEntrep[it]=NumberEntrep*(it+1)/10;
        Tthetae[it]= EntreprenVecPair[decileEntrep[it]].first;
        Tthetaw[it]= WorkersVecPair[decileWorkers[it]].first;
        cout << Tthetae[it] << " Tthetae[it] in obtaining deciles"<< endl;
        
    }
    
    //Vectors are sorted. Now I will obtain the corresponding taxes payed proportionally.
    double TaxesPayedTotal=0;
    //Entrepreneur
    vector<double> WsolTAXES(NumberEntrep);
    vector<double> ProductionTAXES(NumberEntrep);
    vector<double> TaxesAbsoluteTAXES(NumberEntrep);
    vector<double> PropInformalDemandedTAXES(NumberEntrep);
    vector<double> WorkersTotalDemandedTAXES(NumberEntrep);
    vector<double> TaxesPayed(NumberEntrep);
    for(int it=0; it<NumberEntrep; it++){
        
        
        WsolTAXES=profitsFinMaxim(InitProfdouble,  wi,wf,aalpha,ddelta, ggamma,bbeta,
                                  ssigma, EntreprenVecPair[it].first,c);
        
        //Wsol contains: ni, nf, z, prof
        
        //1. Production
        ProductionTAXES[it]=production(WsolTAXES[0],  WsolTAXES[1],  aalpha,  Tthetae[it],c);

        
        //2. Taxes payed;
        TaxesPayed[it]=TcActual(WsolTAXES[2],WsolTAXES[0],WsolTAXES[1],aalpha,EntreprenVecPair[it].first,wi,wf,c);
        TaxesPayedTotal+=TaxesPayed[it];
    }
    
    //Now I will obtain the proportion of taxes payed in each decile
    double taxesCumulative[10];
    double sumaimpuestos=0;
    //Decile 1.
    for(int it=0; it<decileEntrep[0]; it++){
        sumaimpuestos+=TaxesPayed[it]/TaxesPayedTotal;
        taxesCumulative[0]=sumaimpuestos;

    }
    
    //Decile 2.
    for(int it=decileEntrep[0]+1; it<decileEntrep[1]; it++){
        sumaimpuestos+=TaxesPayed[it]/TaxesPayedTotal;
        taxesCumulative[1]=sumaimpuestos;

    }
    
    //Decile 3.
    for(int it=decileEntrep[1]+1; it<decileEntrep[2]; it++){
        sumaimpuestos+=TaxesPayed[it]/TaxesPayedTotal;
        taxesCumulative[2]=sumaimpuestos;

    }
    
    //Decile 4.
    for(int it=decileEntrep[2]+1; it<decileEntrep[3]; it++){
        sumaimpuestos+=TaxesPayed[it]/TaxesPayedTotal;
        taxesCumulative[3]=sumaimpuestos;

    }
    //Decile 5.
    for(int it=decileEntrep[3]+1; it<decileEntrep[4]; it++){
        sumaimpuestos+=TaxesPayed[it]/TaxesPayedTotal;
        taxesCumulative[4]=sumaimpuestos;

    }
    //Decile 6.
    for(int it=decileEntrep[4]+1; it<decileEntrep[5]; it++){
        sumaimpuestos+=TaxesPayed[it]/TaxesPayedTotal;
        taxesCumulative[5]=sumaimpuestos;

    }
    //Decile 7.
    for(int it=decileEntrep[5]+1; it<decileEntrep[6]; it++){
        sumaimpuestos+=TaxesPayed[it]/TaxesPayedTotal;
        taxesCumulative[6]=sumaimpuestos;

    }
    //Decile 8.
    for(int it=decileEntrep[6]+1; it<decileEntrep[7]; it++){
        sumaimpuestos+=TaxesPayed[it]/TaxesPayedTotal;
        taxesCumulative[7]=sumaimpuestos;

    }
    //Decile 9.
    for(int it=decileEntrep[7]+1; it<decileEntrep[8]; it++){
        sumaimpuestos+=TaxesPayed[it]/TaxesPayedTotal;
        taxesCumulative[8]=sumaimpuestos;

    }
    
   
    
    
    
    
    //Entrepreneur
    vector<double> Wsol(4);
    vector<double> Production(9);
    vector<double> TaxesAbsolute(9);
    vector<double> PropInformalDemanded(9);
    vector<double> WorkersTotalDemanded(9);
    
    
    
    for(int it=0;it<9;it++){
        Wsol=profitsFinMaxim(InitProfdouble,  wi,wf,aalpha,ddelta, ggamma,bbeta,
                             ssigma, Tthetae[it],c);
        
  
        
        //Wsol contains: ni, nf, z, prof
        
        //1. Production
        Production[it]=Tthetae[it]*pow((c+Wsol[0]+Wsol[1]),aalpha);
        
        //2. Taxes payed;
        TaxesAbsolute[it]=TcActual(Wsol[2],Wsol[0],Wsol[1],aalpha,Tthetae[it],wi,wf,c);
        
        //3. Workers employed
        WorkersTotalDemanded[it]=Wsol[0]+Wsol[1];

        
        //4. Proportioninformaldemanded
        PropInformalDemanded[it]=Wsol[0]/(max(WorkersTotalDemanded[it],0.001));
        if(WorkersTotalDemanded[it]==0){
            PropInformalDemanded[it]=1;
        }
        
        cout << " ----"<< endl;
        cout <<Wsol[0] << " informaldemanded "<< endl;
        cout <<  WorkersTotalDemanded[it] << " total labor demand  total"<< endl;
        
    }
    //Taxes payed proportionally
    vector<double> TaxesProportionally(9);
    double sumTaxes=0;
    for(int it=0;it<9;it++){
        sumTaxes+=TaxesAbsolute[it];
    }
    
    for(int it=0; it<9; it++){
        TaxesProportionally[it]=TaxesAbsolute[it]/max(sumTaxes,0.01);
    }
    
    
    
    
    //WorkerMoments.
    vector<double> SolWorker(2);
    vector<double> TotalIncome(9);
    vector<double> InformalLaborSupplyProp(9);
    vector<double> TotalLaborSupply(9);
    
    
    for(int it=0; it<9; it++){
        
        //0. Solving the problem of the worker. Returns informal and formal.
        SolWorker=valueWorkerFinMaxim(InitLWorkersDouble,wf,wi,kkappa,rrho, psi,chi,Tthetaw[it]);
        
        //1. TotalIncome
        TotalIncome[it]=SolWorker[0]*Tthetaw[it]*wi+SolWorker[1]*Tthetaw[it]*wf;
        
        //2. Informal labor supply proportion
        
        InformalLaborSupplyProp[it]=SolWorker[0]/max(SolWorker[0]+SolWorker[1],0.00001);
        
        //3. Total Labor Supply
        TotalLaborSupply[it]=SolWorker[0]+SolWorker[1];
        cout << " ----"<< endl;
        cout <<TotalLaborSupply[it] << " total labor supply "<< endl;
        cout <<InformalLaborSupplyProp[it] << " total labor supply  informal"<< endl;
        
        
    }
    
    //We want labor supply total as proportion of the median
    vector<double> LaborSupplyProportion(9);
    for(int it=0;it<9;it++){
        LaborSupplyProportion[it]=TotalLaborSupply[it]/TotalLaborSupply[4];
    }
    
    
    //The answer vector will return:
    //1. Proportion of workers, entrepreneurs
    //2. Production[it]
    //3. Taxes payed TaxesProportionally[it]
    //4. Workers employed WorkersTotalDemanded[it]
    //5. Proportioninformaldemanded PropInformalDemanded[it]
    
    //6. TotalIncome  TotalIncome[it]
    //7. Informal labor supply proportion InformalLaborSupplyProp[it]
    //8. Total Labor Supply TotalLaborSupply[it]
    
    
    vector<vector<double> > answer(10);
    for(int it=0; it<10; it++){
        answer[it].resize(10);
    }
    //Need to transform into double.
    double doubNumberEntrep=NumberEntrep;
    double doubNumberWorkers=NumberWorkers;
    
    for(int it=0; it<9; it++){
        answer[it][0]=Production[it];
        answer[it][1]=TaxesProportionally[it];
        answer[it][2]=WorkersTotalDemanded[it];
        answer[it][3]=PropInformalDemanded[it];
        answer[it][4]=TotalIncome[it];
        answer[it][5]=InformalLaborSupplyProp[it];
        answer[it][6]=LaborSupplyProportion[it];
        answer[it][7]=doubNumberEntrep/(doubNumberWorkers+doubNumberEntrep);
        answer[it][8]=propWages;
        
        
        
    }
    
    
    if(NumberOriginalEntrep==0){
        answer[0][0]=666;
        answer[0][1]=666;
        answer[0][2]=666;
        answer[0][3]=666;
        answer[0][4]=666;
        answer[0][5]=666;
    }
    
    return(answer);
}



//Now combining moments and equilibrium wages finder
vector<vector<double> > EquilibriumMoments(arma::vec Others, arma::vec WagesInit,
                                           arma::vec armaInitLWorkers,arma::vec armaInitProf){
    
    
    //1. Find the equilibrium wages
    //Testing the moment generating function
    arma::vec WagesEquilibrium=EqWagesNumericVector(Others, WagesInit);

    //2. Once we find the equilibrium wages, obtain the moments.
    vector<vector<double> >answ(9);
    for(int it=0;it<9;it++){
        answ[it].resize(9);
        
    }
    answ=TheoMoments(Others, WagesEquilibrium,
                     armaInitLWorkers,armaInitProf);
    return(answ);
}








//Computing distances from theoretical to empirical moments.

// [[Rcpp::export]]
double DistanceEstimator(arma::vec Others, arma::vec WagesInit,
                         arma::vec armaInitLWorkers,arma::vec armaInitProf){
    
    
    
    
    for(int it=0; it<20; it++){
        cout <<Others[it]<< " othersit in distance estimator"<< endl;
        cout << it << " it "<< endl;
    }
    
    //0. Initializing distance
    double distance=0;
    //1. Obtaining theoretical distance from the parameters:
    vector<vector<double> >Theomoments(10);
    for(int it=0;it<10;it++){
        Theomoments[it].resize(10);
    }
    Theomoments=EquilibriumMoments(Others, WagesInit,
                                   armaInitLWorkers,armaInitProf);
    
    //Computing the distances
    //answer[0]=Production[it];
    //answer[it][1]=TaxesProportionally[it];
    //answer[it][2]=WorkersTotalDemanded[it];
    //answer[it][3]=PropInformalDemanded[it];
    //answer[it][4]=TotalIncome[it];
    //answer[it][5]=InformalLaborSupplyProp[it];
    //answer[it][6]=TotalLaborSupply[it];
    //answer[it][7]=doubNumberEntrep/(doubNumberWorkers+doubNumberEntrep);
    
    //Constructing the empirical moments:
    vector<double> PropEntrepEmpirical(10);
    for(int it=0; it<10; it++){
        PropEntrepEmpirical[it]=0.3625;
    }
    
    
    
    //Total income
    vector<double> IncomeEmpirical(9);
    IncomeEmpirical[0]=82.756/225.682;
    IncomeEmpirical[1]=135.012/225.682;
    IncomeEmpirical[2]=164.702/225.682;
    IncomeEmpirical[3]=197.829/225.682;
    IncomeEmpirical[4]=225.682/225.682;
    IncomeEmpirical[5]=268.271/225.682;
    IncomeEmpirical[6]=318.264/225.682;
    IncomeEmpirical[7]=394.568/225.682;
    IncomeEmpirical[8]=595.318/225.682;
    
    //Transforming theoretical income in proportion in terms of median:
    for(int it=0; it<9; it++){
        Theomoments[it][4]=Theomoments[it][4]/max(Theomoments[4][4],0.001);
        Theomoments[it][0]=Theomoments[it][0]/max(Theomoments[4][0],0.001);
    }
    
    //Production
    vector<double> Production(9);
    
    Production[0]=900.900/5246.640;
    Production[1]=1653.750/5246.640;
    Production[2]=2520.000/5246.640;
    Production[3]=3698.415/5246.640;
    Production[4]=5246.640/5246.640;
    Production[5]=7257.600/5246.640;
    Production[6]=10584.000/5246.640;
    Production[7]=17050.950/5246.640;
    Production[8]=38690.505/5246.640;
    
    
    //Taxes
    vector<double> Taxes(10);
    
    Taxes[0]=0.002457482;
    Taxes[1]=0.006255534;
    Taxes[2]=0.12598472;
    Taxes[3]=0.20355296;
    Taxes[4]=0.32547876;
    Taxes[5]=0.47825881;
    Taxes[6]=0.72803947;
    Taxes[7]=0.116116270;
    Taxes[8]=0.233266590;
    Taxes[9]=1;
    
    //Proportion ifnformal supply
    vector<double> InformalSupplyProportion(9);
    InformalSupplyProportion[0]=0.9082840;
    InformalSupplyProportion[1]=0.8318318;
    InformalSupplyProportion[2]=0.7062315;
    InformalSupplyProportion[3]=0.6227545;
    InformalSupplyProportion[4]=0.5970149;
    InformalSupplyProportion[5]=0.5134328;
    InformalSupplyProportion[6]=0.4464286;
    InformalSupplyProportion[7]=0.3731343;
    InformalSupplyProportion[8]=0.2298507;
    
    //Total labor supply with respect to the median.
    //TAKING THE VERSION CORRESPONDING TO
    vector<double> TotalLaborSupply(10);
    TotalLaborSupply[0]=0.748;
    TotalLaborSupply[1]=0.904;
    TotalLaborSupply[2]=0.991;
    TotalLaborSupply[3]=1.02;
    TotalLaborSupply[4]=1.0000000;
    TotalLaborSupply[5]=1.03;
    TotalLaborSupply[6]=1.05;
    TotalLaborSupply[7]=1.04;
    TotalLaborSupply[8]=1.02;
    TotalLaborSupply[9]=0.997;
    double test=1;
    
    //Open file where moments will be stored
    
    
    //Computing the final distance
    for(int it=0; it<9; it++){
        if(test==1){
            cout << " ------ Moments inside distance "<< endl;
            cout << it << " it "<< endl;
            cout << Production[it] << " Production[it] " << endl; //checked
            cout << Theomoments[it][0] << " Theomoments[it][0] " << endl; //checked
            cout << Taxes[it] << " Taxes[it] " << endl; //Checked
            cout << Theomoments[it][1] << " Theomoments[it][1] " << endl;
            cout << IncomeEmpirical[it] << " IncomeEmpirical[it] " << endl; // Checked
            cout << Theomoments[it][3] << " Theomoments[it][3] demand informal prop" << endl;
            cout << Theomoments[it][4] << " Theomoments[it][4] " << endl; //Checked
            cout << InformalSupplyProportion[it] << " InformalSupplyProportion[it] " << endl; //Checked
            cout << Theomoments[it][5] << " Theomoments[it][5] " << endl;//Checked
            cout << TotalLaborSupply[it] << " TotalLaborSupply[it] " << endl; //Checked
            cout << Theomoments[it][6] << " Theomoments[it][6],2 " << endl; //Labor supply theoretical with respect to median
            cout << PropEntrepEmpirical[it] << " PropEntrepEmpirical[it] " << endl;
            cout << Theomoments[it][7] << " Theomoments[it][7] " << endl; //Proportion entrepreneurs
            
            
            
        }
        
        
        
        distance+=pow((Production[it]-Theomoments[it][0])/(max(0.0001,Production[it])),2);
        
        cout << "contribution distance production "<<pow((Production[it]-Theomoments[it][0])/(max(0.0001,Production[it])),2)<< endl;
        
        distance+=pow((Taxes[it]-Theomoments[it][1])/(max(0.0001,Taxes[it])),2);
        
        cout << "contribution taxes  "<<pow((Taxes[it]-Theomoments[it][1])/(max(0.0001,Taxes[it])),2)<< endl;
        
        
        distance+=pow((IncomeEmpirical[it]-Theomoments[it][4])/(max(0.0001,IncomeEmpirical[it])),2);
        cout << "contribution IncomeEmpirical  "<<pow((IncomeEmpirical[it]-Theomoments[it][4])/(max(0.0001,IncomeEmpirical[it])),2)<< endl;
        
        distance+=pow((InformalSupplyProportion[it]-Theomoments[it][5])/(max(0.0001,InformalSupplyProportion[it])),2);
        cout << "contribution InformalSupplyProportion  "<<pow((InformalSupplyProportion[it]-Theomoments[it][5])/(max(0.0001,InformalSupplyProportion[it])),2)<< endl;
        
        distance+=pow((TotalLaborSupply[it]-Theomoments[it][6])/(max(0.0001,TotalLaborSupply[it])),2);
        cout << "contribution TotalLaborSupply  "<<pow((TotalLaborSupply[it]-Theomoments[it][6])/(max(0.0001,TotalLaborSupply[it])),2)<< endl;
        
        distance+=pow((PropEntrepEmpirical[it]-Theomoments[it][7])/(max(0.0001,PropEntrepEmpirical[it])),2);
        
        cout << "contribution PropEntrepEmpirical  "<<pow((PropEntrepEmpirical[it]-Theomoments[it][7])/(max(0.0001,PropEntrepEmpirical[it])),2)<< endl;
        
        //I will add a part corresponding to the 0.8 approximately retribution to work. Need to adapt.
        distance+=pow((0.75-Theomoments[it][8])/0.75,2);
        
        cout << "contribution PropEntrepEmpirical  "<<pow((0.75-Theomoments[it][8])/0.75,2)<< endl;
        
    }
    
    
    
    
    
    //ThMomentsCSV.open("ThMomentsCSV.csv", ios::out | ios::app);
    
    
    //ThMomentsCSV << ", " << endl;
    //ThMomentsCSV << " " << endl;
    
    //Moments relative to the proportion of informality in firm correspond to
    //4th, 7th and 8th decile. These proportions are 99.6, 97, and 91.2%.
    distance+=pow((Theomoments[3][3]-0.996)/(0.996),2);
    distance+=pow((Theomoments[6][3]-0.97)/(0.912),2);
    distance+=pow((Theomoments[7][3]-0.912)/(0.912),2);
    distance+=pow((Theomoments[8][3]-0.001)/(0.001),2);
    
    cout << " pow((Theomoments[3][3]-0.996)/(0.996),2) "<< pow((Theomoments[3][3]-0.996)/(0.996),2) << endl;
    cout << " pow((Theomoments[6][3]-0.996)/(0.996),2) "<< pow((Theomoments[6][3]-0.97)/(0.97),2) << endl;
    cout << " pow((Theomoments[7][3]-0.996)/(0.996),2) "<< pow((Theomoments[7][3]-0.91)/(0.91),2) << endl;
    cout << " pow((Theomoments[8][3]-0.996)/(0.996),2) "<< pow((Theomoments[8][3]-0.01)/(0.01),2) << endl;
    
    
    //Writing the csv file of the moments
    ThMoments0CSV.open("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolsGenerated/ThMoments0CSV.csv", ios::out | ios::app);
    ThMoments1CSV.open("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolsGenerated/ThMoments1CSV.csv", ios::out | ios::app);
    ThMoments2CSV.open("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolsGenerated/ThMoments2CSV.csv", ios::out | ios::app);
    ThMoments3CSV.open("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolsGenerated/ThMoments3CSV.csv", ios::out | ios::app);
    ThMoments4CSV.open("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolsGenerated/ThMoments4CSV.csv", ios::out | ios::app);
    ThMoments5CSV.open("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolsGenerated/ThMoments5CSV.csv", ios::out | ios::app);
    ThMoments6CSV.open("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolsGenerated/ThMoments6CSV.csv", ios::out | ios::app);
    ThMoments7CSV.open("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolsGenerated/ThMoments7CSV.csv", ios::out | ios::app);
    ThMoments8CSV.open("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolsGenerated/ThMoments8CSV.csv", ios::out | ios::app);
    for(int it=0; it<9; it++){
        ThMoments0CSV<<Theomoments[it][0] << " , " ;
        ThMoments1CSV<<Theomoments[it][1] << " , " ;
        ThMoments2CSV<<Theomoments[it][2] << " , " ;
        ThMoments3CSV<<Theomoments[it][3] << " , " ;
        ThMoments4CSV<<Theomoments[it][4] << " , " ;
        ThMoments5CSV<<Theomoments[it][5] << " , " ;
        ThMoments6CSV<<Theomoments[it][6] << " , " ;
        ThMoments7CSV<<Theomoments[it][7] << " , " ;
        ThMoments8CSV<<Theomoments[it][8] << " , " ;
        cout << Theomoments[it][8] << "Theomoments[it][8]"<< endl;
        
    }
    
    ThMoments0CSV << endl;
    ThMoments0CSV.close();
    
    ThMoments1CSV << endl;
    ThMoments1CSV.close();
    
    ThMoments2CSV << endl;
    ThMoments2CSV.close();
    
    ThMoments3CSV << endl;
    ThMoments3CSV.close();
    
    ThMoments4CSV << endl;
    ThMoments4CSV.close();
    
    ThMoments5CSV << endl;
    ThMoments5CSV.close();
    
    ThMoments6CSV << endl;
    ThMoments6CSV.close();
    
    ThMoments7CSV << endl;
    ThMoments7CSV.close();
    
    ThMoments8CSV << endl;
    ThMoments8CSV.close();
    
    
    
    //ThMomentsCSV << " " << endl;
    //ThMomentsCSV << endl << endl;
    //ThMomentsCSV.close();
    
    //Writing the csv file of the parameters
    //ParametersCSV.close();
    
    
    //Loading the parameters into a csv file
    
    
    
    
    ParametersCSV.open("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolsGenerated/ParametersCSV.csv", ios::out | ios::app);
    cout << " loading parameters in the csv file0"<< endl;
    //Printing the parameters
    for(int it=0; it<15; it++){
        cout <<Others[it] << "PARAAAAAAAAAAMAMAMA" << endl;
        ParametersCSV<<Others[it] << " , " ;
    }
    
    ParametersCSV << endl;
    ParametersCSV.close();
    cout << distance << " distance "<< endl;
    
    
    
    
    
    
    DistanceMoments.open("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolsGenerated/DistanceMoments.csv", ios::out | ios::app);
    
    DistanceMoments<<distance << " , " ;
    DistanceMoments << endl;
    DistanceMoments.close();
    
    
    return(distance);
}


//STANDARDIZING THE distance estimator to get inputs as parameters, variables
double StandardizedDistanceEstimator(const double *Parameters, double AdditionalVars[8]){
    
    //First, loading the parameters
    double ggamma=Parameters[0];
    double ddelta=Parameters[1];
    double bbeta=Parameters[2];
    double ssigma=Parameters[3];
    double kkappa=Parameters[4];
    double psi=Parameters[5];
    double chi=Parameters[6];
    double rrho=Parameters[7];
    double mmu1=Parameters[8];
    double mmu2=Parameters[9];
    double ssigma1=Parameters[10];
    double ssigma2=Parameters[11];
    double rho12=Parameters[12];
    double c=Parameters[13];
    double aalpha=Parameters[14];
    
    cout << aalpha << " aalpha in standardized distance"<< endl;
    cout << ggamma << " ggamma in standardized distance"<< endl;
    cout << ddelta << " ddelta in standardized distance"<< endl;
    cout << mmu1 << " mmu1 in standardized distance"<< endl;
    cout << c << " c in standardized distance"<< endl;
    cout << psi << " psi in standardized distance"<< endl;
    cout << chi << " chi in standardized distance"<< endl;
    cout << rho12 << " rho12 in standardized distance"<< endl;
    cout << ssigma1 << " ssigma1 in standardized distance"<< endl;
    cout << kkappa << " kkappa in standardized distance"<< endl;
    
    //Now loading the remaining of the variables necessary for the analysis
    //double aalpha=AdditionalVars[0]; Loading before aalpha was a parameter
    
    //Wages
    arma::vec WagesVector=arma::zeros(2);
    WagesVector[0]=AdditionalVars[0];
    WagesVector[1]=AdditionalVars[1];
    
    //Initializing optimizers of workers and entrepreneurs
    arma::vec InitLWorkers=arma::zeros(2);
    InitLWorkers[0]=AdditionalVars[2];
    InitLWorkers[1]=AdditionalVars[3];
    
    //Initializing entrepreneurs moments
    arma::vec InitProf=arma::zeros(3);
    InitProf[0]=AdditionalVars[4];
    InitProf[1]=AdditionalVars[5];
    InitProf[2]=AdditionalVars[6];
    
    
    //The necessary parameters to be estimated are:
    //Constructing the vector [Others]
    arma::vec Others=arma::zeros(20);
    
    Others[0]=aalpha;
    Others[1]=ggamma;
    Others[2]=ddelta;
    Others[3]=bbeta;
    Others[4]=ssigma;
    Others[5]=kkappa;
    Others[6]=rrho;
    Others[7]=psi;
    Others[8]=chi;
    Others[9]=mmu1;
    Others[10]=mmu2;
    Others[11]=ssigma1;
    Others[12]=ssigma2;
    Others[13]=rho12;
    Others[14]=c;
    Others[15]=InitLWorkers[0];
    Others[16]=InitLWorkers[1];
    Others[17]=InitProf[0];
    Others[18]=InitProf[1];
    Others[19]=InitProf[2];
    
    
    
    
    
    double distance=DistanceEstimator(Others,WagesVector,InitLWorkers,InitProf);
    cout << distance << " distance standardized 1!!"<< endl;
    return(distance);
}







//11. MinimizeExcessDemands
int iteratMinimizeDistance=0;
double MinimizeDistanceTranslate(unsigned n, const double *x, double *grad, void *MinimizeDistance_data){
    //Void to double
    //cout << profitsFinMaxim_data<< " this stuff"<< endl;
    cout << " minimizedistance translate0000 "<< endl;
    double *Params = (double *)MinimizeDistance_data;
    cout << " minimizedistance translate "<< endl;
    
    
    double result=StandardizedDistanceEstimator(x,Params);
    cout << " minimizedistance translate1 "<< endl;
    cout << result << " evalF"<< endl;
    iteratMinimizeDistance=iteratMinimizeDistance+1;
    printf("Iteration=(%d); Feval=%0.10g\n", iteratMinimizeDistance, result);
    return(result);
}






//Finally, function such that running it will simply minimize the situation
arma::vec MinimizingDistance(arma::vec Parameters, arma::vec AdditionalVars){
    
    //First, we need to load the elements in the Parameters vector.
    
    double ggamma=Parameters[0];
    double ddelta=Parameters[1];
    double bbeta=Parameters[2];
    double ssigma=Parameters[3];
    double kkappa=Parameters[4];
    double psi=Parameters[5];
    double chi=Parameters[6];
    double rrho=Parameters[7];
    double mmu1=Parameters[8];
    double mmu2=Parameters[9];
    double ssigma1=Parameters[10];
    double ssigma2=Parameters[11];
    double rho12=Parameters[12];
    double c=Parameters[13];
    double aalpha=Parameters[14];
    
    
    cout << aalpha << " aaalpha in minimizingdistance"<< endl;
    cout << ggamma << " ggamma in minimizingdistance"<< endl;
    cout << c << " c in minimizingdistance"<< endl;
    
    cout << AdditionalVars[0]<< " AdditionalVars[0] "<< endl;
    cout << AdditionalVars[1]<< " AdditionalVars[1] "<< endl;
    cout << AdditionalVars[2]<< " AdditionalVars[2] "<< endl;
    cout << AdditionalVars[3]<< " AdditionalVars[3] "<< endl;
    cout << AdditionalVars[4]<< " AdditionalVars[4] "<< endl;
    cout << AdditionalVars[5]<< " AdditionalVars[5] "<< endl;
    cout << AdditionalVars[6]<< " AdditionalVars[6] "<< endl;
    cout << AdditionalVars[7]<< " AdditionalVars[7] "<< endl;
    
    //We also need to translate it to a double to use the minimizer
    double ParamDoubles[15]={};
    ParamDoubles[0]=ggamma;
    ParamDoubles[1]=ddelta;
    ParamDoubles[2]=bbeta;
    ParamDoubles[3]=ssigma;
    ParamDoubles[4]=kkappa;
    ParamDoubles[5]=psi;
    ParamDoubles[6]=chi;
    ParamDoubles[7]=rrho;
    ParamDoubles[8]=mmu1;
    ParamDoubles[9]=mmu2;
    ParamDoubles[10]=ssigma1;
    ParamDoubles[11]=ssigma2;
    ParamDoubles[12]=rho12;
    ParamDoubles[13]=c;
    ParamDoubles[14]=aalpha;
    
    
    //Establishing upper bounds
    double ubParameters[15];
    ubParameters[0]=10;
    ubParameters[1]=700;
    ubParameters[2]=700;
    ubParameters[3]=10;
    ubParameters[4]=500;
    ubParameters[5]=4;
    ubParameters[6]=600;
    ubParameters[7]=10;
    ubParameters[8]=4;
    ubParameters[9]=4;
    ubParameters[10]=4.5;
    ubParameters[11]=4.5;
    ubParameters[12]=1;
    ubParameters[13]=100;
    ubParameters[14]=1;
    
    //Establishing the lower bounds
    double lbParameters[15];
    lbParameters[0]=0.01;
    lbParameters[1]=0.1;
    lbParameters[2]=0.000000001;
    lbParameters[3]=0.000000001;
    lbParameters[4]=0.000000001;
    lbParameters[5]=0.000000001;
    lbParameters[6]=0.000000001;
    lbParameters[7]=0.00000001;
    lbParameters[8]=0.00000001;
    lbParameters[9]=0.00000001;
    lbParameters[10]=0.00000001;
    lbParameters[11]=0.00000001;
    lbParameters[12]=0.00000001;
    lbParameters[13]=0.00000001;
    lbParameters[14]=0.5;

    
    //Putting the Other parameters in the double format
    double doubAdditionalPar[8];
    doubAdditionalPar[0]=AdditionalVars[0];
    doubAdditionalPar[1]=AdditionalVars[1];
    doubAdditionalPar[2]=AdditionalVars[2];
    doubAdditionalPar[3]=AdditionalVars[3];
    doubAdditionalPar[4]=AdditionalVars[4];
    doubAdditionalPar[5]=AdditionalVars[5];
    doubAdditionalPar[6]=AdditionalVars[6];
    doubAdditionalPar[7]=AdditionalVars[7];


    
    
    for(int it=0; it<9; it++){
        cout <<doubAdditionalPar[it]<< " doubAdditionalPar[it]"<< endl;
    }
    
    //Start creating the optimizer object
    nlopt_opt MinimizeDistance;
    MinimizeDistance = nlopt_create(NLOPT_LN_NELDERMEAD, 15);
    
    nlopt_set_lower_bounds(MinimizeDistance, lbParameters);
    nlopt_set_upper_bounds(MinimizeDistance, ubParameters);
    nlopt_set_min_objective(MinimizeDistance, MinimizeDistanceTranslate,(void *)&doubAdditionalPar);
    
    nlopt_set_ftol_abs(MinimizeDistance,0.1);
    double distanceValue;
    
    cout << " algo mas"<< endl;
    nlopt_optimize(MinimizeDistance, ParamDoubles, &distanceValue); //Not running, why?
    cout << " here trying hard"<< endl;
    if (nlopt_optimize(MinimizeDistance, ParamDoubles, &distanceValue) < 0) {
        printf("nlopt failed!\n");
    }
    else {
        printf("found minimum at f(%g) = %0.10g\n",ParamDoubles[0], distanceValue);
    }
    nlopt_destroy(MinimizeDistance);
    vector<double> Parfinal(15);
    for(int it=0; it<14; it++){
        Parfinal[it]=ParamDoubles[it];
    }
    
    //And loading the distance//Note: why did i do this?
    Parfinal[13]=distanceValue;
    
    
    
    
    return(Parfinal);
}


void SobolRun(arma::vec WagesVectorIn,
              arma::vec InitLWorkersDecision,
              arma::vec InitProfDecision){
    double a=2;
    cout << a << " a "<< endl;
    
    std::ifstream theFile ("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolsGenerated/SobolDim15.csv");
    
    int SIZEOBS=20000;
    int NVAR=15;
    double MYARRAY[SIZEOBS][NVAR];
    // ...
    
    std::string line;
    std::vector<std::vector<std::string> > values;
    int it=0;
    int it2=0;
    std::string line_value;
    std::vector<std::string> line_values;
    std::stringstream ss;
    while(std::getline(theFile, line))
    {
        
        ss<<line;
        //std::stringstream ss(line);
        //std::string item;
        //cout <<  << "linevalprev"<<endl;
        while(std::getline(ss, line_value, ','))
        {
            line_values.push_back(line_value);
            MYARRAY[it][it2] = ::atof(line_value.c_str());
            //cout << MYARRAY[it][it2]<< " hahaha"<< endl;
            //MYARRAY[it][it2]=std::stod (line_value); only for c++11 compi
            it2=it2+1;
            if (it2==NVAR){ //later change 4 for
                it2=0;
            }
        }
        values.push_back(line_values);
        
        //For c++11 used values.emplace_back(line_values);
        
        //cout << line_value<< "line_value2"<< endl;
        it=it+1;
        //Free the string types
        line_value.clear();
        line_values.clear();
        ss.clear();
        
    }
    
    
    
    
    
    
    vector<double> VectorOthers(21);
    VectorOthers[15]=InitLWorkersDecision[0];
    VectorOthers[16]=InitLWorkersDecision[1];
    VectorOthers[17]=InitProfDecision[0];
    VectorOthers[18]=InitProfDecision[1];
    VectorOthers[19]=InitProfDecision[2];
    
    cout << " loading the situation"<< endl;
    vector<double>PARLOAD(16);
#pragma omp parallel for shared(WagesVectorIn, InitLWorkersDecision,InitProfDecision,VectorOthers) private(PARLOAD)
    for(int it=0; it<10000; it++){
        cout << " testing the loading sobol "<< endl;
        for(int par=0; par<15;par++){
            PARLOAD[par]=MYARRAY[it][par];
        }
        
        for(int it=0; it<15; it++){
            
            VectorOthers[it]=PARLOAD[it];
            cout << it << " it "<< endl;
            cout <<VectorOthers[it]<< " VectorOthers[it]"<< endl;
        }
        
        for(int it=0; it<20; it++){
            cout << it << " it "<< endl;
            cout <<VectorOthers[it]<< " VectorOthers[it]INSOBOL"<< endl;
        }
        
        DistanceEstimator(VectorOthers, WagesVectorIn,InitLWorkersDecision,InitProfDecision);
    }
    
}



int main(int argc, const char * argv[]) {
    
    
    
    
    
    
    
    vector<double> A(3);
    vector<double> V(3);
    V[0]=0;
    V[1]=1;
    V[2]=2;
    A[0]=50;
    A[1]=40;
    A[2]=30;
    pair <double, double> make_pair(920, 12); ;
    
    
    int N=10;
    vector<pair<double, double> > myVec (N, std::make_pair(0, 0));
    myVec[0].second=-10;
    myVec[1].second=0.9;
    myVec[2].second=10.8;
    myVec[3].second=1.8;
    
    myVec[0].first=19.2;
    myVec[1].first=-0.2;
    myVec[2].first=-1.2;
    myVec[3].first=2.2;
    
    cout <<myVec[0].second << " myVec[0].second"<< endl;
    cout <<myVec[1].second << " myVec[1].second"<< endl;
    cout <<myVec[2].second << " myVec[2].second"<< endl;
    cout <<myVec[3].second << " myVec[3].second"<< endl;
    
    cout <<myVec[0].first << " myVec[0].first"<< endl;
    cout <<myVec[1].first << " myVec[1].first"<< endl;
    cout <<myVec[2].first << " myVec[2].first"<< endl;
    cout <<myVec[3].first << " myVec[3].first"<< endl;
    
    std::sort(myVec.begin(), myVec.end());
    
    cout <<myVec[0].second << " myVec[0].second"<< endl;
    cout <<myVec[1].second << " myVec[1].second"<< endl;
    cout <<myVec[2].second << " myVec[2].second"<< endl;
    cout <<myVec[3].second << " myVec[3].second"<< endl;
    
    cout <<myVec[0].first << " myVec[0].first"<< endl;
    cout <<myVec[1].first << " myVec[1].first"<< endl;
    cout <<myVec[2].first << " myVec[2].first"<< endl;
    cout <<myVec[3].first << " myVec[3].first"<< endl;
    
    
    int x=0;
    
    //Testing functions
    //Setting parameters to test
    
    
    
    
    
    
    double aalpha=0.8;
    double tthetae=12.84209;
    double tthetaw=3.873585;
    double wi=8.81;
    double wf=8.13;
    
    double ni=10;
    double nf=50;
    double ggamma=0.28;
    double ddelta=0.12;
    double bbeta=0.15;
    double ssigma=0.2;
    double kkappa=0.1;
    double psi=0.4;
    double chi=2.5;
    double rrho=0.9;
    
    
    
    double lf=0.5;
    double li=0.5;
    double z=24;
    double mmu1=0.2;
    double mmu2=1.3;
    double ssigma1=0.1;
    double ssigma2=0.9;
    double rho12=0.08;
    double c=10;
    
    
    
    

    
    //These were the original parameters. Trying different combinations:
    
    aalpha=0.8;
    ggamma=0.198503;
    ddelta=0.180862;
    bbeta=0.247681;
    ssigma=0.468333;
    kkappa=0.0689548;
    psi=1.19302;
    chi=1.09956;
    rrho=2.74542;
    mmu1=0.240092;
    mmu2=2.50935;
    ssigma1=0.113429;
    ssigma2=1.16547;
    rho12=0.10219;
    c=0;
    
    
    //Final definition
    
    aalpha=0.682324;
    ggamma=6.82592;
    ddelta=647.281;
    bbeta=324.671;
    ssigma=8.25678;
    kkappa=263.156;
    psi=3.17603;
    chi=592.294;
    rrho=3.01674;
    mmu1=0.420156;
    mmu2=3.92797;
    ssigma1=1.51879;
    ssigma2=4.14574;
    rho12=0.0519922;
    c=1;
    
    //0. Payroll taxes marginal
    double Payrolltest=Tn(nf);
    cout << Payrolltest<< " Tn'(nf) "<< endl;
    
    //1. Payroll taxes definition
    double PayrolltestActual=TnActual(nf);
    cout << PayrolltestActual<< " Tn(nf)"<< endl;
    
    //2. Pre-tax profits
    double profmTest=profm( ni,  nf,  aalpha,  tthetae,  wi,  wf,c);
    cout << profmTest<< " profmTest"<< endl;
    
    //3. Corporate tax profits, marginal
    double TcTest=Tc(z,ni,nf,aalpha,tthetae,wi,wf,c);
    cout << TcTest << " TcTest " << endl;
    
    
    //4. Tc Actual corporate taxes
    double TcActualTest=TcActual(z,ni,nf,aalpha,tthetae,wi,wf,c);
    cout << TcActualTest << " TcActualtest "<< endl;
    
    //5. PIT
    double PITtest=PIT(tthetaw,wf,lf);
    cout << PITtest << " PITtest " << endl;
    
    //6. PIT Marginal
    double PITtestMarginal=PITM(tthetaw,wf,lf);
    cout << PITtestMarginal<< " PITtestMarginal "<< endl;
    
    //7. FinProfits
    double Params[9]={};
    
    Params[0]=wi;
    Params[1]=wf;
    Params[2]=aalpha;
    Params[3]=ddelta;
    Params[4]=ggamma;
    Params[5]=bbeta;
    Params[6]=ssigma;
    Params[7]=tthetae;
    Params[8]=c;
    //int v=0;
    //void value = *(double *)Params;
    //void *p=&Params;
    //static_cast<void*>(&a)
    
    double Args[3]={};
    Args[0]=ni;
    Args[1]=nf;
    Args[2]=z;
    
    cout << z << " z entering finprofits"<< endl;
    double FinProfitsTest=FinProfits(Args, Params);
    cout << FinProfitsTest << " FinProfitsTest "<< endl;
    
    //8-9. Maximize profits
    double lb[3];
    lb[0]=0;
    lb[1]=0;
    lb[2]=0;
    
    cout << " p0"<< endl;
    nlopt_opt optProfits;
    optProfits = nlopt_create(NLOPT_LN_NELDERMEAD, 3);
    
    nlopt_set_lower_bounds(optProfits, lb);
    nlopt_set_max_objective(optProfits, profitsFinMaxim,(void *)&Params);
    
    
    nlopt_set_xtol_rel(optProfits, 1.0e-8);
    nlopt_set_ftol_abs(optProfits,1.0e-8);
    double minf3; /* the minimum objective value, upon return */
    double x3[3];
    x3[0]=ni;
    x3[1]=nf;
    x3[2]=z;
    
    //memcpy(Params, (x3**)pointer, sizeof Params);
    
    //cout <<FinProfits(x3,Params)<< " test inicial"<< endl;
    //cout << x3[2]<< " x[2]"<< endl;
    nlopt_optimize(optProfits, x3, &minf3);
    if (nlopt_optimize(optProfits, x3, &minf3) < 0) {
        printf("nlopt failed!\n");
    }
    else {
        printf("found minimum at f(%g) = %0.10g\n", x3[0], minf3);
    }
    
    nlopt_destroy(optProfits);
    //cout << x3[0]<< " optimal ni"<< endl;
    //cout << x3[1]<< " optimal nf "<< endl;
    //cout << x3[2]<< " optimal z"<< endl;
    
    //Testing maximum functino of profits
    vector <double> MaxProf;
    MaxProf.resize(4);
    double InitialCond[3];
    InitialCond[0]=ni;
    InitialCond[1]=nf;
    InitialCond[2]=z;
    
    cout << " before maxprof"<< endl;
    
    MaxProf=profitsFinMaxim( InitialCond,  wi,wf,aalpha,ddelta,ggamma, bbeta,ssigma,tthetae,c);
    //Finding the actual optimal values found
    cout << MaxProf[0] << " Maxprof[0]"<< endl;
    cout << MaxProf[1] << " Maxprof[1]"<< endl;
    cout << MaxProf[2] << " Maxprof[2]"<< endl;
    cout << MaxProf[3] << " Maxprof[3]"<< endl;
    
    
    //10. Value of workers
    double ParamWorkers[7]={};
    
    ParamWorkers[0]=wf;
    ParamWorkers[1]=wi;
    ParamWorkers[2]=kkappa;
    ParamWorkers[3]=rrho;
    ParamWorkers[4]=psi;
    ParamWorkers[5]=chi;
    ParamWorkers[6]=tthetaw;
    
    //int v=0;
    //void value = *(double *)Params;
    //void *p=&Params;
    //static_cast<void*>(&a)
    
    double ArgWorkers[2]={};
    ArgWorkers[0]=li;
    ArgWorkers[1]=lf;
    double ValWorkersT=ValueWorkers(ArgWorkers, ParamWorkers);
    cout << ValWorkersT << " ValWorkersT "<< endl;
    
    //Testing maximization of value of workers
    vector <double> ValAnsWorker;
    ValAnsWorker.resize(3);
    double InitialCondWorker[2];
    InitialCond[0]=li;
    InitialCond[1]=lf;
    
    
    ValAnsWorker=valueWorkerFinMaxim(InitialCondWorker,wf,wi, kkappa,rrho,psi, chi, tthetaw);
    cout << ValAnsWorker[0] << " Valansworker0"<< endl;
    cout << ValAnsWorker[1] << " Valansworker1"<< endl;
    cout << ValAnsWorker[2] << " Valansworker2"<< endl;
    
    
    
    
    //Testing the decision function
    vector<vector<double> > Decision;
    Decision.resize(4);
    for (int it=0; it<4; it++){
        Decision[it].resize(4);
    }
    //Defining the inputs
    vector<double> ParamsDecision;
    ParamsDecision.resize(12);
    
    ParamsDecision[0]=wi;
    ParamsDecision[1]=wf;
    ParamsDecision[2]=aalpha;
    ParamsDecision[3]=ggamma;
    ParamsDecision[4]=ddelta;
    ParamsDecision[5]=bbeta;
    ParamsDecision[6]=ssigma;
    ParamsDecision[7]=kkappa;
    ParamsDecision[8]=rrho;
    ParamsDecision[9]=psi;
    ParamsDecision[10]=chi;
    ParamsDecision[11]=c;
    
    
    
    
    
    
    
    
    
    
    
    
    vector<double> InitLWorkersDecision;
    InitLWorkersDecision.resize(2);
    InitLWorkersDecision[0]=li;
    InitLWorkersDecision[1]=lf;
    vector<double> InitProfDecision;
    InitProfDecision.resize(3);
    InitProfDecision[0]=ni;
    InitProfDecision[1]=nf;
    InitProfDecision[2]=z;
    vector<double>TthetaDecision;
    TthetaDecision.resize(2);
    TthetaDecision[0]=tthetaw;
    TthetaDecision[1]=tthetae;
    
    
    vector<double> Wages;
    Wages.resize(2);
    Wages[0]=wi;
    Wages[1]=wf;
    
    
    cout << "-----------" << endl;
    
    //CSV files were initially open before the definition of functions. This is the moment to close them
    ParametersCSV.close();
    ThMoments0CSV.close();
    ThMoments1CSV.close();
    ThMoments2CSV.close();
    ThMoments3CSV.close();
    ThMoments4CSV.close();
    ThMoments5CSV.close();
    ThMoments6CSV.close();
    ThMoments7CSV.close();
    ThMoments8CSV.close();
    DistanceMoments.close();
    EquilibriumValue.close();
    
    
    vector<double> WagesInEqWages;
    WagesInEqWages.resize(2);
    double Winitial[2];
    
    Winitial[0]=wi;
    Winitial[1]=wf;
    
    WagesInEqWages[0]=Winitial[0];
    WagesInEqWages[1]=Winitial[1];
    
    //Define vector for excessDemandFunctions
    vector<double> ParamsDecisionExcessDemand;
    ParamsDecisionExcessDemand.resize(15);
    ParamsDecisionExcessDemand[0]=aalpha;
    ParamsDecisionExcessDemand[1]=ggamma;
    ParamsDecisionExcessDemand[2]=ddelta;
    ParamsDecisionExcessDemand[3]=bbeta;
    ParamsDecisionExcessDemand[4]=ssigma;
    ParamsDecisionExcessDemand[5]=kkappa;
    ParamsDecisionExcessDemand[6]=rrho;
    ParamsDecisionExcessDemand[7]=psi;
    ParamsDecisionExcessDemand[8]=chi;
    ParamsDecisionExcessDemand[9]=mmu1;
    ParamsDecisionExcessDemand[10]=mmu2;
    ParamsDecisionExcessDemand[11]=ssigma1;
    ParamsDecisionExcessDemand[12]=ssigma2;
    ParamsDecisionExcessDemand[13]=rho12;
    ParamsDecisionExcessDemand[14]=c;
    
    double Others[21];
    for(int it=0;it<15;it++){
        Others[it]=ParamsDecisionExcessDemand[it];
    }
    
    cout << InitProfDecision[2]<< "InitProfDecision[2] "<< endl;
    Others[14]=InitLWorkersDecision[0];
    Others[15]=InitLWorkersDecision[1];
    Others[16]=InitProfDecision[0];
    Others[17]=InitProfDecision[1];
    Others[18]=InitProfDecision[2];
    
    //Checking what does other contain
    for(int it=0; it<20; it++){
        cout << " -----"<< endl;
        cout << it << " it"<< endl;
        cout << Others[it]<< " others it"<< endl;
    }
    
    
    
    vector<double> VOthers;
    VOthers.resize(20);
    for(int i=0;i<20;i++){
        VOthers[i]=Others[i];
    }
    
    //Vector of initial Wages
    arma::vec WagesVectorIn=arma::zeros(2);
    WagesVectorIn.resize(2);
    
    WagesVectorIn[0]=wi;
    WagesVectorIn[1]=wf;
    
    WagesVectorIn[0]=8.18169;
    WagesVectorIn[1]=8.53154;
    
    
    //Vector of final answer
    vector<double> InitProf;
    InitProf.resize(3);
    InitProf[0]=Others[11];
    InitProf[1]=Others[12];
    InitProf[2]=Others[13];
    double WageEx[2];
    WageEx[0]=wi;
    WageEx[1]=wf;
    
    
    
    Wages[0]=WageEx[0];
    Wages[1]=WageEx[1];
    
    WageEx[0]=6.9425;
    WageEx[1]=7.69912;
    
    //cout << StandardizedExcessDemands(WageEx,  Others) << " Standardized1" << endl;
    cout << " ---"<< endl;
    WageEx[0]=6.9425;
    WageEx[1]=7.69912;
    //cout << StandardizedExcessDemands(WageEx,  Others) << " modified" << endl;
    arma::vec WagesEquilibrium=arma::zeros(2);
    
    WageEx[0]=6.94398;
    WageEx[1]=7.699953;
    //cout << StandardizedExcessDemands(WageEx,  Others) << " modified2" << endl;
    
    
    //Vector of others
    arma::vec VecOthers=arma::zeros(20);
    for(int it=0; it<20; it++){
        VecOthers[it]=Others[it];
    }
    
    
    //cout << " running distance estimator "<< endl;
    Decision=iDecision(TthetaDecision,ParamsDecision,InitLWorkersDecision,InitProfDecision);
    //DistanceEstimator(VecOthers, WagesVectorIn,InitLWorkersDecision,InitProfDecision);
    
    
    
    //Running the minimizing distance estimator
    arma::vec MinimizeDistanceParameters(15);
    arma::vec AdditionalVAriables(8);
    
    //Loading the parameters first. Modified to 15 elements to include aalpha.
    double ParStandardDistance[15];
    ParStandardDistance[0]=ggamma;
    ParStandardDistance[1]=ddelta;
    ParStandardDistance[2]=bbeta;
    ParStandardDistance[3]=ssigma;
    ParStandardDistance[4]=kkappa;
    ParStandardDistance[5]=psi;
    ParStandardDistance[6]=chi;
    ParStandardDistance[7]=rrho;
    ParStandardDistance[8]=mmu1;
    ParStandardDistance[9]=mmu2;
    ParStandardDistance[10]=ssigma1;
    ParStandardDistance[11]=ssigma2;
    ParStandardDistance[12]=rho12;
    ParStandardDistance[13]=c;
    ParStandardDistance[14]=aalpha;
    
    //Loading additional parameters
    double AddPar[7];
    AddPar[0]=wi;
    AddPar[1]=wf;
    AddPar[2]=InitLWorkersDecision[0];
    AddPar[3]=InitLWorkersDecision[1];
    AddPar[4]=InitProfDecision[0];
    AddPar[5]=InitProfDecision[1];
    AddPar[6]=InitProfDecision[2];
    cout <<InitProfDecision[0] << " InitProfDecision[0]"<< endl;
    cout <<InitProfDecision[1] << " InitProfDecision[1]"<< endl;
    cout <<InitProfDecision[2] << " InitProfDecision[2]"<< endl;
    //Loading the additional vectors

    
    
    
    
    
    
    
    
    
    aalpha=0.682324;
    ggamma=6.82592;
    ddelta=647.281;
    bbeta=324.671;
    ssigma=8.25678;
    kkappa=263.156;
    psi=3.17603;
    chi=532.294;
    rrho=3.01674;
    mmu1=0.420156;
    mmu2=3.92797;
    ssigma1=1.51879;
    ssigma2=4.14574;
    rho12=0.0519922;
    c=1;
    
    
    ParamsDecision[0]=wi;
    ParamsDecision[1]=wf;
    ParamsDecision[2]=aalpha;
    ParamsDecision[3]=ggamma;
    ParamsDecision[4]=ddelta;
    ParamsDecision[5]=bbeta;
    ParamsDecision[6]=ssigma;
    ParamsDecision[7]=kkappa;
    ParamsDecision[8]=rrho;
    ParamsDecision[9]=psi;
    ParamsDecision[10]=chi;
    ParamsDecision[11]=c;
    
    ParamsDecisionExcessDemand.resize(15);
    ParamsDecisionExcessDemand[0]=aalpha;
    ParamsDecisionExcessDemand[1]=ggamma;
    ParamsDecisionExcessDemand[2]=ddelta;
    ParamsDecisionExcessDemand[3]=bbeta;
    ParamsDecisionExcessDemand[4]=ssigma;
    ParamsDecisionExcessDemand[5]=kkappa;
    ParamsDecisionExcessDemand[6]=rrho;
    ParamsDecisionExcessDemand[7]=psi;
    ParamsDecisionExcessDemand[8]=chi;
    ParamsDecisionExcessDemand[9]=mmu1;
    ParamsDecisionExcessDemand[10]=mmu2;
    ParamsDecisionExcessDemand[11]=ssigma1;
    ParamsDecisionExcessDemand[12]=ssigma2;
    ParamsDecisionExcessDemand[13]=rho12;
    ParamsDecisionExcessDemand[14]=c;
    
    ParStandardDistance[0]=ggamma;
    ParStandardDistance[1]=ddelta;
    ParStandardDistance[2]=bbeta;
    ParStandardDistance[3]=ssigma;
    ParStandardDistance[4]=kkappa;
    ParStandardDistance[5]=psi;
    ParStandardDistance[6]=chi;
    ParStandardDistance[7]=rrho;
    ParStandardDistance[8]=mmu1;
    ParStandardDistance[9]=mmu2;
    ParStandardDistance[10]=ssigma1;
    ParStandardDistance[11]=ssigma2;
    ParStandardDistance[12]=rho12;
    ParStandardDistance[13]=c;
    ParStandardDistance[14]=aalpha;
    
    for(int it=0;it<15;it++){
        Others[it]=ParamsDecisionExcessDemand[it];
    }
    
    cout << ParamsDecisionExcessDemand[14] << " ParamsDecisionExcessDemand[14]=c; "<< endl;
    
    cout << InitProfDecision[2]<< "InitProfDecision[2] "<< endl;
    Others[15]=InitLWorkersDecision[0];
    Others[16]=InitLWorkersDecision[1];
    Others[17]=InitProfDecision[0];
    Others[18]=InitProfDecision[1];
    Others[19]=InitProfDecision[2];
    
    
    for(int it=0; it<20; it++){
        VecOthers[it]=Others[it];
    }
    
    for(int it=0; it<15; it++){
        MinimizeDistanceParameters[it]=ParStandardDistance[it];
    }
    
    
    for(int it=0; it<8; it++){
        AdditionalVAriables[it]=AddPar[it];
    }
    
    
    cout << c<< " ccccc"<< endl;
    cout << ParamsDecisionExcessDemand[14]<< " ParamsDecisionExcessDemand[14]"<< endl;
    cout << Others[14] << " others14"<< endl;
    cout << VecOthers[14]<< " vecothers14"<< endl;
    
    cout << "theomomtn"<< endl;
    
    cout << " muajajajajajajaja"<< endl;
    
    //InitLWorkersDecision.resize(2);
    //wi
    //wf
    //InitLWorkersDecision[0]=li;
    //InitLWorkersDecision[1]=lf;
    //InitProfDecision.resize(3);
    //InitProfDecision[0]=ni;
    //InitProfDecision[1]=nf;
    //InitProfDecision[2]=z;
    
    
    //DistanceEstimator(VecOthers, WagesVectorIn,InitLWorkersDecision,InitProfDecision);
    cout << " muajajajajajajaja222222"<< endl;
    MinimizingDistance(MinimizeDistanceParameters,AdditionalVAriables);
    cout << " muajajajajajajaja222222"<< endl;
    //
    //EquilibriumMoments(VecOthers, WagesVectorIn,InitLWorkersDecision,InitProfDecision);
    //EqWagesNumericVector(VecOthers, WagesVectorIn);
    //StandardizedDistanceEstimator(ParStandardDistance,  AddPar);
    
    //StandardizedExcessDemands(WageEx,  Others);
    
    //ExcessDemandFunctions(Wages,ParamsDecisionExcessDemand, InitLWorkersDecision, InitProfDecision);
    
    cout << " changonorrea"<< endl;
    //Decision=iDecision(TthetaDecision,ParamsDecision,InitLWorkersDecision,InitProfDecision);
    
    
    SobolRun(Wages,InitLWorkersDecision,InitProfDecision);
    //WagesEquilibrium=EqWagesNumericVector(VecOthers, WagesVectorIn);
    //DistanceEstimator(VecOthers, WagesVectorIn,InitLWorkersDecision,InitProfDecision);
    //EquilibriumMoments(VecOthers, WagesVectorIn,InitLWorkersDecision,InitProfDecision);
    
    
    cout << "end of theomomtn"<< endl;
    cout << " minimizing the distance!!!!!!!!! finally"<< endl;
    cout << " minimizing the distance!!!!!!!!! finally"<< endl;
    cout << " minimizing the distance!!!!!!!!! finally"<< endl;
    cout << " minimizing the distance!!!!!!!!! finally"<< endl;
    cout << " minimizing the distance!!!!!!!!! finally"<< endl;
    cout << " minimizing the distance!!!!!!!!! finally"<< endl;
    
    
    cout << MinimizeDistanceParameters[0] << " MinimizeDistanceParameters[0]"<< endl;
    
    
    
    
    
    //MinimizingDistance(MinimizeDistanceParameters,AdditionalVAriables);
    //cout << " running sobolrun "<< endl;
    //SobolRun(Wages,InitLWorkersDecision,InitProfDecision);
    
    //Decision=iDecision(TthetaDecision,ParamsDecision,InitLWorkersDecision,InitProfDecision);
    
    cout << Decision[0][0]<< " Decision "<< endl;
    cout << Decision[0][1]<< " ni "<< endl;
    cout << Decision[1][1]<< " nf "<< endl;
    cout << Decision[2][1]<< " z "<< endl;
    cout << Decision[3][1]<< " prof "<< endl;
    cout << Decision[0][2]<< " li "<< endl;
    cout << Decision[1][2]<< " lf "<< endl;
    cout << Decision[2][2]<< " VWorker "<< endl;
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    cout << " NONSTANDARD------------- "<< endl;
    //cout << ExcessDemandFunctions(Wages,ParamsDecisionExcessDemand, InitLWorkersDecision, InitProfDecision)<<  endl;
    
    
    cout << " STANDARD------------- "<< endl;
    //cout << StandardizedExcessDemands(WageEx,  Others) << " Standardized1" << endl;
    //cout << ExcessDemandFunctions(Wages,ParamsDecisionExcessDemand, InitLWorkersDecision, InitProfDecision)<<  endl;
    //cout << ExcessDemandFunctions(Wages,ParamsDecisionExcessDemand, InitLWorkersDecision, InitProfDecision)<<  endl;
    //cout << ExcessDemandFunctions(Wages,ParamsDecisionExcessDemand, InitLWorkersDecision, InitProfDecision)<<  endl;
    
    //Minimizing excess demands squared
    //ExcessDemandsTotal(unsigned n, const double *x, double *grad, void *ExcessDemandsTotal_data)
    double lbExcessDemmand[2];
    lbExcessDemmand[0]=0;
    lbExcessDemmand[1]=0;
    
    
    nlopt_opt ExcessDemand;
    ExcessDemand = nlopt_create(NLOPT_LN_NELDERMEAD, 2);
    
    nlopt_set_lower_bounds(ExcessDemand, lbExcessDemmand);
    nlopt_set_min_objective(ExcessDemand, ExcessDemandsTotal,(void *)&Others);
    
    
    nlopt_set_xtol_rel(ExcessDemand, 1.0e-8);
    nlopt_set_ftol_abs(ExcessDemand,1.0e-8);
    double ExcessDemandValueFinal; /* the minimum objective value, upon return */
    
    nlopt_destroy(ExcessDemand);
    
    
    
    
    
    
    
    //nlopt_optimize(ExcessDemand, Winitial, &ExcessDemandValueFinal);
    //if (nlopt_optimize(ExcessDemand, Winitial, &ExcessDemandValueFinal) < 0) {
    //    printf("nlopt failed!\n");
    //}
    //else {
    //    printf("found minimum at f(%g) = %0.10g\n",Winitial[0], ExcessDemandValueFinal);
    //}
    
    //Test to write a csv file
    
    
    
    //Trying the Excess Demand functions at wi, wf
    
    //Vector de ExcessDemandFunctions es diferente del vector de idecision. Vector de Excesss
    //demand function empieza en aalpha!
    
    
    cout << Winitial[0] << " wi eq "<< endl;
    cout << Winitial[1] << " wf eq" << endl;
    
    
    //Trying the function of optimizers
    vector<double> EqWagesVector;
    EqWagesVector.resize(2);
    
    
    
    
    //Vector de ExcessDemandFunctions es diferente del vector de idecision.
    
    TthetaDecision[0]=2.5;
    TthetaDecision[1]=0.67;
    
    cout << InitProfDecision[0] << "InitProfDecision[0]"<<endl;
    cout << InitProfDecision[1] << "InitProfDecision[1]"<<endl;
    cout << InitProfDecision[2] << "InitProfDecision[2]"<<endl;
    InitProfDecision[2]=2;
    
    
    
    
    cout << Decision[0][0]<< " Decision "<< endl;
    cout << Decision[0][1]<< " ni "<< endl;
    cout << Decision[1][1]<< " nf "<< endl;
    cout << Decision[2][1]<< " z "<< endl;
    cout << Decision[3][1]<< " prof "<< endl;
    cout << Decision[0][2]<< " li "<< endl;
    cout << Decision[1][2]<< " lf "<< endl;
    cout << Decision[2][2]<< " VWorker "<< endl;
    
    
    
    cout << " obtaining equilibrium wages"<< endl;
    
    //WagesEquilibrium=EqWagesNumericVector(VecOthers, WagesVectorIn);
    
    cout << " end of obtaining equilibrium wages"<< endl;
    //Testing the moment generating function
    cout << "Testing theoretical moments "<< endl;
    //TheoMoments(ParamsDecisionExcessDemand, WagesVectorIn,InitLWorkersDecision,InitProfDecision);
    
    cout << " end of theoretical moments test"<< endl;
    
    
    //Testing the distance estimators
    
    
    
    //
    //DistanceEstimator(VecOthers, WagesVectorIn,InitLWorkersDecision,InitProfDecision);
    
    cout << " initializing equilibrium moments "<< endl;
    //EquilibriumMoments(VecOthers, WagesVectorIn,InitLWorkersDecision,InitProfDecision);
    
    cout << " finalizing equilibrium moments "<< endl;
    
    //Testing standardized distance estimtor loading parameters
    cout << " Testing standardized distance estimtor loading parameters"<< endl;
    cout << " HEEEREEEEEE      "<< endl;
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    cout << WagesEquilibrium[0] << " WagesEquilibrium[0]"<< endl;
    cout << WagesEquilibrium[1] << " WagesEquilibrium[1]"<< endl;
    
    
    
    
    
    
    
    Wages[0]=wi;
    Wages[1]=wf;
    
    
    
    
    int n = 2;
    double m2[] = {0.2, 0.02,  0.02,  0.2};
    double *c2 = cholesky(m2, n);
    show_matrix(c2, n);
    
    cout << " -----"<< endl;
    cout << c2[0] << " c[0][0]"<< endl;
    cout << c2[1] << " c[1][0]"<< endl;
    cout << c2[2] << " c[2][0]"<< endl;
    cout << c2[3] << " c[3][0]"<< endl;
    
    
    free(c2);
    
    
    cout << "Hello, World!\n";
    return 0;
}
