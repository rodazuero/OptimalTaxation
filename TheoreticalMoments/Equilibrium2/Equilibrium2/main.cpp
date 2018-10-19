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
#include <iostream>
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
//#include <Rcpp.h> -> Not necessary if rcpparmadillo included
#include <RcppArmadillo.h>
//#include <bits/stdc++.h>

using std::vector;
using namespace std;
using namespace Rcpp;


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

//2. Pre-tax profits.Checked
long double profm(double ni, double nf, double aalpha, double tthetae, double wi, double wf){
    double pi1=tthetae*pow((ni+nf),aalpha)-wi*ni-wf*nf-TnActual(nf*wf);
    return(pi1);
    
}

//3. Corporate tax profits, marginal
long double Tc(double z, double ni, double nf, double aalpha,double tthetae, double wi, double wf){
    double profm=tthetae*pow((ni+nf),aalpha)-wi*ni-wf*nf-TnActual(nf*wf);
    double arg=profm-z;
    double firsterm=-exp(-arg);
    //If negative profits, set zero marginal rate
    if (profm<0){
        firsterm=0;
    }
    return(firsterm);
}

//4. Tc Actual corporate taxes
//4. Tc Actual corporate taxes
long double TcActual(double z, double ni, double nf, double aalpha, double tthetae, double wi, double wf){
    double profm=tthetae*pow((ni+nf),aalpha)-wi*ni-wf*nf-TnActual(nf*wf);
    double arg=profm-z;
    double ans=0;
    double tax=0;
    double prod=tthetae*pow((ni+nf),aalpha);
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
long double FinProfits(const double *Args, double paramvec[8]){
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
    
    
    double ni=Args[0];
    double nf=Args[1];
    double z=Args[2];
    
    //Operational
    double term1=profm(ni, nf, aalpha, tthetae, wi, wf);
    
    //Corporate taxes
    double taxes=TcActual(z,ni,nf,aalpha,tthetae,wi,wf);
    
    //Cost of evacion
    double evcost=pow(z,1+ssigma)*(bbeta/(1+ssigma));
    
    //Evasion costs
    double infcost=pow(ni,1+ggamma)*(ddelta/(1+ggamma));
    
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
                               double tthetae){
    
    double Params[8]={};
    Params[0]=wi;
    Params[1]=wf;
    Params[2]=aalpha;
    Params[3]=ddelta;
    Params[4]=ggamma;
    Params[5]=bbeta;
    Params[6]=ssigma;
    Params[7]=tthetae;
    
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
    double ddelta=Params[3];
    double ggamma=Params[4];
    double bbeta=Params[5];
    double ssigma=Params[6];
    double kkappa=Params[7];
    double rrho=Params[8];
    double psi=Params[9];
    double chi=Params[10];
    
    
    
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
    ansProfits=profitsFinMaxim(InitialVProfits,  wi,wf, aalpha,ddelta,ggamma,bbeta,ssigma,tthetae);
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
    
    
    
    
    
    
    int M=10000;
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
    ParamsDecision.resize(11);
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
    cout << InformalExcessDemand << " InformalExcessDemand "<< endl;
    cout << FormalExcessDemand << " FormalExcessDemand "<< endl;
    cout << ExcessTotal << " ExcessTotal"<< endl;
    return(ExcessTotal);
}





// Standardizing the excess demand to be minimized

double StandardizedExcessDemands(const double *Wages, double Others[20]){
    //Loaging everything from others to parameters
    vector<double>WagesVector;
    WagesVector.resize(2);
    WagesVector[0]=Wages[0];
    WagesVector[1]=Wages[1];
    
    vector<double>Params;
    Params.resize(15);
    for(int it=0;it<15;it++){
        Params[it]=Others[it];
    }
    
    
    vector<double>InitLWorkers;
    InitLWorkers.resize(2);
    InitLWorkers[0]=Others[14];
    InitLWorkers[1]=Others[15];
    
    
    vector<double> InitProf;
    InitProf.resize(3);
    InitProf[0]=Others[16];
    InitProf[1]=Others[17];
    InitProf[2]=Others[18];
    
    
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
    cout << result << " evalF"<< endl;
    iteratExcessDemands=iteratExcessDemands+1;
    printf("Iteration=(%d); Feval=%0.10g\n", iteratExcessDemands, result);
    return(result);
}

//12. Creating function to find equilibrium wages conditional on parameters
// [[Rcpp::export]]
vector<double> EqWages(vector <double> Others, vector<double> WagesInit){
    
    //Others: Vector of parameters used to find equilibrium.
    //        Need to reconvert to double.
    
    double DOthers[20];
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
    
    
    
    
    
    
    
    
    
    
    nlopt_optimize(ExcessDemand, Winitial, &ExcessDemandValueFinal);
    if (nlopt_optimize(ExcessDemand, Winitial, &ExcessDemandValueFinal) < 0) {
        printf("nlopt failed!\n");
    }
    else {
        printf("found minimum at f(%g) = %0.10g\n",Winitial[0], ExcessDemandValueFinal);
    }
    
    vector<double> WagesVector;
    WagesVector.resize(2);
    WagesVector[0]=Winitial[0];
    WagesVector[1]=Winitial[1];
    
    return(WagesVector);
}


// [[Rcpp::export]]
arma::vec EqWagesNumericVector(arma::vec Others, arma::vec WagesInit){
    
    
    //Others: Vector of parameters used to find equilibrium.
    //        Need to reconvert to double.
    cout << Others[15] << " Others[15] in "<< endl;
    double DOthers[19];
    for(int it=0; it<19;it++){
        DOthers[it]=Others[it];
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
    
    
    
    arma::vec WagesVector=arma::zeros(2);
    WagesVector[0]=Winitial[0];
    WagesVector[1]=Winitial[1];
    
    return(WagesVector);
    
}



int main(int argc, const char * argv[]) {
    //Testing functions
    //Setting parameters to test
    double aalpha=0.8;
    double tthetae=12.84209;
    double tthetaw=3.873585;
    double wi=6.81;
    double wf=7.13;

    double ni=2.3;
    double nf=0.53*70;
    double ggamma=0.28;
    double ddelta=0.12;
    double bbeta=0.15;
    double ssigma=0.2;
    double kkappa=0.1;
    double psi=0.4;
    double chi=1.5;
    double rrho=0.9;
    double lf=2.1;
    double li=2.1;
    double z=24;
    double mmu1=0.2;
    double mmu2=1.3;
    double ssigma1=0.1;
    double ssigma2=0.9;
    double rho12=0.08;
    
    
    //0. Payroll taxes marginal
    double Payrolltest=Tn(nf);
    cout << Payrolltest<< " Tn'(nf) "<< endl;
    
    //1. Payroll taxes definition
    double PayrolltestActual=TnActual(nf);
    cout << PayrolltestActual<< " Tn(nf)"<< endl;
    
    //2. Pre-tax profits
    double profmTest=profm( ni,  nf,  aalpha,  tthetae,  wi,  wf);
    cout << profmTest<< " profmTest"<< endl;
    
    //3. Corporate tax profits, marginal
    double TcTest=Tc(z,ni,nf,aalpha,tthetae,wi,wf);
    cout << TcTest << " TcTest " << endl;
    
    
    //4. Tc Actual corporate taxes
    double TcActualTest=TcActual(z,ni,nf,aalpha,tthetae,wi,wf);
    cout << TcActualTest << " TcActualtest "<< endl;
    
    //5. PIT
    double PITtest=PIT(tthetaw,wf,lf);
    cout << PITtest << " PITtest " << endl;
    
    //6. PIT Marginal
    double PITtestMarginal=PITM(tthetaw,wf,lf);
    cout << PITtestMarginal<< " PITtestMarginal "<< endl;
    
    //7. FinProfits
    double Params[8]={};
    
    Params[0]=wi;
    Params[1]=wf;
    Params[2]=aalpha;
    Params[3]=ddelta;
    Params[4]=ggamma;
    Params[5]=bbeta;
    Params[6]=ssigma;
    Params[7]=tthetae;
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
    
    MaxProf=profitsFinMaxim( InitialCond,  wi,wf,aalpha,ddelta,ggamma, bbeta,ssigma,tthetae);
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
    ParamsDecision.resize(11);
    
    ParamsDecision[0]=wi;
    ParamsDecision[1]=wf;
    ParamsDecision[2]=aalpha;
    ParamsDecision[3]=ddelta;
    ParamsDecision[4]=ggamma;
    ParamsDecision[5]=bbeta;
    ParamsDecision[6]=ssigma;
    ParamsDecision[7]=kkappa;
    ParamsDecision[8]=rrho;
    ParamsDecision[9]=psi;
    ParamsDecision[10]=chi;
    
    
    
    
    
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
    
    
    
    Decision=iDecision(TthetaDecision,ParamsDecision,InitLWorkersDecision,InitProfDecision);
    
    cout << Decision[0][0]<< " Decision "<< endl;
    cout << Decision[0][1]<< " ni "<< endl;
    cout << Decision[1][1]<< " nf "<< endl;
    cout << Decision[2][1]<< " z "<< endl;
    cout << Decision[3][1]<< " prof "<< endl;
    cout << Decision[0][2]<< " li "<< endl;
    cout << Decision[1][2]<< " lf "<< endl;
    cout << Decision[2][2]<< " VWorker "<< endl;
    
    vector<double> Wages;
    Wages.resize(2);
    Wages[0]=wi;
    Wages[1]=wf;
    
    //Define vector for excessDemandFunctions
    vector<double> ParamsDecisionExcessDemand;
    ParamsDecisionExcessDemand.resize(14);
    ParamsDecisionExcessDemand[0]=aalpha;
    ParamsDecisionExcessDemand[1]=ddelta;
    ParamsDecisionExcessDemand[2]=ggamma;
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
    
    
    
    //ExcessDemandFunctions(Wages,ParamsDecisionExcessDemand, InitLWorkersDecision, InitProfDecision);
    
    
    
    double Others[20];
    for(int it=0;it<14;it++){
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
    double Winitial[2];
    
    Winitial[0]=wi;
    Winitial[1]=wf;
    
    
    
    
    
    
    
    //nlopt_optimize(ExcessDemand, Winitial, &ExcessDemandValueFinal);
    //if (nlopt_optimize(ExcessDemand, Winitial, &ExcessDemandValueFinal) < 0) {
    //    printf("nlopt failed!\n");
    //}
    //else {
    //    printf("found minimum at f(%g) = %0.10g\n",Winitial[0], ExcessDemandValueFinal);
    //}
    
    
    
    //Trying the Excess Demand functions at wi, wf
    
    //Vector de ExcessDemandFunctions es diferente del vector de idecision. Vector de Excesss
    //demand function empieza en aalpha!
    
    
    cout << Winitial[0] << " wi eq "<< endl;
    cout << Winitial[1] << " wf eq" << endl;
    
    
    //Trying the function of optimizers
    vector<double> EqWagesVector;
    EqWagesVector.resize(2);
    
    //Vector of others
    arma::vec VecOthers=arma::zeros(20);
    for(int it=0; it<20; it++){
        VecOthers[it]=Others[it];
    }
    
    //Vector of initial Wages
    arma::vec WagesVectorIn=arma::zeros(2);
    WagesVectorIn.resize(2);
    
    WagesVectorIn[0]=wi;
    WagesVectorIn[1]=wf;
    
    //Vector de ExcessDemandFunctions es diferente del vector de idecision.
    
    TthetaDecision[0]=2.5;
    TthetaDecision[1]=0.67;
    
    cout << InitProfDecision[0] << "InitProfDecision[0]"<<endl;
    cout << InitProfDecision[1] << "InitProfDecision[1]"<<endl;
    cout << InitProfDecision[2] << "InitProfDecision[2]"<<endl;
    InitProfDecision[2]=2;
    Decision=iDecision(TthetaDecision,ParamsDecision,InitLWorkersDecision,InitProfDecision);
    
    cout << Decision[0][0]<< " Decision "<< endl;
    cout << Decision[0][1]<< " ni "<< endl;
    cout << Decision[1][1]<< " nf "<< endl;
    cout << Decision[2][1]<< " z "<< endl;
    cout << Decision[3][1]<< " prof "<< endl;
    cout << Decision[0][2]<< " li "<< endl;
    cout << Decision[1][2]<< " lf "<< endl;
    cout << Decision[2][2]<< " VWorker "<< endl;
    
    //Vector of final answer
    
    WageEx[0]=6.9425;
    WageEx[1]=7.69912;
    cout << StandardizedExcessDemands(WageEx,  Others) << " Standardized1" << endl;
    cout << " ---"<< endl;
    WageEx[0]=6.9425;
    WageEx[1]=7.69912;
    cout << StandardizedExcessDemands(WageEx,  Others) << " modified" << endl;
    arma::vec WagesEquilibrium=arma::zeros(2);
    
    WageEx[0]=6.94398;
    WageEx[1]=7.699953;
    cout << StandardizedExcessDemands(WageEx,  Others) << " modified2" << endl;
    
    
    vector<double> WagesInEqWages;
    WagesInEqWages.resize(2);
    WagesInEqWages[0]=Winitial[0];
    WagesInEqWages[1]=Winitial[1];
    
    vector<double> VOthers;
    VOthers.resize(20);
    for(int i=0;i<20;i++){
        VOthers[i]=Others[i];
    }
    
    
    //vector<double> WagesEq=EqWages(VOthers, WagesInEqWages);
    //cout << WagesEq[0]<< " wi "<< endl;
    //cout << WagesEq[1]<< " wi "<< endl;
    
    
    WagesEquilibrium=EqWagesNumericVector(VecOthers, WagesVectorIn);
    
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
