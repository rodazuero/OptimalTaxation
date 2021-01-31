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
#include <omp.h>
#include <random>
#include <unistd.h>
#include <nlopt.hpp>
#include <utility>
#include <chrono>
//#include <Rcpp.h> -> Not necessary if rcpparmadillo included
#include <armadillo>
#include <iterator>
//#include <bits/stdc++.h>
#include "interpolation.h"
#include "interpolation.cpp"
#include <cstdio>
#include <cstdlib>
#include "ap.cpp"
#include "alglibinternal.cpp"
#include "alglibmisc.cpp"
//#include "linalg.cpp"
#include "specialfunctions.cpp"
#include "dataanalysis.cpp"
#include "diffequations.cpp"
#include "fasttransforms.cpp"
#include "integration.cpp"
#include "linalg.cpp"
#include "optimization.cpp"
#include "solvers.cpp"
//#include "specialfunctions.cpp"
#include "statistics.cpp"



using std::vector;
using namespace std::chrono;
using namespace std;
//using namespace Rcpp;

ofstream ThMoments0CSV("/home/ec2-user/InAws/ThMoments0CSV.csv");
ofstream ThMoments1CSV("/home/ec2-user/InAws/ThMoments1CSV.csv");
ofstream ThMoments2CSV("/home/ec2-user/InAws/ThMoments2CSV.csv");
ofstream ThMoments3CSV("/home/ec2-user/InAws/ThMoments3CSV.csv");
ofstream ThMoments4CSV("/home/ec2-user/InAws/ThMoments4CSV.csv");
ofstream ThMoments5CSV("/home/ec2-user/InAws/ThMoments5CSV.csv");
ofstream ThMoments6CSV("/home/ec2-user/InAws/ThMoments6CSV.csv");
ofstream ThMoments7CSV("/home/ec2-user/InAws/ThMoments7CSV.csv");
ofstream ThMoments8CSV("/home/ec2-user/InAws/ThMoments8CSV.csv");
ofstream ParametersCSV("/home/ec2-user/InAws/ParametersCSV.csv");
ofstream DistanceMoments("/home/ec2-user/InAws/DistanceMoments.csv");
ofstream EquilibriumValue("/home/ec2-user/InAws/EquilibriumValue.csv");
ofstream AllocationsOP("/home/ec2-user/InAws/AllocationsOP.csv");


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

//Import the planner´s output to check implementability

vector<vector<double>> Planner(int simdata) {

    std::ifstream theFile("/home/ec2-user/InAws/Taxes_info.csv");
    int SIZEPLANNER = 1000;
    int INPUTS = 11;
    //vector<vector<double>> PLANNERINFO(SIZEPLANNER);
    //for (int it = 0; it < SIZEPLANNER; it++) {
    //    PLANNERINFO[it].resize(INPUTS);
    //}
    vector<vector<double>> PLANNERINFO(SIZEPLANNER, vector<double>(INPUTS));
    //vector<vector<double>> PLANNERINFO[SIZEPLANNER][INPUTS];
    // ...

    std::string line;
    std::vector<std::vector<std::string> > values;
    int it = 0;
    int it2 = 0;
    std::string line_value;
    std::vector<std::string> line_values;
    std::stringstream ss;
    while (std::getline(theFile, line))
    {

        ss << line;
        //std::stringstream ss(line);
        //std::string item;
        //cout <<  << "linevalprev"<<endl;
        while (std::getline(ss, line_value, ';'))
        {
            line_values.push_back(line_value);
            PLANNERINFO[it][it2] = ::atof(line_value.c_str());
            //cout << MYARRAY[it][it2]<< " hahaha"<< endl;
            //MYARRAY[it][it2]=std::stod (line_value); only for c++11 compi
            it2 = it2 + 1;
            if (it2 == INPUTS) {
                it2 = 0;
                //cout << "Linea" << it << endl;
            }
        }
        values.push_back(line_values);

        //For c++11 used values.emplace_back(line_values);

        //cout << line_value<< "line_value2"<< endl;
        it = it + 1;
        //Free the string types
        line_value.clear();
        line_values.clear();
        ss.clear();
        if (it == SIZEPLANNER) {
            cout << "Ultima" << PLANNERINFO[SIZEPLANNER - 1][0] << PLANNERINFO[SIZEPLANNER - 1][10] << endl;
        }
    }

    cout << "Theta w in " << SIZEPLANNER << ":" << PLANNERINFO[SIZEPLANNER - 1][0] << endl;
    cout << "Theta e in " << SIZEPLANNER << ":" << PLANNERINFO[SIZEPLANNER - 1][1] << endl;
    cout << "Theta w in " << 1 << ":" << PLANNERINFO[1][0] << endl;
    cout << "Theta e in " << 1 << ":" << PLANNERINFO[1][1] << endl;

    if (1 == 2) {
        int m = SIZEPLANNER;
        int n = INPUTS;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {

                // Prints ' ' if j != n-1 else prints '\n'           
                cout << PLANNERINFO[i][j] << " \n"[j == n - 1];
            }
        }
    }
    return(PLANNERINFO);
}

// A. Create splines for taxes functions

//A.1.1. Set the spline for Tn
alglib::spline1dinterpolant setSplineTn(vector<vector<double>> taxesmatrix) {
    
    std::vector<double> Xtn(999), Ytn(999);

    for (int it = 0; it < 999; it++) {
        Xtn[it] = taxesmatrix[it + 1][2];
        Ytn[it] = taxesmatrix[it + 1][3];
    }

    alglib::real_1d_array AXtn, AYtn;
    AXtn.setcontent(Xtn.size(), &(Xtn[0]));
    AYtn.setcontent(Ytn.size(), &(Ytn[0]));

    alglib::spline1dinterpolant splineTn;
    
    //1. Linear Spline
    //alglib::spline1dbuildlinear(AXtn, AYtn, Xtn.size(), splineTn);

    //2. Cubic spline (natural- second derivatives in bounds = 0)
    //alglib::spline1dbuildcubic(AXtn, AYtn, Xtn.size(), 2, 0.0, 2, 0.0, splineTn);

    //3. Hermite spline (includes first derivative)
    std::vector<double> Dtn(999);

    for (int it = 0; it < 999; it++) {
        Dtn[it] = taxesmatrix[it + 1][4];
    }

    alglib::real_1d_array ADtn;
    ADtn.setcontent(Dtn.size(), &(Dtn[0]));
    
    alglib::spline1dbuildhermite(AXtn, AYtn, ADtn, Xtn.size(), splineTn);

   
    if (1==2) {

        cout << "-------Spline for Tn------" << endl;

        for (size_t i = 0; i < 20; i++) {
            printf("%f %f\n", Xtn[i], Ytn[i]);
        }
        
        printf("\n");
        for (int i = 500; i < 520; i++) {
            double x = i;
            printf("%f %f\n", x, alglib::spline1dcalc(splineTn, x));
        }
        printf("\n");

    }

    cout << "-------End Spline for Tn------" << endl;

    return(splineTn);

}

//A.1.2 Tn Spline
long double TnSpline(alglib::spline1dinterpolant splineTn,double nomina) {
    double myspline = alglib::spline1dcalc(splineTn, nomina);
    return(myspline);
}

//A.2.1. Set the spline for Tc
alglib::spline1dinterpolant setSplineTc(vector<vector<double>> taxesmatrix) {

    std::vector<double> Xtc(999), Ytc(999);

    for (int it = 0; it < 999; it++) {
        Xtc[it] = taxesmatrix[it + 1][5];
        Ytc[it] = taxesmatrix[it + 1][6];
    }

    alglib::real_1d_array AXtc, AYtc;
    AXtc.setcontent(Xtc.size(), &(Xtc[0]));
    AYtc.setcontent(Ytc.size(), &(Ytc[0]));

    alglib::spline1dinterpolant splineTc;

    //1. Linear Spline
    //alglib::spline1dbuildlinear(AXtc, AYtc, Xtc.size(), splineTc);

    //2. Cubic spline (natural- second derivatives in bounds = 0)
    //alglib::spline1dbuildcubic(AXtc, AYtc, Xtc.size(), 2, 0.0, 2, 0.0, splineTc);

    //3. Hermite spline (includes first derivative)
    std::vector<double> Dtc(999);

    for (int it = 0; it < 999; it++) {
        Dtc[it] = taxesmatrix[it + 1][7];
    }

    alglib::real_1d_array ADtc;
    ADtc.setcontent(Dtc.size(), &(Dtc[0]));

    alglib::spline1dbuildhermite(AXtc, AYtc, ADtc, Xtc.size(), splineTc);

    if (1 == 2) {
        cout << "-------Spline for Tc------" << endl;

        for (size_t i = 0; i < 20; i++) {
            printf("%f %f\n", Xtc[i], Ytc[i]);
        }
        
            printf("\n");
        for (int i = 500; i < 520; i++) {
            double x = i;
            printf("%f %f\n", x, alglib::spline1dcalc(splineTc, x));
        }
        printf("\n");

    }

    cout << "-------End Spline for Tc------" << endl;

    return(splineTc);

}

//A.2.2 Tc Spline
long double TcSpline(alglib::spline1dinterpolant splineTc, double pretax) {
    double myspline = alglib::spline1dcalc(splineTc, pretax);
    return(myspline);
}

//A.3.1. Set the spline for Tl
alglib::spline1dinterpolant setSplineTl(vector<vector<double>> taxesmatrix) {

    std::vector<double> Xtl(499), Ytl(499);

    for (int it = 0; it < 499; it++) {
        Xtl[it] = taxesmatrix[it + 1][8];
        Ytl[it] = taxesmatrix[it + 1][9];
    }

    alglib::real_1d_array AXtl, AYtl;
    AXtl.setcontent(Xtl.size(), &(Xtl[0]));
    AYtl.setcontent(Ytl.size(), &(Ytl[0]));

    alglib::spline1dinterpolant splineTl;

    //1. Linear Spline
    //alglib::spline1dbuildlinear(AXtl, AYtl, Xtl.size(), splineTl);

    //2. Cubic spline (natural- second derivatives in bounds = 0)
    //alglib::spline1dbuildcubic(AXtl, AYtl, Xtl.size(), 2, 0.0, 2, 0.0, splineTl);

    //3. Hermite spline (includes first derivative)
    std::vector<double> Dtl(499);

    for (int it = 0; it < 499; it++) {
        Dtl[it] = taxesmatrix[it + 1][10];
    }

    alglib::real_1d_array ADtl;
    ADtl.setcontent(Dtl.size(), &(Dtl[0]));

    alglib::spline1dbuildhermite(AXtl, AYtl, ADtl, Xtl.size(), splineTl);

    if (1 == 2) {
        cout << "-------Spline for Tl------" << endl;

        for (size_t i = 0; i < 20; i++) {
            printf("%f %f\n", Xtl[i], Ytl[i]);
        }
        
            printf("\n");
        for (int i = 500; i < 520; i++) {
            double x = i;
            printf("%f %f\n", x, alglib::spline1dcalc(splineTl, x));
        }
        printf("\n");

    }

    cout << "-------End Spline for Tl------" << endl;

    return(splineTl);

}

//A.2.2 Tl Spline
long double TlSpline(alglib::spline1dinterpolant splineTl, double income) {
    double myspline = alglib::spline1dcalc(splineTl, income);
    return(myspline);
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

//2. Production. C*0 to eliminate c
long double production(double ni, double nf, double aalpha, double tthetae,double c){
    double prod=tthetae*pow((c*1+ni+nf),aalpha);
    
    return(prod);
}

//3. Pre-tax profits.Checked. C*0 to eliminate c
long double profm(double ni, double nf, double aalpha, double tthetae, double wi, double wf, double c){
    double pi1=tthetae*pow((c*1+ni+nf),aalpha)-wi*ni-wf*nf-TnActual(nf*wf);
    return(pi1);
    
}

//3.1. Pre-tax profits.Checked. C*0 to eliminate c + Spline for Tn
long double profmSpline(double ni, double nf, double aalpha, double tthetae, double wi, double wf, double c, alglib::spline1dinterpolant splineTn) {
    long double taxn = TnSpline(splineTn, nf * wf);
    double pi1 = tthetae * pow((c * 1 + ni + nf), aalpha) - wi * ni - wf * nf - taxn;
    return(pi1);

}

//4. Corporate tax profits, marginal. C*0 to eliminate c
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
long double FinProfits(const double *Args, double paramvec[9], alglib::spline1dinterpolant splineTn, alglib::spline1dinterpolant splineTc){
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
    
    
    
    
    
    
    double ni=Args[0];
    double nf=Args[1];
    double z=Args[2];
    
    //Operational
    //double term1=profm(ni, nf, aalpha, tthetae, wi, wf,c);
    double term1 = profmSpline(ni, nf, aalpha, tthetae, wi, wf, c, splineTn);
    
    //Corporate taxes
    //double taxes=TcActual(z,ni,nf,aalpha,tthetae,wi,wf,c);
    long double taxes = TcSpline(splineTc, term1);
    

    //Cost of evasion
    
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

typedef struct {
    double* myParams;
    alglib::spline1dinterpolant mysplineTn;
    alglib::spline1dinterpolant mysplineTc;
}profitsFin_data;

double profitsFinMaxim(unsigned n, const double *x, double *grad, void *profitsFinMaxim_data){
    //Void to double
    //cout << profitsFinMaxim_data<< " this stuff"<< endl;
    
    //double *Params = (double *)profitsFinMaxim_data;
    
    profitsFin_data* d = (profitsFin_data*)profitsFinMaxim_data;

    double* Params = d->myParams;
    alglib::spline1dinterpolant splTc= d->mysplineTc;
    alglib::spline1dinterpolant splTn= d->mysplineTn;
    
    double result=FinProfits(x, Params, splTn, splTc);
    //cout << result << " evalF"<< endl;
    iterat=iterat+1;
    //printf("Iteration=(%d); Feval=%0.10g\n", iterat, result);
    return(result);
}

typedef struct {
    double* myParams;
    alglib::spline1dinterpolant mysplineTn;
}constraint_data;

double profitsConstraint(unsigned n, const double* x, double* grad, void* data)
{
    //double* Params = (double*)profitsFinMaxim_data;
    
    constraint_data* d = (constraint_data*) data;

    alglib::spline1dinterpolant splTn= d->mysplineTn;
    double* Params= d->myParams;
    
    double wi = Params[0];
    double wf = Params[1];
    double aalpha = Params[2];
    double ddelta = Params[3];
    double ggamma = Params[4];
    double bbeta = Params[5];
    double ssigma = Params[6];
    double tthetae = Params[7];
    double c = Params[8];

    //double cons = x[2]- tthetae*pow((c*1 + x[0] + x[1]), aalpha) + wi*x[0] + wf*x[1] + TnActual(wf*x[1]);
    
    //printf("Iteration=(%d); Feval=%0.10g\n", iterat, result);
    return(x[2] - tthetae * pow((c*1 + x[0] + x[1]), aalpha) + wi*x[0] + wf*x[1] + TnSpline(splTn, wf*x[1]));
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
                               double c,
                               alglib::spline1dinterpolant splineTn, 
                               alglib::spline1dinterpolant splineTc){
    
    
    //Modifying parameters to find internal solution to profits
    
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

    //Set the structures to pass in the optimization
    profitsFin_data datap = {Params,  splineTn, splineTc};
     
  
    constraint_data data = {Params, splineTn};

    
    //Defining object to be optimized
    //First. lower bounds for ni, nf, z.
    double lb[3];
    lb[0]=0;
    lb[1]=0;
    lb[2]=0;
         
   // double pre_profits= profm(double ni, double nf, double aalpha, double tthetae, double wi, double wf, double c)

     //Setting up optimization

    nlopt_opt optProfits;
    nlopt_opt localOpt;
    nlopt_opt subsOpt;

    //Define the optimization algorithm

        //1.Nelder-Mead (Local, non-gradient)
        //optProfits = nlopt_create(NLOPT_LN_NELDERMEAD, 3);
        //2.DIRECT -Dividing rectangles- (Global, non-gradient)
        //optProfits = nlopt_create(NLOPT_GN_ORIG_DIRECT, 3);
        //3.COBYLA (Local, non-gradient)
        //optProfits = nlopt_create(NLOPT_LN_COBYLA, 3);
        //4.ISRES (Global, non-gradient)
        //optProfits = nlopt_create(NLOPT_GN_ISRES, 3);

    //Do not support nonlinear constraints 
        //5. CRS (Global, non-gradient)
        //optProfits = nlopt_create(NLOPT_GN_CRS2_LM, 3);
        //6. ESCH (Global, non-gradient)
        //optProfits = nlopt_create(NLOPT_GN_ESCH, 3);
        //7. StoGO (Global, non-gradient)
        //optProfits = nlopt_create(NLOPT_GD_STOGO, 3);

        //8. MLSL -Multi Level Single Linkage- (Global, non-gradient)
        //optProfits = nlopt_create(NLOPT_G_MLSL_LDS, 3);
        //localOpt= nlopt_create(NLOPT_LN_COBYLA, 3);
        //nlopt_set_local_optimizer(optProfits,localOpt);

    //Augmented Lagragian + Subsidiary algorithm
        //9. NLOPT_AUGLAG
        optProfits = nlopt_create(NLOPT_AUGLAG, 3);

        //Set the subsidiary algorithm (local or global)
        subsOpt= nlopt_create(NLOPT_LN_NELDERMEAD, 3);
        //subsOpt= nlopt_create(NLOPT_GN_CRS2_LM, 3);
        //subsOpt= nlopt_create(NLOPT_GN_ESCH, 3);
        //subsOpt= nlopt_create(NLOPT_GD_STOGO, 3);

        nlopt_set_local_optimizer(optProfits,subsOpt);

        nlopt_set_lower_bounds(optProfits, lb);
        nlopt_set_max_objective(optProfits, profitsFinMaxim, (profitsFin_data*)&datap);

    // Upper bound for z (cannot exceed "pre-tax profits")
    nlopt_add_inequality_constraint(optProfits, profitsConstraint, (constraint_data*)&data, 1e-8);

    //Tolerance, abs or rel.
    //General 
    nlopt_set_xtol_rel(optProfits, 1.0e-8);
    nlopt_set_ftol_abs(optProfits, 1.0e-8);
    //Subsidiary algorithm
    nlopt_set_xtol_rel(subsOpt, 1.0e-8);
    nlopt_set_ftol_abs(subsOpt, 1.0e-8);

    //Number of evaluations
    //nlopt_set_population(optProfits, 1000000);
    //nlopt_set_maxeval(optProfits, 1.0e+18);

    double minf3; /* the minimum objective value, upon return */

    //Setting up initial conditions
    double x3[3];
    x3[0] = InitialCond[0];
    x3[1] = InitialCond[1];
    x3[2] = InitialCond[2];

    if(1==2){
        
        //Obtaining the analytical solution to the problem using a linear tax system:
        double d=0.2;
        
        //NI:
        double ninumerator1=(wf*(1+c)-wi*(1-d));
        double nifraction=ninumerator1/ddelta;
        double niexponent=1/ggamma;
        double niInitial=pow(nifraction,niexponent);
        
        
        //NF
        double nfnumerator1=(wf*(1+c)/(aalpha*tthetae));
        double nf1=pow(nfnumerator1,(1/(aalpha-1)));
        
        double nfnumerator2=(wf*(1+c)-wi)*(1-d);
        double nffraction2=nfnumerator2/ddelta;
        double nf2=pow(nffraction2,1/ggamma);
        double nfInitial=max(nf1-nf2,0.001);
        
        double zInitial=pow(d/bbeta,1/ssigma);
        
        x3[0]=niInitial;
        x3[1]=nfInitial;
        x3[2]=zInitial;
    }
    
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
    double Prof=FinProfits(Args, Params, splineTn, splineTc);
    
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
double ValueWorkers(const double *Args, double ParamWorkers[7], alglib::spline1dinterpolant splineTl){
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
    //double Taxes=PIT(tthetaw,wf,lf);
    double base=tthetaw*wf*lf;
    double Taxes=TlSpline(splineTl, base);
    
    
    double ans=term1-term2-term3-Taxes;
    return(ans);
}



//11. Maximize value of workers

typedef struct {
    double* myParamsW;
    alglib::spline1dinterpolant mysplineTl;
}workers_data;

int iteratvalworkers=0;
double valueWorkerMax(unsigned n, const double *x, double *grad, void *valueWorkerMax_data){
    //Void to double
    //cout << profitsFinMaxim_data<< " this stuff"<< endl;
    //double *Params = (double *)valueWorkerMax_data;

    workers_data* d = (workers_data*)valueWorkerMax_data;

    double* Params = d->myParamsW;
    alglib::spline1dinterpolant splTl= d->mysplineTl;

    
    double result=ValueWorkers(x, Params, splTl);
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
                                   double tthetaw,
                                   alglib::spline1dinterpolant splineTl){
    double Params[7]={};
    Params[0]=wf;
    Params[1]=wi;
    Params[2]=kkappa;
    Params[3]=rrho;
    Params[4]=psi;
    Params[5]=chi;
    Params[6]=tthetaw;
    
    
    
    //Set the structures to pass in the optimization
    workers_data dataw = { Params,  splineTl};
    
    //Defining the object to be maximized
    double lbWorker[2];
    lbWorker[0]=0;
    lbWorker[1]=0;
    
    
    
    nlopt_opt optWorkers;
    optWorkers = nlopt_create(NLOPT_LN_NELDERMEAD, 2);
    
    nlopt_set_lower_bounds(optWorkers, lbWorker);
    nlopt_set_max_objective(optWorkers, valueWorkerMax,(workers_data*)&dataw);
    
    
    nlopt_set_xtol_rel(optWorkers, 1.0e-8);
    nlopt_set_ftol_abs(optWorkers,1.0e-8);
    double MaxfWorker; /* the minimum objective value, upon return */
    double x3Worker[2];
    x3Worker[0]=InitialCond[0];
    x3Worker[1]=InitialCond[1];
    
    //Optimal analytical solutions
    double d=0.3;
    
    
    
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
                                  vector<double> InitProf, 
                                  alglib::spline1dinterpolant splineTn,
                                  alglib::spline1dinterpolant splineTc,
                                  alglib::spline1dinterpolant splineTl){
    
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
    
    
    
    
    ansValWorker=valueWorkerFinMaxim(InitialVworkers,wf,wi, kkappa,rrho,psi, chi, tthetaw, splineTl);
    
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
    ansProfits=profitsFinMaxim(InitialVProfits,  wi,wf, aalpha,ddelta,ggamma,bbeta,ssigma,tthetae,c, splineTn, splineTc);
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


int EquilibriumOP(int simdata,double vec1, double vec2, double vec3, double vec4, double vec5, double vec6, double vec7, double vec8, double vec9, double vec10, double vec11, double vec12, double vec13, double vec14, double vec15, arma::vec WagesInit,
        arma::vec armaInitLWorkers, arma::vec armaInitProf, vector<vector<double>> PLANNERINFO) {


        cout << armaInitProf[0] << "initprofdecission0 in distancenonvectorized" << endl;
        cout << armaInitProf[1] << "initprofdecission1 in distancenonvectorized" << endl;
        cout << armaInitProf[2] << "initprofdecission2 in distancenonvectorized" << endl;
        cout << armaInitLWorkers[0] << "InitLWorkersDecision in distancenonvectorized" << endl;
        cout << armaInitLWorkers[1] << "InitLWorkersDecision in distancenonvectorized" << endl;
        vector<double> VectorOthers(20);
        VectorOthers[0] = vec1;
        VectorOthers[1] = vec2;
        VectorOthers[2] = vec3;
        VectorOthers[3] = vec4;
        VectorOthers[4] = vec5;
        VectorOthers[5] = vec6;
        VectorOthers[6] = vec7;
        VectorOthers[7] = vec8;
        VectorOthers[8] = vec9;
        VectorOthers[9] = vec10;
        VectorOthers[10] = vec11;
        VectorOthers[11] = vec12;
        VectorOthers[12] = vec13;
        VectorOthers[13] = vec14;
        VectorOthers[14] = vec15;


        VectorOthers[15] = armaInitLWorkers[0];
        VectorOthers[16] = armaInitLWorkers[1];
        VectorOthers[17] = armaInitProf[0];
        VectorOthers[18] = armaInitProf[1];
        VectorOthers[19] = armaInitProf[2];


        cout << " about to find the wages and obtain the decision matrix estimator" << endl;

        //1. Set the equilibrium wages

        vector<double> WagesEquilibrium;
        WagesEquilibrium.resize(2);
        
        WagesEquilibrium[0]=0.504;
        WagesEquilibrium[1]=1.054;

        cout << WagesEquilibrium[0] << " wi " << endl;
        cout << WagesEquilibrium[1] << " wf " << endl;
        //2. Once we have the equilibrium wages, obtain the decision matrix.

        cout << " Obtain the decision matrix" << endl;
        //Initialize the variables to iterate over
        vector<double> Ttheta;
        Ttheta.resize(2);

        int Decision = 0;

        //Params for decission has the following structure:
        vector<double> ParamsDecision;
        ParamsDecision.resize(12);
        ParamsDecision[0] = WagesEquilibrium[0];
        ParamsDecision[1] = WagesEquilibrium[1];
        ParamsDecision[2] = VectorOthers[0];
        ParamsDecision[3] = VectorOthers[1];
        ParamsDecision[4] = VectorOthers[2];
        ParamsDecision[5] = VectorOthers[3];
        ParamsDecision[6] = VectorOthers[4];
        ParamsDecision[7] = VectorOthers[5];
        ParamsDecision[8] = VectorOthers[6];
        ParamsDecision[9] = VectorOthers[7];
        ParamsDecision[10] = VectorOthers[8];
        ParamsDecision[11] = VectorOthers[14];

        vector<double>InitLWorkers;
        InitLWorkers.resize(2);
        InitLWorkers[0] = armaInitLWorkers[0];
        InitLWorkers[1] = armaInitLWorkers[1];


        vector<double> InitProf;
        InitProf.resize(3);
        InitProf[0] = armaInitProf[0];
        InitProf[1] = armaInitProf[1];
        InitProf[2] = armaInitProf[2];
    
        vector<vector<double> > DecVector;
        DecVector.resize(4);
        for (int it = 0; it < 4; it++) {
            DecVector[it].resize(4);
        }

        vector<double> StorageAllocation;
        StorageAllocation.resize(15);

        //Initializing random number generator for the distribution
        int SEED = 2581633;
        typedef boost::mt19937 RNGType;
        RNGType rng(SEED);
        boost::random::uniform_real_distribution<> unidistrib(1, 2);
        boost::variate_generator<RNGType, boost::normal_distribution<> >
            generator(rng,
                boost::normal_distribution<>());

        double S0rr = gen_normal_3(generator);
        S0rr = gen_normal_3(generator);

        //We need to obtain the cholesky decomposition of the variance-covariance matrix
        int n = 2;


        //double m2[] = {0.01, 0.003,  0.003,  0.01};
        double mmu1 = VectorOthers[9];
        double mmu2 = VectorOthers[10];
        double ssigma1 = VectorOthers[11];
        double ssigma2 = VectorOthers[12];
        double rho12 = VectorOthers[13];


        double m2[] = { ssigma1, rho12,  rho12,  ssigma2 };
        double* c2 = cholesky(m2, n);

        int M = 500;

        //Draws from the normal distribution
        double z1 = 0;
        double z2 = 0;
        //Mean of log-normal thing

        //Writing the csv file of the allocations
        AllocationsOP.open("/home/ec2-user/InAws/AllocationsOP.csv", ios::out | ios::app);

        //Print the headline in the CSV file
        AllocationsOP << "thetaw" << " , " << "thetae" << " , " << "ni" << " , " << "nf" << " , " << "z" << " , " << "li" << " , " << "lf" << " , " << "pretax profits" << " , " << "check z" << " , " << "baseTn" << " , " << "Tn" << " , " << "BaseTc" << " , " << "Tc" << " , " << "BaseTl" << " , " << "Tl";

        AllocationsOP << endl;
        AllocationsOP.close();

        //Loading the Splines
        alglib::spline1dinterpolant ssplineTn = setSplineTn(PLANNERINFO);
        alglib::spline1dinterpolant ssplineTc = setSplineTc(PLANNERINFO);
        alglib::spline1dinterpolant ssplineTl = setSplineTl(PLANNERINFO);
        
        //AllocationsOP<< "Bla" << " , ";
        for (int it = 0; it < M; it++) {
        
          //cout <<"it: " << it+1 << endl;

            if (simdata == 1) {
                //Obtaining the draws from a standard normal distribution
                z1 = gen_normal_3(generator);
                z2 = gen_normal_3(generator);

                //Doing the corresponding transformation to obtain them from bivariate distribution
                Ttheta[0] = exp(c2[0] * z1 + c2[1] * z2 + mmu1);
                Ttheta[1] = exp(c2[2] * z1 + c2[3] * z2 + mmu2);
            }
            else {
                //cout << "ELSE " << endl;
                //cout <<"it: " << it+1 <<"Theta w " << PLANNERINFO[it+1][0] << endl;
                //cout << "it: " << it + 1 << "Theta e " << PLANNERINFO[it+1][1] << endl;
                Ttheta[0] = PLANNERINFO[it+1][0];
                Ttheta[1] = PLANNERINFO[it+1][1];
            }
            //cout << Ttheta[0] << " ttheta0"<< endl;
            //cout << Ttheta[1]<< " ttheta1"<< endl;

            //Ttheta[0]=1.2;
            //Ttheta[1]=5.5;
            //Ttheta[0]=unidistrib(rng);
            //Ttheta[1]=unidistrib(rng);

            //cout << Ttheta[0] << " Ttheta[0]"<< endl;
            //cout << Ttheta[1] << " Ttheta[1]"<< endl;


            //Computing the decissions
            DecVector = iDecision(Ttheta, ParamsDecision, InitLWorkers, InitProf, ssplineTn, ssplineTc,ssplineTl);

            //cout << " finished the decision matrix" << endl;

            //Load allocations into the storage vector (the one that will go to the CSV file)
            //(thetaw, thetae, ni, nf, z, li, lf, pretaxprofits, check, basetn, tn, basetc, tc, basetl, tl)
            StorageAllocation[0] = Ttheta[0];
            StorageAllocation[1] = Ttheta[1];
            StorageAllocation[2] = DecVector[0][1];
            StorageAllocation[3] = DecVector[1][1];
            StorageAllocation[4] = DecVector[2][1];
            StorageAllocation[5] = DecVector[0][2];
            StorageAllocation[6] = DecVector[1][2];
            StorageAllocation[7] = profmSpline(StorageAllocation[2], StorageAllocation[3], ParamsDecision[2], Ttheta[1], WagesEquilibrium[0], WagesEquilibrium[1], ParamsDecision[11], ssplineTn);
            StorageAllocation[8] = (StorageAllocation[4] < StorageAllocation[7]) ? 1 : 0;
            StorageAllocation[9] = StorageAllocation[3] * WagesEquilibrium[1];
            StorageAllocation[10] = TnSpline(ssplineTn ,StorageAllocation[9]);
            StorageAllocation[11] = StorageAllocation[7];
            StorageAllocation[12] = TcSpline(ssplineTc, StorageAllocation[11]);
            StorageAllocation[13] = StorageAllocation[0] * WagesEquilibrium[1] * StorageAllocation[6];
            StorageAllocation[14] = TlSpline(ssplineTl,StorageAllocation[13]);

            if (1==2){
                cout << StorageAllocation[0] << "thetaw" << endl;
                cout << StorageAllocation[1] << "thetae" << endl;
                cout << StorageAllocation[2] << "ni" << endl;
                cout << StorageAllocation[3] << "nf" << endl;
                cout << StorageAllocation[4] << "z" << endl;
                cout << StorageAllocation[5] << "li" << endl;
                cout << StorageAllocation[6] << "lf" << endl;
                cout << StorageAllocation[7] << "pretax" << endl;
                cout << StorageAllocation[8] << "check" << endl;
               }
        
            //Writing the csv file of the allocations
            AllocationsOP.open("/home/ec2-user/InAws/AllocationsOP.csv", ios::out | ios::app);

            for (int j = 0; j < 15; j++) {
                AllocationsOP<< StorageAllocation[j] << " , ";
                        }
            AllocationsOP << endl;
            AllocationsOP.close();
        }
     
    
        return(0);
}


void OptParamsTM(int simdata,arma::vec WagesVectorIn,
    arma::vec InitLWorkersDecision,
    arma::vec InitProfDecision, vector<vector<double>> PLANNERSIZE) {
    double a = 2;
    cout << a << " a " << endl;

    std::ifstream theFile("/home/ec2-user/InAws/OptParams15.csv");
    cout << " hello" << endl;
    int SIZEOBS = 1;
    int NVAR = 15;
    double MYARRAY[SIZEOBS][NVAR];
    // ...

    std::string line;
    std::vector<std::vector<std::string> > values;
    int it = 0;
    int it2 = 0;
    std::string line_value;
    std::vector<std::string> line_values;
    std::stringstream ss;
    while (std::getline(theFile, line))
    {

        ss << line;
        //std::stringstream ss(line);
        //std::string item;
        //cout <<  << "linevalprev"<<endl;
        while (std::getline(ss, line_value, ','))
        {
            line_values.push_back(line_value);
            MYARRAY[it][it2] = ::atof(line_value.c_str());
            //cout << MYARRAY[it][it2]<< " hahaha"<< endl;
            //MYARRAY[it][it2]=std::stod (line_value); only for c++11 compi
            it2 = it2 + 1;
            if (it2 == NVAR) { //later change 4 for
                it2 = 0;
            }
        }
        values.push_back(line_values);

        //For c++11 used values.emplace_back(line_values);

        //cout << line_value<< "line_value2"<< endl;
        it = it + 1;
        //Free the string types
        line_value.clear();
        line_values.clear();
        ss.clear();

    }






    vector<double> VectorOthers(21);
    VectorOthers[15] = InitLWorkersDecision[0];
    VectorOthers[16] = InitLWorkersDecision[1];
    VectorOthers[17] = InitProfDecision[0];
    VectorOthers[18] = InitProfDecision[1];
    VectorOthers[19] = InitProfDecision[2];

    cout << " 00000 in distance 00" << endl;
    cout << VectorOthers[14] << " VectorOthers[14]" << endl;
    cout << VectorOthers[15] << " VectorOthers[15]" << endl;
    cout << VectorOthers[16] << " VectorOthers[16]" << endl;
    cout << VectorOthers[17] << " VectorOthers[17]" << endl;

    cout << " loading the situation" << endl;

    double var1;
    double var2;
    double var3;
    double var4;
    double var5;
    double var6;
    double var7;
    double var8;
    double var9;
    double var10;
    double var11;
    double var12;
    double var13;
    double var14;
    double var15;

    cout << " Loading the parameters " << endl;
 
    var1 = MYARRAY[0][0];
    var2 = MYARRAY[0][1];
    var3 = MYARRAY[0][2];
    var4 = MYARRAY[0][3];
    var5 = MYARRAY[0][4];
    var6 = MYARRAY[0][5];
    var7 = MYARRAY[0][6];
    var8 = MYARRAY[0][7];
    var9 = MYARRAY[0][8];
    var10 = MYARRAY[0][9];
    var11 = MYARRAY[0][10];
    var12 = MYARRAY[0][11];
    var13 = MYARRAY[0][12];
    var14 = MYARRAY[0][13];
    var15 = MYARRAY[0][14];

    printf("insideloop %g %g  %g %g %g %g %g %g %g %g %g %g %g %g %g", var1, var2, var3, var4, var5, var6, var7, var8, var9, var10, var11, var12, var13, var14, var15);
    EquilibriumOP(simdata,var1, var2, var3, var4, var5, var6, var7, var8, var9, var10, var11, var12, var13, var14, var15, WagesVectorIn, InitLWorkersDecision, InitProfDecision,PLANNERSIZE);
}



int main(int argc, const char * argv[]) {
    
    
    
    
    
    
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
    
    
    //Second definition
    
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
    
    
    
    //Third attempt
    
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
    
    //Trying again
    aalpha=0.505664;
    ggamma=8.13167;
    ddelta=127.133;
    bbeta=23.168;
    ssigma=6.08206;
    kkappa=233.343;
    psi=0.366272;
    chi=178.087;
    rrho=8.98245;
    mmu1=0.323633;
    mmu2=3.12637;
    ssigma1=0.911768;
    ssigma2=4.5519;
    rho12=0.285791;
    c=1;
    
    //Initial small parameters
    aalpha=0.6124;
    ggamma=0.7341;
    ddelta=0.12873;
    bbeta=0.2315;
    ssigma=0.1827;
    kkappa=0.1021;
    psi=0.4528;
    chi=2.0192;
    rrho=0.0912;
    mmu1=1.2528;
    mmu2=1.7627;
    ssigma1=1.0921;
    ssigma2=1.1675;
    rho12=0.2782;
    c=0.01;
    
    
    //Initial small parameters
    aalpha=0.6124;
    ggamma=0.05341;
    ddelta=0.12873;
    bbeta=0.2315;
    ssigma=0.1827;
    kkappa=0.31021;
    psi=0.24528;
    chi=2.8192;
    rrho=0.10912;
    mmu1=0.8;
    mmu2=1.9;
    ssigma1=0.37921;
    ssigma2=0.9675;
    rho12=0.42;
    c=8;
    wi=2;
    wf=2;
    ni=10;
    nf=50;
    z=50;
    lf=50;
    li=50;
    
    //Defining tthetae and tthetaw to the mean:
    tthetaw=exp(mmu1+0.5*pow(ssigma1,2));
    tthetae=exp(mmu2+0.5*pow(ssigma2,2));
    
    cout << tthetaw << " mean tthetaw "<< endl;
    cout << tthetae << " mean tthetae "<< endl;
    
    //1. Payroll taxes definition
    double PayrolltestActual=TnActual(nf);
    cout << PayrolltestActual<< " Tn(nf)"<< endl;
    
    //2 prodution
    double productiontest=production( ni,  nf,  aalpha,  tthetae, c);
    cout << productiontest << " productiontest "<< endl;
    
    //3. Pre-tax profits
    cout << " ----testing profm ----"<<endl;
    double profmTest=profm( ni,  nf,  aalpha,  tthetae,  wi,  wf,c);
    cout << profmTest<< " profmTest"<< endl;
    
    //4. Tc Actual corporate taxes
    cout << "----- testing TcActualTest ---"<< endl;
    double TcActualTest=TcActual(z,ni,nf,aalpha,tthetae,wi,wf,c);
    cout << TcActualTest << " TcActualtest "<< endl;
    
    //5. PIT
    cout << " ---- testing PIT "<< endl;
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
    
    
    double Args[3]={};
    Args[0]=ni;
    Args[1]=nf;
    Args[2]=z;
    
    // cout << z << " z entering finprofits"<< endl;
    //double FinProfitsTest=FinProfits(Args, Params);
    //cout << FinProfitsTest << " FinProfitsTest "<< endl;
    
    
    //Value of workers
    double ParamWorkers[7]={};
    
    ParamWorkers[0]=wf;
    ParamWorkers[1]=wi;
    ParamWorkers[2]=kkappa;
    ParamWorkers[3]=rrho;
    ParamWorkers[4]=psi;
    ParamWorkers[5]=chi;
    ParamWorkers[6]=tthetaw;
    
    double ArgWorkers[2]={};
    ArgWorkers[0]=li;
    ArgWorkers[1]=lf;
    
 
    
    //-------------------------------------------------------------------
    cout << " Block of testing the Value of Workers and Profits to be balanced "<< endl;
    // Block PW.
    // Block of testing that VAlue workers and profits are balanced.
    //Testing maximize profits and value of workers. Trying in minimum, median and other value.
    
    
    //mean: tthetae=tthetae, tthetaw=tthetaw
    
    //We compute the variance of the distributions
    
    double VarMmu1=(exp(ssigma1*ssigma1)-1)*exp(2*mmu1+ssigma1*ssigma1);
    double VarMmu2=(exp(ssigma2*ssigma2)-1)*exp(2*mmu2+ssigma2*ssigma2);
    
    
    double tthetawmedian=exp(mmu1);
    double tthetaemedian=exp(mmu2);
    
    double tthetawmedianLow=tthetaw-2*VarMmu1;
    double tthetaemedianLow=tthetae-2*VarMmu2;
    
    
    double tthetawmedianHigh=tthetaw+2*VarMmu1;
    double tthetaemedianHigh=tthetae+2*VarMmu2;
    
    //Defining the corresponding ttheta's that we use to test. tthetae is already set for the mean.
    
    
    //Getting ready profits
    
    
    
    
    vector <double> MaxProf;
    MaxProf.resize(4);
    double InitialCond[3];
    InitialCond[0]=ni;
    InitialCond[1]=nf;
    InitialCond[2]=z;
    
    //Getting ready the value of workers.
    
    
    //Finding the actual optimal values found
    
    
  
    
    
    //Testing maximization of value of workers
    vector <double> ValAnsWorker;
    ValAnsWorker.resize(3);
    double InitialCondWorker[2];
    InitialCondWorker[0]=li;
    InitialCondWorker[1]=lf;
    

    

   //4. Testing at 5.757359322880715
    //cout << " -------------------------------------------------------------------- " << endl;
    //std::cout << " Testing at 5.757359322880715 " << std::endl;
    //cout << " Sin restriccion " << endl;
    //double zmax = profm(8.9796, 100.43, 0.73, 5.757359322880715, 0.504, 1.054, 10.0);
    //double IC[3];
    //IC[0] = ni;
    //IC[1] = nf;
    //IC[2] = zmax;

    //double t_wi = 0.504;
    //double t_wf = 1.054;
    //double t_alpha = 0.73;
    //double t_delta = 0.12873;
    //double t_gamma = 0.7341;
    //double t_beta = 0.2135;
    //double t_sigma = 0.1817;
    //double t_thetae = 5.757359322880715;
    //double t_c = 10.0;

    //auto start = high_resolution_clock::now();
    //MaxProf = profitsFinMaxim(IC, t_wi, t_wf, t_alpha, t_delta, t_gamma, t_beta, t_sigma, t_thetae, t_c);
    //auto stop = high_resolution_clock::now();
    
    //auto duration = duration_cast<microseconds>(stop - start);
    //    cout << "Time taken by profits optimization: "
    //    << duration.count() << " microseconds" << endl;
    //double zub = profm(MaxProf[0], MaxProf[1], t_alpha, t_thetae, t_wi, t_wf, t_c);

    
    //cout << " ni optimal (profits) " << MaxProf[0] << endl;
    //cout << " nf optimal (profits)" << MaxProf[1] << endl;
    //cout << " z optimal (profits)" << MaxProf[2] << endl;
    //cout << " z-max" << zub << endl;
    //cout << " z-max2" << 5.757359322880715 * pow((t_c*1 + MaxProf[0] + MaxProf[1]), 0.73) - 0.504*MaxProf[0] - 1.054*MaxProf[1] - TnActual(1.054* MaxProf[1]) << endl;
    //cout << " profits " << MaxProf[3] << endl;


    //ValAnsWorker = valueWorkerFinMaxim(InitialCondWorker, wf, wi, kkappa, rrho, psi, chi, 5.757359322880715);
    //cout << ValAnsWorker[0] << " li (workers)" << endl;
    //cout << ValAnsWorker[1] << " lf (workers) " << endl;
    //cout << ValAnsWorker[2] << " Value of workers" << endl;

    
    
    
    
    
    
    
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
 
    double WageEx[2];
    WageEx[0]=wi;
    WageEx[1]=wf;
    
    
    
    Wages[0]=WageEx[0];
    Wages[1]=WageEx[1];
    

    
    cout << " Reading inputs from the Planner's problem " << endl;
    int simdata = 0.0;
    vector<vector<double>> matplanner=Planner(simdata);
    cout << " end Reading inputs from the Planner's problem" << endl;

    alglib::spline1dinterpolant splineTn = setSplineTn(matplanner);
    alglib::spline1dinterpolant splineTc = setSplineTc(matplanner);
    alglib::spline1dinterpolant splineTl = setSplineTl(matplanner);

    long double trysplineTn = TnSpline(splineTn, 505);
    cout << " Tn spline in 505: " << trysplineTn <<endl;

    long double trysplineTc = TcSpline(splineTc, 5);
    cout << " Tc spline in 505: " << trysplineTc << endl;

    long double trysplineTl = TlSpline(splineTl, 5);
    cout << " Tl spline in 505: " << trysplineTl << endl;

    cout << " running Theoretical Moments/Allocations for OptParams" << endl;
    OptParamsTM(simdata,Wages,InitLWorkersDecision,InitProfDecision,matplanner);
    cout << " end Theoretical Moments/Allocations for OptParams"<< endl;
        
    //cout << " running sobol!"<< endl;
    //SobolRun(Wages,InitLWorkersDecision,InitProfDecision);
    //cout << " end of running sobol "<< endl;
    
    
    //-------------------------------------------------------
    //Block of initial prediction of moments
    //Obtaining initial theoretical moments with initial parameters specified
    
    //cout << " -----------Obtaining the initial predicted moments--------------------------------- "  << endl;
    //DistanceEstimator(VecOthers, WagesVectorIn,InitLWorkersDecision,InitProfDecision);
    //cout << " ------------- Obtaining the predicted moments done ---------------" << endl;
    
    
    
    //--------------------------------------------------------
    
    
    
    
    // ------------------------------------------
    
    // Minimization block
    // Doing the minimization of variables here:
    // Minimization block
    
    
    //cout << " Minimization block "<< endl;
    //cout << " Initializing distance minimization "<< endl;
    //MinimizingDistance(MinimizeDistanceParameters,AdditionalVAriables);
    //cout << " Distance minimization done"<< endl;
    // ------------------------------------------
    
    
    
    
    
    //EquilibriumMoments(VecOthers, WagesVectorIn,InitLWorkersDecision,InitProfDecision);
    //EqWagesNumericVector(VecOthers, WagesVectorIn);
    //StandardizedDistanceEstimator(ParStandardDistance,  AddPar);
    
    //StandardizedExcessDemands(WageEx,  Others);
    
    //ExcessDemandFunctions(Wages,ParamsDecisionExcessDemand, InitLWorkersDecision, InitProfDecision);
    
    //cout << " changonorrea"<< endl;
    //Decision=iDecision(TthetaDecision,ParamsDecision,InitLWorkersDecision,InitProfDecision);
    
    
    //----->   SobolRun(Wages,InitLWorkersDecision,InitProfDecision);
    
    
    //WagesEquilibrium=EqWagesNumericVector(VecOthers, WagesVectorIn);
    //
    
    
    
    
    
    
    //cout << MinimizeDistanceParameters[0] << " MinimizeDistanceParameters[0]"<< endl;
    
    
    
    
    
    //MinimizingDistance(MinimizeDistanceParameters,AdditionalVAriables);
    //cout << " running sobolrun "<< endl;
    //SobolRun(Wages,InitLWorkersDecision,InitProfDecision);
    
    //Decision=iDecision(TthetaDecision,ParamsDecision,InitLWorkersDecision,InitProfDecision);
    
    //cout << Decision[0][0]<< " Decision "<< endl;
    //cout << Decision[0][1]<< " ni "<< endl;
    //cout << Decision[1][1]<< " nf "<< endl;
    //cout << Decision[2][1]<< " z "<< endl;
    //cout << Decision[3][1]<< " prof "<< endl;
    //cout << Decision[0][2]<< " li "<< endl;
    //cout << Decision[1][2]<< " lf "<< endl;
    //cout << Decision[2][2]<< " VWorker "<< endl;
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    //cout << ExcessDemandFunctions(Wages,ParamsDecisionExcessDemand, InitLWorkersDecision, InitProfDecision)<<  endl;
    
    
    //cout << StandardizedExcessDemands(WageEx,  Others) << " Standardized1" << endl;
    //cout << ExcessDemandFunctions(Wages,ParamsDecisionExcessDemand, InitLWorkersDecision, InitProfDecision)<<  endl;
    //cout << ExcessDemandFunctions(Wages,ParamsDecisionExcessDemand, InitLWorkersDecision, InitProfDecision)<<  endl;
    //cout << ExcessDemandFunctions(Wages,ParamsDecisionExcessDemand, InitLWorkersDecision, InitProfDecision)<<  endl;
    
    //Minimizing excess demands squared
    //ExcessDemandsTotal(unsigned n, const double *x, double *grad, void *ExcessDemandsTotal_data)
    //double lbExcessDemmand[2];
    //lbExcessDemmand[0]=0;
    //lbExcessDemmand[1]=0;
    
    
    //nlopt_opt ExcessDemand;
    //ExcessDemand = nlopt_create(NLOPT_LN_NELDERMEAD, 2);
    
    //nlopt_set_lower_bounds(ExcessDemand, lbExcessDemmand);
    //nlopt_set_min_objective(ExcessDemand, ExcessDemandsTotal,(void *)&Others);
    
    
    //nlopt_set_xtol_rel(ExcessDemand, 1.0e-8);
    //nlopt_set_ftol_abs(ExcessDemand,1.0e-8);
    //double ExcessDemandValueFinal; /* the minimum objective value, upon return */
    
    //nlopt_destroy(ExcessDemand);
    
    
    
    
    
    
    
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
    
    
    //cout << Winitial[0] << " wi eq "<< endl;
    //cout << Winitial[1] << " wf eq" << endl;
    
    
    //Trying the function of optimizers
    //vector<double> EqWagesVector;
    //EqWagesVector.resize(2);
    
    
    
    
    //Vector de ExcessDemandFunctions es diferente del vector de idecision.
    
    //TthetaDecision[0]=2.5;
    //TthetaDecision[1]=0.67;
    
    //cout << InitProfDecision[0] << "InitProfDecision[0]"<<endl;
    //cout << InitProfDecision[1] << "InitProfDecision[1]"<<endl;
    //cout << InitProfDecision[2] << "InitProfDecision[2]"<<endl;
    //InitProfDecision[2]=2;
    
    
    
    
    //cout << Decision[0][0]<< " Decision "<< endl;
    //cout << Decision[0][1]<< " ni "<< endl;
    //cout << Decision[1][1]<< " nf "<< endl;
    //cout << Decision[2][1]<< " z "<< endl;
    //cout << Decision[3][1]<< " prof "<< endl;
    //cout << Decision[0][2]<< " li "<< endl;
    //cout << Decision[1][2]<< " lf "<< endl;
    //cout << Decision[2][2]<< " VWorker "<< endl;
    
    
    
    //cout << " obtaining equilibrium wages"<< endl;
    
    //WagesEquilibrium=EqWagesNumericVector(VecOthers, WagesVectorIn);
    
    //cout << " end of obtaining equilibrium wages"<< endl;
    //Testing the moment generating function
    //cout << "Testing theoretical moments "<< endl;
    //TheoMoments(ParamsDecisionExcessDemand, WagesVectorIn,InitLWorkersDecision,InitProfDecision);
    
    //cout << " end of theoretical moments test"<< endl;
    
    
    //Testing the distance estimators
    //DistanceEstimator(VecOthers, WagesVectorIn,InitLWorkersDecision,InitProfDecision);
    
    //cout << " initializing equilibrium moments "<< endl;
    //EquilibriumMoments(VecOthers, WagesVectorIn,InitLWorkersDecision,InitProfDecision);
    //cout << " finalizing equilibrium moments "<< endl;
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    //cout << WagesEquilibrium[0] << " WagesEquilibrium[0]"<< endl;
    //cout << WagesEquilibrium[1] << " WagesEquilibrium[1]"<< endl;
    
    
    
    
    
    
    
    return 0;
}
