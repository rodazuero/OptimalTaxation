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

ofstream AllocationsOPW("/home/ec2-user/InAws/AllocationsOPW.csv");


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
long double production(double n, double aalpha, double tthetae,double c){
    double prod=tthetae*pow((c*1+n),aalpha);
    
    return(prod);
}

//3. Pre-tax profits.Checked. C*0 to eliminate c
long double profm(double n, double aalpha, double tthetae, double w, double c){
    double pi1=tthetae*pow((c*1+n),aalpha)-w*n-TnActual(n*w);
    return(pi1);
    
}

//3.1. Pre-tax profits.Checked. C*0 to eliminate c + Spline for Tn
long double profmSpline(double n, double aalpha, double tthetae, double w, double c, alglib::spline1dinterpolant splineTn) {
    long double taxn = TnSpline(splineTn, n * w);
    double pi1 = tthetae * pow((c * 1 + n), aalpha) - w * n - taxn;
    return(pi1);

}

//4. Corporate tax profits, marginal. C*0 to eliminate c
long double Tc(double z, double n, double aalpha,double tthetae, double w,double c){
    double profm=tthetae*pow((c*1+n),aalpha)-w*n-TnActual(n*w);
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
long double TcActual(double z, double n, double aalpha, double tthetae, double w,double c){
    
    
    
    double profm=tthetae*pow((c*1+n),aalpha)-w*n-TnActual(n*w);
    double arg=profm-z;
    double ans=0;
    double tax=0;
    double prod=tthetae*pow((c*1+n),aalpha);
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
long double PIT(double tthetaw, double w, double l){
    double x=tthetaw*w*l;
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
long double PITM(double tthetaw, double w, double l){
    return(0.5);
}

//7. Final profits
long double FinProfits(const double *Args, double paramvec[8], alglib::spline1dinterpolant splineTn, alglib::spline1dinterpolant splineTc){
    //vector<double> paramvec=*(vector<double>*)params;
    //Before, it was (const double *Args, void *params){
    
    double w=paramvec[0];
    double aalpha=paramvec[1];
    double ddelta=paramvec[2];
    double ggamma=paramvec[3];
    double bbeta=paramvec[4];
    double ssigma=paramvec[5];
    double tthetae=paramvec[6];
    double c=paramvec[7];//added
    
    
    // cout << ggamma << " ggamma "<< endl;
    // cout << ddelta << " ddelta "<< endl;
    
    
    
    
    
    
    double n=Args[0];
    double z=Args[1];
    
    //Operational
    //double term1=profm(n, aalpha, tthetae, w,c);
    double term1 = profmSpline(n, aalpha, tthetae, w, c, splineTn);
    double Tcbase = term1-z;

    //Corporate taxes
    //double taxes=TcActual(z,n,aalpha,tthetae,w,c);
    long double taxes = TcSpline(splineTc, Tcbase);
    

    //Cost of evasion
    
    double evcost=pow(z,1+ssigma)*(bbeta/(1+ssigma));
    
    //Evasion costs
    double infcost=0.0;
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
    
    double w = Params[0];
    double aalpha = Params[1];
    double ddelta = Params[2];
    double ggamma = Params[3];
    double bbeta = Params[4];
    double ssigma = Params[5];
    double tthetae = Params[6];
    double c = Params[7];

    //printf("Iteration=(%d); Feval=%0.10g\n", iterat, result);
    return(x[1] - tthetae * pow((c*1 + x[0]), aalpha) + w*x[0] + TnSpline(splTn, w*x[0]));
}


//9. Defininf a function such that for given parameters, returns the optimal inputs and the maximum profits
vector<double> profitsFinMaxim(double InitialCond[2],
                               double w,
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
    
    double Params[8]={};
    
    Params[0]=w;
    Params[1]=aalpha;
    Params[2]=ddelta;
    Params[3]=ggamma;
    Params[4]=bbeta;
    Params[5]=ssigma;
    Params[6]=tthetae;
    Params[7]=c;

    //Set the structures to pass in the optimization
    profitsFin_data datap = {Params,  splineTn, splineTc};
     
  
    constraint_data data = {Params, splineTn};

    
    //Defining object to be optimized
    //First. lower bounds for ni, nf, z.
    double lb[2];
    lb[0]=0.0;
    lb[1]=0,0;
         
   

     //Setting up optimization

    nlopt_opt optProfits;
    //nlopt_opt localOpt;
    //nlopt_opt subsOpt;

    //Define the optimization algorithm

        //1.Nelder-Mead (Local, non-gradient)
        optProfits = nlopt_create(NLOPT_LN_NELDERMEAD, 2);
        //2.DIRECT -Dividing rectangles- (Global, non-gradient)
        //optProfits = nlopt_create(NLOPT_GN_ORIG_DIRECT, 2);
        //3.COBYLA (Local, non-gradient)
        //optProfits = nlopt_create(NLOPT_LN_COBYLA, 2);
        //4.ISRES (Global, non-gradient)
        //optProfits = nlopt_create(NLOPT_GN_ISRES, 2);

    //Do not support nonlinear constraints 
        //5. CRS (Global, non-gradient)
        //optProfits = nlopt_create(NLOPT_GN_CRS2_LM, 2);
        //6. ESCH (Global, non-gradient)
        //optProfits = nlopt_create(NLOPT_GN_ESCH, 2);
        //7. StoGO (Global, non-gradient)
        //optProfits = nlopt_create(NLOPT_GD_STOGO, 2);

        //8. MLSL -Multi Level Single Linkage- (Global, non-gradient)
        //optProfits = nlopt_create(NLOPT_G_MLSL_LDS, 2);
        //localOpt= nlopt_create(NLOPT_LN_COBYLA, 2);
        //nlopt_set_local_optimizer(optProfits,localOpt);

    //Augmented Lagragian + Subsidiary algorithm
        //9. NLOPT_AUGLAG
        //optProfits = nlopt_create(NLOPT_AUGLAG, 2);

        //Set the subsidiary algorithm (local or global)
        //subsOpt= nlopt_create(NLOPT_LN_NELDERMEAD, 2);
        //subsOpt= nlopt_create(NLOPT_GN_CRS2_LM, 2);
        //subsOpt= nlopt_create(NLOPT_GN_ESCH, 2);
        //subsOpt= nlopt_create(NLOPT_GD_STOGO, 2);

        //nlopt_set_local_optimizer(optProfits,subsOpt);

        nlopt_set_lower_bounds(optProfits, lb);
        nlopt_set_max_objective(optProfits, profitsFinMaxim, (profitsFin_data*)&datap);

    // Upper bound for z (cannot exceed "pre-tax profits")
    //nlopt_add_inequality_constraint(optProfits, profitsConstraint, (constraint_data*)&data, 1e-8);

    //Tolerance, abs or rel.
    //General 
    nlopt_set_xtol_rel(optProfits, 1.0e-8);
    //nlopt_set_ftol_abs(optProfits, 1.0e-8);
    //Subsidiary algorithm
    //nlopt_set_xtol_rel(subsOpt, 1.0e-8);
    //nlopt_set_ftol_abs(subsOpt, 1.0e-8);

    //Number of evaluations
    //nlopt_set_population(optProfits, 1000000);
    //nlopt_set_maxeval(optProfits, 1.0e+18);

    double minf3; /* the minimum objective value, upon return */

    //Setting up initial conditions
    double x3[2];
    x3[0] = InitialCond[0];
    x3[1] = InitialCond[1];


    //Finding optimal
    nlopt_optimize(optProfits, x3, &minf3);

    
    //if (nlopt_optimize(optProfits, x3, &minf3) < 0) {
    //    printf("nlopt failed!\n");
    //}
    //else {
    //   printf("found minimum at f(%g) = %0.10g\n", x3[0], minf3);
    //}
    
    
    double n=x3[0];
    double z=x3[1];

    double Args[3];
        
    Args[0]=n;
    Args[1]=z;    
    
    
    //Findin the profits
    double Prof=FinProfits(Args, Params, splineTn, splineTc);
    
    vector<double> Ans;
    Ans[0]=n;
    Ans[1]=z;
    Ans[2]=Prof;
    nlopt_destroy(optProfits);
    
    return(Ans);
}





//10. Value of workers
double ValueWorkers(const double *Args, double ParamWorkers[6], alglib::spline1dinterpolant splineTl){
    double w=ParamWorkers[0];
    double kkappa=ParamWorkers[1];
    double rrho=ParamWorkers[2];
    double psi=ParamWorkers[3];
    double chi=ParamWorkers[4];
    double tthetaw=ParamWorkers[5];
    
    double l=Args[0];

 
    //Term1 is total income
    double term1=(w*l)*tthetaw;
    
    //Term2 total disutility from working
    double term2=pow(l,1+psi)*(chi/(1+psi));
    
    //Term 3 is consumption penalty from informality
    double term3=0.0;
    
    //Term4 is net transfers to.from gogernment
    //double Taxes=PIT(tthetaw,wf,lf);
    double base=tthetaw*w*l;
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



vector<double> valueWorkerFinMaxim(double InitialCond[1],
                                   double w,
                                   double kkappa,
                                   double rrho,
                                   double psi,
                                   double chi,
                                   double tthetaw,
                                   alglib::spline1dinterpolant splineTl){
    double Params[6]={};
    Params[0]=w;
    Params[1]=kkappa;
    Params[2]=rrho;
    Params[3]=psi;
    Params[4]=chi;
    Params[5]=tthetaw;
    
    
    
    //Set the structures to pass in the optimization
    workers_data dataw = { Params,  splineTl};
    
    //Defining the object to be maximized
    double lbWorker[1];
    lbWorker[0]=0;
    
    
    nlopt_opt optWorkers;
    optWorkers = nlopt_create(NLOPT_LN_NELDERMEAD, 1);
    
    nlopt_set_lower_bounds(optWorkers, lbWorker);
    nlopt_set_max_objective(optWorkers, valueWorkerMax,(workers_data*)&dataw);
    
    
    nlopt_set_xtol_rel(optWorkers, 1.0e-8);
    nlopt_set_ftol_abs(optWorkers,1.0e-8);
    double MaxfWorker; /* the minimum objective value, upon return */
    double x3Worker[1];
    x3Worker[0]=InitialCond[0];
    
      
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
    ans.resize(2);
    ans[0]=x3Worker[0];
        
    ans[1]=MaxfWorker;
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
    double w=Params[0];
    double aalpha=Params[1];
    double ggamma=Params[2];
    double ddelta=Params[3];
    double bbeta=Params[4];
    double ssigma=Params[5];
    double kkappa=Params[6];
    double rrho=Params[7];
    double psi=Params[8];
    double chi=Params[9];
    double c=Params[10];
    
    
    
    
    
    //Loading ttheta
    double tthetaw=Ttheta[0];
    double tthetae=Ttheta[1];

        
    //Computing the value of the workers
    
    vector<double> ansValWorker;
    ansValWorker.resize(2);
    double InitialVworkers[1];
    InitialVworkers[0]=InitLWorkers[0];
     
    
    
    
    if(1==1){
        cout << " inside of the idecision stuff "<< endl;
        cout <<  w<< " w "<< endl;
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
    }
    
        
    ansValWorker=valueWorkerFinMaxim(InitialVworkers,w, kkappa,rrho,psi, chi, tthetaw, splineTl);
    
    //Loading answers of value of workers
    double l=ansValWorker[0];
    double ValueofWorker=ansValWorker[1];
    
    if(1==1){
        cout << " Worker "<< endl;
        cout <<  l<< "l"<< endl;
    }
    

    
    //Computing the Profits answer
    vector<double> ansProfits;
    ansProfits.resize(3);
    double InitialVProfits[2];
    InitialVProfits[0]=InitProf[0];
    InitialVProfits[1]=InitProf[1];
    ansProfits=profitsFinMaxim(InitialVProfits,  w, aalpha,ddelta,ggamma,bbeta,ssigma,tthetae,c, splineTn, splineTc);
    //Loading the answer of vlaue of profits
    double n=ansProfits[0];
    double z=ansProfits[1];
    double prof=ansProfits[2];
    
    if(1==1){
        cout << " Entrepreneur "<< endl;
        cout << n<< "n"<< endl;
        cout << z<< "z"<< endl;
        cout << prof<< "prof"<< endl;
    }
    
    vector<vector<double> > ans;
    ans.resize(3);
    for (int it=0; it<3; it++){
        ans[it].resize(3);
    }
    //Loading to the matrix ans the solutions
    ans[0][0]=0;
    if(prof>ValueofWorker){
        ans[0][0]=1;
    }
    
    ans[0][1]=n;
    ans[1][1]=z;
    ans[2][1]=prof;
    
    ans[0][2]=l;
    ans[1][2]=ValueofWorker;

    if(1==1){
        cout << " Entrepreneur/Worker Matrix"<< endl;
        cout <<  ans[0][1]<< "n"<< endl;
        cout <<  ans[1][1]<< "z"<< endl;
        cout <<  ans[2][1]<< "prof"<< endl;
        cout <<  ans[0][2]<< "l"<< endl;
        cout <<  ans[1][2]<< "ValueWorker"<< endl;
    }


    return(ans);
}


int EquilibriumOP(int simdata,double vec1, double vec2, double vec3, double vec4, double vec5, double vec6, double vec7, double vec8, double vec9, double vec10, double vec11, double vec12, double vec13, double vec14, double vec15, arma::vec WagesInit,
        arma::vec armaInitLWorkers, arma::vec armaInitProf, vector<vector<double>> PLANNERINFO) {


        cout << armaInitProf[0] << "initprofdecission0 in distancenonvectorized" << endl;
        cout << armaInitProf[1] << "initprofdecission1 in distancenonvectorized" << endl;
        cout << armaInitLWorkers[0] << "InitLWorkersDecision in distancenonvectorized" << endl;
        vector<double> VectorOthers(18);
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
        VectorOthers[16] = armaInitProf[0];
        VectorOthers[17] = armaInitProf[1];
 

        cout << " about to find the wages and obtain the decision matrix estimator" << endl;

        //1. Set the equilibrium wages

        vector<double> WagesEquilibrium;
        WagesEquilibrium.resize(1);
        
        WagesEquilibrium[0]=0.504;

        cout << WagesEquilibrium[0] << " wf " << endl;
        
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
        ParamsDecision[1] = VectorOthers[0];
        ParamsDecision[2] = VectorOthers[1];
        ParamsDecision[3] = VectorOthers[2];
        ParamsDecision[4] = VectorOthers[3];
        ParamsDecision[5] = VectorOthers[4];
        ParamsDecision[6] = VectorOthers[5];
        ParamsDecision[7] = VectorOthers[6];
        ParamsDecision[8] = VectorOthers[7];
        ParamsDecision[9] = VectorOthers[8];
        ParamsDecision[10] = VectorOthers[14];

        vector<double>InitLWorkers;
        InitLWorkers.resize(1);
        InitLWorkers[0] = armaInitLWorkers[0];

        vector<double> InitProf;
        InitProf.resize(2);
        InitProf[0] = armaInitProf[0];
        InitProf[1] = armaInitProf[1];
        
        vector<vector<double> > DecVector;
        DecVector.resize(3);
        for (int it = 0; it < 3; it++) {
            DecVector[it].resize(3);
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
        AllocationsOPW.open("/home/ec2-user/InAws/AllocationsOPW.csv", ios::out | ios::app);

        //Print the headline in the CSV file
        AllocationsOPW << "thetaw" << " , " << "thetae" << " , " << "ni" << " , " << "nf" << " , " << "z" << " , " << "li" << " , " << "lf" << " , " << "pretax profits" << " , " << "check z" << " , " << "baseTn" << " , " << "Tn" << " , " << "BaseTc" << " , " << "Tc" << " , " << "BaseTl" << " , " << "Tl";

        AllocationsOPW << endl;
        AllocationsOPW.close();

        //Loading the Splines
        alglib::spline1dinterpolant ssplineTn = setSplineTn(PLANNERINFO);
        alglib::spline1dinterpolant ssplineTc = setSplineTc(PLANNERINFO);
        alglib::spline1dinterpolant ssplineTl = setSplineTl(PLANNERINFO);
        
        //AllocationsOPW<< "Bla" << " , ";
        for (int it = 0; it < M; it++) {
        
          cout <<"it: " << it+1 << endl;

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
            //(thetaw, thetae, n, 0.0, z, l, 0.0, pretaxprofits, check, basetn, tn, basetc, tc, basetl, tl)
            StorageAllocation[0] = Ttheta[0];
            StorageAllocation[1] = Ttheta[1];
            StorageAllocation[2] = DecVector[0][1];
            StorageAllocation[3] = 0.0;
            StorageAllocation[4] = DecVector[1][1];
            StorageAllocation[5] = DecVector[0][2];
            StorageAllocation[6] = 0.0;
            StorageAllocation[7] = profmSpline(StorageAllocation[2], ParamsDecision[2], Ttheta[1], WagesEquilibrium[0], ParamsDecision[11], ssplineTn);
            StorageAllocation[8] = (StorageAllocation[4] < StorageAllocation[7]) ? 1 : 0;
            StorageAllocation[9] = StorageAllocation[2] * WagesEquilibrium[0];
            StorageAllocation[10] = TnSpline(ssplineTn ,StorageAllocation[9]);
            StorageAllocation[11] = StorageAllocation[7]-StorageAllocation[4];
            StorageAllocation[12] = TcSpline(ssplineTc, StorageAllocation[11]);
            StorageAllocation[13] = StorageAllocation[0] * WagesEquilibrium[0] * StorageAllocation[5];
            StorageAllocation[14] = TlSpline(ssplineTl,StorageAllocation[13]);

            if (1==1){
                cout << StorageAllocation[0] << "thetaw" << endl;
                cout << StorageAllocation[1] << "thetae" << endl;
                cout << StorageAllocation[2] << "n" << endl;
                cout << StorageAllocation[4] << "z" << endl;
                cout << StorageAllocation[5] << "l" << endl;
                cout << StorageAllocation[7] << "pretax" << endl;
                cout << StorageAllocation[8] << "check" << endl;
               }
        
            //Writing the csv file of the allocations
            AllocationsOPW.open("/home/ec2-user/InAws/AllocationsOPW.csv", ios::out | ios::app);

            for (int j = 0; j < 15; j++) {
                AllocationsOPW<< StorageAllocation[j] << " , ";
                        }
            AllocationsOPW << endl;
            AllocationsOPW.close();
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
    double w=8.81;
    
    double n=50;
    double ggamma=0.28;
    double ddelta=0.12;
    double bbeta=0.15;
    double ssigma=0.2;
    double kkappa=0.1;
    double psi=0.4;
    double chi=2.5;
    double rrho=0.9;
    
    
    
    double l=0.5;
    double z=24;
    double mmu1=0.2;
    double mmu2=1.3;
    double ssigma1=0.1;
    double ssigma2=0.9;
    double rho12=0.08;
    double c=10;
    
        
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
    
    z=50;
    w=2;
    n=10;
    l=50;
    

    //Defining tthetae and tthetaw to the mean:
    tthetaw=exp(mmu1+0.5*pow(ssigma1,2));
    tthetae=exp(mmu2+0.5*pow(ssigma2,2));
    
    cout << tthetaw << " mean tthetaw "<< endl;
    cout << tthetae << " mean tthetae "<< endl;
    
    //1. Payroll taxes definition
    double PayrolltestActual=TnActual(n);
    cout << PayrolltestActual<< " Tn(n)"<< endl;
    
    //2 prodution
    double productiontest=production( n,  aalpha,  tthetae, c);
    cout << productiontest << " productiontest "<< endl;
    
    //3. Pre-tax profits
    cout << " ----testing profm ----"<<endl;
    double profmTest=profm( n,  aalpha,  tthetae, w,c);
    cout << profmTest<< " profmTest"<< endl;
    
    //4. Tc Actual corporate taxes
    cout << "----- testing TcActualTest ---"<< endl;
    double TcActualTest=TcActual(z,n,aalpha,tthetae,w,c);
    cout << TcActualTest << " TcActualtest "<< endl;
    
    //5. PIT
    cout << " ---- testing PIT "<< endl;
    double PITtest=PIT(tthetaw,w,l);
    cout << PITtest << " PITtest " << endl;
    
    //6. PIT Marginal
    double PITtestMarginal=PITM(tthetaw,w,l);
    cout << PITtestMarginal<< " PITtestMarginal "<< endl;
    
    //7. FinProfits
    double Params[8]={};
    
    Params[0]=w;
    Params[1]=aalpha;
    Params[2]=ddelta;
    Params[3]=ggamma;
    Params[4]=bbeta;
    Params[5]=ssigma;
    Params[6]=tthetae;
    Params[7]=c;
    
    
    double Args[2]={};
    Args[0]=n;
    Args[1]=z;
   
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

    cout << z << " z entering finprofits"<< endl;
    double FinProfitsTest=FinProfits(Args, Params, splineTn, splineTc);
    cout << FinProfitsTest << " FinProfitsTest "<< endl;
    
    
    //Value of workers
    double ParamWorkers[6]={};
    
    ParamWorkers[0]=w;
    ParamWorkers[1]=kkappa;
    ParamWorkers[2]=rrho;
    ParamWorkers[3]=psi;
    ParamWorkers[4]=chi;
    ParamWorkers[5]=tthetaw;

    
    double ArgWorkers[1]={};
    ArgWorkers[0]=l;
    
    //iDecision

     vector<vector<double> > Decision;
        Decision.resize(3);
        for (int it = 0; it < 3; it++) {
            Decision[it].resize(3);
        }

    vector <double> ThetaDec;
    ThetaDec.resize(2);
    ThetaDec[0]=5.75;
    ThetaDec[1]=9.98;

    vector <double> ParamDec;
    ParamDec.resize(11);
    
    ParamDec[0]=w;
    ParamDec[1]=aalpha;
    ParamDec[2]=ggamma;
    ParamDec[3]=ddelta;
    ParamDec[4]=bbeta;
    ParamDec[5]=ssigma;
    ParamDec[6]=kkappa;
    ParamDec[7]=rrho;
    ParamDec[8]=psi;
    ParamDec[9]=chi;
    ParamDec[10]=c;
    
    vector <double> ArgsDec;
    ArgsDec.resize(2);
    ArgsDec[0]=Args[0];
    ArgsDec[1]=Args[1];

    vector <double> ArgWorkersDec;
    ArgWorkersDec.resize(1);
    ArgWorkersDec[0]=ArgWorkers[0];

    cout << " iDecision"<< endl;
    Decision=iDecision(ThetaDec,ParamDec,ArgWorkersDec,ArgsDec, splineTn, splineTc, splineTl);
    
    Decision[0][1]=0.0001;
    
        cout << " iDecision Matrix"<< endl;
        cout <<  Decision[0][1]<< "n"<< endl;
        cout <<  Decision[1][1]<< "z"<< endl;
        cout <<  Decision[2][1]<< "prof"<< endl;
        cout <<  Decision[0][2]<< "l"<< endl;
        cout <<  Decision[1][2]<< "ValueWorker"<< endl;
    

    
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
    MaxProf.resize(3);
    double InitialCond[2];
    InitialCond[0]=n;
    InitialCond[1]=z;
    
    //Getting ready the value of workers.
        
    //Finding the actual optimal values found
        
    //Testing maximization of value of workers
    vector <double> ValAnsWorker;
    ValAnsWorker.resize(2);
    double InitialCondWorker[1];
    InitialCondWorker[0]=l;    

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
    InitLWorkersDecision.resize(1);
    InitLWorkersDecision[0]=l;
    vector<double> InitProfDecision;
    InitProfDecision.resize(2);
    InitProfDecision[0]=n;;
    InitProfDecision[1]=z;
    vector<double>TthetaDecision;
    TthetaDecision.resize(2);
    TthetaDecision[0]=tthetaw;
    TthetaDecision[1]=tthetae;
    
    
    vector<double> Wages;
    Wages.resize(1);
    Wages[0]=w;
 
    double WageEx[1];
    WageEx[0]=w; 
    
    Wages[0]=WageEx[0];
   

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
