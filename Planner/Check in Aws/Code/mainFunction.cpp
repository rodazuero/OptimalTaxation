//  OPTIMAL TAXATION AND INFORMALITY
//  Code to explore the optimization algorithms available in NLopt
//  mainFunction.cpp
//
//  Created by Rodrigo Azuero Melo on 5/11/18. Modified Jan-2021.
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

//Define the quadratic function

long double Quadratic(const double* Args, double paramvec[3]) {

    double a = paramvec[0];
    double b = paramvec[1];
    double c = paramvec[2];

    double myx = Args[0];

    //Without penalization (restriction)
    double my_function = a*pow(myx, 2) + b*myx + c;
    
    return(my_function);
}



// Minimize the quadratic function 

//Define the function to be maximized by object nlopt
int ite = 0;

double quadMin(unsigned n, const double *x, double *grad, void *cdata){
    
    //Void to double
    //cout << profitsFinMaxim_data<< " this stuff"<< endl;
    double *Params = (double *)cdata;
    
    double a = Params[0];
    double b = Params[1];
    double c = Params[2];

    if (grad) {
        grad[0] = a*2.0*x[0] + b;
    }
    
    double result = Quadratic(x, Params);
    ite= ite+1;
    
    return(result);
}

//Define the constraint
typedef struct {
    double a, b, c;
}my_params;

double constraint(unsigned n, const double* x, double* grad, void* cdata)
{
    //double* Params = (double*)profitsFinMaxim_data;
    my_params* d = (my_params*) cdata;
    
    double aa = d->a;
    double bb = d->b;
    double cc = d->c;
    
    if (grad) {
        grad[0] = 1.0;
    }

   //double rest = -1.5 - x[0];
    return(-0.5 - x[0]);
}

//Definine the Nlopt optimization function
vector<double> Minimization(double InitialCond, double a,double b,double c) {


    //Modify parameters to find internal solution 

    double Params[3] = {};

    Params[0] = a;
    Params[1] = b;
    Params[2] = c;
   
    //Set the structure to pass in the optimization
    my_params data = { Params[0], Params[1], Params[2]};

    cout << " PARa " << data.a << " PARb " << data.b << " PARc " << data.c << endl;

    //Define the object to be optimized
    //First. lower bounds for ni, nf, z.
    //double lb[1];
    //lb[0] = -5.0;
    
    // double pre_profits= profm(double ni, double nf, double aalpha, double tthetae, double wi, double wf, double c)

     //Set up optimization

    nlopt_opt optProfits;
    nlopt_opt localOpt;
    nlopt_opt subsOpt;

    //Define the optimization algorithm

        //1.Nelder-Mead (Local, non-gradient)
        //optProfits = nlopt_create(NLOPT_LN_NELDERMEAD, 1);
        //2.DIRECT -Dividing rectangles- (Global, non-gradient)
        //optProfits = nlopt_create(NLOPT_GN_ORIG_DIRECT, 1);
        //3.COBYLA (Local, non-gradient)
        //optProfits = nlopt_create(NLOPT_LN_COBYLA, 1);
        //4.ISRES (Global, non-gradient)
        optProfits = nlopt_create(NLOPT_GN_ISRES, 1);
        
    //Those who do not support nonlinear constraints 
        //5. CRS (Global, non-gradient)
        //optProfits = nlopt_create(NLOPT_GN_CRS2_LM, 1);
        //6. ESCH (Global, non-gradient)
        //optProfits = nlopt_create(NLOPT_GN_ESCH, 1);
        //7. StoGO (Global, non-gradient)
        //optProfits = nlopt_create(NLOPT_GD_STOGO, 1);

        //8. MLSL -Multi Level Single Linkage- (Global, non-gradient)
        //optProfits = nlopt_create(NLOPT_G_MLSL_LDS, 1);
        //localOpt= nlopt_create(NLOPT_LN_COBYLA, 1);
        //nlopt_set_local_optimizer(optProfits,localOpt);

    //Augmented Lagragian + Subsidiary algorithm
        //9. NLOPT_AUGLAG
        //optProfits = nlopt_create(NLOPT_AUGLAG, 1);

        //Set the subsidiary algorithm (local or global)
        //subsOpt= nlopt_create(NLOPT_LN_NELDERMEAD, 1);
        //subsOpt= nlopt_create(NLOPT_GN_CRS2_LM, 1);
        //subsOpt= nlopt_create(NLOPT_GN_ESCH, 1);
        //subsOpt= nlopt_create(NLOPT_GD_STOGO, 1);
                
        //nlopt_set_local_optimizer(optProfits,subsOpt);


   
    nlopt_set_min_objective(optProfits, quadMin, (void*)&Params);

    // Upper bound for z (cannot exceed "pre-tax profits")
    nlopt_add_inequality_constraint(optProfits, constraint, (my_params*)&data, 1e-8);

    //Tolerance, abs or rel.
    //General 
    nlopt_set_xtol_rel(optProfits, 1.0e-8);
    nlopt_set_ftol_abs(optProfits, 1.0e-8);
    //Subsidiary algotithm
    //nlopt_set_xtol_rel(subsOpt, 1.0e-8);
    //nlopt_set_ftol_abs(subsOpt, 1.0e-8);

    //Number of evaluations
    //nlopt_set_population(optProfits, 1000000);
    nlopt_set_maxeval(optProfits, 1.0e+18);

    //Set up initial conditions
    double x1[1];
    x1[0] = InitialCond;
   

    double minf3; /* the minimum objective value, upon return */
    //Find optimal
    nlopt_optimize(optProfits, x1, &minf3);

    printf("found minimum after %d evaluations\n", ite);
    //if (nlopt_optimize(optProfits, x3, &minf3) < 0) {
    //    printf("nlopt failed!\n");
    //}
    //else {
    //   printf("found minimum at f(%g) = %0.10g\n", x3[0], minf3);
    //}

    double control[1];
    control[0]= x1[0];
    
    //Find the profits
    long double MinFound = Quadratic(control, Params);

    vector<double> Ans;
    Ans.resize(2);
    Ans[0] = control[0];
    Ans[1] = MinFound;

    nlopt_destroy(optProfits);

    return(Ans);

}


int main(int argc, const char* argv[]) {



    //Testing Quadratic Function
    cout << " -------------------------------------------------------------------- " << endl;
    cout << " Testing Quadratic Function " << endl;


    const double IC = -300;


    //double Args[1] = {};
    //Args[0] = IC;

    double t_a = 1.0;
    double t_b = 2.0;
    double t_c = 1.0;

    auto start = high_resolution_clock::now();
    vector<double> MyMin = Minimization(IC, t_a, t_b, t_c);
    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken by profits optimization: "
        << duration.count() << " microseconds" << endl;

    cout << " x optimal: " << MyMin[0] << endl;
    cout << " Min value" << MyMin[1] << endl;

    std::vector<double> X(5), Y(5);
    X[0] = 0.1;
    X[1] = 0.4;
    X[2] = 1.2;
    X[3] = 1.8;
    X[4] = 2.0;
    Y[0] = 0.1;
    Y[1] = 0.7;
    Y[2] = 0.6;
    Y[3] = 1.1;
    Y[4] = 0.9;

    alglib::real_1d_array AX, AY;
    AX.setcontent(X.size(), &(X[0]));
    AY.setcontent(Y.size(), &(Y[0]));

    alglib::spline1dinterpolant spline;
    alglib::spline1dbuildcubic(AX, AY, X.size(), 2, 0.0, 2, 0.0, spline);
    //alglib::spline1dbuildcubic(AX, AY, spline);


    for (size_t i = 0; i < X.size(); i++) {
        printf("%f %f\n", X[i], Y[i]);
    }

    printf("\n");
    for (int i = -50; i < 250; i++) {
        double x = 0.01 * i;
        printf("%f %f\n", x, alglib::spline1dcalc(spline, x));
    }
    printf("\n");


    return 0;

}