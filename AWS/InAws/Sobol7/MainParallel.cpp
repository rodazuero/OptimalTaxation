//
//  main.cpp
//  BEHAVIORAL38
//
//  Created by Rodrigo Azuero on 2/24/16.
//  Copyright (c) 2016 Rodrigo Azuero Melo. All rights reserved.
//



#include <algorithm>
#include <boost/math/distributions/normal.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/random.hpp>
#include <boost/random/variate_generator.hpp>
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
#include <unistd.h>
		
using std::vector;
using namespace std;

// Functions: Start with `F_' and then underscore continues
// Data inputs: All start with `in'
// Greek letters: start with double letters. e.g. aalpha, bbeta and ggamma.
// Likelihood functions: After the usual `F_' for functions we use `like_' to denote them.

//The structure of the program is as follows: The first block will contain the predicting equations. This are endogenous variables
//determined within the model or data that depends on parameters of the model (Wages, Skills, Effort, Consumption). The second part
//will include the likelihood computation for each part (i.e. likelihood of observing such wages, likelihood of observing skills, not
//necessarily of consumption but indeed of effort levels given the observed level). The third part includes the computation of all the likelihood
//function

//This file computes the likelihood as specified in the folder "Likelihoodfunction"


//===================================================================================//
//-1. Before starting with the function definitions I will define the transformations
//that are necessary for each parameter. I will not export them to R and they will all
//start with FT to denote Function and Transformation.
//==================================================================================//

//This function generates a discrete distribution.
//It is used exclusively as an example for how to draw from a
//discrete distribution so that later we get draws for the
//bootstrap filter.

double DISCDIST(double a){
    double resto=a*a;
    //NOW ATEMPT WITH DOUBLE VECTOR
    std::vector<double> probs(5);
    probs[0]=1;
    probs[1]=1;
    probs[2]=1;
    probs[3]=1;
    probs[4]=10;
    boost::random::discrete_distribution<> dist(probs.begin(), probs.end());
    //boost::random::discrete_distribution<> dist(probs);
    boost::mt19937 eng(std::time(0));
    std::cout << dist(eng);
    return(resto);
}

//Define template for random number generator
template<class T>
double gen_normal_3(T &generator)
{
    return generator();
}



//-1.0 Theta0 transformation
//--------------------------
// [[Rcpp::export]]
long double FT_ttheta0(double ttheta0, double ttheta1){
    //If ttheta0 too large we need to fix it
    if (ttheta0>700){
        ttheta0=700;
    }
    double result=exp(ttheta0)/(1+exp(ttheta0)+exp(ttheta1));
    if (result>0.999){
        result=0.99;
    }
    if (result<1.0e-5){
        result=1.0e-5;
    }
    return(result);
}

//-1.1 Ttheta1 transformation
//-------------------------
// [[Rcpp::export]]
long double FT_ttheta1(double ttheta0, double ttheta1){
    if (ttheta1>700){
        ttheta1=700;
    }
    
    double result=exp(ttheta1)/(1+exp(ttheta0)+exp(ttheta1));
    if (result>0.999){
        result=0.99;
    }
    if (result<1.0e-5){
        result=1.0e-5;
    }
    
    return(result);
}

//-1.2 Pphi transformation
//Function checked!
//-------------------------
// [[Rcpp::export]]
long double FT_pphi(double pphi){
    double result=1-exp(-pphi);
    // Fix pphi if it is zero
    if (result==0){
        result=1.0e-20;
    }
    
    // Fix pphi if it is too large
    if (pphi>4){
        result=1-exp(-4);
    }
    
    //Fix pphi if it is too small
    else if (pphi<-2.5){
        result = 1-exp(2.5);
    }
    return(result);
}


//-1.3 Exp transformation
//-------------------------
// [[Rcpp::export]]
long double FT_exp(double logstd){
    double R=exp(logstd);
    if (R==0){
        R=1.0e-17;
    }
    //0.0.1 We also need to take into account the fact that R can be too large
    if (logstd>700){
        R=exp(700);
    }
    return(R);
}


//-1.4 Gf transformation
//-----------------------------------------------------------------
// [[Rcpp::export]]
long double FT_gf(double gf){
    double result=gf;
    //We need to arrange this in case ggammaf or ggammam are small
    if (gf>0.999){
        result=0.99;
    }
    if (gf<1.0e-5){
        result=1.0e-5;
    }
    return(result);
}
//-1.5// aalphas. We need constraints in aalphas. They should be positive and sum up to one
//moreover, we need aalpha4=aalpha40+aalpha41 to be positive because otherwise we might get
//into trouble when trying to find the adequate level of effort for single mothers. The normalization of making them add
//to one needs to be done in place as we need all three arguments to be normalized and divide them by the sum of
//everyhing.
//-1.5.0 aalpha40
//------------------------------------------------------
// [[Rcpp::export]]
long double FT_aalpha40(double aalpha40){
    double a4m=1/(1+exp(aalpha40));
    double result=a4m;
    if (a4m<1.0e-5){
        result=1.0e-5;
    }
    if (a4m>0.9999){
        result=0.999;
    }
    return(result);
}
//-1.5.1 aalpha41
//------------------------------------------------------
// [[Rcpp::export]]
long double FT_aalpha41(double aalpha40, double aalpha41){
    double a40=FT_aalpha40(aalpha40);
    double a41=1/(1+exp(aalpha41));
    double result=a41*a40;
    return(result);
}

//-1.5.2 aalpha1
//-----------------------------------------------------
// [[Rcpp::export]]
long double FT_aalpha1(double aalpha1m){
    double a1m=1/(1+exp(aalpha1m));
    double result=a1m;
    if (a1m<1.0e-5){
        result=1.0e-5;
    }
    if (a1m>0.9999){
        result=0.999;
    }
    return(result);
}


//-1.5.2 aalpha2
//---------------------------------------------------
// [[Rcpp::export]]
long double FT_aalpha2(double aalpha2m){
    double a2m=1/(1+exp(aalpha2m));
    double result=a2m;
    if (a2m<1.0e-5){
        result=1.0e-5;
    }
    if (a2m>0.9999){
        result=0.999;
    }
    return(result);
}

//-1.5. aalpha3
//---------------------------------------------------
// [[Rcpp::export]]
long double FT_aalpha3(double aalpha3m){
    double a3m=1/(1+exp(aalpha3m));
    double result=a3m;
    if (a3m<1.0e-5){
        result=1.0e-5;
    }
    if (a3m>0.9999){
        result=0.999;
    }
    return(result);
}
//-1.6. Bounding transformations. This will serve
//the purpose of putting close to zero
//parameters that are really small or putting
//close to one other ones. Such is the case of the aalpha parameters
//----------------------------------------------------------------------
// [[Rcpp::export]]
long double FT_approach(double aalpha){
    double result=aalpha;
    if (aalpha<1.0e-5){
        result=1.0e-5;
    }
    if (aalpha>0.9999){
        result=0.999;
    }
    return(result);
}

//===========================================
//Second transformations of the aalphas:
//==========================================
//

long double FT_aalpha_trans(double aalpha1, double aalpha2, double aalpha3, double aalpha40, double aalpha41){
    if (aalpha1>700){
        aalpha1=700;
    }
    double result=exp(aalpha1)/(exp(aalpha1)+exp(aalpha2)+exp(aalpha3)+exp(aalpha40)+exp(aalpha41));
    return(result);
}




//===================================================================================//
//0. Block of function definition
//===================================================================================//




//----------
//0.0. Wages
//------------
//===================================================================================//
// [[Rcpp::export]]
long double F_predwage(double bbeta0, double bbeta1, double bbeta2, double bbeta3,  double Schooling, double Age){
    //1. Defining the output
    double wage=bbeta0+bbeta1*Schooling+bbeta2*Age+bbeta3*Age*Age;
    //1.1 If wage>700 we need to impose the max numerically possible
    if (wage>700){
        wage=exp(700);
    }
    //1.2
    //If such is not the case, we can define it as exp of the predicted value
    else{
        wage=exp(wage);
    }
    if(1==2){
        cout << "--------"<< endl;
        cout << bbeta0 << " bbeta0 " << endl;
        cout << bbeta1 << " bbeta1 " << endl;
        cout << bbeta2 << " bbeta2 " << endl;
        cout << bbeta3 << " bbeta3 " << endl;
        cout << Schooling << " Schooling " << endl;
        cout << Age << " Age " << endl;
        cout << "--------"<< endl;
    }
    return(wage);
}





//---------------------------
//0.1 Production of Skills
//---------------------------
//===================================================================================//
// [[Rcpp::export]]
long double F_predskills(double ddelta0, double ddelta1, double ddelta2, double ddelta3,
                         double ddelta4,
                         double Age,  double ttheta0, double ttheta1, double ttheta2,
                         double pphi, double gf, double gm, double Fe, double Me,
                         double I, double S0,
                         double childcare, double PG, double Hmembers){
    
    //0 Compute the level of effort in the CES production of Effort
    //0.0 If effort levels are zero, replace them by something close to zero.
    if (Fe<=0){
        Fe=0.01;
    }
    if (Me<=0){
        Me=0.01;
    }
    
    
    //And compute the effort by the CES production function
    double e=gf*pow(Fe,pphi)+gm*pow(Me,pphi);
    e=pow(e,1/pphi);
    if (e<=0){
        e=1.0e-5;
    }
    
    double R=exp(ddelta0+ddelta1*Age+ddelta2*childcare+ddelta3*PG+ddelta4*Hmembers);
    
    
    //1 If investments are zero or negative, put them close to zero
    if (I<=0){
        I=0.01;
    }
    //2 The production of skills
    
    double skills=R*pow(S0,ttheta0)*pow(I,ttheta1)*pow(e,ttheta2);
    double isf=std::isfinite(skills);
    if (isf==0){
        skills=5.0e+10;
    }
    if(1==2){
        cout << ddelta0 << " ddelta0 " << endl;
        cout << ddelta1*Age << " ddelta1*Age " << endl;
        cout << ddelta2*childcare << " ddelta2*childcare  "<< endl;
        cout << ddelta3*PG << " ddelta3*PG "<< endl;
        cout << ddelta4*Hmembers << " ddelta4*Hmembers "<< endl;
        cout << R << " r in skills " << endl;
        cout << pow(S0,ttheta0) << " pow(S0,ttheta0) " << endl;
        cout << I << " Investment level " << endl;
        cout << pow(I,ttheta1) << " pow(I,ttheta1) "<< endl;
        cout << pow(e,ttheta2) << " pow(e,ttheta2)  " << endl;
    }
    return(skills);
}

//----------------------------------------------
//0.2 Effort
// This function is defined in 09-02-2015, pg7.
//------------------------------------------------

//We will have to define several intermediate functions


//---------------------------------------------------
//0.2.0: Tau_j. This is an intermediate function defined
//in the same page as mentioned before
//---------------------------------------------------
//===================================================================================//
//0.2.0.0: tau_f
//===================================================================================//
// [[Rcpp::export]]
long double F_tau_f(double ggammaf, double aalpha4, double mmu, double hf){
    double tau_f=aalpha4*(mmu)*(1+hf);
    if (tau_f==0){
        tau_f=1.0e-10;
    }
    tau_f=ggammaf/tau_f;
    return(tau_f);
}
//===================================================================================//
//0.2.0.1: tau_m
//===================================================================================//
// [[Rcpp::export]]
long double F_tau_m(double ggammam, double aalpha4, double mmu, double hm){
    double tau_m=aalpha4*(1-mmu)*(1+hm);
    if (tau_m==0){
        tau_m=1.0e-10;
    }
    tau_m=ggammam/tau_m;
    return(tau_m);
}

//---------------------------------------------------
//0.2.1. ces_tau as defined in page 7
//----------------------------------------------------

//===================================================================================//
//0.2.1.0 ces_tau_f
//===================================================================================//
// [[Rcpp::export]]
long double F_ces_tau_f(double tau_f, double tau_m, double ggammaf,  double pphi){
    //0. Define the denominator of the  exponent
    double denexp=1-pphi;
    double ggammam=1-ggammaf;
    //0.0 If the denominator is infinite we need to Fix
    if (pphi==1){
        denexp=1.0e-10;
    }
    //0.1 Defining the exponent
    double exponent=pphi/denexp;
    if (exponent==0){
        exponent=1.0e-5;
    }
    //1. Now define the denominator
    double denominator=ggammaf*pow(tau_f,exponent)+ggammam*pow(tau_m,exponent);
    if (denominator<=0){
        denominator=1.0e-5;
    }
    if (denominator>=1.5e+10){
        denominator=1.5e+10;
    }
    //2. Defining the numerator
    double numerator=pow(tau_f,exponent);
    if (numerator>=1.5e+10){
        numerator=1.5e+10;
    }
    //3. Final definition of the ces_tau_f Functions
    double ces_tau_f=numerator/denominator;
    return(ces_tau_f);
}
//===================================================================================//
//0.2.1.1 ces_tau_m
//===================================================================================//
// [[Rcpp::export]]
long double F_ces_tau_m(double tau_f, double tau_m, double ggammaf,  double pphi){
    //0. Define the denominator of the  exponent
    double denexp=1-pphi;
    double ggammam=1-ggammaf;
    //0.0 If the denominator is infinite we need to Fix
    if (pphi==1){
        denexp=1.0e-10;
    }
    //0.1 Defining the exponent
    double exponent=pphi/denexp;
    if (exponent==0){
        exponent=1.0e-5;
    }
    
    //1. Now define the denominator
    double denominator=ggammaf*pow(tau_f,exponent)+ggammam*pow(tau_m,exponent);
    if (denominator<=0){
        denominator=1.0e-5;
    }
    if (denominator>=1.5e+10){
        denominator=1.5e+10;
    }
    //2. Defining the numerator
    double numerator=pow(tau_m,exponent);
    if (numerator>=1.5e+10){
        numerator=1.5e+10;
    }
    //3. Final definition of the ces_tau_f Functions
    double ces_tau_m=numerator/denominator;
    return(ces_tau_m);
}


//---------------------------------------------------
//0.2.2. kkappa2 weighted utility of skills
//----------------------------------------------------
//===================================================================================//
// [[Rcpp::export]]
long double F_kkappa2(double mmu, double aalpha2f, double aalpha2m){
    double kkappa2=aalpha2f*mmu+aalpha2m*(1-mmu);
    return(kkappa2);
}

//--------------------------------------------
//0.2.3 Effort levels
//--------------------------------------------
//===================================================================================//
//0.2.3.0 Effort of father
//===================================================================================//
// [[Rcpp::export]]
long double F_effort_f(double mmu, double ggammaf,double aalpha2f, double aalpha2m, double aalpha4f,double aalpha4m,double hf, double hm,double pphi,  double ttheta2){
    //0. Defining the tau inputs
    double ggammam=1-ggammaf;
    double tau_f=F_tau_f(ggammaf, aalpha4f, mmu, hf);
    double tau_m=F_tau_m(ggammam, aalpha4m, mmu, hm);
    //1. Defining the ces inputs
    double ces_tau_m=F_ces_tau_m( tau_f,  tau_m,  ggammaf,  pphi);
    
    //2. Defining kkappa2
    double kkappa2=F_kkappa2(mmu, aalpha2f,aalpha2m);
    
    //3. Finally we have the total effort
    double effort_f=tau_f*kkappa2*ttheta2*ces_tau_m;
    if (effort_f>5.0e+10){
        effort_f=5.0e+10;
    }
    if (effort_f<=0){
        effort_f=1.0e-10;
    }
    return(effort_f);
}
//===================================================================================//
//0.2.3.1 Effort of mother
//===================================================================================//
// [[Rcpp::export]]
long double F_effort_m(double mmu, double ggammaf,double aalpha2f, double aalpha2m, double aalpha4f,double aalpha4m,double hf, double hm,double pphi,  double ttheta2){
    //0. Defining the tau inputs
    
    double ggammam=1-ggammaf;
    double tau_f=F_tau_f(ggammaf, aalpha4f, mmu, hf);
    double tau_m=F_tau_m(ggammam, aalpha4m, mmu, hm);
    //1. Defining the ces inputs
    double ces_tau_m=F_ces_tau_f( tau_f,  tau_m,  ggammaf,  pphi);
    //2. Defining kkappa2
    double kkappa2=F_kkappa2(mmu, aalpha2f,aalpha2m);
    //3. Finally we have the total effort
    double effort_m=tau_m*kkappa2*ttheta2*ces_tau_m;
    
    if (effort_m>5.0e+10){
        effort_m=5.0e+10;
    }
    if (effort_m<=0){
        effort_m=1.0e-10;
    }
    return(effort_m);
    
    
}
//======================================
//Effort of father in the first period
//===================================================================================//
//0.2.3.0 Effort of father
//===================================================================================//
// [[Rcpp::export]]
long double F_effort_f10(double mmu, double ggammaf,double aalpha2f, double aalpha2m, double aalpha4f,double aalpha4m,double hf, double hm,double pphi,
                         double ttheta0,  double ttheta2, double BBETA){
    //0. Defining the tau inputs
    double ggammam=1-ggammaf;
    double tau_f=F_tau_f(ggammaf, aalpha4f, mmu, hf);
    double tau_m=F_tau_m(ggammam, aalpha4m, mmu, hm);
    //1. Defining the ces inputs
    double ces_tau_m=F_ces_tau_m( tau_f,  tau_m,  ggammaf,  pphi);
    
    //2. Defining kkappa2
    double kkappa2=F_kkappa2(mmu, aalpha2f,aalpha2m);
    
    //3. Finally we have the total effort
    double effort_f=tau_f*(kkappa2*ttheta2+BBETA*kkappa2*ttheta2*ttheta0)*ces_tau_m;
    if (effort_f>5.0e+10){
        effort_f=5.0e+10;
    }
    if (effort_f<=0){
        effort_f=1.0e-10;
    }
    return(effort_f);
}
//===================================================================================//
//0.2.3.1 Effort of mother
//===================================================================================//
// [[Rcpp::export]]
long double F_effort_m10(double mmu, double ggammaf,double aalpha2f, double aalpha2m, double aalpha4f,double aalpha4m,double hf, double hm,double pphi,
                         double ttheta0,  double ttheta2, double BBETA){
    //0. Defining the tau inputs
    
    double ggammam=1-ggammaf;
    double tau_f=F_tau_f(ggammaf, aalpha4f, mmu, hf);
    double tau_m=F_tau_m(ggammam, aalpha4m, mmu, hm);
    //1. Defining the ces inputs
    double ces_tau_m=F_ces_tau_f( tau_f,  tau_m,  ggammaf,  pphi);
    //2. Defining kkappa2
    double kkappa2=F_kkappa2(mmu, aalpha2f,aalpha2m);
    //3. Finally we have the total effort
    double effort_m=tau_m*((kkappa2*ttheta2+BBETA*kkappa2*ttheta2*ttheta0))*ces_tau_m;
    
    if (effort_m>5.0e+10){
        effort_m=5.0e+10;
    }
    if (effort_m<=0){
        effort_m=1.0e-10;
    }
    return(effort_m);
    
    
}


//===================================================================================//
//0.3 Utility Functions
//===================================================================================//
//---------------------------
// [[Rcpp::export]]
long double F_utility(double aalpha1j, double aalpha2j, double aalpha3j, double aalpha4j,  double cj, double ej, double hj, double S){
    if (S<=0){
        S=1.0e-5;
    }
    
    if (S>=1.0e+10){
        S=5.0e+10;
    }
    
    
    double utility=aalpha1j*log(cj+0.1)+aalpha2j*log(S)-aalpha3j*hj-aalpha4j*(1+hj)*ej;
    if (utility>5.0e+10){
        utility=5.0e+10;
    }
    
    if (utility<-5.0e+10){
        utility=-5.0e+10;
    }
    
    return(utility);
}

//Utility 10
// [[Rcpp::export]]
long double F_utility10(double aalpha1j, double aalpha2j, double aalpha3j, double aalpha4j,
                        double aalpha5j, double cj, double ej, double hj, double S,
                        double a){
    
    
    if (S<=0){
        S=1.0e-5;
    }
    
    if (S>=1.0e+10){
        S=5.0e+10;
    }
    
    double utility=aalpha1j*log(cj+0.1)+aalpha2j*log(S)-aalpha3j*hj-
    aalpha4j*(1+hj)*ej-aalpha5j*hj*(1-a);
    if (utility>5.0e+10){
        utility=5.0e+10;
    }
    if (utility<-5.0e+10){
        utility=-5.0e+10;
    }
    
    return(utility);
}

//===================================================================================//
//0.4 Consumption level
//===================================================================================//
// [[Rcpp::export]]
long double F_consumption(double aalpha1j, double aalpha2j, double ttheta1, double wj, double yj, double hj){
    double num=aalpha1j*(yj+wj*hj);
    //double isf=std::isfinite(num);
    if (num>=5.0e+10){
        num=5.0e+10;
    }
    if (num<=0){
        num=1.0e-10;
    }
    double den=aalpha2j*ttheta1+aalpha1j;
    if (den<=0){
        den=1.0e-5;
    }
    double consumption=num/den;
    if (consumption<=0){
        consumption=1.0e-5;
    }
    if (consumption>5.0e+10){
        consumption=5.0e+10;
    }
    
    return(consumption);
}

//===================================================================================//
//0.4.1 Consumption level IF JOINT Father
//===================================================================================//
// [[Rcpp::export]]
long double F_consumption_TOG(double aalpha1f, double aalpha2f, double aalpha2m,
                              double Investment,double ttheta1, double mmu,
                              double priceINV){
    
    double num=aalpha1f*mmu*Investment*priceINV;
    if (num<=0){
        num=1.0e-10;
    }
    if (num>=1.0e+10){
        num=1.0e+10;
    }
    double den=aalpha2f*mmu+aalpha2m*(1-mmu);
    den=den*ttheta1;
    if (den<=0){
        den=1.0e-5;
    }
    if (den>=1.0e+10){
        den=1.0e+10;
    }
    double consumption=num/den;
    
    if (consumption<=0){
        consumption=1.0e-5;
    }
    if (consumption>5.0e+10){
        consumption=5.0e+10;
    }
    return(consumption);
}

//===================================================================================//
//0.4.1 Consumption level IF JOINT Mother
//===================================================================================//
// [[Rcpp::export]]
long double M_consumption_TOG(double aalpha1m, double aalpha2f, double aalpha2m,
                              double Investment,double ttheta1, double mmu,
                              double priceINV){
    
    double num=aalpha1m*(1-mmu)*Investment*priceINV;
    if (num>=5.0e+10){
        num=5.0e+10;
    }
    if (num<=0){
        num=1.0e-10;
    }
    
    
    double den=aalpha2f*mmu+aalpha2m*(1-mmu);
    
    den=den*ttheta1;
    if (den<=0){
        den=1.0e-5;
    }
    if (den>=1.0e+10){
        den=1.0e+10;
    }
    
    double consumption=num/den;
    if (consumption<=0){
        consumption=1.0e-5;
    }
    if (consumption>5.0e+10){
        consumption=5.0e+10;
    }
    
    return(consumption);
}

//===================================================================================//
//0.4.1 Consumption level IF JOINT father 10
//===================================================================================//
// [[Rcpp::export]]
long double F_consumption_TOG10(double aalpha1f10, double aalpha2f, double aalpha2m,
                                double aalpha2f10, double aalpha2m10,
                                double Investment, double ttheta0,
                                double ttheta1,
                                double mmu, double priceINV){
    
    double num=aalpha1f10*(mmu)*Investment*priceINV;
    if (num>=5.0e+10){
        num=5.0e+10;
    }
    if (num<=0){
        num=1.0e-10;
    }
    
    
    double den=ttheta1*(aalpha2f10*mmu+aalpha2m10*(1-mmu))+
    0.92*ttheta0*ttheta1*(aalpha2f*mmu+aalpha2m*(1-mmu));
    
    if (den<=0){
        den=1.0e-5;
    }
    if (den>=1.0e+10){
        den=1.0e+10;
    }
    
    double consumption=num/den;
    if (consumption<=0){
        consumption=1.0e-5;
    }
    if (consumption>5.0e+10){
        consumption=5.0e+10;
    }
    
    return(consumption);
}

//===================================================================================//
//0.4.1 Consumption level IF JOINT Mother 10
//===================================================================================//
// [[Rcpp::export]]
long double M_consumption_TOG10(double aalpha1m10, double aalpha2f, double aalpha2m,
                                double aalpha2f10, double aalpha2m10,
                                double Investment, double ttheta0,
                                double ttheta1,
                                double mmu, double priceINV){
    
    double num=aalpha1m10*(1-mmu)*Investment*priceINV;
    if (num>=5.0e+10){
        num=5.0e+10;
    }
    if (num<=0){
        num=1.0e-10;
    }
    
    
    double den=ttheta1*(aalpha2f10*mmu+aalpha2m10*(1-mmu))+
    0.92*ttheta0*ttheta1*(aalpha2f*mmu+aalpha2m*(1-mmu));
    
    if (den<=0){
        den=1.0e-5;
    }
    if (den>=1.0e+10){
        den=1.0e+10;
    }
    
    double consumption=num/den;
    if (consumption<=0){
        consumption=1.0e-5;
    }
    if (consumption>5.0e+10){
        consumption=5.0e+10;
    }
    
    return(consumption);
}
//===================================================================================//
//0.5 Investment decision if single
//===================================================================================//
// [[Rcpp::export]]
long double F_investment(double aalpha1j, double aalpha2j, double ttheta1, double wj, double yj, double hj,double priceINV){
    double num=aalpha2j*ttheta1*(yj+wj*hj);
    
    if (num>=5.0e+10){
        num=5.0e+10;
    }
    if (num<=0){
        num=1.0e-10;
    }
    double den=priceINV*(aalpha1j+aalpha2j*ttheta1);
    if (den==0){
        den=1.0e-5;
    }
    
    double consumption=num/den;
    if (consumption<=0){
        consumption=1.0e-5;
    }
    if (consumption>5.0e+10){
        consumption=5.0e+10;
    }
    
    return(consumption);
}



//===================================================================================//
//0.5 Investment decision if single 2010
//===================================================================================//
// [[Rcpp::export]]
long double F_investment10(double aalpha1j10, double aalpha2j,
                           double aalpha2j10,
                           double ttheta1,
                           double ttheta0,
                           double wj, double yj, double hj,double priceINV){
    
    
    double num=(aalpha2j10*ttheta1+0.92*aalpha2j*ttheta0*ttheta1)*(yj+wj*hj);
    
    if (num>=5.0e+10){
        num=5.0e+10;
    }
    if (num<=0){
        num=1.0e-10;
    }
    double den=priceINV*(aalpha1j10+aalpha2j10*ttheta1+aalpha2j*ttheta0*ttheta1*0.92);
    if (den==0){
        den=1.0e-5;
    }
    
    double consumption=num/den;
    if (consumption<=0){
        consumption=1.0e-5;
    }
    if (consumption>5.0e+10){
        consumption=5.0e+10;
    }
    
    return(consumption);
}


//===================================================================================//
//0.5.1 Investment decision if joint
//===================================================================================//
//Verified
// [[Rcpp::export]]
long double F_invcouple(double aalpha1m, double aalpha1f, double aalpha2m, double aalpha2f,
                        double mmu, double hf,double hm,double wf,
                        double wm, double Yf, double Ym, double ttheta1, double priceINV){
    
    
    
    
    double TI=hf*wf+hm*wm+Yf+Ym;
    double k1=aalpha1f*mmu+aalpha1m*(1-mmu);
    double k2=aalpha2f*mmu+aalpha2m*(1-mmu);
    double denominator=(k1+k2*ttheta1)*priceINV;
    if (denominator<=0){
        denominator=1.0e-5;
    }
    if (denominator>=1.0e+10){
        denominator=1.0e+10;
    }
    double Investment=TI*k2*ttheta1/denominator;
    if (Investment>5.0e+10){
        Investment=5.0e+10;
    }
    
    if (Investment<=0){
        Investment=1.0e-5;
    }
    return(Investment);
}

//===================================================================================//
//0.5.1 Investment decision if joint 2010
//===================================================================================//
// [[Rcpp::export]]
long double F_invcouple10(double aalpha1f10, double aalpha1m10, double aalpha2f,
                          double aalpha2m,
                          double aalpha2f10,
                          double aalpha2m10,
                          double mmu,
                          double hf,
                          double hm,
                          double wf,
                          double wm,
                          double Yf, double Ym,
                          double ttheta0,
                          double ttheta1, double priceINV){
    
    double TI=hf*wf+hm*wm+Yf+Ym;
    double k1=aalpha1f10*mmu+aalpha1m10*(1-mmu);
    double k2=aalpha2f*mmu+aalpha2m*(1-mmu);
    double denominator=(k1+k2*ttheta1+
                        0.92*ttheta0*ttheta1*(aalpha2f10*mmu+(1-mmu)*aalpha2m10))*priceINV;
    if (denominator<=0){
        denominator=1.0e-5;
    }
    if (denominator>=1.0e+10){
        denominator=1.0e+10;
    }
    double Investment=TI*(k2*ttheta1+k2*ttheta0*ttheta1*0.92)/denominator;
    if (Investment>5.0e+10){
        Investment=5.0e+10;
    }
    
    if (Investment<=1.0e-5){
        Investment=1.0e-5;
    }
    
    return(Investment);
}



//================================
//0.6 Bargaining power
//================================
// [[Rcpp::export]]
long double F_mmu(double llambda0, double llambda1, double llambda2, double llambda3,
                  double llambda4, double llambda5, double llambda6, double llambda7,
                  double llambda8,
                  double wf, double wm,
                  double Yf, double Ym, double kkappa, double mmuLB,
                  double mmuUB, double Fage12, double Mage12,
                  double Fyrschool12, double Myrschool12,
                  double FMRATIO, double Unemployment, double Wageratio, double Distance){
    //0 Defining ratio of nonlabor income
    
    double ymratio=0;
    if (Yf+Ym==0){
        ymratio=0.5;
    }
    else {
        ymratio=Yf/(Yf+Ym);
    }
    if (wm==0){
        wm=1.0e-5;
    }
    double wmratio=wf/wm;
    //1 Defining what goes inside of the exponent
    double inexp=llambda0+llambda1*(wmratio)+llambda2*ymratio+
    llambda3*(Fage12-Mage12)+llambda4*(Fyrschool12-Myrschool12)+
    llambda5*FMRATIO+llambda6*Unemployment+llambda7*Wageratio+llambda8*Distance
    +kkappa;
    if (inexp>100){
        inexp=100;
    }
    double mmu=mmuLB+(mmuUB-mmuLB)*((exp(inexp))/(1+exp(inexp)));
    return(mmu);
}




//=================================================================================
//=================================================================================
//=================================================================================
//=================================================================================
//=================================================================================
//=================================================================================







//--------------------------------------------
//1. Block of Likelihood Functions definition
//--------------------------------------------
//=================================================================================
//1.0. Likelihood of Wages
//Checked conditional on F_predwage. Works fine
//=================================================================================
// [[Rcpp::export]]
long double F_likelihood_wage(double bbeta0, double bbeta1, double bbeta2, double bbeta3, double stdwage, double Schooling, double Age, double Wage){
    
    //1. Obtaining the predicted Wage:
    double predwage=F_predwage(bbeta0,bbeta1,bbeta2,bbeta3,Schooling,Age);
    //2. Once we obtain the predicted wage we proceed to
    //computing the likelihood Function
    //2.0. The first step is to define the normal distribution
    boost::math::normal_distribution<> normdist(0,1);
    
    //2.1 Now I will compute the actual likelihood of the observed Wage
    double likelihood_wage=log(Wage)-log(predwage);
    likelihood_wage=likelihood_wage/stdwage;
    likelihood_wage=pdf(normdist,likelihood_wage);
    
    likelihood_wage=likelihood_wage/stdwage;
    //2.2 If the pdf given is equal to zero we will approximate it to very close to zero
    if (likelihood_wage==0){
        likelihood_wage=1.0e-100;
    }
    likelihood_wage=log(likelihood_wage);
    return(likelihood_wage);
}



//=================================================================================
//1.1 Likelihood of skills --Checked conditional on F_predskills
//=================================================================================
// [[Rcpp::export]]
long double F_likelihood_skills(double ddelta0, double ddelta1,
                                double ddelta2, double ddelta3,
                                double ddelta4,
                                double ttheta0, double ttheta1,
                                double ttheta2, double pphi,
                                double gf, double gm,
                                double Fe, double Me,
                                double I, double S0,
                                double obS, double stdskills,
                                double Age, double CHILDCARE,
                                double PG,
                                double Hmembers){
    
    
    //1. Getting predicted level of skills
    double predskills=F_predskills(ddelta0, ddelta1, ddelta2, ddelta3, ddelta4,Age, ttheta0, ttheta1, ttheta2, pphi, gf, gm, Fe, Me, I, S0, CHILDCARE, PG,Hmembers);
    
    
    //1. Computing the likelihood of Skills
    boost::math::normal_distribution<> normdist (0,1);
    
    if (predskills==0){
        predskills=1.0e-5;
    }
    double loglik=log(obS)-log(predskills);
    
    loglik=loglik/stdskills;
    loglik=pdf(normdist,loglik);
    loglik=loglik/stdskills;
    if(loglik==0){
        loglik=-28;
    }
    
    else{
        loglik=log(loglik);
    }
    
    return(loglik);
}


//=================================================================================
//1.2 Likelihood of effort observed
//=================================================================================

//------------------------------------------
//1.2.0 Effort observed by father-not single -- Checked, conditional on F_effort_f
//------------------------------------------
// [[Rcpp::export]]
long double F_likelihood_feffort(double Feffort_observed, double mmu, double ggammaf,  double aalpha2f, double aalpha2m, double aalpha4m,  double aalpha4f,   double hf, double hm, double ttheta0, double ttheta1, double ttheta2,double pphi ,double stdeffort, double hchores){
    
    //0. Predicted level of effort of father
    //-----------------------------------------
    
    double effort_f=F_effort_f(mmu,ggammaf,aalpha2f,aalpha2m,aalpha4f,aalpha4m,hf,hm,pphi,ttheta2);
    
    //2. Computing the likelihood
    //-------------------------------
    
    //2.0 Obtaining the normal distribution
    boost::math::normal_distribution<> normdist(0,1);
    //2.1 And computing the likelihood
    double loglik_feffort=Feffort_observed-effort_f;
    
    loglik_feffort=loglik_feffort/stdeffort;
    loglik_feffort=pdf(normdist,loglik_feffort);
    loglik_feffort=loglik_feffort/stdeffort;
    //2.2 If the pdf given is zero, we will approximate it to the smallest possible value
    if(loglik_feffort==0){
        loglik_feffort=1.0e-100;
    }
    loglik_feffort=log(loglik_feffort);
    return(loglik_feffort);
}

//--------------------------------------------------------
// 1.2.1 Likelihood of effort observed by mother-not single -- checked conditional on F_effort_m
// -------------------------------------------------------
// [[Rcpp::export]]
long double F_likelihood_meffort(double Meffort_observed, double mmu, double ggammaf,  double aalpha2f, double aalpha2m,  double aalpha4m,   double aalpha4f,  double hf, double hm, double ttheta0, double ttheta1, double ttheta2, double pphi ,double stdeffort, double hchores){
    
    //0. Predicted level of effort of mother
    //-----------------------------------------
    double effort_m=F_effort_m(mmu,ggammaf,aalpha2f,aalpha2m,aalpha4f,aalpha4m,hf,hm,pphi,ttheta2);
    
    
    //2. Computing the likelihood
    //-------------------------------
    
    //2.0 Obtaining the normal distribution
    boost::math::normal_distribution<> normdist(0,1);
    //2.1 And computing the likelihood
    double loglik_meffort=Meffort_observed-effort_m;
    
    loglik_meffort=loglik_meffort/stdeffort;
    
    loglik_meffort=pdf(normdist,loglik_meffort);
    
    loglik_meffort=loglik_meffort/stdeffort;
    
    //2.2 If the pdf given is zero, we will approximate it to the smallest possible value
    if(loglik_meffort==0){
        loglik_meffort=1.0e-100;
    }
    
    loglik_meffort=log(loglik_meffort);
    return(loglik_meffort);
}



//------------------------------------------------------------
// 1.3 Likelihood of utility observed for discrete decisions
// we will do the four combinations: mother and father if married
// or single--checked conditional on intermediate functions
// ---------------------------------------------------------
// [[Rcpp::export]]
long double F_like_util_momsingle(double aalpha1m,  double aalpha3m, double aalpha4m,  double hm, double wm, double ym, double hchores, double Investment, double em, double S ,  double stdh1 ){
    //1. Defining the non-random components of the utilities
    //1.0 Consumption levels
    double cons_h0=ym+0.1-Investment;
    //1.0.0 If consumption level given by ym-I is smaller or equal to zero we need to fix it:
    if (cons_h0<=0){
        cons_h0=1.0e-5;
    }
    //1.0.1 Consumption level if works defined in reference to non-work Consumption
    double cons_h1=cons_h0+wm;
    
    
    
    
    //2. Defining the non-random components of the utilities:
    
    
    double UNR0=aalpha1m*log(cons_h0)-aalpha4m*em;
    double UNR1=aalpha1m*log(cons_h1)-aalpha3m-aalpha4m*em*2;
    //3.0. And now computing the likelihood. Defining numerator
    double nume=UNR1-UNR0;
    //3.1 denominator
    //double deno=pow(pow(stdh1,2)+1,0.5); //Remember one is normalized to one
    //before I put normalized one sd to be one. Now the problem is that it cannot be as small as we want
    //the whole variance which might be an issue. In this point I will normalize the sum of both
    double deno=stdh1;
    if (deno==0){
        deno=1.0e-10;
    }
    //1.3 Defining what goes inside of the cdf
    double incdf=nume/deno;
    //1.4 Computing the corresponding cdf
    boost::math::normal_distribution<> normdist (0 ,1);
    double cdfT=cdf(normdist,incdf);
    //If one or the other we need to take that into account
    if (hm==0){
        cdfT=1-cdfT;
    }
    if (cdfT==0){
        cdfT=1.0e-10;
    }
    cdfT=log(cdfT);
    return(cdfT);
}

//1.4 Likelihood of Investment levels without using predicted levels function
// [[Rcpp::export]]
long double F_likelihood_investment_single(double OBsinvestment,double priceINV, double aalpha1, double aalpha2, double ttheta0, double  Cfactorinv12, double  Mwage12, double  Mnly12, double  Mfraclabor12, double stdinvestment){
    double Totalincome=Mwage12*Mfraclabor12+Mnly12;
    double den=aalpha2*ttheta0+aalpha1;
    //0. If denominator close to zero re-arrange it:
    if (den==0){
        den=1.0e-5;
    }
    double num=aalpha2*ttheta0;
    double Investment=Totalincome*(num/den);
    Investment=Investment/priceINV;
    //1.0 Obtaining the normal distribution
    boost::math::normal_distribution<> normdist(0,1);
    //1.1 And computing the likelihood
    double loglik_investment=OBsinvestment-Investment;
    loglik_investment=loglik_investment/stdinvestment;
    loglik_investment=pdf(normdist,loglik_investment);
    loglik_investment=loglik_investment/stdinvestment;
    //2.2 If the pdf given is zero, we will approximate it to the smallest possible value
    if(loglik_investment==0){
        loglik_investment=1.0e-100;
    }
    loglik_investment=log(loglik_investment);
    return(loglik_investment);
    
}
//1.4.1 Likelihood of Investment levels  using predicted levels function
// [[Rcpp::export]]
long double F_likelihood_investment(double OBsinvestment,
                                    double Predinvestment,double MEASINV){
    //Getting a normal distribution
    boost::math::normal_distribution<> normdist (0,1);
    
    //Fixing not to take logs of zero
    if (OBsinvestment<=1.0e-5){
        OBsinvestment=1.0e-5;
    }
    if (Predinvestment<=1.0e-5){
        Predinvestment=1.0e-5;
    }
    
    //If variance is too small we need to fix it
    if (MEASINV<=1.0e-5){
        MEASINV=1.0e-5;
    }
    
    double loglikeinv=log(OBsinvestment)-log(Predinvestment);
    loglikeinv=loglikeinv/MEASINV;
    loglikeinv=pdf(normdist,loglikeinv);
    loglikeinv=loglikeinv/MEASINV;
    if (loglikeinv<=1.0e-5){
        loglikeinv=1.0e-5;
    }
    loglikeinv=log(loglikeinv);
    return(loglikeinv);
}



//1.5 Likelihood of mmu

// [[Rcpp::export]]
long double F_likelihood_mmu(double llambda0, double llambda1, double llambda2,
                             double llambda3, double llambda4,
                             double llambda5, double llambda6, double llambda7, double llambda8,
                             double Fwageinlik, double Mwageinlik,
                             double Mnly12, double Fnly12,
                             double kkappa, double mmuLB, double mmuUB,
                             double Fage12, double Mage12,
                             double Fyrschool12, double Myrschool12,
                             double FMRATIO, double Unemployment, double Wageratio, double Distance,
                             double Hbarg,
                             double MEASMMu){
    
    //Getting a normal distribution to compute the pdf
    boost::math::normal_distribution<> normdist (0,1);
    
    //First I will find the predicted level of mmu
    double loglikmmu=F_mmu(llambda0,llambda1,llambda2,llambda3,llambda4,
                           llambda5,llambda6,llambda7,llambda8,
                           Fwageinlik,Mwageinlik,
                           Mnly12,Fnly12,kkappa,mmuLB,mmuUB,Fage12,Mage12,
                           Fyrschool12,Myrschool12,FMRATIO,Unemployment,Wageratio,Distance);
    
    if(loglikmmu<=0){
        loglikmmu=1.0e-5;
    }
    loglikmmu=log(Hbarg)-log(loglikmmu);
    loglikmmu=loglikmmu/MEASMMu;
    loglikmmu=pdf(normdist,loglikmmu);
    loglikmmu=loglikmmu/MEASMMu;
    if (loglikmmu==0){
        loglikmmu=1.0e-10;
    }
    loglikmmu=log(loglikmmu);
    return(loglikmmu);
}


//1.6 I will define a generic loglikelihood for an observed variable
//And an error term that follows a log-normal distribution


// [[Rcpp::export]]
long double F_loglikelihood_generic(double observed, double predicted, double variance){
    
    //Fixing all the parameters
    if (observed<=1.0e-5){
        observed=1.0e-5;
    }
    if (predicted<=1.0e-5){
        predicted=1.0e-5;
    }
    if (variance<=1.0e-5){
        variance=1.0e-5;
    }
    //Getting a normal distribution to compute the pdf
    boost::math::normal_distribution<> normdist (0,1);
    //Getting the estimate
    double likelihood=log(observed)-log(predicted);
    likelihood=likelihood/variance;
    likelihood=pdf(normdist,likelihood);
    likelihood=likelihood/variance;
    if (likelihood<=1.0e-5){
        likelihood=1.0e-5;
    }
    return(likelihood);
}


//================================
//A. Factor likelihood systems
//================================


//Defining maximum between two numbers function

double maximum(double a, double b){
    if (a>=b){
        return a;
    }
    else{
        return b;
    }
}

//A.0 Ordered probit for n=2. I decided to treat all remaining variables for n larger than 2 as continuous.
double ORDPROB(double aalpha,double factor,int z, double sd){
    //Parameters: incdf: aalpha*ttheta, inside of the cdf
    //n: number of categories
    //z: value taken by the measure
    //sd: sd of the epsilon term
    
    //Initialize normal distr
    boost::math::normal_distribution<> normdist (0,1);
    
    
    //Initialize answer
    double ans=0;
    
    
    //Defining the threshold to compare with
    
    double threshold=-aalpha*factor/sd;
    
    
    
    
    
    //And defining the value of the normal at threshold cdf
    double norm=cdf(normdist,threshold);
    
    
    
    //If corner solution (z=0):
    if (z==0){
        ans=norm;
    }else if (z==1){
        ans=1-norm;
    }
    return(ans);
    
    
}

//A.1 Measurement likelihood system for Skills in 2010
double MEASloglikSkills2010(double Cmeasure1,
                            double aalpha1,
                            double Var1,
                            double Cmeasure2,
                            double aalpha2,
                            double Var2,
                            double Cmeasure3,
                            double aalpha3,
                            double Var3,
                            double Cmeasure4,
                            double aalpha4,
                            double Var4,
                            double Cmeasure5,
                            double aalpha5,
                            double Var5,
                            double Cmeasure6,
                            double aalpha6,
                            double Var6,
                            double Cmeasure7,
                            double aalpha7,
                            double Var7,
                            double Cmeasure8,
                            double aalpha8,
                            double Var8,
                            double Cmeasure9,
                            double aalpha9,
                            double Var9,
                            double Cmeasure10,
                            double aalpha10,
                            double Var10,
                            double Cmeasure11,
                            double aalpha11,
                            double Var11,
                            double logSKILL2010){
    
    if (logSKILL2010==0){
        logSKILL2010=1.0e-5;
    }
    logSKILL2010=log(logSKILL2010);
    //1. Compute the z score for each measurement
    double z1=max(abs((Cmeasure1-aalpha1*logSKILL2010)/Var1),1.0e-5);
    double z2=max(abs((Cmeasure2-aalpha2*logSKILL2010)/Var2),1.0e-5);
    double z3=max(abs((Cmeasure3-aalpha3*logSKILL2010)/Var3),1.0e-5);
    double z4=max(abs((Cmeasure4-aalpha4*logSKILL2010)/Var4),1.0e-5);
    double z5=max(abs((Cmeasure5-aalpha5*logSKILL2010)/Var5),1.0e-5);
    double z6=max(abs((Cmeasure6-aalpha6*logSKILL2010)/Var6),1.0e-5);
    double z7=max(abs((Cmeasure7-aalpha7*logSKILL2010)/Var7),1.0e-5);
    double z8=max(abs((Cmeasure8-aalpha8*logSKILL2010)/Var8),1.0e-5);
    double z9=max(abs((Cmeasure9-aalpha9*logSKILL2010)/Var9),1.0e-5);
    double z10=max(abs((Cmeasure10-aalpha10*logSKILL2010)/Var10),1.0e-5);
    double z11=max(abs((Cmeasure11-aalpha11*logSKILL2010)/Var11),1.0e-5);
    //2. //Getting a normal distribution to compute the pdf
    
    if(1==2){
        
        
        cout << Cmeasure1 << " measure1 1" << endl;
        cout << aalpha1*logSKILL2010 << " pred1 1" << endl;
        cout << Cmeasure2 << " measure1 2" << endl;
        cout << aalpha2*logSKILL2010 << " pred1 2" << endl;
        cout << Cmeasure3 << " measure1 3" << endl;
        cout << aalpha3*logSKILL2010 << " pred1 3" << endl;
        cout << Cmeasure4 << " measure1 4" << endl;
        cout << aalpha4*logSKILL2010 << " pred1 4" << endl;
        cout << Cmeasure5 << " measure1 5" << endl;
        cout << aalpha5*logSKILL2010 << " pred1 5" << endl;
        cout << Cmeasure6 << " measure1 6" << endl;
        cout << aalpha6*logSKILL2010 << " pred1 6" << endl;
        cout << Cmeasure7 << " measure1 7" << endl;
        cout << aalpha7*logSKILL2010 << " pred1 7" << endl;
        cout << Cmeasure8 << " measure1 8" << endl;
        cout << aalpha8*logSKILL2010 << " pred1 8" << endl;
        cout << Cmeasure9 << " measure1 9" << endl;
        cout << aalpha9*logSKILL2010 << " pred1 9" << endl;
        cout << Cmeasure10 << " measure10 1" << endl;
        cout << aalpha10*logSKILL2010 << " pred1 10" << endl;
        cout << Cmeasure11 << " measure11 1" << endl;
        cout << aalpha11*logSKILL2010 << " pred1 11" << endl;
        
        
        
        
        cout << Cmeasure1-aalpha1*logSKILL2010 << " error 1" << endl;
        cout << Cmeasure2-aalpha2*logSKILL2010 << " error 2" << endl;
        cout << Cmeasure3-aalpha3*logSKILL2010 << " error 3" << endl;
        cout << Cmeasure4-aalpha4*logSKILL2010 << " error 4" << endl;
        cout << Cmeasure5-aalpha5*logSKILL2010 << " error 5" << endl;
        cout << Cmeasure6-aalpha6*logSKILL2010 << " error 6" << endl;
        cout << Cmeasure7-aalpha7*logSKILL2010 << " error 7" << endl;
        cout << Cmeasure8-aalpha8*logSKILL2010 << " error 8" << endl;
        cout << Cmeasure9-aalpha9*logSKILL2010 << " error 9" << endl;
        cout << Cmeasure10-aalpha10*logSKILL2010 << " error 9" << endl;
        cout << Cmeasure11-aalpha11*logSKILL2010 << " error 10" << endl;
    }
    
    boost::math::normal_distribution<> normdist (0,1);
    
    //3. Compute the likelihood contribution of each measure
    
    
    z1=log(max(abs(pdf(normdist,z1)/Var1),1.0e-5));
    z2=log(max(abs(pdf(normdist,z2)/Var2),1.0e-5));
    z3=log(max(abs(pdf(normdist,z3)/Var3),1.0e-5));
    z4=log(max(abs(pdf(normdist,z4)/Var4),1.0e-5));
    z5=log(max(abs(pdf(normdist,z5)/Var5),1.0e-5));
    z6=log(max(abs(pdf(normdist,z6)/Var6),1.0e-5));
    z7=log(max(abs(pdf(normdist,z7)/Var7),1.0e-5));
    z8=log(max(abs(pdf(normdist,z8)/Var8),1.0e-5));
    z9=log(max(abs(pdf(normdist,z9)/Var9),1.0e-5));
    z10=log(max(abs(pdf(normdist,z10)/Var10),1.0e-5));
    z11=log(max(abs(pdf(normdist,z11)/Var11),1.0e-5));
    
    
    
    //4. Computing the final result
    //I multiply z4 times zero as it is tvip and we took it out.
    double result=z1+z2+z3+0*z4+z5+z6+z7+z8+z9+z10+z11;
    result=exp(result);
    return(result);
}


//--------------------
//A.2 Measurement likelihood system for Skills in 2012
//--------------------
double MEASloglikSkills2012(double Cmeasure1,
                            double aalpha1,
                            double Var1,
                            double Cmeasure2,
                            double aalpha2,
                            double Var2,
                            double Cmeasure3,
                            double aalpha3,
                            double Var3,
                            double Cmeasure4,
                            double aalpha4,
                            double Var4,
                            double Cmeasure5,
                            double aalpha5,
                            double Var5,
                            double Cmeasure6,
                            double aalpha6,
                            double Var6,
                            double Cmeasure7,
                            double aalpha7,
                            double Var7,
                            double Cmeasure8,
                            double aalpha8,
                            double Var8,
                            double Cmeasure9,
                            double aalpha9,
                            double Var9,
                            double Cmeasure10,
                            double aalpha10,
                            double Var10,
                            double Cmeasure11,
                            double aalpha11,
                            double Var11,
                            double Cmeasure12,
                            double aalpha12,
                            double Var12,
                            double Cmeasure13,
                            double aalpha13,
                            double Var13,
                            double logSKILL2012){
    
    
    
    
    if (logSKILL2012==0){
        logSKILL2012=1.0e-5;
    }
    
    
    logSKILL2012=log(logSKILL2012);
    //1. Compute the z score for each measurement
    double z1=max(abs((Cmeasure1-aalpha1*logSKILL2012)/Var1),1.0e-5);
    double z2=max(abs((Cmeasure2-aalpha2*logSKILL2012)/Var2),1.0e-5);
    double z3=max(abs((Cmeasure3-aalpha3*logSKILL2012)/Var3),1.0e-5);
    double z4=max(abs((Cmeasure4-aalpha4*logSKILL2012)/Var4),1.0e-5);
    double z5=max(abs((Cmeasure5-aalpha5*logSKILL2012)/Var5),1.0e-5);
    double z6=max(abs((Cmeasure6-aalpha6*logSKILL2012)/Var6),1.0e-5);
    double z7=max(abs((Cmeasure7-aalpha7*logSKILL2012)/Var7),1.0e-5);
    double z8=max(abs((Cmeasure8-aalpha8*logSKILL2012)/Var8),1.0e-5);
    double z9=max(abs((Cmeasure9-aalpha9*logSKILL2012)/Var9),1.0e-5);
    double z10=max(abs((Cmeasure10-aalpha10*logSKILL2012)/Var10),1.0e-5);
    double z11=max(abs((Cmeasure11-aalpha11*logSKILL2012)/Var11),1.0e-5);
    double z12=max(abs((Cmeasure12-aalpha12*logSKILL2012)/Var12),1.0e-5);
    double z13=max(abs((Cmeasure13-aalpha13*logSKILL2012)/Var13),1.0e-5);
    //2. //Getting a normal distribution to compute the pdf
    
    
    
    boost::math::normal_distribution<> normdist (0,1);
    
    //3. Compute the likelihood contribution of each measure
    
    
    z1=log(max(abs(pdf(normdist,z1)/Var1),1.0e-5));
    z2=log(max(abs(pdf(normdist,z2)/Var2),1.0e-5));
    z3=log(max(abs(pdf(normdist,z3)/Var3),1.0e-5));
    z4=log(max(abs(pdf(normdist,z4)/Var4),1.0e-5));
    z5=log(max(abs(pdf(normdist,z5)/Var5),1.0e-5));
    z6=log(max(abs(pdf(normdist,z6)/Var6),1.0e-5));
    z7=log(max(abs(pdf(normdist,z7)/Var7),1.0e-5));
    z8=log(max(abs(pdf(normdist,z8)/Var8),1.0e-5));
    z9=log(max(abs(pdf(normdist,z9)/Var9),1.0e-5));
    z10=log(max(abs(pdf(normdist,z10)/Var10),1.0e-5));
    z11=log(max(abs(pdf(normdist,z11)/Var11),1.0e-5));
    z12=log(max(abs(pdf(normdist,z12)/Var12),1.0e-5));
    z13=log(max(abs(pdf(normdist,z13)/Var13),1.0e-5));
    
    //4. Computing the final result
    double result=z1+z2+z3+z4+z5+z6+z7+z8+z9+z10+0*z11+0*z12+z13;
    result=exp(result);
    
    
    if(1==2){
        cout << result << " result " << endl;
    }
    return(result);
}




//--------------------
//A.3 Measurement likelihood system for Skills in 2010
//--------------------
double MEASloglikSkillsBIRTH(int Cmeasure1,
                             double aalpha1,
                             double Var1,
                             int Cmeasure2,
                             double aalpha2,
                             double Var2,
                             int Cmeasure3,
                             double aalpha3,
                             double Var3,
                             int Cmeasure4,
                             double aalpha4,
                             double Var4,
                             int Cmeasure5,
                             double aalpha5,
                             double Var5,
                             int Cmeasure6,
                             double aalpha6,
                             double Var6,
                             int Cmeasure7,
                             double aalpha7,
                             double Var7,
                             int Cmeasure8,
                             double aalpha8,
                             double Var8,
                             int Cmeasure9,
                             double aalpha9,
                             double Var9,
                             int Cmeasure10,
                             double aalpha10,
                             double Var10,
                             int Cmeasure11,
                             double aalpha11,
                             double Var11,
                             int Cmeasure12,
                             double aalpha12,
                             double Var12,
                             int Cmeasure13,
                             double aalpha13,
                             double Var13,
                             int Cmeasure14,
                             double aalpha14,
                             double Var14,
                             int Cmeasure15,
                             double aalpha15,
                             double Var15,
                             int Cmeasure16,
                             double aalpha16,
                             double Var16,
                             double Cmeasure17,
                             double aalpha17,
                             double Var17,
                             double Cmeasure18,
                             double aalpha18,
                             double Var18,
                             double Cmeasure19,
                             double aalpha19,
                             double Var19,
                             double Cmeasure20,
                             double aalpha20,
                             double Var20,
                             int Cmeasure21,
                             double aalpha21,
                             double Var21,
                             double Cmeasure22,
                             double aalpha22,
                             double Var22,
                             double Cmeasure23,
                             double aalpha23,
                             double Var23,
                             double logBIRTH){
    
    if (logBIRTH==0){
        logBIRTH=1.0e-5;
    }
    logBIRTH=log(logBIRTH);
    
    //1. Compute the z score for each measurement
    double z1=ORDPROB(aalpha1,logBIRTH,Cmeasure1,Var1);
    double z2=ORDPROB(aalpha2,logBIRTH,Cmeasure2,Var2);
    double z3=ORDPROB(aalpha3,logBIRTH,Cmeasure3,Var3);
    double z4=ORDPROB(aalpha4,logBIRTH,Cmeasure4,Var4);
    double z5=ORDPROB(aalpha5,logBIRTH,Cmeasure5,Var5);
    double z6=ORDPROB(aalpha6,logBIRTH,Cmeasure6,Var6);
    double z7=ORDPROB(aalpha7,logBIRTH,Cmeasure7,Var7);
    double z8=ORDPROB(aalpha8,logBIRTH,Cmeasure8,Var8);
    double z9=ORDPROB(aalpha9,logBIRTH,Cmeasure9,Var9);
    double z10=ORDPROB(aalpha10,logBIRTH,Cmeasure10,Var10);
    double z11=ORDPROB(aalpha11,logBIRTH,Cmeasure11,Var11);
    double z12=ORDPROB(aalpha12,logBIRTH,Cmeasure12,Var12);
    double z13=ORDPROB(aalpha13,logBIRTH,Cmeasure13,Var13);
    double z14=ORDPROB(aalpha14,logBIRTH,Cmeasure14,Var14);
    double z15=ORDPROB(aalpha15,logBIRTH,Cmeasure15,Var15);
    double z16=ORDPROB(aalpha16,logBIRTH,Cmeasure16,Var16);
    double z17=max(abs((Cmeasure17-aalpha17*logBIRTH)/Var17),1.0e-5);
    double z18=max(abs((Cmeasure18-aalpha18*logBIRTH)/Var18),1.0e-5);
    double z19=max(abs((Cmeasure19-aalpha19*logBIRTH)/Var19),1.0e-5);
    double z20=max(abs((Cmeasure20-aalpha20*logBIRTH)/Var20),1.0e-5);
    double z21=ORDPROB(aalpha21,logBIRTH,Cmeasure21,Var21);
    double z22=max(abs((Cmeasure22-aalpha22*logBIRTH)/Var22),1.0e-5);
    double z23=max(abs((Cmeasure23-aalpha23*logBIRTH)/Var23),1.0e-5);
    
    
    //2. //Getting a normal distribution to compute the pdf
    
    boost::math::normal_distribution<> normdist (0,1);
    
    
    //3. Compute the likelihood contribution of each measure
    z1=log(max(z1,1.0e-5));
    z2=log(max(z2,1.0e-5));
    z3=log(max(z3,1.0e-5));
    z4=log(max(z4,1.0e-5));
    z5=log(max(z5,1.0e-5));
    z6=log(max(z6,1.0e-5));
    z7=log(max(z7,1.0e-5));
    z8=log(max(z8,1.0e-5));
    z9=log(max(z9,1.0e-5));
    z10=log(max(z10,1.0e-5));
    z11=log(max(z11,1.0e-5));
    z12=log(max(z12,1.0e-5));
    z13=log(max(z13,1.0e-5));
    z14=log(max(z14,1.0e-5));
    z15=log(max(z15,1.0e-5));
    z16=log(max(z16,1.0e-5));
    z17=log(max(abs(pdf(normdist,z17)/Var17),1.0e-5));
    z18=log(max(abs(pdf(normdist,z18)/Var18),1.0e-5));
    z19=log(max(abs(pdf(normdist,z19)/Var19),1.0e-5));
    z20=log(max(abs(pdf(normdist,z20)/Var20),1.0e-5));
    z21=log(max(z21,1.0e-5));
    z22=log(max(abs(pdf(normdist,z22)/Var22),1.0e-5));
    z23=log(max(abs(pdf(normdist,z23)/Var23),1.0e-5));
    
    
    //cout << z18 << " z18 " << endl;
    
    //4. Computing the final result
    double result=z1+z2+z3+z4+z5+z6+z7+z8+z9+z10+
    z11+z12+z13+z14+z15+z16+z17+z18+z19+z20+z21+z22+z23;
    
    
    
    //cout << result << " result " << endl;
    result=exp(result);
    return(result);
}


//--------------------
//A.5 Measurement likelihood system for effort in 2010
//--------------------
double MEASloglikEFFORT2010(int Cmeasure1,
                            double aalpha1,
                            double Var1,
                            int Cmeasure2,
                            double aalpha2,
                            double Var2,
                            int Cmeasure3,
                            double aalpha3,
                            double Var3,
                            int Cmeasure4,
                            double aalpha4,
                            double Var4,
                            int Cmeasure5,
                            double aalpha5,
                            double Var5,
                            int Cmeasure6,
                            double aalpha6,
                            double Var6,
                            double logSKILL2012){
    if (logSKILL2012<=1.0e-5){
        logSKILL2012=1.0e-5;
    }
    logSKILL2012=log(logSKILL2012);
    //1. Compute the z score for each measurement
    double z1=ORDPROB(aalpha1,logSKILL2012,Cmeasure1,Var1);
    double z2=ORDPROB(aalpha2,logSKILL2012,Cmeasure2,Var2);
    double z3=ORDPROB(aalpha3,logSKILL2012,Cmeasure3,Var3);
    double z4=ORDPROB(aalpha4,logSKILL2012,Cmeasure4,Var4);
    double z5=ORDPROB(aalpha5,logSKILL2012,Cmeasure5,Var5);
    double z6=ORDPROB(aalpha6,logSKILL2012,Cmeasure6,Var6);
    //2. //Getting a normal distribution to compute the pdf
    
    boost::math::normal_distribution<> normdist (0,1);
    
    
    //3. Compute the likelihood contribution of each measure
    z1=log(max(z1,1.0e-5));
    z2=log(max(z2,1.0e-5));
    z3=log(max(z3,1.0e-5));
    z4=log(max(z4,1.0e-5));
    z5=log(max(z5,1.0e-5));
    z6=log(max(z6,1.0e-5));
    
    
    //4. Computing the final result
    double result=z1+z2+z3+z4+z5+z6;
    result=exp(result);
    //cout << result << " result " << endl;
    return(result);
}



//--------------------
//A.6 Measurement likelihood system for Father's effort in 2012
//--------------------
double MEASloglikFEFFORT2012(double Cmeasure1,
                             double aalpha1,
                             double Var1,
                             double Cmeasure2,
                             double aalpha2,
                             double Var2,
                             double Cmeasure3,
                             double aalpha3,
                             double Var3,
                             double Cmeasure4,
                             double aalpha4,
                             double Var4,
                             double Cmeasure5,
                             double aalpha5,
                             double Var5,
                             double Cmeasure6,
                             double aalpha6,
                             double Var6,
                             double Cmeasure7,
                             double aalpha7,
                             double Var7,
                             double Cmeasure8,
                             double aalpha8,
                             double Var8,
                             double Cmeasure9,
                             double aalpha9,
                             double Var9,
                             double Cmeasure10,
                             double aalpha10,
                             double Var10,
                             double Cmeasure11,
                             double aalpha11,
                             double Var11,
                             double Cmeasure12,
                             double aalpha12,
                             double Var12,
                             double Cmeasure13,
                             double aalpha13,
                             double Var13,
                             double Cmeasure14,
                             double aalpha14,
                             double Var14,
                             double logEF){
    if (logEF<=1.0e-5){
        logEF=1.0e-5;
    }
    logEF=log(logEF);
    
    //1. Compute the z score for each measurement
    double z1=max(abs((Cmeasure1-aalpha1*logEF)/Var1),1.0e-5);
    double z2=max(abs((Cmeasure2-aalpha2*logEF)/Var2),1.0e-5);
    double z3=max(abs((Cmeasure3-aalpha3*logEF)/Var3),1.0e-5);
    double z4=max(abs((Cmeasure4-aalpha4*logEF)/Var4),1.0e-5);
    double z5=max(abs((Cmeasure5-aalpha5*logEF)/Var5),1.0e-5);
    double z6=max(abs((Cmeasure6-aalpha6*logEF)/Var6),1.0e-5);
    double z7=max(abs((Cmeasure7-aalpha7*logEF)/Var7),1.0e-5);
    double z8=max(abs((Cmeasure8-aalpha8*logEF)/Var8),1.0e-5);
    double z9=max(abs((Cmeasure9-aalpha9*logEF)/Var9),1.0e-5);
    double z10=max(abs((Cmeasure10-aalpha10*logEF)/Var10),1.0e-5);
    double z11=max(abs((Cmeasure11-aalpha11*logEF)/Var11),1.0e-5);
    double z12=max(abs((Cmeasure12-aalpha12*logEF)/Var12),1.0e-5);
    double z13=max(abs((Cmeasure13-aalpha13*logEF)/Var13),1.0e-5);
    double z14=max(abs((Cmeasure14-aalpha14*logEF)/Var14),1.0e-5);
    //2. //Getting a normal distribution to compute the pdf
    
    boost::math::normal_distribution<> normdist (0,1);
    
    
    //3. Compute the likelihood contribution of each measure
    z1=log(max(abs(pdf(normdist,z1)/Var1),1.0e-5));
    z2=log(max(abs(pdf(normdist,z2)/Var2),1.0e-5));
    z3=log(max(abs(pdf(normdist,z3)/Var3),1.0e-5));
    z4=log(max(abs(pdf(normdist,z4)/Var4),1.0e-5));
    z5=log(max(abs(pdf(normdist,z5)/Var5),1.0e-5));
    z6=log(max(abs(pdf(normdist,z6)/Var6),1.0e-5));
    z7=log(max(abs(pdf(normdist,z7)/Var7),1.0e-5));
    z8=log(max(abs(pdf(normdist,z8)/Var8),1.0e-5));
    z9=log(max(abs(pdf(normdist,z9)/Var9),1.0e-5));
    z10=log(max(abs(pdf(normdist,z10)/Var10),1.0e-5));
    z11=log(max(abs(pdf(normdist,z11)/Var11),1.0e-5));
    z12=log(max(abs(pdf(normdist,z12)/Var12),1.0e-5));
    z13=log(max(abs(pdf(normdist,z13)/Var13),1.0e-5));
    z14=log(max(abs(pdf(normdist,z14)/Var14),1.0e-5));
    
    
    //4. Computing the final result
    double result=z1+z2+z3+z4+z5+z6+z7+z8+z9+z10+z11+z12+z13+z14;
    result=exp(result);
    //cout << result << " result " << endl;
    return(result);
}




//=========
//A.8 Measures for bargaining power
//=========


double MEASloglikBARG(double Cmeasure1,
                      double Var1,
                      double aalpha1,
                      double Cmeasure2,
                      double Var2,
                      double aalpha2,
                      double Cmeasure3,
                      double Var3,
                      double aalpha3,
                      double Cmeasure4,
                      double Var4,
                      double aalpha4,
                      double Cmeasure5,
                      double Var5,
                      double aalpha5,
                      double Cmeasure6,
                      double Var6,
                      double aalpha6,
                      double Cmeasure7,
                      double Var7,
                      double aalpha7,
                      double Cmeasure8,
                      double Var8,
                      double aalpha8,
                      double Cmeasure9,
                      double Var9,
                      double aalpha9,
                      double Cmeasure10,
                      double Var10,
                      double aalpha10,
                      double Cmeasure11,
                      double Var11,
                      double aalpha11,
                      double Cmeasure12,
                      double Var12,
                      double aalpha12,
                      double Cmeasure13,
                      double Var13,
                      double aalpha13,
                      double Cmeasure14,
                      double Var14,
                      double aalpha14,
                      double Cmeasure15,
                      double Var15,
                      double aalpha15,
                      double Cmeasure16,
                      double Var16,
                      double aalpha16,
                      double Cmeasure17,
                      double Var17,
                      double aalpha17,
                      double Cmeasure18,
                      double Var18,
                      double aalpha18,
                      double Cmeasure19,
                      double Var19,
                      double aalpha19,
                      double logHbarg){
    
    //0. In order to standardize we will do it centered in
    // 0 so that we have between -0.3 and 0.3.
    logHbarg=logHbarg-0.5;
    
    //1. Compute the z score for each measurement
    
    
    double z1=max(abs((Cmeasure1-aalpha1*logHbarg)/Var1),1.0e-5);
    double z2=max(abs((Cmeasure2-aalpha2*logHbarg)/Var2),1.0e-5);
    double z3=max(abs((Cmeasure3-aalpha3*logHbarg)/Var3),1.0e-5);
    double z4=max(abs((Cmeasure4-aalpha4*logHbarg)/Var4),1.0e-5);
    double z5=max(abs((Cmeasure5-aalpha5*logHbarg)/Var5),1.0e-5);
    double z6=max(abs((Cmeasure6-aalpha6*logHbarg)/Var6),1.0e-5);
    double z7=max(abs((Cmeasure7-aalpha7*logHbarg)/Var7),1.0e-5);
    double z8=max(abs((Cmeasure8-aalpha8*logHbarg)/Var8),1.0e-5);
    double z9=max(abs((Cmeasure9-aalpha9*logHbarg)/Var9),1.0e-5);
    double z10=max(abs((Cmeasure10-aalpha10*logHbarg)/Var10),1.0e-5);
    
    double z11=ORDPROB(aalpha11,logHbarg,Cmeasure11,Var11);
    double z12=ORDPROB(aalpha12,logHbarg,Cmeasure12,Var12);
    double z13=ORDPROB(aalpha13,logHbarg,Cmeasure13,Var13);
    
    
    
    double z14=max(abs((Cmeasure14-aalpha14*logHbarg)/Var14),1.0e-5);
    double z15=max(abs((Cmeasure15-aalpha15*logHbarg)/Var15),1.0e-5);
    
    
    double z16=ORDPROB(aalpha16,logHbarg,Cmeasure16,Var16);
    double z17=ORDPROB(aalpha17,logHbarg,Cmeasure17,Var17);
    double z18=ORDPROB(aalpha18,logHbarg,Cmeasure18,Var18);
    double z19=ORDPROB(aalpha19,logHbarg,Cmeasure19,Var19);
    //2. //Getting a normal distribution to compute the pdf
    
    boost::math::normal_distribution<> normdist (0,1);
    
    
    //3. Compute the likelihood contribution of each measure
    z1=log(max(abs(pdf(normdist,z1)/Var1),1.0e-5));
    z2=log(max(abs(pdf(normdist,z2)/Var2),1.0e-5));
    z3=log(max(abs(pdf(normdist,z3)/Var3),1.0e-5));
    z4=log(max(abs(pdf(normdist,z4)/Var4),1.0e-5));
    z5=log(max(abs(pdf(normdist,z5)/Var5),1.0e-5));
    z6=log(max(abs(pdf(normdist,z6)/Var6),1.0e-5));
    z7=log(max(abs(pdf(normdist,z7)/Var7),1.0e-5));
    z8=log(max(abs(pdf(normdist,z8)/Var8),1.0e-5));
    z9=log(max(abs(pdf(normdist,z9)/Var9),1.0e-5));
    z10=log(max(abs(pdf(normdist,z10)/Var10),1.0e-5));
    z11=log(max(z11,1.0e-5));
    z12=log(max(z12,1.0e-5));
    z13=log(max(z13,1.0e-5));
    z14=log(max(abs(pdf(normdist,z14)/Var14),1.0e-5));
    z15=log(max(abs(pdf(normdist,z15)/Var15),1.0e-5));
    z16=log(max(z16,1.0e-5));
    z17=log(max(z17,1.0e-5));
    z18=log(max(z18,1.0e-5));
    z19=log(max(z19,1.0e-5));
    //4. Computing the final result
    double result=z1+z2+z3+z4+z5+z6+z7+z8+z9+z10+
    z11+z12+z13+z14+z15+z16+z17+z18+z19;
    result=exp(result);
    //cout << result << " result " << endl;
    return(result);
}



//A.9 Measurement of primary caregiver characteristics
double MEASloglikPG(double Cmeasure1,
                    double aalpha1,
                    double Var1,
                    double Cmeasure2,
                    double aalpha2,
                    double Var2,
                    double Cmeasure3,
                    double aalpha3,
                    double Var3,
                    double Cmeasure4,
                    double aalpha4,
                    double Var4,
                    double Cmeasure5,
                    double aalpha5,
                    double Var5,
                    double Cmeasure6,
                    double aalpha6,
                    double Var6,
                    double Cmeasure7,
                    double aalpha7,
                    double Var7,
                    double Cmeasure8,
                    double aalpha8,
                    double Var8,
                    double logSKPG){
    
    if (logSKPG==0){
        logSKPG=1.0e-5;
    }
    logSKPG=log(logSKPG);
    //1. Compute the z score for each measurement
    double z1=max(abs((Cmeasure1-aalpha1*logSKPG)/Var1),1.0e-5);
    double z2=max(abs((Cmeasure2-aalpha2*logSKPG)/Var2),1.0e-5);
    double z3=max(abs((Cmeasure3-aalpha3*logSKPG)/Var3),1.0e-5);
    double z4=max(abs((Cmeasure4-aalpha4*logSKPG)/Var4),1.0e-5);
    double z5=max(abs((Cmeasure5-aalpha5*logSKPG)/Var5),1.0e-5);
    double z6=max(abs((Cmeasure6-aalpha6*logSKPG)/Var6),1.0e-5);
    double z7=max(abs((Cmeasure7-aalpha7*logSKPG)/Var7),1.0e-5);
    double z8=max(abs((Cmeasure8-aalpha8*logSKPG)/Var8),1.0e-5);
    
    //2. //Getting a normal distribution to compute the pdf
    
    
    
    boost::math::normal_distribution<> normdist (0,1);
    
    //3. Compute the likelihood contribution of each measure
    
    
    z1=log(max(abs(pdf(normdist,z1)/Var1),1.0e-5));
    z2=log(max(abs(pdf(normdist,z2)/Var2),1.0e-5));
    z3=log(max(abs(pdf(normdist,z3)/Var3),1.0e-5));
    z4=log(max(abs(pdf(normdist,z4)/Var4),1.0e-5));
    z5=log(max(abs(pdf(normdist,z5)/Var5),1.0e-5));
    z6=log(max(abs(pdf(normdist,z6)/Var6),1.0e-5));
    z7=log(max(abs(pdf(normdist,z7)/Var7),1.0e-5));
    z8=log(max(abs(pdf(normdist,z8)/Var8),1.0e-5));
    
    
    
    
    //4. Computing the final result
    // PSI is multiplied by zero (z8) as we dropped to include more
    // people in sample
    double result=z1+z2+z3+z4+z5+z6+z7+0*z8;
    result=exp(result);
    return(result);
}




//=========
//A.10 Measures fOR INVESTMENT IN 2010
//=========


double MEASloginv10(int Cmeasure1,
                    double Var1,
                    double aalpha1,
                    int Cmeasure2,
                    double Var2,
                    double aalpha2,
                    int Cmeasure3,
                    double Var3,
                    double aalpha3,
                    int Cmeasure4,
                    double Var4,
                    double aalpha4,
                    int Cmeasure5,
                    double Var5,
                    double aalpha5,
                    int Cmeasure6,
                    double Var6,
                    double aalpha6,
                    int Cmeasure7,
                    double Var7,
                    double aalpha7,
                    int Cmeasure8,
                    double Var8,
                    double aalpha8,
                    double logINV){
    //0. In order to standardize we will do it centered in
    // 0 so that we have between -0.3 and 0.3.
    
    if (logINV<=1.0e-5){
        logINV=1.0e-5;
    }
    logINV=log(logINV);
    
    //1. Compute the z score for each measurement
    
    double z1=ORDPROB(aalpha1,logINV,Cmeasure1,Var1);
    
    double z2=ORDPROB(aalpha2,logINV,Cmeasure2,Var2);
    
    double z3=ORDPROB(aalpha3,logINV,Cmeasure3,Var3);
    double z4=ORDPROB(aalpha4,logINV,Cmeasure4,Var4);
    double z5=ORDPROB(aalpha5,logINV,Cmeasure5,Var5);
    double z6=ORDPROB(aalpha6,logINV,Cmeasure6,Var6);
    double z7=ORDPROB(aalpha7,logINV,Cmeasure7,Var7);
    double z8=ORDPROB(aalpha8,logINV,Cmeasure8,Var8);
    
    //2. //Getting a normal distribution to compute the pdf
    
    boost::math::normal_distribution<> normdist (0,1);
    
    
    //3. Compute the likelihood contribution of each measure
    
    z1=log(max(z1,1.0e-5));
    z2=log(max(z2,1.0e-5));
    z3=log(max(z3,1.0e-5));
    z4=log(max(z4,1.0e-5));
    z5=log(max(z5,1.0e-5));
    z6=log(max(z6,1.0e-5));
    z7=log(max(z7,1.0e-5));
    z8=log(max(z8,1.0e-5));
    
    //4. Computing the final result
    double result=z1+z2+z3+z4+z5+z6+z7+z8;
    result=exp(result);
    //cout << result << " result " << endl;
    return(result);
}


//Invetsment in 2012
//=========
//A.11 Measures fOR INVESTMENT IN 2012
//=========


double MEASloginv12(double Cmeasure1,
                    double Var1,
                    double aalpha1,
                    double Cmeasure2,
                    double Var2,
                    double aalpha2,
                    double Cmeasure3,
                    double Var3,
                    double aalpha3,
                    double Cmeasure4,
                    double Var4,
                    double aalpha4,
                    double Cmeasure5,
                    double Var5,
                    double aalpha5,
                    double Cmeasure6,
                    double Var6,
                    double aalpha6,
                    double Cmeasure7,
                    double Var7,
                    double aalpha7,
                    double Cmeasure8,
                    double Var8,
                    double aalpha8,
                    double Cmeasure9,
                    double Var9,
                    double aalpha9,
                    double Cmeasure10,
                    double Var10,
                    double aalpha10,
                    double Cmeasure11,
                    double Var11,
                    double aalpha11,
                    int Cmeasure12,
                    double Var12,
                    double aalpha12,
                    int Cmeasure13,
                    double Var13,
                    double aalpha13,
                    int Cmeasure14,
                    double Var14,
                    double aalpha14,
                    int Cmeasure15,
                    double Var15,
                    double aalpha15,
                    int Cmeasure16,
                    double Var16,
                    double aalpha16,
                    int Cmeasure17,
                    double Var17,
                    double aalpha17,
                    int Cmeasure18,
                    double Var18,
                    double aalpha18,
                    int Cmeasure19,
                    double Var19,
                    double aalpha19,
                    double Cmeasure20,
                    double Var20,
                    double aalpha20,
                    double Cmeasure21,
                    double Var21,
                    double aalpha21,
                    double logINV){
    //0. In order to standardize we will do it centered in
    // 0 so that we have between -0.3 and 0.3.
    
    if (logINV<=1.0e-5){
        logINV=1.0e-5;
    }
    logINV=log(logINV);
    
    //1. Compute the z score for each measurement
    double z1=max(abs((Cmeasure1-aalpha1*logINV)/Var1),1.0e-5);
    double z2=max(abs((Cmeasure2-aalpha2*logINV)/Var2),1.0e-5);
    double z3=max(abs((Cmeasure3-aalpha3*logINV)/Var3),1.0e-5);
    double z4=max(abs((Cmeasure4-aalpha4*logINV)/Var4),1.0e-5);
    double z5=max(abs((Cmeasure5-aalpha5*logINV)/Var5),1.0e-5);
    double z6=max(abs((Cmeasure6-aalpha6*logINV)/Var6),1.0e-5);
    double z7=max(abs((Cmeasure7-aalpha7*logINV)/Var7),1.0e-5);
    double z8=max(abs((Cmeasure8-aalpha8*logINV)/Var8),1.0e-5);
    double z9=max(abs((Cmeasure9-aalpha9*logINV)/Var9),1.0e-5);
    double z10=max(abs((Cmeasure10-aalpha10*logINV)/Var10),1.0e-5);
    double z11=max(abs((Cmeasure11-aalpha11*logINV)/Var11),1.0e-5);
    
    double z12=ORDPROB(aalpha12,logINV,Cmeasure12,Var12);
    double z13=ORDPROB(aalpha13,logINV,Cmeasure13,Var13);
    double z14=ORDPROB(aalpha14,logINV,Cmeasure14,Var14);
    double z15=ORDPROB(aalpha15,logINV,Cmeasure15,Var15);
    double z16=ORDPROB(aalpha16,logINV,Cmeasure15,Var16);
    double z17=ORDPROB(aalpha17,logINV,Cmeasure16,Var17);
    double z18=ORDPROB(aalpha18,logINV,Cmeasure17,Var18);
    double z19=ORDPROB(aalpha19,logINV,Cmeasure18,Var19);
    double z20=max(abs((Cmeasure20-aalpha20*logINV)/Var20),1.0e-5);
    double z21=max(abs((Cmeasure21-aalpha21*logINV)/Var21),1.0e-5);
    
    
    //2. //Getting a normal distribution to compute the pdf
    
    boost::math::normal_distribution<> normdist (0,1);
    
    
    //3. Compute the likelihood contribution of each measure
    z1=log(max(abs(pdf(normdist,z1)/Var1),1.0e-5));
    z2=log(max(abs(pdf(normdist,z2)/Var2),1.0e-5));
    z3=log(max(abs(pdf(normdist,z3)/Var3),1.0e-5));
    z4=log(max(abs(pdf(normdist,z4)/Var4),1.0e-5));
    z5=log(max(abs(pdf(normdist,z5)/Var5),1.0e-5));
    z6=log(max(abs(pdf(normdist,z6)/Var6),1.0e-5));
    z7=log(max(abs(pdf(normdist,z7)/Var7),1.0e-5));
    z8=log(max(abs(pdf(normdist,z8)/Var8),1.0e-5));
    z9=log(max(abs(pdf(normdist,z9)/Var9),1.0e-5));
    z10=log(max(abs(pdf(normdist,z10)/Var10),1.0e-5));
    z11=log(max(abs(pdf(normdist,z11)/Var11),1.0e-5));
    
    z12=log(max(z12,1.0e-5));
    z13=log(max(z13,1.0e-5));
    z14=log(max(z14,1.0e-5));
    z15=log(max(z15,1.0e-5));
    z16=log(max(z16,1.0e-5));
    z17=log(max(z17,1.0e-5));
    z18=log(max(z18,1.0e-5));
    z19=log(max(z19,1.0e-5));
    z20=log(max(abs(pdf(normdist,z20)/Var20),1.0e-5));
    z21=log(max(abs(pdf(normdist,z21)/Var21),1.0e-5));
    
    //4. Computing the final result
    double result=z1+z2+z3+z4+z5+z6+z7+z8+z9+z10+z11+z12+z13+z14+z15+z16+z17+z18+z19+z20+z21;
    result=exp(result);
    return(result);
}






//This one is not for log-normal but simply normal
// [[Rcpp::export]]
long double F_likgeneric(double observed, double predicted, double variance){
    
    
    if (variance<=1.0e-5){
        variance=1.0e-5;
    }
    //Getting a normal distribution to compute the pdf
    boost::math::normal_distribution<> normdist (0,1);
    //Getting the estimate
    double likelihood=(observed)-(predicted);
    likelihood=likelihood/variance;
    likelihood=pdf(normdist,likelihood);
    likelihood=likelihood/variance;
    if (likelihood<=1.0e-5){
        likelihood=1.0e-5;
    }
    return(likelihood);
}


//----------------------------------------------------
//2. Block of defining likelihood for all the dataset
//---------------------------------------------------
// [[Rcpp::export]]
vector<vector<long double> > F_likelihoodSCALAR(vector<double> par_likelihood,
                                                double  Cliveswithfather12,
                                                double Cliveswithmother12,
                                                double Hhchores12,
                                                double Meffort12,
                                                double Feffort12,
                                                double Mage12,
                                                double Fage12,
                                                double Cedad_meses12,
                                                double Ctestsfactorsss2012,
                                                double Ctestsfactor2ss_10,
                                                double Myrschool12,
                                                double Fyrschool12,
                                                double Ffraclabor12,
                                                double Mfraclabor12,
                                                double Mwage12,
                                                double Fwage12,
                                                double Mnly12,
                                                double Fnly12,
                                                double Cfactorinv12,
                                                double Cchildcare12,
                                                double Ccareskills12,
                                                double Hbarg,
                                                double Feffort10,
                                                double Meffort10,
                                                double Cfactorinv10,
                                                double Cbirthfactor,
                                                double Mwage10,
                                                double Fwage10,
                                                double Mnly10,
                                                double Fnly10,
                                                double Cedad_meses10,
                                                double Cliveswithfather10,
                                                double Cliveswithmother10,
                                                double Mage10,
                                                double Fage10,
                                                double Mfraclabor10,
                                                double Ffraclabor10,
                                                double Cchildcare10,
                                                double Hchildcareobs,
                                                double PG,
                                                double MTJH,
                                                double MRATIO,
                                                double Unemployment,
                                                double Wageratio,
                                                double  Distance,
                                                double Magegroup10,
                                                double  Magegroup12,
                                                double  Hmemberstotal10,
                                                double  Hmemberstotal12,
                                                double CSTDtepsi_pb_coo10,
                                                double CSTDtepsi_pb_len10,
                                                double CSTDtepsi_pb_mot10,
                                                double CSTDtvip_pb10,
                                                double CSTDcbcl1_pb_110,
                                                double CSTDcbcl1_pb_210,
                                                double CSTDcbcl1_pb_310,
                                                double CSTDcbcl1_pb_410,
                                                double CSTDcbcl1_pb_510,
                                                double CSTDcbcl1_pb_610,
                                                double CSTDcbcl1_pb_710,
                                                double CSTDtadi_pb_cog12,
                                                double CSTDtadi_pb_mot12,
                                                double CSTDtadi_pb_len12,
                                                double CSTDtadi_pb_se12,
                                                double CSTDbt_112,
                                                double CSTDbt_212,
                                                double CSTDbt_312,
                                                double CSTDbt_412,
                                                double CSTDbt_512,
                                                double CSTDbt_t12,
                                                double CSTDhtks_st12,
                                                double CSTDbdst_st12,
                                                double CSTDppvt_t12,
                                                int Ccondpregg3_1,
                                                int Ccondpregg3_2,
                                                int Ccondpregg3_3,
                                                int Ccondpregg3_4,
                                                int Ccondpregg3_5,
                                                int Ccondpregg3_6,
                                                int Ccondpregg3_7,
                                                int Ccondpregg3_8,
                                                int Ccondpregg3_9,
                                                int Ccondpregg4a_1,
                                                int Ccondpregg4a_2,
                                                int Ccondpregg4a_3,
                                                int Ccondpregg4a_4,
                                                int Ccondpregg4a_5,
                                                int Ccondpregg4a_6,
                                                int Ccondpregg4a_7,
                                                double Ccondpregg7b,
                                                double Ccondpregg8b,
                                                double Ccondpregg9,
                                                double Ccondpregg11b,
                                                int Ccondpregpreterm,
                                                double Ccondpregg24,
                                                double Ccondpregg23,
                                                double Cwais_pb_num,
                                                double Cwais_pb_vo,
                                                double Cbfi_pb_ama,
                                                double Cbfi_pb_ape,
                                                double Cbfi_pb_ext,
                                                double Cbfi_pb_neu,
                                                double Cbfi_pb_res,
                                                double Cpsi_pb_total,
                                                double Hbargg2a,
                                                double Hbargg2b,
                                                double Hbargg2c,
                                                double Hbargg2d,
                                                double Hbargg2e,
                                                double Hbargg2f,
                                                double Hbargg2g,
                                                double Hbargg2h,
                                                double Hbargg2i,
                                                double Hbargg2j,
                                                int Hwomdecide12,
                                                int Hfathdecides12,
                                                int Hbothdecide,
                                                double Hcaresacwom,
                                                double Hcaresacman,
                                                int Hopofman121,
                                                int Hopofman122,
                                                int Hopofman123,
                                                int Hopofman124,
                                                int Cinvh3,
                                                int Cinvh4,
                                                int Cinvh5,
                                                int Cinvh6,
                                                int Cinvh7,
                                                int Cinvh8,
                                                int Cinvh9,
                                                int Cinvh13,
                                                double Cinvf11a,
                                                double Cinvf11b,
                                                double Cinvf11c,
                                                double Cinvf11d,
                                                double Cinvf11e,
                                                double Cinvf11f,
                                                double Cinvf11g,
                                                double Cinvf11h,
                                                double Cinvf11i,
                                                double Cinvf11j,
                                                double Cinvf11k,
                                                int Cinvhm2_10,
                                                int Cinvhm2_11,
                                                int Cinvhm2_12,
                                                int Cinvhm2_13,
                                                int Cinvhm2_14,
                                                int Cinvhm2_15,
                                                int Cinvhm2_16,
                                                int Cinvhm2_18,
                                                double Csharesbedroomhowmany12,
                                                double Csharesbedhowmany12,
                                                int g42_a2,
                                                int g42_b2,
                                                int g42_c2,
                                                int g42_d2,
                                                int g42_e2,
                                                int g42_f2,
                                                int g42_a1,
                                                int g42_b1,
                                                int g42_c1,
                                                int g42_d1,
                                                int g42_e1,
                                                int g42_f1,
                                                double f21a_p_t,
                                                double f21b_p_t,
                                                double f21c_p_t,
                                                double f21d_p_t,
                                                double f21e_p_t,
                                                double f21f_p_t,
                                                double f21g_p_t,
                                                double f21h_p_t,
                                                double f21i_p_t,
                                                double f21j_p_t,
                                                double f21k_p_t,
                                                double f21l_p_t,
                                                double f21m_p_t,
                                                double f21n_p_t,
                                                double f21a_m_t,
                                                double f21b_m_t,
                                                double f21c_m_t,
                                                double f21d_m_t,
                                                double f21e_m_t,
                                                double f21f_m_t,
                                                double f21g_m_t,
                                                double f21h_m_t,
                                                double f21i_m_t,
                                                double f21j_m_t,
                                                double f21k_m_t,
                                                double f21l_m_t,
                                                double f21m_m_t,
                                                double f21n_m_t,
                                                double HDens,
                                                int SEED){
    
    //Before everything I will specify whether or not
    //we will run the SMOOTHING DISTRIBUTION
    int SMOOTHINGFILTER=0;
    
    
    //------------------------------------------
    //0. Performing the right transformations to
    //guarantee that parameters lie within the
    //specified set
    //-----------------------------------------
    
    //0.1 Aalphas
    //0.1 Aalphas
    double aalpha1m=FT_approach(FT_exp(par_likelihood[0])/
                                (FT_exp(par_likelihood[0])+FT_exp(par_likelihood[1])+FT_exp(par_likelihood[2])+FT_exp(par_likelihood[3])));
    
    double aalpha2m=FT_approach(FT_exp(par_likelihood[1])/
                                (FT_exp(par_likelihood[0])+FT_exp(par_likelihood[1])+FT_exp(par_likelihood[2])+FT_exp(par_likelihood[3])));
    
    double aalpha3m=FT_approach(FT_exp(par_likelihood[2])/
                                (FT_exp(par_likelihood[0])+FT_exp(par_likelihood[1])+FT_exp(par_likelihood[2])+FT_exp(par_likelihood[3])));
    
    double aalpha40m=FT_approach(FT_exp(par_likelihood[3])/
                                 (FT_exp(par_likelihood[0])+FT_exp(par_likelihood[1])+FT_exp(par_likelihood[2])+FT_exp(par_likelihood[3])));
    
    
    
    
    double aalpha41m=FT_approach(FT_exp(par_likelihood[4])/
                                 (1+FT_exp(par_likelihood[4])))*aalpha40m;
    
    double aalpha4m=0;
    
    
    
    
    
    
    //0.2 Skills
    double ttheta0=FT_exp(par_likelihood[5])/
    (FT_exp(par_likelihood[5])+FT_exp(par_likelihood[6])+FT_exp(par_likelihood[7]));
    
    double ttheta1=FT_exp(par_likelihood[6])/
    (FT_exp(par_likelihood[5])+FT_exp(par_likelihood[6])+FT_exp(par_likelihood[7]));
    
    double ttheta2=FT_exp(par_likelihood[7])/
    (FT_exp(par_likelihood[5])+FT_exp(par_likelihood[6])+FT_exp(par_likelihood[7]));
    
    double pphi=FT_pphi(par_likelihood[8]);
    double stdskills=FT_exp(par_likelihood[9]);
    //double stdskills=0.1;
    
    
    double ggammaf=FT_gf(FT_exp(par_likelihood[10])/
                         (FT_exp(par_likelihood[10])+FT_exp(par_likelihood[11])));
    
    
    double ggammam=FT_gf(FT_exp(par_likelihood[11])/
                         (FT_exp(par_likelihood[10])+FT_exp(par_likelihood[11])));
    
    
    //0.3 Wages for mother
    double bbeta0m=par_likelihood[12];
    double bbeta1m=par_likelihood[13];
    double bbeta2m=par_likelihood[14];
    double bbeta3m=par_likelihood[15];
    double stdwm=FT_exp(par_likelihood[16]);
    
    //0.4 Measurement system
    double stdeffortmom=FT_exp(par_likelihood[17]);
    double stdinvestment=FT_exp(par_likelihood[18]);
    
    //0.5 Price of investment
    double priceINV0=FT_exp(par_likelihood[19]);
    
    double priceINV1=FT_exp(par_likelihood[324]);
    
    double priceINV=0;
    
    //0.6 Childcare situation
    double ddelta0=par_likelihood[20]; //Intercept
    //double ddelta1=0.002; //Age shifter
    double ddelta1=par_likelihood[21]; //Age shifter
    double ddelta2=par_likelihood[22]; //Childcare production
    //double ddelta2=0.1; //Childcare production
    double ddelta3=par_likelihood[23]; //PG in 2010
    double ddelta3_12=par_likelihood[24]; //Childcare production PG in 2012
    
    
    
    
    //Preferences for the father
    //0.1 Aalphas
    double aalpha1f=FT_approach(FT_exp(par_likelihood[25])/
                                (FT_exp(par_likelihood[25])+FT_exp(par_likelihood[26])+
                                 FT_exp(par_likelihood[27])+FT_exp(par_likelihood[28])));
    
    double aalpha2f=FT_approach(FT_exp(par_likelihood[26])/
                                (FT_exp(par_likelihood[25])+FT_exp(par_likelihood[26])+
                                 FT_exp(par_likelihood[27])+FT_exp(par_likelihood[28])));
    
    double aalpha3f=FT_approach(FT_exp(par_likelihood[27])/
                                (FT_exp(par_likelihood[25])+FT_exp(par_likelihood[26])+
                                 FT_exp(par_likelihood[27])+FT_exp(par_likelihood[28])));
    
    double aalpha40f=FT_approach(FT_exp(par_likelihood[28])/
                                 (FT_exp(par_likelihood[25])+FT_exp(par_likelihood[26])+
                                  FT_exp(par_likelihood[27])+FT_exp(par_likelihood[28])));
    
    double aalpha41f=FT_approach(FT_exp(par_likelihood[29])/
                                 (1+FT_exp(par_likelihood[29])))*aalpha40f;
    double aalpha4f=0;
    
    
    
    //0.3 Wages
    double bbeta0f=par_likelihood[30];
    double bbeta1f=par_likelihood[31];
    double bbeta2f=par_likelihood[32];
    double bbeta3f=par_likelihood[33];
    double stdwf=FT_exp(par_likelihood[34]);
    
    //0.4 Measurement system
    double stdeffortfat=FT_exp(par_likelihood[35]);
    
    
    //0.6 Bargaining power parameters
    double llambda0=par_likelihood[36]; //Intercept
    double llambda1=par_likelihood[37]; //Wage ratio
    double llambda2=par_likelihood[38]; //Non-labor income ratio
    double llambda3=par_likelihood[39]; //Age difference
    double llambda4=par_likelihood[40]; //Education difference
    double stdmmu=FT_exp(par_likelihood[41]); //STD of barg shock
    double mmuLB=0.2; //Lower bound for bargaining power
    double mmuUB=0.8; //Upper bound for bargaining power
    double mupred=0; //Predicted level of bargaining power
    double mupred10=0; //Predicted level of bargaining power 2010
    //0.7 Measurement system parameters
    
    double MEASSkills=FT_exp(par_likelihood[42]);
    
    
    //0.8 Preferences in the first period
    //Mother
    double aalpha1m10=FT_approach(FT_exp(par_likelihood[47])/
                                  (FT_exp(par_likelihood[47])+
                                   FT_exp(par_likelihood[48])+
                                   FT_exp(par_likelihood[49])+
                                   FT_exp(par_likelihood[50])+
                                   FT_exp(par_likelihood[51])));
    
    double aalpha2m10=FT_approach(FT_exp(par_likelihood[48])/
                                  (FT_exp(par_likelihood[47])+
                                   FT_exp(par_likelihood[48])+
                                   FT_exp(par_likelihood[49])+
                                   FT_exp(par_likelihood[50])+
                                   FT_exp(par_likelihood[51])));
    
    double aalpha3m10=FT_approach(FT_exp(par_likelihood[49])/
                                  (FT_exp(par_likelihood[47])+
                                   FT_exp(par_likelihood[48])+
                                   FT_exp(par_likelihood[49])+
                                   FT_exp(par_likelihood[50])+
                                   FT_exp(par_likelihood[51])));
    
    double aalpha40m10=FT_approach(FT_exp(par_likelihood[50])/
                                   (FT_exp(par_likelihood[47])+
                                    FT_exp(par_likelihood[48])+
                                    FT_exp(par_likelihood[49])+
                                    FT_exp(par_likelihood[50])+
                                    FT_exp(par_likelihood[51])));
    
    double aalpha5m10=FT_approach(FT_exp(par_likelihood[51])/
                                  (FT_exp(par_likelihood[47])+
                                   FT_exp(par_likelihood[48])+
                                   FT_exp(par_likelihood[49])+
                                   FT_exp(par_likelihood[50])+
                                   FT_exp(par_likelihood[51])));
    
    double aalpha41m10=FT_approach(FT_exp(par_likelihood[52])/
                                   (1+FT_exp(par_likelihood[52])))*aalpha40m10;
    double aalpha4m10=0;
    //Father
    double aalpha1f10=FT_approach(FT_exp(par_likelihood[53])/
                                  (FT_exp(par_likelihood[53])+
                                   FT_exp(par_likelihood[54])+
                                   FT_exp(par_likelihood[55])+
                                   FT_exp(par_likelihood[56])+
                                   FT_exp(par_likelihood[57])));
    
    double aalpha2f10=FT_approach(FT_exp(par_likelihood[54])/
                                  (FT_exp(par_likelihood[53])+
                                   FT_exp(par_likelihood[54])+
                                   FT_exp(par_likelihood[55])+
                                   FT_exp(par_likelihood[56])+
                                   FT_exp(par_likelihood[57])));
    
    double aalpha3f10=FT_approach(FT_exp(par_likelihood[55])/
                                  (FT_exp(par_likelihood[53])+
                                   FT_exp(par_likelihood[54])+
                                   FT_exp(par_likelihood[55])+
                                   FT_exp(par_likelihood[56])+
                                   FT_exp(par_likelihood[57])));
    
    
    double aalpha40f10=FT_approach(FT_exp(par_likelihood[56])/
                                   (FT_exp(par_likelihood[53])+
                                    FT_exp(par_likelihood[54])+
                                    FT_exp(par_likelihood[55])+
                                    FT_exp(par_likelihood[56])+
                                    FT_exp(par_likelihood[57])));
    
    double aalpha5f10=FT_approach(FT_exp(par_likelihood[57])/
                                  (FT_exp(par_likelihood[53])+
                                   FT_exp(par_likelihood[54])+
                                   FT_exp(par_likelihood[55])+
                                   FT_exp(par_likelihood[56])+
                                   FT_exp(par_likelihood[57])));
    
    
    double aalpha41f10=FT_approach(FT_exp(par_likelihood[58])/
                                   (1+FT_exp(par_likelihood[58])))*aalpha40m10;
    double aalpha4f10=0;
    
    //Price of childcare
    double pchildcare0=FT_exp(par_likelihood[59]);
    double pchildcare1=FT_exp(par_likelihood[60]); //Regarding the Distance
    
    //Preference shocks
    double MshockWA=FT_exp(par_likelihood[61]);
    double MshockNWA=FT_exp(par_likelihood[62]);
    double MshockWNA=FT_exp(par_likelihood[63]);
    double MshockNWNA=FT_exp(par_likelihood[64]);
    double FshockWA=FT_exp(par_likelihood[65]);
    double FshockNWA=FT_exp(par_likelihood[66]);
    double FshockWNA=FT_exp(par_likelihood[67]);
    double FshockNWNA=FT_exp(par_likelihood[68]);
    
    
    //Final elements to take into account:
    //Distribution factors
    double llambda5=par_likelihood[69];//Male  to male ratio
    double llambda6=par_likelihood[70];//Unemployment ratio
    double llambda7=par_likelihood[71];//Wage ratio
    double llambda8=par_likelihood[72];//Distance a centros de protección
    
    //Preference heterogeneity
    
    //Preference heterogeneity terms will have the same restrictions as aalpha41.
    //Between zero and one.
    
    //Mujer trabajadora y jefa de hogar 2010
    double aalpha3m10_mtjh=FT_approach(FT_exp(par_likelihood[73])/
                                       (1+FT_exp(par_likelihood[73])))*aalpha3m10;
    
    //Mujer trabajadora y jefa de hogar 2012
    double aalpha3m12_mtjh=FT_approach(FT_exp(par_likelihood[74])/
                                       (1+FT_exp(par_likelihood[74])))*aalpha3m;
    
    //Age preferences for leisure
    double aalpha3mage10=FT_approach(FT_exp(par_likelihood[75])/
                                     (1+FT_exp(par_likelihood[75])))*aalpha3m10;
    
    double aalpha3mage12=FT_approach(FT_exp(par_likelihood[76])/
                                     (1+FT_exp(par_likelihood[76])))*aalpha3m;
    
    
    //People with children
    double ddelta4=par_likelihood[77]; //How having one additional person affects exp(R).
    
    //Loading parameters for the measurement system
    //1. Skills in 2010
    //1.1 CSTDtepsi_pb_coo10
    double M_A1_S2010=par_likelihood[78];
    double M_VS1_S2010=max(exp(par_likelihood[79]),0.0001);
    
    //1.2 CSTDtepsi_pb_len10
    double M_A2_S2010=par_likelihood[80];
    double M_VS2_S2010=max(exp(par_likelihood[81]),0.0001);
    
    //1.3 CSTDtepsi_pb_mot10
    double M_A3_S2010=par_likelihood[82];
    double M_VS3_S2010=max(exp(par_likelihood[83]),0.0001);
    
    
    //1.4 CSTDtvip_pb10
    double M_A4_S2010=par_likelihood[84];
    double M_VS4_S2010=max(exp(par_likelihood[85]),0.0001);
    
    
    //1.5 CSTDcbcl1_pb_110
    double M_A5_S2010=par_likelihood[86];
    double M_VS5_S2010=max(exp(par_likelihood[87]),0.0001);
    
    
    //1.6 CSTDcbcl1_pb_210
    double M_A6_S2010=par_likelihood[88];
    double M_VS6_S2010=max(exp(par_likelihood[89]),0.0001);
    
    
    //1.7 CSTDcbcl1_pb_310
    double M_A7_S2010=par_likelihood[90];
    double M_VS7_S2010=max(exp(par_likelihood[91]),0.0001);
    
    
    //1.8 CSTDcbcl1_pb_410
    double M_A8_S2010=par_likelihood[92];
    double M_VS8_S2010=max(exp(par_likelihood[93]),0.0001);
    
    
    
    //1.9 CSTDcbcl1_pb_510
    double M_A9_S2010=par_likelihood[94];
    double M_VS9_S2010=max(exp(par_likelihood[95]),0.0001);
    
    
    
    //1.10 CSTDcbcl1_pb_610
    double M_A10_S2010=par_likelihood[96];
    double M_VS10_S2010=max(exp(par_likelihood[97]),0.0001);
    
    
    
    //1.11 CSTDcbcl1_pb_710
    double M_A11_S2010=par_likelihood[98];
    double M_VS11_S2010=max(exp(par_likelihood[99]),0.0001);
    
    
    
    
    //Loading parameters for the measurement system in skills 2012
    //1. Skills in 2012
    //1.1 CSTDtadi_pb_cog12
    double M_A1_S2012=par_likelihood[100];
    double M_VS1_S2012=max(exp(par_likelihood[101]),0.0001);
    
    
    //1.2 CSTDtadi_pb_mot12
    double M_A2_S2012=par_likelihood[102];
    double M_VS2_S2012=max(exp(par_likelihood[103]),0.0001);
    
    
    //1.3 CSTDtadi_pb_len12
    double M_A3_S2012=par_likelihood[104];
    double M_VS3_S2012=max(exp(par_likelihood[105]),0.0001);
    
    
    //1.4 CSTDtadi_pb_se12
    double M_A4_S2012=par_likelihood[106];
    double M_VS4_S2012=max(exp(par_likelihood[107]),0.0001);
    
    
    //1.5 CSTDbt_112
    double M_A5_S2012=par_likelihood[108];
    double M_VS5_S2012=max(exp(par_likelihood[109]),0.0001);
    
    
    //1.6 CSTDbt_212
    double M_A6_S2012=par_likelihood[110];
    double M_VS6_S2012=max(exp(par_likelihood[111]),0.0001);
    
    
    //1.7 CSTDbt_312
    double M_A7_S2012=par_likelihood[112];
    double M_VS7_S2012=max(exp(par_likelihood[113]),0.0001);
    
    
    //1.8 CSTDbt_412
    double M_A8_S2012=par_likelihood[114];
    double M_VS8_S2012=max(exp(par_likelihood[115]),0.0001);
    
    
    //1.9 CSTDbt_512
    double M_A9_S2012=par_likelihood[116];
    double M_VS9_S2012=max(exp(par_likelihood[117]),0.0001);
    
    
    //1.10 CSTDbt_t12
    double M_A10_S2012=par_likelihood[118];
    double M_VS10_S2012=max(exp(par_likelihood[119]),0.0001);
    
    
    //1.11 CSTDhtks_st12
    double M_A11_S2012=par_likelihood[120];
    double M_VS11_S2012=max(exp(par_likelihood[121]),0.0001);
    
    
    
    //1.12 CSTDbdst_st12
    double M_A12_S2012=par_likelihood[122];
    double M_VS12_S2012=max(exp(par_likelihood[123]),0.0001);
    
    
    //1.13 CSTDppvt_t12
    double M_A13_S2012=par_likelihood[124];
    double M_VS13_S2012=max(exp(par_likelihood[125]),0.0001);
    
    
    
    //Loading parameters for the measurement system in skills at birth
    //2. Skills a t birth
    //-----------------------------------------
    
    //1. CCONDPREGG3_1
    double M_A1_SBIRTH=par_likelihood[126];
    double M_VS1_SBIRTH=max(exp(par_likelihood[127]),0.0001);
    
    
    //2. CCONDPREGG3_2
    double M_A2_SBIRTH=par_likelihood[128];
    double M_VS2_SBIRTH=max(exp(par_likelihood[129]),0.0001);
    
    
    
    //3. Ccondpregg3_3
    double M_A3_SBIRTH=par_likelihood[130];
    double M_VS3_SBIRTH=max(exp(par_likelihood[131]),0.0001);
    
    
    
    //4. Ccondpregg3_4
    double M_A4_SBIRTH=par_likelihood[132];
    double M_VS4_SBIRTH=max(exp(par_likelihood[133]),0.0001);
    
    
    
    //5. Ccondpregg3_5
    double M_A5_SBIRTH=par_likelihood[134];
    double M_VS5_SBIRTH=max(exp(par_likelihood[135]),0.0001);
    
    
    
    //6. Ccondpregg3_6
    double M_A6_SBIRTH=par_likelihood[136];
    double M_VS6_SBIRTH=max(exp(par_likelihood[137]),0.0001);
    
    
    
    //7. Ccondpregg3_7
    double M_A7_SBIRTH=par_likelihood[138];
    double M_VS7_SBIRTH=max(exp(par_likelihood[139]),0.0001);
    
    
    
    //8. Ccondpregg3_8
    double M_A8_SBIRTH=par_likelihood[140];
    double M_VS8_SBIRTH=max(exp(par_likelihood[141]),0.0001);
    
    
    
    //9. Ccondpregg3_9
    double M_A9_SBIRTH=par_likelihood[142];
    double M_VS9_SBIRTH=max(exp(par_likelihood[143]),0.0001);
    
    
    
    //10. Ccondpregg4a_1
    double M_A10_SBIRTH=par_likelihood[144];
    double M_VS10_SBIRTH=max(exp(par_likelihood[145]),0.0001);
    
    
    
    //11. Ccondpregg4a_2
    double M_A11_SBIRTH=par_likelihood[146];
    double M_VS11_SBIRTH=max(exp(par_likelihood[147]),0.0001);
    
    
    
    //12. Ccondpregg4a_3
    double M_A12_SBIRTH=par_likelihood[148];
    double M_VS12_SBIRTH=max(exp(par_likelihood[149]),0.0001);
    
    
    
    //13. Ccondpregg4a_4
    double M_A13_SBIRTH=par_likelihood[150];
    double M_VS13_SBIRTH=max(exp(par_likelihood[151]),0.0001);
    
    
    
    //14. Ccondpregg4a_5
    double M_A14_SBIRTH=par_likelihood[152];
    double M_VS14_SBIRTH=max(exp(par_likelihood[153]),0.0001);
    
    
    
    //15. Ccondpregg4a_6
    double M_A15_SBIRTH=par_likelihood[154];
    double M_VS15_SBIRTH=max(exp(par_likelihood[155]),0.0001);
    
    
    
    //16. Ccondpregg4a_7
    double M_A16_SBIRTH=par_likelihood[156];
    double M_VS16_SBIRTH=max(exp(par_likelihood[157]),0.0001);
    
    
    
    //17. Ccondpregg7b
    double M_A17_SBIRTH=par_likelihood[158];
    double M_VS17_SBIRTH=max(exp(par_likelihood[159]),0.0001);
    
    
    
    //18. Ccondpregg8b
    double M_A18_SBIRTH=par_likelihood[160];
    double M_VS18_SBIRTH=max(exp(par_likelihood[161]),0.0001);
    
    
    
    //19. Ccondpregg9
    double M_A19_SBIRTH=par_likelihood[162];
    double M_VS19_SBIRTH=max(exp(par_likelihood[163]),0.0001);
    
    
    
    //20. Ccondpregg11b
    double M_A20_SBIRTH=par_likelihood[164];
    double M_VS20_SBIRTH=max(exp(par_likelihood[165]),0.0001);
    
    
    
    //21. Ccondpregpreterm
    double M_A21_SBIRTH=par_likelihood[166];
    double M_VS21_SBIRTH=max(exp(par_likelihood[167]),0.0001);
    
    
    //22. Ccondpregg24
    double M_A22_SBIRTH=par_likelihood[168];
    double M_VS22_SBIRTH=max(exp(par_likelihood[169]),0.0001);
    
    
    //23. Ccondpregg23
    double M_A23_SBIRTH=par_likelihood[170];
    double M_VS23_SBIRTH=max(exp(par_likelihood[171]),0.0001);
    
    
    //4. Skills of primary caregiver.
    //4.1 Cwais_pb_num
    double M_A1_SPG=par_likelihood[172];
    double M_V1_SPG=max(exp(par_likelihood[173]),0.0001);
    
    //4.2 Cwais_pb_vo
    double M_A2_SPG=par_likelihood[174];
    double M_V2_SPG=max(exp(par_likelihood[175]),0.0001);
    
    //4.3 Cbfi_pb_ama
    double M_A3_SPG=par_likelihood[176];
    double M_V3_SPG=max(exp(par_likelihood[177]),0.0001);
    
    //4.4 Cbfi_pb_ape
    double M_A4_SPG=par_likelihood[178];
    double M_V4_SPG=max(exp(par_likelihood[179]),0.0001);
    
    //4.5 Cbfi_pb_ext
    double M_A5_SPG=par_likelihood[180];
    double M_V5_SPG=max(exp(par_likelihood[181]),0.0001);
    
    //4.6 Cbfi_pb_neu
    double M_A6_SPG=par_likelihood[182];
    double M_V6_SPG=max(exp(par_likelihood[183]),0.0001);
    
    //4.7 Cbfi_pb_res
    double M_A7_SPG=par_likelihood[184];
    double M_V7_SPG=max(exp(par_likelihood[185]),0.0001);
    
    //4.8 Cpsi_pb_total
    double M_A8_SPG=par_likelihood[186];
    double M_V8_SPG=max(exp(par_likelihood[187]),0.0001);
    
    //------------------------------------------
    //5 BARGAINING POWER
    
    //5.1. Hbargg2a
    double M_A1_BARG=par_likelihood[188];
    double M_V1_BARG=max(exp(par_likelihood[189]),0.0001);
    
    
    //5.2. Hbargg2b
    double M_A2_BARG=par_likelihood[190];
    double M_V2_BARG=max(exp(par_likelihood[191]),0.0001);
    
    
    //5.3. Hbargg2c
    double M_A3_BARG=par_likelihood[192];
    double M_V3_BARG=max(exp(par_likelihood[193]),0.0001);
    
    
    
    
    //5.4. Hbargg2d
    double M_A4_BARG=par_likelihood[194];
    double M_V4_BARG=max(exp(par_likelihood[195]),0.0001);
    
    
    //5.5. Hbargg2e
    double M_A5_BARG=par_likelihood[196];
    double M_V5_BARG=max(exp(par_likelihood[197]),0.0001);
    
    
    //5.6. Hbargg2f
    double M_A6_BARG=par_likelihood[198];
    double M_V6_BARG=max(exp(par_likelihood[199]),0.0001);
    
    
    //5.7. Hbargg2g
    double M_A7_BARG=par_likelihood[200];
    double M_V7_BARG=max(exp(par_likelihood[201]),0.0001);
    
    
    //5.8. Hbargg2h
    double M_A8_BARG=par_likelihood[202];
    double M_V8_BARG=max(exp(par_likelihood[203]),0.0001);
    
    
    //5.9. Hbargg2i
    double M_A9_BARG=par_likelihood[204];
    double M_V9_BARG=max(exp(par_likelihood[205]),0.0001);
    
    
    //5.10. Hbargg2j
    double M_A10_BARG=par_likelihood[206];
    double M_V10_BARG=max(exp(par_likelihood[207]),0.0001);
    
    
    //5.11. Hwomdecide12
    double M_A11_BARG=par_likelihood[208];
    double M_V11_BARG=max(exp(par_likelihood[209]),0.0001);
    
    
    //5.12. Hfathdecides12
    double M_A12_BARG=par_likelihood[210];
    double M_V12_BARG=max(exp(par_likelihood[211]),0.0001);
    
    
    //5.13. Hbothdecide
    double M_A13_BARG=par_likelihood[212];
    double M_V13_BARG=max(exp(par_likelihood[213]),0.0001);
    
    
    //5.14. Hcaresacwom
    double M_A14_BARG=par_likelihood[214];
    double M_V14_BARG=max(exp(par_likelihood[215]),0.0001);
    
    
    //5.15. Hcaresacman
    double M_A15_BARG=par_likelihood[216];
    double M_V15_BARG=max(exp(par_likelihood[217]),0.0001);
    
    
    //5.16. Hopofman121
    double M_A16_BARG=par_likelihood[218];
    double M_V16_BARG=max(exp(par_likelihood[219]),0.0001);
    
    
    //5.17. Hopofman122
    double M_A17_BARG=par_likelihood[220];
    double M_V17_BARG=max(exp(par_likelihood[221]),0.0001);
    
    
    //5.18. Hopofman123
    double M_A18_BARG=par_likelihood[222];
    double M_V18_BARG=max(exp(par_likelihood[223]),0.0001);
    
    
    //5.19. Hopofman124
    double M_A19_BARG=par_likelihood[224];
    double M_V19_BARG=max(exp(par_likelihood[225]),0.0001);
    
    
    
    //6. Investment in 2010 measures
    //6.1 Cinvh3
    double M_A1_INV10=par_likelihood[226];
    double M_V1_INV10=max(exp(par_likelihood[227]),0.0001);
    
    //6.2 Cinvh4
    double M_A2_INV10=par_likelihood[228];
    double M_V2_INV10=max(exp(par_likelihood[229]),0.0001);
    
    //6.3 Cinv5
    double M_A3_INV10=par_likelihood[230];
    double M_V3_INV10=max(exp(par_likelihood[231]),0.0001);
    
    //6.4 Cinvh6
    double M_A4_INV10=par_likelihood[232];
    double M_V4_INV10=max(exp(par_likelihood[233]),0.0001);
    
    //6.5 Cinvh7
    double M_A5_INV10=par_likelihood[234];
    double M_V5_INV10=max(exp(par_likelihood[235]),0.0001);
    
    //6.6 Cinvh8
    double M_A6_INV10=par_likelihood[236];
    double M_V6_INV10=max(exp(par_likelihood[237]),0.0001);
    
    //6.7 Cinvh9
    double M_A7_INV10=par_likelihood[238];
    double M_V7_INV10=max(exp(par_likelihood[239]),0.0001);
    
    //6.8 Cinvh13
    double M_A8_INV10=par_likelihood[240];
    double M_V8_INV10=max(exp(par_likelihood[241]),0.0001);
    
    
    
    
    
    //7. Investments in 2012
    //7.1 Cinvf11a
    double M_A1_INV12=par_likelihood[242];
    double M_V1_INV12=max(exp(par_likelihood[243]),0.0001);
    
    //7.2 Cinvf11b
    double M_A2_INV12=par_likelihood[244];
    double M_V2_INV12=max(exp(par_likelihood[245]),0.0001);
    
    //7.3 Cinvf11c
    double M_A3_INV12=par_likelihood[246];
    double M_V3_INV12=max(exp(par_likelihood[247]),0.0001);
    
    //7.4 Cinvf11d
    double M_A4_INV12=par_likelihood[248];
    double M_V4_INV12=max(exp(par_likelihood[249]),0.0001);
    
    //7.5 Cinvf11e
    double M_A5_INV12=par_likelihood[250];
    double M_V5_INV12=max(exp(par_likelihood[251]),0.0001);
    
    //7.6 Cinvf11f
    double M_A6_INV12=par_likelihood[252];
    double M_V6_INV12=max(exp(par_likelihood[253]),0.0001);
    
    //7.7 Cinvf11g
    double M_A7_INV12=par_likelihood[254];
    double M_V7_INV12=max(exp(par_likelihood[255]),0.0001);
    
    //7.8 Cinvf11h
    double M_A8_INV12=par_likelihood[256];
    double M_V8_INV12=max(exp(par_likelihood[257]),0.0001);
    
    //7.9 Cinvf11i
    double M_A9_INV12=par_likelihood[258];
    double M_V9_INV12=max(exp(par_likelihood[259]),0.0001);
    
    //7.10 Cinvf11j
    double M_A10_INV12=par_likelihood[260];
    double M_V10_INV12=max(exp(par_likelihood[261]),0.0001);
    
    //7.11 Cinvf11k
    double M_A11_INV12=par_likelihood[262];
    double M_V11_INV12=max(exp(par_likelihood[263]),0.0001);
    
    //7.12 Cinvhm2_10
    double M_A12_INV12=par_likelihood[264];
    double M_V12_INV12=max(exp(par_likelihood[265]),0.0001);
    
    //7.13 Cinvhm2_11
    double M_A13_INV12=par_likelihood[266];
    double M_V13_INV12=max(exp(par_likelihood[267]),0.0001);
    
    //7.14 Cinvhm2_12
    double M_A14_INV12=par_likelihood[268];
    double M_V14_INV12=max(exp(par_likelihood[269]),0.0001);
    
    //7.15 Cinvhm2_13
    double M_A15_INV12=par_likelihood[270];
    double M_V15_INV12=max(exp(par_likelihood[271]),0.0001);
    
    //7.16 Cinvhm2_14
    double M_A16_INV12=par_likelihood[272];
    double M_V16_INV12=max(exp(par_likelihood[273]),0.0001);
    
    //7.17 Cinvhm2_15
    double M_A17_INV12=par_likelihood[274];
    double M_V17_INV12=max(exp(par_likelihood[275]),0.0001);
    
    //7.18 Cinvhm2_16
    double M_A18_INV12=par_likelihood[276];
    double M_V18_INV12=max(exp(par_likelihood[277]),0.0001);
    
    //7.19 Cinvhm2_18
    double M_A19_INV12=par_likelihood[278];
    double M_V19_INV12=max(exp(par_likelihood[279]),0.0001);
    
    //7.20 Csharesbedroomhowmany12
    double M_A20_INV12=par_likelihood[280];
    double M_V20_INV12=max(exp(par_likelihood[281]),0.0001);
    
    //7.21 Csharesbedhowmany12
    double M_A21_INV12=par_likelihood[282];
    double M_V21_INV12=max(exp(par_likelihood[283]),0.0001);
    
    
    //8. PARENT's effort in 2010
    //8.1 g42_a
    double M_A1_EFF10=par_likelihood[284];
    double M_V1_EFF10=max(exp(par_likelihood[285]),0.0001);
    
    //8.2 g42_b
    double M_A2_EFF10=par_likelihood[286];
    double M_V2_EFF10=max(exp(par_likelihood[287]),0.0001);
    
    //8.3 g42_c
    double M_A3_EFF10=par_likelihood[288];
    double M_V3_EFF10=max(exp(par_likelihood[289]),0.0001);
    
    //8.4 g42_d
    double M_A4_EFF10=par_likelihood[290];
    double M_V4_EFF10=max(exp(par_likelihood[291]),0.0001);
    
    //8.5 g42_e
    double M_A5_EFF10=par_likelihood[292];
    double M_V5_EFF10=max(exp(par_likelihood[293]),0.0001);
    
    //8.6 g42_f
    double M_A6_EFF10=par_likelihood[294];
    double M_V6_EFF10=max(exp(par_likelihood[295]),0.0001);
    
    
    //10. Effort of parents in 2012
    //10.1 f21a_p_t
    double M_A1_EFF12=par_likelihood[296];
    double M_V1_EFF12=max(exp(par_likelihood[297]),0.0001);
    
    //10.2 f21b_p_t
    double M_A2_EFF12=par_likelihood[298];
    double M_V2_EFF12=max(exp(par_likelihood[299]),0.0001);
    
    //10.3 f21c_p_t
    double M_A3_EFF12=par_likelihood[300];
    double M_V3_EFF12=max(exp(par_likelihood[301]),0.0001);
    
    //10.4 f21d_p_t
    double M_A4_EFF12=par_likelihood[302];
    double M_V4_EFF12=max(exp(par_likelihood[303]),0.0001);
    
    //10.5 f21e_p_t
    double M_A5_EFF12=par_likelihood[304];
    double M_V5_EFF12=max(exp(par_likelihood[305]),0.0001);
    
    //10.6 f21f_p_t
    double M_A6_EFF12=par_likelihood[306];
    double M_V6_EFF12=max(exp(par_likelihood[307]),0.0001);
    
    //10.7 f21g_p_t
    double M_A7_EFF12=par_likelihood[308];
    double M_V7_EFF12=max(exp(par_likelihood[309]),0.0001);
    
    //10.8 f21h_p_t
    double M_A8_EFF12=par_likelihood[310];
    double M_V8_EFF12=max(exp(par_likelihood[311]),0.0001);
    
    //10.9 f21i_p_t
    double M_A9_EFF12=par_likelihood[312];
    double M_V9_EFF12=max(exp(par_likelihood[313]),0.0001);
    
    //10.10 f21j_p_t
    double M_A10_EFF12=par_likelihood[314];
    double M_V10_EFF12=max(exp(par_likelihood[315]),0.0001);
    
    //10.11 f21k_p_t
    double M_A11_EFF12=par_likelihood[316];
    double M_V11_EFF12=max(exp(par_likelihood[317]),0.0001);
    
    //10.12 f21l_p_t
    double M_A12_EFF12=par_likelihood[318];
    double M_V12_EFF12=max(exp(par_likelihood[319]),0.0001);
    
    //10.13 f21m_p_t
    double M_A13_EFF12=par_likelihood[320];
    double M_V13_EFF12=max(exp(par_likelihood[321]),0.0001);
    
    //10.14 f21n_p_t
    double M_A14_EFF12=par_likelihood[322];
    double M_V14_EFF12=max(exp(par_likelihood[323]),0.0001);
    
    
    //Measurement systems for the bootstrap filter
    
    
    //==============================
    
    //1. Initializing Likelihood
    double loglik=0;
    double logcontribobs_ii=0; //Contribution of observation ii
    double logcontribobs_rr=0; //Contribution of simulation rr
    double logcontribobs_rr10=0;
    double logcontribdecision=0; // Sum of contributions of decisions correct
    double logcontribdecision10=0;
    double logcontribdecision_rr=0;// Contribution of decision in simulation rr
    double logcontribdecision_rr10=0;
    double likwageM=0;
    double likwageM10=0;
    double likwageF=0;
    double likwageF10=0;
    double likbehavior=0;
    double likbehavior10=0;
    double loglikskills=0;
    double loglikskills10=0;
    double loglikskills_rr=0;
    double loglikskills_rr10=0;
    double loglikeffm=0;
    double loglikeffm10=0;
    double loglikeffm_rr=0;
    double loglikeffm_rr10=0;
    double loglikefff=0;
    double loglikefff_rr=0;
    double loglikefff10=0;
    double loglikefff_rr10=0;
    double loglikmmu=0;
    double loglikmmu_rr=0;
    double loglikeinv=0;
    double loglikeinv_rr=0;
    double loglikeinv10=0;
    double loglikeinv_rr10=0;
    
    //2. All the initializations that take place in for and while loop are done here
    //NW->No work; W-> Work
    int rr=0; //Initializing the number of times for the integration
    //2.1 We need to initialize all factor variables in the possible alternatives.
    //work-no work, childcare-no childcare, etc.
    // Effort of mother
    
    double Meffort_factor12shock=0; //Shock of effort for mother
    double MeffortMean12NW=0; //Mean of effort of mother initializer if single
    double MeffortMean12NWNW_TOG=0; //Mean of effort of mother initializer if nw with spouse
    double MeffortMean12WNW_TOG=0; //Mean of effort if mother does't work and with spouse
    double MeffortMean12NWW_TOG=0; //Mean of effort if mother does't work and with spouse
    double MeffortMean12WW_TOG=0; //Mean of effort if mother does't work and with spouse
    double MeffortMean12W=0;  //Mean of effort if mother works
    double Meffort_factor12NW=0; // Effort level of mother initializer
    double Meffort_factor12W=0; // Effort level of mother initializer
    double Meffort_factor10shock=0; //Shock of effort for mother
    double MeffortMean10NW=0; //Mean of effort of mother initializer if single
    double MeffortMean10NWNW_TOG=0; //Mean of effort of mother initializer if nw with spouse
    double MeffortMean10WNW_TOG=0; //Mean of effort if mother does't work and with spouse
    double MeffortMean10NWW_TOG=0; //Mean of effort if mother does't work and with spouse
    double MeffortMean10WW_TOG=0; //Mean of effort if mother does't work and with spouse
    double MeffortMean10W=0;  //Mean of effort if mother works
    double Meffort_factor10NW=0; // Effort level of mother initializer
    double Meffort_factor10W=0; // Effort level of mother initializer
    double MEFFORTNWNW10NOSHOCK=0;
    double MEFFORTNWW10NOSHOCK=0;
    double MEFFORTWNW10NOSHOCK=0;
    double MEFFORTWW10NOSHOCK=0;
    
    double MEFFORTNWNW12NOSHOCK=0;
    double MEFFORTNWW12NOSHOCK=0;
    double MEFFORTWNW12NOSHOCK=0;
    double MEFFORTWW12NOSHOCK=0;
    // Effort of father
    double Feffort_factor12shock=0; //Shock of effort for father
    
    double FeffortMean12NWNW_TOG=0; // mean of effort if father nw and with spouse
    double FeffortMean12WNW_TOG=0; // mean of effort if father nw and with spouse
    double FeffortMean12NWW_TOG=0; // mean of effort if father nw and with spouse
    double FeffortMean12WW_TOG=0; // mean of effort if father nw and with spouse
    
    double Feffort_factor10shock=0; //Shock of effort for father
    
    double FeffortMean10NWNW_TOG=0; // mean of effort if father nw and with spouse
    double FeffortMean10WNW_TOG=0; // mean of effort if father nw and with spouse
    double FeffortMean10NWW_TOG=0; // mean of effort if father nw and with spouse
    double FeffortMean10WW_TOG=0; // mean of effort if father nw and with spouse
    double FEFFORTNWNW10NOSHOCK=0;
    double FEFFORTNWW10NOSHOCK=0;
    double FEFFORTWNW10NOSHOCK=0;
    double FEFFORTWW10NOSHOCK=0;
    
    double FEFFORTNWNW12NOSHOCK=0;
    double FEFFORTNWW12NOSHOCK=0;
    double FEFFORTWNW12NOSHOCK=0;
    double FEFFORTWW12NOSHOCK=0;
    
    // Skills
    double CSkills_factor12NW=0; //Skills level initializer
    double CSkills_factor12W=0; //Skills level initializer
    double CSkills_factor12NWNW_TOG=0; //Skills level initializer if tog and nwnw
    double CSkills_factor12NWW_TOG=0; //Skills level initializer if tog and nww
    double CSkills_factor12WNW_TOG=0; //Skills level initializer if tog and wnw
    double CSkills_factor12WW_TOG=0; //Skills level initializer if tog and nwnw
    double CSkills_factor12shock=0;  //Skills shock
    double CSkillsMean12NW=0; // Mean of skills level initializer
    
    double CSkillsMean12W=0; // Mean of skills level initializer
    
    double CSkills_factor10NWA=0; //Skills level initializer
    double CSkills_factor10NWNA=0; //Skills level initializer
    double CSkills_factor10WA=0; //Skills level initializer
    double CSkills_factor10WNA=0; //Skills level initializer
    double CSkills_factor10NWNWA_TOG=0; //Skills level initializer if tog and nwnw
    double CSkills_factor10NWNWNA_TOG=0; //Skills level initializer if tog and nwnw
    double CSkills_factor10NWWA_TOG=0; //Skills level initializer if tog and nww
    double CSkills_factor10NWWNA_TOG=0; //Skills level initializer if tog and nww
    double CSkills_factor10WNWA_TOG=0; //Skills level initializer if tog and wnw
    double CSkills_factor10WNWNA_TOG=0; //Skills level initializer if tog and wnw
    double CSkills_factor10WWA_TOG=0; //Skills level initializer if tog and nwnw
    double CSkills_factor10WWNA_TOG=0; //Skills level initializer if tog and nwnw
    double CSkills_factor10shock=0;  //Skills shock
    double CSkillsMean10NWA=0; // Mean of skills level initializer
    double CSkillsMean10NWNA=0; // Mean of skills level initializer
    double CSkillsMean10WA=0; // Mean of skills level initializer
    double CSkillsMean10WNA=0; // Mean of skills level initializer
    
    //Investment
    double CInvestmentMean12NW=0;  // Investment level initializer
    double CInvestmentMean12W=0;  // Investment level initializer
    double CInvestmentMean12NWNW_TOG=0; //Investment nwnw tog
    double CInvestmentMean12WNW_TOG=0; //Investment nw tog
    double CInvestmentMean12NWW_TOG=0; //Investment nww tog
    double CInvestmentMean12WW_TOG=0; //Investment ww tog
    double CInv_factor12shock=0; // Shock of investment
    double CInv_factor12NW=0;//Mean of investment level
    double CInv_factor12W=0;//Mean of investment level
    
    double CInvestmentMean10NWA=0;  // Investment level initializer
    double CInvestmentMean10NWNA=0;  // Investment level initializer
    double CInvestmentMean10WA=0;  // Investment level initializer
    double CInvestmentMean10WNA=0;  // Investment level initializer
    double CInvestmentMean10NWNWA_TOG=0; //Investment nwnw tog
    double CInvestmentMean10NWNWNA_TOG=0; //Investment nwnw tog
    double CInvestmentMean10WNWA_TOG=0; //Investment nw tog
    double CInvestmentMean10WNWNA_TOG=0; //Investment nw tog
    double CInvestmentMean10NWWA_TOG=0; //Investment nww tog
    double CInvestmentMean10NWWNA_TOG=0; //Investment nww tog
    double CInvestmentMean10WWA_TOG=0; //Investment ww tog
    double CInvestmentMean10WWNA_TOG=0; //Investment ww tog
    double CInv_factor10shock=0; // Shock of investment
    double CInv_factor10NWA=0;//Mean of investment level
    double CInv_factor10NWNA=0;//Mean of investment level
    double CInv_factor10WA=0;//Mean of investment level
    double CInv_factor10WNA=0;//Mean of investment level
    double CINVFACTOR10NWNWNA_noshock=0;
    double CINVFACTOR10NWNWA_noshock=0;
    double CINVFACTOR10WNWNA_noshock=0;
    double CINVFACTOR10WNWA_noshock=0;
    double CINVFACTOR10NWWNA_noshock=0;
    double CINVFACTOR10NWWA_noshock=0;
    double CINVFACTOR10WWNA_noshock=0;
    double CINVFACTOR10WWA_noshock=0;
    
    double CINVFACTOR12NWNW_noshock=0;
    double CINVFACTOR12WNW_noshock=0;
    double CINVFACTOR12NWW_noshock=0;
    double CINVFACTOR12WW_noshock=0;
    
    //Consumption
    //Mother
    double MconsumptionNW=0; //nw single
    double MconsumptionW; //w single
    double MconsumptionNWNW_TOG=0; //nwnw tog
    double MconsumptionWNW_TOG=0; //wnw tog
    double MconsumptionNWW_TOG=0; //nww tog
    double MconsumptionWW_TOG=0; //ww tog
    
    
    double MconsumptionNWA10=0; //nw single
    double MconsumptionNWNA10=0; //nw single
    double MconsumptionWA10=0; //w single
    double MconsumptionWNA10=0; //w single
    double MconsumptionNWNWA_TOG10=0; //nwnw tog
    double MconsumptionNWNWNA_TOG10=0; //nwnw tog
    double MconsumptionWNWNA_TOG10=0; //wnw tog
    double MconsumptionWNWA_TOG10=0; //wnw tog
    double MconsumptionNWWNA_TOG10=0; //nww tog
    double MconsumptionNWWA_TOG10=0; //nww tog
    double MconsumptionWWA_TOG10=0; //ww tog
    double MconsumptionWWNA_TOG10=0; //ww tog
    
    //Father
    
    double FconsumptionNWNW_TOG=0; //nwnw tog
    
    double FconsumptionWNW_TOG=0; //wnw tog
    double FconsumptionNWW_TOG=0; //nww tog
    double FconsumptionNWWNA_TOG=0; //nww tog
    double FconsumptionWW_TOG=0; //ww tog
    
    
    
    double FconsumptionNWNWA_TOG10=0; //nwnw tog
    double FconsumptionNWNWNA_TOG10=0; //nwnw tog
    double FconsumptionWNWA_TOG10=0; //wnw tog
    double FconsumptionWNWNA_TOG10=0; //wnw tog
    double FconsumptionNWWA_TOG10=0; //nww tog
    double FconsumptionNWWNA_TOG10=0; //nww tog
    double FconsumptionWWA_TOG10=0; //ww tog
    double FconsumptionWWNA_TOG10=0; //ww tog
    
    
    //Preference shocks initializer
    //Mother
    double MprefshockW=0; //w single
    double MprefshockNW=0; //nw single
    double MprefshockWA10=0; //w single
    double MprefshockWNA10=0; //w single
    double MprefshockNWA10=0; //nw single
    double MprefshockNWNA10=0; //nw single
    
    //Father
    double FprefshockW=0; //w single
    double FprefshockNW=0; //nw single
    double FprefshockWA10=0; //w single
    double FprefshockWNA10=0; //w single
    double FprefshockNWA10=0; //nw single
    double FprefshockNWNA10=0; //nw single
    
    
    //Utility levels
    //Mother
    double MutilityNW=0; //Utility for mother if she works sing
    double MutilityW=0; //Utility for mother if she doesn't work single
    double MutilityNWNW_TOG=0; //Utility nwnw tog
    double MutilityWNW_TOG=0; //Utility wnw tog
    double MutilityNWW_TOG=0; //Utility nww tog
    double MutilityWW_TOG=0; //Utility ww tog
    
    double MutilityNWA10=0; //Utility for mother if she works sing
    double MutilityNWNA10=0; //Utility for mother if she works sing
    double MutilityWA10=0; //Utility for mother if she doesn't work single
    double MutilityWNA10=0; //Utility for mother if she doesn't work single
    double MutilityNWNWA_TOG10=0; //Utility nwnw tog
    double MutilityNWNWNA_TOG10=0; //Utility nwnw tog
    double MutilityWNWA_TOG10=0; //Utility wnw tog
    double MutilityWNWNA_TOG10=0; //Utility wnw tog
    double MutilityNWWA_TOG10=0; //Utility nww tog
    double MutilityNWWNA_TOG10=0; //Utility nww tog
    double MutilityWWA_TOG10=0; //Utility ww tog
    double MutilityWWNA_TOG10=0; //Utility ww tog
    //Father
    
    double FutilityNWNW_TOG=0; //Utility nwnw tog
    double FutilityWNW_TOG=0; //Utility wnw tog
    double FutilityNWW_TOG=0; //Utility nww tog
    double FutilityWW_TOG=0; //Utility ww tog
    
    
    double FutilityNWNWA_TOG10=0; //Utility nwnw tog
    double FutilityNWNWNA_TOG10=0; //Utility nwnw tog
    double FutilityWNWA_TOG10=0; //Utility wnw tog
    double FutilityWNWNA_TOG10=0; //Utility wnw tog
    double FutilityNWWA_TOG10=0; //Utility nww tog
    double FutilityNWWNA_TOG10=0; //Utility nww tog
    double FutilityWWA_TOG10=0; //Utility ww tog
    double FutilityWWNA_TOG10=0; //Utility ww tog
    
    //Welfare of both
    double WelfareNWNW=0;
    double WelfareWNW=0;
    double WelfareNWW=0;
    double WelfareWW=0;
    
    double WelfareNWNWA10=0;
    double WelfareNWNWNA10=0;
    double WelfareWNWA10=0;
    double WelfareWNWNA10=0;
    double WelfareNWWA10=0;
    double WelfareNWWNA10=0;
    double WelfareWWA10=0;
    double WelfareWWNA10=0;
    
    //Indicator for the decision taken
    double Decision=0; //Predicted decision
    double obDecision=0; //Observed decision
    
    double Decision10=0; //Predicted decision
    double obDecision10=0; //Observed decision
    
    //Wage initializer
    double Mwageinlik=0; //Wage of mother initializer
    double Mwageinlik10=0; //Wage of mother in 2010 initializer
    double Fwageinlik=0; //Wage of mother initializer
    double Fwageinlik10=0;//Wage of father in 2010 initializer
    
    //Bargaining power
    double MMushock=0; //Shock for unobserved heterogeneity of Pareto weight
    double MMushock10=0; //Shock for unobserved heterogeneity of Pareto weight
    
    
    //Final decisions used for the smoothing distribution
    double Meffort12Final=0;
    double Feffort12Final=0;
    double Investment12Final=0;
    
    double Meffort10Final=0;
    double Feffort10Final=0;
    double Investment10Final=0;
    
    double AuxiliarV=0; //AuxiliarV is a variable used for intermediate steps
    double AuxiliarV1=0;
    double AuxiliarV2=0;
    
    //Price of childcare initializer
    double pricechildcare=0;
    
    
    //Set number of simulations to perform
    double RR=100;
    
    //Particles definitions
    //0 Skills at birth
    double S0rr=0; //Skills at birth
    
    vector<double> S0Vector; //Vector storking skills at birth
    S0Vector.resize(RR);
    
    vector<double> S0Vectoraux; //Vector storking skills at t=1 auxiliar
    S0Vectoraux.resize(RR);
    
    
    double loglik_S0rr=0; //Likelihood contribution of skills rr
    double loglik_S0=0; //Total likelihood contribution of sum of rr
    
    //1. Skills of primarycaregiver
    double PGrr=0; //Skills of PG
    vector<double> PGVector; //Vector storking skills PG
    PGVector.resize(RR);
    double loglik_PGrr=0; //Likelihood contribution of particle rr
    double loglik_PG=0; //likelihood contribution of the sum of all particles.
    
    
    //2. Skills at t=1
    vector<double> S1Vector; //Vector storking skills at t=1
    S1Vector.resize(RR);
    double loglik_S1rr=0; //Likelihood contribution of skills rr
    
    vector<double> S1Vectoraux; //Vector storking skills at t=1 auxiliar
    S1Vectoraux.resize(RR);
    
    vector<double> AUXVECTOR; //Vector storking skills at t=1 auxiliar
    AUXVECTOR.resize(RR);
    
    //2. Skills at t=2
    vector<double> S2Vector; //Vector storking skills at t=1
    S2Vector.resize(RR);
    double loglik_S2rr=0; //Likelihood contribution of skills rr
    
    vector<double> S2Vectoraux; //Vector storking skills at t=1
    S2Vectoraux.resize(RR);
    
    //Defining the remaining parameters for the bootstrap filter
    double wsum0=0;
    vector<double> w0_rr_vectorAUX; //Vector storing weights
    w0_rr_vectorAUX.resize(RR);
    
    vector<double> w0_rr_vector; //Vector storing weights
    w0_rr_vector.resize(RR);
    
    vector<double> w1_rr_vector; //Vector storing weights
    w1_rr_vector.resize(RR);
    
    vector<double> wtilde0_rr_vector;
    wtilde0_rr_vector.resize(RR);
    
    
    vector<double> w1_rr_vectorAUX; //Vector storing weights
    w1_rr_vectorAUX.resize(RR);
    
    
    vector<double> what1_rr_vector; //Vector storing weights
    what1_rr_vector.resize(RR);
    
    vector<double> wnormalized1_rr_vector; //Vector storing weights
    wnormalized1_rr_vector.resize(RR);
    
    vector<double> wtilde1_rr_vector; //Vector storing weights
    wtilde1_rr_vector.resize(RR); //Unnormalized weights. They will finally give the likelihood
    
    double wsum1=0;
    vector<double> wsum1_rr_vector; //Vector storing weights
    wsum1_rr_vector.resize(RR);
    
    vector<double> w2_rr_vector; //Vector storing weights
    w2_rr_vector.resize(RR);
    
    vector<double> w2_rr_vectorAUX; //Vector storing weights
    w2_rr_vectorAUX.resize(RR);
    
    vector<double> what2_rr_vector; //Vector storing weights
    what2_rr_vector.resize(RR);//Unnormalized weights. They will finally give the likelihood
    
    vector<double> wtilde2_rr_vector; //Vector storing weights
    wtilde2_rr_vector.resize(RR);
    
    
    double wsum2=0;
    vector<double> wsum2_rr_vector; //Vector storing weights
    wsum2_rr_vector.resize(RR);
    
    //Initializing the elements of the measurement system
    //1. Skills in 2010
    double SKMeasure1=0;
    double SKMeasure2=0;
    double SKMeasure3=0;
    double SKMeasure4=0;
    double SKMeasure5=0;
    double SKMeasure6=0;
    double SKMeasure7=0;
    double SKMeasure8=0;
    double SKMeasure9=0;
    double SKMeasure10=0;
    double SKMeasure11=0;
    
    //2. Skills in 2012
    double SK2012Measure1=0;
    double SK2012Measure2=0;
    double SK2012Measure3=0;
    double SK2012Measure4=0;
    double SK2012Measure5=0;
    double SK2012Measure6=0;
    double SK2012Measure7=0;
    double SK2012Measure8=0;
    double SK2012Measure9=0;
    double SK2012Measure10=0;
    double SK2012Measure11=0;
    double SK2012Measure12=0;
    double SK2012Measure13=0;
    
    //3. SKILLS AT BIRTH
    int SKBIRTHMEASURE1=0;
    int SKBIRTHMEASURE2=0;
    int SKBIRTHMEASURE3=0;
    int SKBIRTHMEASURE4=0;
    int SKBIRTHMEASURE5=0;
    int SKBIRTHMEASURE6=0;
    int SKBIRTHMEASURE7=0;
    int SKBIRTHMEASURE8=0;
    int SKBIRTHMEASURE9=0;
    int SKBIRTHMEASURE10=0;
    int SKBIRTHMEASURE11=0;
    int SKBIRTHMEASURE12=0;
    int SKBIRTHMEASURE13=0;
    int SKBIRTHMEASURE14=0;
    int SKBIRTHMEASURE15=0;
    int SKBIRTHMEASURE16=0;
    double SKBIRTHMEASURE17=0;
    double SKBIRTHMEASURE18=0;
    double SKBIRTHMEASURE19=0;
    double SKBIRTHMEASURE20=0;
    int SKBIRTHMEASURE21=0;
    double SKBIRTHMEASURE22=0;
    double SKBIRTHMEASURE23=0;
    
    //4. Skills OF PRIMARY CAREGIVER
    double SKPG1=0;
    double SKPG2=0;
    double SKPG3=0;
    double SKPG4=0;
    double SKPG5=0;
    double SKPG6=0;
    double SKPG7=0;
    double SKPG8=0;
    
    //5. Measures of bargaining
    double MBARG1=0;
    double MBARG2=0;
    double MBARG3=0;
    double MBARG4=0;
    double MBARG5=0;
    double MBARG6=0;
    double MBARG7=0;
    double MBARG8=0;
    double MBARG9=0;
    double MBARG10=0;
    int MBARG11=0;
    int MBARG12=0;
    int MBARG13=0;
    double MBARG14=0;
    double MBARG15=0;
    int MBARG16=0;
    int MBARG17=0;
    int MBARG18=0;
    int MBARG19=0;
    
    //6. Measures of investment in 2010
    int MINV1_10=0;
    int MINV2_10=0;
    int MINV3_10=0;
    int MINV4_10=0;
    int MINV5_10=0;
    int MINV6_10=0;
    int MINV7_10=0;
    int MINV8_10=0;
    
    //7. MEASURES of investment in 2012
    double MINV1_12=0;
    double MINV2_12=0;
    double MINV3_12=0;
    double MINV4_12=0;
    double MINV5_12=0;
    double MINV6_12=0;
    double MINV7_12=0;
    double MINV8_12=0;
    double MINV9_12=0;
    double MINV10_12=0;
    double MINV11_12=0;
    int MINV12_12=0;
    int MINV13_12=0;
    int MINV14_12=0;
    int MINV15_12=0;
    int MINV16_12=0;
    int MINV17_12=0;
    int MINV18_12=0;
    int MINV19_12=0;
    double MINV20_12=0;
    double MINV21_12=0;
    
    //8. MEASURES OF father's effor in 2010
    int FEFF1_10=0;
    int FEFF2_10=0;
    int FEFF3_10=0;
    int FEFF4_10=0;
    int FEFF5_10=0;
    int FEFF6_10=0;
    
    //9. MEASURES OF mother's effort in 2010
    int MEFF1_10=0;
    int MEFF2_10=0;
    int MEFF3_10=0;
    int MEFF4_10=0;
    int MEFF5_10=0;
    int MEFF6_10=0;
    
    //10. MEASURES OF FATHER'S EFFORT IN 2012
    double FEFF1_12=0;
    double FEFF2_12=0;
    double FEFF3_12=0;
    double FEFF4_12=0;
    double FEFF5_12=0;
    double FEFF6_12=0;
    double FEFF7_12=0;
    double FEFF8_12=0;
    double FEFF9_12=0;
    double FEFF10_12=0;
    double FEFF11_12=0;
    double FEFF12_12=0;
    double FEFF13_12=0;
    double FEFF14_12=0;
    
    //11. MEASURES OF MOTHER'S EFFORT IN 2012
    double MEFF1_12=0;
    double MEFF2_12=0;
    double MEFF3_12=0;
    double MEFF4_12=0;
    double MEFF5_12=0;
    double MEFF6_12=0;
    double MEFF7_12=0;
    double MEFF8_12=0;
    double MEFF9_12=0;
    double MEFF10_12=0;
    double MEFF11_12=0;
    double MEFF12_12=0;
    double MEFF13_12=0;
    double MEFF14_12=0;
    
    
    //Parameters used for the smoothing distribution
    
    //Defining the vectors w_{t|T} as defined in
    //"Fast particle smoothing:If I had a million particles"
    //First generating the elements that will store the bootstrap filter smoother
    
    vector<double>  W0MATRIX;
    W0MATRIX.resize(RR);
    
    vector<double>  W1MATRIX;
    W1MATRIX.resize(RR);
    
    vector<double>  W2MATRIX;
    W2MATRIX.resize(RR);
    
    vector<double>  S0MATRIX;
    S0MATRIX.resize(RR);
    
    vector<double>  S1MATRIX;
    S1MATRIX.resize(RR);
    
    vector<double>  S2MATRIX;
    S2MATRIX.resize(RR);
    
    
    vector<double> VECWEIGHTSSMOOTH0;
    VECWEIGHTSSMOOTH0.resize(RR);
    double weightsmooth0=0;
    double denomL0=0;
    double numL0=0;
    
    vector<double> VECWEIGHTSSMOOTH1;
    VECWEIGHTSSMOOTH1.resize(RR);
    double weightsmooth1=0;
    double denomL1=0;
    double numL1=0;
    
    vector<double> VECWEIGHTSSMOOTH2;
    VECWEIGHTSSMOOTH2.resize(RR);
    double weightsmooth2=0;
    
    
    
    int iterator=0; //Auxiliary iterator
    
    
    //And some other variables necessary
    double VMAX10=0; //Value of observed decision in 10
    double VMAX12=0; //Value of observed decision in 12
    
    double VOBS10=0; //VALUE of observed decision in 10
    double VOBS12=0; //VALUE of observed decision in 12
    
    double VSUM10=0; //Maximum of decission in 10
    double VSUM12=0; //Maximum of decisions in 12
    
    double VDIF10=0; //Difference between observed and maximum in 2010
    double VDIF12=0; //Difference between observed and maximum in 2012
    
    double TAUSMOOTHING=0.4; //Smoothing parameter of the logit simulator.
    
    double VCONTRIB10=0; //Contribution value of observed decision 2010 once smoothed.
    double VCONTRIB12=0; //Contribution value of observed decision 2012 once smoothed.
    
    double Decision10STORED=0; //Stored decision of 2010
    double Decision12STORED=0; //Stored decision of 2012
    
    double logcontribobdecision10SMOOTHED=0; //contribution of observed decision in 2010
    double logcontribobdecision12SMOOTHED=0; //contribution of observed decision in 2012
    
    
    vector<vector<long double> > OUTPUT;
    OUTPUT.resize(RR);
    for (int n=0; n<RR; n++){
        OUTPUT[n].resize(7);
    }
    
    
    //Setting seed to get draws.
    //This will generate the same random number every time.
    typedef boost::mt19937 RNGType;
    RNGType rng(SEED);
    
    boost::variate_generator<RNGType, boost::normal_distribution<> >
    generator(rng,
              boost::normal_distribution<>());
    
    //Boost to define uniform draws between 0 and 1
    boost::random::uniform_real_distribution<> unidistrib(0,1);
    
    //If I want different random numbers each time:
    //iboost::variate_generator<boost::mt19937, boost::normal_distribution<> >
    //    generator(boost::mt19937(time(0)),
    //        boost::normal_distribution<>());
    loglik=0;
    
    //Initializing the likelihood contributions for observation ii
    logcontribobs_rr=0; //Contribution of simulation rr
    likwageM=0;
    likwageM10=0;
    likwageF=0;
    likbehavior=0;
    likbehavior10=0;
    loglikskills=0;
    loglikskills10=0;
    loglikeffm=0;
    loglikeffm10=0;
    loglikefff=0;
    loglikefff10=0;
    loglikmmu=0;
    loglikeinv=0;
    loglikeinv10=0;
    logcontribdecision=0;
    logcontribdecision10=0;
    logcontribobdecision10SMOOTHED=0;
    logcontribobdecision12SMOOTHED=0;
    
    //Measurement system initiazlation.
    //Skills in 2010
    SKMeasure1=CSTDtepsi_pb_coo10;
    SKMeasure2=CSTDtepsi_pb_len10;
    SKMeasure3=CSTDtepsi_pb_mot10;
    SKMeasure4=CSTDtvip_pb10;
    SKMeasure5=CSTDcbcl1_pb_110;
    SKMeasure6=CSTDcbcl1_pb_210;
    SKMeasure7=CSTDcbcl1_pb_310;
    SKMeasure8=CSTDcbcl1_pb_410;
    SKMeasure9=CSTDcbcl1_pb_510;
    SKMeasure10=CSTDcbcl1_pb_610;
    SKMeasure11=CSTDcbcl1_pb_710;
    
    //Skills in 2012
    SK2012Measure1=CSTDtadi_pb_cog12;
    SK2012Measure2=CSTDtadi_pb_mot12;
    SK2012Measure3=CSTDtadi_pb_len12;
    SK2012Measure4=CSTDtadi_pb_se12;
    SK2012Measure5=CSTDbt_112;
    SK2012Measure6=CSTDbt_212;
    SK2012Measure7=CSTDbt_312;
    SK2012Measure8=CSTDbt_412;
    SK2012Measure9=CSTDbt_512;
    SK2012Measure10=CSTDbt_t12;
    SK2012Measure11=CSTDhtks_st12;
    SK2012Measure12=CSTDbdst_st12;
    SK2012Measure13=CSTDppvt_t12;
    
    
    //Skills at birth
    SKBIRTHMEASURE1=Ccondpregg3_1;
    SKBIRTHMEASURE2=Ccondpregg3_2;
    SKBIRTHMEASURE3=Ccondpregg3_3;
    SKBIRTHMEASURE4=Ccondpregg3_4;
    SKBIRTHMEASURE5=Ccondpregg3_5;
    SKBIRTHMEASURE6=Ccondpregg3_6;
    SKBIRTHMEASURE7=Ccondpregg3_7;
    SKBIRTHMEASURE8=Ccondpregg3_8;
    SKBIRTHMEASURE9=Ccondpregg3_9;
    SKBIRTHMEASURE10=Ccondpregg4a_1;
    SKBIRTHMEASURE11=Ccondpregg4a_2;
    SKBIRTHMEASURE12=Ccondpregg4a_3;
    SKBIRTHMEASURE13=Ccondpregg4a_4;
    SKBIRTHMEASURE14=Ccondpregg4a_5;
    SKBIRTHMEASURE15=Ccondpregg4a_6;
    SKBIRTHMEASURE16=Ccondpregg4a_7;
    SKBIRTHMEASURE17=Ccondpregg7b;
    SKBIRTHMEASURE18=Ccondpregg8b;
    SKBIRTHMEASURE19=Ccondpregg9;
    SKBIRTHMEASURE20=Ccondpregg11b;
    SKBIRTHMEASURE21=Ccondpregpreterm;
    SKBIRTHMEASURE22=Ccondpregg24;
    SKBIRTHMEASURE23=Ccondpregg23;
    
    //SKILLS OF PRIMARY CAREGIVER
    SKPG1=Cwais_pb_num;
    SKPG2=Cwais_pb_vo;
    SKPG3=Cbfi_pb_ama;
    SKPG4=Cbfi_pb_ape;
    SKPG5=Cbfi_pb_ext;
    SKPG6=Cbfi_pb_neu;
    SKPG7=Cbfi_pb_res;
    SKPG8=Cpsi_pb_total;
    
    //Measures of bargaining
    MBARG1=Hbargg2a;
    MBARG2=Hbargg2b;
    MBARG3=Hbargg2c;
    MBARG4=Hbargg2d;
    MBARG5=Hbargg2e;
    MBARG6=Hbargg2f;
    MBARG7=Hbargg2g;
    MBARG8=Hbargg2h;
    MBARG9=Hbargg2i;
    MBARG10=Hbargg2j;
    MBARG11=Hwomdecide12;
    MBARG12=Hfathdecides12;
    MBARG13=Hbothdecide;
    MBARG14=Hcaresacwom;
    MBARG15=Hcaresacman;
    MBARG16=Hopofman121;
    MBARG17=Hopofman122;
    MBARG18=Hopofman123;
    MBARG19=Hopofman124;
    
    //Investment in 2010
    MINV1_10=Cinvh3;
    MINV2_10=Cinvh4;
    MINV3_10=Cinvh5;
    MINV4_10=Cinvh6;
    MINV5_10=Cinvh7;
    MINV6_10=Cinvh8;
    MINV7_10=Cinvh9;
    MINV8_10=Cinvh13;
    
    //Investment in 2012
    MINV1_12=Cinvf11a;
    MINV2_12=Cinvf11b;
    MINV3_12=Cinvf11c;
    MINV4_12=Cinvf11d;
    MINV5_12=Cinvf11e;
    MINV6_12=Cinvf11f;
    MINV7_12=Cinvf11g;
    MINV8_12=Cinvf11h;
    MINV9_12=Cinvf11i;
    MINV10_12=Cinvf11j;
    MINV11_12=Cinvf11k;
    MINV12_12=Cinvhm2_10;
    MINV13_12=Cinvhm2_11;
    MINV14_12=Cinvhm2_12;
    MINV15_12=Cinvhm2_13;
    MINV16_12=Cinvhm2_14;
    MINV17_12=Cinvhm2_15;
    MINV18_12=Cinvhm2_16;
    MINV19_12=Cinvhm2_18;
    MINV20_12=Csharesbedroomhowmany12;
    MINV21_12=Csharesbedhowmany12;
    
    //Fathjer's effort in 2010:
    FEFF1_10=g42_a2;
    FEFF2_10=g42_b2;
    FEFF3_10=g42_c2;
    FEFF4_10=g42_d2;
    FEFF5_10=g42_e2;
    FEFF6_10=g42_f2;
    
    //Mothjer's effort in 2010:
    MEFF1_10=g42_a1;
    MEFF2_10=g42_b1;
    MEFF3_10=g42_c1;
    MEFF4_10=g42_d1;
    MEFF5_10=g42_e1;
    MEFF6_10=g42_f1;
    
    //FATHER'S EFFORT IN 2012
    FEFF1_12=f21a_p_t;
    FEFF2_12=f21b_p_t;
    FEFF3_12=f21c_p_t;
    FEFF4_12=f21d_p_t;
    FEFF5_12=f21e_p_t;
    FEFF6_12=f21f_p_t;
    FEFF7_12=f21g_p_t;
    FEFF8_12=f21h_p_t;
    FEFF9_12=f21i_p_t;
    FEFF10_12=f21j_p_t;
    FEFF11_12=f21k_p_t;
    FEFF12_12=f21l_p_t;
    FEFF13_12=f21m_p_t;
    FEFF14_12=f21n_p_t;
    
    //mother's effort in 2012
    MEFF1_12=f21a_p_t;
    MEFF2_12=f21b_p_t;
    MEFF3_12=f21c_p_t;
    MEFF4_12=f21d_p_t;
    MEFF5_12=f21e_p_t;
    MEFF6_12=f21f_p_t;
    MEFF7_12=f21g_p_t;
    MEFF8_12=f21h_p_t;
    MEFF9_12=f21i_p_t;
    MEFF10_12=f21j_p_t;
    MEFF11_12=f21k_p_t;
    MEFF12_12=f21l_p_t;
    MEFF13_12=f21m_p_t;
    MEFF14_12=f21n_p_t;
    
    
    
    //=============================
    //Definitions needed within the household:
    //=============================
    
    //Declaring the aalpha4
    aalpha4m=aalpha40m+aalpha41m*Hhchores12*1;
    aalpha4f=aalpha40f+aalpha41f*Hhchores12*1;
    aalpha4m10=aalpha40m10+aalpha41m10*Hhchores12*1;
    aalpha4f10=aalpha40f10+aalpha41f10*Hhchores12*1;
    //The price of childcare
    pricechildcare=pchildcare0+(1/(1+Hchildcareobs))*pchildcare1;
    
    //Price of investment
    priceINV=priceINV0-HDens*priceINV1;
    //priceINV=priceINV0-(1/(1+Hchildcareobs))*priceINV1;
    
    //Declaring aalpha3:
    
    
    //aalpha3m10=aalpha3m10-aalpha3m10_mtjh*MTJH+aalpha3mage10*Magegroup10;
    //aalpha3m=aalpha3m-aalpha3m12_mtjh*MTJH+aalpha3mage12*Magegroup12;
    //cout << aalpha3m10 << " aalpha3m10 post" << endl;
    
    
    //========================
    //1.1 Wages likelihood
    //========================
    //1.1.1 2010
    //=======================
    //1.1.1.1  For mothers
    if (Mfraclabor10==0){
        likwageM10=0;
        Mwageinlik10=F_predwage(bbeta0m,bbeta1m,bbeta2m,bbeta3m,
                                Myrschool12,Mage10);
    }
    if (Mfraclabor10==1){
        likwageM10=F_likelihood_wage(bbeta0m,bbeta1m,bbeta2m,bbeta3m,
                                     stdwm,Myrschool12,Mage10,Mwage10);
        Mwageinlik10=Mwage10;
    }
    
    //1.1.1.2 For fathers
    if (Ffraclabor10==0){
        likwageF10=0;
        Fwageinlik10=F_predwage(bbeta0f,bbeta1f,bbeta2f,bbeta3f,
                                Fyrschool12,Fage10);
    }
    if (Ffraclabor10==1){
        likwageF10=F_likelihood_wage(bbeta0f,bbeta1f,bbeta2f,bbeta3f,
                                     stdwf,Fyrschool12,Fage10,Fwage10);
        Fwageinlik10=Fwage10;
    }
    if (1==2){
        cout << likwageF10 << " likwageF10 " << endl;
    }
    //=======================
    //2012
    //=======================
    //1.1.1 For mother
    if (Mfraclabor12==0){
        likwageM=0;
        Mwageinlik=F_predwage(bbeta0m,bbeta1m,bbeta2m,bbeta3m,Myrschool12,Mage12);
    }
    else if (Mfraclabor12==1){
        likwageM=F_likelihood_wage(bbeta0m,bbeta1m,bbeta2m,bbeta3m,
                                   stdwm,Myrschool12,Mage12,Mwage12);
        Mwageinlik=Mwage12;
    }
    //1.1.2 For father
    
    if (Ffraclabor12==0){
        likwageF=0;
        Fwageinlik=F_predwage(bbeta0f,bbeta1f,bbeta2f,bbeta3f,Fyrschool12,Fage12);
    }
    else if (Ffraclabor12==1){
        likwageF=F_likelihood_wage(bbeta0f,bbeta1f,bbeta2f,bbeta3f,
                                   stdwf,Fyrschool12,Fage12,Fwage12);
        Fwageinlik=Fwage12;
    }
    
    
    //Likelihood that requires simulation
    rr=0;
    
    //Defining the mean of everything that can be set out of the loop for the
    //RR simulations
    
    //0. Single mothers 2010
    //========================================================
    
    //1.1 Effort levels
    //1.1.1 If doesn't work
    MeffortMean10NW=0;
    //1.1.2 If works
    MeffortMean10W=0;
    //1.2 Investment levels
    //1.2.1 If doesn't work
    CInvestmentMean10NWNA=F_investment10(aalpha1m10,aalpha2m,aalpha2m10,
                                         ttheta1,ttheta0,
                                         Mwageinlik10,Mnly10,0,priceINV);
    
    CInvestmentMean10NWA=F_investment10(aalpha1m10,aalpha2m,aalpha2m10,
                                        ttheta1,ttheta0,
                                        Mwageinlik10,Mnly10-pricechildcare,0,priceINV);
    
    
    //1.2.2 If works
    CInvestmentMean10WA=F_investment10(aalpha1m10,aalpha2m,aalpha2m10,
                                       ttheta1,ttheta0,
                                       Mwageinlik10,Mnly10,1,priceINV);
    CInvestmentMean10WNA=F_investment10(aalpha1m10,aalpha2m,aalpha2m10,
                                        ttheta1,ttheta0,
                                        Mwageinlik10,Mnly10-pricechildcare,1,priceINV);
    //1. Single mothers 2012
    //========================================================
    
    //1.1 Effort levels
    //1.1.1 If doesn't work
    MeffortMean12NW=0;
    //1.1.2 If works
    MeffortMean12W=0;
    //1.2 Investment levels
    //1.2.1 If doesn't work
    CInvestmentMean12NW=F_investment(aalpha1m,aalpha2m,ttheta1,
                                     Mwageinlik,Mnly12,0,priceINV);
    
    
    //1.2.2 If works
    CInvestmentMean12W=F_investment(aalpha1m,aalpha2m,ttheta1,
                                    Mwageinlik,Mnly12,1,priceINV);
    
    
    
    //======================================================
    //Identifying the observed decision in each year
    //======================================================
    
    //2010
    //2010 Single mothers
    if (Cliveswithfather10==0 & Cliveswithmother10==1){
        //Decisions:
        //1. Mother does not work and no childcare
        //2. Mother does not work and childcare
        //3. Mother works and no childcare
        //4. Mother works and childcare
        obDecision10=(1-Mfraclabor10)*(1-Cchildcare10)*1+
        (1-Mfraclabor10)*Cchildcare10*2+
        Mfraclabor10*(1-Cchildcare10)*3+
        Mfraclabor10*Cchildcare10;
    }
    
    //2010 Married couples
    if (Cliveswithfather10==1 & Cliveswithmother10==1){
        //Decisions:
        //1. Both work and childcare
        //2. Both work and no childcare
        //3. Father works, mother doesn't and childcare
        //4. Father works, mother doesn't and no childcare
        //5. Father doesn't work, mother works and childcare
        //6. Father doesn't work, mother works and no childcare
        //7. Father doesn't work, mother doesn't work and childcare
        //8. Father doesn't work, mother doesn't work and no childcare
        
        
        obDecision10=Ffraclabor10*Mfraclabor10*Cchildcare10*1+
        Ffraclabor10*Mfraclabor10*(1-Cchildcare10)*2+
        Ffraclabor10*(1-Mfraclabor10)*Cchildcare10*3+
        Ffraclabor10*(1-Mfraclabor10)*(1-Cchildcare10)*4+
        (1-Ffraclabor10)*Mfraclabor10*Cchildcare10*5+
        (1-Ffraclabor10)*Mfraclabor10*(1-Cchildcare10)*6+
        (1-Ffraclabor10)*(1-Mfraclabor10)*(Cchildcare10)*7+
        (1-Ffraclabor10)*(1-Mfraclabor10)*(1-Cchildcare10)*8;
        
        
    }
    
    //2012
    //2012 Single mothers
    if (Cliveswithfather12==0 & Cliveswithmother12==1){
        //Decisions:
        //1. Mother works
        //0. Mother does not work
        obDecision=(Mfraclabor12==1);
    }
    
    //2012 Married couples
    if (Cliveswithfather12==1 & Cliveswithmother12==1){
        //Decisions:
        //1. Both work
        //2. Father works, mother doesn't
        //3. Mother works, father doesn't
        //4. Neither father or mother work.
        
        obDecision=Ffraclabor12*Mfraclabor12*1+
        Ffraclabor12*(1-Mfraclabor12)*2+
        (1-Ffraclabor12)*Mfraclabor12*3+
        (1-Ffraclabor12)*(1-Mfraclabor12)*4;
    }
    //cout << " here1 "<< endl;
    
    
    
    //Start of simulations. I will start with the particles that need
    // to be drawn in t=0. These are the skills at period zero or at birth,
    //given by S0, and the skills of the primary caregiver. In first stance
    //I will assume just normal distributions for both. Analyzing the distribution of
    //both variables in STATA, the assumption of standard normal does not seem unreasonable.
    
    
    //==========================================
    //Draws of the birth outcomes and PG skills
    //==========================================
    
    rr=0;
    loglik_S0rr=0;
    loglik_PGrr=0;
    loglik_S1rr=0;
    loglik_S2rr=0;
    wsum0=0;
    wsum1=0;
    wsum2=0;
    
    for (int mm=0; mm<RR; mm=mm+1){
        //Skills at birth
        S0rr=gen_normal_3(generator);
        S0Vector[mm]=S0rr;
        
        //Skills of primary caregiver
        PGrr=gen_normal_3(generator);
        PGVector[mm]=PGrr;
        
        //Likelihood of skills at birth and skills of primary caregiver
        //Before it was loglik_S0rr +=....
        w0_rr_vector[mm]=MEASloglikSkillsBIRTH(SKBIRTHMEASURE1, M_A1_SBIRTH, M_VS1_SBIRTH,SKBIRTHMEASURE2, M_A2_SBIRTH, M_VS2_SBIRTH, SKBIRTHMEASURE3, M_A3_SBIRTH, M_VS3_SBIRTH, SKBIRTHMEASURE4, M_A4_SBIRTH, M_VS4_SBIRTH, SKBIRTHMEASURE5, M_A5_SBIRTH, M_VS5_SBIRTH, SKBIRTHMEASURE6, M_A6_SBIRTH, M_VS6_SBIRTH, SKBIRTHMEASURE7, M_A7_SBIRTH, M_VS7_SBIRTH, SKBIRTHMEASURE8, M_A8_SBIRTH, M_VS8_SBIRTH, SKBIRTHMEASURE9, M_A9_SBIRTH, M_VS9_SBIRTH, SKBIRTHMEASURE10, M_A10_SBIRTH, M_VS10_SBIRTH, SKBIRTHMEASURE11, M_A11_SBIRTH, M_VS11_SBIRTH, SKBIRTHMEASURE12, M_A12_SBIRTH, M_VS12_SBIRTH, SKBIRTHMEASURE13, M_A13_SBIRTH, M_VS13_SBIRTH, SKBIRTHMEASURE14, M_A14_SBIRTH, M_VS14_SBIRTH, SKBIRTHMEASURE15, M_A15_SBIRTH, M_VS15_SBIRTH, SKBIRTHMEASURE16, M_A16_SBIRTH, M_VS16_SBIRTH, SKBIRTHMEASURE17, M_A17_SBIRTH, M_VS17_SBIRTH, SKBIRTHMEASURE18, M_A18_SBIRTH, M_VS18_SBIRTH, SKBIRTHMEASURE19, M_A19_SBIRTH, M_VS19_SBIRTH,SKBIRTHMEASURE20, M_A20_SBIRTH, M_VS20_SBIRTH,SKBIRTHMEASURE21, M_A21_SBIRTH, M_VS21_SBIRTH,SKBIRTHMEASURE22, M_A22_SBIRTH, M_VS22_SBIRTH,SKBIRTHMEASURE23, M_A23_SBIRTH, M_VS23_SBIRTH , exp(S0rr));
        
        
        loglik_PGrr+=MEASloglikPG(SKPG1, M_A1_SPG, M_V1_SPG, SKPG2, M_A2_SPG, M_V2_SPG,SKPG3, M_A3_SPG, M_V3_SPG,SKPG4, M_A4_SPG, M_V4_SPG,SKPG5, M_A5_SPG, M_V5_SPG,SKPG6, M_A6_SPG, M_V6_SPG,SKPG7, M_A7_SPG, M_V7_SPG,SKPG8, M_A8_SPG, M_V8_SPG,exp(PGrr));
        
    }
    //Storing the weights that are NOT normalized and that will give the likelihood
    //For the period 0 they need not to be reshufled in order to compute the likelihood of particles in t=0.
    wtilde0_rr_vector=w0_rr_vector;
    
    
    //Normalizing weights so that we have that they all add up to one.
    AuxiliarV=0;
    for (int mm=0; mm<RR; mm=mm+1){
        AuxiliarV+=w0_rr_vector[mm];
    }
    for (int mm=0; mm<RR; mm=mm+1){
        w0_rr_vector[mm]=w0_rr_vector[mm]/AuxiliarV;
    }
    AuxiliarV=0;
    
    
    //And getting the draws of the particles according to the specified weights
    boost::random::discrete_distribution<int> distskills0(w0_rr_vector);
    iterator=0;
    wsum0=0;
    for (int mm=0; mm<RR;mm=mm+1){
        iterator=distskills0(rng);
        S0Vectoraux[mm]=S0Vector[iterator];
        w0_rr_vectorAUX[mm]=w0_rr_vector[iterator];
        wsum0 +=wtilde0_rr_vector[iterator];
        
    }
    
    //Updating the weights and the skills
    for (int mm=0; mm<RR;mm=mm+1){
        S0Vectoraux[mm]=exp(S0Vectoraux[mm]);
        S0Vector[mm]=exp(S0Vector[mm]);
    }
    
    //Adding lines to exp skills as mentioned
    //Redefining vector of skills
    S0Vector=S0Vectoraux;
    w0_rr_vector=w0_rr_vectorAUX;
    
    
    
    //========================================================
    //2010 simulations
    //========================================================
    rr=0;
    rr=0;
    loglik_S1rr=0;
    loglik_S2rr=0;
    wsum1=0;
    wsum2=0;
    
    
    while (rr<RR){
        logcontribobs_rr10=0;
        loglikskills_rr10=0;
        loglikeffm_rr10=0;
        loglikefff_rr10=0;
        loglikeinv_rr10=0;
        Decision10=0;
        logcontribdecision_rr10=0;
        
        //1. Draw shocks from the corresponding distributions
        //===========================================================
        //1.1 Single mothers:
        if (Cliveswithfather10==0 & Cliveswithmother10==1){
            //cout << " here10 "<< endl;
            //========================================================
            //1. Effort levels
            //========================================================
            //Draw a shock for the effort levels
            Meffort_factor10shock=gen_normal_3(generator);
            
            //========================================================
            //1.1 If mom does'nt work
            //========================================================
            // Getting a draw from a standard normal distribution
            
            // Now transforming it as a draw from the given distribution
            Meffort_factor10NW=Meffort_factor10shock*stdeffortmom;
            
            Meffort_factor10NW=exp(Meffort_factor10NW)*MeffortMean10NW;
            
            //Effort can't be negative, we are drawing from a truncated normal distribution
            if (Meffort_factor10NW<0){
                Meffort_factor10NW=0;
            }
            //cout << " here11 "<< endl;
            //We also bound the effort from above
            if (Meffort_factor10NW>=1.0e+6){
                Meffort_factor10NW=1.0e+6;
            }
            //========================================================
            //1.2 If mom works
            // Getting a draw from a standard normal distribution
            
            // Now transforming it as a draw from the given distribution
            
            Meffort_factor10W=Meffort_factor10shock*stdeffortmom;
            
            Meffort_factor10W=exp(Meffort_factor10W)*MeffortMean10W;
            //cout << " here12 "<< endl;
            //Effort can't be negative, we are drawing from a truncated normal distribution
            if (Meffort_factor10W<0){
                Meffort_factor10W=0;
            }
            if (Meffort_factor10W>=1.0e+6){
                Meffort_factor10W=1.0e+6;
            }
            //========================================================
            
            
            //========================================================
            //2. Investment levels
            //========================================================
            //Getting the draw for the investment shock
            CInv_factor10shock=gen_normal_3(generator);
            //========================================================
            //2.1 If mom doesn't work
            //========================================================
            //3. Now transforming it as a draw from the given distribution
            //CHILDCARE
            CInv_factor10NWA=CInv_factor10shock*stdinvestment;
            CInv_factor10NWA=CInv_factor10NWA+CInvestmentMean10NWA;
            //cout << " here13 "<< endl;
            //Investment can't be negative
            if (CInv_factor10NWA<0){
                CInv_factor10NWA=0;
            }
            
            //NO CHILDCARE
            CInv_factor10NWNA=CInv_factor10shock*stdinvestment;
            CInv_factor10NWNA=CInv_factor10NWNA+CInvestmentMean10NWNA;
            //cout << " here13 "<< endl;
            //Investment can't be negative
            if (CInv_factor10NWNA<0){
                CInv_factor10NWNA=0;
            }
            //========================================================
            //2.2 If mom works
            //========================================================
            //3. Now transforming it as a draw from the given distribution
            //CHILDCARE
            CInv_factor10WA=CInv_factor10shock*stdinvestment;
            CInv_factor10WA=CInv_factor10WA+CInvestmentMean10WA;
            //cout << " here14 "<< endl;
            //Investment can't be negative
            if (CInv_factor10WA<0){
                CInv_factor10WA=0;
            }
            
            
            //NO CHILDCARE
            CInv_factor10WNA=CInv_factor10shock*stdinvestment;
            CInv_factor10WNA=CInv_factor10WNA+CInvestmentMean10WNA;
            //cout << " here14 "<< endl;
            //Investment can't be negative
            if (CInv_factor10WNA<0){
                CInv_factor10WNA=0;
            }
            //==========================================================
            //3. Skills levels
            //==========================================================
            //Shock of skills
            CSkills_factor10shock=gen_normal_3(generator);
            //========================================================
            //1. If mom doesn't works
            //========================================================
            //3. Now transforming it as a draw from the given distribution
            //No childcare
            CSkillsMean10NWNA=F_predskills(ddelta0,ddelta1,
                                           ddelta2,ddelta3,ddelta4,
                                           Cedad_meses10,ttheta0,
                                           ttheta1,ttheta2, pphi,ggammaf, ggammam,
                                           Feffort10,Meffort_factor10NW,
                                           CInv_factor10NWNA,
                                           S0Vector[rr],
                                           0,
                                           PGVector[rr],Hmemberstotal10);
            
            //Notice that it is no longer
            //Cbirthfactor and PG as we
            //are doing the bootstrap filter.
            
            if (CSkills_factor10shock*stdskills<=500){
                CSkills_factor10NWNA=CSkills_factor10shock*stdskills;
            }
            if (CSkills_factor10shock*stdskills>=500){
                CSkills_factor10NWNA=500;
            }
            
            CSkills_factor10NWNA=exp(CSkills_factor10NWNA)*CSkillsMean10NWNA;
            //Childcare
            CSkillsMean10NWA=F_predskills(ddelta0,ddelta1,ddelta2,ddelta3,ddelta4,
                                          Cedad_meses10,ttheta0,
                                          ttheta1,ttheta2, pphi,ggammaf, ggammam,
                                          Feffort10,Meffort_factor10NW,
                                          CInv_factor10NWA,S0Vector[rr],1,PGVector[rr],Hmemberstotal10);
            //cout << " here15 "<< endl;
            if (CSkills_factor10shock*stdskills<=500){
                CSkills_factor10NWA=CSkills_factor10shock*stdskills;
            }
            if (CSkills_factor10shock*stdskills>=500){
                CSkills_factor10NWA=500;
            }
            
            CSkills_factor10NWA=exp(CSkills_factor10NWA)*CSkillsMean10NWA;
            //========================================================
            //2. If mom works
            //========================================================
            //3. Now transforming it as a draw from the given distribution
            //no childcare
            CSkillsMean10WNA=F_predskills(ddelta0,ddelta1,ddelta2,ddelta3,ddelta4,
                                          Cedad_meses10,ttheta0,
                                          ttheta1,ttheta2, pphi,ggammaf, ggammam,
                                          Feffort10,Meffort_factor10W,
                                          CInv_factor10WNA,S0Vector[rr],0,PGVector[rr],Hmemberstotal10);
            //cout << " here16 "<< endl;
            if (CSkills_factor10shock*stdskills<=500){
                CSkills_factor10WNA=CSkills_factor10shock*stdskills;
            }
            if (CSkills_factor10shock*stdskills>=500){
                CSkills_factor10WNA=500;
            }
            
            
            CSkills_factor10WNA=exp(CSkills_factor10WNA)*CSkillsMean10WNA;
            //childcare
            CSkillsMean10WA=F_predskills(ddelta0,ddelta1,ddelta2,ddelta3,ddelta4,
                                         Cedad_meses10,ttheta0,
                                         ttheta1,ttheta2, pphi,ggammaf, ggammam,
                                         Feffort10,Meffort_factor10W,
                                         CInv_factor10WA,S0Vector[rr],1,PGVector[rr],Hmemberstotal10);
            //cout << " here16 "<< endl;
            if (CSkills_factor10shock*stdskills<=500){
                CSkills_factor10WA=CSkills_factor10shock*stdskills;
            }
            if (CSkills_factor10shock*stdskills>=500){
                CSkills_factor10WA=500;
            }
            
            
            CSkills_factor10WA=exp(CSkills_factor10WA)*CSkillsMean10WA;
            //cout << " here17 "<< endl;
            //===========================================================
            //===========================================================
            //4 Preference shocks
            //===========================================================
            //1. If mom doesn't work
            //===========================================================
            //No childcare
            MprefshockNWNA10=gen_normal_3(generator);
            MprefshockNWNA10=MprefshockWNA10*MshockNWNA;
            //Childcare
            MprefshockNWA10=gen_normal_3(generator);
            MprefshockNWA10=MprefshockWA10*MshockNWA;
            //===========================================================
            //2. If mom works
            //===========================================================
            //No childcare
            MprefshockWNA10=gen_normal_3(generator);
            MprefshockWNA10=MprefshockWNA10*MshockWNA;
            //childcare
            MprefshockWA10=gen_normal_3(generator);
            MprefshockWA10=MprefshockWA10*MshockWA;
            //===========================================================
            //===========================================================
            //5 Consumption
            //===========================================================
            //1. If mom doesn't work
            //===========================================================
            //No childcare
            MconsumptionNWNA10=Mnly10-priceINV*CInv_factor10NWNA;
            if (MconsumptionNWNA10<0){
                MconsumptionNWNA10=1.0e-10;
            }
            //childcare
            MconsumptionNWA10=Mnly10-priceINV*CInv_factor10NWA-pricechildcare;
            if (MconsumptionNWA10<0){
                MconsumptionNWA10=1.0e-10;
            }
            //===========================================================
            //2. If mom works
            //===========================================================
            //no childcare
            MconsumptionWNA10=Mnly10-priceINV*CInv_factor10NWNA+Mwageinlik10;
            if (MconsumptionWNA10<0){
                MconsumptionWNA10=1.0e-10;
            }
            //childcare
            MconsumptionWA10=Mnly10-priceINV*CInv_factor10NWA+Mwageinlik10-pricechildcare;
            if (MconsumptionWA10<0){
                MconsumptionWA10=1.0e-10;
            }
            
            
            
            //cout << " here18 "<< endl;
            //===========================================================
            //===========================================================
            //6. Behavioral model
            //==========================================================
            //1. If mom doesn't work
            //==========================================================
            //no childcare
            MutilityNWNA10=F_utility10(aalpha1m10,aalpha2m10,aalpha3m10,
                                       aalpha4m10,aalpha5m10,
                                       MconsumptionNWNA10,
                                       Meffort_factor10NW,
                                       0,CSkills_factor10NWNA,0)+MprefshockNWNA10;
            //childcare
            MutilityNWA10=F_utility10(aalpha1m10,aalpha2m10,aalpha3m10,
                                      aalpha4m10,aalpha5m10,
                                      MconsumptionNWA10,
                                      Meffort_factor10NW,
                                      0,CSkills_factor10NWA,1)+MprefshockNWA10;
            //===========================================================
            //2. If mom works
            //===========================================================
            
            //no childcare
            MutilityWNA10=F_utility10(aalpha1m10,aalpha2m10,aalpha3m10,
                                      aalpha4m,aalpha5m10,
                                      MconsumptionWNA10,
                                      Meffort_factor10W,
                                      1,CSkills_factor10WNA,0)+MprefshockWNA10;
            //CHILDCARE
            MutilityWA10=F_utility10(aalpha1m10,aalpha2m10,aalpha3m10,
                                     aalpha4m,aalpha5m10,
                                     MconsumptionWA10,
                                     Meffort_factor10W,1,CSkills_factor10WA,
                                     1)+MprefshockWA10;
            //cout << " here19 "<< endl;
            //===========================================================
            //===========================================================
            //cout << " here20 "<< endl
            //Identifying the decision
            if ((MutilityNWNA10> MutilityNWA10) &&
                (MutilityNWNA10> MutilityWNA10) &&
                (MutilityNWNA10> MutilityWA10)){
                Decision10=1;
            }
            else if ((MutilityNWA10> MutilityWNA10) &&
                     (MutilityNWA10> MutilityWA10)){
                Decision10=2;
            }
            else if (MutilityWNA10> MutilityWA10){
                Decision10=3;
            }
            else{
                Decision10=4;
            }
            
            Decision10=(Decision10==obDecision10);
            
            //Likelihood contributions of the measurement systems
            if (Mfraclabor10==1 && Cchildcare10==0){
                
                loglikskills_rr10=MEASloglikSkills2010(SKMeasure1, M_A1_S2010,M_VS1_S2010,SKMeasure2, M_A2_S2010,M_VS2_S2010,SKMeasure3, M_A3_S2010,M_VS3_S2010,SKMeasure4, M_A4_S2010,M_VS4_S2010,SKMeasure5, M_A5_S2010,M_VS5_S2010,SKMeasure6, M_A6_S2010,M_VS6_S2010,SKMeasure7, M_A7_S2010,M_VS7_S2010,SKMeasure8, M_A8_S2010,M_VS8_S2010,SKMeasure9, M_A9_S2010,M_VS9_S2010,SKMeasure10, M_A10_S2010,M_VS10_S2010,SKMeasure11, M_A11_S2010,M_VS11_S2010,
                                                       CSkills_factor10WNA);
                
                
                
                
                //1. Skills
                //loglikskills_rr10=F_loglikelihood_generic(Ctestsfactor2ss_10,
                //CSkills_factor10WNA,MEASSkills);
                
                //Storing the particle weight and the particle itself
                
                what1_rr_vector[rr]=loglikskills_rr10;
                S1Vector[rr]=CSkills_factor10WNA;
                
                
                //2. Effort
                loglikeffm_rr10=MEASloglikEFFORT2010(MEFF1_10, M_A1_EFF10, M_V1_EFF10,
                                                     MEFF2_10, M_A2_EFF10, M_V2_EFF10,
                                                     MEFF3_10, M_A3_EFF10, M_V3_EFF10,
                                                     MEFF4_10, M_A4_EFF10, M_V4_EFF10,
                                                     MEFF5_10, M_A5_EFF10, M_V5_EFF10,
                                                     MEFF6_10, M_A6_EFF10, M_V6_EFF10,Meffort_factor10W);
                
                loglikefff_rr10=0;
                
                //3. Investment
                loglikeinv_rr10=MEASloginv10(MINV1_10,M_V1_INV10,M_A1_INV10,
                                             MINV2_10,M_V2_INV10,M_A2_INV10,
                                             MINV3_10,M_V3_INV10,M_A3_INV10,
                                             MINV4_10,M_V4_INV10,M_A4_INV10,
                                             MINV5_10,M_V5_INV10,M_A5_INV10,
                                             MINV6_10,M_V6_INV10,M_A6_INV10,
                                             MINV7_10,M_V7_INV10,M_A7_INV10,
                                             MINV8_10,M_V8_INV10,M_A8_INV10,
                                             CInv_factor10WNA);
                
                //Storing the final Decisions
                Feffort10Final=Feffort10;
                Meffort10Final=Meffort_factor10W;
                Investment10Final=CInv_factor10WNA;
                
                
                
            }
            
            if (Mfraclabor10==1 && Cchildcare10==1){
                
                //1. Skills
                loglikskills_rr10=MEASloglikSkills2010(SKMeasure1, M_A1_S2010,M_VS1_S2010,SKMeasure2, M_A2_S2010,M_VS2_S2010,SKMeasure3, M_A3_S2010,M_VS3_S2010,SKMeasure4, M_A4_S2010,M_VS4_S2010,SKMeasure5, M_A5_S2010,M_VS5_S2010,SKMeasure6, M_A6_S2010,M_VS6_S2010,SKMeasure7, M_A7_S2010,M_VS7_S2010,SKMeasure8, M_A8_S2010,M_VS8_S2010,SKMeasure9, M_A9_S2010,M_VS9_S2010,SKMeasure10, M_A10_S2010,M_VS10_S2010,SKMeasure11, M_A11_S2010,M_VS11_S2010,CSkills_factor10WA);
                
                //loglikskills_rr10=F_loglikelihood_generic(Ctestsfactor2ss_10,
                //CSkills_factor10WA,MEASSkills);
                
                //Storing the particle weight and the particle itself
                
                what1_rr_vector[rr]=loglikskills_rr10;
                S1Vector[rr]=CSkills_factor10WA;
                
                
                
                //2. Effort
                loglikeffm_rr10=MEASloglikEFFORT2010(MEFF1_10, M_A1_EFF10, M_V1_EFF10,
                                                     MEFF2_10, M_A2_EFF10, M_V2_EFF10,
                                                     MEFF3_10, M_A3_EFF10, M_V3_EFF10,
                                                     MEFF4_10, M_A4_EFF10, M_V4_EFF10,
                                                     MEFF5_10, M_A5_EFF10, M_V5_EFF10,
                                                     MEFF6_10, M_A6_EFF10, M_V6_EFF10,Meffort_factor10W);
                //loglikeffm_rr10=F_loglikelihood_generic(Meffort10,Meffort_factor10W,MEASMeffort);
                
                loglikefff_rr10=0;
                
                //3. Investment
                
                loglikeinv_rr10=MEASloginv10(MINV1_10,M_V1_INV10,M_A1_INV10,
                                             MINV2_10,M_V2_INV10,M_A2_INV10,
                                             MINV3_10,M_V3_INV10,M_A3_INV10,
                                             MINV4_10,M_V4_INV10,M_A4_INV10,
                                             MINV5_10,M_V5_INV10,M_A5_INV10,
                                             MINV6_10,M_V6_INV10,M_A6_INV10,
                                             MINV7_10,M_V7_INV10,M_A7_INV10,
                                             MINV8_10,M_V8_INV10,M_A8_INV10,
                                             CInv_factor10WA);
                
                
                //Storing the final Decisions
                Feffort10Final=Feffort10;
                Meffort10Final=Meffort_factor10W;
                Investment10Final=CInv_factor10WA;
                
                
            }
            if (Mfraclabor10==0 && Cchildcare10==0){
                
                //1. Skills
                loglikskills_rr10=MEASloglikSkills2010(SKMeasure1, M_A1_S2010,M_VS1_S2010,SKMeasure2, M_A2_S2010,M_VS2_S2010,SKMeasure3, M_A3_S2010,M_VS3_S2010,SKMeasure4, M_A4_S2010,M_VS4_S2010,SKMeasure5, M_A5_S2010,M_VS5_S2010,SKMeasure6, M_A6_S2010,M_VS6_S2010,SKMeasure7, M_A7_S2010,M_VS7_S2010,SKMeasure8, M_A8_S2010,M_VS8_S2010,SKMeasure9, M_A9_S2010,M_VS9_S2010,SKMeasure10, M_A10_S2010,M_VS10_S2010,SKMeasure11, M_A11_S2010,M_VS11_S2010,CSkills_factor10NWNA);
                
                //loglikskills_rr=F_loglikelihood_generic(Ctestsfactor2ss_10,
                //CSkills_factor10NWNA,MEASSkills);
                //Storing the particle weight and the particle itself
                
                what1_rr_vector[rr]=loglikskills_rr10;
                S1Vector[rr]=CSkills_factor10NWNA;
                
                
                
                
                
                //Storing the particle weight
                
                
                
                //2. Effort
                loglikeffm_rr10=MEASloglikEFFORT2010(MEFF1_10, M_A1_EFF10, M_V1_EFF10,
                                                     MEFF2_10, M_A2_EFF10, M_V2_EFF10,
                                                     MEFF3_10, M_A3_EFF10, M_V3_EFF10,
                                                     MEFF4_10, M_A4_EFF10, M_V4_EFF10,
                                                     MEFF5_10, M_A5_EFF10, M_V5_EFF10,
                                                     MEFF6_10, M_A6_EFF10, M_V6_EFF10,Meffort_factor10NW);
                loglikeffm_rr=
                //F_loglikelihood_generic(Meffort10,Meffort_factor10NW,MEASMeffort);
                
                
                loglikefff_rr=0;
                
                //3. Investment
                loglikeinv_rr10=MEASloginv10(MINV1_10,M_V1_INV10,M_A1_INV10,
                                             MINV2_10,M_V2_INV10,M_A2_INV10,
                                             MINV3_10,M_V3_INV10,M_A3_INV10,
                                             MINV4_10,M_V4_INV10,M_A4_INV10,
                                             MINV5_10,M_V5_INV10,M_A5_INV10,
                                             MINV6_10,M_V6_INV10,M_A6_INV10,
                                             MINV7_10,M_V7_INV10,M_A7_INV10,
                                             MINV8_10,M_V8_INV10,M_A8_INV10,
                                             CInv_factor10NWNA);
                //Storing the final Decisions
                Feffort10Final=Feffort10;
                Meffort10Final=Meffort_factor10NW;
                Investment10Final=CInv_factor10NWNA;
                
                
            }
            if (Mfraclabor10==0 && Cchildcare10==1){
                
                //1. Skills
                
                loglikskills_rr10=MEASloglikSkills2010(SKMeasure1, M_A1_S2010,M_VS1_S2010,SKMeasure2, M_A2_S2010,M_VS2_S2010,SKMeasure3, M_A3_S2010,M_VS3_S2010,SKMeasure4, M_A4_S2010,M_VS4_S2010,SKMeasure5, M_A5_S2010,M_VS5_S2010,SKMeasure6, M_A6_S2010,M_VS6_S2010,SKMeasure7, M_A7_S2010,M_VS7_S2010,SKMeasure8, M_A8_S2010,M_VS8_S2010,SKMeasure9, M_A9_S2010,M_VS9_S2010,SKMeasure10, M_A10_S2010,M_VS10_S2010,SKMeasure11, M_A11_S2010,M_VS11_S2010,CSkills_factor10NWA);
                
                
                //loglikskills_rr=F_loglikelihood_generic(Ctestsfactor2ss_10,
                //CSkills_factor10NWA,MEASSkills);
                
                //Storing the particle weight and the particle itself
                
                what1_rr_vector[rr]=loglikskills_rr10;
                S1Vector[rr]=CSkills_factor10NWA;
                
                
                
                //Storing the particle weight
                
                
                
                //2. Effort
                loglikeffm_rr10=MEASloglikEFFORT2010(MEFF1_10, M_A1_EFF10, M_V1_EFF10,
                                                     MEFF2_10, M_A2_EFF10, M_V2_EFF10,
                                                     MEFF3_10, M_A3_EFF10, M_V3_EFF10,
                                                     MEFF4_10, M_A4_EFF10, M_V4_EFF10,
                                                     MEFF5_10, M_A5_EFF10, M_V5_EFF10,
                                                     MEFF6_10, M_A6_EFF10, M_V6_EFF10,Meffort_factor10NW);
                
                //loglikeffm_rr=F_loglikelihood_generic(Meffort10, Meffort_factor10NW,MEASMeffort);
                //cout << " here28 "<< endl;
                
                loglikefff_rr=0;
                
                //3. Investment
                loglikeinv_rr10=MEASloginv10(MINV1_10,M_V1_INV10,M_A1_INV10,
                                             MINV2_10,M_V2_INV10,M_A2_INV10,
                                             MINV3_10,M_V3_INV10,M_A3_INV10,
                                             MINV4_10,M_V4_INV10,M_A4_INV10,
                                             MINV5_10,M_V5_INV10,M_A5_INV10,
                                             MINV6_10,M_V6_INV10,M_A6_INV10,
                                             MINV7_10,M_V7_INV10,M_A7_INV10,
                                             MINV8_10,M_V8_INV10,M_A8_INV10,
                                             CInv_factor10NWA);
                //Storing the final Decisions
                Feffort10Final=Feffort10;
                Meffort10Final=Meffort_factor10NW;
                Investment10Final=CInv_factor10NWA;
                
            }
            
            
            
        }//End of single mothers
        //If both parents are present.
        
        
        if (Cliveswithfather10==1 & Cliveswithmother10==1){
            
            //======================
            //0. Mmu
            //======================
            //Getting the draw for the unobserved heterogeneity of mmy
            MMushock10=gen_normal_3(generator);
            MMushock10=MMushock10*stdmmu;
            //Getting the predicted level of bargaining power
            mupred10=F_mmu(llambda0,llambda1,llambda2,llambda3,llambda4,
                           llambda5,llambda6,llambda7,llambda8,
                           Fwageinlik10,Mwageinlik10,Fnly10,Mnly10,MMushock,mmuLB,mmuUB,
                           Fage10,Mage10,Fyrschool12,Myrschool12,
                           MRATIO,Unemployment,Wageratio,Distance);
            
            
            
            //========================================================
            //1. Effort levels
            //========================================================
            
            //Draw a shock for the effort levels
            Meffort_factor10shock=exp(gen_normal_3(generator)*stdeffortmom);
            Feffort_factor10shock=exp(gen_normal_3(generator)*stdeffortfat);
            //Getting the means
            //1.1 Effort levels
            //1.1.1 NWNW
            MeffortMean10NWNW_TOG=F_effort_m10
            (mupred10,ggammaf,aalpha2f10,aalpha2m10,
             aalpha4f10,aalpha4m10,0,0,pphi,ttheta0,
             ttheta2,0.92)*Meffort_factor10shock;
            MEFFORTNWNW10NOSHOCK=F_effort_m10
            (mupred10,ggammaf,aalpha2f10,aalpha2m10,
             aalpha4f10,aalpha4m10,0,0,pphi,ttheta0,
             ttheta2,0.92);
            
            
            FeffortMean10NWNW_TOG=F_effort_f10
            (mupred10,ggammaf,aalpha2f10,aalpha2m10,
             aalpha4f,aalpha4m,0,0,pphi,ttheta0,
             ttheta2,0.92)*Feffort_factor10shock;
            FEFFORTNWNW10NOSHOCK=F_effort_f10
            (mupred10,ggammaf,aalpha2f10,aalpha2m10,
             aalpha4f,aalpha4m,0,0,pphi,ttheta0,
             ttheta2,0.92);
            
            //1.1.2 Father works mother not
            MeffortMean10WNW_TOG=F_effort_m10
            (mupred10,ggammaf,aalpha2f,aalpha2m,
             aalpha4f,aalpha4m,1,0,pphi,ttheta0,
             ttheta2,0.92)*Meffort_factor10shock;
            MEFFORTWNW10NOSHOCK=F_effort_m10
            (mupred10,ggammaf,aalpha2f,aalpha2m,
             aalpha4f,aalpha4m,1,0,pphi,ttheta0,
             ttheta2,0.92);
            
            FeffortMean10WNW_TOG=F_effort_f10
            (mupred10,ggammaf,aalpha2f10,aalpha2m10,
             aalpha4f,aalpha4m,1,0,pphi,ttheta0,
             ttheta2,0.92)*Feffort_factor10shock;
            FEFFORTWNW10NOSHOCK=F_effort_f10
            (mupred10,ggammaf,aalpha2f10,aalpha2m10,
             aalpha4f,aalpha4m,1,0,pphi,ttheta0,
             ttheta2,0.92);
            
            //1.1.3 Mother works father doesn't
            MeffortMean10NWW_TOG=F_effort_m10
            (mupred10,ggammaf,aalpha2f10,aalpha2m10,
             aalpha4f10,aalpha4m10,0,1,pphi,ttheta0,
             ttheta2,0.92)*Meffort_factor10shock;
            MEFFORTNWW10NOSHOCK=F_effort_m10
            (mupred10,ggammaf,aalpha2f10,aalpha2m10,
             aalpha4f10,aalpha4m10,0,1,pphi,ttheta0,
             ttheta2,0.92);
            
            
            FeffortMean10NWW_TOG=F_effort_f10
            (mupred10,ggammaf,aalpha2f10,aalpha2m10,
             aalpha4f10,aalpha4m10,0,1,pphi,ttheta0,
             ttheta2,0.92)*Feffort_factor10shock;
            FEFFORTNWW10NOSHOCK=F_effort_f10
            (mupred10,ggammaf,aalpha2f10,aalpha2m10,
             aalpha4f10,aalpha4m10,0,1,pphi,ttheta0,
             ttheta2,0.92);
            //1.1.4 Both work
            MeffortMean10WW_TOG=F_effort_m10
            (mupred10,ggammaf,aalpha2f10,aalpha2m10,
             aalpha4f10,aalpha4m10,1,1,pphi,ttheta0,
             ttheta2,0.92)*Meffort_factor10shock;
            MEFFORTWW10NOSHOCK=F_effort_m10
            (mupred10,ggammaf,aalpha2f10,aalpha2m10,
             aalpha4f10,aalpha4m10,1,1,pphi,ttheta0,
             ttheta2,0.92);
            
            FeffortMean10WW_TOG=F_effort_f10
            (mupred10,ggammaf,aalpha2f10,aalpha2m10,
             aalpha4f10,aalpha4m10,1,1,pphi,ttheta0,
             ttheta2,0.92)*Feffort_factor10shock;
            FEFFORTWW10NOSHOCK=F_effort_f10
            (mupred10,ggammaf,aalpha2f10,aalpha2m10,
             aalpha4f10,aalpha4m10,1,1,pphi,ttheta0,
             ttheta2,0.92);
            //Need to fix the effort levels either if they are negative or
            //if they are too large
            //Effort can't be negative, we are drawing from a truncated normal distribution
            //Fix the effort levels
            //Mother fixing
            if (MeffortMean10NWNW_TOG<0){
                MeffortMean10NWNW_TOG=0;
            }
            if (MeffortMean10NWNW_TOG>=1.0e+6){
                MeffortMean10NWNW_TOG=1.0e+6;
            }
            
            
            if (MeffortMean10WNW_TOG<0){
                MeffortMean10WNW_TOG=0;
            }
            if (MeffortMean10WNW_TOG>=1.0e+6){
                MeffortMean10WNW_TOG=1.0e+6;
            }
            
            
            if (MeffortMean10NWW_TOG<0){
                MeffortMean10NWW_TOG=0;
            }
            if (MeffortMean10NWW_TOG>=1.0e+6){
                MeffortMean10NWW_TOG=1.0e+6;
            }
            
            
            if (MeffortMean10WW_TOG<0){
                MeffortMean10WW_TOG=0;
            }
            if (MeffortMean10WW_TOG>=1.0e+6){
                MeffortMean10WW_TOG=1.0e+6;
            }
            //Father fixing
            if (FeffortMean10NWNW_TOG<0){
                FeffortMean10NWNW_TOG=0;
            }
            if (FeffortMean10NWNW_TOG>=1.0e+6){
                FeffortMean10NWNW_TOG=1.0e+6;
            }
            
            
            if (FeffortMean10WNW_TOG<0){
                FeffortMean10WNW_TOG=0;
            }
            if (FeffortMean10WNW_TOG>=1.0e+6){
                FeffortMean10WNW_TOG=1.0e+6;
            }
            
            
            if (FeffortMean10NWW_TOG<0){
                FeffortMean10NWW_TOG=0;
            }
            if (FeffortMean10NWW_TOG>=1.0e+6){
                FeffortMean10NWW_TOG=1.0e+6;
            }
            
            
            if (FeffortMean10WW_TOG<0){
                FeffortMean10WW_TOG=0;
            }
            if (FeffortMean10WW_TOG>=1.0e+6){
                FeffortMean10WW_TOG=1.0e+6;
            }
            
            
            //========================================================
            //2. Investment levels
            //========================================================
            //Getting the draw for the investment shock
            CInv_factor10shock=gen_normal_3(generator);
            //========================================================
            //2.1 If mom doesn't work
            //========================================================
            //3. Now transforming it as a draw from the given distribution
            CInv_factor10shock=CInv_factor10shock*stdinvestment;
            //4. Now get the mean of the investment as a function of mmu
            //for the possible combinations
            //1. NWNWNA
            CInvestmentMean10NWNWNA_TOG=F_invcouple10
            (aalpha1f10,aalpha1m10,aalpha2f,aalpha2m,
             aalpha2f10,aalpha2m10,
             mupred10,0,0,Fwageinlik10,Mwageinlik10,
             Fnly10,Mnly10,ttheta0,
             ttheta1,priceINV)*exp(CInv_factor10shock);
            if (CInvestmentMean10NWNWNA_TOG<=0){
                CInvestmentMean10NWNWNA_TOG=1.0e-5;
            }
            CINVFACTOR10NWNWNA_noshock=F_invcouple10
            (aalpha1f10,aalpha1m10,aalpha2f,aalpha2m,
             aalpha2f10,aalpha2m10,
             mupred10,0,0,Fwageinlik10,Mwageinlik10,
             Fnly10,Mnly10,ttheta0,
             ttheta1,priceINV);
            //2. NWNWA
            CInvestmentMean10NWNWA_TOG=F_invcouple10
            (aalpha1f10,aalpha1m10,aalpha2f,aalpha2m,
             aalpha2f10,aalpha2m10,
             mupred10,0,0,Fwageinlik10,Mwageinlik10,
             Fnly10,Mnly10-pricechildcare,ttheta0,
             ttheta1,priceINV)*exp(CInv_factor10shock);
            if (CInvestmentMean10NWNWA_TOG<=0){
                CInvestmentMean10NWNWA_TOG=1.0e-5;
            }
            CINVFACTOR10NWNWA_noshock=F_invcouple10
            (aalpha1f10,aalpha1m10,aalpha2f,aalpha2m,
             aalpha2f10,aalpha2m10,
             mupred10,0,0,Fwageinlik10,Mwageinlik10,
             Fnly10,Mnly10-pricechildcare,ttheta0,
             ttheta1,priceINV);
            //3. WNWNA
            CInvestmentMean10WNWNA_TOG=F_invcouple10
            (aalpha1f10,aalpha1m10,aalpha2f,aalpha2m,
             aalpha2f10,aalpha2m10,
             mupred10,1,0,Fwageinlik10,Mwageinlik10,
             Fnly10,Mnly10,ttheta0,
             ttheta1,priceINV)*exp(CInv_factor10shock);
            
            CINVFACTOR10WNWNA_noshock=F_invcouple10
            (aalpha1f10,aalpha1m10,aalpha2f,aalpha2m,
             aalpha2f10,aalpha2m10,
             mupred10,1,0,Fwageinlik10,Mwageinlik10,
             Fnly10,Mnly10,ttheta0,
             ttheta1,priceINV);
            if (CInvestmentMean10WNWNA_TOG<=0){
                CInvestmentMean10WNWNA_TOG=1.0e-5;
            }
            
            //4. WNWA
            CInvestmentMean10WNWA_TOG=F_invcouple10
            (aalpha1f10,aalpha1m10,aalpha2f,aalpha2m,
             aalpha2f10,aalpha2m10,
             mupred10,1,0,Fwageinlik10,Mwageinlik10,
             Fnly10,Mnly10-pricechildcare,ttheta0,
             ttheta1,priceINV)*exp(CInv_factor10shock);
            CINVFACTOR10WNWA_noshock=F_invcouple10
            (aalpha1f10,aalpha1m10,aalpha2f,aalpha2m,
             aalpha2f10,aalpha2m10,
             mupred10,1,0,Fwageinlik10,Mwageinlik10,
             Fnly10,Mnly10-pricechildcare,ttheta0,
             ttheta1,priceINV);
            if (CInvestmentMean10WNWA_TOG<=0){
                CInvestmentMean10WNWA_TOG=1.0e-5;
            }
            
            //5. NWWNA
            CInvestmentMean10NWWNA_TOG=F_invcouple10
            (aalpha1f10,aalpha1m10,aalpha2f,aalpha2m,
             aalpha2f10,aalpha2m10,
             mupred10,0,1,Fwageinlik10,Mwageinlik10,
             Fnly10,Mnly10,ttheta0,
             ttheta1,priceINV)*exp(CInv_factor10shock);
            CINVFACTOR10NWWNA_noshock=F_invcouple10
            (aalpha1f10,aalpha1m10,aalpha2f,aalpha2m,
             aalpha2f10,aalpha2m10,
             mupred10,0,1,Fwageinlik10,Mwageinlik10,
             Fnly10,Mnly10,ttheta0,
             ttheta1,priceINV);
            if (CInvestmentMean10NWWNA_TOG<=0){
                CInvestmentMean10NWWNA_TOG=1.0e-5;
            }
            
            //6. NWWA
            CInvestmentMean10NWWA_TOG=F_invcouple10
            (aalpha1f10,aalpha1m10,aalpha2f,aalpha2m,
             aalpha2f10,aalpha2m10,
             mupred10,0,1,Fwageinlik10,Mwageinlik10,
             Fnly10,Mnly10-pricechildcare,ttheta0,
             ttheta1,priceINV)*exp(CInv_factor10shock);
            CINVFACTOR10NWWA_noshock=F_invcouple10
            (aalpha1f10,aalpha1m10,aalpha2f,aalpha2m,
             aalpha2f10,aalpha2m10,
             mupred10,0,1,Fwageinlik10,Mwageinlik10,
             Fnly10,Mnly10-pricechildcare,ttheta0,
             ttheta1,priceINV);
            if (CInvestmentMean10NWWA_TOG<=0){
                CInvestmentMean10NWWA_TOG=1.0e-5;
            }
            
            //7. WWNA
            CInvestmentMean10WWNA_TOG=F_invcouple10
            (aalpha1f10,aalpha1m10,aalpha2f,aalpha2m,
             aalpha2f10,aalpha2m10,
             mupred10,1,1,Fwageinlik10,Mwageinlik10,
             Fnly10,Mnly10,ttheta0,
             ttheta1,priceINV)*exp(CInv_factor10shock);
            CINVFACTOR10WWNA_noshock=F_invcouple10
            (aalpha1f10,aalpha1m10,aalpha2f,aalpha2m,
             aalpha2f10,aalpha2m10,
             mupred10,1,1,Fwageinlik10,Mwageinlik10,
             Fnly10,Mnly10,ttheta0,
             ttheta1,priceINV);
            if (CInvestmentMean10WWNA_TOG<=0){
                CInvestmentMean10WWNA_TOG=1.0e-5;
            }
            
            
            //8. WWA
            CInvestmentMean10WWA_TOG=F_invcouple10
            (aalpha1f10,aalpha1m10,aalpha2f,aalpha2m,
             aalpha2f10,aalpha2m10,
             mupred10,1,1,Fwageinlik10,Mwageinlik10,
             Fnly10,Mnly10-pricechildcare,ttheta0,
             ttheta1,priceINV)*exp(CInv_factor10shock);
            CINVFACTOR10WWA_noshock=F_invcouple10
            (aalpha1f10,aalpha1m10,aalpha2f,aalpha2m,
             aalpha2f10,aalpha2m10,
             mupred10,1,1,Fwageinlik10,Mwageinlik10,
             Fnly10,Mnly10-pricechildcare,ttheta0,
             ttheta1,priceINV);
            if (CInvestmentMean10WWA_TOG<=0){
                CInvestmentMean10WWA_TOG=1.0e-5;
            }
            //==========================================================
            //3. Skills levels
            //==========================================================
            //Shock of skills
            CSkills_factor10shock=gen_normal_3(generator);
            if (CSkills_factor10shock*stdskills<500){
                CSkills_factor10shock=exp(CSkills_factor10shock*stdskills);
            }
            if (CSkills_factor10shock*stdskills>=500){
                CSkills_factor10shock=1000;
            }
            
            //Now getting the different levels of skills for all possi-
            //bilities
            //3.1 NWNWNA
            CSkills_factor10NWNWNA_TOG=F_predskills
            (ddelta0,ddelta1,ddelta2,ddelta3,ddelta4,Cedad_meses10,
             ttheta0,ttheta1,ttheta2,pphi,
             ggammaf,ggammam,FEFFORTNWNW10NOSHOCK,
             MEFFORTNWNW10NOSHOCK,
             CINVFACTOR10NWNWNA_noshock,
             S0Vector[rr],0,PGVector[rr],Hmemberstotal10)*CSkills_factor10shock;
            
            
            
            //3.1 NWNWA
            CSkills_factor10NWNWA_TOG=F_predskills
            (ddelta0,ddelta1,ddelta2,ddelta3,ddelta4,Cedad_meses10,
             ttheta0,ttheta1,ttheta2,pphi,
             ggammaf,ggammam,FEFFORTNWNW10NOSHOCK,
             MEFFORTNWNW10NOSHOCK,
             CINVFACTOR10NWNWA_noshock,
             S0Vector[rr],1,PGVector[rr],Hmemberstotal10)*CSkills_factor10shock;
            
            
            //3.2 WNWNA
            CSkills_factor10WNWNA_TOG=F_predskills
            (ddelta0,ddelta1,ddelta2,ddelta3,ddelta4,Cedad_meses10,
             ttheta0,ttheta1,ttheta2,pphi,
             ggammaf,ggammam,FEFFORTWNW10NOSHOCK,
             MEFFORTWNW10NOSHOCK,
             CINVFACTOR10WNWNA_noshock,
             S0Vector[rr],0,PGVector[rr],Hmemberstotal10)*CSkills_factor10shock;
            
            //3.2 WNWA
            CSkills_factor10WNWA_TOG=F_predskills
            (ddelta0,ddelta1,ddelta2,ddelta3,ddelta4,Cedad_meses10,
             ttheta0,ttheta1,ttheta2,pphi,
             ggammaf,ggammam,FEFFORTWNW10NOSHOCK,
             MEFFORTWNW10NOSHOCK,
             CINVFACTOR10WNWA_noshock,
             S0Vector[rr],1,PGVector[rr],Hmemberstotal10)*CSkills_factor10shock;
            
            
            //3.3 NWWNA
            CSkills_factor10NWWNA_TOG=F_predskills
            (ddelta0,ddelta1,ddelta2,ddelta3,ddelta4,Cedad_meses10,
             ttheta0,ttheta1,ttheta2,pphi,
             ggammaf,ggammam,FEFFORTNWW10NOSHOCK,
             MEFFORTNWW10NOSHOCK,
             CINVFACTOR10NWWNA_noshock,
             S0Vector[rr],0,PGVector[rr],Hmemberstotal10)*CSkills_factor10shock;
            
            //3.3 NWWA
            CSkills_factor10NWWA_TOG=F_predskills
            (ddelta0,ddelta1,ddelta2,ddelta3,ddelta4,Cedad_meses10,
             ttheta0,ttheta1,ttheta2,pphi,
             ggammaf,ggammam,FEFFORTNWW10NOSHOCK,
             MEFFORTNWW10NOSHOCK,
             CINVFACTOR10NWWA_noshock,
             S0Vector[rr],1,PGVector[rr],Hmemberstotal10)*CSkills_factor10shock;
            
            
            //3.4 WWNA
            CSkills_factor10WWNA_TOG=F_predskills
            (ddelta0,ddelta1,ddelta2,ddelta3,ddelta4,Cedad_meses10,
             ttheta0,ttheta1,ttheta2,pphi,
             ggammaf,ggammam,FEFFORTWW10NOSHOCK,
             MEFFORTWW10NOSHOCK,
             CINVFACTOR10WWNA_noshock,
             S0Vector[rr],0,PGVector[rr],Hmemberstotal10)*CSkills_factor10shock;
            
            //3.4 WWA
            
            CSkills_factor10WWA_TOG=F_predskills
            (ddelta0,ddelta1,ddelta2,ddelta3,ddelta4,Cedad_meses10,
             ttheta0,ttheta1,ttheta2,pphi,
             ggammaf,ggammam,FEFFORTWW10NOSHOCK,
             MEFFORTWW10NOSHOCK,
             CINVFACTOR10WWA_noshock,
             S0Vector[rr],1,PGVector[rr],Hmemberstotal10)*CSkills_factor10shock;
            
            //===========================================================
            //===========================================================
            //4 Preference shocks
            
            //1. If mom doesn't work
            MprefshockNWNA10=gen_normal_3(generator);
            MprefshockNWNA10=MprefshockNWNA10*MshockNWNA;
            
            MprefshockNWA10=gen_normal_3(generator);
            MprefshockNWA10=MprefshockNWA10*MshockNWA;
            
            
            //2. If mom works
            MprefshockWNA10=gen_normal_3(generator);
            MprefshockWNA10=MprefshockWNA10*MshockWNA;
            
            MprefshockWA10=gen_normal_3(generator);
            MprefshockWA10=MprefshockWA10*MshockWA;
            
            //3. If FATHER doesn't work
            FprefshockNWNA10=gen_normal_3(generator);
            FprefshockNWNA10=FprefshockWNA10*FshockNWNA;
            
            FprefshockWNA10=gen_normal_3(generator);
            FprefshockWNA10=FprefshockWNA10*FshockWNA;
            
            //4. If FATHER works
            FprefshockWNA10=gen_normal_3(generator);
            FprefshockWNA10=FprefshockWNA10*FshockWNA;
            
            FprefshockWA10=gen_normal_3(generator);
            FprefshockWA10=FprefshockWA10*FshockWA;
            
            
            //===========================================================
            //===========================================================
            //5 Consumption
            //===========================================================
            //1. Mother
            //===========================================================
            
            //1. NWNW
            MconsumptionNWNWNA_TOG10=M_consumption_TOG10
            (aalpha1m10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
             CInvestmentMean10NWNWNA_TOG,ttheta0,ttheta1,mupred10,priceINV);
            
            MconsumptionNWNWA_TOG10=M_consumption_TOG10
            (aalpha1m10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
             CInvestmentMean10NWNWA_TOG,ttheta0,ttheta1,mupred10,priceINV);
            
            //2. WNW
            MconsumptionWNWNA_TOG10=M_consumption_TOG10
            (aalpha1m10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
             CInvestmentMean10WNWNA_TOG,ttheta0,ttheta1,mupred10,priceINV);
            
            MconsumptionWNWA_TOG10=M_consumption_TOG10
            (aalpha1m10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
             CInvestmentMean10WNWA_TOG,ttheta0,ttheta1,mupred10,priceINV);
            
            //3. NWW
            MconsumptionNWWNA_TOG10=M_consumption_TOG10
            (aalpha1m10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
             CInvestmentMean10NWWNA_TOG,ttheta0,ttheta1,mupred10,priceINV);
            
            MconsumptionNWWA_TOG10=M_consumption_TOG10
            (aalpha1m10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
             CInvestmentMean10NWWA_TOG,ttheta0,ttheta1,mupred10,priceINV);
            //4. WW
            MconsumptionWWNA_TOG10=M_consumption_TOG10
            (aalpha1m10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
             CInvestmentMean10WWNA_TOG,ttheta0,ttheta1,mupred10,priceINV);
            
            MconsumptionWWA_TOG10=M_consumption_TOG10
            (aalpha1m10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
             CInvestmentMean10WWA_TOG,ttheta0,ttheta1,mupred10,priceINV);
            
            //===========================================================
            //1. Father
            //===========================================================
            
            //1. NWNW
            FconsumptionNWNWNA_TOG10=F_consumption_TOG10
            (aalpha1f10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
             CInvestmentMean10NWNWNA_TOG,ttheta0,ttheta1,mupred10,priceINV);
            
            FconsumptionNWNWA_TOG10=F_consumption_TOG10
            (aalpha1f10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
             CInvestmentMean10NWNWA_TOG,ttheta0,ttheta1,mupred10,priceINV);
            //2. WNW
            FconsumptionWNWNA_TOG10=F_consumption_TOG10
            (aalpha1f10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
             CInvestmentMean10WNWNA_TOG,ttheta0,ttheta1,mupred10,priceINV);
            
            FconsumptionWNWA_TOG10=F_consumption_TOG10
            (aalpha1f10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
             CInvestmentMean10WNWA_TOG,ttheta0,ttheta1,mupred10,priceINV);
            //3. NWW
            FconsumptionNWWNA_TOG10=F_consumption_TOG10
            (aalpha1f10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
             CInvestmentMean10NWWNA_TOG,ttheta0,ttheta1,mupred10,priceINV);
            
            FconsumptionNWWA_TOG10=F_consumption_TOG10
            (aalpha1f10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
             CInvestmentMean10NWWA_TOG,ttheta0,ttheta1,mupred10,priceINV);
            
            //4 WW
            FconsumptionWWNA_TOG10=F_consumption_TOG10
            (aalpha1f10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
             CInvestmentMean10WWNA_TOG,ttheta0,ttheta1,mupred10,priceINV);
            
            FconsumptionWWA_TOG10=F_consumption_TOG10
            (aalpha1f10,aalpha2f,aalpha2m,aalpha2f10,aalpha2m10,
             CInvestmentMean10WWA_TOG,ttheta0,ttheta1,mupred10,priceINV);
            
            //===========================================================
            //===========================================================
            //6. Behavioral model
            //==========================================================
            //1. Mother
            //==========================================================
            //1. NWNW
            
            MutilityNWNWNA_TOG10=F_utility10
            (aalpha1m10,aalpha2m10,aalpha3m10,
             aalpha4m10,aalpha5m10,MconsumptionNWNWNA_TOG10,MeffortMean10NWNW_TOG,0,CSkills_factor10NWNWNA_TOG,0)
            +MprefshockNWNA10;
            
            MutilityNWNWA_TOG10=F_utility10
            (aalpha1m10,aalpha2m10,aalpha3m10,
             aalpha4m10,aalpha5m10,MconsumptionNWNWA_TOG10,MeffortMean10NWNW_TOG,0,CSkills_factor10NWNWA_TOG,1)
            +MprefshockNWA10;
            
            //2. WNW
            
            MutilityWNWNA_TOG10=F_utility10
            (aalpha1m10,aalpha2m10,aalpha3m10,aalpha4m10,aalpha5m10,
             MconsumptionWNWNA_TOG10,MeffortMean10WNW_TOG,0,
             CSkills_factor10WNWNA_TOG,0)+MprefshockNWNA10;
            
            MutilityWNWA_TOG10=F_utility10
            (aalpha1m10,aalpha2m10,aalpha3m10,aalpha4m10,aalpha5m10,
             MconsumptionWNWA_TOG10,MeffortMean10WNW_TOG,0,
             CSkills_factor10WNWA_TOG,1)+MprefshockNWA10;
            
            
            //3. NWW
            
            MutilityNWWNA_TOG10=F_utility10
            (aalpha1m10,aalpha2m10,aalpha3m10,aalpha4m10,aalpha5m10,
             MconsumptionNWWNA_TOG10,MeffortMean10NWW_TOG,1,
             CSkills_factor10NWWNA_TOG,0)+MprefshockWNA10;
            
            MutilityNWWA_TOG10=F_utility10
            (aalpha1m10,aalpha2m10,aalpha3m10,aalpha4m10,aalpha5m10,
             MconsumptionNWWA_TOG10,MeffortMean10NWW_TOG,1,
             CSkills_factor10NWWA_TOG,1)+MprefshockWA10;
            
            //4. WW
            MutilityWWNA_TOG10=F_utility10
            (aalpha1m10,aalpha2m10,aalpha3m10,aalpha4m10,aalpha5m10,
             MconsumptionWWNA_TOG10, MeffortMean10WW_TOG,1,
             CSkills_factor10WWNA_TOG,0)+MprefshockWNA10;
            
            MutilityWWA_TOG10=F_utility10
            (aalpha1m10,aalpha2m10,aalpha3m10,aalpha4m10,aalpha5m10,
             MconsumptionWWA_TOG10, MeffortMean10WW_TOG,1,
             CSkills_factor10WWA_TOG,1)+MprefshockWA10;
            //==========================================================
            //2. Father
            //==========================================================
            
            //1. NWNW
            FutilityNWNWNA_TOG10=F_utility10
            (aalpha1f10,aalpha2f10,aalpha3f10,aalpha4f10,aalpha5f10,
             FconsumptionNWNWNA_TOG10, FeffortMean10NWNW_TOG,0,
             CSkills_factor10NWNWNA_TOG,0)+ FprefshockNWNA10;
            
            FutilityNWNWA_TOG10=F_utility10
            (aalpha1f10,aalpha2f10,aalpha3f10,aalpha4f10,aalpha5f10,
             FconsumptionNWNWA_TOG10, FeffortMean10NWNW_TOG,0,
             CSkills_factor10NWNWA_TOG,1)+ FprefshockNWA10;
            
            //2. WNW
            
            FutilityWNWNA_TOG10=F_utility10
            (aalpha1f10,aalpha2f10,aalpha3f10,aalpha4f10,aalpha5f10,
             FconsumptionWNWNA_TOG10,FeffortMean10WNW_TOG,1,
             CSkills_factor10WNWNA_TOG,0)+FprefshockWNA10;
            
            FutilityWNWA_TOG10=F_utility10
            (aalpha1f10,aalpha2f10,aalpha3f10,aalpha4f10,aalpha5f10,
             FconsumptionWNWA_TOG10,FeffortMean10WNW_TOG,1,
             CSkills_factor10WNWA_TOG,1)+FprefshockWA10;
            //3.NWW
            FutilityNWWNA_TOG10=F_utility10
            (aalpha1f10,aalpha2f10,aalpha3f10,aalpha4f10,aalpha5f10,
             FconsumptionNWWNA_TOG,FeffortMean10NWW_TOG,0,
             CSkills_factor10NWWNA_TOG,0)+FprefshockNWNA10;
            
            FutilityNWWA_TOG10=F_utility10
            (aalpha1f10,aalpha2f10,aalpha3f10,aalpha4f10,aalpha5f10,
             FconsumptionNWWA_TOG10,FeffortMean10NWW_TOG,0,
             CSkills_factor10NWWA_TOG,1)+FprefshockNWA10;
            
            //4. WW
            FutilityWWNA_TOG10=F_utility10
            (aalpha1f10,aalpha2f10,aalpha3f10,aalpha4f10,aalpha5f10,
             FconsumptionWWNA_TOG10,FeffortMean12WW_TOG,1,
             CSkills_factor10WWNA_TOG,0)+FprefshockWNA10;
            
            FutilityWWA_TOG10=F_utility10
            (aalpha1f10,aalpha2f10,aalpha3f10,aalpha4f10,aalpha5f10,
             FconsumptionWWA_TOG10,FeffortMean12WW_TOG,1,
             CSkills_factor10WWA_TOG,1)+FprefshockWA10;
            
            //Seeing the utility levels for each case:
            
            
            //==========================================================
            // Defining the utilities for the possible combinations
            //==========================================================
            
            WelfareNWNWNA10=mupred10*FutilityNWNWNA_TOG10+
            (1-mupred10)*MutilityNWNWNA_TOG10;
            
            WelfareNWNWA10=mupred10*FutilityNWNWA_TOG10+
            (1-mupred10)*MutilityNWNWA_TOG10;
            
            
            WelfareWNWNA10=mupred10*FutilityWNWNA_TOG10+
            (1-mupred10)*MutilityWNWNA_TOG10;
            
            WelfareWNWA10=mupred10*FutilityWNWA_TOG10+
            (1-mupred10)*MutilityWNWA_TOG10;
            
            WelfareNWWNA10=mupred10*FutilityNWWNA_TOG10+
            (1-mupred10)*MutilityNWWNA_TOG10;
            
            WelfareNWWA10=mupred10*FutilityNWWA_TOG10+
            (1-mupred10)*MutilityNWWA_TOG10;
            
            WelfareWWNA10=mupred10*FutilityWWNA_TOG10+
            (1-mupred10)*MutilityWWNA_TOG10;
            
            WelfareWWA10=mupred10*FutilityWWA_TOG10+
            (1-mupred10)*MutilityWWA_TOG10;
            
            
            
            //==========================================================
            //Identifying the correct decision
            //==========================================================
            if ( (WelfareNWNWNA10>WelfareNWNWA10) &&
                ( WelfareNWNWNA10>WelfareNWWNA10)  &&
                ( WelfareNWNWNA10>WelfareNWWA10)  &&
                ( WelfareNWNWNA10>WelfareWNWNA10)  &&
                ( WelfareNWNWNA10>WelfareWNWA10)  &&
                ( WelfareNWNWNA10>WelfareWWNA10)  &&
                ( WelfareNWNWNA10>WelfareWWA10)) {
                Decision10=8;
            }
            else if ((WelfareNWNWA10>WelfareNWWNA10) &&
                     (WelfareNWNWA10>WelfareNWWA10) &&
                     (WelfareNWNWA10>WelfareWNWNA10) &&
                     (WelfareNWNWA10>WelfareWNWA10) &&
                     (WelfareNWNWA10>WelfareWWNA10) &&
                     (WelfareNWNWA10>WelfareWWA10)) {
                Decision10=7;
            }
            
            else if ((WelfareNWWNA10>WelfareNWWA10) &&
                     (WelfareNWWNA10>WelfareWNWNA10) &&
                     (WelfareNWWNA10>WelfareWNWA10) &&
                     (WelfareNWWNA10>WelfareWWNA10) &&
                     (WelfareNWWNA10>WelfareWWA10)) {
                Decision10=6;
            }
            
            else if ((WelfareNWWA10>WelfareWNWNA10) &&
                     (WelfareNWWA10>WelfareWNWA10) &&
                     (WelfareNWWA10>WelfareWWNA10) &&
                     (WelfareNWWA10>WelfareWWA10)) {
                Decision10=5;
            }
            
            else if ((WelfareWNWNA10>WelfareWNWA10) &&
                     (WelfareWNWNA10>WelfareWWNA10) &&
                     (WelfareWNWNA10>WelfareWWA10)) {
                Decision10=4;
            }
            
            else if ((WelfareWNWA10>WelfareWWNA10) &&
                     (WelfareWNWA10>WelfareWWA10)) {
                Decision10=3;
            }
            
            else if ((WelfareWWNA10>WelfareWWA10)) {
                Decision10=2;
            }
            
            else{
                Decision10=1;
            }
            
            
            
            //Likelihood of observed decisions
            
            Decision10STORED=Decision10;
            
            //Storing the value of the maximum decision in 2010
            VMAX10=(Decision10==1)*WelfareWWA10+
            (Decision10==2)*WelfareWWNA10+
            (Decision10==3)*WelfareWNWA10+
            (Decision10==4)*WelfareWNWNA10+
            (Decision10==5)*WelfareNWWA10+
            (Decision10==6)*WelfareNWWNA10+
            (Decision10==7)*WelfareNWNWA10+
            (Decision10==8)*WelfareNWNWNA10;
            
            
            //Storing the value of the obseved decision in 2010
            VOBS10=
            (obDecision10==1)*WelfareWWA10+
            (obDecision10==2)*WelfareWWNA10+
            (obDecision10==3)*WelfareWNWA10+
            (obDecision10==4)*WelfareWNWNA10+
            (obDecision10==5)*WelfareNWWA10+
            (obDecision10==6)*WelfareNWWNA10+
            (obDecision10==7)*WelfareNWNWA10+
            (obDecision10==8)*WelfareNWNWNA10;
            
            VSUM10=(exp(WelfareWWA10-VMAX10)+
                    exp(WelfareWWNA10-VMAX10)+
                    exp(WelfareWNWA10-VMAX10)+
                    exp(WelfareWNWNA10-VMAX10)+
                    exp(WelfareNWWA10-VMAX10)+
                    exp(WelfareNWWNA10-VMAX10)+
                    exp(WelfareNWNWA10-VMAX10)+
                    exp(WelfareNWNWNA10-VMAX10))/TAUSMOOTHING;
            
            VDIF10=exp(VOBS10-VMAX10)/TAUSMOOTHING;
            
            VCONTRIB10=VDIF10/VSUM10;
            
            Decision10=(Decision10==obDecision10);
            //7. Likelihood of measurement system
            //===============================
            
            //7.1 NWNWNA
            if (Mfraclabor10==0 && Ffraclabor10==0 && Cchildcare10==0){
                //1. Skills
                
                loglikskills_rr10=MEASloglikSkills2010(SKMeasure1, M_A1_S2010,M_VS1_S2010,SKMeasure2, M_A2_S2010,M_VS2_S2010,SKMeasure3, M_A3_S2010,M_VS3_S2010,SKMeasure4, M_A4_S2010,M_VS4_S2010,SKMeasure5, M_A5_S2010,M_VS5_S2010,SKMeasure6, M_A6_S2010,M_VS6_S2010,SKMeasure7, M_A7_S2010,M_VS7_S2010,SKMeasure8, M_A8_S2010,M_VS8_S2010,SKMeasure9, M_A9_S2010,M_VS9_S2010,SKMeasure10, M_A10_S2010,M_VS10_S2010,SKMeasure11, M_A11_S2010,M_VS11_S2010,CSkills_factor10NWNWNA_TOG);
                
                //loglikskills_rr10=F_loglikelihood_generic(Ctestsfactor2ss_10,
                //CSkills_factor10NWNWNA_TOG,MEASSkills);
                
                //Storing the particle weight and the particle itself
                
                what1_rr_vector[rr]=loglikskills_rr10;
                S1Vector[rr]=CSkills_factor10NWNWNA_TOG;
                
                
                //2. Mother's effort
                loglikeffm_rr10=MEASloglikEFFORT2010(MEFF1_10, M_A1_EFF10, M_V1_EFF10,
                                                     MEFF2_10, M_A2_EFF10, M_V2_EFF10,
                                                     MEFF3_10, M_A3_EFF10, M_V3_EFF10,
                                                     MEFF4_10, M_A4_EFF10, M_V4_EFF10,
                                                     MEFF5_10, M_A5_EFF10, M_V5_EFF10,
                                                     MEFF6_10, M_A6_EFF10, M_V6_EFF10,MeffortMean10NWNW_TOG);
                
                //loglikeffm_rr10=F_loglikelihood_generic(Meffort10, MeffortMean10NWNW_TOG,MEASMeffort);
                
                //3. Father's effort
                loglikefff_rr10=MEASloglikEFFORT2010(FEFF1_10, M_A1_EFF10, M_V1_EFF10,
                                                     FEFF2_10, M_A2_EFF10, M_V2_EFF10,
                                                     FEFF3_10, M_A3_EFF10, M_V3_EFF10,
                                                     FEFF4_10, M_A4_EFF10, M_V4_EFF10,
                                                     FEFF5_10, M_A5_EFF10, M_V5_EFF10,
                                                     FEFF6_10, M_A6_EFF10, M_V6_EFF10,FeffortMean10NWNW_TOG);
                
                //loglikefff_rr10=F_loglikelihood_generic(Feffort10,FeffortMean10NWNW_TOG,MEASFeffort);
                
                //4. Investment level
                loglikeinv_rr10=MEASloginv10(MINV1_10,M_V1_INV10,M_A1_INV10,
                                             MINV2_10,M_V2_INV10,M_A2_INV10,
                                             MINV3_10,M_V3_INV10,M_A3_INV10,
                                             MINV4_10,M_V4_INV10,M_A4_INV10,
                                             MINV5_10,M_V5_INV10,M_A5_INV10,
                                             MINV6_10,M_V6_INV10,M_A6_INV10,
                                             MINV7_10,M_V7_INV10,M_A7_INV10,
                                             MINV8_10,M_V8_INV10,M_A8_INV10,
                                             CInvestmentMean10NWNWNA_TOG);
                //Storing the final Decisions
                Feffort10Final=FeffortMean10NWNW_TOG;
                Meffort10Final=MeffortMean10NWNW_TOG;
                Investment10Final=CInvestmentMean10NWNWNA_TOG;
            }
            
            //7.2 NWNWA
            if (Mfraclabor10==0 && Ffraclabor10==0 && Cchildcare10==1){
                loglikskills_rr10=MEASloglikSkills2010(SKMeasure1, M_A1_S2010,M_VS1_S2010,SKMeasure2, M_A2_S2010,M_VS2_S2010,SKMeasure3, M_A3_S2010,M_VS3_S2010,SKMeasure4, M_A4_S2010,M_VS4_S2010,SKMeasure5, M_A5_S2010,M_VS5_S2010,SKMeasure6, M_A6_S2010,M_VS6_S2010,SKMeasure7, M_A7_S2010,M_VS7_S2010,SKMeasure8, M_A8_S2010,M_VS8_S2010,SKMeasure9, M_A9_S2010,M_VS9_S2010,SKMeasure10, M_A10_S2010,M_VS10_S2010,SKMeasure11, M_A11_S2010,M_VS11_S2010,CSkills_factor10NWNWA_TOG);
                //1. Skills
                //loglikskills_rr10=F_loglikelihood_generic(Ctestsfactor2ss_10,
                //CSkills_factor10NWNWA_TOG,MEASSkills);
                
                //Storing the particle weight and the particle itself
                
                what1_rr_vector[rr]=loglikskills_rr10;
                S1Vector[rr]=CSkills_factor10NWNWA_TOG;
                
                
                //2. Mother's effort
                loglikeffm_rr10=MEASloglikEFFORT2010(MEFF1_10, M_A1_EFF10, M_V1_EFF10,
                                                     MEFF2_10, M_A2_EFF10, M_V2_EFF10,
                                                     MEFF3_10, M_A3_EFF10, M_V3_EFF10,
                                                     MEFF4_10, M_A4_EFF10, M_V4_EFF10,
                                                     MEFF5_10, M_A5_EFF10, M_V5_EFF10,
                                                     MEFF6_10, M_A6_EFF10, M_V6_EFF10,MeffortMean10NWNW_TOG);
                
                //loglikeffm_rr10=F_loglikelihood_generic(Meffort10,MeffortMean10NWNW_TOG,MEASMeffort);
                
                //3. Father's effort
                loglikefff_rr10=MEASloglikEFFORT2010(FEFF1_10, M_A1_EFF10, M_V1_EFF10,
                                                     FEFF2_10, M_A2_EFF10, M_V2_EFF10,
                                                     FEFF3_10, M_A3_EFF10, M_V3_EFF10,
                                                     FEFF4_10, M_A4_EFF10, M_V4_EFF10,
                                                     FEFF5_10, M_A5_EFF10, M_V5_EFF10,
                                                     FEFF6_10, M_A6_EFF10, M_V6_EFF10,FeffortMean10NWNW_TOG);
                
                //loglikefff_rr10=F_loglikelihood_generic(Feffort10,FeffortMean10NWNW_TOG,MEASFeffort);
                
                //4. Investment level
                
                
                loglikeinv_rr10=MEASloginv10(MINV1_10,M_V1_INV10,M_A1_INV10,
                                             MINV2_10,M_V2_INV10,M_A2_INV10,
                                             MINV3_10,M_V3_INV10,M_A3_INV10,
                                             MINV4_10,M_V4_INV10,M_A4_INV10,
                                             MINV5_10,M_V5_INV10,M_A5_INV10,
                                             MINV6_10,M_V6_INV10,M_A6_INV10,
                                             MINV7_10,M_V7_INV10,M_A7_INV10,
                                             MINV8_10,M_V8_INV10,M_A8_INV10,
                                             CInvestmentMean10NWNWA_TOG);
                //Storing the final Decisions
                Feffort10Final=FeffortMean10NWNW_TOG;
                Meffort10Final=MeffortMean10NWNW_TOG;
                Investment10Final=CInvestmentMean10NWNWA_TOG;
            }
            
            
            //7.3 WNWNA
            if (Mfraclabor10==0 && Ffraclabor10==1 && Cchildcare10==0){
                loglikskills_rr10=MEASloglikSkills2010(SKMeasure1, M_A1_S2010,M_VS1_S2010,SKMeasure2, M_A2_S2010,M_VS2_S2010,SKMeasure3, M_A3_S2010,M_VS3_S2010,SKMeasure4, M_A4_S2010,M_VS4_S2010,SKMeasure5, M_A5_S2010,M_VS5_S2010,SKMeasure6, M_A6_S2010,M_VS6_S2010,SKMeasure7, M_A7_S2010,M_VS7_S2010,SKMeasure8, M_A8_S2010,M_VS8_S2010,SKMeasure9, M_A9_S2010,M_VS9_S2010,SKMeasure10, M_A10_S2010,M_VS10_S2010,SKMeasure11, M_A11_S2010,M_VS11_S2010,CSkills_factor10WNWNA_TOG);
                //1. Skills
                //loglikskills_rr10=F_loglikelihood_generic(Ctestsfactor2ss_10,
                //CSkills_factor10WNWNA_TOG,MEASSkills);
                
                //Storing the particle weight and the particle itself
                
                what1_rr_vector[rr]=loglikskills_rr10;
                S1Vector[rr]=CSkills_factor10WNWNA_TOG;
                
                
                //2. Mother's effort
                loglikeffm_rr10=MEASloglikEFFORT2010(MEFF1_10, M_A1_EFF10, M_V1_EFF10,
                                                     MEFF2_10, M_A2_EFF10, M_V2_EFF10,
                                                     MEFF3_10, M_A3_EFF10, M_V3_EFF10,
                                                     MEFF4_10, M_A4_EFF10, M_V4_EFF10,
                                                     MEFF5_10, M_A5_EFF10, M_V5_EFF10,
                                                     MEFF6_10, M_A6_EFF10, M_V6_EFF10,MeffortMean10WNW_TOG);
                
                //loglikeffm_rr10=F_loglikelihood_generic(Meffort10,MeffortMean10WNW_TOG,MEASMeffort);
                
                //3. Father's effort
                loglikefff_rr10=MEASloglikEFFORT2010(FEFF1_10, M_A1_EFF10, M_V1_EFF10,
                                                     FEFF2_10, M_A2_EFF10, M_V2_EFF10,
                                                     FEFF3_10, M_A3_EFF10, M_V3_EFF10,
                                                     FEFF4_10, M_A4_EFF10, M_V4_EFF10,
                                                     FEFF5_10, M_A5_EFF10, M_V5_EFF10,
                                                     FEFF6_10, M_A6_EFF10, M_V6_EFF10,FeffortMean10WNW_TOG);
                //loglikefff_rr10=F_loglikelihood_generic(Feffort10,FeffortMean10WNW_TOG,MEASFeffort);
                
                //4. Investment level
                loglikeinv_rr10=MEASloginv10(MINV1_10,M_V1_INV10,M_A1_INV10,
                                             MINV2_10,M_V2_INV10,M_A2_INV10,
                                             MINV3_10,M_V3_INV10,M_A3_INV10,
                                             MINV4_10,M_V4_INV10,M_A4_INV10,
                                             MINV5_10,M_V5_INV10,M_A5_INV10,
                                             MINV6_10,M_V6_INV10,M_A6_INV10,
                                             MINV7_10,M_V7_INV10,M_A7_INV10,
                                             MINV8_10,M_V8_INV10,M_A8_INV10,
                                             CInvestmentMean10WNWNA_TOG);
                
                //Storing the final Decisions
                Feffort10Final=FeffortMean10NWNW_TOG;
                Meffort10Final=MeffortMean10NWNW_TOG;
                Investment10Final=CInvestmentMean10WNWNA_TOG;
            }
            
            
            //7.4 WNWA
            if (Mfraclabor10==0 && Ffraclabor10==1 && Cchildcare10==1){
                
                
                //1. Skills
                loglikskills_rr10=MEASloglikSkills2010(SKMeasure1, M_A1_S2010,M_VS1_S2010,SKMeasure2, M_A2_S2010,M_VS2_S2010,SKMeasure3, M_A3_S2010,M_VS3_S2010,SKMeasure4, M_A4_S2010,M_VS4_S2010,SKMeasure5, M_A5_S2010,M_VS5_S2010,SKMeasure6, M_A6_S2010,M_VS6_S2010,SKMeasure7, M_A7_S2010,M_VS7_S2010,SKMeasure8, M_A8_S2010,M_VS8_S2010,SKMeasure9, M_A9_S2010,M_VS9_S2010,SKMeasure10, M_A10_S2010,M_VS10_S2010,SKMeasure11, M_A11_S2010,M_VS11_S2010,CSkills_factor10WNWA_TOG);
                //loglikskills_rr10=F_loglikelihood_generic(Ctestsfactor2ss_10,
                //CSkills_factor10WNWA_TOG,MEASSkills);
                
                //Storing the particle weight and the particle itself
                
                what1_rr_vector[rr]=loglikskills_rr10;
                S1Vector[rr]=CSkills_factor10WNWA_TOG;
                
                
                //2. Mother's effort
                
                loglikeffm_rr10=MEASloglikEFFORT2010(MEFF1_10, M_A1_EFF10, M_V1_EFF10,
                                                     MEFF2_10, M_A2_EFF10, M_V2_EFF10,
                                                     MEFF3_10, M_A3_EFF10, M_V3_EFF10,
                                                     MEFF4_10, M_A4_EFF10, M_V4_EFF10,
                                                     MEFF5_10, M_A5_EFF10, M_V5_EFF10,
                                                     MEFF6_10, M_A6_EFF10, M_V6_EFF10,MeffortMean10WNW_TOG);
                
                
                //loglikeffm_rr10=F_loglikelihood_generic(Meffort10, MeffortMean10WNW_TOG,MEASMeffort);
                
                //3. Father's effort
                loglikefff_rr10=MEASloglikEFFORT2010(FEFF1_10, M_A1_EFF10, M_V1_EFF10,
                                                     FEFF2_10, M_A2_EFF10, M_V2_EFF10,
                                                     FEFF3_10, M_A3_EFF10, M_V3_EFF10,
                                                     FEFF4_10, M_A4_EFF10, M_V4_EFF10,
                                                     FEFF5_10, M_A5_EFF10, M_V5_EFF10,
                                                     FEFF6_10, M_A6_EFF10, M_V6_EFF10,FeffortMean10WNW_TOG);
                
                //loglikefff_rr10=F_loglikelihood_generic(Feffort10,FeffortMean10WNW_TOG,MEASFeffort);
                
                //4. Investment level
                loglikeinv_rr10=MEASloginv10(MINV1_10,M_V1_INV10,M_A1_INV10,
                                             MINV2_10,M_V2_INV10,M_A2_INV10,
                                             MINV3_10,M_V3_INV10,M_A3_INV10,
                                             MINV4_10,M_V4_INV10,M_A4_INV10,
                                             MINV5_10,M_V5_INV10,M_A5_INV10,
                                             MINV6_10,M_V6_INV10,M_A6_INV10,
                                             MINV7_10,M_V7_INV10,M_A7_INV10,
                                             MINV8_10,M_V8_INV10,M_A8_INV10,
                                             CInvestmentMean10WNWA_TOG);
                //Storing the final Decisions
                Feffort10Final=FeffortMean10WNW_TOG;
                Meffort10Final=MeffortMean10WNW_TOG;
                Investment10Final=CInvestmentMean10WNWA_TOG;
            }
            
            
            //7.5 NWWNA
            if (Mfraclabor10==1 && Ffraclabor10==0 && Cchildcare10==0){
                //1. Skills
                loglikskills_rr10=MEASloglikSkills2010(SKMeasure1, M_A1_S2010,M_VS1_S2010,SKMeasure2, M_A2_S2010,M_VS2_S2010,SKMeasure3, M_A3_S2010,M_VS3_S2010,SKMeasure4, M_A4_S2010,M_VS4_S2010,SKMeasure5, M_A5_S2010,M_VS5_S2010,SKMeasure6, M_A6_S2010,M_VS6_S2010,SKMeasure7, M_A7_S2010,M_VS7_S2010,SKMeasure8, M_A8_S2010,M_VS8_S2010,SKMeasure9, M_A9_S2010,M_VS9_S2010,SKMeasure10, M_A10_S2010,M_VS10_S2010,SKMeasure11, M_A11_S2010,M_VS11_S2010,CSkills_factor10NWWNA_TOG);
                
                //loglikskills_rr10=F_loglikelihood_generic(Ctestsfactor2ss_10,
                //CSkills_factor10NWWNA_TOG,MEASSkills);
                
                
                //Storing the particle weight and the particle itself
                
                what1_rr_vector[rr]=loglikskills_rr10;
                S1Vector[rr]=CSkills_factor10NWWNA_TOG;
                
                
                //2. Mother's effort
                
                loglikeffm_rr10=MEASloglikEFFORT2010(MEFF1_10, M_A1_EFF10, M_V1_EFF10,
                                                     MEFF2_10, M_A2_EFF10, M_V2_EFF10,
                                                     MEFF3_10, M_A3_EFF10, M_V3_EFF10,
                                                     MEFF4_10, M_A4_EFF10, M_V4_EFF10,
                                                     MEFF5_10, M_A5_EFF10, M_V5_EFF10,
                                                     MEFF6_10, M_A6_EFF10, M_V6_EFF10,MeffortMean10NWW_TOG);
                
                
                //loglikeffm_rr10=F_loglikelihood_generic(Meffort10, MeffortMean10NWW_TOG,MEASMeffort);
                
                //3. Father's effort
                
                loglikefff_rr10=MEASloglikEFFORT2010(FEFF1_10, M_A1_EFF10, M_V1_EFF10,
                                                     FEFF2_10, M_A2_EFF10, M_V2_EFF10,
                                                     FEFF3_10, M_A3_EFF10, M_V3_EFF10,
                                                     FEFF4_10, M_A4_EFF10, M_V4_EFF10,
                                                     FEFF5_10, M_A5_EFF10, M_V5_EFF10,
                                                     FEFF6_10, M_A6_EFF10, M_V6_EFF10,FeffortMean10NWW_TOG);
                
                //loglikefff_rr10=F_loglikelihood_generic(Feffort10,FeffortMean10NWW_TOG,MEASFeffort);
                
                //4. Investment level
                loglikeinv_rr10=MEASloginv10(MINV1_10,M_V1_INV10,M_A1_INV10,
                                             MINV2_10,M_V2_INV10,M_A2_INV10,
                                             MINV3_10,M_V3_INV10,M_A3_INV10,
                                             MINV4_10,M_V4_INV10,M_A4_INV10,
                                             MINV5_10,M_V5_INV10,M_A5_INV10,
                                             MINV6_10,M_V6_INV10,M_A6_INV10,
                                             MINV7_10,M_V7_INV10,M_A7_INV10,
                                             MINV8_10,M_V8_INV10,M_A8_INV10,
                                             CInvestmentMean10NWWNA_TOG);
                //Storing the final Decisions
                Feffort10Final=FeffortMean10NWW_TOG;
                Meffort10Final=MeffortMean10NWW_TOG;
                Investment10Final=CInvestmentMean10NWWNA_TOG;
            }
            
            
            //7.6 NWWA
            if (Mfraclabor10==1 && Ffraclabor10==0 && Cchildcare10==1){
                //1. Skills
                loglikskills_rr10=MEASloglikSkills2010(SKMeasure1, M_A1_S2010,M_VS1_S2010,SKMeasure2, M_A2_S2010,M_VS2_S2010,SKMeasure3, M_A3_S2010,M_VS3_S2010,SKMeasure4, M_A4_S2010,M_VS4_S2010,SKMeasure5, M_A5_S2010,M_VS5_S2010,SKMeasure6, M_A6_S2010,M_VS6_S2010,SKMeasure7, M_A7_S2010,M_VS7_S2010,SKMeasure8, M_A8_S2010,M_VS8_S2010,SKMeasure9, M_A9_S2010,M_VS9_S2010,SKMeasure10, M_A10_S2010,M_VS10_S2010,SKMeasure11, M_A11_S2010,M_VS11_S2010,CSkills_factor10NWWA_TOG);
                
                //loglikskills_rr10=F_loglikelihood_generic(Ctestsfactor2ss_10,
                //CSkills_factor10NWWA_TOG,MEASSkills);
                
                //Storing the particle weight and the particle itself
                
                what1_rr_vector[rr]=loglikskills_rr10;
                S1Vector[rr]=CSkills_factor10NWWA_TOG;
                
                
                //2. Mother's effort
                
                loglikeffm_rr10=MEASloglikEFFORT2010(MEFF1_10, M_A1_EFF10, M_V1_EFF10,
                                                     MEFF2_10, M_A2_EFF10, M_V2_EFF10,
                                                     MEFF3_10, M_A3_EFF10, M_V3_EFF10,
                                                     MEFF4_10, M_A4_EFF10, M_V4_EFF10,
                                                     MEFF5_10, M_A5_EFF10, M_V5_EFF10,
                                                     MEFF6_10, M_A6_EFF10, M_V6_EFF10,MeffortMean10NWW_TOG);
                
                
                //loglikeffm_rr10=F_loglikelihood_generic(Meffort10,MeffortMean10NWW_TOG,MEASMeffort);
                
                //3. Father's effort
                loglikefff_rr10=MEASloglikEFFORT2010(FEFF1_10, M_A1_EFF10, M_V1_EFF10,
                                                     FEFF2_10, M_A2_EFF10, M_V2_EFF10,
                                                     FEFF3_10, M_A3_EFF10, M_V3_EFF10,
                                                     FEFF4_10, M_A4_EFF10, M_V4_EFF10,
                                                     FEFF5_10, M_A5_EFF10, M_V5_EFF10,
                                                     FEFF6_10, M_A6_EFF10, M_V6_EFF10,FeffortMean10NWW_TOG);
                
                
                //loglikefff_rr10=F_loglikelihood_generic(Feffort10,FeffortMean10NWW_TOG,MEASFeffort);
                
                //4. Investment level
                loglikeinv_rr10=MEASloginv10(MINV1_10,M_V1_INV10,M_A1_INV10,
                                             MINV2_10,M_V2_INV10,M_A2_INV10,
                                             MINV3_10,M_V3_INV10,M_A3_INV10,
                                             MINV4_10,M_V4_INV10,M_A4_INV10,
                                             MINV5_10,M_V5_INV10,M_A5_INV10,
                                             MINV6_10,M_V6_INV10,M_A6_INV10,
                                             MINV7_10,M_V7_INV10,M_A7_INV10,
                                             MINV8_10,M_V8_INV10,M_A8_INV10,
                                             CInvestmentMean10NWWA_TOG);
                
                //Storing the final Decisions
                Feffort10Final=FeffortMean10NWW_TOG;
                Meffort10Final=MeffortMean10NWW_TOG;
                Investment10Final=CInvestmentMean10NWWA_TOG;
            }
            //7.7 WWNA
            if (Mfraclabor10==1 && Ffraclabor10==1 && Cchildcare10==0){
                //1. Skills
                loglikskills_rr10=MEASloglikSkills2010(SKMeasure1, M_A1_S2010,M_VS1_S2010,SKMeasure2, M_A2_S2010,M_VS2_S2010,SKMeasure3, M_A3_S2010,M_VS3_S2010,SKMeasure4, M_A4_S2010,M_VS4_S2010,SKMeasure5, M_A5_S2010,M_VS5_S2010,SKMeasure6, M_A6_S2010,M_VS6_S2010,SKMeasure7, M_A7_S2010,M_VS7_S2010,SKMeasure8, M_A8_S2010,M_VS8_S2010,SKMeasure9, M_A9_S2010,M_VS9_S2010,SKMeasure10, M_A10_S2010,M_VS10_S2010,SKMeasure11, M_A11_S2010,M_VS11_S2010,CSkills_factor10WWNA_TOG);
                //loglikskills_rr10=F_loglikelihood_generic(Ctestsfactor2ss_10,
                //CSkills_factor10WWNA_TOG,MEASSkills);
                
                //Storing the particle weight and the particle itself
                
                what1_rr_vector[rr]=loglikskills_rr10;
                S1Vector[rr]=CSkills_factor10WWNA_TOG;
                
                
                
                //2. Mother's effort
                loglikeffm_rr10=MEASloglikEFFORT2010(MEFF1_10, M_A1_EFF10, M_V1_EFF10,
                                                     MEFF2_10, M_A2_EFF10, M_V2_EFF10,
                                                     MEFF3_10, M_A3_EFF10, M_V3_EFF10,
                                                     MEFF4_10, M_A4_EFF10, M_V4_EFF10,
                                                     MEFF5_10, M_A5_EFF10, M_V5_EFF10,
                                                     MEFF6_10, M_A6_EFF10, M_V6_EFF10,MeffortMean10WW_TOG);
                
                //loglikeffm_rr10=F_loglikelihood_generic(Meffort10,MeffortMean10WW_TOG,MEASMeffort);
                
                //3. Father's effort
                loglikefff_rr10=MEASloglikEFFORT2010(FEFF1_10, M_A1_EFF10, M_V1_EFF10,
                                                     FEFF2_10, M_A2_EFF10, M_V2_EFF10,
                                                     FEFF3_10, M_A3_EFF10, M_V3_EFF10,
                                                     FEFF4_10, M_A4_EFF10, M_V4_EFF10,
                                                     FEFF5_10, M_A5_EFF10, M_V5_EFF10,
                                                     FEFF6_10, M_A6_EFF10, M_V6_EFF10,FeffortMean10WW_TOG);
                
                //loglikefff_rr10=F_loglikelihood_generic(Feffort10, FeffortMean10WW_TOG,MEASFeffort);
                
                //4. Investment level
                loglikeinv_rr10=MEASloginv10(MINV1_10,M_V1_INV10,M_A1_INV10,
                                             MINV2_10,M_V2_INV10,M_A2_INV10,
                                             MINV3_10,M_V3_INV10,M_A3_INV10,
                                             MINV4_10,M_V4_INV10,M_A4_INV10,
                                             MINV5_10,M_V5_INV10,M_A5_INV10,
                                             MINV6_10,M_V6_INV10,M_A6_INV10,
                                             MINV7_10,M_V7_INV10,M_A7_INV10,
                                             MINV8_10,M_V8_INV10,M_A8_INV10,
                                             CInvestmentMean10WWNA_TOG);
                //Storing the final Decisions
                Feffort10Final=FeffortMean10WNW_TOG;
                Meffort10Final=MeffortMean10WNW_TOG;
                Investment10Final=CInvestmentMean10WWNA_TOG;
                
            }
            
            //7.7 WWA
            if (Mfraclabor10==1 && Ffraclabor10==1 && Cchildcare10==1){
                loglikskills_rr10=MEASloglikSkills2010(SKMeasure1, M_A1_S2010,M_VS1_S2010,SKMeasure2, M_A2_S2010,M_VS2_S2010,SKMeasure3, M_A3_S2010,M_VS3_S2010,SKMeasure4, M_A4_S2010,M_VS4_S2010,SKMeasure5, M_A5_S2010,M_VS5_S2010,SKMeasure6, M_A6_S2010,M_VS6_S2010,SKMeasure7, M_A7_S2010,M_VS7_S2010,SKMeasure8, M_A8_S2010,M_VS8_S2010,SKMeasure9, M_A9_S2010,M_VS9_S2010,SKMeasure10, M_A10_S2010,M_VS10_S2010,SKMeasure11, M_A11_S2010,M_VS11_S2010,CSkills_factor10WWA_TOG);
                //1. Skills
                // loglikskills_rr10=F_loglikelihood_generic(Ctestsfactor2ss_10,
                //CSkills_factor10WWA_TOG,MEASSkills);
                
                //Storing the particle weight and the particle itself
                
                what1_rr_vector[rr]=loglikskills_rr10;
                S1Vector[rr]=CSkills_factor10WWA_TOG;
                
                
                
                //2. Mother's effort
                
                loglikeffm_rr10=MEASloglikEFFORT2010(MEFF1_10, M_A1_EFF10, M_V1_EFF10,
                                                     MEFF2_10, M_A2_EFF10, M_V2_EFF10,
                                                     MEFF3_10, M_A3_EFF10, M_V3_EFF10,
                                                     MEFF4_10, M_A4_EFF10, M_V4_EFF10,
                                                     MEFF5_10, M_A5_EFF10, M_V5_EFF10,
                                                     MEFF6_10, M_A6_EFF10, M_V6_EFF10,MeffortMean10WW_TOG);
                
                //loglikeffm_rr10=F_loglikelihood_generic(Meffort10,MeffortMean10WW_TOG,MEASMeffort);
                
                //3. Father's effort
                loglikefff_rr10=MEASloglikEFFORT2010(FEFF1_10, M_A1_EFF10, M_V1_EFF10,
                                                     FEFF2_10, M_A2_EFF10, M_V2_EFF10,
                                                     FEFF3_10, M_A3_EFF10, M_V3_EFF10,
                                                     FEFF4_10, M_A4_EFF10, M_V4_EFF10,
                                                     FEFF5_10, M_A5_EFF10, M_V5_EFF10,
                                                     FEFF6_10, M_A6_EFF10, M_V6_EFF10,FeffortMean10WW_TOG);
                
                //loglikefff_rr10=F_loglikelihood_generic(Feffort10,FeffortMean10WW_TOG,MEASFeffort);
                
                //4. Investment level
                loglikeinv_rr10=MEASloginv10(MINV1_10,M_V1_INV10,M_A1_INV10,
                                             MINV2_10,M_V2_INV10,M_A2_INV10,
                                             MINV3_10,M_V3_INV10,M_A3_INV10,
                                             MINV4_10,M_V4_INV10,M_A4_INV10,
                                             MINV5_10,M_V5_INV10,M_A5_INV10,
                                             MINV6_10,M_V6_INV10,M_A6_INV10,
                                             MINV7_10,M_V7_INV10,M_A7_INV10,
                                             MINV8_10,M_V8_INV10,M_A8_INV10,
                                             CInvestmentMean10WWA_TOG);
                //Storing the final Decisions
                Feffort10Final=FeffortMean10WW_TOG;
                Meffort10Final=MeffortMean10WW_TOG;
                Investment10Final=CInvestmentMean10WWA_TOG;
                
            }
            
            
            //5. Measurement system for mmu. No in 2010
            
        }//If living together
        
        
        
        //Computing the sum of the weights
        
        
        
        
        
        //0 Decisions (behavioral system)
        logcontribdecision10+=Decision10; //Storing the number of times correct thing
        logcontribobdecision10SMOOTHED+=VCONTRIB10;
        
        
        //2. Mother's effort
        loglikeffm10+=loglikeffm_rr10;
        //3. Father's effort
        loglikefff10+=loglikefff_rr10;
        //4. Investment
        loglikeinv10+=loglikeinv_rr10;
        
        //5. Mmu
        
        rr=rr+1;
    } //Finish the RR simulations
    
    //Storing the weights that are NOT normalized and that will give the likelihood
    wtilde1_rr_vector=what1_rr_vector;
    
    //Block of normalizing the vector what1_rr_vector
    AuxiliarV=0;
    for(int mm=0; mm<RR; mm=mm+1){
        AuxiliarV+=what1_rr_vector[mm];
    }
    for(int mm=0; mm<RR; mm=mm+1){
        what1_rr_vector[mm]=what1_rr_vector[mm]/AuxiliarV;
    }
    AuxiliarV=0;
    //==========new block of taking out particles randomly (START)==========//
    boost::random::discrete_distribution<int> distskills1(what1_rr_vector);
    iterator=0;
    wsum1=0;
    for (int mm=0; mm<RR;mm=mm+1){
        iterator=distskills1(rng);
        S1Vectoraux[mm]=S1Vector[iterator];
        w1_rr_vectorAUX[mm]=what1_rr_vector[iterator];
        wsum1+=wtilde1_rr_vector[iterator]; //Liklihood is adding up the weight (original weight) of those chosen (iterator)
        
        //Re-adapting the skills in t=0. This is done for the smoothing distribution
        S0Vectoraux[mm]=S0Vector[iterator];
        w0_rr_vectorAUX[mm]=w0_rr_vector[iterator];
    }
    //Updating the weights and the skills
    
    //Redefining vector of skills
    S1Vector=S1Vectoraux;
    what1_rr_vector=w1_rr_vectorAUX;
    w1_rr_vector=what1_rr_vector;
    
    S0Vector=S0Vectoraux;
    w0_rr_vector=w0_rr_vectorAUX;
    
    if(1==2){
        for(int mm=0; mm<RR; mm=mm+1){
            cout << S1Vector[mm] << " S1Vector[rr] re-shuffled" << endl;
        }
    }
    
    
    //========================================================
    //2012 simulations
    //========================================================
    rr=0;
    wsum2=0;
    while (rr<RR){
        logcontribobs_rr=0;
        loglikskills_rr=0;
        loglikeffm_rr=0;
        loglikefff_rr=0;
        loglikmmu_rr=0;
        loglikeinv_rr=0;
        Decision=0;
        logcontribdecision_rr=0;
        
        //1. Draw shocks from the corresponding distributions
        //===========================================================
        //1.1 Single mothers:
        if (Cliveswithfather12==0 & Cliveswithmother12==1){
            
            //========================================================
            //1. Effort levels
            //========================================================
            //Draw a shock for the effort levels
            Meffort_factor12shock=gen_normal_3(generator);
            //========================================================
            //1.1 If mom does'nt work
            //========================================================
            // Getting a draw from a standard normal distribution
            
            // Now transforming it as a draw from the given distribution
            Meffort_factor12NW=Meffort_factor12shock*stdeffortmom;
            Meffort_factor12NW=exp(Meffort_factor12NW)*MeffortMean12NW;
            //Effort can't be negative, we are drawing from a truncated normal distribution
            if (Meffort_factor12NW<0){
                Meffort_factor12NW=0;
            }
            //========================================================
            //1.2 If mom works
            // Getting a draw from a standard normal distribution
            
            // Now transforming it as a draw from the given distribution
            Meffort_factor12W=Meffort_factor12shock*stdeffortmom;
            Meffort_factor12W=exp(Meffort_factor12W)*MeffortMean12W;
            
            //Effort can't be negative, we are drawing from a truncated normal distribution
            if (Meffort_factor12W<0){
                Meffort_factor12W=0;
            }
            //========================================================
            
            
            //========================================================
            //2. Investment levels
            //========================================================
            //Getting the draw for the investment shock
            CInv_factor12shock=gen_normal_3(generator);
            //========================================================
            //2.1 If mom doesn't work
            //========================================================
            //3. Now transforming it as a draw from the given distribution
            CInv_factor12NW=CInv_factor12shock*stdinvestment;
            CInv_factor12NW=CInv_factor12NW+CInvestmentMean12NW;
            
            //Investment can't be negative
            if (CInv_factor12NW<0){
                CInv_factor12NW=0;
            }
            //========================================================
            //2.2 If mom works
            //========================================================
            //3. Now transforming it as a draw from the given distribution
            CInv_factor12W=CInv_factor12shock*stdinvestment;
            CInv_factor12W=CInv_factor12W+CInvestmentMean12W;
            //Investment can't be negative
            if (CInv_factor12W<0){
                CInv_factor12W=0;
            }
            
            //==========================================================
            //3. Skills levels
            //==========================================================
            //Shock of skills
            CSkills_factor12shock=gen_normal_3(generator);
            //========================================================
            //1. If mom doesn't works
            //========================================================
            //3. Now transforming it as a draw from the given distribution
            CSkillsMean12NW=F_predskills(ddelta0,ddelta1,ddelta2,ddelta3_12,ddelta4,
                                         Cedad_meses12,ttheta0,
                                         ttheta1,ttheta2, pphi,ggammaf, ggammam,
                                         Feffort12,Meffort_factor12NW,
                                         CInv_factor12NW,S1Vector[rr],
                                         Cchildcare12,PGVector[rr],Hmemberstotal12);
            
            if (CSkills_factor12shock*stdskills<=500){
                CSkills_factor12NW=CSkills_factor12shock*stdskills;
            }
            if (CSkills_factor12shock*stdskills>=500){
                CSkills_factor12NW=500;
            }
            
            
            CSkills_factor12NW=exp(CSkills_factor12NW)*CSkillsMean12NW;
            //========================================================
            //2. If mom works
            //========================================================
            //3. Now transforming it as a draw from the given distribution
            CSkillsMean12W=F_predskills(ddelta0,ddelta1,ddelta2,ddelta3_12,ddelta4,
                                        Cedad_meses12,ttheta0,
                                        ttheta1,ttheta2, pphi,ggammaf, ggammam,
                                        Feffort12,Meffort_factor12W,
                                        CInv_factor12W,
                                        S1Vector[rr],
                                        Cchildcare12,
                                        PGVector[rr],Hmemberstotal12);
            
            if (CSkills_factor12shock*stdskills<=500){
                CSkills_factor12W=CSkills_factor12shock*stdskills;
            }
            if (CSkills_factor12shock*stdskills>=500){
                CSkills_factor12W=500;
            }
            
            
            CSkills_factor12W=exp(CSkills_factor12W)*CSkillsMean12W;
            
            //===========================================================
            //===========================================================
            //4 Preference shocks
            //===========================================================
            //1. If mom doesn't work
            //===========================================================
            MprefshockNW=gen_normal_3(generator);
            MprefshockNW=MprefshockW*MshockNWA;
            //===========================================================
            //2. If mom works
            //===========================================================
            MprefshockW=gen_normal_3(generator);
            MprefshockW=MprefshockW*MshockWA;
            //===========================================================
            //===========================================================
            //5 Consumption
            //===========================================================
            //1. If mom doesn't work
            //===========================================================
            MconsumptionNW=Mnly12-priceINV*CInv_factor12NW;
            if (MconsumptionNW<0){
                MconsumptionNW=1.0e-10;
            }
            //===========================================================
            //2. If mom works
            //===========================================================
            MconsumptionW=Mnly12-priceINV*CInv_factor12NW+Mwageinlik;
            if (MconsumptionW<0){
                MconsumptionW=1.0e-10;
            }
            
            
            
            //===========================================================
            //===========================================================
            //6. Behavioral model
            //==========================================================
            //1. If mom doesn't work
            //==========================================================
            MutilityNW=F_utility(aalpha1m,aalpha2m,aalpha3m,aalpha4m,MconsumptionNW,
                                 Meffort_factor12NW,0,CSkills_factor12NW)+MprefshockNW;
            
            //===========================================================
            //2. If mom works
            //===========================================================
            MutilityW=F_utility(aalpha1m,aalpha2m,aalpha3m,aalpha4m,MconsumptionW,
                                Meffort_factor12W,1,CSkills_factor12W)+MprefshockW;
            
            //===========================================================
            //===========================================================
            //Computing the actual likelihood
            
            Decision=(MutilityW>MutilityNW);
            //We want to call the decision=1 if the decission is correct. Then if person
            //decides not to work the thing will be 1-that.
            Decision=(obDecision==Decision);
            
            //Likelihood contributions of the measurement systems
            if (Mfraclabor12==1){
                
                //1. Skills
                
                loglikskills_rr= MEASloglikSkills2012(SK2012Measure1, M_A1_S2012,M_VS1_S2012,SK2012Measure2, M_A2_S2012,M_VS2_S2012,SK2012Measure3, M_A3_S2012,M_VS3_S2012,SK2012Measure4, M_A4_S2012,M_VS4_S2012,SK2012Measure5, M_A5_S2012,M_VS5_S2012,SK2012Measure6, M_A6_S2012,M_VS6_S2012,SK2012Measure7, M_A7_S2012,M_VS7_S2012,SK2012Measure8, M_A8_S2012,M_VS8_S2012,SK2012Measure9, M_A9_S2012,M_VS9_S2012,SK2012Measure10, M_A10_S2012,M_VS10_S2012,SK2012Measure11, M_A11_S2012,M_VS11_S2012,SK2012Measure12,M_A12_S2012,M_VS12_S2012,SK2012Measure13, M_A13_S2012, M_VS13_S2012,CSkills_factor12W);
                
                
                
                //loglikskills_rr=F_loglikelihood_generic(Ctestsfactorsss2012,
                //                                        CSkills_factor12W,MEASSkills);
                
                //Storing the particle weight and the particle itself
                what2_rr_vector[rr]=loglikskills_rr;
                S2Vector[rr]=CSkills_factor12W;
                
                //2. Effort
                loglikeffm_rr=MEASloglikFEFFORT2012(MEFF1_12, M_A1_EFF12, M_V1_EFF12,
                                                    MEFF2_12, M_A2_EFF12, M_V2_EFF12,
                                                    MEFF3_12, M_A3_EFF12, M_V3_EFF12,
                                                    MEFF4_12, M_A4_EFF12, M_V4_EFF12,
                                                    MEFF5_12, M_A5_EFF12, M_V5_EFF12,
                                                    MEFF6_12, M_A6_EFF12, M_V6_EFF12,
                                                    MEFF7_12, M_A7_EFF12, M_V7_EFF12,
                                                    MEFF8_12, M_A8_EFF12, M_V8_EFF12,
                                                    MEFF9_12, M_A9_EFF12, M_V9_EFF12,
                                                    MEFF10_12, M_A10_EFF12, M_V10_EFF12,
                                                    MEFF11_12, M_A11_EFF12, M_V11_EFF12,
                                                    MEFF12_12, M_A12_EFF12, M_V12_EFF12,
                                                    MEFF13_12, M_A13_EFF12, M_V13_EFF12,
                                                    MEFF14_12, M_A14_EFF12, M_V14_EFF12, Meffort_factor12W);
                //loglikeffm_rr=F_loglikelihood_generic(Meffort12,Meffort_factor12W,MEASMeffort);
                loglikefff_rr=0;
                
                //3. Investment
                loglikeinv_rr=MEASloginv12(MINV1_12,M_V1_INV12,M_A1_INV12,
                                           MINV2_12,M_V2_INV12,M_A2_INV12,
                                           MINV3_12,M_V3_INV12,M_A3_INV12,
                                           MINV4_12,M_V4_INV12,M_A4_INV12,
                                           MINV5_12,M_V5_INV12,M_A5_INV12,
                                           MINV6_12,M_V6_INV12,M_A6_INV12,
                                           MINV7_12,M_V7_INV12,M_A7_INV12,
                                           MINV8_12,M_V8_INV12,M_A8_INV12,
                                           MINV9_12,M_V9_INV12,M_A9_INV12,
                                           MINV10_12,M_V10_INV12,M_A10_INV12,
                                           MINV11_12,M_V11_INV12,M_A11_INV12,
                                           MINV12_12,M_V12_INV12,M_A12_INV12,
                                           MINV13_12,M_V13_INV12,M_A13_INV12,
                                           MINV14_12,M_V14_INV12,M_A14_INV12,
                                           MINV15_12,M_V15_INV12,M_A15_INV12,
                                           MINV16_12,M_V16_INV12,M_A16_INV12,
                                           MINV17_12,M_V17_INV12,M_A17_INV12,
                                           MINV18_12,M_V18_INV12,M_A18_INV12,
                                           MINV19_12,M_V19_INV12,M_A19_INV12,
                                           MINV20_12,M_V20_INV12,M_A20_INV12,
                                           MINV21_12,M_V21_INV12,M_A21_INV12,
                                           CInv_factor12W);
                
                //Storing the final Decisions
                Feffort12Final=Feffort12;
                Meffort12Final=Meffort_factor12W;
                Investment12Final=CInv_factor12W;
                
            }
            if (Mfraclabor12==0){
                
                //1. Skills
                
                loglikskills_rr= MEASloglikSkills2012(SK2012Measure1, M_A1_S2012,M_VS1_S2012,SK2012Measure2, M_A2_S2012,M_VS2_S2012,SK2012Measure3, M_A3_S2012,M_VS3_S2012,SK2012Measure4, M_A4_S2012,M_VS4_S2012,SK2012Measure5, M_A5_S2012,M_VS5_S2012,SK2012Measure6, M_A6_S2012,M_VS6_S2012,SK2012Measure7, M_A7_S2012,M_VS7_S2012,SK2012Measure8, M_A8_S2012,M_VS8_S2012,SK2012Measure9, M_A9_S2012,M_VS9_S2012,SK2012Measure10, M_A10_S2012,M_VS10_S2012,SK2012Measure11, M_A11_S2012,M_VS11_S2012,SK2012Measure12,M_A12_S2012,M_VS12_S2012,SK2012Measure13, M_A13_S2012, M_VS13_S2012,CSkills_factor12NW);
                
                
                //loglikskills_rr=F_loglikelihood_generic(Ctestsfactorsss2012,
                //                                        CSkills_factor12NW,MEASSkills);
                
                //Storing the particle weight and the particle itself
                what2_rr_vector[rr]=loglikskills_rr;
                S2Vector[rr]=CSkills_factor12W;
                
                
                //2. Effort
                
                loglikeffm_rr=MEASloglikFEFFORT2012(MEFF1_12, M_A1_EFF12, M_V1_EFF12,
                                                    MEFF2_12, M_A2_EFF12, M_V2_EFF12,
                                                    MEFF3_12, M_A3_EFF12, M_V3_EFF12,
                                                    MEFF4_12, M_A4_EFF12, M_V4_EFF12,
                                                    MEFF5_12, M_A5_EFF12, M_V5_EFF12,
                                                    MEFF6_12, M_A6_EFF12, M_V6_EFF12,
                                                    MEFF7_12, M_A7_EFF12, M_V7_EFF12,
                                                    MEFF8_12, M_A8_EFF12, M_V8_EFF12,
                                                    MEFF9_12, M_A9_EFF12, M_V9_EFF12,
                                                    MEFF10_12, M_A10_EFF12, M_V10_EFF12,
                                                    MEFF11_12, M_A11_EFF12, M_V11_EFF12,
                                                    MEFF12_12, M_A12_EFF12, M_V12_EFF12,
                                                    MEFF13_12, M_A13_EFF12, M_V13_EFF12,
                                                    MEFF14_12, M_A14_EFF12, M_V14_EFF12, Meffort_factor12NW);
                
                
                //loglikeffm_rr=F_loglikelihood_generic(Meffort12,Meffort_factor12NW,MEASMeffort);
                loglikefff_rr=0;
                
                //3. Investment
                loglikeinv_rr=MEASloginv12(MINV1_12,M_V1_INV12,M_A1_INV12,
                                           MINV2_12,M_V2_INV12,M_A2_INV12,
                                           MINV3_12,M_V3_INV12,M_A3_INV12,
                                           MINV4_12,M_V4_INV12,M_A4_INV12,
                                           MINV5_12,M_V5_INV12,M_A5_INV12,
                                           MINV6_12,M_V6_INV12,M_A6_INV12,
                                           MINV7_12,M_V7_INV12,M_A7_INV12,
                                           MINV8_12,M_V8_INV12,M_A8_INV12,
                                           MINV9_12,M_V9_INV12,M_A9_INV12,
                                           MINV10_12,M_V10_INV12,M_A10_INV12,
                                           MINV11_12,M_V11_INV12,M_A11_INV12,
                                           MINV12_12,M_V12_INV12,M_A12_INV12,
                                           MINV13_12,M_V13_INV12,M_A13_INV12,
                                           MINV14_12,M_V14_INV12,M_A14_INV12,
                                           MINV15_12,M_V15_INV12,M_A15_INV12,
                                           MINV16_12,M_V16_INV12,M_A16_INV12,
                                           MINV17_12,M_V17_INV12,M_A17_INV12,
                                           MINV18_12,M_V18_INV12,M_A18_INV12,
                                           MINV19_12,M_V19_INV12,M_A19_INV12,
                                           MINV20_12,M_V20_INV12,M_A20_INV12,
                                           MINV21_12,M_V21_INV12,M_A21_INV12,
                                           CInv_factor12NW);
                //Storing the final Decisions
                Feffort12Final=Feffort12;
                Meffort12Final=Meffort_factor12NW;
                Investment12Final=CInv_factor12NW;
                
            }
            
            
            
        }//End of single mothers
        //If both parents are present.
        
        if (Cliveswithfather12==1 & Cliveswithmother12==1){
            
            //======================
            //0. Mmu
            //======================
            //Getting the draw for the unobserved heterogeneity of mmy
            MMushock=gen_normal_3(generator);
            MMushock=MMushock*stdmmu;
            //Getting the predicted level of bargaining power
            mupred=F_mmu(llambda0,llambda1,llambda2,llambda3,llambda4,
                         llambda5,llambda6,llambda7,llambda8,
                         Fwageinlik,Mwageinlik,
                         Fnly12,Mnly12,MMushock,mmuLB,mmuUB,
                         Fage12,Mage12,Fyrschool12,Myrschool12,
                         MRATIO,Unemployment,Wageratio,Distance);
            
            
            
            //========================================================
            //1. Effort levels
            //========================================================
            
            //Draw a shock for the effort levels
            Meffort_factor12shock=exp(gen_normal_3(generator)*stdeffortmom);
            Feffort_factor12shock=exp(gen_normal_3(generator)*stdeffortfat);
            
            //Getting the means
            //1.1 Effort levels
            //1.1.1 NWNW
            MeffortMean12NWNW_TOG=F_effort_m(mupred,ggammaf,aalpha2f,aalpha2m,
                                             aalpha4f,aalpha4m,0,0,pphi,
                                             ttheta2)*Meffort_factor12shock;
            MEFFORTNWNW12NOSHOCK=F_effort_m(mupred,ggammaf,aalpha2f,aalpha2m,
                                            aalpha4f,aalpha4m,0,0,pphi,
                                            ttheta2);
            
            FeffortMean12NWNW_TOG=F_effort_f(mupred,ggammaf,aalpha2f,aalpha2m,
                                             aalpha4f,aalpha4m,0,0,pphi,
                                             ttheta2)*Feffort_factor12shock;
            FEFFORTNWNW12NOSHOCK=F_effort_f(mupred,ggammaf,aalpha2f,aalpha2m,
                                            aalpha4f,aalpha4m,0,0,pphi,
                                            ttheta2);
            //1.1.2 Father works mother not
            MeffortMean12WNW_TOG=F_effort_m(mupred,ggammaf,aalpha2f,aalpha2m,
                                            aalpha4f,aalpha4m,1,0,pphi,
                                            ttheta2)*Meffort_factor12shock;
            MEFFORTWNW12NOSHOCK=F_effort_m(mupred,ggammaf,aalpha2f,aalpha2m,
                                           aalpha4f,aalpha4m,1,0,pphi,
                                           ttheta2);
            FeffortMean12WNW_TOG=F_effort_f(mupred,ggammaf,aalpha2f,aalpha2m,
                                            aalpha4f,aalpha4m,1,0,pphi,
                                            ttheta2)*Feffort_factor12shock;
            FEFFORTWNW12NOSHOCK=F_effort_f(mupred,ggammaf,aalpha2f,aalpha2m,
                                           aalpha4f,aalpha4m,1,0,pphi,
                                           ttheta2);
            //1.1.3 Mother works father doesn't
            MeffortMean12NWW_TOG=F_effort_m(mupred,ggammaf,aalpha2f,aalpha2m,
                                            aalpha4f,aalpha4m,0,1,pphi,
                                            ttheta2)*Meffort_factor12shock;
            MEFFORTNWW12NOSHOCK=F_effort_m(mupred,ggammaf,aalpha2f,aalpha2m,
                                           aalpha4f,aalpha4m,0,1,pphi,
                                           ttheta2);
            
            FeffortMean12NWW_TOG=F_effort_f(mupred,ggammaf,aalpha2f,aalpha2m,
                                            aalpha4f,aalpha4m,0,1,pphi,
                                            ttheta2)*Feffort_factor12shock;
            FEFFORTNWW12NOSHOCK=F_effort_f(mupred,ggammaf,aalpha2f,aalpha2m,
                                           aalpha4f,aalpha4m,0,1,pphi,
                                           ttheta2);
            //1.1.4 Both work
            MeffortMean12WW_TOG=F_effort_m(mupred,ggammaf,aalpha2f,aalpha2m,
                                           aalpha4f,aalpha4m,1,1,pphi,
                                           ttheta2)*Meffort_factor12shock;
            MEFFORTWW12NOSHOCK=F_effort_m(mupred,ggammaf,aalpha2f,aalpha2m,
                                          aalpha4f,aalpha4m,1,1,pphi,
                                          ttheta2);
            FeffortMean12WW_TOG=F_effort_f(mupred,ggammaf,aalpha2f,aalpha2m,
                                           aalpha4f,aalpha4m,1,1,pphi,
                                           ttheta2)*Feffort_factor12shock;
            FEFFORTWW12NOSHOCK=F_effort_f(mupred,ggammaf,aalpha2f,aalpha2m,
                                          aalpha4f,aalpha4m,1,1,pphi,
                                          ttheta2);
            
            //========================================================
            //2. Investment levels
            //========================================================
            //Getting the draw for the investment shock
            CInv_factor12shock=gen_normal_3(generator);
            //========================================================
            //2.1 If mom doesn't work
            //========================================================
            //3. Now transforming it as a draw from the given distribution
            CInv_factor12shock=CInv_factor12shock*stdinvestment;
            //4. Now get the mean of the investment as a function of mmu
            //for the possible combinations
            
            CInvestmentMean12NWNW_TOG=F_invcouple(aalpha1m,aalpha1f,
                                                  aalpha2m,aalpha2f,
                                                  mupred,0,0,
                                                  Fwageinlik,Mwageinlik,
                                                  Fnly12,Mnly12,
                                                  ttheta1,priceINV)*exp(CInv_factor12shock);
            CINVFACTOR12NWNW_noshock=F_invcouple(aalpha1m,aalpha1f,
                                                 aalpha2m,aalpha2f,
                                                 mupred,0,0,
                                                 Fwageinlik,Mwageinlik,
                                                 Fnly12,Mnly12,
                                                 ttheta1,priceINV);
            if (CInvestmentMean12NWNW_TOG<=0){
                CInvestmentMean12NWNW_TOG=1.0e-3;
            }
            CInvestmentMean12WNW_TOG=F_invcouple(aalpha1m,aalpha1f,
                                                 aalpha2m,aalpha2f,
                                                 mupred,1,0,
                                                 Fwageinlik,Mwageinlik,
                                                 Fnly12,Mnly12,
                                                 ttheta1,priceINV)*exp(CInv_factor12shock);
            CINVFACTOR12WNW_noshock=F_invcouple(aalpha1m,aalpha1f,
                                                aalpha2m,aalpha2f,
                                                mupred,1,0,
                                                Fwageinlik,Mwageinlik,
                                                Fnly12,Mnly12,
                                                ttheta1,priceINV);
            if (CInvestmentMean12WNW_TOG<=0){
                CInvestmentMean12WNW_TOG=1.0e-3;
            }
            CInvestmentMean12NWW_TOG=F_invcouple(aalpha1m,aalpha1f,
                                                 aalpha2m,aalpha2f,
                                                 mupred,0,1,
                                                 Fwageinlik,Mwageinlik,
                                                 Fnly12,Mnly12,
                                                 ttheta1,priceINV)*exp(CInv_factor12shock);
            CINVFACTOR12NWW_noshock=F_invcouple(aalpha1m,aalpha1f,
                                                aalpha2m,aalpha2f,
                                                mupred,0,1,
                                                Fwageinlik,Mwageinlik,
                                                Fnly12,Mnly12,
                                                ttheta1,priceINV);
            if (CInvestmentMean12NWW_TOG<=0){
                CInvestmentMean12NWW_TOG=1.0e-3;
            }
            CInvestmentMean12WW_TOG=F_invcouple(aalpha1m,aalpha1f,
                                                aalpha2m,aalpha2f,
                                                mupred,1,1,
                                                Fwageinlik,Mwageinlik,
                                                Fnly12,Mnly12,
                                                ttheta1,priceINV)*exp(CInv_factor12shock);
            CINVFACTOR12WW_noshock=F_invcouple(aalpha1m,aalpha1f,
                                               aalpha2m,aalpha2f,
                                               mupred,1,1,
                                               Fwageinlik,Mwageinlik,
                                               Fnly12,Mnly12,
                                               ttheta1,priceINV);
            if (CInvestmentMean12WW_TOG<=0){
                CInvestmentMean12WW_TOG=1.0e-3;
            }
            
            //==========================================================
            //3. Skills levels
            //==========================================================
            //Shock of skills
            CSkills_factor12shock=gen_normal_3(generator);
            if (CSkills_factor12shock*stdskills<=500){
                CSkills_factor12shock=exp(CSkills_factor12shock*stdskills);
            }
            if (CSkills_factor12shock*stdskills>500){
                CSkills_factor12shock=exp(500);
            }
            
            
            //Now getting the different levels of skills for all possi-
            //bilities
            //3.1 NWNW
            CSkills_factor12NWNW_TOG=F_predskills(ddelta0,ddelta1,ddelta2,ddelta3_12,ddelta4,
                                                  Cedad_meses12,
                                                  ttheta0,ttheta1,ttheta2,pphi,
                                                  ggammaf,ggammam,FEFFORTNWNW12NOSHOCK,
                                                  MEFFORTNWNW12NOSHOCK,
                                                  CINVFACTOR12NWNW_noshock,
                                                  S1Vector[rr],
                                                  Cchildcare12,
                                                  PGVector[rr],Hmemberstotal12)*CSkills_factor12shock;
            
            //3.2 WNW
            CSkills_factor12WNW_TOG=F_predskills(ddelta0,ddelta1,ddelta2,ddelta3_12,ddelta4,
                                                 Cedad_meses12,
                                                 ttheta0,ttheta1,ttheta2,pphi,
                                                 ggammaf,ggammam,FEFFORTWNW12NOSHOCK,
                                                 MEFFORTWNW12NOSHOCK,
                                                 CINVFACTOR12WNW_noshock,
                                                 S1Vector[rr],
                                                 Cchildcare12,
                                                 PGVector[rr],Hmemberstotal12)*CSkills_factor12shock;
            
            //3.3 NWW
            CSkills_factor12NWW_TOG=F_predskills(ddelta0,ddelta1,ddelta2,ddelta3_12,ddelta4,
                                                 Cedad_meses12,
                                                 ttheta0,ttheta1,ttheta2,pphi,
                                                 ggammaf,ggammam,FEFFORTNWW12NOSHOCK,
                                                 MEFFORTNWW12NOSHOCK,
                                                 CINVFACTOR12NWW_noshock,
                                                 S1Vector[rr],Cchildcare12,
                                                 PGVector[rr],Hmemberstotal12)*CSkills_factor12shock;
            
            //3.4 WW
            CSkills_factor12WW_TOG=F_predskills(ddelta0,ddelta1,ddelta2,ddelta3_12,ddelta4,
                                                Cedad_meses12,
                                                ttheta0,ttheta1,ttheta2,pphi,
                                                ggammaf,ggammam,FEFFORTWW12NOSHOCK,
                                                MEFFORTWW12NOSHOCK,
                                                CINVFACTOR12WW_noshock,
                                                S1Vector[rr],Cchildcare12,
                                                PGVector[rr],Hmemberstotal12)*CSkills_factor12shock;
            
            //===========================================================
            //===========================================================
            //4 Preference shocks
            
            //1. If mom doesn't work
            MprefshockNW=gen_normal_3(generator);
            MprefshockNW=MprefshockW*MshockNWA;
            //2. If mom works
            MprefshockW=gen_normal_3(generator);
            MprefshockW=MprefshockW*MshockWA;
            //3. If FATHER doesn't work
            FprefshockNW=gen_normal_3(generator);
            FprefshockNW=FprefshockW*FshockNWA;
            //4. If FATHER works
            FprefshockW=gen_normal_3(generator);
            FprefshockW=FprefshockW*FshockWA;
            
            
            
            //===========================================================
            //===========================================================
            //5 Consumption
            //===========================================================
            //1. Mother
            //===========================================================
            MconsumptionNWNW_TOG=M_consumption_TOG(aalpha1m,aalpha2f,aalpha2m,
                                                   CInvestmentMean12NWNW_TOG,
                                                   ttheta1,mupred,priceINV);
            
            
            MconsumptionWNW_TOG=M_consumption_TOG(aalpha1m,aalpha2f,aalpha2m,
                                                  CInvestmentMean12WNW_TOG,
                                                  ttheta1,mupred,priceINV);
            
            MconsumptionNWW_TOG=M_consumption_TOG(aalpha1m,aalpha2f,aalpha2m,
                                                  CInvestmentMean12NWW_TOG,
                                                  ttheta1,mupred,priceINV);
            
            MconsumptionWW_TOG=M_consumption_TOG(aalpha1m,aalpha2f,aalpha2m,
                                                 CInvestmentMean12WW_TOG,
                                                 ttheta1,mupred,priceINV);
            
            //===========================================================
            //1. Father
            //===========================================================
            FconsumptionNWNW_TOG=F_consumption_TOG(aalpha1f,aalpha2f,aalpha2m,
                                                   CInvestmentMean12NWNW_TOG,
                                                   ttheta1,mupred,priceINV);
            
            FconsumptionWNW_TOG=F_consumption_TOG(aalpha1f,aalpha2f,aalpha2m,
                                                  CInvestmentMean12WNW_TOG,
                                                  ttheta1,mupred,priceINV);
            
            FconsumptionNWW_TOG=F_consumption_TOG(aalpha1f,aalpha2f,aalpha2m,
                                                  CInvestmentMean12NWW_TOG,
                                                  ttheta1,mupred,priceINV);
            
            FconsumptionWW_TOG=F_consumption_TOG(aalpha1f,aalpha2f,aalpha2m,
                                                 CInvestmentMean12WW_TOG,
                                                 ttheta1,mupred,priceINV);
            
            
            
            
            //===========================================================
            //===========================================================
            //6. Behavioral model
            //==========================================================
            //1. Mother
            //==========================================================
            
            MutilityNWNW_TOG=F_utility(aalpha1m,aalpha2m,aalpha3m,aalpha4m,
                                       MconsumptionNWNW_TOG,MeffortMean12NWNW_TOG,0,CSkills_factor12NWNW_TOG)+
            MprefshockNW;
            
            MutilityWNW_TOG=F_utility(aalpha1m,aalpha2m,aalpha3m,aalpha4m,
                                      MconsumptionWNW_TOG,MeffortMean12WNW_TOG,0,CSkills_factor12WNW_TOG)+MprefshockNW;
            
            MutilityNWW_TOG=F_utility(aalpha1m,aalpha2m,aalpha3m,aalpha4m,
                                      MconsumptionNWW_TOG,MeffortMean12NWW_TOG,1,CSkills_factor12NWW_TOG)+MprefshockW;
            
            MutilityWW_TOG=F_utility(aalpha1m,aalpha2m,aalpha3m,aalpha4m,
                                     MconsumptionWW_TOG,MeffortMean12WW_TOG,1,CSkills_factor12WW_TOG)+
            MprefshockW;
            
            //==========================================================
            //2. Father
            //==========================================================
            
            FutilityNWNW_TOG=F_utility(aalpha1f,aalpha2f,aalpha3f,aalpha4f
                                       ,FconsumptionNWNW_TOG,
                                       FeffortMean12NWNW_TOG,0,CSkills_factor12NWNW_TOG)+
            FprefshockNW;
            
            FutilityWNW_TOG=F_utility(aalpha1f,aalpha2f,aalpha3f,aalpha4f,
                                      FconsumptionWNW_TOG,
                                      FeffortMean12WNW_TOG,1,CSkills_factor12WNW_TOG)+
            FprefshockW;
            
            FutilityNWW_TOG=F_utility(aalpha1f,aalpha2f,aalpha3f,aalpha4f,
                                      FconsumptionNWW_TOG,
                                      FeffortMean12NWW_TOG,0,CSkills_factor12NWW_TOG)+
            FprefshockNW;
            
            FutilityWW_TOG=F_utility(aalpha1f,aalpha2f,aalpha3f,aalpha4f,
                                     FconsumptionWW_TOG,
                                     FeffortMean12WW_TOG,1,CSkills_factor12WW_TOG)+
            FprefshockW;
            
            //Seeing the utility levels for each case:
            
            
            //==========================================================
            // Defining the utilities for the possible combinations
            //==========================================================
            
            WelfareNWNW=mupred*FutilityNWNW_TOG+(1-mupred)*MutilityNWNW_TOG;
            WelfareWNW=mupred*FutilityWNW_TOG+(1-mupred)*MutilityWNW_TOG;
            WelfareNWW=mupred*FutilityNWW_TOG+(1-mupred)*MutilityNWW_TOG;
            WelfareWW=mupred*FutilityWW_TOG+(1-mupred)*MutilityWW_TOG;
            
            //==========================================================
            //Identifying the correct decision
            //==========================================================
            if ( (WelfareNWNW>WelfareWNW) && ( WelfareNWNW>WelfareNWW) && ( WelfareNWNW>WelfareWW)) {
                Decision=4;
            }
            else if ( (WelfareWNW>WelfareNWW) && ( WelfareNWNW>WelfareWW) ) {
                Decision=2;
            }
            else if ( (WelfareNWW>WelfareWW)){
                Decision=3;
            }
            else{
                Decision=1;
            }
            //Storing the value of the maximum decision in 2010
            VMAX12=(Decision==1)*WelfareWW+
            (Decision==2)*WelfareWNW+
            (Decision==3)*WelfareNWW+
            (Decision==4)*WelfareNWNW;
            
            
            //Storing the value of the obseved decision in 2010
            VOBS12=
            (obDecision==1)*WelfareWW+
            (obDecision==2)*WelfareWNW+
            (obDecision==3)*WelfareNWW+
            (obDecision==4)*WelfareNWNW;
            
            VSUM12=(exp(WelfareWW-VMAX12)+
                    exp(WelfareWNW-VMAX12)+
                    exp(WelfareNWW-VMAX12)+
                    exp(WelfareNWNW-VMAX12))/TAUSMOOTHING;
            
            
            
            VDIF12=exp(VOBS12-VMAX12)/TAUSMOOTHING;
            
            
            if(VDIF12>9.0e10){
                VDIF12=1.0e5;
            }
            VCONTRIB12=VDIF12/VSUM12;
            
            Decision=(Decision==obDecision);
            //7. Likelihood of measurement system
            //===============================
            
            //7.1 NWNW
            if (Mfraclabor12==0 & Ffraclabor12==0){
                
                //1. Skills
                
                loglikskills_rr= MEASloglikSkills2012(SK2012Measure1, M_A1_S2012,M_VS1_S2012,SK2012Measure2, M_A2_S2012,M_VS2_S2012,SK2012Measure3, M_A3_S2012,M_VS3_S2012,SK2012Measure4, M_A4_S2012,M_VS4_S2012,SK2012Measure5, M_A5_S2012,M_VS5_S2012,SK2012Measure6, M_A6_S2012,M_VS6_S2012,SK2012Measure7, M_A7_S2012,M_VS7_S2012,SK2012Measure8, M_A8_S2012,M_VS8_S2012,SK2012Measure9, M_A9_S2012,M_VS9_S2012,SK2012Measure10, M_A10_S2012,M_VS10_S2012,SK2012Measure11, M_A11_S2012,M_VS11_S2012,SK2012Measure12,M_A12_S2012,M_VS12_S2012,SK2012Measure13, M_A13_S2012, M_VS13_S2012,CSkills_factor12NWNW_TOG);
                
                
                //loglikskills_rr=F_loglikelihood_generic(Ctestsfactorsss2012,
                //                                        CSkills_factor12NWNW_TOG,MEASSkills);
                
                //Storing the particle weight and the particle itself
                
                what2_rr_vector[rr]=loglikskills_rr;
                S2Vector[rr]=CSkills_factor12NWNW_TOG;
                
                
                //2. Mother's effort
                
                loglikeffm_rr=MEASloglikFEFFORT2012(MEFF1_12, M_A1_EFF12, M_V1_EFF12,
                                                    MEFF2_12, M_A2_EFF12, M_V2_EFF12,
                                                    MEFF3_12, M_A3_EFF12, M_V3_EFF12,
                                                    MEFF4_12, M_A4_EFF12, M_V4_EFF12,
                                                    MEFF5_12, M_A5_EFF12, M_V5_EFF12,
                                                    MEFF6_12, M_A6_EFF12, M_V6_EFF12,
                                                    MEFF7_12, M_A7_EFF12, M_V7_EFF12,
                                                    MEFF8_12, M_A8_EFF12, M_V8_EFF12,
                                                    MEFF9_12, M_A9_EFF12, M_V9_EFF12,
                                                    MEFF10_12, M_A10_EFF12, M_V10_EFF12,
                                                    MEFF11_12, M_A11_EFF12, M_V11_EFF12,
                                                    MEFF12_12, M_A12_EFF12, M_V12_EFF12,
                                                    MEFF13_12, M_A13_EFF12, M_V13_EFF12,
                                                    MEFF14_12, M_A14_EFF12, M_V14_EFF12, MeffortMean12NWNW_TOG);
                
                
                //loglikeffm_rr=F_loglikelihood_generic(Meffort12,MeffortMean12NWNW_TOG,MEASMeffort);
                
                //3. Father's effort
                
                loglikefff_rr=MEASloglikFEFFORT2012(FEFF1_12, M_A1_EFF12, M_V1_EFF12,
                                                    FEFF2_12, M_A2_EFF12, M_V2_EFF12,
                                                    FEFF3_12, M_A3_EFF12, M_V3_EFF12,
                                                    FEFF4_12, M_A4_EFF12, M_V4_EFF12,
                                                    FEFF5_12, M_A5_EFF12, M_V5_EFF12,
                                                    FEFF6_12, M_A6_EFF12, M_V6_EFF12,
                                                    FEFF7_12, M_A7_EFF12, M_V7_EFF12,
                                                    FEFF8_12, M_A8_EFF12, M_V8_EFF12,
                                                    FEFF9_12, M_A9_EFF12, M_V9_EFF12,
                                                    FEFF10_12, M_A10_EFF12, M_V10_EFF12,
                                                    FEFF11_12, M_A11_EFF12, M_V11_EFF12,
                                                    FEFF12_12, M_A12_EFF12, M_V12_EFF12,
                                                    FEFF13_12, M_A13_EFF12, M_V13_EFF12,
                                                    FEFF14_12, M_A14_EFF12, M_V14_EFF12, FeffortMean12NWNW_TOG);
                
                
                //loglikefff_rr=F_loglikelihood_generic(Feffort12,FeffortMean12NWNW_TOG,MEASFeffort);
                
                //4. Investment level
                loglikeinv_rr=MEASloginv12(MINV1_12,M_V1_INV12,M_A1_INV12,
                                           MINV2_12,M_V2_INV12,M_A2_INV12,
                                           MINV3_12,M_V3_INV12,M_A3_INV12,
                                           MINV4_12,M_V4_INV12,M_A4_INV12,
                                           MINV5_12,M_V5_INV12,M_A5_INV12,
                                           MINV6_12,M_V6_INV12,M_A6_INV12,
                                           MINV7_12,M_V7_INV12,M_A7_INV12,
                                           MINV8_12,M_V8_INV12,M_A8_INV12,
                                           MINV9_12,M_V9_INV12,M_A9_INV12,
                                           MINV10_12,M_V10_INV12,M_A10_INV12,
                                           MINV11_12,M_V11_INV12,M_A11_INV12,
                                           MINV12_12,M_V12_INV12,M_A12_INV12,
                                           MINV13_12,M_V13_INV12,M_A13_INV12,
                                           MINV14_12,M_V14_INV12,M_A14_INV12,
                                           MINV15_12,M_V15_INV12,M_A15_INV12,
                                           MINV16_12,M_V16_INV12,M_A16_INV12,
                                           MINV17_12,M_V17_INV12,M_A17_INV12,
                                           MINV18_12,M_V18_INV12,M_A18_INV12,
                                           MINV19_12,M_V19_INV12,M_A19_INV12,
                                           MINV20_12,M_V20_INV12,M_A20_INV12,
                                           MINV21_12,M_V21_INV12,M_A21_INV12,
                                           CInvestmentMean12NWNW_TOG);
                //Storing the final Decisions
                Feffort12Final=FeffortMean12NWNW_TOG;
                Meffort12Final=MeffortMean12NWNW_TOG;
                Investment12Final=CInvestmentMean12NWNW_TOG;
            }
            
            
            //7.2 WNW
            if (Mfraclabor12==0 & Ffraclabor12==1){
                
                //1. Skills
                loglikskills_rr= MEASloglikSkills2012(SK2012Measure1, M_A1_S2012,M_VS1_S2012,SK2012Measure2, M_A2_S2012,M_VS2_S2012,SK2012Measure3, M_A3_S2012,M_VS3_S2012,SK2012Measure4, M_A4_S2012,M_VS4_S2012,SK2012Measure5, M_A5_S2012,M_VS5_S2012,SK2012Measure6, M_A6_S2012,M_VS6_S2012,SK2012Measure7, M_A7_S2012,M_VS7_S2012,SK2012Measure8, M_A8_S2012,M_VS8_S2012,SK2012Measure9, M_A9_S2012,M_VS9_S2012,SK2012Measure10, M_A10_S2012,M_VS10_S2012,SK2012Measure11, M_A11_S2012,M_VS11_S2012,SK2012Measure12,M_A12_S2012,M_VS12_S2012,SK2012Measure13, M_A13_S2012, M_VS13_S2012,CSkills_factor12WNW_TOG);
                
                //loglikskills_rr=F_loglikelihood_generic(Ctestsfactorsss2012,
                //                                        CSkills_factor12WNW_TOG,MEASSkills);
                
                //Storing the particle weight and the particle itself
                
                what2_rr_vector[rr]=loglikskills_rr;
                S2Vector[rr]=CSkills_factor12WNW_TOG;
                
                
                //2. Mother's effort
                loglikeffm_rr=MEASloglikFEFFORT2012(MEFF1_12, M_A1_EFF12, M_V1_EFF12,
                                                    MEFF2_12, M_A2_EFF12, M_V2_EFF12,
                                                    MEFF3_12, M_A3_EFF12, M_V3_EFF12,
                                                    MEFF4_12, M_A4_EFF12, M_V4_EFF12,
                                                    MEFF5_12, M_A5_EFF12, M_V5_EFF12,
                                                    MEFF6_12, M_A6_EFF12, M_V6_EFF12,
                                                    MEFF7_12, M_A7_EFF12, M_V7_EFF12,
                                                    MEFF8_12, M_A8_EFF12, M_V8_EFF12,
                                                    MEFF9_12, M_A9_EFF12, M_V9_EFF12,
                                                    MEFF10_12, M_A10_EFF12, M_V10_EFF12,
                                                    MEFF11_12, M_A11_EFF12, M_V11_EFF12,
                                                    MEFF12_12, M_A12_EFF12, M_V12_EFF12,
                                                    MEFF13_12, M_A13_EFF12, M_V13_EFF12,
                                                    MEFF14_12, M_A14_EFF12, M_V14_EFF12, MeffortMean12WNW_TOG);
                
                //loglikeffm_rr=F_loglikelihood_generic(Meffort12, MeffortMean12WNW_TOG,MEASMeffort);
                
                //3. Father's effort
                
                loglikefff_rr=MEASloglikFEFFORT2012(FEFF1_12, M_A1_EFF12, M_V1_EFF12,
                                                    FEFF2_12, M_A2_EFF12, M_V2_EFF12,
                                                    FEFF3_12, M_A3_EFF12, M_V3_EFF12,
                                                    FEFF4_12, M_A4_EFF12, M_V4_EFF12,
                                                    FEFF5_12, M_A5_EFF12, M_V5_EFF12,
                                                    FEFF6_12, M_A6_EFF12, M_V6_EFF12,
                                                    FEFF7_12, M_A7_EFF12, M_V7_EFF12,
                                                    FEFF8_12, M_A8_EFF12, M_V8_EFF12,
                                                    FEFF9_12, M_A9_EFF12, M_V9_EFF12,
                                                    FEFF10_12, M_A10_EFF12, M_V10_EFF12,
                                                    FEFF11_12, M_A11_EFF12, M_V11_EFF12,
                                                    FEFF12_12, M_A12_EFF12, M_V12_EFF12,
                                                    FEFF13_12, M_A13_EFF12, M_V13_EFF12,
                                                    FEFF14_12, M_A14_EFF12, M_V14_EFF12, FeffortMean12WNW_TOG);
                
                //loglikefff_rr=F_loglikelihood_generic(Feffort12,FeffortMean12WNW_TOG,MEASFeffort);
                
                //4. Investment level
                
                loglikeinv_rr=MEASloginv12(MINV1_12,M_V1_INV12,M_A1_INV12,
                                           MINV2_12,M_V2_INV12,M_A2_INV12,
                                           MINV3_12,M_V3_INV12,M_A3_INV12,
                                           MINV4_12,M_V4_INV12,M_A4_INV12,
                                           MINV5_12,M_V5_INV12,M_A5_INV12,
                                           MINV6_12,M_V6_INV12,M_A6_INV12,
                                           MINV7_12,M_V7_INV12,M_A7_INV12,
                                           MINV8_12,M_V8_INV12,M_A8_INV12,
                                           MINV9_12,M_V9_INV12,M_A9_INV12,
                                           MINV10_12,M_V10_INV12,M_A10_INV12,
                                           MINV11_12,M_V11_INV12,M_A11_INV12,
                                           MINV12_12,M_V12_INV12,M_A12_INV12,
                                           MINV13_12,M_V13_INV12,M_A13_INV12,
                                           MINV14_12,M_V14_INV12,M_A14_INV12,
                                           MINV15_12,M_V15_INV12,M_A15_INV12,
                                           MINV16_12,M_V16_INV12,M_A16_INV12,
                                           MINV17_12,M_V17_INV12,M_A17_INV12,
                                           MINV18_12,M_V18_INV12,M_A18_INV12,
                                           MINV19_12,M_V19_INV12,M_A19_INV12,
                                           MINV20_12,M_V20_INV12,M_A20_INV12,
                                           MINV21_12,M_V21_INV12,M_A21_INV12,
                                           CInvestmentMean12WNW_TOG);
                
                //Storing the final Decisions
                Feffort12Final=FeffortMean12WNW_TOG;
                Meffort12Final=MeffortMean12WNW_TOG;
                Investment12Final=CInvestmentMean12WNW_TOG;
            }
            //7.3 NWW
            if (Mfraclabor12==1 & Ffraclabor12==0){
                
                //1. Skills
                loglikskills_rr= MEASloglikSkills2012(SK2012Measure1, M_A1_S2012,M_VS1_S2012,SK2012Measure2, M_A2_S2012,M_VS2_S2012,SK2012Measure3, M_A3_S2012,M_VS3_S2012,SK2012Measure4, M_A4_S2012,M_VS4_S2012,SK2012Measure5, M_A5_S2012,M_VS5_S2012,SK2012Measure6, M_A6_S2012,M_VS6_S2012,SK2012Measure7, M_A7_S2012,M_VS7_S2012,SK2012Measure8, M_A8_S2012,M_VS8_S2012,SK2012Measure9, M_A9_S2012,M_VS9_S2012,SK2012Measure10, M_A10_S2012,M_VS10_S2012,SK2012Measure11, M_A11_S2012,M_VS11_S2012,SK2012Measure12,M_A12_S2012,M_VS12_S2012,SK2012Measure13, M_A13_S2012, M_VS13_S2012,CSkills_factor12NWW_TOG);
                
                //loglikskills_rr=F_loglikelihood_generic(Ctestsfactorsss2012,
                //                                        CSkills_factor12NWW_TOG,MEASSkills);
                
                //Storing the particle weight and the particle itself
                
                what2_rr_vector[rr]=loglikskills_rr;
                S2Vector[rr]=CSkills_factor12NWW_TOG;
                
                
                //2. Mother's effort
                
                loglikeffm_rr=MEASloglikFEFFORT2012(MEFF1_12, M_A1_EFF12, M_V1_EFF12,
                                                    MEFF2_12, M_A2_EFF12, M_V2_EFF12,
                                                    MEFF3_12, M_A3_EFF12, M_V3_EFF12,
                                                    MEFF4_12, M_A4_EFF12, M_V4_EFF12,
                                                    MEFF5_12, M_A5_EFF12, M_V5_EFF12,
                                                    MEFF6_12, M_A6_EFF12, M_V6_EFF12,
                                                    MEFF7_12, M_A7_EFF12, M_V7_EFF12,
                                                    MEFF8_12, M_A8_EFF12, M_V8_EFF12,
                                                    MEFF9_12, M_A9_EFF12, M_V9_EFF12,
                                                    MEFF10_12, M_A10_EFF12, M_V10_EFF12,
                                                    MEFF11_12, M_A11_EFF12, M_V11_EFF12,
                                                    MEFF12_12, M_A12_EFF12, M_V12_EFF12,
                                                    MEFF13_12, M_A13_EFF12, M_V13_EFF12,
                                                    MEFF14_12, M_A14_EFF12, M_V14_EFF12, MeffortMean12NWW_TOG);
                
                
                //loglikeffm_rr=F_loglikelihood_generic(Meffort12,MeffortMean12NWW_TOG,MEASMeffort);
                
                //3. Father's effort
                
                loglikefff_rr=MEASloglikFEFFORT2012(FEFF1_12, M_A1_EFF12, M_V1_EFF12,
                                                    FEFF2_12, M_A2_EFF12, M_V2_EFF12,
                                                    FEFF3_12, M_A3_EFF12, M_V3_EFF12,
                                                    FEFF4_12, M_A4_EFF12, M_V4_EFF12,
                                                    FEFF5_12, M_A5_EFF12, M_V5_EFF12,
                                                    FEFF6_12, M_A6_EFF12, M_V6_EFF12,
                                                    FEFF7_12, M_A7_EFF12, M_V7_EFF12,
                                                    FEFF8_12, M_A8_EFF12, M_V8_EFF12,
                                                    FEFF9_12, M_A9_EFF12, M_V9_EFF12,
                                                    FEFF10_12, M_A10_EFF12, M_V10_EFF12,
                                                    FEFF11_12, M_A11_EFF12, M_V11_EFF12,
                                                    FEFF12_12, M_A12_EFF12, M_V12_EFF12,
                                                    FEFF13_12, M_A13_EFF12, M_V13_EFF12,
                                                    FEFF14_12, M_A14_EFF12, M_V14_EFF12, FeffortMean12NWW_TOG);
                
                //loglikefff_rr=F_loglikelihood_generic(Feffort12,FeffortMean12NWW_TOG,MEASFeffort);
                
                //4. Investment level
                
                loglikeinv_rr=MEASloginv12(MINV1_12,M_V1_INV12,M_A1_INV12,
                                           MINV2_12,M_V2_INV12,M_A2_INV12,
                                           MINV3_12,M_V3_INV12,M_A3_INV12,
                                           MINV4_12,M_V4_INV12,M_A4_INV12,
                                           MINV5_12,M_V5_INV12,M_A5_INV12,
                                           MINV6_12,M_V6_INV12,M_A6_INV12,
                                           MINV7_12,M_V7_INV12,M_A7_INV12,
                                           MINV8_12,M_V8_INV12,M_A8_INV12,
                                           MINV9_12,M_V9_INV12,M_A9_INV12,
                                           MINV10_12,M_V10_INV12,M_A10_INV12,
                                           MINV11_12,M_V11_INV12,M_A11_INV12,
                                           MINV12_12,M_V12_INV12,M_A12_INV12,
                                           MINV13_12,M_V13_INV12,M_A13_INV12,
                                           MINV14_12,M_V14_INV12,M_A14_INV12,
                                           MINV15_12,M_V15_INV12,M_A15_INV12,
                                           MINV16_12,M_V16_INV12,M_A16_INV12,
                                           MINV17_12,M_V17_INV12,M_A17_INV12,
                                           MINV18_12,M_V18_INV12,M_A18_INV12,
                                           MINV19_12,M_V19_INV12,M_A19_INV12,
                                           MINV20_12,M_V20_INV12,M_A20_INV12,
                                           MINV21_12,M_V21_INV12,M_A21_INV12,
                                           CInvestmentMean12NWW_TOG);
                
                //Storing the final Decisions
                Feffort12Final=FeffortMean12NWW_TOG;
                Meffort12Final=MeffortMean12NWW_TOG;
                Investment12Final=CInvestmentMean12NWW_TOG;
            }
            //7.4 WW
            if (Mfraclabor12==1 & Ffraclabor12==1){
                
                //1. Skills
                loglikskills_rr= MEASloglikSkills2012(SK2012Measure1, M_A1_S2012,M_VS1_S2012,SK2012Measure2, M_A2_S2012,M_VS2_S2012,SK2012Measure3, M_A3_S2012,M_VS3_S2012,SK2012Measure4, M_A4_S2012,M_VS4_S2012,SK2012Measure5, M_A5_S2012,M_VS5_S2012,SK2012Measure6, M_A6_S2012,M_VS6_S2012,SK2012Measure7, M_A7_S2012,M_VS7_S2012,SK2012Measure8, M_A8_S2012,M_VS8_S2012,SK2012Measure9, M_A9_S2012,M_VS9_S2012,SK2012Measure10, M_A10_S2012,M_VS10_S2012,SK2012Measure11, M_A11_S2012,M_VS11_S2012,SK2012Measure12,M_A12_S2012,M_VS12_S2012,SK2012Measure13, M_A13_S2012, M_VS13_S2012,CSkills_factor12WW_TOG);
                
                //loglikskills_rr=F_loglikelihood_generic(Ctestsfactorsss2012,
                //                                        CSkills_factor12WW_TOG,MEASSkills);
                
                //Storing the particle weight and the particle itself
                
                what2_rr_vector[rr]=loglikskills_rr;
                S2Vector[rr]=CSkills_factor12WW_TOG;
                
                
                //2. Mother's effort
                loglikeffm_rr=MEASloglikFEFFORT2012(MEFF1_12, M_A1_EFF12, M_V1_EFF12,
                                                    MEFF2_12, M_A2_EFF12, M_V2_EFF12,
                                                    MEFF3_12, M_A3_EFF12, M_V3_EFF12,
                                                    MEFF4_12, M_A4_EFF12, M_V4_EFF12,
                                                    MEFF5_12, M_A5_EFF12, M_V5_EFF12,
                                                    MEFF6_12, M_A6_EFF12, M_V6_EFF12,
                                                    MEFF7_12, M_A7_EFF12, M_V7_EFF12,
                                                    MEFF8_12, M_A8_EFF12, M_V8_EFF12,
                                                    MEFF9_12, M_A9_EFF12, M_V9_EFF12,
                                                    MEFF10_12, M_A10_EFF12, M_V10_EFF12,
                                                    MEFF11_12, M_A11_EFF12, M_V11_EFF12,
                                                    MEFF12_12, M_A12_EFF12, M_V12_EFF12,
                                                    MEFF13_12, M_A13_EFF12, M_V13_EFF12,
                                                    MEFF14_12, M_A14_EFF12, M_V14_EFF12, MeffortMean12WW_TOG);
                
                
                //loglikeffm_rr=F_loglikelihood_generic(Meffort12,MeffortMean12WW_TOG,MEASMeffort);
                
                //3. Father's effort
                
                loglikefff_rr=MEASloglikFEFFORT2012(FEFF1_12, M_A1_EFF12, M_V1_EFF12,
                                                    FEFF2_12, M_A2_EFF12, M_V2_EFF12,
                                                    FEFF3_12, M_A3_EFF12, M_V3_EFF12,
                                                    FEFF4_12, M_A4_EFF12, M_V4_EFF12,
                                                    FEFF5_12, M_A5_EFF12, M_V5_EFF12,
                                                    FEFF6_12, M_A6_EFF12, M_V6_EFF12,
                                                    FEFF7_12, M_A7_EFF12, M_V7_EFF12,
                                                    FEFF8_12, M_A8_EFF12, M_V8_EFF12,
                                                    FEFF9_12, M_A9_EFF12, M_V9_EFF12,
                                                    FEFF10_12, M_A10_EFF12, M_V10_EFF12,
                                                    FEFF11_12, M_A11_EFF12, M_V11_EFF12,
                                                    FEFF12_12, M_A12_EFF12, M_V12_EFF12,
                                                    FEFF13_12, M_A13_EFF12, M_V13_EFF12,
                                                    FEFF14_12, M_A14_EFF12, M_V14_EFF12, FeffortMean12WW_TOG);
                
                //loglikefff_rr=F_loglikelihood_generic(Feffort12,FeffortMean12WW_TOG,MEASFeffort);
                
                //4. Investment level
                loglikeinv_rr=MEASloginv12(MINV1_12,M_V1_INV12,M_A1_INV12,
                                           MINV2_12,M_V2_INV12,M_A2_INV12,
                                           MINV3_12,M_V3_INV12,M_A3_INV12,
                                           MINV4_12,M_V4_INV12,M_A4_INV12,
                                           MINV5_12,M_V5_INV12,M_A5_INV12,
                                           MINV6_12,M_V6_INV12,M_A6_INV12,
                                           MINV7_12,M_V7_INV12,M_A7_INV12,
                                           MINV8_12,M_V8_INV12,M_A8_INV12,
                                           MINV9_12,M_V9_INV12,M_A9_INV12,
                                           MINV10_12,M_V10_INV12,M_A10_INV12,
                                           MINV11_12,M_V11_INV12,M_A11_INV12,
                                           MINV12_12,M_V12_INV12,M_A12_INV12,
                                           MINV13_12,M_V13_INV12,M_A13_INV12,
                                           MINV14_12,M_V14_INV12,M_A14_INV12,
                                           MINV15_12,M_V15_INV12,M_A15_INV12,
                                           MINV16_12,M_V16_INV12,M_A16_INV12,
                                           MINV17_12,M_V17_INV12,M_A17_INV12,
                                           MINV18_12,M_V18_INV12,M_A18_INV12,
                                           MINV19_12,M_V19_INV12,M_A19_INV12,
                                           MINV20_12,M_V20_INV12,M_A20_INV12,
                                           MINV21_12,M_V21_INV12,M_A21_INV12,
                                           CInvestmentMean12WW_TOG);
                //Storing the final Decisions
                Feffort12Final=FeffortMean12WW_TOG;
                Meffort12Final=MeffortMean12WW_TOG;
                Investment12Final=CInvestmentMean12WW_TOG;
            }
            
            //5. Measurement system for mmu, doesn't depend on labor supply
            loglikmmu_rr=MEASloglikBARG(MBARG1, M_V1_BARG, M_A1_BARG,
                                        MBARG2, M_V2_BARG, M_A2_BARG,
                                        MBARG3, M_V3_BARG, M_A3_BARG,
                                        MBARG4, M_V4_BARG, M_A4_BARG,
                                        MBARG5, M_V5_BARG, M_A5_BARG,
                                        MBARG6, M_V6_BARG, M_A6_BARG,
                                        MBARG7, M_V7_BARG, M_A7_BARG,
                                        MBARG8, M_V8_BARG, M_A8_BARG,
                                        MBARG9, M_V9_BARG, M_A9_BARG,
                                        MBARG10, M_V10_BARG, M_A10_BARG,
                                        MBARG11, M_V11_BARG, M_A11_BARG,
                                        MBARG12, M_V12_BARG, M_A12_BARG,
                                        MBARG13, M_V13_BARG, M_A13_BARG,
                                        MBARG14, M_V14_BARG, M_A14_BARG,
                                        MBARG15, M_V15_BARG, M_A15_BARG,
                                        MBARG16, M_V16_BARG, M_A16_BARG,
                                        MBARG17, M_V17_BARG, M_A17_BARG,
                                        MBARG18, M_V18_BARG, M_A18_BARG,
                                        MBARG19, M_V19_BARG, M_A19_BARG,
                                        mupred);
            // loglikmmu_rr=F_loglikelihood_generic(Hbarg,mupred,MEASMMu);
            
        }//If living together
        
        
        
        
        //Adding up likelihood contributions of simulations
        
        
        //0 Decisions (behavioral system)
        logcontribdecision+=Decision; //Storing the number of times correct thing
        logcontribobdecision12SMOOTHED+=VCONTRIB12;
        
        
        
        //2. Mother's effort
        loglikeffm+=loglikeffm_rr;
        
        //3. Father's effort
        loglikefff+=loglikefff_rr;
        
        //4. Investment
        loglikeinv+=loglikeinv_rr;
        
        //5. Mmu
        loglikmmu+=loglikmmu_rr;
        
        rr=rr+1;
    } //Finish the RR simulations
    //At the end of the R computations we will perform the fraction of times person
    //decided to work.
    //Block of particle filter weights and re-sample
    //==========new block of taking out particles randomly (START)==========//
    //Storing the weights that are NOT normalized and that will give the likelihood
    wtilde2_rr_vector=what2_rr_vector;
    
    //Normalizing weights of second period
    
    AuxiliarV=0;
    for(int mm=0; mm<RR;mm=mm+1){
        AuxiliarV+=what2_rr_vector[mm];
    }
    for(int mm=0; mm<RR; mm=mm+1){
        what2_rr_vector[mm]=what2_rr_vector[mm]/AuxiliarV;
    }
    AuxiliarV=0;
    
    
    boost::random::discrete_distribution<int> distskills2(what2_rr_vector);
    iterator=0;
    wsum2=0;
    for (int mm=0; mm<RR;mm=mm+1){
        iterator=distskills2(rng);
        S2Vectoraux[mm]=S2Vector[iterator];
        w2_rr_vectorAUX[mm]=what2_rr_vector[iterator];
        wsum2+=wtilde2_rr_vector[iterator];//Liklihood is adding up the weight (original weight) of those chosen (iterator)
        
        //Re-adapting the skills in t=0 and t=1. This is done for the smoothing distribution. Also, re-sampling
        //done in every period guarantees contribution of the likelihood is given as the sum of weights.
        
        
        S0Vectoraux[mm]=S0Vector[iterator];
        w0_rr_vectorAUX[mm]=w0_rr_vector[iterator];
        
        S1Vectoraux[mm]=S1Vector[iterator];
        w1_rr_vectorAUX[mm]=w1_rr_vector[iterator];
        
        
    }
    //Updating the weights and the skills
    
    //Redefining vector of skills
    
    S1Vector=S1Vectoraux;
    what1_rr_vector=w1_rr_vectorAUX;
    w1_rr_vector=what1_rr_vector;
    
    S0Vector=S0Vectoraux;
    w0_rr_vector=w0_rr_vectorAUX;
    
    S2Vector=S2Vectoraux;
    what2_rr_vector=w2_rr_vectorAUX;
    w2_rr_vector=what2_rr_vector;
    
    //Normalizing the weights of skills for probabilities
    AuxiliarV=0;
    AuxiliarV1=0;
    AuxiliarV2=0;
    for (int mm=0; mm<RR;mm=mm+1){
        AuxiliarV+=w0_rr_vector[mm];
        AuxiliarV1+=w1_rr_vector[mm];
        AuxiliarV2+=w2_rr_vector[mm];
    }
    
    for (int mm=0; mm<RR;mm=mm+1){
        w0_rr_vector[mm]=w0_rr_vector[mm]/AuxiliarV;
        w1_rr_vector[mm]=w1_rr_vector[mm]/AuxiliarV1;
        w2_rr_vector[mm]=w2_rr_vector[mm]/AuxiliarV2;
    }
    AuxiliarV=0;
    AuxiliarV1=0;
    AuxiliarV2=0;
    
    
    //==============================================//
    //Smoothing distribution of the particle filter.//
    //==============================================//
    //Need to define predicted investments, effort levels and childcare decissions.
    //So that I don't have to do it by cases
    if (SMOOTHINGFILTER==1){
        //==========BOOTSTRAP SMOOTHER FILTER ===================================//
        
        //First I will re-normalize the weights used as input to add up to one.
        for (int mm=0;mm<RR;mm=mm+1){
            
            weightsmooth0+=w0_rr_vector[mm];
            weightsmooth1+=w1_rr_vector[mm];
            weightsmooth2+=w2_rr_vector[mm];
        }
        
        //Normalizing it
        for (int mm=0; mm<RR; mm=mm+1){
            w0_rr_vector[mm]=w0_rr_vector[mm]/weightsmooth0;
            
            w1_rr_vector[mm]=w1_rr_vector[mm]/weightsmooth1;
            w2_rr_vector[mm]=w2_rr_vector[mm]/weightsmooth2;
            
        }
        
        
        //Start the backwards recursion to find the corresponding weights.
        VECWEIGHTSSMOOTH2=w2_rr_vector;
        
        //Generating matrices in t==2
        for (int indexOUT=0; indexOUT<RR; indexOUT=indexOUT+1){
            S2MATRIX[indexOUT]=S2Vector[indexOUT];
            W2MATRIX[indexOUT]=VECWEIGHTSSMOOTH2[indexOUT];
            if(1==2){
                cout <<  S2MATRIX[indexOUT] << "  S2MATRIX[ii][indexOUT] " << endl;
                cout << indexOUT << " indexout 1 " << endl;
            }
        }
        denomL1=0;
        denomL0=0;
        numL0=0;
        numL1=0;
        //#pragma omp parallel for private(MEASSkills,denomL1,numL1) shared(ii, RR)
        for(int indexOUT=0; indexOUT<RR; indexOUT=indexOUT+1){
            for( int Middle=0; Middle<RR; Middle=Middle+1){
                for(int iNner=0; iNner<RR; iNner=iNner+1){
                    MEASSkills=F_predskills(ddelta0,ddelta1,ddelta2,ddelta3_12,ddelta4,
                                            Cedad_meses12,
                                            ttheta0,ttheta1,ttheta2,pphi,
                                            ggammaf,ggammam,FeffortMean12WW_TOG,
                                            Meffort12Final,
                                            Investment12Final,
                                            S1Vector[iNner],Cchildcare12,
                                            PGVector[rr],Hmemberstotal12);
                    MEASSkills=F_loglikelihood_generic(S2Vector[Middle],MEASSkills,stdskills);
                    denomL1+=w1_rr_vector[iNner]*MEASSkills;
                }
                
                MEASSkills=F_predskills(ddelta0,ddelta1,ddelta2,ddelta3_12,ddelta4,
                                        Cedad_meses12,
                                        ttheta0,ttheta1,ttheta2,pphi,
                                        ggammaf,ggammam,Feffort12Final,
                                        Meffort12Final,
                                        Investment12Final,
                                        S1Vector[indexOUT],Cchildcare12,
                                        PGVector[rr],Hmemberstotal12);
                
                MEASSkills=F_loglikelihood_generic(S2Vector[Middle],MEASSkills,stdskills);
                numL1+=MEASSkills*VECWEIGHTSSMOOTH2[Middle]/denomL1;
            }
            VECWEIGHTSSMOOTH1[indexOUT]=w1_rr_vector[indexOUT]*numL1;
            S1MATRIX[indexOUT]=S1Vector[indexOUT];
            W1MATRIX[indexOUT]=VECWEIGHTSSMOOTH1[indexOUT];
            if(1==2){
                cout << indexOUT << " indexOUT2  " << endl;
            }
        }
        
        denomL1=0;
        denomL0=0;
        numL0=0;
        numL1=0;
        //#pragma omp parallel for private(MEASSkills,denomL0,numL0) shared(ii)
        for (int indexOUT=0; indexOUT<RR; indexOUT=indexOUT+1){ //Looping across ii particles
            for (int Middle=0; Middle<RR; Middle=Middle+1){
                for (int iNner=0; iNner<RR; iNner=iNner+1){
                    //The probability p(x_{t+1}|x_{t})  is given by the
                    //difference between the predicted level of skills in kk and the
                    //one found in vector jj. MEASskills will be used as auxuliar in this part
                    //1. Evaluating p(x_{t+1}|x_t)
                    //1.a Getting the predicted x_{t+1}|x_t
                    //cout << indexOUT << " INDEXOUT3  " << endl;
                    MEASSkills=F_predskills(ddelta0,ddelta1,ddelta2,ddelta3,ddelta4,
                                            Cedad_meses10,
                                            ttheta0,ttheta1,ttheta2,pphi,
                                            ggammaf,ggammam,Feffort10Final,
                                            Meffort10Final,
                                            Investment10Final,
                                            S0Vector[iNner],0,
                                            PGVector[rr],Hmemberstotal10);
                    
                    //b. Getting the p(x_{t+1}|x_t)
                    MEASSkills=F_loglikelihood_generic(S1Vector[Middle],MEASSkills,stdskills);
                    
                    //c. Getting the denominator and adding it up:
                    denomL0+=w0_rr_vector[iNner]*MEASSkills;
                    
                }
                //Computing the numerator
                //1. Evaluating p(x_{t+1}|x_t)
                //1.a Getting the predicted x_{t+1}|x_t
                
                MEASSkills=F_predskills(ddelta0,ddelta1,ddelta2,ddelta3,ddelta4,
                                        Cedad_meses10,
                                        ttheta0,ttheta1,ttheta2,pphi,
                                        ggammaf,ggammam,FeffortMean12WW_TOG,
                                        Meffort10Final,
                                        Investment10Final,
                                        S0Vector[indexOUT],0,
                                        PGVector[rr],Hmemberstotal10);
                
                //b. Getting the p(x_{t+1}|x_t)
                MEASSkills=F_loglikelihood_generic(S1Vector[Middle],MEASSkills,stdskills);
                numL0+=MEASSkills*VECWEIGHTSSMOOTH1[Middle]/denomL0;
                
                //c. Getting the denominator and adding it up:
                //denomL1+=w1_rr_vector[kk]*MEASSkills;
                
                
            }
            VECWEIGHTSSMOOTH0[indexOUT]=w0_rr_vector[indexOUT]*numL0;
            //cout << VECWEIGHTSSMOOTH0[iiii] << " VECWEIGHTSSMOOTH0[iiii]  " << endl;
            //I will now store the smoothed weights and the particles in the corresponding matrices
            
            S0MATRIX[indexOUT]=S0Vector[indexOUT];
            W0MATRIX[indexOUT]=VECWEIGHTSSMOOTH0[indexOUT];
        }//End of looping across iiii particles
        //ofstream optparamS0("S0MATRIX.csv");
        //ofstream optparamS1("S1MATRIX.csv");
        //ofstream optparamS2("S2MATRIX.csv");
        //ofstream optparamW0("W0MATRIX.csv");
        //ofstream optparamW1("W1MATRIX.csv");
        //ofstream optparamW2("W2MATRIX.csv");
        //for (int ii=0;ii<SIZE;ii=ii+1){
        //  for (int iiii=0; iiii<RR; iiii=iiii+1){
        //    optparamS0<<S0MATRIX[ii][iiii] << " , " ;
        //  optparamS1<<S1MATRIX[ii][iiii] << " , " ;
        //optparamS2<<S2MATRIX[ii][iiii] << " , " ;
        // optparamW0<<W0MATRIX[ii][iiii] << " , " ;
        // optparamW1<<W1MATRIX[ii][iiii] << " , " ;
        // optparamW2<<W2MATRIX[ii][iiii] << " , " ;
        // }
        // optparamS0 << " " << endl;
        // optparamS1 << " " << endl;
        // optparamS2 << " " << endl;
        // optparamW0 << " " << endl;
        // optparamW1 << " " << endl;
        // optparamW2 << " " << endl;
        // }
        //  optparamS0.close();
        //  optparamS1.close();
        //  optparamS2.close();
        //  optparamW0.close();
        //  optparamW1.close();
        //  optparamW2.close();
        
        
        //Setting up the output properly:
        
        //Setting up the ouptut properly:
        
        for(int rr=0; rr<RR; rr=rr+1){
            
            OUTPUT[rr][1]=S0MATRIX[rr];
            OUTPUT[rr][2]=S1MATRIX[rr];
            OUTPUT[rr][3]=S2MATRIX[rr];
            OUTPUT[rr][4]=W0MATRIX[rr];
            OUTPUT[rr][5]=W1MATRIX[rr];
            OUTPUT[rr][6]=W2MATRIX[rr];
            if(1==2){
                cout << S2MATRIX[rr] << " S2MATRIX[rr]"<< endl;
                cout << OUTPUT[rr][3]<< " OUTPUT[rr][3]"<< endl;
            }
            
        }
        //==========BOOTSTRAP SMOOTHER FILTER end ===================================//
    }
    
    
    
    //==============================================================//
    // Likelihood Function, if not running the smoothing distribution
    //==============================================================//
    
    
    //0. Behavioral decisions
    
    
    //2010. Without the smoothing logit estimator: logcontribdecision10=logcontribdecision10/RR;
    //logcontribdecision10=logcontribdecision10/RR;
    //Smoothing via logit simulator
    logcontribdecision10=logcontribobdecision10SMOOTHED/RR;
    
    
    //2012. //Without the logit estimator:
    //logcontribdecision=logcontribdecision/RR;
    //Smoothing via the logit simulator
    logcontribdecision=logcontribobdecision12SMOOTHED/RR;
    
    //cout << "----------"<< endl;
    //cout << logcontribdecision10 << " contrib12 not smoothed" << endl;
    //cout << logcontribdecision << " Contrib12 smoothed " << endl;
    //cout << "----------"<< endl;
    //We need to fix in those cases where we didn't catch a single time. The likelihood
    //will then be zero.
    
    
    //==================
    //Decisions
    //==================
    //2010
    if (logcontribdecision10==0){
        logcontribdecision10=1.0e-20;
    }
    logcontribdecision10=log(logcontribdecision10);
    //2012
    if (logcontribdecision==0){
        logcontribdecision=1.0e-20;
    }
    logcontribdecision=log(logcontribdecision);
    
    //=============
    //Skills
    //=============
    //1. Skills 2010-Already taking logs with
    //the function MEASloglikSkills2010
    
    loglikskills10=loglikskills10/RR;
    
    if (loglikskills10==0){
        loglikskills10=1.0e-20;
    }
    loglikskills10=log(loglikskills10);
    
    //1. Skills 2012
    loglikskills=loglikskills/RR;
    if (loglikskills==0){
        loglikskills=1.0e-20;
    }
    loglikskills=log(loglikskills);
    
    //==================
    //2. Mother's effort
    //==================
    //2010
    if (Cliveswithmother10==1){
        loglikeffm10=loglikeffm10/RR;
        if (loglikeffm10==0){
            loglikeffm10=1.0e-20;
        }
        loglikeffm10=log(loglikeffm10);
    }
    
    //2012
    if (Cliveswithmother12==1){
        loglikeffm=loglikeffm/RR;
        if (loglikeffm==0){
            loglikeffm=1.0e-20;
        }
        loglikeffm=log(loglikeffm);
    }
    
    
    
    //3. Father's effort
    
    //2010
    if (Cliveswithfather10==1){
        loglikefff10=loglikefff10/RR;
        if (loglikefff10==0){
            loglikefff10=1.0e-20;
        }
        loglikefff10=log(loglikefff10);
    }
    
    //2012
    if (Cliveswithfather12==1){
        loglikefff=loglikefff/RR;
        if (loglikefff==0){
            loglikefff=1.0e-20;
        }
        loglikefff=log(loglikefff);
    }
    //==================
    //4. Investment
    //=================
    
    //2010
    loglikeinv10=loglikeinv10/RR;
    if (loglikeinv10==0){
        loglikeinv10=1.0e-20;
    }
    loglikeinv10=log(loglikeinv10);
    
    //2012
    loglikeinv=loglikeinv/RR;
    if (loglikeinv==0){
        loglikeinv=1.0e-20;
    }
    loglikeinv=log(loglikeinv);
    
    //5. Mmu
    //2012
    if (Cliveswithfather12==1){
        loglikmmu=loglikmmu/RR;
        if (loglikmmu==0){
            loglikmmu=1.0e-20;
        }
        loglikmmu=log(loglikmmu);
    }
    
    //Skills at birth
    loglik_S0=loglik_S0rr;
    if (loglik_S0==0){
        loglik_S0=1.0e-20;
    }
    loglik_S0=log(loglik_S0);
    
    //Primary caregiver skills
    loglik_PG=loglik_PGrr/RR;
    if (loglik_PG==0){
        loglik_PG=1.0e-20;
    }
    loglik_PG=log(loglik_PG);
    
    //Summing weights of particle filter.
    
    
    wsum0=wsum0/RR;
    
    if (wsum0==0){
        wsum0=1.0e-20;
    }
    
    wsum0=log(wsum0);
    
    wsum1=wsum1/RR;
    
    if (wsum1==0){
        wsum1=1.0e-20;
    }
    
    wsum1=log(wsum1);
    
    wsum2=wsum2/RR;
    
    if (wsum2==0){
        wsum2=1.0e-20;
    }
    wsum2=log(wsum2);
    
    
    logcontribobs_ii+=likwageF10+likwageM10+likwageF+likwageM+
    loglikefff+loglikefff10+
    loglikeffm10+loglikeffm+
    logcontribdecision+logcontribdecision10+
    loglikeinv+loglikeinv10+
    wsum2+wsum1+wsum0+
    loglikmmu+loglik_PG;
    
    if(1==2){
        cout << "--------- " << endl;
        cout << WelfareWWA10 << " WelfareWWA10" << endl;
        cout << WelfareWWNA10 << " WelfareWWNA10" << endl;
        cout << WelfareNWWA10 << " WelfareNWWA10" << endl;
        cout << WelfareNWWNA10 << " WelfareNWWNA10" << endl;
        cout << WelfareWNWA10 << " WelfareWNWA10" << endl;
        cout << WelfareWNWNA10 << " WelfareWNWNA10" << endl;
        cout << WelfareNWNWA10 << " WelfareNWNWA10" << endl;
        cout << WelfareNWNWNA10 << " WelfareNWNWNA10" << endl;
        cout << " -- decisions  " << endl;
        cout << obDecision10 << " observed decision 2010" << endl;
        cout << Decision10 << " predicted decision 2010" << endl;
        cout << exp(logcontribdecision10) << " proportion of times predicted " << endl;
        cout << " ------  " << endl;
        cout << likwageF10 << " likwageF10"<< endl;
        cout << likwageM10 << " likwageM10"<< endl;
        cout << likwageF << " likwageF"<< endl;
        cout << likwageM << " likwageM"<< endl;
        cout << loglikefff << " likwageeff"<< endl;
        cout << loglikeffm << " likwageefmm"<< endl;
        cout << loglikefff10 << " likwageeff10"<< endl;
        cout << loglikeffm10 << " likwageefmm10"<< endl;
        cout << logcontribdecision << " likwagecontribob"<< endl;
        cout << logcontribdecision10 << " likwagecontri10"<< endl;
        cout << loglikeinv << " likwageinv"<< endl;
        cout << loglikeinv10 << " likwageinv10"<< endl;
        cout << loglikmmu << " loglikmmu"<< endl;
        cout << loglik_PG << " loglik_PG"<< endl;
        cout << wsum0 << " wsum0 " << endl;
        cout << wsum1 << " wsum1 " << endl;
        cout << wsum2 << " wsum2 " << endl;
        cout << logcontribobs_ii  <<  " logcontribobs_ii  " << endl;
        cout << "--------- " << endl;
    }
    
    //To the
    loglik=logcontribobs_ii;
    loglik=-loglik;
    
    
    
    OUTPUT[0][0]=loglik;
    return(OUTPUT);
}//End of function



//Define the function likelihood as function of vectors
long double F_likelihood(vector<double> par_likelihood,
                         vector<double>  V_Cliveswithfather12,
                         vector<double> V_Cliveswithmother12,
                         vector<double> V_Hhchores12,
                         vector<double> V_Meffort12,
                         vector<double> V_Feffort12,
                         vector<double> V_Mage12,
                         vector<double> V_Fage12,
                         vector<double> V_Cedad_meses12,
                         vector<double> V_Ctestsfactorsss2012,
                         vector<double> V_Ctestsfactor2ss_10,
                         vector<double> V_Myrschool12,
                         vector<double> V_Fyrschool12,
                         vector<double> V_Ffraclabor12,
                         vector<double> V_Mfraclabor12,
                         vector<double> V_Mwage12,
                         vector<double> V_Fwage12,
                         vector<double> V_Mnly12,
                         vector<double> V_Fnly12,
                         vector<double> V_Cfactorinv12,
                         vector<double> V_Cchildcare12,
                         vector<double> V_Ccareskills12,
                         vector<double> V_Hbarg,
                         vector<double> V_Feffort10,
                         vector<double> V_Meffort10,
                         vector<double> V_Cfactorinv10,
                         vector<double> V_Cbirthfactor,
                         vector<double> V_Mwage10,
                         vector<double> V_Fwage10,
                         vector<double> V_Mnly10,
                         vector<double> V_Fnly10,
                         vector<double> V_Cedad_meses10,
                         vector<double> V_Cliveswithfather10,
                         vector<double> V_Cliveswithmother10,
                         vector<double> V_Mage10,
                         vector<double> V_Fage10,
                         vector<double> V_Mfraclabor10,
                         vector<double> V_Ffraclabor10,
                         vector<double> V_Cchildcare10,
                         vector<double> V_Hchildcareobs,
                         vector<double> V_PG,
                         vector<double> V_MTJH,
                         vector<double> V_MRATIO,
                         vector<double> V_Unemployment,
                         vector<double> V_Wageratio,
                         vector<double>  V_Distance,
                         vector<double> V_Magegroup10,
                         vector<double>  V_Magegroup12,
                         vector<double>  V_Hmemberstotal10,
                         vector<double>  V_Hmemberstotal12,
                         vector <double> V_CSTDtepsi_pb_coo10,
                         vector <double> V_CSTDtepsi_pb_len10,
                         vector <double> V_CSTDtepsi_pb_mot10,
                         vector <double> V_CSTDtvip_pb10,
                         vector <double> V_CSTDcbcl1_pb_110,
                         vector <double> V_CSTDcbcl1_pb_210,
                         vector <double> V_CSTDcbcl1_pb_310,
                         vector <double> V_CSTDcbcl1_pb_410,
                         vector <double> V_CSTDcbcl1_pb_510,
                         vector <double> V_CSTDcbcl1_pb_610,
                         vector <double> V_CSTDcbcl1_pb_710,
                         vector <double> V_CSTDtadi_pb_cog12,
                         vector <double> V_CSTDtadi_pb_mot12,
                         vector <double> V_CSTDtadi_pb_len12,
                         vector <double> V_CSTDtadi_pb_se12,
                         vector <double> V_CSTDbt_112,
                         vector <double> V_CSTDbt_212,
                         vector <double> V_CSTDbt_312,
                         vector <double> V_CSTDbt_412,
                         vector <double> V_CSTDbt_512,
                         vector <double> V_CSTDbt_t12,
                         vector <double> V_CSTDhtks_st12,
                         vector <double> V_CSTDbdst_st12,
                         vector <double> V_CSTDppvt_t12,
                         vector <int> V_Ccondpregg3_1,
                         vector <int> V_Ccondpregg3_2,
                         vector <int> V_Ccondpregg3_3,
                         vector <int> V_Ccondpregg3_4,
                         vector <int> V_Ccondpregg3_5,
                         vector <int> V_Ccondpregg3_6,
                         vector <int> V_Ccondpregg3_7,
                         vector <int> V_Ccondpregg3_8,
                         vector <int> V_Ccondpregg3_9,
                         vector <int> V_Ccondpregg4a_1,
                         vector <int> V_Ccondpregg4a_2,
                         vector <int> V_Ccondpregg4a_3,
                         vector <int> V_Ccondpregg4a_4,
                         vector <int> V_Ccondpregg4a_5,
                         vector <int> V_Ccondpregg4a_6,
                         vector <int> V_Ccondpregg4a_7,
                         vector <double> V_Ccondpregg7b,
                         vector <double> V_Ccondpregg8b,
                         vector <double> V_Ccondpregg9,
                         vector <double> V_Ccondpregg11b,
                         vector <int> V_Ccondpregpreterm,
                         vector <double> V_Ccondpregg24,
                         vector <double> V_Ccondpregg23,
                         vector <double> V_Cwais_pb_num,
                         vector <double> V_Cwais_pb_vo,
                         vector <double> V_Cbfi_pb_ama,
                         vector <double> V_Cbfi_pb_ape,
                         vector <double> V_Cbfi_pb_ext,
                         vector <double> V_Cbfi_pb_neu,
                         vector <double> V_Cbfi_pb_res,
                         vector <double> V_Cpsi_pb_total,
                         vector <double> V_Hbargg2a,
                         vector <double> V_Hbargg2b,
                         vector <double> V_Hbargg2c,
                         vector <double> V_Hbargg2d,
                         vector <double> V_Hbargg2e,
                         vector <double> V_Hbargg2f,
                         vector <double> V_Hbargg2g,
                         vector <double> V_Hbargg2h,
                         vector <double> V_Hbargg2i,
                         vector <double> V_Hbargg2j,
                         vector <int> V_Hwomdecide12,
                         vector <int> V_Hfathdecides12,
                         vector <int> V_Hbothdecide,
                         vector <double> V_Hcaresacwom,
                         vector <double> V_Hcaresacman,
                         vector <int> V_Hopofman121,
                         vector <int> V_Hopofman122,
                         vector <int> V_Hopofman123,
                         vector <int> V_Hopofman124,
                         vector<int> V_Cinvh3,
                         vector<int> V_Cinvh4,
                         vector<int> V_Cinvh5,
                         vector<int> V_Cinvh6,
                         vector<int> V_Cinvh7,
                         vector<int> V_Cinvh8,
                         vector<int> V_Cinvh9,
                         vector<int> V_Cinvh13,
                         vector<double> V_Cinvf11a,
                         vector<double> V_Cinvf11b,
                         vector<double> V_Cinvf11c,
                         vector<double> V_Cinvf11d,
                         vector<double> V_Cinvf11e,
                         vector<double> V_Cinvf11f,
                         vector<double> V_Cinvf11g,
                         vector<double> V_Cinvf11h,
                         vector<double> V_Cinvf11i,
                         vector<double> V_Cinvf11j,
                         vector<double> V_Cinvf11k,
                         vector<int> V_Cinvhm2_10,
                         vector<int> V_Cinvhm2_11,
                         vector<int> V_Cinvhm2_12,
                         vector<int> V_Cinvhm2_13,
                         vector<int> V_Cinvhm2_14,
                         vector<int> V_Cinvhm2_15,
                         vector<int> V_Cinvhm2_16,
                         vector<int> V_Cinvhm2_18,
                         vector<double> V_Csharesbedroomhowmany12,
                         vector<double> V_Csharesbedhowmany12,
                         vector<int> V_g42_a2,
                         vector<int> V_g42_b2,
                         vector<int> V_g42_c2,
                         vector<int> V_g42_d2,
                         vector<int> V_g42_e2,
                         vector<int> V_g42_f2,
                         vector<int> V_g42_a1,
                         vector<int> V_g42_b1,
                         vector<int> V_g42_c1,
                         vector<int> V_g42_d1,
                         vector<int> V_g42_e1,
                         vector<int> V_g42_f1,
                         vector<double> V_f21a_p_t,
                         vector<double> V_f21b_p_t,
                         vector<double> V_f21c_p_t,
                         vector<double> V_f21d_p_t,
                         vector<double> V_f21e_p_t,
                         vector<double> V_f21f_p_t,
                         vector<double> V_f21g_p_t,
                         vector<double> V_f21h_p_t,
                         vector<double> V_f21i_p_t,
                         vector<double> V_f21j_p_t,
                         vector<double> V_f21k_p_t,
                         vector<double> V_f21l_p_t,
                         vector<double> V_f21m_p_t,
                         vector<double> V_f21n_p_t,
                         vector<double> V_f21a_m_t,
                         vector<double> V_f21b_m_t,
                         vector<double> V_f21c_m_t,
                         vector<double> V_f21d_m_t,
                         vector<double> V_f21e_m_t,
                         vector<double> V_f21f_m_t,
                         vector<double> V_f21g_m_t,
                         vector<double> V_f21h_m_t,
                         vector<double> V_f21i_m_t,
                         vector<double> V_f21j_m_t,
                         vector<double> V_f21k_m_t,
                         vector<double> V_f21l_m_t,
                         vector<double> V_f21m_m_t,
                         vector<double> V_f21n_m_t,
                         vector<double> V_Hdens,
                         int SIZE){
    
    //Before everything I will specify whether or not
    //we will run the SMOOTHING DISTRIBUTION
    int SMOOTHINGFILTER=0;
    double RR=100;
    
    vector<double> LOGLIK;
    LOGLIK.resize(SIZE);
    //Defining the variables
    double  Cliveswithfather12=0;
    double Cliveswithmother12=0;
    double Hhchores12=0;
    double Meffort12=0;
    double Feffort12=0;
    double Mage12=0;
    double Fage12=0;
    double Cedad_meses12=0;
    double Ctestsfactorsss2012=0;
    double Ctestsfactor2ss_10=0;
    double Myrschool12=0;
    double Fyrschool12=0;
    double Ffraclabor12=0;
    double Mfraclabor12=0;
    double Mwage12=0;
    double Fwage12=0;
    double Mnly12=0;
    double Fnly12=0;
    double Cfactorinv12=0;
    double Cchildcare12=0;
    double Ccareskills12=0;
    double Hbarg=0;
    double Feffort10=0;
    double Meffort10=0;
    double Cfactorinv10=0;
    double Cbirthfactor=0;
    double Mwage10=0;
    double Fwage10=0;
    double Mnly10=0;
    double Fnly10=0;
    double Cedad_meses10=0;
    double Cliveswithfather10=0;
    double Cliveswithmother10=0;
    double Mage10=0;
    double Fage10=0;
    double Mfraclabor10=0;
    double Ffraclabor10=0;
    double Cchildcare10=0;
    double Hchildcareobs=0;
    double PG=0;
    double MTJH=0;
    double MRATIO=0;
    double Unemployment=0;
    double Wageratio=0;
    double  Distance=0;
    double Magegroup10=0;
    double  Magegroup12=0;
    double  Hmemberstotal10=0;
    double  Hmemberstotal12=0;
    double CSTDtepsi_pb_coo10=0;
    double CSTDtepsi_pb_len10=0;
    double CSTDtepsi_pb_mot10=0;
    double CSTDtvip_pb10=0;
    double CSTDcbcl1_pb_110=0;
    double CSTDcbcl1_pb_210=0;
    double CSTDcbcl1_pb_310=0;
    double CSTDcbcl1_pb_410=0;
    double CSTDcbcl1_pb_510=0;
    double CSTDcbcl1_pb_610=0;
    double CSTDcbcl1_pb_710=0;
    double CSTDtadi_pb_cog12=0;
    double CSTDtadi_pb_mot12=0;
    double CSTDtadi_pb_len12=0;
    double CSTDtadi_pb_se12=0;
    double CSTDbt_112=0;
    double CSTDbt_212=0;
    double CSTDbt_312=0;
    double CSTDbt_412=0;
    double CSTDbt_512=0;
    double CSTDbt_t12=0;
    double CSTDhtks_st12=0;
    double CSTDbdst_st12=0;
    double CSTDppvt_t12=0;
    int Ccondpregg3_1=0;
    int Ccondpregg3_2=0;
    int Ccondpregg3_3=0;
    int Ccondpregg3_4=0;
    int Ccondpregg3_5=0;
    int Ccondpregg3_6=0;
    int Ccondpregg3_7=0;
    int Ccondpregg3_8=0;
    int Ccondpregg3_9=0;
    int Ccondpregg4a_1=0;
    int Ccondpregg4a_2=0;
    int Ccondpregg4a_3=0;
    int Ccondpregg4a_4=0;
    int Ccondpregg4a_5=0;
    int Ccondpregg4a_6=0;
    int Ccondpregg4a_7=0;
    double Ccondpregg7b=0;
    double Ccondpregg8b=0;
    double Ccondpregg9=0;
    double Ccondpregg11b=0;
    int Ccondpregpreterm=0;
    double Ccondpregg24=0;
    double Ccondpregg23=0;
    double Cwais_pb_num=0;
    double Cwais_pb_vo=0;
    double Cbfi_pb_ama=0;
    double Cbfi_pb_ape=0;
    double Cbfi_pb_ext=0;
    double Cbfi_pb_neu=0;
    double Cbfi_pb_res=0;
    double Cpsi_pb_total=0;
    double Hbargg2a=0;
    double Hbargg2b=0;
    double Hbargg2c=0;
    double Hbargg2d=0;
    double Hbargg2e=0;
    double Hbargg2f=0;
    double Hbargg2g=0;
    double Hbargg2h=0;
    double Hbargg2i=0;
    double Hbargg2j=0;
    int Hwomdecide12=0;
    int Hfathdecides12=0;
    int Hbothdecide=0;
    double Hcaresacwom=0;
    double Hcaresacman=0;
    int Hopofman121=0;
    int Hopofman122=0;
    int Hopofman123=0;
    int Hopofman124=0;
    int Cinvh3=0;
    int Cinvh4=0;
    int Cinvh5=0;
    int Cinvh6=0;
    int Cinvh7=0;
    int Cinvh8=0;
    int Cinvh9=0;
    int Cinvh13=0;
    double Cinvf11a=0;
    double Cinvf11b=0;
    double Cinvf11c=0;
    double Cinvf11d=0;
    double Cinvf11e=0;
    double Cinvf11f=0;
    double Cinvf11g=0;
    double Cinvf11h=0;
    double Cinvf11i=0;
    double Cinvf11j=0;
    double Cinvf11k=0;
    int Cinvhm2_10=0;
    int Cinvhm2_11=0;
    int Cinvhm2_12=0;
    int Cinvhm2_13=0;
    int Cinvhm2_14=0;
    int Cinvhm2_15=0;
    int Cinvhm2_16=0;
    int Cinvhm2_18=0;
    double Csharesbedroomhowmany12=0;
    double Csharesbedhowmany12=0;
    int g42_a2=0;
    int g42_b2=0;
    int g42_c2=0;
    int g42_d2=0;
    int g42_e2=0;
    int g42_f2=0;
    int g42_a1=0;
    int g42_b1=0;
    int g42_c1=0;
    int g42_d1=0;
    int g42_e1=0;
    int g42_f1=0;
    double f21a_p_t=0;
    double f21b_p_t=0;
    double f21c_p_t=0;
    double f21d_p_t=0;
    double f21e_p_t=0;
    double f21f_p_t=0;
    double f21g_p_t=0;
    double f21h_p_t=0;
    double f21i_p_t=0;
    double f21j_p_t=0;
    double f21k_p_t=0;
    double f21l_p_t=0;
    double f21m_p_t=0;
    double f21n_p_t=0;
    double f21a_m_t=0;
    double f21b_m_t=0;
    double f21c_m_t=0;
    double f21d_m_t=0;
    double f21e_m_t=0;
    double f21f_m_t=0;
    double f21g_m_t=0;
    double f21h_m_t=0;
    double f21i_m_t=0;
    double f21j_m_t=0;
    double f21k_m_t=0;
    double f21l_m_t=0;
    double f21m_m_t=0;
    double f21n_m_t=0;
    double HDens=0;
    
    //Defining the output of the function scalar defined before
    vector<vector<long double> > OUTPUT;
    OUTPUT.resize(RR);
    for (int n=0; n<RR; n++){
        OUTPUT[n].resize(7);
    }
    
    //If running the smoothing filter we need to define additional stuff:
    
    //DEFINING THE NECESSARY VARIABLES
    
    
    vector<vector<double> > S0MATRIX;
    S0MATRIX.resize(SIZE);
    for (int n=0; n<SIZE; n++){
        S0MATRIX[n].resize(RR);
    }
    
    vector<vector<double> > S1MATRIX;
    S1MATRIX.resize(SIZE);
    for (int n=0; n<SIZE; n++){
        S1MATRIX[n].resize(RR);
    }
    
    vector<vector<double> > S2MATRIX;
    S2MATRIX.resize(SIZE);
    for (int n=0; n<SIZE; n++){
        S2MATRIX[n].resize(RR);
    }
    
    
    vector<vector<double> > W0MATRIX;
    W0MATRIX.resize(SIZE);
    for (int n=0; n<SIZE; n++){
        W0MATRIX[n].resize(RR);
    }
    
    vector<vector<double> > W1MATRIX;
    W1MATRIX.resize(SIZE);
    for (int n=0; n<SIZE; n++){
        W1MATRIX[n].resize(RR);
    }
    
    vector<vector<double> > W2MATRIX;
    W2MATRIX.resize(SIZE);
    for (int n=0; n<SIZE; n++){
        W2MATRIX[n].resize(RR);
    }
    
    
    
    
    
    
    
#pragma omp parallel for shared(LOGLIK, par_likelihood) private(Cliveswithfather12, Cliveswithmother12, Hhchores12, Meffort12, Feffort12, Mage12, Fage12, Cedad_meses12, Ctestsfactorsss2012, Ctestsfactor2ss_10, Myrschool12, Fyrschool12, Ffraclabor12, Mfraclabor12, Mwage12, Fwage12, Mnly12, Fnly12, Cfactorinv12, Cchildcare12, Ccareskills12, Hbarg, Feffort10, Meffort10, Cfactorinv10, Cbirthfactor, Mwage10, Fwage10, Mnly10, Fnly10, Cedad_meses10, Cliveswithfather10, Cliveswithmother10, Mage10, Fage10, Mfraclabor10, Ffraclabor10, Cchildcare10, Hchildcareobs, PG, MTJH,MRATIO,Unemployment, Wageratio,Distance,Magegroup10,Magegroup12,Hmemberstotal10,Hmemberstotal12,CSTDtepsi_pb_coo10,CSTDtepsi_pb_len10,CSTDtepsi_pb_mot10,CSTDtvip_pb10,CSTDcbcl1_pb_110,CSTDcbcl1_pb_210,CSTDcbcl1_pb_310,CSTDcbcl1_pb_410,CSTDcbcl1_pb_510,CSTDcbcl1_pb_610,CSTDcbcl1_pb_710,CSTDtadi_pb_cog12, CSTDtadi_pb_mot12, CSTDtadi_pb_len12, CSTDtadi_pb_se12, CSTDbt_112, CSTDbt_212, CSTDbt_312, CSTDbt_412, CSTDbt_512, CSTDbt_t12,  CSTDhtks_st12, CSTDbdst_st12, CSTDppvt_t12,Ccondpregg3_1, Ccondpregg3_2 ,Ccondpregg3_3, Ccondpregg3_4, Ccondpregg3_5, Ccondpregg3_6, Ccondpregg3_7, Ccondpregg3_8, Ccondpregg3_9, Ccondpregg4a_1, Ccondpregg4a_2, Ccondpregg4a_3, Ccondpregg4a_4, Ccondpregg4a_5, Ccondpregg4a_6, Ccondpregg4a_7, Ccondpregg7b, Ccondpregg8b, Ccondpregg9, Ccondpregg11b, Ccondpregpreterm, Ccondpregg24, Ccondpregg23,Cwais_pb_num, Cwais_pb_vo, Cbfi_pb_ama,  Cbfi_pb_ape, Cbfi_pb_ext, Cbfi_pb_neu, Cbfi_pb_res, Cpsi_pb_total,Hbargg2a, Hbargg2b, Hbargg2c, Hbargg2d, Hbargg2e, Hbargg2f, Hbargg2g, Hbargg2h, Hbargg2i, Hbargg2j, Hwomdecide12, Hfathdecides12, Hbothdecide, Hcaresacwom, Hcaresacman,  Hopofman121,  Hopofman122,  Hopofman123,  Hopofman124, Cinvh3, Cinvh4, Cinvh5, Cinvh6, Cinvh7, Cinvh8, Cinvh9, Cinvh13,Cinvf11a, Cinvf11b,  Cinvf11c, Cinvf11d, Cinvf11e, Cinvf11f, Cinvf11g, Cinvf11h, Cinvf11i, Cinvf11j, Cinvf11k, Cinvhm2_10, Cinvhm2_11, Cinvhm2_12, Cinvhm2_13, Cinvhm2_14, Cinvhm2_15, Cinvhm2_16, Cinvhm2_18, Csharesbedroomhowmany12,  Csharesbedhowmany12,g42_a2, g42_b2, g42_c2, g42_d2, g42_e2, g42_f2,g42_a1, g42_b1, g42_c1, g42_d1, g42_e1, g42_f1,f21a_p_t, f21b_p_t, f21c_p_t, f21d_p_t, f21e_p_t, f21f_p_t, f21g_p_t, f21h_p_t, f21i_p_t, f21j_p_t, f21k_p_t, f21l_p_t, f21m_p_t, f21n_p_t,f21a_m_t, f21b_m_t, f21c_m_t, f21d_m_t, f21e_m_t, f21f_m_t, f21g_m_t, f21h_m_t, f21i_m_t, f21j_m_t, f21k_m_t, f21l_m_t, f21m_m_t, f21n_m_t,HDens,OUTPUT)
    for (int ii=0; ii<950; ii=ii+1){
        if(1==2){
            cout << " ---- " << endl;
            cout << ii << " ii in parallelizedloop"<<endl;
            cout << " ---- " << endl;
        }
        //Loading vectors to usual ones
        Cliveswithfather12=V_Cliveswithfather12[ii];
        Cliveswithmother12=V_Cliveswithmother12[ii];
        Hhchores12=V_Hhchores12[ii];
        Meffort12= V_Meffort12[ii];
        Feffort12=V_Feffort12[ii];
        Mage12=V_Mage12[ii];
        Fage12=V_Fage12[ii];
        
        Cedad_meses12=V_Cedad_meses12[ii];
        Ctestsfactorsss2012=V_Ctestsfactorsss2012[ii];
        Ctestsfactor2ss_10=V_Ctestsfactor2ss_10[ii];
        Myrschool12=V_Myrschool12[ii];
        Fyrschool12=V_Fyrschool12[ii];
        Ffraclabor12=V_Ffraclabor12[ii];
        Mfraclabor12=V_Mfraclabor12[ii];
        Mwage12=V_Mwage12[ii];
        Fwage12=V_Fwage12[ii];
        Mnly12=V_Mnly12[ii];
        Fnly12=V_Fnly12[ii];
        Cfactorinv12=V_Cfactorinv12[ii];
        Cchildcare12=V_Cchildcare12[ii];
        Ccareskills12=V_Ccareskills12[ii];
        Hbarg=V_Hbarg[ii];
        Feffort10=V_Feffort10[ii];
        Meffort10=V_Meffort10[ii];
        Cfactorinv10=V_Cfactorinv10[ii];
        Cbirthfactor=V_Cbirthfactor[ii];
        Mwage10=V_Mwage10[ii];
        Fwage10=V_Fwage10[ii];
        Mnly10=V_Mnly10[ii];
        Fnly10=V_Fnly10[ii];
        Cedad_meses10=V_Cedad_meses10[ii];
        Cliveswithfather10=V_Cliveswithfather10[ii];
        Cliveswithmother10=V_Cliveswithmother10[ii];
        Mage10=V_Mage10[ii];
        Fage10=V_Fage10[ii];
        Mfraclabor10=V_Mfraclabor10[ii];
        Ffraclabor10=V_Ffraclabor10[ii];
        Cchildcare10=V_Cchildcare10[ii];
        Hchildcareobs=V_Hchildcareobs[ii];
        PG=V_PG[ii];
        MTJH=V_MTJH[ii];
        MRATIO=V_MRATIO[ii];
        Unemployment=V_Unemployment[ii];
        Wageratio=V_Wageratio[ii];
        Distance=V_Distance[ii];
        Magegroup10=V_Magegroup10[ii];
        Magegroup12=V_Magegroup12[ii];
        Hmemberstotal10=V_Hmemberstotal10[ii];
        Hmemberstotal12=V_Hmemberstotal12[ii];
        CSTDtepsi_pb_coo10=V_CSTDtepsi_pb_coo10[ii];
        CSTDtepsi_pb_len10=V_CSTDtepsi_pb_len10[ii];
        CSTDtepsi_pb_mot10=V_CSTDtepsi_pb_mot10[ii];
        CSTDtvip_pb10=V_CSTDtvip_pb10[ii];
        CSTDcbcl1_pb_110=V_CSTDcbcl1_pb_110[ii];
        CSTDcbcl1_pb_210=V_CSTDcbcl1_pb_210[ii];
        CSTDcbcl1_pb_310=V_CSTDcbcl1_pb_310[ii];
        CSTDcbcl1_pb_410=V_CSTDcbcl1_pb_410[ii];
        CSTDcbcl1_pb_510=V_CSTDcbcl1_pb_510[ii];
        CSTDcbcl1_pb_610=V_CSTDcbcl1_pb_610[ii];
        CSTDcbcl1_pb_710=V_CSTDcbcl1_pb_710[ii];
        CSTDtadi_pb_cog12=V_CSTDtadi_pb_cog12[ii];
        CSTDtadi_pb_mot12=V_CSTDtadi_pb_mot12[ii];
        CSTDtadi_pb_len12=V_CSTDtadi_pb_len12[ii];
        CSTDtadi_pb_se12=V_CSTDtadi_pb_se12[ii];
        CSTDbt_112=V_CSTDbt_112[ii];
        CSTDbt_212=V_CSTDbt_212[ii];
        CSTDbt_312=V_CSTDbt_312[ii];
        CSTDbt_412=V_CSTDbt_412[ii];
        CSTDbt_512=V_CSTDbt_512[ii];
        CSTDbt_t12=V_CSTDbt_t12[ii];
        CSTDhtks_st12=V_CSTDhtks_st12[ii];
        CSTDbdst_st12=V_CSTDbdst_st12[ii];
        CSTDppvt_t12=V_CSTDppvt_t12[ii];
        Ccondpregg3_1=V_Ccondpregg3_1[ii];
        Ccondpregg3_2=V_Ccondpregg3_2[ii];
        Ccondpregg3_3=V_Ccondpregg3_3[ii];
        Ccondpregg3_4=V_Ccondpregg3_4[ii];
        Ccondpregg3_5=V_Ccondpregg3_5[ii];
        Ccondpregg3_6=V_Ccondpregg3_6[ii];
        Ccondpregg3_7=V_Ccondpregg3_7[ii];
        Ccondpregg3_8=V_Ccondpregg3_8[ii];
        Ccondpregg3_9=V_Ccondpregg3_9[ii];
        Ccondpregg4a_1=V_Ccondpregg4a_1[ii];
        Ccondpregg4a_2=V_Ccondpregg4a_2[ii];
        Ccondpregg4a_3=V_Ccondpregg4a_3[ii];
        Ccondpregg4a_4=V_Ccondpregg4a_4[ii];
        Ccondpregg4a_5=V_Ccondpregg4a_5[ii];
        Ccondpregg4a_6=V_Ccondpregg4a_6[ii];
        Ccondpregg4a_7=V_Ccondpregg4a_7[ii];
        Ccondpregg7b=V_Ccondpregg7b[ii];
        Ccondpregg8b=V_Ccondpregg8b[ii];
        Ccondpregg9=V_Ccondpregg9[ii];
        Ccondpregg11b=V_Ccondpregg11b[ii];
        Ccondpregpreterm=V_Ccondpregpreterm[ii];
        Ccondpregg24=V_Ccondpregg24[ii];
        Ccondpregg23=V_Ccondpregg23[ii];
        Cwais_pb_num=V_Cwais_pb_num[ii];
        Cwais_pb_vo=V_Cwais_pb_vo[ii];
        Cbfi_pb_ama=V_Cbfi_pb_ama[ii];
        Cbfi_pb_ape=V_Cbfi_pb_ape[ii];
        Cbfi_pb_ext=V_Cbfi_pb_ext[ii];
        Cbfi_pb_neu=V_Cbfi_pb_neu[ii];
        Cbfi_pb_res=V_Cbfi_pb_res[ii];
        Cpsi_pb_total=V_Cpsi_pb_total[ii];
        Hbargg2a=V_Hbargg2a[ii];
        Hbargg2b=V_Hbargg2b[ii];
        Hbargg2c=V_Hbargg2c[ii];
        Hbargg2d=V_Hbargg2d[ii];
        Hbargg2e=V_Hbargg2e[ii];
        Hbargg2f=V_Hbargg2f[ii];
        Hbargg2g=V_Hbargg2g[ii];
        Hbargg2h=V_Hbargg2h[ii];
        Hbargg2i=V_Hbargg2i[ii];
        Hbargg2j=V_Hbargg2j[ii];
        Hwomdecide12=V_Hwomdecide12[ii];
        Hfathdecides12=V_Hfathdecides12[ii];
        Hbothdecide=V_Hbothdecide[ii];
        Hcaresacwom=V_Hcaresacwom[ii];
        Hcaresacman=V_Hcaresacman[ii];
        Hopofman121=V_Hopofman121[ii];
        Hopofman122=V_Hopofman122[ii];
        Hopofman123=V_Hopofman123[ii];
        Hopofman124=V_Hopofman124[ii];
        Cinvh3=V_Cinvh3[ii];
        Cinvh4=V_Cinvh4[ii];
        Cinvh5=V_Cinvh5[ii];
        Cinvh6=V_Cinvh6[ii];
        Cinvh7=V_Cinvh7[ii];
        Cinvh8=V_Cinvh8[ii];
        Cinvh9=V_Cinvh9[ii];
        Cinvh13=V_Cinvh13[ii];
        Cinvf11a=V_Cinvf11a[ii];
        Cinvf11b=V_Cinvf11b[ii];
        Cinvf11c=V_Cinvf11c[ii];
        Cinvf11d=V_Cinvf11d[ii];
        Cinvf11e=V_Cinvf11e[ii];
        Cinvf11f=V_Cinvf11f[ii];
        Cinvf11g=V_Cinvf11g[ii];
        Cinvf11h=V_Cinvf11h[ii];
        Cinvf11i=V_Cinvf11i[ii];
        Cinvf11j=V_Cinvf11j[ii];
        Cinvf11k=V_Cinvf11k[ii];
        Cinvhm2_10=V_Cinvhm2_10[ii];
        Cinvhm2_11=V_Cinvhm2_11[ii];
        Cinvhm2_12=V_Cinvhm2_12[ii];
        Cinvhm2_13=V_Cinvhm2_13[ii];
        Cinvhm2_14=V_Cinvhm2_14[ii];
        Cinvhm2_15=V_Cinvhm2_15[ii];
        Cinvhm2_16=V_Cinvhm2_16[ii];
        Cinvhm2_18=V_Cinvhm2_18[ii];
        Csharesbedroomhowmany12=V_Csharesbedroomhowmany12[ii];
        Csharesbedhowmany12=V_Csharesbedhowmany12[ii];
        g42_a2=V_g42_a2[ii];
        g42_b2=V_g42_b2[ii];
        g42_c2=V_g42_c2[ii];
        g42_d2=V_g42_d2[ii];
        g42_e2=V_g42_e2[ii];
        g42_f2=V_g42_f2[ii];
        g42_a1=V_g42_a1[ii];
        g42_b1=V_g42_b1[ii];
        g42_c1=V_g42_c1[ii];
        g42_d1=V_g42_d1[ii];
        g42_e1=V_g42_e1[ii];
        g42_f1=V_g42_f1[ii];
        f21a_p_t=V_f21a_p_t[ii];
        f21b_p_t=V_f21b_p_t[ii];
        f21c_p_t=V_f21c_p_t[ii];
        f21d_p_t=V_f21d_p_t[ii];
        f21e_p_t=V_f21e_p_t[ii];
        f21f_p_t=V_f21f_p_t[ii];
        f21g_p_t=V_f21g_p_t[ii];
        f21h_p_t=V_f21h_p_t[ii];
        f21i_p_t=V_f21i_p_t[ii];
        f21j_p_t=V_f21j_p_t[ii];
        f21k_p_t=V_f21k_p_t[ii];
        f21l_p_t=V_f21l_p_t[ii];
        f21m_p_t=V_f21m_p_t[ii];
        f21n_p_t=V_f21n_p_t[ii];
        f21a_m_t=V_f21a_m_t[ii];
        f21b_m_t=V_f21b_m_t[ii];
        f21c_m_t=V_f21c_m_t[ii];
        f21d_m_t=V_f21d_m_t[ii];
        f21e_m_t=V_f21e_m_t[ii];
        f21f_m_t=V_f21f_m_t[ii];
        f21g_m_t=V_f21g_m_t[ii];
        f21h_m_t=V_f21h_m_t[ii];
        f21i_m_t=V_f21i_m_t[ii];
        f21j_m_t=V_f21j_m_t[ii];
        f21k_m_t=V_f21k_m_t[ii];
        f21l_m_t=V_f21l_m_t[ii];
        f21m_m_t=V_f21m_m_t[ii];
        f21n_m_t=V_f21n_m_t[ii];
        HDens=V_Hdens[ii];
        
        int SEED=2581633+ii;
        OUTPUT=F_likelihoodSCALAR(par_likelihood, Cliveswithfather12, Cliveswithmother12, Hhchores12, Meffort12, Feffort12, Mage12, Fage12, Cedad_meses12, Ctestsfactorsss2012, Ctestsfactor2ss_10, Myrschool12, Fyrschool12, Ffraclabor12, Mfraclabor12, Mwage12, Fwage12, Mnly12, Fnly12, Cfactorinv12, Cchildcare12, Ccareskills12, Hbarg, Feffort10, Meffort10, Cfactorinv10, Cbirthfactor, Mwage10, Fwage10, Mnly10, Fnly10, Cedad_meses10, Cliveswithfather10, Cliveswithmother10, Mage10, Fage10, Mfraclabor10, Ffraclabor10, Cchildcare10, Hchildcareobs, PG, MTJH,MRATIO,Unemployment, Wageratio,Distance,Magegroup10,Magegroup12,Hmemberstotal10,Hmemberstotal12,CSTDtepsi_pb_coo10,CSTDtepsi_pb_len10,CSTDtepsi_pb_mot10,CSTDtvip_pb10,CSTDcbcl1_pb_110,CSTDcbcl1_pb_210,CSTDcbcl1_pb_310,CSTDcbcl1_pb_410,CSTDcbcl1_pb_510,CSTDcbcl1_pb_610,CSTDcbcl1_pb_710,CSTDtadi_pb_cog12, CSTDtadi_pb_mot12, CSTDtadi_pb_len12, CSTDtadi_pb_se12, CSTDbt_112, CSTDbt_212, CSTDbt_312, CSTDbt_412, CSTDbt_512, CSTDbt_t12,  CSTDhtks_st12, CSTDbdst_st12, CSTDppvt_t12,Ccondpregg3_1, Ccondpregg3_2 ,Ccondpregg3_3, Ccondpregg3_4, Ccondpregg3_5, Ccondpregg3_6, Ccondpregg3_7, Ccondpregg3_8, Ccondpregg3_9, Ccondpregg4a_1, Ccondpregg4a_2, Ccondpregg4a_3, Ccondpregg4a_4, Ccondpregg4a_5, Ccondpregg4a_6, Ccondpregg4a_7, Ccondpregg7b, Ccondpregg8b, Ccondpregg9, Ccondpregg11b, Ccondpregpreterm, Ccondpregg24, Ccondpregg23,Cwais_pb_num, Cwais_pb_vo, Cbfi_pb_ama,  Cbfi_pb_ape, Cbfi_pb_ext, Cbfi_pb_neu, Cbfi_pb_res, Cpsi_pb_total,Hbargg2a, Hbargg2b, Hbargg2c, Hbargg2d, Hbargg2e, Hbargg2f, Hbargg2g, Hbargg2h, Hbargg2i, Hbargg2j, Hwomdecide12, Hfathdecides12, Hbothdecide, Hcaresacwom, Hcaresacman,  Hopofman121,  Hopofman122,  Hopofman123,  Hopofman124, Cinvh3, Cinvh4, Cinvh5, Cinvh6, Cinvh7, Cinvh8, Cinvh9, Cinvh13,Cinvf11a, Cinvf11b,  Cinvf11c, Cinvf11d, Cinvf11e, Cinvf11f, Cinvf11g, Cinvf11h, Cinvf11i, Cinvf11j, Cinvf11k, Cinvhm2_10, Cinvhm2_11, Cinvhm2_12, Cinvhm2_13, Cinvhm2_14, Cinvhm2_15, Cinvhm2_16, Cinvhm2_18, Csharesbedroomhowmany12,  Csharesbedhowmany12,g42_a2, g42_b2, g42_c2, g42_d2, g42_e2, g42_f2,g42_a1, g42_b1, g42_c1, g42_d1, g42_e1, g42_f1,f21a_p_t, f21b_p_t, f21c_p_t, f21d_p_t, f21e_p_t, f21f_p_t, f21g_p_t, f21h_p_t, f21i_p_t, f21j_p_t, f21k_p_t, f21l_p_t, f21m_p_t, f21n_p_t,f21a_m_t, f21b_m_t, f21c_m_t, f21d_m_t, f21e_m_t, f21f_m_t, f21g_m_t, f21h_m_t, f21i_m_t, f21j_m_t, f21k_m_t, f21l_m_t, f21m_m_t, f21n_m_t,HDens,SEED);
        LOGLIK[ii]=OUTPUT[0][0];
        
        
        //If I run the smoothing distribution, generate the csv properly:
        if (SMOOTHINGFILTER==1){
            for(int rr=0; rr<RR; rr=rr+1){
                S0MATRIX[ii][rr]=OUTPUT[rr][1];
                S1MATRIX[ii][rr]=OUTPUT[rr][2];
                S2MATRIX[ii][rr]=OUTPUT[rr][3];
                
                W0MATRIX[ii][rr]=OUTPUT[rr][4];
                W1MATRIX[ii][rr]=OUTPUT[rr][5];
                W2MATRIX[ii][rr]=OUTPUT[rr][6];
                if(1==2){
                    cout << S2MATRIX[ii][rr] << " W2MATRIX[ii][rr] "<< endl;
                }
                
                
                
            }// END RR LOOP
            
            
        }//end if smoothing filter distribution
        
        
        
    }// end ii loop
    double FIN=0;
    //And adding up all the vector
    for (int JJ=0; JJ<950; JJ=JJ+1){
        FIN+=LOGLIK[JJ];
        
    }
    
    //And running the loading parameters of smoothing distribution
    //We need to do this in a loop different to the previous one because
    //the previous one was parallelized. We need to do it not in parallel to not generate
    //troubles with storing the csv file.
    
    
    if(SMOOTHINGFILTER==1){
        //DEFINING THE OUTPUT WHERE ALL ELEMENTS WILL BE STORED
        ofstream optparamS0("S0MATRIX.csv");
        ofstream optparamS1("S1MATRIX.csv");
        ofstream optparamS2("S2MATRIX.csv");
        ofstream optparamW0("W0MATRIX.csv");
        ofstream optparamW1("W1MATRIX.csv");
        ofstream optparamW2("W2MATRIX.csv");
        for(int ii=0; ii<950; ii=ii+1){
            for(int rr=0; rr<RR; rr=rr+1){
                
                optparamS0<<S0MATRIX[ii][rr] << " , " ;
                optparamS1<<S1MATRIX[ii][rr] << " , " ;
                optparamS2<<S2MATRIX[ii][rr] << " , " ;
                optparamW0<<W0MATRIX[ii][rr] << " , " ;
                optparamW1<<W1MATRIX[ii][rr] << " , " ;
                optparamW2<<W2MATRIX[ii][rr] << " , " ;
                
            }// END RR LOOP
            optparamS0 << " " << endl;
            optparamS1 << " " << endl;
            optparamS2 << " " << endl;
            optparamW0 << " " << endl;
            optparamW1 << " " << endl;
            optparamW2 << " " << endl;
            
        }//end JJ loop
    }
    
    
    
    return(FIN);
}



//====================
//Trying to define likelihood function exclusively as function of parameters
long double F_likelihood_FIN(const double *PARDOUBLE){
    //0. Define the current directory
    
    //chdir("/Users/rodrigoazuero/Dropbox/BACKUPRODRIGO/Research/Chile/RR/BEHAVIORAL70");
    
    int NVAR=209; //Number of variables in file
    int SIZEOBS=950; //Number of observations in file
    
    //First step making the vector<double> from the const(double *par)
    vector<double> PAR;
    PAR.resize(325);
    cout << "----- "<< endl;
    
    for (int it=0; it<325; it=it+1){
        PAR[it]=PARDOUBLE[324-it];
        //PAR[3]=0;
        //PAR[4]=0;
        //PAR[59]=-100.5;
        //PAR[60]=0;
        //PAR[59]=-10;
        //PAR[5]=-100;
        //PAR[6]=-100;
        cout << PAR[it] << " PARLOADING " << endl;
    }
    
    
    //========================================
    //BLOCK OF GETTING THE DATA INTO VECTORS
    //=========================================
    std::ifstream theFile ("BEHAVIORAL.csv");
    double MYARRAY[SIZEOBS+1][NVAR+1];
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
    
    //Defining the vectors
    vector<double> Hmemberstotal12;
    Hmemberstotal12.resize(SIZEOBS);
    
    vector<double> Hmemberstotal10;
    Hmemberstotal10.resize(SIZEOBS);
    
    vector<double> Cedad_meses10;
    Cedad_meses10.resize(SIZEOBS);
    
    vector<double> Ctestsfactor1_10;
    Ctestsfactor1_10.resize(SIZEOBS);
    
    vector<double> Mfraclabor10;
    Mfraclabor10.resize(SIZEOBS);
    
    vector<double> Mwage10;
    Mwage10.resize(SIZEOBS);
    
    vector<double> Mnlincome10;
    Mnlincome10.resize(SIZEOBS);
    
    vector<double> Ffraclabor10;
    Ffraclabor10.resize(SIZEOBS);
    
    vector<double> Fwage10;
    Fwage10.resize(SIZEOBS);
    
    vector<double> Fnlincome10;
    Fnlincome10.resize(SIZEOBS);
    
    vector<double> Cfactorbirth;
    Cfactorbirth.resize(SIZEOBS);
    
    vector<double> Meffort10;
    Meffort10.resize(SIZEOBS);
    
    vector<double> Feffort10;
    Feffort10.resize(SIZEOBS);
    
    vector<double> CfactorInv10;
    CfactorInv10.resize(SIZEOBS);
    
    vector<double> Cchildcare10;
    Cchildcare10.resize(SIZEOBS);
    
    vector<double> Cliveswithfather10;
    Cliveswithfather10.resize(SIZEOBS);
    
    vector<double> Cliveswithmother10;
    Cliveswithmother10.resize(SIZEOBS);
    
    vector<double> Cedad_meses12;
    Cedad_meses12.resize(SIZEOBS);
    
    vector<double> Ctestsfactor4_2012;
    Ctestsfactor4_2012.resize(SIZEOBS);
    
    vector<double> Cliveswithfather12;
    Cliveswithfather12.resize(SIZEOBS);
    
    vector<double> Cliveswithmother12;
    Cliveswithmother12.resize(SIZEOBS);
    
    vector<double> Mage12;
    Mage12.resize(SIZEOBS);
    
    vector<double> Myrschool12;
    Myrschool12.resize(SIZEOBS);
    
    vector<double> Mfraclabor12;
    Mfraclabor12.resize(SIZEOBS);
    
    vector<double> Fage12;
    Fage12.resize(SIZEOBS);
    
    vector<double> Fyrschool12;
    Fyrschool12.resize(SIZEOBS);
    
    vector<double> Ffraclabor12;
    Ffraclabor12.resize(SIZEOBS);
    
    vector<double> Mwage12;
    Mwage12.resize(SIZEOBS);
    
    vector<double> Mnlincome12;
    Mnlincome12.resize(SIZEOBS);
    
    vector<double> Fwage12;
    Fwage12.resize(SIZEOBS);
    
    vector<double> Fnlincome12;
    Fnlincome12.resize(SIZEOBS);
    
    vector<double> Hbarg;
    Hbarg.resize(SIZEOBS);
    
    vector<double> Hbarg1;
    Hbarg1.resize(SIZEOBS);
    
    vector<double> Hbarg2;
    Hbarg2.resize(SIZEOBS);
    
    vector<double> Hbarg3;
    Hbarg3.resize(SIZEOBS);
    
    vector<double> Hbarg4;
    Hbarg4.resize(SIZEOBS);
    
    vector<double> Meffort12;
    Meffort12.resize(SIZEOBS);
    
    vector<double> Feffort12;
    Feffort12.resize(SIZEOBS);
    
    vector<double> CfactorInv12;
    CfactorInv12.resize(SIZEOBS);
    
    vector<double> Hchores12;
    Hchores12.resize(SIZEOBS);
    
    vector<double> Cchildcare12;
    Cchildcare12.resize(SIZEOBS);
    
    vector<double> Ccareskills;
    Ccareskills.resize(SIZEOBS);
    
    vector<double> Hchildcareobs;
    Hchildcareobs.resize(SIZEOBS);
    
    vector<double> FMRATIO;
    FMRATIO.resize(SIZEOBS);
    
    vector<double> Unemployment;
    Unemployment.resize(SIZEOBS);
    
    vector<double> Wageratio;
    Wageratio.resize(SIZEOBS);
    
    vector<double> Distance;
    Distance.resize(SIZEOBS);
    
    vector<double> MTJH;
    MTJH.resize(SIZEOBS);
    
    vector<double> Magegroup10;
    Magegroup10.resize(SIZEOBS);
    
    vector<double> Magegroup12;
    Magegroup12.resize(SIZEOBS);
    
    vector<double> Fage10;
    Fage10.resize(SIZEOBS);
    
    vector<double> Mage10;
    Mage10.resize(SIZEOBS);
    
    vector<double> Mage2_12;
    Mage2_12.resize(SIZEOBS);
    
    vector<double> Fage2_12;
    Fage2_12.resize(SIZEOBS);
    
    vector<double> Mage2_10;
    Mage2_10.resize(SIZEOBS);
    
    vector<double> Fage2_10;
    Fage2_10.resize(SIZEOBS);
    
    vector<double> Mwagepred;
    Mwagepred.resize(SIZEOBS);
    
    vector<double> Fwagepred;
    Fwagepred.resize(SIZEOBS);
    
    vector<double> Agedif;
    Agedif.resize(SIZEOBS);
    
    vector<double> Edudif;
    Edudif.resize(SIZEOBS);
    
    vector<double> ymratio;
    ymratio.resize(SIZEOBS);
    
    vector<double> ymdif;
    ymdif.resize(SIZEOBS);
    
    vector<double> wageratio;
    wageratio.resize(SIZEOBS);
    
    //The measurements of skills in 2010
    vector<double> CSTDtepsi_pb_coo10;
    CSTDtepsi_pb_coo10.resize(SIZEOBS);
    
    vector<double> CSTDtepsi_pb_len10;
    CSTDtepsi_pb_len10.resize(SIZEOBS);
    
    vector<double> CSTDtepsi_pb_mot10;
    CSTDtepsi_pb_mot10.resize(SIZEOBS);
    
    vector<double> CSTDtvip_pb10;
    CSTDtvip_pb10.resize(SIZEOBS);
    
    vector<double> CSTDcbcl1_pb_110;
    CSTDcbcl1_pb_110.resize(SIZEOBS);
    
    vector<double> CSTDcbcl1_pb_210;
    CSTDcbcl1_pb_210.resize(SIZEOBS);
    
    vector<double> CSTDcbcl1_pb_310;
    CSTDcbcl1_pb_310.resize(SIZEOBS);
    
    vector<double> CSTDcbcl1_pb_410;
    CSTDcbcl1_pb_410.resize(SIZEOBS);
    
    vector<double> CSTDcbcl1_pb_510;
    CSTDcbcl1_pb_510.resize(SIZEOBS);
    
    vector<double> CSTDcbcl1_pb_610;
    CSTDcbcl1_pb_610.resize(SIZEOBS);
    
    vector<double> CSTDcbcl1_pb_710;
    CSTDcbcl1_pb_710.resize(SIZEOBS);
    
    //Measures of skills in 2012
    vector<double> CSTDtadi_pb_cog12;
    CSTDtadi_pb_cog12.resize(SIZEOBS);
    
    vector<double> CSTDtadi_pb_mot12;
    CSTDtadi_pb_mot12.resize(SIZEOBS);
    
    vector<double> CSTDtadi_pb_len12;
    CSTDtadi_pb_len12.resize(SIZEOBS);
    
    vector<double> CSTDtadi_pb_se12;
    CSTDtadi_pb_se12.resize(SIZEOBS);
    
    vector<double> CSTDbt_112;
    CSTDbt_112.resize(SIZEOBS);
    
    vector<double> CSTDbt_212;
    CSTDbt_212.resize(SIZEOBS);
    
    vector<double> CSTDbt_312;
    CSTDbt_312.resize(SIZEOBS);
    
    vector<double> CSTDbt_412;
    CSTDbt_412.resize(SIZEOBS);
    
    
    vector<double> CSTDbt_512;
    CSTDbt_512.resize(SIZEOBS);
    
    vector<double> CSTDbt_t12;
    CSTDbt_t12.resize(SIZEOBS);
    
    vector<double> CSTDhtks_st12;
    CSTDhtks_st12.resize(SIZEOBS);
    
    vector<double> CSTDbdst_st12;
    CSTDbdst_st12.resize(SIZEOBS);
    
    vector<double> CSTDppvt_t12;
    CSTDppvt_t12.resize(SIZEOBS);
    
    //Birth outcomes
    
    vector<int> Ccondpregg3_1;
    Ccondpregg3_1.resize(SIZEOBS);
    
    vector<int> Ccondpregg3_2;
    Ccondpregg3_2.resize(SIZEOBS);
    
    vector<int> Ccondpregg3_3;
    Ccondpregg3_3.resize(SIZEOBS);
    
    vector<int> Ccondpregg3_4;
    Ccondpregg3_4.resize(SIZEOBS);
    
    vector<int> Ccondpregg3_5;
    Ccondpregg3_5.resize(SIZEOBS);
    
    vector<int> Ccondpregg3_6;
    Ccondpregg3_6.resize(SIZEOBS);
    
    vector<int> Ccondpregg3_7;
    Ccondpregg3_7.resize(SIZEOBS);
    
    vector<int> Ccondpregg3_8;
    Ccondpregg3_8.resize(SIZEOBS);
    
    vector<int> Ccondpregg3_9;
    Ccondpregg3_9.resize(SIZEOBS);
    
    vector<int> Ccondpregg4a_1;
    Ccondpregg4a_1.resize(SIZEOBS);
    
    vector<int> Ccondpregg4a_2;
    Ccondpregg4a_2.resize(SIZEOBS);
    
    vector<int> Ccondpregg4a_3;
    Ccondpregg4a_3.resize(SIZEOBS);
    
    vector<int> Ccondpregg4a_4;
    Ccondpregg4a_4.resize(SIZEOBS);
    
    vector<int> Ccondpregg4a_5;
    Ccondpregg4a_5.resize(SIZEOBS);
    
    vector<int> Ccondpregg4a_6;
    Ccondpregg4a_6.resize(SIZEOBS);
    
    vector<int> Ccondpregg4a_7;
    Ccondpregg4a_7.resize(SIZEOBS);
    
    vector<double> Ccondpregg7b;
    Ccondpregg7b.resize(SIZEOBS);
    
    vector<double> Ccondpregg8b;
    Ccondpregg8b.resize(SIZEOBS);
    
    vector<double> Ccondpregg9;
    Ccondpregg9.resize(SIZEOBS);
    
    vector<double> Ccondpregg11b;
    Ccondpregg11b.resize(SIZEOBS);
    
    vector<int> Ccondpregpreterm;
    Ccondpregpreterm.resize(SIZEOBS);
    
    vector<double> Ccondpregg24;
    Ccondpregg24.resize(SIZEOBS);
    
    vector<double> Ccondpregg23;
    Ccondpregg23.resize(SIZEOBS);
    
    
    //Skills of primary caregiver
    vector<double> Cwais_pb_num;
    Cwais_pb_num.resize(SIZEOBS);
    
    vector<double> Cwais_pb_vo;
    Cwais_pb_vo.resize(SIZEOBS);
    
    vector<double> Cbfi_pb_ama;
    Cbfi_pb_ama.resize(SIZEOBS);
    
    vector<double> Cbfi_pb_ape;
    Cbfi_pb_ape.resize(SIZEOBS);
    
    vector<double> Cbfi_pb_ext;
    Cbfi_pb_ext.resize(SIZEOBS);
    
    vector<double> Cbfi_pb_neu;
    Cbfi_pb_neu.resize(SIZEOBS);
    
    vector<double> Cbfi_pb_res;
    Cbfi_pb_res.resize(SIZEOBS);
    
    vector<double> Cpsi_pb_total;
    Cpsi_pb_total.resize(SIZEOBS);
    
    
    //Bargaining power
    vector<double> Hbargg2a;
    Hbargg2a.resize(SIZEOBS);
    
    vector<double> Hbargg2b;
    Hbargg2b.resize(SIZEOBS);
    
    vector<double> Hbargg2c;
    Hbargg2c.resize(SIZEOBS);
    
    vector<double> Hbargg2d;
    Hbargg2d.resize(SIZEOBS);
    
    vector<double> Hbargg2e;
    Hbargg2e.resize(SIZEOBS);
    
    vector<double> Hbargg2f;
    Hbargg2f.resize(SIZEOBS);
    
    vector<double> Hbargg2g;
    Hbargg2g.resize(SIZEOBS);
    
    vector<double> Hbargg2h;
    Hbargg2h.resize(SIZEOBS);
    
    vector<double> Hbargg2i;
    Hbargg2i.resize(SIZEOBS);
    
    vector<double> Hbargg2j;
    Hbargg2j.resize(SIZEOBS);
    
    vector<int> Hwomdecide12;
    Hwomdecide12.resize(SIZEOBS);
    
    vector<int> Hfathdecides12;
    Hfathdecides12.resize(SIZEOBS);
    
    vector<int> Hbothdecide;
    Hbothdecide.resize(SIZEOBS);
    
    vector<double> Hcaresacwom;
    Hcaresacwom.resize(SIZEOBS);
    
    vector<double> Hcaresacman;
    Hcaresacman.resize(SIZEOBS);
    
    vector<int> Hopofman121;
    Hopofman121.resize(SIZEOBS);
    
    vector<int> Hopofman122;
    Hopofman122.resize(SIZEOBS);
    
    vector<int> Hopofman123;
    Hopofman123.resize(SIZEOBS);
    
    vector<int> Hopofman124;
    Hopofman124.resize(SIZEOBS);
    
    
    //6. Measures of investment in 2010. All of them are integers
    
    vector<int> Cinvh3;
    Cinvh3.resize(SIZEOBS);
    
    vector<int> Cinvh4;
    Cinvh4.resize(SIZEOBS);
    
    vector<int> Cinvh5;
    Cinvh5.resize(SIZEOBS);
    
    vector<int> Cinvh6;
    Cinvh6.resize(SIZEOBS);
    
    vector<int> Cinvh7;
    Cinvh7.resize(SIZEOBS);
    
    vector<int> Cinvh8;
    Cinvh8.resize(SIZEOBS);
    
    vector<int> Cinvh9;
    Cinvh9.resize(SIZEOBS);
    
    vector<int> Cinvh13;
    Cinvh13.resize(SIZEOBS);
    
    //7. Measures of investment in 2012. All of them are continuous
    vector<double> Cinvf11a;
    Cinvf11a.resize(SIZEOBS);
    
    vector<double> Cinvf11b;
    Cinvf11b.resize(SIZEOBS);
    
    vector<double> Cinvf11c;
    Cinvf11c.resize(SIZEOBS);
    
    vector<double> Cinvf11d;
    Cinvf11d.resize(SIZEOBS);
    
    vector<double> Cinvf11e;
    Cinvf11e.resize(SIZEOBS);
    
    vector<double> Cinvf11f;
    Cinvf11f.resize(SIZEOBS);
    
    vector<double> Cinvf11g;
    Cinvf11g.resize(SIZEOBS);
    
    vector<double> Cinvf11h;
    Cinvf11h.resize(SIZEOBS);
    
    vector<double> Cinvf11i;
    Cinvf11i.resize(SIZEOBS);
    
    vector<double> Cinvf11j;
    Cinvf11j.resize(SIZEOBS);
    
    vector<double> Cinvf11k;
    Cinvf11k.resize(SIZEOBS);
    
    vector<int> Cinvhm2_10;
    Cinvhm2_10.resize(SIZEOBS);
    
    vector<int> Cinvhm2_11;
    Cinvhm2_11.resize(SIZEOBS);
    
    vector<int> Cinvhm2_12;
    Cinvhm2_12.resize(SIZEOBS);
    
    vector<int> Cinvhm2_13;
    Cinvhm2_13.resize(SIZEOBS);
    
    vector<int> Cinvhm2_14;
    Cinvhm2_14.resize(SIZEOBS);
    
    vector<int> Cinvhm2_15;
    Cinvhm2_15.resize(SIZEOBS);
    
    vector<int> Cinvhm2_16;
    Cinvhm2_16.resize(SIZEOBS);
    
    vector<int> Cinvhm2_18;
    Cinvhm2_18.resize(SIZEOBS);
    
    vector<double> Csharesbedroomhowmany12;
    Csharesbedroomhowmany12.resize(SIZEOBS);
    
    vector<double> Csharesbedhowmany12;
    Csharesbedhowmany12.resize(SIZEOBS);
    
    //8. Time investments of the father. All of them are integers
    vector<int> g42_a2;
    g42_a2.resize(SIZEOBS);
    
    vector<int> g42_b2;
    g42_b2.resize(SIZEOBS);
    
    vector<int> g42_c2;
    g42_c2.resize(SIZEOBS);
    
    vector<int> g42_d2;
    g42_d2.resize(SIZEOBS);
    
    vector<int> g42_e2;
    g42_e2.resize(SIZEOBS);
    
    vector<int> g42_f2;
    g42_f2.resize(SIZEOBS);
    
    
    //9. Time investments of the mother. All of them are integers
    vector<int> g42_a1;
    g42_a1.resize(SIZEOBS);
    
    vector<int> g42_b1;
    g42_b1.resize(SIZEOBS);
    
    vector<int> g42_c1;
    g42_c1.resize(SIZEOBS);
    
    vector<int> g42_d1;
    g42_d1.resize(SIZEOBS);
    
    vector<int> g42_e1;
    g42_e1.resize(SIZEOBS);
    
    vector<int> g42_f1;
    g42_f1.resize(SIZEOBS);
    
    
    //10. Time investments of father in 2012. All of them are continuous.
    vector<double> f21a_p_t;
    f21a_p_t.resize(SIZEOBS);
    
    vector<double> f21b_p_t;
    f21b_p_t.resize(SIZEOBS);
    
    vector<double> f21c_p_t;
    f21c_p_t.resize(SIZEOBS);
    
    vector<double> f21d_p_t;
    f21d_p_t.resize(SIZEOBS);
    
    vector<double> f21e_p_t;
    f21e_p_t.resize(SIZEOBS);
    
    vector<double> f21f_p_t;
    f21f_p_t.resize(SIZEOBS);
    
    vector<double> f21g_p_t;
    f21g_p_t.resize(SIZEOBS);
    
    vector<double> f21h_p_t;
    f21h_p_t.resize(SIZEOBS);
    
    vector<double> f21i_p_t;
    f21i_p_t.resize(SIZEOBS);
    
    vector<double> f21j_p_t;
    f21j_p_t.resize(SIZEOBS);
    
    vector<double> f21k_p_t;
    f21k_p_t.resize(SIZEOBS);
    
    vector<double> f21l_p_t;
    f21l_p_t.resize(SIZEOBS);
    
    vector<double> f21m_p_t;
    f21m_p_t.resize(SIZEOBS);
    
    vector<double> f21n_p_t;
    f21n_p_t.resize(SIZEOBS);
    
    //11. Time investments of father in 2012. All of them are continuous.
    vector<double> f21a_m_t;
    f21a_m_t.resize(SIZEOBS);
    
    vector<double> f21b_m_t;
    f21b_m_t.resize(SIZEOBS);
    
    vector<double> f21c_m_t;
    f21c_m_t.resize(SIZEOBS);
    
    vector<double> f21d_m_t;
    f21d_m_t.resize(SIZEOBS);
    
    vector<double> f21e_m_t;
    f21e_m_t.resize(SIZEOBS);
    
    vector<double> f21f_m_t;
    f21f_m_t.resize(SIZEOBS);
    
    vector<double> f21g_m_t;
    f21g_m_t.resize(SIZEOBS);
    
    vector<double> f21h_m_t;
    f21h_m_t.resize(SIZEOBS);
    
    vector<double> f21i_m_t;
    f21i_m_t.resize(SIZEOBS);
    
    vector<double> f21j_m_t;
    f21j_m_t.resize(SIZEOBS);
    
    vector<double> f21k_m_t;
    f21k_m_t.resize(SIZEOBS);
    
    vector<double> f21l_m_t;
    f21l_m_t.resize(SIZEOBS);
    
    vector<double> f21m_m_t;
    f21m_m_t.resize(SIZEOBS);
    
    vector<double> f21n_m_t;
    f21n_m_t.resize(SIZEOBS);
    
    vector<double> HDens1;
    HDens1.resize(SIZEOBS);
    
    vector<double> HDens5;
    HDens5.resize(SIZEOBS);
    
    vector<double> HDens10;
    HDens10.resize(SIZEOBS);
    
    for (int it=0; it<SIZEOBS;it++){
        Hmemberstotal12[it]=MYARRAY[it][0];
        Hmemberstotal10[it]=MYARRAY[it][1];
        Cedad_meses10[it]=MYARRAY[it][2];
        Ctestsfactor1_10[it]=exp(MYARRAY[it][3]);
        Mfraclabor10[it]=MYARRAY[it][4];
        Mwage10[it]=MYARRAY[it][5];
        Mnlincome10[it]=MYARRAY[it][6];
        Ffraclabor10[it]=MYARRAY[it][7];
        Fwage10[it]=MYARRAY[it][8];
        Fnlincome10[it]=MYARRAY[it][9];
        Cfactorbirth[it]=exp(MYARRAY[it][10]);
        Meffort10[it]=MYARRAY[it][11];
        Feffort10[it]=MYARRAY[it][12];
        CfactorInv10[it]=exp(MYARRAY[it][13]);
        Cchildcare10[it]=MYARRAY[it][14];
        Cliveswithfather10[it]=MYARRAY[it][15];
        Cliveswithmother10[it]=MYARRAY[it][16];
        Cedad_meses12[it]=MYARRAY[it][17];
        Ctestsfactor4_2012[it]=exp(MYARRAY[it][18]);
        Cliveswithfather12[it]=MYARRAY[it][19];
        Cliveswithmother12[it]=MYARRAY[it][20];
        Mage12[it]=MYARRAY[it][21];
        Myrschool12[it]=MYARRAY[it][22];
        Mfraclabor12[it]=MYARRAY[it][23];
        Fage12[it]=MYARRAY[it][24];
        Fyrschool12[it]=MYARRAY[it][25];
        Ffraclabor12[it]=MYARRAY[it][26];
        Mwage12[it]=MYARRAY[it][27];
        Mnlincome12[it]=MYARRAY[it][28];
        Fwage12[it]=MYARRAY[it][29];
        Fnlincome12[it]=MYARRAY[it][30];
        Hbarg[it]=MYARRAY[it][31];
        Hbarg1[it]=MYARRAY[it][32];
        Hbarg2[it]=MYARRAY[it][33];
        Hbarg3[it]=MYARRAY[it][34];
        Hbarg4[it]=MYARRAY[it][35];
        Meffort12[it]=MYARRAY[it][36];
        Feffort12[it]=MYARRAY[it][37];
        CfactorInv12[it]=exp(MYARRAY[it][38]);
        Hchores12[it]=MYARRAY[it][39];
        Cchildcare12[it]=MYARRAY[it][40];
        Ccareskills[it]=MYARRAY[it][41];
        Hchildcareobs[it]=MYARRAY[it][42];
        FMRATIO[it]=MYARRAY[it][43];
        Unemployment[it]=MYARRAY[it][44];
        Wageratio[it]=MYARRAY[it][45];
        Distance[it]=MYARRAY[it][46];
        MTJH[it]=MYARRAY[it][47];
        Magegroup10[it]=MYARRAY[it][48];
        Magegroup12[it]=MYARRAY[it][49];
        Fage10[it]=MYARRAY[it][50];
        Mage10[it]=MYARRAY[it][51];
        Mage2_12[it]=MYARRAY[it][52];
        Fage2_12[it]=MYARRAY[it][53];
        Mage2_10[it]=MYARRAY[it][54];
        Fage2_10[it]=MYARRAY[it][55];
        Mwagepred[it]=MYARRAY[it][56];
        Fwagepred[it]=MYARRAY[it][57];
        Agedif[it]=MYARRAY[it][58];
        Edudif[it]=MYARRAY[it][59];
        ymratio[it]=MYARRAY[it][60];
        ymdif[it]=MYARRAY[it][61];
        wageratio[it]=MYARRAY[it][62];
        CSTDtepsi_pb_coo10[it]=MYARRAY[it][63];
        CSTDtepsi_pb_len10[it]=MYARRAY[it][64];
        CSTDtepsi_pb_mot10[it]=MYARRAY[it][65];
        CSTDtvip_pb10[it]=MYARRAY[it][66];
        CSTDcbcl1_pb_110[it]=MYARRAY[it][67];
        CSTDcbcl1_pb_210[it]=MYARRAY[it][68];
        CSTDcbcl1_pb_310[it]=MYARRAY[it][69];
        CSTDcbcl1_pb_410[it]=MYARRAY[it][70];
        CSTDcbcl1_pb_510[it]=MYARRAY[it][71];
        CSTDcbcl1_pb_610[it]=MYARRAY[it][72];
        CSTDcbcl1_pb_710[it]=MYARRAY[it][73];
        CSTDtadi_pb_cog12[it]=MYARRAY[it][74];
        CSTDtadi_pb_mot12[it]=MYARRAY[it][75];
        CSTDtadi_pb_len12[it]=MYARRAY[it][76];
        CSTDtadi_pb_se12[it]=MYARRAY[it][77];
        CSTDbt_112[it]=MYARRAY[it][78];
        CSTDbt_212[it]=MYARRAY[it][79];
        CSTDbt_312[it]=MYARRAY[it][80];
        CSTDbt_412[it]=MYARRAY[it][81];
        CSTDbt_512[it]=MYARRAY[it][82];
        CSTDbt_t12[it]=MYARRAY[it][83];
        CSTDhtks_st12[it]=MYARRAY[it][84];
        CSTDbdst_st12[it]=MYARRAY[it][85];
        CSTDppvt_t12[it]=MYARRAY[it][86];
        Ccondpregg3_1[it]=MYARRAY[it][87];
        Ccondpregg3_2[it]=MYARRAY[it][88];
        Ccondpregg3_3[it]=MYARRAY[it][89];
        Ccondpregg3_4[it]=MYARRAY[it][90];
        Ccondpregg3_5[it]=MYARRAY[it][91];
        Ccondpregg3_6[it]=MYARRAY[it][92];
        Ccondpregg3_7[it]=MYARRAY[it][93];
        Ccondpregg3_8[it]=MYARRAY[it][94];
        Ccondpregg3_9[it]=MYARRAY[it][95];
        Ccondpregg4a_1[it]=MYARRAY[it][96];
        Ccondpregg4a_2[it]=MYARRAY[it][97];
        Ccondpregg4a_3[it]=MYARRAY[it][98];
        Ccondpregg4a_4[it]=MYARRAY[it][99];
        Ccondpregg4a_5[it]=MYARRAY[it][100];
        Ccondpregg4a_6[it]=MYARRAY[it][101];
        Ccondpregg4a_7[it]=MYARRAY[it][102];
        Ccondpregg7b[it]=MYARRAY[it][103];
        Ccondpregg8b[it]=MYARRAY[it][104];
        Ccondpregg9[it]=MYARRAY[it][105];
        Ccondpregg11b[it]=MYARRAY[it][106];
        Ccondpregpreterm[it]=MYARRAY[it][107];
        Ccondpregg24[it]=MYARRAY[it][108];
        Ccondpregg23[it]=MYARRAY[it][109];
        Cwais_pb_num[it]=MYARRAY[it][110];
        Cwais_pb_vo[it]=MYARRAY[it][111];
        Cbfi_pb_ama[it]=MYARRAY[it][112];
        Cbfi_pb_ape[it]=MYARRAY[it][113];
        Cbfi_pb_ext[it]=MYARRAY[it][114];
        Cbfi_pb_neu[it]=MYARRAY[it][115];
        Cbfi_pb_res[it]=MYARRAY[it][116];
        Cpsi_pb_total[it]=MYARRAY[it][117];
        Hbargg2a[it]=MYARRAY[it][118];
        Hbargg2b[it]=MYARRAY[it][119];
        Hbargg2c[it]=MYARRAY[it][120];
        Hbargg2d[it]=MYARRAY[it][121];
        Hbargg2e[it]=MYARRAY[it][122];
        Hbargg2f[it]=MYARRAY[it][123];
        Hbargg2g[it]=MYARRAY[it][124];
        Hbargg2h[it]=MYARRAY[it][125];
        Hbargg2i[it]=MYARRAY[it][126];
        Hbargg2j[it]=MYARRAY[it][127];
        Hwomdecide12[it]=MYARRAY[it][128];
        Hfathdecides12[it]=MYARRAY[it][129];
        Hbothdecide[it]=MYARRAY[it][130];
        Hcaresacwom[it]=MYARRAY[it][131];
        Hcaresacman[it]=MYARRAY[it][132];
        Hopofman121[it]=MYARRAY[it][133];
        Hopofman122[it]=MYARRAY[it][134];
        Hopofman123[it]=MYARRAY[it][135];
        Hopofman124[it]=MYARRAY[it][136];
        Cinvh3[it]=MYARRAY[it][137];
        Cinvh4[it]=MYARRAY[it][138];
        Cinvh5[it]=MYARRAY[it][139];
        Cinvh6[it]=MYARRAY[it][140];
        Cinvh7[it]=MYARRAY[it][141];
        Cinvh8[it]=MYARRAY[it][142];
        Cinvh9[it]=MYARRAY[it][143];
        Cinvh13[it]=MYARRAY[it][144];
        Cinvf11a[it]=MYARRAY[it][145];
        Cinvf11b[it]=MYARRAY[it][146];
        Cinvf11c[it]=MYARRAY[it][147];
        Cinvf11d[it]=MYARRAY[it][148];
        Cinvf11e[it]=MYARRAY[it][149];
        Cinvf11f[it]=MYARRAY[it][150];
        Cinvf11g[it]=MYARRAY[it][151];
        Cinvf11h[it]=MYARRAY[it][152];
        Cinvf11i[it]=MYARRAY[it][153];
        Cinvf11j[it]=MYARRAY[it][154];
        Cinvf11k[it]=MYARRAY[it][155];
        Cinvhm2_10[it]=MYARRAY[it][156];
        Cinvhm2_11[it]=MYARRAY[it][157];
        Cinvhm2_12[it]=MYARRAY[it][158];
        Cinvhm2_13[it]=MYARRAY[it][159];
        Cinvhm2_14[it]=MYARRAY[it][160];
        Cinvhm2_15[it]=MYARRAY[it][161];
        Cinvhm2_16[it]=MYARRAY[it][162];
        Cinvhm2_18[it]=MYARRAY[it][163];
        Csharesbedroomhowmany12[it]=MYARRAY[it][164];
        Csharesbedhowmany12[it]=MYARRAY[it][165];
        g42_a2[it]=MYARRAY[it][166];
        g42_b2[it]=MYARRAY[it][167];
        g42_c2[it]=MYARRAY[it][168];
        g42_d2[it]=MYARRAY[it][169];
        g42_e2[it]=MYARRAY[it][170];
        g42_f2[it]=MYARRAY[it][171];
        g42_a1[it]=MYARRAY[it][172];
        g42_b1[it]=MYARRAY[it][173];
        g42_c1[it]=MYARRAY[it][174];
        g42_d1[it]=MYARRAY[it][175];
        g42_e1[it]=MYARRAY[it][176];
        g42_f1[it]=MYARRAY[it][177];
        f21a_p_t[it]=MYARRAY[it][178];
        f21b_p_t[it]=MYARRAY[it][179];
        f21c_p_t[it]=MYARRAY[it][180];
        f21d_p_t[it]=MYARRAY[it][181];
        f21e_p_t[it]=MYARRAY[it][182];
        f21f_p_t[it]=MYARRAY[it][183];
        f21g_p_t[it]=MYARRAY[it][184];
        f21h_p_t[it]=MYARRAY[it][185];
        f21i_p_t[it]=MYARRAY[it][186];
        f21j_p_t[it]=MYARRAY[it][187];
        f21k_p_t[it]=MYARRAY[it][188];
        f21l_p_t[it]=MYARRAY[it][189];
        f21m_p_t[it]=MYARRAY[it][190];
        f21n_p_t[it]=MYARRAY[it][191];
        
        f21a_m_t[it]=MYARRAY[it][192];
        f21b_m_t[it]=MYARRAY[it][193];
        f21c_m_t[it]=MYARRAY[it][194];
        f21d_m_t[it]=MYARRAY[it][195];
        f21e_m_t[it]=MYARRAY[it][196];
        f21f_m_t[it]=MYARRAY[it][197];
        f21g_m_t[it]=MYARRAY[it][198];
        f21h_m_t[it]=MYARRAY[it][199];
        f21i_m_t[it]=MYARRAY[it][200];
        f21j_m_t[it]=MYARRAY[it][201];
        f21k_m_t[it]=MYARRAY[it][202];
        f21l_m_t[it]=MYARRAY[it][203];
        f21m_m_t[it]=MYARRAY[it][204];
        f21n_m_t[it]=MYARRAY[it][205];
        HDens1[it]=MYARRAY[it][206];
        HDens5[it]=MYARRAY[it][207];
        HDens10[it]=MYARRAY[it][208];
        
    }//Finishmloading data
    
    //LOADING THE MEASUREMENT SYSTEM OF THE UNOBSERVED FACTORS
    
    
    
    //===============================================
    //Block to get the initial parameters to optimize
    //===============================================
    
    long double test2=F_likelihood(PAR, Cliveswithfather12, Cliveswithmother12, Hchores12, Meffort12, Feffort12, Mage12, Fage12, Cedad_meses12, Ctestsfactor4_2012, Ctestsfactor1_10, Myrschool12, Fyrschool12, Ffraclabor12, Mfraclabor12, Mwage12, Fwage12, Mnlincome12, Fnlincome12, CfactorInv12, Cchildcare12, Ccareskills, Hbarg4, Feffort10, Meffort10, CfactorInv10, Cfactorbirth, Mwage10, Fwage10, Mnlincome10, Fnlincome10, Cedad_meses10, Cliveswithfather10, Cliveswithmother10, Mage10, Fage10, Mfraclabor10, Ffraclabor10, Cchildcare10, Hchildcareobs, Ccareskills, MTJH,FMRATIO,Unemployment, Wageratio,Distance,Magegroup10,Magegroup12,Hmemberstotal10,Hmemberstotal12,CSTDtepsi_pb_coo10,CSTDtepsi_pb_len10,CSTDtepsi_pb_mot10,CSTDtvip_pb10,CSTDcbcl1_pb_110,CSTDcbcl1_pb_210,CSTDcbcl1_pb_310,CSTDcbcl1_pb_410,CSTDcbcl1_pb_510,CSTDcbcl1_pb_610,CSTDcbcl1_pb_710,CSTDtadi_pb_cog12, CSTDtadi_pb_mot12, CSTDtadi_pb_len12, CSTDtadi_pb_se12, CSTDbt_112, CSTDbt_212, CSTDbt_312, CSTDbt_412, CSTDbt_512, CSTDbt_t12,  CSTDhtks_st12, CSTDbdst_st12, CSTDppvt_t12,Ccondpregg3_1, Ccondpregg3_2 ,Ccondpregg3_3, Ccondpregg3_4, Ccondpregg3_5, Ccondpregg3_6, Ccondpregg3_7, Ccondpregg3_8, Ccondpregg3_9, Ccondpregg4a_1, Ccondpregg4a_2, Ccondpregg4a_3, Ccondpregg4a_4, Ccondpregg4a_5, Ccondpregg4a_6, Ccondpregg4a_7, Ccondpregg7b, Ccondpregg8b, Ccondpregg9, Ccondpregg11b, Ccondpregpreterm, Ccondpregg24, Ccondpregg23,Cwais_pb_num, Cwais_pb_vo, Cbfi_pb_ama,  Cbfi_pb_ape, Cbfi_pb_ext, Cbfi_pb_neu, Cbfi_pb_res, Cpsi_pb_total,Hbargg2a, Hbargg2b, Hbargg2c, Hbargg2d, Hbargg2e, Hbargg2f, Hbargg2g, Hbargg2h, Hbargg2i, Hbargg2j, Hwomdecide12, Hfathdecides12, Hbothdecide, Hcaresacwom, Hcaresacman,  Hopofman121,  Hopofman122,  Hopofman123,  Hopofman124, Cinvh3, Cinvh4, Cinvh5, Cinvh6, Cinvh7, Cinvh8, Cinvh9, Cinvh13,Cinvf11a, Cinvf11b,  Cinvf11c, Cinvf11d, Cinvf11e, Cinvf11f, Cinvf11g, Cinvf11h, Cinvf11i, Cinvf11j, Cinvf11k, Cinvhm2_10, Cinvhm2_11, Cinvhm2_12, Cinvhm2_13, Cinvhm2_14, Cinvhm2_15, Cinvhm2_16, Cinvhm2_18, Csharesbedroomhowmany12,  Csharesbedhowmany12,g42_a2, g42_b2, g42_c2, g42_d2, g42_e2, g42_f2,g42_a1, g42_b1, g42_c1, g42_d1, g42_e1, g42_f1,f21a_p_t, f21b_p_t, f21c_p_t, f21d_p_t, f21e_p_t, f21f_p_t, f21g_p_t, f21h_p_t, f21i_p_t, f21j_p_t, f21k_p_t, f21l_p_t, f21m_p_t, f21n_p_t,f21a_m_t, f21b_m_t, f21c_m_t, f21d_m_t, f21e_m_t, f21f_m_t, f21g_m_t, f21h_m_t, f21i_m_t, f21j_m_t, f21k_m_t, f21l_m_t, f21m_m_t, f21n_m_t,HDens5,950);
    
    return(test2);
    
}

//Defining the likelihood function that is going to be minimized finally
int iterat=0;
double FLIKMINIMIZED(unsigned n, const double *x, double *grad, void *FLIKMNIMZED_data){
    double result=F_likelihood_FIN(x);
    cout << result << " evalF"<< endl;
    iterat=iterat+1;
    printf("Iteration=(%d); Feval=%0.10g\n", iterat, result);
    return(result);
}



//FUNCTIONS TO READ AND WRITE HOPSPACK OPTIMIZER
// --------------------------------------------------------------------
//  Internal Function readInputFile
// --------------------------------------------------------------------
static bool  readInputFile (const string &            szExeName,
                            const string &            szFileName,
                            vector< double > &  x)
{
    const int  REQBUFLEN = 10;
    char       reqBuf[REQBUFLEN + 1];
    
    ifstream  fp;
    
    //---- OPEN THE INPUT FILE.
    fp.open (szFileName.c_str());
    if (!fp)
    {
        cerr << szExeName << " - ERROR opening input file '"
        << szFileName << "'." << endl;
        return( false );
    }
    
    //---- READ AND VERIFY THE INPUT REQUEST TYPE.
    fp >> reqBuf;
    for (int  i = 0; i < REQBUFLEN + 1; i++)
    {
        if ((reqBuf[i] == '\n') || (reqBuf[i] == '\r'))
        {
            reqBuf[i] = 0;
            break;
        }
    }
    if ((reqBuf[0] != 'F') || (reqBuf[1] != 0))
    {
        cerr << szExeName << " - ERROR reading request type from '"
        << szFileName << "'." << endl;
        cerr << "  Read '" << reqBuf << "'" << endl;
        fp.close();
        return( false );
    }
    
    //---- READ THE LENGTH OF x.
    int  n;
    fp >> n;
    if (n != 325)
    {
        cerr << szExeName << " - ERROR reading n from '"
        << szFileName << "'." << endl;
        cerr << "  Read " << n << ", but problem size should be 2" << endl;
        fp.close();
        return( false );
    }
    
    x.resize (n);
    
    //---- READ x.
    for (int  i = 0; i < n; i++)
        fp >> x[i];
    
    fp.close();
    
    //---- RETURN SUCCESS.
    return( true );
}


// --------------------------------------------------------------------
//  Internal Function writeOutputFile
// --------------------------------------------------------------------
static bool  writeOutputFile (const string &  szExeName,
                              const string &  szFileName,
                              const double    f)
{
    ofstream  fp;
    
    fp.open (szFileName.c_str(), ios::trunc);
    if (!fp)
    {
        cerr << szExeName << " - ERROR opening output file '"
        << szFileName << "'." << endl;
        return( false );
    }
    
    //-- WRITE ALL DIGITS USING SCIENTIFIC NOTATION.
    fp.setf (ios::scientific);
    fp.precision (15);
    
    //---- WRITE THE NUMBER OF OBJECTIVES AND THEIR VALUES TO THE OUTPUT.
    fp << "1" << endl;
    fp << f << endl;
    
    fp.close();
    
    return( true );
}


//Likelihood function used for hopspack. Input is vector rather than double:
long double F_likelihood_FINvecDOUBLE(const vector<double>PAR){
    double PARDOUBLE[325]={};
    for(int ii=0; ii<325; ii=ii+1){
        PARDOUBLE[ii]=PAR[ii];
    }
    double test=F_likelihood_FIN(PARDOUBLE);
    return(test);
}

int  main (const int           argc,
           const char * const  argv[])

{
    //If want to run HOPSPACK SET CONDITION TO ONE (1)
    int HOPSPACK=0; //If want to run hopspack run 1
    int optimize=0; //If want to do nonlinearoptimization part put 1.
    if (HOPSPACK!=1){
        
        //0. Define the current directory
        //chdir("/Users/rodrigoazuero/Dropbox/BACKUPRODRIGO/Research/Chile/RR/BEHAVIORAL70");
        
        //===============================================
        //Block to get the initial parameters to optimize
        //===============================================
        
        //I will generate two vectors. One of type vector and
        //another one double because it seems that nlopt only allows
        //me to use double types objects for the minimization.
        
        //PAROPTFOUND.csv has length of 325 because in previous file (49) I put it to save up until
        //the number 325. I load the 325 parameters in PAR but PARDOUBLE should contain only 325.
        
        std::string line;
        std::vector<std::vector<std::string> > values;
        vector<double> PAR;
        PAR.resize(325);
        
        
        
        
        std::ifstream PARAMETERS ("PAR32.csv");
        
        std::string linePARAM;
        int itpar=0;
        while(std::getline(PARAMETERS, line))
        {
            PAR[itpar]=::atof (line.c_str());
            //PAR[itpar]=std::stod (line); only C++11s
            
            
            
            itpar=itpar+1;
        }
        //And the ones that were not stored
        
        //Loading The parameters into PARDOUBLE
        double PARDOUBLE[326]={};
        for (int ii=0;ii<325;ii=ii+1){
            PARDOUBLE[ii]=PAR[324-ii];
        }
        //PARDOUBLE[60]=0;
        
        cout << " before evaluating likelihood  " << endl;
        
        cout << F_likelihood_FIN(PARDOUBLE) << " LIKELIHOOD " << endl;
        //cout << F_SMOOTHING_FIN(PARDOUBLE) << "smoothing " << endl;
        
        std::cout << "Hello, World!\n";
        //Now doing for the likelihood function
        cout <<"---------------"<< endl;
        
        
        
        //Setoptimize=1 if I want to perform the optimization
        
        if(optimize==1){
            
            nlopt_opt opt3;
            opt3 = nlopt_create(NLOPT_LN_NELDERMEAD, 325); /* algorithm and dimensionality */
            nlopt_set_min_objective(opt3, FLIKMINIMIZED, NULL);
            nlopt_set_xtol_rel(opt3, 1.0);
            nlopt_set_ftol_abs(opt3,1.0);
	    //nlopt_set_maxtime(opt3,  10000);
            //nlopt_set_maxeval(opt3,100);
            nlopt_set_maxeval(opt3,10000);
            double x3[325]={0};
            for(int it=0;it<325;it=it+1){
                x3[it]=PARDOUBLE[it];
            }
            //PARDOUBLE[60]=;
            double minf3; /* the minimum objective value, upon return */
            nlopt_optimize(opt3, x3, &minf3);
            // cout << minf3 << " minf3 " << endl;
            //cout << nlopt_optimize(opt3, x3, &minf3) << " code " << endl;
            if (nlopt_optimize(opt3, x3, &minf3) < 0) {
                printf("nlopt failed!\n");
            }
            else {
                printf("found minimum at f(%g) = %0.10g\n", x3[0], minf3);
            }
            
            //Saving the vector of optimal parameters in cv file
            ofstream optparam("PAR33.csv");
            for (int it=0;it<325;it=it+1){
                optparam<<x3[324-it] << endl;
            }
            optparam.close();
            
            
            nlopt_destroy(opt3);
        }
    }//END IF HOPSPACK NOT RUN.
    
    
    //Run hopspack if decided
    if (HOPSPACK==1){
        vector< double >  x;
        double            f;
        if (argc != 5)
        {
            cerr << "usage: " << argv[0]
            << " <input file> <output file> <tag> <type>" << endl;
            return( -1 );
        }
        
        if (readInputFile (argv[0], argv[1], x) == false)
            return( -2 );
        
        f = F_likelihood_FINvecDOUBLE (x);
        
        if (writeOutputFile (argv[0], argv[2], f) == false)
            return( -3 );
        
    }
    
    return 0;
    
}
