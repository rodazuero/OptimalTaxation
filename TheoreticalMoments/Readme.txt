

Version 1: Version without utility penalty on providing labor supply. 

Version 2: In this version of the code, we run the new version of the model with lf+li not necessarily bounded by one. Utility penalty from
providing both types of labor and extra consumption utility from informal labor income. 


Version 3: In this version we start to make more complex tax functions. Payroll taxes no longer as they are but aiming to generate the increasing pattern. I also start guess of the equilibrium in values nearby. 


Version 4: We include the fact that corporate tax rates can depend on other elements such as total income and number of workers. Distribution of skills to be log-normal. 


Lower Triangular                     Transpose
     2	     0	     0		     2	     6	    -8	
     6	     1	     0		     0	     1	     5	
    -8	     5	     3		     0	     0	     3	


Version 5: Start to make the code and the calibrated version look similar. The first step is to set up corporate income taxes. For this, note that the production levels are around 0.4 to 0.9. Theoretical distribution in the equilibrium:

In the empirical distribution, we have: (in USD)

10% 913.5
20% 1682.1
30% 2583
40% 3780
50% 5355.315
60% 7560
70% 11088
80% 18614.988
90% 48762




In the theoretical distribution we have approximately
10% 0.4 
50% 0.65
90% 0.8

If we want 0.65 to be 5,355-> Divide by 8000 for the taxes.
If we want 0.8 to be  18614-> Divide by 60952.5
I will divide by 20,000 and then play with mean of log-normal distribution to see. 

(Recall all is annual)

That is,: taxes are given in the following table monthly:
5000-20
8000-50
13000-200
20000-400
30000-600
30% if annual income exceeds 30,000
36% if annual income exceeds 30,000 and labor force is in excess of 20. 
In the empirical distribution, 20 employees correspond to percentile 98 of demand, which is about 0.33 in the theoretical distribution. 



With this in mind, the theoretical tax function should look something like this:


Limit, tax payed

5000*0.31*12/8000-20*12/8000
8000*12/8000-50*12/8000
13000*12/8000-200*12/8000
20000*12/8000-400*12/8000
30000*12/8000-600*12/8000
>30000*12/8000- 0.3 of profits
>30000*12/8000 and more than 0.33 employees, pays 0.35 of taxes. 


5000*12/60000-20*12/8000
8000*12/8000-50*12/8000
13000*12/8000-200*12/8000
20000*12/8000-400*12/8000
30000*12/8000-600*12/8000
>30000*12/8000- 0.3 of profits
>30000*12/8000 and more than 0.33 employees, pays 0.35 of taxes. 





Transforming it to annual *12, and to fir the theoretical situations

First approach, transforming everything to have similar scale to the actual one. 
Modified variance-covariance matrix and mean of tthetas. 

Decidí expandir los rangos de tthetas. Problema con valores iniciales de ni, nf. Los puse en 0.1,0.1. Doing this solved the problem we had in the workers. 

We run into a simliar problem for the demand of informal workers. For this reason, I decreased the expected demand from them, initial point to be 2.3
Evasion same thing, we do it at around 240. 


#Version 6. We modify this version to be able to use as input the variance-covariance matrix of the log-normal distribution as well as its mean. 

Function ExcessDemandFunctions needs additional inputs:
mmu1>- mean of log-tthetaw
mmu2<- mean of log-tthetae
ssigma1-> sd of log-tthetaw
ssigma2-> sd of log-tthetae
rho12-> cov of log-tthetaw, log-tthetae


MAryan found a typo: needed to multiply by ttheta the labor supply of individuals. Fixed this. This will be done in CodeV7.

First attempt to solve for the eq: run the same code fixing the typo. Does't work due to the amount of demand (little informal too formal). Problem was solved.  Trying to make more workers, reduce mmu2 from 4 to 3. Decision to do the production in hundreds of dollar a year. With this in mind, the corresponding taxes should be:



Code v8. Generating the tax breaks of RUS etc. 
(Payments)




20 (6.3*12) -> 75.6 usd year -> 0.756 hundreds of dollar
50 (15.75*12)->189 USD year -> 1.89 hundredths of dollars
200 (756 usd a year) -> 7.56 hundreds of dollars
400 (1512 a year) 15.12 hundreds of dollars
600 (2268 a year) 22.68 hundreds of dollars

Categories of income are:
5,000 (18,900 year) 189 hundreds of dollars
8,000 (30,240 year) 302.4 hundreds of dollars
13,000 (49,140 a year) 491.4 hundreds of dollars
20,000 (75,600 a year) 756 hundreds of dollars
30,000 (113,400 a year) 1,13.4 hundreds of dollars

#Code v9. modification of new equi. V8 included the first version with specific 
details of tax system. Now we will incorporate pit. For workers, I assume it is in hundreds.

Tax brackets in soles:


0   ->[0,24150]      -> [0,7607]     -> Divide by 10 
15  ->[24150,117300] -> [7607,36949] ->
21  ->[117300,210450]-> [36949,66291]
30  ->[210450,+]     -> [66291,+]


#Typo found, alpha should be around 0.8 try to fix this in the code version10. What I did was return usual tax funciton (0.3*Profits) to re-calibrate. 
Config works:



double aalpha=0.5;
    double tthetae=8.600943;
    double tthetaw=5.813159;
    double wi=0.9;
    double wf=1.68;
    double ni=2.3;
    double nf=0.53*70;
    double ggamma=0.7;
    double ddelta=0.7;
    double bbeta=0.2;
    double ssigma=0.2;
    double kkappa=0.02;
    double psi=0.4;
    double chi=15;
    double rrho=0.1;
    double lf=2.1;
    double li=2.1;
    double z=24;
    double mmu1=1.3;
    double mmu2=1.5;
    double ssigma1=1.5;
    double ssigma2=1.3;
    double rho12=0.3;


    #Need to increase tax evasion parameter. 


#Code v14. going back to see which code actually works. 


#Trying to modify the tax function to have a smooth thing that can generate something similar to the one in the legislation. 


Works in C++ but not in R?? 
long double PIT(double tthetaw, double wf, double lf){
    double taxburden=tthetaw*wf*lf;
    double taxrate=pow(taxburden,0.015)-1;
    double ans=taxrate*taxburden;
    return(ans);
}


#Modifying bit the wages does not change that much the overall outcomes. 


#10-01-2018. Need to:

1. Increase proportion of people who are entrepreneurs. Theoretical 3% real 36%. 
2. INcrease inequality in labor income. 

3. Increase Informality of labor supply

#Solutions:
Increase mmu2 to increase the proportion of people who are entrepreneurs. Setting it initially to 2 just lead to 3%. 

Trying to decrease the distribution of mmu so that the demands do not blow out of proportion. This will imply that we have to change the tax settings as well. 

The parameters were effectively decreased to: 


double aalpha=0.8;
    double tthetae=12.84209;
    double tthetaw=3.873585;
    double wi=1.38;
    double wf=1.46;

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
    double mmu1=0.4;
    double mmu2=0.7;
    double ssigma1=0.2;
    double ssigma2=0.3;
    double rho12=0.1;


This works nicely. Lesson learned: if we increase the distribution of ttheta, demands and supplies explode!

This combination works better!!!

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

Problem keeps being the small number of entrepreneurs. Probably increasing chi or other stuff might help


//Running the optimizer of distance estimator in the C++ file. It is necessary to readjust it so that it can compile in R. The difference is due to the function called to sort vectors. 

//It is necessary to compute theoretically the aggregate sum of taxes for each decile and compare that with the corresponding empirical part. 

Modified increasing 100+nf+ni to increase entrepreneurs!


//Tesla: remove everything from RCPP. 
Necesario poner #include <armadillo> y tambi;en en las librerias a buscar incluir el directorio -I/usr/lib64/R/library/RcppArmadillo/include/ que va a encontrar el archivo armadillo.h. 
También, #include <iterator>.
Change std::begin and std::end to nothing! std::begin(algo)-> algo




g++ mainTesla.cpp -std=c++0x -I/usr/include/c++/3.4.6/backward/  -I/usr/lib64/R/library/BH/include/  -I/home/razu/install/include/ -I/usr/lib64/R/library/RcppArmadillo/include/ -I/home/razu/boost_1_61_0/ -L/home/razu/install/lib/ -lnlopt   -fopenmp -Wall   -lgsl -lgslcblas  -o main

g++ testTesla.cpp -std=c++0x -I/usr/include/c++/3.4.6/backward/  -I/usr/lib64/R/library/BH/include/  -I/home/razu/install/include/ -I/usr/lib64/R/library/RcppArmadillo/include/ -I/home/razu/boost_1_61_0/ -L/home/razu/install/lib/ -lnlopt   -fopenmp -Wall   -lgsl -lgslcblas  -o main

*Code to run in parallel:

g++ mainTesla.cpp  -I/usr/include/c++/3.4.6/backward/  -I/usr/lib64/R/library/BH/include/  -I/home/razu/install/include/ -I/usr/lib64/R/library/RcppArmadillo/include/ -I/home/razu/boost_1_61_0/ -L/home/razu/install/lib/ -lnlopt   -fopenmp -Wall   -lgsl -lgslcblas  -o main

g++ mainTesla.cpp -std=c++0x -I/usr/include/c++/3.4.6/backward/  -I/usr/lib64/R/library/BH/include/  -I/home/razu/install/include/ -I/usr/lib64/R/library/RcppArmadillo/include/ -I/home/razu/boost_1_61_0/ -L/home/razu/install/lib/ -lnlopt   -fopenmp -Wall   -lgsl -lgslcblas  -o main



g++ mainTesla.cpp -std=gnu++0x -I/usr/include/c++/3.4.6/backward/  -I/usr/lib64/R/library/BH/include/  -I/home/razu/install/include/ -I/usr/lib64/R/library/RcppArmadillo/include/ -I/home/razu/boost_1_61_0/ -L/home/razu/install/lib/ -lnlopt   -fopenmp -Wall   -lgsl -lgslcblas  -o main


g++ mainTesla.cpp  -I/usr/include/c++/3.4.6/backward/  -I/usr/lib64/R/library/BH/include/  -I/home/razu/install/include/  -I/home/razu/boost_1_61_0/ -L/home/razu/install/lib/ -lnlopt   -fopenmp -Wall   -lgsl -lgslcblas  -o paralleloutput

#Remember to change the size of MYVAR and the limits of the pragma omp for loop when the size of the sobol calibration is changed

#Initial labor supply answers increased to 24.1. 
#I will also increase the initial wages guess because it is not throwing much equilibria from around 8 to something higher////

#When expanding the range of parameters very few equilibria are found with the sobol. need to be careful with this. 

To not see the output: nohup ./main > nohupNLOPT.out&





#Sobol description. 
#The first attempt was to generate a quite wide range of parameters. 

#Problem with lack of equilibria found iseems to be related to initial conditions found. 

Initially, we had:
li=24
lf=24
ni=2.3
nf=0.53*8
z=24. 

New first initial conditions:
li=0.5
lf=0.5
ni=10
nf=50
z=24.

*Once we run the sobol, we find that there is a problem still with the number of entrepreneurs and workers. We need to work on it as the closest numbers are always very far from the situation optimal. Tests for which can be actually close to the real number of entrepreneurs. The ones with high levels of entrepreneurs:


Analysis of other variables and the Proportion of informality. Some bounds on the moments that we obtain. 

The problem is that ddelta goes from 0 to 600 and for a large number of values (up until 50-600) there is no informality. For this reason, it is probably a good idea to re-run everything bounding ddelta to be under 50. 

1. Need to have an upper bound on ddelta at around 50. Values larger than that will have zero informality. Before that, it was up until 600. 

#The best combination found was:
aalpha=0.8;
    ggamma=0.198503;
    ddelta=0.180862;
    bbeta=0.247681;
    ssigma=0.468333;
    kkappa=0.0689548;
    psi=2.19302;
    chi=1.09956;
    rrho=2.74542;
    mmu1=0.240092;
    mmu2=1.50935;
    ssigma1=0.113429;
    ssigma2=1.16547;
    rho12=0.10219;

#Is there anything we can do to increase the number of entrepreneurs with these parameters with slight modifications? I will try the following:

1. lowering aalpha from 0.8 to 0.7. doesn't work. 
to 0.5-> it increases but only to 2% from 0.7%. 
to 0.3->  no longer works. 

2. increasing mmu2:
from 1.50935 to 3-> Doesn't work. 
From 1.50935 to 2-> Doesn't work

3. Increasing mmu1 
0.240092 to 0.5-> Doesn't work
0.240092 to 0.75-> Doesn't work
0.240092 to 1-> Doesn't work
0.240092 to 1.5-> Doesn't work
0.240092 to 5.1-> Doesn't work

4. Lowering ssigma2:
1.16547 to 0.7-> Doesn't work
1.16547 to 0.5-> ligeramente a 1%. 
1.16547 to 0.2-> ligeramente a 1%. 

5. Lowering cchi:
1.09956 to 0.7-> Doesn't work. 
1.09956 to 0.5-> Doesn't work. 
1.09956 to 0.2-> Doesn't work. 
1.09956 to 2.0-> Doesn't work. 

6. Increasing psi 
2.1930 to 10.19302->doesn't work
2.1930 to 20.19302->doesn't work
2.1930 to 1.19302->doesn't work


#Now I will try to give for free some work for the entrepreneurs such that
we count them as self-employed and they don't need to hire someone or else. 

Doing
double profm=tthetae*pow((10+ni+nf),aalpha)-wi*ni-wf*nf-TnActual(nf*wf);

It does increase the number of entrepreneurs but breaks all actually. 

Parece que sí va a tener que ser por este lado. 

50-> leads to entrepreneurs of approximately 14%. 

100-> 29%. 


#Increasing to have parameter “c” in the production function so that 
#we have more entrepreneurs implies:

#1. Having more parameters in all production functions. 
#2. Final Profits should have parameter paramvec of size 9. Following modifications:
	a. FinProfits: c=paramvec[8];
	b. profitsFinMaxim: Params[8]=c;
	c. iDecision: c=Params[11];
	d. ExcessDemandFunctions. Increase one variable in parameters:
		ParamsDecision[11]=Params[9];. All variables from 9 to 13 should then
		be increased by one. 
		double mmu1=Params[10];
    		double mmu2=Params[11];
    		double ssigma1=Params[12];
    		double ssigma2=Params[13];
    		double rho12=Params[14];

	e. StandardizedExcessDemands: Vector Others[20] increased to Others[21]
	f. EqWages DOthers[20] increased to [21]
	g. EqWagesNumericVector [19] to [20] in DOthers. 
	h. TheoMoments: 19 to 20. Others
			ParamsDecision expand to 12 and 11->9. 
			In DOthers squeeze to allow c=DOthers[9];
	i. StandardizedDistanceEstimator-> c=parameters[13]. Others[9] is C 
		and displace the rest. 
	j. MinimizingDistance-> c is PArameters[13];
	k. SobolRun. VecOthers to 21 and VecOthers increased one (15)
	

#Passing parameters:
FinProfits<-
double wi=paramvec[0];
    double wf=paramvec[1];
    double aalpha=paramvec[2];
    double ddelta=paramvec[3];
    double ggamma=paramvec[4];
    double bbeta=paramvec[5];
    double ssigma=paramvec[6];
    double tthetae=paramvec[7];
    double c=paramvec[8];//ADDED


profitsFinMaxim
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

iDecision
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
    double c=Params[11];//ADDED


Exces
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


StandardizedExcessDemands
Params.resize(16);
    for(int it=0;it<16;it++){
        Params[it]=Others[it];
    }


EqWages
double DOthers[21];
    for(int it=0; it<21;it++){
        DOthers[it]=Others[it];
    }

EqWagesNumericVector
double DOthers[20];
    for(int it=0; it<20;it++){
        DOthers[it]=Others[it];
        cout << DOthers[it]<< " DOthers[it] in eqwagesnumeric" << endl;
    }

TheoMoments
double DOthers[20];
    for(int it=0; it<20;it++){
        DOthers[it]=Others[it];
    }
ParamsDecision.resize(12);
ParamsDecision[11]=DOthers[14];
double c=DOthers[14];


#ANALYZING step by step:
SobolRun(Wages,InitLWorkersDecision,InitProfDecision);

SobolRun reads SobolDim.csv and loads it into VectorOthers of size 20. 
VectorOthers is passed on to distanceEstimator. 

VectorOthers[14] is inserted to be c. It is the 15th element in 14 step. 
Then, 15, 16, are initial conditions. 


DISTANCEESTIMATOR->OTHERS OF SIZE 20 PASSED TO EQUILIBRIUMMOMENTS. Exactly as the 
VectorOthers input in SobolRun passed to equilibrium moments. 

EQUILIBRIUMMOMENTS->OTHERS OF SIZE 20 PASSED TO EqWagesNumericVector

Passed as it is to ExcessDemandTotal
Passed to StandardizedExcessDemands. Here is the weird thing. It is passed to ExcessDemandFunctions


#There is a problem with delta and gamma passing them. 

#In tesla: in addition to modifying everything as in the normal one, it use necessary to change DistanceNonVectorized. 





