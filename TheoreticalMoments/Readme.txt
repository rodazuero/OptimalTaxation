

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


