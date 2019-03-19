//Sobol run structure



1. SobolRun(arma::vec WagesVectorIn,
              arma::vec InitLWorkersDecision,
              arma::vec InitProfDecision){


	//Call 
	DistanceNonVectorized( var1,var2,var3,var4,var5,var6,var7,var8,var9,var10,var11,var12,var13,var14,var15,WagesVectorIn,InitLWorkersDecision,InitProfDecision);
}

2. double DistanceNonVectorized(double vec1, double vec2, double vec3, double vec4,double vec5, double vec6, double vec7, double vec8,double vec9, double vec10, double vec11, double vec12,double vec13, double vec14, double vec15, arma::vec WagesVectorIn,
                             arma::vec InitLWorkersDecision,
                             arma::vec InitProfDecision){

vector<double> VectorOthers(20);
VectorOthers[0]=vec1;
    VectorOthers[1]=vec2;
    VectorOthers[2]=vec3;
    VectorOthers[3]=vec4;
    VectorOthers[4]=vec5;
    VectorOthers[5]=vec6;
    VectorOthers[6]=vec7;
    VectorOthers[7]=vec8;
    VectorOthers[8]=vec9;
    VectorOthers[9]=vec10;
    VectorOthers[10]=vec11;
    VectorOthers[11]=vec12;
    VectorOthers[12]=vec13;
    VectorOthers[13]=vec14;
    VectorOthers[14]=vec15;
    
    
    VectorOthers[15]=InitLWorkersDecision[0];
    VectorOthers[16]=InitLWorkersDecision[1];
    VectorOthers[17]=InitProfDecision[0];
    VectorOthers[18]=InitProfDecision[1];
    VectorOthers[19]=InitProfDecision[2];
return(DistanceEstimator(VectorOthers, WagesVectorIn,InitLWorkersDecision,InitProfDecision));
}



3. double DistanceEstimator(arma::vec Others, arma::vec WagesInit,
                         arma::vec armaInitLWorkers,arma::vec armaInitProf){



vector<vector<double> >Theomoments(10);
    for(int it=0;it<10;it++){
        Theomoments[it].resize(10);
    }
    cout << " finding theoretical moments"<< endl;
    Theomoments=EquilibriumMoments(Others, WagesInit,
                                   armaInitLWorkers,armaInitProf);
}


4. vector<vector<double> > EquilibriumMoments(arma::vec Others, arma::vec WagesInit,
                                           arma::vec armaInitLWorkers,arma::vec armaInitProf){



 }



 5. vector<vector<double> > TheoMoments(arma::vec Others, arma::vec WagesEquilibrium,
                                    arma::vec armaInitLWorkers,arma::vec armaInitProf){


 	double DOthers[21];
    for(int it=0; it<20;it++){
        DOthers[it]=Others[it];
    }


 }




 6. vector<vector<double> > EquilibriumMoments(arma::vec Others, arma::vec WagesInit,
                                           arma::vec armaInitLWorkers,arma::vec armaInitProf){
    
    cout << " finding equilibrium moments"<< endl;
    //1. Find the equilibrium wages
    //Testing the moment generating function
    arma::vec WagesEquilibrium=EqWagesNumericVector(Others, WagesInit);
    cout << WagesEquilibrium[0]<< " wi "<< endl;
    cout << WagesEquilibrium[1]<< " wf "<< endl;
    //2. Once we find the equilibrium wages, obtain the moments.
    vector<vector<double> >answ(10);
    for(int it=0;it<9;it++){
        answ[it].resize(10);
    }
    cout << " just about to run theomoments"<< endl;
    answ=TheoMoments(Others, WagesEquilibrium,
                     armaInitLWorkers,armaInitProf);
    cout << " finished running theomoments"<< endl;
    return(answ);
}


7. arma::vec EqWagesNumericVector(arma::vec Others, arma::vec WagesInit){
    
    
    //Others: Vector of parameters used to find equilibrium.
    //        Need to reconvert to double.
    cout << Others[15] << " Others[15] in "<< endl;
    double DOthers[21];
    for(int it=0; it<20;it++){
        DOthers[it]=Others[it];
        cout << it << " it "<< endl;
        cout << DOthers[it]<< " DOthers[it] in eqwagesnumeric" << endl;
    }


    nlopt_set_min_objective(ExcessDemand, ExcessDemandsTotal,(void *)&DOthers);




