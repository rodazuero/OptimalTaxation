rm(list=ls(all=TRUE))

setwd('/Users/razuero/Dropbox/OptmalTaxationShared/Data/git/LastGitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/')

require(RcppArmadillo)



RcppArmadillo.package.skeleton("Calibration",example_code=FALSE)

R CMD build Calibration
R CMD check Calibration
R CMD INSTALL -l /Users/razuero/Dropbox/OptmalTaxationShared/Data/git/LastGitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2 Calibration_1.0.tar.gz


library(Calibration, lib.loc = "/Users/razuero/Dropbox/OptmalTaxationShared/Data/git/LastGitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/")






sourceCpp("main.cpp", verbose=TRUE, rebuild=TRUE)
Rcpp::compileAttributes()  