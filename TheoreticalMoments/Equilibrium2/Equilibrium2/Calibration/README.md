Calibration
=======


***This package computes the estimation and generates the theoretical moments predicted by the optimal parameters in Azuero, Hernandez, and Wills (2020)***


Three functions are created and called with the Calibration:: framework. The package requires nlopt, RcppArmadillo, and Boost. 

```{r}
sudo yum install nlopt Boost
```

Modifications in the locadl Makevars file might be necessary depending on the location of the RcppArmadillo headers. 
The package generates three functions: `DistanceEstimator`, `EqWages`, and `EqWagesNumericVector`. They are called in the following way:


```{r}
# CRAN install.  
library(Calibration, lib.loc = "Local/Dir")

Calibration::DistanceEstimator(Parameters,InitialWages)

```
Parameters are described in the CPP file and InitialWages are the initial suggestion to find the equilibrium. 


### Instructions to re-create the package when the cpp file is modified:

The cpp file needs to be included in the /src/ folder. 

Once everything is ready to be generated, the cpp code needs to be compiled in R. Do this with:

```{r}
Rcpp::compileAttributes()  
```

It is important to note that we need to modify the RcppExports.cpp needs to be modified
manually to include the following line:

```
using std::vector;
```
Afterwards

```
R CMD build Calibration
R CMD check Calibration
R CMD install L -l /local/dir/ Calibration_1.0.tar.gz
```





