Calibration
=======


***This package computes the estimation and generates the theoretical moments predicted by the optimal parameters in Azuero, Hernandez, and Wills (2020)***


Three functions are created and called with the Calibration:: framework. The package requires nlopt, RcppArmadillo, and Boost. 

```{r}
sudo yum install nlopt Boost
```

Modifications in the local Makevars file might be necessary depending on the location of the RcppArmadillo headers. 
The package generates three functions: `DistanceEstimator`, `EqWages`, and `EqWagesNumericVector`. They are called in the following way:


```{r}
# CRAN install.  
library(Calibration, lib.loc = "Local/Dir")

Calibration::DistanceEstimator(Parameters,InitialWages)

```
Parameters are described in the CPP file and InitialWages are the initial suggestion to find the equilibrium. 
