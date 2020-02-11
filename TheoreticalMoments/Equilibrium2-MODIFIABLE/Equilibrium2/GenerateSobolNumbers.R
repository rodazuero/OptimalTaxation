
library(randtoolbox)

#Lessons learned by the naive sobol: 
#1. Decrease the values of ddelta. Have it between 0 and 300 
#2. Decrease chi. Between 0 and 200 


#Naive sobol with huge range. Very few equilibria found
set.seed(2581633)
Rand<-sobol(50000,dim=15,seed=2581633)


#This is the pure g


#We will generate sobols centered around the optimized version of the parameters
#that were last found. 
#Let's start by reading the file storing the parameters
ParamsInitial<-read.csv("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/LastGitVersion/OptimalTaxation/AWS/InAws/Sobol4/Modelfit/Model1/Param.csv", header = T, sep=",")







#The C++ file is actually reading a 15 dimension sobol, we need to put them correctly
#And rho12

Rand<-round(Rand,2)
sobolsdir="/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/AWS/InAws/Sobol3/SobolDim15.csv"

write.table(Rand,file = sobolsdir, sep=",",  col.names=FALSE,row.names = FALSE)


dim(Rand)
1613-1625
48-67
