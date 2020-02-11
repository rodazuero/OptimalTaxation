
library(randtoolbox)

#Lessons learned by the naive sobol: 
#1. Decrease the values of ddelta. Have it between 0 and 300 
#2. Decrease chi. Between 0 and 200 








#We will generate sobols centered around the optimized version of the parameters
#that were last found. 
#Let's start by reading the file storing the parameters
ParamsInitial<-read.csv("/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/LastGitVersion/OptimalTaxation/AWS/InAws/Sobol4/Modelfit/Model2/Param.csv", header = T, sep=",")




#Naive sobol with huge range. Very few equilibria found
set.seed(2581633)
Rand<-sobol(10000,dim=15,seed=2581633)
Rand[,1]<-ParamsInitial$aalpha+(Rand[,1]-0.8)/4  #aalpha between 0.67 and 0.92
Rand[,2]<-ParamsInitial$ggamma+(Rand[,2]-0.5)*6.5 #ggamma between 0.13 and 6.62
Rand[,3]<-Rand[,3]*(700-0.1)+0.01 #ddelta  (Originally 0-700)
Rand[,4]<-Rand[,4]*(700-0.1)+0.01 #bbeta  between 0.1 and 2
Rand[,5]<-Rand[,5]*(10-0.1)+0.01 ###ssigma  between 0.1 and 10
Rand[,6]<-Rand[,6]*(10-0.1)+0.01 #Kappa between 0.1 and 2
Rand[,7]<-Rand[,7]*(700-0.1)+0.01 #359 rho. 
Rand[,8]<-Rand[,8]*(10-0.25)+0.25 ###psi between 0.25 and 4. 
Rand[,9]<-Rand[,9]*(100-0.1)+0.01 #chi between 0 and 100
Rand[,10]<-Rand[,10]*5+0.01 #Mmu1 between 0.5 and 3
Rand[,11]<-Rand[,11]*5 +0.01#mmu3 between 0.5 and 3. 
Rand[,12]<-Rand[,12]*5+0.01 #ssigma1 between 0.5 and 3
Rand[,13]<-Rand[,13]*5+0.01 #ssigma2 between 0.5 and 3
Rand[,14]<-Rand[,14]*1+0.01   #rho12 between 0 and 1
Rand[,15]<-Rand[,15]*30   #c between 0 and 100


#The C++ file is actually reading a 15 dimension sobol, we need to put them correctly
#And rho12

Rand<-round(Rand,2)
sobolsdir="/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/LastGitVersion/OptimalTaxation/AWS/InAws/Sobol5/SobolDim15.csv"
write.table(Rand,file = sobolsdir, sep=",",  col.names=FALSE,row.names = FALSE)

