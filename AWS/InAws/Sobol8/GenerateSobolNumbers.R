
library(randtoolbox)

#Lessons learned by the naive sobol: 
#1. Decrease the values of ddelta. Have it between 0 and 300 
#2. Decrease chi. Between 0 and 200 








#We will generate sobols centered around the optimized version of the parameters
#that were last found. 
#Let's start by reading the file storing the parameters


ParamsInitial<-read.csv('/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/LastGitVersion/OptimalTaxation/AWS/InAws/Sobol7/Modelfit/Model1/Param.csv', header = T, sep=',')

#Naive sobol with huge range. Very few equilibria found
set.seed(2581633)
Rand<-sobol(10000,dim=15,seed=2581633)
Rand[,1]<-ParamsInitial$aalpha+(Rand[,1]-0.5)/5  #aalpha between 0.67 and 0.82
Rand[,2]<-ParamsInitial$ggamma+(Rand[,2]-0.5)*2 #ggamma between 1.69 and 3.69
Rand[,3]<-ParamsInitial$ddelta+(Rand[,3]-0.5)*100 #ddelta  (Originally range of 100)
Rand[,4]<-ParamsInitial$bbeta+(Rand[,4]-0.5)*100 #bbeta  -range of 100
Rand[,5]<-ParamsInitial$ssigma+(Rand[,5]-0.5)*4 ###ssigma  range of 4
Rand[,6]<-ParamsInitial$kkappa+(Rand[,6]-0.5)*0.5 #Kappa range of 0.5
Rand[,7]<-ParamsInitial$rho+(Rand[,7]-0.5)*40 #rho range of 40. 
Rand[,8]<-ParamsInitial$psi+(Rand[,8]-0.5)*40 ###psi range of 40. 
Rand[,9]<-ParamsInitial$chi+(Rand[,9]-0.5)*10 #chi range of 10
Rand[,10]<-ParamsInitial$mmu1+(Rand[,10]-0.5)*4 #Mmu1 range of 2
Rand[,11]<-ParamsInitial$mmu2 +(Rand[,11]-0.0)*5 #mmu2 range of 2
Rand[,12]<-ParamsInitial$ssigma1+(Rand[,12]-0.5)*2 #ssigma1 range of 2
Rand[,13]<-ParamsInitial$ssigma2+(Rand[,13]-0.5)*2 #ssigma2 range of 2
Rand[,14]<-ParamsInitial$rho12+(Rand[,14]-0.5)*0.5   #rho12 range of 0.5
Rand[,15]<-(Rand[,15])*200  #c between 3 and 40
range(Rand[,15])

#The C++ file is actually reading a 15 dimension sobol, we need to put them correctly
#And rho12

Rand<-round(Rand,2)
sobolsdir="/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/LastGitVersion/OptimalTaxation/AWS/InAws/Sobol8/SobolDim15.csv"
write.table(Rand,file = sobolsdir, sep=",",  col.names=FALSE,row.names = FALSE)

