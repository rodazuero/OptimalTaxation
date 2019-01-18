
library(randtoolbox)

#Lessons learned by the naive sobol: 
#1. Decrease the values of ddelta. Have it between 0 and 300 
#2. Decrease chi. Between 0 and 200 


#Naive sobol with huge range. Very few equilibria found

Rand<-sobol(20000,dim=14,seed=2581633)
Rand[,1]<-Rand[,1]*(0.9-0.01)+0.01  #Between 0.1 and 0.9
Rand[,2]<-Rand[,2]*(10-0.1)+0.01 ###ggamma between 0.1 and 10
Rand[,3]<-Rand[,3]*(700-0.1)+0.01 #ddelta between 0.1 and 5 (Originally 0-700)
Rand[,4]<-Rand[,4]*(700-0.1)+0.01 #bbeta  between 0.1 and 2
Rand[,5]<-Rand[,5]*(10-0.1)+0.01 ###ssigma  between 0.1 and 10
Rand[,6]<-Rand[,6]*(700-0.1)+0.01 #Kappa between 0.1 and 2
Rand[,7]<-Rand[,7]*(4-0.25)+0.25 ###psi between 0.25 and 4. 
Rand[,8]<-Rand[,8]*(600-0.1)+0.01 #chi between 0.1 and 4 (Originally 0-700 but no eq>200)
Rand[,9]<-Rand[,9]*(10-0.1)+0.01 #rrho between 0.1 and 10. Update: We haveLarger than psi
Rand[,10]<-Rand[,10]*4+0.01 #Mmu1 between 0.5 and 3
Rand[,11]<-Rand[,11]*4 +0.01#mmu3 between 0.5 and 3. 
Rand[,12]<-Rand[,12]*5+0.01 #ssigma1 between 0.5 and 3
Rand[,13]<-Rand[,13]*5+0.01 #ssigma2 between 0.5 and 3
Rand[,14]<-Rand[,14]*1+0.01   #rho12 between 0 and 1
Rand<-cbind(Rand,Rand[,14])


#The C++ file is actually reading a 15 dimension sobol, we need to put them correctly
#And rho12


write.table(Rand,file = "/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolsGenerated/SobolDim15.csv", sep=",",  col.names=FALSE,row.names = FALSE)


dim(Rand)

