
library(randtoolbox)

#Lessons learned by the naive sobol: 
#1. Decrease the values of ddelta. Have it between 0 and 300 
#2. Decrease chi. Between 0 and 200 


#Naive sobol with huge range. Very few equilibria found
set.seed(2581633)
Rand<-sobol(200000,dim=15,seed=2581633)
Rand[,1]<-Rand[,1]*(0.85-0.6)+0.6  #aalpha between 0.6 and 0.85
Rand[,2]<-Rand[,2]*(1-0.01)+0.01 ###ggamma between 0.01 and 1
Rand[,3]<-Rand[,3]*(10-0.01)+0.01 #ddelta between 0.1 and 10 (Originally 0-700)
Rand[,4]<-Rand[,4]*(700-0.1)+0.01 #bbeta  between 0.1 and 2
Rand[,5]<-Rand[,5]*(1-0.01)+0.01 ###ssigma  between 0.01 and 1
Rand[,6]<-Rand[,6]*(10-0.01)+0.01 #Kappa between 0.01 and 10
Rand[,7]<-Rand[,2] #Gamma = rrho as suggested by juan manuel. 
#Rand[,7]<-Rand[,7]*(1-0.01)+0.01 #rrho between 0.01 and 1. Update: We haveLarger than psi
Rand[,8]<-Rand[,8]*(1-0.01)+0.01  ###psi between 0.01 and 1
Rand[,9]<-Rand[,9]*(10-0.1)+0.01 #chi between 0.1 and 10 (Originally 0-700 but no eq>200)
Rand[,10]<-Rand[,10]*(3-0.3)+0.3 #Mmu1 between 0.3 and 3
Rand[,11]<-Rand[,11]*(3-0.3)+0.3#mmu3 between 0.3 and 3. 
Rand[,12]<-Rand[,12]*(1-0.3)+0.3 #ssigma1 between 0.3 and 1
Rand[,13]<-Rand[,13]*(1-0.3)+0.3 #ssigma2 between 0.5 and 3
Rand[,14]<-Rand[,14]*(1-0.01)+0.01  #rho12 between 0.01 and 1
Rand[,15]<-Rand[,15]*(30-0.01)+0.01 #c between 0.01 and 30

rows <- sample(nrow(Rand))
Rand <- Rand[rows, ]

#The C++ file is actually reading a 15 dimension sobol, we need to put them correctly
#And rho12
Rand<-round(Rand,2)
sobolsdir="/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/LocalCopy/OptimalTaxation/AWS/InAws/SobolDim15.csv"

write.table(Rand,file = sobolsdir, sep=",",  col.names=FALSE,row.names = FALSE)

