#File designed to send files to amazon server
#Define environmental variables

#1. Specify version of model if needed (not needed)
#version=1 

#2. Amazon server name
SERV=ec2-34-222-212-182


#3. Local directory where file is stored
DIR=

#4. Key stored in:
KEYDIR=/Users/rodrigoazuero/Documents/Parallelcomputing/amazoncomputing/

#5. Key name
KEYNAME=rodazuero-key-pair-nvirginia.pem

#6. File to transfer from server to pc
FROMSERVERTOPC=/Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/gitVersion/OptimalTaxation/AWS/InAws/mainTesla.cpp
#scp -i /Users/rodrigoazuero/Documents/Parallelcomputing/amazoncomputing/rodazuero-key-pair-nvirginia.pem /Users/rodrigoazuero/Dropbox/OptmalTaxationShared/Data/git/OptimalTaxation/TheoreticalMoments/Equilibrium2/Equilibrium2/SobolDim14 ec2-user@ec2-52-87-167-229.compute-1.amazonaws.com:~
toserver:
	scp -i $(KEYDIR)$(KEYNAME) $(FROMSERVERTOPC) ec2-user@$(SERV).compute-1.amazonaws.com:~
	

accessserver:
	ssh -i $(KEYDIR)$(KEYNAME) ec2-user@$(SERV).compute-1.amazonaws.com


fromserver:
	scp -i $(KEYDIR)$(KEYNAME) ec2-user@$(SERV).compute-1.amazonaws.com:/home/ec2-user/$(FROMSERVERTOPC) /Users/rodrigoazuero/Downloads

#---------------------
#Housekeeping        #
#---------------------

#1. C++ compiler
CC=g++

#2. Flags
CFLAGS=-fopenmp -Wall  

#3. Name of output
OUTPUT=main
OUTPUTNLOPT=mainnlopt
OUTPUTSMOOTH=smooth
OUTPUTPARALLEL=paralleloutput
#4. Files to run in main:
SOURCE=mainTesla.cpp
SOURCEHOPSPACK=main_hopspack.cpp
SOURCENLOPT=main_nlopt.cpp
SOURCESMOOTH=mainSmoothingDistribution.cpp
SOURCEPARALLEL=mainTesla.cpp
SOURCEOPTIMIZER=mainOptimizer.cpp

#5. Home directory. in TESLA is /home/razu/
HOME=/home/ec2-user/

#6. Dynamic libraries 
DLIB=-lgsl -lgslcblas

#7. Set number of threads	
export OMP_NUM_THREADS=36

#-------------------------------------
#3. Include directories:
INCDIR1=/usr/lib64/R/library/RcppArmadillo/include/
INCDIR2=$(HOME)install/include/
INCDIR3=/home/ec2-user/R/x86_64-koji-linux-gnu-library/3.4/Rcpp/include
INCDIR4=/usr/lib64/R/include

#Generating the INC for all the paths ->this is for tesla
INC=$(INCDIR1) $(INCDIR2) $(INCDIR3) $(INCDIR4)
#INC=$(INCDIR3)
INC_PARAMS=$(foreach d, $(INC), -I$d)
#-------------------------------------

#3.1 Include directories in my pc
INCDIRMYPC1=PATH=/Users/rodrigoazuero/local/lib/gcc/x86_64-apple-darwin13.0.0/4.8.2/include
INCMYPC=$(INCDIRMYPC1) 
INC_PARAMS_MYPC=$(foreach d, $(INCMYPC), -I$d)
#-------------------------------------

#4. Lnlopt library 
LNLOPT=-L$(HOME)install/lib/ -lnlopt




configureserver:
	sudo yum install htop
	sudo yum install boost
	sudo yum install boost-devel
	sudo yum install gsl 
	sudo yum install gsl-devel
	sudo yum install gcc-c++


installR:
	sudo amazon-linux-extras install R3.4

InstallRcppArmadillo:
	sudo Rscript -e 'install.packages("RcppArmadillo", repos = "http://cran.us.r-project.org")'

myprog:
	$(CC) $(SOURCE) $(INC_PARAMS) $(LNLOPT)   $(CFLAGS) $(DLIB)  -o $(OUTPUT)

nlopt:
	$(CC) $(SOURCENLOPT) $(INC_PARAMS) $(LNLOPT)   $(CFLAGS) $(DLIB)  -o $(OUTPUTNLOPT) 

runnlopt: nlopt
	nohup ./$(OUTPUTNLOPT) > nohupNLOPT.out&

myprogmypc:
	$(CC) $(SOURCE) $(INC_PARAMS_MYPC) -lnlopt   $(CFLAGS) $(DLIB)  -o $(OUTPUT)

myproghopspack:
	$(CC) $(SOURCEHOPSPACK) $(INC_PARAMS) $(LNLOPT)   $(CFLAGS) $(DLIB)  -o $(OUTPUT)

nloptinst:
	wget http://ab-initio.mit.edu/nlopt/nlopt-2.4.2.tar.gz -P $(HOME)
	tar  -zxvf $(HOME)nlopt-2.4.2.tar.gz -C$(HOME)
	mkdir $(HOME)/install
	cd $(HOME)nlopt-2.4.2 && ./configure --prefix=$(HOME)install && make && make install

parallel: 
	$(CC) $(SOURCEPARALLEL) $(INC_PARAMS) $(LNLOPT)   $(CFLAGS) $(DLIB)  -o $(OUTPUTPARALLEL)

parallelrun: parallel
	nohup ./$(OUTPUTPARALLEL)  > nohupPARALLR.out&

optimizer: 
	$(CC) $(SOURCEOPTIMIZER) $(INC_PARAMS) $(LNLOPT)   $(CFLAGS) $(DLIB)  -o $(OUTPUTPARALLEL)

optimizerrun: optimizer
	nohup ./$(OUTPUTPARALLEL)  > nohupPARALLR.out&


paralleltesla:
	$(CC) $(SOURCEPARALLEL) -std=gnu++0x $(INC_PARAMS) $(LNLOPT) $(CFLAGS) $(DLIB) -o $(OUTPUTPARALLEL)

parallelruntesla: paralleltesla
	nohup ./$(OUTPUTPARALLEL)  > nohupPARALLR.out&

smoothdist: 
	$(CC) $(SOURCESMOOTH) $(INC_PARAMS) $(LNLOPT)   $(CFLAGS) $(DLIB)  -o $(OUTPUTSMOOTH)

smoothdistrun: smoothdist
	nohup ./$(OUTPUTSMOOTH) > nohupSMOOTHDIST.out&

nloptuninst:
	rm $(HOME)nlopt-2.4.2.tar.gz
	rm -r $(HOME)nlopt-2.4.2





hopspackinstall:
	wget http://www.sandia.gov/hopspack/code_drops/hopspack-2.0.2-src.tar.gz
	tar -zxvf hopspack-2.0.2-src.tar.gz
	cd hopspack-2.0.2-src && mkdir build_mt
	cd hopspack-2.0.2-src && cmake ../hopspack-2.0.2-src -Dmt=yes
	cd hopspack-2.0.2-src && make

hopspackrun: myproghopspack
	nohup ../hopspack-2.0.2-src/HOPSPACK_main_threaded hopspackCHILE.txt> nohuphopspack.out&

clean: 
	rm $(OUTPUT)


#In my pc: omp.h is in /Users/rodrigoazuero/local/lib/gcc/x86_64-apple-darwin13.0.0/4.8.2/include
#Include it to path PATH=$PATH:/Users/rodrigoazuero/local/lib/gcc/x86_64-apple-darwin13.0.0/4.8.2/include
#Also include the executable so that hopspack finds it. PATH=$PATH:/home/razu/BEHAVIORAL55
#Myproghopspack: g++ main.cpp  -I/usr/include/c++/3.4.6/backward/  -I/usr/lib64/R/library/BH/include/  -I/home/razu/install/include/ -L/home/razu/install/lib/ -lnlopt   -fopenmp -Wall   -lgsl -lgslcblas  -o main
#g++ main_hopspack.cpp  -I/usr/include/c++/3.4.6/backward/  -I/usr/lib64/R/library/BH/include/  -I/home/razu/install/include/ -L/home/razu/install/lib/ -lnlopt   -fopenmp -Wall   -lgsl -lgslcblas  -o main

#Locating C++  libraries in TESLA

#locate iostream.h->usr/include/c++/3.4.6/backward
#locate locate discrete_distribution.hpp->/usr/lib64/R/library/BH/include/
# //This is where i find the iostream.h

#After running nloptinst the nlopt library will be installed in /home/razu/install/lib

#In order to run it in the background: execute nohup ./main.out & 

# g++  main.cpp  -I/home/razu/install/include/ -L/home/razu/install/lib/ -lnlopt -Wall -lgsl -lgslcblas   -o main
# g++  main.cpp  -I/home/razu/install/include/ -I/usr/include/c++/3.4.6/backward/  -I/usr/lib64/R/library/BH/include/ -L/home/razu/install/lib/ -lnlopt -Wall -lgsl -lgslcblas   -o main
#g++   -I/usr/include/c++/3.4.6/backward/  -I/usr/lib64/R/library/BH/include/  -I/home/razu/install/include/ -L/home/razu/install/lib/ -lnlopt -Wall -lgsl -lgslcblas  main.cpp -omain


#g++ mainTesla.cpp -std=gnu++0x -I/usr/include/c++/3.4.6/backward/  -I/usr/lib64/R/library/BH/include/  -I/home/razu/install/include/ -I/usr/lib64/R/library/RcppArmadillo/include/ -I/home/razu/boost_1_61_0/ -L/home/razu/install/lib/ -lnlopt   -fopenmp -Wall   -lgsl -lgs

#g++ mainTesla.cpp  -I/usr/lib64/R/library/RcppArmadillo/include/ -I/home/ec2-user/install/include -L/home/ec2-user/install/lib -lnlopt -Wall -lgsl -lgslcblas   -o main

