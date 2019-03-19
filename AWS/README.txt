1. Go to http://mpitutorial.com/tutorials/installing-mpich2/
2. Read the bullshit and download the mphch2:http://www.mpich.org/
3. On the directory where the zip file is stored, run: 
>>> tar -xzf mpich2-1.4.tar.gz
>>> cd mpich2-1.4
Note that versions need to be changed if updated. 
4. Run 
>>> ./configure --disable-fortran --prefix=/Users/rodrigoazuero/Documents/Parallelcomputing/mpich-3.2/mpich-install 2>&1 | tee c.txt

5. Run
>>> make 2>&1 | tee m.txt


6. Run
>>> sudo make install 2>&1 | tee mi.txt

7. Run
>>>  PATH=/Users/rodrigoazuero/Documents/Parallelcomputing/mpich-3.2/mpich-install/bin:$PATH ; export PATH

8. Test versions:
>>> which mpicc
>>> which mpiexec



6. Test version
>>> mpiexec --version
HYDRA build details:
    Version:                         3.1.4
    Release Date:                    Fri Feb 20 15:02:56 CST 2015
    CC:                              gcc
    CXX:                             g++
    F77:
    F90:

7. To run an example do this: (steps from http://mpitutorial.com/tutorials/mpi-hello-world/)

>>> cd "/Users/rodrigoazuero/Documents/Parallelcomputing"
>>> git clone https://github.com/wesleykendall/mpitutorial
>>> cd mpitutorial/tutorials/mpi-hello-world/code
>>> cat makefile
>>> export MPICC=/Users/rodrigoazuero/Documents/Parallelcomputing/mpich-3.2/mpich-install/bin/mpicc
>>> make
>>> export MPIRUN=/Users/rodrigoazuero/Documents/Parallelcomputing/mpich-3.2/mpich-install/bin/mpicc
>>> ./mpi_hello_world

Hello world from processor Rodrigos-MacBook-Pro.local, rank 0 out of 1 processors

Done!


************************************************************
 Now I need a cluster to run this. Amazon EC2 services. Let's see
************************************************************

Read: http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/EC2_GetStarted.html?r=1874


Instance id: i-29023eb1
Public instance id:ec2-54-88-197-63.com
Path to key: /Users/rodrigoazuero/Documents/Parallelcomputing/amazoncomputing

Enable inbound SSH traffic from your IP address to your instance
Ensure that the security group associated with your instance allows incoming SSH traffic from your IP address. For more information, see Authorizing Network Access to Your Instances. This can be checked here: https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/AccessingInstancesLinux.html

Security group. Security group should be created with:
1. HTTP/port range 80/ source 0.0.0.0/0
2. HTTPs/port range 443/ source 0.0.0.0/0
3. SSH/port range 443/ source ipaddress/32


#1. First need to set up the key as private. 
chmod 400 /path/my-key-pair.pem

*Access the server by:
ssh -i /Users/rodrigoazuero/Documents/Parallelcomputing/amazoncomputing/rodazuero-key-pair-nvirginia.pem ec2-user@ec2-52-91-171-108.compute-1.amazonaws.com

ssh -i /Users/rodrigoazuero/Documents/Parallelcomputing/amazoncomputing/rodazuero-key-pair-nvirginia.pem ec2-user@ec2-52-91-171-108.compute-1.amazonaws.com

ssh -i /Users/rodrigoazuero/Documents/Parallelcomputing/amazoncomputing/rodazuero-key-pair-nvirginia.pem ubuntu@ec2-52-207-212-62.compute-1.amazonaws.com


ssh -i /Users/rodrigoazuero/Dropbox/BACKUPRODRIGO/Parallelcomputing/amazoncomputing/RodrigoWashingtonDC.pem ec2-user@ec2-34-211-29-255.us-west-2.compute.amazonaws.com

ssh -i /Users/rodrigoazuero/Dropbox/BACKUPRODRIGO/Parallelcomputing/amazoncomputing/RodrigoWashingtonDC.pem ec2-user@ec2-18-236-221-228.us-west-2.compute.amazonaws.com


ssh -i /Users/rodrigoazuero/Dropbox/cons/attemptconsult.pem ubuntu@54.237.202.145

*Transfer files
scp -i /Users/rodrigoazuero/Dropbox/BACKUPRODRIGO/Parallelcomputing/amazoncomputing/RodrigoWashingtonDC.pem InAws.zip ec2-user@ec2-34-222-212-182.us-west-2.compute.amazonaws.com:~

scp -i /Users/rodrigoazuero/Documents/Parallelcomputing/amazoncomputing/rodazuero-key-pair-nvirginia.pem ATTEMPT6.zip ubuntu@ec2-54-172-36-135.compute-1.amazonaws.com:/srv/shiny-server/

scp -i /Users/rodrigoazuero/Documents/Parallelcomputing/amazoncomputing/rodazuero-key-pair-nvirginia.pem hello_world.cpp ec2-user@ec2-54-172-36-135.compute-1.amazonaws.com:~


scp -i /Users/rodrigoazuero/Documents/Parallelcomputing/amazoncomputing/rodazuero-key-pair-nvirginia.pem /Users/rodrigoazuero/Documents/Research/Chile/RR/BEHAVIORAL64/PARALLELATTEMPT/MainParallel.cpp ec2-user@ec2-52-87-167-229.compute-1.amazonaws.com:/home/ec2-user/BEHAVIORAL64


*INstall compiler shet

sudo yum install gcc-c++
sudo yum groupinstall "Development Tools"

*Preparing to c++* Boost library


wget https://sourceforge.net/projects/boost/files/boost/1.60.0/boost_1_60_0.tar.gz/download

*Transfer the boost library
scp -i /Users/rodrigoazuero/Documents/Parallelcomputing/amazoncomputing/rodazuero-key-pair-nvirginia.pem Filestotransfer/boost_1_60_0.tar.bz2 ec2-user@ec2-54-88-197-63.compute-1.amazonaws.com:~

*Transer nlopt
scp -i /Users/rodrigoazuero/Documents/Parallelcomputing/amazoncomputing/rodazuero-key-pair-nvirginia.pem Filestotransfer/nlopt-2.4.2.tar.gz ec2-user@ec2-54-172-36-135.compute-1.amazonaws.com:~

*Install boost library 
sudo yum install boost-devel
sudo yum install gsl-devel

Directory is:
/home/ec2-user


sudo wget http://download.opensuse.org/repositories/home:zhonghuaren/Fedora_23/home:zhonghuaren.repo
sudo yum install nlopt-c
cd /etc/yum.repos.d/

*nlopt file
tar -zxvf  nlopt-2.4.2.tar.gz
cd nlopt-2.4.2
./configure && make && sudo make install

*------------------------------------*
*Steps for using Amazon cloud server
*------------------------------------*

1. Go to https://console.aws.amazon.com/ec2/v2/home?region=us-east-1#
2. Click on Launch Instance. (On the right it should be a region where I am. e.g. N. Virginia).
3. Chose Free tier only.
4. Chose Amazon Liux AMI. 
5. Step 6. Configure security group: chose rodazuero_SG_nvirginia. 
6. Launch
7. Chose existing key (rodazuero-key-pair-nvirginia)
8. Get the launch id: i-8e79840d
9. On the instance, look for the public DNS: (ec2-52-87-238-57)
10. ssh -i /Users/rodrigoazuero/Documents/Parallelcomputing/amazoncomputing/rodazuero-key-pair-nvirginia.pem ec2-user@ec2-52-201-230-107.compute-1.amazonaws.com

*--------------------------------------------------------------*
*Steps for configuring instance to be able to compile c++ files*
*--------------------------------------------------------------*
0. Update everything
>>> sudo yum update

1. Install c++ compilers
>>> sudo yum install gcc-c++
>>> sudo yum install clang
>>> sudo yum groupinstall "Development Tools"

2. Install boost library and gel library:
>>> sudo yum install boost-devel
>>> sudo yum install gsl-devel

3. Install NLOPT library
Exit of server and run the following in my local terminal:
>>> cd "/Users/rodrigoazuero/Documents/Parallelcomputing/amazoncomputing"

4. Now transfer the file necessary to install nlopt

>>> scp -i /Users/rodrigoazuero/Documents/Parallelcomputing/amazoncomputing/rodazuero-key-pair-nvirginia.pem Filestotransfer/nlopt-2.4.2.tar.gz ec2-user@ec2-54-210-20-199.compute-1.amazonaws.com:~

Enter again to the server and check ls that nlopt is there. 
5. Unzip file:
>>> tar -zxvf  nlopt-2.4.2.tar.gz

6. Install the file:
>>> cd nlopt-2.4.2
>>> ./configure && make && sudo make install

{{{{{
(Optional Step

You can check if the libraries are ready and the compilers also:
>>> scp -i /Users/rodrigoazuero/Documents/Parallelcomputing/amazoncomputing/rodazuero-key-pair-nvirginia.pem hello_world.cpp ec2-user@ec2-52-87-238-57.compute-1.amazonaws.com:~

>>>g++ hello_world.cpp -o hello_world
>>>./hello_world
}}}}}



7. Transfer the necessary files to compile the desired optimization process. 
In my local terminal run:


>>> cd "/Users/rodrigoazuero/Documents/Parallelcomputing/amazoncomputing/Chileattempt"

>>> scp -i /Users/rodrigoazuero/Documents/Parallelcomputing/amazoncomputing/rodazuero-key-pair-nvirginia.pem filestotransfer.zip ec2-user@ec2-52-87-238-57.compute-1.amazonaws.com:~

scp -i /Users/rodrigoazuero/Documents/Parallelcomputing/amazoncomputing/rodazuero-key-pair-nvirginia.pem filetotransfer_debug.zip ec2-user@ec2-54-172-36-135.compute-1.amazonaws.com:~


>>> ssh -i /Users/rodrigoazuero/Documents/Parallelcomputing/amazoncomputing/rodazuero-key-pair-nvirginia.pem ec2-user@ec2-52-87-167-229.compute-1.amazonaws.com




filetotransfer_debug

8. Unzip and everything. 
Go back to server and unzip
>>> unzip filestotransfer.zip

9. Go to the directory and compile
>>> cd filestotransfer
>>>g++ main.cpp -lgsl -lnlopt -o main
OR
>>> g++ main.cpp -lgsl -lgslcblas -lnlopt -o main
>>> ./main


g++ MainParallel.cpp -lgsl -lgslcblas -lnlopt -o main

g++ main.cpp -lgsl -lgslcblas -lnlopt
g++ main.cpp -lgsl  -lnlopt -lgslcblas
g++ main.cpp  -lgslcblas -lnlopt -lgsl
g++ main.cpp  -lgslcblas -lgsl -lnlopt
g++ main.cpp -lnlopt -lgsl -lgslcblas 
g++ main.cpp  -lnlopt -lgsl -lgslcblas 
./a.out

g++ main.cpp -lgsl -lgslcblas -lnlopt


/Users/rodrigoazuero/Documents/Parallelcomputing/amazoncomputing/Chileattempt



*Finally, to transfer data from amazon server to my pc
scp -i /Users/rodrigoazuero/Documents/Parallelcomputing/amazoncomputing/rodazuero-key-pair-nvirginia.pem  ec2-user@ec2-52-201-230-107.compute-1.amazonaws.com:/home/ec2-user/filetotransfer_debug/PAROPTFOUND.csv /Users/rodrigoazuero/Documents/Parallelcomputing/amazoncomputing/downloads



43387.1final attempt 
Hello, World!
---------------
43525.3 evalF
Iteration=(1); Feval=43525.30145
42948.6 evalF
Iteration=(2); Feval=42948.63563



