This repository holds the BS-SOLCTRA project used for visualization purposes.

To compile on Theta it's important to first switch the default Intel compiler suite to GNU:

$ module swap PrgEnv-intel PrgEnv-gnu

Compilation should be done through the provided Makefile by typing: 

$make

To execute, go into resultados/ and configure the different runtime parameters in submit_THETA.sh

To submit your job just type:

$ qsub submit_THETA.sh

To check the status of your job type:

$ qstat -u <username> 
 
or 

$ watch -n 5 qstat -u <username>

To check information on your project's allocation use:

$ sbank



