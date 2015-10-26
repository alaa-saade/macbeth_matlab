MaCBetH : Matrix Completion with the Bethe Hessian

Files included in this package :
    complete.m : The main reconstruction code. 
    macbeth_demo.m : a sample code to test the completion function
    mex_all.m : a file to compile the mex files needed by the code
    subroutines folder : all subroutines needed by the different parts of the algorithm. 
   
IMPORTANT : it is highly recommended that the minFunc software from 
http://www.cs.ubc.ca/~schmidtm/Software/minFunc.html be installed and on the matlab path. Without it, MaCBetH will attempt to perform the local optimization using Matlab's Optimization Toolbox, resulting in very slow and memory intensive processing. 

USAGE : from matlab, run mex_all to compile the mex files. Make sure that the folder containing the minFunc files in on the matlab path. Then :
- macbeth_demo : runs macbeth on a synthetic example
- complete : is main reconstruction code that completes a partially observed matrix (assumed to be centered).
use help macbeth_demo and help complete for usage of these functions. 

Comments and remarks regarding bugs or functionalities are more than welcome :-).
