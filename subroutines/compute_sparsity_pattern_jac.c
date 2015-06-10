/*==========================================================
 *
 * compute_current_estimate.c computes the current estimate 
 * of the revealed entries.
 *
 * Necessary for optimization subroutines
 *
 * This is a MEX-file for MATLAB. Compile with "mex compute_current_estimate.c" 
 *
 *========================================================*/

#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* The computational routine */
 void compute_sparsity_pattern_jac(double *I, double *J, double nnz,double n, double m, double r,double *u,double *v)
 {

    for(int k=0;k<nnz;k++){
        int i = (int)I[k];
        int j = (int)J[k];

        for(int mu=1;mu<r+1;mu++){
            u[(int)(2*r*k+mu-1)] = k+1;
            v[(int)(2*r*k+mu-1)] = i+(mu-1)*n;
            u[(int)(2*r*k+r+mu-1)] = k+1;
            v[(int)(2*r*k+r+mu-1)] = n*r+j+(mu-1)*m;
        }
    }
}
    // The gateway function 
void mexFunction( int nlhs, mxArray *plhs[],
  int nrhs, const mxArray *prhs[])
{


    double *I;
    double *J;

    double n;
    double m;
    double r;
    double nnz;

    double *u;
    double *v;             

    /* check for proper number of arguments */
    if(nrhs!=6) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","6 inputs required: I, J, nnz, n, m, r");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","2 outputs required: u, v");
    }

    I = mxGetPr(prhs[0]);
    J = mxGetPr(prhs[1]);
    nnz = mxGetScalar(prhs[2]);
    n = (int)mxGetScalar(prhs[3]);
    m = (int)mxGetScalar(prhs[4]);
    r = (int)mxGetScalar(prhs[5]);
    
    /* create the output vectors */
    plhs[0] = mxCreateDoubleMatrix(1,nnz*2*r,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,nnz*2*r,mxREAL);

    /* get a pointer to the real data in the output vectors */
    u = mxGetPr(plhs[0]);
    v = mxGetPr(plhs[1]);

    /* call the computational routine */
    compute_sparsity_pattern_jac(I,J,nnz,n,m,r,u,v);
}
