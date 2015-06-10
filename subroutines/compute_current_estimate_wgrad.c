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
 void compute_current_estimate_wgrad(double *x, double *I, double *J,double *VAL,double n, double m, double r,double nnz,double *COST,double *grad)
 {

    for(int foo=0;foo<(n+m)*r;foo++){
        grad[foo] = 0.0;
    }
    COST[0] = 0.0;    

    for(int k=0;k<nnz;k++){

        double xyt = 0.0;
        int i = (int)I[k];
        int j = (int)J[k];

        for(int mu=1;mu<r+1;mu++){
            xyt += x[(int)(i+(mu-1)*n-1)]*x[(int)(n*r+j+(mu-1)*m-1)];
        }

        COST[0] += (VAL[k]-xyt)*(VAL[k]-xyt);

        for(int mu=1;mu<r+1;mu++){
            grad[(int)(i+(mu-1)*n-1)] += 2.0*(-VAL[k] + xyt)*x[(int)(n*r+j+(mu-1)*m-1)];   
            grad[(int)(n*r+j+(mu-1)*m-1)] += 2.0*(-VAL[k] + xyt)*x[(int)(i+(mu-1)*n-1)];
        }

    }
}
    // The gateway function 
void mexFunction( int nlhs, mxArray *plhs[],
  int nrhs, const mxArray *prhs[])
{


    double *x;                  
    double *I;
    double *J;
    double *VAL;

    double n;
    double m;
    double r;
    double nnz;

    double *COST;
    double *grad;             

    /* check for proper number of arguments */
    if(nrhs!=8) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","8 inputs required: x, I, J, VAL, n, m, r, nnz");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","2 outputs required: COST, grad");
    }


    x = mxGetPr(prhs[0]);
    I = mxGetPr(prhs[1]);
    J = mxGetPr(prhs[2]);
    VAL = mxGetPr(prhs[3]);
    n = (int)mxGetScalar(prhs[4]);
    m = (int)mxGetScalar(prhs[5]);
    r = (int)mxGetScalar(prhs[6]);
    nnz = mxGetScalar(prhs[7]);


    /* create the output vectors */
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix((n+m)*r,1,mxREAL);

    /* get a pointer to the real data in the output vectors */
    COST = mxGetPr(plhs[0]);
    grad = mxGetPr(plhs[1]);

    /* call the computational routine */
    compute_current_estimate_wgrad(x,I,J,VAL,n,m,r,nnz,COST,grad);
}
