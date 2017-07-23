#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
     /* Declarations */
     const mxArray *coefsData, *parData, *spInfoData, *knotsData, *indexData;
    
    double *coefs, *params, *spInfo, *knots, *index;
    int N, splpieces, id, vecLength, sn;
    double speed, frR, nf;
    double splstart, splend, offsetjump, splos, tr, curt, t, obj;
    double ax, bx, cx, dx, ay, by, cy, dy;
    int i, ind, cnt, os;
    int knoffset, indexoffset;
    double tmp1, tmp2, tmp3;
    double px_, py_, pxs, pys, tsq,  sqrtpxspys;
    
    int indtmp;
    double brtmp;
    int allcnt;
    
    double *fx, *dfx;
    
    double t1,t2,t3,t4,t5,t6,c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y;
    double cmt1,cmt2,cmt3,cmt4,cmt5,cmt6;
    double t14,t24,t25,t35,t36;
    double c2x_cmt1_c1x_cmt4, c2y_cmt1_c1y_cmt4;
    int *svpos;
     
     /* Copy input pointers */
     coefsData = prhs[0];    /* state vec */
    parData = prhs[1];      /* speed, frame rate */
    spInfoData = prhs[2];   /* pieces, starts, ends */
    knotsData = prhs[3];   /* knots */
    indexData = prhs[4];    /* index */
    
    coefs  = mxGetPr(coefsData);
    params = mxGetPr(parData);
    spInfo = mxGetPr(spInfoData);
    knots = mxGetPr(knotsData);
    index  = mxGetPr(indexData);
    
    /* How many coefficients? */
    vecLength = (int)mxGetM(coefsData);
    
     /* How many splines? */
     N = (int)mxGetM(spInfoData);     
     mxAssert((int)mxGetN(spInfoData)==4, "spInfo must have 4 columns");
     

     
     
     /* // Allocate memory and assign output pointer */
      plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
      plhs[1] = mxCreateDoubleMatrix(vecLength, 1, mxREAL);
     
     /* //Get a pointer to the data space in our newly allocated memory  */
      fx = mxGetPr(plhs[0]);
      dfx = mxGetPr(plhs[1]);
     
     offsetjump=8; splos=0;
     knoffset = 0;
     indexoffset = 0;
     
      *fx=0;
     for (id = 0; id < vecLength; id++) {
         dfx[id]=0;
     }
     
}