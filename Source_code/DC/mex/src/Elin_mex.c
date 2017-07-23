#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
     /* Declarations */
     const mxArray *coefsData, *parData, *spInfoData, *breaksData, *indexData;
    
    double *coefs, *params, *spInfo, *breaks, *index;
    int N, splpieces, id, vecLength;
    double speed, frR, nf;
    double splstart, splend, offsetjump, splos, tr, t, obj;
    double ax, bx, cx, dx, ay, by, cy, dy;
    int i, ind, cnt, os;
    int broffset, indexoffset;
    double tmp1, tmp2, tmp3;
    double px_, py_, pxs, pys, tsq,  sqrtpxspys;
    
    int indtmp;
    double brtmp;
    int allcnt;
    
    double *fx, *dfx;
     
     /* Copy input pointers */
     coefsData = prhs[0];    /* state vec */
    parData = prhs[1];      /* speed, frame rate */
    spInfoData = prhs[2];   /* pieces, starts, ends */
    breaksData = prhs[3];   /* breaks */
    indexData = prhs[4];    /* index */
    
    coefs  = mxGetPr(coefsData);
    params = mxGetPr(parData);
    spInfo = mxGetPr(spInfoData);
    breaks = mxGetPr(breaksData);
    index  = mxGetPr(indexData);
    
    /* How many coefficients? */
    vecLength = (int)mxGetM(coefsData);
    
     /* How many splines? */
     N = (int)mxGetM(spInfoData);     
     mxAssert((int)mxGetN(spInfoData)==3, "spInfo must have 3 columns");
     
     /* Get parameters */
     nf = params[0]; speed = params[1];     frR = params[2];
     
     
     /* // Allocate memory and assign output pointer */
      plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
      plhs[1] = mxCreateDoubleMatrix(vecLength, 1, mxREAL);
     
     /* //Get a pointer to the data space in our newly allocated memory  */
      fx = mxGetPr(plhs[0]);
      dfx = mxGetPr(plhs[1]);
     
     offsetjump=8; splos=0;
     broffset = 0;
     indexoffset = 0;
     
      *fx=0;
     for (id = 0; id < vecLength; id++) {
         dfx[id]=0;
     }
     
     
    /* for each target */
      allcnt=0;
     for (id = 0; id < N; id++) {
        
        /* get spline Info */
        splpieces = (int)spInfo[id];
        splstart =  spInfo[id+N];
        splend =    spInfo[id+2*N];
        
        /* for all frames */
        cnt=0;
        for (tr = splstart; tr <= splend; tr += 1.)
        {
             indtmp=(int)index[cnt+indexoffset]-1;
             brtmp=breaks[indtmp + broffset];
                         
              t = tr - brtmp;
             
              os = ((int)(index[cnt+indexoffset])-1)*offsetjump + splos;
              ax=coefs[0+os];bx=coefs[1+os];cx=coefs[2+os];dx=coefs[3+os];
              ay=coefs[4+os];by=coefs[5+os];cy=coefs[6+os];dy=coefs[7+os];
              
              tsq=t*t;
			  px_=(cx + 2*bx*t + 3*ax*tsq);
			  py_=(cy + 2*by*t + 3*ay*tsq);
			  pxs=px_*px_; pys=py_*py_;
			  
			  sqrtpxspys=sqrt(pxs+pys);
        

			obj=(speed - (frR*sqrtpxspys)/nf)*(speed - (frR*sqrtpxspys)/nf); 
			
			*fx+=obj;
			dfx[0 + os] += -(6*frR*tsq*(speed - (frR*sqrtpxspys)/nf)*px_)/(nf*sqrtpxspys);
			dfx[1 + os] += -(4*frR*t*(speed - (frR*sqrtpxspys)/nf)*px_)/(nf*sqrtpxspys);
			dfx[2 + os] += -(frR*(speed - (frR*sqrtpxspys)/nf)*(2*cx + 4*bx*t + 6*ax*tsq))/(nf*sqrtpxspys);
			/* dfx[3 + os] += 0;*/
			dfx[4 + os] += -(6*frR*tsq*(speed - (frR*sqrtpxspys)/nf)*py_)/(nf*sqrtpxspys);
			dfx[5 + os] += -(4*frR*t*(speed - (frR*sqrtpxspys)/nf)*py_)/(nf*sqrtpxspys);
			dfx[6 + os] += -(frR*(speed - (frR*sqrtpxspys)/nf)*(2*cy + 4*by*t + 6*ay*tsq))/(nf*sqrtpxspys);
			/* dfx[7 + os] += 0; */

              cnt++;
        }
        splos += splpieces*8;
        broffset += (splpieces+1);
        indexoffset += cnt;
    }
    
    
    
}