#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
     /* Declarations */
     const mxArray *coefsData, *parData, *spInfoData, *breaksData, *indexData;
    
    double *coefs, *params, *spInfo, *breaks, *index;
    int N, splpieces, id, vecLength;
    double speed, frR;
    double splstart, splend, offsetjump, splos, tr, t, obj;
    double ax, bx, cx, dx, ay, by, cy, dy;
    int i, ind, cnt, os;
    int broffset, indexoffset;
    double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11;
    double px, py, pxs, pys, tsq,  sqrtpxspys, frRsq, twofrRsq, pxspyssq, pxspyscu;
	double px_, py_, px__, py__, px_spy_ssq, k,ksq,kqu,_1_ksq, epsil;
    
    int indtmp;
    double brtmp;

    
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
     speed = params[0];     frR = params[1]; frRsq=frR*frR; twofrRsq=2*frRsq;
	 k=.1; ksq=k*k; kqu=k*k*k*k; _1_ksq=1.0/ksq;
	 epsil=.1;
     
     
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
     
     //printf("%f %f\n",speed,frR);
     
    /* for each target */
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
              ax=coefs[0+os];bx=coefs[1+os];cx=coefs[2+os];
              ay=coefs[4+os];by=coefs[5+os];cy=coefs[6+os];
              
              tsq=t*t;
			  px_=(cx + 2*bx*t + 3*ax*tsq);
				py_=(cy + 2*by*t + 3*ay*tsq);
				 px__=(2*bx + 6*ax*t);
				 py__=(2*by + 6*ay*t);
				// tmp1=(tmp1);
				// tmp2=tmp2.;
				
              //px_spy_ssq=(px_^2. + py_^2.);
              pxspyssq = pxs+pys; pxspyssq *= pxspyssq; pxspyscu = pxspyssq*(pxs+pys);
              tmp1=py__ * px_ - px__ * py_;
			  tmp2=tmp1*tmp1;
			  
			  
			  tmp3=px_*px_ + py_*py_;
			  tmp4=tmp3*tmp3*tmp3;
			  tmp5=sqrt(tmp2/ksq + 1);
			  tmp6=(tmp5 - 1)*(tmp5 - 1);
			  tmp7=tmp3*tmp3;
			  /*
			  
			obj=-log(1./((frRsq*kqu*tmp6)/tmp7 + 1.));
			*fx+=obj;
			dfx[0 + os] += -((12*frRsq*kqu*tsq*tmp6*px_)/tmp4 + (2*frRsq*ksq*(6*t*py_ - 3*tsq*py__)*(tmp5 - 1)*(tmp1))/(tmp7*tmp5))/((frRsq*kqu*tmp6)/tmp7 + 1);
			dfx[1 + os] += -((8*frRsq*kqu*t*tmp6*px_)/tmp4 + (2*frRsq*ksq*(tmp5 - 1)*(tmp1)*(2*cy - 2*t*py__ + 4*by*t + 6*ay*tsq))/(tmp7*tmp5))/((frRsq*kqu*tmp6)/tmp7 + 1);
			dfx[2 + os] += -((2*frRsq*kqu*tmp6*(2*cx + 4*bx*t + 6*ax*tsq))/tmp4 - (2*frRsq*ksq*py__*(tmp5 - 1)*(tmp1))/(tmp7*tmp5))/((frRsq*kqu*tmp6)/tmp7 + 1);
			dfx[3 + os] += 0;
			dfx[4 + os] += -((12*frRsq*kqu*tsq*tmp6*py_)/tmp4 - (2*frRsq*ksq*(6*t*px_ - 3*tsq*px__)*(tmp5 - 1)*(tmp1))/(tmp7*tmp5))/((frRsq*kqu*tmp6)/tmp7 + 1);
			dfx[5 + os] += -((8*frRsq*kqu*t*tmp6*py_)/tmp4 - (2*frRsq*ksq*(tmp5 - 1)*(tmp1)*(2*cx - 2*t*px__ + 4*bx*t + 6*ax*tsq))/(tmp7*tmp5))/((frRsq*kqu*tmp6)/tmp7 + 1);
			dfx[6 + os] += -((2*frRsq*kqu*tmp6*(2*cy + 4*by*t + 6*ay*tsq))/tmp4 + (2*frRsq*ksq*px__*(tmp5 - 1)*(tmp1))/(tmp7*tmp5))/((frRsq*kqu*tmp6)/tmp7 + 1);
			dfx[7 + os] += 0;
			*/
			obj=-log(1/((frRsq*(epsil + tmp2))/tmp7 + 1));
			*fx+=obj;
			dfx[0 + os] += -((2*frRsq*(6*t*py_ - 3*tsq*py__)*tmp1)/tmp7 + (12*frRsq*tsq*(epsil + tmp2)*px_)/tmp4)/((frRsq*(epsil + tmp2))/tmp7 + 1);
			dfx[1 + os] += -((2*frRsq*tmp1*(2*cy - 2*t*py__ + 4*by*t + 6*ay*tsq))/tmp7 + (8*frRsq*t*(epsil + tmp2)*px_)/tmp4)/((frRsq*(epsil + tmp2))/tmp7 + 1);
			dfx[2 + os] += -((2*frRsq*(epsil + tmp2)*(2*cx + 4*bx*t + 6*ax*tsq))/tmp4 - (2*frRsq*py__*tmp1)/tmp7)/((frRsq*(epsil + tmp2))/tmp7 + 1);
			dfx[4 + os] += ((2*frRsq*(6*t*px_ - 3*tsq*px__)*tmp1)/tmp7 - (12*frRsq*tsq*(epsil + tmp2)*py_)/tmp4)/((frRsq*(epsil + tmp2))/tmp7 + 1);
			dfx[5 + os] += ((2*frRsq*tmp1*(2*cx - 2*t*px__ + 4*bx*t + 6*ax*tsq))/tmp7 - (8*frRsq*t*(epsil + tmp2)*py_)/tmp4)/((frRsq*(epsil + tmp2))/tmp7 + 1);
			dfx[6 + os] += -((2*frRsq*(epsil + tmp2)*(2*cy + 4*by*t + 6*ay*tsq))/tmp4 + (2*frRsq*px__*tmp1)/tmp7)/((frRsq*(epsil + tmp2))/tmp7 + 1);
 

              cnt++;
        }
        splos += splpieces*8;
        broffset += (splpieces+1);
        indexoffset += cnt;
    }
    
    
    
}