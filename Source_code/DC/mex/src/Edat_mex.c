#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    /* Declarations */
    const mxArray *coefsData, *parData, *spInfoData, *breaksData, *indexData, *ptsData;
    
    double *coefs, *params, *spInfo, *breaks, *index, *alldpoints;
    int N, Npts,  id, vecLength;
    double k, ksq, nf, nfsq;
    double offsetjump, splos, t, obj;
    double ax, bx, cx, dx, ay, by, cy, dy;
    int os, det, df;
    int broffset, indexoffset;
    double tmp1;
    double px, py, pxDx, pyDy, tsq, tcu;
    double Dx, Dy, sc;
    
	double spt;
    double ptx,pty,pts;
    int ptt, ptl,  cos;
    	
    
    double *fx, *dfx;
    
    /* Copy input pointers */
    coefsData = prhs[0];    /* state vec */
    parData = prhs[1];      /* speed, frame rate */
    spInfoData = prhs[2];   /* pieces, starts, ends */
    breaksData = prhs[3];   /* breaks */
    indexData = prhs[4];    /* index */
    ptsData = prhs[5];      /* all points */
    
    coefs  = mxGetPr(coefsData);
    params = mxGetPr(parData);
    spInfo = mxGetPr(spInfoData);
    breaks = mxGetPr(breaksData);
    index  = mxGetPr(indexData);
    alldpoints = mxGetPr(ptsData);
    
    /* How many coefficients? */
    vecLength = (int)mxGetM(coefsData);
    
    /* How many splines? */
    N = (int)mxGetM(spInfoData);
    mxAssert((int)mxGetN(spInfoData)==3, "spInfo must have 3 columns");
    
    /* How many detections? */
    Npts = (int)mxGetN(ptsData);
    
    /* Get parameters */
	df = (int)params[0];
    k = params[2]; ksq=k*k;    nf = params[1]; nfsq=nf*nf;
    
    
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
    // for (id = 0; id < N; id++) {
    
    
    for (det = 0; det < Npts; det++) {
        
        ptl=(int)alldpoints[det*7+4];  /* label of detection */
        
        if (ptl > 0) {
            ptx=alldpoints[det*7];         /* x */
            pty=alldpoints[det*7+1];       /* y */
            pts=alldpoints[det*7+2];       /* scores */
            ptt=(int)alldpoints[det*7+3];  /* t of detection */
            spt=alldpoints[det*7+5];  /* t to evaluate spline */
            cos=(int)alldpoints[det*7+6];  /* state vector offset */
            
//             printf("%f %f %f %i %i %i %i\n",ptx,pty,pts,ptt,ptl, spt,cos];
            
            /* only consider labeled points */
            
            
            os = cos;
            ax=coefs[0+os];bx=coefs[1+os];cx=coefs[2+os];dx=coefs[3+os];
            ay=coefs[4+os];by=coefs[5+os];cy=coefs[6+os];dy=coefs[7+os];
            
            Dx=ptx;Dy=pty; sc=pts;
            
            t=spt;
            
            tsq = t*t; tcu=tsq*t;
            //px = ax*tcu + bx*tsq + cx*t + dx;
            //py = ay*tcu + by*tsq + cy*t + dy;
			px=t*(t*(t*ax+bx)+cx)+dx;
			py=t*(t*(t*ay+by)+cy)+dy;
            pxDx = px-Dx;
			pxDx = dx - Dx + cx * t + ax*tcu + bx*tsq;
            pyDy = py-Dy;
            
            
            if (df == 2) {
            obj = (sc*(pxDx*pxDx + pyDy*pyDy))/nfsq;
            *fx += obj;
		      
            
            dfx[0 + os] = dfx[0 + os] + (2*sc*tcu*pxDx)/nfsq;
            dfx[1 + os] = dfx[1 + os] + (2*sc*tsq*pxDx)/nfsq;
            dfx[2 + os] = dfx[2 + os] + (2*sc*t*pxDx)/nfsq;
            dfx[3 + os] = dfx[3 + os] + (sc*(2*dx - 2*Dx + 2*cx*t + 2*ax*tcu + 2*bx*tsq))/nfsq;
            dfx[4 + os] = dfx[4 + os] + (2*sc*tcu*pyDy)/nfsq;
            dfx[5 + os] = dfx[5 + os] + (2*sc*tsq*pyDy)/nfsq;
            dfx[6 + os] = dfx[6 + os] + (2*sc*t*pyDy)/nfsq;
            dfx[7 + os] = dfx[7 + os] + (sc*(2*dy - 2*Dy + 2*cy*t + 2*ay*tcu + 2*by*tsq))/nfsq;
			}
			else if (df == 1 || df == 4) {
			/* Pseudo Huber
			tmp1=sqrt((pxDx*pxDx + pyDy*pyDy)/(ksq*nfsq) + 1);
			obj=k*sc*(tmp1 - 1);
			*fx+=obj;
			dfx[0 + os] += (k*sc*tcu*pxDx)/(ksq*nfsq*tmp1);
			dfx[1 + os] += (k*sc*tsq*pxDx)/(ksq*nfsq*tmp1);
			dfx[2 + os] += (k*sc*t*pxDx)/(ksq*nfsq*tmp1);
			dfx[3 + os] += (k*sc*(2*dx - 2*Dx + 2*cx*t + 2*ax*tcu + 2*bx*tsq))/(2*ksq*nfsq*tmp1);
			dfx[4 + os] += (k*sc*tcu*pyDy)/(ksq*nfsq*tmp1);
			dfx[5 + os] += (k*sc*tsq*pyDy)/(ksq*nfsq*tmp1);
			dfx[6 + os] += (k*sc*t*pyDy)/(ksq*nfsq*tmp1);
			dfx[7 + os] += (k*sc*(2*dy - 2*Dy + 2*cy*t + 2*ay*tcu + 2*by*tsq))/(2*ksq*nfsq*tmp1);
			*/ 
			
			/* simple Charbonnier */
			tmp1=sqrt(k + (pxDx*pxDx + pyDy*pyDy)/nfsq);
			obj=sc*tmp1;
			*fx+=obj;
			dfx[0 + os] += (sc*tcu*(pxDx))/(nfsq*tmp1);
			dfx[1 + os] += (sc*tsq*(pxDx))/(nfsq*tmp1);
			dfx[2 + os] += (sc*t*(pxDx))/(nfsq*tmp1);
			dfx[3 + os] += (sc*(2*dx - 2*Dx + 2*cx*t + 2*ax*tcu + 2*bx*tsq))/(2*nfsq*tmp1);
			dfx[4 + os] += (sc*tcu*(pyDy))/(nfsq*tmp1);
			dfx[5 + os] += (sc*tsq*(pyDy))/(nfsq*tmp1);
			dfx[6 + os] += (sc*t*(pyDy))/(nfsq*tmp1);
			dfx[7 + os] += (sc*(2*dy - 2*Dy + 2*cy*t + 2*ay*tcu + 2*by*tsq))/(2*nfsq*tmp1);
			
			}
        }
        
        /* Where does this point belong to? Spline=?, piece=? */
        
    }
    
    
    
    //}
    
    
    
}