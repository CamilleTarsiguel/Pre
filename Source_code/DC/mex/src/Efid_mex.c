#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    /* Declarations */
    const mxArray *coefsData, *parData, *spInfoData, *breaksData, *indexData, *ptsData, *ptInfoData;
    
    double *coefs, *params, *spInfo, *breaks, *index, *alldpoints, *ptInfo;
    int N, splpieces, Npts,  id, vecLength;
    double siga, sigb, nf, nfsq, k, ksq;
    double splstart, splend, offsetjump, splos, t, obj, dist, disthu, tr;
    double ax, bx, cx, dx, ay, by, cy, dy;
    int i, ind, cnt, co, os, det, detrel, detrel2;
    int broffset, indexoffset;
    double tmp1, tmpprod, fac, tmpexp;
    double px, py,pxs, pys, pxDx, pyDy, tsq, tcu;
    double Dx, Dy, sc;
	int loopused;
	int ptctr, ptsfrom, ptsto, seqLength;
	int *detind, *allptstp;
	double *tmpfx, *tmpdfx;
    
	double spt;
    double ptx,pty,pts;
    int ptt, ptl,  cos;
    int indtmp, inttr;
    double brtmp;
    
    double *fx, *dfx;
	
    
    /* Copy input pointers */
    coefsData = prhs[0];    /* state vec */
    parData = prhs[1];      /* speed, frame rate */
    spInfoData = prhs[2];   /* pieces, starts, ends */
    breaksData = prhs[3];   /* breaks */
    indexData = prhs[4];    /* index */
    ptsData = prhs[5];      /* all points */
	ptInfoData = prhs[6];      /* all points */
    
    coefs  = mxGetPr(coefsData);
    params = mxGetPr(parData);
    spInfo = mxGetPr(spInfoData);
    breaks = mxGetPr(breaksData);
    index  = mxGetPr(indexData);
    alldpoints = mxGetPr(ptsData);
	ptInfo = mxGetPr(ptInfoData);
	
    
    /* How many coefficients? */
    vecLength = (int)mxGetM(coefsData);
    
    /* How many splines? */
    N = (int)mxGetM(spInfoData);
    mxAssert((int)mxGetN(spInfoData)==3, "spInfo must have 3 columns");
	
	mxAssert((int)mxGetN(ptInfoData)==2, "ptInfo must have 2 columns");
    
    /* How many detections? */
    Npts = (int)mxGetN(ptsData);
	seqLength = (int)mxGetM(ptInfoData);
    
    /* Get parameters */
	siga = params[0]; sigb=params[1]*siga; k=params[2]; ksq=k*k;
    
    
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
    detind = (int*)malloc(Npts * sizeof (int));
	
	
	allptstp = (int*)malloc(Npts * sizeof(int));
	for (det = 0; det < Npts; det++) {
		allptstp[det]=(int)alldpoints[det*7+3];
	}

    				tmpfx = (double*)malloc(Npts * sizeof(double));
				tmpdfx = (double*)malloc(8*Npts * sizeof(double));
    
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
              ax=coefs[0+os];bx=coefs[1+os];cx=coefs[2+os];dx=coefs[3+os];
              ay=coefs[4+os];by=coefs[5+os];cy=coefs[6+os];dy=coefs[7+os];
              
              tsq=t*t; tcu=tsq*t;
              px=dx + cx*t + bx*tsq + ax*tcu; pxs = px*px;
              py=dy + cy*t + by*tsq + ay*tcu; pys = py*py;
			  
			  // loopused=0;
			  
			  /* first determine, how many points are relevant */
			  ptctr=0;
			  
			  
			  inttr=(int)tr;
			  /* assuming points are sorted accrd. to time */
			  ptsfrom=ptInfo[inttr-1];
			  ptsto=ptInfo[inttr-1 + seqLength];
			  // for (det = ptsfrom; det <= ptsto; det++) {
			  // ptctr++;
			  // }
			  
			  // for (det = 0; det < Npts && allptstp[det] <= inttr; det++) {				
			  
				// if (allptstp[det] == inttr) {
					// detind[ptctr]=det;
					// ptctr++;
				// }		
			  // }
			   
			  ptctr=ptsto-ptsfrom+1;
			// ptctr=0;
			if (ptctr > 0) {

				// tmpfx=1;
				for (detrel = 0; detrel < ptctr; detrel++) {
					// det=detind[detrel];
					// printf("%i\n",detind[detrel]-(detrel+ptsfrom));
					det=detrel+ptsfrom;

					ptx=alldpoints[det*7];         /* x */
					pty=alldpoints[det*7+1];       /* y */

					
					Dx=ptx;Dy=pty;

					pxDx = px-Dx;
					pyDy = py-Dy;
					//dist=sqrt(pxDx*pxDx+pyDy*pyDy);
					disthu=sqrt((pxDx*pxDx + pyDy*pyDy)/ksq + 1);
					
					if (disthu<params[1]*2) {
						
						// tmpexp=exp(sigb - siga*dist);
						tmpexp=exp(sigb - k * siga*(disthu - 1));
				
						obj=1./(tmpexp + 1);
						tmpfx[detrel] = obj;

						
						
						tmp1=(tmpexp + 1); tmp1=tmp1*tmp1;
						
						/*
						tmpdfx[0+detrel*8] = (siga*tcu*tmpexp*pxDx)/(tmp1*dist);
						tmpdfx[1+detrel*8] = (siga*tsq*tmpexp*pxDx)/(tmp1*dist);
						tmpdfx[2+detrel*8] = (siga*t*tmpexp*pxDx)/(tmp1*dist);
						tmpdfx[3+detrel*8] = (siga*tmpexp*(2*dx - 2*Dx + 2*cx*t + 2*ax*tcu + 2*bx*tsq))/(2*tmp1*dist);
						tmpdfx[4+detrel*8] = (siga*tcu*tmpexp*pyDy)/(tmp1*dist);
						tmpdfx[5+detrel*8] = (siga*tsq*tmpexp*pyDy)/(tmp1*dist);
						tmpdfx[6+detrel*8] = (siga*t*tmpexp*pyDy)/(tmp1*dist);
						tmpdfx[7+detrel*8] = (siga*tmpexp*(2*dy - 2*Dy + 2*cy*t + 2*ay*tcu + 2*by*tsq))/(2*tmp1*dist); */
						
						tmpdfx[0+detrel*8] = (k*siga*tcu*tmpexp*pxDx)/(ksq*disthu*tmp1);
						tmpdfx[1+detrel*8] = (k*siga*tsq*tmpexp*pxDx)/(ksq*disthu*tmp1);
						tmpdfx[2+detrel*8] = (k*siga*t*tmpexp*pxDx)/(ksq*disthu*tmp1);
						tmpdfx[3+detrel*8] = (k*siga*tmpexp*(2*dx - 2*Dx + 2*cx*t + 2*ax*tcu + 2*bx*tsq))/(2*ksq*disthu*tmp1);
						tmpdfx[4+detrel*8] = (k*siga*tcu*tmpexp*pyDy)/(ksq*disthu*tmp1);
						tmpdfx[5+detrel*8] = (k*siga*tsq*tmpexp*pyDy)/(ksq*disthu*tmp1);
						tmpdfx[6+detrel*8] = (k*siga*t*tmpexp*pyDy)/(ksq*disthu*tmp1);
						tmpdfx[7+detrel*8] = (k*siga*tmpexp*(2*dy - 2*Dy + 2*cy*t + 2*ay*tcu + 2*by*tsq))/(2*ksq*disthu*tmp1);
						
						
					} else {
						tmpfx[detrel] = 1;
						tmpdfx[0+detrel*8] = 0;
						tmpdfx[1+detrel*8] = 0;
						tmpdfx[2+detrel*8] = 0;
						tmpdfx[3+detrel*8] = 0;
						tmpdfx[4+detrel*8] = 0;
						tmpdfx[5+detrel*8] = 0;
						tmpdfx[6+detrel*8] = 0;
						tmpdfx[7+detrel*8] = 0;
					}
				}
				
				tmpprod=1.;
				for (detrel = 0; detrel < ptctr; detrel++) {
					tmpprod *= tmpfx[detrel];
				}
				// if (id == 1) {
					// for (detrel = 0; detrel < ptctr; detrel++) {
						// printf("%f ",tmpfx[detrel]);
					// }
					// printf("%f %f\n",tr,tmpprod);
				// }
				// printf("%f %f\n",tr,tmpprod);
				*fx+=tmpprod;
				
				/* derivative */
				for (co=0; co<8; co++) {
					fac=0;
					for (detrel = 0; detrel < ptctr; detrel++) {
					
						tmpprod=1.;
						for (detrel2 = 0; detrel2 < ptctr; detrel2++) {
							if (detrel != detrel2)
								tmpprod *= tmpfx[detrel2];
						}
						
						fac += tmpdfx[co + detrel*8] * tmpprod;
					}
					dfx[co + os] += fac;
					
				}
				
				
				

			} else {
				*fx += 1;
			}
			cnt++;
		}
		// printf("%i %f\n",id+1,*fx);
		splos += splpieces*8;
        broffset += (splpieces+1);
        indexoffset += cnt;
    
    
    
    }
	free(detind);
	free(allptstp);
    				free(tmpfx);
				free(tmpdfx);
    
    
}