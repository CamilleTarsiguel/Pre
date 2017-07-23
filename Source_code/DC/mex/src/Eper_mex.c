#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
     /* Declarations */
     const mxArray *coefsData, *parData, *spInfoData, *breaksData, *indexData, *imOnGPData;
    
    double *coefs, *params, *spInfo, *breaks, *index, *imOnGP;
    int N, splpieces, id, vecLength;
    double k, ksqnfsq, frR, nf, nfsq, T;
    double splstart, splend, offsetjump, splos, tr, t, tsq, tcu, obj;
    double ax, bx, cx, dx, ay, by, cy, dy, px, py;
	double dl, dt, dr, db, ol, ot, or, ob;
	double x0, x1, y0, y1, x0x1sq, y0y1sq;
	double derx, dery, derxl, derxt, derxr, derxb, deryl, deryt, deryr, deryb;
	double *dderl, *ddert, *dderr, *dderb;
    int i, ind, cnt, os;
    int broffset, indexoffset;
    double tmp1, tmp2, tmp3;
    double px_, py_, pxs, pys,  sqrtpxspys;
    
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
	imOnGPData = prhs[5];   /* tracking area corners */
    
    coefs  = mxGetPr(coefsData);
    params = mxGetPr(parData);
    spInfo = mxGetPr(spInfoData);
    breaks = mxGetPr(breaksData);
    index  = mxGetPr(indexData);
	imOnGP = mxGetPr(imOnGPData);
    
    /* How many coefficients? */
    vecLength = (int)mxGetM(coefsData);
    
     /* How many splines? */
     N = (int)mxGetM(spInfoData);     
     mxAssert((int)mxGetN(spInfoData)==3, "spInfo must have 3 columns");
     mxAssert((int)mxGetN(imOnGPData)==8, "imOnGP must have 8 values");
     
     /* Get parameters */
     nf = params[0];  nfsq=nf*nf;
	 k = params[1]; ksqnfsq=k*k*nfsq;
	 T = params[2];
     

     /* // Allocate memory and assign output pointer */
      plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
      plhs[1] = mxCreateDoubleMatrix(vecLength, 1, mxREAL);
      //plhs[2] = mxCreateDoubleMatrix(116, 1, mxREAL);
     
     /* //Get a pointer to the data space in our newly allocated memory  */
      fx = mxGetPr(plhs[0]);
      dfx = mxGetPr(plhs[1]);
      //dfs = mxGetPr(plhs[2]);
     
     offsetjump=8; splos=0;
     broffset = 0;
     indexoffset = 0;
     
      *fx=0;
     for (id = 0; id < vecLength; id++) {
         dfx[id]=0;
     }
     
	 dderl = (double*)malloc(8 * sizeof(double));
	 ddert = (double*)malloc(8 * sizeof(double));
	 dderr = (double*)malloc(8 * sizeof(double));
	 dderb = (double*)malloc(8 * sizeof(double));
				 
     // printf("%f %f %f\n",nf, speed,frR);
     
    /* for each target */
      allcnt=0;
     for (id = 0; id < N; id++) {
        
        /* get spline Info */
        splpieces = (int)spInfo[id];
        splstart =  spInfo[id+N];
        splend =    spInfo[id+2*N];
        
		
        /* first frame */
        cnt=0;
		if (splstart > 1.5) {

			indtmp=(int)index[cnt+indexoffset]-1;
			brtmp=breaks[indtmp + broffset];
			
			t = splstart - brtmp;

			os = ((int)(index[cnt+indexoffset])-1)*offsetjump + splos;
			ax=coefs[0+os];bx=coefs[1+os];cx=coefs[2+os];dx=coefs[3+os];
			ay=coefs[4+os];by=coefs[5+os];cy=coefs[6+os];dy=coefs[7+os];

			tsq=t*t; tcu=tsq*t;
			px=dx + cx*t + bx*tsq + ax*tcu;
			py=dy + cy*t + by*tsq + ay*tcu;

			// left
			
			x0=imOnGP[0];y0=imOnGP[1]; x1=imOnGP[2];y1=imOnGP[3]; 
			x0x1sq=(x1-x0)*(x1-x0); y0y1sq=(y1-y0)*(y1-y0); 
			tmp1=(x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1))*(x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1));
			tmp2=(tmp1/(ksqnfsq*(x0x1sq + y0y1sq)) - 1)*(tmp1/(ksqnfsq*(x0x1sq + y0y1sq)) - 1);
			dl = ((y0-y1)*px + (x1-x0)*py + (x0*y1-x1*y0) ) / (sqrt(x0x1sq + y0y1sq) ); dl=dl/nf;
			ol = (k*k/6)*(1-pow(1-pow((dl/k),2.),3.0));
			derxl=  (6*tmp2*(y0 - y1)*(x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1)))/(k*nfsq*(x0x1sq + y0y1sq));
			deryl= -(6*tmp2*(x0 - x1)*(x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1)))/(k*nfsq*(x0x1sq + y0y1sq));
			if (abs(dl)>=k) { derxl=0; deryl=0; ol=k*k/6; }
			ol=ol*6/k;
			
			
			// top
			x0=imOnGP[2];y0=imOnGP[3]; x1=imOnGP[4];y1=imOnGP[5];
			x0x1sq=(x1-x0)*(x1-x0); y0y1sq=(y1-y0)*(y1-y0);
			tmp1=(x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1))*(x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1));
			tmp2=(tmp1/(ksqnfsq*(x0x1sq + y0y1sq)) - 1)*(tmp1/(ksqnfsq*(x0x1sq + y0y1sq)) - 1);
			dt = ((y0-y1)*px + (x1-x0)*py + (x0*y1-x1*y0) ) / (sqrt(x0x1sq + y0y1sq) ); dt=dt/nf;
			ot = (k*k/6)*(1-pow(1-pow((dt/k),2.),3.0));
			derxt=  (6*tmp2*(y0 - y1)*(x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1)))/(k*nfsq*(x0x1sq + y0y1sq));
			deryt= -(6*tmp2*(x0 - x1)*(x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1)))/(k*nfsq*(x0x1sq + y0y1sq));
			if (abs(dt)>=k) { derxt=0; deryt=0; ot=k*k/6;  }
			ot=ot*6/k;
			
			// right
			x0=imOnGP[4];y0=imOnGP[5]; x1=imOnGP[6];y1=imOnGP[7];
			x0x1sq=(x1-x0)*(x1-x0); y0y1sq=(y1-y0)*(y1-y0);
			tmp1=(x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1))*(x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1));
			tmp2=(tmp1/(ksqnfsq*(x0x1sq + y0y1sq)) - 1)*(tmp1/(ksqnfsq*(x0x1sq + y0y1sq)) - 1);
			dr = ((y0-y1)*px + (x1-x0)*py + (x0*y1-x1*y0) ) / (sqrt(x0x1sq + y0y1sq) ); dr=dr/nf;
			or = (k*k/6)*(1-pow(1-pow((dr/k),2.),3.0));
			derxr=  (6*tmp2*(y0 - y1)*(x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1)))/(k*nfsq*(x0x1sq + y0y1sq));
			deryr= -(6*tmp2*(x0 - x1)*(x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1)))/(k*nfsq*(x0x1sq + y0y1sq));
			if (abs(dr)>=k) { derxr=0; deryr=0; or=k*k/6;  }
			or=or*6/k;
			
			// bottom
			x0=imOnGP[6];y0=imOnGP[7]; x1=imOnGP[0];y1=imOnGP[1];
			x0x1sq=(x1-x0)*(x1-x0); y0y1sq=(y1-y0)*(y1-y0);
			tmp1=(x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1))*(x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1));
			tmp2=(tmp1/(ksqnfsq*(x0x1sq + y0y1sq)) - 1)*(tmp1/(ksqnfsq*(x0x1sq + y0y1sq)) - 1);
			db = ((y0-y1)*px + (x1-x0)*py + (x0*y1-x1*y0) ) / (sqrt(x0x1sq + y0y1sq) ); db=db/nf;
			ob = (k*k/6)*(1-pow(1-pow((db/k),2.),3.0));
			derxb=  (6*tmp2*(y0 - y1)*(x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1)))/(k*nfsq*(x0x1sq + y0y1sq));
			deryb= -(6*tmp2*(x0 - x1)*(x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1)))/(k*nfsq*(x0x1sq + y0y1sq));
			if (abs(db)>=k) { derxb=0; deryb=0; ob=k*k/6;  }
			ob=ob*6/k;
			

			obj=ol*ot*or*ob;
			derx=derxl*ot*or*ob + ol*derxt*or*ob + ol*ot*derxr*ob + ol*ot*or*derxb;
			dery=deryl*ot*or*ob + ol*deryt*or*ob + ol*ot*deryr*ob + ol*ot*or*deryb;
			
			*fx+=obj;
			dfx[3 + os] = derx;
			dfx[7 + os] = dery;	
			
		}
    

		/* persistence end */
		cnt = splend-splstart;


		if (splend < T-.5) {

			indtmp=(int)index[cnt+indexoffset]-1;
			brtmp=breaks[indtmp + broffset];
			
			t = splend - brtmp;

			os = ((int)(index[cnt+indexoffset])-1)*offsetjump + splos;
			ax=coefs[0+os];bx=coefs[1+os];cx=coefs[2+os];dx=coefs[3+os];
			ay=coefs[4+os];by=coefs[5+os];cy=coefs[6+os];dy=coefs[7+os];

			// printf("%f\n",t);
			tsq=t*t; tcu=tsq*t;
			px=dx + cx*t + bx*tsq + ax*tcu;
			py=dy + cy*t + by*tsq + ay*tcu;
			
			// left
			x0=imOnGP[0];y0=imOnGP[1]; x1=imOnGP[2];y1=imOnGP[3]; 
			x0x1sq=(x1-x0)*(x1-x0); y0y1sq=(y1-y0)*(y1-y0); 
			tmp1=((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0);
			tmp2=(tmp1/(ksqnfsq*(x0x1sq + y0y1sq)) - 1)*(tmp1/(ksqnfsq*(x0x1sq + y0y1sq)) - 1);
			dl = ((y0-y1)*px + (x1-x0)*py + (x0*y1-x1*y0) ) / (sqrt(x0x1sq + y0y1sq) ); dl=dl/nf;
			ol = (k*k/6)*(1-pow(1-pow((dl/k),2.),3.0));
			
			dderl[0] = -(6*tcu*(y0 - y1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));
			dderl[1] = -(6*tsq*(y0 - y1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));
			dderl[2] = -(6*t*(y0 - y1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));
			dderl[3] = -(6*(y0 - y1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));
			dderl[4] = (6*tcu*(x0 - x1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));
			dderl[5] = (6*tsq*(x0 - x1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));
			dderl[6] = (6*t*(x0 - x1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));
			dderl[7] = (6*(x0 - x1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));       
			if (abs(dl)>=k) { for (i=0;i<8;i++){dderl[i]=0;} ol=k*k/6; }
			ol=ol*6/k;

			// top
			x0=imOnGP[2];y0=imOnGP[3]; x1=imOnGP[4];y1=imOnGP[5]; 
			x0x1sq=(x1-x0)*(x1-x0); y0y1sq=(y1-y0)*(y1-y0); 
			tmp1=((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0);
			tmp2=(tmp1/(ksqnfsq*(x0x1sq + y0y1sq)) - 1)*(tmp1/(ksqnfsq*(x0x1sq + y0y1sq)) - 1);
			dt = ((y0-y1)*px + (x1-x0)*py + (x0*y1-x1*y0) ) / (sqrt(x0x1sq + y0y1sq) ); dt=dt/nf;
			ot = (k*k/6)*(1-pow(1-pow((dt/k),2.),3.0));
			ddert[0] = -(6*tcu*(y0 - y1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));
			ddert[1] = -(6*tsq*(y0 - y1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));
			ddert[2] = -(6*t*(y0 - y1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));
			ddert[3] = -(6*(y0 - y1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));
			ddert[4] = (6*tcu*(x0 - x1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));
			ddert[5] = (6*tsq*(x0 - x1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));
			ddert[6] = (6*t*(x0 - x1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));
			ddert[7] = (6*(x0 - x1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));       
			if (abs(dt)>=k) { for (i=0;i<8;i++){ddert[i]=0;} ot=k*k/6; }
			ot=ot*6/k;
        
			// right
			x0=imOnGP[4];y0=imOnGP[5]; x1=imOnGP[6];y1=imOnGP[7]; 
			x0x1sq=(x1-x0)*(x1-x0); y0y1sq=(y1-y0)*(y1-y0); 
			tmp1=((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0);
			tmp2=(tmp1/(ksqnfsq*(x0x1sq + y0y1sq)) - 1)*(tmp1/(ksqnfsq*(x0x1sq + y0y1sq)) - 1);
			dr = ((y0-y1)*px + (x1-x0)*py + (x0*y1-x1*y0) ) / (sqrt(x0x1sq + y0y1sq) ); dr=dr/nf;
			or = (k*k/6)*(1-pow(1-pow((dr/k),2.),3.0));
			dderr[0] = -(6*tcu*(y0 - y1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));
			dderr[1] = -(6*tsq*(y0 - y1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));
			dderr[2] = -(6*t*(y0 - y1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));
			dderr[3] = -(6*(y0 - y1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));
			dderr[4] = (6*tcu*(x0 - x1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));
			dderr[5] = (6*tsq*(x0 - x1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));
			dderr[6] = (6*t*(x0 - x1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));
			dderr[7] = (6*(x0 - x1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));       
			if (abs(dr)>=k) { for (i=0;i<8;i++){dderr[i]=0;} or=k*k/6; }
			or=or*6/k;
        
			// bottom
			x0=imOnGP[6];y0=imOnGP[7]; x1=imOnGP[0];y1=imOnGP[1]; 
			x0x1sq=(x1-x0)*(x1-x0); y0y1sq=(y1-y0)*(y1-y0); 
			tmp1=((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0);
			tmp2=(tmp1/(ksqnfsq*(x0x1sq + y0y1sq)) - 1)*(tmp1/(ksqnfsq*(x0x1sq + y0y1sq)) - 1);
			db = ((y0-y1)*px + (x1-x0)*py + (x0*y1-x1*y0) ) / (sqrt(x0x1sq + y0y1sq) ); db=db/nf;
			ob = (k*k/6)*(1-pow(1-pow((db/k),2.),3.0));
			
			dderb[0] = -(6*tcu*(y0 - y1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));
			dderb[1] = -(6*tsq*(y0 - y1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));
			dderb[2] = -(6*t*(y0 - y1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));
			dderb[3] = -(6*(y0 - y1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));
			dderb[4] = (6*tcu*(x0 - x1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));
			dderb[5] = (6*tsq*(x0 - x1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));
			dderb[6] = (6*t*(x0 - x1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));
			dderb[7] = (6*(x0 - x1)*tmp2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nfsq*(x0x1sq + y0y1sq));       
			if (abs(db)>=k) { for (i=0;i<8;i++){dderb[i]=0;} ob=k*k/6; }
			ob=ob*6/k;
        
        
			obj=ol*ot*or*ob;
			*fx += obj;
        
			for (i=0;i<8;i++){
				dfx[os+i] =  dderl[i]*ot*or*ob + ol*ddert[i]*or*ob + ol*ot*dderr[i]*ob + ol*ot*or*dderb[i];
			}
		}

		
		cnt = splend-splstart+1;
			  
        splos += splpieces*8;
        broffset += (splpieces+1);
        indexoffset += cnt;
    }
    
    free(dderl); free(ddert);
	free(dderr); free(dderb);
    
}