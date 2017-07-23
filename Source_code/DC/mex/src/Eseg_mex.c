#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
     /* Declarations */
     const mxArray *coefsData, *parData, *spInfoData, *breaksData, *indexData;
    
    double *coefs,*params, *spInfo, *breaks, *index;
    int N, splpieces, id, vecLength;
    double ax, bx, cx, dx, ay, by, cy, dy, c2x, c2y, d2x, d2y;
    double splstart, splend;
    double offsetjump, splos, tr, t, t2, obj, T, T2;
    double px, py, pxdx, pydy, tsq, tcu, t2sq, p2x, p2y;
    double extraWt;
    double *fx, *dfx;
    
    int seg, broffset, os;
         int i, ind, cnt;
    int indexoffset;
	int indtmp;
	double brtmp;
	
    /* Copy input pointers */
    coefsData = prhs[0];    /* state vec */    
	parData = prhs[1];      /* seg2 weight */
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
     
     
     /* // Allocate memory and assign output pointer */
     plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
     plhs[1] = mxCreateDoubleMatrix(vecLength, 1, mxREAL);
     
     /* //Get a pointer to the data space in our newly allocated memory  */
     fx = mxGetPr(plhs[0]);
     dfx = mxGetPr(plhs[1]);
     
     *fx=0;
     
     //printf("%f %f\n",speed,frR);
     
     broffset=0;
     os=0;
	 	 
	 extraWt = params[0];
     
    /* for each target */
    for (id = 0; id < N; id++) {
        
        /* get spline Info */
        splpieces = (int)spInfo[id];
//         splstart =  spInfo[id+N];
//         splend =    spInfo[id+2*N];
        
        /* only care if more than one segment */
        for (seg = 0; seg < splpieces-1; seg++){
            
            t = breaks[broffset+1] - breaks[broffset];
//             printf("%i %i %i %i   %f %f %f\n",id, seg,os,broffset,breaks[broffset],breaks[broffset+1],t );
            
             ax=coefs[0+os];bx=coefs[1+os];cx=coefs[2+os];dx=coefs[3+os];
             ay=coefs[4+os];by=coefs[5+os];cy=coefs[6+os];dy=coefs[7+os];
             d2x=coefs[11+os]; d2y=coefs[15+os];                    
            
            
            
            tsq = t*t; tcu=tsq*t;
            
            px = ax*tcu + bx*tsq + cx*t + dx;
            py = ay*tcu + by*tsq + cy*t + dy;
            pxdx=px-d2x; pydy=py-d2y;
            
            
            obj=pxdx*pxdx + pydy*pydy;
            *fx += obj;
        dfx[0+os] = dfx[0+os] + 2*tcu*pxdx;
        dfx[1+os] = dfx[1+os] + 2*tsq*pxdx;
        dfx[2+os] = dfx[2+os] + 2*t*pxdx;
        dfx[3+os] = dfx[3+os] + 2*dx - 2*d2x + 2*t*cx + 2*tcu*ax + 2*tsq*bx;
        dfx[4+os] = dfx[4+os] + 2*tcu*pydy;
        dfx[5+os] = dfx[5+os] + 2*tsq*pydy;
        dfx[6+os] = dfx[6+os] + 2*t*pydy;
        dfx[7+os] = dfx[7+os] + 2*dy - 2*d2y + 2*t*cy + 2*tcu*ay + 2*tsq*by;
        dfx[11+os] = dfx[11+os] + 2*d2x - 2*dx - 2*t*cx - 2*tcu*ax - 2*tsq*bx;
        dfx[15+os] = dfx[15+os] + 2*d2y - 2*dy - 2*t*cy - 2*tcu*ay - 2*tsq*by;            
            
           /* 1st derivative */
            c2x=coefs[10+os]; c2y=coefs[14+os];
            
            px = 3*ax*tsq + 2*bx*t + cx;
            py = 3*ay*tsq + 2*by*t + cy;
            pxdx=px-c2x; pydy=py-c2y;

            obj = pxdx*pxdx + pydy*pydy;
            *fx += obj;
        dfx[0 + os] = dfx[0 + os] + 6*tsq*(cx - c2x + 2*t*bx + 3*tsq*ax);
        dfx[1 + os] = dfx[1 + os] + 4*t*(cx - c2x + 2*t*bx + 3*tsq*ax);
        dfx[2 + os] = dfx[2 + os] + 2*cx - 2*c2x + 4*t*bx + 6*tsq*ax;
        dfx[4 + os] = dfx[4 + os] + 6*tsq*(cy - c2y + 2*t*by + 3*tsq*ay);
        dfx[5 + os] = dfx[5 + os] + 4*t*(cy - c2y + 2*t*by + 3*tsq*ay);
        dfx[6 + os] = dfx[6 + os] + 2*cy - 2*c2y + 4*t*by + 6*tsq*ay;
        dfx[10 + os] = dfx[10 + os] + 2*c2x - 2*cx - 4*t*bx - 6*tsq*ax;
        dfx[14 + os] = dfx[14 + os] + 2*c2y - 2*cy - 4*t*by - 6*tsq*ay;
        
        
            os += 8;
            broffset++;
        }
		
		
		
        os+=8;
        broffset+=2;
//         printf("%i %i\n",id, broffset);
        splos += 8;
//         broffset += (splpieces+1);
//         indexoffset += cnt;
    }
    
     broffset=0;
     os=0;
     
    /* for each target */
	/*
    for (id = 0; id < N; id++) {
        
        /* get spline Info 
        splpieces = (int)spInfo[id];
        for (seg = 0; seg < splpieces; seg++){
            
			if (seg==0) {
				T=breaks[broffset];
				T2=T-1;
				
				ax=coefs[0+os];bx=coefs[1+os];cx=coefs[2+os];
				ay=coefs[4+os];by=coefs[5+os];cy=coefs[6+os];
				
				tsq = T*T;
				t2sq = T2*T2;
				
				px = 3*ax*tsq + 2*bx*T +cx;
				py = 3*ay*tsq + 2*by*T +cy;
				p2x = 3*ax*t2sq + 2*bx*T2 +cx;
				p2y = 3*ay*t2sq + 2*by*T2 +cy;
				
				
				pxdx=px-p2x; pydy=py-p2y;
	            obj = pxdx*pxdx + pydy*pydy;
				*fx += extraWt*obj;
				dfx[0 + os] += extraWt*2*(3*tsq - 3*t2sq)*(2*T*bx - 2*T2*bx + 3*tsq*ax - 3*t2sq*ax);
				dfx[1 + os] += extraWt*2*(2*T - 2*T2)*(2*T*bx - 2*T2*bx + 3*tsq*ax - 3*t2sq*ax);
				dfx[4 + os] += extraWt*2*(3*tsq - 3*t2sq)*(2*T*by - 2*T2*by + 3*tsq*ay - 3*t2sq*ay);
				dfx[5 + os] += extraWt*2*(2*T - 2*T2)*(2*T*by - 2*T2*by + 3*tsq*ay - 3*t2sq*ay);

			}
			if (seg==splpieces-1) {
				T=breaks[broffset+1];
				T2=T+1;
				
				ax=coefs[0+os];bx=coefs[1+os];cx=coefs[2+os];
				ay=coefs[4+os];by=coefs[5+os];cy=coefs[6+os];
				
				tsq = T*T;
				t2sq = T2*T2;
				
				px = 3*ax*tsq + 2*bx*T;
				py = 3*ay*tsq + 2*by*T;
				p2x = 3*ax*t2sq + 2*bx*T2;
				p2y = 3*ay*t2sq + 2*by*T2;
				
				
				pxdx=px-p2x; pydy=py-p2y;
	            obj = pxdx*pxdx + pydy*pydy;
				*fx += extraWt*obj;
				dfx[0 + os] += extraWt*2*(3*tsq - 3*t2sq)*(2*T*bx - 2*T2*bx + 3*tsq*ax - 3*t2sq*ax);
				dfx[1 + os] += extraWt*2*(2*T - 2*T2)*(2*T*bx - 2*T2*bx + 3*tsq*ax - 3*t2sq*ax);
				dfx[4 + os] += extraWt*2*(3*tsq - 3*t2sq)*(2*T*by - 2*T2*by + 3*tsq*ay - 3*t2sq*ay);
				dfx[5 + os] += extraWt*2*(2*T - 2*T2)*(2*T*by - 2*T2*by + 3*tsq*ay - 3*t2sq*ay);

			}
			
            os += 8;
            broffset++;
        }
        splos += 8;
		broffset++;
    }
	*/
	
	offsetjump=8; splos=0;
    broffset = 0;
    indexoffset = 0;
	 
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
            if (tr==splstart) {
             indtmp=(int)index[cnt+indexoffset]-1;
             brtmp=breaks[indtmp + broffset];
                         
              t = tr - brtmp;
			  T=t;
			  T2=T-1;
             
              os = ((int)(index[cnt+indexoffset])-1)*offsetjump + splos;
              ax=coefs[0+os];bx=coefs[1+os];cx=coefs[2+os];
              ay=coefs[4+os];by=coefs[5+os];cy=coefs[6+os];
				tsq = T*T;
				t2sq = T2*T2;
				
				px = 3*ax*tsq + 2*bx*T +cx;
				py = 3*ay*tsq + 2*by*T +cy;
				p2x = 3*ax*t2sq + 2*bx*T2 +cx;
				p2y = 3*ay*t2sq + 2*by*T2 +cy;
				
				
				pxdx=px-p2x; pydy=py-p2y;
	            obj = pxdx*pxdx + pydy*pydy;
				*fx += extraWt*obj;
				dfx[0 + os] += extraWt*2*(3*tsq - 3*t2sq)*(2*T*bx - 2*T2*bx + 3*tsq*ax - 3*t2sq*ax);
				dfx[1 + os] += extraWt*2*(2*T - 2*T2)*(2*T*bx - 2*T2*bx + 3*tsq*ax - 3*t2sq*ax);
				dfx[4 + os] += extraWt*2*(3*tsq - 3*t2sq)*(2*T*by - 2*T2*by + 3*tsq*ay - 3*t2sq*ay);
				dfx[5 + os] += extraWt*2*(2*T - 2*T2)*(2*T*by - 2*T2*by + 3*tsq*ay - 3*t2sq*ay);
				
			}
			if (tr==splend) {
             indtmp=(int)index[cnt+indexoffset]-1;
             brtmp=breaks[indtmp + broffset];
                         
              t = tr - brtmp;
			  T=t;
			  T2=T+1;
             
              os = ((int)(index[cnt+indexoffset])-1)*offsetjump + splos;
              ax=coefs[0+os];bx=coefs[1+os];cx=coefs[2+os];
              ay=coefs[4+os];by=coefs[5+os];cy=coefs[6+os];
				tsq = T*T;
				t2sq = T2*T2;
				
				px = 3*ax*tsq + 2*bx*T +cx;
				py = 3*ay*tsq + 2*by*T +cy;
				p2x = 3*ax*t2sq + 2*bx*T2 +cx;
				p2y = 3*ay*t2sq + 2*by*T2 +cy;
				
				
				pxdx=px-p2x; pydy=py-p2y;
	            obj = pxdx*pxdx + pydy*pydy;
				*fx += extraWt*obj;
				dfx[0 + os] += extraWt*2*(3*tsq - 3*t2sq)*(2*T*bx - 2*T2*bx + 3*tsq*ax - 3*t2sq*ax);
				dfx[1 + os] += extraWt*2*(2*T - 2*T2)*(2*T*bx - 2*T2*bx + 3*tsq*ax - 3*t2sq*ax);
				dfx[4 + os] += extraWt*2*(3*tsq - 3*t2sq)*(2*T*by - 2*T2*by + 3*tsq*ay - 3*t2sq*ay);
				dfx[5 + os] += extraWt*2*(2*T - 2*T2)*(2*T*by - 2*T2*by + 3*tsq*ay - 3*t2sq*ay);			
			}
              cnt++;
        }
        splos += splpieces*8;
        broffset += (splpieces+1);
        indexoffset += cnt;
	}
    
}