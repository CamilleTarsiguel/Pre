#include "mex.h"
#include <math.h>

/*
double min(double a, double b) {
	return a < b ? a : b;
}
double max(double a, double b) {
	return a > b ? a : b;
}
*/

int min(int a, int b) {
	return a < b ? a : b;
}
int max(int a, int b) {
	return a > b ? a : b;
}


int ROUND(double number)
{
  return (number >= 0) ? ((int)(number + 0.5)) : ((int)(number - 0.5));
}
	
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
     /* Declarations */
     const mxArray *coefsData, *parData, *spInfoData, *breaksData, *indexData;
    
    double *coefs, *params, *spInfo, *breaks, *index, *brInfo;
    int N, Npieces, spl1pieces, spl2pieces, id1, id2, p1, p2, vecLength;
    double siga, sigb, targetSize;
    double spl1start, spl1end, spl2start, spl2end;
	double offsetjump, splos, tr, t, t1, t2, obj, dist;
	double t1sq, t2sq, t1cu, t2cu;
    double ax, bx, cx, dx, ay, by, cy, dy;
	double a2x, b2x, c2x, d2x, a2y, b2y, c2y, d2y;
    int i, ind, cnt, os1, os2;
    int broffset, indexoffset;
    double tmp1, tmp2, tmp3, expexp;
    double px, py, p2x, p2y, px_p2x, py_p2y, px_p2xsq, py_p2ysq;
	double tot, t1locc, t2locc;
	double s1d, s2d, e1d, e2d;
	int s1,s2,e1,e2;
	double sigscale;
	double br1first, br1last, br2first, br2last;
    
    int indtmp;
    double brtmp;
    int allcnt;
	
	int omit;
	int *ommitting;
	int minstart, maxend, maxlength;
    
    double *fx, *dfx;
     
     /* Copy input pointers */
     coefsData = prhs[0];    /* state vec */
    parData = prhs[1];      /* speed, frame rate */
    spInfoData = prhs[2];   /* pieces, starts, ends */
    breaksData = prhs[4];   /* breaks */
    indexData = prhs[5];    /* index */
    
    coefs  = mxGetPr(coefsData);
    params = mxGetPr(parData);
    spInfo = mxGetPr(spInfoData);
	brInfo = mxGetPr(prhs[3]);
    breaks = mxGetPr(breaksData);
    index  = mxGetPr(indexData); /* pcindex, brindex */
    
    /* How many coefficients? */
    vecLength = (int)mxGetM(coefsData);
    
     /* How many splines? */
     N = (int)mxGetM(spInfoData);     
     mxAssert((int)mxGetN(spInfoData)==3, "spInfo must have 3 columns");
	 mxAssert((int)mxGetN(prhs[3])==2, "brInfo must have 2 columns");
     
     /* Get parameters */
     siga = params[0];     targetSize = params[1]; sigb=(targetSize/2)*siga;
	 sigscale=1.;
     
     
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
     for (i = 0; i < vecLength; i++) {
         dfx[i]=0;
     }
	 Npieces=0;
	 for (p1 = 0; p1<N; p1++) {
		Npieces += (int)spInfo[p1];
	}
     
     //printf("%f %f\n",speed,frR);
     
	 /* what is max length */
	 minstart=spInfo[N]; maxend=spInfo[2*N];
	 for (id1=0; id1<N; id1++) {
		//if (minstart>(int)spInfo[id1+N])
			//minstart=(int)spInfo[id1+N];
		if (maxend<(int)spInfo[id1+2*N])
			maxend=(int)spInfo[id1+2*N];
	 }
	 maxlength=maxend+1;
	 
	  ommitting = (int*)malloc(N*N*maxlength * sizeof(int));
	  memset(ommitting, 0, N*N*maxlength * sizeof(int));
	  // ommitting = (int*)calloc(N*N*maxlength, sizeof(int));
	 // for (p1=0; p1<N*N*maxlength; p1++) {
		// ommitting[p1]=0;		
		// }
	
	// for (p1 = 0; p1 < Npieces; p1++) {
		// id1=(int)index[2*p1];
		// s1 =  (int)spInfo[id1+N];
		// e1 =    (int)spInfo[id1+2*N];		
		// for (p2 = p1+1; p2 < Npieces; p2++) {
			// id2=(int)index[2*p2];
			// s2 =  (int)spInfo[id2+N];
			// e2 =    (int)spInfo[id2+2*N];
			// for (tot = max(s1,s2); tot <= min(e1,e2); tot++){
				// ommitting[N*N*(int)tot + N*id2 + id1]=1;
			// }
		// }
	// }
		
		
    /* for each pair of pieces */
      allcnt=0;
     for (p1 = 0; p1 < Npieces; p1++) {
	 
		id1=(int)index[2*p1];
	 	spl1pieces = (int)spInfo[id1];
		spl1start =  spInfo[id1+N];
		spl1end =    spInfo[id1+2*N];
		br1first = brInfo[id1]; 
		br1last = brInfo[id1+N]; 
		s1d=breaks[(int)index[2*p1+1]]; e1d=breaks[(int)index[2*p1+1]+1];
		
		if (spl1start>s1d) s1d=spl1start;
		if (spl1end<e1d) e1d=spl1end;
		if ((p1>0 && (int)index[2*(p1-1)] != id1) || p1==0)
			if (spl1start<br1first-.5)
				s1d=spl1start;
						
		if ((p1<Npieces-1 && (int)index[2*(p1+1)] != id1) || p1==Npieces-1)
			if (spl1end>br1last+.5)
				e1d=spl1end;
						
		s1=ROUND(s1d); e1=ROUND(e1d);
		os1=(p1)*8;
		ax=coefs[0+os1];bx=coefs[1+os1];cx=coefs[2+os1];dx=coefs[3+os1];
		ay=coefs[4+os1];by=coefs[5+os1];cy=coefs[6+os1];dy=coefs[7+os1];
						
        for (p2 = p1+1; p2 < Npieces; p2++) {
		
			/* what splines are we in? */
			
			id2=(int)index[2*p2];
			
			/* get spline Info */


			spl2pieces = (int)spInfo[id2];
			spl2start =  spInfo[id2+N];
			spl2end =    spInfo[id2+2*N];
			

			br2first = brInfo[id2];
			br2last = brInfo[id2+N];
			

			s2d=breaks[(int)index[2*p2+1]]; e2d=breaks[(int)index[2*p2+1]+1];
			
			
			/*
			s1=max(s1,spl1start); e1=min(e1,spl1end);
			s2=max(s2,spl2start); e2=min(e2,spl2end);
			s1=ROUND(s1); e1=ROUND(e1);
			s2=ROUND(s2); e2=ROUND(e2);
			
			s1=spl1start; s2=spl2start;
			e1=spl1end; e2=spl2end;
			*/
			    
				
				
                if (spl2start>s2d) s2d=spl2start;                
                if (spl2end<e2d) e2d=spl2end;
				

						
				if ((p2>0 && (int)index[2*(p2-1)] != id2) || p2==0)
                 if (spl2start<br2first-.5)
                     s2d=spl2start;

				if ((p2<Npieces-1 && (int)index[2*(p2+1)] != id2) || p2==Npieces-1)
                 if (spl2end>br2last+.5)
                     e2d=spl2end;

				 
			
			s2=ROUND(s2d); e2=ROUND(e2d);
			//s1=s1d; s2=s2d; e1=e1d; e2=e2d;
				// printf("%i %i %i %i %f %f %f %f\n",s1,s2,e1,e2,s1d,s2d,e1d,e2d);
			// if (p1==0 && p2==1){
			// printf("%i %i %i %i %f %f %f %f\n",p1,p2,id1,id2,s1,e1,s2,e2);
			// }
			
			/* for all frames */
			/*
			omit=0;
			if (max(s1,s2) == min(e1,e2))
				omit=1;
				
			if (omit==1)
				if (spl1start==spl2end || spl2start == spl1end)
					omit=0;
					*/
			// printf("%i %i %i %i, %i. %i %i %i %i, %i %i %i %i\n",id1+1,id2+1,p1+1,p2+1,omit,s1,s2,e1,e2,(int)spl1start,(int)spl2start,(int)spl1end,(int)spl2end);

					os2=(p2)*8;

					  a2x=coefs[0+os2];b2x=coefs[1+os2];c2x=coefs[2+os2];d2x=coefs[3+os2];
					  a2y=coefs[4+os2];b2y=coefs[5+os2];c2y=coefs[6+os2];d2y=coefs[7+os2];
					  
			
			if ((s1<=e2) && (e1>=s2) && (id1 != id2)) {
				for (tot = max(s1,s2); tot <= min(e1,e2); tot++){
					//printf("o[%d] %d\n",N*N*(int)tot + N*id2 + id1,ommitting[N*N*(int)tot + N*id2 + id1]);
					if (!ommitting[N*N*(int)tot + N*id2 + id1]) {
					ommitting[N*N*(int)tot + N*id2 + id1]=1;
					ommitting[N*N*(int)tot + N*id1 + id2]=1;
						t1 = tot - breaks[(int)index[2*p1+1]];
						t2 = tot - breaks[(int)index[2*p2+1]];
						


					  t1sq=t1*t1; t1cu=t1sq*t1;
					  t2sq=t2*t2; t2cu=t2sq*t2;
					  
					  px = ax*t1cu + bx*t1sq + cx*t1 + dx; 
					  py = ay*t1cu + by*t1sq + cy*t1 + dy;
					  p2x = a2x*t2cu + b2x*t2sq + c2x*t2 + d2x; 
					  p2y = a2y*t2cu + b2y*t2sq + c2y*t2 + d2y;
					  px_p2x=px-p2x; px_p2xsq = px_p2x*px_p2x;
					  py_p2y=py-p2y; py_p2ysq = py_p2y*py_p2y;
					  dist = sqrt(px_p2xsq + py_p2ysq);	
					  // printf("%i %i %f\n",p1,p2,dist);
					  
					  /* do not bother if too far apart */
					  if (dist>targetSize)
						continue;
					  
					  
					  obj = 1.-1./(1+exp(-siga*dist+sigb));

					   // if (p1==0 && p2==8){
					  
					   // printf("%i %i %f %f %.15f %.15f\n",id1+1, id2+1, t1, t2, dist,obj);
					   // }
					  
					  *fx +=obj;
					  tmp1=(d2x - dx + c2x*t2 - cx*t1 + a2x*t2cu - ax*t1cu + b2x*t2sq - bx*t1sq);
					  tmp2=(d2y - dy + c2y*t2 - cy*t1 + a2y*t2cu - ay*t1cu + b2y*t2sq - by*t1sq);
					  expexp=exp(sigb - siga*dist);
					  tmp3=(expexp + 1)*(expexp + 1);
					  

	dfx[0 + os1] += (siga*t1cu*expexp*tmp1)/(tmp3*dist);
	dfx[1 + os1] += (siga*t1sq*expexp*tmp1)/(tmp3*dist);
	dfx[2 + os1] += (siga*t1*expexp*tmp1)/(tmp3*dist);
	dfx[3 + os1] += (siga*expexp*(2*d2x - 2*dx + 2*c2x*t2 - 2*cx*t1 + 2*a2x*t2cu - 2*ax*t1cu + 2*b2x*t2sq - 2*bx*t1sq))/(2*tmp3*dist);
	dfx[4 + os1] += (siga*t1cu*expexp*tmp2)/(tmp3*dist);
	dfx[5 + os1] += (siga*t1sq*expexp*tmp2)/(tmp3*dist);
	dfx[6 + os1] += (siga*t1*expexp*tmp2)/(tmp3*dist);
	dfx[7 + os1] += (siga*expexp*(2*d2y - 2*dy + 2*c2y*t2 - 2*cy*t1 + 2*a2y*t2cu - 2*ay*t1cu + 2*b2y*t2sq - 2*by*t1sq))/(2*tmp3*dist);
	dfx[0 + os2] += -(siga*t2cu*expexp*tmp1)/(tmp3*dist);
	dfx[1 + os2] += -(siga*t2sq*expexp*tmp1)/(tmp3*dist);
	dfx[2 + os2] += -(siga*t2*expexp*tmp1)/(tmp3*dist);
	dfx[3 + os2] += -(siga*expexp*(2*d2x - 2*dx + 2*c2x*t2 - 2*cx*t1 + 2*a2x*t2cu - 2*ax*t1cu + 2*b2x*t2sq - 2*bx*t1sq))/(2*tmp3*dist);
	dfx[4 + os2] += -(siga*t2cu*expexp*tmp2)/(tmp3*dist);
	dfx[5 + os2] += -(siga*t2sq*expexp*tmp2)/(tmp3*dist);
	dfx[6 + os2] += -(siga*t2*expexp*tmp2)/(tmp3*dist);
	dfx[7 + os2] += -(siga*expexp*(2*d2y - 2*dy + 2*c2y*t2 - 2*cy*t1 + 2*a2y*t2cu - 2*ay*t1cu + 2*b2y*t2sq - 2*by*t1sq))/(2*tmp3*dist);
					  
					}
				}
				
			}
			
			/*
			for (tot = max(s1,s2); tot <= min(e1,e2); tot++){
									ommitting[N*N*(int)tot + N*id2 + id1]=1;
						ommitting[N*N*(int)tot + N*id1 + id2]=1;
						}
						*/


		}
	}    
		 free(ommitting);
}