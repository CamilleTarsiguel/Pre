#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  
  /* Declarations */
  const mxArray *coefsData, *parData, *spInfoData, *knotsData, *indexData, *ptsData;
  
  double *coefs, *params, *spInfo, *knots, *index, *alldpoints;
  int N, Npts,  id, vecLength, cnt;
  double k, ksq, nf, nfsq;
  double offsetjump, splos, t, obj;
  double ax, bx, cx, dx, ay, by, cy, dy;
  int os, det, df;
  int broffset;
  double px, py, pxDx, pyDy, tsq, tcu;
  double Dx, Dy, sc;
  
  double spt;
  double ptx,pty,pts;
  int ptt, ptl,  cos, cos2, kos;
  
  
  double *fx, *dfx;
  
  /* temp vars */    
  double curt,t1,t2,t3,t4,t5,t6,c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y;
  double cmt1,cmt2,cmt3,cmt4,cmt5,cmt6;
  double t14,t24,t25,t35,t34,t36;
  double c2x_cmt1_c1x_cmt4, c2y_cmt1_c1y_cmt4;
  double _1_nfsq, tmp1, tmp2, tmp3, tmp4, tmp5;
  
  
  int *svpos;
  
  /* Copy input pointers */
  coefsData = prhs[0];    /* state vec */
  parData = prhs[1];      /* speed, frame rate */
  spInfoData = prhs[2];   /* pieces, starts, ends */
  knotsData = prhs[3];   /* knots */
  indexData = prhs[4];    /* index */
  ptsData = prhs[5];      /* all points */
  
  coefs  = mxGetPr(coefsData);
  params = mxGetPr(parData);
  spInfo = mxGetPr(spInfoData);
  knots = mxGetPr(knotsData);
  index  = mxGetPr(indexData);
  alldpoints = mxGetPr(ptsData);
  
  /* How many coefficients? */
  vecLength = (int)mxGetM(coefsData);
  
  /* How many splines? */
  N = (int)mxGetM(spInfoData);
  mxAssert((int)mxGetN(spInfoData)==4, "spInfo must have 4 columns");
  
  /* How many detections? */
  Npts = (int)mxGetN(ptsData);
  
  /* Get parameters */
  df = (int)params[0];
  k = params[2]; ksq=k*k;    nf = params[1]; nfsq=nf*nf; _1_nfsq=1.0/nfsq;
  
  
  /* // Allocate memory and assign output pointer */
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(vecLength, 1, mxREAL);
  
  /* //Get a pointer to the data space in our newly allocated memory  */
  fx = mxGetPr(plhs[0]);
  dfx = mxGetPr(plhs[1]);
  
  offsetjump=8; splos=0;
  broffset = 0;
  
  *fx=0;
  for (id = 0; id < vecLength; id++) {
    dfx[id]=0;
  }
  svpos = (int*)malloc(8 * sizeof(int));
  
  
  /* for each target */
  // for (id = 0; id < N; id++) {
  
  cnt=0;
  for (det = 0; det < Npts; det++) {
    
    ptl=(int)alldpoints[det*9+4];  /* label of detection */
    
    if (ptl > 0) {
      ptx=alldpoints[det*9];         /* x */
      pty=alldpoints[det*9+1];       /* y */
      pts=alldpoints[det*9+2];       /* scores */
      ptt=(int)alldpoints[det*9+3];  /* t of detection */
      spt=alldpoints[det*9+5];  /* t to evaluate spline */
      cos=(int)alldpoints[det*9+6];  /* state vector offset */
      cos2=(int)alldpoints[det*9+7];  /* state vector offset */
      kos=(int)alldpoints[det*9+8];  /* knots vector offset */
      
      //             printf("%.1f %.1f %.1f %i %i %.1f %i %i\n",ptx,pty,pts,ptt,ptl, spt,cos,kos);
      Dx=ptx;Dy=pty; sc=pts;
      
      /* only consider labeled points */
      
      curt=spt;
            
      svpos[0]=cos+0;
      svpos[1]=cos+1;
      svpos[2]=cos+2;
      svpos[3]=cos+3;
      svpos[4]=cos2+0;
      svpos[5]=cos2+1;
      svpos[6]=cos2+2;
      svpos[7]=cos2+3;
      
      t1=knots[kos+0];
      t2=knots[kos+1];
      t3=knots[kos+2];
      t4=knots[kos+3];
      t5=knots[kos+4];
      t6=knots[kos+5];
      
      
      c1x=coefs[svpos[0]];c2x=coefs[svpos[1]];c3x=coefs[svpos[2]];c4x=coefs[svpos[3]];
      c1y=coefs[svpos[4]];c2y=coefs[svpos[5]];c3y=coefs[svpos[6]];c4y=coefs[svpos[7]];
      
      cmt1=curt-t1;cmt2=curt-t2;cmt3=curt-t3;cmt4=curt-t4;cmt5=curt-t5;cmt6=curt-t6;
      t14=t1-t4;        t24=t2-t4;        t25=t2-t5;        t35=t3-t5;     t34=t3-t4;   t36=t3-t6;
      c2x_cmt1_c1x_cmt4=(c2x*cmt1 - c1x*cmt4);
      c2y_cmt1_c1y_cmt4=(c2y*cmt1 - c1y*cmt4);
      
      if (df == 2) {
	
	
	mexPrintf("WARNING! L2 DATA NOT IMPLEMENTED\n");
	
      }
      else if (df == 1 || df == 4) {
	

	tmp1=((((c2x*cmt1-c1x*cmt4)*cmt4)/t14-((c3x*cmt2-c2x*cmt5)*cmt2)/t25)*cmt4)/t24;
	tmp2=((((c3x*cmt2-c2x*cmt5)*cmt5)/t25-((c4x*cmt3-c3x*cmt6)*cmt3)/t36)*cmt3)/t35;
	tmp3=((((c2y*cmt1-c1y*cmt4)*cmt4)/t14-((c3y*cmt2-c2y*cmt5)*cmt2)/t25)*cmt4)/t24;
	tmp4=((((c3y*cmt2-c2y*cmt5)*cmt5)/t25-((c4y*cmt3-c3y*cmt6)*cmt3)/t36)*cmt3)/t35;
	tmp5=1.0/sqrt(k+_1_nfsq*(pow(Dx+(tmp1-tmp2)/t34,2.0)+pow(Dy+(tmp3-tmp4)/t34,2.0)));
	
	/* simple Charbonnier */
    obj = sc/tmp5;
    *fx+=obj;
    dfx[svpos[0]] += -(_1_nfsq*sc*(Dx+(tmp1-tmp2)/t34)*pow(curt-t4,3.0)*tmp5)/(t14*t24*t34);
    dfx[svpos[1]] += (_1_nfsq*sc*(Dx+(tmp1-tmp2)/t34)*tmp5*((((cmt1*cmt4)/t14+(cmt2*cmt5)/t25)*cmt4)/t24+(cmt3*pow(curt-t5,2.0))/(t25*t35)))/t34;
    dfx[svpos[2]] += -(_1_nfsq*sc*(Dx+(tmp1-tmp2)/t34)*tmp5*((((cmt2*cmt5)/t25+(cmt3*cmt6)/t36)*cmt3)/t35+(pow(curt-t2,2.0)*cmt4)/(t24*t25)))/t34;
    dfx[svpos[3]] += (_1_nfsq*sc*(Dx+(tmp1-tmp2)/t34)*pow(curt-t3,3.0)*tmp5)/(t34*t35*t36);
    dfx[svpos[4]] += -(_1_nfsq*sc*(Dy+(tmp3-tmp4)/t34)*pow(curt-t4,3.0)*tmp5)/(t14*t24*t34);
    dfx[svpos[5]] += (_1_nfsq*sc*(Dy+(tmp3-tmp4)/t34)*tmp5*((((cmt1*cmt4)/t14+(cmt2*cmt5)/t25)*cmt4)/t24+(cmt3*pow(curt-t5,2.0))/(t25*t35)))/t34;
    dfx[svpos[6]] += -(_1_nfsq*sc*(Dy+(tmp3-tmp4)/t34)*tmp5*((((cmt2*cmt5)/t25+(cmt3*cmt6)/t36)*cmt3)/t35+(pow(curt-t2,2.0)*cmt4)/(t24*t25)))/t34;
    dfx[svpos[7]] += (_1_nfsq*sc*(Dy+(tmp3-tmp4)/t34)*pow(curt-t3,3.0)*tmp5)/(t34*t35*t36);

	
      }
      cnt++;
      
    }
    
    /* Where does this point belong to? Spline=?, piece=? */
    
  }
  
  
  
  //}
  
  free(svpos); 
  
}
