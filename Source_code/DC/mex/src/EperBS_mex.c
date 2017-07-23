#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  
  /* Declarations */
  const mxArray *coefsData, *parData, *spInfoData, *knotsData, *indexData, *imOnGPData;
  
  double *coefs, *params, *spInfo, *knots, *index, *imOnGP;
  int N, splpieces, id, vecLength, sn;
  double k, ksqnfsq, frR, nf, nfsq, T;
  double splstart, splend, offsetjump, splos, tr, t, tsq, tcu, obj;
  double ax, bx, cx, dx, ay, by, cy, dy, px, py;
  double dl, dt, dr, db, ol, ot, or, ob;
  double x0, x1, y0, y1, x0x1sq, y0y1sq;
  double derx, dery, derxl, derxt, derxr, derxb, deryl, deryt, deryr, deryb;
  double *dderl, *ddert, *dderr, *dderb;
  int i, ind, cnt, os;
  int knoffset, indexoffset;
  double tmp1, tmp2, tmp3, tmp4, tmp5;
  double px_, py_, pxs, pys,  sqrtpxspys;
  double _1_ksq, _1_nfsq;
  
  int indtmp;
  double brtmp;
  int allcnt;
  
  double *fx, *dfx;
  
  double t1,t2,t3,t4,t5,t6,c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y,curt;
  double cmt1,cmt2,cmt3,cmt4,cmt5,cmt6;
  double t14,t24,t25,t35,t34,t36;
  double c2x_cmt1_c1x_cmt4, c2y_cmt1_c1y_cmt4;
  int *svpos;
  
  /* Copy input pointers */
  coefsData = prhs[0];    /* state vec */
  parData = prhs[1];      /* speed, frame rate */
  spInfoData = prhs[2];   /* pieces, starts, ends */
  knotsData = prhs[3];   /* knots */
  indexData = prhs[4];    /* index */
  imOnGPData = prhs[5];   /* tracking area corners */
  
  coefs  = mxGetPr(coefsData);
  params = mxGetPr(parData);
  spInfo = mxGetPr(spInfoData);
  knots = mxGetPr(knotsData);
  index  = mxGetPr(indexData);
  
  imOnGP = mxGetPr(imOnGPData);
  
  /* How many coefficients? */
  vecLength = (int)mxGetM(coefsData);
  
  /* How many splines? */
  N = (int)mxGetM(spInfoData);     
  mxAssert((int)mxGetN(spInfoData)==4, "spInfo must have 4 columns");
  mxAssert((int)mxGetN(imOnGPData)==8, "imOnGP must have 8 values");
  
  
  /* Get parameters */
  nf = params[0];  nfsq=nf*nf; _1_nfsq=1.0/nfsq; 
  k = params[1]; ksqnfsq=k*k*nfsq;
  _1_ksq=1.0/(k*k);
  T = params[2];
  
  
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
  
  svpos = (int*)malloc(8 * sizeof(int));
  dderl = (double*)malloc(8 * sizeof(double));
  ddert = (double*)malloc(8 * sizeof(double));
  dderr = (double*)malloc(8 * sizeof(double));
  dderb = (double*)malloc(8 * sizeof(double));
  
  /* for each target */
  allcnt=0;
  for (id = 0; id < N; id++) {
    
    /* get spline Info */
    os = (int)spInfo[id];
    splstart =  spInfo[id+N];
    splend =    spInfo[id+2*N];
    splpieces = spInfo[id+3*N];
    
    sn=splpieces+3;
    
    // 	mexPrintf("%d\n",sn);
    
    
    /* for all frames */
    cnt=0;
    for (tr = splstart; tr <= splend; tr += 1.)
    {
      
      if ((tr==splstart && tr>1.5) || (tr==splend && tr < T-.5)) { 
	
	
	curt=tr;
	indtmp=(int)index[cnt+indexoffset]-1;
	
	// 	    mexPrintf("%.1f %d %d %d\n",tr,indtmp-2+knoffset,indtmp+3+knoffset,knoffset);
	
	svpos[0]=os+indtmp+0;
	svpos[1]=os+indtmp+1;
	svpos[2]=os+indtmp+2;
	svpos[3]=os+indtmp+3;
	svpos[4]=os+indtmp+0+sn;
	svpos[5]=os+indtmp+1+sn;
	svpos[6]=os+indtmp+2+sn;
	svpos[7]=os+indtmp+3+sn;
	
	indtmp += 3;
	t1=knots[indtmp-2+knoffset];
	t2=knots[indtmp-1+knoffset];
	t3=knots[indtmp+0+knoffset];
	t4=knots[indtmp+1+knoffset];
	t5=knots[indtmp+2+knoffset];
	t6=knots[indtmp+3+knoffset];
	
	c1x=coefs[svpos[0]];c2x=coefs[svpos[1]];c3x=coefs[svpos[2]];c4x=coefs[svpos[3]];
	c1y=coefs[svpos[4]];c2y=coefs[svpos[5]];c3y=coefs[svpos[6]];c4y=coefs[svpos[7]];
	
	
	
	cmt1=curt-t1;cmt2=curt-t2;cmt3=curt-t3;cmt4=curt-t4;cmt5=curt-t5;cmt6=curt-t6;
	t14=t1-t4;        t24=t2-t4;        t25=t2-t5;        t34=t3-t4; t35=t3-t5;        t36=t3-t6;
	c2x_cmt1_c1x_cmt4=(c2x*cmt1 - c1x*cmt4);
	c2y_cmt1_c1y_cmt4=(c2y*cmt1 - c1y*cmt4);
	
	px=-(((((c2x*(curt - t1) - c1x*(curt - t4))*(curt - t4))/(t1 - t4) - ((c3x*(curt - t2) - c2x*(curt - t5))*(curt - t2))/(t2 - t5))*(curt - t4))/(t2 - t4) - ((((c3x*(curt - t2) - c2x*(curt - t5))*(curt - t5))/(t2 - t5) - ((c4x*(curt - t3) - c3x*(curt - t6))*(curt - t3))/(t3 - t6))*(curt - t3))/(t3 - t5))/(t3 - t4);
	py=-(((((c2y*(curt - t1) - c1y*(curt - t4))*(curt - t4))/(t1 - t4) - ((c3y*(curt - t2) - c2y*(curt - t5))*(curt - t2))/(t2 - t5))*(curt - t4))/(t2 - t4) - ((((c3y*(curt - t2) - c2y*(curt - t5))*(curt - t5))/(t2 - t5) - ((c4y*(curt - t3) - c3y*(curt - t6))*(curt - t3))/(t3 - t6))*(curt - t3))/(t3 - t5))/(t3 - t4);
	
	// left
	x0=imOnGP[0];y0=imOnGP[1]; x1=imOnGP[2];y1=imOnGP[3]; 
	x0x1sq=(x1-x0)*(x1-x0); y0y1sq=(y1-y0)*(y1-y0); 
	tmp1=((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0);
	tmp2=(tmp1/(ksqnfsq*(x0x1sq + y0y1sq)) - 1)*(tmp1/(ksqnfsq*(x0x1sq + y0y1sq)) - 1);
	dl = ((y0-y1)*px + (x1-x0)*py + (x0*y1-x1*y0) ) / (sqrt(x0x1sq + y0y1sq) ); dl=dl/nf;
	ol = (k*k/6)*(1-pow(1-pow((dl/k),2.),3.0));
	
	tmp1=((((c2y*cmt1-c1y*cmt4)*cmt4)/t14-((c3y*cmt2-c2y*cmt5)*cmt2)/t25)*cmt4)/t24;
	tmp2=((tmp1-((((c3y*cmt2-c2y*cmt5)*cmt5)/t25-((c4y*cmt3-c3y*cmt6)*cmt3)/t36)*cmt3)/t35)*(x0-x1))/t34;
	tmp3=((((c2x*cmt1-c1x*cmt4)*cmt4)/t14-((c3x*cmt2-c2x*cmt5)*cmt2)/t25)*cmt4)/t24;
	tmp4=((tmp3-((((c3x*cmt2-c2x*cmt5)*cmt5)/t25-((c4x*cmt3-c3x*cmt6)*cmt3)/t36)*cmt3)/t35)*(y0-y1))/t34;
	tmp5=_1_nfsq*pow((_1_ksq*_1_nfsq*pow(x0*y1-x1*y0+tmp2-tmp4,2.0))/(x0x1sq+y0y1sq)-1.0,2.0);
	
	dderl[0] = (tmp5*pow(curt-t4,3.0)*(y0-y1)*(x0*y1-x1*y0+tmp2-tmp4)*6.0)/(k*(x0x1sq+y0y1sq)*t14*t24*t34);
	dderl[1] = (tmp5*(y0-y1)*((((cmt1*cmt4)/t14+(cmt2*cmt5)/t25)*cmt4)/t24+(cmt3*pow(curt-t5,2.0))/(t25*t35))*(x0*y1-x1*y0+tmp2-tmp4)*-6.0)/(k*(x0x1sq+y0y1sq)*t34);
	dderl[2] = (tmp5*(y0-y1)*((((cmt2*cmt5)/t25+(cmt3*cmt6)/t36)*cmt3)/t35+(pow(curt-t2,2.0)*cmt4)/(t24*t25))*(x0*y1-x1*y0+tmp2-tmp4)*6.0)/(k*(x0x1sq+y0y1sq)*t34);
	dderl[3] = (tmp5*pow(curt-t3,3.0)*(y0-y1)*(x0*y1-x1*y0+tmp2-tmp4)*-6.0)/(k*(x0x1sq+y0y1sq)*t34*t35*t36);
	dderl[4] = (tmp5*pow(curt-t4,3.0)*(x0-x1)*(x0*y1-x1*y0+tmp2-tmp4)*-6.0)/(k*(x0x1sq+y0y1sq)*t14*t24*t34);
	dderl[5] = (tmp5*(x0-x1)*((((cmt1*cmt4)/t14+(cmt2*cmt5)/t25)*cmt4)/t24+(cmt3*pow(curt-t5,2.0))/(t25*t35))*(x0*y1-x1*y0+tmp2-tmp4)*6.0)/(k*(x0x1sq+y0y1sq)*t34);
	dderl[6] = (tmp5*(x0-x1)*((((cmt2*cmt5)/t25+(cmt3*cmt6)/t36)*cmt3)/t35+(pow(curt-t2,2.0)*cmt4)/(t24*t25))*(x0*y1-x1*y0+tmp2-tmp4)*-6.0)/(k*(x0x1sq+y0y1sq)*t34);
	dderl[7] = (tmp5*pow(curt-t3,3.0)*(x0-x1)*(x0*y1-x1*y0+tmp2-tmp4)*6.0)/(k*(x0x1sq+y0y1sq)*t34*t35*t36);
	if (abs(dl)>=k) { for (i=0;i<8;i++){dderl[i]=0;} ol=k*k/6; }
	ol=ol*6/k;
	
	
	// top
	x0=imOnGP[2];y0=imOnGP[3]; x1=imOnGP[4];y1=imOnGP[5]; 
	x0x1sq=(x1-x0)*(x1-x0); y0y1sq=(y1-y0)*(y1-y0); 
	tmp1=((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0);
	tmp2=(tmp1/(ksqnfsq*(x0x1sq + y0y1sq)) - 1)*(tmp1/(ksqnfsq*(x0x1sq + y0y1sq)) - 1);
	dt = ((y0-y1)*px + (x1-x0)*py + (x0*y1-x1*y0) ) / (sqrt(x0x1sq + y0y1sq) ); dt=dt/nf;
	ot = (k*k/6)*(1-pow(1-pow((dt/k),2.),3.0));
	tmp1=((((c2y*cmt1-c1y*cmt4)*cmt4)/t14-((c3y*cmt2-c2y*cmt5)*cmt2)/t25)*cmt4)/t24;
	tmp2=((tmp1-((((c3y*cmt2-c2y*cmt5)*cmt5)/t25-((c4y*cmt3-c3y*cmt6)*cmt3)/t36)*cmt3)/t35)*(x0-x1))/t34;
	tmp3=((((c2x*cmt1-c1x*cmt4)*cmt4)/t14-((c3x*cmt2-c2x*cmt5)*cmt2)/t25)*cmt4)/t24;
	tmp4=((tmp3-((((c3x*cmt2-c2x*cmt5)*cmt5)/t25-((c4x*cmt3-c3x*cmt6)*cmt3)/t36)*cmt3)/t35)*(y0-y1))/t34;
	tmp5=_1_nfsq*pow((_1_ksq*_1_nfsq*pow(x0*y1-x1*y0+tmp2-tmp4,2.0))/(x0x1sq+y0y1sq)-1.0,2.0);
	
	ddert[0] = (tmp5*pow(curt-t4,3.0)*(y0-y1)*(x0*y1-x1*y0+tmp2-tmp4)*6.0)/(k*(x0x1sq+y0y1sq)*t14*t24*t34);
	ddert[1] = (tmp5*(y0-y1)*((((cmt1*cmt4)/t14+(cmt2*cmt5)/t25)*cmt4)/t24+(cmt3*pow(curt-t5,2.0))/(t25*t35))*(x0*y1-x1*y0+tmp2-tmp4)*-6.0)/(k*(x0x1sq+y0y1sq)*t34);
	ddert[2] = (tmp5*(y0-y1)*((((cmt2*cmt5)/t25+(cmt3*cmt6)/t36)*cmt3)/t35+(pow(curt-t2,2.0)*cmt4)/(t24*t25))*(x0*y1-x1*y0+tmp2-tmp4)*6.0)/(k*(x0x1sq+y0y1sq)*t34);
	ddert[3] = (tmp5*pow(curt-t3,3.0)*(y0-y1)*(x0*y1-x1*y0+tmp2-tmp4)*-6.0)/(k*(x0x1sq+y0y1sq)*t34*t35*t36);
	ddert[4] = (tmp5*pow(curt-t4,3.0)*(x0-x1)*(x0*y1-x1*y0+tmp2-tmp4)*-6.0)/(k*(x0x1sq+y0y1sq)*t14*t24*t34);
	ddert[5] = (tmp5*(x0-x1)*((((cmt1*cmt4)/t14+(cmt2*cmt5)/t25)*cmt4)/t24+(cmt3*pow(curt-t5,2.0))/(t25*t35))*(x0*y1-x1*y0+tmp2-tmp4)*6.0)/(k*(x0x1sq+y0y1sq)*t34);
	ddert[6] = (tmp5*(x0-x1)*((((cmt2*cmt5)/t25+(cmt3*cmt6)/t36)*cmt3)/t35+(pow(curt-t2,2.0)*cmt4)/(t24*t25))*(x0*y1-x1*y0+tmp2-tmp4)*-6.0)/(k*(x0x1sq+y0y1sq)*t34);
	ddert[7] = (tmp5*pow(curt-t3,3.0)*(x0-x1)*(x0*y1-x1*y0+tmp2-tmp4)*6.0)/(k*(x0x1sq+y0y1sq)*t34*t35*t36);
	if (abs(dt)>=k) { for (i=0;i<8;i++){ddert[i]=0;} ot=k*k/6; }
	ot=ot*6/k;
	
	// right
	x0=imOnGP[4];y0=imOnGP[5]; x1=imOnGP[6];y1=imOnGP[7]; 
	x0x1sq=(x1-x0)*(x1-x0); y0y1sq=(y1-y0)*(y1-y0); 
	tmp1=((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0);
	tmp2=(tmp1/(ksqnfsq*(x0x1sq + y0y1sq)) - 1)*(tmp1/(ksqnfsq*(x0x1sq + y0y1sq)) - 1);
	dr = ((y0-y1)*px + (x1-x0)*py + (x0*y1-x1*y0) ) / (sqrt(x0x1sq + y0y1sq) ); dr=dr/nf;
	or = (k*k/6)*(1-pow(1-pow((dr/k),2.),3.0));
	tmp1=((((c2y*cmt1-c1y*cmt4)*cmt4)/t14-((c3y*cmt2-c2y*cmt5)*cmt2)/t25)*cmt4)/t24;
	tmp2=((tmp1-((((c3y*cmt2-c2y*cmt5)*cmt5)/t25-((c4y*cmt3-c3y*cmt6)*cmt3)/t36)*cmt3)/t35)*(x0-x1))/t34;
	tmp3=((((c2x*cmt1-c1x*cmt4)*cmt4)/t14-((c3x*cmt2-c2x*cmt5)*cmt2)/t25)*cmt4)/t24;
	tmp4=((tmp3-((((c3x*cmt2-c2x*cmt5)*cmt5)/t25-((c4x*cmt3-c3x*cmt6)*cmt3)/t36)*cmt3)/t35)*(y0-y1))/t34;
	tmp5=_1_nfsq*pow((_1_ksq*_1_nfsq*pow(x0*y1-x1*y0+tmp2-tmp4,2.0))/(x0x1sq+y0y1sq)-1.0,2.0);
	
	dderr[0] = (tmp5*pow(curt-t4,3.0)*(y0-y1)*(x0*y1-x1*y0+tmp2-tmp4)*6.0)/(k*(x0x1sq+y0y1sq)*t14*t24*t34);
	dderr[1] = (tmp5*(y0-y1)*((((cmt1*cmt4)/t14+(cmt2*cmt5)/t25)*cmt4)/t24+(cmt3*pow(curt-t5,2.0))/(t25*t35))*(x0*y1-x1*y0+tmp2-tmp4)*-6.0)/(k*(x0x1sq+y0y1sq)*t34);
	dderr[2] = (tmp5*(y0-y1)*((((cmt2*cmt5)/t25+(cmt3*cmt6)/t36)*cmt3)/t35+(pow(curt-t2,2.0)*cmt4)/(t24*t25))*(x0*y1-x1*y0+tmp2-tmp4)*6.0)/(k*(x0x1sq+y0y1sq)*t34);
	dderr[3] = (tmp5*pow(curt-t3,3.0)*(y0-y1)*(x0*y1-x1*y0+tmp2-tmp4)*-6.0)/(k*(x0x1sq+y0y1sq)*t34*t35*t36);
	dderr[4] = (tmp5*pow(curt-t4,3.0)*(x0-x1)*(x0*y1-x1*y0+tmp2-tmp4)*-6.0)/(k*(x0x1sq+y0y1sq)*t14*t24*t34);
	dderr[5] = (tmp5*(x0-x1)*((((cmt1*cmt4)/t14+(cmt2*cmt5)/t25)*cmt4)/t24+(cmt3*pow(curt-t5,2.0))/(t25*t35))*(x0*y1-x1*y0+tmp2-tmp4)*6.0)/(k*(x0x1sq+y0y1sq)*t34);
	dderr[6] = (tmp5*(x0-x1)*((((cmt2*cmt5)/t25+(cmt3*cmt6)/t36)*cmt3)/t35+(pow(curt-t2,2.0)*cmt4)/(t24*t25))*(x0*y1-x1*y0+tmp2-tmp4)*-6.0)/(k*(x0x1sq+y0y1sq)*t34);
	dderr[7] = (tmp5*pow(curt-t3,3.0)*(x0-x1)*(x0*y1-x1*y0+tmp2-tmp4)*6.0)/(k*(x0x1sq+y0y1sq)*t34*t35*t36);
	if (abs(dr)>=k) { for (i=0;i<8;i++){dderr[i]=0;} or=k*k/6; }
	or=or*6/k;
	
	// bottom
	x0=imOnGP[6];y0=imOnGP[7]; x1=imOnGP[0];y1=imOnGP[1]; 
	x0x1sq=(x1-x0)*(x1-x0); y0y1sq=(y1-y0)*(y1-y0); 
	tmp1=((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0);
	tmp2=(tmp1/(ksqnfsq*(x0x1sq + y0y1sq)) - 1)*(tmp1/(ksqnfsq*(x0x1sq + y0y1sq)) - 1);
	db = ((y0-y1)*px + (x1-x0)*py + (x0*y1-x1*y0) ) / (sqrt(x0x1sq + y0y1sq) ); db=db/nf;
	ob = (k*k/6)*(1-pow(1-pow((db/k),2.),3.0));
	
	tmp1=((((c2y*cmt1-c1y*cmt4)*cmt4)/t14-((c3y*cmt2-c2y*cmt5)*cmt2)/t25)*cmt4)/t24;
	tmp2=((tmp1-((((c3y*cmt2-c2y*cmt5)*cmt5)/t25-((c4y*cmt3-c3y*cmt6)*cmt3)/t36)*cmt3)/t35)*(x0-x1))/t34;
	tmp3=((((c2x*cmt1-c1x*cmt4)*cmt4)/t14-((c3x*cmt2-c2x*cmt5)*cmt2)/t25)*cmt4)/t24;
	tmp4=((tmp3-((((c3x*cmt2-c2x*cmt5)*cmt5)/t25-((c4x*cmt3-c3x*cmt6)*cmt3)/t36)*cmt3)/t35)*(y0-y1))/t34;
	tmp5=_1_nfsq*pow((_1_ksq*_1_nfsq*pow(x0*y1-x1*y0+tmp2-tmp4,2.0))/(x0x1sq+y0y1sq)-1.0,2.0);
	
	dderb[0] = (tmp5*pow(curt-t4,3.0)*(y0-y1)*(x0*y1-x1*y0+tmp2-tmp4)*6.0)/(k*(x0x1sq+y0y1sq)*t14*t24*t34);
	dderb[1] = (tmp5*(y0-y1)*((((cmt1*cmt4)/t14+(cmt2*cmt5)/t25)*cmt4)/t24+(cmt3*pow(curt-t5,2.0))/(t25*t35))*(x0*y1-x1*y0+tmp2-tmp4)*-6.0)/(k*(x0x1sq+y0y1sq)*t34);
	dderb[2] = (tmp5*(y0-y1)*((((cmt2*cmt5)/t25+(cmt3*cmt6)/t36)*cmt3)/t35+(pow(curt-t2,2.0)*cmt4)/(t24*t25))*(x0*y1-x1*y0+tmp2-tmp4)*6.0)/(k*(x0x1sq+y0y1sq)*t34);
	dderb[3] = (tmp5*pow(curt-t3,3.0)*(y0-y1)*(x0*y1-x1*y0+tmp2-tmp4)*-6.0)/(k*(x0x1sq+y0y1sq)*t34*t35*t36);
	dderb[4] = (tmp5*pow(curt-t4,3.0)*(x0-x1)*(x0*y1-x1*y0+tmp2-tmp4)*-6.0)/(k*(x0x1sq+y0y1sq)*t14*t24*t34);
	dderb[5] = (tmp5*(x0-x1)*((((cmt1*cmt4)/t14+(cmt2*cmt5)/t25)*cmt4)/t24+(cmt3*pow(curt-t5,2.0))/(t25*t35))*(x0*y1-x1*y0+tmp2-tmp4)*6.0)/(k*(x0x1sq+y0y1sq)*t34);
	dderb[6] = (tmp5*(x0-x1)*((((cmt2*cmt5)/t25+(cmt3*cmt6)/t36)*cmt3)/t35+(pow(curt-t2,2.0)*cmt4)/(t24*t25))*(x0*y1-x1*y0+tmp2-tmp4)*-6.0)/(k*(x0x1sq+y0y1sq)*t34);
	dderb[7] = (tmp5*pow(curt-t3,3.0)*(x0-x1)*(x0*y1-x1*y0+tmp2-tmp4)*6.0)/(k*(x0x1sq+y0y1sq)*t34*t35*t36);
	
	if (abs(db)>=k) { for (i=0;i<8;i++){dderb[i]=0;} ob=k*k/6; }
	ob=ob*6/k;
	
	
	obj=ol*ot*or*ob;
	*fx += obj;
	
	for (i=0;i<8;i++){
	  dfx[svpos[i]] = dderl[i]*ot*or*ob + ol*ddert[i]*or*ob + ol*ot*dderr[i]*ob + ol*ot*or*dderb[i];
	}
	
	
	
      }
      
      cnt++;
    }
    splos += splpieces*8;
    knoffset += (splpieces+7);
    indexoffset += cnt;
  }
  
  
  free(svpos);   
  free(dderl); free(ddert);
  free(dderr); free(dderb);
}