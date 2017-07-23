#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* compute angular velocity energy with b-splines */
    
    /* Declarations */
    const mxArray *coefsData, *parData, *spInfoData, *knotsData, *indexData;
    
    double *coefs, *params, *spInfo, *knots, *index;
    int N, splpieces, id, vecLength, sn;
    double speed, frR, nf;
    double splstart, splend, offsetjump, splos, tr, curt, t, obj;
    int i, ind, cnt, os;
    int knoffset, indexoffset;
    double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;
    double tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, tmp17, tmp18, tmp19;
    double tmp20, tmp21, tmp22, tmp23, tmp24, tmp25, tmp26, tmp27, tmp28, tmp29;
    double tmp30, tmp31, tmp32, tmp33, tmp34, tmp35, tmp36, tmp37, tmp38, tmp39;
    double px, py, pxs, pys, tsq,  sqrtpxspys, frRsq, twofrRsq, pxspyssq, pxspyscu;
    double px_, py_, px__, py__, px_spy_ssq, k,ksq,kqu,_1_ksq, epsil, hory;
    
    int indtmp;
    double brtmp;
    int allcnt;
    
    double *fx, *dfx;
    
    double t1,t2,t3,t4,t5,t6,c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y;
    double cmt1,cmt2,cmt3,cmt4,cmt5,cmt6;
    double t14,t24,t25,t35,t36,t34;
    double c2x_cmt1_c1x_cmt4, c2y_cmt1_c1y_cmt4;
    double t34sq, cmt4sq, cmt3sq;
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
    
    
    /* Get parameters */
    speed = params[0];     frR = params[1]; frRsq=frR*frR; twofrRsq=2*frRsq;
    k=params[2]; ksq=k*k; kqu=k*k*k*k; _1_ksq=1.0/ksq;
    epsil=.1;
    hory=params[3];
    
    
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
    
    
    /* for each target */
    allcnt=0;
    for (id = 0; id < N; id++) {
        
        /* get spline Info */
        os = (int)spInfo[id];
        splstart =  spInfo[id+N];
        splend =    spInfo[id+2*N];
        splpieces = spInfo[id+3*N];
        
        // B-splines have one coefficient per segment + 3
        sn=splpieces+3;
        
// 	mexPrintf("%d\n",sn);
        
        
        /* for all frames */
        cnt=0;
        for (tr = splstart; tr <= splend; tr += 1.)
        {
            curt=tr; // current frame to evaluate
            
            indtmp=(int)index[cnt+indexoffset]-1; // which segment are we in?
            
            // where in the state vector are we? x: 0-3,  y: 4-7
            svpos[0]=os+indtmp+0;
            svpos[1]=os+indtmp+1;
            svpos[2]=os+indtmp+2;
            svpos[3]=os+indtmp+3;
            svpos[4]=os+indtmp+0+sn;
            svpos[5]=os+indtmp+1+sn;
            svpos[6]=os+indtmp+2+sn;
            svpos[7]=os+indtmp+3+sn;
            
            // which knots do we use?
            indtmp += 3;
            t1=knots[indtmp-2+knoffset];
            t2=knots[indtmp-1+knoffset];
            t3=knots[indtmp+0+knoffset];
            t4=knots[indtmp+1+knoffset];
            t5=knots[indtmp+2+knoffset];
            t6=knots[indtmp+3+knoffset];
            
            // get the necessary coefficients
            c1x=coefs[svpos[0]];c2x=coefs[svpos[1]];c3x=coefs[svpos[2]];c4x=coefs[svpos[3]];
            c1y=coefs[svpos[4]];c2y=coefs[svpos[5]];c3y=coefs[svpos[6]];c4y=coefs[svpos[7]];
            
            
            
            cmt1=curt-t1; cmt2=curt-t2; cmt3=curt-t3; cmt4=curt-t4; cmt5=curt-t5; cmt6=curt-t6;
            t14=t1-t4;        t24=t2-t4;        t25=t2-t5;        t35=t3-t5;        t36=t3-t6; t34=t3-4;
            c2x_cmt1_c1x_cmt4=(c2x*cmt1 - c1x*cmt4);
            c2y_cmt1_c1y_cmt4=(c2y*cmt1 - c1y*cmt4);
            t34sq=pow(t34,2.0);
            cmt4sq=pow(cmt4,2.0);
            cmt3sq=pow(cmt3,2.0);
            tmp1=((c2y*cmt1*2.0-c1y*cmt4*2.0)/t14-(c3y*cmt2*2.0-c2y*cmt5*2.0)/t25-((c1y-c2y)*cmt4*2.0)/t14+((c2y-c3y)*cmt2*2.0)/t25)/t24;
            tmp2=((c3y*cmt2*2.0-c2y*cmt5*2.0)/t25-(c4y*cmt3*2.0-c3y*cmt6*2.0)/t36-((c2y-c3y)*cmt5*2.0)/t25+((c3y-c4y)*cmt3*2.0)/t36)/t35;
            tmp3=(cmt4*((c1y*2.0-c2y*2.0)/t14-(c2y*2.0-c3y*2.0)/t25))/t24;
            tmp4=(cmt3*((c2y*2.0-c3y*2.0)/t25-(c3y*2.0-c4y*2.0)/t36))/t35;
            tmp5=(((c2y*cmt1-c1y*cmt4)*cmt4)/t14-((c3y*cmt2-c2y*cmt5)*cmt2)/t25)/t24;
            tmp6=(cmt3*((c3y*cmt2-c2y*cmt5)/t25-(c4y*cmt3-c3y*cmt6)/t36-((c2y-c3y)*cmt5)/t25+((c3y-c4y)*cmt3)/t36))/t35;
            tmp7=(((c3y*cmt2-c2y*cmt5)*cmt5)/t25-((c4y*cmt3-c3y*cmt6)*cmt3)/t36)/t35;
            tmp8=(((c3x*cmt2-c2x*cmt5)*cmt5)/t25-((c4x*cmt3-c3x*cmt6)*cmt3)/t36)/t35;
            tmp9=(((cmt1*cmt4)/t14+(cmt2*cmt5)/t25)/t24+(cmt4*(cmt1/t14+cmt2/t25+cmt4/t14+cmt5/t25))/t24+pow(curt-t5,2.0)/(t25*t35)+(cmt3*cmt5*2.0)/(t25*t35));
            tmp10=(((c2x*cmt1-c1x*cmt4)*cmt4)/t14-((c3x*cmt2-c2x*cmt5)*cmt2)/t25)/t24;
            tmp11=((c2x*cmt1*2.0-c1x*cmt4*2.0)/t14-(c3x*cmt2*2.0-c2x*cmt5*2.0)/t25-((c1x-c2x)*cmt4*2.0)/t14+((c2x-c3x)*cmt2*2.0)/t25)/t24;
            tmp12=(cmt4*((c1x*2.0-c2x*2.0)/t14-(c2x*2.0-c3x*2.0)/t25))/t24;
            tmp13=((((c3y*cmt2-c2y*cmt5)*cmt5)/t25-((c4y*cmt3-c3y*cmt6)*cmt3)/t36)*cmt3)/t35;
            tmp14=((c3x*cmt2*2.0-c2x*cmt5*2.0)/t25-(c4x*cmt3*2.0-c3x*cmt6*2.0)/t36-((c2x-c3x)*cmt5*2.0)/t25+((c3x-c4x)*cmt3*2.0)/t36)/t35;
            tmp15=(cmt3*((c2x*2.0-c3x*2.0)/t25-(c3x*2.0-c4x*2.0)/t36))/t35;
            tmp16=((curt*2.0-t1*2.0)/t14+(curt*2.0-t2*2.0)/t25+(curt*2.0-t4*2.0)/t14+(curt*2.0-t5*2.0)/t25)/t24;
            tmp17=((curt*2.0-t2*2.0)/t25+(curt*2.0-t3*2.0)/t36+(curt*2.0-t5*2.0)/t25+(curt*2.0-t6*2.0)/t36)/t35;
            tmp18=(cmt4*((c2y*cmt1-c1y*cmt4)/t14-(c3y*cmt2-c2y*cmt5)/t25-((c1y-c2y)*cmt4)/t14+((c2y-c3y)*cmt2)/t25))/t24;
            tmp19=(cmt4*((c2x*cmt1-c1x*cmt4)/t14-(c3x*cmt2-c2x*cmt5)/t25-((c1x-c2x)*cmt4)/t14+((c2x-c3x)*cmt2)/t25))/t24;
            tmp20=(cmt3*((c3x*cmt2-c2x*cmt5)/t25-(c4x*cmt3-c3x*cmt6)/t36-((c2x-c3x)*cmt5)/t25+((c3x-c4x)*cmt3)/t36))/t35;
            tmp21=(((((c2y*cmt1-c1y*cmt4)*cmt4)/t14-((c3y*cmt2-c2y*cmt5)*cmt2)/t25)*cmt4)/t24-tmp13)/t34;
            tmp22=(tmp1-tmp2-tmp3+tmp4);
            tmp23=(tmp11-tmp14-tmp12+tmp15);
            tmp24=(tmp5-tmp7+tmp18-tmp6);
            tmp25=(tmp10-tmp8+tmp19-tmp20);
            tmp26=pow(tmp5-tmp7+tmp18-tmp6,2.0);
            tmp27=pow(tmp10-tmp8+tmp19-tmp20,2.0);
            tmp28=frRsq*(epsil+pow(1.0/t34sq*tmp23*tmp24-1.0/t34sq*tmp22*tmp25,2.0));
            tmp29=1.0/pow(epsil+1.0/t34sq*tmp27+1.0/t34sq*tmp26,2.0)*2.0;
            tmp30=(tmp28*1.0/pow(epsil+1.0/t34sq*tmp27+1.0/t34sq*tmp26,2.0)+1.0);
            tmp31=(cmt3*(cmt2/t25+cmt3/t36+cmt5/t25+cmt6/t36))/t35;
            tmp32=(((cmt2*cmt5)/t25+(cmt3*cmt6)/t36)/t35+tmp31+pow(curt-t2,2.0)/(t24*t25)+(cmt2*cmt4*2.0)/(t24*t25));
            tmp33=(tmp28*1.0/pow(epsil+1.0/t34sq*tmp27+1.0/t34sq*tmp26,2.0));
            tmp34=1.0/pow(epsil+1.0/t34sq*tmp27+1.0/t34sq*tmp26,3.0);
            
            // if we need to scale according to perspective (ETH sequences)
            if (k>0){
                obj= -log(1.0/((frRsq*(epsil+pow(1.0/t34sq*(((c2x*cmt1*2.0-c1x*cmt4*2.0)/t14-(c3x*(curt-t2)*2.0-c2x*(curt-t5)*2.0)/(t2-t5)-((c1x-c2x)*cmt4*2.0)/t14+((c2x-c3x)*(curt-t2)*2.0)/(t2-t5))/(t2-t4)-((c3x*(curt-t2)*2.0-c2x*(curt-t5)*2.0)/(t2-t5)-(c4x*(cmt3)*2.0-c3x*(curt-t6)*2.0)/t36-((c2x-c3x)*(curt-t5)*2.0)/(t2-t5)+((c3x-c4x)*(cmt3)*2.0)/t36)/t35-(cmt4*((c1x*2.0-c2x*2.0)/t14-(c2x*2.0-c3x*2.0)/(t2-t5)))/(t2-t4)+((cmt3)*((c2x*2.0-c3x*2.0)/(t2-t5)-(c3x*2.0-c4x*2.0)/t36))/t35)*((((c2y*cmt1-c1y*cmt4)*cmt4)/t14-((c3y*(curt-t2)-c2y*(curt-t5))*(curt-t2))/(t2-t5))/(t2-t4)-(((c3y*(curt-t2)-c2y*(curt-t5))*(curt-t5))/(t2-t5)-((c4y*(cmt3)-c3y*(curt-t6))*(cmt3))/t36)/t35+(cmt4*((c2y*cmt1-c1y*cmt4)/t14-(c3y*(curt-t2)-c2y*(curt-t5))/(t2-t5)-((c1y-c2y)*cmt4)/t14+((c2y-c3y)*(curt-t2))/(t2-t5)))/(t2-t4)-((cmt3)*((c3y*(curt-t2)-c2y*(curt-t5))/(t2-t5)-(c4y*(cmt3)-c3y*(curt-t6))/t36-((c2y-c3y)*(curt-t5))/(t2-t5)+((c3y-c4y)*(cmt3))/t36))/t35)-1.0/t34sq*(((c2y*cmt1*2.0-c1y*cmt4*2.0)/t14-(c3y*(curt-t2)*2.0-c2y*(curt-t5)*2.0)/(t2-t5)-((c1y-c2y)*cmt4*2.0)/t14+((c2y-c3y)*(curt-t2)*2.0)/(t2-t5))/(t2-t4)-((c3y*(curt-t2)*2.0-c2y*(curt-t5)*2.0)/(t2-t5)-(c4y*(cmt3)*2.0-c3y*(curt-t6)*2.0)/t36-((c2y-c3y)*(curt-t5)*2.0)/(t2-t5)+((c3y-c4y)*(cmt3)*2.0)/t36)/t35-(cmt4*((c1y*2.0-c2y*2.0)/t14-(c2y*2.0-c3y*2.0)/(t2-t5)))/(t2-t4)+((cmt3)*((c2y*2.0-c3y*2.0)/(t2-t5)-(c3y*2.0-c4y*2.0)/t36))/t35)*((((c2x*cmt1-c1x*cmt4)*cmt4)/t14-((c3x*(curt-t2)-c2x*(curt-t5))*(curt-t2))/(t2-t5))/(t2-t4)-(((c3x*(curt-t2)-c2x*(curt-t5))*(curt-t5))/(t2-t5)-((c4x*(cmt3)-c3x*(curt-t6))*(cmt3))/t36)/t35+(cmt4*((c2x*cmt1-c1x*cmt4)/t14-(c3x*(curt-t2)-c2x*(curt-t5))/(t2-t5)-((c1x-c2x)*cmt4)/t14+((c2x-c3x)*(curt-t2))/(t2-t5)))/(t2-t4)-((cmt3)*((c3x*(curt-t2)-c2x*(curt-t5))/(t2-t5)-(c4x*(cmt3)-c3x*(curt-t6))/t36-((c2x-c3x)*(curt-t5))/(t2-t5)+((c3x-c4x)*(cmt3))/t36))/t35),2.0))*1.0/pow(epsil+1.0/t34sq*pow((((c2x*cmt1-c1x*cmt4)*cmt4)/t14-((c3x*(curt-t2)-c2x*(curt-t5))*(curt-t2))/(t2-t5))/(t2-t4)-(((c3x*(curt-t2)-c2x*(curt-t5))*(curt-t5))/(t2-t5)-((c4x*(cmt3)-c3x*(curt-t6))*(cmt3))/t36)/t35+(cmt4*((c2x*cmt1-c1x*cmt4)/t14-(c3x*(curt-t2)-c2x*(curt-t5))/(t2-t5)-((c1x-c2x)*cmt4)/t14+((c2x-c3x)*(curt-t2))/(t2-t5)))/(t2-t4)-((cmt3)*((c3x*(curt-t2)-c2x*(curt-t5))/(t2-t5)-(c4x*(cmt3)-c3x*(curt-t6))/t36-((c2x-c3x)*(curt-t5))/(t2-t5)+((c3x-c4x)*(cmt3))/t36))/t35,2.0)+1.0/t34sq*pow((((c2y*cmt1-c1y*cmt4)*cmt4)/t14-((c3y*(curt-t2)-c2y*(curt-t5))*(curt-t2))/(t2-t5))/(t2-t4)-(((c3y*(curt-t2)-c2y*(curt-t5))*(curt-t5))/(t2-t5)-((c4y*(cmt3)-c3y*(curt-t6))*(cmt3))/t36)/t35+(cmt4*((c2y*cmt1-c1y*cmt4)/t14-(c3y*(curt-t2)-c2y*(curt-t5))/(t2-t5)-((c1y-c2y)*cmt4)/t14+((c2y-c3y)*(curt-t2))/(t2-t5)))/(t2-t4)-((cmt3)*((c3y*(curt-t2)-c2y*(curt-t5))/(t2-t5)-(c4y*(cmt3)-c3y*(curt-t6))/t36-((c2y-c3y)*(curt-t5))/(t2-t5)+((c3y-c4y)*(cmt3))/t36))/t35,2.0),2.0))/(k+pow(hory+(((((c2y*cmt1-c1y*cmt4)*cmt4)/t14-((c3y*(curt-t2)-c2y*(curt-t5))*(curt-t2))/(t2-t5))*cmt4)/(t2-t4)-((((c3y*(curt-t2)-c2y*(curt-t5))*(curt-t5))/(t2-t5)-((c4y*(cmt3)-c3y*(curt-t6))*(cmt3))/t36)*(cmt3))/t35)/(t34),2.0))+1.0));
                *fx+=obj;
                dfx[svpos[0]] += ((frRsq*((cmt4sq*1.0/t34sq*tmp22*3.0)/(t14*t24)-(cmt4*1.0/t34sq*tmp24*6.0)/(t14*t24))*(1.0/t34sq*tmp23*tmp24-1.0/t34sq*tmp22*tmp25)*tmp29)/(k+pow(hory+tmp21,2.0))+(tmp28*cmt4sq*1.0/t34sq*tmp34*tmp25*1.2E1)/((k+pow(hory+tmp21,2.0))*t14*t24))/(tmp33/(k+pow(hory+tmp21,2.0))+1.0);
                dfx[svpos[1]] += ((frRsq*(1.0/t34sq*(tmp16+(curt*2.0-t3*2.0)/(t25*t35)+(curt*4.0-t5*4.0)/(t25*t35)+((2.0/t14+2.0/t25)*cmt4)/t24)*tmp24-1.0/t34sq*tmp9*tmp22)*(1.0/t34sq*tmp23*tmp24-1.0/t34sq*tmp22*tmp25)*tmp29)/(k+pow(hory+tmp21,2.0))-(tmp28*1.0/t34sq*tmp34*tmp9*tmp25*4.0)/(k+pow(hory+tmp21,2.0)))/(tmp33/(k+pow(hory+tmp21,2.0))+1.0);
                dfx[svpos[2]] += -((frRsq*(1.0/t34sq*(tmp17+(curt*4.0-t2*4.0)/(t24*t25)+(curt*2.0-t4*2.0)/(t24*t25)+((2.0/t25+2.0/t36)*cmt3)/t35)*tmp24-1.0/t34sq*tmp32*tmp22)*(1.0/t34sq*tmp23*tmp24-1.0/t34sq*tmp22*tmp25)*tmp29)/(k+pow(hory+tmp21,2.0))-(tmp28*1.0/t34sq*tmp34*tmp32*tmp25*4.0)/(k+pow(hory+tmp21,2.0)))/(tmp33/(k+pow(hory+tmp21,2.0))+1.0);
                dfx[svpos[3]] += -((frRsq*((cmt3sq*1.0/t34sq*tmp22*3.0)/(t35*t36)-(cmt3*1.0/t34sq*tmp24*6.0)/(t35*t36))*(1.0/t34sq*tmp23*tmp24-1.0/t34sq*tmp22*tmp25)*tmp29)/(k+pow(hory+tmp21,2.0))+(tmp28*cmt3sq*1.0/t34sq*tmp34*tmp25*1.2E1)/((k+pow(hory+tmp21,2.0))*t35*t36))/(tmp33/(k+pow(hory+tmp21,2.0))+1.0);
                dfx[svpos[4]] += ((frRsq*((cmt4sq*1.0/t34sq*tmp23*3.0)/(t14*t24)-(cmt4*1.0/t34sq*tmp25*6.0)/(t14*t24))*(1.0/t34sq*tmp23*tmp24-1.0/t34sq*tmp22*tmp25)*1.0/pow(epsil+1.0/t34sq*tmp27+1.0/t34sq*tmp26,2.0)*-2.0)/(k+pow(hory+tmp21,2.0))+(tmp28*(hory+tmp21)*pow(cmt4,3.0)*1.0/pow(k+pow(hory+tmp21,2.0),2.0)*tmp29)/(t14*t24*t34)+(tmp28*cmt4sq*1.0/t34sq*tmp34*tmp24*1.2E1)/((k+pow(hory+tmp21,2.0))*t14*t24))/(tmp33/(k+pow(hory+tmp21,2.0))+1.0);
                dfx[svpos[5]] += -((frRsq*(1.0/t34sq*(tmp16+(curt*2.0-t3*2.0)/(t25*t35)+(curt*4.0-t5*4.0)/(t25*t35)+((2.0/t14+2.0/t25)*cmt4)/t24)*tmp25-1.0/t34sq*tmp9*tmp23)*(1.0/t34sq*tmp23*tmp24-1.0/t34sq*tmp22*tmp25)*tmp29)/(k+pow(hory+tmp21,2.0))+(tmp28*1.0/t34sq*tmp34*tmp9*tmp24*4.0)/(k+pow(hory+tmp21,2.0))+(tmp28*(hory+tmp21)*1.0/pow(k+pow(hory+tmp21,2.0),2.0)*((((cmt1*cmt4)/t14+(cmt2*cmt5)/t25)*cmt4)/t24+(cmt3*pow(curt-t5,2.0))/(t25*t35))*tmp29)/t34)/(tmp33/(k+pow(hory+tmp21,2.0))+1.0);
                dfx[svpos[6]] += ((frRsq*(1.0/t34sq*(tmp17+(curt*4.0-t2*4.0)/(t24*t25)+(curt*2.0-t4*2.0)/(t24*t25)+((2.0/t25+2.0/t36)*cmt3)/t35)*tmp25-1.0/t34sq*tmp32*tmp23)*(1.0/t34sq*tmp23*tmp24-1.0/t34sq*tmp22*tmp25)*tmp29)/(k+pow(hory+tmp21,2.0))+(tmp28*1.0/t34sq*tmp34*tmp32*tmp24*4.0)/(k+pow(hory+tmp21,2.0))+(tmp28*(hory+tmp21)*1.0/pow(k+pow(hory+tmp21,2.0),2.0)*((((cmt2*cmt5)/t25+(cmt3*cmt6)/t36)*cmt3)/t35+(pow(curt-t2,2.0)*cmt4)/(t24*t25))*tmp29)/t34)/(tmp33/(k+pow(hory+tmp21,2.0))+1.0);
                dfx[svpos[7]] += -((frRsq*((cmt3sq*1.0/t34sq*tmp23*3.0)/(t35*t36)-(cmt3*1.0/t34sq*tmp25*6.0)/(t35*t36))*(1.0/t34sq*tmp23*tmp24-1.0/t34sq*tmp22*tmp25)*1.0/pow(epsil+1.0/t34sq*tmp27+1.0/t34sq*tmp26,2.0)*-2.0)/(k+pow(hory+tmp21,2.0))+(tmp28*(hory+tmp21)*pow(cmt3,3.0)*1.0/pow(k+pow(hory+tmp21,2.0),2.0)*tmp29)/(t34*t35*t36)+(tmp28*cmt3sq*1.0/t34sq*tmp34*tmp24*1.2E1)/((k+pow(hory+tmp21,2.0))*t35*t36))/(tmp33/(k+pow(hory+tmp21,2.0))+1.0);
            }
            else{
                
                obj= -log(1.0/(frRsq*(epsil+pow(1.0/t34sq*(((c2x*cmt1*2.0-c1x*cmt4*2.0)/t14-(c3x*(curt-t2)*2.0-c2x*(curt-t5)*2.0)/(t2-t5)-((c1x-c2x)*cmt4*2.0)/t14+((c2x-c3x)*(curt-t2)*2.0)/(t2-t5))/(t2-t4)-((c3x*(curt-t2)*2.0-c2x*(curt-t5)*2.0)/(t2-t5)-(c4x*(cmt3)*2.0-c3x*(curt-t6)*2.0)/t36-((c2x-c3x)*(curt-t5)*2.0)/(t2-t5)+((c3x-c4x)*(cmt3)*2.0)/t36)/t35-(cmt4*((c1x*2.0-c2x*2.0)/t14-(c2x*2.0-c3x*2.0)/(t2-t5)))/(t2-t4)+((cmt3)*((c2x*2.0-c3x*2.0)/(t2-t5)-(c3x*2.0-c4x*2.0)/t36))/t35)*((((c2y*cmt1-c1y*cmt4)*cmt4)/t14-((c3y*(curt-t2)-c2y*(curt-t5))*(curt-t2))/(t2-t5))/(t2-t4)-(((c3y*(curt-t2)-c2y*(curt-t5))*(curt-t5))/(t2-t5)-((c4y*(cmt3)-c3y*(curt-t6))*(cmt3))/t36)/t35+(cmt4*((c2y*cmt1-c1y*cmt4)/t14-(c3y*(curt-t2)-c2y*(curt-t5))/(t2-t5)-((c1y-c2y)*cmt4)/t14+((c2y-c3y)*(curt-t2))/(t2-t5)))/(t2-t4)-((cmt3)*((c3y*(curt-t2)-c2y*(curt-t5))/(t2-t5)-(c4y*(cmt3)-c3y*(curt-t6))/t36-((c2y-c3y)*(curt-t5))/(t2-t5)+((c3y-c4y)*(cmt3))/t36))/t35)-1.0/t34sq*(((c2y*cmt1*2.0-c1y*cmt4*2.0)/t14-(c3y*(curt-t2)*2.0-c2y*(curt-t5)*2.0)/(t2-t5)-((c1y-c2y)*cmt4*2.0)/t14+((c2y-c3y)*(curt-t2)*2.0)/(t2-t5))/(t2-t4)-((c3y*(curt-t2)*2.0-c2y*(curt-t5)*2.0)/(t2-t5)-(c4y*(cmt3)*2.0-c3y*(curt-t6)*2.0)/t36-((c2y-c3y)*(curt-t5)*2.0)/(t2-t5)+((c3y-c4y)*(cmt3)*2.0)/t36)/t35-(cmt4*((c1y*2.0-c2y*2.0)/t14-(c2y*2.0-c3y*2.0)/(t2-t5)))/(t2-t4)+((cmt3)*((c2y*2.0-c3y*2.0)/(t2-t5)-(c3y*2.0-c4y*2.0)/t36))/t35)*((((c2x*cmt1-c1x*cmt4)*cmt4)/t14-((c3x*(curt-t2)-c2x*(curt-t5))*(curt-t2))/(t2-t5))/(t2-t4)-(((c3x*(curt-t2)-c2x*(curt-t5))*(curt-t5))/(t2-t5)-((c4x*(cmt3)-c3x*(curt-t6))*(cmt3))/t36)/t35+(cmt4*((c2x*cmt1-c1x*cmt4)/t14-(c3x*(curt-t2)-c2x*(curt-t5))/(t2-t5)-((c1x-c2x)*cmt4)/t14+((c2x-c3x)*(curt-t2))/(t2-t5)))/(t2-t4)-((cmt3)*((c3x*(curt-t2)-c2x*(curt-t5))/(t2-t5)-(c4x*(cmt3)-c3x*(curt-t6))/t36-((c2x-c3x)*(curt-t5))/(t2-t5)+((c3x-c4x)*(cmt3))/t36))/t35),2.0))*1.0/pow(epsil+1.0/t34sq*pow((((c2x*cmt1-c1x*cmt4)*cmt4)/t14-((c3x*(curt-t2)-c2x*(curt-t5))*(curt-t2))/(t2-t5))/(t2-t4)-(((c3x*(curt-t2)-c2x*(curt-t5))*(curt-t5))/(t2-t5)-((c4x*(cmt3)-c3x*(curt-t6))*(cmt3))/t36)/t35+(cmt4*((c2x*cmt1-c1x*cmt4)/t14-(c3x*(curt-t2)-c2x*(curt-t5))/(t2-t5)-((c1x-c2x)*cmt4)/t14+((c2x-c3x)*(curt-t2))/(t2-t5)))/(t2-t4)-((cmt3)*((c3x*(curt-t2)-c2x*(curt-t5))/(t2-t5)-(c4x*(cmt3)-c3x*(curt-t6))/t36-((c2x-c3x)*(curt-t5))/(t2-t5)+((c3x-c4x)*(cmt3))/t36))/t35,2.0)+1.0/t34sq*pow((((c2y*cmt1-c1y*cmt4)*cmt4)/t14-((c3y*(curt-t2)-c2y*(curt-t5))*(curt-t2))/(t2-t5))/(t2-t4)-(((c3y*(curt-t2)-c2y*(curt-t5))*(curt-t5))/(t2-t5)-((c4y*(cmt3)-c3y*(curt-t6))*(cmt3))/t36)/t35+(cmt4*((c2y*cmt1-c1y*cmt4)/t14-(c3y*(curt-t2)-c2y*(curt-t5))/(t2-t5)-((c1y-c2y)*cmt4)/t14+((c2y-c3y)*(curt-t2))/(t2-t5)))/(t2-t4)-((cmt3)*((c3y*(curt-t2)-c2y*(curt-t5))/(t2-t5)-(c4y*(cmt3)-c3y*(curt-t6))/t36-((c2y-c3y)*(curt-t5))/(t2-t5)+((c3y-c4y)*(cmt3))/t36))/t35,2.0),2.0)+1.0));
                *fx+=obj;
                dfx[svpos[0]] += (frRsq*((cmt4sq*1.0/t34sq*tmp22*3.0)/(t14*t24)-(cmt4*1.0/t34sq*tmp24*6.0)/(t14*t24))*(1.0/t34sq*tmp23*tmp24-1.0/t34sq*tmp22*tmp25)*tmp29+(tmp28*cmt4sq*1.0/t34sq*tmp34*tmp25*1.2E1)/(t14*t24))/tmp30;
                dfx[svpos[1]] += (frRsq*(1.0/t34sq*(tmp16+(curt*2.0-t3*2.0)/(t25*t35)+(curt*4.0-t5*4.0)/(t25*t35)+((2.0/t14+2.0/t25)*cmt4)/t24)*tmp24-1.0/t34sq*tmp9*tmp22)*(1.0/t34sq*tmp23*tmp24-1.0/t34sq*tmp22*tmp25)*tmp29-tmp28*1.0/t34sq*tmp34*tmp9*tmp25*4.0)/tmp30;
                dfx[svpos[2]] += -(frRsq*(1.0/t34sq*(tmp17+(curt*4.0-t2*4.0)/(t24*t25)+(curt*2.0-t4*2.0)/(t24*t25)+((2.0/t25+2.0/t36)*cmt3)/t35)*tmp24-1.0/t34sq*tmp32*tmp22)*(1.0/t34sq*tmp23*tmp24-1.0/t34sq*tmp22*tmp25)*tmp29-tmp28*1.0/t34sq*tmp34*tmp32*tmp25*4.0)/tmp30;
                dfx[svpos[3]] += -(frRsq*((cmt3sq*1.0/t34sq*tmp22*3.0)/(t35*t36)-(cmt3*1.0/t34sq*tmp24*6.0)/(t35*t36))*(1.0/t34sq*tmp23*tmp24-1.0/t34sq*tmp22*tmp25)*tmp29+(tmp28*cmt3sq*1.0/t34sq*tmp34*tmp25*1.2E1)/(t35*t36))/tmp30;
                dfx[svpos[4]] += -(frRsq*((cmt4sq*1.0/t34sq*tmp23*3.0)/(t14*t24)-(cmt4*1.0/t34sq*tmp25*6.0)/(t14*t24))*(1.0/t34sq*tmp23*tmp24-1.0/t34sq*tmp22*tmp25)*tmp29-(tmp28*cmt4sq*1.0/t34sq*tmp34*tmp24*1.2E1)/(t14*t24))/tmp30;
                dfx[svpos[5]] += -(frRsq*(1.0/t34sq*(tmp16+(curt*2.0-t3*2.0)/(t25*t35)+(curt*4.0-t5*4.0)/(t25*t35)+((2.0/t14+2.0/t25)*cmt4)/t24)*tmp25-1.0/t34sq*tmp9*tmp23)*(1.0/t34sq*tmp23*tmp24-1.0/t34sq*tmp22*tmp25)*tmp29+tmp28*1.0/t34sq*tmp34*tmp9*tmp24*4.0)/tmp30;
                dfx[svpos[6]] += (frRsq*(1.0/t34sq*(tmp17+(curt*4.0-t2*4.0)/(t24*t25)+(curt*2.0-t4*2.0)/(t24*t25)+((2.0/t25+2.0/t36)*cmt3)/t35)*tmp25-1.0/t34sq*tmp32*tmp23)*(1.0/t34sq*tmp23*tmp24-1.0/t34sq*tmp22*tmp25)*tmp29+tmp28*1.0/t34sq*tmp34*tmp32*tmp24*4.0)/tmp30;
                dfx[svpos[7]] += (frRsq*((cmt3sq*1.0/t34sq*tmp23*3.0)/(t35*t36)-(cmt3*1.0/t34sq*tmp25*6.0)/(t35*t36))*(1.0/t34sq*tmp23*tmp24-1.0/t34sq*tmp22*tmp25)*tmp29-(tmp28*cmt3sq*1.0/t34sq*tmp34*tmp24*1.2E1)/(t35*t36))/tmp30;
            }
            
            
            cnt++;
        }
        splos += splpieces*8;
        knoffset += (splpieces+7);
        indexoffset += cnt;
    }
    
    
    free(svpos);
}
