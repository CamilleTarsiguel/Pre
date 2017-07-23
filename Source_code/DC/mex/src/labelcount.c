#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {


	/* //Declarations */
	const mxArray *labelingdata;
    const mxArray *nLabelsData;
    

	double *labeling;
	double *res;
    double *nLabels;
    
	int L,i;
	
	labelingdata = prhs[0];
    nLabelsData = prhs[1];
/*    printf("%d %d\n",nrhs, nlhs);*/
	
	/* //Get matrix x */
	labeling = mxGetPr(labelingdata);
    nLabels = mxGetPr(nLabelsData);
	
	L = mxGetN(labelingdata); 
    
    /* Allocate memory and assign output pointer */
    plhs[0] = mxCreateDoubleMatrix(1,(int)nLabels[0],mxREAL);
	res = mxGetPr(plhs[0]);
	
	for (i=0; i<L;i++) {
        /*printf("%d %d\n",i,(int)labeling[i]);*/
		res[(int)labeling[i]-1]=res[(int)labeling[i]-1]+1;
	}

}