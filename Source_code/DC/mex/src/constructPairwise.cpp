#include "mex.h"
#include <math.h>
#include <vector>

// TODO: Fix data types! doubles to int, etc.
// TODO: sparse matrices instead of full for the Neighborhood
// TODO: memcpy instead of loop at end
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    
    /* //Declarations */
    const mxArray *TNData;
    const mxArray *nAlphaData;
    const mxArray *labelingdata;
    const mxArray *oldata;
    
    
    double *TNeighborsAux;
    double *Dpw1, *Dpw2, *Dpwi;
    double *oLabelA;
    double oLabel;
    int tmp1,tmp2;
//     double *notAlphaInd;
    double *alphaIndicator;
    double *labeling;
    int nAlphaNodes;
    
//     double *nLabels;
    
    int nPoints,i,nNeighbors;
    int a;
    int idx;
    
    TNData = prhs[0];
    nAlphaData = prhs[1];
    labelingdata = prhs[2];
    oldata = prhs[3];
    /*    nLabelsData = prhs[1];
     * printf("%d %d\n",nrhs, nlhs);*/
    
    /* //Get matrix x */
    TNeighborsAux = mxGetPr(TNData);
    alphaIndicator = mxGetPr(nAlphaData);
    labeling = mxGetPr(labelingdata);
    oLabelA = mxGetPr(oldata);
    oLabel = oLabelA[0];
//     notAlphaInd = (int*)(nAlphaData);
//     nLabels = mxGetPr(nLabelsData);
    
    nPoints = mxGetN(TNData);
    nAlphaNodes = mxGetN(nAlphaData);
    
    std::vector<int> N1;
    std::vector<int> N2;
    std::vector<int> Id;
    
//     double UnAux[nPoints];
    
    // count notalphas
    int nNotAlphas=0;
    for (a=0; a<nPoints; a++)
        if (!alphaIndicator[a])
            nNotAlphas++;
        
    double *UnAux;
    plhs[0] = mxCreateDoubleMatrix(1,nNotAlphas,mxREAL);	    UnAux= mxGetPr(plhs[0]);
    
    
    //printf("deb a");
    tmp1=0;
    for (a=0; a<nPoints; a++) {
        if (!alphaIndicator[a]) {
            tmp1++;
            tmp2=0;
            for (i=0; i<nPoints; i++) {
                if (!alphaIndicator[i])
                    tmp2++;
                
                idx=a*nPoints+i;
                if (TNeighborsAux[idx]) {
                    if (alphaIndicator[i]) {
                        UnAux[tmp1-1]++;
                    } else {
                        if (tmp1<tmp2) {
                            N1.push_back(tmp1); N2.push_back(tmp2);
                        }else{
                            N1.push_back(tmp2); N2.push_back(tmp1);
                            
                        }

                        // diff between temp pairwise and exclusions
                        if (oLabel>0) { // exclusions
                            if (labeling[a]==labeling[i] && labeling[a]!=oLabel){
                                Id.push_back(1);
                            } else {
                                Id.push_back(0);
                            }
                        }
                        else
                        {
                            if (labeling[a]==labeling[i]){
                                Id.push_back(1);
                            } else {
                                Id.push_back(0);
                            }
                        }
                            
                        
                        TNeighborsAux[idx]=0;
                        TNeighborsAux[i*nPoints+a]=0;
                    }
                }
            }
        }
    }
    //printf("deb b: %d",N1.size());
//     }
//                    tmp1=0;
//             for nodeId=notAlphasInd
//                 tmp1=tmp1+1;
//                 nodesNeighbors=find(TNeighborsAux(nodeId,:));
//
//                 tmp2=0;
//                 for nnodeId=nodesNeighbors
//                     %                     nnodeId
//                     % neighbor is alpha
//                     if ~isempty(find(alpInds==nnodeId, 1))
//                         UnAux(2,tmp1)=UnAux(2,tmp1)-opt.pairwiseFactor;
//                     else % neighbor is not alpha, keep pairwise connection
//                         tmp2=find(notAlphasInd==nnodeId);
//                         if tmp1<tmp2
//                             Dpw=[Dpw [tmp1;tmp2]];
//                         else
//                             Dpw=[Dpw [tmp2;tmp1]];
//                         end
//
//
//                         if labeling(nodeId)==labeling(nnodeId) % both are equal, so either leave both or flip both
//                             Dpwi=[Dpwi 1];
//                         else
//                             Dpwi=[Dpwi 0];
//                         end
//                         TNeighborsAux(tmp1,tmp2)=0;TNeighborsAux(tmp2,tmp1)=0;
//                     end
//                     %                     UnAux
//                     %                     Dpw
//                     %                     Dpwi
//                     %                     pause
//                 end
    
    //printf("deb c");
    nNeighbors=N1.size();
    //printf("deb %d: ",nNeighbors);
    /* Allocate memory and assign output pointer */
    plhs[1] = mxCreateDoubleMatrix(1,nNeighbors,mxREAL);	Dpw1 = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(1,nNeighbors,mxREAL);	Dpw2 = mxGetPr(plhs[2]);
    plhs[3] = mxCreateDoubleMatrix(1,nNeighbors,mxREAL);	Dpwi = mxGetPr(plhs[3]);
    //printf("deb d");
    
    for (i=0; i<N1.size(); i++){
//             Dpw1[i]=0;Dpw2[i]=0;
        Dpw1[i]=N1[i];
        Dpw2[i]=N2[i];
        Dpwi[i]=Id[i];
    }
    //printf("deb d\n");
// 	for (i=0; i<nPoints;i++) {
//         /*printf("%d %d\n",i,(int)TNeighborsAux[i]);*/
// 		Dpw[(int)TNeighborsAux[i]-1]=Dpw[(int)TNeighborsAux[i]-1]+1;
// 	}
    
}