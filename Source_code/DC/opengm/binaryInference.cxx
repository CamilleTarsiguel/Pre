//include matlab headers
#include "mex.h"

#include <vector>


#include <opengm/graphicalmodel/graphicalmodel.hxx>
// #include <opengm/graphicalmodel/graphicalmodel_hdf5.hxx>
//#include <opengm/utilities/metaprogramming.hxx>
//#include <opengm/graphicalmodel/space/discretespace.hxx>
#include <opengm/graphicalmodel/space/simplediscretespace.hxx>
#include <opengm/functions/explicit_function.hxx>
#include <opengm/functions/potts.hxx>


// inference
#include <opengm/operations/minimizer.hxx>
#include <opengm/operations/adder.hxx>


// TRW
#include <opengm/inference/external/trws.hxx>

// QPBO 
#include <opengm/inference/external/qpbo.hxx>

// MQPBO
// #include <opengm/inference/mqpbo.hxx>

// Gibbs
#include <opengm/inference/gibbs.hxx>

// FastPD
// #include <opengm/inference/external/fastPD.hxx>

#include <opengm/inference/trws/trws_trws.hxx>



typedef double ValueType;
typedef size_t IndexType;
typedef size_t LabelType;
typedef opengm::GraphicalModel<
   ValueType,
   opengm::Adder,
   opengm::meta::TypeListGenerator<opengm::ExplicitFunction<ValueType> >::type,
   opengm::DiscreteSpace<IndexType, LabelType>
   > GraphicalModelType;

typedef opengm::ExplicitFunction<ValueType> ExplicitFunctionType;
typedef opengm::PottsFunction<ValueType> PottsFunctionType;
typedef GraphicalModelType::FunctionIdentifier FunctionIdentifier; 

typedef opengm::SimpleDiscreteSpace<> Space;
//typedef opengm::DiscreteSpace<> DSpace;

//typedef OPENGM_TYPELIST_2(ExplicitFunction<double>, PottsFunction<double>) FunctionTypelist;
//typedef opengm::GraphicalModel<double, opengm::Adder, FunctionTypelist, Space> Model;
typedef opengm::GraphicalModel<double, opengm::Adder, 
   OPENGM_TYPELIST_2(ExplicitFunctionType, PottsFunctionType), Space> Model; 
//typedef opengm::GraphicalModel<double, opengm::Adder, 
//   OPENGM_TYPELIST_2(ExplicitFunctionType, PottsFunctionType), DSpace> Model2; 
//typedef opengm::GraphicalModel<double, opengm::Adder> GmType;

//TRW-S*
typedef opengm::external::TRWS<Model> TRWS;

// QPBO
typedef opengm::external::QPBO<Model> QPBO;


// MQPBO
//typedef opengm::Minimizer gmmin;
//typedef opengm::MQPBO<GmType,opengm::Minimizer> MQPBOType;
// typedef opengm::MQPBO<Model,opengm::Minimizer> MQPBO;

// Gibbs
typedef opengm::Gibbs<Model, opengm::Minimizer> Gibbs;

// FastPD
// typedef opengm::external::FastPD<Model> FastPD;

// ICM
#include <opengm/inference/icm.hxx>
typedef opengm::ICM<Model, opengm::Minimizer> ICM;

typedef opengm::TRWSi<Model, opengm::Minimizer> TRWSi;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  //check if data is in correct format
  if(nlhs != 2) {
     mexErrMsgTxt("Wrong number of output variables specified (two needed)\n");
  }
  if(nrhs != 17) {
     mexErrMsgTxt("Wrong number of input variables specified (17 needed)\n");
  }

   double* unary             = (double*) mxGetData(prhs[0]);
   double* pairwise_node     = (double*) mxGetData(prhs[1]); 
   double* pairwise_fac      = (double*) mxGetData(prhs[2]); 
   double* lcost             = (double*) mxGetData(prhs[3]); 
   double* lcost_ind         = (double*) mxGetData(prhs[4]); 
   double* exclusion_fac     = (double*) mxGetData(prhs[5]);
   double* new_pw_node       = (double*) mxGetData(prhs[7]);
   double* new_pw_node_id    = (double*) mxGetData(prhs[8]);
   double* new_epw_node      = (double*) mxGetData(prhs[9]);
   double* new_epw_node_id   = (double*) mxGetData(prhs[10]);
   double* notAlphasNodes    = (double*) mxGetData(prhs[11]);
   double* pwLcost 			  = (double*) mxGetData(prhs[12]);
   double* nAux				 = (double*) mxGetData(prhs[13]);
   double* inf_params         = (double*) mxGetData(prhs[14]);
   double* auxSwitch         = (double*) mxGetData(prhs[15]);
   double* InfAlg            = (double*) mxGetData(prhs[16]);

   
  
   size_t numberOfVariables = mxGetN(prhs[0]);
   size_t numberOfLabels = mxGetM(prhs[0]);
   size_t outlierLabel = numberOfLabels-1;
   size_t numberOfPairwiseFactors   = mxGetN(prhs[1]);
   size_t numberOfAuxiliaryVars   = (int)nAux[0];
   size_t numberOfNewPairwiseFactors   = mxGetN(prhs[7]);
   size_t numberOfNewEPairwiseFactors  = mxGetN(prhs[9]);
   size_t numberOfPairwiseLcosts  = mxGetN(prhs[12]);

   const double Lconst=1000000;

   std::vector<LabelType> numStates(numberOfVariables,numberOfVariables);
   //GraphicalModelType gm(opengm::DiscreteSpace<IndexType, LabelType >(numStates.begin(), numStates.end()) );
   Space myspace(numberOfVariables, numberOfLabels); 
   
   Model mymodel(myspace);
   //Model mygmmodel(myspace);
   //gmmin mymin();


   //const int aaa = 23;//numberOfLabels;
   //size_t numbersOfLabels[aaa];
   //for(size_t var=0; var<numberOfVariables; ++var) {
	  // numbersOfLabels[var]=numberOfLabels;
   //}

   
   //DSpace mydspace();
   //Model2 mymodel(mydspace);
   //space(numbersOfLabels, numbersOfLabels + 3);
   //opengm::GraphicalModel<double,opengm::Adder, 
	//   OPENGM_TYPELIST_2(ExplicitFunctionType, PottsFunctionType), SS> mymodel;

//   for(size_t var=0; var<numberOfVariables; ++var) {
	//   mymodel.addVariable(
   //}

	// UNARIES
   double zerosonly=0;
   for(size_t var=0; var<numberOfVariables; ++var) {
      LabelType shape[] = {numberOfLabels};
      ExplicitFunctionType func(shape,shape+1); 
      for(size_t i=0; i<numberOfLabels; ++i) {
         func(i) = unary[i+var*numberOfLabels];		 
      }
      zerosonly+=unary[var*numberOfLabels];
      FunctionIdentifier fid = mymodel.addFunction(func);
      mymodel.addFactor(fid, &var, &var+1);
   }

      //ADD Pairwise (potts)
   if (numberOfNewPairwiseFactors>0) {	   
	   LabelType shape[] = {numberOfLabels, numberOfLabels};
	   
		PottsFunctionType pottsf(numberOfLabels, numberOfLabels, 0.0, pairwise_fac[0]); 
		Model::FunctionIdentifier fidp = mymodel.addFunction(pottsf);
		
		for(size_t f=0;f<numberOfNewPairwiseFactors;++f) {			
			IndexType v1 = (IndexType)new_pw_node[2*f];
			IndexType v2 = (IndexType)new_pw_node[2*f+1];		
			
			size_t variableIndices[] = {v1,v2};
			
			if ((int)new_pw_node_id[f]) {
				mymodel.addFactor(fidp,variableIndices, variableIndices+2);
				//printf("p1 %d %d %d %f\n",f,v1,v2,-pairwise_fac[0]);
			} else {
				ExplicitFunctionType exExplf(shape,shape+2,0.0);
				exExplf(1,1)=-pairwise_fac[0];
				Model::FunctionIdentifier fide = mymodel.addFunction(exExplf);		
				mymodel.addFactor(fide,variableIndices, variableIndices+2);
				//printf("p0 %d %d %d %f\n",f,v1,v2,-pairwise_fac[0]);
			}
				
			
		}
   }   

	// Exclusion New Pairwise (reverse potts)
   if (numberOfNewEPairwiseFactors>0) {	   
	   LabelType shape[] = {numberOfLabels, numberOfLabels};
	   
		PottsFunctionType pottsf(numberOfLabels, numberOfLabels, exclusion_fac[0], 0.0); 
		Model::FunctionIdentifier fidp = mymodel.addFunction(pottsf);
		
		for(size_t f=0;f<numberOfNewEPairwiseFactors;++f) {			
			IndexType v1 = (IndexType)new_epw_node[2*f];
			IndexType v2 = (IndexType)new_epw_node[2*f+1];		
			
			size_t variableIndices[] = {v1,v2};
			
			if ((int)new_epw_node_id[f]) {
				mymodel.addFactor(fidp,variableIndices, variableIndices+2);
				//printf("e1 %d %d %d %f\n",f,v1,v2,exclusion_fac[0]);
			} else {
				ExplicitFunctionType exExplf(shape,shape+2,0.0);
				exExplf(1,1)=exclusion_fac[0];
				Model::FunctionIdentifier fide = mymodel.addFunction(exExplf);		
				mymodel.addFactor(fide,variableIndices, variableIndices+2);
				//printf("e0 %d %d %d %f\n",f,v1,v2,exclusion_fac[0]);
			}
				
			
		}
   }   

   // this is for mexCallMATLAB
	mxArray *array_ptr;
	mxArray *input_array[1];
	
	array_ptr = mxCreateDoubleMatrix(1, 1, mxREAL);
	double pr[] = {.5};
	memcpy(mxGetPr(array_ptr), pr, sizeof(pr));
   input_array[0] = array_ptr;
   int      num_out, num_in;
   num_out=0; num_in=1;
  //printf("time - D %d (%f sec)\n", clock()-start_time, (float) (clock() - start_time) / CLOCKS_PER_SEC); //start_time = clock(); 

   // Auxiliary pairwise connections for label cost
   if (numberOfPairwiseFactors>0) {
	   LabelType shape[] = {numberOfLabels, numberOfLabels};
	   
		
		for(size_t f=0;f<numberOfPairwiseFactors;++f) {
			IndexType v1 = (IndexType)pairwise_node[2*f];
			IndexType v2 = (IndexType)pairwise_node[2*f+1];
			
			// THIS IS INEFFICIENT
			// not every pairwise factor has a different function
			// consider reducing number of functions
			ExplicitFunctionType exExplf(shape,shape+2,0.0); 
			
// 			exExplf(0,1)=lcost[(int)lcost_ind[f]]; 
			exExplf(0,0)=lcost[(int)lcost_ind[f]]+Lconst; 
// 			exExplf(1,0)=0;
			
		
			size_t variableIndices[] = {v1,v2};
			//printf("%d %d %d %d %f\n",f,v1,v2,(int)lcost_ind[f],lcost[(int)lcost_ind[f]]);
			//printf("%d %d %d %d %f\n",f,v1,v2,(int)lcost_ind[f],lcost[(int)lcost_ind[f]]+Lconst);
// 			//printf("%f %f %f %f\n",exExplf(0,0),exExplf(0,1),exExplf(1,0),exExplf(1,1));
			
			Model::FunctionIdentifier fide = mymodel.addFunction(exExplf);		
			mymodel.addFactor(fide,variableIndices, variableIndices+2);
		
		}
   }   
   
   if (numberOfPairwiseLcosts>0) {
	   LabelType shape[] = {numberOfLabels, numberOfLabels};
	   for(size_t f=0;f<numberOfPairwiseLcosts;++f) {
		   for(size_t g=f+1;g<numberOfPairwiseLcosts;++g) {
			   if (pwLcost[numberOfPairwiseLcosts*f+g] > 0) {
				IndexType v1 = (IndexType)f;
				IndexType v2 = (IndexType)g;
				
				ExplicitFunctionType exExplf(shape,shape+2,0.0); 
				exExplf(1,1)=pwLcost[numberOfPairwiseLcosts*f+g];
			   
				size_t variableIndices[] = {v1,v2};
				//printf("%d %d %d %d %f\n",f,g,v1,v2,pwLcost[numberOfPairwiseLcosts*f+g]);
// 				//printf("%f %f %f %f\n",exExplf(0,0),exExplf(0,1),exExplf(1,0),exExplf(1,1));
				
				Model::FunctionIdentifier fide = mymodel.addFunction(exExplf);		
				mymodel.addFactor(fide,variableIndices, variableIndices+2);

			   }
			}
		   
	   }
	   
   }

	double* labeling = (double*) mxGetData(prhs[6]);
	size_t labelingsize = mxGetN(prhs[6]);
	Model::ValueType value = 0;	//printf("time - E %d (%f sec)\n", clock()-start_time, (float) (clock() - start_time) / CLOCKS_PER_SEC); //start_time = clock(); 
	
	double *Eval;
	double *labelingoutput;
	
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	Eval = mxGetPr(plhs[0]);
	*Eval=value;

	plhs[1] = mxCreateDoubleMatrix(1, numberOfVariables, mxREAL);
	labelingoutput = mxGetPr(plhs[1]);

	
   	// Inference
	std::vector<Model::LabelType> x;
	double val=0;

			std::vector<LabelType> startingPoint;
		for(size_t var=0; var<numberOfVariables; ++var) {
			startingPoint.push_back(static_cast<LabelType>(labeling[var]));
// 			mexPrintf("%d ",static_cast<LabelType>(labeling[var]));
		}
// 		mexPrintf("\n");

	unsigned int ia = (int)InfAlg[0];
	if (ia==1) {
		if ((int)mxGetN(prhs[14])!=3)
			mexErrMsgTxt("3 Parameters needed for TRWS\n");

		TRWS::Parameter param;
		param.numberOfIterations_ = (size_t)inf_params[0];
		param.useRandomStart_ = (bool)inf_params[1];
		param.doBPS_ = (bool)inf_params[2];		
		param.energyType_ = TRWS::Parameter::VIEW;

		TRWS trws(mymodel,param);
		trws.infer();
		trws.arg(x);
		val=trws.value();
// 		printf("TRWS minvalue %f\n",trws.value());
	} 
	else if (ia==2)
	{
				// QPBO
		if ((int)mxGetN(prhs[14])!=3)
			mexErrMsgTxt("3 Parameters needed for QPBO\n");

 		QPBO::Parameter qparam;
		qparam.useProbeing_=(bool)inf_params[0];
		qparam.strongPersistency_=(bool)inf_params[1];
		qparam.useImproveing_=(bool)inf_params[2];

			
		QPBO qpbo(mymodel,qparam);
		qpbo.setStartingPoint(startingPoint.begin());
// 		opengm::hdf5::save(mymodel, "qpbo_crash.h5", "mymodel");
// 		std::vector<LabelType> startingPoint0;
// 		std::vector<LabelType> startingPoint1;
// 		for(size_t var=0; var<numberOfVariables; ++var) {
// 			startingPoint0.push_back(static_cast<LabelType>(0));
// 			startingPoint1.push_back(static_cast<LabelType>(1));
// 		}
// 		qpbo.setStartingPoint(startingPoint0.begin());
// 		printf("QPBO minvalue %f\n",qpbo.value());
// 		qpbo.setStartingPoint(startingPoint1.begin());
// 		printf("QPBO minvalue %f\n",qpbo.value());
 		qpbo.infer();		
		qpbo.arg(x);
		val=qpbo.value();
		
// 		printf("QPBO minvalue %f\n",qpbo.value());

		//
		//typedef MQPBOType::Parameter mqpboPara;
		//para.useKovtunsMethod_=1;
		//TRWS trws(a);
		//MQPBOType mqpbo(mymodel);

		//printf("%d %d %d\n", mymodel.numberOfFactors(),mymodel.numberOfFunctions(0),mymodel.numberOfVariables());
	} 
	else if (ia==3)
	{
					// QPBO
		// if ((int)mxGetN(prhs[14])!=3)
			// mexErrMsgTxt("3 Parameters needed for QPBO\n");

 		// QPBO::Parameter qparam;
		// qparam.useProbeing_=(bool)inf_params[0];
		// qparam.strongPersistency_=(bool)inf_params[1];
		// qparam.useImproveing_=(bool)inf_params[2];

			
		// QPBO qpbo(mymodel,qparam);
		// qpbo.setStartingPoint(startingPoint.begin());
 		// qpbo.infer();		
		// qpbo.arg(x);
		// val=qpbo.value();
		// if ((int)mxGetN(prhs[14])!=4)
			// mexErrMsgTxt("3 Parameters needed for MQPBO\n");

		// MQPBO::Parameter mqp;
		// mqp.rounds_=(size_t)inf_params[0];
		// mqp.strongPersistency_=(bool)inf_params[1];
		// mqp.useKovtunsMethod_=(bool)inf_params[2];
		// enum PermutationType {NONE, RANDOM, MINMARG};
		// mqp.permutationType_=MQPBO::PermutationType((size_t)inf_params[4]);
		// MQPBO mqpbo(mymodel,mqp);
		
		// mqpbo.setStartingPoint(startingPoint.begin());
		// mqpbo.infer();
		// mqpbo.arg(x);
		// val=mqpbo.value();
		//printf("MQPBO minvalue %f\n",mqpbo.value());
		
		//FastPD fastpd(mymodel);
		//fastpd.infer();
		////fastpd.arg(x);
		//val=fastpd.value();
		//printf("bound %f\n",fastpd.bound());
	 //   printf("FPD  minvalue %f\n",fastpd.value());

	} 
	else if (ia==4)
	{

		ICM icm(mymodel);

		icm.setStartingPoint(startingPoint.begin());
		icm.infer();
		icm.arg(x);
		val=icm.value();
		//printf("ICM minvalue %f\n",icm.value());

		//printf("Saving...\n");
		//opengm::hdf5::save(mymodel, "gm.h5", "mygm");

			      //   Parameter(): useKovtunsMethod_(true), probing_(false),  strongPersistency_(false), rounds_(0), permutationType_(NONE) {};
         //std::vector<LabelType> label_;
         //bool useKovtunsMethod_;
         //const bool probing_; //do not use this!
         //bool strongPersistency_;
         //size_t rounds_;
         //PermutationType permutationType_;

	} 	
	else if (ia==5)
	{
		if ((int)mxGetN(prhs[14])!=3)
			mexErrMsgTxt("3 Parameters needed for TRWSi\n");

		TRWSi::Parameter param;
		param.maxNumberOfIterations_=(size_t)inf_params[0];
		param.fastComputations_=(bool)inf_params[1];
		param.absolutePrecision_=(bool)inf_params[2];

		TRWSi trwsi(mymodel,param);		
		trwsi.setStartingPoint(startingPoint.begin());
		trwsi.infer();
		trwsi.arg(x);
		val=trwsi.value();
	}
	else if (ia==6)
	{
		if ((int)mxGetN(prhs[14])!=3)
			mexErrMsgTxt("3 Parameters needed for TRWSi\n");

		TRWSi::Parameter param;
		param.maxNumberOfIterations_=(size_t)inf_params[0];
		param.fastComputations_=(bool)inf_params[1];
		param.absolutePrecision_=(bool)inf_params[2];

		TRWSi trwsi(mymodel,param);		
		trwsi.setStartingPoint(startingPoint.begin());
		trwsi.infer();
		trwsi.arg(x);
		val=trwsi.value();	
		// if ((int)mxGetN(prhs[14])!=1)
			// mexErrMsgTxt("1 Parameter needed for FastPD\n");

		// FastPD::Parameter param;
		// param.numberOfIterations_ = (size_t)inf_params[0];

		// opengm::hdf5::save(mymodel, "gm.h5", "mygm");
		// FastPD fastpd(mymodel, param);		
		// fastpd.setStartingPoint(startingPoint.begin());
		// fastpd.infer();
		// fastpd.arg(x);
		// val=fastpd.value();
	}
	else
	{
		printf("Unknown inference algorithm");
	}
	// retreive optimal labeling
	for(size_t j = 0; j < x.size(); ++j) {
		labelingoutput[j]=(double)x[j];
	}
	*Eval = val;

	
 	   
}
