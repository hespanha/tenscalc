/**********************************/
/* Definition of Matlab's gateway */
/**********************************/

#ifdef createGateway

MEXfunction tmpC_testSet
Cfunction tmpC_testSet_raw
include testSet_raw.c

inputs
      double A[n1,n2,n3]
      double B[n1,n2,n3]
      double C[n1,n2,n3]

outputs
      double D[n1*2,n2,n3]


preprocess(n1,n2,n3)

      code=csparse();

      Tvariable A [n1,n2,n3];
      Tvariable B [n1,n2,n3];
      Tvariable C [n1,n2,n3];

      Tvariable ABC [3*n1,n2,n3];

      AB=ABC(1:2*n1,:,:);
      BC=ABC(n1+1:end,:,:);

      declareSet(code,ABC(1:n1,:,:),'setA');
      declareSet(code,ABC(n1+1:2*n1,:,:),'setB');
      declareSet(code,ABC(2*n1+1:3*n1,:,:),'setC');
      declareGet(code,-AB+BC,'getSum');

      fprintf('After Set/Get declarations\n');
      disp(code)

      compile2C(code,inf,inf,'tmpC_testSet_CS.c');
      
      fprintf('After Compile2C\n');
      disp(code)
#endif

#include "tmpC_testSet_CS.c"
	
void tmpC_testSet_raw(const double *A,
		      const double *B,
		      const double *C,
		      double *D,
		      mwSize n1,
		      mwSize n2,
		      mwSize n3)
{
  clock_t dt=clock();

  setA(A);
  setB(B);
  setC(C);

  getSum(D);

  dt=clock()-dt;
  mexPrintf("      testset_raw (sets+compute+get): %.1lf us\n",dt/(double)CLOCKS_PER_SEC*1e6);

}
		  
		  
