/**********************************/
/* Definition of Matlab's gateway */
/**********************************/

#ifdef createGateway

MEXfunction tmpC_testLUatomic
Cfunction tmpC_testLUatomic_raw
include testLUatomic_raw.c

inputs
      double A[n,n]
      double B[n,1]

outputs
      double X[n,1]


preprocess(n)

      code=csparse();

      Tvariable A [n,n];
      Tvariable B [n];

      declareSet(code,A,'setA');
      declareSet(code,B,'setB');
      luA=declareAlias(code,lu(A),'luA',true);
      declareGet(code,luA\B,'getX');
      
      fprintf('After Set/Get declarations\n');
      disp(code)

      compile2C(code,inf,inf,'tmpC_testLU_CS.c','tmpC_testLU_CS.h','tmpC_testLU_CS.log','',0);
      
      fprintf('After Compile2C\n');
      disp(code)
#endif

#include "tmpC_testLU_CS.c"
      
void tmpC_testLUatomic_raw(const double *A,
			   const double *B,
			   double *X,
			   mwSize n)
{
  clock_t dt;
  
  setA(A);
  setB(B);

  dt=clock();
  getX(X);
  dt=clock()-dt;
  mexPrintf("  getX: %.1lf us\n",dt/(double)CLOCKS_PER_SEC*1e6);
}
		  
		  
