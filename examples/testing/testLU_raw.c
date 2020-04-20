/**********************************/
/* Definition of Matlab's gateway */
/**********************************/

#ifdef createGateway

MEXfunction tmpC_testLU
Cfunction tmpC_testLU_raw
include testLU_raw.c

inputs
      double A[n,n]
      double B[n,m]

outputs
      double LU[n,n]
      double X[n,m]


preprocess(n,m)

      code=csparse();

      Tvariable A [n,n];
      Tvariable B [n,m];

      declareSet(code,A,'setA');
      declareSet(code,B,'setB');
      declareGet(code,lu(A),'getLU');
      declareGet(code,lu(A)\B,'getX');
      
      fprintf('After Set/Get declarations\n');
      disp(code)

	compile2C(code,inf,inf,'tmpC_testLU_CS.c','tmpC_testLU_CS.h','tmpC_testLU_CS.log','',0);
      
      fprintf('After Compile2C\n');
      disp(code)
#endif

#include "tmpC_testLU_CS.c"
      
void tmpC_testLU_raw(const double *A,
		     const double *B,
		     double *LU,
		     double *X,
		     mwSize n,
		     mwSize m)
{
  clock_t dt;
  
  setA(A);
  setB(B);

  dt=clock();
  getLU(LU);
  dt=clock()-dt;
  mexPrintf("  getLU: %.1lf us\n",dt/(double)CLOCKS_PER_SEC*1e6);

  dt=clock();
  getX(X);
  dt=clock()-dt;
  mexPrintf("  getX: %.1lf us\n",dt/(double)CLOCKS_PER_SEC*1e6);

}
		  
		  
