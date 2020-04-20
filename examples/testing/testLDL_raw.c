/**********************************/
/* Definition of Matlab's gateway */
/**********************************/

#ifdef createGateway

MEXfunction tmpC_testLDL
Cfunction tmpC_testLDL_raw
include testLDL_raw.c

inputs
      double A[n,n]
      double B[n,m]

outputs
      double LDL[n,n]
      double X[n,m]


preprocess(n,m)

      code=csparse();

      Tvariable A [n,n];
      Tvariable B [n,m];

      declareSet(code,A,'setA');
      declareSet(code,B,'setB');
      declareGet(code,ldl(A),'getLDL');
      declareGet(code,ldl(A)\B,'getX');
      
      fprintf('After Set/Get declarations\n');
      disp(code)

      compile2C(code,inf,inf,'tmpC_testLDL_CS.c','tmpC_testLDL_CS.h','tmpC_testLDL_CS.log','',0);
      
      fprintf('After Compile2C\n');
      disp(code)
#endif

#include "tmpC_testLDL_CS.c"

void tmpC_testLDL_raw(const double *A,
		      const double *B,
		      double *LDL,
		      double *X,
		      mwSize n,
		      mwSize m)
{
  clock_t dt;

  setA(A);
  setB(B);

  dt=clock();
  getLDL(LDL);
  dt=clock()-dt;
  mexPrintf("  getLDL: %.1lf us\n",dt/(double)CLOCKS_PER_SEC*1e6);
  
  dt=clock();
  getX(X);
  dt=clock()-dt;
  mexPrintf("  getX: %.1lf us\n",dt/(double)CLOCKS_PER_SEC*1e6);
  
}
		  
		  
