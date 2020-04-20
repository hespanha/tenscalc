/**********************************/
/* Definition of Matlab's gateway */
/**********************************/

#ifdef createGateway

MEXfunction tmpC_testSum
Cfunction tmpC_testSum_raw
include testSum_raw.c

inputs
      double A[n1,n2,n3]
      double B[n1,n2,n3]
      double C[n1,n2,n3]

outputs
      double D[n1,n2,n3]


preprocess(n1,n2,n3)

      code=csparse();

      Tvariable A [n1,n2,n3];
      Tvariable B [n1,n2,n3];
      Tvariable C [n1,n2,n3];

      declareSet(code,A,'setA');
      declareSet(code,B,'setB');
      declareSet(code,C,'setC');
      declareGet(code,-A+B+C,'getSum');

      fprintf('After Set/Get declarations\n');
      disp(code)

	compile2C(code,inf,inf,'tmpC_testSum_CS.c');
      
      fprintf('After Compile2C\n');
      disp(code)
#endif

#include "tmpC_testSum_CS.c"

void tmpC_testSum_raw(const double *A,
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
  mexPrintf("      testsum_raw (sets+compute+get): %.1lf us\n",dt/(double)CLOCKS_PER_SEC*1e6);

}
		  
		  
