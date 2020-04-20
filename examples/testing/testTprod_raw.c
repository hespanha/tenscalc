/**********************************/
/* Definition of Matlab's gateway */
/**********************************/

#ifdef createGateway

MEXfunction tmpC_testTprod
Cfunction tmpC_testTprod_raw
include testTprod_raw.c

inputs
      double A[n1,n2,n3]
      double B[n1,n2,n3]
      double C[n1,n2,n3]

outputs
      double D[n2,n3]


preprocess(n1,n2,n3)

      code=csparse();

      Tvariable A [n1,n2,n3];
      Tvariable B [n1,n2,n3];
      Tvariable C [n1,n2,n3];

      declareSet(code,A,'setA');
      declareSet(code,B,'setB');
      declareSet(code,C,'setC');
      ind=[-1,1,2];
      declareGet(code,tprod(A,ind,B,ind,C,ind),'getTprod');
      
      compile2C(code,inf,inf,'tmpC_testTprod_CS.c');
      
#endif

#include "tmpC_testTprod_CS.c"

void tmpC_testTprod_raw(const double *A,
			const double *B,
			const double *C,
			double *D,
			mwSize n1,
			mwSize n2,
			mwSize n3)
{

  setA(A);
  setB(B);
  setC(C);

  clock_t dt=clock();

  getTprod(D);

  dt=clock()-dt;
  mexPrintf("  testTprod_raw: %.1lf us\n",dt/(double)CLOCKS_PER_SEC*1e6);

}
		  
		  
