/**********************************/
/* Definition of Matlab's gateway */
/**********************************/

#ifdef createGateway

MEXfunction tmpC_testTprod2
Cfunction tmpC_testTprod2_raw
include testTprod2_raw.c

inputs
      double A[N,N]

outputs
      double plus_26[N-1,N]
      double plus_30[N-2,N]
      double tprod_37[N,N]


preprocess(N)

      code=csparse();

      eye_21=Teye([N,N]);
      plus_26=eye_21(2:end,:)-eye_21(1:end-1,:)
      plus_30=plus_26(2:end,:)-plus_26(1:end-1,:)
      tprod_37=tprod(plus_30,[-1,1],plus_30,[-1,2])
      
      declareGet(code,full(plus_26),'get_plus_26');
      declareGet(code,full(plus_30),'get_plus_30');
      declareGet(code,full(tprod_37),'get_tprod_37');
      

      compile2C(code,inf,inf,'tmpC_testTprod2_CS.c');
      code
      
#endif

#include "tmpC_testTprod2_CS.c"

void tmpC_testTprod2_raw(const double *A,
			 double *plus_26,
			 double *plus_30,
			 double *tprod_37,
			 mwSize N)
{

  clock_t dt=clock();

  get_plus_26(plus_26);
  get_plus_30(plus_30);
  get_tprod_37(tprod_37);

  dt=clock()-dt;
  mexPrintf("  testTprod2_raw: %.1lf us\n",dt/(double)CLOCKS_PER_SEC*1e6);

}
		  
		  
