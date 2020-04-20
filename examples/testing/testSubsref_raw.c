/**********************************/
/* Definition of Matlab's gateway */
/**********************************/

#ifdef createGateway

MEXfunction tmpC_testSubsref
Cfunction tmpC_testSubsref_raw
include testSubsref_raw.c

inputs
      double X[N]
      double dt1[N-1]

outputs
      double DX[N-1]
      double DDX[N-2]


preprocess(N)

      code=csparse();

      Tvariable X N;  
      Tvariable dt1 N-1;

      DX=(X(2:end)-X(1:end-1)).*dt1;             % velocity
      DDX=(DX(2:end)-DX(1:end-1)).*dt1(1:end-1); % acceleration

      declareSet(code,X,'setX');
      declareSet(code,dt1,'setDt1');
      declareGet(code,DX,'getDX');
      declareGet(code,DDX,'getDDX');

      fprintf('After Set/Get declarations\n');
      disp(code)

      compile2C(code,inf,inf,'tmpC_testSubsref_CS.c','tmpC_testSubsref_CS.h','tmpC_testSubsref_CS.log','',0);
      
      fprintf('After Compile2C\n');
      disp(code)
#endif

#include "tmpC_testSubsref_CS.c"

void tmpC_testSubsref_raw(const double *X,
			  const double *dt1,
			  double *DX,
			  double *DDX,
			  mwSize N)
{

  setX(X);
  setDt1(dt1);

  clock_t dt=clock();

  getDX(DX);
  getDDX(DDX);

  dt=clock()-dt;
  mexPrintf("  testsum_raw: %.1lf us\n",dt/(double)CLOCKS_PER_SEC*1e6);

}
		  
		  
