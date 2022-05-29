#include <stdio.h>
#include "umfpack.h"

double *null = (double *) NULL ;
void *Symbolic[1],*Numeric[1];
double Ax0[9];
double *Ax[1]={Ax0};
int Ap0[4]={0,3,6,9};
int *Ap[1]={Ap0};
int Ai0[9]={0,1,2,0,1,2,0,1,2};
int *Ai[1]={Ai0};
double b0[3];
double *b[1]={b0};
double x0[3];
double *x[1]={x0};

int main (void)
{
#ifndef BB
    Ax[0][0]=0;
    Ax[0][1]=1;
    Ax[0][2]=2;
    Ax[0][3]=3;
    Ax[0][4]=4;
    Ax[0][5]=5;
    Ax[0][6]=6;
    Ax[0][7]=4;
    Ax[0][8]=5;
    b[0][0]=9;
    b[0][1]=10;
    b[0][2]=12;

#endif
  (void)umfpack_di_symbolic(3,3,Ap[0],Ai[0],Ax[0],&Symbolic[0],null,null);
  (void)umfpack_di_numeric(Ap[0],Ai[0],Ax[0],Symbolic[0],&Numeric[0],null,null);
  umfpack_di_free_symbolic(&Symbolic[0]);

	{int atomicID=0;
	/* b[atomicID][0]=m[9]; */
	/* b[atomicID][1]=m[10]; */
	/* b[atomicID][2]=m[11]; */
	(void)umfpack_di_solve(UMFPACK_A,Ap[atomicID],Ai[atomicID],Ax[atomicID],x[atomicID],b[atomicID],Numeric,null,null);
	}
  (void) umfpack_di_solve (UMFPACK_A, Ap[0], Ai[0], Ax[0], x[0], b[0], Numeric[0], null, null) ;
  umfpack_di_free_numeric (&Numeric[0]) ;

  int n=sizeof(Ap0)/sizeof(int)-1;
  for (int i = 0 ; i < n ; i++) printf ("b [%d] = %g, x [%d] = %g\n", i, b[0][i],i, x[0] [i]) ;
  return (0) ;
}
