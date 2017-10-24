/*
  Copyright 2012-2017 Joao Hespanha

  This file is part of Tencalc.

  TensCalc is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  TensCalc is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with TensCalc.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>

#ifdef __APPLE__
#include <unistd.h>
#elif __linux__
#include <unistd.h>
#elif _WIN32
#include <windows.h>
#endif

/******************************************************************************/
/* ipmPD_CSsolver() - Basic iteration of an interior point method                   */
/******************************************************************************/

/* Functions to do the basic computations, typically generated by ipmPD_CS.m  */
extern void initPrimalDual__();
extern void initDualEq__();
extern void initDualIneq__();
extern void updatePrimalDual__();
extern void setAlphaPrimal__(double *alpha);
extern void setAlphaDualEq__(double *alpha);
extern void setAlphaDualIneq__(double *alpha);
extern void setMu__(double *mu);
extern void getJ__(double *J);
extern void getNorminf_G__(double *norminf_eq);
extern void getNorminf_Grad__(double *norminf_grad);
extern void getGapMinFMinLambda__(double *gap,double *ineq,double *dual);
extern void getMaxAlphas_a__(double *alphaPrimal,double *alphaDual);
extern void getMaxAlphas_s__(double *alphaPrimal,double *alphaDual);
extern void getMinF_a__(double *ineq);
extern void getMinF_s__(double *ineq);
extern void getRho__(double *sigma);
extern void saveWW__(char *filename);

#ifdef DYNAMIC_LIBRARY
#ifdef __APPLE__
#define EXPORT __attribute__((visibility("default")))
#elif __linux__
#define EXPORT __attribute__((visibility("default")))
#elif _WIN32
#define EXPORT __declspec(dllexport)
#endif
#endif

//#define DEBUG

/*****************/
/* Main function */
/*****************/

#if verboseLevel>=2
#include <mex.h>
#define printf2(...) mexPrintf(__VA_ARGS__)
#else
#define printf2(...) 
#endif

#if verboseLevel>=3
#define printf3(...) mexPrintf(__VA_ARGS__)
#else
#define printf3(...) 
#endif

#define MIN(A,B) (((A)<(B))?(A):(B))
#define MAX(A,B) (((A)>(B))?(A):(B))

#ifdef DEBUG
void printMatrix(const char *name,double *mat,int m, int n)
{
  int i,j,nnz=0;
  if (0) {
    printf("%s[%d,%d] =\n",name,m,n);
    for (i=0;i<m;i++) {
      printf("%2d: ",i);
      for (j=0;j<n;j++) {
	printf("%10g ",mat[i+m*j]);
	if (fabs(mat[i+m*j])>1e-7) nnz++;
      }
      printf("\n"); }
  } else {
    printf("\n%s =[\n",name);
    for (i=0;i<m;i++) {
      for (j=0;j<n;j++) {
	printf("%g,",mat[i+m*j]);
	if (fabs(mat[i+m*j])>1e-7) nnz++;
      }
      printf(";"); }
    printf("]; %% (nnz=%d)\n",nnz);
  }
}
#endif

EXPORT void ipmPD_CSsolver(
	       /* inputs */
	       double  *mu0,
	       int32_t *maxIter,
	       int32_t *saveIter,
	       /* outputs */
	       int32_t *status,
	       int32_t *iter,
	       double  *time
               )
{
  *iter=0; // iteration number

#ifdef DEBUG
  double buffer[(nU+nF+nG)*(nU+nF+nG)];
#endif
  
  double norminf_grad,alphaMax_=alphaMax;

#if nF>0
  double mu=*mu0,muMin=desiredDualityGap/nF/2,gap,ineq,ineq1,dual,
    alphaPrimal=0,alphaDualEq=0,alphaDualIneq=0;
#if skipAffine != 1
  double sigma;
#endif
#endif 

#if nG>0
  double norminf_eq;
#endif

#if verboseLevel>=2
  double J;
#endif

  printf2("%s.c (skipAffine=%d,delta=%g,allowSave=%d): %d primal variables, %d eq. constr., %d ineq. constr.\n",__FUNCTION__,skipAffine,(double)delta,allowSave,nU,nG,nF);
  printf3("Iter   cost      |grad|      |eq|    inequal     dual      gap       mu       alphaA    sigma     alphaP     alphaI     alphaE   time[us]\n");
#if nF>0
  printf3("%3d:<-mx tol-> %10.2e%10.2e                    %10.2e%10.2e\n",*maxIter,gradTolerance,equalTolerance,desiredDualityGap,muMin);
#else
  printf3("%3d:<-mx tol-> %10.2e%10.2e\n",*maxIter,gradTolerance,equalTolerance);
#endif
  
#if verboseLevel>=1
  clock_t dt0=clock();
#endif
#if verboseLevel>=3
  clock_t dt1;
#endif

  //initPrimalDual__();
  initPrimal__();

#if nF>0
  //getGap__(&gap);
  //if (mu < gap/nF/10) mu = gap/nF/10;
  setMu__(&mu);
#endif

#if nG>0
  initDualEq__();
#endif
#if nF>0
  initDualIneq__();
#endif

  while (1) {
#if verboseLevel>=3
    dt1=clock();
#endif

    (*iter)++;
    printf3("%3d:",*iter);

    if ((*iter) > (*maxIter)) {
      printf3("maximum # iterations (%d) reached.\n",*maxIter);
      (*status) = 8;
      break; }

    /*************************/
    /* Check exit conditions */
    /*************************/
#if verboseLevel>=3
    getJ__(&J);
    printf3("%11.3e",J);
#endif
    getNorminf_Grad__(&norminf_grad);
    printf3("%10.2e",norminf_grad);

    if (isnan(norminf_grad)) {
	printf3("  -> failed to invert hessian\n");
	(*status) = 4;
#if allowSave==1
	printf("Saving \"" saveNamePrefix "_WW.values\" due to status = 4\n");
	saveWW__(saveNamePrefix "_WW.values");
#endif
	break;
      }

#if nG>0
    getNorminf_G__(&norminf_eq);
#if verboseLevel>=3
    printf3("%10.2e",norminf_eq);
#endif
#else
    printf3("    -eq-  ");
#endif
#if nF>0
    getGapMinFMinLambda__(&gap,&ineq,&dual);
    printf3("%10.2e%10.2e%10.2e",ineq,dual,gap);
    if (ineq<=0) {
        printf3("  -> (primal) variables violate constraints\n");
        (*status) = 1;
        break;
    }
    if (dual<=0) {
        printf3("  -> negative value for dual variables\n");
        (*status) = 2;
        break;
    }
#else
    printf3("   -ineq-    -dual-    -gap-  ");
#endif

    if (norminf_grad<=gradTolerance
#if nF>0
        && gap<=desiredDualityGap
#endif
#if nG>0
        && norminf_eq<=equalTolerance
#endif
         ) {
               printf3("  -> clean exit\n");
               (*status) = 0;
               break;
    }

#if nF>0
    printf3("%10.2e",mu);
#else
    printf3("   -mu-   ");
#endif

#if nF==0
    /****************************************/
    /******  NO INEQUALITY CONSTRAINTS ******/
    /****************************************/
    setAlphaPrimal__(&alphaMax_);
    printf3(" -alphaA-  -sigma- ");
    printf3("%10.2e                   ",alphaMax_);

#if allowSave==1
    if ((*iter)==(*saveIter)) {
      printf("Saving \"" saveNamePrefix "_WW.values\" due to iter = saveIter\n");
      saveWW__(saveNamePrefix "_WW.values");
    }
#endif

    updatePrimalDual__();
#else

    /******************************************/
    /******  WITH INEQUALITY CONSTRAINTS ******/
    /******************************************/

#if skipAffine==1
    printf3(" -alphaA-  -sigma-");
#else
    /*******************************************************************/
    /** Affine search direction                                       **/
    /*******************************************************************/

    getMaxAlphas_a__(&alphaPrimal,&alphaDualIneq);

    alphaMax_ = (alphaPrimal<alphaMax)?alphaPrimal:alphaMax;
    if (alphaDualIneq<alphaMax_)
      alphaMax_=alphaDualIneq;

    if (alphaMax_ >= alphaMin) {
      // try max
      alphaPrimal=alphaMax_;
      setAlphaPrimal__(&alphaPrimal);getMinF_a__(&ineq);
      if (ineq<0) {
	// try min
        alphaPrimal=alphaMin;
	setAlphaPrimal__(&alphaPrimal);getMinF_a__(&ineq);
	if (ineq>0) {
	  // try between min and max
	  for (alphaPrimal = alphaMax_*.95 ; alphaPrimal >= alphaMin ; alphaPrimal /= 2) {
	    setAlphaPrimal__(&alphaPrimal);getMinF_a__(&ineq);
	    if (ineq>=0) {
	      break; }
	  }
	  if (alphaPrimal < alphaMin) {
	    alphaPrimal = 0;
	    setAlphaPrimal__(&alphaPrimal);
	  }
	} else {
	  alphaPrimal = 0;
	  setAlphaPrimal__(&alphaPrimal);
	}
      }
    } else {
      alphaPrimal = 0;
      setAlphaPrimal__(&alphaPrimal);
    }
    setAlphaDualIneq__(&alphaPrimal);
    printf3("%10.2e",alphaPrimal);

    // update mu based on sigma, but this only seems to be safe for:
    // 1) "long" newton steps in the affine direction
    // 2) equality constraints fairly well satisfied (perhaps not very important)
    if (alphaPrimal> alphaMax/2
#if nG>0
	&& norminf_eq<100*equalTolerance
#endif
        ) {
      getRho__(&sigma);
      if (sigma>1) sigma=1;
      if (sigma<0) sigma=0;
#if delta==2
      sigma=sigma*sigma;
#else
      sigma=sigma*sigma*sigma;
#endif
      printf3("%10.2e",sigma);
      mu = sigma*gap/nF;
      if (mu < muMin) mu = muMin;
      setMu__(&mu); 
    } else {
      printf3("  -sigma- ");
    }
#endif  // skipAffine==1

    /*******************************************************************/
    /** Combined search direction                                     **/
    /*******************************************************************/

#ifdef DEBUG
    getWW__(buffer);
    //printMatrix("WW",buffer,nU+nG,nU+nG);
    printMatrix("WW",buffer,nU+nG+nF,nU+nG+nF);
    getb_s__(buffer);
    //printMatrix("b_s",buffer,nU+nG,1);
    printMatrix("b_s",buffer,nU+nG+nF,1);
    getU__(buffer);
    printMatrix("U",buffer,nU,1);
    getLambda__(buffer);
    printMatrix("Lambda",buffer,nF,1);
    getDx_s__(buffer);
    //printMatrix("dx",buffer,nU+nG,1);
    printMatrix("dx",buffer,nU+nG+nF,1);
#endif
    
#if allowSave==1
    if ((*iter)==(*saveIter)) {
      printf("Saving \"" saveNamePrefix "_WW.values\" due to iter = saveIter\n");
      saveWW__(saveNamePrefix "_WW.values");
    }
#endif

    getMaxAlphas_s__(&alphaPrimal,&alphaDualIneq);

#if coupledAlphas==1
    if (alphaDualIneq<alphaPrimal)
      alphaPrimal=alphaDualIneq;
#endif

#define STEPBACK .99

    alphaPrimal *= STEPBACK;

    alphaMax_ = (alphaPrimal<alphaMax)?alphaPrimal:alphaMax;

    if (alphaMax_ >= alphaMin) {
      // try max
      alphaPrimal=alphaMax_/STEPBACK;
      setAlphaPrimal__(&alphaPrimal);getMinF_s__(&ineq);
      if (isnan(ineq)) {
	  printf3("  -> failed to invert hessian\n");
	  (*status) = 4;
#if allowSave==1
	  printf("Saving \"" saveNamePrefix "_WW.values\" due to status = 4\n");
	  saveWW__(saveNamePrefix "_WW.values");
#endif
	  break;
	}
      if (ineq>0) {
	// recheck just to be safe in case not convex
        alphaPrimal *= STEPBACK;
	setAlphaPrimal__(&alphaPrimal);getMinF_s__(&ineq1);
      }
      if (ineq<=0 || ineq1<=ineq/10) {
        // try min
	alphaPrimal=alphaMin/STEPBACK;
	setAlphaPrimal__(&alphaPrimal);getMinF_s__(&ineq);
	if (ineq>0) {
          // try between min and max
	  for (alphaPrimal = alphaMax_*.95;alphaPrimal >= alphaMin;alphaPrimal /= 2) {
            setAlphaPrimal__(&alphaPrimal);getMinF_s__(&ineq);
	    if (ineq>0) {
	      // backtrace just a little
              alphaPrimal *= STEPBACK;
	      // recheck just to be safe in case not convex
	      setAlphaPrimal__(&alphaPrimal);getMinF_s__(&ineq1);
	      if (ineq1>ineq/10) {
		break; }
	    }
	  } 
	  if (alphaPrimal < alphaMin) {
	    alphaPrimal = 0;
	    setAlphaPrimal__(&alphaPrimal);
	  }
	} else {
	  alphaPrimal = 0;
	  setAlphaPrimal__(&alphaPrimal);
	}
      }
    } else {
      alphaPrimal = 0;
      setAlphaPrimal__(&alphaPrimal);
    }

#if coupledAlphas==1
    alphaDualEq=alphaPrimal;
    alphaDualIneq=alphaPrimal;
#else
    alphaDualIneq *= STEPBACK;
    if (alphaDualIneq>alphaMax_)
      alphaDualIneq=alphaMax_;
    alphaDualEq = alphaMax;
#endif
    setAlphaDualEq__(&alphaDualEq);
    setAlphaDualIneq__(&alphaDualIneq);
    updatePrimalDual__();
    
#if nG>0
    printf3("%10.2e %10.2e %10.2e",alphaPrimal,alphaDualIneq,alphaDualEq);
#else
    printf3("%10.2e %10.2e   -eq-    ",alphaPrimal,alphaDualIneq);
#endif

#if skipAffine==1
    // More aggressive if
    // 1) "long" newton steps in the affine direction
    // 2) small gradient
    // 3) equality constraints fairly well satisfied
    // (2+3 mean close to the central path)
    int th_grad=norminf_grad<MAX(1e-1,1e2*gradTolerance);
#if nG>0
    int th_eq=norminf_eq<MAX(1e-3,1e2*equalTolerance);
#endif
    if (alphaPrimal>.5 && th_grad
#if nG>0
	&& th_eq
#endif
	) {
      mu *= muFactorAggressive;
      if (mu < muMin) mu = muMin;
      setMu__(&mu); 
      printf3(" *");
    } else {
      mu *= muFactorConservative;
      if (mu < muMin) mu = muMin;
      setMu__(&mu);
      if (alphaPrimal<.1) {
	mu = MIN(1e2,1.25*mu);
	setMu__(&mu);
	initDualIneq__();
      }
#if verboseLevel>=3
      if (th_grad)
	printf3("g");
      else
	printf3(" ");
#if nG>0
      if (th_eq)
	printf3("e");
      else
#endif
	printf3(" ");
#endif
    }
#endif
    
    // if no motion, slowly increase mu
    if (alphaPrimal<alphaMin && alphaDualIneq<alphaMin && alphaDualEq<alphaMin) {
      mu /= (muFactorConservative*muFactorConservative); // square to compensate for previous decrease
      if (mu < muMin) mu = muMin;
      setMu__(&mu); }

#endif
#if verboseLevel>=3
    dt1=clock()-dt1;
#endif
    printf3("%8.1lfus\n",dt1/(double)CLOCKS_PER_SEC*1e6);

  } // while(1)

#if allowSave==1
  if ((*saveIter)==0 && (*status)==0) {
      printf("  Saving \"" saveNamePrefix "_WW.values\" due to saveIter = 0\n");
      saveWW__(saveNamePrefix "_WW.values");
    }
#endif

  if ((*status) == 8) {
    getNorminf_Grad__(&norminf_grad);
    if (norminf_grad>gradTolerance) {
      (*status) |= 16;
    }
#if nG>0
    getNorminf_G__(&norminf_eq);
    if (norminf_eq>equalTolerance) {
      (*status) |= 32;
    }
#endif
#if nF>0
    getGapMinFMinLambda__(&gap,&ineq,&dual);
    if (gap>desiredDualityGap) {
      (*status) |= 64;
    }
    if (mu>muMin) {
      (*status) |= 128;
    }
    if (alphaPrimal<=alphaMin && alphaDualIneq<alphaMin && alphaDualEq<alphaMin) 
      (*status) |= 1792; // (256|512|1024);
    else if (alphaPrimal<=.1 && alphaDualIneq<.1 && alphaDualEq<.1)
      (*status) |= 1536; // (512|1024);
    else if (alphaPrimal<=.5 && alphaDualIneq<.5 && alphaDualEq<.5)
      (*status) |= 1024;
#endif
  }
  
#if verboseLevel>=1
  (*time)=(clock()-dt0)/(double)CLOCKS_PER_SEC;
#endif
    
#if verboseLevel>=2
  getJ__(&J);
  if ((*status)<8) {
    getNorminf_Grad__(&norminf_grad);
#if nG>0
    getNorminf_G__(&norminf_eq);
#endif
#if nF>0
    getGapMinFMinLambda__(&gap,&ineq,&dual);
#endif
  }
  
  printf2("%3d:status=0x%X, ",(*iter),(*status));
  printf2("cost=%13.5e",J);
#if nG>0
  printf2(", |eq|=%10.2e",norminf_eq);
#endif
#if nF>0
  printf2(", ineq=%10.2e,\n              dual=%10.2e, gap=%10.2e, last alpha=%10.2e",ineq,dual,gap,alphaPrimal);
#endif
  printf2(", |grad|=%10.2e",norminf_grad);
  printf2(" (%.1lfus,%.2lfus/iter)\n",(*time)*1e6,(*time)/(double)(*iter)*1e6);
#endif
}

