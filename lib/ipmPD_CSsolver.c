/*
  This file is part of Tencalc.

  Copyright (C) 2012-21 The Regents of the University of California
  (author: Dr. Joao Hespanha).  All rights reserved.
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
/* ipmPD_CSsolver() - Basic iteration of an interior point method             */
/******************************************************************************/

/* Functions to do the basic computations, typically generated by ipmPD_CS.m  */
extern void initPrimalDual__();
extern void initDualEq__();
extern void initDualIneq__();
extern void setAddEye2Hessian1__(const double *addEye2Hessian);
extern void setAddEye2Hessian2__(const double *addEye2Hessian);
extern void updatePrimalDual__();
extern void scaleIneq__();
extern void scaleCost__();
extern void getScale4Cost__(double *scale4Cost);
extern void setAlphaPrimal__(const double *alpha);
extern void setAlphaDualEq__(const double *alpha);
extern void setAlphaDualIneq__(const double *alpha);
extern void setMu__(const double *mu);
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

#if verboseLevel>=4
#define printf4(...) mexPrintf(__VA_ARGS__)
#else
#define printf4(...)
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

#if (useLDL != 0)
#define useInertia 1 // set to 1 to use inertia
#else
#define useInertia 0 // cannot use inertia without LDL factorization
#endif

EXPORT void ipmPD_CSsolver(
	       /* inputs */
	       double  *mu0,
	       int32_t *maxIter,
	       int32_t *saveIter,
	       double  *addEye2Hessian,
	       /* outputs */
	       int32_t *status,
	       int32_t *iter,
	       double  *time
	       )
{
  *iter=0; // iteration number

#if setAddEye2Hessian != 0
  double addEye2Hessian1=addEye2Hessian[0];
  double addEye2Hessian2=addEye2Hessian[1];
  setAddEye2Hessian1__(&addEye2Hessian1);
  setAddEye2Hessian2__(&addEye2Hessian2);
  double derr,curvature;
  int updateaAddEye2Hessian1=0,updateaAddEye2Hessian2=0;
#if (useInertia != 0)
  double mp,mn;
#endif
#endif

#define mpDesired nU
#if smallerNewtonMatrix != 0
#define mnDesired nG
#else
#define mnDesired (nF+nG)
#endif

  double norminf_grad,alphaMax_=alphaMax;

#if nF>0
  double mu=*mu0,muMin,gap,ineq,ineq1,dual,
    scale4cost,desiredDualityGapVar=desiredDualityGap,alphaPrimal=0,alphaDualEq=0,alphaDualIneq=0;
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

  initPrimal__();

#if scaleCost != 0
  scaleCost__();
#if nF>0
  getScale4Cost__(&scale4cost);
  desiredDualityGapVar *= scale4cost;
#endif
#endif

#if nF>0
#if scaleInequalities != 0
  scaleIneq__();
#endif
  setMu__(&mu);
  muMin=desiredDualityGapVar/nF/2;
#endif

  printf2("%s.c (coupledAlphas=%d,skipAffine=%d,delta=%g,allowSave=%d,addEye2Hessian=%d,adjustAddEye2Hessian=%d,muFactorAggresive=%g,muFactorConservative=%g,LDL=%d,umfpack=%d):\n   %d primal variables, %d eq. constr., %d ineq. constr.\n",
	  __FUNCTION__,coupledAlphas,skipAffine,(double)delta,allowSave,setAddEye2Hessian,adjustAddEye2Hessian,muFactorAggressive,muFactorConservative,useLDL,useUmfpack,nU,nG,nF);
#if verboseLevel>=3
  char *header="Iter      cost   |grad|   |eq|    ineq.    dual    gap   l(mu) "
#if (setAddEye2Hessian != 0)
    "l(Hp) l(Hd) "
#endif
#if (adjustAddEye2Hessian != 0)
#if (useInertia != 0)
    "eig+ eig-  d.err. "
#else
    "curvature  d.err. "
#endif
#endif
    "alphaA  sigma   alphaP  alphaDI alphaDE      time\n";
  printf3(header);
#endif
#if nF>0
  printf3("%4d:<-mx des.->%8.1e%8.1e                %8.1e%6.1f",
	  *maxIter,gradTolerance,equalTolerance,desiredDualityGapVar,log10(muMin));
#else
  printf3("%4d:<-mx tol.->%8.1e%8.1e                              ",*maxIter,gradTolerance,equalTolerance);
#endif
#if (setAddEye2Hessian != 0) && (adjustAddEye2Hessian != 0)
  printf3("%6.1f",log10(addEye2Hessian1tolerance));
#if (useInertia != 0)
  printf3("      %5d%5d",mpDesired,mnDesired);
#endif
#endif
  printf3("\n");

#if verboseLevel>=1
  clock_t dt0=clock();
#endif
#if verboseLevel>=3
  clock_t dt1;
#endif

#if nF>0
  initDualIneq__();
#endif

#if nG>0
  initDualEqX__();
  //initDualEq__();
#endif

  while (1) {
    (*iter)++;
#if verboseLevel>=3
    if ((*iter) % 50 ==0)
      printf3(header);
    dt1=clock();
    printf3("%4d:",*iter);
#endif

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
    printf3("%8.1e",norminf_grad);

    if (isnan(norminf_grad)) {
	printf2("  -> failed to invert hessian\n");
	(*status) = 4;
#if allowSave!=0
	printf("Saving \"" saveNamePrefix "_WW.values\" due to status = 4 (grad)\n");
	saveWW__(saveNamePrefix "_WW.values");
#endif
	break;
      }

#if nG>0
    getNorminf_G__(&norminf_eq);
    printf3("%8.1e",norminf_eq);
#else
    printf3("  -eq-  ");
#endif
#if nF>0
    getGapMinFMinLambda__(&gap,&ineq,&dual);
    printf3("%8.1e%8.1e%8.1e",ineq,dual,gap);
    if (ineq<=0) {
	printf2("  -> (primal) variables violate constraints\n");
	(*status) = 1;
	break;
    }
    if (dual<=0) {
	printf2("  -> negative value for dual variables\n");
	(*status) = 2;
	break;
    }
#else
    printf3(" -ineq-  -dual-   -gap- ");
#endif

    if (norminf_grad<=gradTolerance
#if nF>0
	&& gap<=desiredDualityGapVar
#endif
#if nG>0
	&& norminf_eq<=equalTolerance
#endif
#if (setAddEye2Hessian != 0) && (adjustAddEye2Hessian != 0)
	&& addEye2Hessian1<=addEye2Hessian1tolerance
#endif
	 ) {
	       printf2("  -> clean exit\n");
	       (*status) = 0;
	       break;
    }

#if nF>0
    printf3("%6.1f",log10(mu));
#else
    printf3(" -mu- ");
#endif

    /*************************/
    /* Adjust addEye2Hessian */
    /*************************/

#define addEye2HessianMAX 1e2
#define addEye2HessianMIN 1e-20
#define maxDirectionError 1e-9

#if (setAddEye2Hessian != 0) && (adjustAddEye2Hessian != 0)

    // updates delayed from the last iteration
    if (updateaAddEye2Hessian1) {
      setAddEye2Hessian1__(&addEye2Hessian1);
      updateaAddEye2Hessian1=0;
    }
    if (updateaAddEye2Hessian2) {
      setAddEye2Hessian2__(&addEye2Hessian2);
      updateaAddEye2Hessian2=0;
    }

    getCurvature__(&curvature);
#if (useInertia != 0)
    if (curvature<=0) {
      getHessInertia__(&mp,&mn);
    }
#endif // #if (useInertia != 0)
    getDirectionError__(&derr);
    if (curvature>0
#if (useInertia != 0)
	|| (mp==mpDesired && mn==mnDesired)
#endif // #if (useInertia != 0)
	) {
#if (useInertia == 0)
      printf3("%6.1f%6.1f %+8.1e %8.1e",log10(addEye2Hessian1),log10(addEye2Hessian2),curvature,derr);
#else
#if verboseLevel>=3
      if (curvature>0) {
	printf3("%6.1f%6.1f c=%+6.0e %8.1e",log10(addEye2Hessian1),log10(addEye2Hessian2),curvature,derr);
      } else {
	printf3("%6.1f%6.1f%5.0f%5.0f%8.1e",log10(addEye2Hessian1),log10(addEye2Hessian2),mp,mn,derr);
      }
#endif
#endif // #if (useInertia == 0)
      if (addEye2Hessian1>addEye2HessianMIN && derr<maxDirectionError) {
      //if (addEye2Hessian1>addEye2HessianMIN && norminf_grad<=10*gradTolerance) {
	addEye2Hessian1=MAX(.75*addEye2Hessian1,addEye2HessianMIN);
	updateaAddEye2Hessian1=1; // update at next iteration
      }
      if (addEye2Hessian2>addEye2HessianMIN && derr<maxDirectionError) {
	addEye2Hessian2=MAX(.75*addEye2Hessian2,addEye2HessianMIN);
	updateaAddEye2Hessian2=1; // update at next iteration
      }
    } else {
      for (int ii=0;ii<20;ii++) {
	int change=0;
	if (curvature<0)
#if (useInertia == 0)
	  {
	    // inertia not available, increase both
	    printf4("%6.1f%6.1f %+8.1e %8.1e\n                                                              ",
		    log10(addEye2Hessian1),log10(addEye2Hessian2),curvature,derr);
	    if (addEye2Hessian1<addEye2HessianMAX) {
	      addEye2Hessian1=MIN(10*MAX(addEye2Hessian1,addEye2HessianMIN),addEye2HessianMAX);
	      setAddEye2Hessian1__(&addEye2Hessian1);
	      change=1;
	    }
	    if (addEye2Hessian2<addEye2HessianMAX) {
	      addEye2Hessian2=MIN(10*MAX(addEye2Hessian2,addEye2HessianMIN),addEye2HessianMAX);
	      setAddEye2Hessian2__(&addEye2Hessian2);
	      change=1;
	    }
	  }
#else
	{
	  if (mp<mpDesired) {
	    // not enough positive eigenvales
	    printf4("%6.1f*%5.1f%5.0f*%4.0f%8.1e\n                                                              ",
		    log10(addEye2Hessian1),log10(addEye2Hessian2),mp,mn,derr);
	    if (addEye2Hessian1<addEye2HessianMAX) {
	      addEye2Hessian1=MIN(10*MAX(addEye2Hessian1,addEye2HessianMIN),addEye2HessianMAX);
	      setAddEye2Hessian1__(&addEye2Hessian1);
	      change=1;
	    }
	    if (addEye2Hessian2<addEye2HessianMAX) {
	      addEye2Hessian2=MIN(2*MAX(addEye2Hessian2,addEye2HessianMIN),addEye2HessianMAX);
	      setAddEye2Hessian2__(&addEye2Hessian2);
	      change=1;
	    }
	  } else if (mn<mnDesired) {
	    // not enough negative eigenvales
	    printf4("%6.1f%6.1f*%4.0f%5.0f*%7.1e\n                                                              ",
		    log10(addEye2Hessian1),log10(addEye2Hessian2),mp,mn,derr);
	    if (addEye2Hessian1<addEye2HessianMAX) {
	      addEye2Hessian1=MIN(2*MAX(addEye2Hessian1,addEye2HessianMIN),addEye2HessianMAX);
	      setAddEye2Hessian1__(&addEye2Hessian1);
	      change=1;
	    }
	    if (addEye2Hessian2<addEye2HessianMAX) {
	      addEye2Hessian2=MIN(10*MAX(addEye2Hessian2,addEye2HessianMIN),addEye2HessianMAX);
	      setAddEye2Hessian2__(&addEye2Hessian2);
	      change=1;
	    }
	  }
	}
#endif // #if (useInertia == 0)
	if (! change)
	  break;

	getCurvature__(&curvature);
#if (useInertia != 0)
	if (curvature<=0) {
	  getHessInertia__(&mp,&mn);
	}
#endif // #if (useInertia != 0)
	getDirectionError__(&derr);

      }
#if (useInertia == 0)
      printf3("%6.1f%6.1f %+8.1e %8.1e",log10(addEye2Hessian1),log10(addEye2Hessian2),curvature,derr);
#else
#if verboseLevel>=3
      if (curvature>0)
	printf3("%6.1f%6.1f c=%+6.0e %8.1e",log10(addEye2Hessian1),log10(addEye2Hessian2),curvature,derr);
      else
	printf3("%6.1f%6.1f%5.0f%5.0f%8.1e",log10(addEye2Hessian1),log10(addEye2Hessian2),mp,mn,derr);
#endif
#endif // #if (useInertia != 0)
    }
#else // #if (setAddEye2Hessian != 0) && (adjustAddEye2Hessian != 0)
#if setAddEye2Hessian != 0
    printf3("%6.1f%6.1f",log10(addEye2Hessian1),log10(addEye2Hessian2));
#endif
#endif // #if (setAddEye2Hessian != 0) && (adjustAddEye2Hessian != 0)

#if nF==0
    /****************************************/
    /******  NO INEQUALITY CONSTRAINTS ******/
    /****************************************/
    setAlphaPrimal__(&alphaMax_);
#if nG>0
    setAlphaDualEq__(&alphaMax_);
#endif
    printf3("  -alpA- -sigm- ");
    printf3("%8.1e                ",alphaMax_);

#if allowSave!=0
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

#if skipAffine!=0
    printf3(" -alpA-  -sigm- ");
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
    printf3("%8.1e",alphaPrimal);

    // update mu based on sigma, but this only seems to be safe for:
    // 1) "long" newton steps in the affine direction
    // 2) equality constraints fairly well satisfied (perhaps not very important)
    if (alphaPrimal> alphaMax/2
#if nG>0
	&& (norminf_eq<100*equalTolerance || norminf_eq<1e-3)
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
      printf3("%8.1e",sigma);
      mu=sigma*gap/nF;
      if (mu < muMin) mu=muMin;
      setMu__(&mu);
    } else {
      printf3(" -sigm- ");
    }
#endif  // skipAffine!=0

    /*******************************************************************/
    /** Combined search direction                                     **/
    /*******************************************************************/

#if allowSave!=0
    if ((*iter)==(*saveIter)) {
      printf("Saving \"" saveNamePrefix "_WW.values\" due to iter = saveIter\n");
      saveWW__(saveNamePrefix "_WW.values");
    }
#endif

    getMaxAlphas_s__(&alphaPrimal,&alphaDualIneq);

#if coupledAlphas!=0
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
	  printf2("  -> failed to invert hessian\n");
	  (*status) = 4;
#if allowSave!=0
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

#if coupledAlphas!=0
    alphaDualEq=alphaPrimal;
    alphaDualIneq=alphaPrimal;
#else
    alphaDualIneq *= STEPBACK;
    if (alphaDualIneq>alphaMax_)
      alphaDualIneq=alphaMax_;
    //alphaDualEq = alphaMax;
    alphaDualEq = alphaDualIneq;
#endif
#if nG>0
    setAlphaDualEq__(&alphaDualEq);
#endif
    setAlphaDualIneq__(&alphaDualIneq);
    updatePrimalDual__();

#if nG>0
    printf3("%8.1e%8.1e%8.1e",alphaPrimal,alphaDualIneq,alphaDualEq);
#else
    printf3("%8.1e%8.1e  -eq-  ",alphaPrimal,alphaDualIneq);
#endif

#if skipAffine!=0
    // More aggressive if
    // 1) "long" newton steps in the affine direction
    // 2) small gradient
    // 3) equality constraints fairly well satisfied
    // (2+3 mean close to the central path)
    //int th_grad=norminf_grad<MAX(1e-1,1e2*gradTolerance);
    int th_grad=norminf_grad<MAX(1e-6,1e0*gradTolerance);
#if nG>0
    //int th_eq=norminf_eq<MAX(1e-3,1e2*equalTolerance);
    int th_eq=norminf_eq<MAX(1e-5,1e0*equalTolerance);
#endif
    if (alphaPrimal>alphaMax_/2 && th_grad
#if nG>0
	&& th_eq
#endif
	) {
      //mu *= muFactorAggressive;
      mu *= MIN(muFactorAggressive,pow(mu,.5));
      if (mu < muMin) mu=muMin;
      setMu__(&mu);
      printf3(" * ");
    } else {
      if (alphaPrimal<.1) {
	//mu = MIN(1e2,1.25*mu);
	mu *= 1.1;
	if (mu > *mu0) mu = *mu0;
	setMu__(&mu);
	initDualIneq__();
	printf3("^");
      } else if (alphaPrimal>.99
#if nG>0
	&& th_eq
#endif
		 ) {
	mu *= muFactorConservative;
	if (mu < muMin) mu = muMin;
	setMu__(&mu);
	printf3("v");
      } else {
	printf3("=");
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
#else
	printf3("   ");
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

#if allowSave!=0
  if ((*saveIter)==0 && (*status)==0) {
      printf("  Saving \"" saveNamePrefix "_WW.values\" due to saveIter = 0\n");
      saveWW__(saveNamePrefix "_WW.values");
    }
#endif

  if ((*status)==8) {
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
    if (gap>desiredDualityGapVar) {
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

  if (*status) {
    char sep='(';
    printf2("%3d:status=0x%X ",(*iter),(*status));
    if ((*status) & 16) {
      printf2("%clarge gradient",sep);
      sep=',';
    }
    if ((*status) & 32) {
      printf2("%cbad equality const.",sep);
      sep=',';
    }
    if ((*status) & 64) {
      printf2("%clarge duality gap",sep);
      sep=',';
    }
    if ((*status) & 128) {
      printf2("%clarge mu",sep);
      sep=',';
    }
    if ((*status) & 256) {
      printf2("%calpha negligible",sep);
    } else if ((*status) & 512) {
      printf2("%calpha<.1",sep);
    } else if ((*status) & 1024) {
      printf2("%calpha<.5",sep);
    }
    printf2(")\n                ");
  } else {
    printf2("%4d:status=0x%X, ",(*iter),(*status));
  }
  printf2("cost=%13.5e",J);
  printf2(", |grad|=%10.2e",norminf_grad);
#if setAddEye2Hessian != 0
  printf2(", addEye2Hessian=[%10.2e,%10.2e]",addEye2Hessian1,addEye2Hessian2);
#endif
#if nG>0
  printf2(", |eq|=%10.2e",norminf_eq);
#endif
#if nF>0
  printf2(", ineq=%10.2e,\n                 dual=%10.2e, gap=%10.2e, last alpha=%10.2e",
	  ineq,dual,gap,alphaPrimal);
#endif
  printf2(" (%.1lfus,%.2lfus/iter)\n",(*time)*1e6,(*time)/(double)(*iter)*1e6);
#endif  // verboseLevel>=2
}
