/* Created by script createGateway.m on 22-Jun-2017 21:15:49 */

/* START OF #included "GPL.c" */
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

    
    
/* END OF #included "GPL.c" */
/* mex -largeArrayDims -I/Users/hespanha/matlab_projects/tenscalc/TimDavis/SuiteSparse/include -L/Users/hespanha/matlab_projects/tenscalc/TimDavis/SuiteSparse/lib COPTIMFLAGS="-Ofast -msse -msse2 -msse3 -msse4 -msse4.1 -DNDEBUG" CFLAGS="\$CFLAGS -Wall -Werror -Wno-unused-variable -Wno-unused-result -std=gnu99" appendInstruction.c -lumfpack -outdir  */

#ifdef __linux__
#include <dlfcn.h>
#include <unistd.h>
#endif
#ifdef __APPLE__
#include <dlfcn.h>
#include <unistd.h>
#endif
#ifdef _WIN32
#include <windows.h>
#include <stdint.h>
#endif
#include <fcntl.h>
#include <mex.h>

#ifdef __linux__
void *libHandle=NULL;
#endif
#ifdef __APPLE__
void *libHandle=NULL;
#endif
#ifdef _WIN32
HMODULE libHandle=NULL;
#endif
void (*PappendInstruction4MEX)(
   /* inputs */
   int32_t *type,
   double *parameters,
   int64_t *operands,
   /* outputs */
   int64_t *index,
   /* sizes */
   mwSize mp,
   mwSize np,
   mwSize mo,
   mwSize no)=NULL;

#include "instructionsTableFunctions.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
   /* inputs */
   int32_t *type;
   double *parameters;
   int64_t *operands;
   /* outputs */
   int64_t *index;
   /* sizes */
   mwSize mp;
   mwSize np;
   mwSize mo;
   mwSize no;

   /* Process inputs */

   /* Check # inputs */
   if(nrhs!=3) {
      mexErrMsgIdAndTxt("appendInstruction:nrhs", "3 inputs required, %d found.",nrhs);
      return; }

   /* input type */
   if (mxGetNumberOfDimensions(prhs[0])!=2)
       mexErrMsgIdAndTxt("appendInstruction:prhs","input 1 (type) should have 2 dimensions, %d found.",mxGetNumberOfDimensions(prhs[0]));
   { const mwSize *dims=mxGetDimensions(prhs[0]);
   if (dims[0]!=1)
       mexErrMsgIdAndTxt("appendInstruction:prhs","input 1 (type) should have %d (=1) in dimension 1, %d found.",1,dims[0]);
   if (dims[1]!=1)
       mexErrMsgIdAndTxt("appendInstruction:prhs","input 1 (type) should have %d (=1) in dimension 2, %d found.",1,dims[1]);
   }
   if (!mxIsInt32(prhs[0]))
       mexErrMsgIdAndTxt("appendInstruction:prhs","input 1 (type) should have type int32");
   type=mxGetData(prhs[0]);
   /* input parameters */
   if (mxGetNumberOfDimensions(prhs[1])!=2)
       mexErrMsgIdAndTxt("appendInstruction:prhs","input 2 (parameters) should have 2 dimensions, %d found.",mxGetNumberOfDimensions(prhs[1]));
   { const mwSize *dims=mxGetDimensions(prhs[1]);
   mp=dims[0];
   np=dims[1];
   }
   if (!mxIsDouble(prhs[1]))
       mexErrMsgIdAndTxt("appendInstruction:prhs","input 2 (parameters) should have type double");
   parameters=mxGetData(prhs[1]);
   /* input operands */
   if (mxGetNumberOfDimensions(prhs[2])!=2)
       mexErrMsgIdAndTxt("appendInstruction:prhs","input 3 (operands) should have 2 dimensions, %d found.",mxGetNumberOfDimensions(prhs[2]));
   { const mwSize *dims=mxGetDimensions(prhs[2]);
   mo=dims[0];
   no=dims[1];
   }
   if (!mxIsInt64(prhs[2]))
       mexErrMsgIdAndTxt("appendInstruction:prhs","input 3 (operands) should have type int64");
   operands=mxGetData(prhs[2]);

   /* Process outputs */

   /* Check # outputs */
   if(nlhs!=1) {
      mexErrMsgIdAndTxt("appendInstruction:nrhs", "1 outputs required, %d found.",nlhs);
      return; }

   /* output index */
   { mwSize dims[]={1,1};
     plhs[0] = mxCreateNumericArray(2,dims,mxINT64_CLASS,mxREAL);
     index=mxGetData(plhs[0]); }

   /* Call function */
   if (!PappendInstruction4MEX) {
#ifdef __linux__
     libHandle = dlopen("/Users/hespanha/GitHub/tenscalc/lib/csparse/instructionsTable.so", RTLD_NOW);
     if (!libHandle) { printf("[%s] Unable to open library: %s\n",__FILE__, dlerror());return; }
     PappendInstruction4MEX = dlsym(libHandle, "appendInstruction4MEX");
     if (!PappendInstruction4MEX) { printf("[%s] Unable to get symbol: %s\n",__FILE__, dlerror());return; }
#elif __APPLE__
     libHandle = dlopen("/Users/hespanha/GitHub/tenscalc/lib/csparse/instructionsTable.dylib", RTLD_NOW);
     if (!libHandle) { printf("[%s] Unable to open library: %s\n",__FILE__, dlerror());return; }
     PappendInstruction4MEX = dlsym(libHandle, "appendInstruction4MEX");
     if (!PappendInstruction4MEX) { printf("[%s] Unable to get symbol: %s\n",__FILE__, dlerror());return; }
#elif _WIN32
     libHandle = LoadLibrary("/Users/hespanha/GitHub/tenscalc/lib/csparse/instructionsTable.dll");
     if (!libHandle) { printf("[%s] Unable to open library\n",__FILE__);return; }
     PappendInstruction4MEX = GetProcAddress(libHandle, "appendInstruction4MEX");
     if (!PappendInstruction4MEX) { printf("[%s] Unable to get symbol\n",__FILE__);return; }
#endif // _WIN32
   }
   PappendInstruction4MEX(type,parameters,operands,index,mp,np,mo,no);
}
