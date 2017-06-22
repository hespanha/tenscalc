/* Created by script createGateway.m on 21-Jun-2017 23:54:06 */

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
/* mex -largeArrayDims -I/Users/hespanha/matlab_projects/tenscalc/TimDavis/SuiteSparse/include -L/Users/hespanha/matlab_projects/tenscalc/TimDavis/SuiteSparse/lib COPTIMFLAGS="-Ofast -msse -msse2 -msse3 -msse4 -msse4.1 -DNDEBUG" CFLAGS="\$CFLAGS -Wall -Werror -Wno-unused-variable -Wno-unused-result -std=gnu99" getInstruction.c -lumfpack -outdir  */

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
#include "mex.h"

#ifdef __linux__
void *libHandle=NULL;
#endif
#ifdef __APPLE__
void *libHandle=NULL;
#endif
#ifdef _WIN32
HMODULE libHandle=NULL;
#endif
void (*PgetInstruction4MEX)(
   /* inputs */
   int64_t *index,
   /* outputs */
   int32_t *type,
   mxArray **parameters,
   mxArray **operands)=NULL;

#include "instructionsTableFunctions.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
   /* inputs */
   int64_t *index;
   /* outputs */
   int32_t *type;
   mxArray **parameters;
   mxArray **operands;

   /* Process inputs */

   /* Check # inputs */
   if(nrhs!=1) {
      mexErrMsgIdAndTxt("getInstruction:nrhs", "1 inputs required, %d found.",nrhs);
      return; }

   /* input index */
   if (mxGetNumberOfDimensions(prhs[0])!=2)
       mexErrMsgIdAndTxt("getInstruction:prhs","input 1 (index) should have 2 dimensions, %d found.",mxGetNumberOfDimensions(prhs[0]));
   { const mwSize *dims=mxGetDimensions(prhs[0]);
   if (dims[0]!=1)
       mexErrMsgIdAndTxt("getInstruction:prhs","input 1 (index) should have %d (=1) in dimension 1, %d found.",1,dims[0]);
   if (dims[1]!=1)
       mexErrMsgIdAndTxt("getInstruction:prhs","input 1 (index) should have %d (=1) in dimension 2, %d found.",1,dims[1]);
   }
   if (!mxIsInt64(prhs[0]))
       mexErrMsgIdAndTxt("getInstruction:prhs","input 1 (index) should have type int64");
   index=mxGetData(prhs[0]);

   /* Process outputs */

   /* Check # outputs */
   if(nlhs!=3) {
      mexErrMsgIdAndTxt("getInstruction:nrhs", "3 outputs required, %d found.",nlhs);
      return; }

   /* output type */
   { mwSize dims[]={1,1};
     plhs[0] = mxCreateNumericArray(2,dims,mxINT32_CLASS,mxREAL);
     type=mxGetData(plhs[0]); }
   /* output parameters */
   parameters=plhs+1;
   /* output operands */
   operands=plhs+2;

   /* Call function */
   if (!PgetInstruction4MEX) {
#ifdef __linux__
     libHandle = dlopen("instructionsTable.so", RTLD_NOW);
     if (!libHandle) { printf("[%s] Unable to open library: %s\n",__FILE__, dlerror());return; }
     PgetInstruction4MEX = dlsym(libHandle, "getInstruction4MEX");
     if (!PgetInstruction4MEX) { printf("[%s] Unable to get symbol: %s\n",__FILE__, dlerror());return; }
#elif __APPLE__
     libHandle = dlopen("./instructionsTable.dylib", RTLD_NOW);
     if (!libHandle) { printf("[%s] Unable to open library: %s\n",__FILE__, dlerror());return; }
     PgetInstruction4MEX = dlsym(libHandle, "getInstruction4MEX");
     if (!PgetInstruction4MEX) { printf("[%s] Unable to get symbol: %s\n",__FILE__, dlerror());return; }
#elif _WIN32
     libHandle = LoadLibrary("./instructionsTable.dll");
     if (!libHandle) { printf("[%s] Unable to open library\n",__FILE__);return; }
     PgetInstruction4MEX = GetProcAddress(libHandle, "getInstruction4MEX");
     if (!PgetInstruction4MEX) { printf("[%s] Unable to get symbol\n",__FILE__);return; }
#endif // _WIN32
   }
   PgetInstruction4MEX(index,type,parameters,operands);
}
