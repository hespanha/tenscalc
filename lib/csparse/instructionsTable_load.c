/* Created by script createGateway.m on 09-Jul-2017 16:49:58 */

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
/* mex -largeArrayDims -I/Users/hespanha/matlab_projects/tenscalc/TimDavis/SuiteSparse/include -L/Users/hespanha/matlab_projects/tenscalc/TimDavis/SuiteSparse/lib COPTIMFLAGS="-Ofast -msse -msse2 -msse3 -msse4 -msse4.1 -DNDEBUG" CFLAGS="\$CFLAGS -Wall -Werror -Wno-unused-variable -Wno-unused-result -std=gnu99" instructionsTable_load.c -lumfpack -outdir  */

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
void (*PinitInstructionsTable)();
void (*PinstructionsTableHeight4MEX)();
void (*PappendInstruction4MEX)();
void (*PappendUniqueInstruction4MEX)();
void (*PgetInstruction4MEX)();
void (*PfindInstructionsByType4MEX)();
void (*PgetDependencies4MEX)();
void (*PwriteCinstructionsC)();
void (*PwriteAsmInstructionsC)();

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
   /* inputs */
   double *load;

   /* Process inputs */

   /* Check # inputs */
   if(nrhs!=1) {
      mexErrMsgIdAndTxt("instructionsTable_load:nrhs", "1 inputs required, %d found.",nrhs);
      return; }

   /* input load */
   if (mxGetNumberOfDimensions(prhs[0])!=2)
       mexErrMsgIdAndTxt("instructionsTable_load:prhs","input 1 (load) should have 2 dimensions, %d found.",mxGetNumberOfDimensions(prhs[0]));
   { const mwSize *dims=mxGetDimensions(prhs[0]);
   if (dims[0]!=1)
       mexErrMsgIdAndTxt("instructionsTable_load:prhs","input 1 (load) should have %d (=1) in dimension 1, %d found.",1,dims[0]);
   if (dims[1]!=1)
       mexErrMsgIdAndTxt("instructionsTable_load:prhs","input 1 (load) should have %d (=1) in dimension 2, %d found.",1,dims[1]);
   }
   if (!mxIsDouble(prhs[0]))
       mexErrMsgIdAndTxt("instructionsTable_load:prhs","input 1 (load) should have type double");
   load=mxGetData(prhs[0]);

   /* Process outputs */

   /* Check # outputs */
   if(nlhs!=0) {
      mexErrMsgIdAndTxt("instructionsTable_load:nrhs", "0 outputs required, %d found.",nlhs);
      return; }


   /* Call function */
  if (!libHandle || load[0]) {
#ifdef __linux__
     libHandle = dlopen("/Users/hespanha/GitHub/tenscalc/lib/csparse/instructionsTable.so", RTLD_NOW);
     if (!libHandle) { printf("[%s] Unable to open library: %s\n",__FILE__, dlerror());return; }
     PinitInstructionsTable = dlsym(libHandle, "initInstructionsTable");
     if (!PinitInstructionsTable) { printf("[%s] Unable to get symbol: %s\n",__FILE__, dlerror());return; }
     PinstructionsTableHeight4MEX = dlsym(libHandle, "instructionsTableHeight4MEX");
     if (!PinstructionsTableHeight4MEX) { printf("[%s] Unable to get symbol: %s\n",__FILE__, dlerror());return; }
     PappendInstruction4MEX = dlsym(libHandle, "appendInstruction4MEX");
     if (!PappendInstruction4MEX) { printf("[%s] Unable to get symbol: %s\n",__FILE__, dlerror());return; }
     PappendUniqueInstruction4MEX = dlsym(libHandle, "appendUniqueInstruction4MEX");
     if (!PappendUniqueInstruction4MEX) { printf("[%s] Unable to get symbol: %s\n",__FILE__, dlerror());return; }
     PgetInstruction4MEX = dlsym(libHandle, "getInstruction4MEX");
     if (!PgetInstruction4MEX) { printf("[%s] Unable to get symbol: %s\n",__FILE__, dlerror());return; }
     PfindInstructionsByType4MEX = dlsym(libHandle, "findInstructionsByType4MEX");
     if (!PfindInstructionsByType4MEX) { printf("[%s] Unable to get symbol: %s\n",__FILE__, dlerror());return; }
     PgetDependencies4MEX = dlsym(libHandle, "getDependencies4MEX");
     if (!PgetDependencies4MEX) { printf("[%s] Unable to get symbol: %s\n",__FILE__, dlerror());return; }
     PwriteCinstructionsC = dlsym(libHandle, "writeCinstructionsC");
     if (!PwriteCinstructionsC) { printf("[%s] Unable to get symbol: %s\n",__FILE__, dlerror());return; }
     PwriteAsmInstructionsC = dlsym(libHandle, "writeAsmInstructionsC");
     if (!PwriteAsmInstructionsC) { printf("[%s] Unable to get symbol: %s\n",__FILE__, dlerror());return; }
#elif __APPLE__
     libHandle = dlopen("/Users/hespanha/GitHub/tenscalc/lib/csparse/instructionsTable.dylib", RTLD_NOW);
     if (!libHandle) { printf("[%s] Unable to open library: %s\n",__FILE__, dlerror());return; }
     PinitInstructionsTable = dlsym(libHandle, "initInstructionsTable");
     if (!PinitInstructionsTable) { printf("[%s] Unable to get symbol: %s\n",__FILE__, dlerror());return; }
     PinstructionsTableHeight4MEX = dlsym(libHandle, "instructionsTableHeight4MEX");
     if (!PinstructionsTableHeight4MEX) { printf("[%s] Unable to get symbol: %s\n",__FILE__, dlerror());return; }
     PappendInstruction4MEX = dlsym(libHandle, "appendInstruction4MEX");
     if (!PappendInstruction4MEX) { printf("[%s] Unable to get symbol: %s\n",__FILE__, dlerror());return; }
     PappendUniqueInstruction4MEX = dlsym(libHandle, "appendUniqueInstruction4MEX");
     if (!PappendUniqueInstruction4MEX) { printf("[%s] Unable to get symbol: %s\n",__FILE__, dlerror());return; }
     PgetInstruction4MEX = dlsym(libHandle, "getInstruction4MEX");
     if (!PgetInstruction4MEX) { printf("[%s] Unable to get symbol: %s\n",__FILE__, dlerror());return; }
     PfindInstructionsByType4MEX = dlsym(libHandle, "findInstructionsByType4MEX");
     if (!PfindInstructionsByType4MEX) { printf("[%s] Unable to get symbol: %s\n",__FILE__, dlerror());return; }
     PgetDependencies4MEX = dlsym(libHandle, "getDependencies4MEX");
     if (!PgetDependencies4MEX) { printf("[%s] Unable to get symbol: %s\n",__FILE__, dlerror());return; }
     PwriteCinstructionsC = dlsym(libHandle, "writeCinstructionsC");
     if (!PwriteCinstructionsC) { printf("[%s] Unable to get symbol: %s\n",__FILE__, dlerror());return; }
     PwriteAsmInstructionsC = dlsym(libHandle, "writeAsmInstructionsC");
     if (!PwriteAsmInstructionsC) { printf("[%s] Unable to get symbol: %s\n",__FILE__, dlerror());return; }
#elif _WIN32
     libHandle = LoadLibrary("/Users/hespanha/GitHub/tenscalc/lib/csparse/instructionsTable.dll");
     if (!libHandle) { printf("[%s] Unable to open library\n",__FILE__);return; }
     PinitInstructionsTable = GetProcAddress(libHandle, "initInstructionsTable");
     if (!PinitInstructionsTable) { printf("[%s] Unable to get symbol\n",__FILE__);return; }
     PinstructionsTableHeight4MEX = GetProcAddress(libHandle, "instructionsTableHeight4MEX");
     if (!PinstructionsTableHeight4MEX) { printf("[%s] Unable to get symbol\n",__FILE__);return; }
     PappendInstruction4MEX = GetProcAddress(libHandle, "appendInstruction4MEX");
     if (!PappendInstruction4MEX) { printf("[%s] Unable to get symbol\n",__FILE__);return; }
     PappendUniqueInstruction4MEX = GetProcAddress(libHandle, "appendUniqueInstruction4MEX");
     if (!PappendUniqueInstruction4MEX) { printf("[%s] Unable to get symbol\n",__FILE__);return; }
     PgetInstruction4MEX = GetProcAddress(libHandle, "getInstruction4MEX");
     if (!PgetInstruction4MEX) { printf("[%s] Unable to get symbol\n",__FILE__);return; }
     PfindInstructionsByType4MEX = GetProcAddress(libHandle, "findInstructionsByType4MEX");
     if (!PfindInstructionsByType4MEX) { printf("[%s] Unable to get symbol\n",__FILE__);return; }
     PgetDependencies4MEX = GetProcAddress(libHandle, "getDependencies4MEX");
     if (!PgetDependencies4MEX) { printf("[%s] Unable to get symbol\n",__FILE__);return; }
     PwriteCinstructionsC = GetProcAddress(libHandle, "writeCinstructionsC");
     if (!PwriteCinstructionsC) { printf("[%s] Unable to get symbol\n",__FILE__);return; }
     PwriteAsmInstructionsC = GetProcAddress(libHandle, "writeAsmInstructionsC");
     if (!PwriteAsmInstructionsC) { printf("[%s] Unable to get symbol\n",__FILE__);return; }
#endif // _WIN32
  }
#ifdef __linux__
  if (load[0]==0) { while (!dlclose(libHandle)) printf(".");}
#elif __APPLE__
  if (load[0]==0) { while (!dlclose(libHandle)) printf("."); }
#elif _WIN32
  if (load[0]==0) { while (FreeLibrary(libHandle)) printf("."); }
#endif // _WIN32
}
