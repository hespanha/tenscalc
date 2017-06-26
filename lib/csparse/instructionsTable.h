/* Created by script createGateway.m on 25-Jun-2017 19:01:17 */

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
/* gcc -I/Applications/MATLAB_R2016b.app/include -I/Applications/MATLAB_R2016b.app/extern/include -I/Users/hespanha/matlab_projects/tenscalc/TimDavis/SuiteSparse/include -L/Applications/MATLAB_R2016b.app/bin/maci64 -L/Users/hespanha/matlab_projects/tenscalc/TimDavis/SuiteSparse/lib -Wall -Werror -Wno-unused-variable -Wno-unused-result -std=gnu99  -DDYNAMIC_LIBRARY -DNDEBUG -dynamiclib -fvisibility=hidden -current_version 1.0 -compatibility_version 1.0 -Ofast -msse -msse2 -msse3 -msse4 -msse4.1 instructionsTable.c -lmx -lmex -lmat -lstdc++ -lumfpack -o /Users/hespanha/GitHub/tenscalc/lib/csparse/instructionsTable.dylib */

#include <mex.h>
#include <stdint.h>
#ifdef DYNAMIC_LIBRARY
#ifdef __APPLE__
#define EXPORT __attribute__((visibility("default")))
#elif __linux__
#define EXPORT __attribute__((visibility("default")))
#elif _WIN32
#define EXPORT __declspec(dllexport)
#endif
#else
#define EXPORT 
#endif
EXPORT void initInstructionsTable();
EXPORT void instructionsTableHeight4MEX(   /* outputs */
   int64_t *height);
EXPORT void appendInstruction4MEX(   /* inputs */
   int32_t *type,
   double *parameters,
   int64_t *operands,
   /* outputs */
   int64_t *index,
   /* sizes */
   mwSize mp,
   mwSize np,
   mwSize mo,
   mwSize no);
EXPORT void appendUniqueInstruction4MEX(   /* inputs */
   int32_t *type,
   double *parameters,
   int64_t *operands,
   /* outputs */
   int64_t *index,
   /* sizes */
   mwSize mp,
   mwSize np,
   mwSize mo,
   mwSize no);
EXPORT void getInstruction4MEX(   /* inputs */
   int64_t *index,
   /* outputs */
   int32_t *type,
   mxArray **parameters,
   mxArray **operands);
EXPORT void findInstructionsByType4MEX(   /* inputs */
   int32_t *type,
   /* outputs */
   mxArray **ifOfType);
EXPORT void getDependencies4MEX(   /* outputs */
   mxArray **children,
   mxArray **parents);
EXPORT void writeCinstructionsC(   /* inputs */
   int64_t *indices,
   int64_t *memoryLocations,
   /* outputs */
   int64_t *countFlops,
   /* sizes */
   mwSize nInstructions,
   mwSize NInstructions);
EXPORT void writeAsmInstructionsC(   /* inputs */
   int64_t *indices,
   int64_t *memoryLocations,
   /* sizes */
   mwSize nInstructions,
   mwSize NInstructions);
