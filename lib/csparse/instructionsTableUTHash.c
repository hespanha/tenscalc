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

#ifdef createGateway

MEXfunction  initInstructionsTable

Cfunction initInstructionsTable
include instructionsTableFunctions.h

MEXfunction  instructionsTableHeight

Cfunction instructionsTableHeight4MEX
include instructionsTableFunctions.h

outputs
      int64 height [1]

MEXfunction appendInstruction

Cfunction appendInstruction4MEX
include instructionsTableFunctions.h

inputs
      int32  type [1]
      double parameters [mp,np]
      int64 operands   [mo,no]
outputs
      int64 index [1]

MEXfunction appendUniqueInstruction

Cfunction appendUniqueInstruction4MEX
include instructionsTableFunctions.h

inputs
      int32  type [1]
      double parameters [mp,np]
      int64 operands   [mo,no]
outputs
      int64 index [1]


MEXfunction getInstruction

Cfunction getInstruction4MEX
include instructionsTableFunctions.h

inputs
      int64 index [1]

outputs
      int32 type [1]
      double parameters [~]
      int64 operands [~]

MEXfunction findInstructionsByType

Cfunction findInstructionsByType4MEX
include instructionsTableFunctions.h

inputs
      int32 type [1]

outputs
      int8 ifOfType [~]

MEXfunction getDependencies

Cfunction getDependencies4MEX
include instructionsTableFunctions.h

outputs
      int64 children [~]
      int64 parents [~]

MEXfunction  writeCinstructionsC
Cfunction writeCinstructionsC
include instructionsTableFunctions.h

inputs 
      int64 indices [nInstructions]
      int64 memoryLocations [NInstructions]

outputs
      int64 countFlops[16]

MEXfunction  writeAsmInstructionsC
Cfunction writeAsmInstructionsC
include instructionsTableFunctions.h

inputs 
      int64 indices [nInstructions]
      int64 memoryLocations [NInstructions]


#endif

#include <stdint.h>   /* needed by uint64_t */
#include <inttypes.h> /* needed by PRId64 */
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <time.h>
#include <math.h>

//  generate a fair miss rate may benefit from the built-in Bloom filter support 
#define HASH_BLOOM 28
  
#ifdef __linux__
#include <unistd.h>
#elif __APPLE__
#include <unistd.h>
#elif _WIN32
#include <stdint.h>
#include <windows.h>
#endif
        
// Using UTHash

/* The Bernstein hash function, used in Perl prior to v5.6. Note (x<<5+x)=x*33. */
#define HASH_BER0(key,keylen,hashv)                                              \
do {                                                                             \
  unsigned _hb_keylen = (unsigned)keylen;                                        \
  const unsigned char *_hb_key = (const unsigned char*)(key);                    \
  (hashv) = 0;                                                                   \
  while (_hb_keylen-- != 0U) {                                                   \
    (hashv) = (((hashv) << 5) + (hashv)) + *_hb_key++;                           \
  }                                                                              \
} while (0)
#define HASH_BER1(key,keylen,hashv)                                              \
do {                                                                             \
  unsigned _hb_keylen = (unsigned)keylen;                                        \
  const unsigned char *_hb_key = (const unsigned char*)(key);                    \
  while (_hb_keylen-- != 0U) {                                                   \
    (hashv) = (((hashv) << 5) + (hashv)) + *_hb_key++;                           \
  }                                                                              \
} while (0)

/* The Bernstein hash function, used in Perl prior to v5.6. Note (x<<5+x)=x*33. */
#define HASH_TENSCALC(key,keylen,hashv)                                          \
do {                                                                             \
  instruction_t *ptrI=(instruction_t*)key;                                       \
  int64_t  *_keyop=instructionsTable.operandsBuffer+ptrI->operands;              \
  unsigned _keylenop=sizeof(int64_t)*(ptrI->nOperands);                          \
  double   *_keypar=instructionsTable.parametersBuffer+ptrI->parameters;         \
  unsigned _keylenpar=sizeof(double)*(ptrI->nParameters);                        \
  HASH_BER0(key,keylen,hashv);						         \
  HASH_BER1(_keyop,_keylenop,hashv);                                             \
  HASH_BER1(_keypar,_keylenpar,hashv);                                           \
 } while (0)

#define HASH_FUNCTION HASH_TENSCALC
#include "uthash.h"



#include "instructionsTableFunctions.h"

// Structure used to store instructions

#define MAX_INSTRUCTIONS_PER_TABLE  5000000LL
#define MAX_PARAMETERS_PER_TABLE   20000000LL // lasso2 example needs very large
#define MAX_OPERANDS_PER_TABLE     10000000LL 

#define MAX_TERMS_PERLINE 100

typedef struct instruction_s {
  // KEY STARTS HERE
  instructionType_t type;           // type of instruction
  int64_t nOperands;                // # of operands
  int64_t nParameters;              // # of parameters
  // KEY ENDS HERE
  int64_t operands;                 // position of operands in buffer
  int64_t parameters;               // position of parameters in buffer
  UT_hash_handle hh;                // makes this structure hashable
  int64_t index;
} instruction_t;

typedef struct instructionsTable_s {
  // instructions array
  instruction_t instructionsBuffer[MAX_INSTRUCTIONS_PER_TABLE]; // buffer to store instructions
  int64_t       nInstructions;  // # of instructions currently in the instructions buffer
  
  // parameters buffer 
  double        parametersBuffer[MAX_PARAMETERS_PER_TABLE];    // buffer to store parameters
  int64_t       nParameters;     // # of parameters currently in the parameters buffer

  // operands buffer 
  int64_t       operandsBuffer[MAX_OPERANDS_PER_TABLE];    // buffer to store operands
  int64_t       nOperands;     // # of operands currently in the operands buffer

  // sorted instructions
  boolean_T     isSorted;
  int64_t       sortedIndices[MAX_INSTRUCTIONS_PER_TABLE];

} instructionsTable_t;

#define verboseLevel 0

#ifdef IGNORE_MEX
#undef printf
#define printf0(...) printf(__VA_ARGS__)
#else
#define printf0(...) mexPrintf(__VA_ARGS__)
#endif 

#if verboseLevel>=1
#define printf1(...) printf(__VA_ARGS__)
#else
#define printf1(...) 
#endif

#if verboseLevel>=2
#define printf2(...) printf(__VA_ARGS__)
#else
#define printf2(...) 
#endif

#if verboseLevel>=3
#define printf3(...) printf(__VA_ARGS__)
#else
#define printf3(...) 
#endif

#ifdef DYNAMIC_LIBRARY
#ifdef __linux__
/* Initializer */
__attribute__((constructor))
static void initializer(void) {
  printf0("%s: loading dynamic library (UTHash)\n", __FILE__);
}
/* Finalizer */
__attribute__((destructor))
static void finalizer(void) {
  printf0("%s: unloading dynamic library (UTHash)\n", __FILE__);
  initInstructionsTable(); // free memory
}
#elif __APPLE__
/* Initializer */
__attribute__((constructor))
static void initializer(void) {
  printf0("%s: loading dynamic library (UTHash)\n", __FILE__);
}
/* Finalizer */
__attribute__((destructor))
static void finalizer(void) {
  printf0("%s: unloading dynamic library (UTHash)\n", __FILE__);
  initInstructionsTable(); // free memory
}
#elif _WIN32
#include <windows.h>
BOOL WINAPI DllMain(HINSTANCE hinstDLL, DWORD fdwReason, LPVOID lpvReserved)
{
    if (fdwReason == DLL_PROCESS_ATTACH) {
      printf0("%s: loading dynamic library\n", __FILE__);
      return TRUE; }
    else if (fdwReason == DLL_PROCESS_DETACH) {
      printf0("%s: unloading dynamic library\n", __FILE__);
      initInstructionsTable(); // free memory
      return TRUE; }
}
#endif
#endif


/*******************************************************************************
* Initializes the instructions table.
*******************************************************************************/
EXPORT instructionsTable_t instructionsTable;
EXPORT instruction_t *instructionsTable_hash=NULL; // initialize hash table

EXPORT void initInstructionsTable()
{
  // Currently the instruction table is a global variable, 
  // but it could be dynamically allocated using malloc.

  if (instructionsTable_hash==NULL) {
    printf0("initInstructionsTable: Initializing table... ");
  } else {
    printf0("initInstructionsTable: Initializing table (removing %"PRId64" instructions, %"PRId64" parameters, %"PRId64" operands)... ",
	    instructionsTable.nInstructions,instructionsTable.nParameters,instructionsTable.nOperands);
    HASH_CLEAR(hh,instructionsTable_hash);
  }
  instructionsTable.nInstructions=0;
  instructionsTable.nParameters=0;
  instructionsTable.nOperands=0;
  instructionsTable.isSorted=false;
  printf0("done\n");
}

/*******************************************************************************
* Get number of instructions
*******************************************************************************/
EXPORT int64_t instructionsTableHeight()
{
  // Currently the instruction table is a global variable, 
  // but it could be dynamically allocated using malloc.

  return instructionsTable.nInstructions;
}

#ifndef IGNORE_MEX
EXPORT void instructionsTableHeight4MEX(/* outputs */
				  int64_t *height)
{
  
  (*height)=instructionsTableHeight();
}
#endif

/*******************************************************************************
* Appends instruction at the end of the table, and return index to it.
*
* In the process, computes the hash value that is used to make comparisons
* faster (by first comparing the hash value and only then making a
* full comparison
*
* Returns -1 if any of buffers is full. 
*******************************************************************************/
#if verboseLevel<3
inline
#endif
void instruction2buffer(instruction_t *ptr,
			int64_t index,
			instructionType_t type,
			int64_t nParameters,
			double *parameters,
			int64_t nOperands,
			int64_t *operands)
{
  int64_t j;

  // clear structure
  memset(ptr, 0, sizeof(instruction_t));
  
  // copy instruction index
  ptr->index=index;
  
  // copy instruction type table
  ptr->type=type;
  
  // Save parameters to parametersBuffer
  if (nParameters>0){
    double *destination=instructionsTable.parametersBuffer+instructionsTable.nParameters,
      *source=parameters;
    for (j=nParameters-1;j>=0;j--) {
      *(destination++)=*(source++); }; 
    // save parameters, nParameters
    ptr->parameters=instructionsTable.nParameters;
    ptr->nParameters=nParameters;
    instructionsTable.nParameters += nParameters;
  } else {
    ptr->parameters=0;
    ptr->nParameters=0;
  }
  
  // Save operands to operandsBuffer
  if (nOperands>0) {
    int64_t *destination=instructionsTable.operandsBuffer+instructionsTable.nOperands,
      *source=operands;
    for (j=nOperands-1;j>=0;j--) {
      *(destination++)=*(source++); }; 
    // save operands, nOperands
    ptr->operands=instructionsTable.nOperands;
    ptr->nOperands=nOperands;
    instructionsTable.nOperands += nOperands; 
  } else {
    ptr->operands=0;
    ptr->nOperands=0;
  }

#if verboseLevel>=3
  printf3("  instruction index=%"PRId64": type=%d, nParameters=%"PRId64", nOperands=%"PRId64"",
	  ptr->index,ptr->type,ptr->nParameters,ptr->nOperands);
  if (ptr->nParameters>0) {
    int l;
    printf3(", parameters = [");
    for (l=0;l<ptr->nParameters;l++) {
      printf3("%g,",instructionsTable.parametersBuffer[ptr->parameters+l]); }
    printf3("]");
  }
  if (ptr->nOperands>0) {
    int l;
    printf3(", operands = [");
    for (l=0;l<ptr->nOperands;l++) {
      printf3("%"PRId64",",instructionsTable.operandsBuffer[ptr->operands+l]); }
    printf3("]");
  }
  printf3(": ");
#endif
}


EXPORT int64_t appendInstruction(instructionType_t type,
				 int64_t nParameters,
				 double *parameters,
				 int64_t nOperands,
				 int64_t *operands)
{
  if (instructionsTable.nInstructions+1     > MAX_INSTRUCTIONS_PER_TABLE ||
      instructionsTable.nParameters+nParameters > MAX_PARAMETERS_PER_TABLE ||
      instructionsTable.nOperands+nOperands > MAX_OPERANDS_PER_TABLE) {
    printf0("appendInstruction: table is full (I:%"PRId64">%"PRId64" | P:%"PRId64">%"PRId64" | O:%"PRId64">%"PRId64")\n",
	    instructionsTable.nInstructions+1,MAX_INSTRUCTIONS_PER_TABLE,
	    instructionsTable.nParameters+nParameters,MAX_PARAMETERS_PER_TABLE,
	    instructionsTable.nOperands+nOperands,MAX_OPERANDS_PER_TABLE);
    return -1;
  }
  
  int64_t index=instructionsTable.nInstructions++;
  instruction_t *ptr=instructionsTable.instructionsBuffer+index;

  instruction2buffer(ptr,index,type,nParameters,parameters,nOperands,operands);
  
  HASH_ADD(hh,instructionsTable_hash,type,
	   offsetof(instruction_t, nParameters)  // offset of last key field
	   + sizeof(int64_t)                    // size of last key field
	   - offsetof(instruction_t, type),     // offset of first key field
	   ptr);
  
  printf3("appendInstruction: appended instruction %"PRId64", type %d, nParameters %"PRId64", nOperands %"PRId64" (nInstructions=%"PRId64",nParameters=%"PRId64",nOperands=%"PRId64")\n",
	  index,type,
	  nParameters,
	  nOperands,
	  instructionsTable.nInstructions,instructionsTable.nParameters,instructionsTable.nOperands);

  instructionsTable.sortedIndices[index]=index;
  instructionsTable.isSorted=false;

  return index+1;  // 1-based indexing  
}


#ifndef IGNORE_MEX
EXPORT void appendInstruction4MEX(  /* inputs */
				  int32_t *type,
				  double *parameters,
				  int64_t *operands,
				  /* outputs */
				  int64_t *index,
				  /* sizes */
				  mwSize mp,
				  mwSize np,
				  mwSize mo,
				  mwSize no)
{
  printf3("appendInstruction4MEX: type=%d, Psize=%llux%llu, Osize=%llux%llu\n",
	  *type,mp,np,mo,no);
  (*index)=appendInstruction(*type,mp*np,parameters,mo*no,operands);
}
#endif

/*******************************************************************************
 * Compares two instructions.
 * Returns 0 if equal, -1 if the first is "smaller" than the second, and +1 if "larger"
*******************************************************************************/

EXPORT int compareInstructions(int64_t *index1,
			       int64_t *index2)
{
  instruction_t *ptr1=instructionsTable.instructionsBuffer+(*index1);
  instruction_t *ptr2=instructionsTable.instructionsBuffer+(*index2);
  int64_t j;

  //printf3("Comparing instructions %"PRId64" with %"PRId64" ",*index1,*index2);
  if (ptr1->type<ptr2->type)
    return -1;
  else if (ptr1->type>ptr2->type)
    return 1;
  
  if (ptr1->nOperands<ptr2->nOperands)
    return -1;
  else if (ptr1->nOperands>ptr2->nOperands)
    return 1;
  
  if (ptr1->nParameters<ptr2->nParameters)
    return -1;
  else if (ptr1->nParameters>ptr2->nParameters)
    return 1;
  
  // printf3("same header... ");

  { double
      *pt1=instructionsTable.parametersBuffer+ptr1->parameters,
      *pt2=instructionsTable.parametersBuffer+ptr2->parameters,
      cmp;
    for (j=ptr1->nParameters; j>0 ;j--) {
      cmp=(*(pt1++))-(*(pt2++));
      if (cmp<0)
	return -1;
      else if (cmp>0)
	return 1;
    }
  }
    
  //printf3("same parameters... ");

  { int64_t
      *pt1=instructionsTable.operandsBuffer+ptr1->operands,
      *pt2=instructionsTable.operandsBuffer+ptr2->operands,
      cmp;
    for (j=ptr1->nOperands; j>0 ;j--) {
      cmp=(*(pt1++))-(*(pt2++));
      if (cmp<0)
	return -1;
      else if (cmp>0)
	return 1;
    }
  }
  
  printf3("equal... ");

  return 0;
}


/*******************************************************************************
 * Search for instruction, if already in the table, return index, 
 * otherwise append at the end of the table and return the index of the new instruction. 
 * Returns -1 is either buffer is full.
*******************************************************************************/
EXPORT int64_t appendUniqueInstruction(instructionType_t type,
				       int64_t nParameters,
				       double *parameters,
				       int64_t nOperands,
				       int64_t *operands)
{
  // save values, prior ro adding instruction
  int64_t 
    oldP=instructionsTable.nParameters,
    oldO=instructionsTable.nOperands;

  if (instructionsTable.nInstructions+1     > MAX_INSTRUCTIONS_PER_TABLE ||
      instructionsTable.nParameters+nParameters > MAX_PARAMETERS_PER_TABLE ||
      instructionsTable.nOperands+nOperands > MAX_OPERANDS_PER_TABLE) {
    printf0("appendInstruction: table is full (I:%"PRId64">%"PRId64" | P:%"PRId64">%"PRId64" | O:%"PRId64">%"PRId64")\n",
	    instructionsTable.nInstructions+1,MAX_INSTRUCTIONS_PER_TABLE,
	    instructionsTable.nParameters+nParameters,MAX_PARAMETERS_PER_TABLE,
	    instructionsTable.nOperands+nOperands,MAX_OPERANDS_PER_TABLE);
    return -1;
  }
  
  int64_t index=instructionsTable.nInstructions++;
  instruction_t
    *ptrNew=instructionsTable.instructionsBuffer+index,
    *ptrSearch;

  printf3("appendUniqueInstruction: appending instruction type %d, nParameters %"PRId64", nOperands %"PRId64" (nInstructions=%"PRId64",nParameters=%"PRId64",nOperands=%"PRId64")\n",
	  type,
	  nParameters,
	  nOperands,
	  instructionsTable.nInstructions,instructionsTable.nParameters,instructionsTable.nOperands);
  
  instruction2buffer(ptrNew,index,type,nParameters,parameters,nOperands,operands);
  
  HASH_FIND(hh,instructionsTable_hash,
	    &ptrNew->type,
	    offsetof(instruction_t, nParameters)  // offset of last key field
	    + sizeof(int64_t)                    // size of last key field
	    - offsetof(instruction_t, type),     // offset of first key field
	    ptrSearch);
  
  // Check if HASH_FIND was fooled by same key, but not same parameters/operators
  if (ptrSearch) {
    if (compareInstructions(&index,&ptrSearch->index)!=0) {
      printf3("****Should make sure HASH_FIND looks into parameter & operator values\n");
      ptrSearch=NULL;
    }
  }

  /*
  // DEBUG: override HASH_FIND and search over instructions for match
  if (index>0) {
    int64_t i;
    ptrSearch=NULL;
    for (i=index-1;i>=0;i--)
      if (compareInstructions(&i,&index)==0) {
	// instruction exists, remove the just added instruction, return index
	printf2("FOUND %"PRId64" %"PRId64"\n",i,index);
	ptrSearch=instructionsTable.instructionsBuffer+i;
	break;
    }
    printf("%"PRId64" ",(int64_t)ptrSearch);
  }
  */

  if (ptrSearch) {
    // instruction exists, remove the just added instruction, return index
    instructionsTable.nInstructions--;
    instructionsTable.nParameters=oldP;
    instructionsTable.nOperands=oldO;
    index=ptrSearch->index;
    printf3("reusing (%"PRId64")\n",index);
  } else {
    ptrNew=instructionsTable.instructionsBuffer+index;
    HASH_ADD(hh,instructionsTable_hash,
	     type,
	     offsetof(instruction_t, nParameters)  // offset of last key field
	     + sizeof(int64_t)                    // size of last key field
	     - offsetof(instruction_t, type),     // offset of first key field
	     ptrNew);
    instructionsTable.sortedIndices[index]=index;
    instructionsTable.isSorted=false;
    printf3("new (%"PRId64")\n",index);
  }
  return index+1;  // 1-based indexing  
}

#ifndef IGNORE_MEX
EXPORT void appendUniqueInstruction4MEX(/* inputs */
					int32_t *type,
					double *parameters,
					int64_t *operands,
					/* outputs */
					int64_t *index,
					/* sizes */
					mwSize mp,
					mwSize np,
					mwSize mo,
					mwSize no)
{
  printf3("appendUniqueInstruction4MEX: type=%d, Psize=%llux%llu, Osize=%llux%llu\n",
	  *type,mp,np,mo,no);
  (*index)=appendUniqueInstruction(*type,mp*np,parameters,mo*no,operands);
}
#endif

/*******************************************************************************
 * Sort instructions
 *******************************************************************************/

EXPORT void sortInstructions()
{
  qsort(instructionsTable.sortedIndices,instructionsTable.nInstructions,sizeof(int64_t),(int (*)(const void *, const void *))&compareInstructions);
  //heapsort(instructionsTable.sortedIndices,instructionsTable.nInstructions,sizeof(int64_t),(int (*)(const void *, const void *))&compareInstructions);
  //mergesort(instructionsTable.sortedIndices,instructionsTable.nInstructions,sizeof(int64_t),(int (*)(const void *, const void *))&compareInstructions);
  instructionsTable.isSorted=true;
}

/*******************************************************************************
 * Get instruction, given its index.
 * Returns -1 if the index is not valid
 *******************************************************************************/
EXPORT int getInstruction(int64_t index,
			  instructionType_t *type, // return type by reference
			  int64_t *nParameters,    // return nParameters by reference
			  double  **parameters,    // return pointer to parameters by reference
			  int64_t *nOperands,      // return nOperands by reference
			  int64_t **operands)      // return pointer to operands by reference
{
  index--; // back to 0-based indexing
  if (index<0 || index >=instructionsTable.nInstructions) {
    printf0("getInstruction: invalid index %"PRId64" (expected in [0,%"PRId64"]\n",index,instructionsTable.nInstructions-1);
    return -1; }
  instruction_t *ptr = instructionsTable.instructionsBuffer+index;
  *type=ptr->type;
  *nParameters=ptr->nParameters;
  *parameters=instructionsTable.parametersBuffer+ptr->parameters;
  *nOperands=ptr->nOperands;
  *operands=instructionsTable.operandsBuffer+ptr->operands;
  printf2("getInstruction(%"PRId64"): type=%d, nParameters=%"PRId64", nOperands=%"PRId64"\n",
	  index,*type,
	  *nParameters,
	  *nOperands);

  return 0;
}

#ifndef IGNORE_MEX
EXPORT void getInstruction4MEX(  /* inputs */
			       int64_t *index,
			       /* outputs */
			       int32_t *type,
			       mxArray **parameters,
			       mxArray **operands
)
{
  int64_t *ptrOperands;
  double *ptrParameters;
  int64_t nParameters,nOperands;
  
  int64_t rc=getInstruction(*index,(instructionType_t*)type,
			    &nParameters,&ptrParameters,&nOperands,&ptrOperands);
  printf3("getInstruction4MEX(%"PRId64"): type=%d, nParameters=%"PRId64", nOperands=%"PRId64"\n",
	  *index,*type,
	  nParameters,
	  nOperands);

  if (rc>=0) {
    { mwSize dims[]={1,nParameters};
      (*parameters)=mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
      double *ptr=mxGetData(*parameters); 
      memcpy(ptr,ptrParameters,sizeof(*ptrParameters)*nParameters); }

    { mwSize dims[]={1,nOperands};
      (*operands)=mxCreateNumericArray(2,dims,mxINT64_CLASS,mxREAL);
      int64_t *ptr=mxGetData(*operands); 
      memcpy(ptr,ptrOperands,sizeof(*ptrOperands)*nOperands); }
  } else {
    printf0("getInstruction: invalid index\n");
  }
}
#endif

/*******************************************************************************
 * find instructions of a given type.
 *******************************************************************************/
EXPORT void findInstructionsByType(instructionType_t type, // return type by reference
				   int8_t *isOfType)       // boolean array where true/false will be stored
{
  instruction_t *ptr = instructionsTable.instructionsBuffer;
  for (int64_t i=instructionsTable.nInstructions;i>0;i--)
    *(isOfType++)=((ptr++)->type==type)?1:0;
  return;
}

#ifndef IGNORE_MEX
EXPORT void findInstructionsByType4MEX(/* inputs */
				       int32_t *type,
				       /* outputs */
				       mxArray **isOfType
)
{
  mwSize dims[]={instructionsTable.nInstructions,1};
  (*isOfType)=mxCreateNumericArray(2,dims,mxINT8_CLASS,mxREAL);
  int8_t *ptr=mxGetData(*isOfType); 
  findInstructionsByType(*type,ptr);
}
#endif

/*******************************************************************************
 * get dependencies
 *******************************************************************************/
EXPORT void getDependencies(int64_t *children,       // array of children (dependents)
			    int64_t *parents)        // array of parents
{
  int64_t n=0;
  instruction_t *ptr = instructionsTable.instructionsBuffer;
  int64_t instr=1; // 1-based indexing
  for (int64_t i=instructionsTable.nInstructions;i>0;i--,instr++) {
    int64_t *operands=instructionsTable.operandsBuffer+ptr->operands;
    for (int64_t j=(ptr++)->nOperands;j>0;j--) {
      if (++n>instructionsTable.nOperands) {
	printf0("\ngetDependencies: error too many dependencies!\n");
      } else {
	*(children++)=instr; 
	*(parents++)=*(operands++); }
    }
  }
    
  if (n!=instructionsTable.nOperands) {
    printf0("\ngetDependencies: error wrong number of dependencies (%"PRId64", expected %"PRId64")\n",n,instructionsTable.nOperands);
  }
  
  return;
}

#ifndef IGNORE_MEX
EXPORT void getDependencies4MEX(/* outputs */
				mxArray **children,
				mxArray **parents)
{
  mwSize dims[]={instructionsTable.nOperands,1};
  (*children)=mxCreateNumericArray(2,dims,mxINT64_CLASS,mxREAL);
  int64_t *pChildren=mxGetData(*children); 
  (*parents)=mxCreateNumericArray(2,dims,mxINT64_CLASS,mxREAL);
  int64_t *pParents=mxGetData(*parents); 
  getDependencies(pChildren,pParents);
}
#endif

/*******************************************************************************
 * Write C code corresponding to an array of instructions
 * Returns -1 if any index is not valid, otherwise returns 0
 *******************************************************************************/
EXPORT int writeCinstructionsC(/* inputs */
			       int64_t *indices,          // indices of instructions to write
			       int64_t *memoryLocations,  // memory locations for
			                                  // all instructions
                               /* outputs */
                               int64_t *countFlops,          // array with instruction counts

                               /* sizes */
			       mwSize nInstructions,     // # of instructions to write
			       mwSize NInstructions)     // total # of instructions
{
  instructionType_t type;
  int64_t nParameters; 
  double  *parameters;
  int64_t nOperands;   
  int64_t *operands;
  int64_t odiv;
  int nSums;

#define countFlops_nsum     countFlops[P_nsum-1]
#define countFlops_nprod    countFlops[P_nprod-1]
#define countFlops_ndiv     countFlops[P_ndiv-1]
#define countFlops_nif      countFlops[P_nif-1]
#define countFlops_nclp     countFlops[P_nclp-1]
#define countFlops_nabs     countFlops[P_nabs-1]
#define countFlops_nround   countFlops[P_nround-1]
#define countFlops_nceil    countFlops[P_nceil-1]
#define countFlops_nfloor   countFlops[P_nfloor-1]
#define countFlops_nsign    countFlops[P_nsign-1]
#define countFlops_nsqrt    countFlops[P_nsqrt-1]
#define countFlops_npow     countFlops[P_npow-1]
#define countFlops_ntrig    countFlops[P_ntrig-1]
#define countFlops_nlog     countFlops[P_nlog-1]
#define countFlops_nexp     countFlops[P_nexp-1]
#define countFlops_numfpack countFlops[P_numfpack-1]

#if P_nCountFlops != 16
#error "update outputs field of writeCinstructionsC() (line 102)"
#endif
 
  FILE *fid=fopen("tmp_toremove.c","w");
  if (!fid) {
    printf0("writeCinstructionsC: unable to open file (errno=%d, %s)\n",errno,strerror(errno));
    return -1;
  }    

  printf1("writeCinstructionsC: nInstructions=%lu\n",nInstructions);

  for (int64_t i=0;i<nInstructions;i++,indices++) {
    //printf2("   getInstruction(%d) ... ",indices[0]);
    int rc=getInstruction(indices[0],&type,&nParameters,&parameters,&nOperands,&operands);
    //printf2("   getInstruction(%d)=%d\n",indices[0],rc);
    if (rc<0) {
      return rc; }
    switch (type) {

    case I_set:
      break;

    case I_load:
      fprintf(fid,"\tm[%"PRId64"]=%.20e;//load(%"PRId64")\n",memoryLocations[indices[0]-1]-1,*parameters,indices[0]);
      break;

    case I_sum:
      countFlops_nsum  += (nOperands-1);

      fprintf(fid,"\tm[%"PRId64"]=",memoryLocations[indices[0]-1]-1);
      while (nOperands-->0)
	fprintf(fid,(*parameters++>0)?"+m[%"PRId64"]":"-m[%"PRId64"]",
		memoryLocations[(*(operands++))-1]-1);
      fprintf(fid,";//sum(%"PRId64")\n",indices[0]);
      break;

    case I_sumprod:
      countFlops_nsum  += (parameters[1]-1);
      countFlops_nprod += (parameters[0]-1)*parameters[1];
      
      nSums=0;
      
      fprintf(fid,"\tm[%"PRId64"]=",memoryLocations[indices[0]-1]-1);
      if (parameters[1]<=0) {
	printf0("\nERROR: I_sumprod with 0 sums\n");
	break;
      }
      if (parameters[0]<=0) {
	printf0("\nERROR: I_sumprod with 0 prods\n");
	break;
      }
      int64_t s=parameters[1];
      do {
	int64_t p=parameters[0];
	do {
	  fprintf(fid,"m[%"PRId64"]",memoryLocations[(*(operands++))-1]-1);
	  if (--p>0)
	    fprintf(fid,"*");
	  else break;
	} while (1);
	if (--s>0) {
	    // one large line breaks compiler, so break line every few +'s
	  if ((nSums+=parameters[0])<MAX_TERMS_PERLINE) 
	    fprintf(fid,"+");
	  else {
	    nSums=0;
	    fprintf(fid,";\n\tm[%"PRId64"]+=",memoryLocations[indices[0]-1]-1);
	  }
	}
	else break;
      } while (1);
      fprintf(fid,";//sumprod(%"PRId64")\n",indices[0]);
      break;

    case I_plus_minus_dot_div:
      countFlops_nsum  += (nOperands-1)/2;
      countFlops_nprod += (nOperands-1)/2;
      countFlops_ndiv  ++;
      
      fprintf(fid,"\tm[%"PRId64"]=(m[%"PRId64"]-(",
	      memoryLocations[indices[0]-1]-1,memoryLocations[(*(operands++))-1]-1);
      odiv=memoryLocations[(*(operands++))-1]-1;
      nOperands-=2;
      if (nOperands<=0) {
	printf0("\nERROR: I_plus_minus_dot_div with 0 sums\n");
	break;
      }
      do {
	fprintf(fid,"m[%"PRId64"]*m[%"PRId64"]",
		memoryLocations[operands[0]-1]-1,memoryLocations[operands[1]-1]-1);
	operands+=2;
	nOperands-=2;
	if (nOperands>0) {
	  fprintf(fid,"+");
	} else break;
      } while (1);
      fprintf(fid,"))/m[%"PRId64"];//plus-dot-div(%"PRId64")\n",odiv,indices[0]);
      break;

    case I_minus_dot_div:
      countFlops_nsum  += (nOperands-1)/2;
      countFlops_nprod += (nOperands-1)/2;
      countFlops_ndiv  ++;

      fprintf(fid,"\tm[%"PRId64"]=-(",memoryLocations[indices[0]-1]-1);
      odiv=memoryLocations[(*(operands++))-1]-1;
      nOperands--;
      if (nOperands<=0) {
	printf0("\nERROR: I_minus_dot_div with 0 sums\n");
	break;
      }
      do {
	fprintf(fid,"m[%"PRId64"]*",memoryLocations[(*(operands++))-1]-1);
	fprintf(fid,"m[%"PRId64"]",memoryLocations[(*(operands++))-1]-1);
	nOperands-=2;
	if (nOperands>0) {
	  fprintf(fid,"+");
	} else break;
      } while (1);
      fprintf(fid,")/m[%"PRId64"];//plus-dot-div(%"PRId64")\n",odiv,indices[0]);
      break;
      
    case I_plus_minus_dot:
      countFlops_nsum  += (nOperands-1)/2;
      countFlops_nprod += (nOperands-1)/2;

      fprintf(fid,"\tm[%"PRId64"]=m[%"PRId64"]-(",
	      memoryLocations[indices[0]-1]-1,memoryLocations[(*(operands++))-1]-1);
      nOperands--;
      if (nOperands<=0) {
	printf0("\nERROR: I_plus_minus_dot with 0 sums\n");
	break;
      }
      do {
	fprintf(fid,"m[%"PRId64"]*m[%"PRId64"]",
		memoryLocations[operands[0]-1]-1,memoryLocations[operands[1]-1]-1);
	operands+=2;
	nOperands-=2;
	if (nOperands>0) {
	  fprintf(fid,"+");
	} else break;
      } while (1);
      fprintf(fid,");//plus-dot(%"PRId64")\n",indices[0]);
      break;

    case I_minus_dot:
      countFlops_nsum  += (nOperands-1)/2;
      countFlops_nprod += (nOperands-1)/2;

      fprintf(fid,"\tm[%"PRId64"]=-(",memoryLocations[indices[0]-1]-1);
      if (nOperands<=0) {
	printf0("\nERROR: I_minus_dot with 0 sums\n");
	break;
      }
      do {
	fprintf(fid,"m[%"PRId64"]*",memoryLocations[(*(operands++))-1]-1);
	fprintf(fid,"m[%"PRId64"]",memoryLocations[(*(operands++))-1]-1);
	nOperands-=2;
	if (nOperands>0) {
	  fprintf(fid,"+");
	} else break;
      } while (1);
      fprintf(fid,");//minus-dot(%"PRId64")\n",indices[0]);
      break;

    case I_div:
      countFlops_ndiv  ++;

      fprintf(fid,"\tm[%"PRId64"]=m[%"PRId64"]/m[%"PRId64"];//div\n",
	      memoryLocations[indices[0]-1]-1,
	      memoryLocations[operands[0]-1]-1,memoryLocations[operands[1]-1]-1);
      break;

    case I_plus_sqr:
      countFlops_nsum  += (nOperands-1);
      countFlops_nprod += nOperands;

      fprintf(fid,"\tm[%"PRId64"]=",memoryLocations[indices[0]-1]-1);
      do {
	fprintf(fid,"m[%"PRId64"]*",memoryLocations[operands[0]-1]-1);
	fprintf(fid,"m[%"PRId64"]",memoryLocations[(*(operands++))-1]-1);
	nOperands--;
	if (nOperands>0) {
	  fprintf(fid,"+");
	} else break;
      } while (1);
      fprintf(fid,";//plus-sqr(%"PRId64")\n",indices[0]);
      break;
      
    case I_plus_abs:
      countFlops_nsum  += (nOperands-1);
      countFlops_nabs  += nOperands;

      fprintf(fid,"\tm[%"PRId64"]=",memoryLocations[indices[0]-1]-1);
      do {
	fprintf(fid,"fabs(m[%"PRId64"])",memoryLocations[(*(operands++))-1]-1);
	nOperands--;
	if (nOperands>0) {
	  fprintf(fid,"+");
	} else break;
      } while (1);
      fprintf(fid,";//plus-abs(%"PRId64")\n",indices[0]);
      break;
      
    case I_min:
      countFlops_nsum += (nOperands-1);
      countFlops_nif  += (nOperands-1);
      
      fprintf(fid,"\tm[%"PRId64"]=m[%"PRId64"];\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[(*(operands++))-1]-1);
      nOperands--;
      if (nOperands<=0) {
	printf0("\nERROR: I_min with <=1 value\n");
	break;
      }
      while (nOperands-->0) {
	fprintf(fid,"\tif (m[%"PRId64"]>m[%"PRId64"]) m[%"PRId64"]=m[%"PRId64"];//min(%"PRId64")\n",
		memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,
		memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,indices[0]);
	operands++; }
      break;

    case I_min0:
      countFlops_nsum += nOperands;
      countFlops_nif  += nOperands;
      
      fprintf(fid,"\tm[%"PRId64"]=0;\n",memoryLocations[indices[0]-1]-1);
      if (nOperands<=0) {
	printf0("\nERROR: I_min0 with <=0 value\n");
	break;
      }
      while (nOperands-->0) {
	fprintf(fid,"\tif (m[%"PRId64"]>m[%"PRId64"]) m[%"PRId64"]=m[%"PRId64"];//min0(%"PRId64")\n",
		memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,
		memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,indices[0]);
	operands++; }
      break;

    case I_max:
      countFlops_nsum += (nOperands-1);
      countFlops_nif  += (nOperands-1);
      
      fprintf(fid,"\tm[%"PRId64"]=m[%"PRId64"];\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[(*(operands++))-1]-1);
      nOperands--;
      if (nOperands<=0) {
	printf0("\nERROR: I_max with <=1 value\n");
	break;
      }
      while (nOperands-->0) {
	fprintf(fid,"\tif (m[%"PRId64"]<m[%"PRId64"]) m[%"PRId64"]=m[%"PRId64"];//max(%"PRId64")\n",
		memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,
		memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,indices[0]);
	operands++; }
      break;

    case I_max0:
      countFlops_nsum += nOperands;
      countFlops_nif  += nOperands;
      
      fprintf(fid,"\tm[%"PRId64"]=0;\n",memoryLocations[indices[0]-1]-1);
      if (nOperands<=0) {
	printf0("\nERROR: I_max0 with 0 values\n");
	break;
      }
      while (nOperands-->0) {
	fprintf(fid,"\tif (m[%"PRId64"]<m[%"PRId64"]) m[%"PRId64"]=m[%"PRId64"];//max0(%"PRId64")\n",
		memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,
		memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,indices[0]);
	operands++; }
      break;

    case I_max_abs:
      countFlops_nsum += (nOperands-1);
      countFlops_nif  += (nOperands-1);
      countFlops_nabs  += nOperands;
      
      fprintf(fid,"\tm[%"PRId64"]=fabs(m[%"PRId64"]);\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[(*(operands++))-1]-1);
      nOperands--;
      if (nOperands<=0) {
	printf0("\nERROR: I_max_abs with 1 value\n");
	break;
      }
      while (nOperands-->0) {
	fprintf(fid,"\tif (m[%"PRId64"]<fabs(m[%"PRId64"])) m[%"PRId64"]=fabs(m[%"PRId64"]);//max_abs(%"PRId64")\n",
		memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,
		memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,indices[0]);
	operands++; }
      break;

    case I_clp:
      countFlops_nclp += nOperands;
      
      fprintf(fid,"\tm[%"PRId64"]=DBL_MAX;\n",memoryLocations[indices[0]-1]-1);
      if (nOperands<=0) {
	printf0("\nERROR: I_clp with 0 values\n");
	break;
      }
      while (nOperands>0) {
	fprintf(fid,"\tif (m[%"PRId64"]<0) {SCRATCHBOOK_TYPE x=-m[%"PRId64"]/m[%"PRId64"]; if (m[%"PRId64"]>x) m[%"PRId64"]=x;}//clp(%"PRId64")\n",
		memoryLocations[operands[1]-1]-1,memoryLocations[operands[0]-1]-1,
		memoryLocations[operands[1]-1]-1,
		memoryLocations[indices[0]-1]-1,memoryLocations[indices[0]-1]-1,indices[0]);
	operands+=2;
	nOperands-=2;
      }
      break;
    case I_round:
      countFlops_nround  += nOperands;
      
      fprintf(fid,"\tm[%"PRId64"]=round(m[%"PRId64"]);//round(%"PRId64")\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,indices[0]);
      break;
    case I_ceil:
      countFlops_nceil  += nOperands;
      
      fprintf(fid,"\tm[%"PRId64"]=ceil(m[%"PRId64"]);//ceil(%"PRId64")\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,indices[0]);
      break;
    case I_floor:
      countFlops_nfloor  += nOperands;
      
      fprintf(fid,"\tm[%"PRId64"]=floor(m[%"PRId64"]);//floor(%"PRId64")\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,indices[0]);
      break;
    case I_abs:
      countFlops_nabs  += nOperands;
      
      fprintf(fid,"\tm[%"PRId64"]=fabs(m[%"PRId64"]);//abs(%"PRId64")\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,indices[0]);
      break;
    case I_sign:
      countFlops_nsign  += nOperands;
      
      fprintf(fid,"\tm[%"PRId64"]=(m[%"PRId64"]>=0)?1:-1;//sign(%"PRId64")\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,indices[0]);
      break;
    case I_2:
      fprintf(fid,"\tm[%"PRId64"]=2;//2\n",
	      memoryLocations[indices[0]-1]-1);
      break;
    case I_2times:
      countFlops_nprod++;
      fprintf(fid,"\tm[%"PRId64"]=2*m[%"PRId64"];//2*(%"PRId64")\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,indices[0]);
      break;
    case I_6times:
      countFlops_nprod++;
      fprintf(fid,"\tm[%"PRId64"]=6*m[%"PRId64"];//6*(%"PRId64")\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,indices[0]);
      break;
    case I_sqr:
      countFlops_nprod++;
      fprintf(fid,"\tm[%"PRId64"]=m[%"PRId64"]*m[%"PRId64"];//sqr(%"PRId64")\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,
	      memoryLocations[operands[0]-1]-1,indices[0]);
      break;
    case I_3sqr:
      countFlops_nprod +=2;;
      fprintf(fid,"\tm[%"PRId64"]=3*m[%"PRId64"]*m[%"PRId64"];//sqr(%"PRId64")\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,
	      memoryLocations[operands[0]-1]-1,indices[0]);
      break;
    case I_cube:
      countFlops_nprod +=2;;
      fprintf(fid,"\tm[%"PRId64"]=m[%"PRId64"]*m[%"PRId64"]*m[%"PRId64"];//cube(%"PRId64")\n",
	      memoryLocations[indices[0]-1]-1,
	      memoryLocations[operands[0]-1]-1,memoryLocations[operands[0]-1]-1,
	      memoryLocations[operands[0]-1]-1,indices[0]);
      break;
    case I_sqrt:
      countFlops_nsqrt ++;
      fprintf(fid,"\tm[%"PRId64"]=sqrt(m[%"PRId64"]);//sqrt(%"PRId64")\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,indices[0]);
      break;
    case I_Dsqrt:
      countFlops_nprod++;
      countFlops_npow ++;
      fprintf(fid,"\tm[%"PRId64"]=.5*pow(m[%"PRId64"],-.5);//Dsqrt(%"PRId64")\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,indices[0]);
      break;
    case I_DDsqrt:
      countFlops_nprod++;
      countFlops_npow ++;
      fprintf(fid,"\tm[%"PRId64"]=-.25*pow(m[%"PRId64"],-1.5);//DDsqrt(%"PRId64")\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,indices[0]);
      break;
    case I_inv:
      countFlops_ndiv++;
      fprintf(fid,"\tm[%"PRId64"]=1/m[%"PRId64"];//inv(%"PRId64")\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,indices[0]);
    break;
    case I_minus_inv_sqr:
      countFlops_ndiv ++;
      countFlops_nprod++;      
      fprintf(fid,"\tm[%"PRId64"]=-1/(m[%"PRId64"]*m[%"PRId64"]);//-inv sqr(%"PRId64")\n",
	      memoryLocations[indices[0]-1]-1,
	      memoryLocations[operands[0]-1]-1,memoryLocations[operands[0]-1]-1,indices[0]);
    break;
    case I_2_inv_cube:
      countFlops_ndiv++;
      countFlops_nprod+=2;      
      fprintf(fid,"\tm[%"PRId64"]=2/(m[%"PRId64"]*m[%"PRId64"]*m[%"PRId64"]);//2inv cube(%"PRId64")\n",
	      memoryLocations[indices[0]-1]-1,
	      memoryLocations[operands[0]-1]-1,memoryLocations[operands[0]-1]-1,memoryLocations[operands[0]-1]-1,indices[0]);
    break;
    case I_cos:
      countFlops_ntrig++;
      fprintf(fid,"\tm[%"PRId64"]=cos(m[%"PRId64"]);//cos(%"PRId64")\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,indices[0]);
    break;
    case I_minus_cos:
      countFlops_ntrig++;
      fprintf(fid,"\tm[%"PRId64"]=-cos(m[%"PRId64"]);//-cos(%"PRId64")\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,indices[0]);
    break;
    case I_sin:
      countFlops_ntrig++;
      fprintf(fid,"\tm[%"PRId64"]=sin(m[%"PRId64"]);//sin(%"PRId64")\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,indices[0]);
    break;
    case I_minus_sin:
      countFlops_ntrig++;
      fprintf(fid,"\tm[%"PRId64"]=-sin(m[%"PRId64"]);//-sin(%"PRId64")\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,indices[0]);
    break;
    case I_log:
      countFlops_nlog++;
      fprintf(fid,"\tm[%"PRId64"]=log(m[%"PRId64"]);//log(%"PRId64")\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,indices[0]);
    break;
    case I_exp:
      countFlops_nexp++;
      fprintf(fid,"\tm[%"PRId64"]=exp(m[%"PRId64"]);//exp(%"PRId64")\n",
      	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,indices[0]);
      break;
    /*
    case I_exp:
      fprintf(fid,"\tm[%"PRId64"]=(m[%"PRId64"]<-40)?0:exp(m[%"PRId64"]);//exp\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[indices[0]-1]-1,
	      memoryLocations[operands[0]-1]-1);
    break;
    */
    case I_atan:
      countFlops_ntrig++;
      fprintf(fid,"\tm[%"PRId64"]=atan(m[%"PRId64"]);//atan(%"PRId64")\n",
      	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,indices[0]);
      break;
    case I_Datan:
      countFlops_ndiv++;
      countFlops_nsum++;
      countFlops_nprod++;
      fprintf(fid,"\tm[%"PRId64"]=1/(1+m[%"PRId64"]*m[%"PRId64"]);//Datan(%"PRId64")\n",
      	      memoryLocations[indices[0]-1]-1,
	      memoryLocations[operands[0]-1]-1,memoryLocations[operands[0]-1]-1,indices[0]);
      break;
    case I_DDatan:
      countFlops_ndiv +=2;
      countFlops_nsum  +=2;
      countFlops_nprod +=3;
      fprintf(fid,"\tm[%"PRId64"]=-2*m[%"PRId64"]/(1+m[%"PRId64"]*m[%"PRId64"])/(1+m[%"PRId64"]*m[%"PRId64"]);//DDatan(%"PRId64")\n",
      	      memoryLocations[indices[0]-1]-1,
	      memoryLocations[operands[0]-1]-1,
	      memoryLocations[operands[0]-1]-1,memoryLocations[operands[0]-1]-1,
	      memoryLocations[operands[0]-1]-1,memoryLocations[operands[0]-1]-1,indices[0]);
      break;

    case I_luS2A:
      countFlops_numfpack++;
      for (uint64_t j=0;j<nOperands;j++)
	fprintf(fid,"\tAx[%"PRId64"][%"PRId64"]=m[%"PRId64"];\n",(uint64_t)parameters[0]-1,j,memoryLocations[operands[j]-1]-1);
      fprintf(fid,"\tif (Numeric[%"PRId64"]) umfpack_dl_free_numeric(&Numeric[%"PRId64"]);\n",(uint64_t)parameters[0]-1,(uint64_t)parameters[0]-1);


fprintf(fid,"\t(void)umfpack_dl_symbolic(%"PRId64",%"PRId64",Ap[%"PRId64"],Ai[%"PRId64"],Ax[%"PRId64"],&Symbolic[%"PRId64"],null,null);\n",
	      (uint64_t)parameters[1],(uint64_t)parameters[1],
	      (uint64_t)parameters[0]-1,(uint64_t)parameters[0]-1,(uint64_t)parameters[0]-1,(uint64_t)parameters[0]-1);
      fprintf(fid,"\t(void)umfpack_dl_numeric(Ap[%"PRId64"],Ai[%"PRId64"],Ax[%"PRId64"],Symbolic[%"PRId64"],&Numeric[%"PRId64"],null,null);\n",
	      (uint64_t)parameters[0]-1,(uint64_t)parameters[0]-1,
	      (uint64_t)parameters[0]-1,(uint64_t)parameters[0]-1,(uint64_t)parameters[0]-1);
      fprintf(fid,"\tumfpack_dl_free_symbolic(&Symbolic[%"PRId64"]);\n",(uint64_t)parameters[0]-1);
      fprintf(fid,"\tm[%"PRId64"]=%"PRId64";\n",memoryLocations[indices[0]-1]-1,(uint64_t)parameters[0]-1);
      break;  

    case I_luS2Asym:
      countFlops_numfpack++;
      for (uint64_t j=0;j<nOperands;j++)
	fprintf(fid,"\tAx[%"PRId64"][%"PRId64"]=m[%"PRId64"];\n",(uint64_t)parameters[0]-1,j,memoryLocations[operands[j]-1]-1);
      fprintf(fid,"\tif (Numeric[%"PRId64"]) umfpack_dl_free_numeric(&Numeric[%"PRId64"]);\n",(uint64_t)parameters[0]-1,(uint64_t)parameters[0]-1);


      fprintf(fid,"\tdouble Control [UMFPACK_CONTROL];\n");
      fprintf(fid,"\tumfpack_dl_defaults(Control);\n");
      fprintf(fid,"\tControl[UMFPACK_STRATEGY]=UMFPACK_STRATEGY_SYMMETRIC;\n");
      fprintf(fid,"\t(void)umfpack_dl_symbolic(%"PRId64",%"PRId64",Ap[%"PRId64"],Ai[%"PRId64"],Ax[%"PRId64"],&Symbolic[%"PRId64"],null,null);\n",
	      (uint64_t)parameters[1],(uint64_t)parameters[1],
	      (uint64_t)parameters[0]-1,(uint64_t)parameters[0]-1,(uint64_t)parameters[0]-1,(uint64_t)parameters[0]-1);
      fprintf(fid,"\t(void)umfpack_dl_numeric(Ap[%"PRId64"],Ai[%"PRId64"],Ax[%"PRId64"],Symbolic[%"PRId64"],&Numeric[%"PRId64"],null,null);\n",
	      (uint64_t)parameters[0]-1,(uint64_t)parameters[0]-1,
	      (uint64_t)parameters[0]-1,(uint64_t)parameters[0]-1,(uint64_t)parameters[0]-1);
      fprintf(fid,"\tumfpack_dl_free_symbolic(&Symbolic[%"PRId64"]);\n",(uint64_t)parameters[0]-1);
      fprintf(fid,"\tm[%"PRId64"]=%"PRId64";\n",memoryLocations[indices[0]-1]-1,(uint64_t)parameters[0]-1);
      break;  

    case I_mldivideA2F1:
      countFlops_numfpack++;
      fprintf(fid,"\t{long atomicID=m[%"PRId64"];\n",memoryLocations[operands[0]-1]-1);
      for (uint64_t j=1;j<nOperands;j++)
	fprintf(fid,"\tb[atomicID][%"PRId64"]=m[%"PRId64"];\n",j-1,memoryLocations[operands[j]-1]-1);
      fprintf(fid,"\t(void)umfpack_dl_solve(UMFPACK_A,Ap[atomicID],Ai[atomicID],Ax[atomicID],x[atomicID],b[atomicID],Numeric[atomicID],null,null);\n");
      fprintf(fid,"\tm[%"PRId64"]=x[atomicID][0]; }\n",memoryLocations[indices[0]-1]-1);
      break;
      
    case I_mldivideA2Fn:
      countFlops_numfpack++;
      fprintf(fid,"\tm[%"PRId64"]=x[(long)m[%"PRId64"]][%"PRId64"];\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1,(uint64_t)parameters[0]-1);
      break;
      
    default:
      fprintf(fid,"\t//unknown instructions type %d\n",type);
      printf0("instructionsTable: unknown instructions type %d\n",type);
      break;
    }
    //fprintf(fid,"\tif (isnan(m[%"PRId64"])) printf(\"NAN computing scrapbook entry %"PRId64"\\n\");\n",memoryLocations[indices[0]-1]-1,indices[0]);
    //fprintf(fid,"\tif (isinf(m[%"PRId64"])) printf(\"INFINITY computing scrapbook entry %"PRId64"\\n\");\n",memoryLocations[indices[0]-1]-1,indices[0]);
  }

  printf1("writeCinstructionsC: done\n");

  fclose(fid);

  return 0;
}

/*******************************************************************************
 * Write ASM code corresponding to an array of instructions
 * Returns -1 if any index is not valid, otherwise returns 0
 *******************************************************************************/
EXPORT int writeAsmInstructionsC(/* inputs */
			       int64_t *indices,          // indices of instructions to write
			       int64_t *memoryLocations,  // memory locations for all instructions
			       /* sizes */
			       mwSize nInstructions,     // # of instructions to write
			       mwSize NInstructions)     // total # of instructions
{
  instructionType_t type;
  int64_t nParameters; 
  double  *parameters;
  int64_t nOperands;   
  int64_t *operands;
  int64_t odiv,oplus,ic=0;
  int anyConstant=0;
  int64_t label;


  FILE *fid=fopen("tmp_toremove.c","w");
  if (!fid) {
    printf0("writeAsmInstructionsC: unable to open file (errno=%d, %s)\n",errno,strerror(errno));
    return -1;
  }    

  //printf0("writeAsmInstructionsC: nInstructions=%"PRId64"\n",nInstructions);

  fprintf(fid,"\t{");
  // Store constants in memory
  for (int64_t i=0;i<nInstructions;i++) {
    //printf0("   getInstruction(%d) ... ",indices[0]);
    int rc=getInstruction(indices[i],&type,&nParameters,&parameters,&nOperands,&operands); //1-based operands
    //printf0("   getInstruction(%d)=%d\n",indices[0],rc);
    if (rc<0) {
      return rc; }
    if (type==I_load) {
      if ((anyConstant++)) 
	fprintf(fid,",%.20e",parameters[0]);
      else
	fprintf(fid," double c[]={%.20e",parameters[0]);
    }
  }
  if (anyConstant) {
    fprintf(fid,"};\n");
  }
  fprintf(fid,"\tasm(\n");


#define SIZE 8

#define NEWINSTRUCTION()              fprintf(fid,"\t\"")
#define CLEARREG0()                   fprintf(fid,"subsd %%%%xmm0,%%%%xmm0;")
#define CLEARREG1()                   fprintf(fid,"subsd %%%%xmm1,%%%%xmm1;")
#define CLEARREG2()                   fprintf(fid,"subsd %%%%xmm2,%%%%xmm2;")
#define MOVEvalue2REG0(index)         fprintf(fid,"movsd %"PRId64"(%%[c]),%%%%xmm0;",SIZE*(index))
#define MOVEinfty2REG0()              fprintf(fid,"movsd %%[dbl_max],%%%%xmm0;")
#define MOVESCRAP2REG0(index)         fprintf(fid,"movsd %"PRId64"(%%[m]),%%%%xmm0;",SIZE*(index))
#define MOVESCRAP2REG1(index)         fprintf(fid,"movsd %"PRId64"(%%[m]),%%%%xmm1;",SIZE*(index))
#define MOVEREG02SCRAP(index,comment) fprintf(fid,"movsd %%%%xmm0,%"PRId64"(%%[m]);\"//%s\n",SIZE*(index),comment)
#define MOVEREG12SCRAP(index,comment) fprintf(fid,"movsd %%%%xmm1,%"PRId64"(%%[m]);\"//%s\n",SIZE*(index),comment)

#define MOVEREG22REG0()               fprintf(fid,"movsd %%%%xmm2,%%%%xmm0;")
#define MOVEREG22REG1()               fprintf(fid,"movsd %%%%xmm2,%%%%xmm1;")

#define ADDREG12REG0()                fprintf(fid,"addsd %%%%xmm1,%%%%xmm0;")
#define ADDSCRAP2REG0(index)          fprintf(fid,"addsd %"PRId64"(%%[m]),%%%%xmm0;",SIZE*(index))
#define SUBREG02REG1()                fprintf(fid,"subsd %%%%xmm0,%%%%xmm1;")
#define SUBREG12REG0()                fprintf(fid,"subsd %%%%xmm1,%%%%xmm0;")
#define SUBREG12REG2()                fprintf(fid,"subsd %%%%xmm1,%%%%xmm2;")
#define SUBREG02REG2()                fprintf(fid,"subsd %%%%xmm0,%%%%xmm2;")
#define SUBSCRAP2REG0(index)          fprintf(fid,"subsd %"PRId64"(%%[m]),%%%%xmm0;",SIZE*(index))
#define MULREG02REG0()                fprintf(fid,"mulsd %%%%xmm0,%%%%xmm0;")
#define MULREG12REG1()                fprintf(fid,"mulsd %%%%xmm1,%%%%xmm1;")
#define MULSCRAP2REG0(index)          fprintf(fid,"mulsd %"PRId64"(%%[m]),%%%%xmm0;",SIZE*(index))
#define MULSCRAP2REG1(index)          fprintf(fid,"mulsd %"PRId64"(%%[m]),%%%%xmm1;",SIZE*(index))
#define DIVREG1BYSCRAP2REG1(index)    fprintf(fid,"divsd %"PRId64"(%%[m]),%%%%xmm1;",SIZE*(index))
#define DIVREG1BYREG22REG1()          fprintf(fid,"divsd %%%%xmm2,%%%%xmm1;");

#define MINSCRAP2REG0(index)          fprintf(fid,"minsd %"PRId64"(%%[m]),%%%%xmm0;",SIZE*(index)) 
#define MAXSCRAP2REG0(index)          fprintf(fid,"maxsd %"PRId64"(%%[m]),%%%%xmm0;",SIZE*(index)) 
#define MINREG12REG0()                fprintf(fid,"minsd %%%%xmm1,%%%%xmm0;");
#define MAXREG12REG0()                fprintf(fid,"maxsd %%%%xmm1,%%%%xmm0;") 

#define CMPREG22REG1()                fprintf(fid,"ucomisd %%%%xmm2,%%%%xmm1;")
#define CMPREG22REG0()                fprintf(fid,"ucomisd %%%%xmm2,%%%%xmm0;")



  // Write ASM code
  for (int64_t i=0;i<nInstructions;i++,indices++) {
    //printf0("   getInstruction(%d) ... ",indices[0]);
    int rc=getInstruction(indices[0],&type,&nParameters,&parameters,&nOperands,&operands); //1-based operands
    //printf0("   getInstruction(%d)=%d\n",indices[0],rc);
    if (rc<0) {
      return rc; }
    switch (type) {

    case I_set:
      break;

    case I_load:
      NEWINSTRUCTION();
      MOVEvalue2REG0(ic);
      MOVEREG02SCRAP(memoryLocations[indices[0]-1]-1,"load");
      ic=ic+1;
      break;

    case I_sum:
      NEWINSTRUCTION();
      if (*parameters<0)
	CLEARREG0();
      else {
	MOVESCRAP2REG0(memoryLocations[(*(operands++))-1]-1);
	parameters++;
	nOperands--;
      }
      while (nOperands-->0)
	if (*parameters++>0)
	  ADDSCRAP2REG0(memoryLocations[(*(operands++))-1]-1);
	else
	  SUBSCRAP2REG0(memoryLocations[(*(operands++))-1]-1);
      MOVEREG02SCRAP(memoryLocations[indices[0]-1]-1,"sum");
      break;

    case I_sumprod:
      NEWINSTRUCTION();
      int64_t s=parameters[1];
      int64_t p=parameters[0];
      MOVESCRAP2REG0(memoryLocations[(*(operands++))-1]-1);
      while (--p>0)
	MULSCRAP2REG0(memoryLocations[(*(operands++))-1]-1);
      while (--s>0) {
	int64_t p=parameters[0];
	MOVESCRAP2REG1(memoryLocations[(*(operands++))-1]-1);
	while (--p>0)
	  MULSCRAP2REG1(memoryLocations[(*(operands++))-1]-1);
	ADDREG12REG0();
      };
      MOVEREG02SCRAP(memoryLocations[indices[0]-1]-1,"sumprod");
      break;
      
    case I_plus_minus_dot_div:
      NEWINSTRUCTION();
      oplus=memoryLocations[(*(operands++))-1]-1;
      odiv=memoryLocations[(*(operands++))-1]-1;
      nOperands-=2;
      MOVESCRAP2REG0(memoryLocations[(*(operands++))-1]-1);
      MULSCRAP2REG0(memoryLocations[(*(operands++))-1]-1);
      while ((nOperands-=2)>0) {
	MOVESCRAP2REG1(memoryLocations[(*(operands++))-1]-1);
	MULSCRAP2REG1(memoryLocations[(*(operands++))-1]-1);
	ADDREG12REG0();
      }
      MOVESCRAP2REG1(oplus);
      SUBREG02REG1();
      DIVREG1BYSCRAP2REG1(odiv);
      MOVEREG12SCRAP(memoryLocations[indices[0]-1]-1,"plus_minus_dot_div");
      break;

    case I_minus_dot_div:
      NEWINSTRUCTION();
      odiv=memoryLocations[(*(operands++))-1]-1;
      nOperands--;
      MOVESCRAP2REG0(memoryLocations[(*(operands++))-1]-1);
      MULSCRAP2REG0(memoryLocations[(*(operands++))-1]-1);
      while ((nOperands-=2)>0) {
	MOVESCRAP2REG1(memoryLocations[(*(operands++))-1]-1);
	MULSCRAP2REG1(memoryLocations[(*(operands++))-1]-1);
	ADDREG12REG0(); }
      CLEARREG1();
      SUBREG02REG1();
      DIVREG1BYSCRAP2REG1(odiv);
      MOVEREG12SCRAP(memoryLocations[indices[0]-1]-1,"minus_dot_div");
      break;

    case I_plus_minus_dot:
      NEWINSTRUCTION();
      MOVESCRAP2REG0(memoryLocations[(*(operands++))-1]-1);
      nOperands++; // since -=2 will return value post subtraction
      while ((nOperands-=2)>0) {
	MOVESCRAP2REG1(memoryLocations[(*(operands++))-1]-1);
	MULSCRAP2REG1(memoryLocations[(*(operands++))-1]-1);
	SUBREG12REG0(); }
      MOVEREG02SCRAP(memoryLocations[indices[0]-1]-1,"plus_minus_dot");
      break;

    case I_minus_dot:
      NEWINSTRUCTION();
      CLEARREG0();
      nOperands+=2; // since -=2 will return value post subtraction
      while ((nOperands-=2)>0) {
	MOVESCRAP2REG1(memoryLocations[(*(operands++))-1]-1);
	MULSCRAP2REG1(memoryLocations[(*(operands++))-1]-1);
	SUBREG12REG0(); }
      MOVEREG02SCRAP(memoryLocations[indices[0]-1]-1,"minus_dot");
      break;

    case I_div:
      NEWINSTRUCTION();
      MOVESCRAP2REG1(memoryLocations[(*(operands++))-1]-1);
      DIVREG1BYSCRAP2REG1(memoryLocations[(*(operands++))-1]-1);
      MOVEREG12SCRAP(memoryLocations[indices[0]-1]-1,"div");
      break;

    case I_plus_sqr:
      NEWINSTRUCTION();
      MOVESCRAP2REG0(memoryLocations[(*(operands++))-1]-1);
      MULREG02REG0();
      while (--nOperands>0) {
	MOVESCRAP2REG1(memoryLocations[(*(operands++))-1]-1);
	MULREG12REG1();
	ADDREG12REG0(); }
      MOVEREG02SCRAP(memoryLocations[indices[0]-1]-1,"plus_sqr");
      break;
      
    case I_min:
      NEWINSTRUCTION();
      MOVESCRAP2REG0(memoryLocations[(*(operands++))-1]-1);
      while (--nOperands>0)
	MINSCRAP2REG0(memoryLocations[(*(operands++))-1]-1);
      MOVEREG02SCRAP(memoryLocations[indices[0]-1]-1,"min");
      break;

    case I_min0:
      NEWINSTRUCTION();
      CLEARREG0();
      while (nOperands-- >0)
	MINSCRAP2REG0(memoryLocations[(*(operands++))-1]-1);
      MOVEREG02SCRAP(memoryLocations[indices[0]-1]-1,"min0");
      break;

    case I_max:
      NEWINSTRUCTION();
      MOVESCRAP2REG0(memoryLocations[(*(operands++))-1]-1);
      while (--nOperands>0)
	MAXSCRAP2REG0(memoryLocations[(*(operands++))-1]-1);
      MOVEREG02SCRAP(memoryLocations[indices[0]-1]-1,"max");
      break;

    case I_max_abs:
      NEWINSTRUCTION();
      label=0;
      MOVESCRAP2REG0(memoryLocations[(*(operands++))-1]-1);
      CLEARREG2();
      CMPREG22REG0();  // REG0>=0? jump
      fprintf(fid,"jae %"PRId64"%"PRId64"f;",SIZE*(memoryLocations[indices[0]-1]-1),label);
      SUBREG02REG2();  // REG2=-REG0
      MOVEREG22REG0(); // REG0=REG2
      fprintf(fid,"\\\n\t\t%"PRId64"%"PRId64":",SIZE*(memoryLocations[indices[0]-1]-1),label++);
      while (--nOperands>0) {
	MOVESCRAP2REG1(memoryLocations[(*(operands++))-1]-1);
	CLEARREG2();
	CMPREG22REG1();  // REG1>=0? jump
	fprintf(fid,"jae %"PRId64"%"PRId64"f;",SIZE*(memoryLocations[indices[0]-1]-1),label);
	SUBREG12REG2();  // REG2=-REG1
	MOVEREG22REG1(); // REG1=REG2
	fprintf(fid,"\\\n\t\t%"PRId64"%"PRId64":",SIZE*(memoryLocations[indices[0]-1]-1),label++);
	MAXREG12REG0();
      }
      MOVEREG02SCRAP(memoryLocations[indices[0]-1]-1,"max_abs");
      break;

    case I_max0:
      NEWINSTRUCTION();
      CLEARREG0();
      while (nOperands-- >0)
	MAXSCRAP2REG0(memoryLocations[(*(operands++))-1]-1);
      MOVEREG02SCRAP(memoryLocations[indices[0]-1]-1,"max0");
      break;

    case I_clp:
      NEWINSTRUCTION();
      MOVEinfty2REG0();
      label=0;
      do {
	fprintf(fid,"\\\n\t\t%"PRId64"%"PRId64":",SIZE*(memoryLocations[indices[0]-1]-1),label);
	MOVESCRAP2REG1(memoryLocations[*(operands+1)-1]-1);
	CLEARREG2();
	CMPREG22REG1();  // REG1>=0? jump
	fprintf(fid,"jae %"PRId64"%"PRId64"f;",SIZE*(memoryLocations[indices[0]-1]-1),++label);
	SUBREG12REG2();  // REG2=-REG1
	MOVESCRAP2REG1(memoryLocations[*(operands)-1]-1);
	DIVREG1BYREG22REG1(); 
	MINREG12REG0();
	operands+=2;
      } while ((nOperands-=2) >0);
      fprintf(fid,"\\\n\t\t%"PRId64"%"PRId64":",SIZE*(memoryLocations[indices[0]-1]-1),label);
      MOVEREG02SCRAP(memoryLocations[indices[0]-1]-1,"clp");
      break;
      
    case I_abs:
      fprintf(fid,"\tm[%"PRId64"]=fabs(m[%"PRId64"]);//abs\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1);
      break;
    case I_sqrt:
      fprintf(fid,"\tm[%"PRId64"]=sqrt(m[%"PRId64"]);//sqrt\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1);
      break;
    case I_Dsqrt:
      fprintf(fid,"\tm[%"PRId64"]=.5*pow(m[%"PRId64"],-.5);//Dsqrt\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1);
      break;
    case I_DDsqrt:
      fprintf(fid,"\tm[%"PRId64"]=-.25*pow(m[%"PRId64"],-1.5);//DDsqrt\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1);
      break;
    case I_inv:
      fprintf(fid,"\tm[%"PRId64"]=1/m[%"PRId64"];//inv\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1);
    break;
    case I_minus_inv_sqr:
      fprintf(fid,"\tm[%"PRId64"]=-1/(m[%"PRId64"]*m[%"PRId64"]);//-inv-sqr\n",
	      memoryLocations[indices[0]-1]-1,
	      memoryLocations[operands[0]-1]-1,memoryLocations[operands[0]-1]-1);
    break;
    case I_cos:
      fprintf(fid,"\tm[%"PRId64"]=cos(m[%"PRId64"]);//cos\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1);
    break;
    case I_minus_cos:
      fprintf(fid,"\tm[%"PRId64"]=-cos(m[%"PRId64"]);//-cos\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1);
    break;
    case I_sin:
      fprintf(fid,"\tm[%"PRId64"]=sin(m[%"PRId64"]);//sin\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1);
    break;
    case I_minus_sin:
      fprintf(fid,"\tm[%"PRId64"]=-sin(m[%"PRId64"]);//-sin\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1);
    break;
    case I_log:
      fprintf(fid,"\tm[%"PRId64"]=log(m[%"PRId64"]);//log\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1);
    break;
    case I_exp:
      fprintf(fid,"\tm[%"PRId64"]=exp(m[%"PRId64"]);//exp\n",
	      memoryLocations[indices[0]-1]-1,memoryLocations[operands[0]-1]-1);
    break;
    

    default:
      fprintf(fid,"\t//unknown instructions type %d\n",type);
      printf0("instructionsTable: unknown instructions type %d\n",type);
      break;
    }
  }

  if (anyConstant) 
    fprintf(fid,"\t::[m] \"r\" (scratchbook), [dbl_max] \"m\" (dbl_max), [c] \"r\" (c) :\"xmm0\",\"xmm1\",\"xmm2\");\n\t};\n");
  else
    fprintf(fid,"\t::[m] \"r\" (scratchbook), [dbl_max] \"m\" (dbl_max) :\"xmm0\",\"xmm1\",\"xmm2\");\n\t};\n");
  
  //printf0("writeAsmInstructionsC: done\n");

  fclose(fid);
  
  return 0;
}

/*******************************************************************************
 * Main for testing
 *******************************************************************************/

// ! gcc -I/Applications/MATLAB_R2016b.app/extern/include -Ofast -DIGNORE_MEX instructionsTableUTHash.c
	       
int main()
{
  /*
#define Ncalls 200LL
#define maxNparameters 10
#define maxParameters 1
#define maxNoperands 3
#define maxOperands 1
#define type (instructionType_t)23
  */
#define Ncalls 1000000LL
#define maxNparameters 10
#define maxParameters 1
#define maxNoperands 30
#define maxOperands 1
#define type (instructionType_t)23


  double parameters[maxNparameters];
  int64_t operands[maxNoperands];
  int64_t nParameters;
  int64_t nOperands;
  int64_t i,j;
  clock_t t0;

  /* Time appendInstruction */
  if (1) {
    initInstructionsTable();
    t0=clock();
    for (i=0;i<Ncalls;i++) {
      // generate parameters
      nParameters=floor((double)(maxNparameters+1)*((double)rand()/RAND_MAX));
      for (j=0;j<nParameters;j++)
	parameters[j]=floor((double)(maxParameters+1)*((double)rand()/RAND_MAX));
      // generate operands
      nOperands=floor((double)(maxNoperands+1)*((double)rand()/RAND_MAX));
      for (j=0;j<nOperands;j++)
	operands[j]=floor((double)(maxOperands+1)*((double)rand()/RAND_MAX));
      // 
      j=appendInstruction(type,nParameters,parameters,nOperands,operands);
    }
    printf0("Called appendInstruction       %12"PRId64" times: table size =%12"PRId64", %.3f ms, %.3f us/call\n",
	    Ncalls,instructionsTableHeight(),((double)clock()-t0)/CLOCKS_PER_SEC*1e3,((double)clock()-t0)/CLOCKS_PER_SEC*1e6/Ncalls);
  }
  
  /* Time appendUniqueInstruction */
  if (1) {
    initInstructionsTable();
    t0=clock();
    for (i=0;i<Ncalls;i++) {
      // generate parameters
      nParameters=floor((double)(maxNparameters+1)*((double)rand()/RAND_MAX));
      for (j=0;j<nParameters;j++)
	parameters[j]=floor((double)(maxParameters+1)*((double)rand()/RAND_MAX));
      // generate operands
      nOperands=floor((double)(maxNoperands+1)*((double)rand()/RAND_MAX));
      for (j=0;j<nOperands;j++)
	operands[j]=floor((double)(maxOperands+1)*((double)rand()/RAND_MAX));
      // 
      j=appendUniqueInstruction(type,nParameters,parameters,nOperands,operands);
    }
    printf0("Called appendUniqueInstruction %12"PRId64" times: table size =%12"PRId64", %.3f ms, %.3f us/call\n",
	   Ncalls,instructionsTableHeight(),((double)clock()-t0)/CLOCKS_PER_SEC*1e3,((double)clock()-t0)/CLOCKS_PER_SEC*1e6/Ncalls);
  }

  /* Time sortInstructions */
  printf0("Calling sortInstructions: table size =%12"PRId64" (%"PRId64", %"PRId64", %"PRId64", %"PRId64",...)\n",
	 instructionsTableHeight(),
	 instructionsTable.sortedIndices[0],instructionsTable.sortedIndices[1],
	 instructionsTable.sortedIndices[2],instructionsTable.sortedIndices[3]
	 );
  t0=clock();
  sortInstructions();
  printf0("Called  sortInstructions: table size =%12"PRId64" (%"PRId64", %"PRId64", %"PRId64", %"PRId64",...), %.3f ms, %.3f us/call\n",
	 instructionsTableHeight(),
	 instructionsTable.sortedIndices[0],instructionsTable.sortedIndices[1],
	 instructionsTable.sortedIndices[2],instructionsTable.sortedIndices[3],
	 ((double)clock()-t0)/CLOCKS_PER_SEC*1e3,((double)clock()-t0)/CLOCKS_PER_SEC*1e6/Ncalls);

  initInstructionsTable();
  exit(0);
}
