/*
  This file is part of Tencalc.

  Copyright (C) 2012-21 The Regents of the University of California
  (author: Dr. Joao Hespanha).  All rights reserved.
*/

#include "mex.h"
#include <stdint.h>
#include <string.h>

#include "instructionsTableTypes.h"

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

/*******************************************************************************
 * Initializes the instructions table
 *******************************************************************************/
EXPORT void initInstructionsTable();

/*******************************************************************************
 * Appends instruction at the end of the table, and return index to it. a
 * Returns -1 is either buffer is full.
 *******************************************************************************/
EXPORT int64_t appendInstruction(instructionType_t type,
                                 int64_t nParameters,
                                 double *parameters,
                                 int64_t nOperands,
                                 int64_t *operands);

/*******************************************************************************
 * Search for instruction, if already in the table, return index,
 * otherwise append at the end of the table and return the index of the new instruction.
 * Returns -1 is either buffer is full.
 *******************************************************************************/
EXPORT int64_t appendUniqueInstruction(instructionType_t type,
                                       int64_t nParameters,
                                       double *parameters,
                                       int64_t nOperands,
                                       int64_t *operands);

/*******************************************************************************
 * Get instruction, given its index.
 * Returns -1 if the index is not valid
 *******************************************************************************/
EXPORT int getInstruction(int64_t index,
                          instructionType_t *type, // return type by reference
                          int64_t *nParameters,    // return nParameters by reference
                          double **parameters,     // return pointer to parameters by reference
                          int64_t *nOperands,      // return nOperands by reference
                          int64_t **operands);     // return pointer to operands by reference

/*******************************************************************************
 * get dependencies
 *******************************************************************************/
EXPORT void getDependencies(int64_t *children, // array of children (dependents)
                            int64_t *parents); // array of parents

/*******************************************************************************
 * find instructions of a given type.
 *******************************************************************************/
EXPORT void findInstructionsByType(instructionType_t type, // return type by reference
                                   int8_t *isOfType);      // boolean array where indices will be stored

/*******************************************************************************
 * Write C code corresponding to an array of instructions
 * Returns -1 if any index is not valid
 *******************************************************************************/
EXPORT int writeCinstructionsC(                               /* inputs */
                               char *filename,                // filename
                               int64_t *indices,              // indices of instructions to write
                               int64_t *memoryLocations,      // memory locations for
                                                              // all instructions
                               int64_t *minInstructions4Loop, // minimum number of instructions implemented as a loop

                               /* outputs */
                               int64_t *profile, // array with instruction counts

                               /* sizes */
                               mwSize nameLength,     // filename length
                               mwSize nInstructions,  // # of instructions to write
                               mwSize NInstructions); // total # of instructions

/*******************************************************************************
 * Write ASM code corresponding to an array of instructions
 * Returns -1 if any index is not valid, otherwise returns 0
 *******************************************************************************/
EXPORT int writeAsmInstructionsC(                          /* inputs */
                                 char *filename,           // filename
                                 int64_t *indices,         // indices of instructions to write
                                 int64_t *memoryLocations, // memory locations for all instructions
                                 /* sizes */
                                 mwSize nameLength,     // filename length
                                 mwSize nInstructions,  // # of instructions to write
                                 mwSize NInstructions); // total # of instructions
