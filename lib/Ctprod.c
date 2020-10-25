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

/*=================================================================
 * Takes a list of N-dimensional arrays of doubles and returns their tensor product.
 *
 * This is a MEX-file for MATLAB.
 *============================================================*/

/* $Revision: 1.0 $ */

#include "mex.h"

/* If you are using a compiler that equates NaN to zero, you must
 * compile this example using the flag -DNAN_EQUALS_ZERO. For example:
 *
 *     mex -DNAN_EQUALS_ZERO findnz.c
 *
 * This will correctly define the IsNonZero macro for your Compiler */

#if NAN_EQUALS_ZERO
#define IsNonZero(d) ((d)!=0.0 || mxIsNaN(d))
#else
#define IsNonZero(d) ((d)!=0.0)
#endif

void checkTprodsizes(int nlhs,       mxArray *plhs[],
		     int nrhs, const mxArray *prhs[],
		     mwSize  **tprod_index_sizes,
		     size_t *number_of_dimensions_tprod,
		     size_t *number_of_dimensions_sum,
		     size_t *number_of_dimensions_total,
		     size_t *indices[]
		     )
{
  mwSize number_of_dims_object,length_index;
  const mwSize *dim_array_index, *sizes_object;
  double *data_index;
  unsigned int i,j,k;

  /* Cumpute number_of_dimensions_total, number_of_dimensions_tprod */
  *number_of_dimensions_tprod=0;
  *number_of_dimensions_total=0;
  *number_of_dimensions_sum=0;

  for (i=0;i<nrhs-1;i+=2) {
    /* Get the number of dimensions in the input argument.*/
    number_of_dims_object=mxGetNumberOfDimensions(prhs[i]);
    sizes_object=mxGetDimensions(prhs[i]);
    length_index=mxGetNumberOfElements(prhs[i+1]);
    data_index=mxGetPr(prhs[i+1]);

    /* remove trailing singleton dimensions */
    while (number_of_dims_object>0 &&
	   sizes_object[number_of_dims_object-1]==1) {
      number_of_dims_object--;
    }

    if (number_of_dims_object!=length_index) {
      mexErrMsgIdAndTxt( "MATLAB:Ctprod:invalidInput",
			 "tprod requires each indice to have the same number of entries as the size of the corresponding object.");
    }

    for (j=0;j< length_index;j++) {
      if (data_index[j]>0) {
   	if (data_index[j]>0 && data_index[j] > (*number_of_dimensions_tprod)) {
	  *number_of_dimensions_tprod=data_index[j]; }
      } else {
	if (data_index[j]<0 && -data_index[j] > (*number_of_dimensions_sum)) {
	  *number_of_dimensions_sum=-data_index[j]; }
      }
      // mexPrintf("factor %d, data_index[%d]=%3.0f : number_of_dimensions_tprod=%d, number_of_dimensions_sums=%d\n",i/2,j,data_index[j],*number_of_dimensions_tprod,*number_of_dimensions_sum);

    }
  }

  *number_of_dimensions_total = (*number_of_dimensions_sum) + (*number_of_dimensions_tprod);

  /* allocate tprod_index_sizes[]*/
  if (*number_of_dimensions_total>0) {
    *tprod_index_sizes = mxCalloc(*number_of_dimensions_total,sizeof(mwSize));
  } else {
    *tprod_index_sizes = 0L;
  }

  /* fill in dimensions */
  for (i=0;i<nrhs-1;i+=2) {
    /* Get the number of dimensions in the input argument.*/
    number_of_dims_object=mxGetNumberOfDimensions(prhs[i]);
    sizes_object=mxGetDimensions(prhs[i]);
    length_index=mxGetNumberOfElements(prhs[i+1]);
    data_index=mxGetPr(prhs[i+1]);
    indices[i/2]=mxCalloc(length_index,sizeof(size_t));

    for (j=0;j< length_index;j++) {
      if (data_index[j]<0) {
	indices[i/2][j] = -data_index[j]-1;
      } else {
	indices[i/2][j] = data_index[j]+(*number_of_dimensions_sum)-1;
      }
      if ((*tprod_index_sizes)[indices[i/2][j]]==0) {
	(*tprod_index_sizes)[indices[i/2][j]]=sizes_object[j];
      }
      else {
	if ((*tprod_index_sizes)[indices[i/2][j]]!=sizes_object[j]) {
	  for (k=0;k<*number_of_dimensions_total;k++) {
	    mexPrintf("%d,",(*tprod_index_sizes)[k]);
	  }
	  mexPrintf("\n");
	  mexErrMsgIdAndTxt( "MATLAB:Ctprod:invalidInput",
			     "incompatible size in factor %d, index %d (%d~=%d)",i/2,j,(*tprod_index_sizes)[indices[i/2][j]],sizes_object[j]);
	}
      }
      // mexPrintf("factor %d, data_index[%d]=%3.0f -> %d : sizes_object[%d]=%d\n",i/2,j,data_index[j],indices[i/2][j],indices[i/2][j],sizes_object[j]);
    }
  }

  for (i=0;i<*number_of_dimensions_total;i++) {
    if ((*tprod_index_sizes)[i]==0) {
      mexErrMsgIdAndTxt( "MATLAB:Ctprod:invalidInput",
			 "output has dimension with unkown size");
    }
  }
}

void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
    /* Declare variables */
  mwSize number_of_dims_object,length_index;
  size_t elements;
  mwSize j,cmplx;
  mwSize number_of_dims;
  double *pr, *pi, *pind,*output_p;
  const mwSize *dim_array;
  const mwSize *dim_array_index, *sizes_object;
  mwSize  *tprod_index_sizes,*entry;
  size_t number_of_dimensions_tprod;
  size_t number_of_dimensions_sums;
  size_t number_of_dimensions_total;
  size_t output_k,input_k,L,numel;
  unsigned int k,i;
  double *data_index,*data_object,S;
  size_t *indices[nrhs/2+1];

  /* Check for proper number of input and output arguments */
  if (nrhs % 2 == 1) {
    mexErrMsgIdAndTxt( "MATLAB:Ctprod:invalidNumInputs",
		       "Even number of inputs arguments required.");
  }
  if (nlhs > 1){
    mexErrMsgIdAndTxt( "MATLAB:Ctprod:maxlhs",
		       "Too many output arguments.");
  }

  /* Check data type of input argument */
  if (!(mxIsDouble(prhs[0]))) {
    mexErrMsgIdAndTxt( "MATLAB:Ctprod:invalidInputType",
		       "Input array must be of type double.");
  }

  /* Check argument sizes and get tprod() sizes */
  checkTprodsizes(nlhs,plhs,nrhs,prhs,
		  &tprod_index_sizes,
		  &number_of_dimensions_tprod,
		  &number_of_dimensions_sums,
		  &number_of_dimensions_total,
		  indices);

  /* mexPrintf("number_of_dimensions_tprod=%d\n",number_of_dimensions_tprod); */
  /* mexPrintf("number_of_dimensions_sums =%d\n",number_of_dimensions_sums); */
  /* mexPrintf("number_of_dimensions_total=%d\n",number_of_dimensions_total); */
  /* for (k=0;k<number_of_dimensions_total;k++) { */
  /*   mexPrintf("  size[%d]=%d\n",k,tprod_index_sizes[k]);  */
  /* } */

  /* Create output array */
  if (number_of_dimensions_tprod>0) {
    plhs[0]=mxCreateNumericArray(number_of_dimensions_tprod,tprod_index_sizes+number_of_dimensions_sums,mxDOUBLE_CLASS,0);
  } else {
    // mexPrintf("Scalar output\n");
    plhs[0]=mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,0);
  }
  output_p=mxGetPr(plhs[0]);
  if (number_of_dimensions_total>0) {
    entry = mxCalloc(number_of_dimensions_total,sizeof(mwSize)); }

  /* Compute product */

  /* loop over entries */
  numel=mxGetNumberOfElements(plhs[0]);
  output_k=0;
  while (1) {
    /* display current entry */
    /* mexPrintf("Computing entry["); */
    /* for (k=0;k<number_of_dimensions_total;k++) { */
    /*   mexPrintf("i%d=%d,",k,entry[k]); */
    /* } */
    /* mexPrintf("] @ %d\n",output_k); */
    S=1;

    for (i=0;i<nrhs-1;i+=2) {
      sizes_object=mxGetDimensions(prhs[i]);
      data_object=mxGetPr(prhs[i]);
      length_index=mxGetNumberOfElements(prhs[i+1]);
      data_index=mxGetPr(prhs[i+1]);

      input_k=0;
      /* mexPrintf("  input[%d]: [",i/2); */
      /* for (j=0;j< length_index;j++) {  */
      /* 	mexPrintf("%d,",indices[i/2][j]); */
      /* } */
      /* mexPrintf(" ] -> ");  */
      L=1;
      for (j=0;j< length_index;j++) {
	input_k += L*entry[indices[i/2][j]];
	// mexPrintf(" i%d=%d (L=%d); ", indices[i/2][j],entry[indices[i/2][j]],L);
	L *= sizes_object[j];
      }
      S *= data_object[input_k];
      // mexPrintf("\n    data[input_k=%d] = %3.0f, (S=%8f) \n",input_k,data_object[input_k],S);
    }
    output_p[output_k] += S;
    // mexPrintf("  output[output_k=%d] = %3.0f, (+=%8f)\n",output_k,output_p[output_k],S);

    /* next entry ... */
    k=0;
    while (k<number_of_dimensions_total) {
      entry[k]++;
      if (entry[k]>=tprod_index_sizes[k]) {
	entry[k]=0;
	k++;
	if (k>=number_of_dimensions_total) {
	  break; }
      } else {
	if (k>=number_of_dimensions_sums) {
	  output_k++; }
	break;
      }
    }
    if (k>=number_of_dimensions_total) {
      break ; }

    if (output_k>numel) {
      mexPrintf("Last entry [");
      for (k=0;k<number_of_dimensions_total;k++) {
	mexPrintf("i%d=%d,",k,entry[k]);
      }
      mexPrintf("] @ %d\n",output_k);
      mexErrMsgIdAndTxt( "MATLAB:Ctprod:invalidInput",
			 "Internal error: output has %d elemenys, but trying to write on element %d",numel,output_k);
    }
  }

  // mexPrintf("Freeing indices\n");
  for (i=0;i<nrhs-1;i+=2) {
    mxFree(indices[i/2]);
  }
  if (number_of_dimensions_total>0) {
    // mexPrintf("Freeing entry\n");
    mxFree(entry);
    // mexPrintf("Freeing tprod_index_sizes \n");
    mxFree(tprod_index_sizes);
  }


}
