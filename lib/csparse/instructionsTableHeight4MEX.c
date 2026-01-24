/* mex -largeArrayDims COPTIMFLAGS="-Ofast -DNDEBUG" CFLAGS="\$CFLAGS -Wall" instructionsTableHeight4MEX.c */

#include <dlfcn.h>
#include "mex.h"

#include "instructionsTable.h"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* outputs */
    int64_t *height;

    /* Process inputs */

    /* Check # inputs */
    if (nrhs != 0)
    {
        mexErrMsgIdAndTxt("instructionsTableHeight4MEX:nrhs", "0 inputs required, %d found.", nrhs);
        return;
    }

    /* Process outputs */

    /* Check # outputs */
    if (nlhs != 1)
    {
        mexErrMsgIdAndTxt("instructionsTableHeight4MEX:nrhs", "1 outputs required, %d found.", nlhs);
        return;
    }

    /* output height */
    {
        mwSize dims[] = {1, 1};
        plhs[0] = mxCreateNumericArray(2, dims, mxINT64_CLASS, mxREAL);
        height = mxGetData(plhs[0]);
    }

    /* Call function */
    void *libHandle = dlopen("instructionsTable.dylib", RTLD_NOW);
    if (!libHandle)
    {
        printf("[%s] Unable to open library: %s\n", __FILE__, dlerror());
        return;
    }
    void (*instructionsTableHeight4MEX)(/* outputs */
                                        int64_t *height) = dlsym(libHandle, "instructionsTableHeight4MEX");
    if (!instructionsTableHeight4MEX)
    {
        printf("[%s] Unable to get symbol: %s\n", __FILE__, dlerror());
        return;
    }
    instructionsTableHeight4MEX(height);
    dlclose(libHandle);
}
