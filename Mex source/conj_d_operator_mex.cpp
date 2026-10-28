#include "mex.h"
#include <vector>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check the number of inputs and outputs
    if (nrhs != 1) {
        mexErrMsgIdAndTxt("conj_d_operator_simple:invalidInput", "One input required (a).");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("conj_d_operator_simple:invalidOutput", "One output required (v).");
    }

    // Get the input vector 'a'
    double *a = mxGetPr(prhs[0]);
    int N = mxGetNumberOfElements(prhs[0]); // Number of elements in 'a'

    // Compute the number of elements in the output vector 'v'
    int Nw = N * (N - 1) / 2;

    // Create the output vector 'v'
    plhs[0] = mxCreateDoubleMatrix(Nw, 1, mxREAL);
    double *v = mxGetPr(plhs[0]);

    // Compute the output vector 'v'
    int k = 0;
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            if (i > j) {
                v[k] = a[i] + a[j];
                k++;
            }
        }
    }
}