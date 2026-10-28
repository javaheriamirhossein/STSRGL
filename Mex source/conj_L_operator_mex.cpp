#include "mex.h"
#include <vector>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check the number of inputs and outputs
    if (nrhs != 1) {
        mexErrMsgIdAndTxt("conj_L_operator_simple:invalidInput", "One input required (A).");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("conj_L_operator_simple:invalidOutput", "One output required (v).");
    }

    // Get the input matrix 'A'
    double *A = mxGetPr(prhs[0]);
    int N = mxGetM(prhs[0]); // Number of rows in A

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
                v[k] = A[i + i * N] + A[j + j * N] - A[i + j * N] - A[j + i * N];
                k++;
            }
        }
    }
}