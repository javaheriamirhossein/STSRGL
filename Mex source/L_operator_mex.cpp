#include "mex.h"
#include <cmath>
#include <vector>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check the number of inputs and outputs
    if (nrhs != 1) {
        mexErrMsgIdAndTxt("L_operator:invalidInput", "One input required (w).");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("L_operator:invalidOutput", "One output required (L).");
    }

    // Get the input vector 'w'
    double *w = mxGetPr(prhs[0]);
    int Nw = mxGetNumberOfElements(prhs[0]);

    // Compute the size 'N' of the Laplacian matrix
    int N = static_cast<int>((1 + sqrt(1 + 8 * Nw)) / 2);

    // Create the weight matrix 'W'
    mxArray *W = mxCreateDoubleMatrix(N, N, mxREAL);
    double *W_data = mxGetPr(W);

    // Fill the upper triangular part of 'W' with 'w'
    int cnt = 0;
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            W_data[i + j * N] = w[cnt];
            cnt++;
        }
    }

    // Make 'W' symmetric by copying the upper triangular part to the lower triangular part
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            W_data[j + i * N] = W_data[i + j * N];
        }
    }

    // Compute the Laplacian matrix 'L'
    mxArray *L = mxCreateDoubleMatrix(N, N, mxREAL);
    double *L_data = mxGetPr(L);

    // Compute the diagonal of 'L' as the row sums of 'W'
    for (int i = 0; i < N; i++) {
        double row_sum = 0.0;
        for (int j = 0; j < N; j++) {
            row_sum += W_data[i + j * N];
        }
        L_data[i + i * N] = row_sum;
    }

    // Subtract 'W' from the diagonal matrix
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            L_data[i + j * N] -= W_data[i + j * N];
        }
    }

    // Set the output
    plhs[0] = L;

    // Clean up
    mxDestroyArray(W);
}