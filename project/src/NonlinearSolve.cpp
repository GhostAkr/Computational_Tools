//
// Created by ighos on 08.01.2019.
//

#include "../include/NonlinearSolve.h"

double** f1(double* x, int n) {
    double** S = new double*[2];
    S[0] = x;
    S[1] = new double[n];
    for (int i = 0; i < n; i++) {
        S[1][i] = (S[0][i] - 0.1) * (S[0][i] - 0.22) * (S[0][i] - 0.55) * (S[0][i] - 0.7) * (S[0][i] - 0.75);
    }
    return S;
}

double** f2(double* x, int n) {
    double** S = new double*[2];
    S[0] = x;
    S[1] = new double[n];
    for (int i = 0; i < n; i++) {
        S[1][i] = sqrt(S[0][i] + 1) - 1;
    }
    return S;
}

double** f3(double* x, int n) {
    double** S = new double*[2];
    S[0] = x;
    S[1] = new double[n];
    for (int i = 0; i < n; i++) {
        double x2 = S[0][i] * S[0][i];
        double x3 = x2 * S[0][i];
        S[1][i] = 35 * x3 - 67 * x2 - 3 * S[0][i] + 3;
    }
    return S;
}

double** rootsLocale(double** _funcMesh, int n, int* nOfPairs) {
    double eps = 1e-14;
    auto** tmp = new double* [2];  // Temporary array to save correct points and values
    for (int i = 0; i < 2; ++i) {
        tmp[i] = new double [n];
    }
    int numberOfPairs = 0;  // Real number of elements in tmp
    for (int i = 1; i < n; ++i) {
        if ((_funcMesh[1][i] * _funcMesh[1][i - 1] < 0) || (fabs(_funcMesh[1][i]) < eps) || (fabs(_funcMesh[1][i - 1]) < eps)) {  // If signs are different
            numberOfPairs++;
            tmp[0][numberOfPairs - 1] = _funcMesh[0][i - 1];
            tmp[1][numberOfPairs - 1] = _funcMesh[1][i - 1];
            tmp[0][numberOfPairs] = _funcMesh[0][i];
            tmp[1][numberOfPairs] = _funcMesh[1][i];
            numberOfPairs++;
        }
    }
    auto** res = new double* [2];  // Array with correct points with size = numberOfPairs (real size)
    for (int i = 0; i < 2; ++i) {  // Copying elements to result array
        res[i] = new double [numberOfPairs];
        for (int j = 0; j < numberOfPairs; ++j) {
            res[i][j] = tmp[i][j];
        }
    }
    *nOfPairs = numberOfPairs;  // Saving size of result array
    for (int i = 0; i < 2; ++i) {
        delete [] tmp[i];
    }
    delete [] tmp;
    return res;
}