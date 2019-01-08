//
// Created by ighos on 08.01.2019.
//

#include "../include/NonlinearSolve.h"

double f1(double x) {
    return (x - 0.1) * (x - 0.22) * (x - 0.55) * (x - 0.7) * (x - 0.75);
}

double f2(double x) {
    return sqrt(x + 1) - 1;
}

double f3(double x) {
    double x2 = x * x;
    double x3 = x2 * x;
    return 35 * x3 - 67 * x2 - 3 * x + 3;
}

double** rootsLocale(double** _funcMesh, int n, int* nOfPairs) {
    auto** tmp = new double* [2];
    for (int i = 0; i < 2; ++i) {
        tmp[i] = new double [n];
    }
    int numberOfPairs = 0;
    for (int i = 1; i < n; ++i) {
        if (_funcMesh[1][i] * _funcMesh[1][i - 1] < 0) {
            numberOfPairs++;
            tmp[0][numberOfPairs - 1] = _funcMesh[0][i - 1];
            tmp[1][numberOfPairs - 1] = _funcMesh[1][i - 1];
            tmp[0][numberOfPairs] = _funcMesh[0][i];
            tmp[1][numberOfPairs] = _funcMesh[1][i];
            numberOfPairs++;
        }
    }
    auto** res = new double* [2];
    for (int i = 0; i < 2; ++i) {
        res[i] = new double [numberOfPairs];
        for (int j = 0; j < numberOfPairs; ++j) {
            res[i][j] = tmp[i][j];
        }
    }
    return res;
}