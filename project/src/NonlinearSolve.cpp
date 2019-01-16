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

double ff1(double x) {
    return (x - 0.1) * (x - 0.22) * (x - 0.55) * (x - 0.7) * (x - 0.75);
}

double ff2(double x) {
    return sqrt(x + 1) - 1;
}

double ff3(double x) {
    double x2 = x * x;
    double x3 = x2 * x;
    return 35 * x3 - 67 * x2 - 3 * x + 3;
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

double* Bisection(double** _intervals, int nOfPairs, int* _nOfRoots, double f(double)) {
    auto* result = new double [nOfPairs];
    double nullEps = 1e-14;
    double eps = 1e-5;
    *_nOfRoots = 0;
    for (int i = 0; i < nOfPairs; i += 2) {
        double a = _intervals[0][i];
        double b = _intervals[0][i + 1];
        if (fabs(_intervals[1][i]) < nullEps) {
            addRoot(result, a, _nOfRoots);
            continue;
        }
        if (fabs(_intervals[1][i + 1]) < nullEps) {
            addRoot(result, b, _nOfRoots);
            continue;
        }
        double h = b - a;
        double x = (a + b) / 2.0;
        while (fabs(h) > eps) {
            x = (a + b) / 2.0;
            double y = f(x);
            if ((f(a) * y) <= 0) {
                b = x;
            } else {
                a = x;
            }
            h = b - a;
        }
        addRoot(result, x, _nOfRoots);
    }
    return result;
}

void addRoot(double* _targetArray, double _root, int* _arraySize) {
    double eps = 1e-14;
    // Searching if root is in the _targetArray already
    for (int i = 0; i < *_arraySize; ++i) {
        if (fabs(_targetArray[i] - _root) < eps) {
            return;
        }
    }
    _targetArray[(*_arraySize)++] = _root;
}

double fd1(double x) {
    return 5 * (0.024299 + x * (-0.30238 + x * (1.1907 + (-1.856 + x)*x)));
}

double fd2(double x) {
    return 1 / (2 * (sqrt(1 + x)));
}

double fd3(double x) {
    return 105 * x*x - 134 * x - 3;
}

double* Newton(double** _intervals, int nOfPairs, int* _nOfRoots, double f(double)) {
    auto* result = new double[nOfPairs];
    double nullEps = 1e-14;
    double eps = 1e-5;
    double a, b, x1, x2, fa, fb, fd;
    x1 = 0;
    *_nOfRoots = 0;
    for (int i = 0; i < nOfPairs; i += 2) {
        a = _intervals[0][i];
        b = _intervals[0][i + 1];
        if ((fabs(_intervals[1][i]) < nullEps)||(fabs(_intervals[1][i + 1]) < nullEps)) {
            if (fabs(_intervals[1][i]) < nullEps) {
                addRoot(result, a, _nOfRoots);
            }
            if (fabs(_intervals[1][i + 1]) < nullEps) {
                addRoot(result, b, _nOfRoots);

            }
            continue;
        }
        fa = f(a);
        fb = f(b);
        x2 = (fa*b - fb * a) / (fa - fb);
        while (fabs(x1 - x2) >= eps) {
            x1 = x2;
            fd = (f(x1 + eps) - (x1)) / eps;
            x2 = x1 - f(x1) / fd;

        }
        addRoot(result, x2, _nOfRoots);
    }
    return result;
}

//double** Newtonsys(double n) {
//    double** res = new double*[2];
//    double** mesh = new double*[2];
//    mesh[0] = new double[n];
//    mesh[1] = new double[n];
//    double h = 20 / (n - 1);
//    for (int i = 0; i < n; i++) {
//        mesh[0][i] = -10 + h * i;
//        mesh[1][i] = -10 + h * i;
//    }
//}