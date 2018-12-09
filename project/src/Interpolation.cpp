//
// Created by ighos on 15.11.2018.
//

#include <cmath>
#include "../include/Interpolation.h"
#include "../include/Matrix.h"
#include "../include/LinearSystems.h"

// Printing

void Print(double* M, int n) {
    for (int i = 0; i < n; i++) {
        cout << M[i] << endl;
    }
}

void Print(double** M, int n) {
    cout << "X" << "       " << "Y" << endl;
    for (int i = 0; i < n; i++) {
        cout << M[0][i] <<"   "<< M[1][i] << endl;
    }
}

// Meshes

double* Mesh(double a, double b, int n) {
    double* M = new double[n + 1];
    double h = (b - a) / n;
    for (int i = 0; i < n + 1; i++) {
        M[i] = a + h * i;
    }
    return(M);
}

double* MeshCheb(double a, double b, int n) {
    n++;
    auto* result = new double [n];
    double c1 = (a + b) / 2;
    double c2 = (b - a) / 2;
    double denom = 2 * (n + 1);
    for (int i = 0; i < n; ++i) {
        result[i] = c1 + c2 * cos((2 * i + 1) * Pi / denom);
    }
    // Reversing of mesh
    auto* tmp = new double [n];
    for (int i = 0; i < n; ++i) {
        tmp[i] = result[n - i - 1];
    }
    return tmp;
}

double* MeshChebX(double* Mesh, int n, int nOfParts) {
    n++;
    int m = n * nOfParts - 1;
    auto* result = new double [m];
    for (int i = 1; i < n; ++i) {
        double h = (Mesh[i] - Mesh[i - 1]) / nOfParts;
        for (int k = 0; k < nOfParts; ++k) {
            result[(i - 1) * nOfParts + k] = Mesh[i - 1] + h * k;
        }
    }
    result[m - 1] = Mesh[n - 1];
    return result;
}

// Functions

double** func1(double* M, int n) {
    double** S = new double*[2];
    S[0] = M;
    S[1] = new double[n];
    for (int i = 0; i < n; i++) {
        S[1][i] = S[0][i] * S[0][i];
    }
    return(S);
}

double** func2(double* M, int n) {
    double** S = new double*[2];
    S[0] = M;
    S[1] = new double[n];
    for (int i = 0; i < n; i++) {
        S[1][i] = 1 / (1 + 25 * S[0][i] * S[0][i]);
    }
    return(S);
}

double** func3(double* M, int n) {
    double** S = new double*[2];
    S[0] = M;
    S[1] = new double[n];
    for (int i = 0; i < n; i++) {
        S[1][i] = 1 / atan(1 + 10 * S[0][i] * S[0][i]);
    }
    return(S);
}

double** func4(double* M, int n) {
    double** S = new double*[2];
    S[0] = M;
    S[1] = new double[n];
    for (int i = 0; i < n; i++) {
        double x = S[0][i];
        double x2 = x * x;
        double x3 = x2 * x;
        S[1][i] = pow(4.0 * x3 + 2.0 * x2 - 4.0 * x + 2.0, sqrt(2.0)) + asin(1.0 / (5.0 + x - x2)) - 5.0;
    }
    return(S);
}

double** func5(double* M, int n) {
    double** S = new double*[2];
    S[0] = M;
    S[1] = new double[n];
    for (int i = 0; i < n; i++) {
        S[1][i] = exp(S[0][i]);
    }
    return(S);
}

double** funcexp(double* M, int n) {
    double** S = new double*[2];
    S[0] = M;
    S[1] = new double[n];
    for (int i = 0; i < n; i++) {
        S[1][i] = exp(S[0][i]);
    }
    return(S);
}

double** func6(double* M, int n) {
    double** S = new double*[2];
    S[0] = M;
    S[1] = new double[n];
    for (int i = 0; i < n; i++) {
        S[1][i] = S[0][i] * S[0][i] * S[0][i];
    }
    return(S);
}

// Interpolation methods

double** Polynom(double** S, double* R, int n, int m) {
    double** P = new double* [2];
    P[0] = R;
    P[1] = new double [m];
    double c;
    double sum;
    for (int k = 0; k < m; k++) {
        sum = 0.0;
        for (int i = 0; i < n; i++) {
            c = 1.0;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    c *= (P[0][k] - S[0][j]) / (S[0][i] - S[0][j]);
                }
            }
            sum += c * S[1][i];
        }
        P[1][k] = sum;
    }
    return P;
}

double** Spline(double** S, double* R, int n, int m) {
    n++;
    auto** P = new double* [2];
    P[0] = R;
    P[1] = new double [m];
    Matrix* triDiagC = new Matrix;
    triDiagC->matrixNullSet(n - 2, 4);
    auto* A = new double [n - 1];
    for (int i = 0; i < n - 1; ++i) {
        A[i] = S[1][i];
    }
    auto* H = new double [n - 1];
    auto* G = new double [n - 1];
    for (int i = 1; i < n; ++i) {
        H[i - 1] = S[0][i] - S[0][i - 1];
        G[i - 1] = (S[1][i] - S[1][i - 1]) / H[i - 1];
    }
    for (int i = 0; i < n - 2; ++i) {
        if (i == 0) {
            triDiagC->matrixGet()[i][0] = 0.0;
        } else {
            triDiagC->matrixGet()[i][0] = H[i];
        }
        triDiagC->matrixGet()[i][1] = 2.0 * (H[i] + H[i + 1]);
        if (i == n - 3) {
            triDiagC->matrixGet()[i][2] = 0.0;
        } else {
            triDiagC->matrixGet()[i][2] = H[i + 1];
        }
        triDiagC->matrixGet()[i][3] = 3.0 * (G[i + 1] - G[i]);
    }
    Matrix* CSolve = tridiagonalLinearSolve(triDiagC);
    auto* C = new double [n - 1];  // Real coefficients for spline
    for (int i = 0; i < n - 1; ++i) {
        if (i == 0) {
            C[i] = 0.0;
        } else {
            C[i] = CSolve->matrixGet()[i - 1][0];
        }
    }
    auto* B = new double [n - 1];
    for (int i = 0; i < n - 1; ++i) {
        B[i] = G[i] - (C[i + 1] + 2.0 * C[i]) * H[i] / 3.0;
    }
    auto* D = new double [n - 1];
    for (int i = 0; i < n - 1; ++i) {
        D[i] = (C[i + 1] - C[i]) / (3.0 * H[i]);
    }
    // Building of spline
    int j = 0;  // Counter
    for (int k = 0; k < m; ++k) {
        if (P[0][k] > S[0][j + 1]) {
            j++;
        }
        double delta = P[0][k] - S[0][j];
        double delta2 = delta * delta;
        double delta3 = delta2 * delta;
        P[1][k] = A[j] + B[j] * delta + C[j] * delta2 + D[j] * delta3;
    }
    return P;
}

// Other methods

void Extracttofile(double** P, int m, std::string _pathToFile) {
    std::ofstream fileOut(_pathToFile);
    if (!fileOut) {  // Exception
        cout << "Error while reading file" << endl;
        return;
    }
    for (int i = 0; i < m; i++) {
        fileOut << P[0][i] << " " << P[1][i] << endl;
    }
    fileOut.close();
}

double max(double* M, int n) {
    double res = fabs(M[0]);
    for (int i = 1; i < n; i++) {
        if (fabs(M[i]) > res) {
            res = fabs(M[i]);
        }
    }
    return res;
}

double pogr(double** P, int n) {
    double* pogr = new double[n];
    double** justf = funcexp(P[0], n);
    for (int i = 0; i < n; i++) {
        pogr[i] = fabs(P[1][i] - justf[1][i]);
    }
    return max(pogr, n);
}
