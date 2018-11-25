//
// Created by ighos on 15.11.2018.
//

#include <cmath>
#include "../include/Interpolation.h"

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
    double* M = new double[n];
    double h = (b - a) / n;
    for (int i = 0; i < n; i++) {
        M[i] = a + h * i;
    }
    return(M);
}

double* MeshCheb(double a, double b, int n) {
    auto* result = new double [n];
    double c1 = (a + b) / 2;
    double c2 = (b - a) / 2;
    double denom = 2 * (n + 1);
    for (int i = 0; i < n; ++i) {
        result[i] = c1 + c2 * cos((2 * i + 1) * Pi / denom);
    }
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

// Interpolation methods

double** Polynom(double** S, double* R, int n, int m) {
    double** P = new double*[2];
    P[0] = R;
    P[1] = new double[m];
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
