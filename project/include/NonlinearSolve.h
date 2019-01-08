//
// Created by ighos on 08.01.2019.
//

#ifndef LAB_1_NONLINEARSOLVE_H
#define LAB_1_NONLINEARSOLVE_H

#include <iostream>
#include <cmath>
#include "../include/Interpolation.h"

using std::cout;
using std::endl;

// TODO: finish localization of roots
double** rootsLocale(double* _funcMesh, int n, int* nOfPairs);  // n - size of mesh, nOfPairs - number of points in pairs

// TODO: write rest tests
// Functions
double f1(double x);
double f2(double x);
double f3(double x);

#endif //LAB_1_NONLINEARSOLVE_H
