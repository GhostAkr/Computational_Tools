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

double** rootsLocale(double** _funcMesh, int n, int* nOfPairs);  // n - size of mesh, nOfPairs - number of points in pairs
double* Bisection(double** _intervals, int nOfPairs, int* _nOfRoots, double f(double));
void addRoot(double* _targetArray, double _root, int* _arraySize, int _iteration);
void addSysRoot(double** _targetArray, double* _root, int* _arraySize);
double* Newton(double** _intervals, int nOfPairs, int* _nOfRoots, double f(double));
double** Newtonsys(double n, std::string path, int* nroot);

// TODO: write rest tests
// Functions
double** f1(double* x, int n);
double** f2(double* x, int n);
double** f3(double* x, int n);

double ff1(double x);
double ff2(double x);
double ff3(double x);

double* fsys1(double* x);
double* fsys2(double* x);
double** f1jac(double* x);
double** f2jac(double* x);

#endif //LAB_1_NONLINEARSOLVE_H
