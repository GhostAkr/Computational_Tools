//
// Created by ighos on 15.11.2018.
//

#ifndef LAB_1_INTERPOLATION_H
#define LAB_1_INTERPOLATION_H

#define Pi 3.1415

#include <iostream>
#include <string>
#include <fstream>

using std::cout;
using std::endl;

// Printing
void Print(double* M, int n);
void Print(double** M, int n);

// Meshes
double* Mesh(double a,double b, int n);
double* MeshCheb(double a, double b, int n);

// Functions
double** func1(double* M, int n);  // M - argument for function

// Interpolation methods
double** Polynom(double** S,double* R,int n, int m);  // n - size of mesh S, m - size of mesh R

// Other methods
void Extracttofile(double** P, int m, std::string _pathToFile);  // P - table of arguments and function values; m - dimension

#endif //LAB_1_INTERPOLATION_H
