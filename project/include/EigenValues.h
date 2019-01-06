//
// Created by ighos on 18.12.2018.
//

#ifndef LAB_1_EIGENVALUES_H
#define LAB_1_EIGENVALUES_H

#include <iostream>
#include <cmath>
#include "../include/Matrix.h"
#include "../include/LinearSystems.h"

using std::cout;
using std::endl;

Matrix* HessenbergForm(Matrix* _A);
Matrix* QRDecompositionEigen(Matrix* _A);
Matrix* Reverse(Matrix* _A);
double Rayleigh(Matrix* _A);
#endif //LAB_1_EIGENVALUES_H
