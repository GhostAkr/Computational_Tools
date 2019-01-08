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

Matrix* HessenbergForm(const Matrix* _A);
Matrix* QRDecompositionEigen(Matrix* _A);
Matrix* Reverse(Matrix* _A, double eigen);
double Rayleigh(Matrix* _A);

void shiftMatrix(Matrix* _A, double _shift);  // _A +=_shift * IdentityMatrix
double cosinus(Matrix* vec1, Matrix* vec2);  // Cosinus between vectors matr1 and matr2
double scalarProd(Matrix* vec1, Matrix* vec2);  // Scalar product of vec1 and vec2

#endif //LAB_1_EIGENVALUES_H
