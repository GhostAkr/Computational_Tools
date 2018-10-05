//
// Created by ghostakr on 9/16/18.
//

#ifndef LAB_1_LINEARSYSTEMS_H
#define LAB_1_LINEARSYSTEMS_H

#include "../include/Matrix.h"

#define EPS1 1e-7
#define EPS2 1e-14

// Linear system solving (solves equations A x = B)

Matrix* gaussLinearSolve(Matrix* _A);
Matrix* QRDecompositionSolve(Matrix* _A);

// Other methods

int valuationVector(Matrix* _solution, Matrix* _system);  // Calculating of valuation vector
bool onlyDesitionCheck(Matrix* _matrix);
Matrix* rotationMatrix(Matrix* _matrix, int line, int numOfVar);  // Calculating of rotation matrix for 'nuOfVar' element on 'line' line
double vectorNorm(double* _vector, size_t _rows);
double conditionNumber(Matrix* _A);
double normInf(Matrix* _A);
double normOne(Matrix* _A);

#endif //LAB_1_LINEARSYSTEMS_H
