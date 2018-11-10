//
// Created by ghostakr on 9/16/18.
//

#ifndef LAB_1_LINEARSYSTEMS_H
#define LAB_1_LINEARSYSTEMS_H

#include "../include/Matrix.h"

// Linear system solving (solves equations A x = B)

Matrix* gaussLinearSolve(Matrix* _A);
Matrix* QRDecompositionSolve(Matrix* _A, Matrix* Q, Matrix* R);
Matrix* QRBackTurn(Matrix* _Q, Matrix* _R, Matrix* _b);
void pertrubationSolution(Matrix* _A);
Matrix* fixedPointIterationSolve(Matrix* _A);
Matrix* Jacobi(const Matrix* _matrix);
Matrix* SOR(const Matrix* _matrix);
Matrix* Seidel(const Matrix* _matrix);

// Other methods

type valuationVector(Matrix* _solution, Matrix* _system);  // Calculating of valuation vector
bool onlyDesitionCheck(Matrix* _matrix);
Matrix* rotationMatrix(Matrix* _matrix, int line, int numOfVar);  // Calculating of rotation matrix for 'nuOfVar' element on 'line' line
Matrix* rotateMatrix(Matrix* _targetMatrix, int _line, int _numOfVar, type _c, type _s);  // Creates rotation matrix based on current element and rotates _targetMatrix
type vectorNorm(type* _vector, size_t _rows);
void conditionNumber(Matrix* _A);
type normInf(Matrix* _A);
type normOne(Matrix* _A);
type normInfVect(Matrix* _A);
type normOneVect(Matrix* _A);
Matrix* createTridiagonalMatrix(int variant);
Matrix* C_SOR(Matrix* _matrix, type omega);  // Creating of C-Matrix for SOR method with linear system _matrixMatrix* Seidel(const Matrix* _matrix);
Matrix* CSeidel(const Matrix* _matrix);
Matrix* inverseMatrix(const Matrix* _matrix);  // TODO: Move this to Matrix class
type pogrNorm(Matrix* _solution, Matrix* _realSolution);
#endif //LAB_1_LINEARSYSTEMS_H
