//
// Created by ghostakr on 9/16/18.
//

#include "../include/LinearSystems.h"

Matrix* gaussLinearSolve(Matrix* _A) {
    double eps = 10e-5;  // For comparing with 0.0
    size_t rowsA = _A->rowsGet();
    size_t colsA = _A->colsGet();
    // First part
    for (int k = 0; k < colsA - 1; ++k) {
        int mainElem = _A->mainElement(k);
        _A->matrixRowsChange(mainElem, k);
        for (int i = k + 1; i < rowsA; ++i) {
            double coeff = _A->matrixGet()[i][k] / _A->matrixGet()[k][k];
            for (int j = k; j < colsA; ++j) {
                if (fabs(_A->matrixGet()[i][j] -= _A->matrixGet()[k][j] * coeff) < eps) {  // Applying accuracy
                    _A->matrixGet()[i][j] = 0.0;
                }
            }
        }
    }
    // Second part
    Matrix* result = new Matrix;
    result->matrixNullSet(rowsA, 1);
    for (int i = rowsA - 1; i >= 0; --i) {
        double leftSum = 0.0;
        for (int j = rowsA - 1; j >= i; --j) {
            leftSum += _A->matrixGet()[i][j] * result->matrixGet()[j][0];
        }
        result->matrixGet()[i][0] = (_A->matrixGet()[i][colsA - 1] - leftSum) / _A->matrixGet()[i][i];
        if (fabs(result->matrixGet()[i][0]) < eps) {
            result->matrixGet()[i][0] = 0.0;
        }
    }
    return result;
}