//
// Created by ghostakr on 9/16/18.
//

#include "../include/LinearSystems.h"

Matrix* gaussLinearSolve(Matrix* _A, Matrix* _B) {
    Matrix* newA = new Matrix;
    Matrix* newB = new Matrix;
    size_t rowsA = _A->rowsGet();
    size_t rowsB = _B->rowsGet();
    size_t colsA = _A->colsGet();
    size_t colsB = _B->colsGet();
    newA->matrixNullSet(rowsA, colsA);
    newB->matrixNullSet(rowsB, colsB);
    size_t n = rowsB;
    cout << "n = " << n << endl;
    // TODO: fix algorithm
    for (int k = 1; k <= n - 1; ++k) {
        for (int i = k + 1; i <= n; ++i) {
            for (int j = k; j <= n; ++j) {
                double c = _A->matrixGet()[i][j - 1] / _A->matrixGet()[i - 1][j - 1];
                newA->matrixGet()[i][j] = _A->matrixGet()[i][j] - _A->matrixGet()[i - 1][j] * c;
                newB->matrixGet()[i][0] = _B->matrixGet()[i][0] - _B->matrixGet()[i - 1][0] * c;
            }
        }
    }
    //cout << "New matrix A is " << endl;
    //newA->matrixPrint();
    //cout << "New matrix B is " << endl;
    //newB->matrixPrint();
}