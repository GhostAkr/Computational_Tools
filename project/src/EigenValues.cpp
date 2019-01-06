//
// Created by ighos on 18.12.2018.
//

#include "../include/EigenValues.h"

Matrix* HessenbergForm(const Matrix* _A) {
    int n = _A->rowsGet();
    Matrix* Result = Matrix::getCopy(_A);
    for (int k = 0; k < n - 2; ++k) {
        for (int l = k + 2; l < n; ++l) {
            double sqr = sqrt(Result->matrixGet()[k + 1][k] * Result->matrixGet()[k + 1][k] + \
                    Result->matrixGet()[l][k] *Result->matrixGet()[l][k]);
            double alpha = Result->matrixGet()[k + 1][k] / sqr;
            double beta = Result->matrixGet()[l][k] / sqr;
            Matrix* tmp1 = Matrix::getCopy(Result);
            for (int i = 0; i < n; ++i) {
                Result->matrixGet()[k + 1][i] = alpha * tmp1->matrixGet()[k + 1][i] + beta * tmp1->matrixGet()[l][i];
                Result->matrixGet()[l][i] = -beta * tmp1->matrixGet()[k + 1][i] + alpha * tmp1->matrixGet()[l][i];
            }
            delete tmp1;
            Matrix* tmp2 = Matrix::getCopy(Result);
            for (int i = 0; i < n; ++i) {
                Result->matrixGet()[i][k + 1] = alpha * tmp2->matrixGet()[i][k + 1] + beta * tmp2->matrixGet()[i][l];
                Result->matrixGet()[i][l] = -beta * tmp2->matrixGet()[i][k + 1] + alpha * tmp2->matrixGet()[i][l];
            }
            delete tmp2;
        }
    }
    return Result;
}

Matrix* QRDecompositionEigen(Matrix* _A) {
    double epsilon = 1e-14;  // For comparing with null
    int rows = _A->rowsGet();
    Matrix* Result = new Matrix;
    Result->matrixNullSet(rows, 1);
    Matrix* CurA = Matrix::getCopy(_A);
    int iteration[rows];
    for (int i = 0; i < rows; ++i) {
        iteration[i] = 0;
        Matrix* Buff = Matrix::getCopy(CurA);  // Buffer matrix for creation of new A
        CurA = new Matrix;
        CurA->matrixNullSet(rows - i, rows - i);
        for (int k = 0; k < rows - i; ++k) {  // Creation of new A
            for (int l = 0; l < rows - i; ++l) {
                CurA->matrixGet()[k][l] = Buff->matrixGet()[k][l];
            }
        }
        if (i == rows - 1) {  // Last eigenvalue
            Result->matrixGet()[i][0] = CurA->matrixGet()[0][0];
            break;
        }
        while(true) {  // Cycle for each eigenvalue
            iteration[i]++;
            double shift = CurA->matrixGet()[rows - 1 - i][rows - 1 - i];
            bool IsNextIteration = false;
            shiftMatrix(CurA, -shift);  // Shifting
            Matrix* Q = new Matrix;
            Matrix* R = new Matrix;
            Matrix* copyA = Matrix::getCopy(CurA);  // Copy of A for QR-Decomposition
            QRDecompositionSolve(copyA, Q, R);
            CurA = Matrix::matrixComp(R, Q);
            shiftMatrix(CurA, shift);  // Reverse shifting
            for (int k = 0; k < rows - i - 1; ++k) {  // Checking for break condition
               if (fabs(CurA->matrixGet()[rows - i - 1][k]) > epsilon) {
                   IsNextIteration = true;
               }
            }
            if (!IsNextIteration) {  // Breakpoint
                Result->matrixGet()[i][0] = CurA->matrixGet()[rows - i - 1][rows - i - 1];  // Last eigenvalue
                delete copyA;
                delete Q;
                delete R;
                break;
            }
            delete copyA;
            delete Q;
            delete R;
        }
    }
    return Result;
}

void shiftMatrix(Matrix* _A, double _shift) {
    int rows = _A->rowsGet();
    int cols = _A->colsGet();
    if (rows != cols) {  // Exception
        cout << "Number of rows and number of cols are not equal" << endl;
        return;
    }
    for (int i = 0; i < rows; ++i) {
        _A->matrixGet()[i][i] += _shift;
    }
}
