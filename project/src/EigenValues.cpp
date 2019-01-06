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

Matrix* Reverse(Matrix* _A) {
    double eigen = 1;
    int rows = _A->rowsGet();
    Matrix* res = new Matrix;
    Matrix* Q = new Matrix;
    Matrix* R = new Matrix;
    Matrix* t = new Matrix;
    Matrix* M = new Matrix;
    t->matrixNullSet(rows, 1);
    res->matrixNullSet(rows,1);
    M->matrixNullSet(rows, rows);
    //for (int i = 0; i < rows; i++) {
    //	res->matrixGet()[i][0] = 1;
    //}
    res->matrixGet()[0][0] = 1;

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < rows; j++) {
            M->matrixGet()[i][j] = _A->matrixGet()[i][j];
            if (i == j) {
                M->matrixGet()[i][j] -= eigen;
            }
        }
    }
    QRDecompositionSolve(M,Q,R);

    //cout << endl;
    int n = 0;
    while (n < 6) {
        //Q->matrixPrint();
        //for (int i = 0; i < rows; i++) {
        //	M->matrixGet()[i][rows] = res->matrixGet()[i][0];
        //}

        //temp = Matrix::getCopy(M);
        //temp->matrixPrint();
        //res = gaussLinearSolve(temp);
        //res->matrixPrint();
        //cout << endl;
        t = QRBackTurn(Q,R,res);
        //t->matrixPrint();
        Q->matrixTranspose();
        //t->matrixPrint();
        //cout << endl;
        //res->matrixPrint();
        //cout << endl;
        res = Matrix::getCopy(t);
        //res->matrixPrint();
        //cout << endl;
        double div = norm(res);
        for (int i = 0; i < rows; i++) {
            res->matrixGet()[i][0] /= div;
        }
        //res->matrixPrint();
        //cout << endl;
        n++;
    }


    //Matrix::matrixComp(_A, res)->matrixPrint();
    //cout << endl;
    return res;
}

double Rayleigh(Matrix* _A) {
    double eigen = 0;
    int rows = _A->rowsGet();
    Matrix* res = new Matrix;
    Matrix* temp = new Matrix;
    Matrix* M1 = new Matrix;
    Matrix* M = new Matrix;
    M->matrixNullSet(rows,rows + 1);
    //for (int i = 0; i < rows; i++) {
    //	for (int j = 0; j < rows; j++) {
    //		M->matrixGet()[i][j] = _A->matrixGet()[i][j];
    //	}
    //}
    M = Matrix::getCopy(_A);
    temp->matrixNullSet(rows, rows + 1);
    res->matrixNullSet(rows, 1);
    int n = 0;

    res->matrixGet()[0][0] = -0.87;
    res->matrixGet()[1][0] = 0;
    res->matrixGet()[2][0] = -0.25;
    res->matrixGet()[3][0] = -0.43;


    while(n<3) {
        M = Matrix::getCopy(_A);
        eigen =0;
        M1 = Matrix::matrixComp(_A, res);
        //cout << endl;
        for (int i = 0; i < rows; i++) {
            eigen += M1->matrixGet()[i][0] * res->matrixGet()[i][0];
        }
        for (int i = 0; i < rows; i++) {
            M->matrixGet()[i][i] -= eigen;
        }
        for (int i = 0; i < rows; i++) {
            M->matrixGet()[i][rows] = res->matrixGet()[i][0];
        }

        temp = Matrix::getCopy(M);
        //temp->matrixPrint();
        res = gaussLinearSolve(temp);
        //res->matrixPrint();
        //cout << endl << eigen << endl;
        double div = norm(res);
        for (int i = 0; i < rows; i++) {
            res->matrixGet()[i][0] /= div;
        }
        n++;
    }
    return eigen;
}
