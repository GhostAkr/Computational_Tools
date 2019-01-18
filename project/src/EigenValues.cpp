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

Matrix* Reverse(Matrix* _A, double eigen) {
    double eps = 1e-3;
    int rows = _A->rowsGet();
    Matrix* res = new Matrix;
    Matrix* Q = new Matrix;
    Matrix* R = new Matrix;
    Matrix* M = new Matrix;
    res->matrixNullSet(rows,1);
    M->matrixNullSet(rows, rows);
    for (int i = 0; i < rows; i++) {
    	res->matrixGet()[i][0] = 1.0;
    }
    double div = norm(res);
    for (int i = 0; i < rows; i++) {
        res->matrixGet()[i][0] /= div;
    }
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < rows; j++) {
            M->matrixGet()[i][j] = _A->matrixGet()[i][j];
            if (i == j) {
                M->matrixGet()[i][j] -= eigen;
            }
        }
    }
    QRDecompositionSolve(M,Q,R);
    int n = 0;
    Matrix* prev = new Matrix;
    do {
        n++;
        prev = Matrix::getCopy(res);
        Matrix* t = QRBackTurn(Q,R,res);
        res = Matrix::getCopy(t);
        Q->matrixTranspose();
        double div = norm(res);
        for (int i = 0; i < rows; i++) {
            res->matrixGet()[i][0] /= div;
        }
    } while (1.0 - fabs(cosinus(res, prev)) >= eps);
    delete prev;
    cout << "Number of iterations is " << n << endl;
    return res;
}

double Rayleigh(Matrix* _A, Matrix* res) {
    double eps = 1e-14;
    double eigen = 0;
    int rows = _A->rowsGet();
    Matrix* temp = new Matrix;
    Matrix* M1 = new Matrix;
    Matrix* M = new Matrix;
    M->matrixNullSet(rows,rows + 1);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < rows; j++) {
            M->matrixGet()[i][j] = _A->matrixGet()[i][j];
        }
    }
    temp->matrixNullSet(rows, rows + 1);
    int n = 0;
    double div = norm(res);
    for (int i = 0; i < rows; i++) {
        res->matrixGet()[i][0] /= div;
    }

    double prevEigen = 0.0;

    do {
        n++;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < rows; ++j) {
                M->matrixGet()[i][j] = _A->matrixGet()[i][j];
            }
        }
        prevEigen = eigen;
        eigen = 0;
        M1 = Matrix::matrixComp(_A, res);
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
        Matrix* Q = new Matrix;
        Matrix* R = new Matrix;
        QRDecompositionSolve(temp, Q, R);
        res = QRBackTurn(Q, R, res);
        double div = norm(res);
        for (int i = 0; i < rows; i++) {
            res->matrixGet()[i][0] /= div;
        }
    } while (fabs(prevEigen - eigen) > eps);
    M1 = Matrix::matrixComp(_A, res);
    double div1 = norm(M1);
    for (int i = 0; i < rows; i++) {
        M1->matrixGet()[i][0] /= div1;
    }
    cout << "Eigenvector for eigenvalue = " << eigen << " is (" << n << " iterations)" << endl;
    M1->matrixPrint();
    return eigen;
}

double cosinus(Matrix* vec1, Matrix* vec2) {
    if ((vec1->rowsGet() != vec2->rowsGet()) || (vec1->colsGet() != 1) || (vec2->colsGet() != 1)) {
        cout << "Vectors are incorrect" << endl;
        return -1.0;
    }
    double res = 1.0;
    res *= scalarProd(vec1, vec2);
    res /= norm(vec1);
    res /= norm(vec2);
    return res;
}

double scalarProd(Matrix* vec1, Matrix* vec2) {
    if ((vec1->rowsGet() != vec2->rowsGet()) || (vec1->colsGet() != 1) || (vec2->colsGet() != 1)) {
        cout << "Vectors are incorrect" << endl;
        return -1.0;
    }
    int rows = vec1->rowsGet();
    double sum = 0.0;
    for (int i = 0; i < rows; ++i) {
        sum += vec1->matrixGet()[i][0] * vec2->matrixGet()[i][0];
    }
    return sum;
}

bool eigenvalueCheck(Matrix* _A, Matrix* _eigenVals) {
    if (_A->rowsGet() != _A->colsGet()) {  // Exception
        cout << "Matrix A is incorrect" << endl;
        return false;
    }
    int matrixDimension = _A->rowsGet();
    int nOfVals = _eigenVals->rowsGet();
    auto* E = new Matrix;
    E->matrixOneSet(matrixDimension, matrixDimension);
    for (int i = 0; i < nOfVals; ++i) {
        Matrix* B = Matrix::matrixDiff(_A, Matrix::matrixConstComp(E, _eigenVals->matrixGet()[i][0]));
        if (gaussLinearSolve(B)) {
            cout << "Eigenvalue = " << _eigenVals->matrixGet()[i][0] << " is NOT CORRECT" << endl;
        } else {
            cout << "Eigenvalue = " << _eigenVals->matrixGet()[i][0] << " is CORRECT" << endl;
        }
        delete B;
    }
    delete E;
    return true;
}

void eigenvectorCheck(Matrix* _A, Matrix* _eigenVector, double _eigenValue) {
    if ((_A->rowsGet() != _A->colsGet()) || (_A->rowsGet() != _eigenVector->rowsGet())) {  // Exception
        cout << "Matrix A or eigenvector are incorrect" << endl;
        return;
    }
    int matrixDimension = _A->rowsGet();
    Matrix* b = Matrix::matrixConstComp(_eigenVector, _eigenValue);
    Matrix* B = new Matrix;
    B->matrixNullSet(matrixDimension, matrixDimension + 1);
    for (int i = 0; i < matrixDimension; ++i) {
        for (int j = 0; j < matrixDimension + 1; ++j) {
            if (j == matrixDimension) {
                B->matrixGet()[i][j] = b->matrixGet()[i][0];
            } else {
                B->matrixGet()[i][j] = _A->matrixGet()[i][j];
            }
        }
    }
    Matrix* checkSolution = gaussLinearSolve(B);
    Matrix* error = Matrix::matrixDiff(checkSolution, _eigenVector);
    cout << "Error vector for eigenvalue = " << _eigenValue << " is" << endl;
    error->matrixPrint();
    delete b;
    delete B;
    delete checkSolution;
    delete error;
}

