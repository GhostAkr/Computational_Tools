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
    if (!onlyDesitionCheck(_A)) {
        cout << "Linear system has infinite number of solutions or hasn't it at all" << endl;
        return NULL;
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

bool onlyDesitionCheck(Matrix* _matrix) {
    double eps = 10e-5;  // For comparing with 0
    double comp = 1.0;  // Composition of elements on the main diagonal
    for (int i = 0; i < _matrix->rowsGet(); ++i) {
        comp *= _matrix->matrixGet()[i][i];
    }
    if (fabs(comp) < eps) {  // if (comp == 0.0)
        return false;
    } else {
        return true;
    }
}

int valuationVector(Matrix* _solution, Matrix* _system) {  // TODO: Finish valuationVector
    double eps = 10e-5;  // For comparing with 0
    if (_solution->rowsGet() != _system->rowsGet()) {  // Exception
        cout << "Solution and system are not compatible" << endl;
        return -1;
    }
    size_t rows = _solution->rowsGet();
    float* bFloat = new float [rows];  // For calculations with usual accuracy
    double* bDouble = new double [rows];  // For calculations with high accuracy
    for (int i = 0; i < rows; ++i) {
        bFloat[i] = 0.0;
        bDouble[i] = 0.0;
        for (int j = 0; j < rows; ++j) {
            bFloat[i] += _solution->matrixGet()[j][0] * _system->matrixGet()[i][j];
            bDouble[i] += _solution->matrixGet()[j][0] * _system->matrixGet()[i][j];
        }
        //if (fabs(bFloat[i]) < eps) {
        //    bFloat[i] = 0.0;
        //}
        //if (fabs(bDouble[i]) < eps) {
        //    bDouble[i] = 0.0;
        //}
    }
    cout.precision(5);
    cout << "b1 with usual accuracy" << endl;
    for (int i = 0; i < rows; ++i) {
        cout << bFloat[i] << " ";
    }
    cout << endl;
    cout << "b1 with high accuracy" << endl;
    for (int i = 0; i < rows; ++i) {
        cout << bDouble[i] << " ";
    }
    cout << endl;
}

Matrix* rotationMatrix(Matrix* _matrix, int line, int numOfVar) {
    double a1 = _matrix->matrixGet()[numOfVar][numOfVar];  // Temporary variable
    double a2 = _matrix->matrixGet()[line][numOfVar];  // Temporary variable
    double c = a1 / sqrt(pow(a1, 2) + pow(a2, 2));
    double s = a2 / sqrt(pow(a1, 2) + pow(a2, 2));
    Matrix* result = new Matrix;
    size_t rows = _matrix->rowsGet();
    result->matrixNullSet(rows, rows);
    for (int i = 0; i < rows; ++i) {
        result->matrixGet()[i][i] = 1.0;
    }
    result->matrixGet()[numOfVar][numOfVar] = c;
    result->matrixGet()[numOfVar][line] = s;
    result->matrixGet()[line][numOfVar] = -s;
    result->matrixGet()[line][line] = c;
    return result;
}

Matrix* QRDecompositionSolve(Matrix* _A) {
    size_t rows = _A->rowsGet();
    size_t cols = _A->colsGet();
    Matrix* rotResultMatrix = new Matrix;
    rotResultMatrix->matrixOneSet(rows, rows);
    double eps = 10e-5;
    for (int j = 0; j < rows - 1; ++j) {
        for (int i = j + 1; i < rows; ++i) {
            Matrix* rotMatrix = rotationMatrix(_A, i, j);  // Getting of rotation matrix with current _A
            //DecMatrix* tempMatr = rotResultMatrix;  // For clearing of memory
            rotResultMatrix = Matrix::matrixComp(rotMatrix, rotResultMatrix);
            //delete tempMatr;
            _A = Matrix::matrixComp(rotMatrix, _A);
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < _A->colsGet(); ++j) {
                    if (fabs(_A->matrixGet()[i][j]) < eps) {
                        _A->matrixGet()[i][j] = 0.0;
                    }
                }
            }
            delete rotMatrix;
        }
    }
    for (int i = 0; i < rows; ++i) {  // Setting 0 on almost null elements
        for (int j = 0; j < _A->colsGet(); ++j) {
            if (fabs(_A->matrixGet()[i][j]) < eps) {
                _A->matrixGet()[i][j] = 0.0;
            }
        }
    }
    // Solving final equations
    Matrix* b = new Matrix;
    b->matrixNullSet(rows, 1);
    for (int i = 0; i < rows; ++i) {
        b->matrixGet()[i][0] = _A->matrixGet()[i][rows];
    }
    Matrix* tmp = Matrix::matrixComp(rotResultMatrix, b);
    Matrix* result = new Matrix;
    result->matrixNullSet(rows, 1);
    for (int i = rows - 1; i >= 0; --i) {
        double leftSum = 0.0;
        for (int j = rows - 1; j >= i; --j) {
            leftSum += _A->matrixGet()[i][j] * result->matrixGet()[j][0];
        }
        result->matrixGet()[i][0] = (_A->matrixGet()[i][cols - 1] - leftSum) / _A->matrixGet()[i][i];
        if (fabs(result->matrixGet()[i][0]) < eps) {
            result->matrixGet()[i][0] = 0.0;
        }
    }
    rotResultMatrix->matrixTranspose();
    cout << "Q-Matrix is" << endl;
    rotResultMatrix->matrixPrint();
    cout << "R-Matrix is" << endl;
    _A->matrixPrint();
    delete tmp;
    delete rotResultMatrix;
    delete b;
    return result;
}