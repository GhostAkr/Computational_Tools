//
// Created by ghostakr on 9/16/18.
//

#include "../include/LinearSystems.h"

Matrix* gaussLinearSolve(Matrix* _A) {
    double eps = 1e-14;  // For comparing with 0.0
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
    double eps = 1e-10;  // For comparing with 0
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

int valuationVector(Matrix* _solution, Matrix* _system) {
    //double eps = 10e-5;  // For comparing with 0
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
            bFloat[i] += (float)_solution->matrixGet()[j][0] * (float)_system->matrixGet()[i][j];
            bDouble[i] += _solution->matrixGet()[j][0] * _system->matrixGet()[i][j];
        }
    }
    float* residualVectorFloat = new float [rows];
    double* residualVectorDouble = new double [rows];
    for (int i = 0; i < rows; ++i) {
        residualVectorFloat[i] = (float)_system->matrixGet()[i][rows] - bFloat[i];
        residualVectorDouble[i] = _system->matrixGet()[i][rows] - bDouble[i];
    }
    cout.precision(10);
    cout << "Residual with usual accuracy is " << vectorNormFloat(residualVectorFloat, rows) << endl;
    cout << "Residual with high accuracy is " << vectorNormDouble(residualVectorDouble, rows) << endl;
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
    double eps = 1e-14;
    for (int j = 0; j < rows - 1; ++j) {
        int mainElem = _A->mainElement(j);
        _A->matrixRowsChange(j, mainElem);
        for (int i = j + 1; i < rows; ++i) {
            Matrix* rotMatrix = rotationMatrix(_A, i, j);  // Getting of rotation matrix with current _A
            rotResultMatrix = Matrix::matrixComp(rotMatrix, rotResultMatrix);
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
    if (!onlyDesitionCheck(_A)) {
        cout << "Linear system has infinite number of solutions or hasn't it at all" << endl;
        return NULL;
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
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < rows; ++j) {
            if (fabs(rotResultMatrix->matrixGet()[i][j]) < eps) {
                rotResultMatrix->matrixGet()[i][j] = 0.0;
            }
        }
    }
    rotResultMatrix->matrixTranspose();
    cout << "Q-Matrix is" << endl;
    rotResultMatrix->matrixPrint();
    cout << "R-Matrix is" << endl;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < rows; ++j) {
            cout << _A->matrixGet()[i][j] << " ";
        }
        cout << endl;
    }
    delete tmp;
    delete rotResultMatrix;
    delete b;
    return result;
}

double vectorNormDouble(double* _vector, size_t _rows) {
    double norm = 0.0;
    for (int i = 0; i < _rows; ++i) {
        norm += pow(_vector[i], 2);
    }
    norm = sqrt(norm);
    return norm;
}

float vectorNormFloat(float* _vector, size_t _rows) {
    float norm = 0.0;
    for (int i = 0; i < _rows; ++i) {
        norm += pow(_vector[i], 2);
    }
    norm = sqrt(norm);
    return norm;
}

double conditionNumber(Matrix* _A) {
    size_t rowsA = _A->rowsGet();
    size_t colsA = _A->colsGet();
    Matrix* A1 = new Matrix;
    A1->matrixNullSet(rowsA, colsA);
    Matrix* B = new Matrix;
    B->matrixNullSet(rowsA, 1);  // Old right part of A
    for (int i = 0; i < rowsA; ++i) {
        B->matrixGet()[i][0] = _A->matrixGet()[i][colsA - 1];
    }
    for (int k = 0; k < rowsA; ++k) {
        Matrix* Q = new Matrix;
        Q->matrixNullSet(rowsA, colsA);
        for (int i = 0; i < rowsA; ++i) {
            for (int j = 0; j < rowsA; ++j) {
                Q->matrixGet()[i][j] = _A->matrixGet()[i][j];
            }
        }
        for (int i = 0; i < rowsA; ++i) {
            Q->matrixGet()[i][colsA - 1] = 0.0;
        }
        Q->matrixGet()[k][colsA - 1] = 1.0;
        Matrix* X = gaussLinearSolve(Q);
        for (int j = 0; j < rowsA; ++j) {
            A1->matrixGet()[j][k] = X->matrixGet()[j][0];
        }
        delete Q;
    }
    for (int i = 0; i < rowsA; ++i) {
        _A->matrixGet()[i][colsA - 1] = B->matrixGet()[i][0];
    }
    double condOne = normOne(_A) * normOne(A1);
    double condInf = normInf(_A) * normInf(A1);
    cout << "Condition number (1) = " << condOne << endl;
    cout << "Condition number (inf) = " << condInf << endl;
}

double normInf(Matrix* _A) {
    double max = 0.0;
    for (int j = 0; j < _A->colsGet() - 1; ++j) {
        max += fabs(_A->matrixGet()[0][j]);
    }
    for (int i = 1; i < _A->rowsGet(); ++i) {
        double sum = 0.0;
        for (int j = 0; j < _A->colsGet() - 1; ++j) {
            sum += fabs(_A->matrixGet()[i][j]);
        }
        if (sum > max) {
            max = sum;
        }
    }
    return max;
}

double normOne(Matrix* _A) {
    double max = 0.0;
    for (int i = 0; i < _A->rowsGet(); ++i) {
        max += _A->matrixGet()[i][0];
    }
    for (int j = 0; j < _A->colsGet() - 1; ++j) {
        double sum = 0.0;
        for (int i = 0; i < _A->rowsGet(); ++i) {
            sum += fabs(_A->matrixGet()[i][j]);
        }
        if (sum > max) {
            max = sum;
        }
    }
    return max;
}
