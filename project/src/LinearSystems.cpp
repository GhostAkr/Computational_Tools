//
// Created by ghostakr on 9/16/18.
//

#include "../include/LinearSystems.h"

Matrix* gaussLinearSolve(Matrix* _A) {
	type eps = 0.0;  // For comparing with 0.0
	if (sizeof(type) == 4) {
		eps = 1e-9;
	} else {
		eps = 1e-14;
	}
    size_t rowsA = _A->rowsGet();
    size_t colsA = _A->colsGet();
    // First part
    for (int k = 0; k < colsA - 1; ++k) {
        int mainElem = _A->mainElement(k);
        _A->matrixRowsChange(mainElem, k);
        for (int i = k + 1; i < rowsA; ++i) {
            type coeff = _A->matrixGet()[i][k] / _A->matrixGet()[k][k];
            for (int j = k + 1; j < colsA; ++j) {
				_A->matrixGet()[i][j] -= _A->matrixGet()[k][j] * coeff;
				_A->matrixGet()[i][k] = 0.0;
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
        type leftSum = 0.0;
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
	type eps = 0.0;  // For comparing with 0.0
	if (sizeof(type) == 4) {
		eps = 1e-9;
	}
	else {
		eps = 1e-14;
	}
    type comp = 1.0;  // Composition of elements on the main diagonal
    for (int i = 0; i < _matrix->rowsGet(); ++i) {
        comp *= _matrix->matrixGet()[i][i];
    }
    if (fabs(comp) < eps) {  // if (comp == 0.0)
        return false;
    } else {
        return true;
    }
}

type valuationVector(Matrix* _solution, Matrix* _system) {
    if (_solution->rowsGet() != _system->rowsGet()) {  // Exception
        cout << "Solution and system are not compatible" << endl;
        return -1.0;
    }
    size_t rows = _solution->rowsGet();
    type* b = new type [rows];
    for (int i = 0; i < rows; ++i) {
        b[i] = 0.0;
        for (int j = 0; j < rows; ++j) {
			type k = _solution->matrixGet()[j][0] * _system->matrixGet()[i][j];
			b[i] += k;
        }
    }
    type* residualVector = new type [rows];
    for (int i = 0; i < rows; ++i) {
        residualVector[i] = _system->matrixGet()[i][rows] - b[i];
    }
    type residual = vectorNorm(residualVector, rows);
	printf("Residual is %.20f\n", residual);
	delete b;
	delete residualVector;
	return residual;
}

Matrix* rotationMatrix(Matrix* _matrix, int line, int numOfVar) {
    type a1 = _matrix->matrixGet()[numOfVar][numOfVar];  // Temporary variable
    type a2 = _matrix->matrixGet()[line][numOfVar];  // Temporary variable
	type denom = sqrt(a1 * a1 + a2 * a2);
    type c = a1 / denom;
    type s = a2 / denom;
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
	type eps = 0.0;  // For comparing with 0.0
	if (sizeof(type) == 4) {
		eps = 1e-9;
	}
	else {
		eps = 1e-14;
	}
    for (int j = 0; j < rows - 1; ++j) {
        int mainElem = _A->mainElement(j);
        _A->matrixRowsChange(j, mainElem);
        for (int i = j + 1; i < rows; ++i) {
			type a1 = _A->matrixGet()[j][j];  // Temporary variable
			type a2 = _A->matrixGet()[i][j];  // Temporary variable
			type denominator = sqrt(a1 * a1 + a2 * a2);
			type c = a1 / denominator;
			type s = a2 / denominator;
			rotateMatrix(_A, i, j, c, s);
			rotateMatrix(rotResultMatrix, i, j, c, s);
			_A->matrixGet()[i][j] = 0.0;
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
        type leftSum = 0.0;
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

type vectorNorm(type* _vector, size_t _rows) {
    type norm = 0.0;
    for (int i = 0; i < _rows; ++i) {
		norm += _vector[i] * _vector[i];
    }
    norm = sqrt(norm);
    return norm;
}

void conditionNumber(Matrix* _A) {
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
        Q->matrixGet()[k][colsA - 1] = 1.0;
        Matrix* X = gaussLinearSolve(Q);
        for (int j = 0; j < rowsA; ++j) {
            A1->matrixGet()[j][k] = X->matrixGet()[j][0];
        }
        delete Q;
        delete X;
    }
    for (int i = 0; i < rowsA; ++i) {
        _A->matrixGet()[i][colsA - 1] = B->matrixGet()[i][0];
    }
    type condOne = normOne(_A) * normOne(A1);
    type condInf = normInf(_A) * normInf(A1);
	printf("Condition number (1) = %.14f\n", condOne);
	printf("Condition number (inf) = %.14f\n", condInf);
	delete A1;
	delete B;
}

type normInf(Matrix* _A) {
    type max = 0.0;
    for (int j = 0; j < _A->colsGet() - 1; ++j) {
        max += fabs(_A->matrixGet()[0][j]);
    }
    for (int i = 1; i < _A->rowsGet(); ++i) {
        type sum = 0.0;
        for (int j = 0; j < _A->colsGet() - 1; ++j) {
            sum += fabs(_A->matrixGet()[i][j]);
        }
        if (sum > max) {
            max = sum;
        }
    }
    return max;
}

type normOne(Matrix* _A) {
    type max = 0.0;
    for (int i = 0; i < _A->rowsGet(); ++i) {
        max += fabs(_A->matrixGet()[i][0]);
    }
    for (int j = 0; j < _A->colsGet() - 1; ++j) {
        type sum = 0.0;
        for (int i = 0; i < _A->rowsGet(); ++i) {
            sum += fabs(_A->matrixGet()[i][j]);
        }
        if (sum > max) {
            max = sum;
        }
    }
    return max;
}

Matrix* rotateMatrix(Matrix* _targetMatrix, int _line, int _numOfVar, type _c, type _s) {
	// Calculting coefficients
	for (int k = 0; k < _targetMatrix->colsGet(); ++k) {
		type sum1 = _c * _targetMatrix->matrixGet()[_numOfVar][k];
		sum1 += _s * _targetMatrix->matrixGet()[_line][k];
		type sum2 = -_s * _targetMatrix->matrixGet()[_numOfVar][k];
		sum2 += _c * _targetMatrix->matrixGet()[_line][k];
		_targetMatrix->matrixGet()[_numOfVar][k] = sum1;
		_targetMatrix->matrixGet()[_line][k] = sum2;
	}
	return _targetMatrix;
}

void pertrubationSolution(Matrix* _A) {
	Matrix* A1 = new Matrix;  // Additional matrix for counting delta X
	type pertrubation = 0.01;
	size_t rows = _A->rowsGet();
	size_t cols = _A->colsGet();
	A1->matrixNullSet(rows, cols);
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			A1->matrixGet()[i][j] = _A->matrixGet()[i][j];
		}
	}
	Matrix* b = new Matrix;
	b->matrixNullSet(rows, 1);
	for (int k = 0; k < rows; ++k) {
		b->matrixGet()[k][0] = _A->matrixGet()[k][cols - 1];
	}
	type bNormOne = normOne(b);
	type bNormInf = normInf(b);
	for (int k = 0; k < rows; ++k) {
		_A->matrixGet()[k][cols - 1] += pertrubation;
	}
	for (int k = 0; k < rows; ++k) {
		A1->matrixGet()[k][cols - 1] = pertrubation;
	}
	Matrix* deltaB = new Matrix;
	deltaB->matrixNullSet(rows, 1);
	for (int k = 0; k < rows; ++k) {
		deltaB->matrixGet()[k][0] = pertrubation;
	}
	type deltaBNormOne = normOne(deltaB);
	type deltaBNormInf = normInf(deltaB);
	Matrix* solution = gaussLinearSolve(_A);
	cout << "Petrubation solution is " << endl;
	solution->matrixPrint();
	type xNormOne = normOne(solution);
	type xNormInf = normInf(solution);
	Matrix* deltaX = gaussLinearSolve(A1);
	type deltaXNormInf = normInf(deltaX);
	type deltaXNormOne = normOne(deltaX);
	type evaluationOne = 0.0;
	type evaluationInf = 0.0;
	type dXOne = deltaXNormOne / xNormOne;
	type dBOne = deltaBNormOne / bNormOne;
	type dXInf = deltaXNormOne / xNormInf;
	type dBInf = deltaBNormOne / bNormInf;
	cout << "dXOne = " << dXOne << endl;
	cout << "dBOne = " << dBOne << endl;
	evaluationOne = dXOne / dBOne;
	evaluationInf = dXInf / dBInf;
	cout << "Evaluation (one): " << evaluationOne << endl;
	//cout << "Evaluation (inf): " << evaluationInf << endl;
	delete A1;
	delete b;
	delete deltaB;
	delete solution;
	delete deltaX;
}

Matrix* fixedPointIterationSolve(Matrix* _A) {
    int rowsA = _A->rowsGet();
    int colsA = _A->colsGet();
    Matrix* X = new Matrix;
    X->matrixNullSet(rowsA, 1);
    type tau = 0.5;
    type eps = 0.5;
    Matrix* b = new Matrix;
    b->matrixNullSet(rowsA, 1);
    for (int i = 0; i < rowsA; i++) {
        b->matrixGet()[i][0] = _A->matrixGet()[i][colsA-1];
    }
    Matrix* pureA = new Matrix;
    pureA->matrixNullSet(rowsA, rowsA);
    for (int i = 0; i < rowsA; ++i) {
        for (int j = 0; j < rowsA; ++j) {
            pureA->matrixGet()[i][j] = _A->matrixGet()[i][j];
        }
    }
    Matrix* prevX = new Matrix;
    //Matrix* C = new Matrix;
    Matrix* E = new Matrix;
    E->matrixOneSet(rowsA, rowsA);
    Matrix* C = Matrix::matrixConstComp(Matrix::matrixDiff(Matrix::matrixConstComp(pureA, tau), E), -1);
    do {
        delete prevX;
        prevX = Matrix::getCopy(X);
        X = Matrix::matrixDiff(X, Matrix::matrixDiff(Matrix::matrixConstComp(Matrix::matrixComp(pureA, X), tau), Matrix::matrixConstComp(b, tau)));
    } while (normOne(Matrix::matrixDiff(X, prevX)) > (1.0 - normOne(C)) * eps / normOne(C));
    return X;
}
