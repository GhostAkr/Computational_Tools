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
//    cout << "matrix is" << endl;
//    _matrix->matrixPrint();
//    cout << "comp = " << comp << endl;
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

Matrix* QRDecompositionSolve(Matrix* _A, Matrix* Q, Matrix* R) {  /// This algorithm was changed for QREigen method
    size_t rows = _A->rowsGet();
    size_t cols = _A->colsGet();
    Matrix* rotResultMatrix = new Matrix;
    rotResultMatrix->matrixOneSet(rows, rows);
	type eps = 0.0;  // For comparing with 0.0
	if (sizeof(type) == 4) {  /// Features for Laboratory work (eps = 1e-14 for double type)
		eps = 1e-9;
	}
	else {
		eps = 1e-14;
	}
    for (int j = 0; j < rows - 1; ++j) {
        /// Changing of rows based on finding of main element was commented for QREigen method
        //int mainElem = _A->mainElement(j);
        //_A->matrixRowsChange(j, mainElem);
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
    /// Checking for singular matrix was commented for QREigen method
//    if (!onlyDesitionCheck(_A)) {
//        cout << "Linear system has infinite number of solutions or hasn't it at all" << endl;
//        return NULL;
//    }
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
    for(int i = 0; i < rows; ++i) {  /// Changes were made for QREigen method
        rotResultMatrix->matrixGet()[i][rows - 1] *= -1;
    }
    _A->matrixGet()[rows - 1][rows - 1] *= -1;  /// Changes were made for QREigen method
    Q->matrixNullSet(rows, rows);
    R->matrixNullSet(rows, rows);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < rows; ++j) {
            Q->matrixGet()[i][j] = rotResultMatrix->matrixGet()[i][j];
            R->matrixGet()[i][j] = _A->matrixGet()[i][j];
        }
    }
    delete tmp;
    delete rotResultMatrix;
    delete b;
    return result;
}

Matrix* QRBackTurn(Matrix* _Q, Matrix* _R, Matrix* _b) {
    type eps = 1e-14;
    int rows = _R->rowsGet();
    _Q->matrixTranspose();
    Matrix* bb = Matrix::matrixComp(_Q, _b);
    //cout << "bb is" << endl;
    //bb->matrixPrint();
    Matrix* RR = new Matrix;
    RR->matrixNullSet(rows, rows + 1);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < rows + 1; ++j) {
            if (j != rows + 1) {
                RR->matrixGet()[i][j] = _R->matrixGet()[i][j];
            } else {
                RR->matrixGet()[i][j] = bb->matrixGet()[i][0];
            }
        }
    }
    Matrix* result = new Matrix;
    result->matrixNullSet(rows, 1);
    for (int i = rows - 1; i >= 0; --i) {
        type leftSum = 0.0;
        for (int j = rows - 1; j >= i; --j) {
            leftSum += RR->matrixGet()[i][j] * result->matrixGet()[j][0];
        }
        result->matrixGet()[i][0] = (RR->matrixGet()[i][rows] - leftSum) / RR->matrixGet()[i][i];
        if (fabs(result->matrixGet()[i][0]) < eps) {
            result->matrixGet()[i][0] = 0.0;
        }
    }
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

type norm(Matrix* _A) {
    type norm = 0.0;
    double rows = _A->rowsGet();
    for (int i = 0; i < rows; ++i) {
        norm += _A->matrixGet()[i][0] * _A->matrixGet()[i][0];
    }
    return sqrt(norm);
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
    for (int j = 0; j < _A->colsGet(); ++j) {
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

type normInfVect(Matrix* _A) {
    type max = -1.0;
    for (int i = 0; i < _A->rowsGet(); ++i) {
        if (fabs(_A->matrixGet()[i][0]) > max) {
            max = fabs(_A->matrixGet()[i][0]);
        }
    }
    return max;
}

type normOneVect(Matrix* _A) {
    type sum = 0.0;
    for (int i = 0; i < _A->rowsGet(); ++i) {
        sum += fabs(_A->matrixGet()[i][0]);
    }
    return sum;
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
    for (int i = 0; i < rowsA; i++) {
        if (_A->matrixGet()[i][i] < 0) {
            for (int j = 0; j < colsA; j++) {
                _A->matrixGet()[i][j] *= -1.0;
            }
        }
    }
    Matrix* X = new Matrix;
    X->matrixNullSet(rowsA, 1);
    type tau = 0.0072;
    type eps = 1e-7;
    Matrix* b = new Matrix;
    b->matrixNullSet(rowsA, 1);
    for (int i = 0; i < rowsA; i++) {
        type element =  _A->matrixGet()[i][colsA - 1];
        b->matrixGet()[i][0] = element;
        X->matrixGet()[i][0] = element;
    }
    Matrix* pureA = new Matrix;
    pureA->matrixNullSet(rowsA, rowsA);
    for (int i = 0; i < rowsA; ++i) {
        for (int j = 0; j < rowsA; ++j) {
            pureA->matrixGet()[i][j] = _A->matrixGet()[i][j];
        }
    }
    Matrix* prevX = new Matrix;
    Matrix* E = new Matrix;
    E->matrixOneSet(rowsA, rowsA);
    Matrix* C = Matrix::matrixDiff(E, Matrix::matrixConstComp(pureA, tau));
    size_t iteration = 0;
    type normC = normOne(C);
    do {
        iteration++;
        prevX = Matrix::getCopy(X);
        X = Matrix::matrixSum(Matrix::matrixComp(C, X), Matrix::matrixConstComp(b, tau));
    } while (normOneVect(Matrix::matrixDiff(X, prevX)) > (1.0 - normInf(C)) * eps / normInf(C));
    cout << "Number of iterations is " << iteration << endl;
    return X;
}

Matrix* Jacobi(const Matrix* _matrix) {
    int n = 0;
    double eps = 1e-7;
    size_t rows = _matrix->rowsGet();
    size_t cols = _matrix->colsGet();
    Matrix *X0 = new Matrix;
    X0->matrixNullSet(rows, 1);
    for (int i = 0; i < rows; i++)
        X0->matrixGet()[i][0] = _matrix->matrixGet()[i][cols - 1];
    Matrix *res = new Matrix;
    res->matrixNullSet(rows, 1);
    res->getCopy(X0);
    Matrix *_C = new Matrix;
    _C->matrixNullSet(rows, rows);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < rows; j++)
            if (i != j) {
                _C->matrixGet()[i][j] = -_matrix->matrixGet()[i][j] / _matrix->matrixGet()[i][i];
            }
    }
    Matrix* y = new Matrix;
    y->matrixNullSet(rows, 1);
    for (int i = 0; i < rows; ++i) {
        y->matrixGet()[i][0] = _matrix->matrixGet()[i][cols - 1] / _matrix->matrixGet()[i][i];
    }
    type normC = normInf(_C);
    cout << normInf(_C) << endl;
    do {
        n++;
        X0 = Matrix::getCopy(res);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols - 1; j++)
                for (int i = 0; i < rows; i++) {
                    for (int j = 0; j < cols - 1; j++)
                        if (i != j) {
                            res->matrixGet()[i][0] -= X0->matrixGet()[j][0] * _matrix->matrixGet()[i][j];
                        }
                    res->matrixGet()[i][0] += _matrix->matrixGet()[i][cols - 1];
                    res->matrixGet()[i][0] /= _matrix->matrixGet()[i][i];
                }
        }
    } while (normInfVect((Matrix::matrixDiff(X0, res))) > ((1 - normC) * eps / normC));
    cout << "Number of iterations = " << n << endl;
    return (res);
}

Matrix* SOR(const Matrix* _matrix) {  // TODO: Write break condition
    int cols = _matrix->colsGet();
    int rows = _matrix->rowsGet();
    type omega = 1.001;
    type eps = 1e-4;
    Matrix* result = new Matrix;
    result->matrixNullSet(rows, 1);
    for (int i = 0; i < rows; ++i) {  // Reading of right part of equation
        result->matrixGet()[i][0] = _matrix->matrixGet()[i][cols - 1];
    }
    Matrix* b = Matrix::getCopy(result);
    Matrix* prevX = new Matrix;
    size_t iteration = 0;
    Matrix* A = new Matrix;
    A->matrixNullSet(rows, rows);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < rows; ++j) {
            A->matrixGet()[i][j] = _matrix->matrixGet()[i][j];
        }
    }
    Matrix* C = C_SOR(A, omega);
    type normC = normInf(C);
    cout << "Norm C = " << normC << endl;
    cout << "Norm CU = " << CUNorm(C) << endl;
    cout << "Norm CL = " << CLNorm(C) << endl;
    do {
       iteration++;
       prevX = Matrix::getCopy(result);
       for (int i = 0; i < rows; ++i) {
           type sum1 = 0.0;
           if (i != rows - 1) {
               type coefficient1 = _matrix->matrixGet()[i][2] / _matrix->matrixGet()[i][1];
               sum1 = coefficient1 * result->matrixGet()[i + 1][0];
           }
           type sum2 = 0.0;
           if (i != 0) {
               type coefficient2 = _matrix->matrixGet()[i][0] / _matrix->matrixGet()[i][1];
               sum2 = coefficient2 * result->matrixGet()[i - 1][0];
           }
           type coefficient = b->matrixGet()[i][0] / _matrix->matrixGet()[i][1];
           result->matrixGet()[i][0] = (1.0 - omega) * result->matrixGet()[i][0] - omega * sum1 + \
           omega * coefficient - omega * sum2;
       }
       //if (iteration == 1000) {
       //    break;
       //}
    } while (normInfVect((Matrix::matrixDiff(prevX, result))) > ((1 - normC) * eps / normC));
    delete b;
    return result;
}

Matrix* createTridiagonalMatrix(int variant) {
    int n = 200 + variant;
    Matrix* A = new Matrix;
    A->matrixNullSet(n, 4);
    for (int i = 0; i < n; ++i) {
        if (i == 0) {
            A->matrixGet()[i][0] = 0.0;
        } else {
            A->matrixGet()[i][0] = 1.0;
        }
        A->matrixGet()[i][1] = 4.0;
        if (i == n - 1) {
            A->matrixGet()[i][2] = 0.0;
        } else {
            A->matrixGet()[i][2] = 1.0;
        }
        if (i == 0) {
            A->matrixGet()[i][3] = 6.0;
        } else if (i == n - 1) {
            A->matrixGet()[i][3] = 9 - 3 * (n % 2);
        } else {
            A->matrixGet()[i][3] = 10 - 2 * ((i - 1) % 2);
        }
    }
    return A;
}

Matrix* C_SOR(Matrix* _matrix, type omega){
    int rows = _matrix->rowsGet();
    Matrix* D = new Matrix;
    D->matrixNullSet(rows, rows);
    Matrix* L = new Matrix;
    L->matrixNullSet(rows, rows);
    Matrix* U = new Matrix;
    U->matrixNullSet(rows, rows);
    Matrix* E = new Matrix;
    E->matrixOneSet(rows, rows);
    for (int i = 0; i < rows; ++i) {
        D->matrixGet()[i][i] = _matrix->matrixGet()[i][1];
        if (i != rows - 1) {
            U->matrixGet()[i][i + 1] = _matrix->matrixGet()[i][2];
        }
        if (i != 0) {
            L->matrixGet()[i][i - 1] = _matrix->matrixGet()[i][0];
        }
    }
    Matrix* Q = Matrix::matrixSum(Matrix::matrixConstComp(D, 1 / omega), L);
    Matrix* Q1 = inverseMatrix(Q);
    Matrix* Q2 = Matrix::matrixDiff(Matrix::matrixConstComp(D, (1 / omega - 1)), U);
    Matrix* C = Matrix::matrixComp(Q1, Q2);
    delete Q;
    delete D;
    delete L;
    delete E;
    delete U;
    delete Q1;
    return C;
}

Matrix* inverseMatrix(const Matrix* _matrix) {  // TODO: Make exceptions
    Matrix* Q = new Matrix;
    Matrix* R = new Matrix;
    Matrix* E = new Matrix;  // Identity matrix
    int rows = _matrix->rowsGet();
    E->matrixOneSet(rows, rows);
    Matrix* T = new Matrix;  // Temporary equation
    T->matrixNullSet(rows, rows + 1);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < rows; ++ j) {
            T->matrixGet()[i][j] = _matrix->matrixGet()[i][j];
        }
    }
    Matrix* Result = new Matrix;
    Result->matrixNullSet(rows, rows);
    for (int k = 0; k < rows; ++k) {
        for (int i = 0; i < rows; ++i) {
            if (i != k) {
                T->matrixGet()[i][rows] = 0.0;
            } else {
                T->matrixGet()[i][rows] = 1.0;
            }
        }
        Matrix* b = new Matrix;
        b->matrixNullSet(rows, 1);
        for (int i = 0; i < rows; ++i) {
            b->matrixGet()[i][0] = T->matrixGet()[i][rows];
        }
        Matrix* X = new Matrix;
        if (k == 0) {
            X = QRDecompositionSolve(T, Q, R);
        } else {
            X = QRBackTurn(Q, R, b);
        }
        for (int i = 0; i < rows; ++i) {
            Result->matrixGet()[i][k] = X->matrixGet()[i][0];
        }
    }
    delete E;
    delete T;
    return Result;
}


Matrix* CSeidel(const Matrix* _matrix) {
	int cols = _matrix->colsGet();
	int rows = _matrix->rowsGet();
	Matrix* Q = new Matrix;
	Matrix* D = new Matrix;
	D->matrixNullSet(rows, rows);
	Matrix* L = new Matrix;
	L->matrixNullSet(rows, rows);
	Matrix* U = new Matrix;
	U->matrixNullSet(rows, rows);
	Matrix* E = new Matrix;
	E->matrixOneSet(rows, rows);
	for (int i = 0; i < rows; ++i) {
		D->matrixGet()[i][i] = _matrix->matrixGet()[i][1];
		if (i != rows - 1) {
			U->matrixGet()[i][i + 1] = _matrix->matrixGet()[i][2];
		}
		if (i != 0) {
			L->matrixGet()[i][i - 1] = _matrix->matrixGet()[i][0];
		}
	}
	Q = Matrix::matrixSum(D,L);
	Q = inverseMatrix(Q);
	Q = Matrix::matrixConstComp(Q, -1);
	Q = Matrix::matrixComp(Q,U);
	//Q = -(A - U) ^ -1 * U;
	return(Q);
}

Matrix* Seidel(const Matrix* _matrix) {
    int cols = _matrix->colsGet();
    int rows = _matrix->rowsGet();
    type eps = 1e-7;
    Matrix* result = new Matrix;
    result->matrixNullSet(rows, 1);
    for (int i = 0; i < rows; ++i) {  // Reading of right part of equation
        result->matrixGet()[i][0] = _matrix->matrixGet()[i][cols - 1];
    }
    Matrix* b = Matrix::getCopy(result);
    Matrix* prevX = new Matrix;
    size_t iteration = 0;
    Matrix* C = CSeidel(_matrix);
    type normC = normInf(C);
    cout << "Norm C = " << normC << endl;
    cout << "Norm CU = " << CUNorm(C) << endl;
    cout << "Norm CL = " << CLNorm(C) << endl;
    do {
        iteration++;
        prevX = Matrix::getCopy(result);
        result->matrixGet()[0][0] =-_matrix->matrixGet()[0][2]* result->matrixGet()[1][0]+b->matrixGet()[0][0];
        result->matrixGet()[0][0] /= _matrix->matrixGet()[0][1];
        for (int i = 1; i < rows - 1; i++) {
            result->matrixGet()[i][0] = b->matrixGet()[i][0] - _matrix->matrixGet()[i - 1][0] * result->matrixGet()[i - 1][0] - _matrix->matrixGet()[i][2] * result->matrixGet()[i + 1][0];
            result->matrixGet()[i][0] /= _matrix->matrixGet()[i][1];
        }
        result->matrixGet()[rows - 1][0] = -_matrix->matrixGet()[rows - 1][0] * result->matrixGet()[rows - 2][0] + b->matrixGet()[rows - 1][0];
        result->matrixGet()[rows - 1][0] /= _matrix->matrixGet()[rows - 1][1];
    } while (normInfVect((Matrix::matrixDiff(prevX, result))) > ((1 - normC) * eps / normC));
    cout << "Number of iterations is " << iteration << endl;
    delete b;
    delete prevX;
    delete C;
    return result;
}

type pogrNorm(Matrix* _solution, Matrix* _realSolution){
    Matrix* pogr = Matrix::matrixDiff(_solution, _realSolution);
    type norm = normInfVect(pogr);
    return norm;
}

type CLNorm(Matrix* _C) {
    int rows = _C->rowsGet();
    Matrix* tmp = Matrix::getCopy(_C);
    for (int i = 0; i < rows; ++i) {
        for (int j = i; j < rows; ++j) {
            tmp->matrixGet()[i][j] = 0.0;
        }
    }
    type norm = normInf(tmp);
    delete tmp;
    return norm;
}

type CUNorm(Matrix* _C) {
    int rows = _C->rowsGet();
    Matrix* tmp = Matrix::getCopy(_C);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < i + 1; ++j) {
            tmp->matrixGet()[i][j] = 0.0;
        }
    }
    type norm = normInf(tmp);
    delete tmp;
    return norm;
}

Matrix* tridiagonalLinearSolve(const Matrix* _matrix) {
    int rows = _matrix->rowsGet();
    cout << "rows = " << rows << endl;
    auto* a = new double [rows - 1];
    auto* b = new double [rows - 1];
    a[0] = -_matrix->matrixGet()[0][2] / _matrix->matrixGet()[0][1];
    b[0] = _matrix->matrixGet()[0][3] / _matrix->matrixGet()[0][1];
    Matrix* result = new Matrix;
    result->matrixNullSet(rows, 1);
    for (int i = 1; i < rows - 1; ++i) {
        double A = _matrix->matrixGet()[i][0];
        double C = _matrix->matrixGet()[i][1];
        double B = _matrix->matrixGet()[i][2];
        double F = _matrix->matrixGet()[i][3];
        a[i] = -B / (A * a[i - 1] + C);
        b[i] = (F - A * b[i - 1]) / (A * a[i - 1] + C);
    }
//    cout << "Alpha is" << endl;
//    for (int i = 0; i < rows - 1; ++i) {
//        cout << a[i] << endl;
//    }
//    cout << endl;
    double F = _matrix->matrixGet()[rows - 1][3];
    double A = _matrix->matrixGet()[rows - 1][0];
    double C = _matrix->matrixGet()[rows - 1][1];
    result->matrixGet()[rows - 1][0] = (F - A * b[rows - 2]) / (A * a[rows - 2] + C);
    for (int i = rows - 2; i >= 0; --i) {
        result->matrixGet()[i][0] = a[i] * result->matrixGet()[i + 1][0] + b[i];
    }
    delete[] a;
    delete[] b;
//    cout << "Result is " << endl;
//    result->matrixPrint();
    return result;
}
