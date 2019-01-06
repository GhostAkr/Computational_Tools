//
// Created by ighos on 18.12.2018.
//

#include "../include/EigenValues.h"

Matrix* HessenbergForm(Matrix* _A) {
    int n = _A->rowsGet();
    Matrix* Result = Matrix::getCopy(_A);
    for (int k = 0; k < n - 2; ++k) {
        for (int l = k + 2; l < n; ++l) {
            double sqr = sqrt(_A->matrixGet()[k + 1][k] * _A->matrixGet()[k + 1][k] + \
                    _A->matrixGet()[l][k] *_A->matrixGet()[l][k]);
            double alpha = _A->matrixGet()[k + 1][k] / sqr;
            double beta = _A->matrixGet()[l][k] / sqr;
            // TODO: Make algorithm more effective
            Matrix* T = new Matrix;
            T->matrixOneSet(n, n);
            T->matrixGet()[k + 1][k + 1] = alpha;
            T->matrixGet()[k + 1][l] = beta;
            T->matrixGet()[l][k + 1] = -beta;
            T->matrixGet()[l][l] = alpha;
            Matrix* TT = new Matrix;
            TT->matrixOneSet(n, n);
            TT->matrixGet()[k + 1][k + 1] = alpha;
            TT->matrixGet()[k + 1][l] = -beta;
            TT->matrixGet()[l][k + 1] = beta;
            TT->matrixGet()[l][l] = alpha;
            _A = Matrix::matrixComp(T, _A);
            _A = Matrix::matrixComp(_A, TT);
//            cout << "A2 is " << endl;
//            _A->matrixPrint();
//            for (int i = 0; i < n; ++i) {
//                Result->matrixGet()[k + 1][i] = alpha * _A->matrixGet()[k + 1][i] + beta * _A->matrixGet()[l][i];
//                Result->matrixGet()[l][i] = -beta * _A->matrixGet()[k + 1][i] + alpha * _A->matrixGet()[l][i];
//            }
//            cout << "AA1 is " << endl;
//            Result->matrixPrint();
//            for (int i = 0; i < n; ++i) {
//                Result->matrixGet()[i][k + 1] = alpha * _A->matrixGet()[i][k + 1] + beta * _A->matrixGet()[i][l];
//                Result->matrixGet()[i][l] = -beta * _A->matrixGet()[i][k + 1] + alpha * _A->matrixGet()[i][l];
//            }
//            cout << "AA2 is " << endl;
//            Result->matrixPrint();
        }
    }
    return _A;
}

Matrix* QRDecompositionEigen(Matrix* _A) {
    double epsilon = 1e-14;  // For comparing with null
    int rows = _A->rowsGet();
    Matrix* Result = new Matrix;
    Result->matrixNullSet(rows, 1);
    Matrix* CurA = Matrix::getCopy(_A);
    for (int i = 0; i < rows; ++i) {
        cout << "In 'for' cycle" << endl;
        Matrix* Buff = Matrix::getCopy(CurA);  // Buffer matrix for creation of new A
        CurA = new Matrix;
        CurA->matrixNullSet(rows - i, rows - i);
        for (int k = 0; k < rows - i; ++k) {  // Creation of new A
            for (int l = 0; l < rows - i; ++l) {
                CurA->matrixGet()[k][l] = Buff->matrixGet()[k][l];
            }
        }
        cout << "Created matrix A is " << endl;
        CurA->matrixPrint();
        if (i == rows - 1) {  // TODO: Treat last eigenvalue
            cout << "Last is" << endl;
            CurA->matrixPrint();
            break;
        }
        while(true) {  // Cycle for each eigenvalue
            cout << "Current matrix in 'while' cycle" << endl;
            CurA->matrixPrint();
            double shift = CurA->matrixGet()[rows - 1 - i][rows - 1 - i];
            bool IsNextIteration = false;
            Matrix* E = new Matrix;  // Identity matrix
            E->matrixOneSet(rows - i, rows - i);
            CurA = Matrix::matrixDiff(CurA, Matrix::matrixConstComp(E, shift));  // Shifting
            Matrix* Q = new Matrix;
            Matrix* R = new Matrix;
            Matrix* copyA = Matrix::getCopy(CurA);  // Copy of A for QR-Decomposition
            QRDecompositionSolve(copyA, Q, R);
            cout << "Current A is" << endl;
            CurA->matrixPrint();
            cout << "Q matrix is" << endl;
            Q->matrixPrint();
            cout << "R matrix is" << endl;
            R->matrixPrint();
            cout << "Q * R = " << endl;
            Matrix::matrixComp(Q, R)->matrixPrint();
            //CurA = Matrix::matrixSum(Matrix::matrixComp(R, Q), Matrix::matrixConstComp(E, shift));  // Reverse shifting
            CurA = Matrix::matrixComp(R, Q);
            CurA = Matrix::matrixSum(CurA, Matrix::matrixConstComp(E, shift));
            for (int k = 0; k < rows - i - 1; ++k) {  // Checking for break condition
               if (fabs(CurA->matrixGet()[rows - i - 1][k]) > epsilon) {
                   IsNextIteration = true;
               }
            }
            if (!IsNextIteration) {
                cout << "Last matrix with current size is" << endl;
                CurA->matrixPrint();
               Result->matrixGet()[i][0] = CurA->matrixGet()[rows - i - 1][rows - i - 1];
               break;
            }
        }
    }
    return Result;
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

	cout << endl;
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
		t->matrixPrint();
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
		res->matrixPrint();
		cout << endl;
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
		cout << endl;
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
		temp->matrixPrint();
		res = gaussLinearSolve(temp);
		res->matrixPrint();
		cout << endl << eigen << endl;
		double div = norm(res);
		for (int i = 0; i < rows; i++) {
			res->matrixGet()[i][0] /= div;
		}
		n++;
	}
	return eigen;
}