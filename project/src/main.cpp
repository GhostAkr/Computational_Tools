#include <iostream>
#include "../include/Matrix.h"
#include "../include/LinearSystems.h"

using std::cout;
using std::endl;

int main() {
    std::string path = "../data/D6.TXT";
    Matrix* A = new Matrix;
    A->readLinearSystemFromFile(path);
    Matrix* A1 = createTridiagonalMatrix(1);
    //cout << "A1 is" << endl;
    //A1->matrixPrint();
    //Matrix* result = Jacobi(A);
    //Matrix* result = fixedPointIterationSolve(A);
    //Matrix* result = gaussLinearSolve(A);
    //Matrix* result = SOR(A1);
    Matrix* result = Seidel(A1);
//    Matrix* b = new Matrix;
//    b->matrixNullSet(A->rowsGet(), 1);
//    for (int i = 0; i < A->rowsGet(); ++i) {
//        b->matrixGet()[i][0] = A->matrixGet()[i][A->colsGet() - 1];
//    }
//    cout << "b is" << endl;
//    b->matrixPrint();
    //Matrix* Q = new Matrix;
    //Matrix* R = new Matrix;
    //Matrix* result = QRDecompositionSolve(A, Q, R);
    //QRBackTurn(Q, R, b);
    cout << "Result is" << endl;
    result->matrixPrint();
    return 0;
}