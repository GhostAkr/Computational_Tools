#include <iostream>
#include "../include/Matrix.h"
#include "../include/LinearSystems.h"

using std::cout;
using std::endl;

int main() {
    // TODO: Fix output with wrong result
    Matrix* A = new Matrix;
    A->readMatrixFromFile("../data/matrix_1");
    cout << "Matrix A is " << endl;
    A->matrixPrint();
    Matrix* X = QRDecompositionSolve(A);
    cout << "Solution is " << endl;
    X->matrixPrint();
    delete X;
    delete A;
    return 0;
}