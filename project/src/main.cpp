#include <iostream>
#include "../include/Matrix.h"
#include "../include/LinearSystems.h"

using std::cout;
using std::endl;

int main() {
    Matrix* A = new Matrix;
    Matrix* B = new Matrix;
    A->readMatrixFromFile("../data/matrix_1");
    B->readMatrixFromFile("../data/matrix_2");
    cout << "Matrix A is " << endl;
    A->matrixPrint();
    cout << "Matrix B is " << endl;
    B->matrixPrint();
    gaussLinearSolve(A, B);
    delete A;
    delete B;
    return 0;
}