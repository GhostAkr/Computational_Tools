#include <iostream>
#include "../include/Matrix.h"

using std::cout;
using std::endl;

int main() {
    Matrix* A = new Matrix;
    Matrix* B = new Matrix;
    A->readMatrixFromFile("../data/matrix_1");
    B->readMatrixFromFile("../data/matrix_2");
    Matrix* C = Matrix::matrixComp(A, B);
    cout << "Matrix A is " << endl;
    A->matrixPrint();
    cout << "Matrix B is " << endl;
    B->matrixPrint();
    cout << "Matrix C is " << endl;
    C->matrixPrint();
    return 0;
}