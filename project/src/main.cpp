#include <iostream>
#include "../include/Matrix.h"
#include "../include/LinearSystems.h"

using std::cout;
using std::endl;

int main() {
    Matrix* A = new Matrix;
    A->readMatrixFromFile("../data/matrix_1");
    cout << "Matrix A is " << endl;
    A->matrixPrint();
    Matrix* Result = gaussLinearSolve(A);
    cout << "Result is " << endl;
    Result->matrixPrint();
    delete A;
    return 0;
}