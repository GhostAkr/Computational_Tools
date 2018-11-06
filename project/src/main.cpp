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
    cout << "A is" << endl;
    A->matrixPrint();
    //Matrix* result = Jacobi(A);
    //Matrix* result = fixedPointIterationSolve(A);
    //Matrix* result = gaussLinearSolve(A);
    Matrix* result = SOR(A);
    cout << "Result is" << endl;
    result->matrixPrint();
    return 0;
}