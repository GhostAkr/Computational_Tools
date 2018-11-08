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
    Matrix* result = SOR(A1);
    cout << "Result is" << endl;
    result->matrixPrint();
    return 0;
}