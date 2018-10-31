#include <iostream>
#include "../include/Matrix.h"
#include "../include/LinearSystems.h"

using std::cout;
using std::endl;

int main() {
    std::string path = "../data/D6.txt";
    Matrix* A = new Matrix;
    A->readLinearSystemFromFile(path);
    cout << "A is" << endl;
    A->matrixPrint();
    Matrix* result = fixedPointIterationSolve(A);
    cout << "Result is" << endl;
    result->matrixPrint();
    return 0;
}