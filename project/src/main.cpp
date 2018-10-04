#include <iostream>
#include "../include/Matrix.h"
#include "../include/LinearSystems.h"

using std::cout;
using std::endl;

int main() {
    string PathToTest = "../data/D5.TXT";
    Matrix* A1 = new Matrix;
    Matrix* A2 = new Matrix;
    Matrix* oldA = new Matrix;
    oldA->readLinearSystemFromFile(PathToTest);
    A1->readLinearSystemFromFile(PathToTest);
    A2->readLinearSystemFromFile(PathToTest);
    cout << "Linear system is " << endl;
    oldA->matrixPrint();
    conditionNumber(A1);
    Matrix* X1 = gaussLinearSolve(A1);
    if (X1 == NULL) {
        delete X1;
        delete A1;
        delete A2;
        delete oldA;
        return 0;
    }
    cout << "Solution with Gauss method is " << endl;
    X1->matrixPrint();
    valuationVector(X1, oldA);
    Matrix* X2 = QRDecompositionSolve(A2);
    cout << "Solution with QR-decomposition method is " << endl;
    X2->matrixPrint();
    valuationVector(X2, oldA);
    delete X1;
    delete X2;
    delete A1;
    delete A2;
    delete oldA;
    return 0;
}