#include <iostream>
#include "../include/Matrix.h"
#include "../include/LinearSystems.h"
#include "../include/Interpolation.h"
#include "../include/EigenValues.h"

//#include "../src/EigenvaluesProblem.cpp"

using std::cout;
using std::endl;

int main() {
    std::string path = "../data/EigenTest.txt";
    Matrix* A = new Matrix;
    A->readMatrixFromFile(path);
    cout << "Matrix A is " << endl;
    A->matrixPrint();
    Matrix* HesA = HessenbergForm(A);
    cout << "Hessenberg is" << endl;
    HesA->matrixPrint();
    Matrix* Result = QRDecompositionEigen(HesA);
    cout << "Eigen values are " << endl;
    Result->matrixPrint();
    return 0;
}