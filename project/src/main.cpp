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
    cout << "Eigenvalues are " << endl;
    Result->matrixPrint();
    double res = Rayleigh(A);
    cout << "Rayleigh quotient = " << res << endl;
    Matrix* r = Reverse(A);
    cout << "Eigenvectors are " << endl;
    r->matrixPrint();
    return 0;
}