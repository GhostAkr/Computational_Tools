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
    Matrix* Vals = QRDecompositionEigen(HesA);
    cout << "Eigenvalues are " << endl;
    Vals->matrixPrint();
    double res = Rayleigh(A);
    cout << "Rayleigh quotient = " << res << endl;
    for (int i = 0; i < Vals->rowsGet(); ++i) {
        double eigenvalue = Vals->matrixGet()[i][0];
        Matrix* r = Reverse(A, eigenvalue);
        cout << "Eigenvector for eigenvalue = " << eigenvalue << " is" << endl;
        r->matrixPrint();
    }
    return 0;
}