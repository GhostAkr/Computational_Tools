#include <iostream>
#include "../include/Matrix.h"
#include "../include/LinearSystems.h"
#include "../include/Interpolation.h"
#include "../include/EigenValues.h"

//#include "../src/EigenvaluesProblem.cpp"

using std::cout;
using std::endl;

int main() {
    // Matrix
    std::string path = "../data/EigenTest.txt";
    Matrix* A = new Matrix;
    A->readMatrixFromFile(path);
    cout << "Matrix A is " << endl;
    A->matrixPrint();
    Matrix* HesA = HessenbergForm(A);
    cout << "Hessenberg form of A is" << endl;
    HesA->matrixPrint();
    cout << endl;
    // QR method
    Matrix* Vals = QRDecompositionEigen(HesA);
    cout << "Eigenvalues are (QR method)" << endl;
    Vals->matrixPrint();
    cout << endl;
    // Eigenvalues check
    eigenvalueCheck(HesA, Vals);
    cout << endl;
    // Reverse method
    cout << "Eigenvectors are (Reverse method)" << endl;
    for (int i = 0; i < Vals->rowsGet(); ++i) {
        double eigenvalue = Vals->matrixGet()[i][0];
        Matrix* r = Reverse(A, eigenvalue);
        cout << "Eigenvector for eigenvalue = " << eigenvalue << " is" << endl;
        r->matrixPrint();
        cout << "Check..." << endl;
        eigenvectorCheck(A, r, eigenvalue);
        delete r;
    }
    cout << endl;
    // Start vectors
    string path1 = "../data/v1.txt";
    string path2 = "../data/v2.txt";
    string path3 = "../data/v3.txt";
    string path4 = "../data/v4.txt";
    Matrix* v1 = new Matrix;
    Matrix* v2 = new Matrix;
    Matrix* v3 = new Matrix;
    Matrix* v4 = new Matrix;
    v1->readMatrixFromFile(path1);
    v2->readMatrixFromFile(path2);
    v3->readMatrixFromFile(path3);
    v4->readMatrixFromFile(path4);
    // Rayleigh
    cout << "Rayleigh method" << endl;
    Rayleigh(A, v1);
    Rayleigh(A, v2);
    Rayleigh(A, v3);
    Rayleigh(A, v4);
    return 0;
}