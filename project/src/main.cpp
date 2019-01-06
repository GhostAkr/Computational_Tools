#include <iostream>
#include "../include/Matrix.h"
#include "../include/LinearSystems.h"
#include "../include/Interpolation.h"
#include "../include/EigenValues.h"

using std::cout;
using std::endl;

int main() {
	string path = "C:/Users/User/Desktop/Новая папка (3)/Computational Tools (VS)/Computational Tools/data/EigenTest.txt";
	Matrix* A = new Matrix;
	A->readMatrixFromFile(path);
	A->matrixPrint();
	double res = Rayleigh(A);
	cout << res << endl;
	//Matrix* r = new Matrix;
	//r = Reverse(A);
	//r->matrixPrint();
	system("pause");
    return 0;
}