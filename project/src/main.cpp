#include <iostream>
#include "../include/Matrix.h"
#include "../include/LinearSystems.h"
#include "../include/Interpolation.h"
#include "../include/EigenValues.h"
#include "../include/NonlinearSolve.h"

using std::cout;
using std::endl;

int main() {
//    int numberOfIntervals = 100;
//    double a = 0.0;  // Left border
//    double b = 1.0;  // Right border
//    double *mesh = Mesh(0, 1, numberOfIntervals);
////    cout << "Mesh" << endl;
////    Print(mesh, numberOfIntervals + 1);
//    double **funcMesh = f3(mesh, numberOfIntervals + 1);
//    int *numOfPairs = new int;
//    double **localMesh = rootsLocale(funcMesh, numberOfIntervals + 1, numOfPairs);
////    cout << "Correct points" << endl;
////    Print(localMesh, *numOfPairs);
    auto *numOfRoots = new int;
//    double **roots = Newtonsys
//    cout << "Roots" << endl;
//    Print(roots, *numOfRoots);
    string path = "../data/dat.txt";
    auto** roots = Newtonsys(21, path, numOfRoots);
    for (int i = 0; i < *numOfRoots; ++i) {
        cout << roots[0][i] << " " << roots[1][i] << endl;
    }
    //system("pause");
    return 0;
    return 0;
}