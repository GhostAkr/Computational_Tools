#include <iostream>
#include "../include/Matrix.h"
#include "../include/LinearSystems.h"
#include "../include/Interpolation.h"
#include "../include/EigenValues.h"
#include "../include/NonlinearSolve.h"

using std::cout;
using std::endl;

int main() {
    int numberOfIntervals = 27;
    double a = 0.0;  // Left border
    double b = 2.0;  // Right border
    double *mesh = Mesh(a, b, numberOfIntervals);
    double **funcMesh = f4(mesh, numberOfIntervals + 1);
    auto *numOfPairs = new int;
    *numOfPairs = 0;
    double **localMesh = rootsLocale(funcMesh, numberOfIntervals + 1, numOfPairs);
//    double** localMesh = new double* [2];
//    localMesh[0] = new double [2];
//    localMesh[1] = new double [2];
//    localMesh[0][0] = -1;
//    localMesh[0][1] = 10;
//    localMesh[1][0] = -1;
//    localMesh[1][1] = 2.317;
    //*numOfPairs = 2;
    auto *numOfRoots = new int;
    cout << "Roots with bisection" << endl;
    double *rootsBisec = Bisection(localMesh, *numOfPairs,numOfRoots, ff4);
//    cout << "Error for root #1 is " << rootsBisec[0] - 0.1 << endl;
//    cout << "Error for root #2 is " << rootsBisec[1] - 0.22 << endl;
//    cout << "Error for root #3 is " << rootsBisec[2] - 0.55 << endl;
//    cout << "Error for root #4 is " << rootsBisec[3] - 0.7 << endl;
//    cout << "Error for root #5 is " << rootsBisec[4] - 0.75 << endl;

    //Print(rootsBisec, *numOfRoots);
    string path = "../data/grid.dat";
    *numOfRoots = 0;
    cout << "Roots with newton" << endl;
    auto* rootsNewt = Newton(localMesh, *numOfPairs, numOfRoots, ff4);
//    cout << "Error for root #1 is " << rootsNewt[0] - 0.1 << endl;
//    cout << "Error for root #2 is " << rootsNewt[1] - 0.22 << endl;
//    cout << "Error for root #3 is " << rootsNewt[2] - 0.55 << endl;
//    cout << "Error for root #4 is " << rootsNewt[3] - 0.7 << endl;
//    cout << "Error for root #5 is " << rootsNewt[4] - 0.75 << endl;

    //Print(rootsNewt, *numOfRoots);
    *numOfRoots = 0;
    auto** rootsSys = Newtonsys(41, path, numOfRoots);
    cout << "Root in system" << endl;
    for (int i = 0; i < *numOfRoots; ++i) {
        cout << rootsSys[0][i] << " " << rootsSys[1][i] << endl;
    }
    return 0;
}