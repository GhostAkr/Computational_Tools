#include <iostream>
#include "../include/Matrix.h"
#include "../include/LinearSystems.h"
#include "../include/Interpolation.h"
#include "../include/EigenValues.h"
#include "../include/NonlinearSolve.h"

using std::cout;
using std::endl;

int main() {
    int numberOfIntervals = 100;
    double a = 0.0;  // Left border
    double b = 1.0;  // Right border
    double *mesh = Mesh(0, 1, numberOfIntervals);
    double **funcMesh = f1(mesh, numberOfIntervals + 1);
    auto *numOfPairs = new int;
    *numOfPairs = 0;
    double **localMesh = rootsLocale(funcMesh, numberOfIntervals + 1, numOfPairs);
    auto *numOfRoots = new int;
    double *rootsBisec = Bisection(localMesh, *numOfPairs,numOfRoots, ff1);
    cout << "Roots with bisection" << endl;
    Print(rootsBisec, *numOfRoots);
    string path = "../data/dat.txt";
    *numOfRoots = 0;
    auto* rootsNewt = Newton(localMesh, *numOfPairs, numOfRoots, ff1);
    cout << "Roots with newton" << endl;
    Print(rootsNewt, *numOfRoots);
    *numOfRoots = 0;
    auto** rootsSys = Newtonsys(21, path, numOfRoots);
    cout << "Root in system" << endl;
    for (int i = 0; i < *numOfRoots; ++i) {
        cout << rootsSys[0][i] << " " << rootsSys[1][i] << endl;
    }
    return 0;
}