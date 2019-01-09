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
    double* mesh = Mesh(0, 1, numberOfIntervals);
    cout << "Mesh" << endl;
    Print(mesh, numberOfIntervals + 1);
    double** funcMesh = f1(mesh, numberOfIntervals + 1);
    int* numOfPairs = new int;
    double** localMesh = rootsLocale(funcMesh, numberOfIntervals + 1, numOfPairs);
    cout << "Correct points" << endl;
    Print(localMesh, *numOfPairs);
    return 0;
}