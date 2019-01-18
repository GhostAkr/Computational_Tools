#include <iostream>
#include "../include/Matrix.h"
#include "../include/LinearSystems.h"
#include "../include/Interpolation.h"

using std::cout;
using std::endl;

int main() {
    //// Parameters
    std::string path = "../data/SUPERTEST.txt";
    double leftBorder = -5.0;
    double rightBorder = 5.0;
    int numberOfParts = 1280;
    int split = 2;
    int numberOfX = numberOfParts * split;

    //// Choosing of mesh
    double* uniformMesh = Mesh(leftBorder, rightBorder, numberOfParts);
    //double* uniformMesh = MeshCheb(leftBorder, rightBorder, numberOfParts);

    double* functionMesh = Mesh(leftBorder, rightBorder, numberOfX);
    //double* functionMesh = MeshChebX(uniformMesh, numberOfParts, split);

    double** functionValues = func2(uniformMesh, numberOfParts + 1);

    //// Choosing of way to interpolate
    double** interpolant = Spline(functionValues, functionMesh, numberOfParts, numberOfX);
    //double** interpolant = Polynom(functionValues,functionMesh,numberOfParts,numberOfX);

    cout << "Error = " << pogr(interpolant, numberOfX) << endl;

    //// Output
    Extracttofile(interpolant, numberOfX + 1, path);
    return 0;
}