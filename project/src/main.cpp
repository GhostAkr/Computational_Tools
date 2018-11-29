#include <iostream>
#include "../include/Matrix.h"
#include "../include/LinearSystems.h"
#include "../include/Interpolation.h"

using std::cout;
using std::endl;

int main() {
    //// Parameters
    std::string path = "../data/SUPERTEST.txt";
    double leftBorder = -3.0;
    double rightBorder = 3.0;
    int numberOfParts = 100;
    int split = 10;
    int numberOfX = numberOfParts * split - 1;

    //// Choosing of mesh
    //double* uniformMesh = Mesh(leftBorder, rightBorder, numberOfParts);
    double* uniformMesh = MeshCheb(leftBorder, rightBorder, numberOfParts);
    //double* functionMesh = Mesh(leftBorder, rightBorder, numberOfX);
    double* functionMesh = MeshChebX(uniformMesh, numberOfParts, split);

    double** functionValues = func3(uniformMesh, numberOfParts + 1);

    //// Choosing of way to interpolate
    double** interpolant = Spline(functionValues, functionMesh, numberOfParts, numberOfX);
    //double** interpolant = Polynom(functionValues, functionMesh, numberOfParts, numberOfX);

    //// Output
    Print(interpolant, numberOfX);
    Extracttofile(interpolant, numberOfX, path);
    return 0;
}