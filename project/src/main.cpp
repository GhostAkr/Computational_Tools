#include <iostream>
#include "../include/Matrix.h"
#include "../include/LinearSystems.h"
#include "../include/Interpolation.h"

using std::cout;
using std::endl;

int main() {
    //// Parameters
    std::string path = "../data/SUPERTEST.txt";
    double leftBorder = -1.0;
    double rightBorder = 1.0;
    int numberOfParts = 2;
    int split = 2;
    int numberOfX = numberOfParts * split;

    //// Choosing of mesh
    double* uniformMesh = Mesh(leftBorder, rightBorder, numberOfParts);
    //double* uniformMesh = MeshCheb(leftBorder, rightBorder, numberOfParts);
    double* functionMesh = Mesh(leftBorder, rightBorder, numberOfX);
    //double* functionMesh = MeshChebX(uniformMesh, numberOfParts, split);
    cout << "Bulding mesh" << endl;
    Print(uniformMesh, numberOfParts + 1);
    cout << "Calculating mesh" << endl;
    Print(functionMesh, numberOfX + 1);
    double** functionValues = func1(uniformMesh, numberOfParts + 1);
    //double* qw = new double[1];
    //qw[0] = 2.2;
    //// Choosing of way to interpolate
    double** interpolant = Spline(functionValues, functionMesh, numberOfParts, numberOfX);
    //double** interpolant = Polynom(functionValues, qw, numberOfParts, 1);
    //double** interpolant = Polynom(functionValues,functionMesh,numberOfParts,numberOfX);
    //cout<<printf("%.9f",interpolant[1][0] - exp(2.2))<<endl;
    //// Output
    //Print(interpolant, numberOfX);
    Extracttofile(interpolant, numberOfX + 1, path);
    return 0;
}