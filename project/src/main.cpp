#include <iostream>
#include "../include/Matrix.h"
#include "../include/LinearSystems.h"
#include "../include/Interpolation.h"

using std::cout;
using std::endl;

int main() {
    std::string path = "../data/SUPERTEST.txt";
    double leftBorder = -1.0;
    double rightBorder = 1.0;
    int numberOfParts = 10;
    double* uniformMesh = Mesh(leftBorder, rightBorder, numberOfParts);
    double** functionValues = func2(uniformMesh, numberOfParts);
    double* functionMesh = MeshCheb(leftBorder, rightBorder, numberOfParts * 10);
    double** polynom = Polynom(functionValues, functionMesh, numberOfParts, numberOfParts * 10);
    Extracttofile(polynom, numberOfParts * 10, path);
    return 0;
}