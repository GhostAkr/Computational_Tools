#include <iostream>
#include "../include/Matrix.h"
#include "../include/LinearSystems.h"

using std::cout;
using std::endl;

int main() {
    std::string path = "../data/D7.TXT";
    Matrix* A = new Matrix;
    A->readLinearSystemFromFile(path);
    delete A;
    return 0;
}