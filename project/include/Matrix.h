//
// Created by ghostakr on 9/14/18.
//

#ifndef LAB_1_MATRIX_H
#define LAB_1_MATRIX_H

#include <iostream>
#include <string>
#include <fstream>

#define TEST cout << "Test" << endl;  // Macros for debugging

using std::cout;
using std::endl;
using std::string;

class Matrix {
private:
    // Data
    double** data;
    size_t rows, cols;

    // Private methods
    void dataDelete();  // Clearing data pointer
    void dataInit(size_t _rows, size_t _cols);  // Initiating pointer

public:
    // Constructors
    Matrix();

    // Destructors
    ~Matrix();

    // Access methods
    void matrixSet(double** _data, size_t _rows, size_t _cols);
    double** matrixGet();
    void matrixPrint();

    // Setting methods
    void readMatrixFromFile(string _pathToFile);
};


#endif //LAB_1_MATRIX_H
