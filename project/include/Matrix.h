//
// Created by ghostakr on 9/14/18.
//

#ifndef LAB_1_MATRIX_H
#define LAB_1_MATRIX_H

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>

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
    Matrix(size_t _rows, size_t _cols);  // Creating of null-matrix
    // TODO: Repair null constructor

    // Destructors
    ~Matrix();

    // Access methods
    void matrixSet(double** _data, size_t _rows, size_t _cols);
    void matrixNullSet(size_t _rows, size_t _cols);
    double** matrixGet() const;
    void matrixPrint() const;
    size_t rowsGet() const;
    size_t colsGet() const;
    void rowsSet(size_t _rows);
    void colsSet(size_t _cols);

    // Setting methods
    void readMatrixFromFile(string _pathToFile);

    // Operations methods
    static Matrix* matrixSum(const Matrix* _matrix1, const Matrix* _matrix2);
    static Matrix* matrixComp(const Matrix* _matrix1, const Matrix* _matrix2);
    static Matrix* matrixDiff(const Matrix* _matrix1, const Matrix* _matrix2);
    void matrixRowsChange(int firstRow, int secondRow);

    // Other methods
    int mainElement(int iteration);  // Searches for main element in matrix by lines

};

#endif //LAB_1_MATRIX_H
