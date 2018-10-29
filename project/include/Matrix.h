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

typedef double type;  // For calculations with high accuracy
//typedef float type;  // For calculations with usual accuracy

using std::cout;
using std::endl;
using std::string;

class Matrix {
private:
    // Data
    type** data;
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
    void matrixSet(type** _data, size_t _rows, size_t _cols);
    type** matrixGet() const;
    void matrixPrint() const;
    size_t rowsGet() const;
    size_t colsGet() const;
    void rowsSet(size_t _rows);
    void colsSet(size_t _cols);

    // Setting methods
    void readMatrixFromFile(string _pathToFile);
    void readLinearSystemFromFile(string _pathToFile);

    // Operations methods
    static Matrix* matrixSum(const Matrix* _matrix1, const Matrix* _matrix2);
    static Matrix* matrixComp(const Matrix* _matrix1, const Matrix* _matrix2);
    static Matrix* matrixDiff(const Matrix* _matrix1, const Matrix* _matrix2);
    static Matrix* matrixConstComp(const Matrix* _matrix, const type _const);
    static Matrix* getCopy(const Matrix* _matrix);
    void matrixRowsChange(int firstRow, int secondRow);
    void matrixTranspose();

    // Special matrixes
    void matrixNullSet(size_t _rows, size_t _cols);
    void matrixOneSet(size_t _rows, size_t _cols);

    // Other methods
    int mainElement(int iteration);  // Searches for main element in matrix by lines
};

#endif //LAB_1_MATRIX_H
