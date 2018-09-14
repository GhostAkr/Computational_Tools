//
// Created by ghostakr on 9/14/18.
//

#include "../include/Matrix.h"

//  Constructors

Matrix::Matrix() {
    data = NULL;
    rows = 0;
    cols = 0;
}

// Destructors

Matrix::~Matrix() {
    this->dataDelete();
}

// Access methods

double** Matrix::matrixGet() {
    return data;
}

void Matrix::matrixSet(double** _data, size_t _rows, size_t _cols) {
    this->dataDelete();
    cols = _cols;
    rows = _rows;
    data = _data;
}

void Matrix::matrixPrint() {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            cout << data[i][j] << " ";
        }
        cout << endl;
    }
}

// Setting methods

void Matrix::readMatrixFromFile(string _pathToFile) {
    std::ifstream fileIn(_pathToFile);
    if (!fileIn) {  // Exception
        cout << "Error while reading file" << endl;
        return;
    }
    // Creating new matrix
    size_t newCols = 0, newRows = 0;
    fileIn >> newRows;
    fileIn >> newCols;
    this->dataInit(newRows, newCols);
    rows = newRows;
    cols = newCols;
    double cellValue = 0.0;
    int i = 0, j = 0;  // Current position in the matrix
    // Reading data from file
    while (fileIn >> cellValue) {
        data[i][j] = cellValue;
        if (++j == cols) {
            i++;
            j = 0;
        }
    }
    fileIn.close();
}

// Private methods

void Matrix::dataDelete() {
    for (int i = 0; i < cols; i++) {
        delete [] data[i];
    }
    delete [] data;
}

void Matrix::dataInit(size_t _rows, size_t _cols) {
    if (rows != 0 && cols != 0) {
        this->dataDelete();
    }
    data = new double* [_rows];
    for (int i = 0; i < _rows; i++) {
        data[i] = new double [_cols];
    }
}
