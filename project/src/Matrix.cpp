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

Matrix::Matrix(size_t _rows, size_t _cols) {
    rows = _rows;
    cols = _cols;
    dataInit(_rows, _cols);
}

// Destructors

Matrix::~Matrix() {
    this->dataDelete();
}

// Access methods

type** Matrix::matrixGet() const {
    return data;
}

void Matrix::matrixSet(type** _data, size_t _rows, size_t _cols) {
    this->dataDelete();
    cols = _cols;
    rows = _rows;
    data = _data;
}

void Matrix::matrixPrint() const {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            //cout.precision(4);
            //cout << data[i][j] << " ";
			printf("%.16f ", data[i][j]);
        }
        cout << endl;
    }
}

size_t Matrix::rowsGet() const {
    return rows;
}

size_t Matrix::colsGet() const {
    return cols;
}

void Matrix::rowsSet(size_t _rows) {
    rows = _rows;
}

void Matrix::colsSet(size_t _cols) {
    cols = _cols;
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
    type cellValue = 0.0;
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

void Matrix::readLinearSystemFromFile(string _pathToFile) {
    std::ifstream fileIn(_pathToFile);
    if (!fileIn) {  // Exception
        cout << "Error while reading file" << endl;
        return;
    }
    // Creating new matrix
    size_t newCols = 0, newRows = 0;
    fileIn >> newRows;
    newCols = newRows + 1;
    this->dataInit(newRows, newCols);
    rows = newRows;
    cols = newCols;
    type cellValue = 0.0;
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
    for (int i = 0; i < rows; i++) {
        delete [] data[i];
    }
    delete [] data;
}

void Matrix::dataInit(size_t _rows, size_t _cols) {
    if (rows != 0 || cols != 0) {
        this->dataDelete();
    }
    data = new type* [_rows];
    for (int i = 0; i < _rows; i++) {
        data[i] = new type [_cols];
    }
    for (int i = 0; i < _rows; ++i) {
        for (int j = 0; j < _cols; ++j) {
            data[i][j] = 0.0;
        }
    }
}

// Operations methods

type Matrix::norm() {
	if (cols > 1) {
		cout << "It's not a vector" << endl;
		return -1.0;
	}
	type sum = 0.0;
	for (int k = 0; k < rows; ++k) {
		sum += data[k][0] * data[k][0];
	}
	return sqrt(sum);
}

Matrix* Matrix::matrixSum(const Matrix* _matrix1, const Matrix* _matrix2) {
    // Checking for matrix compatibility
    size_t rows1 = 0, cols1 = 0;  // Rows and cols of matrix 1
    size_t rows2 = 0, cols2 = 0;  // Rows and cols of matrix 2
    rows1 = _matrix1->rowsGet();
    cols1 = _matrix1->colsGet();
    rows2 = _matrix2->rowsGet();
    cols2 = _matrix2->colsGet();
    if ((rows1 != rows2) || (cols1 != cols2)) {
        cout << "Matrixes are not compatible" << endl;
        return NULL;
    }
    // Summing
    Matrix* outMatrix = new Matrix;
    type** data1 = _matrix1->matrixGet();
    type** data2 = _matrix2->matrixGet();
    type** newData = new type* [rows1];
    for (int i = 0; i < rows1; ++i) {
        newData[i] = new type [cols1];
    }
    for (int i = 0; i < rows1; ++i) {
        for (int j = 0; j < cols1; ++j) {
            newData[i][j] = data1[i][j] + data2[i][j];
        }

    }
    outMatrix->matrixSet(newData, rows1, cols1);
    return outMatrix;
}

Matrix* Matrix::matrixComp(const Matrix* _matrix1, const Matrix* _matrix2) {
    // TODO: Fix leak in composition
    // Checking for matrix compatibility
    size_t rows1 = 0, cols1 = 0;  // Rows and cols of matrix 1
    size_t rows2 = 0, cols2 = 0;  // Rows and cols of matrix 2
    rows1 = _matrix1->rowsGet();
    cols1 = _matrix1->colsGet();
    rows2 = _matrix2->rowsGet();
    cols2 = _matrix2->colsGet();
    if (cols1 != rows2) {
        cout << "Matrixes are not compatible" << endl;
        return NULL;
    }
    // Composition
    type** data1 = _matrix1->matrixGet();
    type** data2 = _matrix2->matrixGet();
    type** newData = new type* [rows1];
    for (int i = 0; i < rows1; ++i) {
        newData[i] = new type [cols2];
    }
	for (int i = 0; i < rows1; ++i) {
		for (int j = 0; j < cols2; ++j) {
			newData[i][j] = 0.0;
		}
	}
    for (int i = 0; i < rows1; ++i) {
        for (int j = 0; j < cols2; ++j) {
            for (int k = 0; k < cols1; ++k) {
                newData[i][j] += data1[i][k] * data2[k][j];
            }
        }
    }
    Matrix* outMatrix = new Matrix;
    outMatrix->matrixSet(newData, rows1, cols2);
    return outMatrix;
}

Matrix* Matrix::matrixDiff(const Matrix* _matrix1, const Matrix* _matrix2) {
    // Checking for matrix compatibility
    size_t rows1 = 0, cols1 = 0;  // Rows and cols of matrix 1
    size_t rows2 = 0, cols2 = 0;  // Rows and cols of matrix 2
    rows1 = _matrix1->rowsGet();
    cols1 = _matrix1->colsGet();
    rows2 = _matrix2->rowsGet();
    cols2 = _matrix2->colsGet();
    if ((rows1 != rows2) || (cols1 != cols2)) {
        cout << "Matrixes are not compatible" << endl;
        return NULL;
    }
    // Summing
    Matrix* outMatrix = new Matrix;
    type** data1 = _matrix1->matrixGet();
    type** data2 = _matrix2->matrixGet();
    type** newData = new type* [rows1];
    for (int i = 0; i < rows1; ++i) {
        newData[i] = new type [cols1];
    }
    for (int i = 0; i < rows1; ++i) {
        for (int j = 0; j < cols1; ++j) {
            newData[i][j] = data1[i][j] - data2[i][j];
        }

    }
    outMatrix->matrixSet(newData, rows1, cols1);
    return outMatrix;
}

void Matrix::matrixRowsChange(int firstRow, int secondRow) {
    type* buffRow = data[secondRow];
    data[secondRow] = data[firstRow];
    data[firstRow] = buffRow;
}

void Matrix::matrixTranspose() {
    size_t newRows = cols;
    size_t newCols = rows;
    type** newData = new type* [newRows];
    for (int i = 0; i < newRows; ++i) {
        newData[i] = new type [newCols];
    }
    for (int i = 0; i < newRows; ++i) {
        for (int j = 0; j < newCols; ++j) {
            newData[i][j] = data[j][i];
        }
    }
    this->dataDelete();
    data = newData;
    rows = newRows;
    cols = newCols;
}

// Special matrixes

void Matrix::matrixNullSet(size_t _rows, size_t _cols) {
    this->dataInit(_rows, _cols);
    rows = _rows;
    cols = _cols;
}

void Matrix::matrixOneSet(size_t _rows, size_t _cols) {
    if ( _rows != _cols) {
        cout << "It's not possible to make identity matrix with such size" << endl;
        return;
    }
    this->matrixNullSet(_rows, _cols);
    for (int i = 0; i < _rows; ++i) {
        data[i][i] = 1.0;
    }
}



// Other methods

int Matrix::mainElement(int iteration) {
    int indexOfMainElem = iteration;
    type mainElem = data[iteration][iteration];
    for (int k = iteration; k < rows; ++k) {
        if (data[k][iteration] > mainElem) {
            mainElem = data[k][iteration];
            indexOfMainElem = k;
        }
    }
    return indexOfMainElem;
}
