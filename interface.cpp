#include <iostream>
#include <stdexcept>
#include <cmath>
#include <limits>
#include "matrix.h"
#include "vector.h"

using namespace std;

// Function prototypes
template <class T>
void Print(Matrix2<T> matrix);

template <class T>
void Print(Vector<T> vec);

void PerformMatrixOperations();
void PerformVectorOperations();

int main() {
    char choice;
    do {
        cout << "\nOptions:\n";
        cout << "1. Print a Matrix\n";
        cout << "2. Print a Vector\n";
        cout << "3. Perform Matrix Operations\n";
        cout << "4. Perform Vector Operations\n";
        cout << "5. Exit\n";
        cout << "Enter your choice: ";
        cin >> choice;

        switch (choice) {
            case '1': {
                 try {
                    cout << "Enter matrix dimensions (rows cols): ";
                    int rows, cols;
                    cin >> rows >> cols;
                    if (cin.fail()) {
                        throw invalid_argument("Invalid input for dimensions.");
                    }
                    cout << "Enter matrix elements:\n";
                    double* data = new double[rows * cols];
                    for (int i = 0; i < rows * cols; ++i) {
                        cin >> data[i];
                        if (cin.fail()) {
                            delete[] data;
                            throw invalid_argument("Invalid input for matrix element.");
                        }
                    }
                    Matrix2<double> matrix(rows, cols, data);
                    cout << "Matrix:\n";
                    Print(matrix);
                    delete[] data;
                } catch (const invalid_argument& e) {
                    cerr << "Error: " << e.what() << endl;
                    cin.clear();
                    cin.ignore(numeric_limits<streamsize>::max(), '\n');
                }
                break;
            }
            case '2': {
                try {
                    cout << "Enter vector dimensions: ";
                    int dims;
                    cin >> dims;
                    if (cin.fail()) {
                        throw invalid_argument("Invalid input for dimensions.");
                    }
                    cout << "Enter vector elements:\n";
                    vector<double> vecData(dims);
                    for (int i = 0; i < dims; ++i) {
                        cin >> vecData[i];
                        if (cin.fail()) {
                            throw invalid_argument("Invalid input for vector element.");
                        }
                    }
                    Vector<double> vec(vecData);
                    cout << "Vector:\n";
                    Print(vec);
                } catch (const invalid_argument& e) {
                    cerr << "Error: " << e.what() << endl;
                    cin.clear();
                    cin.ignore(numeric_limits<streamsize>::max(), '\n');
                }
                break;
            }
            case '3':
                // Perform Matrix Operations
                PerformMatrixOperations();
                break;
            case '4':
                // Perform Vector Operations
                PerformVectorOperations();
                break;
            case '5':
                cout << "Exiting...\n";
                break;
            default:
                cout << "Invalid choice! Please try again.\n";
        }
    } while (choice != '5');

    return 0;
}

// Function to print a matrix
template <class T>
void Print(Matrix2<T> matrix) {
    int nRows = matrix.GetNumRows();
    int nCols = matrix.GetNumCols();
    for (int row = 0; row < nRows; ++row) {
        for (int col = 0; col < nCols; ++col) {
            cout << matrix.GetElement(row, col) << " ";
        }
        cout << endl;
    }
    cout << endl;
}

// Function to print a vector
template <class T>
void Print(Vector<T> vec) {
    int nDims = vec.GetNumDims();
    for (int i = 0; i < nDims; ++i) {
        cout << vec.GetElement(i) << " ";
    }
    cout << endl << endl;
}

// Function to perform matrix operations
void PerformMatrixOperations() {
    char choice;
    do {
        cout << "\nMatrix Operations:\n";
        cout << "1. Inverse of a Matrix\n";
        cout << "2. Addition of Matrices\n";
        cout << "3. Scalar Addition and Subtraction of a Matrix\n";
        cout << "4. Multiplication of Matrices\n";
        cout << "5. Return to Main Menu\n";
        cout << "Enter your choice: ";
        cin >> choice;

        switch (choice) {
            case '1':
                // Inverse of a Matrix
                try {
                    int rows, cols;
                    cout << "Enter dimensions of the matrix (rows cols): ";
                    cin >> rows >> cols;
                    if (cin.fail() || rows <= 0 || cols <= 0) {
                        throw invalid_argument("Invalid input for dimensions.");
                    }
                    cout << "Enter matrix elements:\n";
                    double* data = new double[rows * cols];
                    for (int i = 0; i < rows * cols; ++i) {
                        cin >> data[i];
                        if (cin.fail()) {
                            delete[] data;
                            throw invalid_argument("Invalid input for matrix element.");
                        }
                    }
                    Matrix2<double> matrix(rows, cols, data);
                    cout << "Original Matrix:\n";
                    Print(matrix);
                    Matrix2<double> inverseMatrix = matrix.inverse();
                    cout << "Inverse of the Matrix:\n";
                    Print(inverseMatrix);
                    delete[] data;
                } catch (const invalid_argument& e) {
                    cerr << "Error: " << e.what() << endl;
                    cin.clear();
                    cin.ignore(numeric_limits<streamsize>::max(), '\n');
                }
                break;
            case '2':
                // Addition of Matrices
                try {
                    int rows, cols;
                    cout << "Enter dimensions of matrices (rows cols): ";
                    cin >> rows >> cols;
                    if (cin.fail() || rows <= 0 || cols <= 0) {
                        throw invalid_argument("Invalid input for dimensions.");
                    }
                    cout << "Enter the number of matrices to add (2 or 3): ";
                    int numMatrices;
                    cin >> numMatrices;
                    if (cin.fail() || (numMatrices != 2 && numMatrices != 3)) {
                        throw invalid_argument("Invalid input for the number of matrices.");
                    }
                    vector<Matrix2<double>> matrices;
                    for (int i = 0; i < numMatrices; ++i) {
                        cout << "Enter elements for matrix " << i + 1 << ":\n";
                        double* data = new double[rows * cols];
                        for (int j = 0; j < rows * cols; ++j) {
                            cin >> data[j];
                            if (cin.fail()) {
                                delete[] data;
                                throw invalid_argument("Invalid input for matrix element.");
                            }
                        }
                        matrices.push_back(Matrix2<double>(rows, cols, data));
                        delete[] data;
                    }
                    if (numMatrices == 2) {
                        Matrix2<double> result = addMatrices(matrices[0], matrices[1]);
                        cout << "Result of addition:\n";
                        Print(result);
                    } else {
                        Matrix2<double> result = addMatrices(matrices[0], matrices[1], matrices[2]);
                        cout << "Result of addition:\n";
                        Print(result);
                    }
                } catch (const invalid_argument& e) {
                    cerr << "Error: " << e.what() << endl;
                    cin.clear();
                    cin.ignore(numeric_limits<streamsize>::max(), '\n');
                }
                break;
            case '3':
                try {
                    int rows, cols;
                    cout << "Enter matrix dimensions (rows cols): ";
                    cin >> rows >> cols;
                    if (cin.fail() || rows <= 0 || cols <= 0) {
                        throw invalid_argument("Invalid input for dimensions.");
                    }
                    cout << "Enter scalar value: ";
                    double scalar;
                    cin >> scalar;
                    if (cin.fail()) {
                        throw invalid_argument("Invalid input for scalar.");
                    }
                    cout << "Enter matrix elements:\n";
                    double* data = new double[rows * cols];
                    for (int i = 0; i < rows * cols; ++i) {
                        cin >> data[i];
                        if (cin.fail()) {
                            delete[] data;
                            throw invalid_argument("Invalid input for matrix element.");
                        }
                    }
                    Matrix2<double> matrix(rows, cols, data);
                    cout << "Original Matrix:\n";
                    Print(matrix);
                    Matrix2<double> result1 = matrix + scalar;
                    cout << "Matrix after scalar addition:\n";
                    Print(result1);
                    Matrix2<double> result2 = matrix - scalar;
                    cout << "Matrix after scalar subtraction:\n";
                    Print(result2);
                    delete[] data;
                } catch (const invalid_argument& e) {
                    cerr << "Error: " << e.what() << endl;
                    cin.clear();
                    cin.ignore(numeric_limits<streamsize>::max(), '\n');
                }
                break;
                
            case '4':
                            
                // Matrix Multiplication
                try {
                    int rows1, cols1, rows2, cols2;
                    cout << "Enter dimensions of the first matrix (rows cols): ";
                    cin >> rows1 >> cols1;
                    if (cin.fail() || rows1 <= 0 || cols1 <= 0) {
                        throw invalid_argument("Invalid input for dimensions of the first matrix.");
                    }

                    cout << "Enter elements for the first matrix:\n";
                    double* data1 = new double[rows1 * cols1];
                    for (int i = 0; i < rows1 * cols1; ++i) {
                        cin >> data1[i];
                        if (cin.fail()) {
                            delete[] data1;
                            throw invalid_argument("Invalid input for matrix element.");
                        }
                    }
                    Matrix2<double> matrix1(rows1, cols1, data1);

                    cout << "Enter dimensions of the second matrix (rows cols): ";
                    cin >> rows2 >> cols2;
                    if (cin.fail() || rows2 <= 0 || cols2 <= 0 || cols1 != rows2) {
                        delete[] data1;
                        throw invalid_argument("Invalid input for dimensions of the second matrix or incompatible dimensions for multiplication.");
                    }

                    cout << "Enter elements for the second matrix:\n";
                    double* data2 = new double[rows2 * cols2];
                    for (int i = 0; i < rows2 * cols2; ++i) {
                        cin >> data2[i];
                        if (cin.fail()) {
                            delete[] data1;
                            delete[] data2;
                            throw invalid_argument("Invalid input for matrix element.");
                        }
                    }
                    Matrix2<double> matrix2(rows2, cols2, data2);

                    cout << "First Matrix:\n";
                    Print(matrix1);
                    cout << "Second Matrix:\n";
                    Print(matrix2);

                    Matrix2<double> result = matrix1 * matrix2;
                    cout << "Result of matrix multiplication:\n";
                    Print(result);

                    delete[] data1;
                    delete[] data2;
                } catch (const invalid_argument& e) {
                    cerr << "Error: " << e.what() << endl;
                    cin.clear();
                    cin.ignore(numeric_limits<streamsize>::max(), '\n');
                }
                break;
            default:
                cout << "Invalid choice! Please try again.\n";
        }
    } while (choice != '5');
}

    
// Function to perform vector operations
// Function to perform vector operations
void PerformVectorOperations() {
    char choice;
    do {
        cout << "\nVector Operations:\n";
        cout << "1. Addition of Vectors\n";
        cout << "2. Dot Product of Vectors\n";
        cout << "3. Scalar Multiplication of a Vector\n";
        cout << "4. Cross Product of 3D Vectors\n";
        cout << "5. Return to Main Menu\n";
        cout << "Enter your choice: ";
        cin >> choice;

        switch (choice) {
            case '1':
                // Addition of Vectors
                try {
                    int dims;
                    cout << "Enter the number of dimensions: ";
                    cin >> dims;
                    if (cin.fail() || dims <= 0) {
                        throw invalid_argument("Invalid input for the number of dimensions.");
                    }

                    cout << "Enter the number of vectors to add (2 or more): ";
                    int numVectors;
                    cin >> numVectors;
                    if (cin.fail() || numVectors < 2) {
                        throw invalid_argument("Invalid input for the number of vectors.");
                    }

                    vector<Vector<double>> vectors;
                    for (int i = 0; i < numVectors; ++i) {
                        cout << "Enter elements for vector " << i + 1 << ":\n";
                        vector<double> vecData(dims);
                        for (int j = 0; j < dims; ++j) {
                            cin >> vecData[j];
                            if (cin.fail()) {
                                throw invalid_argument("Invalid input for vector element.");
                            }
                        }
                        vectors.push_back(Vector<double>(vecData));
                    }

                    // Compute the result of vector addition
                    Vector<double> result = vectors[0];
                    for (int i = 1; i < numVectors; ++i) {
                        result = result + vectors[i];
                    }

                    cout << "Result of vector addition:\n";
                    Print(result);

                } catch (const invalid_argument& e) {
                    cerr << "Error: " << e.what() << endl;
                    cin.clear();
                    cin.ignore(numeric_limits<streamsize>::max(), '\n');
                }
                break;
            case '2':
                // Dot Product of Vectors
                try {
                    int dims;
                    cout << "Enter the number of dimensions: ";
                    cin >> dims;
                    if (cin.fail() || dims <= 0) {
                        throw invalid_argument("Invalid input for the number of dimensions.");
                    }

                    cout << "Enter the number of vectors:(2 or 3) ";
                    int numVectors;
                    cin >> numVectors;
                    if (cin.fail() || numVectors < 2 || numVectors > 3) {
                        throw invalid_argument("Invalid input for the number of vectors.");
                    }

                    vector<Vector<double>> vectors;
                    for (int i = 0; i < numVectors; ++i) {
                        cout << "Enter elements for vector " << i + 1 << ":\n";
                        vector<double> vecData(dims);
                        for (int j = 0; j < dims; ++j) {
                            cin >> vecData[j];
                            if (cin.fail()) {
                                throw invalid_argument("Invalid input for vector element.");
                            }
                        }
                        vectors.push_back(Vector<double>(vecData));
                       
                    }
                    double result;
                    if(numVectors==2)
                        result = dot(vectors[0], vectors[1]); // Assuming at least two vectors are provided
                    else
                        result = dot(vectors[0], vectors[1], vectors[2]); 
                    cout << "Dot product of vectors: " << result << endl;
                    
                } catch (const invalid_argument& e) {
                    cerr << "Error: " << e.what() << endl;
                    cin.clear();
                    cin.ignore(numeric_limits<streamsize>::max(), '\n');
                }
                break;
                        
            case '3':
                // Scalar Multiplication of a Vector
                try {
                    int dims;
                    cout << "Enter the number of dimensions: ";
                    cin >> dims;
                    if (cin.fail() || dims <= 0) {
                        throw invalid_argument("Invalid input for the number of dimensions.");
                    }

                    cout << "Enter scalar value: ";
                    double scalar;
                    cin >> scalar;
                    if (cin.fail()) {
                        throw invalid_argument("Invalid input for scalar.");
                    }

                    cout << "Enter elements for the vector:\n";
                    vector<double> vecData(dims);
                    for (int i = 0; i < dims; ++i) {
                        cin >> vecData[i];
                        if (cin.fail()) {
                            throw invalid_argument("Invalid input for vector element.");
                        }
                    }
                    Vector<double> vec(vecData);

                    // Compute the result of scalar multiplication
                    Vector<double> result = vec * scalar;

                    cout << "Result of scalar multiplication:\n";
                    Print(result);

                } catch (const invalid_argument& e) {
                    cerr << "Error: " << e.what() << endl;
                    cin.clear();
                    cin.ignore(numeric_limits<streamsize>::max(), '\n');
                }
                break;
            case '4':
                // Cross Product of 3D Vectors
                try {
                    cout << "Enter elements for the first 3D vector:\n";
                    vector<double> vecData1(3);
                    for (int i = 0; i < 3; ++i) {
                        cin >> vecData1[i];
                        if (cin.fail()) {
                            throw invalid_argument("Invalid input for vector element.");
                        }
                    }
                    Vector<double> vec1(vecData1);

                    cout << "Enter elements for the second 3D vector:\n";
                    vector<double> vecData2(3);
                    for (int i = 0; i < 3; ++i) {
                        cin >> vecData2[i];
                        if (cin.fail()) {
                            throw invalid_argument("Invalid input for vector element.");
                        }
                    }
                    Vector<double> vec2(vecData2);

                    // Compute the cross product
                    Vector<double> result = Vector<double>::cross(vec1, vec2);

                    cout << "Result of cross product:\n";
                    Print(result);

                } catch (const invalid_argument& e) {
                    cerr << "Error: " << e.what() << endl;
                    cin.clear();
                    cin.ignore(numeric_limits<streamsize>::max(), '\n');
                }
                break;
        }
    } while (choice != '5');
}
