#include <iostream>
#include <stdexcept>
#include <cmath>
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
            case '1':
                // Print a Matrix
                // Implementation for printing a matrix is provided in case 1
                break;
            case '2':
                // Print a Vector
                // Implementation for printing a vector is provided in case 2
                break;
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
void PerformMatrixOperations() {
    char choice;
    do {
        cout << "\nMatrix Operations:\n";
        cout << "1. Inverse of a Matrix\n";
        cout << "2. Addition of Matrices\n";
        cout << "3. Scalar Addition and Subtraction of a Matrix\n";
        cout << "4. Return to Main Menu\n";
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
                    int numMatrices;
                    cout << "Enter the number of matrices to add (2 or 3): ";
                    cin >> numMatrices;
                    if (cin.fail() || (numMatrices != 2 && numMatrices != 3)) {
                        throw invalid_argument("Invalid input for number of matrices.");
                    }
                    int rows, cols;
                    cout << "Enter matrix dimensions (rows cols): ";
                    cin >> rows >> cols;
                    if (cin.fail()) {
                        throw invalid_argument("Invalid input for dimensions.");
                    }

                    // Vector to store matrices
                    vector<Matrix2<double>> matrices;
                    matrices.reserve(numMatrices);

                    // Input matrices
                    for (int i = 0; i < numMatrices; ++i) {
                        cout << "Enter elements for Matrix " << i + 1 << ":\n";
                        double* data = new double[rows * cols];
                        for (int j = 0; j < rows * cols; ++j) {
                            cin >> data[j];
                            if (cin.fail()) {
                                delete[] data;
                                throw invalid_argument("Invalid input for matrix element.");
                            }
                        }
                        matrices.emplace_back(rows, cols, data);
                        delete[] data;
                    }

                    // Perform matrix addition
                    Matrix2<double> result(rows, cols);
                    if (numMatrices == 2) {
                        result = addMatrices(matrices[0], matrices[1]);
                    } else { // numMatrices == 3
                        result = addMatrices(matrices[0], matrices[1], matrices[2]);
                    }

                    cout << "Result Matrix:\n";
                    Print(result);

                } catch (const invalid_argument& e) {
                    cerr << "Error: " << e.what() << endl;
                    cin.clear();
                    cin.ignore(numeric_limits<streamsize>::max(), '\n');
                }
                break;
            }
            case '3':
                // Scalar Addition and Subtraction of a Matrix
                try {
                    // Code to perform scalar addition and subtraction of a matrix
                } catch (const invalid_argument& e) {
                    cerr << "Error: " << e.what() << endl;
                    cin.clear();
                    cin.ignore(numeric_limits<streamsize>::max(), '\n');
                }
                break;
            case '4':
                // Return to Main Menu
                break;
            default:
                cout << "Invalid choice! Please try again.\n";
        }
    } while (choice != '4');
}

// Function to perform vector operations
// void PerformVectorOperations() {
//     char choice;
//     do {
//         cout << "\nVector Operations:\n";
//         cout << "1. Addition of Vectors\n";
//         cout << "2. Dot Product of Vectors\n";
//         cout << "3. Scalar Multiplication of a Vector\n";
//         cout << "4. Cross Product of 3D Vectors\n";
//         cout << "5. Return to Main Menu\n";
//         cout << "Enter your choice: ";
//         cin >> choice;

//         switch (choice) {
//             case '1':
//                 // Addition of Vectors
//                 try {
//                     // Code to perform vector addition
//                 } catch (const invalid_argument& e) {
//                     cerr << "Error: " << e.what() << endl;
//                     cin.clear();
//                     cin.ignore(numeric_limits<streamsize>::max(), '\n');
//                 }
//                 break;
//             case '2':
//                 // Dot Product of Vectors
//                 try {
//                     int numVectors;
//                     cout << "Enter the number of vectors to compute dot product: ";
//                     cin >> numVectors;
//                     if (cin.fail() || numVectors < 2) {
//                         throw invalid_argument("Invalid input for number of vectors.");
//                     }
//                     int dims;
//                     cout << "Enter vector dimensions: ";
//                     cin >> dims;
//                     if (cin.fail()) {
//                         throw invalid_argument("Invalid input for dimensions.");
//                     }
//                     Vector<double> result;
//                     for (int i = 0; i < numVectors; ++i) {
//                         cout << "Enter elements for Vector " << i + 1 << ":\n";
//                         vector<double> vecData(dims);
//                         for (int j = 0; j < dims; ++j) {
//                             cin >> vecData[j];
//                             if (cin.fail()) {
//                                 throw invalid_argument("Invalid input for vector element.");
//                             }
//                         }
//                         Vector<double> vec(vecData);
//                         result = (i == 0) ? vec : Vector<double>::dot(result, vec);
//                     }
//                     cout << "Dot Product: " << result << endl;
//                 } catch (const invalid_argument& e) {
//                     cerr << "Error: " << e.what() << endl;
//                     cin.clear();
//                     cin.ignore(numeric_limits<streamsize>::max(), '\n');
//                 }
//                 break;
//             case '3':
//                 // Scalar Multiplication of a Vector
//                 try {
//                     // Code to perform scalar multiplication of a vector
//                 } catch (const invalid_argument& e) {
//                     cerr << "Error: " << e.what() << endl;
//                     cin.clear();
//                     cin.ignore(numeric_limits<streamsize>::max(), '\n');
//                 }
//                 break;
//             case '4':
//                 // Cross Product of 3D Vectors
//                 try {
//                     // Code to compute and print cross product of 3D vectors
//                 } catch (const invalid_argument& e) {
//                     cerr << "Error: " << e.what() << endl;
//                     cin.clear();
//                     cin.ignore(numeric_limits<streamsize>::max(), '\n');
//                 }
//                 break;
//             case '5':
//                 // Return to Main Menu
//                 break;
//             default:
//                 cout << "Invalid choice! Please try again.\n";
//         }
//     } while (choice != '5');
// }