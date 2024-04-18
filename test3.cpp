#include <iostream>
#include <cmath>
#include "matrix.h"
#include "vector.h"

using namespace std;

// A simple function to print a matrix to stdout.
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

// Function to print a vector to stdout.
template <class T>
void Print(Vector<T> vec) {
    int nDims = vec.GetNumDims();
    for (int i = 0; i < nDims; ++i) {
        cout << vec.GetElement(i) << " ";
    }
    cout << endl;
}

int main() {
    // Test Matrix operations
    std::cout << "Testing Matrix operations...\n";

    // Create a 2x2 matrix
    Matrix2<double> matrix(2, 2, new double[4]{1, 2, 3, 4});

    // Print the original matrix
    std::cout << "Original Matrix:" << std::endl;
    Print(matrix);

    // Test matrix inversion
    Matrix2<double> inverseMatrix = matrix.inverse();

    // Print the inverse matrix
    std::cout << "\nInverse Matrix:" << std::endl;
    Print(inverseMatrix);

    // Test Matrix addition and subtraction
    Matrix2<double> matrix2(2, 2, new double[4]{1, 1, 1, 1});
    Matrix2<double> additionResult = matrix + matrix2;
    Matrix2<double> subtractionResult = matrix - matrix2;

    std::cout << "\nMatrix addition result:" << std::endl;
    Print(additionResult);

    std::cout << "\nMatrix subtraction result:" << std::endl;
    Print(subtractionResult);

    // Test Matrix scalar addition and subtraction
    Matrix2<double> scalarAdditionResult = matrix + 2.0;
    Matrix2<double> scalarSubtractionResult = matrix - 2.0;

    std::cout << "\nMatrix scalar addition result:" << std::endl;
    Print(scalarAdditionResult);

    std::cout << "\nMatrix scalar subtraction result:" << std::endl;
    Print(scalarSubtractionResult);

    // Test == operator
    bool res = matrix==matrix2;
    std::cout << "== operator " << res <<endl;
    // Test Matrix multiplication
    std::cout << "\nTesting Matrix multiplication...\n";
    Matrix2<double> matrix3(2, 2, new double[4]{1, 2, 3, 4});
    Matrix2<double> multiplicationResult = matrix * matrix3;

    std::cout << "\nMatrix multiplication result:" << std::endl;
    Print(multiplicationResult);

    // Test Vector operations
    std::cout << "\nTesting Vector operations...\n";

    // Test Vector initialization
    std::cout << "\nTesting Vector initialization...\n";
    std::vector<double> vecData = {1.0, 2.0, 3.0};
    Vector<double> vec(vecData);

    // Test Vector dimension retrieval
    std::cout << "Vector size: " << vec.GetNumDims() << std::endl;

    // Test Vector element retrieval
    std::cout << "Element at index 1: " << vec.GetElement(1) << std::endl;

    // Test Vector addition
    std::cout << "\nTesting Vector addition...\n";
    Vector<double> vec2(vecData);
    Vector<double> vecAdditionResult = vec + vec2;
    std::cout << "Vector addition result:" << std::endl;
    Print(vecAdditionResult);

    // Test Vector dot product
    std::cout << "Dot product of vec and vec2: " << Vector<double>::dot(vec, vec2) << std::endl;

    // Test Vector scalar multiplication
    std::cout << "\nTesting Vector scalar multiplication...\n";
    Vector<double> scalarMultResult = vec * 2.0;
    std::cout << "Vector scalar multiplication result:" << std::endl;
    Print(scalarMultResult);

    // Test Vector cross product (for 3D vectors)
    std::cout << "\nTesting Vector cross product...\n";
    Vector<double> vec3D1({1.0, 0.0, 0.0});
    Vector<double> vec3D2({0.0, 1.0, 0.0});
    Vector<double> crossProductResult = Vector<double>::cross(vec3D1, vec3D2);
    std::cout << "Cross product of vec3D1 and vec3D2:" << std::endl;
    Print(crossProductResult);

    return 0;
}

