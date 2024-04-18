#ifndef MATRIX2_H
#define MATRIX_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>

template <class T>
class Matrix2
{
public:
    // Constructors
    static_assert(std::is_arithmetic<T>::value, "Matrix2 can only be instantiated with numeric types.");
    Matrix2(); // Default constructor
    Matrix2(int nRows, int nCols); // Construct empty matrix (all elements 0)
    Matrix2(int nRows, int nCols, const T *inputData); // Constructor from input array
    Matrix2(const Matrix2<T> &inputMatrix); // Copy constructor
    Matrix2(int nRows, int nCols, const std::vector<T> *inputData); // Constructor from vector
    // Destructor
    ~Matrix2();

    // Resize the array
    bool resize(int numRows, int nCols);

    // Element access methods
    T GetElement(int row, int col); // Get value at specific position
    bool SetElement(int row, int col, T elementVal); // Set value at specific position
    int GetNumRows(); // Get number of rows
    int GetNumCols(); // Get number of columns

    // Inverse of the matrix
    Matrix2<T> inverse(); // Computes the inverse of the matrix

    // Swap rows
    void swapRows(int row1, int row2); // Swaps two rows of the matrix

    // Overloaded operators
    bool operator== (const Matrix2<T>& rhs); // Equality operator
    template <class... Args>
    friend Matrix2<T> addMatrices(const Matrix2<T>& matrix, Args... args); // Add matrices

    // Addition operators
    template <class U> friend Matrix2<U> operator+ (const Matrix2<U>& lhs, const Matrix2<U>& rhs); // Matrix + Matrix
    template <class U> friend Matrix2<U> operator+ (const U& lhs, const Matrix2<U>& rhs); // Scalar + Matrix
    template <class U> friend Matrix2<U> operator+ (const Matrix2<U>& lhs, const U& rhs); // Matrix + Scalar

    // Subtraction operators
    template <class U> friend Matrix2<U> operator- (const Matrix2<U>& lhs, const Matrix2<U>& rhs); // Matrix - Matrix
    template <class U> friend Matrix2<U> operator- (const U& lhs, const Matrix2<U>& rhs); // Scalar - Matrix
    template <class U> friend Matrix2<U> operator- (const Matrix2<U>& lhs, const U& rhs); // Matrix - Scalar

    // Multiplication operators
    template <class U> friend Matrix2<U> operator* (const Matrix2<U>& lhs, const Matrix2<U>& rhs); // Matrix * Matrix
    template <class U> friend Matrix2<U> operator* (const U& lhs, const Matrix2<U>& rhs); // Scalar * Matrix
    template <class U> friend Matrix2<U> operator* (const Matrix2<U>& lhs, const U& rhs); // Matrix * Scalar

private:
    int Sub2Ind(int row, int col); // Convert 2D index to 1D index
    T *m_matrixData; // Pointer to matrix data
    int m_nRows, m_nCols, m_nElements; // Number of rows, columns, and elements
};

// Default constructor
template <class T>
Matrix2<T>::Matrix2()
{
    m_nRows = 1;
    m_nCols = 1;
    m_nElements = 1;
    m_matrixData = new T[m_nElements];
    m_matrixData[0] = 0.0;
}

// Construct empty matrix (all elements 0)
template <class T>
Matrix2<T>::Matrix2(int nRows, int nCols)
{
    m_nRows = nRows;
    m_nCols = nCols;
    m_nElements = m_nRows * m_nCols;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; i++)
        m_matrixData[i] = 0.0;
}

// Constructor from input array
template<class T>
Matrix2<T>::Matrix2(int nRows, int nCols, const T *inputData)
{
    m_nRows = nRows;
    m_nCols = nCols;
    m_nElements = m_nRows * m_nCols;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; i++)
        m_matrixData[i] = inputData[i];
}

// Constructor from vector
template<class T>
Matrix2<T>::Matrix2(int nRows, int nCols, const std::vector<T> *inputData)
{
    m_nRows = nRows;
    m_nCols = nCols;
    m_nElements = m_nRows * m_nCols;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; i++)
        m_matrixData[i] = inputData->at(i);
}

// Copy constructor
template<class T>
Matrix2<T>::Matrix2(const Matrix2<T>& inputMatrix)
{
    m_nRows = inputMatrix.m_nRows;
    m_nCols = inputMatrix.m_nCols;
    m_nElements = inputMatrix.m_nElements;
    m_matrixData = new T[m_nElements];
    for (int i=0; i<m_nElements; i++)
        m_matrixData[i] = inputMatrix.m_matrixData[i];
}

// Destructor
template<class T>
Matrix2<T>::~Matrix2()
{
    if (m_matrixData != nullptr)
        delete[] m_matrixData;
}

// Resize the array
template <class T>
bool Matrix2<T>::resize(int numRows, int numCols)
{
    m_nRows = numRows;
    m_nCols = numCols;
    m_nElements = (m_nRows * m_nCols);
    delete[] m_matrixData;
    m_matrixData = new T[m_nElements];
    if (m_matrixData != nullptr)
    {
        for (int i = 0; i < m_nElements; i++)
            m_matrixData[i] = 0.0;
        return true;
    }
    else
    {
        return false;
    }
}

// Get value at specific position
template <class T>
T Matrix2<T>::GetElement(int row, int col)
{
    int linearIndex = Sub2Ind(row, col);
    if (linearIndex >= 0)
        return m_matrixData[linearIndex];
    else    
        return 0.0;
}

// Set value at specific position
template <class T>
bool Matrix2<T>::SetElement(int row, int col, T elementValue)
{
    int linearIndex = Sub2Ind(row, col);
    if (linearIndex >= 0)
    {
        m_matrixData[linearIndex] = elementValue;
        return true;
    }
    else
    {
        return false;
    }
}

// Get number of rows
template <class T>
int Matrix2<T>::GetNumRows()
{
    return m_nRows;
}

// Get number of columns
template <class T>
int Matrix2<T>::GetNumCols()
{
    return m_nCols;
}

// Matrix + Matrix
template <class T>
Matrix2<T> operator+(const Matrix2<T>& lhs, const Matrix2<T>& rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i)
        tempResult[i] = lhs.m_matrixData[i] + rhs.m_matrixData[i];
    Matrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// Specialization
template <>
Matrix2<int> operator+(const Matrix2<int>& lhs, const Matrix2<int>& rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = numRows * numCols;
    int *tempResult = new int[numElements];
    for (int i = 0; i < numElements; ++i)
        tempResult[i] = lhs.m_matrixData[i] + rhs.m_matrixData[i];
    Matrix2<int> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// Scalar + Matrix
template <class T>
Matrix2<T> operator+(const T& lhs, const Matrix2<T>& rhs)
{
    int numRows = rhs.m_nRows;
    int numCols = rhs.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i)
        tempResult[i] = lhs + rhs.m_matrixData[i];
    Matrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// Matrix + Scalar
template <class T>
Matrix2<T> operator+(const Matrix2<T>& lhs, const T& rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i)
        tempResult[i] = lhs.m_matrixData[i] + rhs;
    Matrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// Scalar - Matrix
template <class T>
Matrix2<T> operator-(const T& lhs, const Matrix2<T>& rhs)
{
    int numRows = rhs.m_nRows;
    int numCols = rhs.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i)
        tempResult[i] = lhs - rhs.m_matrixData[i];
    Matrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// Matrix - Scalar
template <class T>
Matrix2<T> operator-(const Matrix2<T>& rhs, const T& lhs)
{
    int numRows = rhs.m_nRows;
    int numCols = rhs.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i)
        tempResult[i] = rhs.m_matrixData[i] - lhs;
    Matrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// Matrix - Matrix
template <class T>
Matrix2<T> operator-(const Matrix2<T>& lhs, const Matrix2<T>& rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i)
        tempResult[i] = lhs.m_matrixData[i] - rhs.m_matrixData[i];
    Matrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// Scalar * Matrix
template <class T>
Matrix2<T> operator* (const T& lhs, const Matrix2<T>& rhs)
{
    int numRows = rhs.m_nRows;
    int numCols = rhs.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i)
        tempResult[i] = lhs * rhs.m_matrixData[i];
    Matrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// Matrix * Scalar
template <class T>
Matrix2<T> operator* (const Matrix2<T>& lhs, const T& rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i)
        tempResult[i] = lhs.m_matrixData[i] * rhs;
    Matrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// Matrix * Matrix
template <class T>
Matrix2<T> operator* (const Matrix2<T>& lhs, const Matrix2<T>& rhs)
{
    int r_numRows = rhs.m_nRows;
    int r_numCols = rhs.m_nCols;
    int l_numCols = lhs.m_nCols;

    if (l_numCols == r_numRows)
    {
        // The output will be the same size as the RHS.
        T *tempResult = new T[lhs.m_nRows * rhs.m_nCols];

        // Loop through each row of the LHS.
        for (int lhsRow = 0; lhsRow < lhs.m_nRows; ++lhsRow)
        {
            // Loop through each column on the RHS row.
            for (int rhsCol = 0; rhsCol < r_numCols; ++rhsCol)
            {
                T elementResult = 0.0;

                // Loop through each element of this LHS row.
                for (int lhsCol = 0; lhsCol < l_numCols; ++lhsCol)
                {
                    // Compute the LHS linear index.
                    int lhsLinearIndex = (lhsRow * l_numCols) + lhsCol;

                    // Compute the RHS linear index (based on LHS col).
                    // rhs row number equal to lhs column number.
                    int rhsLinearIndex = (lhsCol * r_numCols) + rhsCol;

                    // Perform the calculation on these elements.
                    elementResult += (lhs.m_matrixData[lhsLinearIndex] * rhs.m_matrixData[rhsLinearIndex]);
                }

                // Store the result.
                int resultLinearIndex = (lhsRow * r_numCols) + rhsCol;
                tempResult[resultLinearIndex] = elementResult;
            }
        }

        Matrix2<T> result(lhs.m_nRows, r_numCols, tempResult);
        delete[] tempResult;
        return result;
    }
    else
    {
        Matrix2<T> result(1, 1);
        return result;
    }
}

// Equality operator
template <class T>
bool Matrix2<T>::operator== (const Matrix2<T>& rhs)
{
    // Check if the matrices are the same size, if not return false.
    if ((this->m_nRows != rhs.m_nRows) && (this->m_nCols != rhs.m_nCols))
        return false;

    // Check if the elements are equal.
    bool flag = true;
    for (int i = 0; i < this->m_nElements; ++i)
    {
        if (this->m_matrixData[i] != rhs.m_matrixData[i])
            flag = false;
    }
    return flag;
}

// Convert 2D index to 1D index
template <class T>
int Matrix2<T>::Sub2Ind(int row, int col)
{
    if((row < m_nRows) && (row >= 0) && (col < m_nCols) && (col >= 0))
        return row*m_nCols + col;
    else
        return -1;
}

// Compute the inverse of the matrix
template <class T>
Matrix2<T> Matrix2<T>::inverse() {
    // Check if the matrix is square
    if (m_nRows != m_nCols) {
        // Not a square matrix, return an empty matrix
        return Matrix2<T>(0, 0);
    }

    // Create a copy of this matrix
    Matrix2<T> result(*this);

    // Create an identity matrix of the same size as the original matrix
    Matrix2<T> identity(m_nRows, m_nCols);
    for (int i = 0; i < m_nRows; ++i) {
        identity.SetElement(i, i, 1);
    }

    // Perform Gaussian elimination with partial pivoting to compute the inverse
    for (int i = 0; i < m_nRows; ++i) {
        // Find the pivot element
        int maxRow = i;
        for (int j = i + 1; j < m_nRows; ++j) {
            if (std::abs(result.GetElement(j, i)) > std::abs(result.GetElement(maxRow, i))) {
                maxRow = j;
            }
        }

        // Swap rows if necessary
        if (maxRow != i) {
            result.swapRows(i, maxRow);
            identity.swapRows(i, maxRow);
        }

        // Divide the current row by the pivot element
        T pivot = result.GetElement(i, i);
        if (pivot == 0) {
            // Matrix is singular, return an empty matrix
            return Matrix2<T>(0, 0);
        }
        for (int j = 0; j < m_nCols; ++j) {
            result.SetElement(i, j, result.GetElement(i, j) / pivot);
            identity.SetElement(i, j, identity.GetElement(i, j) / pivot);
        }

        // Eliminate all other elements in the current column
        for (int j = 0; j < m_nRows; ++j) {
            if (j != i) {
                T factor = result.GetElement(j, i);
                for (int k = 0; k < m_nCols; ++k) {
                    result.SetElement(j, k, result.GetElement(j, k) - factor * result.GetElement(i, k));
                    identity.SetElement(j, k, identity.GetElement(j, k) - factor * identity.GetElement(i, k));
                }
            }
        }
    }

    return identity;
}

// Swap rows of the matrix
template <class T>
void Matrix2<T>::swapRows(int row1, int row2) {
    if (row1 < 0 || row1 >= m_nRows || row2 < 0 || row2 >= m_nRows) {
        // Invalid row indices
        return;
    }

    for (int j = 0; j < m_nCols; ++j) {
        T temp = m_matrixData[row1 * m_nCols + j];
        m_matrixData[row1 * m_nCols + j] = m_matrixData[row2 * m_nCols + j];
        m_matrixData[row2 * m_nCols + j] = temp;
    }
}

// Base case for adding matrices
template <class T>
Matrix2<T> addMatrices(Matrix2<T>& matrix) {
    return matrix; // Base case: return the matrix itself
}

// Recursive case for adding matrices
template <class T, class... Args>
Matrix2<T> addMatrices(Matrix2<T>& matrix, Args... args) {
    return matrix + addMatrices(args...); // Recursive case: add the matrix to the result of adding the rest
}

#endif
