#ifndef MATRIX2_H
#define MATRIX_H

template <class T>
class Matrix2
{
    public:
        //constructors
        Matrix2();
        Matrix2(int nRows, int nCols);
        Matrix2(int nRows, int nCols, const T *inputData);
        Matrix2(const Matrix2<T> &inputMatrix);

        //destructor
        ~Matrix2();

        //to resize the array
        bool resize(int numRows, int nCols);


        //Element access methods
        T GetElement(int row, int col);
        bool SetElement(int row, int col, T elementVal);
        int GetNumRows();
        int GetNumCols();

        //overloading == operator to check equality bw matrix
        bool operator== (const Matrix2<T>& rhs);

        template <class U> friend Matrix2<U> operator+ (const Matrix<U>& lhs, const Matrix2<U>& rhs);
        template <class U> friend Matrix2<U> operator+ (const U& lhs, const Matrix2<U>& rhs);
        template <class U> friend Matrix2<U> operator+ (const Matrix<U>& lhs, const U& rhs);
        
        template <class U> friend Matrix2<U> operator- (const Matrix<U>& lhs, const Matrix2<U>& rhs);
        template <class U> friend Matrix2<U> operator- (const U& lhs, const Matrix2<U>& rhs);
        template <class U> friend Matrix2<U> operator- (const Matrix<U>& lhs, const U& rhs);

        template <class U> friend Matrix2<U> operator* (const Matrix<U>& lhs, const Matrix2<U>& rhs);
        template <class U> friend Matrix2<U> operator* (const U& lhs, const Matrix2<U>& rhs);
        template <class U> friend Matrix2<U> operator* (const Matrix<U>& lhs, const U& rhs);

    private:
        int Sub2Ind(int row, int col);
        T *m_matrixData;
        int m_nRows, m_nCols, m_nElements;

};

//constructors 
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
//Constructor from input array
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
//Copy constructor

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

//Destructor
template<class T>
Matrix2<T>::~Matrix2()
{
    if (m_matrixData != nullptr)
        delete[] m_matrixData;
}


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

// Element functions
template <class T>
T Matrix2<T>::GetElement(int row, int col)
{
    int linearIndex = Sub2Ind(row, col);
    if (linearIndex >= 0)
    {
        if (linearIndex >= 0)
            return m_matrixData[linearIndex];
    }
}

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

template <class T>
int Matrix2<T>::GetNumRows()
{
    return m_nRows;
}

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

template <class T>
Matrix2<T> operator+(const Matrix2<T>& rhs, const T& lhs)
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

// Scalar * Matrix
template <class T>
Matrix2<T> operator*(const T& lhs, const Matrix2<T>& rhs)
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
Matrix2<T> operator*(const Matrix2<T>& lhs, const T& rhs)
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
// matrix * matrix
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
        
        qbMatrix2<T> result(lhs.m_nRows, r_numCols, tempResult);
        delete[] tempResult;
        return result;
    }
    else
    {
        qbMatrix2<T> result(1, 1);
        return result;
    }
}


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

template <class T>
int Matrix2<T>::Sub2Ind(int row, int col)
{
    if((row < M_nRows) && (row >= 0) && (col < m_nCols) && (col >= 0))
        return row*m_nCol + col;
    else
        return -1;
}

#endif