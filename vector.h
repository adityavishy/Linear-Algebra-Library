#ifndef VECTOR_H
#define VECTOR_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

template <class T>
class Vector {
public:
    // Constructors
    static_assert(std::is_arithmetic<T>::value, "Vector can only be instantiated with numeric types.");
    Vector(); // Default constructor
    Vector(std::vector<T> inputData); // Constructor with vector data

    // Destructor
    ~Vector();

    // Function to return number of dimensions
    int GetNumDims() const;

    // Function to get element at index
    T GetElement(int index) const;

    // Overloaded operators
    Vector<T> operator+ (const Vector<T>& rhs) const; // Addition
    Vector<T> operator- (const Vector<T>& rhs) const; // Subtraction
    Vector<T> operator* (const T& rhs) const; // Scalar multiplication

    // Friend function for scalar multiplication
    template <class U> friend Vector<U> operator* (const U &lhs, const Vector<U> &rhs);

    // Static functions
    static T dot(const Vector<T> &a, const Vector<T> &b); // Dot product
    static Vector<T> cross(const Vector<T> &a, const Vector<T> &b); // Cross product

private:
    std::vector<T> m_vectorData; // Vector data
    int m_nDims; // Number of dimensions
};

// Constructors
template <class T>
Vector<T>::Vector() {
    m_nDims = 0;
    m_vectorData = std::vector<T>(); // Initialize with an empty std::vector
}

template <class T>
Vector<T>::Vector(std::vector<T> inputData) {
    m_nDims = inputData.size();
    m_vectorData = inputData; // Directly assign the input std::vector
}

// Destructor
template <class T>
Vector<T>::~Vector() {
}

// Function to return number of dimensions
template <class T>
int Vector<T>::GetNumDims() const {
    return m_nDims;
}

// Function to get element at index
template <class T>
T Vector<T>::GetElement(int index) const {
    if (index < 0 || index >= m_nDims) {
        throw std::out_of_range("Index out of range");
    }
    return m_vectorData.at(index);
}

// Operator overloading: Addition
template <class T>
Vector<T> Vector<T>::operator+(const Vector<T> &rhs) const {
    if (m_nDims != rhs.m_nDims)
        throw std::invalid_argument("Vector dimensions do not match.");

    std::vector<T> resultData;
    for (int i = 0; i < m_nDims; ++i)
        resultData.push_back(m_vectorData.at(i) + rhs.m_vectorData.at(i));

    Vector<T> result(resultData);
    return result;
}

// Operator overloading: Subtraction
template <class T>
Vector<T> Vector<T>::operator-(const Vector<T> &rhs) const {
    if (m_nDims != rhs.m_nDims)
        throw std::invalid_argument("Vector dimensions do not match.");

    std::vector<T> resultData;
    for (int i = 0; i < m_nDims; ++i)
        resultData.push_back(m_vectorData.at(i) - rhs.m_vectorData.at(i));

    Vector<T> result(resultData);
    return result;
}

// Operator overloading: Scalar multiplication
template <class T>
Vector<T> Vector<T>::operator*(const T &rhs) const {
    std::vector<T> resultData;
    for (int i = 0; i < m_nDims; ++i)
        resultData.push_back(m_vectorData.at(i) * rhs);

    Vector<T> result(resultData);
    return result;
}
 
// Friend function for scalar multiplication
template <class T>
Vector<T> operator*(const T &lhs, const Vector<T> &rhs) {
    std::vector<T> resultData;
    for (int i = 0; i < rhs.m_nDims; ++i)
        resultData.push_back(lhs * rhs.m_vectorData.at(i));

    Vector<T> result(resultData);
    return result;
}

// Helper function for recursive dot product calculation
template <class T>
T dot_helper(const Vector<T>& a, const Vector<T>& b) {
    if (a.GetNumDims() != b.GetNumDims()) {
        throw std::invalid_argument("All vectors must have the same number of dimensions.");
    }
    T result = 0;
    for (int i = 0; i < a.GetNumDims(); ++i) {
        result += a.GetElement(i) * b.GetElement(i);
    }
    return result;
}

// Variadic template dot function
template <class T, class... Args>
T dot(const Vector<T>& first, const Args&... args) {
    return (dot_helper(first, args) + ...); // Fold expression to sum dot products
}

// Base case for the recursive template
template <class T>
T dot(const Vector<T>& first, const Vector<T>& second) {
    return dot_helper(first, second);
}

// Static function for dot product
template <class T>
T Vector<T>::dot(const Vector<T> &a, const Vector<T> &b) {
    return dot_helper(a, b);
}

// Static function for cross product
template <class T>
Vector<T> Vector<T>::cross(const Vector<T> &a, const Vector<T> &b) {
    if (a.m_nDims != b.m_nDims)
        throw std::invalid_argument("Vector dimensions must match for the cross-product to be computed.");
    if (a.m_nDims != 3)
        throw std::invalid_argument("The cross-product can only be computed for three-dimensional vectors.");

    std::vector<T> resultData;
    resultData.push_back((a.m_vectorData.at(1) * b.m_vectorData.at(2)) - (a.m_vectorData.at(2) * b.m_vectorData.at(1)));
    resultData.push_back((a.m_vectorData.at(2) * b.m_vectorData.at(0)) - (a.m_vectorData.at(0) * b.m_vectorData.at(2)));
    resultData.push_back((a.m_vectorData.at(0) * b.m_vectorData.at(1)) - (a.m_vectorData.at(1) * b.m_vectorData.at(0)));

    Vector<T> result(resultData);
    return result;
}

#endif
