#ifndef MY_VECTOR_H
#define MY_VECTOR_H

#include <stdexcept>

template <class T>
class MyVector {
public:
    // Constructors
    MyVector();
    MyVector(size_t initialCapacity);
    MyVector(const MyVector& other);

    // Destructor
    ~MyVector();

    // Capacity
    size_t size() const;
    size_t capacity() const;
    bool empty() const;

    // Element access
    T& at(size_t index);
    const T& at(size_t index) const;

    // Modifiers
    void push_back(const T& value);

private:
    T* m_data;
    size_t m_size;
    size_t m_capacity;
};

template <class T>
MyVector<T>::MyVector() : m_data(nullptr), m_size(0), m_capacity(0) {}

template <class T>
MyVector<T>::MyVector(size_t initialCapacity) : m_size(0), m_capacity(initialCapacity) {
    m_data = new T[m_capacity];
}

template <class T>
MyVector<T>::MyVector(const MyVector& other) : m_size(other.m_size), m_capacity(other.m_capacity) {
    m_data = new T[m_capacity];
    for (size_t i = 0; i < m_size; ++i) {
        m_data[i] = other.m_data[i];
    }
}

template <class T>
MyVector<T>::~MyVector() {
    delete[] m_data;
}

template <class T>
size_t MyVector<T>::size() const {
    return m_size;
}

template <class T>
size_t MyVector<T>::capacity() const {
    return m_capacity;
}

template <class T>
bool MyVector<T>::empty() const {
    return m_size == 0;
}

template <class T>
T& MyVector<T>::at(size_t index) {
    if (index >= m_size) {
        throw std::out_of_range("Index out of range");
    }
    return m_data[index];
}

template <class T>
const T& MyVector<T>::at(size_t index) const {
    if (index >= m_size) {
        throw std::out_of_range("Index out of range");
    }
    return m_data[index];
}

template <class T>
void MyVector<T>::push_back(const T& value) {
    if (m_size >= m_capacity) {
        size_t newCapacity = (m_capacity == 0) ? 1 : m_capacity * 2;
        T* newData = new T[newCapacity];
        for (size_t i = 0; i < m_size; ++i) {
            newData[i] = m_data[i];
        }
        delete[] m_data;
        m_data = newData;
        m_capacity = newCapacity;
    }
    m_data[m_size++] = value;
}

#endif // MY_VECTOR_H
