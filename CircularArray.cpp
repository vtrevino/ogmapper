/* 
 * File:   CircularArray.cpp
 * Author: victortrevino
 * 
 * Created on June 22, 2022, 6:47 PM
 */

#include "CircularArray.hpp"

template <class T>
CircularArray<T>::CircularArray() {
    blockSize = length = 10;
    sizeofdata = sizeof(T);
    utilized = 0;
    first = 0;
    last = 9;
    pMem = NULL;
    allocate();
}

template <class T>
CircularArray<T>::CircularArray(uint32_t initialSize, uint32_t newBlocks) {
    blockSize = newBlocks;
    length = initialSize;
    sizeofdata = sizeof(T);
    utilized = 0;
    first = 0;
    last = length - 1;
    pMem = NULL;
    allocate();
}

template <class T>
CircularArray<T>::~CircularArray() {
    free(pMem);
}

template <class T>
void CircularArray<T>::reset() {
    utilized = 0;
    first = 0;
    last = length - 1;
}

template <class T>
void CircularArray<T>::allocate() {
    pMem = (T*) realloc(pMem, length * sizeofdata);
}

template <class T>
uint32_t CircularArray<T>::push(T value) {
    uint32_t added = 0;
    if (++utilized > length) {
        added = blockSize;
        length = length + added;
        allocate();
    }
    if (++last >= length) last = 0;
    pMem[last] = value;
    return added;
}

template <class T>
T CircularArray<T>::pull() {
    if (utilized == 0) return NULL;
    T p = pMem[first];
    if (++first >= length) first = 0;
    utilized--;
    return p;
}

template <class T>
T CircularArray<T>::checkNextValue() {
    if (utilized == 0) return NULL;
    return pMem[first];
}

template <class T>
uint32_t CircularArray<T>::size() {
    return length;
}

template <class T>
uint32_t CircularArray<T>::used() {
    return utilized;
}
