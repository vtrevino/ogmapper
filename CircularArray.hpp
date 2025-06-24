#include <stdio.h>
#include <string.h>
//#include "testFastFiles.h"
#include <unistd.h>
/* 
 * File:   CircularArray.hpp
 * Author: victortrevino
 *
 * Created on June 22, 2022, 6:47 PM
 */

#ifndef CIRCULARARRAY_HPP
#define CIRCULARARRAY_HPP

#include <stdint.h>

template <class T>

class CircularArray {
    uint32_t    length;
    T           *pMem;
    uint32_t    sizeofdata;
    uint32_t    blockSize;
    uint32_t    first, last;
    uint32_t    utilized;
    
public:
                    CircularArray();
                    CircularArray(uint32_t initialSize, uint32_t newBlocks);
    virtual        ~CircularArray();
    void            allocate();
    uint32_t        push(T value); // regresa cuántos nuevos elementos se agregaron
    T               pull();
    T               checkNextValue();
    uint32_t        size();
    uint32_t        used();
    void            reset();
    
private:

};

#endif /* CIRCULARARRAY_HPP */

