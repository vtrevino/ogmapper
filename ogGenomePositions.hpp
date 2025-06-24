
/* 
 * File:   ogGenomePositions.hpp
 * Author: victortrevino
 *
 * Created on May 4, 2022, 6:28 PM
 */

#ifndef OGGENOMEPOSITIONS_HPP
#define OGGENOMEPOSITIONS_HPP

#include <stdlib.h>
#include <stdio.h>
#include "ogGenome.hpp"

class ogGenomePositions {
    uint32_t    nPositions;
    uint64_t    memoryUsed;    
    uint32_t    *gPositions;
    uint32_t    keySize;
    
public:
    ogGenomePositions(uint32_t totalPositions, uint16_t userKeySize);
    virtual    ~ogGenomePositions();
    uint32_t    getNPositions();
    uint32_t    getArraySize();
    void        load(FILE *pFile);
    uint64_t    save(FILE *pFile);
    void        setPosition(uint32_t kPos, uint32_t genPos);
    uint32_t    getPosition(uint32_t kPos);
    uint32_t   *getPointerPosition(uint32_t kPos);
    uint32_t    getKeySize();
    void        printGenomePositionsInfo(uint32_t fromGenPos, uint32_t size, ogGenome *pGenome);
    
private:

};

#endif /* OGGENOMEPOSITIONS_HPP */

