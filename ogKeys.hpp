/* 
 * File:   keys.hpp
 * Author: victortrevino
 *
 * Created on May 4, 2022, 4:03 PM
 */

#ifndef OGKEYS_HPP
#define OGKEYS_HPP

#include <stdio.h>
#include "ogGenomePositions.hpp"


typedef struct ogKey {
    uint32_t        size;
    uint32_t        offset;
} OG_KEY;

        
class ogKeys {
    uint16_t            sizeInChars;
    uint8_t             sizeInBits;         // max 32 bits
    uint32_t            arraySize;          // length of the memory asked in ogOneKey units
    ogKey               *keys;
    ogGenomePositions   *gPos;
    uint64_t            memUsed;
    
public:
                    ogKeys(uint8_t size, uint16_t sizeChars, ogGenomePositions *gp, uint32_t totalKeys);
    virtual         ~ogKeys();
    void            setSizeInChars(uint16_t size);
    void            setSizeInBits(uint8_t size);
    void            setGenomePos(ogGenomePositions *gp);
    uint8_t         getSizeInBits();
    uint16_t        getSizeInChars();
    uint32_t        getArraySize();
    void            addPosToKey(uint32_t genPos, uint32_t key);
    ogKey*          getInfoForKey(uint32_t key);
    void            load(FILE *pFile);
    uint64_t        save(FILE *pFile);
    void            resetKeySizesAndOffSets();
    void            resetKeySizes();
    void            estimateOffSets();
    void            allocate();
    void            deallocate();
    void            freeArray();
    void            incKey(uint32_t key);
    void            loadGenomePositions(FILE *pFile);
    uint64_t        saveGenomePositions(FILE *pFile);
    void            saveSizes(char *filename);
    uint32_t        getSizeForKey(uint32_t key);
    uint32_t*       getPointerPositionForKey(uint32_t key);
    void            printKeyInfo(uint32_t key);
    uint32_t        removeKeySizesLargerThan(uint32_t moreThan);
    
};


#endif /* OGKEYS_HPP */

