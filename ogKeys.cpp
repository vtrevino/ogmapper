#include <cstddef>
#include <stdint.h>
#include <cstdlib>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ogDefinitions.hpp"

#include "ogKeys.hpp"


ogKeys::ogKeys(uint8_t size, uint16_t sizeChars, ogGenomePositions *gp, uint32_t totalKeys) {
    sizeInBits = 0;
    sizeInChars = sizeChars;
    arraySize = totalKeys;
    keys = NULL;
    setSizeInBits(size);
    setGenomePos(gp);
    allocate();
}

ogKeys::~ogKeys() {
    //fprintf(stderr, "<Deallocating ogKeys:");
    //fprintf(stderr, "<~ogK:"); fflush(stderr);
    deallocate();
    //fprintf(stderr, ":ogK>");
}

void ogKeys::freeArray() {
    if (keys != NULL) {
        free(keys);
        keys = NULL;
    }
    arraySize = 0;
}

void ogKeys::setSizeInBits(uint8_t size) {
    sizeInBits = size;
    freeArray();
}

void ogKeys::setGenomePos(ogGenomePositions *gp) {
    gPos = gp;
}

uint8_t ogKeys::getSizeInBits() {
    return sizeInBits;
}

uint32_t ogKeys::getArraySize() {
    return arraySize;
}

void ogKeys::addPosToKey(uint32_t genPos, uint32_t key) {
    if (key >= arraySize) {
        fprintf(stderr, "\n**********************************************\n");
        fprintf(stderr, "**** Problem accessing: %u of %u\n", key, arraySize);
        fprintf(stderr, "**********************************************\n");
    }
    ogKey *pK = keys+key;
    gPos->setPosition(pK->offset + pK->size++, genPos);
}

void ogKeys::incKey(uint32_t key) {
    if (key >= arraySize) {
        fprintf(stderr, "\n**********************************************\n");
        fprintf(stderr, "**** Problem accessing: %u of %u\n", key, arraySize);
        fprintf(stderr, "**********************************************\n");
    }
    keys[key].size++; // used in counting mode (not in storing mode)
}

uint32_t ogKeys::getSizeForKey(uint32_t key) {
    return keys[key].size;
}

void ogKeys::resetKeySizesAndOffSets() {
    memset(keys, 0, memUsed);
}

void ogKeys::resetKeySizes() {
    uint32_t i;
    uint32_t maxOS = 0;
    ogKey *pK = keys;
    for (i=0; i < arraySize; i++) {
        pK->size = 0;
        if (pK->offset > maxOS) maxOS = pK->offset;
        pK++;
    }
    //fprintf(stderr, "max offset=%u\n", maxOS);
}

void ogKeys::estimateOffSets() {
    // after counting mode, update offsets
    uint32_t    i;
    uint32_t    os = 0;
    ogKey *pK = keys;
    for (i=0; i < arraySize; i++) {
        pK->offset = os;
        os += pK->size;
        pK++;
    }
    fprintf(stderr, "Maximum offset = %u\n", os);
}

void ogKeys::allocate() {
    freeArray();
    if (arraySize == 0) arraySize = pow(2, sizeInBits);
    if (arraySize == 0) arraySize = 0xFFFFFFFF;
    memUsed = sizeof(ogKey) * arraySize;
    fprintf(stderr, "Allocating %.1f Mb (%llu bytes) of memory for %u keys having %hhu bits.\n", (float) memUsed/(1024*1024), memUsed, arraySize, sizeInBits);
    keys = (ogKey *) calloc(memUsed, 1);
}

void ogKeys::deallocate() {
    freeArray();
}


void ogKeys::load(FILE *pFile) {
    //FILE           *pfile;
    //char            xfile[1000];
    //setPureFileName(xfile, filename, "ogk", 1000);
    //printf_FileOperation("Loading keys");
    //pfile = fopen(xfile, "r");
    fread(keys, memUsed, 1, pFile);
    //fclose(pfile);
    //printf_FileOperationDone();
}

uint64_t ogKeys::save(FILE *pFile) {
    //FILE           *pfile;
    //char            xfile[1000];
    //setFileName(xfile, filename, sizeInChars, "ogk", 1000);
    printf_FileOperation("Saving keys");
    //pfile = fopen(xfile, "wb");
    uint64_t bytes = memUsed * fwrite(keys, memUsed, 1, pFile);
    //fclose(pfile);
    printf_FileOperationDone();
    return bytes;
}

void ogKeys::loadGenomePositions(FILE *pFile) {
    gPos->load(pFile);
}

uint64_t ogKeys::saveGenomePositions(FILE *pFile) {
    return gPos->save(pFile);
}

void ogKeys::saveSizes(char *filename) {
    FILE           *pfile;
    char            xfile[MAX_FILENAME];
    setFileName(xfile, filename, sizeInChars, "TO_REMOVE.ogX", MAX_FILENAME);
    printf_FileOperationFile("Saving sizes per key", xfile);
    pfile = fopen(xfile, "w");
    uint32_t i;
    ogKey *pK = keys;
    fprintf(pfile, "GenPosPerKey\n");
    for (i=0; i < arraySize; i++) {
        fprintf(pfile, "%u\n",pK->size);
        pK++;
    }    
    fclose(pfile);
    printf_FileOperationDone();
    printf_FileOperationFile("(Informative file, can be deleted)", xfile);
    printf_FileOperationDone();
    //fprintf(stderr, "====================================================\n");
    //fprintf(stderr, "This file is informative, it can be safetly deleted.\n");
    //fprintf(stderr, "====================================================\n");
}

uint32_t* ogKeys::getPointerPositionForKey(uint32_t key) {
    ogKey *pK = keys+key;
    return gPos->getPointerPosition(pK->offset);
}

ogKey* ogKeys::getInfoForKey(uint32_t key) {
    return keys+key;
}

void ogKeys::printKeyInfo(uint32_t key) {
    ogKey *pK = keys+key;
    fprintf(stderr, "Key: %u, Offset: %u, Size: %u\n", key, pK->offset, pK->size);
}

uint32_t ogKeys::removeKeySizesLargerThan(uint32_t moreThan) {
    uint32_t i;
    uint32_t n=0;
    ogKey *pK = keys;
    for (i=0; i < arraySize; i++) {
        if (pK->size >= moreThan) {
            pK->size = 0;
            n++;
        }
        pK++;
    }
    return n;
}