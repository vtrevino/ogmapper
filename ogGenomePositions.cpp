/* 
 * File: ogGenomePositions.cpp
 * Author: victortrevino
 * 
 * Created on May 4, 2022, 6:28 PM
 */

#include <stdio.h>
#include "ogDefinitions.hpp"

#include "ogGenomePositions.hpp"

ogGenomePositions::ogGenomePositions(uint32_t totalPositions, uint16_t userKeySize) {
    nPositions = totalPositions;
    keySize = userKeySize;
    memoryUsed = nPositions * sizeof(uint32_t);
    fprintf(stderr, "Allocating %.1f MB (%llu bytes) of memory for %u genomic positions.\n", (float) memoryUsed/(1024*1024), memoryUsed, totalPositions);
    gPositions = (uint32_t *) malloc(memoryUsed);
}

ogGenomePositions::~ogGenomePositions() {
    //fprintf(stderr, "<Deallocating ogGenomePositions:");
    //fprintf(stderr, "<~ogGP:"); fflush(stderr);
    if (gPositions != NULL) free(gPositions);
    //fprintf(stderr, ":ogGP>");
}

uint32_t ogGenomePositions::getNPositions() {
    return nPositions;
}

uint32_t ogGenomePositions::getArraySize() {
    return nPositions;
}

uint32_t ogGenomePositions::getKeySize() {
    return keySize;
}

void ogGenomePositions::load(FILE *pFile) {
    //FILE *pfile;
    //char xfile[1000];
    //setPureFileName(xfile, filename, "ogo", 1000);
    //printf_FileOperation("Loading genome positions");
    //pfile = fopen(xfile, "r");
    fread(gPositions, sizeof(uint32_t), nPositions, pFile);
    //fclose(pfile);    
    //printf_FileOperationDone();
}

uint64_t ogGenomePositions::save(FILE *pFile) {
    //FILE *pfile;
    //char xfile[1000];
    //setFileName(xfile, filename, keySize, "ogo", 1000);
    printf_FileOperation("Saving genome positions");
    //pfile = fopen(xfile, "wb");
    uint64_t bytes = fwrite(gPositions, 1, memoryUsed, pFile);
    //fclose(pfile);    
    printf_FileOperationDone();
    return bytes;
}

void ogGenomePositions::setPosition(uint32_t kPos, uint32_t genPos) {
    if (kPos >= nPositions) {
        fprintf(stderr, "\n**********************************************\n");
        fprintf(stderr, "**** Problem accessing: %u of %u\n", kPos, nPositions);
        fprintf(stderr, "**********************************************\n");
    }
    gPositions[kPos] = genPos;
}

uint32_t ogGenomePositions::getPosition(uint32_t kPos) {
    return gPositions[kPos];
}

uint32_t * ogGenomePositions::getPointerPosition(uint32_t kPos) {
    return gPositions + kPos;
}

void ogGenomePositions::printGenomePositionsInfo(uint32_t fromGenPos, uint32_t size, ogGenome *pGenome) {
    uint32_t *p = gPositions + fromGenPos;
    uint32_t i;
    ogChromosome *pChr;
    char     dna[10000];
    dna[keySize] = 0;
    fprintf(stderr, "\n...........................................................................................\n");
    for (i=0; i < size; p++, i++) {
        pChr = pGenome->getGenomicCoordinate(*p);
        pGenome->extractPackedGenome(*p, keySize, dna, 0, 0);
        fprintf(stderr, "%u, index:%u, GenPos:%u, Chr:%s, ChrPos:%u, DNA:%s\n", i+1, fromGenPos+i, *p, pChr->name, *p - pChr->start, dna);
    }
}