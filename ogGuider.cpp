/* 
 * File:   ogGuider.cpp
 * Author: victortrevino
 * 
 * Created on July 1, 2022, 8:00 PM
 */

#include <string>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include "ogDefinitions.hpp"
#include "ogGuider.hpp"

ogGuider::ogGuider() {
    //memset(ascii2RevComp, 0, sizeof(ascii2RevComp));
    isSymmetric = 1;
    fixedGuideLength = 1;
    complementSize = 0;
    pComplement = NULL;
    //allocateComplement(10*1024);
    isComplement = 0;
    genomeMode = 0;
    configurationFile[0] = 0;
}

ogGuider::~ogGuider() {
    //fprintf(stderr, "<Deallocating ogGuider:"); fflush(stderr);
    //fprintf(stderr, "<~ogG:"); fflush(stderr);
    // deallocation if needed...
    //free(pComplement); // pComplement now uses the sequenceRead structure
    //fprintf(stderr, ":ogG>");
}

void ogGuider::allocateComplement(uint32_t size) {
    if (size > complementSize) {
        complementSize = size; //((float) size / 1024.0 + 1.0) * 1024;
        pComplement = (char *) realloc(pComplement, complementSize * sizeof(char));
    }
}

ogGuider *ogGuider::clone() {
    //fprintf(stderr, "ogGuider.clone()\n");
    ogGuider *ogG = new ogGuider();
    ogG->isSymmetric = isSymmetric;
    ogG->fixedGuideLength = fixedGuideLength;
    ogG->keySizeInChars = keySizeInChars;
    return ogG;
}

void ogGuider::initializeForIndexing() {
    // Prepare for loading data if needed, default: do nothing.
}

void ogGuider::setSequence(char *pSeq, uint32_t sLen) {
    pSequence = pI = pCurrent = pSeq;
    seqLen = sLen;
    count = countOverLen = 0;
    countTotal = 0;
    //pSeqLimit = pSeq + sLen - 1;
    pSeqLimit = pLimFwd = pSequence + seqLen - keySizeInChars - 1;
    pLimRev = pSequence + keySizeInChars - 1;
    finish = 0;
    isComplement = 0;
}

void ogGuider::setSequenceAsReverseComplement(char *pSeqComp) {
// Old Version
//    //static int pelos = 0;
//    //pelos++;
//    char *pF = pSequence - 1;
//    char *pR = pSequence + seqLen;
//    char w;
//    while (++pF < --pR) {
//        //if (ascii2RevComp[*pF] == 0 || ascii2RevComp[*pR] == 0) {
//        //    fprintf(stderr, "*** PROBLEM with sequence %d *** F=%hhu, R=%hhu, f=%c, r=%c, fPos=%ld, rPos=%ld, Len=%u\n", pelos, *pF, *pR, *pF, *pR, pF-pSequence, pR-pSequence, seqLen);
//        //}
//        w = ascii2RevComp[*pF];
//        *pF = ascii2RevComp[*pR];
//        *pR = w;
//    }
//    if (pF == pR) *pF = ascii2RevComp[*pF];

    //allocateComplement(seqLen);
    pComplement = pSeqComp;

////// OLD METHOD     
//    char *pF = pComplement;
//    char *pK = pSequence + seqLen - 1;
//    uint32_t n;
//    for (n=seqLen; n; n--) {
//        *pF++ = ascii2Complement[*pK--];
//    }
//    *pF = 0;
/////////////////////////////////////////

    count = countOverLen = 0;
    countTotal = 0;
    finish = 0;
    pI = pCurrent = pComplement;
    pSeqLimit = pLimFwd = pComplement + seqLen - keySizeInChars - 1;
    pLimRev = pComplement + keySizeInChars - 1;
    isComplement = 1;
}

void ogGuider::prepareAfterNewSequence() {    
}


void ogGuider::setConfigFile(char *pFileName) {
    strncpy(configurationFile, pFileName, MAX_FILENAME);
}

char *ogGuider::getConfigFile() {
    return configurationFile;
}

char ogGuider::nextGuide() {
    if (finish) return 0;
    if (*pI) {
        pCurrent = pI;
        pI += fixedGuideLength;
        if (pI > pSeqLimit) {
            finish = 1;
            return 0;
        }
        count++;
        return 1;
    } else {
        finish = 1;
        return 0;
    }
}


uint32_t ogGuider::countGuides() {
    if (countTotal == 0) {
        // save state
        char *pIpI = pI;
        char *pCur = pCurrent;
        char end = finish;
        uint32_t guides = count;

        //reset
        pCurrent = pI = (isComplement ? pComplement : pSequence);
        finish = 0;
        count = 0;
        while (nextGuide());    
        countTotal = count;

        // restore state
        pI = pIpI;
        pCurrent = pCur;
        finish = end;
        count = guides;

        //
    }
    return countTotal;
}

char * ogGuider::getCurrentGuide() {
    return pCurrent;
}

uint32_t ogGuider::getCurrentGuidePosition() {
    return pCurrent - (isComplement ? pComplement : pSequence);
}

const char *ogGuider::getShortExtensionName() {
    return "_gDef";
}

uint32_t ogGuider::load(FILE *pFile) {
    // initialize object from loading parameters or files starting with name pFileName
    // default: do nothing.
    // return the size of loaded bytes
    return 0;
}

uint32_t ogGuider::save(FILE *pFile) {
    // save file if needed after indexing preparing for loading in mapping
    // default: do nothing.
    // return the size of written bytes
    return 0;
}

const char name[] = "DefaultGuider";

const char * ogGuider::getName() {
    return name;
}

void ogGuider::setKeySizeInChars(uint16_t keySize) {
    keySizeInChars = keySize;
}


uint16_t ogGuider::getKeySizeInChars() {
    return keySizeInChars;
}


void ogGuider::fixSequenceCase() {
    char *p = pSequence;
    uint32_t l = seqLen;
    //uint32_t f = 0;
    const char delta_A = 'a' - 'A';
    for (; l && *p; l--, p++) {
        if (*p >= 'a') {
            *p -= delta_A;
        }
    }    
}

void ogGuider::setIsSymmetric(char sym) {
    isSymmetric = sym;
}

char ogGuider::getIsSymmetric() {
    return isSymmetric;
}

void ogGuider::setGenomeMode(char mode) {
    genomeMode = mode;
}

