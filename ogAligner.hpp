/* 
 * File:   ogAligner.hpp
 * Author: victortrevino
 *
 * Created on July 15, 2022, 12:04 AM
 */

#ifndef OGALIGNER_HPP
#define OGALIGNER_HPP

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "ogReadKeyMapping.hpp"
#include "ogGenome.hpp"
#include "bindings/cpp/WFAligner.hpp"
//#include "sequenceRead.h"
#include "ogSingleRead.hpp"
#include "ogCigarOperations.hpp"

using namespace wfa;

class ogAligner {
    
    ogGenome                *pGenome;
    
                            // This object contains memory STATE so it is not multi-threding safe
                            // Would need to allocate 1 object per thread and dedicate the object to the threa
    int16_t                 *pMatScore;
    int8_t                  *pMatDecision;
    int16_t                 *pLeftMax;
    int16_t                 *pUpMax;
    int16_t                 *pUpPos, *pLeftPos;
    uint16_t                allocated;
    ogCigarOperations       *pCigarOperations;
    //WFAlignerGapAffine      *pWFAligner;
    //char                    *DNA;
    //uint32_t                dnaAlloc;

public:
    WFAlignerGapAffine2Pieces   *pWFAligner;
    uint64_t                fullAligments;
                            ogAligner(ogGenome *ptrGenome);
    uint64_t                getFullAlignments();
    virtual                ~ogAligner();
    void                    allocate(uint16_t maxLen);
    //void                    allocDNA(uint32_t len);
    uint32_t                checkKeySeq(char *pSeq, uint32_t seqLen, ogReadKeyMapping *pMap);
    uint32_t                align(ogSingleRead *read, ogReadKeyMapping *pMap, char checkKeySeq);
    uint32_t                simpleAlign(char *pSeq, uint32_t seqLen, ogReadKeyMapping *pMap, char show, char dynAlign);
    uint16_t                globalAlignment(char *pA, uint32_t lenA, char *pB, uint32_t lenB);
    uint32_t                complexAlign(char *pSeq, uint32_t seqLen, ogReadKeyMapping *pMap);
    uint32_t                completeAlignment(ogSingleRead *read, ogReadKeyMapping *pMap, char checkKeySeq);
    uint32_t                fastAlignment(ogSingleRead *read, ogReadKeyMapping *pMap, char checkKeySeq);
    //uint32_t                compareSeqAndDNABuildingCIGAR(char *pSeq, char *dna, ogReadKeyMapping *pRKM, ogSingleRead *read);
    char                    buildCIGARfromWFA(char *pCIGARchars, ogSingleRead *read, uint32_t leftPos, uint32_t rightPos);
    uint32_t                windowMatches(char *pStart, uint32_t maxLen, signed char direction, uint16_t maxMismatches, uint16_t window);
    uint32_t                windowMismatches(char *pStart, uint32_t maxLen, signed char direction, uint16_t minMismatches, uint16_t window, uint16_t breakMatches);
    void                    writeOnlyMatchesCigar(ogSingleRead *read, uint32_t leftPos, uint32_t rightPos);

private:

};

#endif /* OGALIGNER_HPP */

