/* 
 * File:   ogGuider.hpp
 * Author: victortrevino
 *
 * Created on July 1, 2022, 8:00 PM
 */

#ifndef OGGUIDER_HPP
#define OGGUIDER_HPP

#include <unistd.h>
#include <stdint.h>
#include <stdio.h>
#include "ogDefinitions.hpp"

class ogGuider {
    
public:
    char                ascii2RevComp[256];    
    char                configurationFile[MAX_FILENAME];
    char               *pSequence;
    uint32_t            seqLen;
    char               *pSeqLimit;
    char               *pLimFwd;
    char               *pLimRev;
    char               *pI;
    char               *pCurrent;
    char                finish;
    uint32_t            count;
    uint32_t            countOverLen;
    uint32_t            countTotal;
    uint16_t            keySizeInChars;
    char                isSymmetric;
    uint8_t             fixedGuideLength;
    char               *pComplement;
    uint32_t            complementSize;
    char                isComplement;
    char                genomeMode;         // is the encoding in genome mode (this is only for indexing).
    uint16_t            maxLenNoPattern;
    uint16_t            fixedStep;          // = maxLenNoPattern / fixedGuideLength
    
                             ogGuider();
    virtual ogGuider        *clone();
    virtual void             setKeySizeInChars(uint16_t keySize);
    virtual uint16_t         getKeySizeInChars();
    virtual void             initializeForIndexing();
    virtual uint32_t         load(FILE *pFile);
    virtual uint32_t         save(FILE *pFile);
    virtual void             setSequence(char *pSeq, uint32_t sLen);
    virtual void             fixSequenceCase();
    virtual void             setSequenceAsReverseComplement(char *pSeqComp);
    virtual void             prepareAfterNewSequence();
    virtual char             nextGuide();
    virtual char            *getCurrentGuide();
    virtual uint32_t         getCurrentGuidePosition();
    virtual uint32_t         countGuides();
    virtual                 ~ogGuider();
    virtual const char      *getShortExtensionName();
    virtual const char      *getName();
    virtual char             getIsSymmetric();
    virtual void             setIsSymmetric(char sym);
    virtual void             allocateComplement(uint32_t size);
    virtual void             setGenomeMode(char mode);
    virtual void             setConfigFile(char *pFileName);
    virtual char            *getConfigFile();
    
private:

};

#endif /* OGGUIDER_HPP */

