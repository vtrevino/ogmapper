/* 
 * File:   ogReadKeyMapping.cpp
 * Author: victortrevino
 * 
 * Created on May 21, 2022, 12:52 PM
 */

#include <string.h>
#include <math.h>
#include "ogDefinitions.hpp"
#include "ogReadKeyMapping.hpp"
#include "ogCandidatePosManager.hpp"
#include "ogKeys.hpp"

ogReadKeyMapping::ogReadKeyMapping(uint32_t maxReadLen, uint16_t threadNum) {
    /*
    pKeys = pMapParameters->pKeys;
    pKeyEncoding = pMapParameters->pKeyEncoding;
    pGenome = pMapParameters->pGenome;
    pGenPos = pMapParameters->pGenPos;
    */
    maxReadLength = maxReadLen;
    /**
    keyForward = (uint32_t *) malloc(maxReadLen * sizeof(uint32_t));
    keyReverse = (uint32_t *) malloc(maxReadLen * sizeof(uint32_t));
    keyPosFwd = (uint32_t *) malloc(maxReadLen * sizeof(uint32_t));
    keyPosRev = (uint32_t *) malloc(maxReadLen * sizeof(uint32_t));
    keyPosFwd = (uint32_t *) malloc(maxReadLen * sizeof(uint32_t));
    keyPosRev = (uint32_t *) malloc(maxReadLen * sizeof(uint32_t));
    keyInfoFwd = (ogKey **) malloc(maxReadLen * sizeof(ogKey *));
    keyInfoRev = (ogKey **) malloc(maxReadLen * sizeof(ogKey *));
    activeKeyFwd  = (char *) malloc(maxReadLen * sizeof(char));
    activeKeyRev  = (char *) malloc(maxReadLen * sizeof(char));
     **/
    pForwardKeys = (ogKeyMap *) malloc(maxReadLen * 2 * sizeof(ogKeyMap));
    pReverseKeys = pForwardKeys + maxReadLen;
    
    pCandPosMan = new ogCandidatePosManager(0);
    forceUpperCase = 0;
    thread = threadNum;
    //relativeFwdPos = relativeRevPos = 0;
    debug = 0;
    sharedMemoryBytes = 0;
    pMemory = NULL;
    pLinkedRKM = NULL;
    isMemoryReleasable = true;
    
    //packLenAllocated = 0;
    //packedReadSeqFwd[0] = NULL;
    //packedReadSeqRev[0] = NULL;
    
    isRead2 = 0; // default readd 1
    DNA = NULL;
    dnaAlloc = 0;
    allocDNA(BUFFER_FOR_UNPACKING_DNA);
}

ogReadKeyMapping::~ogReadKeyMapping() {
    //fprintf(stderr, "<Deallocating ogReadKeyMapping:");
    //fprintf(stderr, "<~ogRKM:"); fflush(stderr);
    if (pCandPosMan != NULL) delete pCandPosMan;
    /**
    free(keyForward);
    free(keyReverse);
    free(keyPosFwd);
    free(keyPosRev);
    free(keyInfoFwd);
    free(keyInfoRev);
    free(activeKeyFwd);
    free(activeKeyRev);
     **/
    free(pForwardKeys); // pReverseKeys is only a offset
    if (pMemory != NULL && isMemoryReleasable) {
        free(pMemory);
    }
    if (dnaAlloc > 0) free(DNA);
    //if (packLenAllocated > 0) free(packedReadSeqFwd[0]);
    //fprintf(stderr, ":ogRKM>");
}

void ogReadKeyMapping::allocDNA(uint32_t len) {
    if (len > dnaAlloc) {
        dnaAlloc    = len * 2;
        DNA         = (char *) realloc(DNA, sizeof(char) * dnaAlloc);
    }    
}

void ogReadKeyMapping::setMasterObjects(ogKeys *pK, ogGenome *pG, ogKeyEncoding *pKE,ogGenomePositions *pGP, ogGuider *pGuide, char forceUpper, ogSamWriter *pSamW, ogMappingParameters *pMapPar) {
    pKeys = pK;
    pKeyEncoding = pKE;
    pGenome = pG;
    pGenPos = pGP;
    pGuider = pGuide;
    keySize = pKE->getSizeInChars();
    forceUpperCase = forceUpper;
    pSamWriter = pSamW;
    pMapParams = pMapPar;
}

void ogReadKeyMapping::setStatistics(float meanLog10, float sdLog10, uint32_t q95) {
    meanKeys_log10 = meanLog10;
    sdKeys_log10 = sdLog10;
    uint32_t msd = (uint32_t) pow(10, meanLog10 + sdLog10*4);
    //q95 *= 2;
    if (msd > q95) {
        nKeysLow = q95;
        nKeysHigh = msd;
    } else {
        nKeysLow = msd;
        nKeysHigh = q95;        
    }
}

uint32_t ogReadKeyMapping::getLowKeyCountLimit() {
    return nKeysLow;
}

uint32_t ogReadKeyMapping::getHighKeyCountLimit() {
    return nKeysHigh;
}

uint32_t ogReadKeyMapping::getMeanKeyPlusStdDev(float nStdDev) {
    double n_log10 = meanKeys_log10 + nStdDev * sdKeys_log10;
    double n = pow(10, n_log10);
    return (uint32_t) n;
}

void ogReadKeyMapping::setRead(ogSingleRead *theRead) {
    
    //seq = p;
    //seqLen = len;
    completed = 0;
    read = theRead;
    mapped = 0;
    //relativeFwdPos = relativeRevPos = 0;

    /**
    uint32_t toAlloc = len * 4 + 4;
    if (toAlloc > packLenAllocated) {
        uint32_t len2 = len >> 1;
        packLenAllocated = toAlloc;
        packedReadSeqFwd[0] = (char *) realloc(packedReadSeqFwd[0], packLenAllocated);
        packedReadSeqFwd[1] = packedReadSeqFwd[0] + len2 + 1;
        packedReadSeqFwd[2] = packedReadSeqFwd[1] + len2 + 1;
        packedReadSeqFwd[3] = packedReadSeqFwd[2] + len2 + 1;
        packedReadSeqRev[0] = packedReadSeqFwd[3] + len2 + 1;
        packedReadSeqRev[1] = packedReadSeqRev[0] + len2 + 1;
        packedReadSeqRev[2] = packedReadSeqRev[1] + len2 + 1;
        packedReadSeqRev[3] = packedReadSeqRev[2] + len2 + 1;
    }
    
    pGenome->packSequence(seq, packedReadSeqFwd[0], 0);
    pGenome->packSequence(seq, packedReadSeqFwd[1], 1);
    pGenome->packSequence(seq, packedReadSeqFwd[2], 2);
    pGenome->packSequence(seq, packedReadSeqFwd[3], 3);
    pGenome->packSequence(pRevComp, packedReadSeqRev[0], 0);
    pGenome->packSequence(pRevComp, packedReadSeqRev[1], 1);
    pGenome->packSequence(pRevComp, packedReadSeqRev[2], 2);
    pGenome->packSequence(pRevComp, packedReadSeqRev[3], 3);
     **/
    
    pGuider->setSequence(read->pSeq, read->lenSeq);
    if (forceUpperCase) pGuider->fixSequenceCase();
    uint32_t i = 0, k;
    uint32_t pos;
    uint32_t size;
    char *gF, *gR;
    nFwdPos = nRevPos = 0;
    /**
    if (0 && pGuider->getIsSymmetric()) { // this can not be done because CACAC does not generate the same behaviour from 5'->3' than to 3'->5'
        while (pGuider->nextGuide()) {
            gF = pGuider->getCurrentGuide();
            i = pGuider->getCurrentGuidePosition();
            if (gF <= pGuider->pLimFwd) {
                keyForward[nFwdPos]  = pKeyEncoding->getFwdKey(gF);
                keyPosFwd[nFwdPos++] = i;
            }
            gR = gF - pGuider->fixedGuideLength - 1;
            if (gR >= pGuider->pLimRev) {
                keyReverse[nRevPos]  = pKeyEncoding->getRevKey(gR);
                keyPosRev[nRevPos++] = i - pGuider->fixedGuideLength - 1;
            }
        }
    } else {
     **/
    
    ogKeyMap *pKMap = pForwardKeys, *pKMi;
    pFwdFirst = pFwdLast = pForwardKeys;
    minFwdKeySize = 0x7FFFFFFF;
    maxFwdKeySize = 0;
    minFwdKeySizePos = 0;
    maxFwdKeySizePos = 0;
    while (pGuider->nextGuide()) {
        //fprintf(stderr, "Current Guide:%s",pGuider->getCurrentGuide());
        k = pKeyEncoding->getFwdKey(pGuider->getCurrentGuide());
        //fprintf(stderr, ", Pos=%u, isValid=%c\n",pGuider->getCurrentGuidePosition(), (pKeyEncoding->isValidKey() ? '1' : '0'));
        if (pKeyEncoding->isValidKey()) {
            pos = pGuider->getCurrentGuidePosition();
            pKMap->position = pos;
            pKMap->active = 1;
            pKMap->keyInfo = *(pKeys->getInfoForKey(k));
            pKMap->next = NULL; // no next
            size = pKMap->keyInfo.size;
            if (size < minFwdKeySize && size > 0) {
                minFwdKeySize = pKMap->keyInfo.size;
                minFwdKeySizePos = nFwdPos;
            }
            if (size > maxFwdKeySize) {
                maxFwdKeySize = pKMap->keyInfo.size;
                maxFwdKeySizePos = nFwdPos;            
            }
            pKMap->index = nFwdPos++;
    //        keyForward  [nFwdPos] = k;
    //        keyPosFwd   [nFwdPos] = pos;
    //        keyInfoFwd  [nFwdPos] = pKeys->getInfoForKey(k);
    //        activeKeyFwd[nFwdPos] = 1;
            
            // Search & Sort
            if (size > 0 && pKMap != pForwardKeys) {
                if (size < pFwdFirst->keyInfo.size) {
                    pKMap->next = pFwdFirst;
                    pFwdFirst = pKMap;
                } else {
                    for (pKMi = pFwdFirst; pKMi->next != NULL && pKMi->next->keyInfo.size <= size; pKMi = pKMi->next);
                    pKMap->next = pKMi->next;
                    pKMi->next = pKMap;
                }
                if (pFwdLast->keyInfo.size < size) {
                    pFwdLast = pKMap;
                }
                //fprintf(stderr, "================   NUEVO ORDEN  FORWARD  ========================\n");
                //for (pKMi = pFwdFirst; pKMi != NULL; pKMi = pKMi->next) {
                //    fprintf(stderr, "Curr=%p, Pos=%u, Size=%u, Next=%p\n", pKMi, pKMi->position, pKMi->keyInfo.size, pKMi->next);
                //}
            }
            //
            pKMap++;
            if (nFwdPos >= maxReadLength || pos >= read->lenSeq) {
                fprintf(stderr, "****** PROBLEMS with Fwd Guides %u / %u / %u / %c / %u *****\n", nFwdPos, maxReadLength, read->lenSeq, *pGuider->getCurrentGuide(), pGuider->getCurrentGuidePosition());
                fprintf(stderr, "%s\n", read->pSeq);
            }
        }
    }
    nActiveFwd = nFwdPos;
    //memset(activeKeyFwd, 1, nFwdPos);
    pGuider->setSequenceAsReverseComplement(read->pSeqRevComp);
    pKMap = pReverseKeys;
    pRevFirst = pRevLast = pReverseKeys;
    minRevKeySize = 0x7FFFFFFF;
    maxRevKeySize = 0;
    minRevKeySizePos = 0;
    maxRevKeySizePos = 0;
    while (pGuider->nextGuide()) {
        // Fwd here is correct because sequence was reversed+complemented
        k = pKeyEncoding->getFwdKey(pGuider->getCurrentGuide()); 
        if (pKeyEncoding->isValidKey()) {
            pos = pGuider->getCurrentGuidePosition();
            pKMap->position = pos;
            pKMap->active = 1;
            pKMap->keyInfo = *(pKeys->getInfoForKey(k));
            pKMap->next = NULL; // no next
            size = pKMap->keyInfo.size;
            if (size < minRevKeySize && size > 0) {
                minRevKeySize = pKMap->keyInfo.size;
                minRevKeySizePos = nRevPos;
            }
            if (size > maxRevKeySize) {
                maxRevKeySize = pKMap->keyInfo.size;
                maxRevKeySizePos = nRevPos;            
            }
            pKMap->index = nRevPos++;
    //        keyReverse  [nRevPos] = k;
    //        keyPosRev   [nRevPos] = pos;
    //        keyInfoRev  [nRevPos] = pKeys->getInfoForKey(k);
    //        activeKeyRev[nRevPos] = 1;
            // Search & Sort
            if (size > 0 && pKMap != pReverseKeys) {
                if (size < pRevFirst->keyInfo.size) {
                    pKMap->next = pRevFirst;
                    pRevFirst = pKMap;
                } else {
                    for (pKMi = pRevFirst; pKMi->next != NULL && pKMi->next->keyInfo.size <= size; pKMi = pKMi->next);
                    pKMap->next = pKMi->next;
                    pKMi->next = pKMap;
                }
                if (pRevLast->keyInfo.size < size) {
                    pRevLast = pKMap;
                }
                //fprintf(stderr, "================   NUEVO ORDEN  REVERSE  ========================\n");
                //for (pKMi = pRevFirst; pKMi != NULL; pKMi = pKMi->next) {
                //    fprintf(stderr, "Curr=%p, Pos=%u, Size=%u, Next=%p\n", pKMi, pKMi->position, pKMi->keyInfo.size, pKMi->next);
                //}
            }
            //
            pKMap++;
            if (nRevPos >= maxReadLength || pos >= read->lenSeq) {
                fprintf(stderr, "****** PROBLEMS with Rev Guides %u / %u / %u / %c / %u *****\n", nRevPos, maxReadLength, read->lenSeq, *pGuider->getCurrentGuide(), pGuider->getCurrentGuidePosition());
                fprintf(stderr, "%s\n", read->pSeq);
            }
        }
    }
    nActiveRev = nRevPos;
    //pGuider->setSequenceAsReverseComplement();        

    completed = 100;
}


void ogReadKeyMapping::clearFwdKeyablePosition(uint32_t pos) {
    // depracated
    //keyPosFwd[pos] = 0;
}

void ogReadKeyMapping::clearRevKeyablePosition(uint32_t pos) {
    // depracated
    //keyPosRev[pos] = 0;
}

uint32_t ogReadKeyMapping::getNFwdPos() {
    return nFwdPos;
}

uint32_t ogReadKeyMapping::getNRevPos() {
    return nRevPos;
}

uint32_t ogReadKeyMapping::getNKeyPos(char isReverse) {
    return (isReverse ? nRevPos : nFwdPos);
}


uint32_t ogReadKeyMapping::getFwdKeyablePosition(uint32_t iKey) {
    //return keyPosFwd[iKey];
    return pForwardKeys[iKey].position;
}

uint32_t ogReadKeyMapping::getRevKeyablePosition(uint32_t iKey) {
    //return keyPosRev[iKey];
    return pReverseKeys[iKey].position;
}

uint32_t ogReadKeyMapping::getKeyablePosition(char RC, uint32_t iKey) {
    if (RC) 
        //return keyPosRev[iKey];
        return pReverseKeys[iKey].position;
    else 
        //return keyPosFwd[iKey];
        return pForwardKeys[iKey].position;
}


ogKeyMap *ogReadKeyMapping::getKeyMapForPosition(char isRev, uint32_t pos) {
    if (isRev) 
        return pReverseKeys+pos;
    else 
        return pForwardKeys+pos;
}


ogKey* ogReadKeyMapping::getFwdInfoForKeyablePosition(uint32_t iKey) {
    //return keyInfoFwd[iKey]; //pKeys->getInfoForKey(keyForward[iKey]);
    return &(pForwardKeys[iKey].keyInfo);
}

ogKey* ogReadKeyMapping::getRevInfoForKeyablePosition(uint32_t iKey) {
    //return keyInfoRev[iKey]; // pKeys->getInfoForKey(keyReverse[iKey]);
    return &(pReverseKeys[iKey].keyInfo);
}

ogKey* ogReadKeyMapping::getInfoForPositionKey(char isReverse, uint32_t p1) {
    if (isReverse)
        //return keyInfoRev[p1];//pKeys->getInfoForKey(keyReverse[p1]);
        return &(pReverseKeys[p1].keyInfo);
    else 
        //return keyInfoFwd[p1];//pKeys->getInfoForKey(keyForward[p1]);        
        return &(pForwardKeys[p1].keyInfo);
}


void ogReadKeyMapping::setPosReadKey(uint32_t pos) {
    // DEPRACATED
    if (completed != 100) {
        // pKeyEncoding->getFRKeys(seq + pos, keyForward[pos], keyReverse[pos]);
        // MAPPING NOT YET IMPLEMENTED IN NextKey VERSION
        // DEPRACATED: nKeys++;
    }
}

void ogReadKeyMapping::setFirstReadKey() {
    // DEPRACATED
    if (completed != 100) {
        setPosReadKey(0);
    }
}

uint32_t ogReadKeyMapping::setLastReadKey() {
    // DEPRACATED
    uint32_t pos = getLastKeyablePosition();
    if (completed != 100) {
        setPosReadKey(pos);
    }
    return pos;
}

uint32_t ogReadKeyMapping::getLastKeyablePosition() {
    return read->lenSeq - pKeyEncoding->getSizeInChars();
}

void ogReadKeyMapping::setAllReadKeys() {
    // DEPRACATED:
    return;
    /**
    if (completed != 100) {
        setFirstReadKey();
        uint32_t kF = keyForward[0], kR = keyReverse[0];
        char *p = seq + pKeyEncoding->getSizeInChars();
        uint16_t k;
        uint16_t kk = (seqLen > maxReadLength ? maxReadLength : seqLen);
        for (k=1; k < kk; k++) {
            //pKeyEncoding->shiftFRKeys(*p++, kF, kR);
            // Mapping Not Implemented in NextKey version
            keyForward[k] = kF;
            keyReverse[k] = kR;
        }
        // DEPRACATED: nKeys += kk;
        nFwdPos = kk;
        nRevPos = kk;
        completed = 100;
    }
     **/
}

void ogReadKeyMapping::setReadKeysFromToBy(uint32_t from, uint32_t to, uint32_t by) {
    // DEPRACATED
    return;
    /**
    if (completed != 100) {
        for (; from < to; from += by) {
            setPosReadKey(from);
            // DEPRACATED: nKeys++;
        }
        completed = 100;
    }
     **/
}

ogCandidatePosManager* ogReadKeyMapping::getCandPosMan() {
    return pCandPosMan;
}

ogKey* ogReadKeyMapping::getInfoForKey(uint32_t k1) {
    return pKeys->getInfoForKey(k1);
}

//uint32_t ogReadKeyMapping::intersectByTwoReadPositionKeys(char isReverse, uint32_t p1, uint32_t p2, uint16_t maxDif) {
//    if (isReverse)
//        return intersectByTwoReadKeys(keyReverse[p1], keyReverse[p2], maxDif);
//    else
//        return intersectByTwoReadKeys(keyForward[p1], keyForward[p2], maxDif);
//}
//
uint32_t ogReadKeyMapping::intersectByTwoReadPositionKeys(char isReverse, uint32_t p1, uint32_t p2, uint16_t maxDif) {
    /**
    ogKey *pK1 = (isReverse ? keyInfoRev[p1] : keyInfoFwd[p1]);
    ogKey *pK2 = (isReverse ? keyInfoRev[p2] : keyInfoFwd[p2]);
    if (isReverse) {
        relativeRevPos = keyPosRev[p1 < p2 ? p1 : p2];
    } else {
        relativeFwdPos = keyPosFwd[p1 < p2 ? p1 : p2];
    }
    return pCandPosMan->setByInsersectingKeys(
                    pGenPos->getPointerPosition(pK1->offset), pK1->size, 
                    pGenPos->getPointerPosition(pK2->offset), pK2->size,
                    maxDif);
     **/
    ogKeyMap *pkm1 = (isReverse ? pReverseKeys : pForwardKeys) + p1;
    ogKeyMap *pkm2 = (isReverse ? pReverseKeys : pForwardKeys) + p2;
    uint32_t r;
    r = pCandPosMan->setByInsersectingKeys(
            pGenPos->getPointerPosition(pkm1->keyInfo.offset), pkm1->keyInfo.size, pkm1->position,
            pGenPos->getPointerPosition(pkm2->keyInfo.offset), pkm2->keyInfo.size, pkm2->position, 
            maxDif, isReverse);
    return r;
}

uint32_t ogReadKeyMapping::intersectByTwoReadKeys(uint32_t k1, uint32_t k2, uint16_t maxDif) {
    // depracated 
    // but kept for compatibility before release    
    ogKey *pK1 = pKeys->getInfoForKey(k1);
    ogKey *pK2 = pKeys->getInfoForKey(k2);
    return pCandPosMan->setByInsersectingKeys(
                    pGenPos->getPointerPosition(pK1->offset), pK1->size, 0, 
                    pGenPos->getPointerPosition(pK2->offset), pK2->size, 0,
                    maxDif, 0);
}


uint32_t ogReadKeyMapping::intersectByTwoKeys(ogKey *pK1, ogKey *pK2, uint16_t maxDif, char isReverse) {
    return pCandPosMan->setByInsersectingKeys(
                    pGenPos->getPointerPosition(pK1->offset), pK1->size, 0, 
                    pGenPos->getPointerPosition(pK2->offset), pK2->size, 0,
                    maxDif, isReverse);
}

uint32_t ogReadKeyMapping::intersectByTwoKeys(ogKey *pK1, ogKey *pK2, uint32_t rdPos1, uint32_t rdPos2, uint16_t maxDif, char isReverse) {
    return pCandPosMan->setByInsersectingKeys(
                    pGenPos->getPointerPosition(pK1->offset), pK1->size, rdPos1, 
                    pGenPos->getPointerPosition(pK2->offset), pK2->size, rdPos2,
                    maxDif, isReverse);
}

uint32_t ogReadKeyMapping::intersectByTwoKeysLeft(ogKey *pK1, ogKey *pK2, uint32_t rdPos1, uint16_t maxDif, char isReverse, uint16_t maxAdd) {
    return pCandPosMan->setByInsersectingKeysLeftKey(
                    pGenPos->getPointerPosition(pK1->offset), pK1->size, rdPos1, 
                    pGenPos->getPointerPosition(pK2->offset), pK2->size,
                    maxDif, isReverse, maxAdd);
}

void ogReadKeyMapping::printKeyInfo(uint32_t k) {
    ogKey *pK = pKeys->getInfoForKey(k);
    pGenPos->printGenomePositionsInfo(pK->offset, pK->size, pGenome);
}

uint32_t* ogReadKeyMapping::getGenomePositionFromReadPositionKey(char isReverse, uint32_t p1, uint32_t iPos) {
    //ogKey *pK = pKeys->getInfoForKey(isReverse ? keyReverse[p1] : keyForward[p1]);
    //return pGenPos->getPointerPosition(pK->offset) + iPos;
    return pGenPos->getPointerPosition(((isReverse ? pReverseKeys : pForwardKeys) + p1)->keyInfo.offset) + iPos;
}

uint32_t ogReadKeyMapping::getKeyForward(uint32_t pos) {
    //return keyForward[pos];
    // depracated
    return 0;
}

uint32_t ogReadKeyMapping::getKeyReverse(uint32_t pos) {
    //return keyReverse[pos];
    // depracated
    return 0;
}

uint32_t ogReadKeyMapping::getKeyFR(uint32_t pos, char isRev) {
    // depracated
    // return (isRev ? keyReverse[pos] : keyForward[pos]);
    return 0;
}

uint32_t ogReadKeyMapping::intersectByAddingReadPositionKey(char isReverse, uint32_t p1, uint16_t maxDif) {
    //uint32_t newPos = *((isReverse ? keyPosRev : keyPosFwd) + p1);
//    if (isReverse)
//        return intersectByAddingReadKey(keyReverse[p1], maxDif);
//    else
//        return intersectByAddingReadKey(keyForward[p1], maxDif);
    //ogKey *pKey = &((isReverse ? pReverseKeys[p1] : pForwardKeys[p1]).keyInfo);
    //return pCandPosMan->intersectAddingKey(pGenPos->getPointerPosition(pKey->offset), pKey->size, maxDif);
    ogKeyMap *pKM = (isReverse ? pReverseKeys : pForwardKeys) + p1;
    return pCandPosMan->intersectAddingKey(pGenPos->getPointerPosition(pKM->keyInfo.offset), pKM->keyInfo.size, pKM->position, maxDif, isReverse);
}

uint32_t ogReadKeyMapping::intersectByAddingReadKey(uint32_t k1, uint16_t maxDif) {
    // depracated 
    // but kept for compatibility before release    
//    ogKey *pK1 = pKeys->getInfoForKey(k1);
//    return pCandPosMan->intersectAddingKey(pGenPos->getPointerPosition(pK1->offset), pK1->size, maxDif);
    return 0;
}

//uint32_t ogReadKeyMapping::intersectByAddingReadPositionKeyPositiveDistance(char isReverse, uint32_t p1, uint16_t maxDif) {
//   if (isReverse)
//       return intersectByAddingReadKeyPositiveDistance(keyReverse[p1], maxDif);
//   else
//       return intersectByAddingReadKeyPositiveDistance(keyForward[p1], maxDif);    
//}

//uint32_t ogReadKeyMapping::intersectByAddingReadKeyPositiveDistance(uint32_t k1, uint16_t maxDif) {
//    ogKey *pK1 = pKeys->getInfoForKey(k1);
//    return pCandPosMan->intersectAddingKeyPositiveDistance(pGenPos->getPointerPosition(pK1->offset), pK1->size, maxDif);
//}


void ogReadKeyMapping::inactivePositionFwd(uint32_t p) {
//    if (activeKeyFwd[p]) {
//        activeKeyFwd[p] = 0;
//        nActiveFwd--;
//    }
    if (pForwardKeys[p].active) {
        pForwardKeys[p].active = 0;
        nActiveFwd--;
    }
}

void ogReadKeyMapping::inactivePositionRev(uint32_t p) {
//    if (activeKeyRev[p]) {
//        activeKeyRev[p] = 0;
//        nActiveRev--;
//    }
    if (pReverseKeys[p].active) {
        pReverseKeys[p].active = 0;
        nActiveRev--;
    }
}


void ogReadKeyMapping::inactivePositionsAroundKeyFwd(uint32_t p) {
//    if (activeKeyFwd[p] == 0) {
//        return;
//    }
////    if (activeKeyFwd[p]) {
////        activeKeyFwd[p] = 0;
////        nActiveFwd--;
////    }
////    return;
//    uint16_t from = keyPosFwd[p];
//    //uint16_t to   = from + keySize;
//    uint16_t i, j;
//    for (i=0; i < nFwdPos; i++) {
//        if (activeKeyFwd[i]) {
//            j = keyPosFwd[i];
//            if ((j > from ? j - from : from - j) < keySize) {
//                activeKeyFwd[i] = 0;
//                nActiveFwd--;
//            }
//        }
//    }
    if (pForwardKeys[p].active) {
        uint16_t from = pForwardKeys[p].position;
        //uint16_t to   = from + keySize;
        uint16_t i, j;
        for (i=0; i < nFwdPos; i++) {
            if (pForwardKeys[i].active) {
                j = pForwardKeys[i].position;
                if ((j > from ? j - from : from - j) < keySize) {
                    pForwardKeys[i].active = 0;
                    nActiveFwd--;
                }
            }
        }        
    }
}

void ogReadKeyMapping::inactivePositionsAroundKeyRev(uint32_t p) {
//    if (activeKeyRev[p] == 0) {
//        return;
//    }
////    if (activeKeyRev[p]) {
////        activeKeyRev[p] = 0;
////        nActiveRev--;
////    }
////    return;
//    uint16_t from = keyPosRev[p];
//    //uint16_t to   = from + keySize - 1;
//    uint16_t i, j;
//    for (i=0; i < nRevPos; i++) {
//        if (activeKeyRev[i]) {
//            j = keyPosRev[i];
//            if ((j > from ? j - from : from - j) < keySize) {
//                activeKeyRev[i] = 0;
//                nActiveRev--;
//            }
//        }
//    }
    
    if (pReverseKeys[p].active) {
        uint16_t from = pReverseKeys[p].position;
        //uint16_t to   = from + keySize;
        uint16_t i, j;
        for (i=0; i < nRevPos; i++) {
            if (pReverseKeys[i].active) {
                j = pReverseKeys[i].position;
                if ((j > from ? j - from : from - j) < keySize) {
                    pReverseKeys[i].active = 0;
                    nActiveRev--;
                }
            }
        }        
    }
    
}

void ogReadKeyMapping::reactivateKeys() {
    uint32_t i;
    if (nActiveFwd < nFwdPos) {
        //memset(activeKeyFwd, 1, nFwdPos);
        for (i=0; i < nFwdPos; i++) pForwardKeys[i].active = 1;
        nActiveFwd = nFwdPos;
    }
    if (nActiveRev < nRevPos) {
        //memset(activeKeyRev, 1, nRevPos);
        for (i=0; i < nRevPos; i++) pReverseKeys[i].active = 1;
        nActiveRev = nRevPos;
    }
}

char ogReadKeyMapping::isFwdKeyActive(uint32_t p) {
    //return activeKeyFwd[p];
    return pForwardKeys[p].active;
}

char ogReadKeyMapping::isRevKeyActive(uint32_t p) {
    //return activeKeyRev[p];
    return pReverseKeys[p].active;
}

// uint32_t ogReadKeyMapping::getRelativePos(char isReverse) {
//     if (isReverse) {
//        return relativeRevPos;
//     } else {
//        return relativeFwdPos;
//     }
// }

void ogReadKeyMapping::printGenomeInfoFromKeyPos(uint32_t pos) {
    pGenPos->printGenomePositionsInfo(getFwdInfoForKeyablePosition(pos)->offset, getFwdInfoForKeyablePosition(pos)->size, pGenome);
}

uint32_t ogReadKeyMapping::getFwdTargetSize(uint32_t iKey) {
    return pForwardKeys[iKey].keyInfo.size;
}

uint32_t ogReadKeyMapping::getRevTargetSize(uint32_t iKey) {
    return pReverseKeys[iKey].keyInfo.size;    
}

uint32_t ogReadKeyMapping::getTargetSize(char RC, uint32_t iKey) {
    if (RC) return pReverseKeys[iKey].keyInfo.size;
    return pForwardKeys[iKey].keyInfo.size;
}


void ogReadKeyMapping::printNotFound(uint64_t readNum) {
    fprintf(stderr, "----------- NOT FOUND %llu ------------\n%s\n", readNum, read->pSeq);
    fprintf(stderr, "nFwdPos=%u, nRevPos=%u\n", nFwdPos, nRevPos);
    fprintf(stderr, "minFwdKeySize=%u, maxFwdKeySize=%u\n", minFwdKeySize, maxFwdKeySize);
    fprintf(stderr, "minFwdKeySizePos=%u, maxFwdKeySizePos=%u\n", minFwdKeySizePos, maxFwdKeySizePos);
    fprintf(stderr, "minRevKeySize=%u, maxRevKeySize=%u\n", minRevKeySize, maxRevKeySize);
    fprintf(stderr, "minRevKeySizePos=%u, maxRevKeySizePos=%u\n", minRevKeySizePos, maxRevKeySizePos);

    int i, j, k;
    
    fprintf(stderr, "Forward Key sizes:");
    for (i=0; i < nFwdPos; i++) {
        fprintf(stderr, "%u ", pForwardKeys[i].keyInfo.size);
    }
    fprintf(stderr, "\n");
    fprintf(stderr, "Reverse Key sizes:");
    for (i=0; i < nRevPos; i++) {
        fprintf(stderr, "%u ", pReverseKeys[i].keyInfo.size);
    }
    fprintf(stderr, "\n");
    fprintf(stderr, "Forward Key pos  :");
    for (i=0; i < nFwdPos; i++) {
        fprintf(stderr, "%u ", pForwardKeys[i].position);
    }
    fprintf(stderr, "\n");
    fprintf(stderr, "Reverse Key pos  :");
    for (i=0; i < nRevPos; i++) {
        fprintf(stderr, "%u ", pReverseKeys[i].position);
    }
    fprintf(stderr, "\n");
}


uint32_t ogReadKeyMapping::getCloserGenomePosition(uint32_t targetGPos, uint32_t *pa, uint32_t aSize) {

    if (aSize < 2) return *pa;
    
    uint32_t *p = pa; // always set to lowest explored value
    uint32_t ga = *p;
    uint32_t aMx = aSize;
    uint32_t i;
    uint32_t d1, d2;

    while (aMx > 2) { // ga != targetGPos && 
        i = aMx >> 1;
        ga = *(p+i);
        if (ga < targetGPos) {
            p += i;
            aMx -= i;
        } else if (ga > targetGPos) {
            aMx = i+1; // +1 para que este disponible de 0 hasta aMx
        } else {
            // Hit
            return ga;
        }
    }
    d1 = (targetGPos > *p ? targetGPos - *p : *p - targetGPos);
    pa = p + 1;
    d2 = (targetGPos > *pa ? targetGPos - *pa : *pa - targetGPos);
    return (d1 <= d2 ? *p : *pa);
    //if ((targetGPos-*p) < (*(p+1)-targetGPos)) {
    //    return *p;
    //} else {
    //    return *(p+1);
    //}
}


// return the position in 0-based index
uint32_t ogReadKeyMapping::getCloserGenomeRelativePosition(uint32_t targetGPos, uint32_t *pa, uint32_t aSize) {
    if (aSize < 2) return 0;
    
    uint32_t *p = pa, *pb; // always set to lowest explored value
    uint32_t ga = *p;
    uint32_t aMx = aSize;
    uint32_t i;
    uint32_t d1, d2;

    while (aMx > 2) { // ga != targetGPos && 
        i = aMx >> 1;
        ga = *(p+i);
        if (ga < targetGPos) {
            p += i;
            aMx -= i;
        } else if (ga > targetGPos) {
            aMx = i+1; // +1 para que este disponible de 0 hasta aMx
        } else {
            // Hit
            return (p+i - pa);
        }
    }
    d1 = (targetGPos > *p ? targetGPos - *p : *p - targetGPos);
    pb = p + 1;
    d2 = (targetGPos > *pb ? targetGPos - *pb : *pb - targetGPos);
    return ((d1 <= d2 ? p : pb) - pa) ;
}

void ogReadKeyMapping::resetByOrderSize(char RC) {
    if (RC) pRevK = pRevFirst;
    else pFwdK = pFwdFirst;
}

uint32_t ogReadKeyMapping::getIndexByOrderSize(char RC) {
    return (RC ? pRevK->index : pFwdK->index);
}

char ogReadKeyMapping::nextByOrderSize(char RC) {
    if (RC) {
        pRevK = pRevK->next;
        if (pRevK == NULL) return 0;
    } else {
        pFwdK = pFwdK->next;
        if (pFwdK == NULL) return 0;
    }
    return 1;
}

char ogReadKeyMapping::finishByOrderSize(char RC) {
    if (RC) return (pRevK == NULL ? 1 : 0);
    else return (pFwdK == NULL ? 1 : 0);
}

ogKey *ogReadKeyMapping::keyByOrderSize(char RC) {
    return (RC ? (pRevK == NULL ? NULL : &pRevK->keyInfo) : (pFwdK == NULL ? NULL : &pFwdK->keyInfo));
}

void *ogReadKeyMapping::getUsableMemory(uint64_t bytes) {
    if (bytes > sharedMemoryBytes) {
        uint64_t n = bytes - sharedMemoryBytes;
        //fprintf(stderr, "<ogRKM: Allocating %.1f Mb>", (float) n / (1024*1024)); fflush(stderr);
        if (pMemory == NULL) {
            // First time, make sure it if full of zeros
            pMemory = calloc(bytes, 1);
        } else {
            pMemory = realloc(pMemory, bytes);
        }
        sharedMemoryBytes = bytes;
        if (pLinkedRKM != NULL) {
            // Estrategia para ahorrar memoria, 
            // Hace que los 2 RdKyMap compartan la misma memory (supone que solo 1 thread opera los 2 RdKyMap)
            pLinkedRKM->pMemory = pMemory;
            pLinkedRKM->sharedMemoryBytes = sharedMemoryBytes;
        }
    }
    return pMemory;
}