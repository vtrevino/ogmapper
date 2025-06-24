
#include <stdio.h>
#include <map>
#include "ogDefinitions.hpp"
#include "ogIndex.hpp"
#include "ogKeyEncoding.hpp"
#include "ogReadKeyMapping.hpp"
#include "ogReusableBinaryTree.hpp"
#include "ogReadKeyMapping.hpp"
#include "ogCandidatePosManager.hpp"
#include "binnedBigCounters.hpp"



// BLAST LIKE (ORIGINAL)
uint16_t _original_blastLike(ogReadKeyMapping *pMap, uint32_t maxTargets) {

    char        x, acc = 0;
    uint32_t    nKeys, i, iRdPos;
    ogKey       *pKey;
    uint32_t    *pGenPos, bestGenPos, bestKeyPos, minError;
    uint32_t    size;
    uint32_t    rdLen = pMap->read->lenSeq;
    uint32_t    errores, maxErr = MAX_SIMPLE_ALIGN_ERR(rdLen);
    
    pMap->pCandPosMan->reset();
    pMap->read->packIfNeeded(pMap->pGenome);
    for (x = 0; x < 2; x++) {
        nKeys = (x == 0 ? pMap->getNFwdPos() : pMap->getNRevPos());
        if (nKeys > 0) {
            //char **pPacked = (x == 0 ? pMap->read->packedReadSeqFwd : pMap->read->packedReadSeqRev);
            pMap->resetByOrderSize(x);
            minError = maxErr+1;
            for (; pMap->finishByOrderSize(x) == 0 && acc < 99; pMap->nextByOrderSize(x)) {
                i = pMap->getIndexByOrderSize(x);
                pKey = pMap->getInfoForPositionKey(x, i);
                if ((size = pKey->size) < maxTargets) {                
                    iRdPos = pMap->getKeyablePosition(x, i);
                    pGenPos = pMap->getGenomePositionFromReadPositionKey(x, i, 0);
                    acc = 0;
                    char offset = -1;
                    //fprintf(stderr, "ReadNum=%llu, len=%u, RC=%c, size=%u, ", pMap->seqNum, rdLen, x+48, size);
                    while (size) {
                        if (iRdPos < *pGenPos) errores = pMap->pGenome->comparePackedSequences(pMap->read->packReadSeqIfNeeded(pMap->pGenome->getPackedOffset(*pGenPos - iRdPos), x), *pGenPos - iRdPos, rdLen, minError);
                        if (errores <= minError) {
                            if (errores < minError) {
                                minError = errores;
                                acc = (rdLen - minError) * 100 / rdLen;
                                pMap->pCandPosMan->restoreFwdCandidatePositions(); // reset counter for fwd or rev
                                pMap->pCandPosMan->checkToAddGenomeAndReadPosition(*pGenPos, iRdPos, x, acc);
                            } else {
                                // (errores == minError)
                                if (pMap->pCandPosMan->checkToAddGenomeAndReadPosition(*pGenPos, iRdPos, x, acc) == 2) break;
                            }
                        }
                        pGenPos++;
                        size--;
                    }
                    //fprintf(stderr, "offset=%c, acc=%u\n", offset+48, acc);
                } else {
                    break;
                }
            }
        }
        if (x == 0) {
            pMap->pCandPosMan->saveFwdCandidatePositions();
        }
    }
    //pMap->pCandPosMan->printWithGenomePositions(pMap->pGenome);
    //if (pMap->read->readNumber == 231) {
    //    fprintf(stderr, "This.count=%d\n",pMap->pCandPosMan->getCount()); fflush(stderr);
    //}

    return pMap->pCandPosMan->getCount();
}

// BLAST LIKE : alternating between Fwd/Rev by minimum key size (it is 33% faster than the "original" algorithm)
uint16_t _bad_blastLike(ogReadKeyMapping *pMap, uint32_t maxTargets) {

    char        x, acc = 0;
    uint32_t    nKeys, i, iRdPos;
    ogKey       *pKeyF, *pKeyR, *pKey;
    uint32_t    *pGenPos, bestGenPos, bestKeyPos, minError;
    uint32_t    size,j=0;
    uint32_t    rdLen = pMap->read->lenSeq;
    uint32_t    errores, maxErr = MAX_ERROR_ACCEPT1_ALIGN(rdLen); //MAX_ERROR_ACCEPT1_PACKED
    
    pMap->pCandPosMan->reset();
    pMap->read->packIfNeeded(pMap->pGenome);
    for (x=0; x < 2; x++) {
        nKeys = (x == 0 ? pMap->getNFwdPos() : pMap->getNRevPos());
        minError = maxErr+1;
        // 99 tiene mucho impacto, si se cambia a 95 se tarda menos de la mitad del tiempo
        for (i=0; acc < 99 && i < nKeys; i++) {
            pKey = pMap->getInfoForPositionKey(x, i);
            if ((size = pKey->size) < maxTargets) {
                iRdPos = pMap->getKeyablePosition(x, i);
                pGenPos = pMap->getGenomePositionFromReadPositionKey(x, i, 0);
                acc = 0;
                // Este es solo para no comparar coordenadas problematicas en ~0
                for(; *pGenPos < iRdPos && size; pGenPos++, size--);
                while (size) {
                    errores = pMap->pGenome->comparePackedSequences(pMap->read->packReadSeqIfNeeded(pMap->pGenome->getPackedOffset(*pGenPos - iRdPos), x), *pGenPos - iRdPos, rdLen, minError);
                    if (errores <= minError) {
                        if (errores < minError) {
                            minError = errores;
                            acc = (rdLen - minError) * 100 / rdLen;
                            pMap->pCandPosMan->reset();
                            pMap->pCandPosMan->addGenomeAndReadPosition(*pGenPos, iRdPos, x, acc);
                        } else {
                            // (errores == minError)
                            pMap->pCandPosMan->checkToAddGenomeAndReadPosition(*pGenPos, iRdPos, x, acc);
                        }
                    }
                    pGenPos++;
                    size--;
                }
                //fprintf(stderr, "offset=%c, acc=%u\n", offset+48, acc);
            }
        }
    }
    return pMap->pCandPosMan->getCount();
}

uint16_t blastLike(ogReadKeyMapping *pMap, uint32_t maxTargets, uint32_t minTargets, char minAcc) {

    char        x, acc = 0;
    uint32_t    nKeys, i, iRdPos;
    ogKey       *pKeyF, *pKeyR, *pKey;
    uint32_t    *pGenPos, bestGenPos, bestKeyPos, minError;
    uint32_t    size,j=0;
    uint32_t    rdLen = pMap->read->lenSeq;
    uint32_t    errores, maxErr = MAX_ERROR_ACCEPT1_ALIGN(rdLen); //MAX_ERROR_ACCEPT1_PACKED //MAX_ERROR_ACCEPT1_ALIGN
    ogSingleRead *pRead = pMap->read;
    ogGenome     *pGenome = pMap->pGenome;
    ogCandidatePosManager *pCandPosMan = pMap->pCandPosMan;
    
    pCandPosMan->reset();
    pRead->packIfNeeded(pGenome);
    pMap->resetByOrderSize(0);
    pMap->resetByOrderSize(1);
    minError = maxErr+1;
    // 99 tiene mucho impacto, si se cambia a 95 se tarda menos de la mitad del tiempo
    while (acc < minAcc && (pMap->finishByOrderSize(0) == 0 || pMap->finishByOrderSize(1) == 0) ) {
        pKeyF = pMap->keyByOrderSize(0);
        pKeyR = pMap->keyByOrderSize(1);
        if (pKeyF == NULL) x = 1;
        else if (pKeyR == NULL) x = 0;
        else x = (pKeyF->size <= pKeyR->size ? 0 : 1);
        //char **pPacked = (x == 0 ? pMap->read->packedReadSeqFwd : pMap->read->packedReadSeqRev);
        pKey = pMap->keyByOrderSize(x);
        if ((size = pKey->size) < maxTargets) {
            if (size >= minTargets) {
                i = pMap->getIndexByOrderSize(x);        
                iRdPos = pMap->getKeyablePosition(x, i);
                pGenPos = pMap->getGenomePositionFromReadPositionKey(x, i, 0);
                acc = 0;
                // Este es solo para no comparar coordenadas problematicas en ~0
                for(; *pGenPos < iRdPos && size; pGenPos++, size--);
                while (size) {
                    //errores = pMap->pGenome->comparePackedSequences(pMap->read->packReadSeqIfNeeded(pMap->pGenome->getPackedOffset(*pGenPos - iRdPos), x), *pGenPos - iRdPos, rdLen, minError);
                    errores = pGenome->comparePackedSequences(pRead->packReadSeqIfNeeded(pGenome->getPackedOffset(*pGenPos - iRdPos), x), *pGenPos - iRdPos, rdLen, minError);
                    if (errores <= minError) {
                        if (errores < minError) {
                            minError = errores;
                            acc = (rdLen - minError) * 100 / rdLen;
                            pCandPosMan->reset();
                            pCandPosMan->addGenomeAndReadPosition(*pGenPos, iRdPos, x, acc);
                        } else {
                            // (errores == minError)
                            pCandPosMan->checkToAddGenomeAndReadPosition(*pGenPos, iRdPos, x, acc);
                        }
                    }
                    pGenPos++;
                    size--;
                }
            }
            //fprintf(stderr, "offset=%c, acc=%u\n", offset+48, acc);
        } else {
            break;
        }
        pMap->nextByOrderSize(x);
    }
    return pCandPosMan->getCount();
}

uint16_t blastLikeSmall(ogReadKeyMapping *pMap) {
    return blastLike(pMap,  pMap->getLowKeyCountLimit() >> 1, 0, 98);
}

uint16_t blastLikeLarge(ogReadKeyMapping *pMap) {
    return blastLike(pMap,  pMap->getLowKeyCountLimit(), 0, 99);
}

uint16_t blastLikeHuge(ogReadKeyMapping *pMap) {
    return blastLike(pMap, (pMap->getHighKeyCountLimit()), pMap->getLowKeyCountLimit(), 100);
}

uint16_t blastLikeAll(ogReadKeyMapping *pMap) {
    return blastLike(pMap, 0xFFFFFFFF, 0, 100);
}


// intersect keys from 2 reads
uint16_t intersectTwoReads(ogReadKeyMapping *pMap, uint32_t maxTargets) {

    ogReadKeyMapping *pMap2 = pMap->pLinkedRKM;
    char        x, y, maxAcc = 0, xBest = 'Z';
    uint32_t    nKeys, i, j, iRdPos, jRdPos, c, k, jBest=0, iBest=0;
    ogKey       *pKey, *pKey2;
    uint32_t    *pGenPos, bestGenPos, bestKeyPos, minError;
    uint32_t    size;
    uint32_t    rdLen = pMap->read->lenSeq;
    uint32_t    errores, maxErr = MAX_SIMPLE_ALIGN_ERR(rdLen);
    uint16_t    MAX_INTERSECT = 64;
    uint32_t    lePos;
    ogCandidatePosManager *pCand = pMap->pCandPosMan;
    
    pCand->reset();
    //pMap2->pCandPosMan->reset();

    // Relativeto pMap1
    minError = rdLen;
    
    for (x = 0; x < 2; x++) {
        y = 1 - x;
        nKeys = (x == 0 ? pMap->getNFwdPos() : pMap->getNRevPos());
        if (nKeys > 0) {
            //char **pPacked = (x == 0 ? pMap->read->packedReadSeqFwd : pMap->read->packedReadSeqRev);
            pMap->resetByOrderSize(x);
            for (; pMap->finishByOrderSize(x) == 0 && maxAcc < 99; pMap->nextByOrderSize(x)) {
                i = pMap->getIndexByOrderSize(x);
                pKey = pMap->getInfoForPositionKey(x, i);
                if (pKey->size < maxTargets) {
                    iRdPos = pMap->getKeyablePosition(x, i);
                    //fprintf(stderr, "/i=%u",i);
                    pMap2->resetByOrderSize(y);                    
                    for (; pMap2->finishByOrderSize(y) == 0 && maxAcc < 99; pMap2->nextByOrderSize(y)) {
                        j = pMap2->getIndexByOrderSize(y);
                        //fprintf(stderr, "/j=%u",j);
                        pKey2 = pMap2->getInfoForPositionKey(y, j);
                        if (pKey->size < maxTargets) {
                            jRdPos = pMap2->getKeyablePosition(y, j);
                            pCand->rewindCountsAfterAddingKeyCountZero();
                            c = pCand->count;
                            if (pMap->intersectByTwoKeysLeft(pKey, pKey2,iRdPos, MAX_PAIR_INSERT_SIZE-iRdPos+jRdPos, x, MAX_INTERSECT)) {
                                if ((pCand->count - c) < MAX_INTERSECT) {
                                    for (k=c; k < pCand->count && maxAcc < 99; k++) {
                                        lePos = pCand->getkPos(k)->genomePosition - iRdPos;
                                        errores = pMap->pGenome->comparePackedSequences(pMap->read->packReadSeqIfNeeded(pMap->pGenome->getPackedOffset(lePos), x), lePos, rdLen, minError);
                                        if (errores < minError) {
                                            maxAcc = (rdLen-errores)*100/rdLen;
                                            minError = errores;
                                            iBest = i;
                                            jBest = j;
                                            xBest = x;
                                            //fprintf(stderr, "/i=%u,j=%u,C=%u,e=%u", i,j,pCand->count,minError);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        pCand->rewindCountsAfterAddingKeyCountZero();
        if (x == xBest) {
            pKey = pMap->getInfoForPositionKey(xBest, iBest);
            iRdPos = pMap->getKeyablePosition(xBest, iBest);
            pKey2 = pMap2->getInfoForPositionKey(1-xBest, jBest);
            jRdPos = pMap2->getKeyablePosition(1-xBest, jBest);
            pMap->intersectByTwoKeysLeft(pKey, pKey2, iRdPos, MAX_PAIR_INSERT_SIZE-iRdPos+jRdPos, xBest, MAX_INTERSECT);
            //fprintf(stderr, "/Count=%u",pCand->count);
            //pCand->printWithGenomePositions(pMap->pGenome);
        }
        if (x == 0) {
            pCand->saveFwdCandidatePositions();
        }
    }
    //pCand->printWithGenomePositions(pMap->pGenome);
    return pCand->getCount();
}


uint16_t Intersect2R_Small(ogReadKeyMapping *pMap) {
    return intersectTwoReads(pMap,  pMap->getLowKeyCountLimit());
}

uint16_t Intersect2R_Large(ogReadKeyMapping *pMap) {
    return intersectTwoReads(pMap,  pMap->getHighKeyCountLimit());
}

uint16_t Intersect2R_Huge(ogReadKeyMapping *pMap) {
    return intersectTwoReads(pMap, 10*pMap->getHighKeyCountLimit());
}

uint16_t Intersect2R_All(ogReadKeyMapping *pMap) {
    return intersectTwoReads(pMap, 0xFFFFFFFF);
}



/// ESTE ES MUY LENTO EN EL LONG RUN
uint16_t blastLikeBinary(ogReadKeyMapping *pMap) {

    char        x, acc;
    uint32_t    nKeys, i, iRdPos, k;
    int32_t     jNode;
    ogKey       *pKey;
    uint32_t    *pGenPos, bestGenPos, bestKeyPos, minError;
    uint32_t    size;
    uint32_t    rdLen = pMap->read->lenSeq;
    uint32_t    errores, maxErr = MAX_SIMPLE_ALIGN_ERR(rdLen);
    ogReusableBinaryTree *pTree = pMap->pCandPosMan->pTree; // memAlloc[pMap->thread];
    
    
    pMap->pCandPosMan->reset();
    for (x = 0; x < 2; x++) {
        nKeys = (x == 0 ? pMap->getNFwdPos() : pMap->getNRevPos());
        if (nKeys > 0) {
            
            // Build tree from last largest node
            pTree->clear();
            pTree->setMaximumDistanceForKey(3);
            ogKeyMap *pLast = (x == 0 ? pMap->pFwdLast : pMap->pRevLast);            
            uint32_t *pGP = pMap->getGenomePositionFromReadPositionKey(x,pLast->index,0);
            uint32_t L = pLast->keyInfo.size;
            iRdPos = pMap->getKeyablePosition(x, pLast->index);
            pTree->allocateIfNeeded(L);
            for (i=L; i; ) {
                i--;
                pTree->setNodeTo(i, pGP[i]-iRdPos, -1);
            }
            pTree->buildTreeFromOrderedNodesSet(L);
            
            //char **pPacked = (x == 0 ? pMap->read->packedReadSeqFwd : pMap->read->packedReadSeqRev);
            pMap->resetByOrderSize(x);
            minError = maxErr;
            for (; pMap->finishByOrderSize(x) == 0; pMap->nextByOrderSize(x)) {
                i = pMap->getIndexByOrderSize(x);
                iRdPos = pMap->getKeyablePosition(x, i);
                pKey = pMap->getInfoForPositionKey(x, i);
                size = pKey->size;
                pGenPos = pMap->getGenomePositionFromReadPositionKey(x, i, 0);
                acc = 0;
                char offset = -1;
                fprintf(stderr, "ReadNum=%llu, len=%u, RC=%c, size=%u, pTree.Size=%d, ", pMap->read->readNumber, rdLen, x+48, size, pTree->nUsed);
                while (size) {
                    k = *pGenPos - iRdPos;
                    jNode = pTree->searchForKey(k);
                    //fprintf(stderr, "/%d ",jNode);
                    if (jNode < 0 || pTree->getNodeContent(jNode) < 0) {
                        errores = pMap->pGenome->comparePackedSequences(pMap->read->packReadSeqIfNeeded(pMap->pGenome->getPackedOffset(k), x), k, rdLen, minError+1);
                        if (jNode < 0) pTree->insertNode(k, errores); 
                        else pTree->setNodeContent(jNode, errores);
                    } else {
                        errores = pTree->getNodeContent(jNode);
                    }
                    if (errores <= minError) {
                        if (errores < minError) {
                            offset = k & 0x00000003;
                            minError = errores;
                            acc = (rdLen - minError) * 100 / rdLen;
                            pMap->pCandPosMan->restoreFwdCandidatePositions(); // reset counter for fwd or rev
                        }
                        if (errores <= minError) {
                            if (pMap->pCandPosMan->checkToAddGenomeAndReadPosition(*pGenPos, iRdPos, x, acc) == 2) break;
                        }
                    }
                    pGenPos++;
                    size--;
                }
                fprintf(stderr, "pTree.Size=%d, offset=%c, acc=%u\n", pTree->nUsed, offset+48, acc);
                if (acc > 98) {
                    break;
                }
            }
        }
        if (x == 0) {
            pMap->pCandPosMan->saveFwdCandidatePositions();
        }
    }
    pMap->pCandPosMan->printWithGenomePositions(pMap->pGenome);
    return pMap->pCandPosMan->getCount();
}



// Algorithm - Divide in 5 pieces, take 2nd and 4th
uint16_t ExactInApartSides(ogReadKeyMapping *pMap) {    
    uint32_t l = pMap->getLastKeyablePosition();
    uint32_t p1 = l/5;                                                          // Parameter 1: 5=How many parts should the read be split into?
    uint32_t p2 = l - p1;
    uint32_t d = p2 - p1 + 2;                                                   // Parameter 2: 5=How many nt of theoretical difference should be considered as small indels?
    uint32_t p3 = 0;
    char RC = 0;
    pMap->setPosReadKey(p1);
    pMap->setPosReadKey(p2);
    uint32_t intCount = pMap->intersectByTwoReadPositionKeys(0, p1, p2, d);
    uint32_t oriCount = intCount;
    if (intCount == 0) {
        intCount = pMap->intersectByTwoReadPositionKeys(1, p1, p2, d);
        oriCount = intCount;
        RC = 1;
    }
    pMap->pCandPosMan->saveFwdCandidatePositions();
    if (intCount > 10) {
        p3 = p1 + ((p2 - p1) >> 1);
        d = p3 - p1 + 5;
        pMap->setPosReadKey(p3);
        intCount = pMap->intersectByAddingReadPositionKey(RC, p3, d);
        if (intCount > 10) {
            p3 = p2 - ((p2 - p1) >> 1);
            d = p3 - p1 + 5;
            pMap->setPosReadKey(p3);
            intCount = pMap->intersectByAddingReadPositionKey(RC, p3, d);            
        }
    }
    /**
    if (intCount > 0 && 0) {
        fprintf(stderr, "\n====================================================\n");
        fprintf(stderr, "Read: %s\n", pMap->seq);
        fprintf(stderr, "%u Genome Positions.\nRead KeyPos 1 %u [%*s]\nRead KeyPos 2 %u [%*s]\n", 
                intCount, p1, pMap->pKeyEncoding->getSizeInChars(), pMap->seq + p1, p2, pMap->pKeyEncoding->getSizeInChars(), pMap->seq + p2);
        int i;
        for (i=0; i < intCount; i++) {
            fprintf(stderr, "Pos %d, Gen Pos %u, ",i+1, pMap->pCandPosMan->getkPos(i));
            ogChromosome *pChr = pMap->pGenome->getGenomicCoordinate(pMap->pCandPosMan->getkPos(i));
            fprintf(stderr, "Chr %s Pos %u\n", pChr->name, pMap->pCandPosMan->getkPos(i) - pChr->start);
        }
    }
    if (intCount > 0) {
        int i = 0;
            fprintf(stderr, "Pos %d, Gen Pos %u, ",i+1, pMap->pCandPosMan->getkPos(i));
            ogChromosome *pChr = pMap->pGenome->getGenomicCoordinate(pMap->pCandPosMan->getkPos(i));
            fprintf(stderr, "Chr %s Pos %u, [original counts %u, RC=%hhu] ", pChr->name, pMap->pCandPosMan->getkPos(i) - pChr->start, oriCount, RC);
    }
    if (intCount == 0) {
        fprintf(stderr, "\n====================================================\n");
        fprintf(stderr, "Read: %s\n", pMap->seq);
        fprintf(stderr, "%u Genome Positions.\nRead KeyPos 1 %u [%*s]\nRead KeyPos 2 %u [%*s]\n", 
                intCount, p1, pMap->pKeyEncoding->getSizeInChars(), pMap->seq + p1, p2, pMap->pKeyEncoding->getSizeInChars(), pMap->seq + p2);
        ogKey *pK = pMap->pKeys->getInfoForKey(pMap->getKeyForward(p1));
        fprintf(stderr, "Key 1 Forward: offset=%u, size=%u\n",pK->offset,pK->size);
        pK = pMap->pKeys->getInfoForKey(pMap->getKeyForward(p2));
        fprintf(stderr, "Key 2 Forward: offset=%u, size=%u\n",pK->offset,pK->size);
        pK = pMap->pKeys->getInfoForKey(pMap->getKeyReverse(p1));
        fprintf(stderr, "Key 1 Reverse: offset=%u, size=%u\n",pK->offset,pK->size);
        pK = pMap->pKeys->getInfoForKey(pMap->getKeyReverse(p2));
        fprintf(stderr, "Key 2 Reverse: offset=%u, size=%u\n",pK->offset,pK->size);
        fprintf(stderr, "Original intersection count:%u, RC=%hhu\n", oriCount, RC);
        if (p3 > 0) {
            fprintf(stderr, "Read KeyPos 3 %u [%s]\n", p3, pMap->seq + p3);
            pK = pMap->pKeys->getInfoForKey(pMap->getKeyForward(p3));
            fprintf(stderr, "Key 3 Forward: offset=%u, size=%u\n",pK->offset,pK->size);            
            pK = pMap->pKeys->getInfoForKey(pMap->getKeyReverse(p3));
            fprintf(stderr, "Key 3 Reverse: offset=%u, size=%u\n",pK->offset,pK->size);            
        }
    }
     **/
    return (uint16_t) intCount;
}

// Algorithm - QUICKEST! Exact Contiguous in extremes
uint16_t ExactInExtremes(ogReadKeyMapping *pMap) {
    uint32_t l = pMap->setLastReadKey();
    uint32_t d = l + 5;                                                   // Parameter 2: 5=How many nt of theoretical difference should be considered as small indels?
    uint32_t p3 = 0;
    char RC = 0;
    pMap->setFirstReadKey();
    uint32_t intCount = pMap->intersectByTwoReadKeys(pMap->getKeyForward(0), pMap->getKeyForward(l), d);
    uint32_t oriCount = intCount;
    if (intCount == 0) {
        intCount = pMap->intersectByTwoReadKeys(pMap->getKeyReverse(0), pMap->getKeyReverse(l), d);
        oriCount = intCount;
        RC = 1;
    }
    pMap->pCandPosMan->saveFwdCandidatePositions();
    if (intCount > 10) {
        p3 = (l >> 1);
        d = p3 + 5;
        pMap->setPosReadKey(p3);
        if (RC) 
            intCount = pMap->intersectByAddingReadKey(pMap->getKeyReverse(p3), d);
        else
            intCount = pMap->intersectByAddingReadKey(pMap->getKeyForward(p3), d);
    }

    return (uint16_t) intCount;
}



// Algorithm - Let 5p Keys
uint16_t Two5pKeys(ogReadKeyMapping *pMap) {
    uint32_t keysize = pMap->pKeyEncoding->getSizeInChars(); // se puede usar pMap->keySize
    uint32_t maxl = pMap->getLastKeyablePosition();
    uint32_t p2 = keysize << 1;                            // Parameter: offset to enlarge key
    uint32_t p1 = 0;
    uint32_t d = p2 - p1 + 5;                                                   // Parameter 2: 5=How many nt of theoretical difference should be considered as small indels?
    uint32_t p3 = 0;
    if (p2 > maxl) p2 = maxl;
    char RC = 0;
    pMap->setPosReadKey(p1);
    pMap->setPosReadKey(p2);
    uint32_t intCount = pMap->intersectByTwoReadKeys(pMap->getKeyForward(p1), pMap->getKeyForward(p2), d);
    if (intCount > 10) {
        p3 = p1 + ((p2 - p1) >> 1);
        d = p3 - p1 + 5;
        pMap->setPosReadKey(p3);
        intCount = pMap->intersectByAddingReadKey(pMap->getKeyForward(p3), d);
    }
    uint32_t intCountR = pMap->intersectByTwoReadKeys(pMap->getKeyReverse(p1), pMap->getKeyReverse(p2), d);
    if (intCountR > 10) {
        p3 = p1 + ((p2 - p1) >> 1);
        d = p3 - p1 + 5;
        pMap->setPosReadKey(p3);
        intCountR = pMap->intersectByAddingReadKey(pMap->getKeyReverse(p3), d);
    }
    RC = intCountR >= intCount ? 1 : 0;
    pMap->pCandPosMan->saveFwdCandidatePositions();
    return (RC ? intCountR : intCount);
}


uint16_t SweepKeys(ogReadKeyMapping *pMap, uint32_t startPos, uint32_t step) {
    //uint32_t keysize = pMap->pKeyEncoding->getSizeInChars();
    uint32_t maxl = pMap->getLastKeyablePosition();
    uint32_t p1,p2,p0;
    uint32_t d;
    uint32_t intCount;
    uint32_t maxForwardCounts = 0;
    uint32_t maxReverseCounts = 0;
    uint32_t maxCount;
    char RC;
    char first;
    if (step < 5) pMap->setAllReadKeys();
    else pMap->setPosReadKey(startPos);
    for (RC=0; RC < 2; RC++) {
        p0 = p1 = startPos;
        intCount = 2;
        maxCount = 0;
        for (first = 1; intCount > 1 && p1 < maxl; p1 = p2, first = 0) {
            p2 = p1 + step;
            if (p2 > maxl) p2 = maxl;
            if (RC == 0) pMap->setPosReadKey(p2);
            d = p2 - p0 + 2;
            if (first) 
                intCount = pMap->intersectByTwoReadPositionKeys(RC, p1, p2, d);
            else 
                intCount = pMap->intersectByAddingReadPositionKey(RC, p2, d);
            if (intCount > 0) maxCount = intCount;
        }
        if (RC) {
            maxReverseCounts = maxCount;
        } else {
            maxForwardCounts = maxCount;
        }
        
    }
    pMap->pCandPosMan->saveFwdCandidatePositions();
    return (maxReverseCounts > maxForwardCounts ? maxReverseCounts : maxForwardCounts);
}

uint16_t NonOverlappingKeys(ogReadKeyMapping *pMap) {
    return SweepKeys(pMap, 0, pMap->pKeyEncoding->getSizeInChars());
}


uint16_t OverlappingBy3nt(ogReadKeyMapping *pMap) {
    return SweepKeys(pMap, 0, 3);
}

uint16_t OverlappingBy6nt(ogReadKeyMapping *pMap) {
    return SweepKeys(pMap, 0, 6);
}

uint16_t OverlappingBy9nt(ogReadKeyMapping *pMap) {
    return SweepKeys(pMap, 0, 9);
}

uint16_t OverlappingBy12nt(ogReadKeyMapping *pMap) {
    return SweepKeys(pMap, 0, 12);
}

uint16_t OverlappingBy15nt(ogReadKeyMapping *pMap) {
    return SweepKeys(pMap, 0, 15);
}

uint16_t OverlappingBy20nt(ogReadKeyMapping *pMap) {
    return SweepKeys(pMap, 0, 20);
}



// Use an histogram scheme but then splitting positions into more specific, thus only 1 Key round is needed
// Genomic position is rooted to first NT of the read
// Counts : BIN array having 0-255 contain the # of keys pointing to that bin, INITIALIZED TO 0
// PosRef : BIN array of Genomic Position of the first Position of this bin to speed up calculations
// ConRef : BIN array of counts of keys "matching" +/- 10 nt to that position to speed up calc
// KeyRef : BIN array of the iKey for the PosRef
// MorRef : BIN array of counts for other positions not being PosRef
uint16_t HistogramMapping(ogReadKeyMapping *pMap, uint32_t maxTargets) {
    struct OtherGenPos {
        uint32_t    genPos;
        uint16_t    iRdPos;
        uint16_t    count;
        uint32_t    nextPosInBin; // 0 means no more
    };
    uint32_t    BINS = 1000, ibin;
    float       fBINS = BINS;
    float       fMaxGenomePos = pMap->pGenome->getGenomeSize();
    float       factorMultiplier =  fBINS / fMaxGenomePos;
    uint16_t    *binCounts; //binCounts[BINS+1];
    uint32_t    *binRefPos; // binRefPos[BINS+1]; // position to real other gen pos
    uint32_t    maxPositions = 65535 * 3;
    OtherGenPos *allGenPos; //allGenPos[maxPositions];
    OtherGenPos *kGenPos;
    uint32_t    allGenPosK;// Next position usable of allGenPos
    
    uint32_t bestKeyPos, bestGenPos, bestCount;
    uint32_t nReadPos, currGenPos;
    uint16_t i, ties, nUsedKeys;
    uint32_t j, k, lastk, kpos, *keyGenPos, *keyGenFinalPos, nTargets, iPos;
    ogChromosome *pChr;
    char    x, bug = pMap->debug, match;
    void    *pV;

    // Get Memory
    uint64_t bytes = (BINS+1)*(sizeof(uint16_t)+sizeof(uint32_t)) + maxPositions * sizeof(OtherGenPos);
    pV = pMap->getUsableMemory(bytes);
    // Redirect pointers
    binCounts = (uint16_t *) pV;
    binRefPos = (uint32_t *) (binCounts + (BINS+1));
    allGenPos = (OtherGenPos *) (binRefPos + (BINS+1));
    //

    pMap->pCandPosMan->reset();
    
    for (x = 0; x < 2; x++) {
        nReadPos = (x == 0 ? pMap->getNFwdPos() : pMap->getNRevPos());
        
        if (nReadPos > 0) {
            // Initialize
            memset(binCounts, 0, (BINS+1)*sizeof(uint16_t));        
            bestKeyPos = 0;
            bestGenPos = 0;
            bestCount = 0;
            ties = 0;
            nUsedKeys = 0;
            allGenPosK = 1; // in 1 to avoid using the 0 and comparisons chains with 0

            for (i=0; i < nReadPos; i++) {
                if ((nTargets=pMap->getTargetSize(x, i)) < maxTargets) {
                    nUsedKeys++;
                    iPos = pMap->getKeyablePosition(x, i);
                    keyGenPos = pMap->getGenomePositionFromReadPositionKey(x, i, 0); // pMap->getInfoForPositionKey(x,i)->offset
                    for (keyGenFinalPos = keyGenPos + nTargets; keyGenPos < keyGenFinalPos; keyGenPos++) {
                        currGenPos = *keyGenPos - iPos; // shift to start of the read (using -iPos) to facilitate referencing and speed up
                        ibin = currGenPos * factorMultiplier;
                        if (ibin < 0) ibin = 0;
                        else if (ibin > BINS) ibin = BINS;

                        if ((++binCounts[ibin]) == 1) {
                            if (allGenPosK < maxPositions) {
                                // FIRST TIME
                                kpos = allGenPosK++;
                                binRefPos[ibin] = kpos;
                                allGenPos[kpos] = { currGenPos, i, 1, 0 };
                            } else {
                                fprintf(stderr, "*** Problems (1), positions overpassed %u *** ", allGenPosK);
                            }
                        } else {
                            if (false) {
                                // New version : do a sorted chain
                                k = binRefPos[ibin];
                                match = 0;
                                lastk = 0;
                                do {
                                    kGenPos = allGenPos + k;
                                    if (currGenPos <= kGenPos->genPos) {
                                        break;
                                    }
                                    lastk = k;
                                    k = kGenPos->nextPosInBin;
                                } while (k != 0);
                                if (k > 0) {
                                    if (currGenPos == kGenPos->genPos || 
                                            (currGenPos > kGenPos->genPos ? currGenPos-kGenPos->genPos < 10 : kGenPos->genPos-currGenPos < 10 )) {
                                        // MATCH!
                                        if (++kGenPos->count > bestCount) {
                                            bestCount = kGenPos->count;
                                            bestKeyPos = kGenPos->iRdPos;
                                            bestGenPos = kGenPos->genPos;
                                        }
                                        match = 1;
                                    } else {
                                        // Chain in the middle or beginning ... k is larger so the new points to k and previous should point to new kpos
                                        if (allGenPosK < maxPositions) {
                                            // CREATE ANOTHER
                                            kpos = allGenPosK++;
                                            allGenPos[kpos] = { currGenPos, i, 1, k };
                                            // Chain
                                            if (lastk == 0) {
                                                binRefPos[ibin] = kpos; // Chain in the beginning 
                                            } else {
                                                allGenPos[lastk].nextPosInBin = kpos; // Chain
                                            }
                                        } else {
                                            fprintf(stderr, "*** Problems (3), positions overpassed %u *** ", allGenPosK);
                                        }
                                    }
                                } else {
                                    // Chain to the end ... 
                                    if (allGenPosK < maxPositions) {
                                        // CREATE ANOTHER
                                        kpos = allGenPosK++;
                                        allGenPos[kpos] = { currGenPos, i, 1, 0 };
                                        allGenPos[lastk].nextPosInBin = kpos; // Chain
                                    } else {
                                        fprintf(stderr, "*** Problems (4), positions overpassed %u *** ", allGenPosK);
                                    }                            
                                }
                            }

                            if (true) {
                                // Check in all positions
                                k = binRefPos[ibin];
                                match = 0;
                                do {
                                    lastk = k;
                                    kGenPos = allGenPos + k;
                                    if (currGenPos == kGenPos->genPos || 
                                            (currGenPos > kGenPos->genPos ? currGenPos-kGenPos->genPos < 10 : kGenPos->genPos-currGenPos < 10 )) {
                                        // MATCH!
                                        if (++kGenPos->count >= bestCount) {
                                            if (kGenPos->count > bestCount) {
                                                bestCount = kGenPos->count;
                                                bestKeyPos = kGenPos->iRdPos;
                                                bestGenPos = kGenPos->genPos;
                                                ties = 0;
                                            } else {
                                                ties++;
                                            }
                                        }
                                        match = 1;
                                        break;
                                    }
                                    k = kGenPos->nextPosInBin;
                                } while (k != 0);
                                if (match == 0) {
                                    if (allGenPosK < maxPositions) {
                                        // CREATE ANOTHER
                                        kpos = allGenPosK++;
                                        allGenPos[kpos] = { currGenPos, i, 1, 0 };
                                        allGenPos[lastk].nextPosInBin = kpos; // Chain
                                    } else {
                                        fprintf(stderr, "*** Problems (2), positions overpassed %u *** ", allGenPosK);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (bestCount > 0 && ties > 0) {
                //fprintf(stderr, "********** TIES %u, BestCount %u, Used Keys %u ***********", ties, bestCount, nUsedKeys);
                if (ties > 3 || (ties > 0 && bestCount * 3 / nUsedKeys < 1)) {
    //                fprintf(stderr, "CLEARED.\n%s\n",pMap->seq);
    //                for (i=0; i < BINS; i++) {
    //                    if (binCounts[i] >= bestCount) {
    //                        // Check
    //                        k = binRefPos[i];
    //                        do {
    //                            kGenPos = allGenPos + k;
    //                            if (kGenPos->count == bestCount) {
    //                                pChr = pMap->pGenome->getGenomicCoordinate(kGenPos->genPos);
    //                                fprintf(stderr, "RC=%hhd, Bin=%hu, Count=%hu, GenPos %u, Chr %u:%u \n", x, i, kGenPos->count, kGenPos->genPos, pChr->number, kGenPos->genPos - pChr->start);
    //                            }
    //                            k = kGenPos->nextPosInBin;
    //                        } while (k != 0);
    //                    }
    //                }
                    bestCount = 0;
                //} else {
                //    fprintf(stderr, "Ok!.\n%s\n",pMap->seq);
                }
            }
            if (bestCount > 0) {
                // If ties > 0 ???
                if (ties > 0) {
                    for (i = 0; i <= BINS; i++) {
                        if (binCounts[i] >= bestCount) {
                            // Check
                            k = binRefPos[i];
                            do {
                                kGenPos = allGenPos + k;
                                if (kGenPos->count == bestCount) {
                                    pMap->pCandPosMan->addGenomeAndReadPosition(kGenPos->genPos + pMap->getKeyablePosition(x, kGenPos->iRdPos), pMap->getKeyablePosition(x, kGenPos->iRdPos), x, bestCount * 100 / (nUsedKeys + ties));
                                    //pChr = pMap->pGenome->getGenomicCoordinate(kGenPos->genPos);
                                    //fprintf(stderr, "RC=%hhd, Bin=%hu, Count=%hu, GenPos %u, Chr %u:%u \n", x, i, kGenPos->count, kGenPos->genPos, pChr->number, kGenPos->genPos - pChr->start);
                                }
                                k = kGenPos->nextPosInBin;
                            } while (k != 0);
                        }
                    }
                } else {
                    pMap->pCandPosMan->addGenomeAndReadPosition(bestGenPos + pMap->getKeyablePosition(x, bestKeyPos), pMap->getKeyablePosition(x, bestKeyPos), x, bestCount * 100 / (nUsedKeys));
                }
            }
        }
        if (x == 0) {
            pMap->pCandPosMan->saveFwdCandidatePositions();
        }
    }

    //fprintf(stderr, "Counts: %u, Fwd:%u, Rev:%u\n", pMap->pCandPosMan->getCount(), pMap->pCandPosMan->getFwdCount(), pMap->pCandPosMan->getRevCount());
    return pMap->pCandPosMan->getCount();

}



// Use an histogram scheme Binned by different levels depending on counts, it uses binnedBigCounts class
uint16_t BinnedHistogramMapping(ogReadKeyMapping *pMap, uint32_t maxTargets) {
    uint32_t    BINS = 1024*8;
    uint32_t    genomeSize = pMap->pGenome->getGenomeSize();
    binnedBigCounters bbc[] = { binnedBigCounters(0, pMap->pGenome->getGenomeSize(), BINS), binnedBigCounters(0, pMap->pGenome->getGenomeSize(), BINS) };
    
    uint32_t bestKeyPos, bestGenPos, bestCount;
    uint32_t nPossibleKeys, currGenPos;
    uint16_t i, ties, nUsedKeys;
    uint32_t j, k, lastk, kpos, *keyGenPos, *keyGenFinalPos, nTargets, iPos;
    ogChromosome *pChr;
    char    x, bug = pMap->debug, match;
    void    *pV;
    unsigned char score;

    // Get Memory
    //

    pMap->pCandPosMan->reset();
    
    for (x = 0; x < 2; x++) {
        nPossibleKeys = (x == 0 ? pMap->getNFwdPos() : pMap->getNRevPos());
        
        if (nPossibleKeys > 0) {
            // Initialize
            nUsedKeys = 0;

            for (i=0; i < nPossibleKeys; i++) {
                if ((nTargets=pMap->getTargetSize(x, i)) < maxTargets) {
                    nUsedKeys++;
                    iPos = pMap->getKeyablePosition(x, i);
                    keyGenPos = pMap->getGenomePositionFromReadPositionKey(x, i, 0); // pMap->getInfoForPositionKey(x,i)->offset
                    for (keyGenFinalPos = keyGenPos + nTargets; keyGenPos < keyGenFinalPos; keyGenPos++) {
                        currGenPos = *keyGenPos - iPos; // shift to start of the read (using -iPos) to facilitate referencing and speed up
                        if (currGenPos < genomeSize) { // iPos is > than *keyGenPos and currGenPos goes "negative" (not handled).
                            bbc[x].incCounter(currGenPos);
                            //fprintf(stderr, "ReadNum=%llu\n", pMap->seqNum);
                        }
                    }
                }
            }
            bestCount = bbc[x].maxCount;
            ties = bbc[x].ties;
            bestGenPos = bbc[x].maxCountItem;
            if (bestCount > 0 && (nUsedKeys > 1 || nUsedKeys == 1 && ties <= ACEPTABLE_CANDIDATE_REGIONS)) {
                //if (bestCount > 0 && ties > 0) {
                //  if (ties > 3 || (ties > 0 && bestCount * 3 / nUsedKeys < 1)) {
                //      bestCount = 0;
                //  }
                //}
                score = bestCount * 80 / (nUsedKeys < 5 ? 5 : nUsedKeys) + nUsedKeys * 20 / nPossibleKeys;
                if (ties > TOO_MANY_CANDIDATE_REGIONS && score < 50 || bestCount < 4 && ties > ACEPTABLE_CANDIDATE_REGIONS || ties > 0 && bestCount < 3 && nUsedKeys >= 6) {
                    bestCount = 0;
                    ties = 0;
                }
                if (bestCount > 0) {
                    if (ties > 0) {
                        //fprintf(stderr, "_______________________________________________________________________________\n");
                        //fprintf(stderr, "Ties=%d, bestCount=%d, KeysInRead=%d, RC=%c, ReadNum=%llu\n%s\n",ties,bestCount,nReadPos,x+48,pMap->seqNum, pMap->seq);
                        bbc[x].startTiedProcess();
                        uint32_t gpos = bbc[x].getNextTiedCounterValue();
                        while (gpos != 0) {
                            pMap->pCandPosMan->addGenomeAndReadPosition(gpos + pMap->getKeyablePosition(x, 0), pMap->getKeyablePosition(x, 0), x, score);
                            gpos = bbc[x].getNextTiedCounterValue();
                        }
                    } else {
                        pMap->pCandPosMan->addGenomeAndReadPosition(bestGenPos + pMap->getKeyablePosition(x, 0), pMap->getKeyablePosition(x, 0), x, score);
                    }
                }
            }
        }
        if (x == 0) {
            pMap->pCandPosMan->saveFwdCandidatePositions();
        }
    }

    //fprintf(stderr, "Counts: %u, Fwd:%u, Rev:%u\n", pMap->pCandPosMan->getCount(), pMap->pCandPosMan->getFwdCount(), pMap->pCandPosMan->getRevCount());
    return pMap->pCandPosMan->getCount();

}




// Use an histogram scheme but then splitting positions into more specific, thus only 1 Key round is needed
// Genomic position is rooted to first NT of the read
// Counts : BIN array having 0-255 contain the # of keys pointing to that bin, INITIALIZED TO 0
// PosRef : BIN array of Genomic Position of the first Position of this bin to speed up calculations
// ConRef : BIN array of counts of keys "matching" +/- 10 nt to that position to speed up calc
// KeyRef : BIN array of the iKey for the PosRef
// MorRef : BIN array of counts for other positions not being PosRef
uint16_t SortedSizesHistogramMapping(ogReadKeyMapping *pMap, uint32_t maxTargets) {
    struct OtherGenPos {
        uint32_t    genPos;
        uint16_t    iRdPos;
        uint16_t    count;
        uint32_t    nextPosInBin; // 0 means no more
    };
    //fprintf(stderr, "**** SortedSizesHistogramMapping ****\n"); fflush(stderr);
    uint32_t    BINS = 1000, ibin;
    float       fBINS = BINS;
    float       fMaxGenomePos = pMap->pGenome->getGenomeSize();
    float       factorMultiplier =  fBINS / fMaxGenomePos;
    uint16_t    *binCounts; //binCounts[BINS+1];
    uint32_t    *binRefPos; //binRefPos[BINS+1]; // position to real other gen pos
    uint32_t    maxPositions = 65536 * 3; // esto aguanta cuando son threads pero luego hay errores: = 8192*1 * 3; // 65536 * 3 esto es el original
    OtherGenPos *allGenPos;//allGenPos[maxPositions];
    OtherGenPos *kGenPos;
    uint32_t    allGenPosK;// Next position usable of allGenPos
    
    uint32_t bestKeyPos, bestGenPos, bestCount;
    uint32_t nReadPos, currGenPos, compGenPos;
    uint16_t i, ties, nUsedKeys;
    uint32_t j, k, lastk, kpos, *keyGenPos, *keyGenFinalPos, nTargets, iPos;
    ogChromosome *pChr;
    char    x, bug = pMap->debug, match;
    void    *pV;

    // Get Memory
    uint64_t bytes = (BINS+1)*(sizeof(uint16_t)+sizeof(uint32_t)) + maxPositions * sizeof(OtherGenPos);
    pV = pMap->getUsableMemory(bytes);
    // Redirect pointers
    binCounts = (uint16_t *) pV;
    binRefPos = (uint32_t *) (binCounts + (BINS+1));
    allGenPos = (OtherGenPos *) (binRefPos + (BINS+1));
    //
    
    pMap->pCandPosMan->reset();
    
    //if (pMap->read->readNumber == 42) {
    //    fprintf(stderr, "Seq 42\n");
    //}
    //if (bug) {
    //    fprintf(stderr, "*** BUG, SeqNum=%llu ***\n", pMap->seqNum);
    //    fprintf(stderr, "%s\n",pMap->seq);
    //}
    
//    if (maxTargets < pMap->minFwdKeySize || maxTargets < pMap->minRevKeySize) {
//        //if (bug) {
//            if ((maxTargets = pMap->minFwdKeySize) < pMap->minRevKeySize) maxTargets = pMap->minRevKeySize;
//            fprintf(stderr, "MinTargets=%u, SeqNum=%llu \n",maxTargets, pMap->seqNum);
//            if (maxTargets > 20000) bug = 1;
//            if (bug) fprintf(stderr, "%s\n",pMap->seq);
//            maxTargets *= 2;
//        //}
//    }
    
    for (x = 0; x < 2; x++) {
        nReadPos = (x == 0 ? pMap->getNFwdPos() : pMap->getNRevPos());
        
        if (nReadPos > 0) {
            // Initialize
            //memset(binCounts, 0, sizeof(binCounts));        
            for (i=0; i <= BINS; i++) binCounts[i] = 0;
            bestKeyPos = 0;
            bestGenPos = 0;
            bestCount = 0;
            ties = 0;
            nUsedKeys = 0;
            allGenPosK = 1; // in 1 to avoid using the 0 and comparisons chains with 0

            pMap->resetByOrderSize(x);

            for (; pMap->finishByOrderSize(x) == 0; pMap->nextByOrderSize(x)) {
                i = pMap->getIndexByOrderSize(x);
                if ((nTargets=pMap->getTargetSize(x, i)) < maxTargets) {
                    nUsedKeys++;
                    iPos = pMap->getKeyablePosition(x, i);
                    keyGenPos = pMap->getGenomePositionFromReadPositionKey(x, i, 0); // pMap->getInfoForPositionKey(x,i)->offset
                    //if (bug) fprintf(stderr, "i=%u, MaxTarg=%u, nUsedKeys=%u, iPos=%u, nTargets=%u\n",i,maxTargets,nUsedKeys,iPos,nTargets);
                    for (keyGenFinalPos = keyGenPos + nTargets; keyGenPos < keyGenFinalPos; keyGenPos++) {
                        currGenPos = *keyGenPos - iPos; // shift to start of the read (using -iPos) to facilitate referencing and speed up
                        compGenPos = currGenPos - 10;
                        ibin = currGenPos * factorMultiplier;
                        if (ibin < 0) ibin = 0;
                        else if (ibin > BINS) ibin = BINS;

                        if ((++binCounts[ibin]) == 1) {
                            if (allGenPosK < maxPositions) {
                                // FIRST TIME
                                kpos = allGenPosK++;
                                binRefPos[ibin] = kpos;
                                allGenPos[kpos] = { currGenPos, i, 1, 0 };
                            } else {
                                fprintf(stderr, "*** Problems (1), positions overpassed %u *** ", allGenPosK);
                            }
                        } else {
                            if (true) {
                                // New version : do a sorted chain
                                k = binRefPos[ibin];
                                match = 0;
                                lastk = 0;
                                do {
                                    kGenPos = allGenPos + k;
                                    if (compGenPos <= kGenPos->genPos) {
                                        break;
                                    }
                                    lastk = k;
                                    k = kGenPos->nextPosInBin;
                                } while (k != 0);
                                if (k > 0) {
                                    if (currGenPos == kGenPos->genPos || 
                                            (currGenPos > kGenPos->genPos ? currGenPos-kGenPos->genPos < 10 : kGenPos->genPos-currGenPos < 10 )) {
                                        // MATCH!
                                        if (++kGenPos->count >= bestCount) {
                                            if (kGenPos->count > bestCount) {
                                                bestCount = kGenPos->count;
                                                bestKeyPos = kGenPos->iRdPos;
                                                bestGenPos = kGenPos->genPos;
                                                ties = 0;
                                            } else {
                                                ties++;
                                                //if (bug) fprintf(stderr, "ties=%u ",ties);
                                            }
                                        }
                                        match = 1;
                                    } else {
                                        // Chain in the middle or beginning ... k is larger so the new points to k and previous should point to new kpos
                                        if (allGenPosK < maxPositions) {
                                            // CREATE ANOTHER
                                            kpos = allGenPosK++;
                                            allGenPos[kpos] = { currGenPos, i, 1, k };
                                            // Chain
                                            if (lastk == 0) {
                                                binRefPos[ibin] = kpos; // Chain in the beginning 
                                            } else {
                                                allGenPos[lastk].nextPosInBin = kpos; // Chain
                                            }
                                        } else {
                                            fprintf(stderr, "*** Problems (3), positions overpassed %u *** ", allGenPosK);
                                        }
                                    }
                                } else {
                                    // Chain to the end ... 
                                    if (allGenPosK < maxPositions) {
                                        // CREATE ANOTHER
                                        kpos = allGenPosK++;
                                        allGenPos[kpos] = { currGenPos, i, 1, 0 };
                                        allGenPos[lastk].nextPosInBin = kpos; // Chain
                                    } else {
                                        fprintf(stderr, "*** Problems (4), positions overpassed %u *** ", allGenPosK);
                                    }                            
                                }
                            }

                            if (false) {
                                // Check in all positions
                                k = binRefPos[ibin];
                                match = 0;
                                do {
                                    lastk = k;
                                    kGenPos = allGenPos + k;
                                    if (currGenPos == kGenPos->genPos || 
                                            (currGenPos > kGenPos->genPos ? currGenPos-kGenPos->genPos < 10 : kGenPos->genPos-currGenPos < 10 )) {
                                        // MATCH!
                                        if (++kGenPos->count >= bestCount) {
                                            if (kGenPos->count > bestCount) {
                                                bestCount = kGenPos->count;
                                                bestKeyPos = kGenPos->iRdPos;
                                                bestGenPos = kGenPos->genPos;
                                                ties = 0;
                                            } else {
                                                ties++;
                                            }
                                        }
                                        match = 1;
                                        break;
                                    }
                                    k = kGenPos->nextPosInBin;
                                } while (k != 0);
                                if (match == 0) {
                                    if (allGenPosK < maxPositions) {
                                        // CREATE ANOTHER
                                        kpos = allGenPosK++;
                                        allGenPos[kpos] = { currGenPos, i, 1, 0 };
                                        allGenPos[lastk].nextPosInBin = kpos; // Chain
                                    } else {
                                        fprintf(stderr, "*** Problems (2), positions overpassed %u *** ", allGenPosK);
                                    }
                                }
                            }
                        }
                    }
                    if (nUsedKeys > 3 && ties == 0 && ((float) bestCount / (float) nUsedKeys) >= 0.74 ) {
                        break;
                    }
                } else {
                    break;
                }
            }
    //        if (bug) {
    //            fprintf(stderr, ">> nUsedKeys=%u, BestCount=%u, ties=%u\n", nUsedKeys, bestCount, ties);
    //            for (i = 0; i <= BINS; i++) {
    //                if (binCounts[i] > 0) {
    //                    // Check
    //                    k = binRefPos[i];
    //                    do {
    //                        kGenPos = allGenPos + k;
    //                        pChr = pMap->pGenome->getGenomicCoordinate(kGenPos->genPos);
    //                        fprintf(stderr, "RC=%hhd, Bin=%hu, BinCount=%u, Node Count=%hu, GenPos %u, Chr %u:%u \n", x, i, binCounts[i], kGenPos->count, kGenPos->genPos, pChr->number, kGenPos->genPos - pChr->start);
    //                        k = kGenPos->nextPosInBin;
    //                    } while (k != 0);
    //                }
    //            }
    //        }
            if (bestCount > 0 && ties > 0) {
                //fprintf(stderr, "::::::: SeqNum %llu, RC: %c, TIES %u, BestCount %u, Used Keys %u ::::::::::\n", pMap->seqNum, x+48, ties, bestCount, nUsedKeys);
                if (ties > TOO_MANY_CANDIDATE_REGIONS || bestCount < 4 && ties > ACEPTABLE_CANDIDATE_REGIONS) {
                    bestCount = 0;
                }
                //if ((ties > 0 && bestCount * 3 / nUsedKeys < 1)) {
    //                fprintf(stderr, "CLEARED.\n%s\n",pMap->seq);
    //                for (i=0; i < BINS; i++) {
    //                    if (binCounts[i] >= bestCount) {
    //                        // Check
    //                        k = binRefPos[i];
    //                        do {
    //                            kGenPos = allGenPos + k;
    //                            if (kGenPos->count == bestCount) {
    //                                pChr = pMap->pGenome->getGenomicCoordinate(kGenPos->genPos);
    //                                fprintf(stderr, "RC=%hhd, Bin=%hu, Count=%hu, GenPos %u, Chr %u:%u \n", x, i, kGenPos->count, kGenPos->genPos, pChr->number, kGenPos->genPos - pChr->start);
    //                            }
    //                            k = kGenPos->nextPosInBin;
    //                        } while (k != 0);
    //                    }
    //                }
                //    bestCount = 0;
                //} else {
                //    fprintf(stderr, "Ok!.\n%s\n",pMap->seq);
                //}
            }
            if (bestCount > 0) {
                // If ties > 0 ???
                if (ties > 0) {
                    for (i = 0; i <= BINS; i++) {
                        if (binCounts[i] >= bestCount) {
                            // Check
                            k = binRefPos[i];
                            do {
                                kGenPos = allGenPos + k;
                                if (kGenPos->count == bestCount) {
                                    pMap->pCandPosMan->addGenomeAndReadPosition(kGenPos->genPos + pMap->getKeyablePosition(x, kGenPos->iRdPos), pMap->getKeyablePosition(x, kGenPos->iRdPos), x, bestCount * 100 / (nUsedKeys + ties));
                                    //pChr = pMap->pGenome->getGenomicCoordinate(kGenPos->genPos);
                                    //fprintf(stderr, "RC=%hhd, Bin=%hu, Count=%hu, GenPos %u, Chr %u:%u \n", x, i, kGenPos->count, kGenPos->genPos, pChr->number, kGenPos->genPos - pChr->start);
                                }
                                k = kGenPos->nextPosInBin;
                            } while (k != 0);
                        }
                    }
                } else {
                    pMap->pCandPosMan->addGenomeAndReadPosition(bestGenPos + pMap->getKeyablePosition(x, bestKeyPos), pMap->getKeyablePosition(x, bestKeyPos), x, bestCount * 100 / (nUsedKeys));
                }
            }
        }
        if (x == 0) {
            pMap->pCandPosMan->saveFwdCandidatePositions();
        }
    }

    //if (pMap->pCandPosMan->getCount() == 0) {
    //    fprintf(stderr, "0000000000000 Counts, ReadNum:%llu\n%s\n",pMap->seqNum, pMap->seq);
    //}
    //fprintf(stderr, "Counts: %u, Fwd:%u, Rev:%u\n", pMap->pCandPosMan->getCount(), pMap->pCandPosMan->getFwdCount(), pMap->pCandPosMan->getRevCount());
    return pMap->pCandPosMan->getCount();

}




// Intersect iterating by key sizes
uint16_t Intersect(ogReadKeyMapping *pMap, uint32_t maxTargets) {
    struct OtherGenPos {
        uint32_t    genPos;
        uint16_t    iRdPos;
        uint16_t    count;
        uint32_t    nextPosInBin; // 0 means no more
    };
    uint32_t    BINS = 1000, ibin;
    float       fBINS = BINS;
    float       fMaxGenomePos = pMap->pGenome->getGenomeSize();
    float       factorMultiplier =  fBINS / fMaxGenomePos;
    uint16_t    binCounts[BINS+1];
    uint32_t    binRefPos[BINS+1]; // position to real other gen pos
    uint32_t    maxPositions = 65535 * 3;
    OtherGenPos allGenPos[maxPositions], *kGenPos;
    uint32_t    allGenPosK;// Next position usable of allGenPos
    
    uint32_t bestKeyPos, bestGenPos, bestCount;
    uint32_t nReadPos, currGenPos;
    uint16_t i, ties, nUsedKeys;
    uint32_t j, k, lastk, kpos, *keyGenPos, *keyGenFinalPos, nTargets, iPos;
    ogChromosome *pChr;
    char    x, bug = pMap->debug, match;

    pMap->pCandPosMan->reset();
    
    //if (bug) {
    //    fprintf(stderr, "*** BUG, SeqNum=%llu ***\n", pMap->seqNum);
    //    fprintf(stderr, "%s\n",pMap->seq);
    //}
    
//    if (maxTargets < pMap->minFwdKeySize || maxTargets < pMap->minRevKeySize) {
//        //if (bug) {
//            if ((maxTargets = pMap->minFwdKeySize) < pMap->minRevKeySize) maxTargets = pMap->minRevKeySize;
//            fprintf(stderr, "MinTargets=%u, SeqNum=%llu \n",maxTargets, pMap->seqNum);
//            if (maxTargets > 20000) bug = 1;
//            if (bug) fprintf(stderr, "%s\n",pMap->seq);
//            maxTargets *= 2;
//        //}
//    }
    
    for (x = 0; x < 2; x++) {
        nReadPos = (x == 0 ? pMap->getNFwdPos() : pMap->getNRevPos());
        
        if (nReadPos > 0) {
            // Initialize
            memset(binCounts, 0, sizeof(binCounts));        
            bestKeyPos = 0;
            bestGenPos = 0;
            bestCount = 0;
            ties = 0;
            nUsedKeys = 0;
            allGenPosK = 1; // in 1 to avoid using the 0 and comparisons chains with 0

            pMap->resetByOrderSize(x);

            for (; pMap->finishByOrderSize(x) == 0; pMap->nextByOrderSize(x)) {
                i = pMap->getIndexByOrderSize(x);
                if ((nTargets=pMap->getTargetSize(x, i)) < maxTargets) {
                    nUsedKeys++;
                    iPos = pMap->getKeyablePosition(x, i);
                    keyGenPos = pMap->getGenomePositionFromReadPositionKey(x, i, 0); // pMap->getInfoForPositionKey(x,i)->offset
                    //if (bug) fprintf(stderr, "i=%u, MaxTarg=%u, nUsedKeys=%u, iPos=%u, nTargets=%u\n",i,maxTargets,nUsedKeys,iPos,nTargets);
                    for (keyGenFinalPos = keyGenPos + nTargets; keyGenPos < keyGenFinalPos; keyGenPos++) {
                        currGenPos = *keyGenPos - iPos; // shift to start of the read (using -iPos) to facilitate referencing and speed up
                        ibin = currGenPos * factorMultiplier;
                        if (ibin < 0) ibin = 0;
                        else if (ibin > BINS) ibin = BINS;

                        if ((++binCounts[ibin]) == 1) {
                            if (allGenPosK < maxPositions) {
                                // FIRST TIME
                                kpos = allGenPosK++;
                                binRefPos[ibin] = kpos;
                                allGenPos[kpos] = { currGenPos, i, 1, 0 };
                            } else {
                                fprintf(stderr, "*** Problems (1), positions overpassed %u *** ", allGenPosK);
                            }
                        } else {
                            if (true) {
                                // New version : do a sorted chain
                                k = binRefPos[ibin];
                                match = 0;
                                lastk = 0;
                                do {
                                    kGenPos = allGenPos + k;
                                    if (currGenPos <= kGenPos->genPos) {
                                        break;
                                    }
                                    lastk = k;
                                    k = kGenPos->nextPosInBin;
                                } while (k != 0);
                                if (k > 0) {
                                    if (currGenPos == kGenPos->genPos || 
                                            (currGenPos > kGenPos->genPos ? currGenPos-kGenPos->genPos < 10 : kGenPos->genPos-currGenPos < 10 )) {
                                        // MATCH!
                                        if (++kGenPos->count >= bestCount) {
                                            if (kGenPos->count > bestCount) {
                                                bestCount = kGenPos->count;
                                                bestKeyPos = kGenPos->iRdPos;
                                                bestGenPos = kGenPos->genPos;
                                                ties = 0;
                                            } else {
                                                ties++;
                                                //if (bug) fprintf(stderr, "ties=%u ",ties);
                                            }
                                        }
                                        match = 1;
                                    } else {
                                        // Chain in the middle or beginning ... k is larger so the new points to k and previous should point to new kpos
                                        if (allGenPosK < maxPositions) {
                                            // CREATE ANOTHER
                                            kpos = allGenPosK++;
                                            allGenPos[kpos] = { currGenPos, i, 1, k };
                                            // Chain
                                            if (lastk == 0) {
                                                binRefPos[ibin] = kpos; // Chain in the beginning 
                                            } else {
                                                allGenPos[lastk].nextPosInBin = kpos; // Chain
                                            }
                                        } else {
                                            fprintf(stderr, "*** Problems (3), positions overpassed %u *** ", allGenPosK);
                                        }
                                    }
                                } else {
                                    // Chain to the end ... 
                                    if (allGenPosK < maxPositions) {
                                        // CREATE ANOTHER
                                        kpos = allGenPosK++;
                                        allGenPos[kpos] = { currGenPos, i, 1, 0 };
                                        allGenPos[lastk].nextPosInBin = kpos; // Chain
                                    } else {
                                        fprintf(stderr, "*** Problems (4), positions overpassed %u *** ", allGenPosK);
                                    }                            
                                }
                            }

                            if (false) {
                                // Check in all positions
                                k = binRefPos[ibin];
                                match = 0;
                                do {
                                    lastk = k;
                                    kGenPos = allGenPos + k;
                                    if (currGenPos == kGenPos->genPos || 
                                            (currGenPos > kGenPos->genPos ? currGenPos-kGenPos->genPos < 10 : kGenPos->genPos-currGenPos < 10 )) {
                                        // MATCH!
                                        if (++kGenPos->count >= bestCount) {
                                            if (kGenPos->count > bestCount) {
                                                bestCount = kGenPos->count;
                                                bestKeyPos = kGenPos->iRdPos;
                                                bestGenPos = kGenPos->genPos;
                                                ties = 0;
                                            } else {
                                                ties++;
                                            }
                                        }
                                        match = 1;
                                        break;
                                    }
                                    k = kGenPos->nextPosInBin;
                                } while (k != 0);
                                if (match == 0) {
                                    if (allGenPosK < maxPositions) {
                                        // CREATE ANOTHER
                                        kpos = allGenPosK++;
                                        allGenPos[kpos] = { currGenPos, i, 1, 0 };
                                        allGenPos[lastk].nextPosInBin = kpos; // Chain
                                    } else {
                                        fprintf(stderr, "*** Problems (2), positions overpassed %u *** ", allGenPosK);
                                    }
                                }
                            }
                        }
                    }
                    if (nUsedKeys > 3 && ties == 0 && ((float) bestCount / (float) nUsedKeys) >= 0.74 ) {
                        break;
                    }
                }
            }
    //        if (bug) {
    //            fprintf(stderr, ">> nUsedKeys=%u, BestCount=%u, ties=%u\n", nUsedKeys, bestCount, ties);
    //            for (i = 0; i <= BINS; i++) {
    //                if (binCounts[i] > 0) {
    //                    // Check
    //                    k = binRefPos[i];
    //                    do {
    //                        kGenPos = allGenPos + k;
    //                        pChr = pMap->pGenome->getGenomicCoordinate(kGenPos->genPos);
    //                        fprintf(stderr, "RC=%hhd, Bin=%hu, BinCount=%u, Node Count=%hu, GenPos %u, Chr %u:%u \n", x, i, binCounts[i], kGenPos->count, kGenPos->genPos, pChr->number, kGenPos->genPos - pChr->start);
    //                        k = kGenPos->nextPosInBin;
    //                    } while (k != 0);
    //                }
    //            }
    //        }
            if (bestCount > 0 && ties > 0) {
                //fprintf(stderr, "********** TIES %u, BestCount %u, Used Keys %u ***********", ties, bestCount, nUsedKeys);
                if (ties > 3 || (ties > 0 && bestCount * 3 / nUsedKeys < 1)) {
    //                fprintf(stderr, "CLEARED.\n%s\n",pMap->seq);
    //                for (i=0; i < BINS; i++) {
    //                    if (binCounts[i] >= bestCount) {
    //                        // Check
    //                        k = binRefPos[i];
    //                        do {
    //                            kGenPos = allGenPos + k;
    //                            if (kGenPos->count == bestCount) {
    //                                pChr = pMap->pGenome->getGenomicCoordinate(kGenPos->genPos);
    //                                fprintf(stderr, "RC=%hhd, Bin=%hu, Count=%hu, GenPos %u, Chr %u:%u \n", x, i, kGenPos->count, kGenPos->genPos, pChr->number, kGenPos->genPos - pChr->start);
    //                            }
    //                            k = kGenPos->nextPosInBin;
    //                        } while (k != 0);
    //                    }
    //                }
                    bestCount = 0;
                //} else {
                //    fprintf(stderr, "Ok!.\n%s\n",pMap->seq);
                }
            }
            if (bestCount > 0) {
                // If ties > 0 ???
                if (ties > 0) {
                    for (i = 0; i <= BINS; i++) {
                        if (binCounts[i] >= bestCount) {
                            // Check
                            k = binRefPos[i];
                            do {
                                kGenPos = allGenPos + k;
                                if (kGenPos->count == bestCount) {
                                    pMap->pCandPosMan->addGenomeAndReadPosition(kGenPos->genPos + pMap->getKeyablePosition(x, kGenPos->iRdPos), pMap->getKeyablePosition(x, kGenPos->iRdPos), x, bestCount * 100 / (nUsedKeys + ties));
                                    //pChr = pMap->pGenome->getGenomicCoordinate(kGenPos->genPos);
                                    //fprintf(stderr, "RC=%hhd, Bin=%hu, Count=%hu, GenPos %u, Chr %u:%u \n", x, i, kGenPos->count, kGenPos->genPos, pChr->number, kGenPos->genPos - pChr->start);
                                }
                                k = kGenPos->nextPosInBin;
                            } while (k != 0);
                        }
                    }
                } else {
                    pMap->pCandPosMan->addGenomeAndReadPosition(bestGenPos + pMap->getKeyablePosition(x, bestKeyPos), pMap->getKeyablePosition(x, bestKeyPos), x, bestCount * 100 / (nUsedKeys));
                }
            }
        }
        if (x == 0) {
            pMap->pCandPosMan->saveFwdCandidatePositions();
        }
    }

    //fprintf(stderr, "Counts: %u, Fwd:%u, Rev:%u\n", pMap->pCandPosMan->getCount(), pMap->pCandPosMan->getFwdCount(), pMap->pCandPosMan->getRevCount());
    return pMap->pCandPosMan->getCount();

}




// Half Apart, left keys are "paired" with right keys in order
uint16_t HalfApartMapping(ogReadKeyMapping *pMap, uint32_t maxTargets) {

    uint32_t nReadPos, halfPos;
    uint32_t i, iPos, intC;
    uint32_t j;
    uint8_t  intraFails;
    uint32_t keyDist = 3; // = 3;
    char    x, bug = pMap->debug;
    
    pMap->pCandPosMan->reset();
    
    for (x = 0; x < 2; x++) {
        if ((nReadPos = (x == 0 ? pMap->getNFwdPos() : pMap->getNRevPos())) > 1) {
            //fprintf(stderr, "nKeys=%d\n",nReadPos); fflush(stderr);
            halfPos = (nReadPos >> 1);

            intC = 0;

            for (i=0; i < halfPos; i++) {
                if (pMap->getTargetSize(x, i) < maxTargets) {
                    iPos = pMap->getKeyablePosition(x, i);
                    intC = 0;
                    intraFails = 0;
                    for (j=halfPos+i; j < nReadPos; j++) {
                        if (pMap->getTargetSize(x, j) < maxTargets) {
                            if (intC == 0) {
                                intC = pMap->intersectByTwoReadPositionKeys(x, i, j, pMap->getKeyablePosition(x,j) - iPos + keyDist);
                                if (intC > 10 || intC == 0) {
                                    intC = 0;
                                    if (++intraFails > 3) break;                                
                                }
                            } else {
                                intC = pMap->intersectByAddingReadPositionKey(x, j, pMap->getKeyablePosition(x,j) - iPos + keyDist);
                                if (intC == 0) {
                                    intC = pMap->pCandPosMan->rewindCountsAfterAddingKeyCountZero();
                                    if (++intraFails > 3) break;
                                }
                            }
                            if (intC > 0 && intC < 3) {
                                break;
                            }
                        }
                    }
                    if (intC > 0 && intC < 3) {
                        break;
                    }
                }
            }
        }
        if (x == 0) {
            if (intC == 0) { // (pMap->pCandPosMan->getCount() > 2)
                pMap->pCandPosMan->reset();
            }
            pMap->pCandPosMan->saveFwdCandidatePositions();
        } else {
            if (intC == 0) { // (pMap->pCandPosMan->getRevCount() > 2)
                pMap->pCandPosMan->restoreFwdCandidatePositions();
            }
            //if (intC > 2) {
            //    pMap->pCandPosMan->removeReverseTargets();
            //}
        }
    }

    return pMap->pCandPosMan->getCount();
}

uint16_t HalfMatchingSmall(ogReadKeyMapping *pMap) {
    return HalfApartMapping(pMap,  pMap->getLowKeyCountLimit());
}

uint16_t HalfMatchingLarge(ogReadKeyMapping *pMap) {
    return HalfApartMapping(pMap,  pMap->getHighKeyCountLimit());
}

uint16_t HalfMatchingHuge(ogReadKeyMapping *pMap) {
    return HalfApartMapping(pMap, 10*pMap->getHighKeyCountLimit());
}



// Generic Algorithm for a read segment and then finding the position by histogram. 
uint16_t SearchMatching(ogReadKeyMapping *pMap, uint32_t maxTargets) {
    uint16_t BINS = 1000;
    float fBINS = BINS;
    uint32_t minGenomePos = 0;
    //uint32_t minGenomePosCorrected = 0;
    uint32_t maxGenomePos = pMap->pGenome->getGenomeSize();
    uint32_t deltaGenomePos;
    float    dGenome;
    float    *binCounts; //binCounts[BINS+1];
    float    *binSums; // binSums[BINS+1];
    float    *binSqSums; //binSqSums[BINS+1];
    uint32_t maxPositions = 10000;
    uint32_t **keyInitialPos; //*keyInitialPos[maxPositions]; // pointers to keys , must be higher than readLength o keep it static per thread
    uint32_t **keyFinalPos; //*keyFinalPos[maxPositions];
    float    sizePos;
    uint32_t *p, *pp;
    uint32_t maxForwardPos = 0;
    uint32_t maxForwardCounts = 0;
    uint32_t maxReversePos = 0;
    uint32_t maxReverseCounts = 0;
    uint16_t startReadPos = 0, endReadPos, nActive, nFwdUsedKeys, nRevUsedKeys;
    //uint32_t lastKeyPos = pMap->getLastKeyablePosition();
    char x;
    float bin;
    //float red;
    uint32_t bestPos, bestGenPos;
    float    maxCount;
    uint32_t i;
    int32_t ibin, pibin, posCorrection;
    float   factor = 1;
    uint32_t last_i = 0, best_i;
    uint32_t last_i_GenPos = 0;
    uint32_t last_i_InPosFwd = 0;
    uint32_t last_i_InPosRev = 0;
    uint32_t currPos;
    ogKey   *pKey;
    ogChromosome *pChr;
    char    ties = 0;
    char    bug = pMap->debug;

    // Get Memory
    uint64_t bytes = (BINS+1)*(sizeof(float)*3) + maxPositions * sizeof(uint32_t);
    void *pV = pMap->getUsableMemory(bytes);
    // Redirect pointers
    binCounts = (float *) pV;
    binSums = (float *) (binCounts + (BINS+1));
    binSqSums = (float *) (binSums + (BINS+1));
    keyInitialPos = (uint32_t **) (binSqSums + (BINS+1));
    keyFinalPos = (uint32_t **) (keyInitialPos + maxPositions);
    //
    
//    if (endReadPos > lastKeyPos) endReadPos = lastKeyPos; 
//    //fprintf(stderr, "Rd=%s\n",pMap->seq);
//    if ((endReadPos - startReadPos)*4 > pMap->seqLen) {
//        pMap->setAllReadKeys();
//    } else {
//        for (i=startReadPos; i <= endReadPos; i++) pMap->setPosReadKey(i);
//    }
    for (x = 0; x < 2; x++) {
        endReadPos = (x == 0 ? pMap->getNFwdPos() : pMap->getNRevPos());
        factor = 0.1;
        memset(binCounts, 0, (BINS+1)*sizeof(float));
        
        minGenomePos = 0;
        maxGenomePos = pMap->pGenome->getGenomeSize();
        
        // Initialize
        nActive = 0;
        //while (true) {
            for (i=startReadPos; i < endReadPos; i++) {
                //if (x == 0) pMap->setPosReadKey(i);
                pKey = pMap->getInfoForPositionKey(x, i);
                // SAVE Starting and End Possible Genome positions
                keyInitialPos[i] = pMap->getGenomePositionFromReadPositionKey(x, i, 0);
                if (pKey->size < maxTargets) {
                    keyFinalPos[i] = keyInitialPos[i] + pKey->size - 1;
                    nActive++;
                } else {
                    keyFinalPos[i] = keyInitialPos[i] - 1;
                }
                //if (bug) fprintf(stderr, "RC:%hhd, iKey:%u, Size:%u, Pos:%u, Passed:%c\n", x, i, pKey->size, pMap->getKeyablePosition(x, i), pKey->size < maxTargets ? '1' : '0');
                //sizePos[i] = (pKey->size == 0 ? 1e-30 : (float) 1 / (float) pKey->size);
            }
        //    if (nActive < (endReadPos >> 1)) {
        //        nActive = 0;
        //        maxTargets <<= 1;
        //    } else {
        //        break;
        //    }
        //}

        while (true) {
            factor *= 10;
            deltaGenomePos = (maxGenomePos - minGenomePos);
            dGenome = deltaGenomePos;
            bin = dGenome / fBINS;
            if (bin < 1) bin = 1;
            //if (bug) fprintf(stderr, "RC:%hhd, MinG:%u, MaxG:%u, bin:%.1f\n", x, minGenomePos, maxGenomePos,bin);
            //red = ((float) (minGenomePos) / bin);
            // Set counts to 0
            //memset(binCounts, 0, sizeof(binCounts));
            bestPos = 0;
            bestGenPos = 0;
            maxCount = 0;
            best_i = 0 ;
            ties = 0;
            //minGenomePosCorrected = (minGenomePos > lastKeyPos ? minGenomePos - lastKeyPos : minGenomePos);
            
            //fprintf(stderr, "\n%u~%u/%.2f-",minGenomePos,maxGenomePos,bin);
            //posCorrection = (x == 0 ? -1 : endReadPos - startReadPos + 1);
            for (i=startReadPos; i < endReadPos; i++) {
                posCorrection = pMap->getKeyablePosition(x, i);
                //if (x) posCorrection--; else posCorrection++;
                // Advance keyInitialPos
                while (keyInitialPos[i] <= keyFinalPos[i] && *keyInitialPos[i]-posCorrection < minGenomePos) keyInitialPos[i]++;
                // Reduce keyFinalPos
                while (keyInitialPos[i] <= keyFinalPos[i] && *keyFinalPos[i]-posCorrection > maxGenomePos) keyFinalPos[i]--;
                sizePos = factor;// * 1 / (float) (keyFinalPos[i] - keyInitialPos[i] + 1);
                pibin = -1; // previous ibin
                for (p=keyInitialPos[i], pp=keyFinalPos[i]; p <= pp; p++) {
                    //last_i = i; // last valid position counting in the correct range
                    //last_i_GenPos = *p;
                    //ibin = (float) (*p - i) / bin - red; // real op: (*p - i - minGenomePos) * fBINS / deltaGenome
                    currPos = *p;
                    ibin = ((float) (currPos - posCorrection - minGenomePos)) / bin;
                    //if (bin < 1000) fprintf(stderr, "%hhd, %u, %u, ibin=%d\n",x,*p,*p-posCorrection,ibin);
                    if (ibin != pibin) {
                        if (ibin < 0) {
                            //fprintf(stderr, "ibin=%d,",ibin);
                            ibin = 0;
                        }
                        if (ibin > BINS) {
                            //fprintf(stderr, "ibin=%d,",ibin);
                            ibin = BINS;
                            //fprintf(stderr, "Problemas!!!\n");
                        }
                        //if (bug) {
                        //    pChr = pMap->pGenome->getGenomicCoordinate(currPos);
                        //    fprintf(stderr, "(%u:%d:Chr%u:%u)",i,ibin,pChr->number,currPos-pChr->start);
                        //}
                        if (binCounts[ibin] < sizePos) {
                            binCounts[ibin] = 0; // reset
                            binSums[ibin] = 0;
                            binSqSums[ibin] = 0;
                        }
                        binCounts[ibin] += sizePos;
                        currPos -= (ibin * bin + minGenomePos); // relative to expected
                        binSums[ibin] += currPos;
                        binSqSums[ibin] += (float) currPos * (float) currPos;
                        if (binCounts[ibin] >= maxCount) {
                            if (binCounts[ibin] == maxCount) {
                                ties = 1;
                            } else {
                                bestPos = ibin;
                                bestGenPos = *p - posCorrection;
                                maxCount = binCounts[bestPos];
                                last_i_GenPos = *p;
                                best_i = i;
                                ties = 0;
                            }
                        }
                        pibin = ibin;
                    }
                }
            }
            if (ties) {
                int correctCounts = binCounts[bestPos] / sizePos; // count is scaled to factor [sizePos]
                float bestMean = binSums[bestPos]  / correctCounts;
                float bestDev = binSqSums[bestPos] / correctCounts - bestMean*bestMean;
                float iMean;
                float iDev;
                for (i=0; i <= BINS; i++) {
                    if (binCounts[i] == maxCount && i != bestPos) {
                        iMean = binSums[i]  / correctCounts;
                        iDev = binSqSums[i] / correctCounts - iMean*iMean;
                        if (iDev < bestDev) {
                            bestPos = i;
                            bestMean = iMean;
                            bestDev = iDev;
                        }
                    }
                }
            }
            //if (bug) fprintf(stderr, "|bin=%.4f|best=%u|c=%.4f\n",bin,bestPos,maxCount);
            if (bestPos > 0) minGenomePos += (uint32_t) (bin * bestPos - bin / 2);
            maxGenomePos  = minGenomePos + (uint32_t) (bin * 1.5);
            if (maxCount < .000001 || bin <= 1 || bin < pMap->read->lenSeq) break; // || bin < BINS
            //if (minGenomePos > lastKeyPos) minGenomePos -= lastKeyPos; else minGenomePos = 0;
            //maxGenomePos  = minGenomePos + bin * (bestPos+1) + (bestPos < BINS ? (bin * (float) binCounts[bestPos+1] / (float) maxCount) : 0);
            //minGenomePos += (uint32_t) (bin * bestPos - (bestPos > 0 ? (bin * (float) binCounts[bestPos-1] / (float) maxCount) : 0));
        }
        minGenomePos = bestGenPos;
        //if (bug) {
        //    fprintf(stderr, "[RC=%hhu, GP:%u, C:%d, bin=%.3f]\n", x, minGenomePos, (int) maxCount, bin);
        //    pChr = pMap->pGenome->getGenomicCoordinate(minGenomePos);
        //    fprintf(stderr, "Chr %u (%s) Pos %u.\n", pChr->number, pChr->name, minGenomePos - pChr->start);
        //}
        if (x) {
            last_i_InPosRev = best_i; //last_i;
            maxReversePos = last_i_GenPos; //minGenomePos;
            maxReverseCounts = maxCount;
            nRevUsedKeys = nActive;
        } else {
            last_i_InPosFwd = best_i; // last_i;
            maxForwardPos = last_i_GenPos; //minGenomePos;
            maxForwardCounts = maxCount;            
            nFwdUsedKeys = nActive;
        }
        
    }

    pMap->pCandPosMan->reset();
    if (maxReverseCounts > maxForwardCounts  && maxReverseCounts > 0) {
        pMap->pCandPosMan->saveFwdCandidatePositions();
        pMap->pCandPosMan->addGenomeAndReadPosition(maxReversePos, pMap->getRevKeyablePosition(last_i_InPosRev), 1, maxReverseCounts * 100 / (nRevUsedKeys)); //pMap->pCandPosMan->setByFixedPosition(1, maxReversePos);
    } else if (maxForwardCounts > 0) {
        pMap->pCandPosMan->addGenomeAndReadPosition(maxForwardPos, pMap->getFwdKeyablePosition(last_i_InPosFwd), 0, maxForwardCounts * 100 / (nFwdUsedKeys)); //pMap->pCandPosMan->setByFixedPosition(0, maxForwardPos);
        pMap->pCandPosMan->saveFwdCandidatePositions();
    } else {
        return 0;
    }
    return 1;        
}


uint16_t DefaultMatchingSmall(ogReadKeyMapping *pMap) {
    return SearchMatching(pMap,  pMap->getLowKeyCountLimit());
}

uint16_t DefaultMatchingLarge(ogReadKeyMapping *pMap) {
    return SearchMatching(pMap,  pMap->getHighKeyCountLimit());
}

uint16_t DefaultMatchingHuge(ogReadKeyMapping *pMap) {
    return SearchMatching(pMap,  10*pMap->getHighKeyCountLimit());
}



uint16_t HistogramMatchingSmall(ogReadKeyMapping *pMap) {
    //return HistogramMapping(pMap,  pMap->getLowKeyCountLimit());
    return SortedSizesHistogramMapping(pMap,  pMap->getLowKeyCountLimit());    
}

uint16_t HistogramMatchingLarge(ogReadKeyMapping *pMap) {
    //return HistogramMapping(pMap,  pMap->getHighKeyCountLimit());
    return SortedSizesHistogramMapping(pMap,  pMap->getHighKeyCountLimit());    
}

uint16_t HistogramMatchingHuge(ogReadKeyMapping *pMap) {
    //return HistogramMapping(pMap,  10*pMap->getHighKeyCountLimit());
    return SortedSizesHistogramMapping(pMap,  10*pMap->getHighKeyCountLimit());    
}

uint16_t HistogramMatchingAll(ogReadKeyMapping *pMap) {
    //return HistogramMapping(pMap,  10*pMap->getHighKeyCountLimit());
    return SortedSizesHistogramMapping(pMap, 0xFFFFFFFF);    
}


uint16_t BinnedHistogramMatchingSmall(ogReadKeyMapping *pMap) {
    //return HistogramMapping(pMap,  pMap->getLowKeyCountLimit());
    return BinnedHistogramMapping(pMap,  pMap->getLowKeyCountLimit());    
}

uint16_t BinnedHistogramMatchingLarge(ogReadKeyMapping *pMap) {
    //return HistogramMapping(pMap,  pMap->getHighKeyCountLimit());
    return BinnedHistogramMapping(pMap,  pMap->getHighKeyCountLimit());    
}

uint16_t BinnedHistogramMatchingHuge(ogReadKeyMapping *pMap) {
    //return HistogramMapping(pMap,  10*pMap->getHighKeyCountLimit());
    return BinnedHistogramMapping(pMap,  10*pMap->getHighKeyCountLimit());    
}

uint16_t BinnedHistogramMatchingAll(ogReadKeyMapping *pMap) {
    //return HistogramMapping(pMap,  10*pMap->getHighKeyCountLimit());
    return BinnedHistogramMapping(pMap,  0xFFFFFFFF);    
}

// Algorithm X - All Keys are checked, target position is "search" by cycles of reducing ranges
uint16_t MapPosMatching(ogReadKeyMapping *pMap) {
    std::map<uint32_t, uint32_t> mapFwd;
    std::map<uint32_t, uint32_t> mapRev;
    std::map<uint32_t, uint32_t> mapKeyFwd;
    std::map<uint32_t, uint32_t> mapKeyRev;
    uint32_t lastKeyPos = pMap->getLastKeyablePosition();
    uint32_t i, *p, *pp, *q;
    uint32_t pos, max, maxFwdPos, maxFwd, maxRevPos, maxRev, maxEqFwd, maxEqRev, maxEq;
    ogKey *pKey;
    char RC = 0;

    pMap->setAllReadKeys();
    maxFwd = 0;
    maxFwdPos = -1;
    maxRev = 0;
    maxEqFwd = 0;
    maxEqRev = 0;
    maxRevPos = -1;
    /**
    fprintf(stderr, "Size/Key Fwd:");
    for (i=0; i <= lastKeyPos; i++) {
        pKey = pMap->getInfoForPositionKey(0, i);
        fprintf(stderr, "%u,", pKey->size);
    }
    fprintf(stderr, "\nSize/Key Rev:");
    for (i=0; i <= lastKeyPos; i++) {
        pKey = pMap->getInfoForPositionKey(1, i);
        fprintf(stderr, "%u,", pKey->size);
    }
    fprintf(stderr, "\n");
     **/
    uint32_t prevFwdKey = -1;
    uint32_t prevRevKey = -1;
    for (i=0; i <= lastKeyPos; i++) {
        pKey = pMap->getInfoForPositionKey(0, i);
        if (pKey->size < 100000 && prevFwdKey != pMap->getKeyForward(i)) {
            prevFwdKey = pMap->getKeyForward(i);
            if (mapKeyFwd[prevFwdKey] == 0) {
                mapKeyFwd[prevFwdKey]++;
                p = pMap->getGenomePositionFromReadPositionKey(0, i, 0);
                pp = p + pKey->size;
                for (q=p; q < pp; q++) {
                    pos = *q - i;
                    if (++mapFwd[pos] >= maxFwd) {
                        if (mapFwd[pos] == maxFwd) {
                            maxEqFwd++;
                        } else {
                            maxFwd = mapFwd[pos];
                            maxFwdPos = pos;
                            maxEqFwd = 0;
                        }
                    }
                }
            }
        }
        pKey = pMap->getInfoForPositionKey(1, i);
        if (pKey->size < 100000 && prevRevKey != pMap->getKeyReverse(i)) {
            prevRevKey = pMap->getKeyReverse(i);
            if (mapKeyRev[prevRevKey] == 0) {
                mapKeyRev[prevRevKey]++;
                p = pMap->getGenomePositionFromReadPositionKey(1, i, 0);
                pp = p + pKey->size;
                for (q=p; q < pp; q++) {
                    pos = *q + i;
                    if (++mapRev[pos] >= maxRev) {
                        if (mapRev[pos] == maxRev) {
                            maxEqRev++;
                        } else {
                            maxRev = mapRev[pos];
                            maxRevPos = pos;
                            maxEqRev = 0;
                        }
                    }
                }
            }
        }
    }

    pMap->pCandPosMan->reset();
    if (maxRev > maxFwd) {
        pMap->pCandPosMan->addGenomeAndReadPosition(maxRevPos, 0, 1, 0);//pMap->pCandPosMan->setByFixedPosition(1, maxRevPos);
        pos = maxRevPos;
        max = maxRev;
        maxEq = maxEqRev;
        RC = 1;
    } else {
        pMap->pCandPosMan->addGenomeAndReadPosition(maxFwdPos, 0, 0, 0); //pMap->pCandPosMan->setByFixedPosition(0, maxFwdPos);
        pos = maxFwdPos;
        max = maxFwd;
        maxEq = maxEqFwd;
        RC = 0;
    }

    /**
    i = 0;
    fprintf(stderr, "%u:%hhu:%u:%s\n", max, RC, maxEq, pMap->seq);
    fprintf(stderr, "Pos %d, Gen Pos %u, ",i+1, pMap->pCandPosMan->getkPos(i));
    ogChromosome *pChr = pMap->pGenome->getGenomicCoordinate(pMap->pCandPosMan->getkPos(i));
    fprintf(stderr, "Chr %s Pos %u.\n", pChr->name, pMap->pCandPosMan->getkPos(i) - pChr->start);
    
    char DNA[101];
    DNA[100] = 0;
    pMap->pGenome->extractPackedGenome(pMap->pCandPosMan->getkPos(i), 100, DNA, 0);
    fprintf(stderr, "GP:%s\n", DNA);
    
     **/
    
    return 1;
}



// Extend the 
uint16_t MovingWindowIntersect(ogReadKeyMapping *pMap, uint16_t start, uint16_t steps, uint16_t window, uint16_t maxIntersection) {
    uint32_t i, j, k;
    ogKey *pKey;
    char RC = 0;
    uint32_t maxExtFR[2] = {0, 0}, maxExtPosFR[2] = {0,0}, maxExtFREq[2] = {0,0}, maxGenPosFR[2] = {0,0}, maxCountsFR[2]={0,0};
    uint32_t maxExt = 0, maxPos = 0, maxEq = 0, maxGenPos, maxCounts;
    uint32_t intC;
    pMap->setAllReadKeys();
    //fprintf(stderr, "Rd=%s\n",pMap->seq);
    char flagTwoPos = 0;
    uint16_t searching = 0;
    for (RC = 0; RC < 2; RC++) {
        uint32_t lastKeyPos = (RC == 0 ? pMap->getNFwdPos() : pMap->getNRevPos());
        for (i=start; i < lastKeyPos - maxExtFR[RC]; i += steps) {
            flagTwoPos = 0;
            searching = window;
            maxExt = 0;
            maxCounts = 0;
            maxGenPos = 0;
            j = i;
            while (j < lastKeyPos) {
                k = j + searching;
                if (k > lastKeyPos) {
                    k = lastKeyPos;
                    searching = k - j;
                }
                if (flagTwoPos) {
                    intC = pMap->intersectByAddingReadPositionKey(RC, k, k - i + 100);
                } else {
                    intC = pMap->intersectByTwoReadPositionKeys(RC, i, k, k - i + 1);
                    if (intC > maxIntersection) {
                        break;
                    }
                }
                if (intC > 0) {
                    // intersection works, so
                    j = k;
                    flagTwoPos = 1;
                    maxExt = k - i;
                    maxCounts = pMap->pCandPosMan->getCount();
                    maxGenPos = pMap->pCandPosMan->getkPos(0)->genomePosition;
                    if (searching < window && searching > 1) searching >>= 1;
                } else {
                    if (searching > 1) {
                        // No ha terminado de encontrar el punto.
                        searching >>= 1;
                        pMap->pCandPosMan->rewindCountsAfterAddingKeyCountZero();
                    } else {
                        // Terminó
                        break;
                    }
                }
            }
            if (maxExt >= maxExtFR[RC]) {
                if (maxExt == maxExtFR[RC]) {
                    maxExtFREq[RC]++;
                } else {
                    maxCountsFR[RC] = maxCounts;
                    maxGenPosFR[RC] = maxGenPos;
                    maxExtFR[RC] = maxExt;
                    maxExtPosFR[RC] = i;
                    maxExtFREq[RC] = 0;
                }
            }
        }
    }
    
    RC = (maxExtFR[0] >= maxExtFR[1] ? 0 : 1);
    pMap->pCandPosMan->reset();
    pMap->pCandPosMan->addGenomeAndReadPosition(maxGenPosFR[RC],0,RC,0); //pMap->pCandPosMan->setByFixedPosition(RC, maxGenPosFR[RC]);

    /*
    char DNA[1001];
    DNA[1000] = 0;
    i = 0;
    fprintf(stderr, "Pos=%u, Ext=%u, Counts=%u, RC=%hhu, Eq=%u, Rd=%s\n", maxPos, maxExt, maxCounts, RC, maxEq, pMap->seq);
    pMap->pGenome->extractPackedGenome(pMap->pCandPosMan->getkPos(i), pMap->keySize+maxExt, DNA, 0);
    DNA[pMap->keySize+maxExt] = 0;
    fprintf(stderr, "Genome Seq (+) =%s\n", DNA);
    pMap->pGenome->extractPackedGenome(pMap->pCandPosMan->getkPos(i), pMap->keySize+maxExt, DNA, 1);
    fprintf(stderr, "Genome Seq (-) =%s\n", DNA);
    fprintf(stderr, "Read   Seq Pos =%s\n", pMap->seq + maxPos);
    fprintf(stderr, "Pos %d, Gen Pos %u, ",i+1, pMap->pCandPosMan->getkPos(i));
    ogChromosome *pChr = pMap->pGenome->getGenomicCoordinate(pMap->pCandPosMan->getkPos(i));
    fprintf(stderr, "Chr %s Pos %u.\n", pChr->name, pMap->pCandPosMan->getkPos(i) - pChr->start);
    */
    
    return maxCountsFR[RC];
}


// Algorithm X - All Keys are checked until count = 1, target position is "search" by cycles of reducing ranges
uint16_t NestedIntersectGeneric(ogReadKeyMapping *pMap, uint32_t maxTargets) {
    //return MovingWindowIntersect(pMap, 0, 5, 20, 1000);
    //return MovingWindowIntersect(pMap, 0, 1, 3, 1000);
    uint32_t i, i2, j;
    char RC = 0;
    const uint32_t NOTFOUND = 0x7FFFFFFF;
    uint32_t maxExtPosFwd, maxCountsFwd = NOTFOUND;
    uint32_t maxExtPosRev, maxCountsRev = NOTFOUND;
    uint32_t maxCounts, zeros;
    uint32_t intC;
    uint32_t lastKeyPos;
    //fprintf(stderr, "Rd=%s\n",pMap->seq);
    pMap->pCandPosMan->reset();
    for (RC = 0; RC < 2; RC++) {
        lastKeyPos = (RC ? pMap->getNRevPos() : pMap->getNFwdPos());
        for (i=0; i < lastKeyPos-1; i++) {
            if (pMap->getTargetSize(RC, i) < maxTargets) {
                i2 = (i < lastKeyPos-2 ? i+2 : i + 1);
                intC = pMap->intersectByTwoReadPositionKeys(RC, i, i2, pMap->getKeyablePosition(RC, i2)-pMap->getKeyablePosition(RC, i)+3);
                for (j=i2+1, zeros=0; intC > 1 && j < lastKeyPos && zeros < 2; j++) {
                    if (pMap->getTargetSize(RC, j) < maxTargets) {
                        maxCounts = intC;
                        intC = pMap->intersectByAddingReadPositionKey(RC, j, pMap->getKeyablePosition(RC, j)-pMap->getKeyablePosition(RC, i)+3);
                        if (intC == 0) {
                            pMap->pCandPosMan->rewindCountsAfterAddingKeyCountZero();
                            intC = maxCounts;
                            zeros++; // to optimize, if two additional keys has no intersection is likely that keys i & i+1 are wrong
                        }
                    }
                }
                if (RC == 0) {
                    if (intC > 0 && intC < maxCountsFwd) {
                        maxCountsFwd = intC;
                        maxExtPosFwd = i; // genomePosition is relative to first key
                    }
                } else {
                    if (intC > 0 && intC < maxCountsRev) {
                        maxCountsRev = intC;
                        maxExtPosRev = i; // genomePosition is relative to first key
                    }
                }
                if (intC == 1) {
                    break;
                }
            }
        }
    }

    if (maxCountsRev < NOTFOUND && maxCountsRev > 0 && maxCountsRev < maxCountsFwd) {
        RC = 1;
        i = maxExtPosRev;
    } else if (maxCountsFwd < NOTFOUND && maxCountsFwd > 0) {
        RC = 0;
        i = maxExtPosFwd;
    } else {
        return 0;
    }

    // Reproduce best config
    lastKeyPos = (RC ? pMap->getNRevPos() : pMap->getNFwdPos());
    pMap->pCandPosMan->reset();
    if (RC == 1) pMap->pCandPosMan->saveFwdCandidatePositions();
    i2 = (i < lastKeyPos-2 ? i+2 : i + 1);
    intC = pMap->intersectByTwoReadPositionKeys(RC, i, i2, pMap->getKeyablePosition(RC, i2)-pMap->getKeyablePosition(RC, i)+3);
    for (j=i2+1, zeros=0; intC > 1 && j < lastKeyPos && zeros < 2; j++) {
        if (pMap->getTargetSize(RC, j) < maxTargets) {
            maxCounts = intC;
            intC = pMap->intersectByAddingReadPositionKey(RC, j, pMap->getKeyablePosition(RC, j)-pMap->getKeyablePosition(RC, i)+3);
            if (intC == 0) {
                pMap->pCandPosMan->rewindCountsAfterAddingKeyCountZero();
                intC = maxCounts;
                zeros++; // to optimize, if two additional keys has no intersection is likely that keys i & i+1 are wrong
            }
        }
    }
    if (RC == 0) pMap->pCandPosMan->saveFwdCandidatePositions();
    if (maxCountsRev == maxCountsFwd) {
        lastKeyPos = (RC ? pMap->getNFwdPos() : pMap->getNRevPos());
        i = maxExtPosRev;
        RC = 1-RC;
        i2 = (i < lastKeyPos-2 ? i+2 : i + 1);
        intC = pMap->intersectByTwoReadPositionKeys(RC, i, i2, pMap->getKeyablePosition(RC, i2)-pMap->getKeyablePosition(RC, i)+3);
        for (j=i2+1, zeros=0; intC > 1 && j < lastKeyPos && zeros < 2; j++) {
            if (pMap->getTargetSize(RC, j) < maxTargets) {
                maxCounts = intC;
                intC = pMap->intersectByAddingReadPositionKey(RC, j, pMap->getKeyablePosition(RC, j)-pMap->getKeyablePosition(RC, i)+3);
                if (intC == 0) {
                    pMap->pCandPosMan->rewindCountsAfterAddingKeyCountZero();
                    intC = maxCounts;
                    zeros++; // to optimize, if two additional keys has no intersection is likely that keys i & i+1 are wrong
                }
            }
        }
        intC += maxCountsFwd; // Fwd because intC ahorita trae Rev y entonces hay que sumar el Fwd
    }
    return intC;
}


uint16_t NestedIntersectSmall(ogReadKeyMapping *pMap) {
    return NestedIntersectGeneric(pMap, pMap->getLowKeyCountLimit());
}

uint16_t NestedIntersectLarge(ogReadKeyMapping *pMap) {
    return NestedIntersectGeneric(pMap, pMap->getHighKeyCountLimit());
}


uint16_t LeftMost2KeysLimitingTargets(ogReadKeyMapping *pMap, uint16_t maxTargets) {
    
    pMap->pCandPosMan->reset();
    
    uint32_t k1 = 0, k2 = 1, k3;
    uint32_t nKF = pMap->getNFwdPos();
    uint32_t nKR = pMap->getNRevPos();
    uint32_t p1,p2,p3;
    uint32_t d ;
    uint16_t fwdCount=0, revCount=0;
    uint16_t i;
    
    for (k1=0,k2=0; fwdCount == 0 && k1 < nKF-1 && k2 < nKF; k1++) {
        if (pMap->getFwdTargetSize(k1) <= maxTargets) {
            p1 = pMap->getFwdKeyablePosition(k1);
            for (k2=k1+1; k2 < nKF && ( (p2=pMap->getFwdKeyablePosition(k2)) - p1 < 10  || pMap->getFwdTargetSize(k2) > maxTargets); k2++);
            if (k2 < nKF) {
                //p2 = pMap->getFwdKeyablePosition(k2);
                // This assumes that LEFT KEYS GENERATING HITS ARE TRUE and remove left key
                fwdCount = pMap->intersectByTwoReadPositionKeys(0, k1, k2, p2-p1+3);
            }
        }
    }
    for(k3=k2+1; fwdCount > 5 && k3 < nKF; k3++) {
        if (pMap->getFwdTargetSize(k3) <= maxTargets) {
            p3 = pMap->getFwdKeyablePosition(k3);
            d = p3 - p1 + 2;
            fwdCount = pMap->intersectByAddingReadPositionKey(0,k3,d);
            if (fwdCount == 0) fwdCount = pMap->pCandPosMan->rewindCountsAfterAddingKeyCountZero();
        }
    }
    
    pMap->pCandPosMan->saveFwdCandidatePositions();
    
    for (k1=0, k2=0; revCount == 0 && k1 < nKR-1 && k2 < nKR; k1++) {
        if (pMap->getRevTargetSize(k1) <= maxTargets) {
            p1 = pMap->getRevKeyablePosition(k1);
            for (k2=k1+1; k2 < nKR && ( (p2=pMap->getRevKeyablePosition(k2)) - p1 < 10 || pMap->getRevTargetSize(k2) > maxTargets); k2++);
            if (k2 < nKR) {
                //p2 = pMap->getRevKeyablePosition(k2);
                // This assumes that LEFT KEYS GENERATING HITS ARE TRUE and remove left key
                revCount = pMap->intersectByTwoReadPositionKeys(1, k1, k2, p2-p1+3);
            }
        }
    }
    for(k3=k2+1; revCount > 5 && k3 < nKR; k3++) {
        if (pMap->getRevTargetSize(k3) <= maxTargets) {
            p3 = pMap->getRevKeyablePosition(k3);
            d = p3 - p1 + 2;
            revCount = pMap->intersectByAddingReadPositionKey(1,k3,d);
            if (revCount == 0) revCount = pMap->pCandPosMan->rewindCountsAfterAddingKeyCountZero();
        }
    }
    //pMap->pCandPosMan->saveFwdCandidatePositions();//pMap->pCandPosMan->setIsReverseComplement(revCount > fwdCount || (revCount == fwdCount && nKR > nKF)? 1 : 0);
    //if (fwdCount + revCount > 100) {
        //fprintf(stderr, "%llu, %u targets, %u fwd, %u rev:\n%s\n", iProc, fwdCount+revCount, fwdCount, revCount, pMap->seq);
    //}
    return fwdCount + revCount;
}


uint16_t LeftMost2KeysQuick(ogReadKeyMapping *pMap) {

    /**
Falta poner aqui el median size porque en HPCEncoding 100 es muy poco.
Distribution of positions per key:
// HPCEncoding
<=0:16069055	<=1:1025	<=2:1636	<=3:2162	<=4:2656
<=5:3396	<=6:3885	<=7:4281	<=8:4864	<=9:5213
<=10:5249	<=20:45592	<=30:25983	<=40:18667	<=50:18603
<=60:21718	<=70:24176	<=80:25186	<=90:24607	<=100:23674
<=200:129959	<=300:33813	<=400:12943	<=500:10621	<=600:13954
<=700:17918	<=800:21013	<=900:22557	<=1000:22457	<=2000:120889
<=3000:23863	<=4000:6572	<=5000:2848	<=6000:1704	<=7000:1090
<=8000:673	<=9000:408	<=10000:324	<=20000:1258	<=30000:285
<=40000:121	<=50000:83	<=60000:54	<=70000:42	<=80000:36
<=90000:18	<=100000:15	<=200000:51	<=300000:8	<=400000:3
<=500000:1	<=600000:1	<=700000:0	<=800000:1	<=900000:1
<=1000000:0	<=4294967295:4	
4/4 Unpacking Genome: Positioning keys.

// Bitwise Encoding
<=0:160803	<=1:387041	<=2:577699	<=3:705631	<=4:781036
<=5:812730	<=6:815588	<=7:799106	<=8:768318	<=9:731294
<=10:689511	<=20:4621462	<=30:2007353	<=40:955634	<=50:535070
<=60:337711	<=70:212208	<=80:145774	<=90:115715	<=100:95766
<=200:331787	<=300:88040	<=400:29333	<=500:26148	<=600:9405
<=700:4032	<=800:2783	<=900:4042	<=1000:4655	<=2000:11850
<=3000:3680	<=4000:1668	<=5000:1007	<=6000:834	<=7000:548
<=8000:420	<=9000:241	<=10000:189	<=20000:727	<=30000:150
<=40000:57	<=50000:58	<=60000:40	<=70000:29	<=80000:12
<=90000:6	<=100000:1	<=200000:16	<=300000:2	<=400000:0
<=500000:0	<=600000:0	<=700000:0	<=800000:0	<=900000:0
<=1000000:0	<=4294967295:6	

     */
    return LeftMost2KeysLimitingTargets(pMap, pMap->getLowKeyCountLimit());
}


uint16_t LeftMost2KeysLong(ogReadKeyMapping *pMap) {
    return LeftMost2KeysLimitingTargets(pMap, pMap->getHighKeyCountLimit());
}


uint16_t LeftMost2KeysHuge(ogReadKeyMapping *pMap) {
    return LeftMost2KeysLimitingTargets(pMap, 10*pMap->getHighKeyCountLimit());
}


uint16_t Min2KeysLimitingTargets(ogReadKeyMapping *pMap, uint16_t maxTargets) {
    
    const uint32_t NOTFOUND = 0x7FFFFFFF;
    //uint16_t fp1, fp2, rp1, rp2;
    uint16_t fno = 65535, rno = 65535, k, i;
    uint16_t fwdCount = 0;
    uint16_t revCount = 0;
    uint32_t a, r1=NOTFOUND, f1=NOTFOUND;
    uint32_t r2=r1, f2=f1, f01flag = 0, r01flag = 0;
    uint32_t fi1=0, ri1=0, fi2=0, ri2=0, fp1, fp2, p, rp1, rp2;
    uint32_t nKF = pMap->getNFwdPos();
    uint32_t nKR = pMap->getNRevPos();
    uint16_t minKeyDist = pMap->keySize / 2;
    //uint16_t maxTargets = 200;
    
    pMap->pCandPosMan->reset();
    
    if (nKF >= 2 && pMap->minFwdKeySize <= maxTargets) {
        for (k=0; pMap->nActiveFwd > 1; k++) { // prev condition k < 10
            // Default: make assumptions
            f1 = NOTFOUND;
            fp1 = NOTFOUND;
            fi1 = 0;
            f2 = NOTFOUND;
            fp2 = NOTFOUND;
            fi2 = 1;
            if (pMap->nActiveFwd > 0) {
                for (i=0; i < nKF; i++) {
                    if (pMap->isFwdKeyActive(i)) { //actF[i]
                        a = pMap->getFwdTargetSize(i);
                        if (a > 0 && a < maxTargets) {
                            p = pMap->getFwdKeyablePosition(i);
                            if (a < f1) {
                                if (p - fp1 >= minKeyDist) {
                                    f2 = f1;
                                    fi2 = fi1;
                                    fp2 = fp1;
                                }
                                f1 = a;
                                fi1 = i;
                                fp1 = p;
                            } else if (a < f2 && p - fp1 >= minKeyDist) {
                                f2 = a;
                                fi2 = i;
                                fp2 = p;
                            }
                        }
                    }
                }
                if (fi2 < fi1) { 
                    a=fi1;
                    fi1=fi2;
                    fi2=a;
                    a=f1;
                    f1=f2;
                    f2=a;
                }
                if (fwdCount == 0) {
                    if (f2 == NOTFOUND || f1 == NOTFOUND) break;
                    fwdCount = pMap->intersectByTwoReadPositionKeys(0, fi1, fi2, pMap->getFwdKeyablePosition(fi2) - pMap->getFwdKeyablePosition(fi1) + 3);
                    if (pMap->getFwdTargetSize(fi1) < pMap->getFwdTargetSize(fi2)) {
                        pMap->inactivePositionFwd(fi1); // NOT AROUND: inactivePositionsAroundKeyRev
                    } else {
                        pMap->inactivePositionFwd(fi2); // because this search for minimum which may be associated to seq errors and divergent keys
                    }
                } else {
                    if (f1 != NOTFOUND) {
                        fwdCount = pMap->intersectByAddingReadPositionKey(0, fi1, pMap->read->lenSeq);
                        if (fwdCount == 0) {
                            fwdCount = pMap->pCandPosMan->rewindCountsAfterAddingKeyCountZero();
                        }
                        pMap->inactivePositionFwd(fi1); // NOT AROUND: inactivePositionsAroundKeyRev
                    } else if (f2 != NOTFOUND) {
                        fwdCount = pMap->intersectByAddingReadPositionKey(0, fi2, pMap->read->lenSeq);
                        if (fwdCount == 0) {
                            fwdCount = pMap->pCandPosMan->rewindCountsAfterAddingKeyCountZero();
                        }
                        pMap->inactivePositionFwd(fi2); // because this search for minimum which may be associated to seq errors and divergent keys
                    } else break;
                }
                if (fwdCount == 1) {
                    break;
                }
                if (pMap->nActiveFwd == 0) break;
            }
        }
    }
    
    pMap->pCandPosMan->saveFwdCandidatePositions();

    if (nKR >= 2 && pMap->minRevKeySize <= maxTargets) {
        for (k=0; pMap->nActiveRev > 1; k++) {  // prev condition k < 10
            // Default: make assumptions
            r1 = NOTFOUND;
            rp1 = NOTFOUND;
            ri1 = 0;
            r2 = NOTFOUND;
            rp2 = NOTFOUND;
            ri2 = 1;
            if (pMap->nActiveRev > 0) {
                for (i=0; i < nKR; i++) {
                    if (pMap->isRevKeyActive(i)) { //actF[i]
                        a = pMap->getRevTargetSize(i);
                        if (a > 0 && a < maxTargets) {
                            p = pMap->getRevKeyablePosition(i);
                            if (a < r1) {
                                if (p - rp1 >= minKeyDist) {
                                    r2 = r1;
                                    ri2 = ri1;
                                    rp2 = rp1;
                                }
                                r1 = a;
                                ri1 = i;
                                rp1 = p;
                            } else if (a < r2 && p - rp1 >= minKeyDist) {
                                r2 = a;
                                ri2 = i;
                                rp2 = p;
                            }
                        }
                    }
                }
                if (ri2 < ri1) { 
                    a=ri1;
                    ri1=ri2;
                    ri2=a;
                    a=r1;
                    r1=r2;
                    r2=a;
                }
                if (revCount == 0) {
                    if (r2 == NOTFOUND || r1 == NOTFOUND) break;
                    revCount = pMap->intersectByTwoReadPositionKeys(1, ri1, ri2, pMap->getRevKeyablePosition(ri2) - pMap->getRevKeyablePosition(ri1) + 3);
                    if (pMap->getRevTargetSize(ri1) < pMap->getFwdTargetSize(ri2)) {
                        pMap->inactivePositionRev(ri1); // NOT AROUND: inactivePositionsAroundKeyRev ... Around is faster but lose 2%, better to be faster
                    } else {
                        pMap->inactivePositionRev(ri2); // because this search for minimum which may be associated to seq errors and divergent keys
                    }
                } else {
                    if (r1 != NOTFOUND) {
                        revCount = pMap->intersectByAddingReadPositionKey(1, ri1, pMap->read->lenSeq);
                        if (revCount == 0) {
                            revCount = pMap->pCandPosMan->rewindCountsAfterAddingKeyCountZero();
                        }
                        pMap->inactivePositionRev(ri1); // NOT AROUND: inactivePositionsAroundKeyRev ... Around is faster but lose 2%, better to be faster
                    } else if (r2 != NOTFOUND) {
                        revCount = pMap->intersectByAddingReadPositionKey(1, ri2, pMap->read->lenSeq);
                        if (revCount == 0) {
                            revCount = pMap->pCandPosMan->rewindCountsAfterAddingKeyCountZero();
                        }
                        pMap->inactivePositionRev(ri2); // because this search for minimum which may be associated to seq errors and divergent keys
                    } else break;
                }
                if (revCount == 1) {
                    break;
                }
                if (pMap->nActiveRev == 0) break;
            }
        }
    }
    
//        if (revCount == 0 && fwdCount == 0) {
////            fprintf(stderr, "\n%s\n%u Keys, f1=%u, f2=%u, r1=%u, r2=%u\n", pMap->seq, nK, fi1, fi2, ri1, ri2);
////            for (i=0; i < nK; i++) {
////                fprintf(stderr, "%u : fPos=%u : rPos=%u : F=%u : R=%u\n", i, pMap->getFwdKeyablePosition(i),  pMap->getRevKeyablePosition(i), pMap->getFwdInfoForKeyablePosition(i)->size, pMap->getRevInfoForKeyablePosition(i)->size);
////            }
//        }
    
    if (revCount == 0 && fwdCount == 0) {
        //fprintf(stderr, "** PROBLEM Zero Counts**\n");
        //fprintf(stderr, "%s\n", pMap->seq);
    }
    if (fwdCount > 5 || revCount > 5) {
        //fprintf(stderr, "%d",k);
    }
    //pMap->pCandPosMan->setIsReverseComplement(revCount >= fwdCount ? 1 : 0);
    //return  pMap->pCandPosMan->getFwdRevCount();
    return pMap->pCandPosMan->getCount();
}



uint16_t Min2KeysSize(ogReadKeyMapping *pMap) {
    return Min2KeysLimitingTargets(pMap, pMap->getLowKeyCountLimit()); // ,200
}


uint16_t Min2KeysLong(ogReadKeyMapping *pMap) {
    return Min2KeysLimitingTargets(pMap, pMap->getHighKeyCountLimit()); // , 65535
}


uint16_t Min2KeysHuge(ogReadKeyMapping *pMap) {
    return Min2KeysLimitingTargets(pMap, 10*pMap->getHighKeyCountLimit()); // , 65535
}



uint16_t MinApart2Keys(ogReadKeyMapping *pMap) {
    
    const uint32_t NOTFOUND = 100000000;
    uint16_t fno = 65535, rno = 65535, k, i, j;
    uint16_t fwdCount = 0;
    uint16_t revCount = 0;
    uint32_t a, r1=NOTFOUND, f1=NOTFOUND;
    uint32_t r2=r1, f2=f1, r3=r1, f3=f1;
    uint16_t fi1=0, ri1=0, fi2=0, ri2=0, fi3=0, ri3=0;
    uint32_t nK;
    uint32_t nK2;
    uint16_t maxTargets = 200; 
    uint16_t unused = 0;
    

    
    
    pMap->pCandPosMan->reset();
    
    nK = pMap->getNFwdPos(); 
    nK2 = nK >> 1;
    if (pMap->minFwdKeySize <= maxTargets) {
        for (k=0; pMap->nActiveFwd > unused; k++) { // prev condition && k < 10 && 
            f1 = f2 = NOTFOUND;
            fi1= 0;
            fi2= nK2;
            unused = 0;
            for (i=0; i < nK2; i++) {
                // left part
                if (pMap->isFwdKeyActive(i)) {
                    a = pMap->getFwdTargetSize(i);
                    if (a < maxTargets && a > 0 && a < f1) {
                        f1 = a;
                        fi1 = i;
                    } else {
                        unused++;
                    }
                }
            }
            for (j=nK2; j < nK; j++) {
                // right part
                //j = i + nK2;
                if (pMap->isFwdKeyActive(j)) {
                    a = pMap->getFwdTargetSize(j);
                    if (a < maxTargets && a > 0 && a < f2 && a != f1) {
                        f2 = a;
                        fi2 = j;
                    } else {
                        unused++;                    
                    }
                }
            }
            if (fwdCount == 0) {
                if (f1 == NOTFOUND || f2 == NOTFOUND) break;
                fwdCount = pMap->intersectByTwoReadPositionKeys(0, fi1, fi2, pMap->getFwdKeyablePosition(fi2) - pMap->getFwdKeyablePosition(fi1) + 3);
                if (fwdCount == 0) {
                    if (pMap->getFwdTargetSize(fi1) < pMap->getFwdTargetSize(fi2)) {
                        pMap->inactivePositionFwd(fi1);
                    } else {
                        pMap->inactivePositionFwd(fi2);
                    }
                } else {
                    pMap->inactivePositionFwd(fi1);
                    pMap->inactivePositionFwd(fi2);
                }
            } else {
                if (f1 != NOTFOUND) {
                    fwdCount = pMap->intersectByAddingReadPositionKey(0, fi1, 1000);
                    pMap->inactivePositionFwd(fi1);
                    if (fwdCount == 0) {
                        fwdCount = pMap->pCandPosMan->rewindCountsAfterAddingKeyCountZero();
                    }
                }
                if (f2 != NOTFOUND) {
                    fwdCount = pMap->intersectByAddingReadPositionKey(0, fi2, 1000);
                    pMap->inactivePositionFwd(fi2);
                    if (fwdCount == 0) {
                        fwdCount = pMap->pCandPosMan->rewindCountsAfterAddingKeyCountZero();
                    }
                }
                if (f1 == NOTFOUND && f2 == NOTFOUND) break;
            }
            if (fwdCount > 0 && fwdCount < 5) {
                break;
            }
        }
    }

    pMap->pCandPosMan->saveFwdCandidatePositions();

    if (pMap->minRevKeySize <= maxTargets) {
        unused = 0;
        for (k=0; pMap->nActiveRev > unused; k++) { // prev cond k < 10 && 
            r1 = r2 = NOTFOUND;
            ri1= 0;
            ri2= nK2;
            nK = pMap->getNRevPos(); 
            nK2 = nK >> 1;
            unused = 0;
            for (i=0; i < nK2; i++) {
                // Left part
                if (pMap->isRevKeyActive(i)) {
                    a = pMap->getRevTargetSize(i);
                    if (a < maxTargets && a > 0 && a < r1) {
                        r1 = a;
                        ri1 = i;
                    } else {
                        unused++;                    
                    }
                }
                //j = i + nK2;
            }
            for (j=nK2; j < nK; j++) {
                // Right part
                if (pMap->isRevKeyActive(j)) {
                    a = pMap->getRevTargetSize(j);
                    if (a < maxTargets && a > 0 && a < r2 && a != r1) {
                        r2 = a;
                        ri2 = j;
                    } else {
                        unused++;
                    }
                }
            }
            if (revCount == 0) {
                if (r1 == NOTFOUND || r2 == NOTFOUND) break;
                revCount = pMap->intersectByTwoReadPositionKeys(1, ri1, ri2, pMap->getRevKeyablePosition(ri2) - pMap->getRevKeyablePosition(ri1) + 3);
                if (revCount == 0) {
                    if (pMap->getRevTargetSize(ri1) < pMap->getRevTargetSize(ri2)) {
                        pMap->inactivePositionRev(ri1);
                    } else {
                        pMap->inactivePositionRev(ri2);
                    }
                } else {
                    pMap->inactivePositionRev(ri1);
                    pMap->inactivePositionRev(ri2);
                }
            } else {
                if (r1 != NOTFOUND) {
                    revCount = pMap->intersectByAddingReadPositionKey(1, ri1, pMap->read->lenSeq);
                    pMap->inactivePositionRev(ri1);
                    if (revCount == 0) {
                        revCount = pMap->pCandPosMan->rewindCountsAfterAddingKeyCountZero();
                    }
                }
                if (r2 != NOTFOUND) {
                    revCount = pMap->intersectByAddingReadPositionKey(1, ri2, pMap->read->lenSeq);
                    pMap->inactivePositionRev(ri2);
                    if (revCount == 0) {
                        revCount = pMap->pCandPosMan->rewindCountsAfterAddingKeyCountZero();
                    }
                }
                if (r1 == NOTFOUND && r2 == NOTFOUND) break;
            }
            if (revCount > 0 && revCount < 5) {
                break;
            }
        }
    }
    
    if (revCount == 0 && fwdCount == 0) {
        //fprintf(stderr, "** PROBLEM Zero Counts iProc=%llu **\n", iProc);
        //fprintf(stderr, "%s\n", pMap->seq);
    }
    if (revCount > 20 || fwdCount > 20) {
        //fprintf(stderr, "** PROBLEM Large Counts iProc=%llu, f=%u, r=%u **\n", iProc, fwdCount, revCount);
        //fprintf(stderr, "%s\n", pMap->seq);        
    }
    //pMap->pCandPosMan->setIsReverseComplement(revCount >= fwdCount ? 1 : 0);
    //return  fwdCount + revCount;
    return pMap->pCandPosMan->getCount();
}


uint16_t AnalyzeAllPos(ogReadKeyMapping *pMap) {
//    This method for allocating specific space is NICE!!!! but now, CandidatePositionManager contains a BinaryTreeObject pTree
//    static ogReusableBinaryTree **memAlloc = (ogReusableBinaryTree **) calloc(512, sizeof(ogReusableBinaryTree)); // 512 threads maximum
    ogReusableBinaryTree *pTree = pMap->pCandPosMan->pTree; // memAlloc[pMap->thread];
    
//    if (pTree == NULL) {
//        // this should not run now that pCandPosMan allocate the tree, but for security ... 
//        pTree = new ogReusableBinaryTree();
//        //memAlloc[pMap->thread] = pTree;
//        pMap->pCandPosMan->pTree = pTree;
//    }
    
    pTree->clear();
    pTree->setMaximumDistanceForKey(10);
    
    const uint16_t NOTFOUND = 65535;
    int32_t  node, val;
    int32_t  max_F, maxNode_F, max2_F, max2Node_F;
    int32_t  max_R, maxNode_R, max2_R, max2Node_R;
    uint32_t maxFwdSize = 100; //pMap->minFwdKeySize * 2; // Parametro (maximo tamaño de una llave para hacer el conteo)
    uint32_t maxRevSize = 100; //pMap->minRevKeySize * 2;
    uint32_t k, i, j, L;
    uint16_t fwdCount = 0, f1=0, fi1=NOTFOUND, r1=0, ri1=NOTFOUND;
    uint16_t revCount = 0;
    uint32_t nKF = pMap->getNFwdPos(), a; // Possible problem with rev pos
    uint32_t nK2 = nKF >> 1;
    uint32_t nKR;
    uint32_t fwdPos=0, revPos=0, fwdPos2=0, revPos2=0;
    ogKey *pK;
    uint32_t *pGP;
    char flag;
    if (maxFwdSize < 100) maxFwdSize = 100;
    if (maxRevSize < 100) maxRevSize = 100;
    if (maxFwdSize > 1000 || maxRevSize > 1000) {
        //fprintf(stderr, "[%u %u iP=%llu] ", maxFwdSize, maxRevSize, iProc);
    }

//    pMap->pCandPosMan->reset();
    
    // Get max positions
    max_F = max2_F = maxNode_F = max2Node_F = 0;    
    val = 0x7FFFFFFF; //val = pMap->getFwdInfoForKeyablePosition(0)->size;
    for (i=0; i < nKF; i++) {
        a = pMap->getFwdTargetSize(i);
        if (a < maxFwdSize && a > f1) {
            f1 = a;
            fi1 = i;
        } if (a < val) {
            val = a;
            node = i;
        }
    }
//    if (fi1 == NOTFOUND && val < maxFwdSize * 2) {
//        maxFwdSize = val * 2;        
//        fi1 = node;
//    }
//    if (fi1 == NOTFOUND) {
//        f1 = NOTFOUND;
//        for (i=0; i < nK; i++) {
//            a = pMap->getFwdInfoForKeyablePosition(i)->size;
//            if (a > 0 && a < f1) {
//                f1 = a;
//                fi1 = i;
//            }
//        }
//        maxFwdSize = f1 * 10;
//    }
    if (fi1 != NOTFOUND) {
        //pK = pMap->getFwdInfoForKeyablePosition(fi1);
        //pGP = pMap->pKeys->getPointerPositionForKey(pMap->getKeyForward(fi1));
        L = pMap->getFwdTargetSize(fi1);
        pGP = pMap->getGenomePositionFromReadPositionKey(0,fi1,0);
        k = pMap->getFwdKeyablePosition(fi1);
        for (i=L; i; ) {
            i--;
            pTree->setNodeTo(i, pGP[i]-k, 1);
        }
        pTree->buildTreeFromOrderedNodesSet(L);
        for (i=0; i < nKF; i++) {
            if (i != fi1) {
                //pK = pMap->getFwdInfoForKeyablePosition(i);
                L = pMap->getFwdTargetSize(i);
                if (L < maxFwdSize && L > 0) {
                    //pGP = pMap->pKeys->getPointerPositionForKey(pMap->getKeyForward(i));
                    pGP = pMap->getGenomePositionFromReadPositionKey(0,i,0);
                    k = pMap->getFwdKeyablePosition(i);
                    flag = 0;
                    for (j=L-1; j; j--) {
                        if ((node = pTree->insertNode(pGP[j]-k, 1)) >= 0) {
                            val = pTree->incNodeContent(node);
                            if (val >= max_F) {
                                if (node != maxNode_F) {
                                    max2_F = max_F;
                                    max2Node_F = maxNode_F;
                                    maxNode_F = node;
                                }
                                max_F = val;
                            } else if (val > max2_F) {
                                max2_F = val;
                                max2Node_F = node;
                            }
                            if (val > 2) {
                                flag = 1;
                            }
                        }
                    }
                    if (flag) {
                        // Se encontraron 3 keys con match, salir agiliza pero no es exahustivo
                        break;
                    }
                }
            }
        }
        fwdCount = max_F;
        fwdPos = pTree->getNodeKey(maxNode_F);
        fwdPos2 = pTree->getNodeKey(max2Node_F);
    } else {
        //fprintf(stderr, " Min Size=%d ",val);
    }


    pTree->clear();
    // Get max positions in REVERSE
    nKR = pMap->getNRevPos();
    nK2 = nKR >> 1;
    max_R = max2_R = maxNode_R = max2Node_R = 0;    
    val = 0x7FFFFFFF;//pMap->getRevInfoForKeyablePosition(0)->size;
    for (i=0; i < nKR; i++) {
        a = pMap->getRevInfoForKeyablePosition(i)->size;
        if (a < maxRevSize && a > r1) {
            r1 = a;
            ri1 = i;
        } else if (a < val) {
            val = a;
            node = i;
        }
    }
//    if (ri1 == NOTFOUND && val < maxRevSize * 2) {
//        maxRevSize = val * 2;
//        ri1 = node;
//    }
//    if (ri1 == NOTFOUND) {
//        r1 = NOTFOUND;
//        for (i=0; i < nKR; i++) {
//            a = pMap->getRevInfoForKeyablePosition(i)->size;
//            if (a > 0 && a < r1) {
//                r1 = a;
//                ri1 = i;
//            }
//        }
//        maxRevSize = r1 * 10;
//    }
    if (ri1 != NOTFOUND) {
        //pK = pMap->getRevInfoForKeyablePosition(ri1);
        //pGP = pMap->pKeys->getPointerPositionForKey(pMap->getKeyReverse(ri1));
        L = pMap->getRevTargetSize(ri1);
        pGP = pMap->getGenomePositionFromReadPositionKey(1,ri1,0);
        k = pMap->getRevKeyablePosition(ri1);
        for (i=L; i; ) {
            i--;
            pTree->setNodeTo(i, pGP[i]-k, 1);
        }
        pTree->buildTreeFromOrderedNodesSet(L);
        for (i=0; i < nKR; i++) {
            if (i != ri1) {
                //pK = pMap->getRevInfoForKeyablePosition(i);
                L = pMap->getRevTargetSize(i);
                if (L < maxRevSize) {
                    //pGP = pMap->pKeys->getPointerPositionForKey(pMap->getKeyReverse(i));
                    pGP = pMap->getGenomePositionFromReadPositionKey(1,i,0);
                    k = pMap->getRevKeyablePosition(i);
                    flag = 0;
                    for (j=L; j; ) {
                        j--;
                        if ((node = pTree->insertNode(pGP[j]-k, 1)) >= 0) {
                            val = pTree->incNodeContent(node);
                            if (val >= max_R) {
                                if (node != maxNode_R) {
                                    max2_R = max_R;
                                    max2Node_R = maxNode_R;
                                    maxNode_R = node;
                                }
                                max_R = val;
                            } else if (val > max2_R) {
                                max2_R = val;
                                max2Node_R = node;
                            }
                            if (val > 2) {
                                flag = 1;
                            }
                        }
                    }
                    if (flag) {
                        // Se encontraron 3 keys con match, salir agiliza pero no es exahustivo
                        break;
                    }
                }
            }
        }
        revCount = max_R;
        revPos = pTree->getNodeKey(maxNode_R);
        revPos2 = pTree->getNodeKey(max2Node_R);
    }
    
    // Here revCount & fwdCount refers to the total number of 
    // keys referring to the position
    // (instead of total number of possible positions)
    if (revCount == 0 && fwdCount == 0) {
        //fprintf(stderr, "** PROBLEM Zero Counts iProc=%llu **\n", iProc);
        //fprintf(stderr, "%s\n", pMap->seq);
        return 0;
    }
    /**
    if (revCount > 20 || fwdCount > 20) {
        fprintf(stderr, "** PROBLEM Large Counts iProc=%llu, f=%u, r=%u **\n", iProc, fwdCount, revCount);
        fprintf(stderr, "%s\n", pMap->seq);        
    }
     **/
    pMap->pCandPosMan->reset();
    if (fwdCount > 2 || revCount > 2) {
        if (revCount > fwdCount || (revCount == fwdCount && nKF > nKR)) { // if == take that of less number of keys
            pMap->pCandPosMan->saveFwdCandidatePositions();
            pMap->pCandPosMan->addGenomeAndReadPosition(revPos, 0, 1, 0);//pMap->pCandPosMan->setByFixedPosition(1, revPos);
            if ((revCount >> 1) < max2_R  && max2_R > 4) {
                // Possible 2 hits
                pMap->pCandPosMan->addGenomeAndReadPosition(revPos2, 0, 1, 0); //pMap->pCandPosMan->addFixedPosition(1, revPos2);
                return 2;
            }
            return  1;
        } else {
            pMap->pCandPosMan->addGenomeAndReadPosition(fwdPos, 0, 0, 0); //pMap->pCandPosMan->setByFixedPosition(0, fwdPos);
            if ((fwdCount >> 1) < max2_F  && max2_F > 4) {
                // Possible 2 hits
                pMap->pCandPosMan->addGenomeAndReadPosition(fwdPos2, 0, 0, 0); //pMap->pCandPosMan->addFixedPosition(0, fwdPos2);
                pMap->pCandPosMan->saveFwdCandidatePositions();
                //fprintf(stderr, "** 2 Hits **");
                return 2;
            }
            pMap->pCandPosMan->saveFwdCandidatePositions();
            return  1;
        }
    } else {
        return 0;
    }
}





// Extend the 
uint16_t MovWinIntersect(ogReadKeyMapping *pMap, uint16_t window, uint16_t steps, uint32_t maxTargets, char doDynamicSteps) {
    uint32_t i, j, k;
    ogKey *pKey;
    char RC = 0;
    uint32_t maxExtFR[2] = {0, 0}, maxExtPosFR[2] = {0,0}, maxExtFREq[2] = {0,0}, maxGenPosFR[2] = {0,0}, maxCountsFR[2]={0,0}, maxRelPosFR[2] = {0,0};
    uint32_t maxExt = 0, maxPos = 0, maxEq = 0, maxGenPos, maxCounts, maxRelPos;
    uint32_t intC, iPos, jPos, limit, next_i;
    char flagTwoPos = 0;
    uint16_t searching = 0;
    
    for (RC = 0; RC < 2; RC++) {
        pMap->pCandPosMan->reset();
        uint32_t lastKeyPos = (RC == 0 ? pMap->getNFwdPos() : pMap->getNRevPos());
        if (lastKeyPos > 1) {
            for (i=0; i < lastKeyPos-1; ) {
                next_i = 0;
                if (pMap->getTargetSize(RC, i) < maxTargets) {
                    maxExt = 0;
                    maxCounts = 0;
                    maxGenPos = 0;
                    maxRelPos = 0;
                    intC = 0;
                    iPos = pMap->getKeyablePosition(RC, i);
                    limit = i + window;
                    if (limit > lastKeyPos) limit = lastKeyPos;
                    for (j = i+1; j < lastKeyPos; j++) {
                        if (pMap->getTargetSize(RC, j) < maxTargets) {
                            jPos = pMap->getKeyablePosition(RC, j);
                            if (intC == 0) {
                                intC = pMap->intersectByTwoReadPositionKeys(RC, i, j, jPos - iPos + 5);
                                if (intC == 0 && next_i == 0) {
                                    next_i = j;
                                }
                            } else {
                                intC = pMap->intersectByAddingReadPositionKey(RC, j, jPos - iPos + 5);
                                if (intC == 0) {
                                    pMap->pCandPosMan->rewindCountsAfterAddingKeyCountZero();
                                    if (next_i == 0) {
                                        next_i = j;
                                    }
                                }
                            }
                            if (intC > 0) {
                                // intersection works, so
                                if (jPos - iPos > maxExt) {
                                    maxExt = jPos - iPos;
                                    maxCounts = intC;
                                    maxGenPos = pMap->pCandPosMan->getkPos(0)->genomePosition - iPos; // aqui supone el mejor hit es el primero
                                    maxRelPos = iPos;
                                }
                            }
                        }
                    }
                    if (maxExt >= maxExtFR[RC]) {
                        if (maxExt == maxExtFR[RC]) {
                            maxExtFREq[RC]++;
                        } else {
                            maxCountsFR[RC] = maxCounts;
                            maxGenPosFR[RC] = maxGenPos;
                            maxExtFR[RC] = maxExt;
                            maxExtPosFR[RC] = i;
                            maxExtFREq[RC] = 0;
                            maxRelPosFR[RC] = maxRelPos;
                        }
                    }
                }
                if (doDynamicSteps && next_i > 0) {
                    i = next_i;
                } else {
                    i += steps;
                }
            }
        }  
    }

    if (maxExtFR[0] == maxExtFR[1] && maxExtFR[0] == 0) {
        return 0;
    }
    
    // Este algoritmo hay que ajustarlo al nuevo pCandPosMan porque supone que el 1er hit es el mejor y por supuesto le falla mucho

    RC = (maxExtFR[0] >= maxExtFR[1] ? 0 : 1);
    pMap->pCandPosMan->reset();
    if (RC) pMap->pCandPosMan->saveFwdCandidatePositions();
    pMap->pCandPosMan->addGenomeAndReadPosition(maxGenPosFR[RC], 0, RC, 0); //pMap->pCandPosMan->setByFixedPosition(RC, maxGenPosFR[RC]);
    if (RC==0) pMap->pCandPosMan->saveFwdCandidatePositions();
    
    return maxCountsFR[RC];
}

uint16_t MovWinIntersectMicro(ogReadKeyMapping *pMap) {
    return MovWinIntersect(pMap, 5, 3, pMap->getLowKeyCountLimit(), 0);
}

uint16_t MovWinIntersectSmall(ogReadKeyMapping *pMap) {
    return MovWinIntersect(pMap, 7, 2, pMap->getLowKeyCountLimit(), 0);
}

uint16_t MovWinIntersectLarge(ogReadKeyMapping *pMap) {
    return MovWinIntersect(pMap, 7, 5, pMap->getHighKeyCountLimit(), 0);
}

uint16_t MovWinIntersectHuge(ogReadKeyMapping *pMap) {
    return MovWinIntersect(pMap, 20, 5, 10*pMap->getHighKeyCountLimit(), 1);
}



// leftmost Intersection with rightmost keys ....
uint16_t ExtremesMapping(ogReadKeyMapping *pMap, uint32_t maxTargets) {

    uint32_t nReadPos, halfPos;
    uint32_t i, iPos, intC;
    uint32_t j;
    uint8_t  intraFails;
    char    x, bug = pMap->debug;
    
    pMap->pCandPosMan->reset();
    
    for (x = 0; x < 2; x++) {
        nReadPos = (x == 0 ? pMap->getNFwdPos() : pMap->getNRevPos());
        halfPos = (nReadPos >> 1);

        intC = 0;
        
        for (i=0; i < halfPos; i++) {
            if (pMap->getTargetSize(x, i) < maxTargets) {
                iPos = pMap->getKeyablePosition(x, i);
                intC = 0;
                for (j=nReadPos-1; j >= halfPos; j--) {
                    if (pMap->getTargetSize(x, j) < maxTargets) {
                        if (intC == 0) {
                            intC = pMap->intersectByTwoReadPositionKeys(x, i, j, pMap->getKeyablePosition(x,j) - iPos + 3);
                            if (intC > 10 || intC == 0) {
                                intC = 0;
                                if (++intraFails > 3) break;                                
                            }
                        } else {
                            intC = pMap->intersectByAddingReadPositionKey(x, j, pMap->getKeyablePosition(x,j) - iPos + 3);
                            if (intC == 0) {
                                intC = pMap->pCandPosMan->rewindCountsAfterAddingKeyCountZero();
                                if (++intraFails > 3) break;                                
                            }
                        }
                        if (intC > 0 && intC < 3) {
                            break;
                        }
                    }
                }
                if (intC > 0 && intC < 3) {
                    break;
                }
            }
        }
        if (x == 0) {
            if (intC == 0) { // (pMap->pCandPosMan->getCount() > 2)
                pMap->pCandPosMan->reset();
            }
            pMap->pCandPosMan->saveFwdCandidatePositions();
        } else {
            if (intC == 0) { // (pMap->pCandPosMan->getRevCount() > 2)
                pMap->pCandPosMan->restoreFwdCandidatePositions();
            }
            //if (intC > 2) {
            //    pMap->pCandPosMan->removeReverseTargets();
            //}
        }
    }

    return pMap->pCandPosMan->getCount();
}

uint16_t ExtremeMatchingSmall(ogReadKeyMapping *pMap) {
    return ExtremesMapping(pMap,  pMap->getLowKeyCountLimit());
}

uint16_t ExtremeMatchingLarge(ogReadKeyMapping *pMap) {
    return ExtremesMapping(pMap,  pMap->getHighKeyCountLimit());
}

uint16_t ExtremeMatchingHuge(ogReadKeyMapping *pMap) {
    return ExtremesMapping(pMap, 10*pMap->getHighKeyCountLimit());
}

typedef struct frCount {
    uint16_t     counts[2];
} FR_COUNT;



// count all keys
uint16_t countGeneKeys(ogReadKeyMapping *pRdMap, uint32_t maxTargets) {

    ogReadKeyMapping    *pMap = pRdMap;
    char            x, y, z, w, RCkeys;
    uint16_t        bestTie, bestFwdTie, bestRevTie;
    ogChromosome *pChr;
    ogGenome    *pGenome = pMap->pGenome;
    uint32_t    nGenes = pGenome->getNChromosomes();
    uint32_t    nKeys, i, nValidKeys, nTotalKeys;
    ogCandidatePosManager *pCand = pMap->pCandPosMan;
    ogKey       *pKey;
    uint32_t    *pGenPos, *pLimGenPos, bc, bestCount, bestGene, prevIndex, prevGene, bestFwdCount, bestFwdGene, fwdKeys, revKeys, bestRevGene, bestRevCount, bestTieGene, bestTieGeneFwd, bestTieGeneRev;
    //uint32_t    nBlocks = nGenes >> 7, iBlock;
    uint32_t    bytes = nGenes * sizeof(FR_COUNT); // + nBlocks; ? // + 128;
    FR_COUNT    *pMem = (FR_COUNT *) pMap->getUsableMemory(bytes);
    Exon        *pXONS = pMap->pCountInfo->pExons;
    uint32_t    iGene, iRdPos;
    
    pCand->reset();
    
    if ((pMap->read->readNumber & 0xFFF) == 0) memset(pMem, 0, bytes); // For security
    //memset(pMem, 0, bytes);
    //memset(pMem, 0, 128);
    
    w = (pRdMap->pLinkedRKM == NULL ? 1 : 2);
    
    for (z = 0; z < 2; z++) {
        // z == 0 : Count
        // z == 1 : Erase counts
        for (x = 0; x < 2; x++) {
            // x == 0 - Para Fwd read del Primer read
            // x == 1 - Para Rev read del Primer read
            //memset(pMem, 0, bytes); // no usar memset indica que Asumimos que los conteos de Fwd y Rev no se mezclan, para ahorrar tiempo
            bestCount = 0;
            bestGene = 0;
            nValidKeys = 0;
            nTotalKeys = 0;
            bestTie = 0;
            bestTieGene = 0;
            for (y = 0; y < w; y++) {
                // x == 0 && y == 0 Maps RKM1_Fwd 
                // x == 0 && y == 1 Maps RKM2_Rev
                // x == 1 && y == 0 Maps RKM1_Rev
                // x == 1 && y == 1 Maps RKM2_Fwd

                RCkeys = (x == y ? 0                     : 1);
                pMap   = (y      ? pRdMap->pLinkedRKM : pRdMap);
                nKeys  = (RCkeys ? pMap->getNRevPos() : pMap->getNFwdPos());
                if (nKeys > 0) {
                    nTotalKeys += nKeys;
                    for (i=0; i < nKeys; i++) {
                        pKey = pMap->getInfoForPositionKey(RCkeys, i);
                        if (pKey->size < maxTargets) {
                            nValidKeys++;
                            pGenPos = pMap->getGenomePositionFromReadPositionKey(RCkeys, i, 0);
                            prevIndex = pMap->pCountInfo->nExons + 1;
                            prevGene  = nGenes + 1;
                            for (pLimGenPos = pGenPos + pKey->size; pGenPos < pLimGenPos; pGenPos++) {
                                if (z) {
                                    pMem[pXONS[*pGenPos].iGene].counts[x] = 0;
                                } else {
                                    if (*pGenPos != prevIndex) {
                                        prevIndex = *pGenPos;
                                        iGene = pXONS[prevIndex].iGene;
                                        if (iGene != prevGene) {
                                            prevGene = iGene;
                                            if ((bc=++pMem[iGene].counts[x]) >= bestCount) {
                                                if (bc > bestCount) {
                                                    bestCount = pMem[iGene].counts[x];
                                                    bestGene = iGene;
                                                    bestTie = 0;
                                                } else {
                                                    bestTieGene = iGene;
                                                    bestTie++;
                                                }
                                                //if (iGene == 15022) {
                                                //    char *pS = (RCkeys ? pMap->read->pSeqRevComp : pMap->read->pSeq);
                                                //    fprintf(stderr, "RN:%llu,i=%u,exon=%u||||||||||||||%s\n", pMap->read->readNumber, i, prevIndex, pS+pMap->getKeyablePosition(RCkeys, i));
                                                //}
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (z == 0) {
                if (x == 0) {
                    bestFwdCount = bestCount;
                    bestFwdGene = bestGene;
                    fwdKeys = nValidKeys;
                    bestFwdTie = bestTie;
                    bestTieGeneFwd = bestTieGene;
                } else {
                    bestRevCount = bestCount;
                    bestRevGene = bestGene;
                    revKeys = nValidKeys;
                    bestRevTie = bestTie;
                    bestTieGeneRev = bestTieGene;
                }
            }
            //if (bestCount > 0) {
            //    pChr = pGenome->getChromosome(bestGene);
            //    //pCand->addGenomeAndReadPosition(pChr->start+1, 0, x, bestCount * 100 / nValidKeys);
            //    fprintf(stderr, "ReadNumber=%llu, Len=%u, Map=%c, RC=%c, i=%u, Gene=%s, Count=%u of %u from %u keys.\n", pMap->read->readNumber, pMap->read->lenSeq, y+48, x+48, pChr->index, pChr->name, bestCount, nValidKeys, nTotalKeys); fflush(stderr);
            //}
        }
    }
    // FALTA MANEJAR TIES
    // FALTA CONSIDERAR QUE EL RATIO DE COUNTS / KEYS > 50%
//    if (bestRevGene == 15022 && bestRevCount > bestFwdCount || bestFwdGene == 15022 && bestFwdCount > bestRevCount) { // bestFwdCount <= 2 && bestCount <= 2 // bestFwdCount == bestCount && bestCount > 4 || 
//        pChr = pGenome->getChromosome(bestFwdGene);
//        fprintf(stderr, "Read=%llu, L=%u, RC=%c, %u=%s, C=%u/%u:%u K %s  %s\n", pMap->read->readNumber, pMap->read->lenSeq, '0', pChr->index, pChr->name, bestFwdCount, fwdKeys, nTotalKeys, pMap->read->pSeq, pMap->pLinkedRKM->read->pSeqRevComp); fflush(stderr);
//        pChr = pGenome->getChromosome(bestRevGene);
//        fprintf(stderr, "Read=%llu, L=%u, RC=%c, %u=%s, C=%u/%u:%u K %s  %s\n", pMap->read->readNumber, pMap->read->lenSeq, '1', pChr->index, pChr->name, bestRevCount, revKeys, nTotalKeys, pMap->read->pSeqRevComp, pMap->pLinkedRKM->read->pSeq); fflush(stderr);
//    }
    if (bestFwdCount > bestRevCount) {
        // Fwd
        z = bestFwdCount * 100 / fwdKeys;
        if (bestFwdTie == 1 && pRdMap->pMapParams->BREAK_TIE_COUNT1) {
            // La mayoría de los ties son con solo 1 gen, entonces se puede romper fácil tomandolo aleatorio, esto aumenta 3% de los reads (recupera el 15% de los reads)
            if (pRdMap->read->readNumber & 0x1) bestFwdGene = bestTieGeneFwd;
            bestFwdTie = 0;
        }
        if (bestFwdTie == 0 && bestFwdCount >= pRdMap->pMapParams->MIN_KEY_COUNTS && z >= pRdMap->pMapParams->MIN_SCORE_COUNT) {
            pCand->addGenomeAndReadPosition(bestFwdGene, 0, 0, z);
        }
        //if (bestFwdTie > 0) {
        //    fprintf(stderr, "[%u:%hu:%u]", bestFwdCount, bestFwdTie, maxTargets);
        //}
    } else if (bestRevCount > bestFwdCount) {
        // Reverse
        z = bestRevCount * 100 / revKeys;
        if (bestRevTie == 1 && pRdMap->pMapParams->BREAK_TIE_COUNT1) {
            // La mayoría de los ties son con solo 1 gen, entonces se puede romper fácil tomandolo aleatorio, esto aumenta 3% de los reads (recupera el 15% de los reads)
            if (pRdMap->read->readNumber & 0x1) bestRevGene = bestTieGeneRev;
            bestRevTie = 0;
        }
        if (bestRevTie == 0 && bestRevCount >= pRdMap->pMapParams->MIN_KEY_COUNTS && z >= pRdMap->pMapParams->MIN_SCORE_COUNT) {
            pCand->addGenomeAndReadPosition(bestRevGene, 0, 1, z);
        }
        //if (bestRevTie > 0) {
        //    fprintf(stderr, "[%u:%hu:%u]", bestRevCount, bestRevTie, maxTargets);
        //}
//    } else {
//        pChr = pGenome->getChromosome(bestFwdGene);
//        fprintf(stderr, "Read=%llu, L=%u, RC=%c, %u=%s, C=%u/%u:%u K %s  %s\n", pMap->read->readNumber, pMap->read->lenSeq, '0', pChr->index, pChr->name, bestFwdCount, fwdKeys, nTotalKeys, pMap->read->pSeq, pMap->pLinkedRKM->read->pSeqRevComp); fflush(stderr);
//        pChr = pGenome->getChromosome(bestRevGene);
//        fprintf(stderr, "Read=%llu, L=%u, RC=%c, %u=%s, C=%u/%u:%u K %s  %s\n", pMap->read->readNumber, pMap->read->lenSeq, '1', pChr->index, pChr->name, bestRevCount, revKeys, nTotalKeys, pMap->read->pSeqRevComp, pMap->pLinkedRKM->read->pSeq); fflush(stderr);
    }
//    if (pCand->getCount() == 0) {
//        pChr = pGenome->getChromosome(bestFwdGene);
//        fprintf(stderr, "Read=%llu, L=%u, RC=%c, %u=%s, C=%u/%u:%u K Tie=%u %s  %s\n", pMap->read->readNumber, pMap->read->lenSeq, '0', pChr->index, pChr->name, bestFwdCount, fwdKeys, nTotalKeys, bestFwdTie, pMap->read->pSeq, pMap->pLinkedRKM->read->pSeqRevComp); fflush(stderr);
//        pChr = pGenome->getChromosome(bestRevGene);
//        fprintf(stderr, "Read=%llu, L=%u, RC=%c, %u=%s, C=%u/%u:%u K Tie=%u %s  %s\n", pMap->read->readNumber, pMap->read->lenSeq, '1', pChr->index, pChr->name, bestRevCount, revKeys, nTotalKeys, bestRevTie, pMap->read->pSeqRevComp, pMap->pLinkedRKM->read->pSeq); fflush(stderr);
//    }
    //pCand->printWithGenomePositions(pMap->pGenome);
    return pCand->getCount();
}


uint16_t CountGeneReadSmall(ogReadKeyMapping *pMap) {
    return countGeneKeys(pMap,  pMap->getLowKeyCountLimit());
}

uint16_t CountGeneReadLarge(ogReadKeyMapping *pMap) {
    return countGeneKeys(pMap,  pMap->getHighKeyCountLimit());
}

uint16_t CountGeneReadHuge(ogReadKeyMapping *pMap) {
    return countGeneKeys(pMap, 10*pMap->getHighKeyCountLimit());
}





// count all matching 
// NO FUNCIONA PORQUE EL GENOME POSITION ESTA INDEXADO A LA SECUENCIA SINO AL EXON/GEN
// Y POR LO TANTO NO REFERENCIA A LA SECUENCIA
uint16_t countGeneKeysBlasting(ogReadKeyMapping *pRdMap, uint32_t maxTargets) {

    ogReadKeyMapping    *pMap = pRdMap;
    char            x, y, z, w, RCkeys;
    uint16_t        bestTie, bestFwdTie, bestRevTie;
    ogChromosome *pChr;
    ogGenome    *pGenome = pMap->pGenome;
    uint32_t    nGenes = pGenome->getNChromosomes();
    uint32_t    nKeys, i, nValidKeys, nTotalKeys;
    ogCandidatePosManager *pCand = pMap->pCandPosMan;
    ogKey       *pKey;
    ogSingleRead *pRead;
    uint32_t    *pGenPos, *pLimGenPos, bc, bestCount, bestGene, prevIndex, prevGene, bestFwdCount, bestFwdGene, fwdKeys, revKeys, bestRevGene, bestRevCount, bestTieGene, bestTieGeneFwd, bestTieGeneRev;
    Exon        *pXONS = pMap->pCountInfo->pExons;
    uint32_t    iGene, iRdPos;
    uint32_t    rdLen = pMap->read->lenSeq;
    uint32_t    errores, matches;
    uint32_t    bytes = nGenes * sizeof(FR_COUNT); // + nBlocks; ? // + 128;
    FR_COUNT    *pMem = (FR_COUNT *) pMap->getUsableMemory(bytes);
    
    pCand->reset();

    if ((pMap->read->readNumber & 0xFFF) == 0) memset(pMem, 0, bytes); // For security
       
    w = (pRdMap->pLinkedRKM == NULL ? 1 : 2);
    for (z = 0; z < 2; z++) {
        // z == 0 : Count
        // z == 1 : Erase counts
        for (x = 0; x < 2; x++) {
            // x == 0 - Para Fwd read del Primer read
            // x == 1 - Para Rev read del Primer read
            bestCount = 0; // Errors
            bestGene = 0;
            nValidKeys = 0;
            nTotalKeys = 0;
            bestTie = 0;
            bestTieGene = 0;
            for (y = 0; y < w; y++) {
                // x == 0 && y == 0 Maps RKM1_Fwd 
                // x == 0 && y == 1 Maps RKM2_Rev
                // x == 1 && y == 0 Maps RKM1_Rev
                // x == 1 && y == 1 Maps RKM2_Fwd

                RCkeys = (x == y ? 0                     : 1);
                pMap   = (y      ? pRdMap->pLinkedRKM : pRdMap);
                nKeys  = (RCkeys ? pMap->getNRevPos() : pMap->getNFwdPos());
                if (nKeys > 0) {
                    pRead = pMap->read;
                    nTotalKeys += nKeys;
                    for (i=0; i < nKeys; i++) {
                        pKey = pMap->getInfoForPositionKey(RCkeys, i);
                        if (pKey->size < maxTargets) {
                            nValidKeys++;
                            iRdPos = pMap->getKeyablePosition(RCkeys, i);
                            pGenPos = pMap->getGenomePositionFromReadPositionKey(RCkeys, i, 0);
                            prevIndex = pMap->pCountInfo->nExons + 1;
                            prevGene  = nGenes + 1;
                            for (pLimGenPos = pGenPos + pKey->size; pGenPos < pLimGenPos; pGenPos++) {
                                if (z) {
                                    pMem[pXONS[*pGenPos].iGene].counts[x] = 0;
                                } else {
                                    if (*pGenPos != prevIndex) {
                                        prevIndex = *pGenPos;
                                        if (iRdPos < prevIndex) {
        //                                    fprintf(stderr, "Read=%llu, L=%u, RC=%c, K %s  %s\n", pMap->read->readNumber, pMap->read->lenSeq, RCkeys, pMap->read->pSeq, pMap->pLinkedRKM->read->pSeqRevComp); fflush(stderr);
        //                                    fprintf(stderr, "prevIndex=%u, pos=%u\n", prevIndex, iRdPos); fflush(stderr);
                                            errores = pGenome->comparePackedSequences(pRead->packReadSeqIfNeeded(pGenome->getPackedOffset(prevIndex - iRdPos), RCkeys), prevIndex - iRdPos, rdLen, rdLen);
                                            matches = rdLen - errores;
                                            matches = (pMem[pXONS[prevIndex].iGene].counts[x] = matches);
                                            if (matches >= bestCount) { 
                                                iGene = pXONS[prevIndex].iGene;
                                                if (matches > bestCount) { 
                                                    bestCount = matches;
                                                    bestGene = iGene;
                                                    bestTie = 0;
                                                } else if (iGene != bestGene) {
                                                    bestTieGene = iGene;
                                                    bestTie++;                                        
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (z == 0) {
                if (x == 0) {
                    bestFwdCount = bestCount;
                    bestFwdGene = bestGene;
                    fwdKeys = nValidKeys;
                    bestFwdTie = bestTie;
                    bestTieGeneFwd = bestTieGene;
                } else {
                    bestRevCount = bestCount;
                    bestRevGene = bestGene;
                    revKeys = nValidKeys;
                    bestRevTie = bestTie;
                    bestTieGeneRev = bestTieGene;
                }
            }
        }
    }

    if (bestFwdCount > bestRevCount) {
        // Fwd
        //z = bestFwdCount * 100 / rdLen; //bestFwdCount * 100 / fwdKeys;
        if (bestFwdTie == 1 && pRdMap->pMapParams->BREAK_TIE_COUNT1) {
            // La mayoría de los ties son con solo 1 gen, entonces se puede romper fácil tomandolo aleatorio, esto aumenta 3% de los reads (recupera el 15% de los reads)
            if (pRdMap->read->readNumber & 0x1) bestFwdGene = bestTieGeneFwd;
            bestFwdTie = 0;
        }
        if (bestFwdTie == 0 && bestFwdCount >= pRdMap->pMapParams->MIN_KEY_COUNTS && bestFwdCount >= pRdMap->pMapParams->MIN_MATCHES) {
            pCand->addGenomeAndReadPosition(bestFwdGene, 0, 0, z);
        }
        //if (bestFwdTie > 0) {
        //    fprintf(stderr, "[%u:%hu:%u]", bestFwdCount, bestFwdTie, maxTargets);
        //}
    } else if (bestRevCount > bestFwdCount) {
        // Reverse
        //z = bestRevCount * 100 / rdLen; //bestRevCount * 100 / revKeys;
        if (bestRevTie == 1 && pRdMap->pMapParams->BREAK_TIE_COUNT1) {
            // La mayoría de los ties son con solo 1 gen, entonces se puede romper fácil tomandolo aleatorio, esto aumenta 3% de los reads (recupera el 15% de los reads)
            if (pRdMap->read->readNumber & 0x1) bestRevGene = bestTieGeneRev;
            bestRevTie = 0;
        }
        if (bestRevTie == 0 && bestRevCount >= pRdMap->pMapParams->MIN_KEY_COUNTS && bestRevCount >= pRdMap->pMapParams->MIN_MATCHES) {
            pCand->addGenomeAndReadPosition(bestRevGene, 0, 1, z);
        }
        //if (bestRevTie > 0) {
        //    fprintf(stderr, "[%u:%hu:%u]", bestRevCount, bestRevTie, maxTargets);
        //}
    }
    if (pCand->getCount() == 0) {
        pChr = pGenome->getChromosome(bestFwdGene);
        fprintf(stderr, "Read=%llu, L=%u, RC=%c, %u=%s, C=%u/%u:%u K %s  %s\n", pMap->read->readNumber, pMap->read->lenSeq, '0', pChr->index, pChr->name, bestFwdCount, fwdKeys, nTotalKeys, pMap->read->pSeq, pMap->pLinkedRKM->read->pSeqRevComp); fflush(stderr);
        pChr = pGenome->getChromosome(bestRevGene);
        fprintf(stderr, "Read=%llu, L=%u, RC=%c, %u=%s, C=%u/%u:%u K %s  %s\n", pMap->read->readNumber, pMap->read->lenSeq, '1', pChr->index, pChr->name, bestRevCount, revKeys, nTotalKeys, pMap->read->pSeqRevComp, pMap->pLinkedRKM->read->pSeq); fflush(stderr);
    }
    //pCand->printWithGenomePositions(pMap->pGenome);
    return pCand->getCount();
}

uint16_t CountGeneReadBlastSmall(ogReadKeyMapping *pMap) {
    return countGeneKeysBlasting(pMap,  pMap->getLowKeyCountLimit());
}

uint16_t CountGeneReadBlastLarge(ogReadKeyMapping *pMap) {
    return countGeneKeysBlasting(pMap,  pMap->getHighKeyCountLimit());
}

uint16_t CountGeneReadBlastHuge(ogReadKeyMapping *pMap) {
    return countGeneKeysBlasting(pMap, 10*pMap->getHighKeyCountLimit());
}






#define SWEEP_MAX_GEN_POS   100


uint16_t SweepAllGP(ogReadKeyMapping *pMap, uint16_t maxTargets) {
    
    pMap->pCandPosMan->reset();
    
    uint32_t d, n;
    uint16_t i, j, minPosKey, anteMinPosKey, nKeys;
    uint32_t minPos, anteMinPos;
    uint32_t hitLen = pMap->read->lenSeq << 1;
    char FR, searching;
    
    // Active GP
    uint16_t    nGenPos, minKeys;
    uint32_t    *genPos[SWEEP_MAX_GEN_POS];
    uint32_t    nActiveGP[SWEEP_MAX_GEN_POS];
    uint16_t    readPos[SWEEP_MAX_GEN_POS];
    
    
    for (FR = 0; FR < 2; FR++) {
        // Reset internal structures
        nGenPos = 0;
        
        // Add all keys
        n = pMap->getNKeyPos(FR);
        for (i=0; nGenPos < SWEEP_MAX_GEN_POS && i < n; i++) {
            if ((d=pMap->getTargetSize(FR, i)) <= maxTargets) {
                // Add Key
                genPos[nGenPos] = pMap->getGenomePositionFromReadPositionKey(FR, i, 0);
                nActiveGP[nGenPos] = d;
                readPos[nGenPos] = pMap->getKeyablePosition(FR, i);
                nGenPos++;
            }
        }
//        BELOW CODE is an alternative by order size and using top 10
//        pMap->resetByOrderSize(FR);
//        n = pMap->getNKeyPos(FR);
//        for (j=0; j < 10 && nGenPos < SWEEP_MAX_GEN_POS && !pMap->finishByOrderSize(FR); j++) {
//            i = pMap->getIndexByOrderSize(FR);
//            if ((d=pMap->getTargetSize(FR, i)) <= maxTargets) {
//                // Add Key
//                genPos[nGenPos] = pMap->getGenomePositionFromReadPositionKey(FR, i, 0);
//                nActiveGP[nGenPos] = d;
//                readPos[nGenPos] = pMap->getKeyablePosition(FR, i);
//                nGenPos++;
//            }
//            pMap->nextByOrderSize(FR);
//        }
        // Search intersections
        searching = 1;
        minKeys = 2;
        while (searching) {
            searching = 0;
            // Buscar menor "actual" entre las posiciones genomicas
            // 
            minPos = 0xFFFFFFFF;
            minPosKey = 0xFFFF;
            for (i=0; i < nGenPos; i++) {
                if (nActiveGP[i]) {
                    if (*(genPos[i]) < minPos) {
                        anteMinPosKey = minPosKey;
                        anteMinPos = minPos;
                        minPosKey = i;
                        minPos = *(genPos[i]);
                        searching = 1;
                    } else if (*(genPos[i]) < anteMinPos) {
                        anteMinPos = *(genPos[i]);
                        anteMinPosKey = i;
                    }
                }
            }
            if (searching) {
                // Buscar entre los menores actuales cuántos matchean siempre y cuando el anteMin esta en rango, de otra forma es "aislado"
                if (anteMinPos - minPos <= hitLen) {
                    // Al menos hay 2 llaves con hitLen buscar mas
                    for (nKeys=0, i=0; i < nGenPos; i++) if (*(genPos[i]) - minPos <= hitLen && nActiveGP[i]) nKeys++;
                    if (nKeys >= minKeys) {
                        if (nKeys > minKeys) {
                            if (FR) {
                                pMap->pCandPosMan->restoreFwdCandidatePositions();
                            } else {
                                pMap->pCandPosMan->reset();
                            }
                            minKeys = nKeys;
                        }
                        pMap->pCandPosMan->addGenomeAndReadPosition(minPos, readPos[minPosKey], FR, nKeys * 100 / nGenPos);
                    }
                    genPos[minPosKey]++;
                    nActiveGP[minPosKey]--;
                } else {
                    genPos[minPosKey]++;
                    nActiveGP[minPosKey]--;
                    while (*(genPos[minPosKey]) < anteMinPos && nActiveGP[minPosKey]) {
                        if (anteMinPos - *(genPos[minPosKey]) <= hitLen) {
                            break;
                        } else {
                            genPos[minPosKey]++;
                            nActiveGP[minPosKey]--;                            
                        }
                    }
                }
            }
        }
        if (FR == 0) pMap->pCandPosMan->saveFwdCandidatePositions();
    }
       
    //pMap->pCandPosMan->saveFwdCandidatePositions();//pMap->pCandPosMan->setIsReverseComplement(revCount > fwdCount || (revCount == fwdCount && nKR > nKF)? 1 : 0);
    //if (fwdCount + revCount > 100) {
        //fprintf(stderr, "%llu, %u targets, %u fwd, %u rev:\n%s\n", iProc, fwdCount+revCount, fwdCount, revCount, pMap->seq);
    //}
    return pMap->pCandPosMan->getCount();
}

uint16_t _with_ante_prev_SweepAllGP(ogReadKeyMapping *pMap, uint16_t maxTargets) {
    
    pMap->pCandPosMan->reset();
    
    uint32_t d, n;
    uint16_t i, j, minPosKey, nPrevKeys, prevPosKey, anteMinPosKey;
    uint32_t prevPos, minPos, anteMinPos;
    uint32_t hitLen = pMap->read->lenSeq << 1;
    char FR, searching;
    
    // Active GP
    uint16_t    nGenPos, minKeys;
    uint32_t    *genPos[SWEEP_MAX_GEN_POS];
    uint32_t    nActiveGP[SWEEP_MAX_GEN_POS];
    uint16_t    readPos[SWEEP_MAX_GEN_POS];
    
    
    for (FR = 0; FR < 2; FR++) {
        // Reset internal structures
        nGenPos = 0;
        
        // Add all keys
        n = pMap->getNKeyPos(FR);
        for (i=0; nGenPos < SWEEP_MAX_GEN_POS && i < n; i++) {
            if ((d=pMap->getTargetSize(FR, i)) <= maxTargets) {
                // Add Key
                genPos[nGenPos] = pMap->getGenomePositionFromReadPositionKey(FR, i, 0);
                nActiveGP[nGenPos] = d;
                readPos[nGenPos] = pMap->getKeyablePosition(FR, i);
                nGenPos++;
            }
        }
        
//        pMap->resetByOrderSize(FR);
//        n = pMap->getNKeyPos(FR);
//        for (j=0; j < 10 && nGenPos < SWEEP_MAX_GEN_POS && !pMap->finishByOrderSize(FR); j++) {
//            i = pMap->getIndexByOrderSize(FR);
//            if ((d=pMap->getTargetSize(FR, i)) <= maxTargets) {
//                // Add Key
//                genPos[nGenPos] = pMap->getGenomePositionFromReadPositionKey(FR, i, 0);
//                nActiveGP[nGenPos] = d;
//                readPos[nGenPos] = pMap->getKeyablePosition(FR, i);
//                nGenPos++;
//            }
//            pMap->nextByOrderSize(FR);
//        }
        
        
        // Search intersections
        prevPosKey = SWEEP_MAX_GEN_POS;
        searching = 1;
        prevPos = 0;
        nPrevKeys = 0;
        minKeys = 2;
        while (searching) {
            searching = 0;
            
            // Buscar menor "actual" entre las posiciones genomicas
            // 
            anteMinPos = minPos = 0xFFFFFFFF;
            anteMinPosKey = minPosKey = 0xFFFF;
            for (i=0; i < nGenPos; i++) {
                if (nActiveGP[i]) {
                    if (*(genPos[i]) < minPos) {
                        anteMinPosKey = minPosKey;
                        anteMinPos = minPos;
                        minPosKey = i;
                        minPos = *(genPos[i]);
                        searching = 1;
                    } else if (*(genPos[i]) < anteMinPos) {
                        anteMinPos = *(genPos[i]);
                        anteMinPosKey = i;
                    }
                }
            }
            if (minPosKey < 0xFFFF) {
                genPos[minPosKey]++;
                nActiveGP[minPosKey]--;
                if (minPos - prevPos <= hitLen) {
                    // Esta dentro del rango del read de la posición anterior
                    if (minPosKey != prevPosKey) nPrevKeys++;
                } else {
                    // Registrar
                    if (nPrevKeys >= minKeys) {
                        if (nPrevKeys > minKeys) {
                            if (FR) {
                                pMap->pCandPosMan->restoreFwdCandidatePositions();
                            } else {
                                pMap->pCandPosMan->reset();
                            }
                            minKeys = nPrevKeys;
                        }
                        pMap->pCandPosMan->addGenomeAndReadPosition(prevPos, readPos[prevPosKey], FR, nPrevKeys * 100 / nGenPos);
                    }
                    
                    if (anteMinPos - minPos <= hitLen) {
                        // means that next minPos is also a candidate for overlap
                        //Actualizar
                        prevPos = minPos;
                        prevPosKey = minPosKey;
                        nPrevKeys = 1;
                    } else {
                        // Means that next minPos IS NOT A CANDIDATE for overlap
                        // Increase until anteMinPos
                        anteMinPos -= hitLen; // this is the possible minimum target
                        while (nActiveGP[minPosKey] && *(genPos[minPosKey]) < anteMinPos) {
                            genPos[minPosKey]++;
                            nActiveGP[minPosKey]--;
                        }
                        if (*(genPos[minPosKey]) < anteMinPos+hitLen && nActiveGP[minPosKey]) {
                            prevPos = *(genPos[minPosKey]);
                            prevPosKey = minPosKey;
                        } else if (anteMinPosKey < 0xFFFF) {
                            prevPos = anteMinPos+hitLen;
                            prevPosKey = anteMinPosKey;
                            genPos[anteMinPosKey]++;
                            nActiveGP[anteMinPosKey]--;
                        } else break;
                        nPrevKeys = 1;
                    }
                }
            } else break;
        }
        if (nPrevKeys >= minKeys) {
            pMap->pCandPosMan->addGenomeAndReadPosition(prevPos, readPos[prevPosKey], FR, nPrevKeys * 100 / nGenPos);            
        }
        if (FR == 0) pMap->pCandPosMan->saveFwdCandidatePositions();
    }
   
    pMap->pCandPosMan->saveFwdCandidatePositions();
    
    //pMap->pCandPosMan->saveFwdCandidatePositions();//pMap->pCandPosMan->setIsReverseComplement(revCount > fwdCount || (revCount == fwdCount && nKR > nKF)? 1 : 0);
    //if (fwdCount + revCount > 100) {
        //fprintf(stderr, "%llu, %u targets, %u fwd, %u rev:\n%s\n", iProc, fwdCount+revCount, fwdCount, revCount, pMap->seq);
    //}
    return pMap->pCandPosMan->getCount();
}

uint16_t SweepAllPosSmall(ogReadKeyMapping *pMap) {
    return SweepAllGP(pMap,  1 + (pMap->getLowKeyCountLimit() >>  2));
}

uint16_t SweepAllPosLarge(ogReadKeyMapping *pMap) {
    return SweepAllGP(pMap,  pMap->getLowKeyCountLimit());
}

uint16_t SweepAllPosHuge(ogReadKeyMapping *pMap) {
    return SweepAllGP(pMap,  pMap->getHighKeyCountLimit());
}
