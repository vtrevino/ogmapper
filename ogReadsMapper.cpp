/* 
 * File:   ogReadsMapper.cpp
 * Author: victortrevino
 * 
 * Created on May 20, 2022, 1:25 PM
 */
#include <stdio.h>
//#include <time.h>
#include <chrono>
#define clock_t     std::chrono::steady_clock::time_point
#define clock()     std::chrono::steady_clock::now()
#define duration_microsec(end, start)  std::chrono::duration_cast<chrono::microseconds>(end - start).count()
#define duration_nanosec(end, start)  std::chrono::duration_cast<chrono::nanoseconds>(end - start).count()
#include "ogDefinitions.hpp"
#include "ogReadsMapper.hpp"
#include "ogReadKeyMapping.hpp"
#include "CircularArray.hpp"
#include "CircularArray.cpp"
#include "ogSingleRead.hpp"

using namespace std;

ogReadsMapper::ogReadsMapper(ogMappingParameters *pMapParameters, ogScheduler *pScheduling, ogAligner *pAlign, ogSamWriter *pSamWri, CountingInfo *pCountInfo, uint32_t maxReadLen, uint16_t threadNum, float meanLog10, float sdLog10, uint32_t q95, char prodMode, char isPairedReads) {
    pMapParams = pMapParameters;
    pAligner = pAlign;
    nMatches = 0;
    nReads = 0;
    nMaps = 0;
    nKeys = 0;
    nFoundFwd = 0;
    nFoundRev = 0;
    elapsed = 0;
    preparation = 0;
    pairing = 0;
    nFound = 0;
    nReadsMappedR1R2 = 0;
    nReadsMappedR2R1 = 0;
    nReadsMappedR1 = 0;
    nReadsMappedR2 = 0;
    nReadsMappedTrans = 0;
    nReadsMappedOther = 0;
    nReadsUnmapped = 0;
    nReadsAlternMaps = 0;
    productionMode = prodMode;
    pClonedGuider = pMapParams->pGuiderDontUseIt->clone();
    isPaired = isPairedReads;
    pRdKeyMap[0] = new ogReadKeyMapping(maxReadLen < MAX_DEF_READ_LEN ? MAX_DEF_READ_LEN : maxReadLen, threadNum);
    pRdKeyMap[0]->setMasterObjects(pMapParams->pKeys,pMapParams->pGenome,pMapParams->pKeyEncoding,pMapParams->pGenPos, pClonedGuider, pMapParams->forceUpperCase, pSamWri, pMapParams);
    pRdKeyMap[0]->setStatistics(meanLog10, sdLog10, q95);
    pRdKeyMap[0]->pCountInfo = pCountInfo;
    if (isPaired) {
        pRdKeyMap[1] = new ogReadKeyMapping(maxReadLen < MAX_DEF_READ_LEN ? MAX_DEF_READ_LEN : maxReadLen, threadNum);
        pRdKeyMap[1]->setMasterObjects(pMapParams->pKeys,pMapParams->pGenome,pMapParams->pKeyEncoding,pMapParams->pGenPos, pClonedGuider, pMapParams->forceUpperCase, pSamWri, pMapParams);
        pRdKeyMap[1]->setStatistics(meanLog10, sdLog10, q95);
        pRdKeyMap[1]->pCountInfo = pCountInfo;
        pRdKeyMap[1]->pLinkedRKM = pRdKeyMap[0];
        pRdKeyMap[0]->pLinkedRKM = pRdKeyMap[1];
        pRdKeyMap[0]->isMemoryReleasable = true;
        pRdKeyMap[1]->isMemoryReleasable = false;
        pRdKeyMap[1]->isRead2 = 1;
    } else {
        pRdKeyMap[1] = NULL;
    }
    pRdKeyMapRd1 = pRdKeyMap[0];
    pRdKeyMapPair = pRdKeyMap[1];
    pScheduler = pScheduling;
    nFunc = pScheduler->nScheduleFunc;
    pCounts = (ogReadMapperCounts *) calloc(nFunc, sizeof(ogReadMapperCounts));
//    tCalls = (uint64_t *) calloc(nFunc, sizeof(uint64_t));
//    tMatches = (uint64_t *)calloc(nFunc, sizeof(uint64_t));
//    tFoundFwd = (uint64_t *)calloc(nFunc, sizeof(uint64_t));
//    tFoundRev = (uint64_t *)calloc(nFunc, sizeof(uint64_t));
//    tElapsed = (uint64_t *)calloc(nFunc, sizeof(uint64_t));
    busy = 0;
//    canWork = 0;
    finish = 0;
    wastingCycles = 0;
    wastingCycles2 = 0;
    wastingEvents = 0;
    maxQsize = 10;
    thread = threadNum;
    //qSeq = CircularArray<sequenceRead *>(maxQsize, 10);
    checkMatch = (pRdKeyMapRd1->pKeyEncoding->isKeyMatchExactSequence() == 0);
    compareKey = (pRdKeyMapRd1->pKeyEncoding->canKeyBeComparedInSequence() != 0);
    dnaAlloc = 0;
    DNA = NULL;
    allocDNA(BUFFER_FOR_UNPACKING_DNA);
    processingFunc = &ogReadsMapper::mapRead;
}

ogReadsMapper::~ogReadsMapper() {
    //fprintf(stderr, "<Deallocating ogReadsMapper %hu:", thread);
    //fprintf(stderr, "<~ogRM %hu:",thread); fflush(stderr);
//    free(tCalls);
//    free(tMatches);
//    free(tFoundFwd);
//    free(tFoundRev);
//    free(tElapsed);
    if (DNA != NULL) free(DNA);
    if (pCounts != NULL) free(pCounts);
    delete pRdKeyMapRd1;
    if (pRdKeyMapPair != NULL) delete pRdKeyMapPair;
    delete pClonedGuider;
    //fprintf(stderr, ":ogRM %hu>", thread);
}

void ogReadsMapper::allocDNA(uint32_t len) {
    if (len > dnaAlloc) {
        dnaAlloc    = len << 1;
        DNA         = (char *) realloc(DNA, sizeof(char) * dnaAlloc);
    }    
}

void ogReadsMapper::setMode(char mode) {
    if (mode == 'M') { 
        // mapping
        processingFunc = &ogReadsMapper::mapRead;
    } else if (mode == 'C') {
        // count
        processingFunc = &ogReadsMapper::countRead;
    }
}

void ogReadsMapper::countRead(ogSingleRead *theRead) {
    char aCigar[] = "M";
    clock_t start_t = clock();
    //clock_t curr_t;
    pRdKeyMapRd1->setRead(theRead); // Calculates Keys
    if (theRead->isPaired) {
        pRdKeyMapPair->setRead(theRead->pairedRead);
    }
    preparation += duration_microsec(clock(), start_t);

    uint32_t k1 = countTheReadWithKeyMap(pRdKeyMapRd1);
    uint32_t iGene;
    //uint32_t k2, k3, k4, gPos;
    nReads++;
    nKeys += pRdKeyMapRd1->nFwdPos;

    ogChromosome                *pChr;
    ogCandidatePosManager       *pCand = pRdKeyMapRd1->pCandPosMan;
    ogGenomeAndReadPosition     *pGenRdPos;
    
    //fprintf(stderr, "CCCC Counting k1=%u\n", k1); fflush(stderr);
    if (k1 == 1) {
        // just 1 hit, count!
        pGenRdPos = pCand->getkPos(0);
        // pGenRdPos->genomePosition contains the Gene index designated
        
        //pChr = pRdKeyMapRd1->pGenome->getGenomicCoordinate(pGenRdPos->genomePosition);
        //pMapParams->pGeneCounts[pChr->number-1]++;
                
        pMapParams->pGeneCounts[iGene = pGenRdPos->genomePosition]++;
        
        if (pMapParams->countExons) {
            ogReadKeyMapping    *pMap;
            char        x, y, w, RCkeys, f;
            uint32_t    nKeys, i, j, antExon, kGeneExon = 0;
            ogKey       *pKey;
            uint32_t    *pGenPos, *pLimGenPos, *pGeneExon = pMapParams->pGeneExon, *pX;
            Exon        *pXONS = pRdKeyMapRd1->pCountInfo->pExons;

            w = (pRdKeyMapRd1->pLinkedRKM == NULL ? 1 : 2);

            for (x = 0; x < 2; x++) {
                RCkeys = (pGenRdPos->isReverse ? 1-x : x);
                pMap = (x ? pRdKeyMapRd1->pLinkedRKM : pRdKeyMapRd1);
                nKeys = (RCkeys == 0 ? pMap->getNFwdPos() : pMap->getNRevPos());
                if (nKeys > 0) {
                    for (i=0; i < nKeys; i++) {
                        pKey = pMap->getInfoForPositionKey(RCkeys, i);
                        pGenPos = pMap->getGenomePositionFromReadPositionKey(RCkeys, i, 0);
                        antExon = 0xFFFFFFFF;
                        for (pLimGenPos = pGenPos + pKey->size; pGenPos < pLimGenPos && pXONS[*pGenPos].iGene != iGene; pGenPos++);
                        for (; pGenPos < pLimGenPos && pXONS[*pGenPos].iGene == iGene; pGenPos++) {
                            if (*pGenPos != antExon) {
                                antExon = *pGenPos; // pXONS[*pGenPos].iExon
                                for (f=0,j=0,pX=pGeneExon; j < kGeneExon; j++, pX++) {
                                    if (*pX == antExon) {
                                        f = 1;
                                        break;
                                    }
                                }
                                // Add to the queue if not found
                                if (f == 0) pGeneExon[kGeneExon++] = antExon;
                            }
                        }
                    }
                }
            }
            for (j=0, pX=pGeneExon; j < kGeneExon; j++, pX++) {
                pMapParams->pExonCounts[*pX]++;
                //if (iGene == 24419-1) {
                //    fprintf(stderr, "Rd=%llu, Exon=%i\n", theRead->readNumber, *pX); fflush(stderr);
                //}
            }
        }

        //fprintf(stderr, "Counting Gene [%u] [%s] [%s] = [%u]\n", pChr->number, pChr->name, pChr->comment, pMapParams->pGeneCounts[pChr->number-1]); fflush(stderr);        
        // Contar exones
        
        if (pMapParams->outputSAMcount) {
            pChr = pRdKeyMapRd1->pGenome->getChromosome(pGenRdPos->genomePosition);
            pRdKeyMapPair->pSamWriter->writeSAMInfo(
                theRead->pId,
                (theRead->isPaired ? SAMFLAG_READ_PAIRED : 0) | SAMFLAG_READS_MAPPED | SAMFLAG_READ1 | (pGenRdPos->isReverse ? SAMFLAG_READ_REVCOMP : 0) | (pGenRdPos->isReverse ? SAMFLAG_MATE_REVCOMP : 0),
                pChr->name,
                0,
                pGenRdPos->score,
                aCigar,
                pChr->name,
                0,
                0,
                (pGenRdPos->isReverse ? theRead->pSeqRevComp : theRead->pSeq),
                (pGenRdPos->isReverse ? theRead->pQualRev    : theRead->pQual),
                NULL,
                theRead->readNumber,
                pGenRdPos->function,
                    theRead
            );
            if (theRead->isPaired) {
                pRdKeyMapPair->pSamWriter->writeSAMInfo(
                    theRead->pairedRead->pId,
                    (theRead->pairedRead->isPaired ? SAMFLAG_READ_PAIRED : 0) | SAMFLAG_READS_MAPPED | SAMFLAG_READ1 | (pGenRdPos->isReverse ? SAMFLAG_READ_REVCOMP : 0) | (pGenRdPos->isReverse ? SAMFLAG_MATE_REVCOMP : 0),
                    pChr->name,
                    0,
                    pGenRdPos->score,
                    aCigar,
                    pChr->name,
                    0,
                    0,
                    (pGenRdPos->isReverse ? theRead->pairedRead->pSeq  : theRead->pairedRead->pSeqRevComp),
                    (pGenRdPos->isReverse ? theRead->pairedRead->pQual : theRead->pairedRead->pQualRev),
                    NULL,
                    theRead->pairedRead->readNumber,
                    pGenRdPos->function,
                        theRead
                );
            }
        }
    } else { // if (k1 == 0)
        nReadsUnmapped++;
    }
    
    elapsed += duration_microsec(clock(), start_t);
    
}


void ogReadsMapper::mapRead(ogSingleRead *theRead) {
    clock_t start_t = clock();
    clock_t curr_t;
    pRdKeyMapRd1->setRead(theRead); // Calculates Keys
    if (theRead->isPaired) {
        pRdKeyMapPair->setRead(theRead->pairedRead);
    }
    curr_t = clock();
    preparation += duration_microsec(curr_t, start_t);
        
    // Map Read 1st Read
    //fprintf(stderr, "<--->"); fflush(stderr);

    // First check no N's
    char read1_OnlyNs = theRead->lenSeq - theRead->Ncount < pRdKeyMapRd1->pKeyEncoding->getSizeInChars();
    char read2_OnlyNs = 0;

    nReads++;
    
    uint32_t k1 = (read1_OnlyNs ? 0 : mapTheReadWithKeyMap(pRdKeyMapRd1));
    uint32_t k2, k3, k4, gPos;
    ogGenomeAndReadPosition *pGenRdPos;
    ogChromosome            *pChr;

    nKeys += pRdKeyMapRd1->nFwdPos;
    
    //fprintf(stderr, "[---]"); fflush(stderr);
    if (theRead->isPaired) {
        read2_OnlyNs = theRead->pairedRead->lenSeq - theRead->pairedRead->Ncount < pRdKeyMapRd1->pKeyEncoding->getSizeInChars();
        if (k1 > 0) {
            // Check intersection with Rd2
            k2 = candidatesIntersectPair(pRdKeyMapRd1, pRdKeyMapPair, 0);
    //fprintf(stderr, "-RN=%llu,",pRdKeyMapRd1->read->readNumber); fflush(stderr);
            if (k2 > 0) {
                // Both Reads MAPPED, k2 - 1 : which is the best candidate position
                // written by candidates Intersect by default
                nReadsMappedR1R2++;
            } else {
                // Print UNPAIRED (call to candidatesIntersectPair is for printing)
                // fprintf(stderr, "Unpaired (R2 does not intersect) ReadNum=%llu\nRd1=%s\nRd2=%s\n", pRdKeyMapRd1->read->readNumber, pRdKeyMapRd1->read->pSeq, pRdKeyMapPair->read->pSeq);
                // Now map first Rd2
                k3 = (read2_OnlyNs ? 0 : mapTheReadWithKeyMap(pRdKeyMapPair));
                if (k3 > 0) {
                    // Check intersection with Rd1
                    k4 = candidatesIntersectPair(pRdKeyMapPair, pRdKeyMapRd1, 0);
                    if (k4 > 0) {
                        // Both Reads MAPPED, Fine, but mapped in the other way, k4 - 1 : which is the best candidate position
                        // fprintf(stderr, "----------------------------- Rd1<>Rd2, Rd2:Rd1 ------------------------------------- ReadNum=%llu\n", pRdKeyMapRd1->read->readNumber);
                        // written by candidate instersect by default
                        nReadsMappedR2R1++;
                    } else {
                        // UNMAPPED or translocation
                        // Seems Translocation since :
                        // - Rd1 MAPPED but does not intersect with Rd2 
                        // - And Rd2 MAPPED but does not intersect with Rd1
                        //fprintf(stderr, "Translocation (R1 & R2 Mapped but not intersect), ReadNum=%llu\nRd1=%s\nRd2=%s\n", pRdKeyMapRd1->read->readNumber, pRdKeyMapRd1->read->pSeq, pRdKeyMapPair->read->pSeq);
                        //fprintf(stderr, "============================= TRANSLOCATION Rd1,Rd2 =================================== ReadNum=%llu\n", pRdKeyMapRd1->read->readNumber);
                        //candidatesIntersectPair(pRdKeyMapRd1, pRdKeyMapPair, 1);
                        //candidatesIntersectPair(pRdKeyMapPair, pRdKeyMapRd1, 1);
                        writeSamFromReadKeyMapsTranslocated(pRdKeyMapRd1, pRdKeyMapPair, pRdKeyMapRd1->pSamWriter);
                        nReadsMappedTrans++; // nReadsMappedOthers??? aqui checar otros tipos
                    }
                } else {
                  // Rd1 Mapped, Rd2 UNMAPPED
                  //fprintf(stderr, "............................ Rd1:Mapped, Rd2:Unmapped ...................................... ReadNum=%llu\n", pRdKeyMapRd1->read->readNumber);
                  writeSamFromGenomeAndReadPositionsMapUnmap(pRdKeyMapRd1, pRdKeyMapPair, pRdKeyMapRd1->pCandPosMan->getkPos(pRdKeyMapRd1->maxScorePos), pRdKeyMapRd1->pSamWriter);
                  nReadsMappedR1++;
                }
            }
        } else { 
            // Not found by first searching Forward Read
            // Now map first Rd2
            k1 = (read2_OnlyNs ? 0 : mapTheReadWithKeyMap(pRdKeyMapPair));
            if (k1 > 0) {
                // Check intersection with Rd1
                k2 = candidatesIntersectPair(pRdKeyMapPair, pRdKeyMapRd1, 0);
                if (k2 > 0) {
                    // Both Reads MAPPED, Fine, but mapped in the other way, k2 - 1 : which is the best candidate position
                    // fprintf(stderr, "----------------------------- Rd2:Rd1 ----------------------------------------------- ReadNum=%llu\n", pRdKeyMapRd1->read->readNumber);
                    // Written by default by candidatesIntersect
                    nReadsMappedR2R1++;
                } else {
                    // Rd2 Mapped, Rd1 UNMAPPED (or translocation?)
                    // fprintf(stderr, "Unpaired (R1 does not intersect), ReadNum=%llu\nRd1=%s\nRd2=%s\n", pRdKeyMapRd1->read->readNumber, pRdKeyMapRd1->read->pSeq, pRdKeyMapPair->read->pSeq);
                    // fprintf(stderr, "============================= Rd2 mapped : Rd1 Unmapped ===================================== ReadNum=%llu\n", pRdKeyMapRd1->read->readNumber);
                    writeSamFromGenomeAndReadPositionsMapUnmap(pRdKeyMapPair, pRdKeyMapRd1, pRdKeyMapPair->pCandPosMan->getkPos(pRdKeyMapPair->maxScorePos), pRdKeyMapPair->pSamWriter);
                    nReadsMappedR2++;
                }
            } else {
                // UNMAPPED
                // fprintf(stderr, "UnMapped ReadNum=%llu\nRd1=%s\nRd2=%s\n", pRdKeyMapRd1->read->readNumber, pRdKeyMapRd1->read->pSeq, pRdKeyMapPair->read->pSeq);
                // fprintf(stderr, "********************************* Rd1 Unmapped & Rd2 Unmapped ********************************* ReadNum=%llu\n", pRdKeyMapRd1->read->readNumber);
                writeSamFromReadKeyMapsUnmapped(pRdKeyMapRd1, pRdKeyMapPair, pRdKeyMapPair->pSamWriter);
                nReadsUnmapped++;
            }
        }
    } else {
        if (k1 > 0) {
            writeSamFrom1GenomeAndReadPosition(pRdKeyMapRd1, pRdKeyMapRd1->pCandPosMan->getkPos(pRdKeyMapRd1->maxScorePos), pRdKeyMapRd1->pSamWriter);
            nReadsMappedR1++;
        } else {
            writeSamFrom1UnmappedRead(pRdKeyMapRd1, pRdKeyMapRd1->pSamWriter);
            nReadsUnmapped++;
        }
    }
    // fprintf(stderr, "\\RN=%llu,",pRdKeyMapRd1->read->readNumber); fflush(stderr);
}

uint16_t ogReadsMapper::candidatesIntersectPair(ogReadKeyMapping *pRdKeyMap1, ogReadKeyMapping *pRdKeyMap2, char print) {
    
    //FALTA MEJORAR EL PAIRING
    //Tal vez poniendo restricciones de cuando, si hay muchos candidatos mejor hacer el map del 2 read e intersección, o limitando la cantidad de opciones de pares solo incluyendo los ... 
    
    clock_t start_t = clock();
    clock_t curr_t;
    uint32_t    i, j;
    uint32_t    len1 = pRdKeyMap1->read->lenSeq, len2 = pRdKeyMap2->read->lenSeq;
    uint32_t    genPos, genPos2, genPos2raw, genDist, rdPos, genPosEnd, genPos2End, left, right, genPos2Ant = 0xFFFFFFFF;
    uint32_t    nCandidates, nKeys, *pInfoOffset, *pIO;
    uint32_t    isRev;
    ogGenomeAndReadPosition *pGenRdPos, *pGenRdPos2;
    ogGenomePositions *pRK2gp = pRdKeyMap2->pGenPos;
    ogKeyMap    *pKinfo;
    uint32_t    iBest = 0, jBest = 0, k, bestTies = 0, maxTies, maxI, maxRdPos = 0; // Max Refer to maximum for a candidate, best refer to a maximum of all candidates
    char        newMinError, score, iScore, maxScore = 0, best = 0, searchFlag, valid, localScore;
    uint32_t    genPosMinError, minError, error, minErrorTies, searchK, localK;
    ogCandidatePosManager *pCand = pRdKeyMap1->pCandPosMan;
    uint32_t    defineAlignMaxErr = MAX_ERROR_ACCEPT1_PACKED(len2);
    uint32_t    acceptAlignMaxErr = MAX_ERROR_CHECK3_ALIGN(len2);
    uint32_t    maxDif2GP = MAX_PAIR_INSERT_SIZE << 1;
    char        **pPacked;
    
    nCandidates = pCand->getCount();

    // General Score Para el Par = (len1+len2 - Error Rd 1 - Error Rd 2)/(len1+len2)
    
    //iBest = nCandidates;
    //if (nCandidates > 20) {
    //    fprintf(stderr, ">>>>> RN:%llu, Candidates:%u\n",pRdKeyMap1->read->readNumber,nCandidates); fflush(stderr);
    //}
    for (i=0; i < nCandidates; i++) {
        pGenRdPos = pCand->getkPos(i);
        if (pGenRdPos->status == 'A' || pGenRdPos->status == 'G') {
            //if (nCandidates > 20) fprintf(stderr, "%c",pGenRdPos->function); fflush(stderr);
            pGenRdPos->pairedKeys = 0;
            pGenRdPos->maxPairedIndex = 0;
            iScore = pGenRdPos->score;
            genPos = pGenRdPos->genomePosition - pGenRdPos->readPosition;
            genPosEnd = genPos + len1;
            isRev = 1 - pGenRdPos->isReverse; // Necesitamos el otro
            nKeys = pRdKeyMap2->getNKeyPos(isRev);
            genPos2Ant = 0xFFFFFFFF;
            pCand->initPairedGenomeAndReadPosition(i);
            minError = len2;
            minErrorTies = 0;
            maxScore = 0;
            maxTies = 0;
            maxI = 0;
            //pPacked = (isRev ? pRdKeyMap2->read->packedReadSeqRev : pRdKeyMap2->read->packedReadSeqFwd);
            
            for (j=0; j < nKeys; j++) {
                pKinfo = pRdKeyMap2->getKeyMapForPosition(isRev, j);
                if (pKinfo->keyInfo.size > 0) {
                    //fprintf(stderr, "(((( Key=%u / %u, Key.Size=%u\n",j,nKeys,pKinfo->keyInfo.size);  fflush(stderr);
                    valid = 0;
                    error = 99;
                    //minError = len2;
                    pInfoOffset = pRK2gp->getPointerPosition(pKinfo->keyInfo.offset);
                    searchK = pRdKeyMap2->getCloserGenomeRelativePosition(genPos, pInfoOffset, pKinfo->keyInfo.size);
                    //fprintf(stderr, "searchK in=%d\n",searchK);  fflush(stderr);
                    if (searchK > 0) {
                        pIO = pInfoOffset + searchK;
                        if ((*pIO - *(pIO - 1)) <= maxDif2GP) {
                            searchK--;
                        }
                    }
                    localScore = 0;
                    localK = 0;
                    for (; searchK < pKinfo->keyInfo.size; searchK++) {
                        //fprintf(stderr, "searchK=%d / %d\n",searchK, pKinfo->keyInfo.size); fflush(stderr);
                        genPos2raw = genPos2 = *(pInfoOffset + searchK);
                        if (genPos2 > pKinfo->position) genPos2 -= pKinfo->position;
                        //if (genPos2 != genPos2Ant) {
                            // To optimize calculations
                            genPos2End = genPos2 + len2;
                            left = genPos2 < genPos ? genPos2 : genPos;
                            right = genPos2End > genPosEnd ? genPos2End : genPosEnd;
                            genDist = (left > right ? left-right : right-left);
                            genPos2Ant = genPos2;
                        //}
                        if (genDist <= MAX_PAIR_INSERT_SIZE) {
                            // 2nd Read Seems to Match
                            error = pRdKeyMap2->pGenome->comparePackedSequences(pRdKeyMap2->read->packReadSeqIfNeeded(pRdKeyMap2->pGenome->getPackedOffset(genPos2), isRev), genPos2, len2, minError);
                            if (error <= minError) {
                                if (error < minError) {
                                    score = (len2 - error) * 100 / len2;
                                    minError = error;
                                    //pCand->resetPairedGenomeAndReadPosition(pGenRdPos); // This reset is wrong because if there are more than 1 candidate there could be problems of indexes
                                }
                                k = pCand->checkPairedGenomeAndReadPosition(pGenRdPos, genPos2, 0, isRev);
                                if (k == 0) {
                                    if (error < acceptAlignMaxErr) {
                                        k = pCand->addPairedGenomeAndReadPosition(pGenRdPos, genPos2raw, pKinfo->position, isRev, score); // genPos2raw, pKinfo->position
                                        valid = 1;
                                        if (score > localScore) {
                                            localScore = score;
                                            localK = k;
                                        }
                                    }
                                } else {
                                    score = (pCand->getkPairedPos(k-1))->score;
                                    error = len2 - (score * len2) / 100;
                                    valid = 1; // si usar para que cuente como pairedKey
                                    if (score > localScore) {
                                        localScore = score;
                                        localK = k;
                                    }
                                }
                            }
                            //fprintf(stderr, "err=%u , maxErr=%u ",error, acceptAlignMaxErr); fflush(stderr);
                            //fprintf(stderr, "sc=%hhu, err=%d, k=%u\n", score, error, k); fflush(stderr);
                        } else if (genPos2 > genPos) {
                            break;
                        }
                    }
                    //fprintf(stderr, "localScore=%hhu, RN=%llu\n",localScore,pRdKeyMap1->read->readNumber); fflush(stderr);
                    if (valid) {
                        score = ((iScore + localScore) >> 1) + ((++pGenRdPos->pairedKeys) * 10 / nKeys);
                        if (score >= maxScore) {
                            if (score > maxScore) {
                                maxScore = score;
                                maxI = localK-1;
                                maxTies = 0;
                                //maxRdPos= pKinfo->position;
                                //if (k==0) fprintf(stderr, "-------------------Pinche 0----------------------\n"); fflush(stderr);
                            } else {
                                maxTies++;
                            }
                            if (score >= best) {
                                if (score > best) {
                                    best = score;
                                    iBest = i + 1;
                                    bestTies = 0;
                                } else {
                                    bestTies++;
                                }
                            }
                        }
                        //if (error < defineAlignMaxErr) break;
                    }
                    //fprintf(stderr, "<<<<\n"); fflush(stderr);
                }
                
                /*****   
                    
                    searchFlag = 1;
                    valid = 1;
                    error = 99;
                    pInfoOffset = pRdKeyMap2->pGenPos->getPointerPosition(pKinfo->keyInfo.offset);
                    genPos2 = pRdKeyMap2->getCloserGenomePosition(genPos, pInfoOffset, pKinfo->keyInfo.size);
                    if (genPos2 > pKinfo->position) genPos2 -= pKinfo->position;
                    while (1) {
                        //fprintf(stderr, "[[[[[[[ RN=%llu, genPos=%u, genPos2=%u, searchK=%u\n", pRdKeyMap1->read->readNumber, genPos, genPos2, searchK);
                        if (genPos2 != genPos2Ant) {
                            // To optimize calculations
                            genPos2End = genPos2 + len2;
                            left = genPos2 < genPos ? genPos2 : genPos;
                            right = genPos2End > genPosEnd ? genPos2End : genPosEnd;
                            genDist = (left > right ? left-right : right-left);
                            genPos2Ant = genPos2;
                        }
                        if (genDist <= MAX_PAIR_INSERT_SIZE) {
                            // 2nd Read Seems to Match
                            k = pCand->checkPairedGenomeAndReadPosition(i, genPos2, 0, isRev);
                            if (k == 0) {
                                char **pPacked = (isRev ? pRdKeyMap2->read->packedReadSeqRev : pRdKeyMap2->read->packedReadSeqFwd);
                                error = pRdKeyMap2->pGenome->comparePackedSequences(pPacked, genPos2, len2, len2);
                                if (error <= minError) {
                                    if (error < minError) {
                                        genPosMinError = genPos2;
                                        minError = error;
                                        minErrorTies = 0;
                                    } else {
                                        minErrorTies++;
                                    }
                                }
                                score = (len2 - error) * 100 / len2;
                                if (score > 70) {
                                } else {
                                    //if (searchFlag == 0) fprintf(stderr, "Found at genPos2=%u, score=%hu\n", genPos2, score);
                                    k = pCand->addPairedGenomeAndReadPosition(i, genPos2, 0, isRev, score);
                                    if (searchFlag || score > 90) break;
                                }
                            } else {
                                score = (pCand->getkPairedPos(k-1))->score;
                                //fprintf(stderr, "Repeated. Bye.\n");
                                break;
                            }
                        } else {
                            // ya se salió de rango
                            if (genPos2 > genPos || searchFlag) {
                                valid = 0;
                                break;
                            }
                        }
                        if (++searchK >= pKinfo->keyInfo.size) {
                            valid = 0;
                            //fprintf(stderr, "bye ]]]]]]]\n");
                            break;
                        }
                        genPos2 = *(pInfoOffset + searchK);
                        if (genPos2 > pKinfo->position) genPos2 -= pKinfo->position;
                        //fprintf(stderr, "new genPos2=%u, searchK=%u ]]]]]]]\n", genPos2, searchK);
                    }
                    if (valid) {
                        score = ((iScore + score) >> 1) + ((++pGenRdPos->pairedKeys) * 10 / nKeys);
                        if (score >= maxScore) {
                            if (score > maxScore) {
                                maxScore = score;
                                maxI = k-1;
                                maxTies = 0;
                            } else {
                                maxTies++;
                            }
                            if (score >= best) {
                                if (score > best) {
                                    best = score;
                                    iBest = i + 1;
                                    bestTies = 0;
                                } else {
                                    bestTies++;
                                }
                            }
                        }
                    }
                    ****/
                if (print) {
                    ogChromosome *pChr = pRdKeyMap1->pGenome->getGenomicCoordinate(genPos);
                    //fprintf(stderr, "ReadNum=%llu | Cand=%u, RC=%hhu, St=%c, Sc=%hhu, Chr %u:%u, GenPos1 L=%u, R=%u | Key=%u, GenPos2 L=%u, R=%u | Dif=%u, Err=%u\n", seqRead->readNumber, i, pGenRdPos->isReverse, pGenRdPos->status, pGenRdPos->score, pChr->number, genPos - pChr->start, genPos, genPosEnd, j, genPos2, genPos2End, genDist, error);
                    fprintf(stderr, "ReadNum=%llu | Cand=%u, RC=%hhu, St=%c, Sc=%hhu, Chr %u:%u, GenPos1 L=%u, R=%u | Key=%u, Valid=%c, GenPos2 L=%u, R=%u | Dif=%u, Sc=%u, Err=%u\n", pRdKeyMap1->read->readNumber, i, pGenRdPos->isReverse, pGenRdPos->status, pGenRdPos->score, pChr->number, genPos - pChr->start, genPos, genPosEnd, j, valid+48, genPos2, genPos2End, genDist, score, error);
                }
            }
            if (maxScore > 0) {
                pGenRdPos->score = maxScore;
                pGenRdPos->status = (maxTies ? 'p' : 'P');
                pGenRdPos->maxPairedIndex = maxI;
                pGenRdPos2 = pCand->getkPairedPos(pGenRdPos->maxPairedIndex);
                //pGenRdPos2->readPosition = maxRdPos;
                //if (pGenRdPos->pairedCounts > 20) {
                //    fprintf(stderr, "RN:%llu Cand:%u/%u, Paired: Counts=%u, Max=%u, Score=%hhu\n", pRdKeyMap1->read->readNumber, i, nCandidates, pGenRdPos->pairedCounts, pGenRdPos->maxPairedIndex - pGenRdPos->pairedIndex, pGenRdPos->score); fflush(stderr);
                //}
            }
        } else if (print) {
            //fprintf(stderr, "[%c]",pGenRdPos->status);
            genPos = pGenRdPos->genomePosition - pGenRdPos->readPosition;
            ogChromosome *pChr = pRdKeyMap1->pGenome->getGenomicCoordinate(genPos);
            //fprintf(stderr, "ReadNum=%llu | Cand=%u, RC=%hhu, St=%c, Sc=%hhu, Chr %u:%u, GenPos1 (L=%u, R=%u)\n", seqRead->readNumber, i, pGenRdPos->isReverse, pGenRdPos->status, pGenRdPos->score, pChr->number, genPos - pChr->start, genPos, genPosEnd);
            fprintf(stderr, "ReadNum=%llu | Cand=%u, RC=%hhu, St=%c, Sc=%hhu, Chr %u:%u, GenPos1 (L=%u, R=%u)\n", pRdKeyMap1->read->readNumber, i, pGenRdPos->isReverse, pGenRdPos->status, pGenRdPos->score, pChr->number, genPos - pChr->start, genPos, genPosEnd);
        }
    }
    if (iBest > 0) {
        pGenRdPos = (pCand->getkPos(iBest-1));
        if (pGenRdPos->status == 'P' || pGenRdPos->status == 'p') {
            //
            //if (pRdKeyMap1->read->readNumber > 9000 || pGenRdPos->maxPairedIndex > 1000) {
            //    fprintf(stderr, "RN=%llu, mpi=%u\n",pRdKeyMap1->read->readNumber,pGenRdPos->maxPairedIndex);
            //}
            //fprintf(stderr, "p1=%p,p2=%p, mxPI=%u\n",pRdKeyMap1,pRdKeyMap2,pGenRdPos->maxPairedIndex); fflush(stderr);
            writeSamFromGenomeAndReadPositions(pRdKeyMap1, pRdKeyMap2, pGenRdPos, pCand->getkPairedPos(pGenRdPos->maxPairedIndex), pRdKeyMap1->pSamWriter);
            if (pGenRdPos->status == 'p') {
                for (i=0; i < pGenRdPos->pairedCounts; i++) {
                    if (pGenRdPos->maxPairedIndex != i+pGenRdPos->pairedIndex) {
                        ogGenomeAndReadPosition *grp2 =pCand->getkPairedPos(i + pGenRdPos->pairedIndex);
                        if (grp2->score == maxScore) {
                            // No es el mejor, escribir alternativo
                            // Hacer el write?
                            nReadsAlternMaps++;
                            //fprintf(stderr, "Writing Alternate %s\n", pRdKeyMap1->read->pId);
                            writeAltRd2SamFromGenomeAndReadPosition(pRdKeyMap1, pRdKeyMap2, pGenRdPos, grp2, pRdKeyMap1->pSamWriter);
                        }
                    }
                }
            }
        } else {
            iBest = 0;
        }
    }
    //fprintf(stderr, "Best=%u, iBest=%u\n", best, iBest);
    // if (bestTies == 0) // Es único
    // if (bestTies > 0) // se pueden checar todos los status P/p para decidir el mas corto, el mas largo, o todos
    curr_t = clock();
    pairing += duration_microsec(curr_t, start_t);
    return iBest;
}

/***
                                // To optimize, only update when necessary
                                seqRead->genPos1_mapScore = pGenRdPos->score;
                                seqRead->genPos2_mapScore = pGenRdPos->score; // Asumimos que es el mismo score?
                                seqRead->genPos1_isRd2 = isSecondRead2;
                                seqRead->posIsReverse1 = pGenRdPos->isReverse;
                                seqRead->posIsReverse2 = isRev;
                                seqRead->genPos1_Left = genPos;
                                seqRead->genPos1_Right = genPosEnd;
                                seqRead->genPos2_Left = genPos2;
                                seqRead->genPos2_Right = genPos2End;
                                seqRead->genPos_leftMax = left;
                                seqRead->genPos_leftMin = (genPos < genPos2 ? genPos : genPos2);
                                seqRead->genPos_rightMin = right;
                                seqRead->genPos_rightMax = (genPosEnd > genPos2End ? genPosEnd : genPos2End);
                                seqRead->genPos_difference = genDist;
                                seqRead->genPos_insertSize = seqRead->genPos_rightMax - seqRead->genPos_leftMin;
                                //if (seqRead->genPos_insertSize > 1000) {
                                //    fprintf(stderr, "ReadNum=%llu | Cand=%u, St=%c, GenPos1 L=%u, R=%u | Key=%u, GenPos2 L=%u, R=%u | Dif=%u | Insert Size=%u\n", seqRead->readNumber, i, pGenRdPos->status, genPos, genPosEnd, j, genPos2, genPos2End, genDist, seqRead->genPos_insertSize);
                                //}
**/


#define min(a,b)    (a < b ? a : b)
#define max(a,b)    (a > b ? a : b)


void ogReadsMapper::buildPlainCIGARorWFA(ogReadKeyMapping *pRKM1, ogSingleRead *pR1, uint32_t genomicPos, char isRev) {
    uint32_t        seqLen = pR1->lenSeq;
    uint32_t        maxErrCheck3 = MAX_ERROR_CHECK3_ALIGN(seqLen);
    //char            **pPacked = (isRev ? pR1->packedReadSeqRev : pR1->packedReadSeqFwd);
    uint32_t        errPack = pRKM1->pGenome->comparePackedSequences(pRKM1->read->packReadSeqIfNeeded(pRKM1->pGenome->getPackedOffset(genomicPos), isRev), genomicPos, seqLen, maxErrCheck3);
    uint32_t        maxErrAcceptPacked = MAX_ERROR_ACCEPT1_PACKED(seqLen);
    
    if (errPack <= maxErrAcceptPacked) {
        pAligner->writeOnlyMatchesCigar(pR1, pR1->cigarLeftPos = genomicPos, pR1->cigarRightPos = genomicPos+seqLen);
    } else {
        allocDNA(seqLen);
        pRKM1->pGenome->extractFromGenomicPos(genomicPos, seqLen, DNA, 0, 0);
        //pAligner->compareSeqAndDNABuildingCIGAR(isRev ? pR1->pSeqRevComp : pR1->pSeq, DNA, pRKM1, pR1);
        //if (pR1->cigarSoft > 8) {
            char *pCIGAR;
            int cigarLen;
            int free = 0; //pR1->lenSeq >> 5;
            pAligner->fullAligments++;
            // WFA2 (https://github.com/smarco/WFA2-lib)
            pAligner->pWFAligner->alignEndsFree(DNA, seqLen, free, free, isRev ? pR1->pSeqRevComp : pR1->pSeq, seqLen, free, free);
            pAligner->pWFAligner->getAlignmentCigar(&pCIGAR, &cigarLen);
            //pGenRdPos->status = 'G'; 
            //pGenRdPos->score = buildCIGARfromWFA(pCIGAR, read); //(seqLen-errPack)*100/seqLen;
            pAligner->buildCIGARfromWFA(pCIGAR, pR1, genomicPos, genomicPos+seqLen);
        //}
    }
}


void ogReadsMapper::writeSamFromGenomeAndReadPositions(ogReadKeyMapping *pRKM1, ogReadKeyMapping *pRKM2, ogGenomeAndReadPosition *grp1, ogGenomeAndReadPosition *grp2, ogSamWriter *pSamWri) {

    //fprintf(stderr, "w1 p1=%p,p2=%p,p3=%p,p5=%p\n",pRKM1,pRKM2,grp1,grp2); fflush(stderr);

    ogSingleRead    *pR1 = pRKM1->read;
    ogSingleRead    *pR2 = pRKM2->read;
    uint32_t        left_r1 = grp1->genomePosition - grp1->readPosition; // - (pR1->cigarLeftType == 'I' || pR1->cigarLeftType == 'D' ? pR1->cigarLeft : 0) + (pR1->cigarLeftType == 'S' ? pR1->cigarLeft : 0) + pR1->cigarIns - pR2->cigarDel;
    uint32_t        right_r1 = left_r1 + pR1->lenSeq;

    if (*pR1->pCigar == 0) {
        buildPlainCIGARorWFA(pRKM1, pR1, left_r1, grp1->isReverse);
        //left_r1 = grp1->genomePosition - grp1->readPosition; // - (pR1->cigarLeftType == 'I' || pR1->cigarLeftType == 'D' ? pR1->cigarLeft : 0) + (pR1->cigarLeftType == 'S' ? pR1->cigarLeft : 0) + pR1->cigarIns - pR2->cigarDel;
    }
    left_r1 = pR1->cigarLeftPos;
    right_r1 = pR1->cigarRightPos;
    //if (pR1->cigarMatches + pR1->cigarNoMatches != pR1->lenSeq || pR1->cigarLeftType != '-' || pR1->cigarIns > 0 || pR1->cigarDel > 0 || pR1->cigarSoft > 0) {
        // no es solo xxxM sino que hay S/D/I y por lo tanto la posición izq hay que corregirla
        //fprintf(stderr,"Correcting Read ID:%s, rd pos:%u, CIGAR:%s, left_r1:%u, right_r1:%u\n", pR1->pId, grp1->readPosition, pR1->pCigar, left_r1, right_r1);
    //    correctLeftRighPositionFromCIGAR(pR1->pCigar, grp1->readPosition, &left_r1, &right_r1, pR1->readIndex == 2);
        //fprintf(stderr,"Post left_r1:%u, right_r1:%u\n", left_r1, right_r1);
    //}
    
    uint32_t        left_r2 = grp2->genomePosition - grp2->readPosition; // - (pR2->cigarLeftType == 'I' || pR2->cigarLeftType == 'D' ? pR2->cigarLeft : 0) + (pR2->cigarLeftType == 'S' ? pR2->cigarLeft : 0) + pR2->cigarIns - pR2->cigarDel;
    uint32_t        right_r2 = left_r2 + pR2->lenSeq;
    if (*pR2->pCigar == 0) {
        buildPlainCIGARorWFA(pRKM2, pR2, left_r2, grp2->isReverse);
        //left_r2 = grp2->genomePosition - grp2->readPosition; // - (pR2->cigarLeftType == 'I' || pR2->cigarLeftType == 'D' ? pR2->cigarLeft : 0) + (pR2->cigarLeftType == 'S' ? pR2->cigarLeft : 0) + pR2->cigarIns - pR2->cigarDel;
    }
    left_r2 = pR2->cigarLeftPos;
    right_r2 = pR2->cigarRightPos;
    //if (pR2->cigarMatches + pR2->cigarNoMatches != pR2->lenSeq || pR2->cigarLeftType != '-' || pR2->cigarIns > 0 || pR2->cigarDel > 0 || pR2->cigarSoft > 0) {
        // no es solo xxxM sino que hay S/D/I y por lo tanto la posición izq hay que corregirla
        //fprintf(stderr,"Correcting Read ID:%s, rd pos:%u, CIGAR:%s, left_r2:%u, right_r2:%u\n", pR2->pId, grp2->readPosition, pR2->pCigar, left_r2, right_r2);
    //    correctLeftRighPositionFromCIGAR(pR2->pCigar, grp2->readPosition, &left_r2, &right_r2, pR2->readIndex == 2);
        //fprintf(stderr,"Post left_r2:%u, right_r2:%u\n", left_r2, right_r2);
    //}
    //fprintf(stderr, "wx pr1=%p,pr2=%p\n",pR1,pR2); fflush(stderr);
    //fprintf(stderr, "w2 p1=%p,p2=%p,p3=%p,p5=%p\n",pRKM1,pRKM2,grp1,grp2); fflush(stderr);
    //fprintf(stderr, "w3 p1=%p,p2=%p,p3=%p,p5=%p\n",pRKM1,pRKM2,grp1,grp2); fflush(stderr);
    uint32_t        left = left_r1 < left_r2 ? left_r1 : left_r2;
    uint32_t        right = right_r1 > right_r2 ? right_r1 : right_r2;
    
    uint32_t        genDist1 = (left_r1 > right_r1 ? left_r1 - right_r1 : right_r1 - left_r1);
    uint32_t        genDist2 = (left_r2 > right_r2 ? left_r2 - right_r2 : right_r2 - left_r2);
    //uint32_t        genDist = (left > right ? left-right : right-left);
    ogChromosome   *pChr1 = pRKM1->pGenome->getGenomicCoordinate(left);
    //char            possibleReverse2 = 0;
    
    uint32_t        L1 = left_r1 - pChr1->start + 1;
    uint32_t        L2 = left_r2 - pChr1->start + 1;
    
    /***
    if ((left_r1 <= left_r2 && right_r1 > right_r2)) { // && pR1->cigarLeft == 0 && pR1->cigarSoft == 0 && pR2->cigarLeft == 0 && pR2->cigarSoft == 0
        //fprintf(stderr, "*** Overlap 1: %s, left_r1:%u, left_r2:%u, right_r1:%u, right_r2:%u, CHR=%s:%u\n", pR1->pId, left_r1, left_r2, right_r1, right_r2, pChr1->name, L1);
        //fprintf(stderr, "Prev CIGAR: %s, ", pR1->pCigar);
        cutCIGARtoSeqLen(pR1->pCigar, right_r2-left_r2+1);
        //fprintf(stderr, "New CIGAR: %s\n", pR1->pCigar);
        genDist = right_r2-left_r2;
    }
    if ((left_r1 > left_r2 && right_r1 <= right_r2)) { // && pR1->cigarLeft == 0 && pR1->cigarSoft == 0 && pR2->cigarLeft == 0 && pR2->cigarSoft == 0
        //fprintf(stderr, "*** Overlap 2: %s, left_r1:%u, left_r2:%u, right_r1:%u, right_r2:%u, CHR=%s:%u\n", pR1->pId, left_r1, left_r2, right_r1, right_r2, pChr1->name, L1);
        //fprintf(stderr, "Prev CIGAR: %s, ", pR2->pCigar);
        cutCIGARtoSeqLen(pR2->pCigar, right_r1-left_r1+1);
        //fprintf(stderr, "New CIGAR: %s\n", pR2->pCigar);
        genDist = right_r1-left_r1;
    }
    ****/
    
    pSamWri->writeSAMInfo(
        pR1->pId,
        SAMFLAG_READ_PAIRED | SAMFLAG_READS_MAPPED | (pR1->readIndex == 1 ? SAMFLAG_READ1 : SAMFLAG_READ2) | (grp1->isReverse ? SAMFLAG_READ_REVCOMP : 0) | (grp2->isReverse ? SAMFLAG_MATE_REVCOMP : 0),
        pChr1->name,
        L1,
        grp1->score,
        pR1->pCigar, // (grp1->isReverse && (pR1->cigarDel > 0 || pR1->cigarIns > 0) ? revertCIGAR(pR1->pCigar, reverseCIGAR) : pR1->pCigar)
        (char *) EQUAL,
        L2,
        (left_r1 < left_r2 ? genDist1 : -genDist1),
        (grp1->isReverse ? pR1->pSeqRevComp : pR1->pSeq),
        (grp1->isReverse ? pR1->pQualRev    : pR1->pQual),
        NULL,
        pR1->readNumber,
        grp1->function,
            pR1
    );
    pSamWri->writeSAMInfo(
        pR2->pId,
        SAMFLAG_READ_PAIRED | SAMFLAG_READS_MAPPED | (pR1->readIndex == 1 ? SAMFLAG_READ2 : SAMFLAG_READ1)  | (grp2->isReverse ? SAMFLAG_READ_REVCOMP : 0) | (grp1->isReverse ? SAMFLAG_MATE_REVCOMP : 0),
        pChr1->name,
        L2,
        grp2->score,
        pR2->pCigar, // (possibleReverse2 && grp2->isReverse && (pR2->cigarDel > 0 || pR2->cigarIns > 0) ? revertCIGAR(pR2->pCigar, reverseCIGAR) : pR2->pCigar),
        (char *) EQUAL,
        L1,
        (left_r1 < left_r2 ? -genDist2 : genDist2),
        (grp2->isReverse ? pR2->pSeqRevComp : pR2->pSeq),
        (grp2->isReverse ? pR2->pQualRev    : pR2->pQual),
        NULL,
        pR2->readNumber,
        grp2->function,
            pR2
    );
}


void ogReadsMapper::writeSamFrom1GenomeAndReadPosition(ogReadKeyMapping *pRKM1, ogGenomeAndReadPosition *grp1, ogSamWriter *pSamWri) {
    ogSingleRead    *pR1 = pRKM1->read;
    uint32_t        left_r1 = grp1->genomePosition - grp1->readPosition;
    uint32_t        right_r1 = left_r1 + pR1->lenSeq;
    ogChromosome   *pChr1 = pRKM1->pGenome->getGenomicCoordinate(left_r1);
    
    if (*pR1->pCigar == 0) {
        //fprintf(stderr, "{ R1 %u:%s:%u: }\n", left_r1, pChr1->name, left_r1 - pChr1->start);
        buildPlainCIGARorWFA(pRKM1, pR1, left_r1, grp1->isReverse);
        //left_r1 = grp1->genomePosition - grp1->readPosition - (pR1->cigarLeftType == 'I' || pR1->cigarLeftType == 'D' ? pR1->cigarLeft : 0);
        left_r1 = pR1->cigarLeftPos;
        right_r1 = pR1->cigarRightPos;
        //allocDNA(pR1->lenSeq);
        //pRKM1->pGenome->extractFromGenomicPos(left_r1, pR1->lenSeq, DNA, 0, 0);
        //pAligner->compareSeqAndDNABuildingCIGAR(grp1->isReverse ? pR1->pSeqRevComp : pR1->pSeq, DNA, pRKM1, pR1);
    }
    pSamWri->writeSAMInfo(
        pR1->pId,
        SAMFLAG_READS_MAPPED | SAMFLAG_READ1 | (grp1->isReverse ? SAMFLAG_READ_REVCOMP : 0),
        pChr1->name,
        left_r1 - pChr1->start + 1, // + (pR1->cigarLeftType == 'S' ? pR1->cigarLeft : 0),
        grp1->score,
        pR1->pCigar,
        (char *) ASTERISK,
        0,
        0,
        (grp1->isReverse ? pR1->pSeqRevComp : pR1->pSeq),
        (grp1->isReverse ? pR1->pQualRev    : pR1->pQual),
        NULL,
        pR1->readNumber,
        grp1->function,
            pR1
    );
}


void ogReadsMapper::writeSamFrom1UnmappedRead(ogReadKeyMapping *pRKM1, ogSamWriter *pSamWri) {
    ogSingleRead    *pR1 = pRKM1->read;
    
    pSamWri->writeSAMInfo(
        pR1->pId,
        SAMFLAG_READ1 | SAMFLAG_READ_UNMAPPED,
        (char *) ASTERISK,
        0,
        0,
        (char *) ASTERISK,
        (char *) ASTERISK,
        0,
        0,
        pR1->pSeq,
        pR1->pQual,
        NULL,
        pR1->readNumber,
        '-',
            pR1
    );
}



void ogReadsMapper::writeSamFromGenomeAndReadPositionsMapUnmap(ogReadKeyMapping *pRKM1, ogReadKeyMapping *pRKM2, ogGenomeAndReadPosition *grp1, ogSamWriter *pSamWri) {
    ogSingleRead    *pR1 = pRKM1->read;
    ogSingleRead    *pR2 = pRKM2->read;
    uint32_t        left_r1 = grp1->genomePosition - grp1->readPosition; // - (pR1->cigarLeftType == 'I'  || pR1->cigarLeftType == 'D' ? pR1->cigarLeft : 0);
    uint32_t        right_r1 = left_r1 + pR1->lenSeq;

    if (*pR1->pCigar == 0) {
        //fprintf(stderr, "{ R1 %u:%s:%u: }\n", left_r1, pChr1->name, left_r1 - pChr1->start);
        buildPlainCIGARorWFA(pRKM1, pR1, left_r1, grp1->isReverse);
        //allocDNA(pR1->lenSeq);
        //pRKM1->pGenome->extractFromGenomicPos(left_r1, pR1->lenSeq, DNA, 0, 0);
        //pAligner->compareSeqAndDNABuildingCIGAR(grp1->isReverse ? pR1->pSeqRevComp : pR1->pSeq, DNA, pRKM1, pR1);
        //left_r1 = grp1->genomePosition - grp1->readPosition - (pR1->cigarLeftType == 'I' || pR1->cigarLeftType == 'D' ? pR1->cigarLeft : 0);
    }
    left_r1 = pR1->cigarLeftPos;
    right_r1 = pR1->cigarRightPos;
    
    ogChromosome   *pChr1 = pRKM1->pGenome->getGenomicCoordinate(left_r1);
    
    pSamWri->writeSAMInfo(
        pR1->pId,
        SAMFLAG_READ_PAIRED | SAMFLAG_MATE_UNMAPPED | (pR1->readIndex == 1 ? SAMFLAG_READ1 : SAMFLAG_READ2) | (grp1->isReverse ? SAMFLAG_READ_REVCOMP : 0),
        pChr1->name,
        left_r1 - pChr1->start + 1, // + (pR1->cigarLeftType == 'S' ? pR1->cigarLeft : 0),
        grp1->score,
        pR1->pCigar,
        (char *) ASTERISK,
        0,
        0,
        (grp1->isReverse ? pR1->pSeqRevComp : pR1->pSeq),
        (grp1->isReverse ? pR1->pQualRev    : pR1->pQual),
        NULL,
        pR1->readNumber,
        grp1->function,
            pR1
    );
    *pR2->pCigar = 0;
    pSamWri->writeSAMInfo(
        pR2->pId,
        SAMFLAG_READ_PAIRED | SAMFLAG_READ_UNMAPPED | (pR1->readIndex == 1 ? SAMFLAG_READ2 : SAMFLAG_READ1)  | (grp1->isReverse ? SAMFLAG_MATE_REVCOMP : 0),
        (char *) ASTERISK,
        0,
        0,
        (char *) ASTERISK,
        pChr1->name,
        left_r1 - pChr1->start + 1, //+ (pR1->cigarLeftType == 'S' ? pR1->cigarLeft : 0),
        0,
        pR2->pSeq,
        pR2->pQual,
        NULL,
        pR2->readNumber,
        '-',
    pR2
    );
}


void ogReadsMapper::writeSamFromReadKeyMapsUnmapped(ogReadKeyMapping *pRKM1, ogReadKeyMapping *pRKM2, ogSamWriter *pSamWri) {
    ogSingleRead    *pR1 = pRKM1->read;
    ogSingleRead    *pR2 = pRKM2->read;
    
    *pR1->pCigar = 0;
    *pR2->pCigar = 0;
    pSamWri->writeSAMInfo(
        pR1->pId,
        SAMFLAG_READ_PAIRED | SAMFLAG_READ_UNMAPPED | SAMFLAG_MATE_UNMAPPED | (pR1->readIndex == 1 ? SAMFLAG_READ1 : SAMFLAG_READ2) ,
        (char *) ASTERISK,
        0,
        0,
        (char *) ASTERISK,
        (char *) ASTERISK,
        0,
        0,
        pR1->pSeq,
        pR1->pQual,
        NULL,
        pR1->readNumber,
        '-',
            pR1
    );
    *pR2->pCigar = 0;
    pSamWri->writeSAMInfo(
        pR2->pId,
        SAMFLAG_READ_PAIRED | SAMFLAG_READ_UNMAPPED | SAMFLAG_MATE_UNMAPPED | (pR1->readIndex == 1 ? SAMFLAG_READ2 : SAMFLAG_READ1),
        (char *) ASTERISK,
        0,
        0,
        (char *) ASTERISK,
        (char *) ASTERISK,
        0,
        0,
        pR2->pSeq,
        pR2->pQual,
        NULL,
        pR2->readNumber,
        '-',
            pR2
    );
}

void ogReadsMapper::writeSamFromReadKeyMapsTranslocated(ogReadKeyMapping *pRKM1, ogReadKeyMapping *pRKM2, ogSamWriter *pSamWri) {
    ogGenomeAndReadPosition *grp1 = pRKM1->pCandPosMan->getkPos(pRKM1->maxScorePos);
    ogGenomeAndReadPosition *grp2 = pRKM2->pCandPosMan->getkPos(pRKM2->maxScorePos);
    ogSingleRead    *pR1 = pRKM1->read;
    ogSingleRead    *pR2 = pRKM2->read;
    uint32_t        left_r1 = grp1->genomePosition - grp1->readPosition; // - (pR1->cigarLeftType == 'I'  || pR1->cigarLeftType == 'D' ? pR1->cigarLeft : 0);
    uint32_t        left_r2 = grp2->genomePosition - grp2->readPosition; // - (pR2->cigarLeftType == 'I'  || pR2->cigarLeftType == 'D' ? pR2->cigarLeft : 0);
    uint32_t        right_r1 = left_r1 + pR1->lenSeq;
    uint32_t        right_r2 = left_r2 + pR2->lenSeq;

    if (*pR1->pCigar == 0) {
        //fprintf(stderr, "{ R1 %u:%s:%u: }\n", left_r1, pChr1->name, left_r1 - pChr1->start);
        buildPlainCIGARorWFA(pRKM1, pR1, left_r1, grp1->isReverse);
        //allocDNA(pR1->lenSeq);
        //pRKM1->pGenome->extractFromGenomicPos(left_r1, pR1->lenSeq, DNA, 0, 0);
        //pAligner->compareSeqAndDNABuildingCIGAR(grp1->isReverse ? pR1->pSeqRevComp : pR1->pSeq, DNA, pRKM1, pR1);
        //left_r1 = grp1->genomePosition - grp1->readPosition - (pR1->cigarLeftType == 'I' || pR1->cigarLeftType == 'D' ? pR1->cigarLeft : 0);
    }
    left_r1 = pR1->cigarLeftPos;
    right_r1 = pR1->cigarRightPos;
    if (*pR2->pCigar == 0) {
        //fprintf(stderr, "{ R2 %u:%s:%u: }\n", left_r2, pChr1->name, left_r2 - pChr1->start);
        buildPlainCIGARorWFA(pRKM2, pR2, left_r2, grp2->isReverse);
        //allocDNA(pR2->lenSeq);
        //pRKM2->pGenome->extractFromGenomicPos(left_r2, pR2->lenSeq, DNA, 0, 0);
        //pAligner->compareSeqAndDNABuildingCIGAR(grp2->isReverse ? pR2->pSeqRevComp : pR2->pSeq, DNA, pRKM2, pR2);
        //left_r2 = grp2->genomePosition - grp2->readPosition - (pR2->cigarLeftType == 'I' || pR2->cigarLeftType == 'D' ? pR2->cigarLeft : 0);
    }
    left_r2 = pR2->cigarLeftPos;
    right_r2 = pR2->cigarRightPos;
    
    ogChromosome   *pChr1 = pRKM1->pGenome->getGenomicCoordinate(left_r1);
    ogChromosome   *pChr2 = pRKM2->pGenome->getGenomicCoordinate(left_r2);
    
    pSamWri->writeSAMInfo(
        pR1->pId,
        SAMFLAG_READ_PAIRED | SAMFLAG_READS_MAPPED | (pR1->readIndex == 1 ? SAMFLAG_READ1 : SAMFLAG_READ2) | (grp1->isReverse ? SAMFLAG_READ_REVCOMP : 0) | (grp2->isReverse ? SAMFLAG_MATE_REVCOMP : 0),
        pChr1->name,
        left_r1 - pChr1->start + 1, // + (pR1->cigarLeftType == 'S' ? pR1->cigarLeft : 0),
        grp1->score,
        pR1->pCigar,
        pChr2->name,
        left_r2 - pChr2->start + 1, // + (pR2->cigarLeftType == 'S' ? pR2->cigarLeft : 0),
        0,
        (grp1->isReverse ? pR1->pSeqRevComp : pR1->pSeq),
        (grp1->isReverse ? pR1->pQualRev    : pR1->pQual),
        NULL,
        pR1->readNumber,
        grp1->function,
            pR1
    );
    pSamWri->writeSAMInfo(
        pR2->pId,
        SAMFLAG_READ_PAIRED | SAMFLAG_READS_MAPPED | (pR1->readIndex == 1 ? SAMFLAG_READ2 : SAMFLAG_READ1)  | (grp2->isReverse ? SAMFLAG_READ_REVCOMP : 0) | (grp1->isReverse ? SAMFLAG_MATE_REVCOMP : 0),
        pChr2->name,
        left_r2 - pChr2->start + 1, // + (pR2->cigarLeftType == 'S' ? pR2->cigarLeft : 0),
        grp2->score,
        pR2->pCigar,
        pChr1->name,
        left_r1 - pChr1->start + 1, // + (pR1->cigarLeftType == 'S' ? pR1->cigarLeft : 0),
        0,
        (grp2->isReverse ? pR2->pSeqRevComp : pR2->pSeq),
        (grp2->isReverse ? pR2->pQualRev    : pR2->pQual),
        NULL,
        pR2->readNumber,
        grp2->function,
            pR2
    );
}


// Se asume que el que se escribe es el "2" [grp2] (usando el "1" como referencia [grp1])
void ogReadsMapper::writeAltRd2SamFromGenomeAndReadPosition(ogReadKeyMapping *pRKM1, ogReadKeyMapping *pRKM2, ogGenomeAndReadPosition *grp1, ogGenomeAndReadPosition *grp2, ogSamWriter *pSamWri) {
    ogSingleRead    *pR1 = pRKM1->read;
    ogSingleRead    *pR2 = pRKM2->read;
    
    uint32_t        left_r1 = grp1->genomePosition - grp1->readPosition;// - (pR1->cigarLeftType == 'I'  || pR1->cigarLeftType == 'D' ? pR1->cigarLeft : 0);
    uint32_t        left_r2 = grp2->genomePosition - grp2->readPosition;// - (pR2->cigarLeftType == 'I'  || pR2->cigarLeftType == 'D' ? pR2->cigarLeft : 0);
    uint32_t        right_r1 = left_r1 + pR1->lenSeq;
    uint32_t        right_r2 = left_r2 + pR2->lenSeq;

    //if (*pR2->pCigar == 0) {
    // El cigar siempre se tiene que recalcular ya que es alternativo
        //fprintf(stderr, "{ R2 %u:%s:%u: }\n", left_r2, pChr1->name, left_r2 - pChr1->start);
        buildPlainCIGARorWFA(pRKM2, pR2, left_r2, grp2->isReverse);
        //left_r2 = grp2->genomePosition - grp2->readPosition - (pR2->cigarLeftType == 'I' || pR2->cigarLeftType == 'D' ? pR2->cigarLeft : 0);
        left_r2 = pR2->cigarLeftPos;
        right_r2 = pR2->cigarRightPos;
        //allocDNA(pR2->lenSeq);
        //pRKM1->pGenome->extractFromGenomicPos(left_r2, pR2->lenSeq, DNA, 0, 0);
        //Aligner->compareSeqAndDNABuildingCIGAR(grp2->isReverse ? pR2->pSeqRevComp : pR2->pSeq, DNA, pRKM2, pR2);
    //}
        
    uint32_t        left = left_r1 < left_r2 ? left_r1 : left_r2;
    //uint32_t        right = right_r1 > right_r2 ? right_r1 : right_r2;
    //uint32_t        genDist = (left > right ? left-right : right-left);
    //uint32_t        genDist1 = (left_r1 > right_r1 ? left_r1 - right_r1 : right_r1 - left_r1);
    uint32_t        genDist2 = (left_r2 > right_r2 ? left_r2 - right_r2 : right_r2 - left_r2);
    ogChromosome   *pChr1 = pRKM1->pGenome->getGenomicCoordinate(left);
    
    pSamWri->writeSAMInfo(
        pR2->pId,
        SAMFLAG_READ_PAIRED | SAMFLAG_READS_MAPPED | (pR1->readIndex == 1 ? SAMFLAG_READ1 : SAMFLAG_READ2) | (grp2->isReverse ? SAMFLAG_READ_REVCOMP : 0) | (grp1->isReverse ? SAMFLAG_MATE_REVCOMP : 0) |
        SAMFLAG_NOT_PRIMARY_ALIGNMENT,
        pChr1->name,
        left_r2 - pChr1->start + 1, //  + (pR2->cigarLeftType == 'S' ? pR2->cigarLeft : 0),
        grp2->score,
        pR2->pCigar,
        (char *) EQUAL,
        left_r1 - pChr1->start + 1, // + (pR2->cigarLeftType == 'S' ? pR2->cigarLeft : 0),
        genDist2,
        (char *) ASTERISK, // to reduce file size as specified in SAM format
        (char *) ASTERISK,  // to reduce file size as specified in SAM format
        NULL,
        pR2->readNumber,
        grp2->function,
            pR2
    );
}


// Se asume que el que se escribe es el "1" [grp2] (usando el "1" como referencia [grp1])
void ogReadsMapper::writeAltRd1SamFromGenomeAndReadPosition(ogReadKeyMapping *pRKM1, ogReadKeyMapping *pRKM2, ogGenomeAndReadPosition *grp1, ogGenomeAndReadPosition *grp2, ogSamWriter *pSamWri) {
    ogSingleRead    *pR1 = pRKM1->read;
    ogSingleRead    *pR2 = pRKM2->read;
    
    uint32_t        left_r1 = grp1->genomePosition - grp1->readPosition; // - (pR1->cigarLeftType == 'I'  || pR1->cigarLeftType == 'D' ? pR1->cigarLeft : 0);
    uint32_t        left_r2 = grp2->genomePosition - grp2->readPosition; // - (pR2->cigarLeftType == 'I'  || pR2->cigarLeftType == 'D' ? pR2->cigarLeft : 0);
    uint32_t        right_r1 = left_r1 + pR1->lenSeq;
    uint32_t        right_r2 = left_r2 + pR2->lenSeq;
    //if (*pR2->pCigar == 0) {
    // El cigar siempre se tiene que recalcular ya que es alternativo
        //fprintf(stderr, "{ R1 %u:%s:%u: }\n", left_r1, pChr1->name, left_r1 - pChr1->start);
        buildPlainCIGARorWFA(pRKM1, pR1, left_r1, grp1->isReverse);
        //left_r1 = grp1->genomePosition - grp1->readPosition - (pR1->cigarLeftType == 'I' || pR1->cigarLeftType == 'D' ? pR1->cigarLeft : 0);
        left_r1 = pR1->cigarLeftPos;
        right_r1 = pR1->cigarRightPos;
        //allocDNA(pR1->lenSeq);
        //pRKM1->pGenome->extractFromGenomicPos(left_r1, pR1->lenSeq, DNA, 0, 0);
        //pAligner->compareSeqAndDNABuildingCIGAR(grp1->isReverse ? pR1->pSeqRevComp : pR1->pSeq, DNA, pRKM1, pR1);
    //}
        
    uint32_t        left = left_r1 < left_r2 ? left_r1 : left_r2;
    //uint32_t        right = right_r1 > right_r2 ? right_r1 : right_r2;
    //uint32_t        genDist = (left > right ? left-right : right-left);
    uint32_t        genDist1 = (left_r1 > right_r1 ? left_r1 - right_r1 : right_r1 - left_r1);
    //uint32_t        genDist2 = (left_r2 > right_r2 ? left_r2 - right_r2 : right_r2 - left_r2);
    ogChromosome   *pChr1 = pRKM1->pGenome->getGenomicCoordinate(left);
    
    pSamWri->writeSAMInfo(
        pR1->pId,
        SAMFLAG_READ_PAIRED | SAMFLAG_READS_MAPPED | (pR1->readIndex == 1 ? SAMFLAG_READ1 : SAMFLAG_READ2)  | (grp1->isReverse ? SAMFLAG_READ_REVCOMP : 0) | (grp2->isReverse ? SAMFLAG_MATE_REVCOMP : 0) |
        SAMFLAG_NOT_PRIMARY_ALIGNMENT,
        pChr1->name,
        left_r1 - pChr1->start + 1, // + (pR1->cigarLeftType == 'S' ? pR1->cigarLeft : 0),
        grp1->score,
        pR1->pCigar,
        (char *) EQUAL,
        left_r2 - pChr1->start + 1, // + (pR2->cigarLeftType == 'S' ? pR2->cigarLeft : 0),
        genDist1,
        (char *) ASTERISK, // to reduce file size as specified in SAM format
        (char *) ASTERISK,  // to reduce file size as specified in SAM format
        NULL,
        pR1->readNumber,
        grp1->function,
            pR1
    );
}

/**
void ogReadsMapper::writeSamFromBothReadsInfo(sequenceRead *seqRead, ogSamWriter *pSamWri) {
    char            isRd2 = seqRead->genPos1_isRd2;
    ogChromosome   *pChr1 = pRdKeyMapRd1->pGenome->getGenomicCoordinate(seqRead->genPos1_Left);
    uint32_t        gPos1 = seqRead->genPos1_Left - pChr1->start + 1;
    uint32_t        gPos2 = seqRead->genPos2_Left - pChr1->start + 1;
    //ogChromosome   *pChr2 = pRdKeyMapRd1->pGenome->getGenomicCoordinate(gPos2);
    if (isRd2) {
        // Primero escribir Rd2 que es el "primario"
        if (*seqRead->pCigar2 == 0) {
            // En teoria no debería entrar aquí
            allocDNA(seqRead->seqLen2);
            fprintf(stderr, "[ %u:%s:%u: ]\n", seqRead->genPos1_Left, pChr1->name, gPos1);
            pRdKeyMapRd1->pGenome->extractFromGenomicPos(pChr1->start+gPos1, seqRead->seqLen2, DNA, 0, 0);
            pAligner->compareSeqAndDNABuildingCIGAR((seqRead->posIsReverse1 ? seqRead->pSeqRevComp2: seqRead->pSeq2), DNA, pRdKeyMapRd1, seqRead);
        }
        pSamWri->writeSAMInfo(
            seqRead->pId2,
            SAMFLAG_READ_PAIRED | SAMFLAG_READS_MAPPED | SAMFLAG_READ2 | (seqRead->posIsReverse1 ? SAMFLAG_READ_REVCOMP : 0) | (seqRead->posIsReverse2 ? SAMFLAG_MATE_REVCOMP : 0),
            pChr1->name,
            gPos1,
            seqRead->genPos1_mapScore,
            seqRead->pCigar2,
            (char *) EQUAL,
            gPos2,
            (gPos1 < gPos2 ? 1 : -1) * seqRead->genPos_insertSize,
            (seqRead->posIsReverse1 ? seqRead->pSeqRevComp2 : seqRead->pSeq2),
            (seqRead->posIsReverse1 ? seqRead->pQualRev2    : seqRead->pQual2),
            NULL,
            seqRead->readNumber,
            '-'
        );
        if (*seqRead->pCigar1 == 0) {
            //fprintf(stderr, "[[ %u:%s:%u: ]]\n", seqRead->genPos2_Left, pChr1->name, gPos2);
            pRdKeyMapRd1->pGenome->extractFromGenomicPos(pChr1->start+gPos2-1, seqRead->seqLen1, DNA, 0, 0);
            pAligner->compareSeqAndDNABuildingCIGAR((seqRead->posIsReverse2 ? seqRead->pSeqRevComp1: seqRead->pSeq1), DNA, seqRead->seqLen1, seqRead->pCigar1, pRdKeyMapRd1, &seqRead->leftSoftR1);
        }
        pSamWri->writeSAMInfo(
            seqRead->pId1,
            SAMFLAG_READ_PAIRED | SAMFLAG_READS_MAPPED | SAMFLAG_READ1 | (seqRead->posIsReverse2 ? SAMFLAG_READ_REVCOMP : 0) | (seqRead->posIsReverse1 ? SAMFLAG_MATE_REVCOMP : 0),
            pChr1->name,
            gPos2,
            seqRead->genPos2_mapScore,
            seqRead->pCigar1,
            (char *) EQUAL,
            gPos1,
            (gPos2 < gPos1 ? 1 : -1) * seqRead->genPos_insertSize,
            (seqRead->posIsReverse2 ? seqRead->pSeqRevComp1 : seqRead->pSeq1),
            (seqRead->posIsReverse2 ? seqRead->pQualRev1    : seqRead->pQual1),
            NULL,
            seqRead->readNumber,
            '-'
        );
    } else {
        if (*seqRead->pCigar1 == 0) {
            fprintf(stderr, "{ %u:%s:%u: }\n", seqRead->genPos1_Left, pChr1->name, gPos1);
            pRdKeyMapRd1->pGenome->extractFromGenomicPos(pChr1->start+gPos1-1, seqRead->seqLen1, DNA, 0, 0);
            pAligner->compareSeqAndDNABuildingCIGAR((seqRead->posIsReverse1 ? seqRead->pSeqRevComp1: seqRead->pSeq1), DNA, seqRead->seqLen1, seqRead->pCigar1, pRdKeyMapRd1, &seqRead->leftSoftR1);
        }
        pSamWri->writeSAMInfo(
            seqRead->pId1,
            SAMFLAG_READ_PAIRED | SAMFLAG_READS_MAPPED | SAMFLAG_READ1 | (seqRead->posIsReverse1 ? SAMFLAG_READ_REVCOMP : 0) | (seqRead->posIsReverse2 ? SAMFLAG_MATE_REVCOMP : 0),
            pChr1->name,
            gPos1,
            seqRead->genPos1_mapScore,
            seqRead->pCigar1,
            (char *) EQUAL,
            gPos2,
            (gPos1 < gPos2 ? 1 : -1) * seqRead->genPos_insertSize,
            (seqRead->posIsReverse1 ? seqRead->pSeqRevComp1 : seqRead->pSeq1),
            (seqRead->posIsReverse1 ? seqRead->pQualRev1    : seqRead->pQual1),
            NULL,
            seqRead->readNumber,
            '-'
        );
        if (*seqRead->pCigar2 == 0) {
            // calculate cigar
            // extract DNA
            //fprintf(stderr, "{{ %u:%s:%u: }}\n", seqRead->genPos1_Left, pChr1->name, gPos1);
            pRdKeyMapRd1->pGenome->extractFromGenomicPos(pChr1->start+gPos2-1, seqRead->seqLen2, DNA, 0, 0);
            pAligner->compareSeqAndDNABuildingCIGAR((seqRead->posIsReverse2 ? seqRead->pSeqRevComp2: seqRead->pSeq2), DNA, seqRead->seqLen2, seqRead->pCigar2, pRdKeyMapRd1, &seqRead->leftSoftR2);
        }
        pSamWri->writeSAMInfo(
            seqRead->pId2,
            SAMFLAG_READ_PAIRED | SAMFLAG_READS_MAPPED | SAMFLAG_READ2 | (seqRead->posIsReverse2 ? SAMFLAG_READ_REVCOMP : 0) | (seqRead->posIsReverse1 ? SAMFLAG_MATE_REVCOMP : 0),
            pChr1->name,
            gPos2,
            seqRead->genPos2_mapScore,
            seqRead->pCigar2,
            (char *) EQUAL,
            gPos1,
            (gPos2 < gPos1 ? 1 : -1) * seqRead->genPos_insertSize,
            (seqRead->posIsReverse2 ? seqRead->pSeqRevComp2 : seqRead->pSeq2),
            (seqRead->posIsReverse2 ? seqRead->pQualRev2    : seqRead->pQual2),
            NULL,
            seqRead->readNumber,
            '-'
        );
    }
}

void ogReadsMapper::writeSamFromRead1MapRead2Unmapped(sequenceRead *seqRead, ogSamWriter *pSamWri) {
    ogChromosome   *pChr1 = pRdKeyMapRd1->pGenome->getGenomicCoordinate(seqRead->genPos1_Left);
    uint32_t        gPos1 = seqRead->genPos1_Left - pChr1->start + 1 + seqRead->leftSoftR1;
    uint32_t        gPos2 = seqRead->genPos2_Left - pChr1->start + 1 + seqRead->leftSoftR2;
    //ogChromosome   *pChr2 = pRdKeyMapRd1->pGenome->getGenomicCoordinate(gPos2);

    pSamWri->writeSAMInfo(
        seqRead->pId1,
        SAMFLAG_READ_PAIRED | SAMFLAG_MATE_UNMAPPED | SAMFLAG_READ1 | (seqRead->posIsReverse1 ? SAMFLAG_READ_REVCOMP : 0) | (seqRead->posIsReverse2 ? SAMFLAG_MATE_REVCOMP : 0),
        pChr1->name,
        gPos1,
        seqRead->genPos1_mapScore,
        seqRead->pCigar1,
        (char *) EQUAL,
        gPos1,
        0,
        (seqRead->posIsReverse1 ? seqRead->pSeqRevComp1 : seqRead->pSeq1),
        (seqRead->posIsReverse1 ? seqRead->pQualRev1    : seqRead->pQual1),
        NULL,
        seqRead->readNumber,
        '-'
    );
    pSamWri->writeSAMInfo(
        seqRead->pId2,
        SAMFLAG_READ_PAIRED | SAMFLAG_READ_UNMAPPED | SAMFLAG_READ2 | (seqRead->posIsReverse2 ? SAMFLAG_READ_REVCOMP : 0) | (seqRead->posIsReverse1 ? SAMFLAG_MATE_REVCOMP : 0),
        pChr1->name,
        gPos1,
        0,
        (char *) ASTERISK,
        (char *) EQUAL,
        gPos1,
        0,
        (seqRead->posIsReverse2 ? seqRead->pSeqRevComp2 : seqRead->pSeq2),
        (seqRead->posIsReverse2 ? seqRead->pQualRev2    : seqRead->pQual2),
        NULL,
        seqRead->readNumber,
        '-'
    );
}


void ogReadsMapper::writeSamFromRead2MapRead1Unmapped(sequenceRead *seqRead, ogSamWriter *pSamWri) {
    ogChromosome   *pChr1 = pRdKeyMapRd1->pGenome->getGenomicCoordinate(seqRead->genPos1_Left);
    uint32_t        gPos1 = seqRead->genPos1_Left - pChr1->start + 1 + seqRead->leftSoftR1;
    uint32_t        gPos2 = seqRead->genPos2_Left - pChr1->start + 1 + seqRead->leftSoftR2;
    //ogChromosome   *pChr2 = pRdKeyMapRd1->pGenome->getGenomicCoordinate(gPos2);
    pSamWri->writeSAMInfo(
        seqRead->pId1,
        SAMFLAG_READ_PAIRED | SAMFLAG_READ_UNMAPPED | SAMFLAG_READ1 | (seqRead->posIsReverse2 ? SAMFLAG_READ_REVCOMP : 0) | (seqRead->posIsReverse1 ? SAMFLAG_MATE_REVCOMP : 0),
        pChr1->name,
        gPos1,
        0,
        (char *) ASTERISK,
        (char *) EQUAL,
        gPos1,
        0,
        (seqRead->posIsReverse1 ? seqRead->pSeqRevComp1 : seqRead->pSeq1),
        (seqRead->posIsReverse1 ? seqRead->pQualRev1    : seqRead->pQual1),
        NULL,
        seqRead->readNumber,
        '-'
    );
    // Luego escribir Rd 2 aunque en mapping fue Rd2 y Rd1
    pSamWri->writeSAMInfo(
        seqRead->pId2,
        SAMFLAG_READ_PAIRED | SAMFLAG_MATE_UNMAPPED | SAMFLAG_READ2 | (seqRead->posIsReverse1 ? SAMFLAG_READ_REVCOMP : 0) | (seqRead->posIsReverse2 ? SAMFLAG_MATE_REVCOMP : 0),
        pChr1->name,
        gPos1,
        seqRead->genPos2_mapScore,
        seqRead->pCigar2,
        (char *) EQUAL,
        gPos1,
        0,
        (seqRead->posIsReverse2 ? seqRead->pSeqRevComp2 : seqRead->pSeq2),
        (seqRead->posIsReverse2 ? seqRead->pQualRev2    : seqRead->pQual2),
        NULL,
        seqRead->readNumber,
        '-'
    );
}

 **/

uint32_t ogReadsMapper::mapTheReadWithKeyMap(ogReadKeyMapping *pRdKeyMap) {
    
    if (pRdKeyMap->mapped) return pRdKeyMap->pCandPosMan->getCount();
    
    ////////////////////////////////////////////////////////////////
    /* PARAMETERS */
    ////////////////////////////////////////////////////////////////
    
    clock_t prev_t, start_t = clock();
    clock_t curr_t = start_t;
    busy = 1;
    uint32_t    j, k=0, n = 0, fwdK, revK;
    char f = 0;
    ogReadMapperCounts *pC;
    uint64_t readNum = pRdKeyMap->read->readNumber; //seqRead->readNumber;
    pRdKeyMap->debug = 0;
    
    nMaps++;
    pC = pCounts;
    curr_t = clock();
    preparation += duration_microsec(curr_t, start_t);
    for (j=0; j < pScheduler->nScheduleFunc; j++) {
        prev_t = curr_t;
        if (j > 0) pRdKeyMap->reactivateKeys();
        pRdKeyMap->pCandPosMan->currentFunction = pScheduler->funcKey[j];
        k = pScheduler->funcSchedule[j](pRdKeyMap);
        n += k;
        pC->calls++;
        if (k) {
            pC->hitting++;
            pC->hits += k;
            //
            pRdKeyMap->pCandPosMan->cleanCandidatePositions(); // Currently, it does nothing.
            //
            
            if (true || k < TOO_MANY_CANDIDATE_REGIONS) {
                if (false) {
                    // Pre-filter : remove those hitting very many, but why?
                    if (k > ACEPTABLE_CANDIDATE_REGIONS) {
                        fwdK = pRdKeyMap->pCandPosMan->getFwdCount();
                        revK = pRdKeyMap->pCandPosMan->getRevCount();
                        if (revK > 0 && fwdK > 0) {
                            // Hits in both strands
                            if (revK <= 4) {
                                if (fwdK > revK << 5) {
                                    // complementary sequence may have large unspecific targets
                                    pRdKeyMap->pCandPosMan->removeForwardTargets();
                                }
                            } else if (fwdK <= 4) {
                                if (revK > fwdK << 5) {
                                    pRdKeyMap->pCandPosMan->removeReverseTargets();
                                }
                            }
                        } else {
                            // Can the hits be reduced?                    
                        }
                        k = pRdKeyMap->pCandPosMan->getCount();
                    }
                }
                k = pAligner->align(pRdKeyMap->read, pRdKeyMap, compareKey);
                if (k) {
                    pC->matches += k;
                    pRdKeyMap->pCandPosMan->estimateFwdRevMatches(0);
                    pC->foundFwd += pRdKeyMap->pCandPosMan->getFwdMatches();
                    nFoundFwd += pRdKeyMap->pCandPosMan->getFwdMatches();
                    pC->foundRev += pRdKeyMap->pCandPosMan->getRevMatches();
                    nFoundRev += pRdKeyMap->pCandPosMan->getRevMatches();
                    f = 1;
                } else {
                    //fprintf(stderr, "*** FAIL ***\n");
                    //pRdKeyMap->printNotFound(readNum);
                    pC->failures++;
                    //fprintf(stderr, "Failure #1:%lld:%s\n",readNum,seq);
                    //pRdKeyMap->pCandPosMan->printWithGenomePositions(pRdKeyMap->pGenome);
                }
                
                
                
                ///////////////////////////////////////////////////////////////////////
                //////  OLD CODE BELOW
                ///////////////////////////////////////////////////////////////////////
                /******
                if (false) {
                    // OLD CODE
                    if (checkMatch) {
                        //fprintf(stderr, "readNum=%llu\n",readNum);
                        k = pAligner->checkKeySeq(seq, seqLen, pRdKeyMap);
                    } else {
                        pRdKeyMap->pCandPosMan->setAllStatus('K'); // Set Validated Key for alignments
                    }
                    if (k) {
                        k = pAligner->align(seqRead, pRdKeyMap, compareKey);
                        if (k) {
                            pC->matches += k;
                            pRdKeyMap->pCandPosMan->estimateFwdRevMatches();
                            pC->foundFwd += pRdKeyMap->pCandPosMan->getFwdMatches();
                            nFoundFwd += pRdKeyMap->pCandPosMan->getFwdMatches();
                            pC->foundRev += pRdKeyMap->pCandPosMan->getRevMatches();
                            nFoundRev += pRdKeyMap->pCandPosMan->getRevMatches();
                            f = 1;
                        } else {
                            //fprintf(stderr, "*** FAIL ***\n");
                            //pRdKeyMap->printNotFound(readNum);
                            pC->failures++;
                            //fprintf(stderr, "Failure #1:%lld:%s\n",readNum,seq);
                            //pRdKeyMap->pCandPosMan->printWithGenomePositions(pRdKeyMap->pGenome);
                        }
                    } else {
                        pC->failures++;
                        //fprintf(stderr, "Failure #2:RdNum=%lld:%s\n",readNum,seq);
                        //pRdKeyMap->pCandPosMan->printWithGenomePositions(pRdKeyMap->pGenome);
                        if (readNum == 56) {
                            fprintf(stderr, "Activating debug\n");
                            pRdKeyMap->debug = 1;
                            pMapParams->funcSchedule[j](pRdKeyMap);
                            pRdKeyMap->debug = 0;
                            fprintf(stderr, "Deactivating debug\n");
                        }
                    }
                }
                //////  OLD CODE ABOVE
                *****/
            } else {
                fprintf(stderr, "Too many candidate regions[%u]. Not mapped. ReadNum=%llu\n%s\n", k, pRdKeyMap->read->readNumber, pRdKeyMap->read->pSeq);
                k = 0;
            }
        }
        //else {
        //    fprintf(stderr, "Not mapped. ReadNum=%llu\n%s\n", pRdKeyMap->read->readNumber, pRdKeyMap->read->pSeq);
        //}
        curr_t = clock();
        pC->elapsed += duration_microsec(curr_t, prev_t);
        if (k && productionMode) {
            break;
        }
        pC++;
    }
    nFound += f;
    nMatches += n;
    elapsed += duration_microsec(clock(), start_t);
    busy = 0;
    //fprintf(stderr, "done mapRead %llu %p\n",readNum,seqRead); fflush(stderr);
//    canWork = 0;
    pRdKeyMap->mapped = 1;
    return k;
}




uint32_t ogReadsMapper::countTheReadWithKeyMap(ogReadKeyMapping *pRdKeyMap) {
    
    if (pRdKeyMap->mapped) return pRdKeyMap->pCandPosMan->getCount();
    
    clock_t prev_t, start_t = clock();
    clock_t curr_t = start_t;
    busy = 1;
    uint32_t    j, k=0, n = 0;
    char f = 0;
    ogReadMapperCounts *pC;
    
    nMaps++;
    pC = pCounts;
    curr_t = clock();
    preparation += duration_microsec(curr_t, start_t);
    for (j=0; j < pScheduler->nScheduleFunc; j++) {
        prev_t = curr_t;
        //if (j > 0) pRdKeyMap->reactivateKeys();
        pRdKeyMap->pCandPosMan->currentFunction = pScheduler->funcKey[j];
        k = pScheduler->funcSchedule[j](pRdKeyMap);
        n += k;
        pC->calls++;
        if (k) {
            f = 1;
            pC->hitting++;
            pC->hits += k;
                pC->matches += k;
                pRdKeyMap->pCandPosMan->estimateFwdRevMatches(1);
                pC->foundFwd += pRdKeyMap->pCandPosMan->getFwdMatches();
                nFoundFwd += pRdKeyMap->pCandPosMan->getFwdMatches();
                pC->foundRev += pRdKeyMap->pCandPosMan->getRevMatches();
                nFoundRev += pRdKeyMap->pCandPosMan->getRevMatches();
        }
        curr_t = clock();
        pC->elapsed += duration_microsec(curr_t, prev_t);
        if (k && productionMode) {
            break;
        }
        pC++;
    }
    nFound += f;
    nMatches += n;
    elapsed += duration_microsec(clock(), start_t);
    busy = 0;
    pRdKeyMap->mapped = 1;
    return k;
}


//void ogReadsMapper::mapReadPairs(char *seq1, uint32_t seq1Len, char *seq2, uint32_t seq2Len) {
//    
//}

void ogReadsMapper::printSummary(uint64_t totalReads) {
    
    fprintf(stderr, "Keys/Read=%.2f, k=%llu, r=%llu, Maps=%llu (%.1f%%), Found=%llu (%.3f%%), Not Found=%llu (%.3f%%), Alignments=%llu (%.3f%%), Matches=%llu, fFwd=%llu, fRev=%llu, Fwd/Rev=%.2f, Matches/Call=%.1f\n",
        (double) nKeys / (double) (nReads), nKeys, nReads,
        nMaps, (float) nMaps*100 / (float) totalReads,
        nFound, //nFoundFwd+nFoundRev,
        100 * (float) (nFound) / ((float) (nMaps) + 0.00001), // nFoundFwd+nFoundRev
        nMaps - nFound,
        100 * (float) (nMaps - nFound) / nMaps,
        pAligner->getFullAlignments(), (100 * ((float) (pAligner->getFullAlignments()) / nMaps)),
        nMatches, 
        nFoundFwd, nFoundRev, (float) nFoundFwd / (float) (nFoundRev),
        (float) nMatches / (float) nMaps
        );
    fprintf(stderr, "Elapsed=%4.3f sec (includes Preparation=%4.3f sec), Mate Pairing=%4.3f sec\n",
        (float) elapsed / (float) 1e6, 
        (float) preparation / (float) 1e6,
        (float) pairing / (float) 1e6
    );
    fprintf(stderr, "Wasting Events = %llu : Cycles = %llu * %llu\n", wastingEvents, wastingCycles2+1, wastingCycles);
    uint64_t r =     nReadsMappedR1R2+nReadsMappedR2R1+nReadsMappedR1+nReadsMappedR2+nReadsMappedTrans+nReadsMappedOther; // +nReadsUnmapped
    if (r > 0) {
        fprintf(stderr, "%llu Reads: R1R2=%llu, R2R1=%llu, R1+R2-=%llu, R1-R2+=%llu, Trans=%llu, Other=%llu, Unmapped=%llu. (Altern=%llu)\n",
        nReads,
        nReadsMappedR1R2,
        nReadsMappedR2R1,
        nReadsMappedR1,
        nReadsMappedR2,
        nReadsMappedTrans,
        nReadsMappedOther,
        nReadsUnmapped,
        nReadsAlternMaps);
        fprintf(stderr, "%% Reads: R1R2=%.3f%%, R2R1=%.3f%%, R1+R2-=%.3f%%, R1-R2+=%.3f%%, Trans=%.3f%%, Other=%.3f%%, Unmapped=%.3f%%, (Altern=%.3f%%)\n",
        (float) (nReadsMappedR1R2*100)/nReads,
        (float) (nReadsMappedR2R1*100)/nReads,
        (float) (nReadsMappedR1*100)/nReads,
        (float) (nReadsMappedR2*100)/nReads,
        (float) (nReadsMappedTrans*100)/nReads,
        (float) (nReadsMappedOther*100)/nReads,
        (float) (nReadsUnmapped*100)/nReads,
        (float) (nReadsAlternMaps*100)/nReads);
    }
    
    int i;
    // FROM: https://stackoverflow.com/questions/24281603/c-underline-output
    // https://misc.flogisoft.com/bash/tip_colors_and_formatting
    char normal[]={0x1b,'[','0',';','3','9','m',0};
    char Uyellow[]={0x1b,'[','4',';','9', '3', 'm',0};
    char yellow[]={0x1b,'[','0',';','9', '3', 'm',0};
    char bold[]={0x1b,'[','1',';','m',0};
    char Ubold[]={0x1b,'[','4',';','1', 'm',0};
    char UXbold[]={0x1b,'[','4',';','1',';','3', '3', 'm',0};
    char UX[]={0x1b,'[','4',';','3', '3', 'm',0};
    char X[]={0x1b,'[','0',';','3', '3', 'm',0};
    fprintf(stderr, "%sFun|      Calls|  %%Hits| %%Fails|%%RelFnd|Time (sec)|  %% | (ms)/Call| Hits/Fd|Mtchs/Fd|Fwd/Rev|%s\n",UXbold,normal);
    ogReadMapperCounts *pC;
    uint64_t totElapsed = 0;
    for (pC = pCounts, i=0; i < nFunc; i++, pC++) totElapsed += pC->elapsed;
    for (pC = pCounts, i=0; i < nFunc; i++, pC++) {
        uint64_t found = (pC->foundFwd + pC->foundRev);
        //fprintf(stderr, "%s[%c] |%11llu|%11llu|%6.2f%%|%6.2f%%|%11llu|%6.2f%%|%13llu|%10.3f|%14.7f|%10llu|%10llu|%8.2f|%7.2f|\n",
        fprintf(stderr, "%s %c |%11llu|%6.2f%%|%6.2f%%|%6.2f%%|%10.3f|%3.0f%%|%10.7f|%8.2f|%8.2f|%7.4f|\n",
                (i == nFunc-1 ? UX : X),    // Color
                pScheduler->funcKey[i],     // Func Character
                pC->calls,                  // Calls
                100 * (float) (pC->hitting) / ((float) pC->calls + 0.00001),        // Hits
                100 * (float) (pC->failures) / ((float) pC->hitting + 0.00001),     // Failures
                100 * (float) (pC->hitting - pC->failures) / ((float) nFound + 0.00001),     // RelFnd
                
                //100 * (float) (found) / ((float) pC->calls + 0.00001),              // RelFnd
                //100 * (float) (found) / ((float) (nFoundFwd+nFoundRev) + 0.00001),  // AbsFnd
                //pC->matches,
                (float) pC->elapsed / (float) 1e6,
                (float) pC->elapsed * 100 / (float) totElapsed,
                ((float) pC->elapsed * 1000) / (float) (1e6 * pC->calls),
                //pC->foundFwd,
                //pC->foundRev,
                (float) pC->hits / (float) (found == 0 ? 1 : found),
                (float) pC->matches / (float) (found == 0 ? 1 : found),
                (float) pC->foundFwd / ((float) pC->foundRev + 0.01)
                );
    }
    fprintf(stderr, "%s\n" ,normal);
}

char ogReadsMapper::isBusy() {
    return (busy != 0);
}

char ogReadsMapper::isAvailable() {
    return (busy == 0);
}

char ogReadsMapper::hasBeenBusy() {
    return (nMaps > 0);
}

void ogReadsMapper::process1Read() {
    if (! isQueueEmpty()) {
        ogQueueMapper *pqm;
        mtx.lock();
            pqm = qSeq.front();
            qSeq.pop();
            //fprintf(stderr, "Queue %p size %lu\n", &qSeq, qSeq.size());
            //sequenceRead *seqR = qSeq.pull();
        mtx.unlock();
        ogSingleRead *read1;
        if (pqm->pRead == NULL) {
            ogFastAQsequence *pFAS1 = pqm->pFAS1;
            ogFastAQsequence *pFAS2 = pqm->pFAS2;
            read1 = new ogSingleRead();
            read1->setRead(pFAS1->name.s, pFAS1->name.l, pFAS1->seq.s, pFAS1->seq.l, pFAS1->qual.s, pFAS1->qual.l, pqm->readNum, 1);
            read1->packIfNeeded(pMapParams->pGenome);
            if (pFAS2 != NULL) {
                ogSingleRead *read2 = new ogSingleRead();
                read2->setRead(pFAS2->name.s, pFAS2->name.l, pFAS2->seq.s, pFAS2->seq.l, pFAS2->qual.s, pFAS2->qual.l, pqm->readNum, 1);
                read2->packIfNeeded(pMapParams->pGenome);
                read1->setPairedRead(read2);
                //pFAS2->pBuffer->readsPending[thread]--;
                (*pFAS2->pRelease)++;
                free(pFAS2);
            }
            //pFAS1->pBuffer->readsPending[thread]--;
            (*pFAS1->pRelease)++;
            free(pFAS1);
        } else {
            // es del tipo ya preparado tipo KSEQ
            read1 = pqm->pRead;
        }
        free(pqm);
        
        //mapRead((char *) seqR.sequence.c_str(), seqR.sequence.length());
        //mapRead(seqR);
        (this->*processingFunc)(read1);
        delete read1;
    }
}

void ogReadsMapper::processTheRead(ogSingleRead *seqR) {
    (this->*processingFunc)(seqR);
}


void ogReadsMapper::waste(char newEvent) {
    if (++wastingCycles == 0) wastingCycles2++;
    wastingEvents += newEvent;
}

//void ogReadsMapper::pushReadToQueue(char *pid, uint32_t pidLen, char *pseq, uint32_t pseqLen, char *pqual, uint32_t pqualLen, uint64_t readNum, char pack) {
//    ogSingleRead *seqR = new ogSingleRead(pid, pidLen, pseq, pseqLen, pqual, pqualLen, readNum, pack);
//    seqR->packIfNeeded(pRdKeyMapRd1->pGenome);
//    pushReadToQueue(seqR);
//}


/**
void ogReadsMapper::pushReadToQueue(kstring_t *name1, kstring_t *qual1, kstring_t *seq1, kstring_t *name2, kstring_t *qual2, kstring_t *seq2, uint16_t *pThreader, uint64_t readNum) {
    ogQueueMapper *pqm = (ogQueueMapper *) malloc(sizeof(ogQueueMapper));
    pqm->pFAS = (ogFastAQGZreader *) malloc(sizeof(ogFastAQGZreader));
    pqm->pFAS->name.l = name1->l;
    pqm->pSeq1 = seq1;
    pqm->pSeq2 = seq2;
    pqm->pQual1 = qual1;
    pqm->pQual2 = qual2;
    pqm->pName1 = name1;
    pqm->pName2 = name2;
    pqm->pThreadUsing = pThreader;
    pqm->readNum = readNum;
    mtx.lock();
        qSeq.push(pqm);
        //qSeq.push(theRead); // seqR will be "deleted" (as object) in this class
        //fprintf(stderr, "Queue %p size %lu\n", &qSeq, qSeq.size());
    mtx.unlock();
}
 **/

void ogReadsMapper::pushReadToQueue(ogSingleRead *pRead, uint64_t readNum) {
    ogQueueMapper *pqm = (ogQueueMapper *) malloc(sizeof(ogQueueMapper));
    pqm->pRead = pRead;
    pqm->pFAS1 = pqm->pFAS2 = NULL;
    pqm->readNum = readNum;
    mtx.lock();
        qSeq.push(pqm);
        //qSeq.push(theRead); // seqR will be "deleted" (as object) in this class
        //fprintf(stderr, "Queue %p size %lu\n", &qSeq, qSeq.size());
    mtx.unlock();
}

void ogReadsMapper::pushReadToQueue(ogFastAQsequence *pFAS1, ogFastAQsequence *pFAS2, uint64_t readNum) {
    ogQueueMapper *pqm = (ogQueueMapper *) malloc(sizeof(ogQueueMapper));
    pqm->pFAS1 = pFAS1;
    pqm->pFAS2 = pFAS2;
    pqm->pRead = NULL;
    pqm->readNum = readNum;
    mtx.lock();
        qSeq.push(pqm);
        //qSeq.push(theRead); // seqR will be "deleted" (as object) in this class
        //fprintf(stderr, "Queue %p size %lu\n", &qSeq, qSeq.size());
    mtx.unlock();
}

char ogReadsMapper::isQueueFull(char secureCheck) {
    //if (qSeq.size() >= maxQsize) { fprintf(stderr, "Q"); }
    if (secureCheck) {
        mtx.lock();
        char res = (qSeq.size() >= maxQsize);
        mtx.unlock();
        return res;
    } else {
        return (qSeq.size() >= maxQsize);        
    }
}

char ogReadsMapper::isQueueEmpty() {
    return (qSeq.size() == 0);
}

char ogReadsMapper::isQueueCritic() {
    return (qSeq.size() >= (maxQsize >> 2));
}

void ogReadsMapper::setMaxQsize(uint16_t mxq) {
    maxQsize = mxq;
}

void ogReadsMapper::setProductionMode(char isProduction) {
    productionMode = isProduction;
}

void ogReadsMapper::setSAMfromRead(ogSAM *pSAM, ogSingleRead *read) {
    pSAM->qname = read->pId;
    pSAM->seq = read->pSeq;
    pSAM->qual = read->pQual;
    pSAM->cigar = read->pCigar;
}

char * ogReadsMapper::revertCIGAR(char *pCigar, char *pDest) {
    uint16_t len = strlen(pCigar);
    char *pLast, *pI, *pCpy, *pStart;
    // Input:  140M2I8M13D10M  (len=14)
    // Output: 10M13D8M2I140M
    pCpy = pDest;
    pLast = pCigar + len;
    pI = pLast - 1;
    while (pI >= pCigar) {
        if (*pI >= 'A') {
            pStart = pI + 1;
            while (pLast > pStart) {
                *pCpy++ = *pStart++;
            }
            pLast = pI + 1;
        } else if (pI == pCigar) {
            while (pLast > pI) {
                *pCpy++ = *pI++;
            }
            break;
        }
        pI--;
    }
    *pCpy = 0;
    return pDest;
}


void ogReadsMapper::cutCIGARtoSeqLen(char *pCigar, uint32_t maxLen) {
    uint32_t cigLen = 0;
    uint32_t lenToUpdate = 0;
    uint32_t l = 0;
    char cigarOpToUpdate = '-';
    char *pC = pCigar;
    char *pA = pCigar;
    char *pT = NULL; // target
    if (*pC == '*') return;
    
    for(; *pC != 0; pC++) {
        if (*pC > '9' || *pC < '0') {
            if (*pC == 'M' || *pC == 'I' || *pC == 'S' || *pC == '=' || *pC == 'X') {
                if ((cigLen += l) > maxLen && pT == NULL) {
                    cigarOpToUpdate = *pC;
                    lenToUpdate = l - (cigLen - maxLen);
                    pT = pA;
                }
            }
            l = 0;
            pA = pC + 1;
        } else {
            l = l * 10 + *pC - 48;
        }
    }
    if (pT != NULL) {
        if (cigarOpToUpdate == 'S' || lenToUpdate == 0) snprintf(pT, pC-pT+10, "%uS", lenToUpdate + cigLen-maxLen);
        else snprintf(pT, pC-pT+10, "%u%c%uS", lenToUpdate, cigarOpToUpdate, cigLen-maxLen);
    }
}


void ogReadsMapper::correctLeftRighPositionFromCIGAR(char *pCigar, uint16_t readPos, uint32_t *pLeft, uint32_t *pRight, char isRd2) {
    uint32_t cigLen = 0;
    uint32_t l = 0;
    char *pC = pCigar;
    if (*pC == '*') return;
    
    for(; *pC != 0 && cigLen < readPos; pC++) {
        if (*pC > '9' || *pC < '0') {
            if (*pC == 'M' || *pC == 'I' || *pC == 'S' || *pC == '=' || *pC == 'X' || *pC == 'D') {
                //if (isRd2) {
                //    if (*pC == 'I') *pRight -= l; else 
                //    if (*pC == 'D') *pRight += l; else
                //    if (*pC == 'S') *pRight -= l;
                //} else {
                    if (*pC == 'I') *pLeft += l; else 
                    if (*pC == 'D') *pLeft -= l; else
                    if (*pC == 'S') *pLeft += l;
                //}
                cigLen += l;
                //if ((cigLen += l) >= readPos) break;
            }
            l = 0;
        } else {
            l = l * 10 + *pC - 48;
        }
    }
    l = 0;
    for(; *pC != 0; pC++) {
        if (*pC > '9' || *pC < '0') {
            if (*pC == 'M' || *pC == 'I' || *pC == 'S' || *pC == '=' || *pC == 'X' || *pC == 'D') {
                cigLen += l;
                //if (isRd2) {
                //    if (*pC == 'I') *pLeft += l; else 
                //    if (*pC == 'D') *pLeft -= l; else
                //    if (*pC == 'S') *pLeft += l;
                //} else {
                    if (*pC == 'I') *pRight -= l; else 
                    if (*pC == 'D') *pRight += l; else
                    if (*pC == 'S') *pRight -= l;
                //}
            }
            l = 0;
        } else {
            l = l * 10 + *pC - 48;
        }
    }
}
