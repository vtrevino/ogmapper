/* 
 * File:   ogAligner.cpp
 * Author: victortrevino
 * 
 * Created on July 15, 2022, 12:04 AM
 */

//#include <stdio>
//#include <string>
#include "ogAligner.hpp"
#include "ogDefinitions.hpp"
//#include "sequenceRead.h"
#include "ogSingleRead.hpp"
//#include "edlib.h"

//static uint64_t iProc = 0;    

ogAligner::ogAligner(ogGenome *ptrGenome) {
    pGenome = ptrGenome;
    allocated = 0; 
    fullAligments = 0;
    pMatScore = pLeftMax = pUpMax = NULL;
    pMatDecision = NULL;
    //allocate(1000);
    //pWFAligner = new WFAlignerGapAffine(4,6,2,WFAligner::Alignment,WFAligner::MemoryHigh);
    //pWFAligner = new WFAlignerGapAffine2Pieces(4,6,2,12,1,WFAligner::Alignment,WFAligner::MemoryHigh);
    pWFAligner = new WFAlignerGapAffine2Pieces(6,4,2,12,1,WFAligner::Alignment,WFAligner::MemoryHigh);
    pCigarOperations = new ogCigarOperations();
    //
    //dnaAlloc = 0;
    //DNA = NULL;
    //allocDNA(10000);
}

ogAligner::~ogAligner() {
    if (allocated > 0) free(pMatScore);
    //free(pLeftMax);
    //free(pMatDecision);
    // Free
    //if (dnaAlloc > 0) free(DNA);
    delete pCigarOperations;
    delete pWFAligner;
}

void ogAligner::allocate(uint16_t maxLen) {
    if (maxLen > allocated) {
        allocated = maxLen + 10;
        uint64_t oneMat = allocated * allocated;
        uint64_t oneVec = allocated;
        pMatScore    = (int16_t *) realloc(pMatScore, sizeof(int16_t) * (oneMat + 4 * oneVec) + sizeof(int8_t) * oneMat);
        pLeftMax     = pMatScore + oneMat;
        pLeftPos     = pLeftMax  + oneVec;
        pUpMax       = pLeftPos  + oneVec;
        pUpPos       = pUpMax    + oneVec;
        pMatDecision = (int8_t *) (pUpPos + oneVec);
    }
}
//
//void ogAligner::allocDNA(uint32_t len) {
//    if (len > dnaAlloc) {
//        dnaAlloc    = len * 2;
//        DNA         = (char *) realloc(DNA, sizeof(char) * dnaAlloc);
//    }    
//}

uint32_t ogAligner::align(ogSingleRead *read, ogReadKeyMapping *pMap, char checkKeySeq) {
    //return completeAlignment(read, pMap, checkKeySeq);
    return fastAlignment(read, pMap, checkKeySeq);
//    uint32_t res = simpleAlign(pSeq, seqLen, pMap, 0, 0);
//    if (res == 0) {
//        //doAlign(pSeq, seqLen, pMap, 1);
//        //fprintf(stderr, "^^^^^^^ Aligning: Id=%llu, Len=%u\n>%s\n", iProc, seqLen, pSeq);
//        //res = simpleAlign(pSeq, seqLen, pMap, 1, 0);
//        res = complexAlign(pSeq, seqLen, pMap);
//    }
//    return res;
}


uint32_t ogAligner::fastAlignment(ogSingleRead *read, ogReadKeyMapping *pMap, char checkKeySeq) {

    pMap->allocDNA(read->lenSeq + (read->isPaired ? read->pairedRead->lenSeq : 0)); //seqRead->seqLen1+seqRead->seqLen2);
    
    ogGenome    *pGenome = pMap->pGenome;
    ogCandidatePosManager *pCandPos = pMap->pCandPosMan;
    ogGenomeAndReadPosition *pGenRdPos;
    char        *pSeq = pMap->read->pSeq; 
    char        *pSeqRev = pMap->read->pSeqRevComp;
    uint32_t    seqLen = pMap->read->lenSeq; // pMap->seqLen; //seqRead->seqLen1;
    uint32_t    extLen = seqLen + NT_EXTENSION + NT_EXTENSION;
    uint32_t    nGenPos = pCandPos->getCount();
    char        maxScore = 0;
    uint32_t    maxScorePos = 0;
    uint32_t    maxScoreTies = 0;
    char        cycle, statusCycle;
    char        *dna = pMap->DNA;
    char        *pCIGAR;
    uint32_t    errPack;
    int         cigarLen = 500;
    char        isRevComp;
    uint32_t    Gpos, leftGPOS, pos;
    uint32_t    i;
    char        **pPacked;

    ///////////////////////////////////////////////
    /* PARAMETERS */
    uint32_t     maxErrAcceptPacked = MAX_ERROR_ACCEPT1_PACKED(seqLen);          // In Cycle 1
    uint32_t     maxErrAcceptAlign = MAX_ERROR_ACCEPT1_ALIGN(seqLen);            // In Cycle 1
    uint32_t     maxErrCheck2 = MAX_ERROR_CHECK2_ALIGN(seqLen);
    uint32_t     maxErrCheck3 = MAX_ERROR_CHECK3_ALIGN(seqLen);
    uint16_t     accepted = 0;
    uint16_t     cigarsMade = 0;
    ///////////////////////////////////////////////

    
    for (i=0; i < nGenPos; i++) {
        pGenRdPos = pCandPos->getkPos(i);
        Gpos = pGenRdPos->genomePosition;       //pCandPos->getkPos(i);
        isRevComp = pGenRdPos->isReverse;       //pCandPos->getIsReverseComplement(i);
        pos = pGenRdPos->readPosition;          //pMap->getRelativePos(isRevComp);
        leftGPOS = Gpos - (pos > Gpos ? 0 : pos);
        //pPacked = (isRevComp ? pMap->read->packedReadSeqRev : pMap->read->packedReadSeqFwd);
        //errPack = pGenome->comparePackedSequences(pPacked, leftGPOS, seqLen, maxErrCheck3);
        errPack = pGenome->comparePackedSequences(pMap->read->packReadSeqIfNeeded(pGenome->getPackedOffset(leftGPOS), isRevComp), leftGPOS, seqLen, maxErrCheck3);
        pGenRdPos->score = (seqLen-errPack)*100/seqLen;
        if (errPack <= maxErrAcceptPacked) {
            pGenRdPos->status = 'A'; 
            accepted++;
        } else if (errPack <= maxErrAcceptAlign) {
            pGenome->extractFromGenomicPos(leftGPOS > NT_EXTENSION ? leftGPOS-NT_EXTENSION : 0, extLen, dna, 0, 0); // Always at reference sequence
            fullAligments++;
            // WFA2 (https://github.com/smarco/WFA2-lib)
            pWFAligner->alignEndsFree(dna, extLen, NT_EXTENSION, NT_EXTENSION, isRevComp ? pSeqRev : pSeq, seqLen, NT_EXTENSION, NT_EXTENSION);
            pWFAligner->getAlignmentCigar(&pCIGAR, &cigarLen);
            pGenRdPos->status = 'G'; 
            pGenRdPos->score = buildCIGARfromWFA(pCIGAR, read, leftGPOS, leftGPOS+seqLen); //(seqLen-errPack)*100/seqLen;
            accepted++;
            cigarsMade++;
        } else if (errPack <= maxErrCheck2) {
            pGenRdPos->status = '2';
        } else if (errPack <= maxErrCheck3) {
            pGenRdPos->status = '3';
        } else {
            pGenRdPos->status = 'g';
        }
        if (pGenRdPos->score >= maxScore) {
            if (pGenRdPos->score > maxScore) {
                maxScore = pGenRdPos->score;
                maxScorePos = i;
                maxScoreTies = 0;
            } else {
                maxScoreTies++;
            }
        }
    }
    for (cycle=2; accepted == 0 && cycle < 4; cycle++) {
        statusCycle = '0'+cycle;
        for (i=0; i < nGenPos; i++) {
            pGenRdPos = pCandPos->getkPos(i);
            if (pGenRdPos->status == statusCycle) {
                Gpos = pGenRdPos->genomePosition;       //pCandPos->getkPos(i);
                isRevComp = pGenRdPos->isReverse;       //pCandPos->getIsReverseComplement(i);
                pos = pGenRdPos->readPosition;          //pMap->getRelativePos(isRevComp);
                leftGPOS = Gpos - (pos > Gpos ? 0 : pos);
                errPack = seqLen - (pGenRdPos->score * seqLen / 100);
                pGenome->extractFromGenomicPos(leftGPOS > NT_EXTENSION ? leftGPOS-NT_EXTENSION : 0, extLen, dna, 0, 0); // Always at reference sequence
                fullAligments++;
                // WFA2 (https://github.com/smarco/WFA2-lib)
                pWFAligner->alignEndsFree(dna, extLen, NT_EXTENSION, NT_EXTENSION, isRevComp ? pSeqRev : pSeq, seqLen, NT_EXTENSION, NT_EXTENSION);
                pWFAligner->getAlignmentCigar(&pCIGAR, &cigarLen);
                pGenRdPos->score = buildCIGARfromWFA(pCIGAR, read, leftGPOS, leftGPOS+seqLen); //(seqLen-errPack)*100/seqLen;
                cigarsMade++;
                if (pGenRdPos->score >= maxScore) {
                    if (pGenRdPos->score > maxScore) {
                        maxScore = pGenRdPos->score;
                        maxScorePos = i;
                        maxScoreTies = 0;
                    } else {
                        maxScoreTies++;
                    }
                }
                errPack = seqLen - (pGenRdPos->score * seqLen / 100);
                if (errPack <= maxErrCheck3) { // maxErrAcceptAlign
                    pGenRdPos->status = 'G';
                    accepted++;
                }
            }
        }
    }
    pMap->maxScore = maxScore;
    pMap->maxScorePos = maxScorePos;
    pMap->maxTies = maxScoreTies;
    
    if (accepted > 1 || cigarsMade > 1) *read->pCigar = 0; // delete cigar becuse it can be wrong, lets recalculate when writing sam
    *read->pCigar = 0;

    return accepted;
}


uint32_t ogAligner::completeAlignment(ogSingleRead *read, ogReadKeyMapping *pMap, char checkKeySeq) {
    //static EdlibEqualityPair    edlEq[4] = {{'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}};
    
    pMap->allocDNA(read->lenSeq+read->pairedRead->lenSeq); //seqRead->seqLen1+seqRead->seqLen2);
    
    char *pSeq = pMap->read->pSeq; // pMap->seq; //seqRead->pSeq1;
    uint32_t seqLen = pMap->read->lenSeq; // pMap->seqLen; //seqRead->seqLen1;

    ogCandidatePosManager *pCandPos = pMap->pCandPosMan;
    ogGenomeAndReadPosition *pGenRdPos;
    uint32_t    nGenPos = pCandPos->getCount();
    uint32_t    n,eK,eL,eR,e,nK, nR, nL, errKey, errSeq;
    char        *dna = pMap->DNA;
    char        cCigar;
    char        *cigarLet;          // cigarLet[500]
    int         cigarLetLen;
    char        *pCIGAR, *pReadSeq, *pDNA;
    int         cigarLen = 500;
    char        isRevComp;
    uint16_t    pos, posDNA;
    ogChromosome *pChr;
    uint32_t    res = 0;
    uint32_t    Gpos, expectedGPos, closerGenPos;
    uint32_t    i,j,L,k,nKeys, keySize = pMap->pKeyEncoding->getSizeInChars();
    uint32_t    lastAlignment_i = 0;
    uint32_t    pI;
    ogKey       *pK;
    uint32_t    minAlignedPos, maxAlignedPos;
    uint32_t    leftAddNt, rightAddNt;
    uint32_t    rightExactPos, leftExactPos;
    char        *pSeqN, *pDnaN, flag, fine, score, status, *pC;
    float       ovER; // 1/16
    uint32_t    M,X,I,D,matches, nomatches, ins, del, gapsOpen, maxScore=0, maxScorePos=0, maxTies=0, nosoft, soft, ms, x, xm, m;

    ///////////////////////////////////////////////
    /* PARAMETERS */
    uint32_t     maxKeyErr = MAX_KEY_ERR(keySize);
    uint32_t     okSimpleAlignThreshold = MAX_SIMPLE_ALIGN_ERR(seqLen);
    uint32_t     okEndSimpleAlignThreshold = MIN_SIMPLE_ALIGN_ERR(seqLen);
    ///////////////////////////////////////////////

    //fprintf(stderr, "RdNum=%llu\n", pMap->seqNum);
    //if (pMap->seqNum == 4388) {
    //    fprintf(stderr, "----------------Found---------------\n");
    //}
    /**
    fprintf(stderr, "|||| Aligning %llu : %s : %u targets\n", alignNum, pSeq, nGenPos);
    for (i=0; i < nGenPos; i++) {
        garp = pCandPos->getkPos(i);
        pChr = pMap->pGenome->getGenomicCoordinate(garp->genomePosition);
        fprintf(stderr, "Rd Pos = %u, Gen Pos = %u:%u, RC = %u, Status = %c, Score = %u\n",garp->readPosition, pChr->number, garp->genomePosition - pChr->start, (garp->isReverse ? 1 : 0), garp->status, garp->score);
    }
     **/
    
    cigarLet = read->pCigar; //(pMap->isRead2 ? seqRead->pCigar2 : seqRead->pCigar1);
    
    for (i=0; i < nGenPos; i++) {
        pGenRdPos = pCandPos->getkPos(i);
        Gpos = pGenRdPos->genomePosition;       //pCandPos->getkPos(i);
        isRevComp = pGenRdPos->isReverse;       //pCandPos->getIsReverseComplement(i);
        pos = pGenRdPos->readPosition;          //pMap->getRelativePos(isRevComp);
        L = seqLen - pos;
        if (pos > Gpos) pos = Gpos;             // can happen only in telomere sequences
        
        //fprintf(stderr, "RdNum=%llu | Cand=%u, RC=%c, St=%c, Sc=%hhu, |Gpos|=%u, GenPos=%u, ReadPos=%u\n", seqRead->readNumber, i, 48+isRevComp, pGenRdPos->status, pGenRdPos->score, Gpos-pos, Gpos, pos);

        fine = 1;
        errKey = 99;
        // Check Key First
        if (checkKeySeq) {
            pGenome->extractFromGenomicPos(Gpos, keySize, dna, isRevComp, 0);
            pReadSeq = pSeq + (isRevComp ? seqLen - pos - keySize : pos);
            e = 0;
            pDNA = dna;
            for (j=0; pReadSeq && pDNA && j < keySize; j++) {
                if ((*pReadSeq++ & 0x5F) != *pDNA++) {
                    e++;
                }
            }
            pCandPos->setPositionKScore(i, errKey = (keySize - e) * 100 / keySize);
            if (e <= maxKeyErr) {
                pCandPos->setPositionKStatus(i, 'K');
            } else {
                pCandPos->setPositionKStatus(i, 'k');
                fine = 0;
            }
        }
        
        if (fine) {
            // Key is fine
            // Now SIMPLE ALIGN
            pGenome->extractFromGenomicPos(Gpos-pos, seqLen, dna, isRevComp, 0);
            errSeq = 0;
            for (n=0,pSeqN=pSeq,pDnaN=dna; *pSeqN && *pDnaN && n < seqLen; n++, pSeqN++, pDnaN++) {
                if ((*pSeqN & 0x5F) != *pDnaN) {
                    errSeq++;
                }
            }
            pCandPos->setPositionKScore(i, score = (seqLen - e) * 100 / seqLen);
            if (errSeq < okSimpleAlignThreshold) {
                fine = 0;   // do not continue for long alignment
                res++;      // but this is a OK result
                lastAlignment_i = i;
                pCandPos->setPositionKStatus(i, 'A'); // Simple alignment is ok
            } else if (errSeq <= okEndSimpleAlignThreshold) {
                pCandPos->setPositionKStatus(i, 'a'); // Simple Align not passed
                fine = 0;
            } else {
                // Must be ALIGNED COMPLETLY
            }
            
            
            fine = 1; // force all alignments.
            /**
            if (false) {
                seqRead->samRecord1.qname = seqRead->pId1;
                seqRead->samRecord1.flag = SAMFLAG_READ1;
                if (isRevComp) seqRead->samRecord1.flag |= SAMFLAG_READ_REVCOMP;
                pChr = pGenome->getGenomicCoordinate(Gpos-pos);
                seqRead->samRecord1.rname = pChr->name;
                seqRead->samRecord1.pos   = Gpos-pos - pChr->start;
                seqRead->samRecord1.mapq  = score; //-10 * log10 (probability is correct) ... instead, we will write % of nt errors vs ref
                cigarLetLen = sprintf(cigarLet, "%uM", seqLen);
                cigarLet[cigarLetLen] = 0;
                seqRead->samRecord1.cigar = cigarLet;
                seqRead->samRecord1.rnext = NULL;
                seqRead->samRecord1.pnext = 0;
                seqRead->samRecord1.tlen = seqLen;
                seqRead->samRecord1.seq  = (isRevComp ? seqRead->pSeqRevComp1 : seqRead->pSeq1);
                seqRead->samRecord1.qual = (isRevComp ? seqRead->pQualRev1 : seqRead->pQual1);
                pMap->pSamWriter->writeSAMRecord(&(seqRead->samRecord1));
            }
             **/
        }

        if (fine) {
            //fprintf(stderr, "Performing Complete Alignment=%llu\n", fullAligments); fflush(stderr);
            // Simple alignment was unable to define mapping
            // PERFORM COMPLETE ALIGNMENT
            
            leftAddNt = rightAddNt = 0;
            nKeys = pMap->getNKeyPos(isRevComp);
            // Find left key position closer to the Gpos
            minAlignedPos = pos;
            for (k=0; leftAddNt == 0 || leftAddNt > 10 && k < nKeys; k++) {
                pI = pMap->getKeyablePosition(isRevComp,k);
                if (pI < pos) {
                    expectedGPos = Gpos - (pos - pI);
                    pK = pMap->getInfoForPositionKey(isRevComp, k);
                    closerGenPos = pMap->getCloserGenomePosition(expectedGPos, pMap->pGenPos->getPointerPosition(pK->offset), pK->size); 
                    if (closerGenPos == expectedGPos) {
                        minAlignedPos = pI;
                        leftAddNt = 0;
                        break;
                    } else if (closerGenPos < expectedGPos && expectedGPos-closerGenPos < MAX_GAP_1READ_ALIGMENT) {
                        if ((X = expectedGPos-closerGenPos) < leftAddNt || leftAddNt == 0) {
                            minAlignedPos = pI;
                            leftAddNt = X;
                        }
                    } else if (closerGenPos > expectedGPos && closerGenPos-expectedGPos < MAX_GAP_1READ_ALIGMENT) {
                        if ((X = closerGenPos-expectedGPos) < leftAddNt || leftAddNt == 0) {
                            minAlignedPos = pI;
                            leftAddNt = X;
                        }
                    }
                } else {
                    break;
                }
            }
            // Find right key position closer to the Gpos
            maxAlignedPos = pos;
            for (k=nKeys-1; rightAddNt == 0 || rightAddNt > 10 && k > 0; k--) {
                pI = pMap->getKeyablePosition(isRevComp,k);
                if (pI > pos) {
                    expectedGPos = Gpos - (pos - pI);
                    pK = pMap->getInfoForPositionKey(isRevComp, k);
                    closerGenPos = pMap->getCloserGenomePosition(expectedGPos, pMap->pGenPos->getPointerPosition(pK->offset), pK->size); 
                    if (closerGenPos == expectedGPos) {
                        maxAlignedPos = pI;
                        rightAddNt = 0;
                        break;
                    } else if (closerGenPos < expectedGPos && expectedGPos-closerGenPos < MAX_GAP_1READ_ALIGMENT) {
                        if ((X = expectedGPos-closerGenPos) < rightAddNt || rightAddNt == 0) {
                            maxAlignedPos = pI;
                            rightAddNt = X;
                        }
                    } else if (closerGenPos > expectedGPos && closerGenPos-expectedGPos < MAX_GAP_1READ_ALIGMENT) {
                        if ((X = closerGenPos-expectedGPos) < rightAddNt || rightAddNt == 0) {
                            maxAlignedPos = pI;
                            rightAddNt = X;
                        }
                    }
                } else {
                    break;
                }
            }

            if (leftAddNt != 0 || rightAddNt != 0) {
                pGenome->extractFromGenomicPos(Gpos-pos-leftAddNt, seqLen+leftAddNt+rightAddNt, dna, isRevComp, 0);
            }
            // Error in Key
            if (isRevComp) {
               pos = seqLen-pos-keySize;
               posDNA = pos + rightAddNt;
            } else {
                posDNA = pos + leftAddNt;
            }
            for (e=0,n=0,pSeqN=pSeq+pos,pDnaN=dna+posDNA; *pSeqN && *pDnaN && n < keySize; n++, pSeqN++, pDnaN++) {
                if ((*pSeqN & 0x5F) != *pDnaN) {
                    e++;
                }
            }
            eK = e;
            nK = n;
            if (eK > maxKeyErr) {
                // error key not aligned not supposed to happen
                //**
                fprintf(stderr, "***************************************\n");
                fprintf(stderr, "** Error in Key , eK=%u, eStatus=%u\n", eK, errKey);
                fprintf(stderr, "***************************************\n");
                fprintf(stderr, "RC:%c, pos=%u, Add[L=%u, R=%u], ReadNum:%llu\n",isRevComp ? '1' : '0', pos, leftAddNt, rightAddNt, pMap->read->readNumber);
                fprintf(stderr, "Key Start in Read:%s\n",pSeq+pos);
                fprintf(stderr, "Key Start in DNA :%s\n", dna+posDNA);
                // **/
            }
            // Error to the Right (Extending)
            for (rightExactPos=pos+n+leftAddNt, flag=1; (float)e/n <= MAX_ERROR_RATE_WFA && *pSeqN && *pDnaN; pSeqN++, pDnaN++, n++) {
                if ((*pSeqN & 0x5F) == *pDnaN) {
                    if (flag) rightExactPos++;
                } else {
                    e++;
                    flag=0;
                }
            }
            eR = (e-eK);
            nR = (n-nK);
            e  = eK;
            n  = nK;
            // Error to the Left (Extending)
            for (leftExactPos=pos+rightAddNt, flag=1, j=pos, pSeqN=pSeq+pos-1, pDnaN=dna+pos-1; 
                j > 0 && (float)e/n <= MAX_ERROR_RATE_WFA && *pSeqN && *pDnaN; 
                pSeqN--, pDnaN--, n++, j--) {
                if ((*pSeqN & 0x5F) == *pDnaN) {
                    if (flag) leftExactPos--;
                } else {
                    e++;
                    flag=0;
                }
            }
            eL = (e - eK);
            nL = (n - nK);
            ovER = (float) (eL+eK+eR) / (nL+nK+nR);
            //if (ovER < maxErrRate) res++;
            /**
            fprintf(stderr, "%u/%u:%c, Err[L=%u, K=%u, R=%u, Rate=%3.1f%%], Len Chck=%u/%u+(%u), qExact [L=%d, R=%d, Size=%u, from=%u, to=%u], Extended [L=%u, R=%u]\n", 
                    i, nGenPos, pCandPos->getPositionKStatus(i), eL, eK, eR, 100*ovER,  // ...ErrRate
                    nL + nK + nR, seqLen, leftAddNt+rightAddNt,
                    pos-leftExactPos, rightExactPos-pos-keySize, rightExactPos-leftExactPos, leftExactPos, rightExactPos, // ... qExact
                    leftAddNt, rightAddNt
                    );
            fprintf(stderr, "DNA: %s\n",dna);
             **/

            ++fullAligments;
            //fprintf(stderr, "FullAlignment=%llu\n", fullAligments); fflush(stderr);
            
            // WFA2 (https://github.com/smarco/WFA2-lib)
            pWFAligner->alignEndsFree(pSeq, seqLen, 0, 0, dna, seqLen+leftAddNt+rightAddNt, 0, 0); // Pattern = Read, Text = Reference Genome
            pWFAligner->getAlignmentCigar(&pCIGAR, &cigarLen);
            //fprintf(stderr, "CIGAR: %s. Aligment score: %d.\n", pCIGAR, pWFAligner->getAlignmentScore());            

            matches = nomatches = ins = del = gapsOpen = 0;
            pDnaN = pCIGAR;
            cigarLetLen = soft = 0;
            //fprintf(stderr, "CIGAR:");
            while ((cCigar = *pDnaN)) {
                for (M=X=0; cCigar == 'M' || cCigar == 'X'; ) {  // soft=nosoft=ms=0
                    //fprintf(stderr, "%c", cCigar);
                    if (cCigar == 'M') {
                        M++;
                    } else {
                        X++;
                    }
                    cCigar = *++pDnaN;
                }
                if (M+X > 0) {
                    if (cigarLetLen == 0 && X > 1) {
                        // Es al inicio y hay mismatches
                        if (M > 0) {
                            // Tambien hay matches ... Search desde el inicio hasta encontrar
                            pC = pCIGAR;
                            x = m = xm = 0;
                            while (*pC == 'M' | *pC == 'X') {
                                if (*pC == 'M') {
                                    if (++m > 4) break;
                                } else {
                                    x++;
                                    xm += m;
                                    m = 0;
                                }
                                pC++;
                            }
                            if (x > 1) { 
                                cigarLetLen += snprintf(cigarLet+cigarLetLen, MAX_CIGAR_SIZE - cigarLetLen, "%uS", x+xm); 
                                X -= x; 
                                M -= xm;
                            }
                        } else {
                            // Todos son mismatches ... Just add soft
                            cigarLetLen += snprintf(cigarLet+cigarLetLen, MAX_CIGAR_SIZE - cigarLetLen, "%uS", X); 
                            X = 0; 
                        }
                    }
                    if (cCigar == 0 && X > 1) {
                        pC = pDnaN - 1;
                        x = m = xm = 0;
                        while (*pC == 'M' | *pC == 'X') {
                            if (*pC == 'M') {
                                if (++m > 4) break;                                    
                            } else {
                                x++;
                                xm += m;
                                m = 0;
                            }
                            if (--pC <= pCIGAR) break;
                        }
                        if (x > 1) {
                            soft = x + xm;
                            X -= x;
                            M -= xm;
                        }
                    }
                    cigarLetLen += snprintf(cigarLet+cigarLetLen, MAX_CIGAR_SIZE - cigarLetLen, "%uM", M+X); 
                    matches += M; 
                    nomatches += X;
                }
                for (I=0; cCigar == 'I'; I++) cCigar = *++pDnaN;
                if (I > 0) { gapsOpen++; cigarLetLen += snprintf(cigarLet+cigarLetLen, MAX_CIGAR_SIZE-cigarLetLen, "%uI", I); ins += I; }
                for (D=0; cCigar == 'D'; D++) cCigar = *++pDnaN;
                if (D > 0) { gapsOpen++; cigarLetLen += snprintf(cigarLet+cigarLetLen, MAX_CIGAR_SIZE-cigarLetLen, "%uD", D); del += D; }
            }
            if (soft > 0) cigarLetLen += snprintf(cigarLet+cigarLetLen, MAX_CIGAR_SIZE-cigarLetLen, "%uS", soft);
            cigarLet[cigarLetLen] = 0;
            ovER = ((float) (nomatches+(gapsOpen * 5)) / (float) (matches+nomatches));
            //fprintf(stderr, "CIGAR let:%s | M=%u,X=%u,D=%u,I=%u | Error=%f\n", cigarLet, matches, nomatches, del, ins, ovER);
            if ((matches+nomatches+del+ins) > seqLen*2) {
                static uint64_t large = 0;
                large++;
                if (false) {
                    fprintf(stderr, "\n[[[FullAlignNum:%llu, LargeNum:%llu, CIGAR let:%s | M=%u,X=%u,D=%u,I=%u,Gaps=%u\nAcc: Key=%u, Seq=%u, WFA=%f | Prev Score=%hhu | Left+=%u | Right+=%u\n", fullAligments, large, cigarLet, matches, nomatches, del, ins, gapsOpen, errKey, errSeq, 1-ovER, pCandPos->getPositionKScore(i), leftAddNt, rightAddNt);
                    pChr = pGenome->getGenomicCoordinate(Gpos); 
                    fprintf(stderr, "LARGE SEQ IN COMPLEX ALIGNMENT. GenPos=%u, Chr %u:%u, isRev=%hhu, pos=%u, RdNum=%llu\n%s\n", Gpos, pChr->number, Gpos - pChr->start, isRevComp, pos, pMap->read->readNumber, pSeq);
                    fprintf(stderr, "ALIGMENT: Score=%d, Status=%d, Nt x Nt:\n", pWFAligner->getAlignmentScore(), pWFAligner->getAlignmentStatus());
                    for (j=0; j < 3; j++) {
                        char *s = (j == 0 ? pSeq : (j == 1 ? pCIGAR : dna));
                        pDnaN = pCIGAR;
                        k = 0;
                        while ((cCigar = *pDnaN)) {
                            while (cCigar == 'M' || cCigar == 'X') { 
                                fprintf(stderr, "%c", *s++);
                                cCigar = *++pDnaN;
                            }
                            while (cCigar == 'I') { 
                                cCigar = *++pDnaN;
                                if (j == 0) fprintf(stderr, "-");
                                if (j != 0) fprintf(stderr, "%c", *s++);
                            }
                            while (cCigar == 'D') {
                                cCigar = *++pDnaN;
                                if (j != 2) fprintf(stderr, "%c", *s++);
                                if (j == 2) fprintf(stderr, "-");
                            }
                        }
                        fprintf(stderr, "\n");
                    }
                    fprintf(stderr, "]]]");
                }
            }
            if (ovER < MAX_ERROR_RATE_WFA) {
                res++;
                lastAlignment_i = i;
                pCandPos->setPositionKStatus(i, 'G'); // Global Alignment OK!
            } else {
                pCandPos->setPositionKStatus(i, 'g'); // Global Alignment not working
            }
            pCandPos->setPositionKScore(i, 100*(1-ovER));

        }
        status = pCandPos->getPositionKStatus(i);
        score = pCandPos->getPositionKScore(i);
        if (status == 'A' || status == 'G') {
            if (score >= maxScore) {
                if (score == maxScore) {
                    maxTies++;
                } else {
                    maxScore = score;
                    maxScorePos = i;
                    maxTies = 0;
                }
            }
        }
    }
    //fprintf(stderr, "RdNum=%llu\n", pMap->seqNum);
    pMap->maxScore = maxScore;
    pMap->maxScorePos = maxScorePos;
    pMap->maxTies = maxTies;
    
    return res;
}

// SE ELIMINO ESTO AL CAMBIAR CIGAR LeftType por cigarLeftPos y cigarRightPos
//uint32_t ogAligner::compareSeqAndDNABuildingCIGAR(char *pSeq, char *dna, ogReadKeyMapping *pRKM, ogSingleRead *read) {
//    // uint32_t seqLen, char *pCIGAR, uint16_t *pLeftSoft
//    uint32_t seqLen = read->lenSeq;
//    char *pCIGAR = read->pCigar;
//    
//    // Esta rutina no contempla Ins ni Del supone que es toda la region
//    // X,M conteos globals // m conteo continuos // x conteo semi-continuo
//    uint32_t n=0, X=0, M=0, m=0, x=0, cigarLen = 0, soft = 0, softStart = 1;
//    char *pSeqN, *pDnaN;
//    //fprintf(stderr, "\nSEQ  :%s\n",pSeq);
//    //fprintf(stderr, "DNA  :%s\n",dna);
//    //fprintf(stderr, "CIGAR:");
//    read->cigarLeft = 0;
//    read->cigarLeftType = '-';
//    for (pSeqN=pSeq,pDnaN=dna; *pSeqN && *pDnaN && n < seqLen; n++, pSeqN++, pDnaN++) {
//        if ((*pSeqN & 0x5F) != *pDnaN) {
//            //fprintf(stderr, "X");
//            X++;
//            if (x > 0 || m < 4) x += m;
//            x++;
//            m = 0;
//        } else {
//            //fprintf(stderr, "M");
//            M++;
//            if (++m > 4) {
//                if (softStart && x > 1) {
////                    if (x == 126) {
////                        fprintf(stderr, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>11111111\n%s\n%s\nx=%u,m=%u,X=%u,M=%u,n1=%ld,n2=%ld,thread=%u,dna*=%p,maxScorePos=%u\n",pSeq,dna,x,m,X,M,pSeqN-pSeq,pDnaN-dna,pRKM->thread,dna,pRKM->maxScorePos);
////                    }
//                    cigarLen += snprintf(pCIGAR + cigarLen, MAX_CIGAR_SIZE-cigarLen, "%uS", x);
//                    soft += x;
//                    read->cigarLeft = x;
//                    read->cigarLeftType = 'S';
//                }
//                softStart = 0;
//                x = 0;
//            }
//        }
//    }
//    if (x > 1) {
//        x += m;
//        m = 0;
//        if (softStart) {
//            // Very rare
////            if (x == 126) {
////                fprintf(stderr, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>22222222\n%s\n%s\nx=%u,m=%u,X=%u,M=%u,n1=%ld,n2=%ld,thread=%u,dna*=%p,maxScorePos=%u,isRd2=%c\n",pSeq,dna,x,m,X,M,pSeqN-pSeq,pDnaN-dna,pRKM->thread,dna,pRKM->maxScorePos,'0'+pRKM->isRead2);
////            }
//            cigarLen += snprintf(pCIGAR + cigarLen, MAX_CIGAR_SIZE-cigarLen, "%uS", x);
//            read->cigarLeft = x;
//            read->cigarLeftType = 'S';
//            soft += x;
//            x = 0;
//        }
//    } else {
//        x = 0;
//    }
//    if (M+X > soft+x) {
//        cigarLen += snprintf(pCIGAR + cigarLen, MAX_CIGAR_SIZE-cigarLen, "%uM", M+X-soft-x);
//    }
//    if (x > 1) {
//        // Soft at the end
////        if (x == 126) {
////            fprintf(stderr, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>333333\n%s\n%s\nx=%u,m=%u,X=%u,M=%u,n1=%ld,n2=%ld,thread=%u,dna*=%p,maxScorePos=%u\n",pSeq,dna,x,m,X,M,pSeqN-pSeq,pDnaN-dna,pRKM->thread,dna,pRKM->maxScorePos);
////        }
//        cigarLen += snprintf(pCIGAR + cigarLen, MAX_CIGAR_SIZE-cigarLen, "%uS", x);
//        soft += x;
//    }
//    pCIGAR[cigarLen] = 0;
//    read->cigarMatches = M;
//    read->cigarNoMatches = X;
//    read->cigarIns = 0;
//    read->cigarDel = 0;
//    read->cigarSoft = soft + (read->cigarLeftType == 'S' ? read->cigarLeft : 0);
//    //fprintf(stderr, ":%s:%uX\n",pCIGAR,X);
//    return X; // X = mismatches totales
//}



uint32_t ogAligner::checkKeySeq(char *pSeq, uint32_t seqLen, ogReadKeyMapping *pMap) {
    
    pMap->allocDNA(seqLen);
    
    ogCandidatePosManager *pCandPos = pMap->pCandPosMan;
    ogGenomeAndReadPosition *pGenRdPos;
    uint32_t    nGenPos = pCandPos->getCount();
    uint32_t    i, j, k, n;
    char        *dna = pMap->DNA;
    char        *pDNA;
    char        *pReadSeq;
    char        isRevComp;
    uint32_t    Gpos, pos, L;
    ogChromosome *pChr;
    uint32_t    keySize = pMap->pKeyEncoding->getSizeInChars();
    uint32_t    maxErr = keySize / 10;
    int32_t     e;
    uint32_t    nOk = 0;
    
    k = 0;
    for (i=0; i < nGenPos; i++) {
        pGenRdPos = pCandPos->getkPos(i);
        isRevComp = pGenRdPos->isReverse; //pCandPos->getIsReverseComplement(i);
        Gpos = pGenRdPos->genomePosition; //pCandPos->getkPos(i);
        //pChr = pGenome->getGenomicCoordinate(Gpos); 
        //pos  = pMap->getRelativePos(isRevComp);     // Relative to read
        pos = pGenRdPos->readPosition;
        pGenome->extractFromGenomicPos(Gpos, keySize, dna, isRevComp, 0);
        if (isRevComp) {
            //n = seqLen-pos-keySize;
            pReadSeq = pSeq + seqLen-pos-keySize;
        } else {
            //n = pos;
            pReadSeq = pSeq + pos;
        }
        e=0;
        pDNA = dna;
        for (j=0; pReadSeq && pDNA && j < keySize; j++) {
            if ((*pReadSeq++ & 0x5F) != *pDNA++) {
                e++;
            }
        }
        pCandPos->setPositionKScore(i, (keySize - e) * 100 / keySize);
        if (e <= maxErr) {
            nOk++;
            pCandPos->setPositionKStatus(i, 'K');
        } else {
            pCandPos->setPositionKStatus(i, 'k');
        }
        // printing
        // pChr = pGenome->getGenomicCoordinate(Gpos); 
        // fprintf(stderr, "==> %u/%u | RC=%c | e=%d | pos=%u | Chr=%u:%u | DNA=%s | RD[key]=%s | len=%u\n", i, nGenPos, (isRevComp ? '1' : '0'), e, pos, pChr->number, Gpos-pChr->start, dna, pSeq+(isRevComp ? seqLen-pos-keySize : pos), seqLen);
    }
    return nOk;
}

uint32_t ogAligner::simpleAlign(char *pSeq, uint32_t seqLen, ogReadKeyMapping *pMap, char show, char dynAlign) {
    
    pMap->allocDNA(seqLen*2);
    
    ogCandidatePosManager *pCandPos = pMap->pCandPosMan;
    ogGenomeAndReadPosition *pGenRdPos;
    uint32_t    nGenPos = pCandPos->getCount();
    int32_t     n,e;
    char        *dna = pMap->DNA;
    char        isRevComp;
    char        lowcase = (dynAlign ? 0 : 1);
    uint16_t    pos;
    ogChromosome *pChr;
    uint32_t    okSimpleAlignThreshold = seqLen * 5 / 100; //seqLen >> 4;
    uint32_t    okEndSimpleAlignThreshold = seqLen >> 1;
    uint32_t    res = 0;
    uint32_t    Gpos;
    uint32_t    i,j,L;
    char        *pSeqN, *pDnaN;

    if (show) { 
        fprintf(stderr, "/////// Aligning: SeqNum=%llu, Len=%u, nTargets=%u\n>%s\n",pMap->read->readNumber, seqLen, nGenPos, pSeq); fflush(stderr); 
    }
    for (i=0; i < nGenPos && i < 255; i++) {
        if (pCandPos->getPositionKStatus(i) == 'K') {
            pGenRdPos = pCandPos->getkPos(i);
            Gpos = pGenRdPos->genomePosition; //pCandPos->getkPos(i);
            isRevComp = pGenRdPos->isReverse; //pCandPos->getIsReverseComplement(i);
            pos = pGenRdPos->readPosition; //pMap->getRelativePos(isRevComp);
            L = (seqLen <= pos ? pMap->keySize : seqLen-pos+1);
            if (pos > Gpos) pos = Gpos;
            if (isRevComp) {
                pGenome->extractFromGenomicPos(Gpos-1, L, dna, 1, 0);
                pGenome->extractFromGenomicPos(Gpos-pos-1, pos, dna+L, 1, lowcase);
            } else {
                pGenome->extractFromGenomicPos(Gpos-pos, pos, dna, 0, lowcase);
                pGenome->extractFromGenomicPos(Gpos, L, dna+pos, 0, 0);
            }
            if (dynAlign) {
                // Esto va a haber que quitarlo para que se haga una rutina igual a esta pero especifica para globalAlignment checando el estatus correcto 'e' 
                e = globalAlignment(pSeq, seqLen, dna, seqLen);
            } else {
                if (show) fprintf(stderr, "_");
                e = 0;
                for (n=0,pSeqN=pSeq,pDnaN=dna; *pSeqN && *pDnaN && n < seqLen; n++, pSeqN++, pDnaN++) {
                    if ((*pSeqN & 0x5F) != *pDnaN) {
                        e++;
                        if (show) fprintf(stderr, "^");
        //                if (n <= left) el++;
        //                else if (n >= right) er++;
                    } else {
                        if (show) fprintf(stderr, " ");
                    }
                }
            }
            if (show) { 
                pChr = pGenome->getGenomicCoordinate(Gpos);
                fprintf(stderr, "\n~%s ", dna); fflush(stderr); 
                fprintf(stderr, "i=%u,e=%d,kPos=%u,3'=%hhu,Chr=%d:%u\n", i+1,  e, Gpos, isRevComp, pChr->number, Gpos-pChr->start); fflush(stderr); 
            }
            pCandPos->setPositionKScore(i, (seqLen - e) * 100 / seqLen);
            if (e < okSimpleAlignThreshold) {
                res++;
                pCandPos->setPositionKStatus(i, 'A'); // Simple alignment is ok (e < 1/16)
                //if (e == 0) break;
            } else if (e <= okEndSimpleAlignThreshold) {
                pCandPos->setPositionKStatus(i, 'e'); // Needs extension (e < 1/4)
                //res++;
    //            // soft clipping en los extremos???
    //            fprintf(stderr, "#### ERRORS AT ENDS: %u:%d:%.2f%%:Left=%u:Right=%u:MxErr=%u\n", n, e, (float) 100*e/n, el, er, okSimpleAlignThreshold);
    //            fprintf(stderr, "R>%s\nG>%s %u,%d,%u,3'=%hhu,Chr=%d:%u\n", pSeq, dna, i+1,  e, pCandPos->getkPos(i), isRevComp, pChr->number, pCandPos->getkPos(i)-pChr->start);
    //            res++;
            } else {
                pCandPos->setPositionKStatus(i, 'a'); // Simple Align not passed
    //            fprintf(stderr, "%u:%d:%.2f%%\n", n, e, (float) 100*e/n);
    //            if ((float) e/n > 0.7) {
    //                fprintf(stderr, "R>%s\nG>%s %u,%d,%u,3'=%hhu,Chr=%d:%u\n", pSeq, dna, i+1,  e, pCandPos->getkPos(i), isRevComp, pChr->number, pCandPos->getkPos(i)-pChr->start);
    //            }
            }
        }
    }
    
    return res;
}


uint32_t ogAligner::complexAlign(char *pSeq, uint32_t seqLen, ogReadKeyMapping *pMap) {
    static uint32_t             MAX_GAP  = 1000;
    static uint64_t             alignNum = 0;
    //static EdlibEqualityPair    edlEq[4] = {{'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}};
    
    pMap->allocDNA(seqLen*2);

    ogCandidatePosManager *pCandPos = pMap->pCandPosMan;
    ogGenomeAndReadPosition *pGenRdPos;
    uint32_t    nGenPos = pCandPos->getCount();
    int32_t     n,eK,eL,eR,e,nK, nR, nL;
    char        *dna = pMap->DNA;
    char        cCigar;
    char        cigarLet[500];
    int         cigarLetLen;
    char        *pCIGAR;
    int         cigarLen = 500;
    char        isRevComp;
    uint16_t    pos, posDNA;
    ogChromosome *pChr;
    //uint32_t    okSimpleAlignThreshold = seqLen * 3 / 100; //>> 4;
    uint32_t    okAlignThreshold = seqLen / 10;
    uint32_t    res = 0;
    uint32_t    Gpos, expectedGPos, closerGenPos;
    uint32_t    i,j,L,k,keyPos, nKeys, keySize = pMap->pKeyEncoding->getSizeInChars();
    uint32_t    maxKeyError = keySize / 10;
    uint32_t    pI;
    ogKey       *pK;
    uint32_t    minAlignedPos, maxAlignedPos;
    uint32_t    leftAddNt, rightAddNt;
    uint32_t    rightExactPos, leftExactPos;
    char        *pSeqN, *pDnaN, flag;
    float       maxErrRate = 0.0625, ovER; // 1/16
    ogGenomeAndReadPosition *garp;
    uint32_t    M,X,I,D,matches, nomatches, ins, del;

    alignNum++;
    /**
    fprintf(stderr, "|||| Aligning %llu : %s : %u targets\n", alignNum, pSeq, nGenPos);
    for (i=0; i < nGenPos && i < 255; i++) {
        garp = pCandPos->getkPos(i);
        pChr = pMap->pGenome->getGenomicCoordinate(garp->genomePosition);
        fprintf(stderr, "Rd Pos = %u, Gen Pos = %u:%u, RC = %u, Status = %c, Score = %u\n",garp->readPosition, pChr->number, garp->genomePosition - pChr->start, (garp->isReverse ? 1 : 0), garp->status, garp->score);
    }
     **/
    for (i=0; i < nGenPos && i < 255; i++) {
        if (pCandPos->getPositionKStatus(i) == 'e') {
            pGenRdPos = pCandPos->getkPos(i);
            Gpos = pGenRdPos->genomePosition; //pCandPos->getkPos(i);
            isRevComp = pGenRdPos->isReverse; //pCandPos->getIsReverseComplement(i);
            pos = pGenRdPos->readPosition; //pMap->getRelativePos(isRevComp);
            L = seqLen - pos;
            if (pos > Gpos) pos = Gpos;

            leftAddNt = rightAddNt = 0;
            nKeys = pMap->getNKeyPos(isRevComp);
            // Find left key position closer to the Gpos
            minAlignedPos = pos;
            for (k=0; leftAddNt == 0 || leftAddNt > 10 && k < nKeys; k++) {
                pI = pMap->getKeyablePosition(isRevComp,k);
                if (pI < pos) {
                    expectedGPos = Gpos - (pos - pI);
                    pK = pMap->getInfoForPositionKey(isRevComp, k);
                    closerGenPos = pMap->getCloserGenomePosition(expectedGPos, pMap->pGenPos->getPointerPosition(pK->offset), pK->size); 
                    if (closerGenPos == expectedGPos) {
                        minAlignedPos = pI;
                        leftAddNt = 0;
                        break;
                    } else if (closerGenPos < expectedGPos && expectedGPos-closerGenPos < MAX_GAP) {
                        if ((X = expectedGPos-closerGenPos) < leftAddNt || leftAddNt == 0) {
                            minAlignedPos = pI;
                            leftAddNt = X;
                        }
                    } else if (closerGenPos > expectedGPos && closerGenPos-expectedGPos < MAX_GAP) {
                        if ((X = closerGenPos-expectedGPos) < leftAddNt || leftAddNt == 0) {
                            minAlignedPos = pI;
                            leftAddNt = X;
                        }
                    }
                } else {
                    break;
                }
            }
            // Find right key position closer to the Gpos
            maxAlignedPos = pos;
            for (k=nKeys-1; rightAddNt == 0 || rightAddNt > 10 && k > 0; k--) {
                pI = pMap->getKeyablePosition(isRevComp,k);
                if (pI > pos) {
                    expectedGPos = Gpos - (pos - pI);
                    pK = pMap->getInfoForPositionKey(isRevComp, k);
                    closerGenPos = pMap->getCloserGenomePosition(expectedGPos, pMap->pGenPos->getPointerPosition(pK->offset), pK->size); 
                    if (closerGenPos == expectedGPos) {
                        maxAlignedPos = pI;
                        rightAddNt = 0;
                        break;
                    } else if (closerGenPos < expectedGPos && expectedGPos-closerGenPos < MAX_GAP) {
                        if ((X = expectedGPos-closerGenPos) < rightAddNt || rightAddNt == 0) {
                            maxAlignedPos = pI;
                            rightAddNt = X;
                        }
                    } else if (closerGenPos > expectedGPos && closerGenPos-expectedGPos < MAX_GAP) {
                        if ((X = closerGenPos-expectedGPos) < rightAddNt || rightAddNt == 0) {
                            maxAlignedPos = pI;
                            rightAddNt = X;
                        }
                    }
                } else {
                    break;
                }
            }

            pGenome->extractFromGenomicPos(Gpos-pos-leftAddNt, seqLen+leftAddNt+rightAddNt, dna, isRevComp, 0);
            // Error in Key
            if (isRevComp) {
               pos = seqLen-pos-keySize;
               posDNA = pos + rightAddNt;
            } else {
                posDNA = pos + leftAddNt;
            }
            for (e=0,n=0,pSeqN=pSeq+pos,pDnaN=dna+posDNA; *pSeqN && *pDnaN && n < keySize; n++, pSeqN++, pDnaN++) {
                if ((*pSeqN & 0x5F) != *pDnaN) {
                    e++;
                }
            }
            eK = e;
            nK = n;
            if (eK > maxKeyError) {
                // error key not aligned not supposed to happen
                /**
                fprintf(stderr, "***************************************\n");
                fprintf(stderr, "** Error in Key                      **\n");
                fprintf(stderr, "***************************************\n");
                fprintf(stderr, "RC:%c, pos=%u, Add[L=%u, R=%u]\n",isRevComp ? '1' : '0', pos, leftAddNt, rightAddNt);
                fprintf(stderr, "Key Start in Read:%s\n",pSeq+pos);
                fprintf(stderr, "Key Start in DNA :%s\n", dna+posDNA);
                 **/
            }
            // Error to the Right (Extending)
            for (rightExactPos=pos+n+leftAddNt, flag=1; (float)e/n <= maxErrRate && *pSeqN && *pDnaN; pSeqN++, pDnaN++, n++) {
                if ((*pSeqN & 0x5F) == *pDnaN) {
                    if (flag) rightExactPos++;
                } else {
                    e++;
                    flag=0;
                }
            }
            eR = (e-eK);
            nR = (n-nK);
            e  = eK;
            n  = nK;
            // Error to the Left (Extending)
            for (leftExactPos=pos+rightAddNt, flag=1, j=pos, pSeqN=pSeq+pos-1, pDnaN=dna+pos-1; 
                j > 0 && (float)e/n <= maxErrRate && *pSeqN && *pDnaN; 
                pSeqN--, pDnaN--, n++, j--) {
                if ((*pSeqN & 0x5F) == *pDnaN) {
                    if (flag) leftExactPos--;
                } else {
                    e++;
                    flag=0;
                }
            }
            eL = (e - eK);
            nL = (n - nK);
            ovER = (float) (eL+eK+eR) / (nL+nK+nR);
            //if (ovER < maxErrRate) res++;
            /**
            fprintf(stderr, "%u/%u:%c, Err[L=%u, K=%u, R=%u, Rate=%3.1f%%], Len Chck=%u/%u+(%u), qExact [L=%d, R=%d, Size=%u, from=%u, to=%u], Extended [L=%u, R=%u]\n", 
                    i, nGenPos, pCandPos->getPositionKStatus(i), eL, eK, eR, 100*ovER,  // ...ErrRate
                    nL + nK + nR, seqLen, leftAddNt+rightAddNt,
                    pos-leftExactPos, rightExactPos-pos-keySize, rightExactPos-leftExactPos, leftExactPos, rightExactPos, // ... qExact
                    leftAddNt, rightAddNt
                    );
            fprintf(stderr, "DNA: %s\n",dna);
             **/

            
            if (true) {
                // WFA2 (https://github.com/smarco/WFA2-lib)

                pWFAligner->alignEndsFree(pSeq, seqLen, 0, 0, dna, seqLen+leftAddNt+rightAddNt, 0, 0); // Pattern = Read, Text = Reference Genome
                pWFAligner->getAlignmentCigar(&pCIGAR, &cigarLen);
                //fprintf(stderr, "CIGAR: %s. Aligment score: %d.\n", pCIGAR, pWFAligner->getAlignmentScore());            

                matches = nomatches = ins = del = 0;
                pDnaN = pCIGAR;
                cigarLetLen = 0;
                while ((cCigar = *pDnaN)) {
                    for (M=X=0; cCigar == 'M' || cCigar == 'X'; ) { 
                        if (cCigar == 'M') M++; else X++; 
                        cCigar = *++pDnaN;
                    }
                    if (M+X > 0) { cigarLetLen += snprintf(cigarLet+cigarLetLen, MAX_CIGAR_SIZE-cigarLetLen, "%uM", M+X); matches += M; nomatches += X; }
                    for (I=0; cCigar == 'I'; I++) cCigar = *++pDnaN;
                    if (I > 0) { cigarLetLen += snprintf(cigarLet+cigarLetLen, MAX_CIGAR_SIZE-cigarLetLen, "%uI", I); ins += I; }
                    for (D=0; cCigar == 'D'; D++) cCigar = *++pDnaN;
                    if (D > 0) { cigarLetLen += snprintf(cigarLet+cigarLetLen, MAX_CIGAR_SIZE-cigarLetLen, "%uD", D); del += D; }
                }
                cigarLet[cigarLetLen] = 0;
                ovER = ((float) (nomatches+((I+D) >> 1)) / (float) (matches+nomatches));
                //fprintf(stderr, "CIGAR let:%s | M=%u,X=%u,D=%u,I=%u | Error=%f\n", cigarLet, matches, nomatches, del, ins, ovER);
                if ((matches+nomatches+del+ins) > seqLen*2) {
                    fprintf(stderr, "CIGAR let:%s | M=%u,X=%u,D=%u,I=%u | Error=%f | Prev Score=%hhu | Left+=%u | Right+=%u\n", cigarLet, matches, nomatches, del, ins, ovER, pCandPos->getPositionKScore(i), leftAddNt, rightAddNt);
                    pChr = pGenome->getGenomicCoordinate(Gpos); 
                    fprintf(stderr, "LARGE SEQ IN COMPLEX ALIGNMENT. GenPos=%u, Chr %u:%u, isRev=%hhu, pos=%u, RdNum=%llu\n%s\n", Gpos, pChr->number, Gpos - pChr->start, isRevComp, pos, pMap->read->readNumber, pSeq);
                }
                if (ovER < maxErrRate) {
                    res++;
                    pCandPos->setPositionKStatus(i, 'G'); // Global Alignment OK!
                } else {
                    pCandPos->setPositionKStatus(i, 'g'); // Global Alignment not working
                }
                pCandPos->setPositionKScore(i, 100*(1-ovER));

                /***
                EdlibAlignConfig     edlConf = edlibNewAlignConfig(
                                -1,                     // k
                                EDLIB_MODE_NW,          // 
                                EDLIB_TASK_PATH,
                                NULL, 0                 // No additional equalities
                            );
                EdlibAlignResult edResult = edlibAlign(pSeq, seqLen, dna, seqLen, edlConf);
                fprintf(stderr, "edLib Status: %d, distance: %d ", edResult.status, edResult.editDistance);
                if (edResult.status == EDLIB_STATUS_OK) {
                    char* cigar = edlibAlignmentToCigar(edResult.alignment, edResult.alignmentLength, EDLIB_CIGAR_STANDARD);
                    fprintf(stderr, "CIGAR: [%s]\n", cigar);
                    free(cigar);
                } else {
                    fprintf(stderr, "no cigar.\n");
                }
                edlibFreeAlignResult(edResult);
                ***/
            }


            if (false) {
                // First find keyIndex of seed position
                keyPos = 0x7FFFFFFF;
                for (i=(isRevComp ? pMap->getNRevPos() : pMap->getNFwdPos()); i > 0; i--) {
                    if (pMap->getKeyablePosition(isRevComp,i) == pos) {
                        keyPos = i;
                        break;
                    }
                }
                // Second find left and right keys containing exact relative positions
                i = keyPos;
                minAlignedPos = pos;
                while (i > 0) {
                    i--;
                    pI = pMap->getKeyablePosition(isRevComp,i);
                    expectedGPos = Gpos - (pos - pI);
                    pK = pMap->getInfoForPositionKey(isRevComp, i);
                    closerGenPos = pMap->getCloserGenomePosition(expectedGPos, pMap->pGenPos->getPointerPosition(pK->offset), pK->size); 
                    if (closerGenPos == expectedGPos) {
                        minAlignedPos = pI;
                    } else {
                        break;
                    }
                }
                i = keyPos;
                maxAlignedPos = pos;
                nKeys = (isRevComp ? pMap->getNRevPos() : pMap->getNFwdPos());
                while (++i < nKeys) {
                    pI = pMap->getKeyablePosition(isRevComp,i);
                    expectedGPos = Gpos - (pI - pos);
                    pK = pMap->getInfoForPositionKey(isRevComp, i);
                    closerGenPos = pMap->getCloserGenomePosition(expectedGPos, pMap->pGenPos->getPointerPosition(pK->offset), pK->size); 
                    if (closerGenPos == expectedGPos) {
                        maxAlignedPos = pI;
                    } else {
                        break;
                    }
                }
                // Listo

                e = globalAlignment(pSeq, seqLen, dna, seqLen);
                if (e < okAlignThreshold) {
                    res++;
                    pCandPos->setPositionKStatus(i, 'G'); // Global Alignment OK!
                } else {
                    pCandPos->setPositionKStatus(i, 'g'); // Global Alignment not working
                }
            }
        }
    }
    
    return res;
}


uint16_t ogAligner::globalAlignment(char *pA, uint32_t lenA, char *pB, uint32_t lenB) {
    
    allocate(lenA);
    allocate(lenB);
    
    // contrary to R like matrices, a row is continous and columns are discontinous
    // Matrix: A VERTICAL memory discontinous and do not depend on the length of B
    //         B HORIZONTAL memory continous, depends on B so it can be extended to allocate insertions
    //          [.,0] [.,1] [.,2] [.,3] [.,4]
    // [0,.]     0      -1    -1    -1    -1
    // [1,.]    -1
    // [2,.]    -1
    // [3,.]    -1
 
    // Parameters
    int16_t penaltyNT=1, match=1, penaltyOpen=4;
    int16_t s, d, l, u, x;
    int8_t  y;
    uint32_t nA = lenA+1;
    uint32_t nB = lenB+1;
    uint32_t ll = lenA + lenB - 1;
    
    uint32_t i,j,k,J1,J2,ij,i_1j,ij_1,i_1j_1,idx,_i,lim;

    // initialize
    //memset(pMatScore, -penaltyOpen, nA);
    //memset(pLeftMax, -penaltyOpen, nA);
    //memset(pUpMax, -penaltyOpen, nB);
    //memset(pLeftPos, 0, nA*sizeof(int16_t));
    //memset(pUpPos, 0, nB*sizeof(int16_t));
    x = -penaltyOpen;

    for (i=0; i <= lenA; i++) {
        pMatScore[i] = x;
        pLeftMax[i] = x;
        pLeftPos[i] = 0;
    }
    for (i=0; i <= lenB; i++) {
        pMatScore[i*nB] = x;
        pUpMax[i] = x;
        pUpPos[i] = 0;
    }
    pMatScore[0] = 0;
    pMatDecision[0] = 0;

    fprintf(stderr, "\n>%s:%d\n]%s:%d\n", pA, lenA, pB, lenB);
    
    pMatScore[nB] = (pA[0] == pB[0] ? match : 0);
    pMatDecision[nB] = 0;
    
    // align first part u: vertical control of A
    
    lim = nA * nB;
    for (i=1; i < ll; i++) {
        if (i < lenA) {
            J1 = 0;
            _i = i + 1;
        } else {
            J1 = i-lenA+1;
            _i = lenA;
        }
        for (j=J1; _i > 0 && j < lenB; j++) {
            _i--;
            // i - 1, j - 1
            i_1j_1 = _i*nB + j;
            // i, j
            ij = i_1j_1 + nB + 1;
            if (ij >= lim) {
                fprintf(stderr, "*Epale* i=%d, j=%d, ij=%d\n",i,j,ij);
            }
            
            // Computing [i,j]
            s = (pA[_i] == pB[j] ? match : 0);
            
            // Diagonal = s + M[i-1,j-1]
            d = s + pMatScore[i_1j_1];
            
            // Left = M[i-1..., j] - penalty
            l = pLeftMax[_i] - (penaltyOpen + (j-pLeftPos[_i])*penaltyNT);
            
            // Up = M[i, j-1...] - penalty
            u = pUpMax[j] - (penaltyOpen + (_i-pUpPos[j])*penaltyNT);
            
            // Max
            if (d >= l) {
                if (d >= u) {
                    x = d;
                    y = 0; // DIAG
                } else {
                    x = u;
                    y = 1; // UP
                }
            } else {
                if (l > u) {
                    x = l;
                    y = -1; // LEFT
                } else {
                    x = u;
                    y = 1; // UP
                }
            }
            
            
            // Store in [i,j]
            pMatScore[ij] = x;
            // Decision
            pMatDecision[ij] = y;
            
            // Update Left
            if (x >= pLeftMax[_i]) {
                pLeftMax[_i] = x;
                pLeftPos[_i] = j;
            }
            
            // Update Up
            if (x >= pUpMax[j]) {
                pUpMax[j] = x;
                pUpPos[j] = _i;
            }
        }
    }
    
    uint16_t mxl = (lenA > lenB ? lenA : lenB);

    FILE *f = fopen("alignment.txt", "w");
    fprintf(f,"       ");
    for (j=0; j < nB; j++) fprintf(f,"%5d", j % 10);
    fprintf(f,"\n       ");
    for (j=0; j < nB; j++) fprintf(f,"%5c", pB[j]);
    fprintf(f,"\n-- ");
    for (i=0; i < nA; i++) {
        idx = i * nB;
        if (i > 0) fprintf(f,"%d%c ", (i-1) % 10, pA[i-1]);
        for (j=0; j < nB; j++) {
            y = pMatDecision[idx+j];
            fprintf(f,"%4d%c",pMatScore[idx+j],(y == 0 ? '\\' : (y > 0 ? '^' : '<')));
        }
        fprintf(f,"\n");
    }
    fclose(f);
    fprintf(stderr, "** GLOBAL ALIGNMENT ** %d\n", x);
    return mxl - x;
}

uint64_t ogAligner::getFullAlignments() {
    return fullAligments;
}

char ogAligner::buildCIGARfromWFA(char *pCIGARchars, ogSingleRead *read, uint32_t leftPos, uint32_t rightPos) {
    char *cigarLet = read->pCigar;
    uint32_t matches = 0, nomatches = 0, ins = 0, del = 0, gapsOpen = 0;
    uint32_t M, X, soft = 0, x, m, xm, I=0, D=0, L, leftM = 0, maxMX = 0, posMX = 0;
    uint16_t transitions = 0;
    char cCigar, *pC, *pCIGARstart = pCIGARchars, *lastMX=cigarLet, *lastOP=cigarLet;
    uint32_t cigarLetLen = 0;
    uint16_t maxTransitions = (read->lenSeq >> 5);
    char print = 0;
    //*pLeftSoft = 0;
    //read->cigarLeft = 0;
    //read->cigarLeftType = '-';
    //read->cigarRight = 0;
    pCigarOperations->reset(leftPos, rightPos);
    while ((cCigar = *pCIGARchars)) {
        for (M=X=0; cCigar == 'M' || cCigar == 'X'; cCigar = *++pCIGARchars) if (cCigar == 'M') M++; else X++;
        if (M+X > 0) pCigarOperations->push((M >= X ? 'M' : 'X'), M+X);
        matches += M;
        nomatches += X;
        for (D=0; cCigar == 'D'; D++) cCigar = *++pCIGARchars;
        if (D > 0) pCigarOperations->push('D', D);
        del += D;
        for (I=0; cCigar == 'I'; I++) cCigar = *++pCIGARchars;
        if (I > 0) pCigarOperations->push('I', I);
        ins += I;
    }
    pCigarOperations->buildCIGAR(read->pCigarOriginal, MAX_CIGAR_SIZE);
    pCigarOperations->reduceOperations();
    pCigarOperations->buildCIGAR(cigarLet, MAX_CIGAR_SIZE);
    read->cigarLeftPos = pCigarOperations->leftGenPos;
    read->cigarRightPos = pCigarOperations->rightGenPos;
    read->cigarMatches = matches;
    read->cigarNoMatches = nomatches;
    read->cigarIns = ins;
    read->cigarDel = del;
    read->cigarSoft = 0;
    //if ((M = matches+nomatches+ins+del) == 0) M = 1;
    return (matches*100) / read->lenSeq;
    
    // old implementation
//    //fprintf(stderr, "CIGAR:");
//    while ((cCigar = *pCIGARchars)) {
//        for (M=X=0; cCigar == 'M' || cCigar == 'X'; cCigar = *++pCIGARchars) if (cCigar == 'M') M++; else X++;
//        if (M+X > 0) {
//            if (cigarLetLen == 0 && X > 2) {
//                // Es al inicio y hay mismatches
//                L = windowMismatches(pCIGARstart, M+X, 1, 3, 8, 5);
//                if (L > 0) {
//                    // hay S y M
//                    read->cigarLeft = L;
//                    read->cigarLeftType = 'S';
//                    cigarLetLen += snprintf(cigarLet+cigarLetLen, MAX_CIGAR_SIZE-cigarLetLen, "%uS", read->cigarLeft);
//                    transitions++;
//                    for (pC=pCIGARstart; L > 0; L--,pC++) if (*pC++ == 'M') M--; else X--;
//                }
//            }
//            if (cCigar == 0 && X > 2 && cigarLetLen > 0) { // Llego al final
//                L = windowMismatches(pCIGARchars-1, M+X, -1, 3, 8, 5);
//                if (L > 0) {
//                    soft = L;
//                    for (pC=pCIGARchars-soft; L > 0; L--) if (*pC++ == 'M') M--; else X--;
//                }
//                
//            }
//        }
//        if (M+X+leftM > 0) {
//            lastMX = cigarLet+cigarLetLen;
//            M += leftM;
//            cigarLetLen += snprintf(cigarLet+cigarLetLen, MAX_CIGAR_SIZE-cigarLetLen, "%uM", M+X); 
//            matches += M; 
//            nomatches += X;
//            leftM = 0;
//            transitions++;
//            if (M+X > maxMX) {
//                maxMX = M+X;
//                posMX = matches + nomatches + ins + (read->cigarLeftType == 'S' ? read->cigarLeft : 0);
//            }
//        }
//        
//        for (D=0; cCigar == 'D'; D++) cCigar = *++pCIGARchars;
//        if (D > 0) {
//            if (cigarLetLen == 0) {
//                // Delecion al inicio, ignorar
//                //read->cigarLeft = D;
//                //read->cigarLeftType = 'D';
//                //leftM += D;
//            } else if (matches+nomatches < 3 && D < 4) {
//                // Delecion "casi al inicio", 
//                // supones muchos matches adelante y si dejas la delecion, se recorre
//                // es decir es error de alineamiento
//                // ignorar
//            } else if (I == D && M+X == 0) {
//                // hubo una inserción (antes) y luego una delección con la misma cantidad, eliminar
//                cigarLetLen = lastOP - cigarLet;
//                if (cCigar == 0) {
//                    soft += I;
//                } else {
//                    cigarLetLen += snprintf(cigarLet+cigarLetLen, MAX_CIGAR_SIZE-cigarLetLen, "%uM", I);
//                    nomatches += I;
//                }
//                ins -= I; // quitar lo sumado
//                print = 1;
//            } else {
//                if (cCigar == 0) {
//                    // D al final, ignorar
//                } else {
//                    gapsOpen++; cigarLetLen += snprintf(lastOP=cigarLet+cigarLetLen, MAX_CIGAR_SIZE-cigarLetLen, "%uD", D); del += D;
//                    transitions++;
//                }
//            }
//        }
//        
//        for (I=0; cCigar == 'I'; I++) cCigar = *++pCIGARchars;
//        if (I > 0) { 
//            if (D == 0 && cCigar == 0 && M+X > 0 && I <= del) {
//                // llegó al final & no viene de una delecion sino de un Match y el # de ins es "comparable" al de deleciones 
//                // ==> No es inserción, es que le falto alinear el pedazo de la delecion al final, y por lo tanto dejar como match.
//                cigarLetLen = snprintf(lastMX, MAX_CIGAR_SIZE-cigarLetLen, "%uM", M+X+I) + (lastMX - read->pCigar);
//                matches += I;
//                read->cigarRight = I; // se agregaron, se requiere corregir posteriormente
//                //transitions++;
//            } else {
//                if (cigarLetLen == 0) {
//                //    // Inserción al inicio cambiar a soft clipping
//                //    cigarLetLen += snprintf(cigarLet, MAX_CIGAR_SIZE, "%uS", I); 
//                //    read->cigarLeftSoft = I; 
//                    read->cigarLeft = I;
//                    read->cigarLeftType = 'I';
//                    leftM += I;
//                    // se agregaron I a la izquierda, requiere corregir posición
//                } else if (cCigar == 0 && M+X > 0) {
//                    // Llegó al final y es insercion, dejar como matches
//                    cigarLetLen = snprintf(lastMX, MAX_CIGAR_SIZE-cigarLetLen, "%uM", M+X+I) + (lastMX - read->pCigar);
//                    matches += I;
//                    transitions++;
//                } else if (I == D) {
//                    // hubo una delección y luego una inserción con la misma cantidad, eliminar
//                    cigarLetLen = lastOP - cigarLet;
//                    if (cCigar == 0) {
//                        soft += I;
//                    } else {
//                        cigarLetLen += snprintf(cigarLet+cigarLetLen, MAX_CIGAR_SIZE-cigarLetLen, "%uM", I);
//                        nomatches += I;
//                    }
//                    del -= D; // quitar lo sumado
//                    print = 1;
//                } else {
//                    gapsOpen++; cigarLetLen += snprintf(lastOP=cigarLet+cigarLetLen, MAX_CIGAR_SIZE-cigarLetLen, "%uI", I); ins += I;
//                    transitions++;
//                }
//            }
//        }
//    }
//    if (soft > 0) {
//        cigarLetLen += snprintf(cigarLet+cigarLetLen, MAX_CIGAR_SIZE-cigarLetLen, "%uS", soft);
//        transitions++;
//    }
//    cigarLet[cigarLetLen] = 0;
//    //fprintf(stderr, "\nALIGN:%s\n",pCIGARstart);
//    //fprintf(stderr, "CIGAR:%s\n",cigarLet);
//    if (transitions > maxTransitions) {
//        ///if (posMX > read->lenSeq) {
//        ///    fprintf(stderr, "+++++++ CIGAR:%s\n", cigarLet);
//        ///}
//        cigarLetLen = 0;
//        read->cigarLeft = 0;
//        read->cigarLeftType = '-';
//        soft = 0;
//        if (posMX > maxMX) {
//            read->cigarLeft = posMX-maxMX;
//            read->cigarLeftType = 'S';
//            cigarLetLen += snprintf(cigarLet+cigarLetLen, MAX_CIGAR_SIZE-cigarLetLen, "%uS", posMX-maxMX);
//        }
//        cigarLetLen += snprintf(cigarLet+cigarLetLen, MAX_CIGAR_SIZE-cigarLetLen, "%uM", maxMX);
//        if (posMX < read->lenSeq) {
//            soft = read->lenSeq - posMX;
//            cigarLetLen += snprintf(cigarLet+cigarLetLen, MAX_CIGAR_SIZE-cigarLetLen, "%uS", soft);
//        }
//        cigarLet[cigarLetLen] = 0;
//        nomatches = del = ins = 0;
//        matches = maxMX;
//        ///if (posMX > read->lenSeq) {
//        ///    fprintf(stderr, "-------- CIGAR:%s\n", cigarLet);
//        ///}
//    }
//    
//    read->cigarMatches = matches;
//    read->cigarNoMatches = nomatches;
//    read->cigarIns = ins;
//    read->cigarDel = del;
//    read->cigarSoft = soft + (read->cigarLeftType == 'S' ? read->cigarLeft : 0);
//    //return (matches*100) / (matches+nomatches+gapsOpen);
//    if ((M = matches+nomatches+ins+del) == 0) M = 1;
//    
//    if (print) {
//        fprintf(stderr, "////// Seq ID: %s, CIGARS: Long:%s, Short:%s\n", read->pId, pCIGARstart, read->pCigar);
//    }
//    return (matches*100) / M;
}


void ogAligner::writeOnlyMatchesCigar(ogSingleRead *read, uint32_t leftPos, uint32_t rightPos) {
    snprintf(read->pCigar, MAX_CIGAR_SIZE, "%dM", read->lenSeq);
    //read->cigarDel = read->cigarIns = read->cigarLeft = read->cigarNoMatches = read->cigarRight = read->cigarSoft = 0;
    //read->cigarLeftType = '-';
    read->cigarDel = read->cigarIns = read->cigarNoMatches = read->cigarSoft = 0;
    read->cigarMatches = read->lenSeq;
    read->cigarLeftPos = leftPos;
    read->cigarRightPos = rightPos;
}

/// CIGAR CASOS
// soft al inicio de los que no matchean
// soft al final de los que no matchean
// I al final pegados con M, dejar como M (se asume si del > 0)
// X=M

/**
char ogAligner::buildCIGARfromWFA(char *pCIGARchars, ogSingleRead *read) {
    char *cigarLet = read->pCigar;
    uint32_t matches = 0, nomatches = 0, ins = 0, del = 0, gapsOpen = 0;
    uint32_t M, X, soft = 0, x, m, xm, I, D;
    char cCigar, *pC, *pCIGARstart = pCIGARchars, *lastMX=cigarLet;
    uint32_t cigarLetLen = 0;
    // *pLeftSoft = 0;
    read->cigarLeftSoft = 0;
    //fprintf(stderr, "CIGAR:");
    while ((cCigar = *pCIGARchars)) {
        for (M=X=0; cCigar == 'M' || cCigar == 'X'; cCigar = *++pCIGARchars) if (cCigar == 'M') M++; else X++;
        if (M+X > 0) {
            if (cigarLetLen == 0 && X > 2) {
                // Es al inicio y hay mismatches
                if (M > 0) {
                    // Tambien hay matches ... Search desde el inicio hasta encontrar
                    pC = pCIGARstart; //cigarLet;
                    x = m = xm = 0;
                    while (*pC == 'M' || *pC == 'X') {
                        if (*pC == 'M') {
                            if (++m > 4) break;
                        } else {
                            x++;
                            xm += m;
                            m = 0;
                        }
                        pC++;
                    }
                    if (x > 1) { 
                        cigarLetLen += snprintf(cigarLet+cigarLetLen, MAX_CIGAR_SIZE-cigarLetLen, "%uS", x+xm); 
                        X -= x; 
                        M -= xm;
                        read->cigarLeftSoft = x+xm;
                    }
                } else {
                    // Todos son mismatches ... Just add soft
                    cigarLetLen += snprintf(cigarLet+cigarLetLen, MAX_CIGAR_SIZE-cigarLetLen, "%uS", X); 
                    read->cigarLeftSoft = X;
                    X = 0; 
                }
            }
            if (cCigar == 0 && X > 2) { // Llego al final
                pC = pCIGARchars - 1;
                x = m = xm = 0;
                while (*pC == 'M' || *pC == 'X') {
                    if (*pC == 'M') {
                        if (++m > 4) break;                                    
                    } else {
                        x++;
                        xm += m;
                        m = 0;
                    }
                    if (--pC < pCIGARstart) break;
                }
                if (x > 1) {
                    soft = x + xm;
                    X -= x;
                    M -= xm;
                }
            }
            if (M+X > 0) {
                lastMX = cigarLet+cigarLetLen;
                cigarLetLen += snprintf(cigarLet+cigarLetLen, MAX_CIGAR_SIZE-cigarLetLen, "%uM", M+X); 
            }
            matches += M; 
            nomatches += X;
        }
        for (D=0; cCigar == 'D'; D++) cCigar = *++pCIGARchars;
        if (D > 0) { gapsOpen++; cigarLetLen += snprintf(cigarLet+cigarLetLen, MAX_CIGAR_SIZE-cigarLetLen, "%uD", D); del += D; }
        for (I=0; cCigar == 'I'; I++) cCigar = *++pCIGARchars;
        if (I > 0) { 
            if (cCigar == 0 && D == 0 && M+X > 0 && I <= del) {
                // llegó al final & no viene de una delecion sino de un Match y el # de ins es "comparable" al de deleciones 
                // ==> No es inserción, es que le falto alinear el pedazo de la delecion al final, y por lo tanto dejar como match.
                snprintf(lastMX, MAX_CIGAR_SIZE-cigarLetLen, "%uM", M+X+I);
                cigarLetLen += I;
            } else {
                gapsOpen++; cigarLetLen += snprintf(cigarLet+cigarLetLen, MAX_CIGAR_SIZE-cigarLetLen, "%uI", I); ins += I; 
            }
        }
    }
    if (soft > 0) cigarLetLen += snprintf(cigarLet+cigarLetLen, MAX_CIGAR_SIZE-cigarLetLen, "%uS", soft);
    cigarLet[cigarLetLen] = 0;
    
    read->cigarMatches = matches;
    read->cigarNoMatches = nomatches;
    read->cigarIns = ins;
    read->cigarDel = del;
    read->cigarSoft = read->cigarLeftSoft + soft;
    //return (matches*100) / (matches+nomatches+gapsOpen);
    if ((M = matches+nomatches+ins+del) == 0) M = 1;
    return (matches*100) / M;
}
*/


// regresa el # de caracteres que hay que avanzar desde pStart en direccion direction que no sobrepasen maxMismatches
uint32_t ogAligner::windowMatches(char *pStart, uint32_t maxLen, signed char direction, uint16_t maxMismatches, uint16_t window) {
    
    uint32_t matches = 0, mismatches = 0, len = 0;
    uint32_t p = 0;
    char *pI = pStart;
    int16_t offset = (direction > 0 ? -window :    window);
    
    for(; len < maxLen; len++, pI += direction) {
        if (*pI == 'M') matches++; else mismatches++;
        if (len > window) { if (*(pI+offset) == 'M') matches--; else mismatches--; }
        if (mismatches >= maxMismatches) {
            // terminar, pero primero encontrar el mejor punto desde el inicio de la ventana hasta el 1er mismatch
            for(pI += offset; *pI == 'M'; pI += direction, len--);
            return len;
        }
    }
    return len;

}

// regresa el # de caracteres que hay que avanzar desde pStart en direccion direction que si sobrepasen maxMismatches
uint32_t ogAligner::windowMismatches(char *pStart, uint32_t maxLen, signed char direction, uint16_t minMismatches, uint16_t window, uint16_t breakMatches) {
    
    uint32_t  len = 0;
    uint16_t mismatches = 0, matches = 0;
    char *pI = pStart;
    int16_t offset = (direction > 0 ? -window :    window);
    
    for(; len < maxLen && len < window; len++, pI += direction) {
        if (*pI == 'M') { 
            if (++matches >= breakMatches) return mismatches;
        } else {
            mismatches++;
            matches = 0;
        }
    }
    
    // si no cumple la primera ventana con mismatches, pos no hay nada que hacer.
    if (mismatches < minMismatches) return 0;
    
    for(; mismatches >= minMismatches && len < maxLen && matches < breakMatches; len++, pI += direction) {
        if (*pI == 'M') ++matches;
        else {
            mismatches++;
            matches = 0;
        }
        if (*(pI+offset) != 'M') mismatches--;
    }
    for (pI -= direction; *pI == 'M' && len > 0; pI -= direction, len--);
    
    return len;

}
    