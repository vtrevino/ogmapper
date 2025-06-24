
/* 
 * File:   ogCandidatePosManager.cpp
 * Author: victortrevino
 *
 * Created on May 4, 2022, 4:27 PM
 */

#include <string>
#include <stdint.h>
#include <stdlib.h>
#include "ogDefinitions.hpp"
#include "ogCandidatePosManager.hpp"


ogCandidatePosManager::ogCandidatePosManager(uint32_t maxArraySize) {
    positions = NULL;
    pairedPositions = NULL;
    arraySize = 0;
    pairedArraySize = 0;
    allocateForSize(maxArraySize);
    allocatePairedForSize(maxArraySize >> 4);
    count = 0;
    savedFwdCount = 0;
    pTree = new ogReusableBinaryTree(); // This is free to use for mapping functions
    changed = 1;
    currentFunction = '-';
}

ogCandidatePosManager::~ogCandidatePosManager() {
    //fprintf(stderr, "<Deallocating ogCandidatePosManager:"); fflush(stderr);
    //fprintf(stderr, "<~ogCPM:"); fflush(stderr);
    free(positions);
    if (pairedPositions != NULL) free(pairedPositions);
    if (pTree != NULL) delete pTree;
    //fprintf(stderr, ":ogCPM>");
}

void ogCandidatePosManager::allocateForSize(uint32_t size) {
    if (size >= arraySize) {
        arraySize = size + size/10;
        if (arraySize < MIN_CANDIDATE_POSITIONS) { arraySize = MIN_CANDIDATE_POSITIONS; }
        //fprintf(stderr, ":ogCPM:arraySize=%u\n",arraySize);
        positions = (ogGenomeAndReadPosition *) realloc(positions, arraySize * sizeof(ogGenomeAndReadPosition));
    }
}

char ogCandidatePosManager::allocatePairedForSize(uint32_t size) {
    if (size >= pairedArraySize) {
        pairedArraySize = size + size/10;
        if (pairedArraySize == 0) pairedArraySize = 128;
        if (pairedArraySize < (MIN_CANDIDATE_POSITIONS >> 4)) { pairedArraySize = (MIN_CANDIDATE_POSITIONS >> 4); }
        //fprintf(stderr, ":ogCPM:pairedArraySize=%u,size=%hu\n",pairedArraySize,size);
        pairedPositions = (ogGenomeAndReadPosition *) realloc(pairedPositions, pairedArraySize * sizeof(ogGenomeAndReadPosition));
        return 1;
    }
    return 0;
}

void ogCandidatePosManager::reset() {
    prevCount = count = savedFwdCount = 0;
    changed = 1;
    pairedArrayCount =  0;
}

uint32_t ogCandidatePosManager::getArraySize() {
    return arraySize;
}


uint32_t ogCandidatePosManager::setByInsersectingKeys(uint32_t *pa, uint32_t aSize, uint16_t aRdPos, uint32_t *pb, uint32_t bSize, uint16_t bRdPos, uint16_t maxDiff, char isRev) {
    uint32_t ai, bi;
    uint32_t ga, gb, kount = count;
    
    // arraySize is not a problem now because addGenome func calls allocate and is handled

    //if (count >= arraySize) {
    //    fprintf(stderr, "***** PROBLEM possible overpass (0) Count=%u, aSize=%u, bSize=%u, maxDiff=%u ****\n", count, aSize, bSize, maxDiff);
    //    return 0; // possible overpass
    //}
    ga = *pa;
    gb = *pb;
    changed = 1;

    ai = bi = 0;
    while (1) {
        //for (ai=bi=0; ai < aSize && bi < bSize && count < arraySize; ) 
        if (ga <= gb) {
            if (gb - ga <= maxDiff) {
                // THIS IS A HIT
                addGenomeAndReadPosition(ga, aRdPos, isRev, 0);
                gb = *++pb;
                //bi++;
                //if (count >= arraySize || bi >= bSize) break;
                if (++bi >= bSize) break;
            }
            ga = *++pa;
            if (++ai >= aSize) break;            
        } else {
            if (ga - gb <= maxDiff) {
                // THIS IS A HIT
                addGenomeAndReadPosition(gb, bRdPos, isRev, 0);
                ga = *++pa;
                //ai++;
                //if (count >= arraySize || ai >= aSize) break;
                if (++ai >= aSize) break;
            }            
            gb = *++pb;
            if (++bi >= bSize) break;
        }
    }
    
    //if (count >= arraySize) {
    //    fprintf(stderr, "***** PROBLEM possible overpass (1) Count=%u, Kount=%u, aSize=%u, bSize=%u, maxDiff=%u ****\n", count, kount, aSize, bSize, maxDiff);
    //}
    if (count-kount > 256) { // DEBE SER PARAMETRO
        /// Demasiado intersect, no vale la pena go deeper, regresar a 0
        count = kount;        
    }
        
    return count - kount;
}


// Similar a setByInsersectingKeys pero siempre tomando de base la posicion genómica izquierda (a)
uint32_t ogCandidatePosManager::setByInsersectingKeysLeftKey(uint32_t *pa, uint32_t aSize, uint16_t aRdPos, uint32_t *pb, uint32_t bSize, uint16_t maxDiff, char isRev, uint16_t maxAdd) {
    uint32_t ai, bi;
    uint32_t ga, gb, kount = count;
    uint32_t maxArray = count + maxAdd;
    
    if (maxArray > arraySize) maxArray = arraySize;
    
    // arraySize is not a problem now because addGenome func calls allocate and is handled
    //if (count >= arraySize) {
    //    fprintf(stderr, "***** PROBLEM possible overpass (0) Count=%u, aSize=%u, bSize=%u, maxDiff=%u ****\n", count, aSize, bSize, maxDiff);
    //    return 0; // possible overpass
    //}
    ga = *pa;
    gb = *pb;
    changed = 1;

    ai = bi = 0;
    while (1) {
        //for (ai=bi=0; ai < aSize && bi < bSize && count < arraySize; ) 
        if (ga <= gb) {
            if (gb - ga <= maxDiff) {
                // THIS IS A HIT
                addGenomeAndReadPosition(ga, aRdPos, isRev, 0);
                gb = *++pb;
                //bi++;
                //if (count >= maxArray || bi >= bSize) break;
                if (++bi >= bSize) break;
            }
            ga = *++pa;
            if (++ai >= aSize) break;            
        } else {
            if (ga - gb <= maxDiff) {
                // THIS IS A HIT
                addGenomeAndReadPosition(ga, aRdPos, isRev, 0); // use "A" because this is for "Left Keys" only
                ga = *++pa;
                //ai++;
                //if (count >= maxArray || ai >= aSize) break;
                if (++ai >= aSize) break;
            }            
            gb = *++pb;
            if (++bi >= bSize) break;
        }
    }
    
    //if (count >= arraySize) {
    //    fprintf(stderr, "***** PROBLEM possible overpass (1) Count=%u, Kount=%u, aSize=%u, bSize=%u, maxDiff=%u ****\n", count, kount, aSize, bSize, maxDiff);
    //}
    if (count-kount > 256) { // DEBE SER PARAMETRO
        /// Demasiado intersect, no vale la pena go deeper, regresar a 0
        count = kount;        
    }
        
    return count - kount;
}

uint32_t ogCandidatePosManager::intersectAddingKey(uint32_t *pa, uint32_t aSize, uint16_t aRdPos, uint16_t maxDiff, char isRev) {
    uint32_t ai, j, k;
    uint32_t ga, gb;
    uint32_t limit = count;

    if (limit > arraySize) limit = arraySize;
    j = savedFwdCount;      // current positions in b to be compared
    k  = savedFwdCount;     // new position to be replaced/saved
    ga = *pa;               // a is the new key
    gb = getkPos(j)->genomePosition;      // b is the internal positions
    ai = 0;
    changed = 1;

    while (1) {
        if (ga <= gb) {
            if (gb - ga <= maxDiff) {
                // THIS IS A HIT, keep ga Position in k
                setGenomeAndReadPositionK(k++, ga, aRdPos, isRev, 0); // positions[k++] = ga;
                if (++j >= limit) break;
                gb = getkPos(j)->genomePosition; // positions[j];
            }
            ga = *++pa;
            if (++ai >= aSize) break;
        } else {
            if (ga - gb <= maxDiff) {
                // THIS IS A HIT, keep k position
                uint16_t rbpos = getkPos(j)->readPosition;
                setGenomeAndReadPositionK(k++, gb, rbpos, isRev, 0); // positions[k++] = gb;
                if (++j >= limit) break;
                gb = getkPos(j)->genomePosition; //positions[j];
                ga = *++pa;
                if (++ai >= aSize) break;
            } else {
                if (++j >= limit) break;
                gb = getkPos(j)->genomePosition; //positions[j];
            }
        }
    }
    if (j > arraySize) {
        fprintf(stderr, "***** PROBLEM possible overpass (2) i=%u, k=%u, savedFwdCount=%u, arraySize=%u ****\n", j, k, savedFwdCount, arraySize);
    }    
    prevCount = count;
    count = k;
    return k - savedFwdCount;
}

uint32_t ogCandidatePosManager::getCount() {
    return count;
}

uint32_t ogCandidatePosManager::getFwdCount() {
    return savedFwdCount;
}

uint32_t ogCandidatePosManager::getRevCount() {
    return count-savedFwdCount;
}

ogGenomeAndReadPosition * ogCandidatePosManager::getkPos(uint32_t k) {
    return positions + k;
}

// void ogCandidatePosManager::setByFixedPosition(uint32_t gpos, uint32_t rdpos, char isRev) {
//     reset();
//     addGenomeAndReadPosition(gpos, rdpos, isRev);
// }

// // Ahora usar addGenomeAndReadPosition 
// void ogCandidatePosManager::addFixedPosition(uint32_t gpos, uint32_t rdpos, char isRev) {
//     if (isRevComp) {
//         positions[savedFwdCount + count++] = pos;
//     } else {
//         positions[savedFwdCount++ + count] = pos;
//     }
//     if (savedFwdCount + count >= arraySize -1 ) {
//         fprintf(stderr, "***** PROBLEM possible overpass (3) ****\n");
//     }
// }

uint32_t ogCandidatePosManager::rewindCountsAfterAddingKeyCountZero() {
    // if (count != 0) PROBLEMS, because positions have been overwritten
    count = prevCount;
    changed = 1;
    return count;
}

void ogCandidatePosManager::saveFwdCandidatePositions() {
    prevCount = savedFwdCount = count;
}

void ogCandidatePosManager::restoreFwdCandidatePositions() {
    changed = 1;
    count = savedFwdCount;
}

//void ogCandidatePosManager::saveCandidatePositions() {
//    uint32_t j = arraySize;
//    uint32_t i = 0;
//    while (i < count) {
//        positions[--j] = positions[i++];
//    }
//    posSaved = j;
//}
//
//uint32_t ogCandidatePosManager::restoreSavedCandidatePositions() {
//    count = 0;
//    return copySavedCandidatePositionsTo(0);
//}
//
//uint32_t ogCandidatePosManager::addSavedCandidatePositions() {
//    return copySavedCandidatePositionsTo(count);
//}
//
//uint32_t ogCandidatePosManager::copySavedCandidatePositionsTo(uint32_t i) {
//    uint32_t j = arraySize;
//    while (j > posSaved) {
//        positions[i++] = positions[--j];
//    }
//    count += (arraySize - j);
//    return count;
//}

char ogCandidatePosManager::checkToAddGenomeAndReadPosition(uint32_t genPos, uint16_t readPos, char isReverse, unsigned char score) {
    //if (count+1 >= arraySize) return 2;
    uint32_t GP = genPos - readPos;
    ogGenomeAndReadPosition *p = positions, *pMax = positions + count;
    for (; p < pMax; p++) {
        if (((p->genomePosition-p->readPosition) == GP) && (isReverse == p->isReverse)) {
            return 0;
        }
    }
    addGenomeAndReadPosition(genPos, readPos, isReverse, score);
    return 1;
}

void ogCandidatePosManager::addGenomeAndReadPosition(uint32_t genPos, uint16_t readPos, char isReverse, unsigned char score) {
    allocateForSize(count);
    ogGenomeAndReadPosition *p = positions + count;
    p->genomePosition = genPos;
    p->readPosition = readPos;
    p->isReverse = isReverse;
    p->status = (score  > 0 ? 's' : 'u');
    p->score = score;
    p->function = currentFunction;
    count++;
    changed = 1;
}

void ogCandidatePosManager::setGenomeAndReadPositionK(uint32_t k, uint32_t genPos, uint16_t readPos, char isReverse, unsigned char score) {
    allocateForSize(k);
    ogGenomeAndReadPosition *p = positions + k;
    p->genomePosition = genPos;
    p->readPosition = readPos;
    p->isReverse = isReverse;
    p->status = (score  > 0 ? 's' : 'u');
    p->score = score;
    p->function = currentFunction;
    changed = 1;
}

void ogCandidatePosManager::initPairedGenomeAndReadPosition(ogGenomeAndReadPosition *p) {
    p->pairedCounts = 0;
    p->pairedIndex = 0;
    p->pairedKeys = 0;
    p->maxPairedIndex = 0;
}

void ogCandidatePosManager::initPairedGenomeAndReadPosition(uint32_t k) {
    initPairedGenomeAndReadPosition(positions + k);
}

uint32_t ogCandidatePosManager::checkPairedGenomeAndReadPosition(uint32_t k, uint32_t genPos, uint16_t readPos, char isReverse) {
    return checkPairedGenomeAndReadPosition(positions + k, genPos, readPos, isReverse);
}

uint32_t ogCandidatePosManager::checkPairedGenomeAndReadPosition(ogGenomeAndReadPosition *gp, uint32_t genPos, uint16_t readPos, char isReverse) {
    uint16_t i;
    uint32_t GP = genPos - readPos;
    ogGenomeAndReadPosition *p = pairedPositions + gp->pairedIndex;
    uint16_t iMax = gp->pairedCounts;
    for (i=0; i < iMax; p++, i++) {
        if (((p->genomePosition-p->readPosition) == GP) && (isReverse == p->isReverse)) {
            return gp->pairedIndex+i+1;
        }
    }
    return 0;
}

char ogCandidatePosManager::checkToAddPairedGenomeAndReadPosition(ogGenomeAndReadPosition *gp, uint32_t genPos, uint16_t readPos, char isReverse, unsigned char score) {
    uint32_t pos = checkPairedGenomeAndReadPosition(gp, genPos, readPos, isReverse);
    if (pos == 0) {
        addPairedGenomeAndReadPosition(gp, genPos, readPos, isReverse, score);
        return 1;
    }
    return pos;
}

char ogCandidatePosManager::checkToAddPairedGenomeAndReadPosition(uint32_t k, uint32_t genPos, uint16_t readPos, char isReverse, unsigned char score) {
    return checkToAddPairedGenomeAndReadPosition(positions + k, genPos, readPos, isReverse, score);
}

uint32_t ogCandidatePosManager::addPairedGenomeAndReadPosition(uint32_t k, uint32_t genPos, uint16_t readPos, char isReverse, unsigned char score) {
    return addPairedGenomeAndReadPosition(positions + k, genPos, readPos, isReverse, score);
}

uint32_t ogCandidatePosManager::addPairedGenomeAndReadPosition(ogGenomeAndReadPosition *gp, uint32_t genPos, uint16_t readPos, char isReverse, unsigned char score) {
    if (gp->pairedCounts == 0) gp->pairedIndex = pairedArrayCount;
    if (allocatePairedForSize(++pairedArrayCount)) { // General count
        //fprintf(stderr, "Paired: Counts=%u, Index=%u, Max=%u\n", gp->pairedCounts, gp->pairedIndex, gp->maxPairedIndex);
    }
    ogGenomeAndReadPosition *p = pairedPositions + gp->pairedIndex + gp->pairedCounts;
    p->genomePosition = genPos;
    p->readPosition = readPos;
    p->isReverse = isReverse;
    p->status = 'P';
    p->score = score;
    return (gp->pairedIndex + ++gp->pairedCounts); // to use this returned index use: idx-1 relative to pairedPositions as in getkPairedPos
}

void ogCandidatePosManager::resetPairedGenomeAndReadPosition(ogGenomeAndReadPosition *gp) {
    // supone que gp es el último de la cola del arreglo de pairedArray
    if (gp->pairedCounts > 0) {
        pairedArrayCount -= gp->pairedCounts;
        gp->pairedCounts = 0;
        gp->pairedIndex = 0;
        gp->pairedKeys = 0;
    }
}

ogGenomeAndReadPosition * ogCandidatePosManager::getkPairedPos(uint32_t k) {
    return (pairedPositions + k);
}

char ogCandidatePosManager::areThereMoreForwardCounts() {
    uint32_t fwd = savedFwdCount;
    uint32_t rev = count - savedFwdCount;
    return (fwd >= rev);
}

void ogCandidatePosManager::setAllStatus(char stat) {
    uint32_t i;
    for (i=0; i < count; i++) {
        positions[i].status = stat;
    }
    changed = 1;
}

void ogCandidatePosManager::setPositionKStatus(uint32_t k, char stat) {
    positions[k].status = stat;
    changed = 1;
}

char ogCandidatePosManager::getPositionKStatus(uint32_t k) {
    return positions[k].status;
}

void ogCandidatePosManager::setPositionKScore(uint32_t k, char score) {
    positions[k].score = score;
    changed = 1;
}

char unsigned ogCandidatePosManager::getPositionKScore(uint32_t k) {
    return positions[k].score;
}

void ogCandidatePosManager::removeForwardTargets() {
    if (savedFwdCount > 0) {
        uint32_t n = count-savedFwdCount;
        memmove(positions, positions+savedFwdCount, n*sizeof(ogGenomeAndReadPosition));
        count = n;
        savedFwdCount = 0;
    }
}

void ogCandidatePosManager::removeReverseTargets() {
    count = savedFwdCount;
    changed = 1;
}

void ogCandidatePosManager::print() {
    int i;
    fprintf(stderr, "----- Printing %u Candidate Positions Fwd:%u, Rev:%u ------\n", count, savedFwdCount, count - savedFwdCount);
    for (i=0; i < count; i++) {
        fprintf(stderr, "i=%u, RdPos=%u, GenPos=%u, RC=%c, status=%c, score=%d\n", i, 
                positions[i].readPosition,
                positions[i].genomePosition,
                positions[i].isReverse+48,
                positions[i].status,
                positions[i].score
                );
    }
}

void ogCandidatePosManager::printWithGenomePositions(ogGenome *pGenome) {
    int i;
    ogChromosome *pChr;
    fprintf(stderr, "----- Printing %u Candidate Positions Fwd:%u, Rev:%u ------\n", count, savedFwdCount, count - savedFwdCount);
    for (i=0; i < count; i++) {
        pChr = pGenome->getGenomicCoordinate(positions[i].genomePosition);
        fprintf(stderr, "i=%u, RdPos=%u, GenPos=%u, Chr=%u:%u, RC=%c, status=%c, score=%d\n", i, 
                positions[i].readPosition,
                positions[i].genomePosition,
                pChr->number, positions[i].genomePosition - pChr->start,
                positions[i].isReverse+48,
                positions[i].status,
                positions[i].score
                );
    }
}

void ogCandidatePosManager::estimateFwdRevMatches(char all) {
    if (changed) {
        uint32_t i;
        ogGenomeAndReadPosition *p = positions;
        char rev;
        fwdMatches = revMatches = 0;
        maxScore[0] = 0;
        maxScore[1] = 0;
        minScore[0] = 100;
        minScore[1] = 100;
        maxScorePos[0] = 0;
        maxScorePos[1] = 0;
        maxScoreTies[0] = 0;
        maxScoreTies[1] = 0;
        for (i=0; i < count; i++, p++) {
            if (all || p->status == 'A' || p->status == 'G') {
                if (p->isReverse) {
                    rev = 1;
                    revMatches++;
                } else {
                    rev = 0;
                    fwdMatches++;
                }
                if (p->score >= maxScore[rev]) {
                    if (p->score == maxScore[rev]) {
                        maxScoreTies[rev]++;
                    } else {
                        maxScore[rev] = p->score;
                        maxScorePos[rev] = i;
                        maxScoreTies[rev] = 0;
                    }
                }
                if (p->score < minScore[rev]) {
                    minScore[rev] = p->score;
                }
            }
        }
        changed = 0;
    }
}

uint32_t ogCandidatePosManager::getFwdMatches() {
    return fwdMatches;
}

uint32_t ogCandidatePosManager::getRevMatches() {
    return revMatches;
}

void ogCandidatePosManager::cleanCandidatePositions() {
    //No hay un criterio claro de porque quitar regiones candidatas sin mapear el read 
    return;
    //uint32_t i;
    //ogGenomeAndReadPosition *p = positions;
    //estimateFwdRevMatches();
    //if (fwdMatches > 0 && maxScore[0]-minScore[0] > 10) {        
    //}
}

