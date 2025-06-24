    /* 
 * File:   ogCigarOperations.hpp
 * Author: victortrevino
 *
 * Created on March 20, 2024, 05:02 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "ogDefinitions.hpp"
#include "ogCigarOperations.hpp"


ogCigarOperations::ogCigarOperations() {
    cigarOp = (OG_CIGAR_OPERATION *) malloc(sizeof(OG_CIGAR_OPERATION) * MAX_CIGAR_SIZE * 2);
    reset(0,0);
}

ogCigarOperations::~ogCigarOperations() {
    free(cigarOp);
}

void ogCigarOperations::reset(uint32_t leftPos, uint32_t rightPos) {
    nCigarOperations = 0;
    sizeSum = sizeSumConsuming = 0;
    leftGenPosSet = leftGenPos = leftPos;
    rightGenPosSet = rightGenPos = rightPos;
}

uint16_t ogCigarOperations::pushOperationsFromCIGAR(char *pCIGARchars) {
    // reset must be called before
    uint32_t n;
    while (*pCIGARchars) {
        for (n = 0; *pCIGARchars >= '0' && *pCIGARchars <= '9'; n = n*10 + *pCIGARchars - 48, pCIGARchars++);
        //fprintf(stderr, "%c", *pCIGARchars);
        if (*pCIGARchars) {
            push(*pCIGARchars++, n);
        }
    }
    return nCigarOperations;
}

char isSizeConsuming(char op) {
    return (op == 'M' || op == 'I' || op == 'X' || op == 'S');
}

void ogCigarOperations::push(char op, uint32_t len) {
    cigarOp[nCigarOperations].operation = op;
    cigarOp[nCigarOperations].size = len;
    sizeSum += len;
    if (isSizeConsuming(op)) sizeSumConsuming += len;
    nCigarOperations++;
}

OG_CIGAR_OPERATION *ogCigarOperations::pull() {
    return (cigarOp + nCigarOperations--);
}

OG_CIGAR_OPERATION *ogCigarOperations::last() {
    return (cigarOp + nCigarOperations);
}

OG_CIGAR_OPERATION *ogCigarOperations::getOperation(uint16_t x) {
    if (x < nCigarOperations ) return cigarOp + x;
    return NULL;
}


uint32_t ogCigarOperations::totalSize() {
    return sizeSum;
}

uint32_t ogCigarOperations::totalSizeConsuming() {
    return sizeSumConsuming;
}

uint32_t ogCigarOperations::buildCIGAR(char *pCigar, uint32_t maxCigarLen) {
    if (nCigarOperations == 0) {
        *pCigar = '*';
        *(pCigar+1) = 0;
        return 1;
    }
    int i;
    char prevOp = '-';
    uint32_t size = 0;
    uint32_t len = 0;
    OG_CIGAR_OPERATION *cop = cigarOp;

    for (i=0; len < maxCigarLen && i < nCigarOperations; i++, cop++) {
        len += snprintf(pCigar+len, maxCigarLen-len, "%d%c", cop->size, cop->operation);
    }
    return len;
    
    /**
    for (i=0; i < nCigarOperations; i++, cop++) {
        if (cop->operation == prevOp) {
            size += cop->size;
        } else {
            len += snprintf(pCigar+len, maxCigarLen-len, "%d%c", size, prevOp);
            size = cop->size;
            prevOp = cop->operation;
        }
    }
    if (prevOp != '-') {
       len += snprintf(pCigar+len, maxCigarLen-len, "%d%c", size, prevOp);
    }
    pCigar[len] = 0;
    return len;
     **/
}



uint32_t ogCigarOperations::sumSizeConsuming(uint16_t from, uint16_t to) {
    uint16_t            i;
    uint32_t            sum = 0;
    OG_CIGAR_OPERATION *cop = cigarOp + from;

    for (i=from; i < to; i++, cop++) {
        if (isSizeConsuming(cop->operation)) sum += cop->size;
    }
    return sum;
}
    

uint32_t ogCigarOperations::sumOperationsType(char op) {
    uint16_t            i;
    uint32_t            sum = 0;
    OG_CIGAR_OPERATION *cop = cigarOp;

    for (i=0; i < nCigarOperations; i++, cop++) {
        if (cop->operation == op) sum += cop->size;
    }
    return sum;
}
    



void ogCigarOperations::removeOperation(uint16_t x) {
    if (x < nCigarOperations) {
        sizeSum -= cigarOp[x].size;
        if (isSizeConsuming(cigarOp[x].operation)) sizeSumConsuming -= cigarOp[x].size;
        if (x < nCigarOperations-1) memmove(cigarOp+x, cigarOp+(x+1), sizeof(OG_CIGAR_OPERATION) * (nCigarOperations-x));
        nCigarOperations--;
    }
}

void ogCigarOperations::setOperation(uint16_t x, char op, uint32_t len) {
    if (x < nCigarOperations) {
        sizeSum -= cigarOp[x].size;
        if (isSizeConsuming(cigarOp[x].operation)) sizeSumConsuming -= cigarOp[x].size;
        cigarOp[x].operation = op;
        cigarOp[x].size = len;
        sizeSum += len;
        if (isSizeConsuming(op)) sizeSumConsuming += len;
    }
}


char ogCigarOperations::reduceOperations() {
    if (nCigarOperations < 2) return nCigarOperations;
    
    char        res = 1;
    char        prevOp = '-', nextOp = '-';
    char        isLast;
    char        isFirst;
    uint32_t    prevSize = 0, prevI, prevD, nextSize = 0;
    uint16_t    i;
    OG_CIGAR_OPERATION *op;
    uint16_t    lastM = 0;
    uint32_t    lastMsize = 0;
    uint16_t    upperM = 0;
    //uint32_t    ins = sumOperationsType('I');
    //uint32_t    del = sumOperationsType('D');
    uint32_t    length = sizeSumConsuming;
    //char        cigarro[10000];
    char shortify = 0;

    uint16_t maxTransitions = (sizeSumConsuming >> 4);
    if (maxTransitions < 3) maxTransitions = 3;
    //if (maxTransitions > 9) {
    //    fprintf(stderr,"------------ EPALE ------------ maxTransitions=%d, sizeSumConsuming=%d\n",maxTransitions, sizeSumConsuming);
    //}
    
    //fprintf(stderr, "--------------------------------------------\n"); fflush(stderr);
    
    //14M 85I 85D 27M 19D 2M 19I 4M 
    //32M 68D 113I 5M 45D 1X
    //3M 102D 6M 142I 40D
    //1M 3D 3M 2I 5M 1I 138M
    //1D 43M 1I 30M 2I 2M 2D 13M 5I 5D 11M 1D 4M 1I 3M 1D 3M 1I 32M ===> 43M 1I 30M 2I 2M 2D 29M 1D 4M 1I 3M 1D 3M 1I 32M

    lastM = lastMsize = upperM = 0;
    for (i=0; i < nCigarOperations; i++) {
        op = getOperation(i);
        if (op->operation == 'M') {
            if (lastMsize < op->size) {
                lastM = i;
                lastMsize = op->size;
            }
            if (op->size >= MIN_CIGAR_OPERATION) upperM++;
        }
    }
    if (upperM > 1) {
        while (res) {
            res = 0;
            lastM = lastMsize = upperM = 0;
            prevI = 0;
            prevD = 0;
            upperM = 0;
            for (i=0; i < nCigarOperations; i++) {
                //buildCIGAR(cigarro, 10000);
                //fprintf(stderr, "CIGAR:%s\n",cigarro); fflush(stderr);
                isLast = (i+1 == nCigarOperations);
                //fprintf(stderr, "Op %d/%d %c ==> ",i,nCigarOperations,isLast ? ' ' : ' '); fflush(stderr);
                if (isLast) {
                    nextOp = '-';
                    nextSize = 0;
                } else {
                    op = getOperation(i+1);
                    nextOp = op->operation;
                    nextSize = op->size;
                }
                if (i == 0) {
                    isFirst = 1;
                    prevOp = '-';
                    prevSize = 0;
                } else {
                    isFirst = 0;
                    op = getOperation(i-1);
                    prevOp = op->operation;
                    prevSize = op->size;
                }
                op = getOperation(i);
                //fprintf(stderr, "Prev <%c:%u> |", prevOp, prevSize); fflush(stderr);
                //fprintf(stderr, "Curr <%c:%u> |", op->operation, op->size); fflush(stderr);
                //fprintf(stderr, "Next <%c:%u>\n", nextOp, nextSize); fflush(stderr);

                if (op->operation == prevOp) {
                    // two equivalent operations, simplify to 1.
                    // asume i > 0 por el uso de prevOp
                    setOperation(i-1, prevOp, prevSize += op->size);
                    removeOperation(i);
                    res = 1;
                    i--;
                } else if ((op->operation == 'I' && prevOp == 'D' || op->operation == 'D' && prevOp == 'I') && prevSize == op->size) {
                    setOperation(i-1, 'X', prevSize);
                    removeOperation(i);
                    i--;
                    res = 1;
                } else if (op->operation == 'X' && prevOp == 'S'  ||  op->operation == 'S' && prevOp == 'X') {
                        if (prevOp == 'S') leftGenPos += op->size; // should be first
                        else if (prevOp == 'X') rightGenPos -= op->size; // should be last
                        setOperation(i-1, 'S', prevSize + op->size);
                        removeOperation(i);
                        i--;
                        res = 1;
                } else if (op->operation == 'X' && isFirst) {
                        leftGenPos += op->size;
                        setOperation(i, 'S', op->size);
                        i--;
                        res = 1;
                } else if (op->operation == 'X' && isLast) {
                        rightGenPos -= op->size;
                        setOperation(i, 'S', op->size);
                        i--;
                        res = 1;
                } else if (op->operation == 'I' && i == 1 && prevSize < MIN_CIGAR_OPERATION && (prevOp == 'M' || prevOp == 'X') && op->size < MIN_CIGAR_OPERATION) {
                    // insercion 'chica' "casi" al inicio
                    // Convertir en "match" porque es una correccion del WFA debido a un posible delete en medio
                    if (leftGenPos >= op->size) {
                        leftGenPos -= op->size;
                        setOperation(i, prevOp, op->size);
                        i--;
                        res = 1;
                    }
                } else if (op->operation == 'D' && i == 1 && op->size < MIN_CIGAR_OPERATION && prevSize < MIN_CIGAR_OPERATION && (prevOp == 'M' || prevOp == 'X')) {
                    // delecion 'chica' "casi" al inicio
                    leftGenPos += op->size;
                    removeOperation(i);
                    i--;
                    res = 1;
                } else if (op->operation == 'I') {
                    prevI = op->size;                        
                    if (prevI >= MIN_CIGAR_OPERATION && (prevD > prevI ? prevD-prevI : prevI-prevD) < 3) {
                        // hay inserción grande compensando deleción grande, ==> S
                        shortify = 1;
                    }
                    if (shortify == 0) {
                        if (isFirst) {
                            //setOperation(i, 'S', op->size);
                            //i--;
                            if (leftGenPos >= op->size) {
                                leftGenPos -= op->size;
                                setOperation(i, op->size > MIN_CIGAR_OPERATION ? 'X' : 'M', op->size);
                            } else {
                                // Esta al inicio del cromosoma (rarísimo)
                                setOperation(i, 'S', op->size);
                                i--;
                            }
                            res = 1;
                        } else if (isLast) {
                            //setOperation(i, 'S', op->size);
                            //i--;
                            rightGenPos += op->size;
                            setOperation(i, op->size > MIN_CIGAR_OPERATION ? 'X' : 'M', op->size);
                            i--;
                            res = 1;                    
                        }
                    }
                } else if (op->operation == 'D') {
                    prevD = op->size;
                    if (prevD >= MIN_CIGAR_OPERATION && (prevD > prevI ? prevD-prevI : prevI-prevD) < 3) {
                        // hay delecion grande compensando insercion grande, ==> S
                        shortify = 1;
                    }
                    if (shortify == 0) {
                        if (isLast) {
                            rightGenPos -= op->size;
                            removeOperation(i);
                            i--;
                            res = 1;
                        } else if (isFirst) {
                            leftGenPos += op->size;
                            removeOperation(i);
                            i--;
                            res = 1;
                        }
                    }
                } else if ((op->operation == 'M' || op->operation == 'X') && isLast && op->size < MIN_CIGAR_OPERATION && (prevOp == 'I' || prevOp == 'D') && prevSize < MIN_CIGAR_OPERATION) {
                    if (prevOp == 'I') {
                        // es una inserción 'chica' "casi" al final
                        rightGenPos += prevSize;
                        setOperation(i-1, 'X', prevSize);
                        i--;
                        res = 1;
                    } else if (prevOp == 'D') {
                        rightGenPos -= prevSize;
                        removeOperation(i-1);
                        i--;
                        res = 1;
                    }
                } else if ((op->operation == 'X' && prevOp == 'M' || op->operation == 'M' && prevOp == 'X') && op->size < MIN_CIGAR_OPERATION && prevSize < MIN_CIGAR_OPERATION) {
                    // asume i > 0 por el uso de prevOp
                    setOperation(i-1, prevSize >= op->size ? prevOp : op->operation, prevSize + op->size);
                    removeOperation(i);
                    res = 1;
                    i--;
                } else if ((op->operation == 'M' || op->operation == 'X') && prevSize == nextSize && (prevOp == 'I' && nextOp == 'D' || prevOp == 'D' && nextOp == 'I') &&  op->size+prevSize < MIN_CIGAR_OPERATION) {
                    // ... I M D ...
                    // ... D M I ...
                    setOperation(i-1, 'M', prevSize + op->size);
                    removeOperation(i+1);
                    removeOperation(i);
                    i--;
                    res = 1;
                } else if (op->operation == 'M') {
                    if (lastMsize < op->size) {
                        lastM = i;
                        lastMsize = op->size;
                    }
                    if (op->size >= MIN_CIGAR_OPERATION) upperM++;
                }
            }
        }
    } else {
        shortify = 1;
    }

    uint32_t left = 0;
    uint32_t right = 0;
    if (nCigarOperations > maxTransitions && (upperM == 1 || lastMsize > (sizeSumConsuming >> 1))) {
        shortify = 1;
    }
    if (upperM == 0 || upperM == 1 && nCigarOperations > 3) shortify = 1;

    if (shortify) {
        if (lastM > 0) left = sumSizeConsuming(0, lastM);
        if (lastM < nCigarOperations-1) right = sumSizeConsuming(lastM+1, nCigarOperations);
        reset(leftGenPosSet, rightGenPosSet);
        if (left+lastMsize+right == length) {
            if (left > 0) { 
                push('S', left); 
                leftGenPos += left;
            }
            push('M', lastMsize);
            if (right > 0) {
                push('S', right);
                rightGenPos -= right;
            }
        } else {
            push('M', length); // problems, mismatches or rare things, then leave them as is
        }
    }
    //fprintf(stderr, "++++++++\n"); fflush(stderr);
    
    return nCigarOperations;
}
