/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.cc to edit this template
 */

/* 
 * File:   ogSingleRead.cpp
 * Author: victortrevino
 * 
 * Created on February 11, 2023, 6:08 PM
 */

#include "ogSingleRead.hpp"

ogSingleRead::ogSingleRead() {
    // char *pid, uint32_t pidLen, char *pseq, uint32_t pseqLen, char *pqual, uint32_t pqualLen, uint64_t readNum, char pack
    memoryAllocated = 0;
    pId = NULL;
    pGenome = NULL;
}

ogSingleRead::~ogSingleRead() {
    // This part of the code is obsolete, so both reads must be deleted.
    // NOPE, siempre si
    if (isPaired) {
        pairedRead->isPaired =  0;
        delete pairedRead; // this adds the functionality of automatically deleting the paired read, so only 1 read must be deleted
    }
    
    free(pId);
}

void ogSingleRead::setRead(char *pid, uint32_t pidLen, char *pseq, uint32_t pseqLen, char *pqual, uint32_t pqualLen, uint64_t readNum, char pack) {
    uint32_t LEN2 = (pseqLen << 1); // seqlen should be the same than pqualen
    uint64_t sum = 
            pidLen      // id just stored as is
            + LEN2      // Fwd & Rev
            + LEN2      // Quality Fwd & Rev
            + 6         // Delimiters (SeqFwd,SeqRev,QualFwd,QualRev,Id,Cigar)
            + MAX_CIGAR_SIZE // CIGAR data
            + MAX_CIGAR_SIZE + 1 // original
            + LEN2 + 16 // Packed sequence in the 4 possible byte frames needs twice seqLen and 2 bytes for possible offsets
            ;

    char *p;
    if (sum > memoryAllocated) {
        memoryAllocated = sum + (sum >> 3);   // 12% more
        p = (char *) realloc(pId, memoryAllocated * sizeof(char));
    } else {
        p = pId;
    }
    lenId = pidLen;
    lenSeq = pseqLen;
    lenQual = pqualLen;
    readNumber = readNum;
    pId   = p;
    pQual = pId + pidLen + 1;
    pSeq  = pQual + pqualLen + 1;
    pCigar = pSeq + pseqLen + 1;
    pCigarOriginal = pCigar + MAX_CIGAR_SIZE + 1;
    pSeqRevComp = pCigarOriginal + MAX_CIGAR_SIZE + 1;
    pQualRev = pSeqRevComp + pseqLen + 1;
    //cigarLeftType = '-';
    //cigarLeft = 0;
    //cigarRight = 0;
    cigarLeftPos = cigarRightPos = 0;
    cigarMatches = cigarNoMatches = cigarSoft = cigarIns = cigarDel = 0;

    pairedRead = NULL;
    
    packNeeded = pack;

    /**
    int _pidlen = strlen(pid);
    int _psqlen = strlen(pseq);
    int _pqulen = strlen(pqual);
    if (_pidlen != pidLen || _psqlen != pseqLen || _pqulen != pqualLen ) {
        fprintf(stderr, "LONGITUDES PROBLEMAS id %d vs %d, seq %d vs %d, qual %d vs %d\n ", _pidlen, pidLen, _psqlen, pseqLen, _pqulen, pqualLen); fflush(stderr);
    }
    **/
    
    strncpy(pId, pid, pidLen);
    pId[pidLen] = 0;
    //Copy is better done in the for below
    //strncpy(pQual, pqual, pqualLen);
    //pQual[pqualLen] = 0;
    //strncpy(pSeq, pseq, pseqLen);
    //pSeq[pseqLen] = 0;
    *pCigar = 0;
    *pCigarOriginal = 0;
    uint32_t i;
    char *pseqP = pseq, *pqualP = pqual;
    char *p1 = pSeq;
    char *p2 = pQualRev - 1;
    char *p3 = pQualRev + pqualLen;
    char *p4 = pQual;
    *p2-- = 0;
    *p3-- = 0;
    Ncount = 0;
    for (i=0; i < lenSeq; i++) {
        //Copy
        if ((*p1 = *pseqP++) == 'N') Ncount++;
        *p4 = *pqualP++;
        *p2-- = ascii2Complement[*p1++];
        *p3-- = *p4++;
    }
    *p1 = 0; // finish strings
    *p4 = 0;
    isPaired = 0;
    //leftSoft = 0;
    readIndex = 1;
    packedInitial = pQualRev + lenQual + 1;
    packedReadSeq[0] = NULL;
    packedReadSeq[1] = NULL;
    packedReadSeq[2] = NULL;
    packedReadSeq[3] = NULL;
    packedReadSeq[4] = NULL;
    packedReadSeq[5] = NULL;
    packedReadSeq[6] = NULL;
    packedReadSeq[7] = NULL;
    packLen = (lenSeq >> 2) + 2;
}

void ogSingleRead::packIfNeeded(ogGenome *pGenomeObject) {
    pGenome = pGenomeObject;
    if (packNeeded) {
        packNeeded = 0;
//        uint32_t packLen = (lenSeq >> 2) + 2;
//        packedReadSeqFwd[0] = packedInitial;
//        packedReadSeqFwd[1] = packedReadSeqFwd[0] + packLen;
//        packedReadSeqFwd[2] = packedReadSeqFwd[1] + packLen;
//        packedReadSeqFwd[3] = packedReadSeqFwd[2] + packLen;
//        packedReadSeqRev[0] = packedReadSeqFwd[3] + packLen;
//        packedReadSeqRev[1] = packedReadSeqRev[0] + packLen;
//        packedReadSeqRev[2] = packedReadSeqRev[1] + packLen;
//        packedReadSeqRev[3] = packedReadSeqRev[2] + packLen;
//        pGenome->packSequence(pSeq, packedReadSeqFwd[0], 0);
//        pGenome->packSequence(pSeq, packedReadSeqFwd[1], 1);
//        pGenome->packSequence(pSeq, packedReadSeqFwd[2], 2);
//        pGenome->packSequence(pSeq, packedReadSeqFwd[3], 3);
//        pGenome->packSequence(pSeqRevComp, packedReadSeqRev[0], 0);
//        pGenome->packSequence(pSeqRevComp, packedReadSeqRev[1], 1);
//        pGenome->packSequence(pSeqRevComp, packedReadSeqRev[2], 2);
//        pGenome->packSequence(pSeqRevComp, packedReadSeqRev[3], 3);
        packedReadSeq[0] = packedInitial;
        packedReadSeq[1] = packedReadSeq[0] + packLen;
        packedReadSeq[2] = packedReadSeq[1] + packLen;
        packedReadSeq[3] = packedReadSeq[2] + packLen;
        pGenome->packSequences(pSeq, packedReadSeq[0], packedReadSeq[1], packedReadSeq[2], packedReadSeq[3]);
        packedReadSeq[4] = packedReadSeq[3] + packLen;
        packedReadSeq[5] = packedReadSeq[4] + packLen;
        packedReadSeq[6] = packedReadSeq[5] + packLen;
        packedReadSeq[7] = packedReadSeq[6] + packLen;
        pGenome->packSequences(pSeqRevComp, packedReadSeq[4], packedReadSeq[5], packedReadSeq[6], packedReadSeq[7]);            

    }
}

void ogSingleRead::setPairedRead(ogSingleRead *theRead) {
    isPaired = 1;
    pairedRead = theRead;
    readIndex = 1;
    theRead->isPaired = 1;
    theRead->pairedRead = this;
    theRead->readIndex = 2;
}

char *ogSingleRead::packReadSeqIfNeeded(char offset, char isRev) {
    char i = (isRev ? offset + 4 : offset);
//    if (packedReadSeq[i] == NULL) {
//        //packedReadSeq[i] = packedInitial + (packLen * i);
//        //pGenome->packSequence(isRev ? pSeqRevComp : pSeq, packedReadSeq[i], offset);
//        packedReadSeq[0] = packedInitial;
//        packedReadSeq[1] = packedReadSeq[0] + packLen;
//        packedReadSeq[2] = packedReadSeq[1] + packLen;
//        packedReadSeq[3] = packedReadSeq[2] + packLen;
//        pGenome->packSequences(pSeq, packedReadSeq[0], packedReadSeq[1], packedReadSeq[2], packedReadSeq[3]);
//        packedReadSeq[4] = packedReadSeq[3] + packLen;
//        packedReadSeq[5] = packedReadSeq[4] + packLen;
//        packedReadSeq[6] = packedReadSeq[5] + packLen;
//        packedReadSeq[7] = packedReadSeq[6] + packLen;
//        pGenome->packSequences(pSeqRevComp, packedReadSeq[4], packedReadSeq[5], packedReadSeq[6], packedReadSeq[7]);            
//    }
    return packedReadSeq[i];
}
