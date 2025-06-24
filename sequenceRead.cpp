/* 
 * File:   sequenceRead.cpp
 * Author: victortrevino
 * 
 * Created on June 9, 2022, 6:21 PM
 */

#include "sequenceRead.h"
#include "ogDefinitions.hpp"

sequenceRead::sequenceRead(char *pid, uint32_t pidLen, char *pseq, uint32_t pseqLen, char *pqual, uint32_t pqualLen, uint64_t readNum) {
    setSequenceRead(pid, pidLen, pseq, pseqLen, pqual, pqualLen, readNum,
                 NULL, 0,     NULL, 0,       NULL,  0);
}

sequenceRead::sequenceRead(
                       char *pid1, uint32_t pidLen1, char *pseq1, uint32_t pseqLen1, char *pqual1, uint32_t pqualLen1, uint64_t readNum,
                       char *pid2, uint32_t pidLen2, char *pseq2, uint32_t pseqLen2, char *pqual2, uint32_t pqualLen2) {
    setSequenceRead(pid1, pidLen1, pseq1, pseqLen1, pqual1, pqualLen1, readNum,
                    pid2, pidLen2, pseq2, pseqLen2, pqual2, pqualLen2);
}

void sequenceRead::setSequenceRead(
                       char *pid1, uint32_t pidLen1, char *pseq1, uint32_t pseqLen1, char *pqual1, uint32_t pqualLen1, uint64_t readNum,
                       char *pid2, uint32_t pidLen2, char *pseq2, uint32_t pseqLen2, char *pqual2, uint32_t pqualLen2) {

    //    id = std::string(pid, pidLen);
//    sequence = std::string(pseq, pseqLen);
//    quality = std::string(pqual, pqualLen);
    isPaired = (pid2 != NULL);
    idLen1 = pidLen1;
    seqLen1 = pseqLen1;
    quaLen1 = pqualLen1;
    readNumber = readNum;
    if (*(pid1+pidLen1) != 0 || *(pseq1+pseqLen1) != 0 || *(pqual1+pqualLen1) != 0) {
        fprintf(stderr, "**** PROBLEMAS CON LOS FINES DE STRING READ single ***\n");
    }
    uint64_t sum = pidLen1 + (pseqLen1 << 1) + (pqualLen1 << 1) + 6 + MAX_CIGAR_SIZE;
    if (isPaired) {
        idLen2 = pidLen2;
        seqLen2 = pseqLen2;
        quaLen2 = pqualLen2;
        sum += pidLen2 + (pseqLen2 << 1) + (pqualLen2 << 1) + 6 + MAX_CIGAR_SIZE;
        if (*(pid2+pidLen2) != 0 || *(pseq2+pseqLen2) != 0 || *(pqual2+pqualLen2) != 0) {
            fprintf(stderr, "**** PROBLEMAS CON LOS FINES DE STRING READ PAIRED ***\n");
        }
    }

    if (sum < 10 || sum > 100000) {
        fprintf(stderr, "**** PROBLEMAS CON LAS LONGITUDES %llu ***\n", sum);
    }
//    if (idLen != strlen(pId) || seqLen != strlen(pseq) || pqualLen != strlen(pqual)) {
//        fprintf(stderr, "**** PROBLEMAS CON LAS LONGITUDES ***\n");
//    }
    p = (char *) malloc(sum * sizeof(char));
    pId1   = p;
    pQual1 = pId1 + pidLen1 + 1;
    pSeq1  = pQual1 + pqualLen1 + 1;
    pCigar1 = pSeq1 + pseqLen1 + 1;
    pSeqRevComp1 = pCigar1 + MAX_CIGAR_SIZE + 1;
    pQualRev1 = pSeqRevComp1 + pseqLen1 + 1;
    
    strncpy(pId1,   pid1, pidLen1);
    pId1[pidLen1] = 0;
    strncpy(pQual1, pqual1, pqualLen1);
    pQual1[pqualLen1] = 0;
    strncpy(pSeq1,  pseq1, pseqLen1);
    pSeq1[pseqLen1] = 0;
    *pCigar1 = 0;
    uint32_t i;
    char *p1 = pSeq1;
    char *p2 = pQualRev1 - 1;
    char *p3 = pQualRev1 + pseqLen1;
    char *p4 = pQual1;
    *p2-- = 0;    
    *p3-- = 0;
    for (i=0; i < seqLen1; i++) {
        *p2-- = ascii2Complement[*p1++];
        *p3-- = *p4++;
    }

    if (isPaired) {
        pId2 = pQualRev1 + pqualLen1 + 1;
        pQual2 = pId2 + pidLen2 + 1;
        pSeq2 = pQual2 + pqualLen2 + 1;
        pCigar2 = pSeq2 + pseqLen2 + 1;
        pSeqRevComp2 = pCigar2 + MAX_CIGAR_SIZE + 1;
        pQualRev2 = pSeqRevComp2 + pseqLen2 + 1;
        
        strncpy(pId2,   pid2, pidLen2);
        pId2[pidLen2] = 0;
        strncpy(pQual2, pqual2, pqualLen2);
        pQual2[pqualLen2] = 0;
        strncpy(pSeq2,  pseq2, pseqLen2);
        pSeq2[pseqLen2] = 0;
        *pCigar2 = 0;
        p1 = pSeq2;
        p2 = pQualRev2 - 1;
        p3 = pQualRev2 + pseqLen2;
        p4 = pQual2;
        *p2-- = 0;    
        *p3-- = 0;
        for (i=0; i < seqLen2; i++) {
            *p2-- = ascii2Complement[*p1++];
            *p3-- = *p4++;
        }
    }

}

sequenceRead::~sequenceRead() {
    //fprintf(stderr, "/");fflush(stderr);
    free(p);
    //fprintf(stderr, "|");fflush(stderr);
}


void sequenceRead::setSamInfo(char isRead2, char *_qname, int _flag, char *_rname, uint64_t _pos, int _mapq, char *_cigar, char *_rnext, uint64_t _posnext, int _tlen, char *_seq, char *_qual) {
    ogSAM *pSAM = (isRead2 ? &samRecord2 : &samRecord1);
    pSAM->qname     = _qname;
    pSAM->flag      = _flag;
    pSAM->rname     = _rname;
    pSAM->pos       = _pos;
    pSAM->mapq      = _mapq;
    pSAM->cigar     = _cigar;
    pSAM->rnext     = _rnext;
    pSAM->pnext     = _posnext;
    pSAM->tlen      = _tlen;
    pSAM->seq       = _seq;
    pSAM->qual      = _qual;
}

/* 
    char        *qname;     // Query template NAME
    int         flag;       // bitwise FLAG
    char        *rname;     // References sequence NAME
    uint64_t    pos;        // 1- based leftmost mapping POSition
    int         mapq;       // MAPping Quality
    char        *cigar;     // CIGAR string
    char        *rnext;     // Ref. name of the mate/next read
    uint64_t    pnext;      // Position of the mate/next read
    int         tlen;       // observed Template LENgth
    char        *seq;       // segment SEQuence
    char        *qual;      // ASCII of Phred-scaled base QUALity+33
 */