/* 
 * File:   sequenceRead.h
 * Author: victortrevino
 *
 * Created on June 9, 2022, 6:21 PM
 */

#ifndef SEQUENCEREAD_H
#define SEQUENCEREAD_H

//#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "ogDefinitions.hpp"



class sequenceRead {
    
public:
    char isPaired;
    char *p;                                    // ALL DATA :: std::string id, sequence, quality;
    char *pId1, *pSeq1, *pQual1, *pSeqRevComp1, *pQualRev1, *pCigar1;
    char *pId2, *pSeq2, *pQual2, *pSeqRevComp2, *pQualRev2, *pCigar2;
    uint32_t quaLen1, seqLen1, idLen1;
    uint32_t quaLen2, seqLen2, idLen2;
    uint64_t readNumber;
    uint32_t genPos1_Left, genPos1_Right;
    uint32_t genPos2_Left, genPos2_Right;
    uint32_t genPos_leftMax, genPos_rightMin;   // if leftMost > rightMost : Overlap between reads
    uint32_t genPos_leftMin, genPos_rightMax;   // if leftMost > rightMost : Overlap between reads
    uint32_t genPos_difference;                 // if leftMost < rightMost : Gap between reads
    uint32_t genPos_insertSize;                 //
    char     genPos1_mapScore;
    char     genPos2_mapScore;
    char     genPos1_isRd2;
    char     posIsReverse1, posIsReverse2;
    uint16_t leftSoftR1; // se agrega para que compile, esta clase esta depracada como sea
    uint16_t leftSoftR2; // se agrega para que compile, esta clase esta depracada como sea
    
    ogSAM    samRecord1;
    ogSAM    samRecord2;
    sequenceRead(char *pid, uint32_t pidLen, char *pseq, uint32_t pseqLen, char *pqual, uint32_t pqualLen, uint64_t readNum);
    sequenceRead(char *pid1, uint32_t pidLen1, char *pseq1, uint32_t pseqLen1, char *pqual1, uint32_t pqualLen1, uint64_t readNum,
                 char *pid2, uint32_t pidLen2, char *pseq2, uint32_t pseqLen2, char *pqual2, uint32_t pqualLen2);
    virtual ~sequenceRead();
    void        setSequenceRead(
                       char *pid1, uint32_t pidLen1, char *pseq1, uint32_t pseqLen1, char *pqual1, uint32_t pqualLen1, uint64_t readNum,
                       char *pid2, uint32_t pidLen2, char *pseq2, uint32_t pseqLen2, char *pqual2, uint32_t pqualLen2);
    void        setSamInfo(char isRead2, char *_qname, int _flag, char *_rname, uint64_t _pos, int _mapq, char *_cigar, char *_rnext, uint64_t _posnext, int _tlen, char *_seq, char *_qual);
    
private:

};

#endif /* SEQUENCEREAD_H */

