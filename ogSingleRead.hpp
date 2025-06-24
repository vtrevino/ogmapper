/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.h to edit this template
 */

/* 
 * File:   singleRead.hpp
 * Author: victortrevino
 *
 * Created on February 11, 2023, 6:08 PM
 */

#ifndef OGSINGLEREAD_HPP
#define OGSINGLEREAD_HPP

//#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "ogDefinitions.hpp"
#include "ogGenome.hpp"


class ogSingleRead {
    
public:

    uint64_t        readNumber;
    char            *pId, *pSeq, *pQual, *pSeqRevComp, *pQualRev, *pCigar;
    uint32_t        lenQual, lenSeq, lenId;
    uint32_t        Ncount;
    uint32_t        memoryAllocated;
    //char            *packedReadSeqFwd[4];
    //char            *packedReadSeqRev[4];
    char            *packedInitial;
    char            *packedReadSeq[8]; // 4 Fwd + 4 Rev
    char            packNeeded;
    char            isPaired;
    char            readIndex; // 1, 2
    //char            cigarLeftType;
    //uint16_t        cigarLeft;
    //uint16_t        cigarRight;
    uint32_t        cigarLeftPos;
    uint32_t        cigarRightPos;
    uint32_t        cigarMatches;
    uint32_t        cigarNoMatches;
    uint32_t        cigarSoft;
    uint32_t        cigarIns;
    uint32_t        cigarDel;
    uint32_t        packLen;
    ogSingleRead    *pairedRead;
    char            *pCigarOriginal;
    ogGenome        *pGenome;
    
                ogSingleRead();
    virtual     ~ogSingleRead();
    void        setPairedRead(ogSingleRead *theRead);
    void        packIfNeeded(ogGenome *pGenomeObject);
    void        setRead(char *pid, uint32_t pidLen, char *pseq, uint32_t pseqLen, char *pqual, uint32_t pqualLen, uint64_t readNum, char pack);
    char        *packReadSeqIfNeeded(char offset, char isRev);
    
private:

};

#endif /* OGSINGLEREAD_HPP */
