
/* 
 * File:   ogGenome.hpp
 * Author: victortrevino
 *
 * Created on May 4, 2022, 7:03 PM
 */

#ifndef OGGENOME_HPP
#define OGGENOME_HPP

#include <stdio.h>
#include <zlib.h>

#include "ogGuider.hpp"
#include "ogKeyEncoding.hpp"

#define ALL_ASCII           256
#define MAX_CHROMOSOME_NAME 128
#define MAX_CHROMOSOME_COMMENT 512

typedef struct ogChrPosWithN {
    uint32_t        position;
    uint8_t         length;
} OG_CHR_POS_WITH_N;


typedef struct ogChromosome {
    uint32_t        index;
    uint32_t        number;
    char            name[MAX_CHROMOSOME_NAME];
    char            comment[MAX_CHROMOSOME_COMMENT];
    uint32_t        size;
    uint32_t        cummulative;
    uint32_t        start;
    uint32_t        initialPack;
    char            initialBytePos;
    uint32_t        validPositions;
    uint32_t        validPositionsOverLen;
    uint32_t        nPositionsWithN;
    ogChrPosWithN   *pPosWithN;
    uint32_t        effectiveSize; // used only in counting
    uint32_t        nExons;        // used only in counting
} OG_INDEX_CHROMOSOME;

class ogGenome {
    
    uint16_t            keySize;
    ogChromosome       *allChromosomes;
    uint32_t            nChromosomes;
    uint32_t            chr_i;
    char               *sourceFileName;
    char               *fileName;
    unsigned char      *packed;
    unsigned char      *packLimit;
    uint32_t            genomeSize;
    uint32_t            memoryPacked;
    gzFile              fastaGZ;
    unsigned char       ascii2pack[ALL_ASCII];
    unsigned char       ascii2pack1[ALL_ASCII];
    unsigned char       ascii2pack2[ALL_ASCII];
    unsigned char       ascii2pack3[ALL_ASCII];
    void               *seqPointer;

    //unsigned char      *pack_i;
    //unsigned char       pack_s;
    uint32_t            validPositions;
    uint32_t            validPositionsOverLen;
    
public:
                    ogGenome(char *sourceFileName, uint16_t kSize);
    virtual        ~ogGenome();
    void            setSourceFile(char *filename);
    void            setKeySize(uint16_t kSize);
    void            estimateChromosomeSizes(char verbose);
    void            openSourceFile();
    void            closeSourceFile();
    char            *readSourceFile();
    uint32_t        getGenomeSize();
    uint32_t        getGenomeValidPositions();
    uint32_t        getGenomeValidPositionsOverLen();
    uint32_t        getNChromosomes();
    uint32_t        getChromosomeSize(uint32_t iChromosome);
    ogChromosome   *getLargestChromosome();
    uint32_t        getChromosomeStart(uint32_t iChromosome);
    uint32_t        getChromosomeValidPositions(uint32_t iChromosome);
    ogChromosome   *getChromosome(uint32_t i); // 0 - based array
    ogChromosome   *getGenomicCoordinate(uint32_t genPos);
    void            loadChromosomesInfo(FILE *pFile);
    void            loadPackedGenome(FILE *pFile);
    uint32_t        saveChromosomesInfo(FILE *pFile);
    uint32_t        savePackedGenome(FILE *pFile);
    void            allocate();
    //char            extractNextNtFromPackedChromosome();
    void            extractPackedGenome(uint32_t genpos, uint32_t size, char *str, char revComp, char lowerCase);
    void            extractPackedChromosome(uint32_t chr, uint32_t pos, uint32_t size, char *str, char revComp, char lowerCase);
    void            extractFromGenomicPos(uint32_t gPos, uint32_t size, char *str, char revComp, char lowerCase);
    void            extractFromPack_i_s(uint32_t size, char *str, char revComp, char lowerCase, unsigned char *pack_i, unsigned char pack_s);
    void            packGenome(uint32_t nNtInform, ogKeyEncoding *pEncoding, ogGuider *pGuider);
    void            printGenomePositionInfo(uint32_t gp);
    void            packSequence(char *seq, char *pack, char offset);
    void            packSequences(char *seq, char *pack0, char *pack1, char *pack2, char*pack3);
    unsigned char   getPackedOffset(uint32_t gPos);
    uint32_t        comparePackedSequences(char *pack, uint32_t gPos, uint32_t len, uint32_t maxErrors);
    void            setNPosToUnpackChr(uint32_t chr, uint32_t chrPos, uint32_t size, char *unpackstr);
    void            printFromGenomicPos(uint32_t gPos, uint32_t size);

    uint32_t        nTotalPositionsWithN;
    
private:

};

#endif /* OGGENOME_HPP */

