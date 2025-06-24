/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#ifndef _OG_INDEX_H
#define _OG_INDEX_H

#include <stdint.h>
#include <cstdlib>
#include <zlib.h>
//#include "kseq.h"
//KSEQ_INIT(gzFile, gzread)

const uint16_t MAX_OGINDEX_FILENAME = 512;
const uint16_t ALL_ASCII = 256;
const uint16_t MAX_SIMILAR_POSITIONS = 128;
const uint16_t MAX_CHROMOSOME_NAME = 128;

typedef struct ogTupleChunk {
    uint32_t        size;
    uint32_t        offset;
} OG_TUPLE_CHUNK;


typedef struct ogIndexOriginalHeader {
    char            indexVersion[32];
    char            fileName[MAX_OGINDEX_FILENAME];      // Source File name
    char            ascii2code[ALL_ASCII];               // Index code values
    uint32_t        indexSize;                           // Maximum Genome Size: 4,294,967,296
    uint8_t         tupleLength;                         // in bytes, default to 3 
    uint8_t         tupleSizeInChars;                    // number of chars per tuple, just informative
    uint16_t        tupleChunkSize;                      // default allocation size for ogTupleChunks
    uint8_t         indexType;                           // Index Type, default to 0
    uint32_t        tupleArraySize;                      // Size of all tuples
    uint64_t        memoryPositions;                     // memory in bytes needed for positions table
    uint64_t        memoryTuples;                        // memory in bytes needed for tuples
    uint64_t        memoryPacked;                        // memory in bytes needed for packed genome
    uint32_t        largestSequence;                     // size of the largest chromosome
    OG_TUPLE_CHUNK  (*pChunks)[];                        // Pointer to an array of chunks (saved but unreal)
} OG_INDEX_HEADER;


typedef struct ogChromosome {
    int             number;
    char            seqName[MAX_CHROMOSOME_NAME];
    uint32_t        size;
    uint32_t        cummulative;
    uint32_t        start;
    ogChromosome    *next;
    ogChromosome    *prev;
} OG_INDEX_CHROMOSOME;

typedef struct ogSequenceKeys {
    uint32_t        *idxForward;        // Array of Forward Keys [0..MAX_READ_LEN-1]
    uint32_t        *idxReverse;        // Array of Reverse Keys [0..MAX_READ_LEN-1]
    OG_TUPLE_CHUNK  *chunkRoot;         // Pointer to initial Tuple Chunks
    char            completed;          // 100 is completed, 1 is partially set to positions
    OG_TUPLE_CHUNK  *tupTabForward;     // Arrays of Forward Tuples corresponding to Keys [0...MAX_READ_LEN-1]
    OG_TUPLE_CHUNK  *tupTabReverse;     // Arrays of Reverse Tuples corresponding to Keys [0...MAX_READ_LEN-1]
    uint16_t        nKeys;
    ogIndexOriginalHeader   header;// esta se debe de quitar
} OG_SEQUENCE_KEYS;

#define MAX_OG_MATCHES  32

typedef struct ogMatches {
    //kseq_t      *kseq;
    uint8_t     nMatches;
    uint8_t     maxMatch;
    uint32_t    position[MAX_OG_MATCHES];
    char        isRevComp[MAX_OG_MATCHES];
    uint16_t    score[MAX_OG_MATCHES];
    char        status; // 0 : not fully matched, 1 : fully matched
} OG_MATCHES;

#define initMatches(ogM) ogM.nMatches=ogM.maxMatch=0; ogM.status = 0;
#define pushMatch(ogM, pos, revcomp, lescore) if (ogM.nMatches < MAX_OG_MATCHES) { ogM.position[ogM.nMatches]=pos; ogM.isRevComp[ogM.nMatches]=revcomp; ogM.score[ogM.nMatches]=lescore; ogM.nMatches++; }
#define setMatchStatus(ogM, lestatus) ogM.status = lestatus;

class ogIndexOriginal {
    OG_INDEX_HEADER     header;
public:
    ogIndexOriginal();
    ~ogIndexOriginal();
    void        setMemoryMode(char lowMemMode);
    int         generateIndex(char *filename, int keySize);
    uint32_t    getKey(char *p);
    void        shiftKey(char * &p, uint32_t &key);
    int         loadIndex(char *filename, int keySize);
    int         mapSequences(char *filename);
    void        intersectPositions(uint32_t *simPositions, uint16_t &simPosK, OG_TUPLE_CHUNK &a, OG_TUPLE_CHUNK &b, uint32_t diff);
    void        setCodeMasks();
    void        setSchedule(char *s);
    void        unPackGenomePosition(char *p, uint32_t pos, uint32_t size);
    void        setMaxReadLength(uint16_t len);
    void        setMaxReads(uint64_t mxreads);
    void        allocateKeys(ogSequenceKeys &readKeys, ogIndexOriginalHeader &header);
    void        deallocateKeys(ogSequenceKeys &readKeys);
    void        getKeysFR(char *p, uint32_t &key1, uint32_t &key2);
    void        shiftKeysFR(char * &p, uint32_t &key1, uint32_t &key2);
    OG_INDEX_CHROMOSOME *findChromosome(uint32_t &p);
    OG_INDEX_CHROMOSOME *findChromosomeSequential(uint32_t &p);
};

#endif

