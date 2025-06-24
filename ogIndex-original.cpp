/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <ctype.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <cstdlib>
#include <stdlib.h>
#include <math.h>
#include <zlib.h>
#include <time.h>
#include "ogIndex-original.hpp"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define MAX_SCHEDULE    100

char xfile[MAX_OGINDEX_FILENAME];
FILE *pfile;
char *ALLPOS = NULL;
char *PACKED = NULL;
char ascii2pack[ALL_ASCII];
char schedule[MAX_SCHEDULE] = "ED";
OG_INDEX_CHROMOSOME     *firstChromosome = NULL;
OG_INDEX_CHROMOSOME*     *allChromosomes;
int  nChromosomes;


uint32_t simPositions[MAX_SIMILAR_POSITIONS];
uint16_t simPosK;
uint32_t keyCode   = 0x00800000; // default
uint32_t keyCodeRC = 0x00000001;
uint32_t maskCode = 0x00FFFFFF;
const uint32_t code1 = 0x00000001;

char lowMemoryMode = 0;

uint16_t    MAX_READ_LEN = 1000;
uint64_t    MAX_READS_MAP = 0xFFFFFFFFFFFFFFFF;


// Estimate the direct (5'->3') key in key1
// and the reverse complement key (3'->5') in key2
void getKeysFR(char *p, uint32_t &key1, uint32_t &key2, ogIndexOriginalHeader &header) {
    key1 = key2 = 0;
    uint32_t    code1 = keyCode;
    uint32_t    code2 = keyCodeRC;
    uint8_t     i;
    for (i=header.tupleSizeInChars; i; i--, p++) {
        if (header.ascii2code[*p])  {
            // taking advantage that code at *p is the same that the complement
            key1 |= code1;
            key2 |= code2;
        }
        code1 >>= 1;
        code2 <<= 1;        
    }
}


void shiftKeysFR(char * &p, uint32_t &key1, uint32_t &key2, ogIndexOriginalHeader &header) {
    key1 <<= 1;
    key2 >>= 1;
    //if (header.ascii2code[*(p+header.tupleSizeInChars_1)]) key |= code1;
    if (header.ascii2code[*p]) {
        key1 |= keyCodeRC; // At 3'
        key2 |= keyCode;   // At 5'
    }
    key1 &= maskCode;
    key2 &= maskCode;
}

void setPosReadKey(ogSequenceKeys &keys, kseq_t *seq, uint16_t pos) {
    if (keys.completed != 100) {
        getKeysFR(seq->seq.s + pos, keys.idxForward[pos], keys.idxReverse[pos], keys.header);
        keys.tupTabForward[pos] = keys.chunkRoot[keys.idxForward[pos]];
        keys.tupTabReverse[pos] = keys.chunkRoot[keys.idxReverse[pos]];
        keys.completed = 1;
        keys.nKeys = 1;
    }
}


void setFirstReadKey(ogSequenceKeys &keys, kseq_t *seq) {
    if (keys.completed != 100) {
        setPosReadKey(keys, seq, 0);
    }
}


uint16_t setLastReadKey(ogSequenceKeys &keys, kseq_t *seq) {
    if (keys.completed != 100) {
        setPosReadKey(keys, seq, seq->seq.l - keys.header.tupleSizeInChars - 2);
    }
    return seq->seq.l - keys.header.tupleSizeInChars - 2;
}



void setAllReadKeys(ogSequenceKeys &keys, kseq_t *seq) {
    if (keys.completed != 100) {
        setFirstReadKey(keys, seq);
        char *p = seq->seq.s + keys.header.tupleSizeInChars - 1;
        uint32_t idxF = keys.idxForward[0];
        uint32_t idxR = keys.idxReverse[0];
        uint16_t k;
        uint16_t kk = (seq->seq.l > MAX_READ_LEN ? MAX_READ_LEN : seq->seq.l);
        for (k=1; k < kk; k++) {
            shiftKeysFR(++p, idxF, idxR, keys.header);
            keys.idxForward[k] = idxF;
            keys.idxReverse[k] = idxR;
            keys.tupTabForward[k] = keys.chunkRoot[keys.idxForward[k]];
            keys.tupTabReverse[k] = keys.chunkRoot[keys.idxReverse[k]];
        }
        keys.completed = 100;
        keys.nKeys = kk;
    }
}

void resetReadKeys(ogSequenceKeys &keys) {
    keys.completed = 0;
    keys.nKeys = 0;
}




// Constructor to create an index file
ogIndexOriginal::ogIndexOriginal() {

    header.pChunks = NULL;
    
    int i;
    for (i=0; i < ALL_ASCII; i++) ascii2pack[i] = 0;
    ascii2pack['A'] = 0x00; 
    ascii2pack['a'] = 0x00; 
    ascii2pack['C'] = 0x40; 
    ascii2pack['c'] = 0x40; 
    ascii2pack['G'] = 0x80; 
    ascii2pack['g'] = 0x80; 
    ascii2pack['T'] = 0xC0; 
    ascii2pack['t'] = 0xC0; 
    ascii2pack['U'] = 0xC0; 
    ascii2pack['u'] = 0xC0; 
    
}


// Destructor
ogIndexOriginal::~ogIndexOriginal() {
    if (header.pChunks != NULL) { free(header.pChunks); header.pChunks = NULL; }
    if (PACKED != NULL) { free(PACKED); PACKED = NULL; }
    if (ALLPOS != NULL) { free(ALLPOS); ALLPOS = NULL; }
    if (allChromosomes != NULL) { free(allChromosomes); allChromosomes = NULL; }
    OG_INDEX_CHROMOSOME *chr = firstChromosome, *chr2;
    while (chr && chr->next) chr = chr->next; // Go to last
    while (chr && chr->prev) { chr2 = chr->prev; free(chr); chr = chr2; }
    if (chr) free(chr);
            
    fprintf(stderr, "memory free.\n");
}

void ogIndexOriginal::setSchedule(char *s) {
    strncpy(schedule,s,MAX_SCHEDULE);
}

void ogIndexOriginal::setCodeMasks() {
    keyCode = 0x00000001 << (header.tupleSizeInChars - 1);
    keyCodeRC = 0x00000001;
    maskCode = 0xFFFFFFFF >> (32 - header.tupleSizeInChars);    
}

void ogIndexOriginal::setMemoryMode(char lowMemMode) {
    lowMemoryMode = lowMemMode;
}

void ogIndexOriginal::setMaxReadLength(uint16_t len) {
    MAX_READ_LEN = len;
}

int ogIndexOriginal::generateIndex(char *filename, int keySize) {

    uint32_t genomeSize;
    uint32_t largest;
    uint32_t slen;
    uint32_t n, m;
    uint32_t idx;
    gzFile fp;
    kseq_t *seq;
    int i;

    strncpy(header.indexVersion, "oriGen Index v0.1;", 32);
    strncpy(header.fileName, filename, MAX_OGINDEX_FILENAME);
    fprintf(stderr, "Reading Genome 1/3: Estimating chromosome sizes.\n");    
    snprintf(xfile, MAX_OGINDEX_FILENAME, "%s.%d.ogs", filename, keySize);
    pfile = fopen(xfile, "w");
    fp = gzopen(filename, "r");
    seq = kseq_init(fp);
    for (largest=genomeSize=idx=0; kseq_read(seq) >= 0; genomeSize += seq->seq.l, idx++) {
        if (seq->seq.l > largest) largest = seq->seq.l;
        fprintf(stderr, "Sequence:%s, length:%ld, cummulative:%ld\n", seq->name.s, seq->seq.l, genomeSize + seq->seq.l);
        fprintf(pfile, "%d\t%s\t%ld\t%ld\n%s\n", idx+1, seq->name.s, seq->seq.l, genomeSize + seq->seq.l, seq->comment.s);
    }
    header.memoryPacked = genomeSize / 4 + 1;
    kseq_destroy(seq);
    gzclose(fp);
    fclose(pfile);
    
    header.indexSize = genomeSize;
    header.tupleLength = 3; // by default
    header.tupleSizeInChars = keySize; // old: header.tupleLength * 8;
    header.tupleArraySize = pow(2, header.tupleSizeInChars);
    header.tupleChunkSize = (genomeSize / header.tupleArraySize); // Estimated number of positions per key
    header.largestSequence = largest;
    header.indexType = 0;
    for (i=0; i < ALL_ASCII; i++) header.ascii2code[i] = 0;
    header.ascii2code['A'] = 1; 
    header.ascii2code['a'] = 1; 
    header.ascii2code['T'] = 1; 
    header.ascii2code['t'] = 1; 
    header.ascii2code['U'] = 1; 
    header.ascii2code['u'] = 1; 
    setCodeMasks();
    
    fprintf(stderr, "(genomeSize) indexSize=%u\n",header.indexSize);
    fprintf(stderr, "tupleSizeInChars=%d\n",header.tupleSizeInChars);
    fprintf(stderr, "ChunkSize=%d\n",header.tupleChunkSize);
    fprintf(stderr, "tupleArraySize=%u\n",header.tupleArraySize);
    
    header.memoryTuples = header.tupleArraySize * sizeof(OG_TUPLE_CHUNK);
    fprintf(stderr, "Allocating memory for keys 1/3=%llu, (%llu MB)\n", header.memoryTuples, header.memoryTuples >> 20);
    // calloc fills 0's
    //header.pChunks = (OG_TUPLE_CHUNK (**)) calloc(header.tupleArraySize, sizeof(OG_TUPLE_CHUNK));
    header.pChunks = (OG_TUPLE_CHUNK (*)[]) calloc(header.tupleArraySize, sizeof(OG_TUPLE_CHUNK));
    
    fprintf(stderr, "Allocating memory for packed genome 2/3=%llu, (%llu MB)\n", header.memoryPacked, header.memoryPacked >> 20);
    PACKED = (char *) malloc(header.memoryPacked);

    header.memoryPositions = header.indexSize * sizeof(uint32_t);
    fprintf(stderr, "Allocating memory genome positions 3/3=%llu, (%llu MB)\n",header.memoryPositions, header.memoryPositions >> 20);
    ALLPOS = (char *) malloc(header.memoryPositions);
    
    // Write Header File
    snprintf(xfile, MAX_OGINDEX_FILENAME, "%s.%d.ogg", filename, keySize);
    fprintf(stderr, "Saving [%s].\n",xfile);
    pfile = fopen(xfile, "wb");
    fwrite(&header, sizeof(header), 1, pfile);
    fclose(pfile);

    fprintf(stderr, "Reading Genome 2/3: Estimating structure sizes.\n");
    snprintf(xfile, MAX_OGINDEX_FILENAME, "%s.%d.ogp", filename, keySize);
    pfile = fopen(xfile, "wb");
    fp = gzopen(filename, "r");
    seq = kseq_init(fp);
    char pack = 0;
    char *packPos = PACKED;
    char a;
    char *p2;
    const char  shift[] = {    0,    2,    4,    6 };
    const char  nt[]    = {  'A',  'C',  'G',  'T' };
    for (slen=0; kseq_read(seq) >= 0; slen += seq->seq.l) {
        char *p = seq->seq.s;
        if (seq->seq.l) {
            fprintf(stderr, "Sequence:%s, length:%ld, cummulative:%ld\n", seq->name.s, seq->seq.l, slen + seq->seq.l);
            m = seq->seq.l - header.tupleSizeInChars;            
            idx = getKey(p);
            p2 = p;
            p += header.tupleSizeInChars - 1;
            for (n=0; n < m; n++) {
                // Pack
                a = ascii2pack[*p2++];
                *packPos |= (a >> shift[pack]);
                if (pack == 3) {
                    packPos++;
                    pack = 0;
                } else {
                    pack++;
                }
                
                // Re-Key
                (*header.pChunks)[idx].size++;
                shiftKey(++p, idx);
                
                // Info
                if (n % 2000000 == 0) {
                    fprintf(stderr, "."); fflush(stderr);
                    if (n > 0 && n % 100000000 == 0) fprintf(stderr, "/\n");
                }
            }
            // Finish packing last segment
            for (; n < seq->seq.l; n++) {
                a = ascii2pack[*p2++];
                *packPos |= (a >> shift[pack]);
                if (pack == 3) {
                    packPos++;
                    pack = 0;
                } else {
                    pack++;
                }                
            }
            fprintf(stderr, "|\n");
        }
    }
    fwrite(PACKED, header.memoryPacked, 1, pfile);
    kseq_destroy(seq);
    gzclose(fp);
    fclose(pfile);

    
    //Check most counts and 0 counts
    
    uint32_t k = 0;
    uint64_t max = 0;
    uint64_t v;
    uint32_t cuts[56] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 
                        10, 20, 30, 40, 50, 60, 70, 80, 90, 
                        100, 200, 300, 400, 500, 600, 700, 800, 900,
                        1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000,
                        10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000,
                        100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000,
                        1000000, (uint32_t) 0xFFFFFFFF };
    uint32_t counts[56] = { 0, 0, 0, 0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0, 0, 0, 0,
                         0, 0 };
    for (n=0; n < header.tupleArraySize; n++) {
        v = (*header.pChunks)[n].size; //(uint64_t) header.pChunks[n];
        for (k=0; k < 56; k++) if (v < cuts[k]) { counts[k]++; break; }
        if (v > max) {
            max = v;
            k = n;
            fprintf(stderr, "Max value = %llu, \tMax position = %u 0x%x\n",max,k,k);
        }
    }
    
    fprintf(stderr, "Max value = %llu, \tMax position = %u 0x%x\n",max,k,k);

    fprintf(stderr, "Distribution of tuple's positions:\n");
    for (k=0; k < 56; k++) fprintf(stderr, "< %u = %u\n", cuts[k], counts[k]);
    
    // assign offset pointers
    fprintf(stderr, "Assigning Memory Pointers.\n");
    uint32_t *currpoint = (uint32_t *) ALLPOS;
    OG_TUPLE_CHUNK *currtuple = (*header.pChunks);
    for (idx=0, n=0; n < header.tupleArraySize; n++) {
        currtuple->offset = idx;
        idx += currtuple->size;
        //currtuple->size = 0; // not restart now to save tuples
        currtuple++;
    }
    
    // Write Tuples
    snprintf(xfile, MAX_OGINDEX_FILENAME, "%s.%d.ogk", filename, keySize);
    fprintf(stderr, "Saving [%s].\n",xfile);
    pfile = fopen(xfile, "wb");
    fwrite(header.pChunks, sizeof(OG_TUPLE_CHUNK), header.tupleArraySize, pfile);
    fclose(pfile);

    
    
    char qq;
    
    if (lowMemoryMode) {
        qq='1';
        fprintf(stderr, "____ Using HALF genome memory access mode ____\n");
    } else {
        qq='-'; // Full memory use
        fprintf(stderr, "^^^^ Using FULL genome memory access mode ^^^^\n");
    }
    
    for (; qq < '3'; qq = (qq == '-' ? '9' : qq+1)) {
        // Restart sizes to save positions orderly
        fprintf(stderr, "Restarting internal sizes (%c).\n", qq);
        currtuple = (*header.pChunks);
        for (n=0; n < header.tupleArraySize; n++) {
            currtuple->size = 0; // restart to add positions 
            currtuple++;
        }


        // Recalculate
        currpoint = (uint32_t *) ALLPOS;
        uint32_t maxpoint = header.memoryPositions / 2, point;
        uint32_t *curroffset;
        fprintf(stderr, "Reading Genome 3/3 (%c): Saving positions into structures.\n", qq);
        fp = gzopen(filename, "r");
        seq = kseq_init(fp);
        char x;
        for (slen=0; kseq_read(seq) >= 0; slen += seq->seq.l) {
            char *p = seq->seq.s;
            if (seq->seq.l) {
                fprintf(stderr, "Sequence:%s, length:%ld, cummulative:%ld\n", seq->name.s, seq->seq.l, slen + seq->seq.l);
                //m = seq->seq.l - header.tupleSizeInChars;
                m = slen + seq->seq.l - header.tupleSizeInChars;
                idx = getKey(p);
                p += header.tupleSizeInChars - 1;
                for (x=0, n = slen; n < m; n++) {
                    if (n % 2000000 == 0) {
                        fprintf(stderr, "-"); fflush(stderr);
                        if (++x % 50 == 0) fprintf(stderr, "/\n");
                    }
                    currtuple = (*header.pChunks)+idx;
                    point = currtuple->offset + currtuple->size;
                    //currpoint[currtuple->offset + currtuple->size] = n;
                    //curroffset = currpoint + point;
                    //*curroffset = n;
                    if (qq == '-') {
                        currpoint[point] = n;
                    } else if (qq == '1') {
                        // First round, store all low half memory
                        if (point < maxpoint) {
                            currpoint[point] = n;
                        }
                    } else { // qq must be '2'
                        // Second round, store all high half memory
                        if (point >= maxpoint) {
                            currpoint[point] = n;
                        }
                    }
                    currtuple->size++;
                    shiftKey(++p, idx);
                }
                fprintf(stderr, "|\n");
            }
        }
        kseq_destroy(seq);
        gzclose(fp);
    }
    
    // Write OffSets
    snprintf(xfile, MAX_OGINDEX_FILENAME, "%s.%d.ogo", filename, keySize);
    fprintf(stderr, "Saving [%s].\n",xfile);
    pfile = fopen(xfile, "wb");
    fwrite(ALLPOS, 1, header.memoryPositions, pfile);
    fclose(pfile);
    
    return 1;
}


uint32_t ogIndexOriginal::getKey(char *p) {
    uint32_t    key = 0;
    uint32_t    code = keyCode;
    uint8_t     i;
    for (i=0; i < header.tupleSizeInChars; i++, p++) {
        if (header.ascii2code[*p]) key |= code;
        code >>= 1;
    }
    
    return key;
}


void ogIndexOriginal::shiftKey(char * &p, uint32_t &key) {
    key <<= 1;
    //if (header.ascii2code[*(p+header.tupleSizeInChars_1)]) key |= code1;
    if (header.ascii2code[*p]) key |= keyCodeRC; // new code is always at 3'
    key &= maskCode;
}


int ogIndexOriginal::loadIndex(char *filename, int keySize) {

    snprintf(xfile, MAX_OGINDEX_FILENAME, "%s.%d.ogs", filename, keySize); // .ogs = Genome & Chromosome SIZES
    fprintf(stderr, "Reading Genome Sizes...\n[%s]\n", xfile);
    pfile = fopen(xfile, "r");
    OG_INDEX_CHROMOSOME *chr = (OG_INDEX_CHROMOSOME *) malloc(sizeof(OG_INDEX_CHROMOSOME));
    OG_INDEX_CHROMOSOME *lastChr = NULL;
    chr->next = chr->prev = NULL;
    firstChromosome = NULL;
    int nchr = 0;
    while (fscanf(pfile, "%d\t%s\t%u\t%u\n", &chr->number, chr->seqName, &chr->size, &chr->cummulative) > 0) {
        chr->start = chr->cummulative - chr->size;
        if (lastChr != NULL) chr->prev = lastChr;
        if (firstChromosome == NULL) firstChromosome = chr;
        else lastChr->next = chr;
        lastChr = chr;
        chr = (OG_INDEX_CHROMOSOME *) malloc(sizeof(OG_INDEX_CHROMOSOME));
        while (fgetc(pfile) != '\n');
        nchr++;
        //fprintf(stderr, "%d\t%s\t%u\t%u\n", lastChr->number, lastChr->seqName, lastChr->size, lastChr->cummulative);
    }
    fprintf(stderr, "N=%d, Chromosomes=%d, Size=%u\n", nchr, lastChr->number, lastChr->cummulative);
    allChromosomes = (OG_INDEX_CHROMOSOME **) malloc(sizeof(OG_INDEX_CHROMOSOME *) * nchr);
    free(chr); // last one not used.    
    for(nchr=0, chr = firstChromosome; chr->next != NULL; nchr++, chr=chr->next) allChromosomes[nchr] = chr;
    nChromosomes = nchr;
    //fseek(fp, 0L, SEEK_END);
    //uint32_t fileSize = ftell(fp);
    //fseek(f, 0, SEEK_SET);
    fclose(pfile);

    snprintf(xfile, MAX_OGINDEX_FILENAME, "%s.%d.ogg", filename, keySize); // .ogg = Genome 
    fprintf(stderr, "Done.\nReading Genome Structure...\n[%s]\n", xfile);
    pfile = fopen(xfile, "r");
    fread(&header, sizeof(header), 1, pfile);
    fclose(pfile);
    setCodeMasks();

    fprintf(stderr, "Done.\nAllocating memory for keys 1/3=%llu (%llu MB)\n", header.memoryTuples, header.memoryTuples >> 20);
    header.pChunks = (OG_TUPLE_CHUNK (*)[]) calloc(header.tupleArraySize, sizeof(OG_TUPLE_CHUNK));

    fprintf(stderr, "Allocating memory for packed genome 2/3=%llu (%llu MB)\n", header.memoryPacked, header.memoryPacked >> 20);
    PACKED = (char *) malloc(header.memoryPacked);
 
    header.memoryPositions = header.indexSize * sizeof(uint32_t);
    fprintf(stderr, "Allocating memory for genome positions 3/3=%llu (%llu MB)\n",header.memoryPositions, header.memoryPositions >> 20);
    ALLPOS = (char *) malloc(header.memoryPositions);

    snprintf(xfile, MAX_OGINDEX_FILENAME, "%s.%d.ogk", filename, keySize); // .ogk = KEYS 
    fprintf(stderr, "Reading Tuple Keys...\n[%s]\n", xfile);
    pfile = fopen(xfile, "r");
    fread(header.pChunks, header.tupleArraySize, sizeof(OG_TUPLE_CHUNK), pfile);
    fclose(pfile);
    
    snprintf(xfile, MAX_OGINDEX_FILENAME, "%s.%d.ogo", filename, keySize); // .ogo = Offsets
    fprintf(stderr, "Done.\nReading Offsets...\n[%s]\n", xfile);
    pfile = fopen(xfile, "r");
    fread(ALLPOS, 1, header.memoryPositions, pfile);
    fclose(pfile);
    
    snprintf(xfile, MAX_OGINDEX_FILENAME, "%s.%d.ogp", filename, keySize); // .ogp = Packed
    fprintf(stderr, "Done.\nReading Packed genome...\n[%s]\n", xfile);
    pfile = fopen(xfile, "r");
    fread(PACKED, 1, header.memoryPacked, pfile);
    fclose(pfile);
    
    return 1;
}


void ogIndexOriginal::allocateKeys(ogSequenceKeys &readKeys, ogIndexOriginalHeader &header) {
    readKeys.idxForward = (uint32_t *) malloc(MAX_READ_LEN * sizeof(uint32_t));
    readKeys.idxReverse = (uint32_t *) malloc(MAX_READ_LEN * sizeof(uint32_t));
    readKeys.tupTabForward = (OG_TUPLE_CHUNK *) malloc(MAX_READ_LEN * sizeof(OG_TUPLE_CHUNK *));
    readKeys.tupTabReverse = (OG_TUPLE_CHUNK *) malloc(MAX_READ_LEN * sizeof(OG_TUPLE_CHUNK *));
    readKeys.chunkRoot = (*header.pChunks);
    readKeys.header = header;
}

void ogIndexOriginal::deallocateKeys(ogSequenceKeys &readKeys) {
    if (readKeys.idxForward != NULL) free(readKeys.idxForward);
    if (readKeys.idxReverse != NULL) free(readKeys.idxReverse);
    if (readKeys.tupTabForward != NULL) free(readKeys.tupTabForward);
    if (readKeys.tupTabReverse != NULL) free(readKeys.tupTabReverse);
    readKeys.idxForward = NULL;
    readKeys.idxReverse = NULL;
    readKeys.tupTabForward = NULL;
    readKeys.tupTabReverse = NULL;
}


void printKeyPositions(ogIndexOriginal &ogIdx, uint32_t offset, uint32_t size) {
    uint32_t i;
    OG_INDEX_CHROMOSOME *chr;
    uint32_t *p = (uint32_t *) ALLPOS + offset;
    fprintf(stderr, "//// Positions Root *: %p, offset: %u, size: %u. Starting *: %p.\n", ALLPOS, offset, size, p);
    for (i=0; i < size; i++, p++) {
        chr = ogIdx.findChromosome(*p);
        fprintf(stderr, "%4u:%10u [%-5s %10u]%c", i+1, *p, chr->seqName, *p - chr->start, (i % 4 == 3 ? '\n' : ' '));
    }
    fprintf(stderr, "\n");
}


// Algorithm - QUICKEST! Exact Contiguous in extremes
uint16_t ExactInExtremes(ogIndexOriginal &ogIdx, kseq_t *seq, ogSequenceKeys &keys, ogMatches &matches) {
    uint16_t pos = setLastReadKey(keys, seq);
    char RC = 0;
    setFirstReadKey(keys, seq);
    //printKeyPositions(ogIdx, keys.tupTabForward[0].offset, keys.tupTabForward[0].size);
    //printKeyPositions(ogIdx, keys.tupTabForward[pos].offset, keys.tupTabForward[pos].size);
    ogIdx.intersectPositions(simPositions, simPosK, keys.tupTabForward[0], keys.tupTabForward[pos], pos+5);
    if (simPosK == 0) {
        ogIdx.intersectPositions(simPositions, simPosK, keys.tupTabReverse[0], keys.tupTabReverse[pos], pos+5);
        RC = 1;
    }
    if (simPosK > 0) {
        uint16_t i;
        for (i=0; i < simPosK; i++) {
            pushMatch(matches, simPositions[i], RC, 100);
        }
        setMatchStatus(matches, 1);
    }
    return simPosK;
}


// Algorithm - Divide in 5 pieces, take 2nd and 4th
uint16_t ExactInApartSides(ogIndexOriginal &ogIdx, kseq_t *seq, ogSequenceKeys &keys, ogMatches &matches) {
    uint32_t l = seq->seq.l - 24;//header.tupleSizeInChars;     // Parameter 0: (this should be taken from tupleSizeInChars)
    uint32_t p1 = l/5;                                          // Parameter 1: 5=How many parts should the read be split into?
    uint32_t p2 = l - p1;
    uint32_t d = p2 - p1 + 5;                                   // Parameter 2: 5=How many nt of theoretical difference should be considered as small indels?
    char RC = 0;
    //ogIdx.setPosReadKey(keys, seq, p1);
    //ogIdx.setPosReadKey(keys, seq, p2);
    setAllReadKeys(keys, seq);
    ogIdx.intersectPositions(simPositions, simPosK, keys.tupTabForward[p1], keys.tupTabForward[p2], d);
    if (simPosK == 0) {
        ogIdx.intersectPositions(simPositions, simPosK, keys.tupTabReverse[p1], keys.tupTabReverse[p2], d);
        RC = 1;
    }
    if (simPosK > 2) {
        p2 = (p2 - p1) >> 1;
        d = p2 - p1 + 3;
        //ogIdx.setPosReadKey(keys, seq, p2);
        if (RC) 
            ogIdx.intersectPositions(simPositions, simPosK, keys.tupTabReverse[p2], keys.tupTabReverse[p2], d);
        else
            ogIdx.intersectPositions(simPositions, simPosK, keys.tupTabForward[p2], keys.tupTabForward[p2], d);
    }
    if (simPosK > 0) {
        uint16_t i;
        for (i=0; i < simPosK; i++) {
            pushMatch(matches, simPositions[i], RC, 100);
        }
        setMatchStatus(matches, 1);
    }
    return simPosK;
}


// Algorithm - Let 5p Keys
uint16_t Two5pKeys(ogIndexOriginal &ogIdx, kseq_t *seq, ogSequenceKeys &keys, ogMatches &matches) {
    uint32_t keysize = 24;                                  // Parameter: header.tupleSizeInChars
    uint32_t maxl = seq->seq.l - keysize;
    uint32_t pos2 = keysize + 5;                            // Parameter: offset to enlarge key
    if (pos2 > maxl) pos2 = maxl;
    char RC = 0;
    setPosReadKey(keys, seq, 0);
    setPosReadKey(keys, seq, pos2);
    ogIdx.intersectPositions(simPositions, simPosK, keys.tupTabForward[0], keys.tupTabForward[pos2], pos2+5);
    if (simPosK == 0) {
        ogIdx.intersectPositions(simPositions, simPosK, keys.tupTabReverse[0], keys.tupTabReverse[pos2], pos2+5);
        RC = 1;
    }
    if (simPosK > 2) {
        pos2 = pos2 >> 1;
        setPosReadKey(keys, seq, pos2);
        if (RC) 
            ogIdx.intersectPositions(simPositions, simPosK, keys.tupTabReverse[pos2], keys.tupTabReverse[pos2], pos2+5);
        else
            ogIdx.intersectPositions(simPositions, simPosK, keys.tupTabForward[pos2], keys.tupTabForward[pos2], pos2+5);
    }
    if (simPosK > 0) {
        uint16_t i;
        for (i=0; i < simPosK; i++) {
            pushMatch(matches, simPositions[i], RC, 100);
        }
        setMatchStatus(matches, 1);
    }
    return simPosK;
}


// Algorithm 3 - MIXTURE at both ends
uint16_t MixtureInBothSides(ogIndexOriginal &ogIdx, kseq_t *seq, ogSequenceKeys &keys, ogMatches &matches) {
    uint32_t l = seq->seq.l - 24;//header.tupleSizeInChars;     // Parameter 0: (this should be taken from tupleSizeInChars)
    uint32_t p1 = l/5;                                          // Parameter 1: 5=How many parts should the read be split into?
    uint32_t p2 = l - p1;
    uint32_t d = p2 - p1 + 5;                                   // Parameter 2: 5=How many nt of theoretical difference should be considered as small indels?
    char RC = 0;
    setAllReadKeys(keys, seq);
    ogIdx.intersectPositions(simPositions, simPosK, keys.tupTabForward[p1], keys.tupTabForward[p2], d);
    if (simPosK == 0) {
        ogIdx.intersectPositions(simPositions, simPosK, keys.tupTabReverse[p1], keys.tupTabReverse[p2], d);
        RC = 1;
    }
    if (simPosK > 2) {
        p2 = (p2 - p1) >> 1;
        d = p2 - p1 + 3;
        //ogIdx.setPosReadKey(keys, seq, p2);
        if (RC) 
            ogIdx.intersectPositions(simPositions, simPosK, keys.tupTabReverse[p2], keys.tupTabReverse[p2], d);
        else
            ogIdx.intersectPositions(simPositions, simPosK, keys.tupTabForward[p2], keys.tupTabForward[p2], d);
    }
    if (simPosK > 0) {
        uint16_t i;
        for (i=0; i < simPosK; i++) {
            pushMatch(matches, simPositions[i], RC, 100);
        }
        setMatchStatus(matches, 1);
    }
    return simPosK;
}



// Algorithm x - default
uint16_t DefaultMatching(ogIndexOriginal &ogIdx, kseq_t *seq, ogSequenceKeys &keys, ogMatches &matches) {
    setAllReadKeys(keys, seq);
   
    uint32_t minGenomePos = 0;
    uint32_t maxGenomePos = 4294967295;
    uint32_t deltaGenomePos;
    uint32_t maxForwardPos = 0;
    uint32_t maxForwardCounts = 0;
    uint32_t maxReversePos = 0;
    uint32_t maxReverseCounts = 0;
    uint32_t i, j, k, bestSP;
    uint32_t inc;
    OG_TUPLE_CHUNK &t = keys.tupTabForward[0];
    uint32_t *initalPOS = (uint32_t *) ALLPOS;
    uint32_t *posK;
    uint32_t gPos;
    uint16_t nPos = 0;
    float bin;
    uint16_t bestPos;
    OG_TUPLE_CHUNK *tuptab;
    
    for (k=0; k < 2; k++) {
        minGenomePos = 0;
        maxGenomePos = 4294967295;

        while (true) {
            // Set counts to 0
            deltaGenomePos=(maxGenomePos-minGenomePos); // >= MAX_SIMILAR_POSITIONS
            memset(simPositions, 0, MAX_SIMILAR_POSITIONS);
            bestSP = 0;
            bin = ((float) deltaGenomePos / (float) MAX_SIMILAR_POSITIONS);
            //fprintf(stderr, "Mode %u, %u ~ %u : %u : %.2f ==> ", k, minGenomePos, maxGenomePos, deltaGenomePos, bin);

            // Create the histogram of positions "matched" bins of size deltaGenomePos / MAX_SIMILAR_POSITIONS
            inc = keys.nKeys / 10;
            if (inc < 1) inc = 1;
            nPos = 0;
            for (i=0; i < keys.nKeys; i += inc) {
                // forward
                //fprintf(stderr, "Key=%d, n=%d\n",i,keys.nKeys);
                t = (k == 0 ? keys.tupTabForward[i] : keys.tupTabReverse[i]);
                posK = initalPOS + t.offset;
                for (j=0; j< t.size; j++, posK++) {
                    //fprintf(stderr, "j=%d, pos=%u\n",j,*posK);
                    gPos = *posK - j; // -j to put relative to the start of the sequence and make counts properly
                    if (gPos >= minGenomePos && gPos <= maxGenomePos) {
                        simPosK = (gPos - minGenomePos) * MAX_SIMILAR_POSITIONS / deltaGenomePos;
                        if (simPositions[simPosK] == 0) nPos++;
                        if (++simPositions[simPosK] > simPositions[bestSP]) bestSP = simPosK;
                    }
                }
            }
            //fprintf(stderr, "Best Pos: %u=%u\n", bestSP, simPositions[bestSP]);
            if (simPositions[bestSP] < 2 || bin <= 1) break;

            minGenomePos += (uint32_t) (bin * bestSP);
            maxGenomePos  = minGenomePos + bin + 1;
            
        }

        if (k == 0) {
            maxForwardPos   = minGenomePos + bin * bestSP;
            maxForwardCounts = simPositions[bestSP];
        } else {
            maxReversePos = minGenomePos + bin * bestSP;            
            maxReverseCounts = simPositions[bestSP];
        }
        
    }
    
    if (maxForwardCounts > 1 || maxReverseCounts > 1) {
        if (maxForwardCounts >= maxReverseCounts) {
            pushMatch(matches, maxForwardPos, 0, 100); // forward
            setMatchStatus(matches, 1);            
        } else {
            pushMatch(matches, maxReversePos, 1, 100); // Reverse
            setMatchStatus(matches, 1);                        
        }
    } else {
        setMatchStatus(matches, 0);
    }
    
    return nPos;
}

void ogIndexOriginal::setMaxReads(uint64_t mxreads) {
    MAX_READS_MAP = mxreads;
}

int ogIndexOriginal::mapSequences(char *filename) {
    const int MAX_IDX = 100;
    uint16_t scheduLen = strlen(schedule);
    uint64_t *nChecked = (uint64_t *) calloc(scheduLen, sizeof(uint64_t));
    uint64_t *nMatched = (uint64_t *) calloc(scheduLen, sizeof(uint64_t));
    uint32_t n, m;
    uint32_t idx[MAX_IDX];
    uint32_t leidx;
    OG_TUPLE_CHUNK *tuptab[MAX_IDX];
    gzFile fp;
    kseq_t *seq;
    uint64_t i;
    char *p;
    char *p2;
    uint16_t k, j, l;
    uint32_t *gPos;
    uint32_t *initalPOS = (uint32_t *) ALLPOS;
    OG_INDEX_CHROMOSOME *iChr;
    char revcomp[50000];
    
    fprintf(stderr, "Reading...\n");
    fp = gzopen(filename, "r");
    seq = kseq_init(fp);
    uint32_t d = 10;
    uint64_t found=0, notfound = 0, foundcomp = 0, notchecked = 0;
    clock_t start_t, end_t;
    double total_t;
    start_t = clock();
    uint64_t cumSimPos=0;
    
    ogSequenceKeys      readKeys;
    allocateKeys(readKeys, header);
    
    ogMatches   matches;

    // Set Schedule
    uint16_t (*funcSchedule[MAX_SCHEDULE])(ogIndexOriginal &ogIdx, kseq_t *seq, ogSequenceKeys &keys, ogMatches &matches);
    for (j=0; j < scheduLen; j++) {
        switch (schedule[j]) {
            case 'E' : funcSchedule[j] = &ExactInExtremes; break;
            case 'D' : funcSchedule[j] = &DefaultMatching; break;
            case 'A' : funcSchedule[j] = &ExactInApartSides; break;
            case '5' : funcSchedule[j] = &Two5pKeys; break;
        }
    }
    
    uint64_t nMatches = 0;
    for (i=1; kseq_read(seq) >= 0 && i < MAX_READS_MAP; i++) {
        //if (i < 10) fprintf(stderr, "%llu %s %s %s\n", i, seq->name.s, seq->qual.s, seq->seq.s);
        if (seq->seq.l > header.tupleSizeInChars && seq->qual.l > header.tupleSizeInChars && seq->name.l > 0) {
            resetReadKeys(readKeys);
            initMatches(matches);
            for (j=0; j < scheduLen; j++) {
                nChecked[j]++;
                nMatches += funcSchedule[j](*this, seq, readKeys, matches);
                if (matches.status) {
                    nMatched[j]++;
                    if (matches.isRevComp[0]) foundcomp++; else found++;
                    break;
                }
            }
            if (i <= 10) {
                iChr = findChromosome(matches.position[matches.maxMatch]);
                fprintf(stderr, "%llu %s %s M:%u P:%u Chr:%d %u\n", i, seq->name.s, seq->seq.s, matches.nMatches, matches.position[matches.maxMatch], iChr->number, matches.position[matches.maxMatch] - iChr->start); // seq->qual.s, 
            }
            if (matches.nMatches == 0) {
                notfound++;
            } else {
                // for now, take the maxMatch
                // matches.maxMatch index point to maximum match
                if (0) {
                    p = seq->seq.s;
                    idx[0] = getKey(p);
                    tuptab[0] = (*header.pChunks) + idx[0];
                    p += header.tupleSizeInChars - 1;
                    leidx = idx[0];
                    for (k=1; k < MAX_IDX && k < seq->seq.l - header.tupleSizeInChars; k++) {
                        shiftKey(++p, leidx);
                        idx[k] = leidx;
                        tuptab[k] = (*header.pChunks) + idx[k];
                    }
                    if (k > 5) intersectPositions(simPositions, simPosK, *tuptab[5], *tuptab[k-1], k - 5 + 5);
                    else intersectPositions(simPositions, simPosK, *tuptab[0], *tuptab[1], 5);
                    if (simPosK > 0) {
                        found++;
                        cumSimPos += simPosK;
                    } else {
                        // no es necesario el complemento dado que generan el mismo codigo
                        char *x = seq->seq.s, *y = revcomp+seq->seq.l;
                        *y-- = 0;
                        for (m=seq->seq.l; m > 0; m--) *y-- = *x++;
                        //for (m=1; m <= seq->seq.l; m++) revcomp[seq->seq.l - m] = seq->seq.s[m-1];
                        //revcomp[seq->seq.l] = 0;
                        p = revcomp;
                        idx[0] = getKey(p);
                        tuptab[0] = (*header.pChunks) + idx[0];
                        p += header.tupleSizeInChars - 1;
                        leidx = idx[0];
                        for (k=1; k < MAX_IDX && k < seq->seq.l - header.tupleSizeInChars; k++) {
                            shiftKey(++p, leidx);
                            idx[k] = leidx;
                            tuptab[k] = (*header.pChunks) + idx[k];
                        }
                        if (k > 5) intersectPositions(simPositions, simPosK, *tuptab[5], *tuptab[k-1], k - 5 + 5);
                        else intersectPositions(simPositions, simPosK, *tuptab[0], *tuptab[1], 5);
                        if (simPosK > 0) {
                            foundcomp++; 
                            cumSimPos += simPosK;
                        } else notfound++;
                    }
                }
            }
        } else notchecked++;
        
        if (i % 1000 == 0) { //  // i % 1000 == 0 // (i & 0x000003FF) == 0
            end_t = clock();
            total_t = (double)(end_t - start_t) / (double) CLOCKS_PER_SEC;
            fprintf(stderr, "%llu: %llu + %llu Ok + compl (%.1f), %llu not found, %llu not checked. Time=%.1f, Speed=%.1f/s Matches/Read=%.1f.\n", i, found, foundcomp, (double) cumSimPos / (double) (found+foundcomp), notfound, notchecked, total_t, (double) i / total_t, (double) nMatches / (double) (found+foundcomp));
        }
    }
    end_t = clock();
    total_t = (double)(end_t - start_t) / (double) CLOCKS_PER_SEC;
    fprintf(stderr, "%llu: %llu + %llu Ok + compl (%.1f), %llu not found, %llu not checked. Time=%.1f, Speed=%.1f/s Matches/Read=%.1f.\n", i, found, foundcomp, (double) cumSimPos / (double) (found+foundcomp), notfound, notchecked, total_t, (double) i / total_t, (double) nMatches / (double) (found+foundcomp));
    fprintf(stderr, "\n %llu sequences keyed.\n", i);
    kseq_destroy(seq);
    gzclose(fp);
    
    for (j=0; j < scheduLen; j++) {
        fprintf(stderr, ">>>> Schedule option '%c' : %llu calls, %llu matches, then %llu\n",schedule[j],nChecked[j],nMatched[j],nChecked[j]-nMatched[j]);
    }
    free(nChecked);
    free(nMatched);
    deallocateKeys(readKeys);
    
    return 1;
}

void ogIndexOriginal::intersectPositions(uint32_t *simPositions, uint16_t &simPosK, OG_TUPLE_CHUNK &a, OG_TUPLE_CHUNK &b, uint32_t diff) {
    uint32_t ai, bi, j, k;
    uint32_t *pa = (uint32_t *) ALLPOS;
    uint32_t *pb = pa;
    uint32_t ga, gb;
    
    if (&b != &a) {
        // starting mode where a & b are DIFFERENT and be the first two hits.
        pa += a.offset;
        pb += b.offset;
        ga = *pa;
        gb = *pb;

        simPosK = 0;
        for (ai=bi=0; ai < a.size && bi < b.size && simPosK < MAX_SIMILAR_POSITIONS; ) {
            if (ga <= gb) {
                if (gb - ga <= diff) {
                    // THIS IS A HIT
                    simPositions[simPosK++] = ga;
                    gb = *++pb;
                    bi++;                
                }
                ga = *++pa;
                ai++;
            } else {
                if (ga - gb <= diff) {
                    // THIS IS A HIT
                    simPositions[simPosK++] = gb;
                    ga = *++pa;
                    ai++;
                }            
                gb = *++pb;
                bi++;
            }
        }
    } else {
        // Following mode where a records are "intersected" with previous simPos
        j = 0; // current simPositions value
        k  = 0; // new position to be replaced/saved
        pa += a.offset;
        ga = *pa;
        gb =  simPositions[j];
        for (ai=0; ai < a.size && j < simPosK; ) {
            if (ga <= gb) {
                if (gb - ga <= diff) {
                    // THIS IS A HIT, keep simPosition
                    simPositions[k++] = simPositions[j++];
                    gb = simPositions[j];
                }
                ga = *++pa;
                ai++;
            } else {
                if (ga - gb <= diff) {
                    // THIS IS A HIT, keep simPosition
                    simPositions[k++] = simPositions[j++];
                    gb = simPositions[j];
                    ga = *++pa;
                    ai++;
                } else {
                    j++;
                }
            }
        }
        simPosK = k;
    }
    
}

OG_INDEX_CHROMOSOME  *ogIndexOriginal::findChromosomeSequential(uint32_t &p) {
    OG_INDEX_CHROMOSOME *chrI = firstChromosome;
    while (chrI->cummulative < p && chrI->next != NULL) {
        chrI = chrI->next;
    }
    return chrI;
}


OG_INDEX_CHROMOSOME  *ogIndexOriginal::findChromosome(uint32_t &p) {
    int min = 0;
    int max = nChromosomes-1;
    int mid;
    OG_INDEX_CHROMOSOME *chrI;
    while (min < max) {
        //fprintf(stderr, "P:%u,min:%d,max:%d\n",p,min,max);
        mid = (max + min) >> 1;
        chrI = allChromosomes[mid];
        if (p > chrI->cummulative) { min = mid+1; }
        else { max = mid - (p <= (chrI->cummulative - chrI->size) ? 1 : 0); }
    }
    return allChromosomes[min];
}


void ogIndexOriginal::unPackGenomePosition(char *p, uint32_t pos, uint32_t size) {
    // pos is absolute position
    
    uint32_t    bytePos = pos / 4;
    char        pack = pos % 4;
    char        *packed = PACKED + bytePos;
    const unsigned char  masks[] = { 0xC0, 0x30, 0x0C, 0x03 };
    const unsigned char  shift[] = {    6,    4,    2,    0 };
    const unsigned char  nt[]    = {  'A',  'C',  'G',  'T' };
    
    for (; size; size--) {
        *p++ = nt[ (*packed & masks[pack]) >> shift[pack] ];
        if (pack == 3) {
            packed++;
            pack = 0;
        } else {
            pack++;
        }
    }
    *p = 0;
    
}



