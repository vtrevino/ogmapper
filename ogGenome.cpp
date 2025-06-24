/* 
 * File:   ogGenome.cpp
 * Author: victortrevino
 * 
 * Created on May 4, 2022, 7:03 PM
 */

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <cstdlib>
#include <zlib.h>
#include "ogDefinitions.hpp"
#include "ogGenome.hpp"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)


const char  shift[]   = {    0,    2,    4,    6 };
const char  unshift[] = {    6,    4,    2,    0 };
const char  ACGT[] = { 'A', 'C', 'G', 'T', 'T', 'G', 'C', 'A', 'a', 'c', 'g', 't', 't', 'g', 'c', 'a' };
const char  ERRORPACKED[] = { 0, 1, 1, 1, 1, 2, 2, 2, 1, 2, 2, 2, 1, 2, 2, 2, 1, 2, 2, 2, 2, 3, 3, 3, 2, 3, 3, 3, 2, 3, 3, 3, 1, 2, 2, 2, 2, 3, 3, 3, 2, 3, 3, 3, 2, 3, 3, 3, 1, 2, 2, 2, 2, 3, 3, 3, 2, 3, 3, 3, 2, 3, 3, 3, 1, 2, 2, 2, 2, 3, 3, 3, 2, 3, 3, 3, 2, 3, 3, 3, 2, 3, 3, 3, 3, 4, 4, 4, 3, 4, 4, 4, 3, 4, 4, 4, 2, 3, 3, 3, 3, 4, 4, 4, 3, 4, 4, 4, 3, 4, 4, 4, 2, 3, 3, 3, 3, 4, 4, 4, 3, 4, 4, 4, 3, 4, 4, 4, 1, 2, 2, 2, 2, 3, 3, 3, 2, 3, 3, 3, 2, 3, 3, 3, 2, 3, 3, 3, 3, 4, 4, 4, 3, 4, 4, 4, 3, 4, 4, 4, 2, 3, 3, 3, 3, 4, 4, 4, 3, 4, 4, 4, 3, 4, 4, 4, 2, 3, 3, 3, 3, 4, 4, 4, 3, 4, 4, 4, 3, 4, 4, 4, 1, 2, 2, 2, 2, 3, 3, 3, 2, 3, 3, 3, 2, 3, 3, 3, 2, 3, 3, 3, 3, 4, 4, 4, 3, 4, 4, 4, 3, 4, 4, 4, 2, 3, 3, 3, 3, 4, 4, 4, 3, 4, 4, 4, 3, 4, 4, 4, 2, 3, 3, 3, 3, 4, 4, 4, 3, 4, 4, 4, 3, 4, 4, 4 };
const unsigned char  LEFTOFFSETMASK[] = { 0xFF, 0x3F, 0x0F, 0x03 };
const unsigned char  RIGHTOFFSETMASK[] = { 0x00, 0xC0, 0xF0, 0xFC };

ogGenome::ogGenome(char *sourceFileName, uint16_t kSize) {
    packed = NULL;
    genomeSize = 0;
    keySize = 0;
    fastaGZ = NULL;
    packed = NULL;
    allChromosomes = NULL;

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
    for (i=0; i < ALL_ASCII; i++) ascii2pack1[i] = ascii2pack[i] >> 2;
    for (i=0; i < ALL_ASCII; i++) ascii2pack2[i] = ascii2pack[i] >> 4;
    for (i=0; i < ALL_ASCII; i++) ascii2pack3[i] = ascii2pack[i] >> 6;
    
    setSourceFile(sourceFileName);
    setKeySize(kSize);
    validPositions = 0;
    validPositionsOverLen = 0;
    nTotalPositionsWithN = 0;
}

ogGenome::~ogGenome() {
    //fprintf(stderr, "<Deallocating ogGenome:"); fflush(stderr);
    //fprintf(stderr, "<~ogGnm:"); fflush(stderr);
    if (packed != NULL) free(packed);
    if (allChromosomes != NULL) {
        if (allChromosomes[0].pPosWithN != NULL) free(allChromosomes[0].pPosWithN);
        free(allChromosomes);
    }
    //fprintf(stderr, ":ogGnm>");
}


void ogGenome::setSourceFile(char *filename) {
    closeSourceFile();
    sourceFileName = filename;
}

void ogGenome::setKeySize(uint16_t kSize) {
    closeSourceFile();
    keySize = kSize;
}

void ogGenome::estimateChromosomeSizes(char verbose) {
    openSourceFile();
    kseq_t     *seq = (kseq_t *) seqPointer;
    if (verbose) {
        fprintf(stderr, "Seq\t      Size\t Cummulative\tName\n");
    }
    for (genomeSize=nChromosomes=0; kseq_read(seq) >= 0; genomeSize += seq->seq.l, nChromosomes++) {
        if (verbose && (nChromosomes % 1000 == 0 || nChromosomes < 50)) fprintf(stderr, "%3d\t%10ld\t%12ld\t%s\n", nChromosomes+1, seq->seq.l, genomeSize + seq->seq.l, seq->name.s);
    }
    memoryPacked = genomeSize / 4 + 1;
    closeSourceFile();
    allocate();
    
}

void ogGenome::openSourceFile() {
    if (fastaGZ == NULL) {
        kseq_t     *seq;
        fastaGZ = gzopen(sourceFileName, "r");
        seq = kseq_init(fastaGZ);
        seqPointer = seq;
    }
}


char *ogGenome::readSourceFile() {
    if (fastaGZ != NULL) {
        kseq_t     *seq = (kseq_t *) seqPointer;
        kseq_read(seq);
        return seq->seq.s;
    }
    return NULL;
}


void ogGenome::closeSourceFile() {
    if (fastaGZ != NULL) {
        kseq_destroy((kseq_t *) seqPointer);
        gzclose(fastaGZ);
    }
    fastaGZ = NULL;
    seqPointer = NULL;
}



void ogGenome::packGenome(uint32_t nNtInform, ogKeyEncoding *pEncoding, ogGuider *pGuider) {
    allocate();
    openSourceFile();
    kseq_t     *seq = (kseq_t *) seqPointer;
    uint32_t    m;
    uint32_t    n;
    uint32_t    mK;
    unsigned char        a;
    char                *p;
    unsigned char       *pack_i = packed;
    unsigned char        pack_s = 0;
    uint32_t             withN;
    char                 informing = 0;
    uint32_t             chrMod = (uint32_t) pow(10,fmax(1,trunc(log10(nChromosomes))-1));
    ogChromosome *chrI;
    if (nNtInform > 0) fprintf(stderr, "Seq\t      Size\tName\n");
    nTotalPositionsWithN = 0;
    for (chr_i=0; chr_i < nChromosomes && kseq_read(seq) >= 0; chr_i++) {
        m = seq->seq.l;
        mK = m - pEncoding->getSizeInChars();
        p = seq->seq.s;
        chrI = allChromosomes + chr_i;
        // Info
        chrI->initialPack = pack_i - packed;
        chrI->initialBytePos = pack_s;
        chrI->start = (chr_i == 0 ? 0 : allChromosomes[chr_i-1].cummulative);
        chrI->index = chr_i;
        chrI->number = chr_i + 1;
        chrI->size  = seq->seq.l;
        chrI->cummulative = chrI->start + chrI->size;
        pGuider->setSequence(seq->seq.s, seq->seq.l);
        pGuider->fixSequenceCase();
        chrI->validPositions = pGuider->countGuides();
        chrI->validPositionsOverLen = pGuider->countOverLen;
        chrI->nPositionsWithN = 0;
        chrI->pPosWithN = NULL;
        withN = 0;
        if (seq->comment.s != NULL) strncpy(chrI->comment, seq->comment.s, MAX_CHROMOSOME_COMMENT); else chrI->comment[0] = 0;
        strncpy(chrI->name, seq->name.s, MAX_CHROMOSOME_NAME); 
        informing = 0;
        if (nNtInform > 0 && (chr_i % chrMod == 0  || chr_i < 50 || chr_i == nChromosomes-1)) { 
            fprintf(stderr, "%3d\t%10u\t%s ", chrI->number, chrI->size, chrI->name); informing = 1; 
        }
        for (n=0; n < m; n++, p++) {
            //a = ascii2pack[*p] >> pack_s;
            //if (n < 20) fprintf(stderr, "%d, *pack_i=%hhX=%hhu, *p=%c=%hhX,pack_s=%hhX, a2p=%hhX, |=%hhX ----- ",n, *pack_i, *pack_i,*p,*p,pack_s,ascii2pack[*p],a);
            // Pack
            *pack_i |= (ascii2pack[*p] >> shift[pack_s]);
            //if (n < 20) fprintf(stderr, "%d, *pack_i=%hhX=%hhu, *p=%c=%hhX,pack_s=%hhX\n",n, *pack_i,*pack_i,*p,*p,pack_s);
            if (pack_s == 3) {
                pack_i++;
                pack_s = 0;
            } else {
                pack_s++;
            }
            // Count N positions, if any
            if (ascii2pack[*p] == 0 && *p != 'A' && *p != 'a') {
                withN++;
            } else {
                if (withN > 0) {
                    chrI->nPositionsWithN += (withN >> 8) + 1;
                    withN = 0;
                }
            }
            // Info
            if (informing && (nNtInform > 0 && n % nNtInform == 0)) {
                fprintf(stderr, "."); fflush(stderr);
                if (n > 0 && n % (nNtInform * 50) == 0) fprintf(stderr, "/\n");
            }
        }
        if (withN > 0) {
            chrI->nPositionsWithN += (withN >> 8) + 1;
            withN = 0;
        }
        validPositions += chrI->validPositions;
        validPositionsOverLen += chrI->validPositionsOverLen;
        nTotalPositionsWithN += chrI->nPositionsWithN;
        if (informing && nNtInform > 0) fprintf(stderr, " (%u Nb)\n", chrI->nPositionsWithN);
        
        //////////////////////////
        /**
        fprintf(stderr, "Checking...\n");
        p = seq->seq.s;
        char check[100];
        char copy[100];
        check[10] = 0;
        copy[10] = 0;
        int k;
        for (n=0; n < m; n++, p++) {
            if (n % 1000000 == 0) {
                extractPackedChromosome(chr_i, n, 10, check, 0, 0);
                for (k=0; k < 10; k++) {
                    copy[k] = *(p+k);
                }
                fprintf(stderr, "%12u: %s == ? %s %d\n",n,copy,check,strcmp(copy,check));
            }
        }
        fprintf(stderr, "\n");
         **/
        //////////////////////////
        
    }
    closeSourceFile();
    
    if (nTotalPositionsWithN > 0) {
        if (informing) {
            fprintf(stderr, "Setting N blocks\n");
            fprintf(stderr, "Seq\t      Size\tName\tN blocks\n");
        }
        chrI = allChromosomes;
        chrI->pPosWithN = (ogChrPosWithN *) malloc(nTotalPositionsWithN * sizeof(ogChrPosWithN));
        ogChrPosWithN *pChrPos = chrI->pPosWithN;
        openSourceFile();
        seq = (kseq_t *) seqPointer;
        for (chr_i=0; chr_i < nChromosomes && kseq_read(seq) >= 0; chr_i++) {
            m = seq->seq.l;
            chrI = allChromosomes + chr_i;
            chrI->pPosWithN = pChrPos; // Even with 0 N pos, this is a valid pointer (not null)
            if (informing && (chr_i <= 50 || chr_i % 1000 == 0)) {
                fprintf(stderr, "%3d\t%10u\t%s\t%u\n", chrI->number, chrI->size, chrI->name, chrI->nPositionsWithN);
            }
            if (chrI->nPositionsWithN > 0) {
                for (n=0, withN=0, p=seq->seq.s; n < m; n++, p++) {
                    if (ascii2pack[*p] == 0 && *p != 'A' && *p != 'a') {
                        withN++;
                    } else {
                        while (withN > 0) {
                            pChrPos->length = (withN > 256 ? 255 : withN-1); // the number of N's is n-1, so this value must be added to 1 to recreate the correct N's
                            pChrPos->position = n - withN;
                            //fprintf(stderr, "Detectando %u N's en chr %u iniciando en %u y hasta %u\n", withN, chr_i, pChrPos->position, pChrPos->position+pChrPos->length);
                            pChrPos++;
                            withN = (withN > 256 ? withN - 256 : 0);
                        }
                        withN = 0;
                    }
                }
                while (withN > 0) {
                    pChrPos->length = (withN > 256 ? 255 : withN-1); // the number of N's is n-1, so this value must be added to 1 to recreate the correct N's
                    pChrPos->position = n - withN;
                    pChrPos++;
                    withN = (withN > 256 ? withN - 256 : 0);
                }
            }
        }
        closeSourceFile();
    }
}



void ogGenome::packSequence(char *seq, char *pack, char offset) { // offset is 0, 1, 2, or 3
    unsigned char        a;
    unsigned char       *p = (unsigned char *) seq;
    unsigned char       *pack_i = (unsigned char *) pack;
    unsigned char        pack_s = offset;
    *pack_i = 0;
    while ((a=*p)) {
        *pack_i |= (ascii2pack[a] >> shift[pack_s]);
        if (pack_s == 3) {
            *++pack_i = 0;
            pack_s = 0;
        } else {
            pack_s++;
        }
        p++;
    }
}

void ogGenome::packSequences(char *seq, char *pack0, char *pack1, char *pack2, char*pack3) {
    unsigned char        a, a2p;
    unsigned char       *p = (unsigned char *) seq;
    unsigned char       *pack_i0 = (unsigned char *) pack0;
    unsigned char       *pack_i1 = (unsigned char *) pack1;
    unsigned char       *pack_i2 = (unsigned char *) pack2;
    unsigned char       *pack_i3 = (unsigned char *) pack3;
    unsigned char        off0 = 0;
    unsigned char        off1 = 1;
    unsigned char        off2 = 2;
    unsigned char        off3 = 3;
    *pack_i0 = *pack_i1 = *pack_i2 = *pack_i3 = 0;
    while ((a=*p)) {
        a2p = ascii2pack[a];
        *pack_i0 |= (a2p >> shift[off0]);
        *pack_i1 |= (a2p >> shift[off1]);
        *pack_i2 |= (a2p >> shift[off2]);
        *pack_i3 |= (a2p >> shift[off3]);
        if (off0 == 3) { *++pack_i0 = 0; off0 = 0; } else { off0++; }
        if (off1 == 3) { *++pack_i1 = 0; off1 = 0; } else { off1++; }
        if (off2 == 3) { *++pack_i2 = 0; off2 = 0; } else { off2++; }
        if (off3 == 3) { *++pack_i3 = 0; off3 = 0; } else { off3++; }
        p++;
    }
}

unsigned char ogGenome::getPackedOffset(uint32_t gPos) { 
    return (gPos & 0x00000003);
}


uint32_t ogGenome::comparePackedSequences(char *pack, uint32_t gPos, uint32_t len, uint32_t maxErrors) { // pack is a pointer to 4 pointer of packed sequences at 0, 1, 2, and 3 offsets
    unsigned char       offset = getPackedOffset(gPos);
    //unsigned char      *pSeq = (unsigned char *) pack[offset];
    unsigned char      *pSeq = (unsigned char *) pack;
    unsigned char      *pGen = packed + (gPos >> 2);
    //unsigned char       xorC;
    uint32_t            errors = 0;

    //for (int i=0; i < 10; i++) {
    //    fprintf(stderr, "%d : %hhu, %hhu, %hhu, %hhu, %hhu :\n", i, *(pGen+i), *(pack[0]+i), *(pack[1]+i), *(pack[2]+i), *(pack[3]+i));
    //}
    
    
    if (len < 4) return 0;
    if (offset > 0) {
        errors += ERRORPACKED[(*pSeq++ ^ *pGen++) & LEFTOFFSETMASK[offset]];
        len = len - 4 + offset;
    }
    if (len >= 4) {
        do { //  len >= 4 && errors <= maxErrors
            if ((errors += ERRORPACKED[(*pSeq++ ^ *pGen++)]) > maxErrors) { len -= 4; break; }
        } while ((len -= 4) >= 4);
    }
    if (len > 0 && errors <= maxErrors) { // 
        errors += ERRORPACKED[(*pSeq ^ *pGen) & RIGHTOFFSETMASK[len]];
    }
    return errors;
}


void ogGenome::setNPosToUnpackChr(uint32_t chr, uint32_t chrPos, uint32_t size, char *unpackstr) {
    uint32_t nPos = allChromosomes[chr].nPositionsWithN;
    if (nPos > 0) {
        uint32_t finChrPos = chrPos + size - 1;
        ogChrPosWithN *pWN = allChromosomes[chr].pPosWithN;
        uint32_t ini = 0;
        uint32_t end = allChromosomes[chr].nPositionsWithN-1;
        uint32_t mid;
        if (chrPos > 0) {
            while (end > ini) {
                mid = (ini + end) >> 1;
                if (pWN[mid].position > chrPos) {
                    if (mid > 0) {
                        end = mid - 1;
                    } else {
                        break;
                    }
                } else if (pWN[mid].position+pWN[mid].length < chrPos) {
                    ini = mid + 1;
                } else {
                    break;
                }
            }
        } else {
            end = 0; // start from first one
        }
        for (mid = (ini + end) >> 1; mid < nPos && pWN[mid].position <= finChrPos; mid++) {
            ini = pWN[mid].position;
            end = ini + pWN[mid].length;
            if (end >= chrPos && ini <= finChrPos) {
                uint32_t i = (ini > chrPos ? ini - chrPos : 0);
                uint32_t f = (end < finChrPos ? end : finChrPos) - chrPos;
                //fprintf(stderr, "Haciendo %u N's en chr %u iniciando en %u y hasta %u\n", f-i+1, chr, chrPos+i, chrPos+f);
                while (i <= f) unpackstr[i++] = 'N';
            }
        }
    }
}


void ogGenome::extractPackedChromosome(uint32_t chr, uint32_t pos, uint32_t size, char *str, char revComp, char lowerCase) {
    if (chr > nChromosomes) {
        fprintf(stderr, "\n****************************************************\n");
        fprintf(stderr, "** Error, accessing chr:%u of max:%u\n", chr, nChromosomes);
        fprintf(stderr, "****************************************************\n");
    }
    int8_t strSign;
    uint8_t agctOffSet;
    unsigned char      *pack_i;
    unsigned char       pack_s;
    pack_i = packed + allChromosomes[chr].initialPack + (pos >> 2);   // eq. pos / 4
    pack_s = allChromosomes[chr].initialBytePos + (pos & 0x00000003);  // eq. pos % 4
    if (pack_s > 3) {
        pack_i++;
        pack_s -= 4;
    }
    extractFromPack_i_s(size, str, revComp, lowerCase, pack_i, pack_s);
    setNPosToUnpackChr(chr, pos, size, str);
}

void ogGenome::printFromGenomicPos(uint32_t gPos, uint32_t size) {
    char *str = (char *) malloc(size+1);
    ogChromosome *chr = getGenomicCoordinate(gPos);
    extractFromGenomicPos(gPos, size, str, 0, 0);
    fprintf(stderr, "Genomic Pos:%u | Chr: %u | %s | Start: %u | Pos: %u | Fragment-Size:%u\n%s\n", gPos, chr->number, chr->name, chr->start, gPos - chr->start, size, str);
    free(str);
}


void ogGenome::extractFromGenomicPos(uint32_t gPos, uint32_t size, char *str, char revComp, char lowerCase) {
    unsigned char      *pack_i;
    unsigned char       pack_s;
    pack_i = packed + (gPos >> 2);
    pack_s = gPos & 0x00000003; 
    extractFromPack_i_s(size, str, revComp, lowerCase, pack_i, pack_s);
    if (nTotalPositionsWithN > 0) {
        ogChromosome *chr = getGenomicCoordinate(gPos);
        setNPosToUnpackChr(chr->number - 1, gPos - chr->start, size, str);
    }
}


void ogGenome::extractFromPack_i_s(uint32_t size, char *str, char revComp, char lowerCase, unsigned char *pack_i, unsigned char pack_s) {
    int8_t strSign;
    uint8_t agctOffSet;
    *(str+size) = 0; // put end of string
    if (revComp) {
        //if (pack_i - packed <= (size >> 2)) { size = ((pack_i - packed + 1) << 2) - 3; } // for protection
        str += size-1;
        agctOffSet = 4;
        strSign = -1;
    } else {
        //if (pack_i + (size >> 2) >= packLimit) { size = ((packLimit - pack_i + 1) << 2) - 3; } // for protection
        agctOffSet = 0;
        strSign = 1;
    }
    if (lowerCase) agctOffSet += 8;
    for (; size > 0; size--) {
        *str = ACGT[agctOffSet + ((*pack_i >> unshift[pack_s]) & 0x03)];
        str += strSign;
        if (pack_s == 3) {
            pack_i++;
            pack_s = 0;
        } else {
            pack_s++;
        }        
    }
}


void ogGenome::extractPackedGenome(uint32_t genpos, uint32_t size, char *str, char revComp, char lowerCase) {
    ogChromosome *chr = getGenomicCoordinate(genpos);
    extractPackedChromosome(chr->number-1, genpos - chr->start, size, str, revComp, lowerCase);
}

//char ogGenome::extractNextNtFromPackedChromosome() {
//    char c = ACGT[(*pack_i >> unshift[pack_s]) & 0x03];
//    if (pack_s == 3) {
//        pack_i++;
//        pack_s = 0;
//    } else {
//        pack_s++;
//    }
//    return c;
//}

uint32_t ogGenome::getGenomeSize() {
    return genomeSize;
}


uint32_t ogGenome::getGenomeValidPositions() {
    return validPositions;
}

uint32_t ogGenome::getGenomeValidPositionsOverLen() {
    return validPositionsOverLen;
}

ogChromosome  *ogGenome::getChromosome(uint32_t i) {
    return allChromosomes + i;
}

ogChromosome  *ogGenome::getGenomicCoordinate(uint32_t genPos) {
    uint32_t min = 0;
    uint32_t max = nChromosomes-1;
    uint32_t mid;
    ogChromosome *chrI;
    while (min < max) {
        //fprintf(stderr, "P:%u,min:%d,max:%d\n",p,min,max);
        mid = (max + min) >> 1;
        chrI = allChromosomes + mid;
        if (genPos > chrI->cummulative) { min = mid+1; }
        else if (genPos < chrI->start) { 
            max = mid - 1; 
        } else {            
            min = mid;
            break;
        }
        //else { max = mid - (genPos <= (chrI->cummulative - chrI->size) ? 1 : 0); }
    }
    return allChromosomes + min;
}

void ogGenome::allocate() {
    if (packed == NULL) {
        fprintf(stderr, "Allocating %.1f MB (%u bytes) of Memory for Packed Genome.\n", (float) memoryPacked / (1024*1024), memoryPacked);
        packed = (unsigned char *) calloc(memoryPacked, 1);
        packLimit = packed + memoryPacked;
    }
    if (allChromosomes == NULL) {
        uint32_t memChr = sizeof(ogChromosome) * nChromosomes;
        fprintf(stderr, "Allocating %u bytes of Memory for Chromosome Info.\n", memChr);
        allChromosomes = (ogChromosome *) malloc(memChr);
    }
}

void ogGenome::loadChromosomesInfo(FILE *pFile) {
    //FILE           *pfile;
    //char            xfile[1000];
    //char            xline[1000];
    int             i;

    //setPureFileName(xfile, sourceFileName, "ogs", 1000);
    //printf_FileOperation("Loading genome structure");
    //fprintf(stderr, "\n");
    //pfile = fopen(xfile, "r");
    //fscanf(pfile, "%[^\n]s ", &xline); // read the line with comment
    fscanf(pFile, "%u\t%u\n", &nChromosomes, &genomeSize);
    memoryPacked = genomeSize / 4 + 1;
    allocate();
    validPositions = 0;
    validPositionsOverLen = 0;
    nTotalPositionsWithN = 0;
    for (i = 0; i < nChromosomes; i++) {
        fscanf(pFile, "%d %[^\t]%*c %u %u %u %hhu %u %u %u %[^\n]s ", 
            &allChromosomes[i].number, allChromosomes[i].name, 
            &allChromosomes[i].size, &allChromosomes[i].cummulative, 
            &allChromosomes[i].initialPack, &allChromosomes[i].initialBytePos,
            &allChromosomes[i].validPositions,  &allChromosomes[i].validPositionsOverLen,
            &allChromosomes[i].nPositionsWithN, allChromosomes[i].comment);
        allChromosomes[i].index = allChromosomes[i].number-1;
        allChromosomes[i].start = (i == 0 ? 0 : allChromosomes[i-1].cummulative);
        validPositions += allChromosomes[i].validPositions;
        validPositionsOverLen += allChromosomes[i].validPositionsOverLen;
        nTotalPositionsWithN += allChromosomes[i].nPositionsWithN;
        /*
        fprintf(stderr, "num:%d, name:%s, size:%u, cum:%u, start:%u.\ncomment:%s\n", 
                allChromosomes[i].number, allChromosomes[i].name, allChromosomes[i].size, 
                allChromosomes[i].cummulative, allChromosomes[i].start, allChromosomes[i].comment);
        */
    }
    //fclose(pfile);
    //printf_FileOperationMsgDone("Loading genome");
    //printf_FileOperationDone();
}
    
void ogGenome::loadPackedGenome(FILE *pFile) {
    int i;
    //setPureFileName(xfile, sourceFileName, "ogp", 1000);
    //printf_FileOperation("Loading packed genome");
    //pfile = fopen(xfile, "r");
    fread(packed, 1, memoryPacked, pFile);
    // Read N positions
    if (nTotalPositionsWithN > 0) {
        allChromosomes[0].pPosWithN = (ogChrPosWithN *) malloc(nTotalPositionsWithN * sizeof(ogChrPosWithN));
        fread(allChromosomes[0].pPosWithN, nTotalPositionsWithN, sizeof(ogChrPosWithN), pFile);
        for (i = 1; i < nChromosomes; i++) {
            allChromosomes[i].pPosWithN = allChromosomes[i-1].pPosWithN + allChromosomes[i-1].nPositionsWithN;
        }
    }
    //fclose(pfile);
    //printf_FileOperationDone();
    
}


uint32_t ogGenome::saveChromosomesInfo(FILE *pFile) {
    //FILE           *pfile;
    char            xfile[MAX_FILENAME];
    kseq_t         *seq;
    int             i;
    ogChromosome   *chrI;
    uint32_t        bytesSaved = 0;

    //setFileName(xfile, pSaveFileName, keySize, "ogs", 1000);
    printf_FileOperation("Saving genome structure");
    //pfile = fopen(xfile, "w");
    bytesSaved += fprintf(pFile, "%d\t%u\n", nChromosomes, genomeSize);
    for (i=0; i < nChromosomes; i++) {
        chrI = allChromosomes+i;
        bytesSaved += fprintf(pFile, "%d\t%s\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n%s\n", 
                chrI->number, chrI->name, chrI->size, 
                chrI->cummulative, chrI->initialPack, 
                chrI->initialBytePos, chrI->validPositions, chrI->validPositionsOverLen,
                chrI->nPositionsWithN, chrI->comment);
    }
    //fclose(pfile);
    printf_FileOperationDone();
    return bytesSaved;    
}

uint32_t ogGenome::savePackedGenome(FILE *pFile) {
    uint32_t bytesSaved = 0;
    //setFileName(xfile, pSaveFileName, keySize, "ogp", 1000);
    //sprintf(xfile, "%s.%d.ogp", sourceFileName, keySize);
    printf_FileOperation("Saving packed genome");
    //pfile = fopen(xfile, "wb");
    //fwrite(packed, memoryPacked, 1, pfile);
    bytesSaved += memoryPacked * fwrite(packed, memoryPacked, 1, pFile);
    if (nTotalPositionsWithN > 0) {
        bytesSaved += nTotalPositionsWithN * fwrite(allChromosomes[0].pPosWithN, nTotalPositionsWithN, sizeof(ogChrPosWithN), pFile);
    }
    //fclose(pfile);
    printf_FileOperationDone();

    return bytesSaved;
}

uint32_t ogGenome::getNChromosomes() {
    return nChromosomes;
}


uint32_t ogGenome::getChromosomeSize(uint32_t iChromosome) {
    return allChromosomes[iChromosome].size;
}

ogChromosome *ogGenome::getLargestChromosome() {
    uint32_t     maxSeqLen = 0, mx=0;
    uint32_t     n;
    for (n=0; n < nChromosomes; n++) {
        if (maxSeqLen < allChromosomes[n].size) {
            maxSeqLen = allChromosomes[n].size;
            mx = n;
        }
    }
    return allChromosomes+mx;
}

uint32_t ogGenome::getChromosomeStart(uint32_t iChromosome) {
    return allChromosomes[iChromosome].start;
}

uint32_t ogGenome::getChromosomeValidPositions(uint32_t iChromosome) {
   return allChromosomes[iChromosome].validPositions;
}

void ogGenome::printGenomePositionInfo(uint32_t gp) {
    uint32_t i;
    ogChromosome *pChr;
    char     dna[10000];
    dna[keySize] = 0;
    fprintf(stderr, "\n...........................................................................................\n");
    pChr = getGenomicCoordinate(gp);
    extractPackedGenome(gp, keySize, dna, 0, 0);
    fprintf(stderr, "GenPos:%u, Chr:%s, ChrPos:%u, DNA:%s\n", gp, pChr->name, gp - pChr->start, dna);
}