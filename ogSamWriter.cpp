/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.cc to edit this template
 */

/* 
 * File:   ogSamWriter.cpp
 * Author: victortrevino
 * 
 * Created on December 21, 2022, 1:54 PM
 */

#include "ogSamWriter.hpp"
#include <stdarg.h>

uint32_t    ogSamWriter::samFileBuffer = 33554432; // 33554432=32 Mb , 134217728=128 Mb default

ogSamWriter::ogSamWriter() {
    pFile = NULL;
    pBamFile = NULL;
    bamUsed = 0;
    isBAM = 0;
    *outFileName = 0;
    tempId = -1;
    unmapMode = 0;
    *readGroup = 0;
}

ogSamWriter::~ogSamWriter() {
    closeFile();
}

char ogSamWriter::openFile(char *filename) {
    closeFile();
    bamUsed = 0;
    if (filename == NULL  ||  strcmp(filename, "stdout") == 0) {
        pFile = stdout;
        *outFileName = 0;
    } else if (strcmp(filename,"NULL") == 0) {
        pFile = NULL;
        *outFileName = 0;
    } else {
        strncpy(outFileName, filename, MAX_SAM_FILENAME);
        int l = strlen(filename);
        if (l > 4) {
            /**
            if ((filename[l-3] == 'B' || filename[l-3] == 'b') && 
                    (filename[l-2] == 'A' || filename[l-2] == 'a') &&
                    (filename[l-1] == 'M' || filename[l-1] == 'm') && 
                    (filename[l-4] == '.') ) {
                isBAM = 1;
            }
             **/
            if ((filename[l-2] == 'G' || filename[l-2] == 'g') &&
                    (filename[l-1] == 'Z' || filename[l-1] == 'z') && 
                    (filename[l-3] == '.') ) {
                isBAM = 1;
            }
        }
        if (isBAM) {
            pBamFile = new BamTools::BgzfData();
            pBamFile->Open(filename, "wb", false);
        } else {
            pFile = fopen(filename, "w");
            setvbuf(pFile, NULL, _IOFBF, ogSamWriter::samFileBuffer); // 128 Mb default
        }
    }
    return (pFile != NULL);
}


void ogSamWriter::assignTemporalFile(char *filename, int tempNumber, char isBAMfile) {
    char newFileName[MAX_SAM_FILENAME];
    closeFile();
    isBAM = isBAMfile;
    tempId = tempNumber;
    strncpy(newFileName, filename, MAX_SAM_FILENAME);
    int n = strlen(filename);
    if (tempId > 0) snprintf(newFileName + n, MAX_SAM_FILENAME - n, "_%d", tempId);
    openFile(newFileName);
}

char ogSamWriter::closeFile() {
    if (pFile != NULL) {
        //sam_print_sam("<<<< CLOSING THE FUCKING FILE >>>>\n");
        fflush(pFile);
        if (pFile != stdout) fclose(pFile);
        pFile = NULL;
    } else if (isBAM && pBamFile != NULL) {
        pBamFile->Close();
    }
    return 1;
}

void ogSamWriter::setIsBAM(char bam) {
    isBAM = bam;
}

void ogSamWriter::sam_print_sam(int flags, const char *format, ...) {
    va_list args;
    va_start(args, format);
    if ((unmapMode) && ((flags & 0x0C) != 0) && (unmapMode != 2 || (flags & 0x0C) == 0x0C)) {
        // flags 0x4 = unmapped rd 1, 0x8 = unmapped rd 2
        // The read/pair is unmapped
        if (pUnmapFile1 != NULL || pUnmapFile2 != NULL) {
            const char *name, *seq, *qual;
            name = va_arg(args, const char *);
            va_arg(args, int); // flag
            va_arg(args, const char *); // ref name
            va_arg(args, int); // pos
            va_arg(args, int); // mapq
            va_arg(args, const char *); // cigar
            va_arg(args, const char *); // ref name next
            va_arg(args, int); // pos next
            va_arg(args, int); // length
            seq = va_arg(args, const char *); // sequence
            qual = va_arg(args, const char *); // quality
            if ((flags & 0x40)) {
                // First in Pair: Rd 1
                fprintf(pUnmapFile1, "@%s\n%s\n+\n%s\n", name, seq, qual);
            } else {
                // Second in Pair: Rd 2
                fprintf(pUnmapFile2, "@%s\n%s\n+\n%s\n", name, seq, qual);
            }
        }
    } else if (isBAM) {
        if (tempId < 0) mtx.lock();
            uint32_t bufSize = vsnprintf(bamBuffer+bamUsed, MAX_BAM_BUFFER, format, args);
            //fprintf(stderr, "[%d]:[%s]\n", buf, bamBuffer); fflush(stderr);
            bamUsed += bufSize;            
            if (bamUsed >= MAX_BAM_USED) {
                pBamFile->Write(bamBuffer, bamUsed);
                bamUsed = 0;
            }
        if (tempId < 0) mtx.unlock();
    } else {
        vfprintf(pFile, format, args);
    }
    va_end(args);
}

char ogSamWriter::writeHeaderChromosomes(ogGenome *pGenome) {
    uint32_t nchr = pGenome->getNChromosomes();
    uint32_t i;
    ogChromosome *pChr;
    for (i=0; i < nchr; i++) {
        pChr = pGenome->getChromosome(i);
        sam_print_sam(0, "@SQ\tSN:%s\tLN:%u\n", pChr->name, pChr->size);
    }
    return 1;
}

char ogSamWriter::writeHeaderProgram(ogGuider *pGuider, ogKeyEncoding *pEncoding, char *commandLine, char *sourceFile1, char *sourceFile2, char *rg) {
    
    if (pGuider != NULL) sam_print_sam(0, "@CO\togMapperGUIDER:%s + %s\n", pGuider->getName(), pGuider->getConfigFile());
    if (pEncoding != NULL) sam_print_sam(0, "@CO\togMapperENCODING:%s\n", pEncoding->getName());
    if (sourceFile1 != NULL) sam_print_sam(0, "@CO\togMapperSOURCEFILE:%s\n", sourceFile1);
    if (sourceFile2 != NULL) sam_print_sam(0, "@CO\togMapperSOURCEFILE:%s\n", sourceFile2);
    if (rg != NULL) sam_print_sam(0, "%s\n", rg);
    //if (commandLine != NULL) sam_print_sam("@PG CL:%s\n", commandLine);
    sam_print_sam(0, "@PG\tID:ogMapper\tPN:ogMapper\tVN:%s\tCL:%s\n", OGMAPPER_VERSION, commandLine);
    
    return 1;
}

char ogSamWriter::writeSAMRecord(ogSAM *pSAM) {
    
    if (pSAM != NULL) {
        sam_print_sam(pSAM->flag,
                "%s\t%u\t%s\t%llu\t%d\t%s\t%s\t%llu\t%d\t%s\t%s\n",
                pSAM->qname,
                pSAM->flag,
                pSAM->rname,
                pSAM->pos,
                pSAM->mapq,
                pSAM->cigar,
                pSAM->rnext == NULL ? ASTERISK : pSAM->rnext,
                pSAM->pnext,
                pSAM->tlen,
                pSAM->seq,
                pSAM->qual  == NULL ? ASTERISK : pSAM->qual
                );
    }
    
    return 1;
}

char ogSamWriter::writeSAMInfo(char *qname, int flags, char *rname, uint64_t pos, int mapq, char *cigar, char *rnext, uint64_t posNext, uint32_t tlen, char *seq, char *qual, char *annot, uint64_t readNum, char aliFunc, ogSingleRead *rd) {

    //if (posNext > 1000000000) {
        // 1 253 408 782
    //    sam_print_sam(stderr, "Pelossssssss\n");
    //}
    
    /**
    uint32_t cigLen = 0;
    uint32_t seqLen = strlen(seq);
    long int l = 0;
    char *pC = cigar;
    char *pA = cigar;
    if (*seq == '*' || *pC == '*') 
        cigLen = seqLen; 
    else
        for(; *pC != 0; pC++) {
            if (*pC > '9' || *pC < '0') {
                if (*pC == 'M' || *pC == 'I' || *pC == 'S' || *pC == '=' || *pC == 'X') {
                    cigLen += l; //strtol(pA, &pA, 10);
                }
                //pA = pC + 1;
                l = 0;
            } else {
                l = l * 10 + *pC - '0';
            }
        }
    if (cigLen != seqLen) {
        sam_print_sam(pFile,"**** NEXT LINE CIGAR LEN= %u, SEQ LEN=%u ****\n", cigLen, seqLen);
    }
     **/
    
    
    
    sam_print_sam(flags,
            "%s\t%u\t%s\t%llu\t%d\t%s\t%s\t%llu\t%d\t%s\t%s\t%srn:i:%llu\taf:A:%c\tcr:Z:%uM,%uX,%uS,%uI,%uD,%soc\t%s\n",
            qname,
            flags,
            rname,
            pos,
            mapq,
            cigar,
            rnext == NULL ? ASTERISK : rnext,
            posNext,
            tlen,
            seq,
            qual  == NULL ? ASTERISK : qual,
            //annot == NULL ? EMPTY : TAB,
            readGroup,
            readNum+1, // Starts on 1
            aliFunc < 32 ? '*' : aliFunc,
            //rd->cigarLeft,rd->cigarLeftType, : %uL%c
            rd->cigarMatches,
            rd->cigarNoMatches,
            rd->cigarSoft,
            rd->cigarIns,
            rd->cigarDel,
            rd->pCigarOriginal,
            annot == NULL ? EMPTY : annot
            );
    return 1;
}

void ogSamWriter::setUnmapMode(FILE *pFileUnmap1, FILE *pFileUnmap2) {
    unmapMode = (pFileUnmap1 != NULL || pFileUnmap2 != NULL);
    pUnmapFile1 = pFileUnmap1;
    pUnmapFile2 = pFileUnmap2;
}
