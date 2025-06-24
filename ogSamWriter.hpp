/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.h to edit this template
 */

/* 
 * File:   ogSamWriter.hpp
 * Author: victortrevino
 *
 * Created on December 21, 2022, 1:54 PM
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "ogGenome.hpp"
#include "ogDefinitions.hpp"
#include "ogGuider.hpp"
#include "ogKeyEncoding.hpp"
#include "ogSingleRead.hpp"
#include "BGZF.h"
#include <mutex>


#ifndef OGSAMWRITER_HPP
#define OGSAMWRITER_HPP

#define MAX_BAM_BUFFER  16777216     //1024*1024*16
#define MAX_BAM_USED    14680064    //*(128-16)/128 = 87.5% of MAX_BAM_BUFFER
#define MAX_SAM_FILENAME 512

class ogSamWriter {
    
    FILE                    *pFile;
    BamTools::BgzfData      *pBamFile;
    char                    bamBuffer[MAX_BAM_BUFFER];
    uint32_t                bamUsed;
    std::mutex              mtx;

public:
    char                    outFileName[MAX_SAM_FILENAME];
    int                     tempId;
    char                    isBAM;
    char                    unmapMode;
    FILE                    *pUnmapFile1;
    FILE                    *pUnmapFile2;
    char                    readGroup[MAX_SAM_FILENAME];
    static uint32_t         samFileBuffer;

    ogSamWriter();
    virtual ~ogSamWriter();
    
    char    openFile(char *filename);
    char    writeHeaderChromosomes(ogGenome *pGenome);
    char    writeHeaderProgram(ogGuider *pGuider, ogKeyEncoding *pEncoding, char *commandLine, char *sourceFile1, char *sourceFile2, char *rg);
    char    writeSAMRecord(ogSAM *pSAM);
    char    writeSAMInfo(char *qname, int flags, char *rname, uint64_t pos, int mapq, char *cigar, char *rnext, uint64_t posNext, uint32_t tlen, char *seq, char *qual, char *annot, uint64_t readNum, char aliFunc, ogSingleRead *rd);
    char    closeFile();
    void    setIsBAM(char bam);
    void    sam_print_sam(int flags, const char *format, ...);
    void    assignTemporalFile(char *filename, int tempNumber, char isBAMfile);
    void    setUnmapMode(FILE *pFileUnmap1, FILE *pFileUnmap2);

private:

};

#endif /* OGSAMWRITER_HPP */

