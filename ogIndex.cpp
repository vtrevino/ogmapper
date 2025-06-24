/* 
 * File:   ogIndex.cpp
 * Author: victortrevino
 * 
 * Created on May 7, 2022, 1:59 PM
 */
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <zlib.h>
//#include <thread>
//#include <chrono>
//#include <sstream>
#include <iostream>
//#include <map>
#include "ogDefinitions.hpp"

#include "ogIndex.hpp"
#include "ogGenome.hpp"
#include "ogKeyEncoding.hpp"
#include "ogGenomePositions.hpp"
#include "ogKeys.hpp"
#include "ogReadsMapper.hpp"
#include "ogPlainEncoding.hpp"
#include "ogStateMachineGuider.hpp"
#include "SimpleStatisticsCalculator.hpp"
#include "ogSamWriter.hpp"

#include "CircularArray.hpp"
#include "CircularArray.cpp"    // Necesita el .CPP para reconstruir el template con ogSingleRead*

#include "kseq.h"
#include "ogStateMachineGuider.hpp"
KSEQ_INIT(gzFile, gzread)

//using namespace std;

// taken from: 
// https://stackoverflow.com/questions/47116974/remove-a-substring-from-a-string-in-c
char *strremove(char *str, const char *sub) {
    char *p, *q, *r;
    if (*sub && (q = r = strstr(str, sub)) != NULL) {
        size_t len = strlen(sub);
        while ((r = strstr(p = r + len, sub)) != NULL) {
            while (p < r)
                *q++ = *p++;
        }
        while ((*q++ = *p++) != '\0')
            continue;
    }
    return str;
}

void copyAndRemoveClassicExtensions(char *pSourceFileName, char *purifiedSourceFileName, int maxLen, const char *pExt, const char *pExt2) {
    strncpy(purifiedSourceFileName, pSourceFileName, maxLen);
    strremove(purifiedSourceFileName, ".FASTQ");
    strremove(purifiedSourceFileName, ".fastq");
    strremove(purifiedSourceFileName, ".FASTA");
    strremove(purifiedSourceFileName, ".fasta");
    strremove(purifiedSourceFileName, ".FQ");
    strremove(purifiedSourceFileName, ".fq");
    strremove(purifiedSourceFileName, ".FA");
    strremove(purifiedSourceFileName, ".fa");
    strremove(purifiedSourceFileName, ".GZIP");
    strremove(purifiedSourceFileName, ".gzip");
    strremove(purifiedSourceFileName, ".GZ");
    strremove(purifiedSourceFileName, ".gz");
    strremove(purifiedSourceFileName, ".GTF");
    strremove(purifiedSourceFileName, ".gtf");
    strremove(purifiedSourceFileName, ".GFF");
    strremove(purifiedSourceFileName, ".gff");
    if (pExt != NULL) strncat(purifiedSourceFileName, pExt, maxLen);
    if (pExt2 != NULL) strncat(purifiedSourceFileName, pExt2, maxLen);
}


ogIndex::ogIndex() {
    doPrepare = 1;
    doProcess = 1;
    kseqReadingMode = 1;
    selEncoding[0] = 0;
    selGuider[0] = 0;
    lowMemoryMode = 0;
    startingRead = 0;
    threads = 1;
    pKeys = NULL;
    pGenome = NULL;
    pEncoding = NULL;
    pGenomePositions = NULL;
    pGuider = NULL;
    strncpy(schedule, "ED", MAX_SCHEDULE_SIZE);
    maxReadLen = MAX_DEF_READ_LEN;
    maxReads = 0xFFFFFFFFFFFFFFFF;
    pMapParams = (ogMappingParameters *) malloc(sizeof(ogMappingParameters));
    //
    // Hay que pasar estos a una estructura nueva USER_PARAMETERS
    pMapParams->MAX_KEY_DISTANCE = 3; // 3 DNA, 100000 RNA
    //
    pMapParams->countExons = 0;
    pMapParams->countTranscripts = 0;
    pMapParams->outputSAMcount = 0;
    pMapParams->MIN_KEY_COUNTS = 4;
    pMapParams->MIN_SCORE_COUNT = 10;
    pMapParams->BREAK_TIE_COUNT1 = 1;
    pMapParams->MIN_MATCHES = 50;
    //
    nThreads = 0;
    nQueue = 10000;
    forceUpperCase = 0;
    productionMode = 1;
    pSamSam = NULL;
    countInfo.GTFfile[0] = 0;
    countInfo.OGXfile[0] = 0;
    countInfo.pExonReads = NULL;
    countInfo.pExons = NULL;
    countInfo.pGeneReads = NULL;
    countInfo.pTranscripts = NULL;
    countInfo.pGeneExonStart = NULL;
    countInfo.pGeneExon = NULL;
    strcpy(genesCountFile, "gene-counts.txt");
    strcpy(exonsCountFile, "exon-counts.txt");
    strcpy(transcriptsCountFile, "transcript-counts.txt");
    isExonMode = 0;
    maxReadsBuffer = 1024 * 16;
    filePerThread = 0;
    unmappedMode = 0;
    *readGroupHeader = 0;
}

void freef(void *p, char msg) {
    fprintf(stderr, "<%c", msg); fflush(stderr);
    free(p);
    fprintf(stderr, ">"); fflush(stderr);
}

ogIndex::~ogIndex() {
    //fprintf(stderr, "<Deallocating ogIndex:");
    //fprintf(stderr, "<~ogIdx:"); fflush(stderr);
    // Deallocate ??
    if (pKeys != NULL) delete pKeys;
    if (pGenome != NULL) delete pGenome;
    if (pGenomePositions != NULL) delete pGenomePositions;
    //if (pGuider != NULL) free(pGuider);
    //if (pEncoding != NULL) delete pEncoding;
    if (pMapParams != NULL) free(pMapParams);
    deleteKeyEncodings();
    deleteGuiders();
    if (countInfo.pExons != NULL) freef(countInfo.pExons, 'X');
    if (countInfo.pTranscripts != NULL) freef(countInfo.pTranscripts, 'T');
    if (countInfo.pExonReads != NULL) freef(countInfo.pExonReads, 'x');
    if (countInfo.pGeneReads != NULL) freef(countInfo.pGeneReads, 'g');
    if (countInfo.pGeneExonStart != NULL) freef(countInfo.pGeneExonStart, 's');
    if (countInfo.pGeneExon != NULL) freef(countInfo.pGeneExon, 'G');
    //fprintf(stderr, ":ogIdx>\n"); fflush(stderr);
}

void ogIndex::setThreads(uint16_t maxThreads) {
    nThreads = maxThreads;
}

void ogIndex::setQueue(uint16_t maxQueue) {
    nQueue = maxQueue;
}

void ogIndex::registerKeyEncoding(ogKeyEncoding *pEncodingClass) {
    std::string s (pEncodingClass->getName());
    protoKEtable[s] = pEncodingClass;
//    fprintf(stderr, "Registered Encodings:\n");
//    std::map<std::string, ogKeyEncoding*>::iterator it = protoKEtable.begin();
//    while (it != protoKEtable.end()) {
//        std::cout << it->first << ':' << it->second << '\n';
//        it++;
//    }
}

void ogIndex::registerGuider(ogGuider *pGuiderClass) {
    std::string s (pGuiderClass->getName());
    protoGuiderTable[s] = pGuiderClass;
//    fprintf(stderr, "Registered Guiders:\n");
//    std::map<std::string, ogGuider*>::iterator it = protoGuiderTable.begin();
//    while (it != protoGuiderTable.end()) {
//        std::cout << it->first << ':' << it->second << '\n';
//        it++;
//    }
}

void ogIndex::deleteKeyEncodings() {
    uint16_t nKE = protoKEtable.size();
    //fprintf(stderr, "<Deallocating %hu key encodings:", nKE);
    //fprintf(stderr, "<~ %hu ke:", nKE); fflush(stderr);
    std::map<std::string, ogKeyEncoding*>::iterator it = protoKEtable.begin();
    while (it != protoKEtable.end()) {
        //fprintf(stderr, "Deallocating encoding, key [%s], name [%s]...\n", (it->first).c_str(),  it->second->getName());
        delete it->second;
        it++;
    }
    //fprintf(stderr, "Ending deallocating %lu key encodings.\n", protoKEtable.size());
    //fprintf(stderr, ":%hu ke>",nKE);
}

void ogIndex::deleteGuiders() {
    uint16_t nG = protoGuiderTable.size();
    //fprintf(stderr, "<Deallocating %hu guiders:", nG);
    //fprintf(stderr, "<~ %hu g:", nG); fflush(stderr);
    std::map<std::string, ogGuider*>::iterator it = protoGuiderTable.begin();
    while (it != protoGuiderTable.end()) {
        //fprintf(stderr, "Deallocating encoding, key [%s], name [%s]...\n", (it->first).c_str(),  it->second->getName());
        delete it->second;
        it++;
    }
    //fprintf(stderr, "Ending deallocating %lu guiders.\n", protoGuiderTable.size());
    //fprintf(stderr, ":%hu g>",nG);
}

ogKeyEncoding* ogIndex::getKeyEncoding(char *pEncoding) {
    std::string encoding = "";
    if (pEncoding != NULL) encoding = std::string(pEncoding);
    std::map<std::string, ogKeyEncoding*>::iterator it = protoKEtable.begin();
    if (encoding.size() >  0) {
        it = protoKEtable.find(encoding);
    }
    if (encoding.size() == 0 && selEncoding[0] > 0) {
        it = protoKEtable.find(std::string(selEncoding));
    }
    if (it == protoKEtable.end()) {
        return NULL;
    }
    return it->second; // field second is the value of the map (first is the key)
}

ogGuider* ogIndex::getGuider(char *pGuider) {
    std::string guider = "";
    if (pGuider != NULL) guider = std::string(pGuider);
    std::map<std::string, ogGuider*>::iterator it = protoGuiderTable.begin();
    if (guider.size() >  0) {
        it = protoGuiderTable.find(guider);
    }
    if (guider.size() == 0 && selGuider[0] > 0) {
        it = protoGuiderTable.find(std::string(selGuider));
    }
    if (it == protoGuiderTable.end()) {
        return NULL;
    }
    ogGuider *pGd = it->second;
    pGd->setConfigFile(guideFile);
    return pGd; // field second is the value of the map (first is the key)
}

void ogIndex::selectKeyEncoding(char *pEncoding) {
    strncpy(selEncoding, pEncoding, MAX_CHAR_NAMES);
}

void ogIndex::selectGuider(char *pGuider) {
    strncpy(selGuider, "", MAX_CHAR_NAMES); // fill zeros
    char *p = pGuider;
    while (*p != 0 && *p != ':') p++;
    strncpy(selGuider, pGuider, (p-pGuider < MAX_CHAR_NAMES ? p-pGuider : MAX_CHAR_NAMES));// MAX_CHAR_NAMES
    if (*p == ':') {
        strncpy(guideFile, p+1, MAX_FILENAME);
    } else {
        strncpy(guideFile, "NO-GUIDE-FILE", MAX_FILENAME);
    }
}

void ogIndex::setMaxReadLength(uint32_t mxRdLen) {
    maxReadLen = mxRdLen;
}

void ogIndex::setStartingRead(uint64_t stRead) {
    startingRead = stRead;
}

void ogIndex::setMaxReads(uint64_t mxRds) {
    maxReads = mxRds;
}

uint32_t ogIndex::getMeanKeyPlusStdDev(float nStdDev) {
    double n_log10 = header.meanKeys_log10 + nStdDev * header.sdKeys_log10;
    double n = pow(10, n_log10);
    return (uint32_t) n;
}

//// v0.5
////#define _FILE_HEAD_POINT     2048
////#define _FILE_START_POINT    4096
// v0.6
#define FILE_HEAD_POINT     6144
#define FILE_START_POINT    8192

FILE *ogIndex::openIndexFileAtSection(char section) {
/////////////// v0.05
// Structure of the OGMAPPER INDEX FILE    
// Bytes    0 to 2047 : Info about the Index for the user (user may type more <index> file.
// Bytes 2048 to 4096 : Info about the Index in the following format:
//                      "8-bytes" of two uint32_t for Size in bytes and "Starting" byte 
//                      (sizeof operation is used for uint32_t)
//                      Starting byte is relative to 0
//               2048 : For the old .ogs file having list of all chromosomes.
//               2056 : For the old .ogp file having the packed genome.
//               2064 : For the old .ogk file having the keys.
//               2072 : For the old .ogo file having the key positions.
//               2080 : For the old .ogg file having the header of the Index.
//               2088 : For the old .ogt file having the tuple info.
//
/////////////// v0.06
// Structure of the OGMAPPER INDEX FILE    
// Bytes    0 to 6144 : Info about the Index for the user (user may type more <index> file.
// Bytes 6144 to 6182 : Info about the Index in the following format:
//                      "8-bytes" of two uint32_t for Size in bytes and "Starting" byte 
//                      (sizeof operation is used for uint32_t)
//                      Starting byte is relative to 0
//               6144 : For the old .ogs file having list of all chromosomes.
//               6152 : For the old .ogp file having the packed genome.
//               6160 : For the old .ogk file having the keys.
//               6168 : For the old .ogo file having the key positions.
//               6176 : For the old .ogg file having the header of the Index.
//               6182 : For the old .ogt file having the tuple info.
//
/////////////// v0.07
// Structure of the OGMAPPER INDEX FILE    
// Bytes    0 to 6144 : Info about the Index for the user (user may type more <index> file.
// Bytes 6144 to 6224 : Info about the Index in the following format:
//                      "16-bytes" of two uint64_t for Size in bytes and "Starting" byte 
//                      (sizeof operation is used for uint64_t)
//                      Starting byte is relative to 0
//               6144 : For the old .ogs file having list of all chromosomes.
//               6160 : For the old .ogp file having the packed genome.
//               6176 : For the old .ogk file having the keys.
//               6192 : For the old .ogo file having the key positions.
//               6208 : For the old .ogg file having the header of the Index.
//               6224 : For the old .ogt file having the tuple info.
//
    FILE *pfile;    
    uint64_t    uis = sizeof(uint64_t);
    uint64_t    readPositions[2];
    
    if (section == 0) {
        // Create
        int i;
        pfile = fopen(indexFileName, "w+b");
        for (i=0; i < FILE_HEAD_POINT/80; i++) fwrite("         .         .         .         .         .         .         .         \n", 80, 1, pfile);
        fseek(pfile, FILE_HEAD_POINT, SEEK_SET);
        readPositions[0] = 0;
        readPositions[1] = FILE_START_POINT;
        for (i=0; i < 32; i++) fwrite(&readPositions, uis, 2, pfile);
        fseek(pfile, 0, SEEK_SET);
    } else {
        pfile = fopen(indexFileName, "r+b");
        fseek(pfile, FILE_HEAD_POINT + uis * 2 * (section-1), SEEK_SET);
        fread(&readPositions, uis, 2, pfile);
        fseek(pfile, readPositions[1], SEEK_SET);
        //fprintf(stderr, ")) Opening Index at Section %d. Physical Position = %u, Size =%u ((\n", (uint32_t) section, readPositions[1], readPositions[0]);
    }
    return pfile;
}

void ogIndex::setSectionSizeInIndexFile(char section, uint64_t size) {
    FILE *pfile;    
    uint64_t    uis = sizeof(uint64_t);
    
    pfile = fopen(indexFileName, "r+b");
    fseek(pfile, FILE_HEAD_POINT + 2*uis*(section-1), SEEK_SET);
    fwrite(&size, uis, 1, pfile);
    fclose(pfile);
}

void ogIndex::setSectionStartInIndexFile(char section, uint64_t start) {
    FILE *pfile;    
    uint64_t    uis = sizeof(uint64_t);
    uint64_t    begin = start + FILE_START_POINT;
    
    pfile = fopen(indexFileName, "r+b");
    fseek(pfile, FILE_HEAD_POINT + 2*uis*(section-1) + uis, SEEK_SET);
    fwrite(&begin, uis, 1, pfile);
    fclose(pfile);
}

void ogIndex::setSectionSizeStartInIndexFile(char section, uint64_t size, uint64_t start) {
    FILE *pfile;    
    uint64_t    uis = sizeof(uint64_t);
    uint64_t    begin = start + FILE_START_POINT;
    
    pfile = fopen(indexFileName, "r+b");
    fseek(pfile, FILE_HEAD_POINT + 2*uis*(section-1), SEEK_SET);
    fwrite(&size, uis, 1, pfile);
    fwrite(&begin, uis, 1, pfile);
    fclose(pfile);
}


void ogIndex::generate(char *pSourceFileName, uint16_t userKeySize, char *pDestFileName, int argc, char *argv[]) {


    uint32_t genomeSize;
    uint32_t largest;
    uint32_t slen;
    uint32_t n, m;
    uint32_t idx;
    gzFile fp;
    int i;
    char    purifiedSourceFileName[MAX_FILENAME];
    FILE    *indexFile;
    
    
    pEncoding = getKeyEncoding(NULL);
    //pEncoding = new ogPlainEncoding();
    if (pEncoding == NULL) {
        fprintf(stderr, "NO KEY ENCODING SPECIFIED, NOT FOUND, OR NONE REGISTERED.\n");
        return;
    }
    pGuider = getGuider(NULL);
    //pGuider = new ogStateMachineGuider();
    if (pGuider == NULL) {
        fprintf(stderr, "NO GUIDER SPECIFIED, NOT FOUND, OR NONE REGISTERED.\n");
        return;
    }
    pGuider->setGenomeMode(1);

    pEncoding->setSizeInChars(userKeySize);
    if (pEncoding->getSizeInBits() < 1) {
        fprintf(stderr, "**** ERROR: Encoding reports 0 bits. ****\n");
        fflush(stderr);
        return;
    }
    if (pEncoding->getSizeInBits() > 27) {
        const char asterisks[] = "**********************************************************************";
        fprintf(stderr, "%s\n", asterisks);
        fprintf(stderr, "%s\n", asterisks);
        fprintf(stderr, "%s\n", asterisks);
        if (pEncoding->getSizeInBits() > 31) {
            fprintf(stderr, "**** ERROR : Encoding reports %d bits MAY REQUIRE A LOT OF MEMORY ****\n", pEncoding->getSizeInBits());
            fprintf(stderr, "%s\n", asterisks);
            fprintf(stderr, "%s\n", asterisks);
            fprintf(stderr, "%s\n", asterisks);
            return;            
        }
        fprintf(stderr, "**** WARNING : Encoding reports %d bits. ****\n", pEncoding->getSizeInBits());
        fprintf(stderr, "%s\n", asterisks);
        fprintf(stderr, "%s\n", asterisks);
        fprintf(stderr, "%s\n", asterisks);
    }
    pGuider->initializeForIndexing();
    pGuider->setKeySizeInChars(pEncoding->getSizeInChars());

    copyAndRemoveClassicExtensions(pSourceFileName, purifiedSourceFileName, MAX_FILENAME, pGuider->getShortExtensionName(), pEncoding->getShortExtensionName());
    strncpy(indexFileName, (*pDestFileName == 0 ? purifiedSourceFileName : pDestFileName), 2048);
    if (*pDestFileName == 0) strncat(indexFileName, ".ogi", 10);
    fprintf(stderr, "Input Genome File: %s\nGuider           : %s:%s\nEncoding         : %s\nUser Key Size    : %d\nKey Size in Bits : %d\nOutput file(s)   : %s\n", 
            pSourceFileName, pGuider->getName(), pGuider->getConfigFile(), pEncoding->getName(), userKeySize, pEncoding->getSizeInBits(), indexFileName);
    
    fprintf(stderr, "==== Index process started ====\n");

    header.nKeyCountsAtQ95 = 200; // Default
    strncpy(header.indexVersion, "oriGen Index v0.8;", MAX_VERSION_NAME);
    strncpy(header.sourceFileName, pSourceFileName, MAX_OGINDEX_FILENAME);
    strncpy(header.encodingName, pEncoding->getName(), MAX_VERSION_NAME);
    strncpy(header.guiderName, pGuider->getName(), MAX_VERSION_NAME);
    header.encodingVersion = pEncoding->version;
    header.userKeySize = userKeySize;
    
    uint64_t subTotalFileSize = 0;

    indexFile = openIndexFileAtSection(0);
    fprintf(indexFile, "** oriGen Project **\nogMapper Index   : %s\n", header.indexVersion);
    fprintf(indexFile, "Input Genome File: %s\nGuider           : %s:%s\nEncoding         : %s\nUser Key Size    : %d\nKey Size in Bits : %d\nOutput file(s)   : %s\n", 
            pSourceFileName, pGuider->getName(), pGuider->getConfigFile(), pEncoding->getName(), userKeySize, pEncoding->getSizeInBits(), indexFileName);
    fprintf(indexFile, "This file is being generated ...\n");
    fclose(indexFile);
    //fprintf(stderr, "S0 Index File [%s] Size = %ld\n", indexFileName, get_file_size(indexFileName));
    
    // **** SAVE GUIDER, section = 1
    setSectionSizeStartInIndexFile(1, 1024, subTotalFileSize);
    indexFile = openIndexFileAtSection(1);                            // 1 = Guider
    uint32_t guiderFileSize = pGuider->save(indexFile);                   // it may save a ogm file
    fclose(indexFile);
    //fprintf(stderr, "S1 Index File [%s] Size = %ld\n", indexFileName, get_file_size(indexFileName));
    setSectionSizeInIndexFile(1, guiderFileSize);
    subTotalFileSize += guiderFileSize;
    
    // **** SAVE ENCODING, section = 2
    setSectionSizeStartInIndexFile(2, 1024, subTotalFileSize);
    indexFile = openIndexFileAtSection(2);                          // 2 = Encoding
    uint32_t encodingFileSize = pEncoding->save(indexFile);              // it may save a oge/ogt
    fclose(indexFile);
    //fprintf(stderr, "S2 Index File [%s] Size = %ld\n", indexFileName, get_file_size(indexFileName));
    setSectionSizeInIndexFile(2, encodingFileSize);
    subTotalFileSize += encodingFileSize;
    
    // **** SAVE GENOME, section = 3
    fprintf(stderr, "1/4 - Reading Genome : Estimating genome/chromosome sizes.\n");
    pGenome = new ogGenome(pSourceFileName, userKeySize);
    pGenome->estimateChromosomeSizes(1); // 1 = verbose
    fprintf(stderr, "2/4 - Reading Genome: Packing genome.\n");
    pGenome->packGenome(pGenome->getGenomeSize() / 500, pEncoding, pGuider);    
    setSectionSizeStartInIndexFile(3, 1024, subTotalFileSize);
    indexFile = openIndexFileAtSection(3);                            // 3 = Genome Structure
    uint32_t chromosomesFileSize = pGenome->saveChromosomesInfo(indexFile);
    fclose(indexFile);
    //fprintf(stderr, "S3 Index File [%s] Size = %ld\n", indexFileName, get_file_size(indexFileName));
    subTotalFileSize += chromosomesFileSize;
    setSectionSizeInIndexFile(3, chromosomesFileSize);

    setSectionSizeStartInIndexFile(4, 1024, subTotalFileSize);
    indexFile = openIndexFileAtSection(4);                            // 4 = Genome Packed
    uint32_t genomeFileSize = pGenome->savePackedGenome(indexFile);
    fclose(indexFile);
    subTotalFileSize += genomeFileSize;
    uint32_t totalKeys = pGenome->getGenomeValidPositions();
    fprintf(stderr, "Valid key positions in Genome: %u.\n", totalKeys);
    
    pGenomePositions = new ogGenomePositions(totalKeys, userKeySize);
    
    pKeys = new ogKeys(pEncoding->getSizeInBits(), userKeySize, pGenomePositions, pEncoding->getTotalKeys());
    pKeys->resetKeySizesAndOffSets();
        

    //char    keyFragment[1000];
    int x = 0;

    char        *pSequence = NULL, *pLargestSeq = NULL;
    uint32_t     maxSeqLen = pGenome->getLargestChromosome()->size+1;
    uint32_t     chrMod = (uint32_t) pow(10,fmax(1,trunc(log10(pGenome->getNChromosomes()))-1));

    fprintf(stderr, "Allocating %.1f MB (%lu bytes) of memory for largest chromosome.\n", ((float) maxSeqLen * sizeof(char))/(1024.0 * 1024.0), maxSeqLen * sizeof(char));
    pLargestSeq = pSequence = (char *) malloc(maxSeqLen * sizeof(char));
    fprintf(stderr, "3/4 - Unpacking Genome: Estimating size for all keys.\n");
    uint32_t mod = totalKeys / 500; //(pGenome->getNChromosomes() * 20);
    uint32_t theKey;
    fprintf(stderr, "Seq\t      Size\tGuides/Keys\n");
    char informing = 0;
    pGenome->openSourceFile();
    for (n=0; n < pGenome->getNChromosomes(); n++) {
        informing = 0;
        if (n < 50 || n % chrMod == 0 || n == pGenome->getNChromosomes()-1) { 
            fprintf(stderr, "%3d\t%10u\t%11u ", n+1, pGenome->getChromosomeSize(n), pGenome->getChromosomeValidPositions(n)); 
            informing = 1; 
        }
        fflush(stderr);
        //slen = pGenome->getChromosomeSize(n) - userKeySize;
        ///// Ya no es necesario extraer porque se vuelve a leer debido a que si hay N, el extractPacked no funciona bien y los guias pueden estar mal estimados.
        ///// En lugar de extractPacked, se lee de pGenome
        pGenome->extractPackedChromosome(n, 0, pGenome->getChromosomeSize(n), pSequence, 0, 0);
        ///// pSequence = pGenome->readSourceFile(); // Lo correcto sería hacer un strncpy hacia pSequence pero ... inge su
        pGuider->setSequence(pSequence, pGenome->getChromosomeSize(n)); // esta siempre esta en upper case porque el extract asi lo hace no necesita fixSequenceCase
        pGuider->fixSequenceCase();
        //fprintf(stderr, "%u ", pEncoding->getValidKeys());
        for (x=0, m=0; pGuider->nextGuide(); m++) {  // && m < slen
            if (informing && ((m % mod) == 0)) {
                fprintf(stderr, "."); fflush(stderr);
                if (++x % 50 == 0) fprintf(stderr, "/\n");
            }
            theKey = pEncoding->getFwdKey(pGuider->getCurrentGuide());
            if (pEncoding->isValidKey()) pKeys->incKey(theKey);
        }
        if (m != pGenome->getChromosomeValidPositions(n)) { 
            fprintf(stderr, " %u vs %u ************** ERROR ************ ",m, pGenome->getChromosomeValidPositions(n)); 
        }
        if (informing) fprintf(stderr, "\n");
    }
    pGenome->closeSourceFile();
    pKeys->estimateOffSets();
    
    // **** SAVE KEYS, section = 5
    setSectionSizeStartInIndexFile(5, 1024, subTotalFileSize);
    indexFile = openIndexFileAtSection(5);                            // 5 = Keys
    uint64_t keysFileSize = pKeys->save(indexFile);
    fclose(indexFile);
    //fprintf(stderr, "S4 Index File [%s] Size = %ld\n", indexFileName, get_file_size(indexFileName));
    setSectionSizeInIndexFile(5, keysFileSize);
    subTotalFileSize += keysFileSize;
    
    pKeys->saveSizes(purifiedSourceFileName);

    
    ///////////////////
    //Check most counts and 0 counts
    uint32_t k = 0, mx = 0;
    uint64_t max = 0;
    uint64_t v;
    uint64_t arrsize = pKeys->getArraySize();
    const uint16_t HKEYS = 138;
    uint32_t cuts[HKEYS] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99,
                        100, 200, 300, 400, 500, 600, 700, 800, 900,
                        1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000,
                        10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000,
                        100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000,
                        1000000, (uint32_t) 0xFFFFFFFF };
    uint32_t counts[HKEYS] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0, 0, 0, 0,
                         0, 0 };
    double sl = 0;
    double l10;
    SimpleStatisticsCalculator stat;
    for (n=0; n < arrsize; n++) {
        v = pKeys->getSizeForKey(n);
        for (k=0; k < HKEYS; k++) if (v <= cuts[k]) { counts[k]++; break; }
        if (v > max) {
            max = v;
            mx = n;
        }
        if (v > 0) stat.enter(log10(v));
    }
    header.meanKeys_log10 = stat.getMean();
    header.sdKeys_log10 = stat.getStandardDeviation();
    fprintf(stderr, "Mean key count (log10) = %f; Std Dev = %f\n",header.meanKeys_log10,header.sdKeys_log10);
    
    // Get Q95 not counting 0
    long tCounts = 0, kCounts = 0;
    for (k=1; k < HKEYS; k++) tCounts += counts[k];
    tCounts = 95 * tCounts / 100;
    for (k=1; k < HKEYS; k++) {
        kCounts += counts[k];
        if (kCounts >= tCounts) {
            header.nKeyCountsAtQ95 = cuts[k];
            break;
        }
    }    
    fprintf(stderr, "Key Count at Q95 = %d\n",header.nKeyCountsAtQ95);

    // Write Header File
    //char xfile[1000];
    //FILE *pfile;
    //setFileName(xfile, purifiedSourceFileName, userKeySize, "ogg", 1000);
    printf_FileOperation("Saving index header");
    //pfile = fopen(xfile, "wb");
    //fwrite(&header, sizeof(header), 1, pfile);
    //fclose(pfile);

    // **** SAVE INDEX HEADER, section = 6
    uint32_t headerFileSize = sizeof(header);
    setSectionSizeStartInIndexFile(6, headerFileSize, subTotalFileSize);
    indexFile = openIndexFileAtSection(6);                            // 6 = Index Header
    fwrite(&header, headerFileSize, 1, indexFile);
    fclose(indexFile);
    //fprintf(stderr, "S5 Index File [%s] Size = %ld\n", indexFileName, get_file_size(indexFileName));
    subTotalFileSize += headerFileSize;

    printf_FileOperationDone();
    
    
    fprintf(stderr, "Maximum genome positions in a key:%llu, Index:%u (0x%x) [",max,mx,mx);
    
    for (k=pEncoding->getSizeInBits(); k > 0; k--) fprintf(stderr, "%u", (mx >> (k-1)) & 0x00000001);

    fprintf(stderr, "]\nDistribution of positions per key:\n");
    for (k=0; k < HKEYS; k++) fprintf(stderr, "<=%u:%u%c", cuts[k], counts[k], (k+1) % 5 == 0 ? '\n' : '\t');
    fprintf(stderr, "\n");
    ////////////////
    
    
    fprintf(stderr, "4/4 Unpacking Genome: Positioning keys.\n");
    char qq;
    
    if (lowMemoryMode) {
        qq='1';
        fprintf(stderr, "____ Using HALF genome memory access mode ____\n");
    } else {
        qq='0'; // Full memory use
        fprintf(stderr, "^^^^ Using FULL genome memory access mode ^^^^\n");
    }

    uint32_t halfKeys = pKeys->getArraySize() >> 1;
    uint32_t genomicPos;
    uint32_t iExon;
    Exon     *pExon;
    if (isExonMode) loadExonAndTranscriptFiles(&countInfo);        
    
    pKeys->resetKeySizes(); // offsets must remain
    fprintf(stderr, "Cycle\tSeq\t      Size\tGuides/Keys\n");
    for (; qq < '3'; qq++) {
        pGenome->openSourceFile();
        uint32_t    gpStart;
        for (n=0; n < pGenome->getNChromosomes(); n++) {
            informing = 0;
            if (n < 50 || n % chrMod == 0 || n == pGenome->getNChromosomes()-1) informing = 1;
            gpStart = pGenome->getChromosomeStart(n);
            if (informing) fprintf(stderr, "  %c  \t%d\t%10u\t%11u ", qq, n+1, pGenome->getChromosomeSize(n), pGenome->getChromosomeValidPositions(n));
            fflush(stderr);
            //slen = gpStart + pGenome->getChromosomeSize(n) - userKeySize;

            ///// Ya no es necesario extraer porque se vuelve a leer debido a que si hay N, el extractPacked no funciona bien y los guias pueden estar mal estimados.
            ///// En lugar de extractPacked, se lee de pGenome
            pGenome->extractPackedChromosome( n, 0, pGenome->getChromosomeSize(n), pSequence, 0, 0);
            ///// pSequence = pGenome->readSourceFile(); // Lo correcto sería hacer un strncpy pero ... inge su
            pGuider->setSequence(pSequence, pGenome->getChromosomeSize(n));
            //fprintf(stderr, "%u ",pEncoding->getValidKeys());
            
            for (x=0, m=0; pGuider->nextGuide(); m++) {
                if (informing && ((m % mod) == 0)) {  // & 0x001FFFFF
                    fprintf(stderr, "."); fflush(stderr);
                    if (++x % 50 == 0) fprintf(stderr, "/\n");
                }
                idx = pEncoding->getFwdKey(pGuider->getCurrentGuide());
                if (pEncoding->isValidKey()) {
                    if (qq == '0' || (qq == '1' && idx < halfKeys) || (qq == '2' && idx >= halfKeys)) {
                        genomicPos = pGuider->getCurrentGuidePosition();
                        if (isExonMode) {
                            iExon = findExonIndexForPseudoPosition(n, genomicPos);
                            pKeys->addPosToKey(iExon, idx);
                            //
                            pExon = countInfo.pExons + iExon;
                            if (pExon->iGene != n || pExon->phyStart > genomicPos || pExon->phyEnd < genomicPos) {
                                fprintf(stderr, "======= Exon does no correspond: Target Chr:Pos=%u:%u, iExon:%u, ChrPos=%u:%u-%u, St:%c, PhyPos:%u-%u, Gene:%u, Exon:%u\n", n, genomicPos, iExon, pExon->chr, pExon->start, pExon->end, pExon->strand, pExon->phyStart, pExon->phyEnd, pExon->iGene, pExon->iExon); fflush(stderr);
                                pExon--;
                                fprintf(stderr, "======= Exon-1 does no correspond: Target Chr:Pos=%u:%u, iExon:%u, ChrPos=%u:%u-%u, St:%c, PhyPos:%u-%u, Gene:%u, Exon:%u\n", n, genomicPos, iExon, pExon->chr, pExon->start, pExon->end, pExon->strand, pExon->phyStart, pExon->phyEnd, pExon->iGene, pExon->iExon); fflush(stderr);
                            }
                        } else {
                            pKeys->addPosToKey(gpStart + genomicPos, idx);
                        }
                    }
                }
            }
            if (m != pGenome->getChromosomeValidPositions(n)) { 
                fprintf(stderr, " %u vs %u *** ERROR ***\n",m, pGenome->getChromosomeValidPositions(n)); 
            } else if (informing) { 
                fprintf(stderr, "\n"); 
            }
        }
        pGenome->closeSourceFile();
        if (qq == '0') break;
    }
    //pKeys->saveGenomePositions(purifiedSourceFileName);
    // **** SAVE GENOME POSITIONS, section = 7
    setSectionSizeStartInIndexFile(7, 1024, subTotalFileSize);
    indexFile = openIndexFileAtSection(7);                            // 7 = Genome Positions
    uint64_t genPosFileSize = pKeys->saveGenomePositions(indexFile);
    //fclose(indexFile);
    //fprintf(stderr, "S6 Index File [%s] Size = %ld\n", indexFileName, get_file_size(indexFileName));
    setSectionSizeInIndexFile(7, genPosFileSize);
    subTotalFileSize += genPosFileSize;
    
    //indexFile = openIndexFileAtSection(7);
    fseek(indexFile, 0, SEEK_SET);
    // Update header
    fprintf(indexFile, "This is a binary file headed with text info - do not modify it.\n");
    fprintf(indexFile, "********************\n");
    fprintf(indexFile, "** oriGen Project **\n");
    fprintf(indexFile, "********************\n");
    fprintf(indexFile, "ogMapper Index   : %s\n", header.indexVersion);
    fprintf(indexFile, "Input Genome File: %s\nGuider           : %s:%s\nEncoding         : %s\nUser Key Size    : %d\nKey Size in Bits : %d\nOutput file(s)   : %s\n", 
            pSourceFileName, pGuider->getName(), pGuider->getConfigFile(), pEncoding->getName(), userKeySize, pEncoding->getSizeInBits(), indexFileName);
    fprintf(indexFile, "== FILE SIZES ==\n");
    fprintf(indexFile, "Guider          : %u (%.1f MB)\n", guiderFileSize, (float) guiderFileSize / (1024.0*1024.0));
    fprintf(indexFile, "Encoding        : %u (%.1f MB)\n", encodingFileSize, (float) encodingFileSize / (1024.0*1024.0));
    fprintf(indexFile, "Genome (packed) : %u (%.1f MB)\n", genomeFileSize, (float) genomeFileSize / (1024.0*1024.0));
    fprintf(indexFile, "Keys            : %llu (%.1f MB)\n", keysFileSize, (float) keysFileSize / (1024.0*1024.0));
    fprintf(indexFile, "Index Spec      : %u (%.1f MB)\n", headerFileSize, (float) headerFileSize / (1024.0*1024.0));
    fprintf(indexFile, "Genome Positions: %llu (%.1f MB)\n", genPosFileSize, (float) genPosFileSize / (1024.0*1024.0));
    fprintf(indexFile, "Header (this)   : %u (%.1f MB)\n", FILE_START_POINT, (float) FILE_START_POINT / (1024.0*1024.0));
    fprintf(indexFile, "Total           : %llu (%.1f MB)\n", FILE_START_POINT+subTotalFileSize, (float) (FILE_START_POINT+subTotalFileSize) / (1024.0*1024.0));
    fprintf(indexFile, "== GENOME ==\n");
    fprintf(indexFile, "Sequences    : %u\n", pGenome->getNChromosomes());
    fprintf(indexFile, "Size         : %u (%.1f M)\n", pGenome->getGenomeSize(), (float) pGenome->getGenomeSize() / (1000.0*1000.0));
    fprintf(indexFile, "Keys (GenPos): %u (%.1f M)\n", pGenome->getGenomeValidPositions(), (float) pGenome->getGenomeValidPositions() / (1000.0*1000.0));
    fprintf(indexFile, "    By Guider: %u (%.1f M)\n", (pGenome->getGenomeValidPositions()-pGenome->getGenomeValidPositionsOverLen()), (float) (pGenome->getGenomeValidPositions()-pGenome->getGenomeValidPositionsOverLen()) / (1000.0*1000.0));
    fprintf(indexFile, "   By OverLen: %u (%.1f M)\n", pGenome->getGenomeValidPositionsOverLen(), (float) pGenome->getGenomeValidPositionsOverLen() / (1000.0*1000.0));
    fprintf(indexFile, "N's          : %u (%.1f M)\n", pGenome->nTotalPositionsWithN, (float) pGenome->nTotalPositionsWithN / (1000.0*1000.0));
    fprintf(indexFile, "== COMMAND LINE USED ==\n");
    for (k=0; k < argc; k++) 
        fprintf(indexFile, "%s ", argv[k]);
    fprintf(indexFile, "\n");
    fprintf(indexFile, "== KEYS ==\n");
    fprintf(indexFile, "Distribution of genome positions per key-size:\n");
    fprintf(indexFile, "Up to size\tCount\n");
    for (k=0; k < HKEYS; k++) 
        fprintf(indexFile, "%u\t%u\n", cuts[k], counts[k]);
    fprintf(indexFile, "Maximum genome positions in a key:%llu, Index:%u\n",max,mx);
    fprintf(indexFile, "== END OF INDEX INFO ==\n");
    fclose(indexFile);
    //fprintf(stderr, "S0 Index File [%s] Size = %ld\n", indexFileName, get_file_size(indexFileName));

    if (pLargestSeq != NULL) free(pLargestSeq);
    
    fprintf(stderr, "==== Index process finished ====\n");
        
}

void ogIndex::setMemoryMode(char lowMemMode) {
    lowMemoryMode = lowMemMode;
}


char ogIndex::load(char *pSourceFileName, uint16_t *userKeySize) {
    char xfile[MAX_FILENAME];
    FILE *pfile;
    
    //copyAndRemoveClassicExtensions(pSourceFileName, purifiedSourceFileName, MAX_FILENAME, pEncoding->getShortExtensionName());
    strncpy(indexFileName, pSourceFileName, 2048);

    setPureFileName(xfile, pSourceFileName, "ogg", MAX_FILENAME);
    //printf_FileOperation("Reading Index Header");
    
    //pfile = fopen(xfile, "r");
//    char c;
//    for (c=1; c < 10; c++) {
//        pfile = openIndexFileAtSection(c);
//        fclose(pfile);        
//    }
    pfile = openIndexFileAtSection(6);
    fread(&header, sizeof(header), 1, pfile);
    fclose(pfile);
    //printf_FileOperationDone();
    if (header.userKeySize != *userKeySize) {
        //fprintf(stderr, "Header file reports %u as key size instead of %u. Check your -k argument.\n", header.userKeySize, userKeySize);
        fprintf(stderr, "Header file reports %u as key size. Reseting key size.\n", header.userKeySize);
        *userKeySize = header.userKeySize;
        //return 0;
    }
    uint32_t msd = (uint32_t) pow(10, header.meanKeys_log10 + header.sdKeys_log10*2);
    if (msd > header.nKeyCountsAtQ95) {
        nKeysLow = header.nKeyCountsAtQ95;
        nKeysHigh = msd;
    } else {
        nKeysLow = msd;
        nKeysHigh = header.nKeyCountsAtQ95;        
    }
    fprintf(stderr, ">> Keys: Low=%u, High=%u >> ", nKeysLow, nKeysHigh);
    
    pEncoding = getKeyEncoding(header.encodingName);
    //fprintf(stderr, ">> Header encoding [%s] version [%f]\n", header.encodingName, header.encodingVersion);
    if (pEncoding == NULL) {
        fprintf(stderr, "NO KEY ENCODING SPECIFIED OR NONE REGISTERED.\n");
        return 0;
    }
    pEncoding->setSizeInChars(*userKeySize);
    pfile = openIndexFileAtSection(2);
    pEncoding->load(pfile);
    fclose(pfile);
    if (pEncoding->getSizeInBits() < 1) {
        fprintf(stderr, "**** ERROR: Encoding reports 0 bits. ****\n");
        return 0;
    }
    fprintf(stderr, "Encoding [%sv%1.1f], KeySize [%u chars] >> ", pEncoding->getName(), pEncoding->version, pEncoding->getSizeInChars());
    
    pGuider = getGuider(header.guiderName);
    //fprintf(stderr, ">> Header guider [%s]\n", header.guiderName);
    if (pGuider == NULL) {
        fprintf(stderr, "NO GUIDER SPECIFIED OR NONE REGISTERED.\n");
        return 0;
    }
    fprintf(stderr, "Guider [%s]\n", pGuider->getName());
    pfile = openIndexFileAtSection(1);
    pGuider->load(pfile);
    fclose(pfile);
    pGuider->setKeySizeInChars(pEncoding->getSizeInChars());
    //fprintf(stderr, ">> Guider loaded [%s]\n", pGuider->getName());

    pGenome = new ogGenome(pSourceFileName, *userKeySize);
    pfile = openIndexFileAtSection(3);
    pGenome->loadChromosomesInfo(pfile);
    fclose(pfile);
    pfile = openIndexFileAtSection(4);
    pGenome->loadPackedGenome(pfile);
    fclose(pfile);
    
    pGenomePositions = new ogGenomePositions(pGenome->getGenomeValidPositions(), *userKeySize);
    pfile = openIndexFileAtSection(7);
    pGenomePositions->load(pfile);
    fclose(pfile);
    
    pKeys = new ogKeys(pEncoding->getSizeInBits(), *userKeySize, pGenomePositions, pEncoding->getTotalKeys());
    pfile = openIndexFileAtSection(5);
    pKeys->load(pfile);
    fclose(pfile);
        
    pMapParams->pKeys = pKeys;
    pMapParams->pGenome = pGenome;
    pMapParams->pGenPos = pGenomePositions;
    pMapParams->pKeyEncoding = pEncoding;
    pMapParams->pGuiderDontUseIt = pGuider;
    pMapParams->keySize = pEncoding->getSizeInChars();
    pMapParams->forceUpperCase = forceUpperCase;
        
    return 1;
}

void ogIndex::setSchedule(char *pSchedule) {
    char *pS = pSchedule;
    int i;
    for (i=0; *pS; pS++, i++) {
        //fprintf(stderr, "Setting Schedule for key [%c] function address %p\n", *pS, protoMapFuncTable[*pS]);
        scheduler.funcSchedule[i] = protoMapFuncTable[*pS]; // (uint16_t (*)(ogReadsMapper *))
        if (scheduler.funcSchedule[i] == NULL) {
            fprintf(stderr, "ERROR: SCHEDULE FUNCTION FOR <%c> in [%s] NOT DEFINED. ABORTING.\n", *pS, pSchedule);
            exit(1);
        }
        scheduler.funcKey[i] = *pS;
        scheduler.nMatches[i] = 0;
        scheduler.nCalls[i] = 0;
        scheduler.nFoundFwd[i] = 0;
        scheduler.nFoundRev[i] = 0;
    }
    scheduler.nScheduleFunc = i;
}

void ogIndex::threadingReads(uint16_t iThread) {
    //fprintf(stderr, "[Th %u rdy] ", iThread+1); fflush(stderr);
    ogReadsMapper *rdM = rdMapr[iThread];
    char newEvent = 0;
    while (rdM->finish == 0 || !rdM->isQueueEmpty()) {
        if (rdM->isQueueEmpty()) {
            rdM->waste(newEvent);
            newEvent = 0;
        } else {
            rdM->process1Read();
            newEvent = 1;
        }
        /*
        //fprintf(stderr, "W%u", iThread);  fflush(stderr);
        if (rdM->canWork) {
            //fprintf(stderr, "M%u", iThread); fflush(stderr);
            rdM->canWork = 0;
            rdM->map1Read();
        } else {
            rdM->waste();
        }
         */
    }
}

void ogIndex::openOutputSAM(char *pFileName, char *fname1, char *fname2) {
    closeOutputSAM();
    pSamSam = new ogSamWriter();
    fprintf(stderr, ">>>> Creating Output File : %s <<<<\n", pFileName);
    pSamSam->openFile(pFileName); // "samout.sam"
    pSamSam->writeHeaderChromosomes(pGenome);
    pSamSam->writeHeaderProgram(pMapParams->pGuiderDontUseIt, pMapParams->pKeyEncoding, pMapParams->cmdLine, fname1, fname2, (*readGroupHeader == 0 ? NULL : readGroupHeader));
}

void ogIndex::closeOutputSAM() {
    if (pSamSam != NULL) pSamSam->closeFile();
    pSamSam = NULL;
}


#define BUF_SIZE        100000

int charPosInStr(char *p, char c) {
    char *p2 = strchr(p, c);
    if (p2 == NULL) return -1;
    return p2 - p;
}
typedef struct exonsStr {
    std::map<std::string, Exon *> exons;
    std::map<std::string, Exon *>::iterator iter;
} ExonsStr;

typedef struct transcriptsStr {
    std::map<std::string, ExonsStr *> transcripts;
    std::map<std::string, ExonsStr *>::iterator iter;
} TranscriptsStr;

#define MAX_COMMENT 500

typedef struct Gene {
    uint32_t    chr;
    uint32_t    start;
    uint32_t    end;
    char        strand;
    char        comment[MAX_COMMENT];
} Gene;

typedef struct genesStr {
    std::map<std::string, TranscriptsStr *> geneTranscripts;
    std::map<std::string, TranscriptsStr *>::iterator transcIter;
    std::map<std::string, ExonsStr *> geneExons;
    std::map<std::string, ExonsStr *>::iterator exonIter;
    std::map<std::string, Gene *> genes;
    std::map<std::string, Gene *>::iterator geneIter;
} GenesStr;

void ogIndex::indexFromGTF(char *pGtfFileName, char *pGenomeFileName, char *outFileName) {
    
    gzFile      gtfFile, gnmFile, outFile;
    char        destFileName[MAX_FILENAME];
    char        buf[BUF_SIZE];
    char        *line = (char *) malloc(10000000);
    char        *pLine;
    char        lastChr[100];
    char        field[1000], field2[1000];
    int         inBuf = 0, newLinePos, tabPos, stPos;
    std::string exonKey, transcriptKey, geneKey;
    kseq_t     *seq;
    uint32_t    nChromosomes, genomeSize;
    uint32_t    nChrGTF, nExons, nSumExons, nAllExons = 0, nAllTranscripts = 0, nTranscripts;
    uint32_t    nAllGenes = 0;
    FILE        *fileExons;
    FILE        *fileTranscripts;

    gnmFile = gzopen(pGenomeFileName, "r");
    seq = kseq_init(gnmFile);
    
    copyAndRemoveClassicExtensions(outFileName == NULL ? pGenomeFileName : outFileName, destFileName, MAX_FILENAME, "-GTF", ".fq.gz");
    outFile = gzopen(destFileName,"w");
    
    for (genomeSize=nChromosomes=0; kseq_read(seq) >= 0; genomeSize += seq->seq.l, nChromosomes++) {
        fprintf(stderr, "%3d\t%10ld\t%12ld\t%s\t%s\n", nChromosomes+1, seq->seq.l, genomeSize + seq->seq.l, seq->name.s, seq->comment.s);
        gtfFile = gzopen(pGtfFileName, "r");
        buf[0] = 0;
        inBuf = 0;
        nChrGTF = 0;
        lastChr[0] = 0;
        nExons = 0;    // exons records reads
        nSumExons = 0; // different exons
        nTranscripts = 0;
        
        GenesStr        allGenes;
        ExonsStr        allExons;
        TranscriptsStr  allTranscripts;
        GenesStr        *gs;
        TranscriptsStr  *ts;
        ExonsStr        *es;
        Exon            *xon;
        Gene            *gn;
        uint32_t        startPos, endPos;
        char            strand;
        char            *gid, *transcid;
        
        while (!gzeof(gtfFile)) {
            // Read a line into line
            newLinePos = charPosInStr(buf, '\n');
            if (newLinePos < 0 || newLinePos >= inBuf) {
                //fprintf(stderr, "Reading GZ %d\n", inBuf); fflush(stderr);
                inBuf += gzread(gtfFile, buf+inBuf, BUF_SIZE-inBuf-1);
                buf[BUF_SIZE-1] = 0;
                newLinePos = charPosInStr(buf, '\n');
                if (newLinePos < 0) {
                    fprintf(stderr, "**** PROBLEMS HOUSTON, NO LF ****\n");
                }
            }
            strncpy(line, buf, newLinePos);
            line[newLinePos] = 0;
            //fprintf(stderr, "%s, %d, %d\n", line, newLinePos, inBuf);
            memmove(buf, buf+newLinePos+1, inBuf - newLinePos - 1);
            inBuf -= (newLinePos + 1);
            //fprintf(stderr, "%d ", inBuf); fflush(stderr);
            //
            
            if (line[0] != '#') {
                tabPos = charPosInStr(line, '\t');
                if (tabPos > 0) {
                    strncpy(field, line, tabPos); // chromosome ID
                    field[tabPos] = 0;
                    if (strcmp(field, lastChr) != 0) {
                        //strlcpy(lastChr, field, 100);
                        strncpy(lastChr, field, 100);
                        nChrGTF++;
                    }
                    if (nChrGTF == (nChromosomes+1)) {
                        // this is a line of the same chromosome
                        pLine = line + tabPos + 1;
                        tabPos = charPosInStr(pLine, '\t'); // source
                        if (tabPos >= 0) {
                            //fprintf(stderr, "pLine0=[%s]\n", pLine);
                            pLine += tabPos+1;
                            tabPos = charPosInStr(pLine, '\t'); // type
                            if (tabPos > 0) {
                                strncpy(field, pLine, tabPos); // type
                                field[tabPos] = 0;
                                if (strcmp(field,"exon") == 0) {
                                    // allGenes.genes[];
                                    pLine += tabPos + 1;
                            //fprintf(stderr, "pLine1=[%s]\n", pLine);
                                    // start position
                                    tabPos = charPosInStr(pLine, '\t');
                                    strncpy(field, pLine, tabPos);
                                    field[tabPos] = 0;
                                    sscanf(field, "%u", &startPos);
                                    pLine += tabPos + 1;
                                    // end position
                            //fprintf(stderr, "pLine2=[%s]\n", pLine);
                                    tabPos = charPosInStr(pLine, '\t');
                                    strncpy(field, pLine, tabPos);
                                    field[tabPos] = 0;
                                    sscanf(field, "%u", &endPos);
                                    pLine += tabPos + 1;
                                    // skip score
                            //fprintf(stderr, "pLine3=[%s]\n", pLine);
                                    pLine += charPosInStr(pLine, '\t') + 1; 
                            //fprintf(stderr, "pLine4=[%s]\n", pLine);
                                    // strand
                                    tabPos = charPosInStr(pLine, '\t');
                                    strncpy(field, pLine, tabPos);
                                    field[tabPos] = 0;
                                    sscanf(field, "%c", &strand);
                                    pLine += tabPos + 1;
                                    // skip codon frame
                            //fprintf(stderr, "pLine5=[%s]\n", pLine);
                                    pLine += charPosInStr(pLine, '\t') + 1; 
                            //fprintf(stderr, "pLine6=[%s]\n", pLine);
                                    // attributes
                                    // get gene_id 
                                    gid = strstr(pLine, "gene_id ");
                                    transcid = strstr(pLine, "transcript_id ");
                                    if (gid != NULL && transcid != NULL) {
                                        gid+=8;
                                        transcid+=14;
                                        tabPos = charPosInStr(gid, ';');                                            
                                        if (*gid == '"' || *gid == '\'') {
                                            stPos = charPosInStr(gid+1, *gid);
                                            gid++;
                                            if (stPos > 0 && stPos < tabPos) tabPos = stPos; else tabPos--;
                                        }
                                        gid[tabPos] = 0;
                                        tabPos = charPosInStr(transcid, ';');
                                        if (*transcid == '"' || *transcid == '\'') {
                                            stPos = charPosInStr(transcid+1, *transcid);
                                            transcid++;
                                            if (stPos > 0 && stPos < tabPos) tabPos = stPos; else tabPos--;
                                        }
                                        transcid[tabPos] = 0;
                                        
                                        //fprintf(stderr, "Gene Id = [%s], Transcript Id = [%s]\n", gid, transcid);
                                        
                                        snprintf(field2, 100, "%09u:%010u-%010u:%c", nChromosomes+1, startPos, endPos, strand);
                                        //fprintf(stderr, "%s\n", field2);
                                        exonKey = field2;
                                        allExons.iter = allExons.exons.find(exonKey);
                                        if (allExons.iter == allExons.exons.end()) {
                                            xon = new Exon();
                                            xon->chr = nChromosomes+1;
                                            xon->start = startPos;
                                            xon->end = endPos;
                                            xon->strand = strand;
                                            xon->iExon = nSumExons++;
                                            allExons.exons[exonKey] = xon;
                                            //strncpy(xon->geneName, gid, 100);
                                            //fprintf(stderr, "%s, nExones=%lu\n", field2, allExons.exons.size());
                                        } else {
                                            xon = allExons.iter->second;
                                        }
                                        
                                        transcriptKey = transcid;
                                        allTranscripts.iter = allTranscripts.transcripts.find(transcriptKey);
                                        if (allTranscripts.iter == allTranscripts.transcripts.end()) {
                                            es = new ExonsStr();
                                            allTranscripts.transcripts[transcriptKey] = es;
                                            nTranscripts++;
                                            //fprintf(stderr, "%s, nTranscripts=%lu\n", transcid, allTranscripts.transcripts.size());
                                        } else {
                                            es = allTranscripts.iter->second;
                                        }
                                        es->exons[exonKey] = xon;
                                        
                                        //sprintf(field2, "%s\tChr=%u (%s), ", gid, nChromosomes + 1, lastChr);
                                        snprintf(field2, 1000, "%s", gid);
                                        geneKey = field2;
                                        allGenes.transcIter = allGenes.geneTranscripts.find(geneKey);
                                        if (allGenes.transcIter == allGenes.geneTranscripts.end()) {
                                            ts = new TranscriptsStr();
                                            allGenes.geneTranscripts[geneKey] = ts;
                                            //fprintf(stderr, "%s, nGeneTranscripts=%lu\n", gid, allGenes.geneTranscripts.size());
                                        } else {
                                            ts = allGenes.transcIter->second;
                                        }
                                        ts->iter = ts->transcripts.find(transcriptKey);
                                        if (ts->iter == ts->transcripts.end()) {
                                            ts->transcripts[transcriptKey] = allTranscripts.transcripts[transcriptKey];
                                        }
                                        
                                        allGenes.exonIter = allGenes.geneExons.find(geneKey);
                                        if (allGenes.exonIter == allGenes.geneExons.end()) {
                                            es = new ExonsStr();
                                            allGenes.geneExons[geneKey] = es;
                                            //fprintf(stderr, "%s, nGeneExons=%lu\n", gid, allGenes.geneExons.size());
                                        } else {
                                            es = allGenes.exonIter->second;
                                        }
                                        es->exons[exonKey] = xon;
                                        
                                        allGenes.geneIter = allGenes.genes.find(geneKey);
                                        if (allGenes.geneIter == allGenes.genes.end()) {
                                            gn = new Gene();
                                            gn->chr = nChromosomes;
                                            gn->strand = strand;
                                            gn->start = startPos;
                                            gn->end = endPos;
                                            snprintf(gn->comment, MAX_COMMENT, "Chr=%u|%s", nChromosomes+1, lastChr);
                                            allGenes.genes[geneKey] = gn;
                                        } else {
                                            gn = allGenes.genes[geneKey];
                                            if (startPos < gn->start) gn->start = startPos;
                                            if (endPos > gn->end) gn->end = endPos;
                                        }
                                        
                                    }
                                    nExons++;
                                }
                            }
                        }
                    }
                }
            }
        }
        fprintf(stderr, "Chr %d, %u exons records, %u non identical exons in total\n", nChromosomes+1, nExons, nSumExons);
        gzclose(gtfFile);
        uint32_t nGenes = 0, nSmallExons = 0, nTotalExons = 0, nOverlappedExons = 0, chrSize = 0;
        for (allGenes.exonIter = allGenes.geneExons.begin(); allGenes.exonIter != allGenes.geneExons.end(); allGenes.exonIter++) {
            geneKey = allGenes.exonIter->first;
            gn = allGenes.genes[geneKey];
            //strncpy(field2, allGenes.exonIter->first.c_str(), 500);
            //fprintf(stderr, "%s {", allGenes.exonIter->first.c_str());
            //for (allTranscripts.iter = allGenes.geneTranscripts[allGenes.exonIter->first]->transcripts.begin(); allTranscripts.iter != allGenes.geneTranscripts[allGenes.exonIter->first]->transcripts.end(); allTranscripts.iter++) {
            //    fprintf(stderr, " %s", allTranscripts.iter->first.c_str());
            //}
            //fprintf(stderr, " }, nExons=%lu\n",allGenes.geneExons[allGenes.exonIter->first]->exons.size());
            //fprintf(stderr, " }, nTranscripts=%lu, nExons=%lu",allGenes.geneTranscripts[allGenes.exonIter->first]->transcripts.size(), allGenes.geneExons[allGenes.exonIter->first]->exons.size());
            es = allGenes.geneExons[geneKey]; // allGenes.exonIter->first
            ts = allGenes.geneTranscripts[geneKey]; // allGenes.exonIter->first
            gzprintf(outFile, ">%s\t%s|%u|%u|%u|%c|%u ex|%u tr\n", geneKey.c_str(), gn->comment, gn->chr+1, gn->start, gn->end, gn->strand, es->exons.size(), ts->transcripts.size());
            uint32_t ntSize = 0, exonSize;
            uint32_t endAnt=0, startAnt=0;
            uint32_t genePos = 0, nExon=0;
            for (es->iter = es->exons.begin(); es->iter != es->exons.end(); es->iter++) {
                if (genePos > 0) { genePos +=3; gzprintf(outFile, "NNN\n"); } // --- ? 
                xon = es->iter->second;
                //if (strcmp(xon->geneName, field2) != 0) {
                //    fprintf(stderr, "Exon compartido, Gene actual [%s], Exon Gene [%s]\n", field2, xon->geneName);
                //      Esto demuestra que si hay exones comapartidos
                //} 
                ntSize += (exonSize = xon->end - xon->start + 1);
                xon->phyStart = genePos;
                strncpy(line, seq->seq.s + xon->start-1, exonSize); // start is 1-based
                line[exonSize] = 0;
                if (gn->strand == '-') {
                    // Reverse complement
                    uint32_t i, j = exonSize-1, k=exonSize >> 1, ch, rc;
                    
                    // reverse...
                    for (i=0; i < k; i++, j--) {
                        ch = line[i];
                        line[i] = line[j];
                        line[j] = ch;
                    }
                    // complement...
                    for (i=0; i < exonSize; i++) {
                        ch = line[i];
                        switch(ch) {
                            case 'A': rc = 'T'; break;
                            case 'a': rc = 't'; break;
                            case 'C': rc = 'G'; break;
                            case 'c': rc = 'g'; break;
                            case 'G': rc = 'C'; break;
                            case 'g': rc = 'c'; break;
                            case 'T': rc = 'A'; break;
                            case 't': rc = 'a'; break;
                            default: rc = ch; break;
                        }
                        line[i] = rc;
                    }
                }
                for (uint32_t i=0; i*80 < exonSize; i++) {
                    strncpy(field, line+i*80, 80);
                    field[80] = 0;
                    gzprintf(outFile, "%s\n",field);
                }
                genePos += exonSize;
                xon->phyEnd = genePos-1; // los limites son incluídos y que sea consistente con la distancia original
                nExon++;
                //printf("// Exon=%d, Key=%s, Start=%u, End=%u, phyStart=%u, phyEnd=%u, iExon=%u\n", nExon, es->iter->first.c_str(),xon->start, xon->end, xon->phyStart, xon->phyEnd, xon->iExon);
                if (exonSize < 48) {
                    //fprintf(stderr, "Exon pequeño [%u nt] gen [%s], de %u a %u\n",exonSize, allGenes.exonIter->first.c_str(),xon->start, xon->end); fflush(stderr);
                    nSmallExons++;
                }
                if (xon->start < endAnt) {
                    //fprintf(stderr, "Inicio=%d antes que el endAnt=%d para el gen [%s]\n", xon->start, endAnt, allGenes.exonIter->first.c_str()); fflush(stderr);
                    //printf("////////////////////////////// OVERLAP\n");
                    nOverlappedExons++;
                    xon->overlapsPrev = 1;
                } else {
                    xon->overlapsPrev = 0;
                }
                endAnt = xon->end;
                startAnt = xon->start;
                nTotalExons++;
            }
            //fprintf(stderr, ", ntSize=%u\n", ntSize);
            chrSize += ntSize;
            nGenes++;
        }
        fprintf(stderr, "ChrExonSize=%u, Genes=%u, Different Exons=%u, Small Exons=%u, Overlapped Exons=%u\n", chrSize, nGenes, nTotalExons, nSmallExons, nOverlappedExons);
        
        // WRITE EXONS
        copyAndRemoveClassicExtensions(outFileName == NULL ? pGtfFileName : outFileName , field, MAX_FILENAME, ".exons", ".ogx");
        fileExons = fopen(field, (nChromosomes == 0 ? "w" : "a"));
        if (nChromosomes == 0) {
            fwrite(&nChromosomes, sizeof(uint32_t), 1, fileExons);   // first 4 bytes chromosomes (genes)
            fwrite(&nAllExons, sizeof(uint32_t), 1, fileExons);         // next  4 bytes should be all exons, add that space.
        }
        Exon *pExon = (Exon *) malloc(nTotalExons * sizeof(Exon));
        nTotalExons = 0;
        nGenes = 0;
        uint32_t exonSize;
        uint32_t genePos;
        for (allGenes.exonIter = allGenes.geneExons.begin(); allGenes.exonIter != allGenes.geneExons.end(); allGenes.exonIter++) {
            es = allGenes.geneExons[allGenes.exonIter->first];
            genePos = 0;
            for (es->iter = es->exons.begin(); es->iter != es->exons.end(); es->iter++) {
                if (genePos > 0) genePos +=3;
                xon = es->iter->second;
                exonSize = xon->end - xon->start + 1;
                xon->phyStart = genePos;
                genePos += exonSize;
                xon->phyEnd = genePos-1;
                xon->iExon = nAllExons + nTotalExons; // Redefine iExon ??? posibles problemas porque se necesita fijo en transcripts
                xon->iGene = nAllGenes + nGenes;
                pExon[nTotalExons++] = *xon;
            }
            nGenes++;
        }
        nAllGenes += nGenes;
        fwrite(pExon, sizeof(Exon), nTotalExons, fileExons);
        fclose(fileExons);
        free(pExon);
        
        // WRITE TRANSCRIPTS
        copyAndRemoveClassicExtensions(outFileName == NULL ? pGtfFileName : outFileName, field2, MAX_FILENAME, ".transcripts", ".ogx");
        fileTranscripts = fopen(field2, (nChromosomes == 0 ? "w" : "a"));
        if (nChromosomes == 0) {
            fwrite(&nChromosomes, sizeof(uint32_t), 1, fileTranscripts);           // first 4 bytes chromosomes (genes)
            fwrite(&nTranscripts, sizeof(uint32_t), 1, fileTranscripts);        // next 4 bytes transcripts
        }
        Transcript *pTranscript = (Transcript *) malloc(nTranscripts * sizeof(Transcript));
        Transcript *pT;
        nTranscripts = 0;
        for (allGenes.transcIter = allGenes.geneTranscripts.begin(); allGenes.transcIter != allGenes.geneTranscripts.end(); allGenes.transcIter++) {
            ts = allGenes.geneTranscripts[allGenes.transcIter->first];
            for (ts->iter = ts->transcripts.begin(); ts->iter != ts->transcripts.end(); ts->iter++) {
                pT = pTranscript + nTranscripts;
                pT->iTranscript = nTranscripts;
                pT->nExons = ts->iter->second->exons.size();
                strncpy(pT->transcriptID, ts->iter->first.c_str(), 31);
                pT->transcriptID[31] = 0;
                nTranscripts++;
            }
        }
        fwrite(pTranscript, sizeof(Transcript), nTranscripts, fileTranscripts);
        fclose(fileTranscripts);
        free(pTranscript);
                
        // FREE EXONS & TRANSCRIPTS
        nAllExons += nTotalExons;
        for (allExons.iter = allExons.exons.begin(); allExons.iter != allExons.exons.end(); allExons.iter++) {
            delete allExons.iter->second;
        }
        for (allTranscripts.iter = allTranscripts.transcripts.begin(); allTranscripts.iter != allTranscripts.transcripts.end(); allTranscripts.iter++) {
            delete allTranscripts.iter->second;
        }
        nAllTranscripts += nTranscripts;
    }
    fileExons= fopen(field, "r+");
    fwrite(&nAllGenes, sizeof(uint32_t), 1, fileExons); // all genes
    fwrite(&nAllExons, sizeof(uint32_t), 1, fileExons); // all exons
    fclose(fileExons);

    fileTranscripts= fopen(field2, "r+");
    fwrite(&nAllGenes, sizeof(uint32_t), 1, fileTranscripts);       // all genes
    fwrite(&nAllTranscripts, sizeof(uint32_t), 1, fileTranscripts); // all Transcripts
    fclose(fileTranscripts);

    gzclose(outFile);
    kseq_destroy((kseq_t *) seq);
    gzclose(gnmFile);
    
    free(line);
    
    // The next step is that the user must run: 
    // ogMapper index <file name generated here>.fq.gz <all parameters>
    // ogMapper count <file name generated here>.fq.gz <all parameters as in map>
    
}

void ogIndex::loadExonAndTranscriptFiles(CountingInfo *pCountingInfo) {
    FILE        *fileExons;
    FILE        *fileTranscripts;
    uint32_t    i, lastGene;

    copyAndRemoveClassicExtensions(pCountingInfo->OGXfile[0] == 0 ? pCountingInfo->GTFfile : pCountingInfo->OGXfile, pCountingInfo->destFileName, MAX_FILENAME, ".exons", ".ogx");

    printf_FileOperationFile("Loading Exon Info", pCountingInfo->destFileName);
    fileExons = fopen(pCountingInfo->destFileName,"r");
    fread(&pCountingInfo->nGenes, 1, sizeof(uint32_t), fileExons);         // all chromosomes=genes
    fread(&pCountingInfo->nExons, 1, sizeof(uint32_t), fileExons);         // all exons
    
    pCountingInfo->pExons = (Exon *) malloc(pCountingInfo->nExons * sizeof(Exon));
    pCountingInfo->pExonReads = (uint32_t *) calloc(pCountingInfo->nExons+1, sizeof(uint32_t));
    pCountingInfo->pGeneReads = (uint32_t *) calloc(pCountingInfo->nGenes+1, sizeof(uint32_t));
    pCountingInfo->pGeneExonStart = (uint32_t *) calloc(pCountingInfo->nGenes+1, sizeof(uint32_t));
    pCountingInfo->pGeneExon = (uint32_t *) calloc(10000, sizeof(uint32_t)); // max of 10000 exons per gene
    
    fread(pCountingInfo->pExons, sizeof(Exon), pCountingInfo->nExons, fileExons);
    fclose(fileExons);
    
    // Define the 1st Exon for each gene into pCountingInfo->pGeneExonStart
    lastGene = 0xFFFFFFFF;
    Exon        *pX = pCountingInfo->pExons;
    uint32_t    *pS = pCountingInfo->pGeneExonStart;
    for (i=0; i < pCountingInfo->nExons; i++, pX++) {
        if (pX->iGene != lastGene) {
            lastGene = pX->iGene;
            pS[lastGene] = i;
        }
        //if (pX->iGene == 0) {
        //    fprintf(stderr, "&&&&&&&& i=%u, Exon: chr=%u, start=%u, end=%u, strand=%c, phyStart=%u, phyEnd=%u, iExon=%u, iGene=%u\n", i, pX->chr, pX->start, pX->end, pX->strand, pX->phyStart, pX->phyEnd, pX->iExon, pX->iGene); fflush(stderr);
        //}
    }
    //
    
    printf_FileOperationDone();
    
    copyAndRemoveClassicExtensions(pCountingInfo->OGXfile[0] == 0 ? pCountingInfo->GTFfile : pCountingInfo->OGXfile, pCountingInfo->destFileName, MAX_FILENAME, ".transcripts", ".ogx");

    printf_FileOperationFile("Loading Transcripts Info", pCountingInfo->destFileName);
    fileTranscripts = fopen(pCountingInfo->destFileName,"r");
    fread(&pCountingInfo->nChromTrans, 1, sizeof(uint32_t), fileTranscripts);      // all "genes" chromosomes
    fread(&pCountingInfo->nTranscripts, 1, sizeof(uint32_t), fileTranscripts);      // all transcripts

    pCountingInfo->pTranscripts = (Transcript *) malloc(pCountingInfo->nTranscripts * sizeof(Transcript));

    fread(pCountingInfo->pTranscripts, sizeof(Transcript), pCountingInfo->nTranscripts, fileTranscripts);
    fclose(fileTranscripts);
    printf_FileOperationDone();
    
}

uint32_t ogIndex::findExonIndexForPseudoPosition(uint32_t xgene, uint32_t xpos) {
    uint32_t min = 0;
    uint32_t max = countInfo.nExons;
    uint32_t mid;
    Exon *exn;
    while (min < max) {
        mid = (max + min) >> 1;
        exn = countInfo.pExons + mid;
        //fprintf(stderr, "Input[Chr:%u,Pos:%u],Curr[Chr:%u,Pos:%u-%u],min:%d,max:%d\n",xchr,xpos,exn->chr,exn->phyStart,exn->phyEnd,min,max); fflush(stderr);
        if (xgene > exn->iGene) { min = mid+1; }
        else if (xgene < exn->iGene) { 
            max = mid - 1;
        } else {
            // Same chromosome, check position
            if (xpos > exn->phyEnd) { min = mid+1; }
            else if (xpos < exn->phyStart) { max = mid-1; }
            else { min = mid; break; }
        }
    }
    return min;
}

Exon* ogIndex::findExonForPseudoPosition(uint32_t xgene, uint32_t xpos) {
    return countInfo.pExons + findExonIndexForPseudoPosition(xgene, xpos);
}


void ogIndex::count(char *pGTFfileName, char *pSourceFileName1, char *pSourceFileName2) {

    FILE            *ofile;
    ogChromosome    *pChr;
    
    strncpy(countInfo.GTFfile, pGTFfileName, 1000);
    loadExonAndTranscriptFiles(&countInfo);
            
    pMapParams->pExonCounts = countInfo.pExonReads;
    pMapParams->pGeneCounts = countInfo.pGeneReads;
    pMapParams->pTranscriptCounts = NULL;
    pMapParams->nGenes = countInfo.nGenes;
    pMapParams->nExons = countInfo.nExons;
    //pMapParams->countExons = 1;       // defined by parameters ???
    //pMapParams->countTranscripts = 0; // defined by parameters ???
    pMapParams->pGeneExonStart = countInfo.pGeneExonStart;
    pMapParams->pGeneExon = countInfo.pGeneExon;

    char        *pSequence = NULL, *pS, *pL;
    uint32_t     maxSeqLen = pGenome->getLargestChromosome()->size+1, effLen, n, i, nEx;
    double       sumReadByLen = 0, readByLen;
    fprintf(stderr, "Allocating %.1f Mb (%lu bytes) of memory for largest chromosome.\n", ((float) maxSeqLen * sizeof(char))/(1024.0 * 1024.0), maxSeqLen * sizeof(char));
    pSequence = (char *) malloc(maxSeqLen * sizeof(char));
    Exon *pX = countInfo.pExons;
    n = 0;
    for (i=0; i < countInfo.nGenes; i++) {
        pChr = pGenome->getChromosome(i);
        // Estimate effective gene size
        memset(pSequence, 0, pChr->size);
        for(; pX->iGene < i && n < countInfo.nExons; pX++, n++);
        for(nEx=0; pX->iGene == i && n < countInfo.nExons; pX++, n++, nEx++) memset(pSequence + pX->phyStart, 1, pX->phyEnd - pX->phyStart + 1);
        for (effLen=0, pS = pSequence, pL = pS + pChr->size; pS < pL; pS++) effLen += *pS;
        pChr->effectiveSize = effLen;
        pChr->nExons = nEx;
    }
    
    fprintf(stderr, "Count Starting...\n");
    mapOrCount(pSourceFileName1, pSourceFileName2, 'C');
        
    uint32_t    geneCount = 0;
    uint64_t    counts = 0;
    uint64_t    exonCount = 0;
    double      tpm;

    
    fprintf(stderr, "\n");
    printf_FileOperationFile("Saving gene counts", genesCountFile);
    ofile = fopen(genesCountFile, "w");
    fprintf(ofile, "#Seq\tComment\tName\tEff Size\tReads\tTPM (pseudo)\n"); fflush(ofile);
    //for (i=0; i < countInfo.nGenes; i++) counts += countInfo.pGeneReads[i];
    for (i=0; i < countInfo.nGenes; i++) {
        pChr = pGenome->getChromosome(i);
        readByLen = ((double) countInfo.pGeneReads[i] / (double) pChr->effectiveSize);
        sumReadByLen += readByLen;
    }
    for (i=0; i < countInfo.nGenes; i++) {
        //if (countInfo.pGeneReads[i] > 0) {
        pChr = pGenome->getChromosome(i);
        
        readByLen = ((double) countInfo.pGeneReads[i] / (double) pChr->effectiveSize);
        tpm = (double) (readByLen * 1000000) / sumReadByLen;
        // Write
            fprintf(ofile, "%u\t%s\t%s\t%u\t%u\t%.4f\n", i+1, pChr->comment, pChr->name, pChr->effectiveSize, countInfo.pGeneReads[i], tpm);
            geneCount++;
            counts += countInfo.pGeneReads[i];
        //}
    }
    fclose(ofile);
    printf_FileOperationDone();
    
    if (pMapParams->countExons) {
        // Suma de Reads By Len en formula de TPM, pero, si los exones se traslapan, solo sumar el mayor
        double sbl = 0, sblant = 0;
        for (sumReadByLen=0, pX = countInfo.pExons, i=0; i < countInfo.nExons; i++, pX++) {
            sbl = (double) countInfo.pExonReads[i] / (double) (pX->end-pX->start+1);
            if (countInfo.pExons[i].overlapsPrev) {
                if (sbl > sblant) sblant = sbl;
            } else {
                sumReadByLen += sblant;
                sblant = sbl;
            }
        }
        sumReadByLen += sblant;
        
        printf_FileOperationFile("Saving exon counts", exonsCountFile);
        ofile = fopen(exonsCountFile, "w");
        fprintf(ofile, "Row\tSeq\tComment\tName\tExon\tCoordinate\tStrand\tSize\tReads\tTPM (pseudo)\tGeneReads\n");
        Exon *pX = countInfo.pExons;
        for (i=0; i < countInfo.nExons; i++, pX++) {
            //if (countInfo.pExonReads[i] > 0) {
                pChr = pGenome->getChromosome(pX->iGene);
                tpm =  ((double) countInfo.pExonReads[i] * 1000000 / (double) (pX->end-pX->start+1)) / sumReadByLen;
                fprintf(ofile, "%u\t%u\t%s\t%s\t%u\t%u:%u-%u\t(%c)\t%u\t%u\t%.4f\t%u\n", i+1, pX->iGene+1, pChr->comment, pChr->name, pX->iExon, pX->chr, pX->start, pX->end, pX->strand, pX->end-pX->start+1, countInfo.pExonReads[i], tpm, countInfo.pGeneReads[pX->iGene]);
            //}
                exonCount++;
        }
        fclose(ofile);
        printf_FileOperationDone();
    }

    if (pMapParams->countTranscripts) {
        printf_FileOperationFile("Saving transcript counts", transcriptsCountFile);
        ofile = fopen(transcriptsCountFile, "w");
        fprintf(ofile, "Row\tSeq\tComment\tName\tExon\tCoordinate\tStrand\tSize\tReads\n");
        fclose(ofile);
        printf_FileOperationDone();        
    }
    
    fprintf(stderr, "\nK=K=K=K=K=K=K=K %u Genes Summing %llu reads, Average=%llu/gene, Exon Counts=%llu\n", geneCount, counts, counts/(geneCount == 0 ? 1 : geneCount), exonCount); fflush(stderr);
    if (pSequence != NULL) free(pSequence);
}


void ogIndex::map(char *pSourceFileName1, char *pSourceFileName2) {
    mapOrCount(pSourceFileName1, pSourceFileName2, 'M');
}

long __get_file_size(char *filename) {
    FILE *fp = fopen(filename, "r");

    if (fp==NULL)
        return -1;

    if (fseek(fp, 0, SEEK_END) < 0) {
        fclose(fp);
        return -1;
    }

    long size = ftell(fp);
    // release the resources when not required
    fclose(fp);
    return size;
}

void ogIndex::mapOrCount(char *pSourceFileName1, char *pSourceFileName2, char mode) {
    uint64_t    r, rAssigned=0;
    kseq_t     *seq1, *seq2 = NULL;
    ogFastAQGZreader *pFAR1, *pFAR2 = NULL;
    ogFastAQsequence *pfSeq1, *pfSeq2 = NULL;
    gzFile     fastaGZ_1, fastaGZ_2;
    char       isPaired = (pSourceFileName2 != NULL);
    uint16_t    iTh;
    uint64_t    fs1, fs2;
    
    fs1 = __get_file_size(pSourceFileName1);
    fastaGZ_1 = gzopen(pSourceFileName1, "r");
    if (isPaired) {
        fs2 = __get_file_size(pSourceFileName1);
        fastaGZ_2 = gzopen(pSourceFileName2, "r");
    }
    
    FILE *pFU1 = NULL;
    FILE *pFU2 = NULL;
    if (unmappedMode) {
        pSamSam->unmapMode = unmappedMode;
        if (unmappedMode == 1 || unmappedMode == 2) {
            char *pF1 = (char *) calloc(MAX_FILENAME, sizeof(char));
            strncpy(pF1, pSourceFileName1, MAX_FILENAME);
            strncat(pF1, ".unmap", MAX_FILENAME);
            fprintf(stderr, "UNMAPPED FILENAME 1:%s\n",pF1); fflush(stderr);
            pFU1 = fopen(pF1, "w");
            if (pSourceFileName2 != NULL) {
                char *pF2 = (char *) calloc(MAX_FILENAME, sizeof(char));
                strncpy(pF2, pSourceFileName2, MAX_FILENAME);
                strncat(pF2, ".unmap", MAX_FILENAME);
                fprintf(stderr, "UNMAPPED FILENAME 2:%s\n",pF2);  fflush(stderr);
                pFU2 = fopen(pF2, "w");
                free(pF2);
                pSamSam->pUnmapFile2 = pFU2;
            }
            free(pF1);
            pSamSam->pUnmapFile1 = pFU1;
        } else {
            fprintf(stderr, "UNMAPPED FILES NOT GENERATED.\n");  fflush(stderr);
        }
    }
    
    // ReadGroup
    if (*readGroupHeader != 0) {
        // Read Group specified
        char *rgid = strstr(readGroupHeader, "ID:");
        if (rgid != NULL) {
            strncat(pSamSam->readGroup,"RG:Z:",MAX_SAM_FILENAME);
            char *dest = pSamSam->readGroup + 5;
            char isNum = 1;
            for(rgid+=3; (*rgid >= ' '); *dest++ = *rgid++) isNum = (isNum && *rgid >= '0' && *rgid <= '9');
            *dest++ = '\t';
            *dest++ = 0;
            if (isNum) *(pSamSam->readGroup+3) = 'i'; // RG:i
        }
    }
    
    ogAligner *pAlign; // = new ogAligner(pGenome);
        
    //ogReadsMapper rdMapr(pMapParams, maxReadLen);
    rdMapr = (ogReadsMapper **) malloc((nThreads+1) * sizeof(ogReadsMapper *));
    
    //pSamSam->writeHeaderProgram(pMapParams->pGuiderDontUseIt, pMapParams->pKeyEncoding, pMapParams->cmdLine, pSourceFileName1, pSourceFileName2);

    
    if(kseqReadingMode) {
        seq1 = kseq_init(fastaGZ_1);
        if (isPaired) seq2 = kseq_init(fastaGZ_2);
    } else {
        pFAR1 = new ogFastAQGZreader();
        pFAR1->setFileReader(&fastaGZ_1, 0, nThreads);
        pFAR1->setFileSize(fs1);
        if (isPaired) {
            pFAR2 = new ogFastAQGZreader();
            pFAR2->setFileReader(&fastaGZ_2, 1, nThreads);
            pFAR2->setFileSize(fs2);
        }
    }
    
    for (iTh=0; iTh <= nThreads; iTh++) {
        
        ogSamWriter *pSam = pSamSam;
        if (filePerThread) {
            if (iTh == 0) {
                pSamSam->tempId = 0;
            } else {
                pSam = new ogSamWriter();
                pSam->assignTemporalFile(pSamSam->outFileName, iTh, pSamSam->isBAM);
                if (unmappedMode) {
                    pSam->pUnmapFile1 = pSamSam->pUnmapFile1; // unmapped goes to same file
                    pSam->pUnmapFile2 = pSamSam->pUnmapFile2;
                    pSam->unmapMode = unmappedMode;
                    strncat(pSam->readGroup,pSamSam->readGroup,MAX_SAM_FILENAME);
                }
            }
        }
        rdMapr[iTh] = new ogReadsMapper(pMapParams, &scheduler, new ogAligner(pGenome), pSam, &countInfo, maxReadLen, iTh, header.meanKeys_log10, header.sdKeys_log10, header.nKeyCountsAtQ95, productionMode, isPaired);
        //rdMapr[iTh]->canWork = 0;
        if (nQueue > 0) rdMapr[iTh]->setMaxQsize(nQueue);
        rdMapr[iTh]->setMode(mode);
    }
    
    uint32_t theKey = 1000000;
    //pKeys->printKeyInfo(theKey);
    //pGenomePositions->printGenomePositionsInfo(pKeys->getInfoForKey(theKey)->offset, pKeys->getInfoForKey(theKey)->size, pGenome);
    
    //clock_t start_t = clock(), prev_t = clock(), curr_t, reading_t = 0;

    uint64_t wastingCycles = 0;
    uint64_t wastingCycles2 = 0;
    uint16_t kTh;
    //CircularArray<ogSingleRead *> *caReads = new CircularArray<ogSingleRead *>(maxReadsBuffer, maxReadsBuffer >> 4);
    ogSingleRead *read1;
    ogSingleRead *read2;
    //uint32_t    nReadObj = 0;
    //for (nReadObj = 0; nReadObj < caReads->size(); nReadObj++) {
    //    caReads->push(new ogSingleRead());
    //}

    fprintf(stderr, "|1 Million Reads: .=25,000 %2u+1 threads|\n", nThreads); fflush(stderr); // mode == 'M' ? "mapping" : "counting"
    fprintf(stderr, "|______________________________________|Reads/s|Mreads|Elap t|Tot r/s|%%Done|Left t|");fflush(stderr); // Map %%|
    if (nThreads > 0) {
        hilos = (thread **) malloc(nThreads * sizeof(thread *));
        for (iTh=0; iTh < nThreads; iTh++) {
            //fprintf(stderr, "%u ", iTh); flush(std::cout);
            hilos[iTh] = new thread(&ogIndex::threadingReads, this, iTh);
        }        
    }
    
    uint64_t reading = 0;
    uint64_t preparing = 0;
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point block = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point prev = std::chrono::steady_clock::now(), curr;
    std::chrono::steady_clock::time_point end;
    iTh = 0;
    int mapToRoot = 0; // Number of reads to be mapped intoRoot
    //if (nThreads > 0) fprintf(stderr, "\n");
    for (r=0; r < maxReads; r++) {
        
        if (kseqReadingMode) {
            if (kseq_read(seq1) < 0) break; // No more reads
            if (isPaired) { kseq_read(seq2); }
        } else {
            pfSeq1 = (ogFastAQsequence *) malloc(sizeof(ogFastAQsequence));
            if (pFAR1->getNextSequenceByType(pfSeq1) == 0) {
                free(pfSeq1);
                break;
            }
            if (isPaired) {
                pfSeq2 = (ogFastAQsequence *) malloc(sizeof(ogFastAQsequence));
                pFAR2->getNextSequenceByType(pfSeq2);
            }
        }
        
        curr = std::chrono::steady_clock::now();
        reading += chrono::duration_cast<chrono::nanoseconds>(curr - prev).count();
        if (r >= startingRead) {
            if (doPrepare) {
                if (isPaired) {
                    if (kseqReadingMode) {
                        if (strcmp(seq1->name.s, seq2->name.s) != 0) {
                            fprintf(stderr, "//////////////////// Different sequences? 1:%s 2:%s\n", seq1->name.s, seq2->name.s);
                        }
                    } else {
                        if (strcmp(pfSeq1->name.s, pfSeq2->name.s) != 0) {
                            fprintf(stderr, "//////////////////// Different sequences? 1:%s 2:%s\n", pfSeq1->name.s, pfSeq2->name.s);
                            break;
                        }
                    }
                }
                preparing += chrono::duration_cast<chrono::nanoseconds>(std::chrono::steady_clock::now() - curr).count();
                rAssigned++;
                // fprintf(stderr, "Mapping read %llu...\n",r);fflush(stderr);
                if (doProcess) {
                    if (nThreads > 0 && mapToRoot <= 0) {
                        kTh = 0;
                        while (mapToRoot <= 0 && rdMapr[iTh]->isQueueFull(0)) {
                            if (++iTh >= nThreads) iTh = 0;
                            if (++kTh % nThreads == 0) {
                                if (++wastingCycles == 0) wastingCycles2++;
                                //if (kTh / nThreads > 1) {
                                    // DEFINITIVELY IS FASTER avoiding the condition
                                    kTh = 0;
                                    mapToRoot = 1; // how many in main thread ? default to 1
                                    //break;
                                //}
                            }
                        }
                        if (mapToRoot <= 0) {
                            // fprintf(stderr, "pushing to >%d< %s\n", iTh, seq1->name.s);fflush(stderr);
                            if (kseqReadingMode) {
                                read1 = new ogSingleRead();
                                read1->setRead(seq1->name.s, seq1->name.l, seq1->seq.s, seq1->seq.l, seq1->qual.s, seq1->qual.l, r, 1);
                                read1->packIfNeeded(pGenome);
                                if (isPaired) {
                                    ogSingleRead *read2 = new ogSingleRead();
                                    read2->setRead(seq2->name.s, seq2->name.l, seq2->seq.s, seq2->seq.l, seq2->qual.s, seq2->qual.l, r, 1);
                                    read2->packIfNeeded(pGenome);
                                    read1->setPairedRead(read2);
                                }
                                rdMapr[iTh]->pushReadToQueue(read1, r);
                            } else {
                                // Process to push ogFastAQRead
                                pfSeq1->pRelease = pFAR1->assignReadToThread(pfSeq1->pBuffer, iTh);
                                if (isPaired) pfSeq2->pRelease = pFAR2->assignReadToThread(pfSeq2->pBuffer, iTh);
                                rdMapr[iTh]->pushReadToQueue(pfSeq1, pfSeq2, r);
                            }
                            if (++iTh >= nThreads) iTh = 0;
                            // fprintf(stderr, "pushed\n");fflush(stderr);
                        }
                    }
                    if (mapToRoot > 0 || nThreads == 0) {
                        // This is the root thread
                        // fprintf(stderr, "[RootThread]");fflush(stderr);
                        //rdMapr[nThreads]->mapRead(read1);
                        if (kseqReadingMode) {
                            read1 = new ogSingleRead();
                            read1->setRead(seq1->name.s, seq1->name.l, seq1->seq.s, seq1->seq.l, seq1->qual.s, seq1->qual.l, r, 1);
                            read1->packIfNeeded(pGenome);
                            if (isPaired) {
                                ogSingleRead *read2 = new ogSingleRead();
                                read2->setRead(seq2->name.s, seq2->name.l, seq2->seq.s, seq2->seq.l, seq2->qual.s, seq2->qual.l, r, 1);
                                read2->packIfNeeded(pGenome);
                                read1->setPairedRead(read2);
                            }
                            rdMapr[nThreads]->processTheRead(read1);
                        } else {
                            pfSeq1->pRelease = pFAR1->assignReadToThread(pfSeq1->pBuffer, nThreads);
                            read1 = new ogSingleRead();
                            read1->setRead(pfSeq1->name.s, pfSeq1->name.l, pfSeq1->seq.s, pfSeq1->seq.l, pfSeq1->qual.s, pfSeq1->qual.l, r, 1);
                            read1->packIfNeeded(pGenome);
                            if (isPaired) {
                                pfSeq2->pRelease = pFAR2->assignReadToThread(pfSeq2->pBuffer, nThreads);
                                ogSingleRead *read2 = new ogSingleRead();
                                read2->setRead(pfSeq2->name.s, pfSeq2->name.l, pfSeq2->seq.s, pfSeq2->seq.l, pfSeq2->qual.s, pfSeq2->qual.l, r, 1);
                                read2->packIfNeeded(pGenome);
                                read1->setPairedRead(read2);
                                //pfSeq2->pBuffer->readsPending[nThreads]--;
                                (*pfSeq2->pRelease)++;
                            }
                            //pfSeq1->pBuffer->readsPending[nThreads]--;
                            (*pfSeq1->pRelease)++;
                            rdMapr[nThreads]->processTheRead(read1);
                            free(pfSeq1);
                            if (isPaired) free(pfSeq2);
                        }
                        //rdMapr[nThreads]->pushReadToQueue(seq->name.s, seq->name.l, seq->seq.s, seq->seq.l, seq->qual.s, seq->qual.l);
                        mapToRoot--;
                        delete read1;
                    }
                } else {
                    if (kseqReadingMode == 0) {
                        free(pfSeq1);
                        if (isPaired) free(pfSeq2);
                    }                    
                }
            } else {
                if (kseqReadingMode == 0) {
                    free(pfSeq1);
                    if (isPaired) free(pfSeq2);
                }
            }
            
            if (r % 25000 == 0) {
                if (r % 1000000 == 0) { 
                    end = std::chrono::steady_clock::now();
                    uint64_t curElapsed = chrono::duration_cast<chrono::milliseconds>(end - start).count();
                    uint64_t blockElapsed = chrono::duration_cast<chrono::milliseconds>(end - block).count();
                    if (r > 0) {
                        float fp = (float) curElapsed / 1000.0;
                        fprintf(stderr, "%7d|%6.0f|%3d:%02d|%7d|", 
                                (int) ((float) 1000000.0 / (float) (blockElapsed / 1000.0)), 
                                (double) r / 1000000.0, 
                                (int) trunc(fp / 60), (int) (fp - 60* trunc((fp / 60.0))), ///(curElapsed / 1000), 
                                (int) ((float) r / ((float) curElapsed / 1000.0)));
                        if (kseqReadingMode) {
                            fp = (float) gzoffset(fastaGZ_1) * 100 / (float) fs1;
                            fprintf(stderr, "%4.1f%%|", fp);
                        } else {
                            fp = pFAR1->estimatedFilePosition();
                            fprintf(stderr, "%4.1f%%|", fp);
                        }
                        fp = ((float) curElapsed / 10.0) / fp - (curElapsed / 1000);
                        //fprintf(stderr, "%6.0f|", fp);
                        fprintf(stderr, "%3d:%02d|", (int) trunc(fp / 60), (int) (fp - 60* trunc((fp / 60.0))));
                        //uint64_t nUnMap = 0;
                        //for (int ith=0; ith <= nThreads; ith++) nUnMap += rdMapr[ith]->nReadsUnmapped;
                        //fprintf(stderr, "%4.1f%%|", (float) (100.0- ((float) nUnMap*100.0)/(float) rAssigned));
                    }
                    fprintf(stderr, "\n"); 
                    block = end;
                }
                fprintf(stderr, "%c", (r % 250000 == 0 ? '|' : '.'));
                //fprintf(stderr, "%c",  (char) ( ( ((uint32_t) (r / 10000)) % 10) + 48));
                fflush(stderr);
            }
            //fprintf(stderr, "%u positions.\n", rdMapr.pRdKeyMap->pCandPosMan->getCount());
        }
        prev = std::chrono::steady_clock::now();
    }
    for (iTh=0; iTh <= nThreads; iTh++) {
        //rdMapr[iTh]->canWork = 0;
        rdMapr[iTh]->finish = 1;
        if (nThreads > 0 && iTh < nThreads) hilos[iTh]->join();
    }
    fprintf(stderr, "\n%llu Reads. %s Done.\n", r, mode == 'M' ? "Mapping" : "Counting"); fflush(stderr);

    end = std::chrono::steady_clock::now();
    
    //fprintf(stderr, "Cleaning read buffer (%u)...", caReads->used());
    //while (caReads->used() > 0) delete caReads->pull();    
    
    if (kseqReadingMode) {
        kseq_destroy(seq1);
        if (isPaired) kseq_destroy(seq2);
        //free(seq1);
        //if (isPaired) free(seq2);
    } else {
        delete pFAR1;
        if (isPaired) delete pFAR2;
    }

    gzclose(fastaGZ_1);
    if (isPaired) gzclose(fastaGZ_2);        
    
    //rdMapr.printSummary();

    uint64_t aligns = 0;
    uint64_t elapsedThreads = 0;
    uint64_t nReadsMappedR1R2 = 0;
    uint64_t nReadsMappedR2R1 = 0;
    uint64_t nReadsMappedR1 = 0;
    uint64_t nReadsMappedR2 = 0;
    uint64_t nReadsMappedTrans = 0;
    uint64_t nReadsMappedOther = 0;
    uint64_t nReadsUnmapped = 0;

    for (iTh=0; iTh <= nThreads; iTh++) {
        if (iTh < nThreads) fprintf(stderr, "[[[ Thread %u ]]]\n", iTh+1); else fprintf(stderr, "[[[ Thread ROOT ]]]\n");
        rdMapr[iTh]->printSummary(rAssigned);
        aligns += rdMapr[iTh]->pAligner->getFullAlignments();
        elapsedThreads += rdMapr[iTh]->elapsed + rdMapr[iTh]->preparation + rdMapr[iTh]->pairing;
        nReadsMappedR1R2 += rdMapr[iTh]->nReadsMappedR1R2;
        nReadsMappedR2R1 += rdMapr[iTh]->nReadsMappedR2R1;
        nReadsMappedR1 += rdMapr[iTh]->nReadsMappedR1;
        nReadsMappedR2 += rdMapr[iTh]->nReadsMappedR2;
        nReadsMappedTrans += rdMapr[iTh]->nReadsMappedTrans;
        nReadsMappedOther += rdMapr[iTh]->nReadsMappedOther;
        nReadsUnmapped += rdMapr[iTh]->nReadsUnmapped;
    }
    uint64_t totalElapsed = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    uint64_t rr =     nReadsMappedR1R2+nReadsMappedR2R1+nReadsMappedR1+nReadsMappedR2+nReadsMappedTrans+nReadsMappedOther; //+nReadsUnmapped
    if (rr > 0) {
        fprintf(stderr, "%llu Reads: R1R2=%llu, R2R1=%llu, R1+R2-=%llu, R1-R2+=%llu, Trans=%llu, Other=%llu, Unmapped=%llu\n",
        rAssigned,
        nReadsMappedR1R2,
        nReadsMappedR2R1,
        nReadsMappedR1,
        nReadsMappedR2,
        nReadsMappedTrans,
        nReadsMappedOther,
        nReadsUnmapped);
        fprintf(stderr, "%% Reads: R1R2=%.3f%%, R2R1=%.3f%%, R1+R2-=%.3f%%, R1-R2+=%.3f%%, Trans=%.3f%%, Other=%.3f%%, Unmapped=%.3f%%\n",
        (float) (nReadsMappedR1R2*100)/rAssigned,
        (float) (nReadsMappedR2R1*100)/rAssigned,
        (float) (nReadsMappedR1*100)/rAssigned,
        (float) (nReadsMappedR2*100)/rAssigned,
        (float) (nReadsMappedTrans*100)/rAssigned,
        (float) (nReadsMappedOther*100)/rAssigned,
        (float) (nReadsUnmapped*100)/rAssigned);
    }
    fprintf(stderr, "Total Reads Processed = %lld\n", rAssigned);
    fprintf(stderr, "Total Reads Mapped    = %lld (%.3f%%)\n", (rAssigned-nReadsUnmapped), (float) (rAssigned-nReadsUnmapped)*100/rAssigned);
    fprintf(stderr, "Total Reads Unmapped  = %lld (%.3f%%)\n", nReadsUnmapped, (float) nReadsUnmapped*100/rAssigned);
    fprintf(stderr, "Elapsed Time          = %.2f s\n", (float) totalElapsed / 1000);
    fprintf(stderr, "Reading Time          = %.2f s\n", (float) reading / (float) 1000000000);
    fprintf(stderr, "Preparing Time        = %.2f s\n", (float) preparing / (float) 1000000000);
    fprintf(stderr, "Map Thread Time (Avg) = %.2f s\n", (float) elapsedThreads / (float) ((nThreads+1) * 1000000));
    fprintf(stderr, "Queue and Other Times = %.2f s\n", (float) (totalElapsed-elapsedThreads/((nThreads+1)*1000.0)-reading/1000000.0-preparing/1000000.0) / (float) 1000);
    fprintf(stderr, "Wasting Cycles        = %lld * %lld\n", wastingCycles2+1, wastingCycles);
    fprintf(stderr, "Avg Mapping All       = %.2f reads/s\n", (float) rAssigned / ((float) chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000));
    fprintf(stderr, "Avg Mapping/Thread    = %.2f reads/s (%u threads)\n", (float) rAssigned / ((float) (nThreads+1) * (float) chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000), nThreads+1);
    fprintf(stderr, "Total Alignments      = %llu (%.3f%%)\n", aligns, ((float) aligns * 100) / rAssigned);

    for (iTh=0; iTh <= nThreads; iTh++) {
        delete rdMapr[iTh];
    }

    for (iTh=0; iTh < nThreads; iTh++) {
        delete hilos[iTh];
    }
    
    free(rdMapr);
    if (nThreads > 0) free(hilos);
    
    if (unmappedMode) {
        if (pSamSam->pUnmapFile2 != NULL) fclose(pSamSam->pUnmapFile2);
        if (pSamSam->pUnmapFile1 != NULL) fclose(pSamSam->pUnmapFile1);
    }
}

void ogIndex::registerMappingFunction(char key, uint16_t (*pFunc)(ogReadKeyMapping *)) {
    protoMapFuncTable[key] = pFunc;
}

// Functions for testing things
// AGCTCGCTAGAAAGCTAGCTAGCTAGCGACTACAGCTACGACTAAA
//           1                                2
// 0123456789012345678901234567890123456789012345
//           111111111122222222223333333333444444
//              L                             I
//              123456789012345678901234
void ogIndex::testRepes(char *pSourceFileName, uint16_t repes) {
    
    uint64_t    r, i, j, n, k, kk, o, l, m;
    char        x;
    kseq_t     *seq;
    gzFile     fastaGZ;
    const uint64_t  maxSpace = 24;
    uint8_t     repI;
    uint8_t     maxRepes = repes;
    uint8_t     minRepes = repes;
    fastaGZ = gzopen(pSourceFileName, "r");
    seq = kseq_init(fastaGZ);

    kk = 0;
    for (r=0; r < maxReads && kseq_read(seq) >= 0; r++) {
        //seq->seq.s, seq->seq.l;
        k = o = l = 0;
        for (i=0; i < seq->seq.l; i++) {
            x = seq->seq.s[i];
            for (j=i+1, n=1; j < seq->seq.l && x == seq->seq.s[j] && n < maxRepes; j++, n++);
            //for (repI=maxRepes; repI >= minRepes; repI--, n=0) {
            //    for (j=i+1, n=1; j < seq->seq.l && x == seq->seq.s[j] && n < repI; j++, n++);
            //    if (n >= repI) break;
            //}
            if (n >= minRepes) {
                if (i > l && i - l >= maxSpace) {
                    o += i - l - maxSpace + 1;
                }
                k++;
                l = i + n;
                i = j - 1;
            }
        }
        if (k == 0) {
            if (i - l >= maxSpace) {
                o += i - l - maxSpace + 1;
            }
        }
        fprintf(stderr, "Seq:%llu\tSize:%lu\tCount:%llu\tOverPassed:%llu\tSum:%llu\tRatio:%.3f\n",r+1,seq->seq.l,k,o,k+o,(float) k/(float) (k+o));
        if (k == 0) {
            fprintf(stderr, "#%s\n",seq->seq.s);
        }
        kk += k + o;
    }
    fprintf(stderr, "Total Count+Overpassed:%llu\n",kk);
    kseq_destroy(seq);
    gzclose(fastaGZ);
 }

void ogIndex::getChr20_22(char *pSourceFileName) {
    
    uint64_t    r;
    char        *p;
    kseq_t     *seq;
    gzFile     fastaGZ;
    const uint64_t  maxSpace = 24;
    uint32_t     k;
    fastaGZ = gzopen(pSourceFileName, "r");
    seq = kseq_init(fastaGZ);

    for (r=0; r < maxReads && kseq_read(seq) >= 0; r++) {
        //seq->seq.s, seq->seq.l;
        if (r >= 19 && r <= 21) {
            fprintf(stderr, ">%s %s\n", seq->name.s, seq->comment.s);
            p = seq->seq.s;
            for (k=1; *p; p++, k++) {
                fprintf(stderr, "%c", *p);
                if (k % 80 == 0) fprintf(stderr, "\n");
            }
            fprintf(stderr, "\n");
        }
    }
    kseq_destroy(seq);
    gzclose(fastaGZ);
 }

void ogIndex::setUpperCase(char upcase) {
    forceUpperCase = upcase;
}


void ogIndex::testStateMachine(char *pSourceFileName, char *pStateFileName, uint16_t maxLenNoResponse) {

    char        states[256][8]; // 8 bytes seems better
    char        ascii2col[256];
    uint8_t     currState;
    uint64_t    r, i, j, n, k, kk, o, oo, l, m, e, ee, delta, preI;
    char        x;
    kseq_t     *seq;
    gzFile     fastaGZ;
    
    fastaGZ = gzopen(pSourceFileName, "r");
    seq = kseq_init(fastaGZ);

    memset(states, 0, sizeof(states));
    
    fprintf(stderr, "Reading file %s\n", pStateFileName);
    
    FILE *pStates = fopen(pStateFileName, "r");
    char line[1000];
    
    for (i=0; fscanf(pStates, "%[^\n]s", line) != EOF; ) {
        if (line[0] != '#') {
            if (line[0] == 0) break;
            sscanf(line, "%hhu ", &currState);
            i = currState;
            sscanf(line, "%hhu %hhu %hhu %hhu %hhu %hhu %hhu %[^\n]s", &currState, &states[i][0], &states[i][1], &states[i][2], &states[i][3], &states[i][4], &states[i][5], (char *)ascii2col);
            fprintf(stderr, "St:%hhu\tA:%hhu\tC:%hhu\tG:%hhu\tT:%hhu\tN:%hhu\t*:%hhu\t%s\n", currState, states[i][0], states[i][1], states[i][2], states[i][3], states[i][4], states[i][5], ascii2col);
        } else {
            fprintf(stderr, "%s\n", line);
        }
        fscanf(pStates,"%c",&x);
        line[0] = 0;
    }
    
    fclose(pStates);

    memset(ascii2col, 5, sizeof(ascii2col));
    ascii2col['A'] = ascii2col['a'] = 0;
    ascii2col['C'] = ascii2col['c'] = 1;
    ascii2col['G'] = ascii2col['g'] = 2;
    ascii2col['T'] = ascii2col['t'] = 3;
    ascii2col['N'] = ascii2col['n'] = 4;

    kk = oo = ee = 0;
    for (r=0; r < maxReads && kseq_read(seq) >= 0; r++) {
        //seq->seq.s, seq->seq.l;
        k = o = e = preI = 0;
        delta = 0;
        for (i=0; i < seq->seq.l; ) {
            currState = 0;
            l = i + maxLenNoResponse;
            if (l >= seq->seq.l) l = seq->seq.l-1;
            for (j=i; j <= l; j++) {
                if ((currState = states[currState][ascii2col[seq->seq.s[j]]]) >= 250) break;
            }
            if (currState == 250) {
                i = j;
                k++;
                delta = 1;
                preI = i;
            } else if (currState == 251 && i - preI < maxLenNoResponse) {
                // failed, move to next nt and go on
                i++;
            } else {
                preI++;
                e+=delta;
                o++;
                i++;
                delta = 0;
            }
        }
        fprintf(stderr, "Seq:%llu\tSize:%lu\tCount:%llu\tOverPassed:%llu\tSum:%llu\tRatio:%.3f\tOvP Events:%llu\n",r+1,seq->seq.l,k,o,k+o,(float) k/(float) (k+o), e);
        if (k == 0) {
            fprintf(stderr, "#%s\n",seq->seq.s);
        }
        kk += k;
        oo += o;
        ee += e;
    }
    fprintf(stderr, "TOTALS... Count:%llu, Overpassed:%llu, Sum:%llu, Ratio:%.3f, OvP Events:%llu\n",kk,oo,kk+oo,(float) kk/(float) (kk+oo),ee);
    kseq_destroy(seq);
    gzclose(fastaGZ);
 }


void ogIndex::setProductionMode(char isProduction) {
    productionMode = isProduction;
}

uint32_t ogIndex::getLowKeyCountLimit() {
    return nKeysLow;
}

uint32_t ogIndex::getHighKeyCountLimit() {
    return nKeysHigh;
}

long ogIndex::get_file_size(char *filename) {
    FILE *fp = fopen(filename, "r");
    if (fp==NULL) return -1;
    if (fseek(fp, 0, SEEK_END) < 0) {
        fclose(fp);
        return -1;
    }
    long size = ftell(fp);
    // release the resources when not required
    fclose(fp);
    return size;
}

void ogIndex::checkKSEQ(char *filename) {
    gzFile zFile = gzopen(filename, "r");
    kseq_t     *seq;
    seq = kseq_init(zFile);
    uint64_t nr = 0;
    while (kseq_read(seq) >= 0  && nr < 10000000) {
        if (++nr % 100000 == 0) fprintf(stderr, "%llu%c", nr, (nr % 500000 == 0 ? '\n' : ' '));
        //nr++;
        fprintf(stderr, "Read=%llu; Lengths: Name=%lu [%s], Comment=%lu, Sequence=%lu, Quality=%lu\n", nr, seq->name.l, seq->name.s, seq->comment.l, seq->seq.l, seq->qual.l);
    }
    kseq_destroy((kseq_t *) seq);
    gzclose(zFile);
}

void ogIndex::setKSEQmode(char ks) {
    kseqReadingMode = ks;
}