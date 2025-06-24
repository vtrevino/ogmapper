/* 
 * File:   ogIndex.hpp
 * Author: victortrevino
 *
 * Created on May 7, 2022, 1:59 PM
 */

#ifndef OGINDEX_HPP
#define OGINDEX_HPP

#include <pthread.h>
#include <map>
#include "ogDefinitions.hpp"
#include "ogKeys.hpp"
#include "ogKeyEncoding.hpp"
#include "ogGenome.hpp"
#include "ogGenomePositions.hpp"
#include "ogReadsMapper.hpp"
#include "ogAligner.hpp"
#include "ogScheduler.hpp"

#include "CircularArray.hpp"

#include "ogKeyEncoding.hpp"

typedef struct ogIndexHeader {
    char            indexVersion[MAX_VERSION_NAME];             // Version of this structure
    char            sourceFileName[MAX_OGINDEX_FILENAME];       // Source File name
    char            encodingName[MAX_VERSION_NAME];             // Which of the registered Encodings this index refer to
    float           encodingVersion;                            // Encoding version
    char            guiderName[MAX_VERSION_NAME];               //
    uint32_t        genomeSize;                                 // Maximum Genome Size: 4,294,967,296, informative only
    uint16_t        userKeySize;
    float           meanKeys_log10;
    float           sdKeys_log10;
    uint32_t        nKeyCountsAtQ95;
} OG_INDEX_HEADER;



class ogIndex {

    char            indexFileName[2048];
    // For Index Creation & Loading
    ogIndexHeader   header;
    char            *pEncodingNames[MAX_ENCODINGS];
    ogKeyEncoding   *pEncodingClasses[MAX_ENCODINGS];
    int             nEncodings;
    char            selEncoding[MAX_CHAR_NAMES];
    char            selGuider[MAX_CHAR_NAMES];
    char            guideFile[MAX_FILENAME];
    char            lowMemoryMode;
    char            forceUpperCase;
    char            kseqReadingMode;
    uint32_t        maxReadLen;
    uint64_t        maxReads;
    uint64_t        startingRead;
    uint32_t        nKeysLow;
    uint32_t        nKeysHigh;
    // For Mapping
    char            schedule[MAX_SCHEDULE_SIZE];
    uint16_t        threads;
    std::map<std::string, ogKeyEncoding*> protoKEtable;                             // taken from http://www.cs.sjsu.edu/~pearce/modules/lectures/oop/types/reflection/prototype.htm
    std::map<char, uint16_t (*)(ogReadKeyMapping *)> protoMapFuncTable;             // taken from http://www.cs.sjsu.edu/~pearce/modules/lectures/oop/types/reflection/prototype.htm
    std::map<std::string, ogGuider*> protoGuiderTable;
    
    uint16_t        nThreads;
    int             nQueue;
    ogReadsMapper   **rdMapr;
    std::thread     **hilos;
    char            productionMode;
    
public:

    ogKeys              *pKeys;
    ogGenome            *pGenome;
    ogKeyEncoding       *pEncoding;
    ogGenomePositions   *pGenomePositions;
    ogMappingParameters *pMapParams;
    ogScheduler         scheduler;
    ogGuider            *pGuider;
    ogSamWriter         *pSamSam;
    char                isExonMode;
    CountingInfo        countInfo;
    char                exonsCountFile[MAX_FILENAME];
    char                genesCountFile[MAX_FILENAME];
    char                transcriptsCountFile[MAX_FILENAME];
    char                doPrepare; // Prepare reads : this is for checking reading speed
    char                doProcess; // Process reads : this is for checking read & preparing speed
    uint32_t            maxReadsBuffer;
    char                filePerThread;
    char                unmappedMode;
    char                readGroupHeader[MAX_FILENAME];

                        ogIndex();
    virtual            ~ogIndex();
    void                registerKeyEncoding(ogKeyEncoding *pEncodingClass);
    void                registerGuider(ogGuider *pGuiderClass);
    void                selectKeyEncoding(char *pEncoding);
    void                selectGuider(char *pGuider);
    void                generate(char *pSourceFileName, uint16_t userKeySize, char *pDestFileName, int argc, char *argv[]);
    char                load(char *pSourceFileName, uint16_t *userKeySize);
    void                setMemoryMode(char lowMemMode);
    ogKeyEncoding*      getKeyEncoding(char *pEncoding);
    ogGuider*           getGuider(char *pGuider);
    uint32_t            getMeanKeyPlusStdDev(float nStdDev);
    void                setThreads(uint16_t maxThreads);
    void                setQueue(uint16_t maxQueue);
    void                registerMap();
    void                setMaxReadLength(uint32_t mxRdLen);
    void                setMaxReads(uint64_t mxRds);
    void                setSchedule(char *pSchedule);
    void                setKSEQmode(char ks);
    void                map(char *pSourceFileName1, char *pSourceFileName2);
    void                registerMappingFunction(char key, uint16_t (*pFunc)(ogReadKeyMapping *));
    //int                 mapSequences(char *filename, uint16_t  threads);
    void                deleteKeyEncodings();
    void                deleteGuiders();
    void                setStartingRead(uint64_t stRead);
    void                threadingReads(uint16_t iThread);
    void                testRepes(char *pSourceFileName, uint16_t repes);
    void                getChr20_22(char *pSourceFileName);
    void                setUpperCase(char upcase);
    void                testStateMachine(char *pSourceFileName, char *pStateFileName, uint16_t maxLenNoResponse);
    void                setProductionMode(char isProduction);
    uint32_t            getLowKeyCountLimit();
    uint32_t            getHighKeyCountLimit();
    void                openOutputSAM(char *pFileName, char *fname1, char *fname2);
    void                closeOutputSAM();
    void                indexFromGTF(char *pGtfFileName, char *pGenomeFileName, char *outFileName);
    void                count(char *pGTFfileName, char *pSourceFileName1, char *pSourceFileName2);
    void                mapOrCount(char *pSourceFileName1, char *pSourceFileName2, char mode);
    void                loadExonAndTranscriptFiles(CountingInfo *pCountingInfo);
    Exon*               findExonForPseudoPosition(uint32_t xgene, uint32_t xpos);
    uint32_t            findExonIndexForPseudoPosition(uint32_t xgene, uint32_t xpos);
    FILE                *openIndexFileAtSection(char section);
    void                setSectionSizeInIndexFile(char section, uint64_t size);
    void                setSectionStartInIndexFile(char section, uint64_t start);
    void                setSectionSizeStartInIndexFile(char section, uint64_t size, uint64_t start);
    long                get_file_size(char *filename);
    void                checkKSEQ(char *filename);
    
private:

};

#endif /* OGINDEX_HPP */

