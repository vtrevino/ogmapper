/* 
 * File:   ogReadsMapper.hpp
 * Author: victortrevino
 *
 * Created on May 20, 2022, 1:25 PM
 */

#ifndef OGREADSMAPPER_HPP
#define OGREADSMAPPER_HPP

#include <string>
#include <cstring>
#include <math.h>
#include <time.h>
#include <thread>
#include <queue>
#include <mutex>
#include <zlib.h>
#include "ogDefinitions.hpp"
#include "sequenceRead.h"
#include "CircularArray.hpp"
#include "ogGuider.hpp"
#include "ogAligner.hpp"
#include "ogMappingParameters.hpp"
#include "ogSingleRead.hpp"
#include "ogScheduler.hpp"
#include "ogFastAQGZreader.hpp"
#include "kseq.h"


typedef struct OG_READ_MAPPER_COUNTS {
    uint64_t                calls;
    uint64_t                hitting;
    uint64_t                hits;
    uint64_t                matches;
    uint64_t                foundFwd;
    uint64_t                foundRev;
    uint64_t                elapsed;
    uint64_t                failures;
} ogReadMapperCounts;

typedef struct OG_QUEUE_MAPPER {
    //kstring_t       *pSeq1, *pSeq2, *pQual1, *pQual2, *pName1, *pName2;
    ogSingleRead        *pRead;
    ogFastAQsequence    *pFAS1, *pFAS2;
    uint64_t            readNum;
} ogQueueMapper;

class ogReadsMapper {
    
public:
    ogMappingParameters     *pMapParams; // Shared
    ogGuider                *pClonedGuider;
    ogAligner               *pAligner;
    ogScheduler             *pScheduler;
    uint64_t                nKeys;
    uint64_t                nMatches;
    uint64_t                nReads;
    uint64_t                nMaps;
    uint64_t                nFound;
    uint64_t                nFoundFwd;
    uint64_t                nFoundRev;
    uint64_t                nReadsMappedR1R2;
    uint64_t                nReadsMappedR2R1;
    uint64_t                nReadsMappedR1;
    uint64_t                nReadsMappedR2;
    uint64_t                nReadsMappedTrans;
    uint64_t                nReadsMappedOther;
    uint64_t                nReadsUnmapped;
    uint64_t                nReadsAlternMaps; // no se suman a todos
    ogReadKeyMapping        *pRdKeyMap[2]; // Exclusive to this object
    ogReadKeyMapping        *pRdKeyMapRd1, *pRdKeyMapPair;
    ogSAM                   sam[2];
    uint64_t                elapsed;
    uint64_t                preparation;
    uint64_t                pairing;
    ogReadMapperCounts      *pCounts;
//    uint64_t                *tCalls;
//    uint64_t                *tMatches;
//    uint64_t                *tFoundFwd;
//    uint64_t                *tFoundRev;
//    uint64_t                *tElapsed;
    uint16_t                nFunc;
    char                    busy;
//    char                    canWork;
    char                    finish;
    //char                    *seq1;
    //uint32_t                seq1Len;
    //char                    *seq2;
    //uint32_t                seq2Len;
    uint64_t                wastingCycles;
    uint64_t                wastingCycles2;
    uint64_t                wastingEvents;
    queue<ogQueueMapper *>   qSeq;
    //CircularArray<sequenceRead *>   qSeq;
    uint16_t                maxQsize;
    std::mutex              mtx;
    uint16_t                thread;
    char                    productionMode;
    uint32_t        nKeysLow;
    uint32_t        nKeysHigh;
    char            checkMatch;
    char            compareKey;
    char            isPaired;
    char            *DNA;       // used to unpack genome DNA sequences for "alignment"....
    uint32_t        dnaAlloc;
    void            (ogReadsMapper::*processingFunc)(ogSingleRead *theRead);
    char            reverseCIGAR[MAX_CIGAR_SIZE];
    
                    ogReadsMapper(ogMappingParameters *pMapParameters, ogScheduler *pScheduling, ogAligner *pAlign, ogSamWriter *pSamWri,  CountingInfo *pCountInfo, uint32_t maxReadLen, uint16_t threadNum, float meanLog10, float sdLog10, uint32_t q95, char prodMode, char isPairedReads);
    virtual        ~ogReadsMapper();
    
    void            process1Read();                     // Take 1 read from the queue and call mapRead/countRead. (for multi-threading)
    void            processTheRead(ogSingleRead *seqR);
    void            mapRead(ogSingleRead *theRead);     // Procedure for mapping
    void            countRead(ogSingleRead *theRead);   // Procedure for counting
    uint32_t        mapTheReadWithKeyMap(ogReadKeyMapping *pRdKeyMap);
    uint32_t        countTheReadWithKeyMap(ogReadKeyMapping *pRdKeyMap);
    //void            mapReadPairs(char *seq1, uint32_t seq1Len, char *seq2, uint32_t seq2Len);
    void            printSummary(uint64_t totalReads);
    char            isBusy();
    char            isAvailable();
    char            hasBeenBusy();
    void            waste(char newEvent);
    //void            pushReadToQueue(char *pid, uint32_t pidLen, char *pseq, uint32_t pseqLen, char *pqual, uint32_t pqualLen, uint64_t readNum, char pack);
    //void            pushReadToQueue(kstring_t *name1, kstring_t *qual1, kstring_t *seq1, kstring_t *name2, kstring_t *qual2, kstring_t *seq2, uint16_t *pThreader, uint64_t readNum); //ogSingleRead *theRead
    void            pushReadToQueue(ogSingleRead *pRead, uint64_t readNum); //ogSingleRead *theRead
    void            pushReadToQueue(ogFastAQsequence *pFAS1, ogFastAQsequence *pFAS2, uint64_t readNum); //ogSingleRead *theRead
    void            setMaxQsize(uint16_t mxq);
    void            setProductionMode(char isProduction);
    char            isQueueFull(char secureCheck);
    char            isQueueEmpty();
    char            isQueueCritic();
    uint16_t        candidatesIntersectPair(ogReadKeyMapping *pRdKeyMap1, ogReadKeyMapping *pRdKeyMap2, char print);
    void            writeSamFromBothReadsInfo(sequenceRead *seqRead, ogSamWriter *pSamWri);
    void            writeSamFromRead1MapRead2Unmapped(sequenceRead *seqRead, ogSamWriter *pSamWri);
    void            writeSamFromRead2MapRead1Unmapped(sequenceRead *seqRead, ogSamWriter *pSamWri);
    void            setSAMfromRead(ogSAM *pSAM, ogSingleRead *read);
    void            writeSamFromGenomeAndReadPositions(ogReadKeyMapping *pRKM1, ogReadKeyMapping *pRKM2, ogGenomeAndReadPosition *grp1, ogGenomeAndReadPosition *grp2, ogSamWriter *pSamWri);
    void            writeSamFromGenomeAndReadPositionsMapUnmap(ogReadKeyMapping *pRKM1, ogReadKeyMapping *pRKM2, ogGenomeAndReadPosition *grp1, ogSamWriter *pSamWri);
    void            writeSamFromReadKeyMapsUnmapped(ogReadKeyMapping *pRKM1, ogReadKeyMapping *pRKM2, ogSamWriter *pSamWri);
    void            writeSamFromReadKeyMapsTranslocated(ogReadKeyMapping *pRKM1, ogReadKeyMapping *pRKM2, ogSamWriter *pSamWri);
    void            writeAltRd2SamFromGenomeAndReadPosition(ogReadKeyMapping *pRKM1, ogReadKeyMapping *pRKM2, ogGenomeAndReadPosition *grp1, ogGenomeAndReadPosition *grp2, ogSamWriter *pSamWri);
    void            writeAltRd1SamFromGenomeAndReadPosition(ogReadKeyMapping *pRKM1, ogReadKeyMapping *pRKM2, ogGenomeAndReadPosition *grp1, ogGenomeAndReadPosition *grp2, ogSamWriter *pSamWri);
    void            writeSamFrom1UnmappedRead(ogReadKeyMapping *pRKM1, ogSamWriter *pSamWri);
    void            writeSamFrom1GenomeAndReadPosition(ogReadKeyMapping *pRKM1, ogGenomeAndReadPosition *grp1, ogSamWriter *pSamWri);
    void            allocDNA(uint32_t len);
    void            buildPlainCIGARorWFA(ogReadKeyMapping *pRKM1, ogSingleRead *pR1, uint32_t genomicPos, char isRev);
    void            setMode(char mode);
    char            *revertCIGAR(char *pCigar, char *pDest);
    void            cutCIGARtoSeqLen(char *pCigar, uint32_t maxLen);
    void            correctLeftRighPositionFromCIGAR(char *pCigar, uint16_t readPos, uint32_t *pLeft, uint32_t *pRight, char isRd2);
    
private:

};

#endif /* OGREADSMAPPER_HPP */

