/* 
 * File:   ogReadKeyMapping.hpp
 * Author: victortrevino
 *
 * Created on May 21, 2022, 12:52 PM
 */

#ifndef OGREADKEYMAPPING_HPP
#define OGREADKEYMAPPING_HPP

#include <stdint.h>
#include "ogKeys.hpp"
#include "ogGenome.hpp"
#include "ogGenomePositions.hpp"
#include "ogKeyEncoding.hpp"
#include "ogCandidatePosManager.hpp"
#include "ogSingleRead.hpp"
#include "ogSamWriter.hpp"
#include "ogMappingParameters.hpp"

typedef struct OG_KEY_MAP {
    uint32_t    position; // relative to read
    char        active;
    ogKey       keyInfo;
    OG_KEY_MAP    *next;
    uint32_t    index;
} ogKeyMap;


class ogReadKeyMapping {
        
public:
    
    ogKeyMap                *pForwardKeys, *pFwdFirst, *pFwdK, *pFwdLast;
    ogKeyMap                *pReverseKeys, *pRevFirst, *pRevK, *pRevLast;
    
    char                     completed;          // 100 is completed, 1 is partially set to positions
    uint32_t                 nFwdPos;
    uint32_t                 nRevPos;
    uint32_t                 nActiveFwd;
    uint32_t                 nActiveRev;
    uint32_t                 minFwdKeySize, maxFwdKeySize;
    uint32_t                 minFwdKeySizePos, maxFwdKeySizePos;
    uint32_t                 minRevKeySize, maxRevKeySize;
    uint32_t                 minRevKeySizePos, maxRevKeySizePos;
    //char                    *seq;
    //uint32_t                 seqLen;
    //uint64_t                 seqNum;
    ogSingleRead            *read;
    char                    mapped; // has the read being mapped
    char                     isRead2;
    uint32_t                 maxReadLength;
    uint8_t                  keySize;
    //char                     *packedReadSeqFwd[4];
    //char                     *packedReadSeqRev[4];
    //uint32_t                 packLenAllocated;
    ogCandidatePosManager   *pCandPosMan;
    uint32_t                maxScore, maxScorePos, maxTies;

    // copy from MapParams to access them a bit quicker
    ogKeys                  *pKeys;
    ogGenome                *pGenome;
    ogKeyEncoding           *pKeyEncoding;
    ogGenomePositions       *pGenPos;
    ogGuider                *pGuider;
    ogSamWriter             *pSamWriter;
    char                     forceUpperCase;
    uint16_t                 thread;
    float                    meanKeys_log10;
    float                    sdKeys_log10;
    uint32_t                nKeysLow;
    uint32_t                nKeysHigh;
    char                    debug;
    uint64_t                sharedMemoryBytes;
    void                    *pMemory;
    char                    isMemoryReleasable;
    ogReadKeyMapping        *pLinkedRKM;
    char                    *DNA;
    CountingInfo            *pCountInfo;
    uint32_t                dnaAlloc;
    ogMappingParameters     *pMapParams;
    
                        ogReadKeyMapping(uint32_t maxReadLen, uint16_t threadNum);
    virtual            ~ogReadKeyMapping();
    void                setPosReadKey(uint32_t pos);
    void                setFirstReadKey();
    uint32_t            setLastReadKey();
    void                setAllReadKeys();
    void                setReadKeysFromToBy(uint32_t from, uint32_t to, uint32_t by);
    void                setRead(ogSingleRead *theRead);
    void                setMasterObjects(ogKeys *pK, ogGenome *pG, ogKeyEncoding *pKE,ogGenomePositions *pGP, ogGuider *pGuide, char forceUpper, ogSamWriter *pSamW, ogMappingParameters *pMapPar);
    void                setStatistics(float meanLog10, float sdLog10, uint32_t q95);
    uint32_t            getMeanKeyPlusStdDev(float nStdDev);

    ogCandidatePosManager* getCandPosMan();
    uint32_t            getLastKeyablePosition();
    uint32_t            getKeyForward(uint32_t pos);
    uint32_t            getKeyReverse(uint32_t pos);
    uint32_t            intersectByTwoReadKeys(uint32_t k1, uint32_t k2, uint16_t maxDif);
    uint32_t            intersectByTwoReadPositionKeys(char isReverse, uint32_t p1, uint32_t p2, uint16_t maxDif);
    uint32_t            intersectByAddingReadKey(uint32_t k1, uint16_t maxDif);
    uint32_t            intersectByAddingReadPositionKey(char isReverse, uint32_t p1, uint16_t maxDif);
    uint32_t            intersectByAddingReadPositionKeyPositiveDistance(char isReverse, uint32_t p1, uint16_t maxDif);
    uint32_t            intersectByAddingReadKeyPositiveDistance(uint32_t k1, uint16_t maxDif);
    uint32_t            intersectByTwoKeys(ogKey *pK1, ogKey *pK2, uint16_t maxDif, char isReverse);
    uint32_t            intersectByTwoKeys(ogKey *pK1, ogKey *pK2, uint32_t rdPos1, uint32_t rdPos2, uint16_t maxDif, char isReverse);
    uint32_t            intersectByTwoKeysLeft(ogKey *pK1, ogKey *pK2, uint32_t rdPos1, uint16_t maxDif, char isReverse, uint16_t maxAdd);
    uint32_t            getKeyFR(uint32_t pos, char isRev);
    ogKeyMap            *getKeyMapForPosition(char isRev, uint32_t pos);
    ogKey*              getInfoForPositionKey(char isReverse, uint32_t p1);
    ogKey*              getInfoForKey(uint32_t k1);
    ogKey*              getFwdInfoForKeyablePosition(uint32_t iKey);
    ogKey*              getRevInfoForKeyablePosition(uint32_t iKey);
    uint32_t*           getGenomePositionFromReadPositionKey(char isReverse, uint32_t p1, uint32_t iPos);
    uint32_t            getNFwdPos();
    uint32_t            getNRevPos();
    uint32_t            getNKeyPos(char isReverse);
    uint32_t            getFwdKeyablePosition(uint32_t iKey);
    uint32_t            getRevKeyablePosition(uint32_t iKey);
    void                printKeyInfo(uint32_t k);
    void                clearFwdKeyablePosition(uint32_t pos);
    void                clearRevKeyablePosition(uint32_t pos);
    uint32_t            getKeyablePosition(char RC, uint32_t iKey);
    void                reactivateKeys();
    void                inactivePositionsAroundKeyRev(uint32_t p);
    void                inactivePositionsAroundKeyFwd(uint32_t p);
    void                inactivePositionFwd(uint32_t p);
    void                inactivePositionRev(uint32_t p);
    char                isFwdKeyActive(uint32_t p);
    char                isRevKeyActive(uint32_t p);
    //uint32_t            getRelativePos(char isReverse);
    void                printGenomeInfoFromKeyPos(uint32_t pos);
    
    uint32_t            getFwdTargetSize(uint32_t iKey);
    uint32_t            getRevTargetSize(uint32_t iKey);
    uint32_t            getTargetSize(char RC, uint32_t iKey);

    void                printNotFound(uint64_t readNum);
    uint32_t            getLowKeyCountLimit();
    uint32_t            getHighKeyCountLimit();
    uint32_t            getCloserGenomePosition(uint32_t targetGPos, uint32_t *pa, uint32_t aSize);
    uint32_t            getCloserGenomeRelativePosition(uint32_t targetGPos, uint32_t *pa, uint32_t aSize);
    void                resetByOrderSize(char RC);
    uint32_t            getIndexByOrderSize(char RC);
    char                nextByOrderSize(char RC);
    char                finishByOrderSize(char RC);
    ogKey              *keyByOrderSize(char RC);
    
    void               *getUsableMemory(uint64_t bytes);
    void                allocDNA(uint32_t len);
    
    
private:

};

#endif /* OGREADKEYMAPPING_HPP */

