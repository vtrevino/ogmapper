
/* 
 * File:   ogCandidatePosManager.hpp
 * Author: victortrevino
 *
 * Created on May 4, 2022, 4:27 PM
 */

#ifndef OGCANDIDATEPOSMANAGER_HPP
#define OGCANDIDATEPOSMANAGER_HPP

#include <stdint.h>
#include "ogReusableBinaryTree.hpp"
#include "ogGenome.hpp"
#include <cstring>

typedef struct ogGenomeAndReadPosition {
    uint32_t        genomePosition;         // global position in the genome of the readPosition (starting position is genomePosition-readPosition)
    uint16_t        readPosition;           // position within the read
    char            isReverse;
    char            status;     // u - unchecked : default
                                // s - unchecked but added with a valid score
                                // k - key was not validated 
                                // K - key validated
                                // A - Simple Alignment passed
                                // a - Simple Alignment not passed
                                // e - Simple Alignment but needs extended for global, 
                                // G - Global Alignment passed
                                // g - Global Alignment not passed
                                // P - Paired with Read 2 (score is a composite of G/A score + G/A score of the pair)
                                // (only created by the other RdKeyMapping indicating that 
                                //     this record is a candidate for record referred by the "pairedKeys" index in the the other RdKeyMapping)
    char unsigned   score;
    char            function;      // function code used to generate the Gen&RdPos
    uint16_t        pairedKeys;    // Count of the number of keys "in agreement" with this position (or the index of the paired position in the other RdKeyMapping)
    uint16_t        pairedCounts;  // Count of paired positions passing insert size filter
    uint32_t        pairedIndex;   // starting position in pairedPosition structure
    uint32_t        maxPairedIndex;// maximum paired position index 
} OG_GENOME_AND_READ_POSITION;



class ogCandidatePosManager {
public:
    ogGenomeAndReadPosition     *positions;         // array for positions, these are added by Forward Operations and Reverse Operations
    ogGenomeAndReadPosition     *pairedPositions;   // smaller array for paired positions
    uint32_t        arraySize;          // memory limit for *positions
    uint32_t        pairedArraySize;
    uint32_t        pairedArrayCount;        
    uint32_t        count;              // number of positions correctly known as positions
                                        // Forward Positions are designated from 0 to savedFwdCount, if savedFwdCount == 0 there are no fwd positions
                                        // Reverse Positions are designated from savedFwdCount to count
    uint32_t        prevCount;          // For rewind function, saves count, thus must be run before adding reverse positions
    uint32_t        savedFwdCount;      // For IntersectAdding
    uint32_t        fwdMatches;
    uint32_t        revMatches;
    unsigned char   maxScore[2];        // maxScore for 0=forward, 1=reverse;
    unsigned char   minScore[2];        // minScore for 0=forward, 1=reverse;
    uint32_t        maxScorePos[2];
    uint32_t        maxScoreTies[2];
    unsigned char   changed;
    char            currentFunction;
    
    
    ogReusableBinaryTree    *pTree;
    
                ogCandidatePosManager(uint32_t maxArraySize);
    virtual    ~ogCandidatePosManager();
    void        allocateForSize(uint32_t size);
    char        allocatePairedForSize(uint32_t size);
    uint32_t    getArraySize();
    uint32_t    setByInsersectingKeys(uint32_t *pa, uint32_t aSize, uint16_t aRdPos, uint32_t *pb, uint32_t bSize, uint16_t bRdPos, uint16_t maxDiff, char isRev);
    uint32_t    setByInsersectingKeysLeftKey(uint32_t *pa, uint32_t aSize, uint16_t aRdPos, uint32_t *pb, uint32_t bSize, uint16_t maxDiff, char isRev, uint16_t maxAdd);
    uint32_t    intersectAddingKey(uint32_t *pa, uint32_t aSize, uint16_t aRdPos, uint16_t maxDiff, char isRev);
    uint32_t    getCount();
    uint32_t    getFwdCount();
    uint32_t    getRevCount();
    ogGenomeAndReadPosition    *getkPos(uint32_t k);
    void        saveFwdCandidatePositions();
    void        restoreFwdCandidatePositions();
    uint32_t    rewindCountsAfterAddingKeyCountZero();
    void        reset();
    void        addGenomeAndReadPosition(uint32_t genPos, uint16_t readPos, char isReverse, unsigned char score);
    char        checkToAddGenomeAndReadPosition(uint32_t genPos, uint16_t readPos, char isReverse, unsigned char score);
    void        setGenomeAndReadPositionK(uint32_t k, uint32_t genPos, uint16_t readPos, char isReverse, unsigned char score);
    char        areThereMoreForwardCounts();
    void        setAllStatus(char stat);
    void        setPositionKStatus(uint32_t k, char stat);
    char        getPositionKStatus(uint32_t k);
    void        setPositionKScore(uint32_t k, char score);
    char unsigned getPositionKScore(uint32_t k);
    void        removeForwardTargets();
    void        removeReverseTargets();
    void        print();
    void        printWithGenomePositions(ogGenome *pGenome);
    void        estimateFwdRevMatches(char all);
    uint32_t    getFwdMatches();
    uint32_t    getRevMatches();
    void        cleanCandidatePositions();
//    void        saveCandidatePositions();
//    uint32_t    restoreSavedCandidatePositions();
//    uint32_t    addSavedCandidatePositions();
//    uint32_t    copySavedCandidatePositionsTo(uint32_t i);
    void        initPairedGenomeAndReadPosition(uint32_t k);
    void        initPairedGenomeAndReadPosition(ogGenomeAndReadPosition *p);
    char        checkToAddPairedGenomeAndReadPosition(ogGenomeAndReadPosition *gp, uint32_t genPos, uint16_t readPos, char isReverse, unsigned char score);
    uint32_t    addPairedGenomeAndReadPosition(ogGenomeAndReadPosition *gp, uint32_t genPos, uint16_t readPos, char isReverse, unsigned char score);
    uint32_t    addPairedGenomeAndReadPosition(uint32_t k, uint32_t genPos, uint16_t readPos, char isReverse, unsigned char score);
    char        checkToAddPairedGenomeAndReadPosition(uint32_t k, uint32_t genPos, uint16_t readPos, char isReverse, unsigned char score);
    uint32_t    checkPairedGenomeAndReadPosition(uint32_t k, uint32_t genPos, uint16_t readPos, char isReverse);
    uint32_t    checkPairedGenomeAndReadPosition(ogGenomeAndReadPosition *gp, uint32_t genPos, uint16_t readPos, char isReverse);
    ogGenomeAndReadPosition    *getkPairedPos(uint32_t k);
    void        resetPairedGenomeAndReadPosition(ogGenomeAndReadPosition *gp);

private:
    
};



#endif /* OGCANDIDATEPOSMANAGER_HPP */

