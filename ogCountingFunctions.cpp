
#include <stdio.h>
#include <map>
#include "ogDefinitions.hpp"
#include "ogIndex.hpp"
#include "ogKeyEncoding.hpp"
#include "ogReadKeyMapping.hpp"
#include "ogReusableBinaryTree.hpp"
#include "ogReadKeyMapping.hpp"
#include "ogCandidatePosManager.hpp"
#include "binnedBigCounters.hpp"



// Half Apart, left keys are "paired" with right keys in order
uint16_t HalfApartCounting(ogReadKeyMapping *pMap, uint32_t maxTargets) {

    uint32_t nReadPos, halfPos;
    uint32_t i, iPos, intC;
    uint32_t j;
    uint8_t  intraFails;
    char    x, bug = pMap->debug;
    
    pMap->pCandPosMan->reset();
    
    for (x = 0; x < 2; x++) {
        if ((nReadPos = (x == 0 ? pMap->getNFwdPos() : pMap->getNRevPos())) > 1) {
            //fprintf(stderr, "nKeys=%d\n",nReadPos); fflush(stderr);
            halfPos = (nReadPos >> 1);

            intC = 0;

            for (i=0; i < halfPos; i++) {
                if (pMap->getTargetSize(x, i) < maxTargets) {
                    iPos = pMap->getKeyablePosition(x, i);
                    intC = 0;
                    intraFails = 0;
                    for (j=halfPos+i; j < nReadPos; j++) {
                        if (pMap->getTargetSize(x, j) < maxTargets) {
                            if (intC == 0) {
                                intC = pMap->intersectByTwoReadPositionKeys(x, i, j, 10000); // pMap->getKeyablePosition(x,j) - iPos + 3
                                if (intC > 10 || intC == 0) {
                                    intC = 0;
                                    if (++intraFails > 3) break;                                
                                }
                            } else {
                                intC = pMap->intersectByAddingReadPositionKey(x, j, 10000); // pMap->getKeyablePosition(x,j) - iPos + 3
                                if (intC == 0) {
                                    intC = pMap->pCandPosMan->rewindCountsAfterAddingKeyCountZero();
                                    if (++intraFails > 3) break;
                                }
                            }
                            if (intC > 0 && intC < 3) {
                                break;
                            }
                        }
                    }
                    if (intC > 0 && intC < 3) {
                        break;
                    }
                }
            }
        }
        if (x == 0) {
            if (intC == 0) { // (pMap->pCandPosMan->getCount() > 2)
                pMap->pCandPosMan->reset();
            }
            pMap->pCandPosMan->saveFwdCandidatePositions();
        } else {
            if (intC == 0) { // (pMap->pCandPosMan->getRevCount() > 2)
                pMap->pCandPosMan->restoreFwdCandidatePositions();
            }
        }
    }

    return pMap->pCandPosMan->getCount();
}

uint16_t HalfCountingSmall(ogReadKeyMapping *pMap) {
    return HalfApartCounting(pMap,  pMap->getLowKeyCountLimit());
}

uint16_t HalfCountingLarge(ogReadKeyMapping *pMap) {
    return HalfApartCounting(pMap,  pMap->getHighKeyCountLimit());
}

uint16_t HalfCountingHuge(ogReadKeyMapping *pMap) {
    return HalfApartCounting(pMap, 10*pMap->getHighKeyCountLimit());
}
