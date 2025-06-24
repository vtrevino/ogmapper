/* 
 * File:   ogStateMachineGuider.hpp
 * Author: victortrevino
 *
 * Created on July 1, 2022, 8:10 PM
 */

#ifndef OGSTATEMACHINEGUIDER_HPP
#define OGSTATEMACHINEGUIDER_HPP

#include <unistd.h>
#include <stdio.h>
#include "ogGuider.hpp"

#define OGSTATE_MACHINE_CONTROL_FILE  "StateMachineControlFile.og"

typedef int16_t  TransitionVector[6]; // A, C, G, T, N, *

class ogStateMachineGuider : public ogGuider {
    char                extName[20];
    int16_t             maxState;
    TransitionVector   *pTransitionMatrix;
    char                ascii2col[256];
    char                stateFileName[MAX_FILENAME];
    
    
public:
                            ogStateMachineGuider();
    virtual                ~ogStateMachineGuider();
    ogStateMachineGuider   *clone();
    void                    initializeForIndexing();
    uint32_t                load(FILE *pFile);
    uint32_t                save(FILE *pFile);
    char                    nextGuide();
    const char             *getName();
    const char             *getShortExtensionName();

private:

};

#endif /* OGSTATEMACHINEGUIDER_HPP */

