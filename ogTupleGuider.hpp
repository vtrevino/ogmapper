/* 
 * File:   ogTupleGuider.hpp
 * Author: victortrevino
 *
 * Created on July 13, 2022, 1:52 AM
 */

#ifndef OGTUPLEGUIDER_HPP
#define OGTUPLEGUIDER_HPP

#include <unistd.h>
#include <stdio.h>
#include "ogGuider.hpp"

#define OGTUPLE_GUIDER_CONTROL_FILE  "TupleGuiderControlFile.og"

class ogTupleGuider  : public ogGuider {
    char                extName[20];
    int16_t             maxState;
    char               *isGuider;
    char               *ntAdvances;
    char                ascii2n[256];
    char                tupleFileName[MAX_FILENAME];
    int16_t             tupleLen;
    int16_t             tupleMask;
    
public:
                            ogTupleGuider();
    virtual                ~ogTupleGuider();
    ogTupleGuider          *clone();
    void                    initializeForIndexing();
    uint32_t                load(FILE *pFile);
    uint32_t                save(FILE *pFile);
    char                    nextGuide();
    const char             *getName();
    const char             *getShortExtensionName();
private:

};

#endif /* OGTUPLEGUIDER_HPP */

