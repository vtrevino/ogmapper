
/* 
 * File:   ogKeyEncoding.hpp
 * Author: victortrevino
 *
 * Created on May 7, 2022, 2:37 PM
 */

#ifndef OGKEYENCODING_HPP
#define OGKEYENCODING_HPP

#include <stdint.h>
#include <stdio.h>
#include "ogGuider.hpp"

#define INVALID_KEY     0xFFFFFFFF

class ogKeyEncoding {
    
public:
    char            name[100];
    float           version;
    uint16_t        sizeInChars;                        // set by user
    uint8_t         sizeInBits;                         // Calculated
    char            validKey;
    
                            ogKeyEncoding();
    virtual                ~ogKeyEncoding();
    virtual void            setSizeInChars(uint16_t  size);
    virtual uint16_t        getSizeInChars();
    virtual uint8_t         getSizeInBits();
    virtual char           *getName();
    virtual const char     *getShortExtensionName();
    virtual uint32_t        save(FILE *pFile);
    virtual void            load(FILE *pFile);
    virtual uint32_t        getFwdKey(char *pSeq); // 5'->3'
    virtual uint32_t        getRevKey(char *pSeq); // 3'->5'
    virtual char            isValidKey();
    virtual char            isKeyMatchExactSequence(); // This determine whether the Key represent the sequence EXACTLY
    virtual uint32_t        getTotalKeys(); // 2^getSizeInBits()
    virtual char            canKeyBeComparedInSequence(); // This determine whether the Key represent the sequence EXACTLY and can be compared directly
    
    
private:

};

#endif /* OGKEYENCODING_HPP */

