/* 
 * File:   ogGappedBitwiseAT1GC0Encoding.hpp
 * Author: victortrevino
 *
 * Created on July 3, 2022, 12:37 PM
 */

#ifndef OGGAPPEDBITWISEAT1GC0ENCODING_HPP
#define OGGAPPEDBITWISEAT1GC0ENCODING_HPP

#include "ogKeyEncoding.hpp"
#define OGBWENC_ALL_ASCII   256


class ogGappedBitwiseAT1GC0Encoding : public ogKeyEncoding {
    char                ascii2bits[OGBWENC_ALL_ASCII];
    uint32_t            keyCode;
    uint32_t            keyCodeRC;
    uint32_t            maskCode;
    char                firstKey;
    uint16_t            nLeft, nGap, nRight;

public:
        
                        ogGappedBitwiseAT1GC0Encoding();
    void                setSizeInChars(uint16_t size);
    const char         *getShortExtensionName();
    uint32_t            getFwdKey(char *pSeq);
    uint32_t            getRevKey(char *pSeq);
    char                isKeyMatchExactSequence();
    char                canKeyBeComparedInSequence();

private:

};

#endif /* OGGAPPEDBITWISEAT1GC0ENCODING_HPP */

