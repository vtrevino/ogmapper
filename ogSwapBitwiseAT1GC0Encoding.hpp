/* 
 * File:   ogSwapBitwiseAT1GC0Encoding.hpp
 * Author: victortrevino
 *
 * Created on July 5, 2022, 6:38 PM
 */

#ifndef OGSWAPBITWISEAT1GC0ENCODING_HPP
#define OGSWAPBITWISEAT1GC0ENCODING_HPP

#include "ogKeyEncoding.hpp"
#define OGBWENC_ALL_ASCII   256


class ogSwapBitwiseAT1GC0Encoding : public ogKeyEncoding {
    char                ascii2bits[OGBWENC_ALL_ASCII];
    uint32_t            keyCode;
    uint32_t            keyCodeRC;
    uint32_t            maskCode;

public:
    ogSwapBitwiseAT1GC0Encoding();
    virtual ~ogSwapBitwiseAT1GC0Encoding();
        
    void                setSizeInChars(uint16_t size);
    const char         *getShortExtensionName();
    uint32_t            getFwdKey(char *pSeq);
    uint32_t            getRevKey(char *pSeq);
    char                isKeyMatchExactSequence();
    char                canKeyBeComparedInSequence();

private:

};

#endif /* OGSWAPBITWISEAT1GC0ENCODING_HPP */

