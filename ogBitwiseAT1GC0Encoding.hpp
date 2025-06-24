
/* 
 * File:   ogBitwiseAT1GC0Encoding.hpp
 * Author: victortrevino
 *
 * Created on May 7, 2022, 3:10 PM
 */

#ifndef OGBITWISEAT1GC0ENCODING_HPP
#define OGBITWISEAT1GC0ENCODING_HPP

#define OGBWENC_ALL_ASCII   256

#include "ogDefinitions.hpp"
#include "ogKeyEncoding.hpp"

class ogBitwiseAT1GC0Encoding : public ogKeyEncoding {
    char                ascii2bits[OGBWENC_ALL_ASCII];
    uint32_t            keyCode;
    uint32_t            keyCodeRC;
    uint32_t            maskCode;
    char                firstKey;
    
public:
                        ogBitwiseAT1GC0Encoding();
    void                setSizeInChars(uint16_t size);
    const char         *getShortExtensionName();
    uint32_t            getFwdKey(char *pSeq);
    uint32_t            getRevKey(char *pSeq);
    char                 isKeyMatchExactSequence();
    
private:
    //void                getLocalFRKeys(char *p);
    //void                shiftLocalFRKeys(char c);

};

#endif /* OGBITWISEAT1GC0ENCODING_HPP */

