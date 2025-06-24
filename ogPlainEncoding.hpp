/* 
 * File:   ogPlainEncoding.hpp
 * Author: victortrevino
 *
 * Created on July 3, 2022, 12:02 PM
 */

#ifndef OGPLAINENCODING_HPP
#define OGPLAINENCODING_HPP

#include <cstring>
#include "ogKeyEncoding.hpp"

#define OGKREPENC_ALL_ASCII   256


class ogPlainEncoding : public ogKeyEncoding {
    uint32_t           ascii2bits[OGKREPENC_ALL_ASCII];
    uint32_t           ascii2bitsRevComp[OGKREPENC_ALL_ASCII];
    
public:
                         ogPlainEncoding();
                         ogPlainEncoding(const ogPlainEncoding& orig);
    virtual             ~ogPlainEncoding();
    const char          *getShortExtensionName();
    void                 setSizeInChars(uint16_t size);
    uint32_t             getFwdKey(char *pSeq);
    uint32_t             getRevKey(char *pSeq);
    char                 isKeyMatchExactSequence();
    char                 canKeyBeComparedInSequence();

private:

};

#endif /* OGPLAINENCODING_HPP */

