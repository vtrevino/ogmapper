/* 
 * File:   ogHomoPolymerCompressed.hpp
 * Author: victortrevino
 *
 * Created on August 19, 2022, 11:16 PM
 */

#ifndef OGHOMOPOLYMERCOMPRESSED_HPP
#define OGHOMOPOLYMERCOMPRESSED_HPP

#include <cstring>

#define OGKREPENC_ALL_ASCII   256


class ogHomoPolymerCompressedEnc : public ogKeyEncoding {
    uint32_t           ascii2bits[OGKREPENC_ALL_ASCII];
    uint32_t           ascii2bitsRevComp[OGKREPENC_ALL_ASCII];
    
public:
                         ogHomoPolymerCompressedEnc();
    virtual             ~ogHomoPolymerCompressedEnc();
    const char          *getShortExtensionName();
    void                 setSizeInChars(uint16_t size);
    uint32_t             getFwdKey(char *pSeq);
    uint32_t             getRevKey(char *pSeq);
    char                 isKeyMatchExactSequence();
    uint32_t             getTotalKeys();
    char                 canKeyBeComparedInSequence();

private:

};

#endif /* OGHOMOPOLYMERCOMPRESSED_HPP */

