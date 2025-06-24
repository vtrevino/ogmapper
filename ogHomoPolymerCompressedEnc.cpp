/* 
 * File:   ogHomoPolymerCompressed.cpp
 * Author: victortrevino
 * 
 * Created on August 19, 2022, 11:16 PM
 */

#include <string>
#include <math.h>
#include "ogDefinitions.hpp"
#include "ogKeyEncoding.hpp"
#include "ogHomoPolymerCompressedEnc.hpp"

const char myname[] = "HPCEncoding\0";

ogHomoPolymerCompressedEnc::ogHomoPolymerCompressedEnc() {
    memset(ascii2bits, 4, sizeof(ascii2bits));
    ascii2bits['A'] = 0x00;
    ascii2bits['a'] = 0x00;
    ascii2bits['C'] = 0x01;
    ascii2bits['c'] = 0x01;
    ascii2bits['G'] = 0x02;
    ascii2bits['G'] = 0x02;
    ascii2bits['T'] = 0x03;
    ascii2bits['t'] = 0x03;
    ascii2bits['U'] = 0x03;
    ascii2bits['u'] = 0x03;
    memset(ascii2bitsRevComp, 4, sizeof(ascii2bitsRevComp));
    ascii2bitsRevComp['A'] = 0x03;
    ascii2bitsRevComp['a'] = 0x03;
    ascii2bitsRevComp['C'] = 0x02;
    ascii2bitsRevComp['c'] = 0x02;
    ascii2bitsRevComp['G'] = 0x01;
    ascii2bitsRevComp['G'] = 0x01;
    ascii2bitsRevComp['T'] = 0x00;
    ascii2bitsRevComp['t'] = 0x00;
    ascii2bitsRevComp['U'] = 0x00;
    ascii2bitsRevComp['u'] = 0x00;
    strncpy(name, myname, MAX_CHAR_NAMES);
    version = 1.0;
}

ogHomoPolymerCompressedEnc::~ogHomoPolymerCompressedEnc() {
}

const char *ogHomoPolymerCompressedEnc::getShortExtensionName() {
    return "_eHPC";
}

void ogHomoPolymerCompressedEnc::setSizeInChars(uint16_t size) {
    sizeInChars = size;
    sizeInBits = sizeInChars * 2;
}

uint32_t ogHomoPolymerCompressedEnc::getFwdKey(char *pSeq) {
    uint32_t    key = 0;
    uint8_t     i;
    char        prev = 'x';
    char        a, nt;
    validKey    = 1;
    for (i=0; *pSeq && i < sizeInChars; pSeq++) {
        if (*pSeq != prev) {
            if ((a = ascii2bits[*pSeq]) == 4) {
                validKey = 0;
                break;
            }
            nt = *pSeq & 0x5F;
            if (prev == 'x') {
                key = ascii2bits[nt];
            } else if (prev == 'A') {
                // C == 1
                if (nt == 'G') key *= 2;
                else if (nt == 'T') key *= 3;
            } else if (prev == 'C') {
                // A == 1
                if (nt == 'G') key *= 2;
                else if (nt == 'T') key *= 3;                
            } else if (prev == 'G') {
                // A == 1
                if (nt == 'C') key *= 2;
                else if (nt == 'T') key *= 3;
            } else if (prev == 'T') {
                // A == 1
                if (nt == 'C') key *= 2;
                else if (nt == 'G') key *= 3;                
            }
            prev = nt;
            i++;
        }
    }
    if (i < sizeInChars) validKey = 0;
    return key;
}

uint32_t ogHomoPolymerCompressedEnc::getRevKey(char *pSeq) {
    uint32_t    key = 0;
    uint8_t     i;
    char        prev = 'x';
    char        a, nt;
    validKey = 1;
    for (i=0; *pSeq && i < sizeInChars; pSeq--) {
        if (*pSeq != prev) {
            if ((a = ascii2bits[*pSeq]) == 4) {
                validKey = 0;
                break;
            }
            nt = *pSeq & 0x5F;
            if (prev == 'x') {
                key = ascii2bits[nt];
            } else if (prev == 'A') {
                // C == 1
                if (nt == 'G') key *= 2;
                else if (nt == 'T') key *= 3;
            } else if (prev == 'C') {
                // A == 1
                if (nt == 'G') key *= 2;
                else if (nt == 'T') key *= 3;                
            } else if (prev == 'G') {
                // A == 1
                if (nt == 'C') key *= 2;
                else if (nt == 'T') key *= 3;
            } else if (prev == 'T') {
                // A == 1
                if (nt == 'C') key *= 2;
                else if (nt == 'G') key *= 3;                
            }
            prev = nt;
            i++;
        }
        /**
        if (*pSeq != prev) {
            prev = *pSeq;
            if ((a = ascii2bitsRevComp[prev]) == 4) {
                validKey = 0;
                break;
            }
            key |= (a << sF);
            if (sF == 0) sF = sizeInBits - 2; else sF -= 2;
            i++;
        }
        **/
    }
    if (i < sizeInChars) validKey = 0;
    return key;
}

char ogHomoPolymerCompressedEnc::isKeyMatchExactSequence() {
    // This determine whether the Key represent the sequence EXACTLY
    return 0;
}

char ogHomoPolymerCompressedEnc::canKeyBeComparedInSequence() {
    // This determine whether the Key represent the sequence EXACTLY
    return 0;
}

uint32_t ogHomoPolymerCompressedEnc::getTotalKeys() {
    return 4*pow(3, sizeInChars-1);
}

