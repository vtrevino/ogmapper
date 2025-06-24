/* 
 * File:   ogPlainEncoding.cpp
 * Author: victortrevino
 * 
 * Created on July 3, 2022, 12:02 PM
 */

#include <string>
#include "ogDefinitions.hpp"
#include "ogKeyEncoding.hpp"
#include "ogPlainEncoding.hpp"

const char myname[] = "PlainEncoding\0";


ogPlainEncoding::ogPlainEncoding() {
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

ogPlainEncoding::~ogPlainEncoding() {
}

const char *ogPlainEncoding::getShortExtensionName() {
    return "_ePlain";
}

void ogPlainEncoding::setSizeInChars(uint16_t size) {
    sizeInChars = size;
    sizeInBits = sizeInChars * 2;
}

uint32_t ogPlainEncoding::getFwdKey(char *pSeq) {
    uint32_t    key = 0;
    uint8_t     i;
    uint8_t     sF = sizeInBits - 2;
    uint32_t    nt;
    validKey = 1;
    for (i=sizeInChars; i; i--) {
        if ((nt = ascii2bits[*pSeq++]) > 3) { validKey = 0; nt = 0; }
        key |= (nt << sF);
        if (sF == 0) sF = sizeInBits - 2; else sF -= 2;
    }
    return key;
}

uint32_t ogPlainEncoding::getRevKey(char *pSeq) {
    uint32_t    key = 0;
    uint8_t     i;
    uint8_t     sF = sizeInBits - 2;
    uint32_t    nt;
    validKey = 1;
    for (i=sizeInChars; i; i--) {
        if ((nt = ascii2bitsRevComp[*pSeq--]) > 3) { validKey = 0; nt = 0; }
        key |= (nt << sF);
        if (sF == 0) sF = sizeInBits - 2; else sF -= 2;
    }
    return key;
}

char ogPlainEncoding::isKeyMatchExactSequence() {
    // This determine whether the Key represent the sequence EXACTLY
    return 1;
}

char ogPlainEncoding::canKeyBeComparedInSequence() {
    // This determine whether the Key represent the sequence EXACTLY
    return 1;
}

