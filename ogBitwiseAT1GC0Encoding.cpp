
/* 
 * File:   ogBitwiseAT1GC0Encoding.cpp
 * Author: victortrevino
 * 
 * Created on May 7, 2022, 3:10 PM
 */

#include <stdio.h>
#include <string.h>
#include "ogKeyEncoding.hpp"
#include "ogDefinitions.hpp"
#include "ogBitwiseAT1GC0Encoding.hpp"

const char myname[] = "BitwiseAT1GC0Encoding\0";

ogBitwiseAT1GC0Encoding::ogBitwiseAT1GC0Encoding() {
    memset(ascii2bits, 2, OGBWENC_ALL_ASCII);
    ascii2bits['A'] = 1;
    ascii2bits['a'] = 1;
    ascii2bits['T'] = 1;
    ascii2bits['t'] = 1;
    ascii2bits['U'] = 1;
    ascii2bits['u'] = 1;
    ascii2bits['G'] = 0;
    ascii2bits['g'] = 0;
    ascii2bits['C'] = 0;
    ascii2bits['c'] = 0;
    strncpy(name, myname, MAX_CHAR_NAMES);
    version = 1.0;
}

const char *ogBitwiseAT1GC0Encoding::getShortExtensionName() {
    return "_eBWat1gc0";
}

void ogBitwiseAT1GC0Encoding::setSizeInChars(uint16_t size) {
    sizeInChars = size;
    sizeInBits = sizeInChars;
    keyCode = 0x00000001 << (sizeInBits - 1);
    keyCodeRC = 0x00000001;
    maskCode = 0xFFFFFFFF >> (32 - sizeInBits);    
}

uint32_t ogBitwiseAT1GC0Encoding::getFwdKey(char *pSeq) {
    uint32_t key = 0;
    uint32_t code1 = keyCode;
    uint8_t  i;
    char     a;
    validKey = 1;
    for (i=sizeInChars; i; i--) {
        if ((a = ascii2bits[*pSeq++]))  {
            if (a == 2) {
                validKey = 0;
                break;
            }
            key |= code1;
        }
        code1 >>= 1;
    }
    return key;
}

uint32_t ogBitwiseAT1GC0Encoding::getRevKey(char *pSeq) {
    // Ahora este encoding es igual porque la secuencia se recorre de 3'->5'
    uint32_t key = 0;
    uint32_t code1 = keyCode;
    uint8_t  i;
    char     a;
    validKey = 1;
    for (i=sizeInChars; i; i--) { // pSeq-- para que sea de 3'->5'
        if ((a = ascii2bits[*pSeq--]))  {
            if (a == 2) {
                validKey = 0;
                break;
            }
            key |= code1;
        }
        code1 >>= 1;
    }
    return key;
}

/**

 * copiada de otro solo para respaldo
void ogBitwiseAT1GC0Encoding::getFRKeys(char *p, uint32_t &keyF, uint32_t &keyR) {
    keyF = keyR = 0;
    uint32_t    code1 = keyCode;
    uint32_t    code2 = keyCodeRC;
    uint8_t     i;
    for (i=sizeInChars; i; i--, p++) {
        if (ascii2bits[*p])  {
            // taking advantage that code at *p is the same that the complement
            keyF |= code1;
            keyR |= code2;
        }
        code1 >>= 1;
        code2 <<= 1;        
    }
}


void ogBitwiseAT1GC0Encoding::shiftLocalFRKeys(char c) {
    currKeyFwd <<= 1;
    currKeyRev >>= 1;
    currKeyFwdPos++;
    currKeyRevPos++;
    //if (header.ascii2code[*(p+header.tupleSizeInChars_1)]) key |= code1;
    if (ascii2bits[c]) {
        currKeyFwd |= keyCodeRC; // At 3'
        currKeyRev |= keyCode;   // At 5'
    }
    currKeyFwd &= maskCode;
    currKeyRev &= maskCode;
}
**/

/**
uint32_t ogBitwiseAT1GC0Encoding::getKey(char *p) {
    uint32_t    key = 0;
    uint32_t    code = keyCode;
    uint8_t     i;
    for (i=0; i < sizeInChars; i++) {
        if (ascii2bits[*p++]) key |= code;
        code >>= 1;
    }
    
    return key;
}

void ogBitwiseAT1GC0Encoding::shiftKey(char c, uint32_t &key) {
    key <<= 1;
    if (ascii2bits[c]) key |= keyCodeRC; // new code is always at 3'
    key &= maskCode;
}

void ogBitwiseAT1GC0Encoding::getFRKeys(char *p, uint32_t &keyF, uint32_t &keyR) {
    keyF = keyR = 0;
    uint32_t    code1 = keyCode;
    uint32_t    code2 = keyCodeRC;
    uint8_t     i;
    for (i=sizeInChars; i; i--, p++) {
        if (ascii2bits[*p])  {
            // taking advantage that code at *p is the same that the complement
            keyF |= code1;
            keyR |= code2;
        }
        code1 >>= 1;
        code2 <<= 1;        
    }
}

void ogBitwiseAT1GC0Encoding::shiftFRKeys(char c, uint32_t &keyF, uint32_t &keyR) {
    keyF <<= 1;
    keyR >>= 1;
    //if (header.ascii2code[*(p+header.tupleSizeInChars_1)]) key |= code1;
    if (ascii2bits[c]) {
        keyF |= keyCodeRC; // At 3'
        keyR |= keyCode;   // At 5'
    }
    keyF &= maskCode;
    keyR &= maskCode;
}
**/

char ogBitwiseAT1GC0Encoding::isKeyMatchExactSequence() {
    // This determine whether the Key represent the sequence EXACTLY
    return 0;
}
