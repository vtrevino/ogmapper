/* 
 * File:   ogSwapBitwiseAT1GC0Encoding.cpp
 * Author: victortrevino
 * 
 * Created on July 5, 2022, 6:38 PM
 */
#include <stdio.h>
#include <string.h>
#include "ogDefinitions.hpp"

#include "ogSwapBitwiseAT1GC0Encoding.hpp"

const char myname[] = "SwapBitwiseAT1GC0Encoding\0";


ogSwapBitwiseAT1GC0Encoding::ogSwapBitwiseAT1GC0Encoding() {
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

ogSwapBitwiseAT1GC0Encoding::~ogSwapBitwiseAT1GC0Encoding() {
}


const char *ogSwapBitwiseAT1GC0Encoding::getShortExtensionName() {
    return "_eSwpBW10";
}


void ogSwapBitwiseAT1GC0Encoding::setSizeInChars(uint16_t size) {
    sizeInChars = size;
    sizeInBits = sizeInChars / 2;
    keyCode = 0x00000001 << (sizeInBits - 1);
    keyCodeRC = 0x00000001;
    maskCode = 0xFFFFFFFF >> (32 - sizeInBits);    
}

uint32_t ogSwapBitwiseAT1GC0Encoding::getFwdKey(char *pSeq) {
    uint32_t key = 0;
    uint32_t code1 = keyCode;
    uint8_t  i;
    char     a;
    validKey = 1;
    for (i=sizeInBits; i; i--, pSeq += 2) {
        if ((a=ascii2bits[*pSeq]))  {
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

uint32_t ogSwapBitwiseAT1GC0Encoding::getRevKey(char *pSeq) {
    // Ahora este encoding es igual porque la secuencia se recorre de 3'->5'
    uint32_t key = 0;
    uint32_t code1 = keyCode;
    uint8_t  i;
    char     a;
    validKey = 1;
    for (i=sizeInBits; i; i--, pSeq -= 2) { // pSeq-- para que sea de 3'->5'
        if ((a=ascii2bits[*pSeq]))  {
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


char ogSwapBitwiseAT1GC0Encoding::isKeyMatchExactSequence() {
    // This determine whether the Key represent the sequence EXACTLY
    return 0;
}

char ogSwapBitwiseAT1GC0Encoding::canKeyBeComparedInSequence() {
    // This determine whether the Key represent the sequence EXACTLY
    return 0;
}
