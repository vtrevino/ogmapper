/* 
 * File:   ogKeyEncoding.cpp
 * Author: victortrevino
 * 
 * Created on May 7, 2022, 2:37 PM
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ogDefinitions.hpp"
#include "ogKeyEncoding.hpp"

ogKeyEncoding::ogKeyEncoding() {
    strncpy(name, "NoEncodingNameSet", MAX_CHAR_NAMES);
    version = 0.0;
}

const char *ogKeyEncoding::getShortExtensionName() {
    return "_NoEnc";
}

ogKeyEncoding::~ogKeyEncoding() {
    //fprintf(stderr, "<Deallocating ogKeyEncoding [%s] version [%f]:>", name, version);
    //if (strlen(name) > 2) fprintf(stderr, "<~ogKE %c%c%c:",name[0],name[1],name[2]); 
    //else fprintf(stderr, "<~ogKE:");
    fflush(stderr);
}

char *ogKeyEncoding::getName() {
    return name;
}

void ogKeyEncoding::setSizeInChars(uint16_t size) {
    sizeInChars = size;
    sizeInBits = 0;
}

uint16_t ogKeyEncoding::getSizeInChars() {
    return sizeInChars;
}

uint8_t ogKeyEncoding::getSizeInBits() {
    return sizeInBits;
}

void ogKeyEncoding::load(FILE *pFile) {
    // do nothing, may be a subclass
}

uint32_t ogKeyEncoding::save(FILE *pFile) {
    // do nothing, may be a subclass
    return 0;
}

uint32_t ogKeyEncoding::getFwdKey(char *pSeq) {
    validKey = 0;
    return 0;
}

uint32_t ogKeyEncoding::getRevKey(char *pSeq) {
    validKey = 0;
    return 0;
}

char ogKeyEncoding::isKeyMatchExactSequence() {
    // This determine whether the Key represent the sequence EXACTLY
    return (sizeInBits == sizeInChars*2);
}

char ogKeyEncoding::isValidKey() {
    return validKey;  // 0 = not valid
}

uint32_t ogKeyEncoding::getTotalKeys() {
    return pow(2, getSizeInBits());
}

char ogKeyEncoding::canKeyBeComparedInSequence() {
    return 1;
}
