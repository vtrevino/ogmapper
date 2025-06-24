/* 
 * File:   ogTupleGuider.cpp
 * Author: victortrevino
 * 
 * Created on July 13, 2022, 1:52 AM
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include "ogDefinitions.hpp"
#include "ogTupleGuider.hpp"
#include <string.h>

ogTupleGuider::ogTupleGuider() {
    strcpy(extName, "_gTupl");
    memset(ascii2n, 4, sizeof(ascii2n));
    ascii2n['A'] = ascii2n['a'] = 0;
    ascii2n['C'] = ascii2n['c'] = 1;
    ascii2n['G'] = ascii2n['g'] = 2;
    ascii2n['T'] = ascii2n['t'] = 3;
    ascii2n['U'] = ascii2n['u'] = 3;
    ascii2n['N'] = ascii2n['n'] = 4; // ERROR
    tupleLen = 1;
    tupleMask = 0;
    genomeMode = 0;
}

ogTupleGuider::~ogTupleGuider() {
    //fprintf(stderr, "<Deallocating guider [ogTupleGuider]:");
    //fprintf(stderr, "<~ogTG:"); fflush(stderr);
    if (isGuider != NULL) free(isGuider);
    if (ntAdvances != NULL) free(ntAdvances);
    //fprintf(stderr, ":ogTG>");
}

ogTupleGuider *ogTupleGuider::clone() {
    ogTupleGuider *ogG = new ogTupleGuider(); 
    ogG->isSymmetric = isSymmetric;
    ogG->fixedGuideLength = fixedGuideLength;
    ogG->keySizeInChars = keySizeInChars;
    ogG->maxLenNoPattern = maxLenNoPattern;
    ogG->tupleLen = tupleLen;
    ogG->tupleMask = tupleMask;
    ogG->fixedStep = fixedStep;
    strncpy(ogG->extName, extName, sizeof(extName));
    strncpy(ogG->tupleFileName, tupleFileName, sizeof(tupleFileName));
    ogG->isGuider = (char *) calloc(tupleLen, sizeof(char));
    memmove(ogG->isGuider, isGuider, tupleLen * sizeof(char));
    ogG->ntAdvances = (char *) calloc(tupleLen, sizeof(char));
    memmove(ogG->ntAdvances, ntAdvances, tupleLen * sizeof(char));
    return ogG;
}

void ogTupleGuider::initializeForIndexing() {
    
    FILE *pFile;
    char line[1000];
    char comments[1000];
    char tuple[100];
    char        x;
    uint16_t    a, b;
    uint16_t    t, i;
    
    strncpy(tupleFileName, configurationFile, MAX_FILENAME);
    // FILE CONTENT FORMAT:
    // # Optional Comments, as many lines as needed
    // # ...
    // Extension: 4 letters (not used as only 1 file is now in index)
    // MAX_LEN_NO_PATTERN : integer (supposed to be short, around 24)
    // isSymmetric : 0=no, 1=yes
    // guideLen : size of the guide in nucleotides, integer. All guideLen-combinations are expected
    // # Optional Comments, as many lines as needed
    // # ...
    // tuple isGuider ntadvance comment
    // ..
    
    //fprintf(stderr, "$$$$ Accessing Tuple control file : %s\n", OGTUPLE_GUIDER_CONTROL_FILE);
    //pFile = fopen(OGTUPLE_GUIDER_CONTROL_FILE, "r");
    fprintf(stderr, "$$$$ Accessing Tuple file : %s\n", tupleFileName);
    pFile = fopen(tupleFileName, "r");
    // next line is "name", at most, 4 chars are used
    fscanf(pFile, "%[^\n]s", line); 
    while (line[0] == '#') {
        fscanf(pFile,"%c",&x);
        fscanf(pFile, "%[^\n]s", line); 
    }
    
    strcpy(extName, "_gTpl");
    strncat(extName, line, 4);
    fscanf(pFile,"%c",&x);
     
    // removed: second line is the path to the state machine file
    // fscanf(pFile, "%[^\n]s", tupleFileName); 

    // next line is maximum length of no pattern
    maxLenNoPattern = 24;
    fscanf(pFile, "%hu", &maxLenNoPattern); 
    // next line is symmetric pattern =1 or =0 no symmetric
    fscanf(pFile, "%hhu", &isSymmetric);
    // next line is pattern length (fixed)
    fscanf(pFile, "%hhu", &fixedGuideLength);
    // fclose(pFile);
    fixedStep = maxLenNoPattern / (fixedGuideLength < 1 ? 1 : fixedGuideLength);
    if (fixedStep < 1) fixedStep = 1;
    
    fprintf(stderr, "%%%%%%%% Tuple Parameters %%%%%%%%\n");
    fprintf(stderr, "Name         : %s\n", line);
    fprintf(stderr, "Extension    : %s\n", extName);
    fprintf(stderr, "Tuple File   : %s\n", tupleFileName);
    fprintf(stderr, "Max Len      : %hu (maximum length without a guide)\n", maxLenNoPattern);
    fprintf(stderr, "Is Symmetric : %hhu\n", isSymmetric);
    fprintf(stderr, "Fixed Length : %hhu\n", fixedGuideLength);
    fprintf(stderr, "Fixed Step   : %hu\n", fixedStep);
    tupleLen = pow(4, fixedGuideLength);
    tupleMask = tupleLen-1;
    fprintf(stderr, "Tuple comb n : %hd\n", tupleLen);    
    fprintf(stderr, "Allocating %lu bytes of memory for tuples.\n", (tupleLen*2) * sizeof(char));
    isGuider = (char *) calloc(tupleLen, sizeof(char));
    ntAdvances = (char *) calloc(tupleLen, sizeof(char));
    memset(isGuider, 0, tupleLen);
    memset(ntAdvances, 1, tupleLen);

    // READING TUPLE DATA
    //fprintf(stderr, "$$$$ Accessing Tuple file : %s\n", tupleFileName);
    //pFile = fopen(tupleFileName, "r");
    for (; fscanf(pFile, "%[^\n]s", line) != EOF; ) {
        if (line[0] != '#' && line[0] != '\n' && line[0] != '\r') {
            if (line[0] == 0) break;
            sscanf(line, "%s %hu %hu %[^\n]s", (char *)tuple, &a, &b, (char *)comments);
            if (strlen(tuple) == fixedGuideLength && (a == 0 || a == 1) && (b > 0)) {
                for (t=0, i=0; i < strlen(tuple); i++) {
                    t <<= 2;
                    t |= ascii2n[tuple[i]];
                }
                isGuider[t] = a;
                ntAdvances[t] = b;
                //fprintf(stderr, "DNA=%s Tuple=%u Guider=%hhu NtAdvance=%hhu Comment=%s\n", tuple, t, isGuider[t], ntAdvances[t], comments);
                if (t != 0 && t % 4 == 0) fprintf(stderr, "\n");
                fprintf(stderr, "[%u]:%s:%hhu:+%hhu|", t, tuple, isGuider[t], ntAdvances[t]);
            } else {
                fprintf(stderr, "line ignored:%hhu\n[%s]\n", *line, line);
                fprintf(stderr, "Tuple Len=%zu, isGuider=%d, ntAdv=%d\n", strlen(tuple), a, b);
            }
        } else {
            //fprintf(stderr, "%s\n", line);
        }
        fscanf(pFile,"%c",&x); // Enter ?
        line[0] = 0;
    }    
    fclose(pFile);
    fprintf(stderr, "\n");
}


uint32_t ogTupleGuider::save(FILE *pFile) {
    uint32_t bytes = 0;
    //char xfile[1000];
    //setFileName(xfile, pFileName, getKeySizeInChars(), "ogt", 1000);
    //FILE *pFile = fopen(xfile, "w");
  //fwrite("StatesMachineGuider (this is a binary file)      \n", sizeof(char), 50, pFile);
    bytes += sizeof(char)*fwrite("TupleGuider         (this is a binary file)      \n", sizeof(char), 50, pFile);
    bytes += fprintf(pFile, "%s\n",tupleFileName);
    bytes += fwrite(&tupleLen, 1, sizeof(tupleLen), pFile);
    bytes += fwrite(&maxLenNoPattern, 1, sizeof(maxLenNoPattern), pFile);    
    bytes += fwrite(&isSymmetric, 1, sizeof(isSymmetric), pFile);    
    bytes += fwrite(&fixedGuideLength, 1, sizeof(fixedGuideLength), pFile);    
    bytes += fwrite(extName, 1, sizeof(extName), pFile);
    bytes += sizeof(char)*fwrite(isGuider, sizeof(char), tupleLen, pFile);
    bytes += sizeof(char)*fwrite(ntAdvances, sizeof(char), tupleLen, pFile);
    //fclose(pFile);
    return bytes;
}

uint32_t ogTupleGuider::load(FILE *pFile) {
    uint32_t bytes = 0;
    uint32_t nbyname = 0;
    //char xfile[1000];
    char dummyname[100];
    char x;
    //setPureFileName(xfile, pFileName, "ogt", 1000);
    //fprintf(stderr, "%%%%%%%% Accessing Guider File.\n");
    //FILE *pFile = fopen(xfile, "r");
    ////fread(dummyname, sizeof(char), 50, pFile);
    bytes += fscanf(pFile, "%[^\n]s", dummyname); bytes += fscanf(pFile,"%c",&x); // Comment
    bytes += fscanf(pFile, "%[^\n]s", tupleFileName); bytes += fscanf(pFile,"%c",&x); // Filename
    bytes += strlen(dummyname);
    bytes += strlen(tupleFileName);
    strncpy(configurationFile, tupleFileName, MAX_FILENAME);
    bytes += fread(&tupleLen, sizeof(tupleLen), 1, pFile);
    bytes += fread(&maxLenNoPattern, sizeof(maxLenNoPattern), 1, pFile);
    bytes += fread(&isSymmetric, sizeof(isSymmetric), 1, pFile);
    bytes += fread(&fixedGuideLength, sizeof(fixedGuideLength), 1, pFile);
    bytes += fread(extName, sizeof(extName), 1, pFile);
    fprintf(stderr, "Allocating %lu bytes of memory for tuple guides.\n", (tupleLen*2) * sizeof(char));
    isGuider = (char *) calloc(tupleLen, sizeof(char));    
    ntAdvances = (char *) calloc(tupleLen, sizeof(char));
    tupleMask = tupleLen - 1;
    bytes += fread(isGuider, sizeof(char), tupleLen, pFile);
    bytes += fread(ntAdvances, sizeof(char), tupleLen, pFile);
    fixedStep = maxLenNoPattern / (fixedGuideLength < 1 ? 1 : fixedGuideLength);
    if (fixedStep < 1) fixedStep = 1;
    //fclose(pFile);
    /**
    fprintf(stderr, "%% Tuple Parameters ");
    fprintf(stderr, "{ Ext: '%s',", extName);
    fprintf(stderr, " File: '%s',", tupleFileName);
    fprintf(stderr, " MaxLenNoGuide: %hu,", maxLenNoPattern);
    fprintf(stderr, " IsSymm: %hhu,", isSymmetric);
    fprintf(stderr, " FixLen: %hhu,", fixedGuideLength);
    fprintf(stderr, " FixStep: %hu,", fixedStep);
    fprintf(stderr, " Comb n: %hd }\n", tupleLen);    
    int t;
    //fprintf(stderr, "[Tuple, Guider, Nt Advance]:");
    //for (t=0; t < tupleLen; t++) {
    //    fprintf(stderr, "[%u,%hhu,%hhu]  ", t, isGuider[t], ntAdvances[t]);
    //    if ((t+1) % 8 == 0 && (t+1 < tupleLen)) fprintf(stderr, "\n[Tuple, Guider, Nt Advance]:");
    //}
    fprintf(stderr, "%% Guides:");
    for (t=0; t < tupleLen; t++) {
        if (isGuider[t]) fprintf(stderr, "%u(+%hhu),", t, ntAdvances[t]);
    }
    //fprintf(stderr, "\n%% No-Guides:");
    //for (t=0; t < tupleLen; t++) {
    //    if (isGuider[t] == 0) fprintf(stderr, "%u(+%hhu),", t, ntAdvances[t]);
    //}
    fprintf(stderr, "\n");
    **/
    return bytes;
}

char ogTupleGuider::nextGuide() {
    if (finish) return 0;
    while (ascii2n[*pI] == 4 && *pI && pI < pSeqLimit) {
        pI++;
    }
    uint16_t maxLen = maxLenNoPattern;
    char *pStart = pI;
    char *pJ = pI;
    uint16_t len, i, t;
    char a;
    while (finish==0) {
        for (t=len=i=0; i < fixedGuideLength && *pJ && len < maxLen; ) {
            t <<= 2;
            t |= (a=ascii2n[*pJ++]);
            if (a==4) { // N or any unrecognized char
                t=len=i=0;
            } else {
                len++, i++; // 
            }
        }
        while (*pJ) {
            if (isGuider[t]) {
                pI = pJ - fixedGuideLength + ntAdvances[t];
                if (pI > pSeqLimit) {
                    finish = 1;
                    return 0;
                }
                pCurrent = pI;
                count++;
                return 1;
            }
            if (len >= maxLen) {
                pI = pStart + fixedStep; // (genomeMode ? 1 : fixedGuideLength); // antes estaba en +1 // ESTE DEBERIA SER +1 cuando es Genomic y fixGuideLen cuando es Read
                if (pI > pSeqLimit) {
                    finish = 1;
                    return 0;
                }
                count++;
                countOverLen++;
                pCurrent = pI;
                return 1;
            }
            for (i=ntAdvances[t]; i && *pJ; i--) {
                t = (t << 2) & tupleMask | (a=ascii2n[*pJ++]);
                len++;
                if (a == 4) break; // N or any unrecognized char
            }
            if (a==4) break; // N or any unrecognized char
        }
        if (*pJ == 0) finish = 1;
    }
    return 0;
}

const char myname[] = "TupleGuider";

const char * ogTupleGuider::getName() {
    return myname;
}

const char * ogTupleGuider::getShortExtensionName() {
    return extName;
}
