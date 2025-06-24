/* 
 * File:   ogStateMachineGuider.cpp
 * Author: victortrevino
 * 
 * Created on July 1, 2022, 8:10 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "ogDefinitions.hpp"
#include "ogStateMachineGuider.hpp"
#include <string.h>

ogStateMachineGuider::ogStateMachineGuider() {
    strcpy(extName, "_gStMch");
    //int i;
    //for (i=0; i < sizeof(ascii2col); i++) {
    //    fprintf(stderr, "%d:%hhu ", i, ascii2col[i]);
    //}
    //fprintf(stderr, "\n");
    memset(ascii2col, 5, sizeof(ascii2col));
    ascii2col['A'] = ascii2col['a'] = 0;
    ascii2col['C'] = ascii2col['c'] = 1;
    ascii2col['G'] = ascii2col['g'] = 2;
    ascii2col['T'] = ascii2col['t'] = 3;
    ascii2col['N'] = ascii2col['n'] = 4;    
    //for (i=0; i < sizeof(ascii2col); i++) {
    //    fprintf(stderr, "%d:%hhu ", i, ascii2col[i]);
    //}
    //fprintf(stderr, "\n");
    pTransitionMatrix = NULL;
}

ogStateMachineGuider::~ogStateMachineGuider() {
    //fprintf(stderr, "<Deallocating guider [ogStateMachineGuider]:");
    //fprintf(stderr, "<~ogSMG:"); fflush(stderr);
    if (pTransitionMatrix != NULL) free(pTransitionMatrix);
    //fprintf(stderr, ":ogSMG>");
}

ogStateMachineGuider *ogStateMachineGuider::clone() {
    //fprintf(stderr, "ogStateMachineGuider.clone()\n");
    ogStateMachineGuider *ogG = new ogStateMachineGuider();
    ogG->isSymmetric = isSymmetric;
    ogG->fixedGuideLength = fixedGuideLength;
    ogG->keySizeInChars = keySizeInChars;
    ogG->maxState = maxState;
    ogG->maxLenNoPattern = maxLenNoPattern;
    ogG->fixedStep = fixedStep;
    strncpy(ogG->extName, extName, sizeof(extName));
    strncpy(ogG->stateFileName, stateFileName, sizeof(stateFileName));
    //strncpy(ogG->ascii2col, ascii2col, sizeof(ascii2col));
    ogG->pTransitionMatrix = (TransitionVector *) calloc(maxState+1, sizeof(TransitionVector));
    memmove(ogG->pTransitionMatrix, pTransitionMatrix, (maxState+1) * sizeof(TransitionVector));
    return ogG;
}

void ogStateMachineGuider::initializeForIndexing() {
    
    FILE *pFile;
    char line[1000];
    char comments[1000];
    char        x;
    
    strncpy(stateFileName, configurationFile, MAX_FILENAME);
    // FILE CONTENT FORMAT:
    // # Optional Comments, as many lines as needed
    // # ...
    // Extension: 4 letters (not used as only 1 file is now in index)
    // MAX_LEN_NO_PATTERN : integer (supposed to be short, around 24)
    // isSymmetric : 0=no, 1=yes
    // guideLen : size of the guide in nucleotides, integer. All guideLen-combinations are expected
    // # Optional Comments, as many lines as needed
    // # ...
    // State     A       C       G       T       N       X       Comment
    // ..

    // NO Control-File ANY MORE
    //fprintf(stderr, "$$$$ Accessing State machine control file : %s\n", OGSTATE_MACHINE_CONTROL_FILE);
    //pFile = fopen(OGSTATE_MACHINE_CONTROL_FILE, "r");
    // first line is "name", at most, 4 chars are used
    fprintf(stderr, "$$$$ Accessing State machine control file : %s\n", stateFileName);
    pFile = fopen(stateFileName, "r+b");
    fscanf(pFile, "%[^\n]s", line); 
    while (line[0] == '#') {
        fscanf(pFile,"%c",&x);
        fscanf(pFile, "%[^\n]s", line); 
    }
    
    strcpy(extName, "_gSM");
    strncat(extName, line, 4);
    fscanf(pFile,"%c",&x);
     
    // second line is the path to the state machine file
    // fscanf(pFile, "%[^\n]s", stateFileName); 

    // third line is maximum length of no pattern
    maxLenNoPattern = 24;
    fscanf(pFile, "%hu", &maxLenNoPattern); 
    // 4th line is symmetric pattern =1 or =0 no symmetric
    fscanf(pFile, "%hhu", &isSymmetric);
    // 5th line is pattern length (fixed)
    fscanf(pFile, "%hhu", &fixedGuideLength);
    // fclose(pFile);
    fixedStep = maxLenNoPattern / (fixedGuideLength < 1 ? 1 : fixedGuideLength);
    if (fixedStep < 1) fixedStep = 1;
    fprintf(stderr, "%%%% State Machine Parameters %%%%\n");
    fprintf(stderr, "Name         : %s\n", line);
    fprintf(stderr, "Extension    : %s\n", extName);
    fprintf(stderr, "StMach File  : %s\n", stateFileName);
    fprintf(stderr, "Max Len      : %hu (maximum length without a guide)\n", maxLenNoPattern);
    fprintf(stderr, "Is Symmetric : %hhu\n", isSymmetric);
    fprintf(stderr, "Fixed Length : %hhu\n", fixedGuideLength);
    fprintf(stderr, "Fixed Step   : %hu\n", fixedStep);
    
    // STATE FILE, CHECKING MAXIMUM STATE NUMBER
    maxState = 100;
    // pFile = fopen(stateFileName, "r");
    int16_t    i, iState;
            
    pTransitionMatrix = (TransitionVector *) calloc(maxState+1, sizeof(TransitionVector));
    // READING STATES FROM FILE
    //pFile = fopen(stateFileName, "r");
    for (; fscanf(pFile, "%[^\n]s", line) != EOF; ) {
        if (line[0] != '#') {
            if (line[0] == 0) break;
            sscanf(line, "%hd ", &i);
            if (i >= 0) {
                if (i > maxState) {
                    maxState = i + 10;
                    pTransitionMatrix = (TransitionVector *) realloc(pTransitionMatrix, (maxState+1)*sizeof(TransitionVector));
                }
                TransitionVector *pTV = &pTransitionMatrix[i];
                sscanf(line, "%hd %hd %hd %hd %hd %hd %hd %[^\n]s", &iState, &pTransitionMatrix[i][0], &pTransitionMatrix[i][1], &pTransitionMatrix[i][2], &pTransitionMatrix[i][3], &pTransitionMatrix[i][4], &pTransitionMatrix[i][5], (char *) comments);
                fprintf(stderr, "St:%hd\tA:%hd\tC:%hd\tG:%hd\tT:%hd\tN:%hd\t*:%hd\t%s\n", iState, pTransitionMatrix[i][0], pTransitionMatrix[i][1], pTransitionMatrix[i][2], pTransitionMatrix[i][3], pTransitionMatrix[i][4], pTransitionMatrix[i][5], comments);
            } else {
                fprintf(stderr, "Negative state, line ignored:\n%s\n", line);                
            }
        } else {
            fprintf(stderr, "%s\n", line);
        }
        fscanf(pFile,"%c",&x); // Enter ?
        line[0] = 0;
    }    
    fprintf(stderr, "Max State    : %hu\n", maxState);
    fprintf(stderr, "Allocated %lu bytes of memory for State Machine.\n", (maxState+1) * sizeof(TransitionVector));
    fclose(pFile);    
}

uint32_t ogStateMachineGuider::save(FILE *pFile) {
    uint32_t bytes = 0;
    //char xfile[1000];
    //setFileName(xfile, pFileName, getKeySizeInChars(), "ogm", 1000);
    //FILE *pFile = fopen(xfile, "w");
    bytes += sizeof(char)*fwrite("StatesMachineGuider (this is a binary file)      \n", sizeof(char), 50, pFile);
    bytes += fprintf(pFile, "%s\n",stateFileName);
    bytes += fwrite(&maxState, 1, sizeof(maxState), pFile);
    bytes += fwrite(&maxLenNoPattern, 1, sizeof(maxLenNoPattern), pFile);    
    bytes += fwrite(&isSymmetric, 1, sizeof(isSymmetric), pFile);    
    bytes += fwrite(&fixedGuideLength, 1, sizeof(fixedGuideLength), pFile);    
    bytes += fwrite(extName, 1, sizeof(extName), pFile);
    bytes += sizeof(TransitionVector) * fwrite(pTransitionMatrix, sizeof(TransitionVector), maxState+1, pFile);
    //fclose(pFile);
    return bytes;
}

uint32_t ogStateMachineGuider::load(FILE *pFile) {
    //char xfile[1000];
    char dummyname[100];
    uint32_t bytes = 0;
    char x;
    //setPureFileName(xfile, pFileName, "ogm", 1000);
    fprintf(stderr, "%%%%%%%% Accessing File [%s].\n.", "");
    //FILE *pFile = fopen(xfile, "r");
    ////fread(dummyname, sizeof(char), 50, pFile);
    bytes += fscanf(pFile, "%[^\n]s", dummyname); bytes += fscanf(pFile,"%c",&x); // Comment
    bytes += fscanf(pFile, "%[^\n]s", stateFileName); bytes += fscanf(pFile,"%c",&x); // Filename
    bytes += strlen(dummyname);
    bytes += strlen(stateFileName);
    strncpy(configurationFile, stateFileName, MAX_FILENAME);    
    bytes += fread(&maxState, sizeof(maxState), 1, pFile);
    bytes += fread(&maxLenNoPattern, sizeof(maxLenNoPattern), 1, pFile);
    bytes += fread(&isSymmetric, sizeof(isSymmetric), 1, pFile);
    bytes += fread(&fixedGuideLength, sizeof(fixedGuideLength), 1, pFile);
    bytes += fread(extName, sizeof(extName), 1, pFile);
    fixedStep = maxLenNoPattern / (fixedGuideLength < 1 ? 1 : fixedGuideLength);
    if (fixedStep < 1) fixedStep = 1;
    fprintf(stderr, "Allocating %lu bytes of memory for State Machine.\n", (maxState+1) * sizeof(TransitionVector));
    pTransitionMatrix = (TransitionVector *) calloc(maxState+1, sizeof(TransitionVector));    
    bytes += fread(pTransitionMatrix, sizeof(TransitionVector), maxState+1, pFile);
    //fclose(pFile);
    fprintf(stderr, "%%%% State Machine Parameters %%%%\n");
    fprintf(stderr, "Extension    : %s\n", extName);
    fprintf(stderr, "StMach File  : %s\n", stateFileName);
    fprintf(stderr, "Max State    : %hu\n", maxState);
    fprintf(stderr, "Max Len      : %hu (maximum length without a guide)\n", maxLenNoPattern);
    fprintf(stderr, "Is Symmetric : %hhu\n", isSymmetric);
    fprintf(stderr, "Fixed Length : %hhu\n", fixedGuideLength);
    fprintf(stderr, "Fixed Step   : %hu\n", fixedStep);
    int i, j, k;
    for (i=0; i <= maxState; i++) {
        TransitionVector *pTV = &pTransitionMatrix[i];
        j = 0;
        for (k=0; k < 6; k++) { if (pTransitionMatrix[i][k]) { j = 1; break; } }
        if (j) {
            fprintf(stderr, "St:%d\tA:%hd\tC:%hd\tG:%hd\tT:%hd\tN:%hd\t*:%hd\n", i, pTransitionMatrix[i][0], pTransitionMatrix[i][1], pTransitionMatrix[i][2], pTransitionMatrix[i][3], pTransitionMatrix[i][4], pTransitionMatrix[i][5]);
        }
    }
    return bytes;
}

char ogStateMachineGuider::nextGuide() {
    if (finish) return 0;
    uint16_t maxLen = maxLenNoPattern;
    char *pStart = pI;
    if (pStart + maxLen > pSeqLimit) maxLen = pSeqLimit - pStart + 1;
    while (true) {
        char *pJ = pI;
        int16_t  state = 0;
        uint16_t len = pJ - pStart;
        //printf("State-Machine Start:[");
        for (; *pJ && len < maxLen; len++) {
            //printf("%c,%d| ",*pJ, state);
            if ((state = pTransitionMatrix[state][ascii2col[*pJ]]) < 0) break;
            pJ++;
        }
        //printf("], State:%d\n",state);
        char x;
        //scanf("%c", &x);
        if (state == -1) {
            // Terminó correctamente, pointer actual inicia el key encoding
            count++;
            pCurrent = pI = pJ;
            return 1;
        } else if (state == -22) {
            // Terminó en un estado fallido, solo incrementar pointer y continuar.
            if (++pI > pSeqLimit) {
                finish = 1;
                return 0;
            }
        } else {
            // Terminó por longitud maxLenNoPattern
            if (len >= maxLenNoPattern) {
                pCurrent = pI = pStart + fixedStep; // (genomeMode ? 1 : fixedGuideLength); // ESTE DEBERIA SER +1 cuando es Genomic y fixGuideLen cuando es Read
                if (pI > pSeqLimit) {
                    finish = 1;
                    return 0;
                }
                count++;
                countOverLen++;
                return 1;
            } else {
                // si llegó aquí es que pasó maxLen pero no maxLenNoPattern y por lo tanto sobrepaso el limite
                finish = 1;
                return 0;
            }
            // Otro estado desconocido, ignorar si no ha terminado
            if (*pJ == 0 || ++pI > pSeqLimit) {
                finish = 1;
                return 0;
            }            
        }
    }
    // In theory this should not be run:
    pCurrent = pI;
    count++;
    return 1;    
}

const char myname[] = "StateMachineGuider";

const char * ogStateMachineGuider::getName() {
    return myname;
}

const char * ogStateMachineGuider::getShortExtensionName() {
    return extName;
}


// ACACACA
//    1111
//       2
//  rrrrrr
//     111