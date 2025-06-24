/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.h to edit this template
 */

/* 
 * File:   ogFastAQGZreader.hpp
 * Author: victortrevino
 *
 * Created on November 30, 2024, 5:45 PM
 */

#ifndef OGFASTAQGZREADER_HPP
#define OGFASTAQGZREADER_HPP

#include <stdio.h>
#include <time.h>
#include <string.h>
#include <zlib.h>
//#include <thread>
//#include <chrono>
//#include <sstream>
#include <iostream>

typedef struct OG_SEQ_STRING {
	uint32_t l;
	char *s;
} ogSeqString;

typedef struct OG_FASTAQ_BUFFER {
    char            *buffer;
    uint32_t        buffSize;
    uint32_t        buffNum;
    OG_FASTAQ_BUFFER  *nextBuffer;
    int             *readsDelivered; // array of total reads delivered
    int             *readsProcessed;   // array of reads pending to be released
} ogFastAQBuffer;

typedef struct OG_FASTAQ_SEQUENCE {
    ogSeqString     name, comment, seq, qual;
    int            *pRelease;
    ogFastAQBuffer *pBuffer;
} ogFastAQsequence;

class ogFastAQGZreader {
    
    char            *buffer;
    uint32_t        buffSize;
    uint32_t        mark;       // begin of the current sequence
    uint32_t        current;
    uint32_t        end;
    char            eof;
    char            eofReach;
    gzFile          *file;
    char            readIdx;       // read 1 or 2
    uint32_t        heap;       // position to write to
    uint32_t        copyFrom, copyTo;
    uint16_t        threads;
    uint64_t        filePos;
    uint64_t        fileSize;
    ogFastAQBuffer  *fasBuffer;
    ogFastAQBuffer  *curBuffer;
    ogFastAQBuffer  *lastBuffer;
    //uint32_t        moveFrom;
    //uint32_t        moveTo;
    
    
    
    
public:
    ogFastAQGZreader();
    ogFastAQGZreader(const ogFastAQGZreader& orig);
    
    void                setFileReader(gzFile *gzF, char xRead, uint16_t nThreads);
    char                isEOF();
    void                read();
    char                get1char();
    void                consumeForNewLines(char c);
    void                reconfigureFAS(ogFastAQsequence *fas, uint32_t nam2m, uint32_t com2m, uint32_t qual2m, uint32_t seq2m);
    void                checkToMove();
    char                getNextSequenceOriginal(ogFastAQsequence *seq); 
    char                getNextSequenceState(ogFastAQsequence *fas);
    char                getNextSequenceNaive(ogFastAQsequence *fas);
    char                getNextSequenceByType(ogFastAQsequence *fas);
    void                cleanBufferAndRead();
    void                setFileSize(uint64_t size);
    float               estimatedFilePosition();

    char                __ReadWhileTextCharsNoSpace();
    char                __ReadUntilNewLines();
    char                __ReadWhileNewLines();

    //static void         assignReadToThread(ogFastAQBuffer *pBuf, uint16_t iThread);
    //static void         releaseReadFromThread(ogFastAQBuffer *pBuf, uint16_t iThread);

    virtual ~ogFastAQGZreader();

    
    static int  *assignReadToThread(ogFastAQBuffer *pBuf, uint16_t iThread) {
        //if (threads > 0) {
            pBuf->readsDelivered[iThread]++;
            int *p = pBuf->readsProcessed + iThread;
            //(*p)++; // ya no incrementa aquí sino que incrementa en el release
            return p;
        //}
    }

    //static void  releaseReadFromThread(ogFastAQBuffer *pBuf, uint16_t iThread) {
        //if (threads > 0) {
    //        pBuf->readsProcessed[iThread]--;
        //}
    //}

private:

};

#endif /* OGFASTAQGZREADER_HPP */

