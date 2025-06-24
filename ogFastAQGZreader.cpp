/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.cc to edit this template
 */

/* 
 * File:   ogFastAQGZreader.cpp
 * Author: victortrevino
 * 
 * Created on November 30, 2024, 5:45 PM
 */

#include "ogFastAQGZreader.hpp"

ogFastAQGZreader::ogFastAQGZreader() {
    file = NULL;
    mark = current = end = heap = copyFrom = copyTo = 0; // moveFrom = moveTo =
    eof = eofReach = 0;
    readIdx = 0;
    buffSize = 4 * 1024 * 1024;
    buffer = (char *) malloc( buffSize );
    fasBuffer = (ogFastAQBuffer *) malloc(sizeof(ogFastAQBuffer));
    fasBuffer->buffer = buffer;
    fasBuffer->nextBuffer = NULL;
    fasBuffer->buffSize = buffSize;
    fasBuffer->readsProcessed = NULL;
    fasBuffer->readsDelivered = NULL;
    fasBuffer->buffNum = 1;
    curBuffer = lastBuffer = fasBuffer;
    threads = 0;
    fileSize = 1;
    filePos = 0;
}

ogFastAQGZreader::ogFastAQGZreader(const ogFastAQGZreader& orig) {
}

ogFastAQGZreader::~ogFastAQGZreader() {
    ogFastAQBuffer *fbuf = fasBuffer, *nextBuf;
    do {
        //if (threads > 0) 
        //fprintf(stderr, "-- Freeing --\n"); fflush(stderr);
        free(fbuf->readsProcessed);
        free(fbuf->buffer);
        nextBuf = fbuf->nextBuffer;
        free(fbuf);
        fbuf = nextBuf;
        //fprintf(stderr, "-- Done --\n");fflush(stderr);
    } while (fbuf != NULL);
}

void ogFastAQGZreader::setFileSize(uint64_t size) {
    fileSize = size;
}

float ogFastAQGZreader::estimatedFilePosition()  {
    return ((float) filePos * 100 / (float) fileSize);
}

void ogFastAQGZreader::read() {
    if (eofReach != 0) { eof = 1; return; }
    if (end >= buffSize) return; // eof != 0 || 
    int tor = buffSize-end;
    int n = gzread(*file, buffer+end, tor);
    if (n < tor) {
        eofReach = 1;
        if (n == 0) eof = 1;
    }
    end += n;
    filePos = gzoffset(*file);
////
    
    /**
    uint32_t ngt = 0, nlf = 0, nff = 0, ncr = 0, no = 0, nat= 0, i;
    char c;
    for (i=end - n; i < end; i++) {
        c = buffer[i];
        switch (c) {
            case '>': ngt++; break;
            case '\n': nlf++; break;
            case '\f': nff++; break;
            case '\r': ncr++; break;
            case '@': nat++; break;
            default :
                if (c < 32) no++;
        }
        if (c == '>') ngt++;
    }
    fprintf(stderr, "::gzread:: tor=%d, n=%d, eof=%c, >=%d, @=%d, \\n=%d, \\f=%d, \\r=%d, \\o=%d\n", tor, n, eof+48, ngt, nat, nlf, nff, ncr, no);
    **/
    
////
}

void ogFastAQGZreader::setFileReader(gzFile *gzF, char xRead, uint16_t nThreads) {
    file = gzF;
    eof = eofReach = 0;
    readIdx = xRead;
    threads = nThreads;
    fasBuffer->readsProcessed = (int *) calloc((threads+1)*2, sizeof(int));
    fasBuffer->readsDelivered = fasBuffer->readsProcessed + threads + 1; // (int *) calloc((threads+1)*1, sizeof(int)); // fasBuffer->readsPending + threads + 1;
    read();
}

char ogFastAQGZreader::isEOF() {
    return eof == 1;
}

char inline ogFastAQGZreader::get1char() {
    if (current < end) return (buffer[current++]);

    cleanBufferAndRead();
    
    return buffer[current++];
}

void inline ogFastAQGZreader::consumeForNewLines(char c) {
    // '\r' || c == '\f' || c == '\n'
    if (eof == 0 && c <= 32) {
        while (eof == 0 && c <= 32) c = get1char();
        current--;
    }
    checkToMove();
    
}

void inline ogFastAQGZreader::reconfigureFAS(ogFastAQsequence *fas, uint32_t nam2m, uint32_t com2m, uint32_t seq2m, uint32_t qual2m) {

    if (fas->name.l > 0) {
        fas->name.s = buffer + mark + nam2m;
        buffer[mark+nam2m+fas->name.l] = 0;
    }
    if (fas->comment.l > 0) {
        fas->comment.s = buffer + mark + com2m;
        buffer[mark+com2m+fas->comment.l] = 0;
    }
    if (fas->seq.l > 0) {
        fas->seq.s = buffer + mark + seq2m;
        buffer[mark+seq2m+fas->seq.l] = 0;
    }
    if (fas->qual.l > 0) {
        fas->qual.s = buffer + mark + qual2m;
        buffer[mark+qual2m+fas->qual.l] = 0;
    }
    
    fas->pBuffer = curBuffer;
    //if (curBuffer->buffer != buffer) fprintf(stdout, "NOT SAME BUFFER");
    
    /**
    fprintf(stderr,"[NAME %u @%s %s]\n[SEQ %u %s]\n[QUAL %u %s]\n",fas->name.l,fas->name.s,fas->comment.s,fas->seq.l,fas->seq.s,fas->qual.l,fas->qual.s); fflush(stderr);
    if (fas->qual.l == 0) {
        int i;
        char *p = fas->seq.s + fas->seq.l - 5;
        for (i=0; i < 15; i++) {
            fprintf(stderr, "{%hu %c%c}", *(p+i), *(p+i), (p+i == fas->qual.s ? '<' : 0)); 
        }
        fprintf(stderr, "\n"); fflush(stderr);
    }
     **/
    //// fprintf(stdout,"@%s %s\n%s\n+\n%s\n",fas->name.s,fas->comment.s,fas->seq.s,fas->qual.s);
    
////
    /**
    if (fas->name.l > 0 || fas->seq.l > 0) {
        fprintf(stderr, "::reconfigureFAS:: SeqLen=%d, QuaLen=%d, Name=[%s], Comment=[%s]\n", fas->seq.l, fas->qual.l, fas->name.s, fas->comment.s);
        fprintf(stderr, "::reconfigureFAS:: Next in buffer %d=[",current);
        int i;
        for (i=current; i < current + 10; i++) fprintf(stderr, "%c",buffer[i]);
        fprintf(stderr, "]\n");
    }
     */
////

}

void inline ogFastAQGZreader::checkToMove() {
    if (copyFrom > 0 && copyFrom > 0) {
        uint32_t d = copyTo - copyFrom;
////
    //// fprintf(stderr, "::checkToMove:: moving memory to=%d (ascii=%hu), from=%d (ascii=%hu='%c'), len=%d\n", heap, buffer[heap], copyFrom, (unsigned char) buffer[copyFrom], buffer[copyFrom], d);
////
        memmove(buffer + heap, buffer + copyFrom, d);
        heap += d;
        copyFrom = copyTo = 0;
    }
}

char ogFastAQGZreader::getNextSequenceOriginal(ogFastAQsequence *fas) {
    //ogFastAQsequence *fas = (ogFastAQsequence *) calloc(1, sizeof(ogFastAQsequence));
    uint32_t    nam2m=0, com2m=0, qual2m=0, seq2m=0; // distances to mark
    uint32_t    posMove = 0;

    char c, multiline=0;
    
    copyFrom = copyTo = 0; // there is nothing to copy
    
    // find sequence starting
    while (eof == 0 && (c = get1char()) != '>' && c != '@');
    if (eof != 0 && current >= end) return 0; 
    mark = current; // set mark = current
    nam2m = 0;
    
    // Find space or delimiter
    while (eof == 0 && ((c = get1char()) > ' '));
    fas->name.l = (current - mark) - 1;
    
    // Is there a "comment" ?
    if (c != '\n' && c != '\f' && c != '\r') {
        com2m = current - mark;
        while (eof == 0 && (c=get1char()) >= ' ');
        fas->comment.l = (current - mark) - com2m - 1;
    }
    consumeForNewLines(c);
    
    // Take next as sequence until blank line or eof or '>'
    seq2m = current - mark; // relative to mark because it may change
    while (eof == 0 && ((c = get1char()) >= ' '));
        
    heap = current - 1;
    // this is a newline or delimiter
    consumeForNewLines(c);

    while (eof == 0 && (c = get1char()) != '>' && c != '+') {
        // sigue siendo secuencia "FASTA/FASTQ"
        //
        copyFrom = current - 1;
        while (eof == 0 && ((c = get1char()) >= ' '));
        copyTo = current - 1;
        // this is a newline or delimiter
        consumeForNewLines(c);
        copyFrom = copyTo = 0;
        multiline = 1;
    }
        
    // Finish sequence
    fas->seq.l = (heap - mark) - seq2m;

    if (c == '>' || eof != 0) {
        current--;
    } else {
        // has to be '+'
        // it was '+' and therefore FASTQ
        // advance until next line
        while (eof == 0 && ((c = get1char()) >= ' '));
        consumeForNewLines(c);

        // FASTQ: Quality measures
        qual2m = current - mark;

        c = 0;
        // guess qual final position 1 line
        if (multiline != 0 && current + fas->seq.l < end) {
            // final position is within the buffer, we can check if there is a new line
            c = buffer[current+fas->seq.l];
            if (c == '\n' || c == '\f' || c == '\r') {
                fas->qual.l = fas->seq.l;
                c = 1;
                current += fas->seq.l;
                consumeForNewLines('\n');
            } else {
                c = 0;
            }
        }
        if (c == 0) {
            while (eof == 0 && ((c = get1char()) >= ' '));
            heap = current - 1;
            consumeForNewLines(c);
            while (eof == 0 && (c = get1char()) != '>' && c != '@') {
                // sigue siendo secuencia "FASTA/FASTQ"
                copyFrom = current - 1;
                while (eof == 0 && ((c = get1char()) >= ' '));
                copyTo = current - 1;
                // this is a newline or delimiter
                consumeForNewLines(c);
                copyFrom = copyTo = 0;
            }
            fas->qual.l = (heap - mark) - qual2m;
            if (c == '>' || c == '@') current--;
        }
    }
    
    reconfigureFAS(fas, nam2m, com2m, seq2m, qual2m);
    return (fas->seq.l ? 1 : 0);
}

char ogFastAQGZreader::getNextSequenceState(ogFastAQsequence *fas) {

    uint32_t    nam2m=0, com2m=0, qual2m=0, seq2m=0; // distances to mark
    char multiline=0;
    
    fas->comment.l = 0;
    fas->qual.l = 0;
    fas->name.l = 0;
    fas->seq.l = 0;
    
    uint32_t i;
    char c;
    char state = 10;
    
    mark = current;
    do {
        for (; current < end; current++) {
            c = buffer[current];
            switch (state) {
                case 10:
                    // Starting, expecting sequence name
                    if (c == '>' || c == '@') {
                        state = 11;
                        mark = current + 1;
                        nam2m = 0; // next char
                        break;
                    }
                    if (c > 32) state = 0;
                    break; // c == '\n' || c == '\f' || c == '\r' or ' ' or tab or so
                case 11 :
                    // Capturing sequence name
                    while (c > ' ' && current < end) c = buffer[++current];
                    if (current < end) {
                        fas->name.l = (current - mark);
                        if (c == ' ' || c == '\t') {
                            state = 12;
                            com2m = current + 1;
                        } else {
                            state = 13;
                        }
                    }
                    break;
                case 12:
                    // Capturing name's comments 
                    while (c >= ' ' && current < end) c = buffer[++current];
                    if (current < end) {
                        fas->comment.l = (current - com2m);
                        state = 13;
                        seq2m = current + 1;
                    }
                    break;
                case 13:
                    // Capturing sequence
                    while (c >= ' ' && current < end) c = buffer[++current];
                    break;
                case 127 : 
                    // Error
                default:
                    //
                    ;
            }
        }
        if (current >= end) read();
    } while (eof == 0);
    return (fas->seq.l ? 1 : 0);
}




char ogFastAQGZreader::getNextSequenceNaive(ogFastAQsequence *fas) {

    uint32_t    nam2m=0, com2m=0, qual2m=0, seq2m=0; // distances to mark

    char consuming = 0; // 0 = none, 1 = name, 2 = comment, 3 = sequence, 4 = ignored comment in +, 5 = quality
    char multiLine = 0;
    
    fas->comment.l = 0;
    fas->qual.l = 0;
    fas->name.l = 0;
    fas->seq.l = 0;
    
    uint32_t i;
    char c;
    char last = '\n';
    
    mark = current;
    heap = 0;
    do {
        for (; current < end; current++) {
            c = buffer[current];
            /**
            if (current % 100000 == 0) { 
                fprintf(stderr, "Current = %u, consuming=%d, buffer=%p; end=%u\n", current, consuming, buffer, end); 
                fflush(stderr); 
            }
             **/
            switch (c) {
                case '>':
                case '@':
                    // 
                    if (last == '\n') {
                        // starting or finishing
                        if (consuming == 3) {
                            // finishing
                            reconfigureFAS(fas, nam2m, com2m, seq2m, qual2m);                 
                            return (fas->seq.l ? 1 : 0);            
                        } else if (consuming >= 4) {
                            last = 'X';
                            if (qual2m == 0) {
                                heap = 0;
                                qual2m = current - mark;
                                if (multiLine == 0 && current + fas->seq.l < end && buffer[current + fas->seq.l] == '\n') {
                                    // uint32_t cu = current;
                                    fas->qual.l = fas->seq.l;
                                    current += fas->seq.l + 1;
                                    reconfigureFAS(fas, nam2m, com2m, seq2m, qual2m);
                                    // printf("Éxito: -3: [%d %d %d] : Qual:[%s] : +3:[%d %d %d]\n", buffer[cu-3], buffer[cu-2], buffer[cu-1], buffer+cu, buffer[current],buffer[current+1],buffer[current+2]);
                                    return (fas->seq.l ? 1 : 0);
                                }
                                while (++current < end && buffer[current] > ' ');
                                current--;
                            } else {
                                // Tiene que terminar porque @ no puede estar al inicio del quality de una 2a linea (tipo multiline)
                                reconfigureFAS(fas, nam2m, com2m, seq2m, qual2m);                 
                                return (fas->seq.l ? 1 : 0);            
                            }
                        } else {
                            // starting
                            mark = current;
                            nam2m = 1;          // this char has to be dismiss
                            last = '>';
                            consuming = 1;
                        }
                    }
                    break;

                case '\n': 
                case '\f': 
                case '\r':
                    if (last != '\n') {
                        // first end line
                        switch (consuming) {
                            case 0:
                                // nothing
                                break;
                            case 1:
                                // name
                                fas->name.l = current - mark - nam2m;
                                break;
                            case 2:
                                // comment
                                fas->comment.l = current - mark - com2m;
                                consuming = 3;
                                break;
                            case 3:
                                // sequence
                                copyTo = current;
                                checkToMove();
                                fas->seq.l = (heap == 0 ? current : heap) - mark - seq2m; // update
                                break;
                            case 4:
                                // Comment in +
                                // ignore text but change to quality
                                consuming = 5;
                                break;
                            case 5:
                                // quality
                                copyTo = current;
                                checkToMove();
                                fas->qual.l = (heap == 0 ? current : heap) - mark - qual2m; // update
                                if (multiLine == 0) {
                                    reconfigureFAS(fas, nam2m, com2m, seq2m, qual2m);
                                    // printf("Éxito: -3: [%d %d %d] : Qual:[%s] : +3:[%d %d %d]\n", buffer[cu-3], buffer[cu-2], buffer[cu-1], buffer+cu, buffer[current],buffer[current+1],buffer[current+2]);
                                    return (fas->seq.l ? 1 : 0);
                                }
                                while (++current < end && buffer[current] > ' ');
                                current--;
                                break;
                        }
                        last = '\n';
                        if (heap == 0) heap = current;
                    }
                    break;

                case '\t':
                case ' ':
                    if (consuming == 1) {
                        // in name
                        fas->name.l = current - mark - nam2m;
                        com2m = current - mark + 1;
                        consuming = 2;
                    }
                    last = '\t';
                    break;

                case '+':
                    if (last == '\n') {
                        // switch to quality, 
                        if (consuming == 3) {
                            consuming = 4; // but before consume until line in case we have comment like in nanopore sequences
                        } else if (consuming == 4) {
                            // also in case of + starting in quality line
                            consuming = 5;
                            heap = 0;
                            qual2m = current - mark;
                            if (multiLine == 0 && current + fas->seq.l < end && buffer[current + fas->seq.l] == '\n') {
                                // uint32_t cu = current;
                                fas->qual.l = fas->seq.l;
                                current += fas->seq.l + 1;
                                reconfigureFAS(fas, nam2m, com2m, seq2m, qual2m);
                                // printf("Éxito: -3: [%d %d %d] : Qual:[%s] : +3:[%d %d %d]\n", buffer[cu-3], buffer[cu-2], buffer[cu-1], buffer+cu, buffer[current],buffer[current+1],buffer[current+2]);
                                return (fas->seq.l ? 1 : 0);
                            }
                            while (++current < end && buffer[current] > ' ');
                            current--;
                        }
                    }
                    last = '+';
                    break;
                case 0:
                    break;
                default :
                    if (last == '\n') {
                        switch (consuming) {
                            case 0:
                                // nothing
                                break;
                            case 1:
                                // name
                                break;
                            case 2:
                                // comment
                                break;
                            case 3:
                                // sequence
                                if (seq2m == 0) {
                                    heap = 0;
                                    seq2m = current - mark;
                                } else {
                                    multiLine = 1;
                                    copyFrom = current;
                                }
                                break;
                            case 4:
                                break;
                            case 5:
                                // quality
                                if (qual2m == 0) {
                                    heap = 0;
                                    qual2m = current - mark;
                                    if (multiLine == 0 && current + fas->seq.l < end && buffer[current + fas->seq.l] == '\n') {
                                        // uint32_t cu = current;
                                        fas->qual.l = fas->seq.l;
                                        current += fas->seq.l + 1;
                                        reconfigureFAS(fas, nam2m, com2m, seq2m, qual2m);
                                        // printf("Éxito: -3: [%d %d %d] : Qual:[%s] : +3:[%d %d %d]\n", buffer[cu-3], buffer[cu-2], buffer[cu-1], buffer+cu, buffer[current],buffer[current+1],buffer[current+2]);
                                        return (fas->seq.l ? 1 : 0);
                                    }
                                } else {
                                    copyFrom = current;
                                }
                                break;
                        }
                        last = 'X';
                    } else {
                        //optimization
                        while (++current < end && buffer[current] > ' ');
                        current--;
                    }
                    break;
            }
        }
        if (current >= end && eof == 0) cleanBufferAndRead();
    } while (eof == 0);
    
    if (last != '\n') {
        // eof and finish without new line
        if (consuming == 3) {
            copyTo = current;
            checkToMove();
            fas->seq.l = heap - mark - seq2m;
        } else if (consuming == 5) {
            // quality
            copyTo = current;
            checkToMove();
            fas->qual.l = heap - mark - qual2m; // update
        }
    }
    reconfigureFAS(fas, nam2m, com2m, seq2m, qual2m);                 
    return (fas->seq.l ? 1 : 0);            
}


void inline ogFastAQGZreader::cleanBufferAndRead() {
    
    // read from file
    if (mark == 0) {
        // make bigger buffer
        buffSize <<= 1;
        buffer = (char *) realloc(buffer, buffSize);
        curBuffer->buffSize = buffSize;
        
        curBuffer->buffer = buffer;
        //fprintf(stderr, "buffsize:%u ",buffSize);
        
    } else {
        // just move the rest of the buffer pending

        uint32_t i, j;
        //fprintf(stderr, "buffsize:%u ",buffSize);
        /**
        fprintf(stderr, "___");
        for (j = mark; j < end; j++) {
            fprintf(stderr, "%c", buffer[j] == 0 ? '|' : buffer[j]);
        }
        fprintf(stderr, "___");
         **/
        ogFastAQBuffer *bufLibre = NULL;

        bufLibre = curBuffer;

        for (i=0; i < threads; i++) {
            if (curBuffer->readsProcessed[i] < curBuffer->readsDelivered[i]) {
                bufLibre = NULL;
                break;
            }
        }
        
        uint32_t pending = buffSize - mark;
        
        if (threads == 0 || bufLibre != NULL) {
            // fprintf(stderr, "R%c:%u ", (eof ? '1' : '0'), mark);
            memmove(buffer, buffer+mark, pending);
            current -= mark;
            end = pending;
            if (heap > 0) heap -= mark;
            if (copyFrom  > 0) copyFrom -= mark;
            if (copyTo  > 0) copyTo -= mark;
            mark = 0;
            //if (moveFrom > pending) moveFrom -= pending;
            //if (moveTo >= pending)  moveTo   -= pending;
        } else {

            // then check if there are no pendings
            // Look for another one
            ogFastAQBuffer *iBuf = fasBuffer, *lastBuf;
            int b = 0;
            int minPend = 0;
            while (iBuf != NULL && bufLibre == NULL) {
                bufLibre = iBuf;
                lastBuf = iBuf;
                for (i=0; i < threads; i++) { // the last thread is used for the root thread so do not need to be checked
                    if (iBuf->readsProcessed[i] < iBuf->readsDelivered[i]) {
                        minPend += iBuf->readsProcessed[i];
                        bufLibre = NULL;
                        iBuf = iBuf->nextBuffer;
                        b++;
                        break;
                    }
                }
            }
            if (bufLibre == NULL) {
                // Hay que crear otro
                curBuffer = (ogFastAQBuffer *) malloc(sizeof(ogFastAQBuffer));
                //fprintf(stderr, "\nr:%c, b:%u, d:%d, p:%d, mtotp:%d", readIdx+'0', lastBuf->buffNum, lastBuf->readsDelivered[i], lastBuf->readsPending[i], minPend);
                curBuffer->buffer = (char *) malloc( buffSize );
                curBuffer->nextBuffer = NULL;
                curBuffer->buffSize = buffSize;
                curBuffer->readsProcessed = (int *) calloc((threads+1)*2, sizeof(int));
                curBuffer->readsDelivered = curBuffer->readsProcessed + threads + 1; // (int *) calloc((threads+1)*1, sizeof(int)); //curBuffer->readsPending + threads + 1;
                curBuffer->buffNum = lastBuf->buffNum + 1;
                lastBuf->nextBuffer = curBuffer;
                lastBuffer = lastBuf;
                /***
                //if (lastBuf->buffNum % 10 == 0) {
                    fprintf(stderr, "\n");
                    iBuf = fasBuffer;
                    int b = 0;
                    while (iBuf != NULL) {
                        fprintf(stderr, "r%c Buffer:%d ==>", readIdx+'0', iBuf->buffNum);
                        for (i=0; i <= threads; i++) {
                            fprintf(stderr, "t%d:(%d)%d/%d ", i, iBuf->readsDelivered[i]-iBuf->readsProcessed[i], iBuf->readsProcessed[i],iBuf->readsDelivered[i]);
                        }
                        fprintf(stderr, "(%d)\n", b);
                        iBuf = iBuf->nextBuffer;
                        b++;
                    }
                //}
                **/
            } else {
                curBuffer = bufLibre;
                //fprintf(stderr, "r");
            }
            
            memcpy(curBuffer->buffer, buffer+mark, pending);
            current -= mark;
            end = pending;
            if (heap > 0) heap -= mark;
            if (copyFrom  > 0) copyFrom -= mark;
            if (copyTo  > 0) copyTo -= mark;
            mark = 0;
            buffer = curBuffer->buffer;
            
        }
    }
    
    // then read the rest of buffer
    read();
    
}

char ogFastAQGZreader::__ReadWhileTextCharsNoSpace() {
    do {
        while (current < end && buffer[current] > ' ') current++;
        if (current >= end) cleanBufferAndRead(); else break;
    } while (eof == 0); 
    return buffer[current];
}

char ogFastAQGZreader::__ReadUntilNewLines() {
    do { 
        while (current < end && buffer[current] >= ' ') current++;
        if (current >= end) cleanBufferAndRead(); else break;
    } while (eof == 0);
    return buffer[current];
}

// aqui lo mas lógico es que solo entre 1 vez (1 new line, entonces es mejor así)

char ogFastAQGZreader::__ReadWhileNewLines() {
    do {
        while (current < end && buffer[current] < ' ') current++;
        if (current >= end) cleanBufferAndRead(); else break;
    } while (eof == 0);
    return buffer[current];
}

            

char ogFastAQGZreader::getNextSequenceByType(ogFastAQsequence *fas) {

    uint32_t    nam2m=0, com2m=0, qual2m=0, seq2m=0; // distances to mark

    char consuming = 0; // 0 = none, 1 = name, 2 = comment, 3 = sequence, 4 = ignored comment in +, 5 = quality
    char multiLine = 0;
    
    fas->comment.l = 0;
    fas->qual.l = 0;
    fas->name.l = 0;
    fas->seq.l = 0;
    
    uint32_t i;
    char c;
    char last = '\n';
    
    mark = current;
    heap = 0;
    
    // identify first character
    do {
        c = __ReadWhileNewLines();
        if (eof) { 
            return 0;
        } else if (c == '>') {
            // FASTA
            // Read Name and Comment
            nam2m = current - mark + 1;
            c = __ReadWhileTextCharsNoSpace();
            fas->name.l = current - mark - nam2m;
            if (c == ' ' || c == '\t') {
                com2m = current - mark;
                c = __ReadUntilNewLines();
                fas->comment.l = current - mark - com2m;
            }
            c = __ReadWhileNewLines();
            // Next Get 1st line of sequence
            seq2m = current - mark;
            c =  __ReadWhileTextCharsNoSpace();
            heap = current;
            c = __ReadWhileNewLines();
            // read next lines if any
            while (c != '>' && eof == 0) {
                copyFrom = current;
                c = __ReadWhileTextCharsNoSpace();
                copyTo = current;
                checkToMove();
                c = __ReadWhileNewLines();
            }
            fas->seq.l = heap - mark - seq2m;
            // finish
            reconfigureFAS(fas, nam2m, com2m, seq2m, qual2m);                 
            return 1;
            
        } else if (c == '@') {
            
            // FASTQ
            // Read Name and Comment
            nam2m = current - mark + 1;
            c = __ReadWhileTextCharsNoSpace();
            //fprintf(stderr, "@%.100s@", buffer+mark+nam2m); fflush(stderr);
            fas->name.l = current - mark - nam2m;
            if (c == ' ' || c == '\t') {
                com2m = current - mark;
                __ReadUntilNewLines();
                fas->comment.l = current - mark - com2m;
            }
            c = __ReadWhileNewLines();
            // Next Get 1st line of sequence
            seq2m = current - mark;
            c = __ReadWhileTextCharsNoSpace();
            //fprintf(stderr, "/%.100s/", buffer+mark+seq2m); fflush(stderr);
            heap = current;
            c = __ReadWhileNewLines();
            // read next lines if any
            while (c != '+' && eof == 0) {
                //fprintf(stderr, "%c", c); fflush(stderr);
                multiLine = 1;
                copyFrom = current;
                c = __ReadWhileTextCharsNoSpace();
                copyTo = current;
                checkToMove();
                c = __ReadWhileNewLines();
            }
            fas->seq.l = heap - mark - seq2m;
            
            // + detected, now consume line/coment (nanopore)
            c = __ReadUntilNewLines();
            c = __ReadWhileNewLines();
            
            // Now read quality
            qual2m = current - mark;
            c = __ReadWhileTextCharsNoSpace();
            heap = current;
            c = __ReadWhileNewLines();
            if (multiLine) {
                // read next lines if any
                while (c != '@' && eof == 0) {
                    copyFrom = current;
                    c = __ReadWhileTextCharsNoSpace();
                    copyTo = current;
                    checkToMove();
                    c = __ReadWhileNewLines();
                }
            }
            fas->qual.l = heap - mark - qual2m;
            reconfigureFAS(fas, nam2m, com2m, seq2m, qual2m);
            
            return 1;
            
        } else {
            fprintf(stderr, "Error parsing FASTQ/FASTA: '%c'\n", c);
            //skip line
            while (current < end && buffer[current] >= ' ') current++;
        }
        if (current >= end) cleanBufferAndRead();
    } while (eof == 0);
    
    return 0;
}
