/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.cc to edit this template
 */

/* 
 * File:   binnedBigCounters.cpp
 * Author: victortrevino
 * 
 * Created on January 9, 2023, 9:58 PM
 */
#include "binnedBigCounters.hpp"
#include <math.h>


binnedBigCounters::binnedBigCounters(uint32_t from, uint32_t to, uint32_t nBins) {
    first = from;
    last = to;
    bins = nBins;
    delta = last - first + 1;
    maxCount = 0;
    ties = 0;
    bin = (double) (delta) / (double) bins;
    factor = (double) bins / (double) (delta);
    /**
    shifts = (int) (log2f((float) (last - first + 1) / (float) bins));
    uint32_t div = 2 << (shifts-1);
    uint16_t binsByShifts = (last - first + 1) / ((float) div);
    if (binsByShifts > 2*bins) {
        shifts = 0;
    } else {
        bins = (nBins < binsByShifts ? binsByShifts : nBins);
    }
     **/
    shifts = 0; // force not to 
    exact = (delta <= bins ? 1 : 0);
    if (exact) {
        bins = delta;
        factor = 1;
        bin = 1;
        shifts = 0;
    }
    /*****
    switch(bins) {
        case 2 : shifts = 1; break;
        case 4 : shifts = 2; break;
        case 8 : shifts = 3; break;
        case 16 : shifts = 4; break;
        case 32 : shifts = 5; break;
        case 64 : shifts = 6; break;
        case 128 : shifts = 7; break;
        case 256 : shifts = 8; break;
        case 512 : shifts = 9; break;
        case 1024 : shifts = 10; break;
        case 2048 : shifts = 11; break;
        case 4096 : shifts = 12; break;
        case 8192 : shifts = 13; break;
        case 16384 : shifts = 14; break;
        case 32768 : shifts = 15; break;
        case 65536 : shifts = 16; break;
        case 131072 : shifts = 17; break;
        case 262144 : shifts = 18; break;
        case 524288 : shifts = 19; break;
        case 1048576 : shifts = 20; break;
        case 2097152 : shifts = 21; break;
        case 4194304 : shifts = 22; break;
        case 8388608 : shifts = 23; break;
        case 16777216 : shifts = 24; break;
        case 33554432 : shifts = 25; break;
        default : shifts = 0; 
    }
     *****/
    next = NULL;
    back = NULL;
    items = calloc(bins, sizeof(itemBigCounters));
    chainedObject = 0;
    root = this;
    tail = this;
    calls = 0;
}

binnedBigCounters::~binnedBigCounters() {
    
    if (back == NULL) {
        uint64_t totalBins = bins;
        uint32_t totalObjects = 1;
        binnedBigCounters *n,*n2;
        //for (n = (binnedBigCounters *) next; n != NULL; totalBins += n->bins, totalObjects++, n = (binnedBigCounters *) n->next);
        //fprintf(stderr, "TotalBins=%llu, TotalObjects=%u, calls=%llu\n", totalBins, totalObjects, calls);
        for (n = (binnedBigCounters *) next; n != NULL; n = n2) {
            n2 = (binnedBigCounters *) n->next;
            delete n; // this is a "lineal" deletion of all objects
        }
    }
    //fprintf(stderr, "object=%u\n", chainedObject);
    
    //if (next != NULL) delete ((binnedBigCounters *) next); // this creates a problem in the stack for more than 64K objects
    free(items);
}

uint32_t binnedBigCounters::getPosition(uint32_t item) {
    //if (shifts > 0) {
    //  return (item - first) >> shifts;
    //} else {
        return (item - first) * factor;
    //}
}

void binnedBigCounters::incCounter(uint32_t item) {
    calls++;
    //fprintf(stderr, "%p\n", items);
    itemBigCounters *p = (itemBigCounters *) items;
    uint32_t pos = getPosition(item);

    if (pos >= bins) {
        // Problemas, pos=11286, bins=8192, first=0, last=3117292070, item=4294967216, bin=380528.817261
        fprintf(stderr, "Problemas, pos=%u, bins=%u, first=%u, last=%u, item=%u, bin=%f\n", pos, bins, first, last, item, bin);
        if (back != NULL) {
            fprintf(stderr, "Encadenado de :\n");
            binnedBigCounters *b = (binnedBigCounters *) back;
            while (b != NULL) {
                fprintf(stderr, "bins=%u, first=%u, last=%u, bin=%f\n", b->bins, b->first, b->last, b->bin);
                b = (binnedBigCounters *) b->back;
            }
        }
    }
    p += pos;
    if (p->mode) {
        if (p->mode == 1) {
            if (exact) {
                incCountForItem(p, item);
            } else {
                if (p->content.data.value != item) {
                    // Create annother Binned Object
                    uint32_t preCount = p->content.data.count;
                    uint32_t preValue = p->content.data.value;
                    uint32_t f = round(bin * pos) + first;
                    //while (f > 0 && getPosition(f) > pos) { fprintf(stderr, "f--\n"); f--; }
                    //while (f > 0 && getPosition(f-1) == pos) { fprintf(stderr, "f--\n"); f--; }
                    uint32_t l = round(bin * (pos+1)) + first - 1;
                    while (getPosition(l+1) == pos) l++;
                    if (l > last) l = last;
                    //fprintf(stderr, "Creando objeto por item=%u (pos=%u), first=%u, last=%u\n", item, pos, f, l);
                    binnedBigCounters *p1 = new binnedBigCounters(f, l, bins);
                    p1->root = this->root;
                    p->mode = 2;
                    p->content.pBV = p1;
                    p1->back = this;
                    while (preCount-- > 0) p1->incCounter(preValue);
                    p1->incCounter(item);

                    // chain
                    // OLD : binnedBigCounters *n = this;
                    binnedBigCounters *n = (binnedBigCounters *) ((binnedBigCounters *) this->root)->tail;
                    while (n->next != NULL) n = (binnedBigCounters *) n->next;
                    n->next = p1;
                    p1->chainedObject = n->chainedObject + 1;
                    //if (p1->chainedObject > 9 && p1->chainedObject % 1000 == 0) {
                    //    fprintf(stderr, "%u ",p1->chainedObject);
                    //}
                    ((binnedBigCounters *) this->root)->tail = p1;
                } else {
                    incCountForItem(p, item);
                }
            }
        } else {
            binnedBigCounters *p2 = (binnedBigCounters *) p->content.pBV;
            p2->incCounter(item);
        }
    } else {
        p->mode = 1;
        p->content.data.count = 1;
        p->content.data.value = item;
    }
}



void binnedBigCounters::incCountForItem(itemBigCounters *p, uint32_t item) {
    if (++p->content.data.count >= maxCount) {
        if (p->content.data.count == maxCount) {
            ties++;
        } else {
            maxCount = p->content.data.count;
            maxCountItem = item;
            ties = 0;
        }
        // Propagate
        binnedBigCounters *b = (binnedBigCounters *) back;
        while (b != NULL) {
            if (b->maxCount <= maxCount) {
                if (b->maxCount == maxCount) {
                    b->ties++;
                } else {
                    b->maxCount = maxCount;
                    b->maxCountItem = maxCountItem;
                    b->ties = 0;
                }
                b = (binnedBigCounters *) b->back;
            } else {
                break;
            }
        }
    }
}

void binnedBigCounters::startTiedProcess() {
    tieItemsBBC = this;
    tiePos = 0;
    //fprintf(stderr, "Starting tied Process, maxCount=%u, ties=%u\n", maxCount, ties);
}


uint32_t binnedBigCounters::getNextTiedCounterValue() {
    // where it was suspended
    binnedBigCounters *pCurr = (binnedBigCounters *) tieItemsBBC;
    uint32_t curBins;
    uint32_t pos = tiePos;
    while (pCurr != NULL) {
        curBins = pCurr->bins;
        itemBigCounters *pCurrIBC = (itemBigCounters *) pCurr->items;
        if (pCurr->maxCount == maxCount) {
            //fprintf(stderr, "pCurr=%p, pCurrIBC=%p, pos=%u, maxCount=%u, ties=%u, bins=%u, first=%u, last=%u, mode=%c\n", pCurr, pCurrIBC, pos, pCurr->maxCount, pCurr->ties, pCurr->bins, pCurr->first, pCurr->last, pCurrIBC->mode+48);
            // Yes, search
            itemBigCounters *p = pCurrIBC + pos;
            while (pos < curBins) {
                if (p->mode == 1 && p->content.data.count == maxCount) {
                    //fprintf(stderr, "Fount at pCurr=%p, pCurrIBC=%p, pos=%u, maxCount=%u, ties=%u, bins=%u, first=%u, last=%u\n", pCurr, pCurrIBC, pos, pCurr->maxCount, pCurr->ties, pCurr->bins, pCurr->first, pCurr->last);
                    tieItemsBBC = pCurr;
                    tiePos = pos+1;
                    return p->content.data.value;
                }
                p++;
                pos++;
            }
        }
        pCurr = (binnedBigCounters *) pCurr->next;
        pos = 0;
    }
    return 0;
}
