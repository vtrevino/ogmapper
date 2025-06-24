/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.h to edit this template
 */

/* 
 * File:   binnedBigCounters.hpp
 * Author: victortrevino
 *
 * Created on January 9, 2023, 9:58 PM
 */

#ifndef BINNEDBIGCOUNTERS_HPP
#define BINNEDBIGCOUNTERS_HPP


#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include "ogDefinitions.hpp"

typedef struct countvalue {
    uint32_t	count;
    uint32_t    value;
} countvalue;

typedef struct itemBigCounters {
    char	mode;		// 0 : nothing, 1: value, 2: pointer to partBigVector
    union	content {
        countvalue      data;
        void           	*pBV;
    } content;
} itemBigCounters;


class binnedBigCounters {
public:
    uint32_t	first;
    uint32_t	last;
    uint32_t	delta;
    uint32_t	bins;    
    char 	shifts;
    double	factor;     // multiplier to estimate position
    double       bin;        // width of each bin
    uint32_t	maxCount;
    uint32_t    maxCountItem;
    uint32_t    ties;
    char	exact; 	// this structure cannot be splitable further
    void	*next; 	// partBigVector
    void        *root;  // root object
    void        *tail;  // pointer to last Object only created by the first root object
    uint32_t    chainedObject;
    void	*back; 	// partBigVector
    void	*items;	// itemBigCounters
    void        *tieItemsBBC; // *binnedBigCounters;
    uint32_t    tiePos;
    uint64_t    calls;
    
                    binnedBigCounters(uint32_t from, uint32_t to, uint32_t nBins);
    virtual        ~binnedBigCounters();
    uint32_t        getPosition(uint32_t item);
    void            incCounter(uint32_t item);
    void            incCountForItem(itemBigCounters *p, uint32_t item);
    void            startTiedProcess();
    uint32_t        getNextTiedCounterValue();
    
private:

};

#endif /* BINNEDBIGCOUNTERS_HPP */

