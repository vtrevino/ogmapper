/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/file.h to edit this template
 */

/* 
 * File:   ogReadManager.hpp
 * Author: victortrevino
 *
 * Created on March 27, 2024, 7:18 PM
 */

#ifndef OGREADMANAGER_HPP
#define OGREADMANAGER_HPP

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "ogDefinitions.hpp"
#include "ogGenome.hpp"

class ogReadManager {
    
public:

    uint32_t        maxReads;
    
    void            **theMallocs;           // mallocs calls saved
    uint32_t        *sizeOfMallocs;         // size of each malloc call
    
    void            **availableMallocs;
    uint32_t        availMallocNext;
    
    void            **busyMallocs;
    uint32_t        busyMallocNext;
    
    ogSingleRead    **availableReads;
    uint32_t        availReadsNext;
    
    ogSingleRead    **busyReads;
    uint32_t        busyReadsNext;
    

private:
    
}

#endif /* OGREADMANAGER_HPP */

