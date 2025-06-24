/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/file.h to edit this template
 */

/* 
 * File:   ogCigarOperations.hpp
 * Author: victortrevino
 *
 * Created on March 20, 2024, 4:55 PM
 */

#ifndef OGCIGAROPERATIONS_HPP
#define OGCIGAROPERATIONS_HPP


#define MIN_CIGAR_OPERATION 12

typedef struct cigarOperation {
    uint32_t        size;
    char            operation;
} OG_CIGAR_OPERATION;


class ogCigarOperations {
    
    OG_CIGAR_OPERATION      *cigarOp;
    uint16_t                nCigarOperations;
    uint32_t                sizeSum;
    uint32_t                sizeSumConsuming;

public:
    ogCigarOperations();
    
    uint32_t                leftGenPos;
    uint32_t                rightGenPos;
    uint32_t                leftGenPosSet;
    uint32_t                rightGenPosSet;
    
    void                    reset(uint32_t leftPos, uint32_t rightPos);
    void                    push(char op, uint32_t len);
    OG_CIGAR_OPERATION      *pull();
    OG_CIGAR_OPERATION      *last();
    uint32_t                buildCIGAR(char *pCigar, uint32_t maxCigarLen);
    uint32_t                totalSize();
    uint32_t                totalSizeConsuming();
    uint32_t                sumSizeConsuming(uint16_t from, uint16_t to);
    uint32_t                sumOperationsType(char op);
    virtual                ~ogCigarOperations();
    OG_CIGAR_OPERATION     *getOperation(uint16_t x);
    void                    removeOperation(uint16_t x);
    void                    setOperation(uint16_t x, char op, uint32_t len);
    char                    reduceOperations();
    uint16_t                pushOperationsFromCIGAR(char *pCIGARchars);
    
private:
    
};


#endif /* OGCIGAROPERATIONS_HPP */

