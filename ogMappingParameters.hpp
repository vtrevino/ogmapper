/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/file.h to edit this template
 */

/* 
 * File:   ogMappingParameters.hpp
 * Author: victortrevino
 *
 * Created on December 22, 2022, 8:10 PM
 */

#ifndef OGMAPPINGPARAMETERS_HPP
#define OGMAPPINGPARAMETERS_HPP

#include "ogKeys.hpp"
#include "ogGenome.hpp"
#include "ogKeyEncoding.hpp"
#include "ogGenomePositions.hpp"
#include "ogCandidatePosManager.hpp"
#include "ogReadKeyMapping.hpp"

typedef struct ogMappingParameters {
    ogKeys              *pKeys;
    ogGenome            *pGenome;
    ogKeyEncoding       *pKeyEncoding;
    ogGenomePositions   *pGenPos;
    ogGuider            *pGuiderDontUseIt;
    char                *cmdLine;
    //std::mutex           mtx;
    
    uint8_t              keySize;
    char                 forceUpperCase;
    //char                *seq1, *seq2;
    //uint32_t             seq1Len, seq2Len;
    //ogCandidatePosManager    *pCandPosMan;
    
    // For counting functions
    uint32_t            *pExonCounts;
    uint32_t            *pGeneCounts;
    uint32_t            *pTranscriptCounts;
    uint32_t            nGenes;
    uint32_t            nExons;
    uint32_t            nTranscripts;
    char                countExons;
    char                countTranscripts;
    char                outputSAMcount;
    Exon                *pExons;
    uint32_t            *pGeneExonStart;
    uint32_t            *pGeneExon;
    uint32_t            nGeneExon;
    uint16_t            MIN_KEY_COUNTS;     // Default: 3
    char                MIN_SCORE_COUNT;    // Score is matchingKeys * 100 / totalKeys, default: 25
    char                BREAK_TIE_COUNT1;   // Default: 1
    uint16_t            MIN_MATCHES;        // number of nt matches in blast like
    //
    
    uint32_t            MAX_KEY_DISTANCE;
} OG_MAPPING_PARAMETERS;


#endif /* OGMAPPINGPARAMETERS_HPP */

