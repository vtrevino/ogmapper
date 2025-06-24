/* 
 * File:   ogDefinitions.hpp
 * Author: victortrevino
 *
 * Created on May 15, 2022, 12:54 PM
 */

#ifndef OGDEFINITIONS_HPP
#define OGDEFINITIONS_HPP

#ifdef __cplusplus
extern "C" {
#endif

#include <time.h>
#include <stdint.h>

#define OGMAPPER_VERSION    "v0.5.0-20-Jan-2025"

#define getAlignmentCigar   getAlignment
    
#define setFileName(str,sourceFileName,keySize,ext,bufsize)         snprintf(str, bufsize, "%s-%d.%s", sourceFileName, keySize, ext);
#define setPureFileName(str,sourceFileName,ext,bufsize)             snprintf(str, bufsize, "%s.%s", sourceFileName, ext);
#define printf_FileOperationFile(msg, xfile)                        fprintf(stderr, ":: %s [%s] ::",msg, xfile); fflush(stderr);
#define printf_FileOperation(msg)                                   fprintf(stderr, ":: %s ::",msg); fflush(stderr);
#define printf_FileOperationDone()                                  fprintf(stderr, "\t[Ok]\n"); fflush(stderr);
#define printf_FileOperationMsgDone(msg)                            fprintf(stderr, ":: %s :: [Ok]\n", msg);  fflush(stderr);

// SIZES
#define MAX_FILENAME                512
#define MAX_OGINDEX_FILENAME        512
#define MAX_VERSION_NAME            32
#define MAX_ENCODINGS               16
#define MAX_SCHEDULE_SIZE           100
#define MIN_CANDIDATE_POSITIONS     8192    //4096    //2048    //65536   //1048576
#define MAX_DEF_READ_LEN            1000
#define MAX_CHAR_NAMES              64
#define MAX_CIGAR_SIZE              255
#define TOO_MANY_CANDIDATE_REGIONS  128
#define ACEPTABLE_CANDIDATE_REGIONS 16
#define MAX_PAIR_INSERT_SIZE        1500
#define BUFFER_FOR_UNPACKING_DNA    10000
// Falta definir un max_read_len


// CALCULATIONS    
#define MAX_GAP_1READ_ALIGMENT      1000
#define MAX_ERROR_RATE_WFA          0.0625              // ~ 1/16
#define MAX_KEY_ERR(keysize)        (keysize / 10)
    // DEPRACATED::
#define MAX_SIMPLE_ALIGN_ERR(seqLen)    (seqLen >> 4)
#define MIN_SIMPLE_ALIGN_ERR(seqLen)    (seqLen >> 1)
    // ::DEPRACATED
#define MAX_ERROR_ACCEPT1_PACKED(seqLen)    (seqLen >> 4)           // ~1/16 = 0.0625 = 93.5% ok ==> Cycle (1) : FINE ACCEPT IMMEDIATLY
#define MAX_ERROR_ACCEPT1_ALIGN(seqLen)     (seqLen >> 3)           // ~1/8 = 0.125 = 87.5% ok  ==> Cycle (1) : PERFORM INS/DEL BUT ACCEPT ANYWAY
#define MAX_ERROR_CHECK2_ALIGN(seqLen)     (seqLen >> 2)            // ~1/4 = 0.25 = 75% ok ==> Cycle (2) : TRY INS/DEL ONLY IF NO ACCEPT RECORDS, ACCEPT IF > ACCEPT_ALIGN
#define MAX_ERROR_CHECK3_ALIGN(seqLen)     (seqLen >> 1)            // WORTH TO BE CHECKED ==>  Cycle (3) : TRY INS/DEL ONLY IF NO ACCEPTED RECORDS, ACCEPT IF > ACCEPT_ALIGN
#define NT_EXTENSION                       00                       // Nucleotides to be extended at extremes for the alignment
    

// Nomenclature BROAD
// https://broadinstitute.github.io/picard/explain-flags.html
#define SAMFLAG_READ_PAIRED             1
#define SAMFLAG_READS_MAPPED            2
#define SAMFLAG_READ_UNMAPPED           4
#define SAMFLAG_MATE_UNMAPPED           8
#define SAMFLAG_READ_REVCOMP            16
#define SAMFLAG_MATE_REVCOMP            32
#define SAMFLAG_READ1                   64
#define SAMFLAG_READ2                   128
#define SAMFLAG_NOT_PRIMARY_ALIGNMENT   256
#define SAMFLAG_ALIGNMENT_FAILS         512
#define SAMFLAG_PCR_OPTICAL_DUPLICATE   1024
#define SAMFLAG_SUPPLEMENTARY_ALIGNMENT 2048

const char EQUAL[] = "=";
const char ASTERISK[] = "*";
const char TAB[] = "\t";
const char EMPTY[] = "";
    
    
using namespace std;


typedef struct ogSAM {
    char        *qname;     // Query template NAME
    int         flag;       // bitwise FLAG
    char        *rname;     // References sequence NAME
    uint64_t    pos;        // 1- based leftmost mapping POSition
    int         mapq;       // MAPping Quality
    char        *cigar;     // CIGAR string
    char        *rnext;     // Ref. name of the mate/next read
    uint64_t    pnext;      // Position of the mate/next read
    int         tlen;       // observed Template LENgth
    char        *seq;       // segment SEQuence
    char        *qual;      // ASCII of Phred-scaled base QUALity+33
    
    int         editScore;  // for NM:i annotation
    int         alignScore; // for AS:i annotation
    uint32_t    genomicPos; // absolute position reported by ogMapper
} OG_SAM;

const unsigned char ascii2Complement[256] = {
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
 60, 61, 62, 63, 64, 'T', 66, 'G', 68, 69,
 70, 'C', 72, 73, 74, 75, 76, 77, 78, 79,
 80, 81, 82, 83, 'A', 'A', 86, 87, 88, 89,
 90, 91, 92, 93, 94, 95, 96, 't', 98, 'g',
 100, 101, 102, 'c', 104, 105, 106, 107, 108, 109,
 110, 111, 112, 113, 114, 115, 'a', 'a', 118, 119,
 120, 121, 122, 123, 124, 125, 126, 127, 128, 129,
 130, 131, 132, 133, 134, 135, 136, 137, 138, 139,
 140, 141, 142, 143, 144, 145, 146, 147, 148, 149,
 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
 160, 161, 162, 163, 164, 165, 166, 167, 168, 169,
 170, 171, 172, 173, 174, 175, 176, 177, 178, 179,
 180, 181, 182, 183, 184, 185, 186, 187, 188, 189,
 190, 191, 192, 193, 194, 195, 196, 197, 198, 199,
 200, 201, 202, 203, 204, 205, 206, 207, 208, 209,
 210, 211, 212, 213, 214, 215, 216, 217, 218, 219,
 220, 221, 222, 223, 224, 225, 226, 227, 228, 229,
 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
 240, 241, 242, 243, 244, 245, 246, 247, 248, 249,
 250, 251, 252, 253, 254, 255 };


typedef struct exon {
    uint32_t    chr;
    uint32_t    start;
    uint32_t    end;
    char        strand;
    uint32_t    phyStart; // physical start in the pseudo fasta sequence
    uint32_t    phyEnd;   // physical end   in the pseudo fasta sequence
    uint32_t    iExon;
    uint32_t    iGene;
    char        overlapsPrev;
} Exon;

typedef struct transcript {
    char        transcriptID[32];
    uint32_t    iTranscript;
    uint32_t    nExons;
    Exon        *pExons;
} Transcript;

typedef struct CountingInfo {
    char        GTFfile[MAX_FILENAME];
    char        OGXfile[MAX_FILENAME];
    char        destFileName[MAX_FILENAME];
    uint32_t    nGenes, nExons, nTranscripts, nChromTrans;
    
    // Memory use
    Exon        *pExons;
    uint32_t    *pExonReads;
    uint32_t    *pGeneReads;
    Transcript  *pTranscripts;
    uint32_t    *pGeneExonStart; // direct access (index) where the first exon of the gene i begins
    uint32_t    *pGeneExon; // to define which exons needs to be counted
    uint32_t     nGeneExon;
    
    // 
} COUNTING_INFO;


#ifdef __cplusplus
}
#endif

#endif /* OGDEFINITIONS_HPP */

