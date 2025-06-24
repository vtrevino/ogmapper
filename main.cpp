#include <stdio.h>
#include <string.h>
#include <zlib.h>
//#include "testFastFiles.h"
#include "ogDefinitions.hpp"
#include "ogIndex.hpp"
#include <unistd.h>
#include "ogKeys.hpp"
#include "ogKeyEncoding.hpp"
#include "ogGenome.hpp"
#include "ogGenomePositions.hpp"
#include "ogBitwiseAT1GC0Encoding.hpp"
#include "ogGappedBitwiseAT1GC0Encoding.hpp"
#include "ogSwapBitwiseAT1GC0Encoding.hpp"
#include "ogHomoPolymerCompressedEnc.hpp"
#include "ogPlainEncoding.hpp"
#include "ogGuider.hpp"
#include "ogStateMachineGuider.hpp"
#include "ogTupleGuider.hpp"
#include "ogReadKeyMapping.hpp"
#include "ogFastAQGZreader.hpp"

extern uint16_t ExactInApartSides(ogReadKeyMapping *pMap) ;
extern uint16_t ExactInExtremes(ogReadKeyMapping *pMap) ;
extern uint16_t Two5pKeys(ogReadKeyMapping *pMap) ;
extern uint16_t DefaultMatchingSmall(ogReadKeyMapping *pMap) ;
extern uint16_t DefaultMatchingLarge(ogReadKeyMapping *pMap) ;
extern uint16_t DefaultMatchingHuge(ogReadKeyMapping *pMap) ;
extern uint16_t HistogramMatchingSmall(ogReadKeyMapping *pMap) ;
extern uint16_t HistogramMatchingLarge(ogReadKeyMapping *pMap) ;
extern uint16_t HistogramMatchingHuge(ogReadKeyMapping *pMap) ;
extern uint16_t HistogramMatchingAll(ogReadKeyMapping *pMap) ;
extern uint16_t MapPosMatching(ogReadKeyMapping *pMap);
extern uint16_t NestedIntersectSmall(ogReadKeyMapping *pMap);
extern uint16_t NestedIntersectLarge(ogReadKeyMapping *pMap);
extern uint16_t NonOverlappingKeys(ogReadKeyMapping *pMap);
extern uint16_t OverlappingBy3nt(ogReadKeyMapping *pMap);
extern uint16_t OverlappingBy6nt(ogReadKeyMapping *pMap);
extern uint16_t OverlappingBy9nt(ogReadKeyMapping *pMap);
extern uint16_t OverlappingBy12nt(ogReadKeyMapping *pMap);
extern uint16_t OverlappingBy15nt(ogReadKeyMapping *pMap);
extern uint16_t OverlappingBy20nt(ogReadKeyMapping *pMap);
extern uint16_t LeftMost2KeysQuick(ogReadKeyMapping *pMap);
extern uint16_t LeftMost2KeysLong(ogReadKeyMapping *pMap);
extern uint16_t LeftMost2KeysHuge(ogReadKeyMapping *pMap);
extern uint16_t Min2KeysSize(ogReadKeyMapping *pMap);
extern uint16_t Min2KeysLong(ogReadKeyMapping *pMap);
extern uint16_t Min2KeysHuge(ogReadKeyMapping *pMap);
extern uint16_t MinApart2Keys(ogReadKeyMapping *pMap);
extern uint16_t AnalyzeAllPos(ogReadKeyMapping *pMap);
extern uint16_t MovWinIntersectSmall(ogReadKeyMapping *pMap);
extern uint16_t MovWinIntersectLarge(ogReadKeyMapping *pMap);
extern uint16_t MovWinIntersectMicro(ogReadKeyMapping *pMap);
extern uint16_t MovWinIntersectHuge(ogReadKeyMapping *pMap);
extern uint16_t HalfMatchingSmall(ogReadKeyMapping *pMap);
extern uint16_t HalfMatchingLarge(ogReadKeyMapping *pMap);
extern uint16_t HalfMatchingHuge(ogReadKeyMapping *pMap);
extern uint16_t ExtremeMatchingSmall(ogReadKeyMapping *pMap);
extern uint16_t ExtremeMatchingLarge(ogReadKeyMapping *pMap);
extern uint16_t ExtremeMatchingHuge(ogReadKeyMapping *pMap);
extern uint16_t BinnedHistogramMatchingSmall(ogReadKeyMapping *pMap);
extern uint16_t BinnedHistogramMatchingLarge(ogReadKeyMapping *pMap);
extern uint16_t BinnedHistogramMatchingHuge(ogReadKeyMapping *pMap);
extern uint16_t BinnedHistogramMatchingAll(ogReadKeyMapping *pMap);
extern uint16_t blastLikeSmall(ogReadKeyMapping *pMap);
extern uint16_t blastLikeLarge(ogReadKeyMapping *pMap);
extern uint16_t blastLikeHuge(ogReadKeyMapping *pMap);
extern uint16_t blastLikeAll(ogReadKeyMapping *pMap);
extern uint16_t Intersect2R_Small(ogReadKeyMapping *pMap);
extern uint16_t Intersect2R_Large(ogReadKeyMapping *pMap);
extern uint16_t Intersect2R_Huge(ogReadKeyMapping *pMap);
extern uint16_t Intersect2R_All(ogReadKeyMapping *pMap);
//extern uint16_t blastLikeBinary(ogReadKeyMapping *pMap);
extern uint16_t CountGeneReadSmall(ogReadKeyMapping *pMap);
extern uint16_t CountGeneReadLarge(ogReadKeyMapping *pMap);
extern uint16_t CountGeneReadHuge(ogReadKeyMapping *pMap);
extern uint16_t CountGeneReadBlastSmall(ogReadKeyMapping *pMap);
extern uint16_t CountGeneReadBlastLarge(ogReadKeyMapping *pMap);
extern uint16_t CountGeneReadBlastHuge(ogReadKeyMapping *pMap);
extern uint16_t SweepAllPosSmall(ogReadKeyMapping *pMap);
extern uint16_t SweepAllPosLarge(ogReadKeyMapping *pMap);
extern uint16_t SweepAllPosHuge(ogReadKeyMapping *pMap);



unsigned long long getTotalSystemMemory()
{
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
}


int main(int argc, char *argv[]) {

    //testFastFiles(argv[1], 0);
    
    char samout1[512] = "samout1.sam";
    char samout2[512] = "samout2.sam";
    char cmdLine[4096]; // somewhere I read this is the maximum in linux
    
    long long totalRam = getTotalSystemMemory();
    fprintf(stderr, "ogMapper version %s\n", OGMAPPER_VERSION);
    fprintf(stderr, "Reported Memory=%lld, argc=%d\n", totalRam, argc);
    
    if (argc < 3) {  
        fprintf(stderr, "Usage: Indexing, Mapping, and Counting sequencing reads in genomes.\n");
        fprintf(stderr,"\n"); 
        fprintf(stderr,"INDEXING:\n"); 
        fprintf(stderr,"\n"); 
        fprintf(stderr, "    ogmapper index [-k <keysize>] [-g <guider>] [-e <encoding>]\n");
        fprintf(stderr, "        [-m 0/1] [-gtf <file>] [-ogx <ogx-prefix-file-name>]\n");
        fprintf(stderr, "        [-o <index-file-no-ext>] <genome-fasta.gz>\n");
        fprintf(stderr, "\n"); 
        fprintf(stderr, "        Prepare files for mapping (step 1) or counting (step 2).\n"); 
        fprintf(stderr, "           Encodings (-e option):\n");
        fprintf(stderr, "               BitwiseAT1GC0Encoding\n");
        fprintf(stderr, "               PlainEncoding\n");
        fprintf(stderr, "               GappedBitwiseAT1GC0Encoding\n");
        fprintf(stderr, "               SwapBitwiseAT1GC0Encoding\n");
        fprintf(stderr, "               HPCEncoding\n");
        fprintf(stderr, "           Guiders (-g option):\n");
        fprintf(stderr, "               StateMachineGuider:<state-file>\n");
        fprintf(stderr, "               DefaultGuider\n");
        fprintf(stderr, "               TupleGuider:<tuple-file>\n");
        fprintf(stderr, "           GTF/Counting (-gtx option used after index-GTF):\n");
        fprintf(stderr, "               Enables index for pseudo-counting mode.\n");
        fprintf(stderr, "               File name used to generate index-GTF.\n");
        fprintf(stderr, "               Must be provided for counting.\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "    ogmapper index-GTF [-o <index-names>] <.gtf[.gz]|.gff[.gz]> <genome.fa.gz>\n");
        fprintf(stderr,"\n"); 
        fprintf(stderr,"         Prepare files for pseudo-counting (step 1).\n"); 
        fprintf(stderr,"         Output:\n"); 
        fprintf(stderr,"         - <genome>-GTF.fq.gz file.\n"); 
        fprintf(stderr,"         - <gtf>.genes.ogx file.\n"); 
        fprintf(stderr,"         - <gtf>.exons.ogx file.\n"); 
        fprintf(stderr,"         - <gtf>.transcripts.ogx file.\n"); 
        fprintf(stderr, "\n");
        fprintf(stderr,"MAPPING:\n"); 
        fprintf(stderr,"\n"); 
        fprintf(stderr,"     ogmapper map [-k <keysize>] [-s <sched>] [-q <nQueue>] [-UC <0/1>]\n");  
        fprintf(stderr,"         [-t <nThreads>] [-p <0/1>] -i <index-file> [-d <n>]\n");  
        fprintf(stderr,"         [-maxreadlen <length>] [-maxreads <n>] [-startread <n>]\n");  
        fprintf(stderr,"         [-o <out.sam>|stdout] [-kseq <0/1>] [-unmapped 0/1/2/3/4] [-R <str>]\n"); 
        fprintf(stderr,"         [-1 <reads.gz>] [-2 <read-1.gz> <read-2.gz>]\n"); 
        fprintf(stderr,"\n"); 
        fprintf(stderr,"         Perform read mapping (step 2).\n"); 
        fprintf(stderr,"         NEEDS an index generated with ogmapper index (step 1).\n"); 
        fprintf(stderr, "\n");
        fprintf(stderr,"COUNTING:\n"); 
        fprintf(stderr,"\n"); 
        fprintf(stderr,"     ogmapper count [-k <keysize>] [-s <sched>] [-q <nQueue>] [-UC <0/1>]\n");  
        fprintf(stderr,"         [-t <nThreads>] [-p <0/1>] -i <complete-index-file> [-d <n>]\n");  
        fprintf(stderr,"         [-maxreadlen <length>] [-maxreads <n>] [-startread <n>]\n");  
        fprintf(stderr,"         [-genes <outfile>] [-exons <outfile>] [-transcripts <outfile>]\n");  
        fprintf(stderr,"         [-o <out.sam>|stdout] [-kseq <0/1>] [-unmapped 0/1/2/3/4] [-R <str>]\n"); 
        fprintf(stderr,"         [-1 <ogx> <reads.gz>] [-2 <ogx> <read-1.gz> <read-2.gz>]\n"); 
        fprintf(stderr,"\n"); 
        fprintf(stderr,"         Perform read pseudo-counting (step 3).\n"); 
        fprintf(stderr,"         NEEDS a indexes generated with ogmapper index-GTF (step 1).\n"); 
        fprintf(stderr,"         Also NEEDS an index generated with ogmapper index with -gtf option (step 2).\n"); 
        fprintf(stderr, "\n");
        fprintf(stderr, "OPTIONS:\n");
        fprintf(stderr, "        -o      Output file. Should appear before -1/-2 parameter.\n");
        fprintf(stderr, "                File may end with .gz but is quite slower.\n");
        fprintf(stderr, "        -k      Key size, in nt depending on the encoding.\n");
        fprintf(stderr, "        -m      Enables low memory access (-m 1) to save memory when indexing.\nDefault 0.\n");
        fprintf(stderr, "                Default to 0.\n");
        fprintf(stderr, "        -s      'schedule' functions to call for read mappings.\n");
        fprintf(stderr, "        -UC     Force uppercase reads.\n");
        fprintf(stderr, "        -kd     Specifies key distance.\n");
        fprintf(stderr, "        -q      Queue size for reads. Default 10000. Removal 0. Recommended >= 1000.\n");
        fprintf(stderr, "        -t      Threads used for processing.\n");
        fprintf(stderr, "        -f      Generate one output file per thread.\n");
        fprintf(stderr, "        -p      'Production' mode -p 1 (faster, default). Production mode 0\n");
        fprintf(stderr, "                designed to test scheduled mapping functions and times.\n");
        fprintf(stderr, "        -kseq n Use of the kseq.h library  (n=1, default) or customized (n=0) for\n");
        fprintf(stderr, "                reading reads. The customized reader can be faster in some systems.\n");
        fprintf(stderr, "        -gtf    Specify file for regions 'pseudo'-counting.\n");
        fprintf(stderr, "        -genes  Output file for gene counting.\n");
        fprintf(stderr, "        -exons  Output file for exon counting.\n");
        fprintf(stderr, "        -transcripts Output file for transcript counting.\n");
        fprintf(stderr, "        -keycount Sets minimum key counts.\n");
        fprintf(stderr, "        -scorecount Sets minimum score counts.\n");
        fprintf(stderr, "        -tie1   Sets tie 1 parameter.\n");
        fprintf(stderr, "        -prepare Do not prepare nor process reads. Used for checking reading times.\n");
        fprintf(stderr, "        -process Do not process reads. Used for checking reading and preparing times.\n");
        //fprintf(stderr, "        -b      Sets buffer size for reads (with 6%% auto growing).\n");
        fprintf(stderr, "        -d <n>   Delete key positions more frequent than <n>. This helps to avoid\n");
        fprintf(stderr, "                 wasting time searching over large populated keys.\n");
        fprintf(stderr, "        -unmapped Activate unmapping mode in which unmapped reads are stored in\n");
        fprintf(stderr, "                 corresponding files adding '.unmap' in RAW/TEXT/PLAIN format.\n");
        fprintf(stderr, "                 These reads will not be included in the .sam output.\n");
        fprintf(stderr, "                 This helps to use a second tool for unmapped reads.\n");
        fprintf(stderr, "                 The mode '-unmapped 1' specify that consider pair unmapped if any read\n");
        fprintf(stderr, "                 is unmapped, and both reads are not included in sam output.\n");
        fprintf(stderr, "                 The mode '-unmapped 2' needs both reads unmapped, otherwise both reads\n");
        fprintf(stderr, "                 will be included in sam output.\n");
        fprintf(stderr, "                 The mode '-unmapped 3' is similar to 1 but not generating .unmap file.\n");
        fprintf(stderr, "                 The mode '-unmapped 4' is similar to 2 but not generating .unmap files.\n");
        fprintf(stderr, "        -maxreads Maximum number of reads to be processed.\n");
        fprintf(stderr, "        -startread Starting read to be processed.\n");
        fprintf(stderr, "        -maxreadlen Read length used to estimate initial memory allocations.\n");
        fprintf(stderr, "        -R       Read group header to add in sam output. Within this, the ID:<id>\n");
        fprintf(stderr, "                 will be used in every read. Example: ’@RG\tID:foo\tSM:bar\tPL:ILLUMINA’.\n");
        fprintf(stderr, "                 This will add the above as header and the tag 'RG:Z:foo' to every read.\n");
        fprintf(stderr, "                 In the case that ID is numeric, 'RG:i:0' will be added instead.\n");
        fprintf(stderr, "        -bs <n>  Sets buffer size sam output file.\n");
        fprintf(stderr, "\n");
        return 1;  
    }

    ogIndex idx;
    uint16_t keySize = 24;
    int val;
    int k = 0;
    int i;
    char schedule[100];
    unsigned int threads;
    unsigned int queue;
    char upper = 0;
    char *pC = cmdLine;
    
    for (i=0; i < argc; i++) {
        pC += snprintf(pC, 4096 - (pC-cmdLine), "%s ", argv[i]);
        //fprintf(stderr, "Argument %d = [%s]\n", i, argv[i]);
    }
    *pC = 0;
    fprintf(stderr, "Command line [%s]\n", cmdLine);
   
    idx.registerKeyEncoding(new ogBitwiseAT1GC0Encoding());
    idx.registerKeyEncoding(new ogPlainEncoding());
    idx.registerKeyEncoding(new ogGappedBitwiseAT1GC0Encoding());
    idx.registerKeyEncoding(new ogSwapBitwiseAT1GC0Encoding());
    idx.registerKeyEncoding(new ogHomoPolymerCompressedEnc());
    
    idx.registerGuider(new ogStateMachineGuider());
    idx.registerGuider(new ogGuider());
    idx.registerGuider(new ogTupleGuider());
    
//    idx.registerMappingFunction('A', &ExactInApartSides);
//    idx.registerMappingFunction('E', &ExactInExtremes);
//    idx.registerMappingFunction('5', &Two5pKeys);    
//    idx.registerMappingFunction('d', &DefaultMatching5p);
//    idx.registerMappingFunction('M', &MapPosMatching);    
//    idx.registerMappingFunction('O', &NonOverlappingKeys);
//    idx.registerMappingFunction('3', &OverlappingBy3nt);
//    idx.registerMappingFunction('6', &OverlappingBy6nt);
//    idx.registerMappingFunction('9', &OverlappingBy9nt);
//    idx.registerMappingFunction('C', &OverlappingBy12nt);
//    idx.registerMappingFunction('F', &OverlappingBy15nt);
//    idx.registerMappingFunction('G', &OverlappingBy20nt);
    idx.registerMappingFunction('l', &LeftMost2KeysQuick);
    idx.registerMappingFunction('L', &LeftMost2KeysLong);
    idx.registerMappingFunction('Y', &LeftMost2KeysHuge);
    idx.registerMappingFunction('s', &Min2KeysSize);
    idx.registerMappingFunction('S', &Min2KeysLong);
    idx.registerMappingFunction('Z', &Min2KeysHuge);
    
    //Not working
    idx.registerMappingFunction('d', &DefaultMatchingSmall);   // Es la misma que D
    idx.registerMappingFunction('D', &DefaultMatchingLarge);   
    idx.registerMappingFunction('E', &DefaultMatchingHuge);   
    //Not working
    
    idx.registerMappingFunction('h', &HistogramMatchingSmall); 
    idx.registerMappingFunction('H', &HistogramMatchingLarge);   
    idx.registerMappingFunction('J', &HistogramMatchingHuge);   
    idx.registerMappingFunction('j', &HistogramMatchingAll);
    idx.registerMappingFunction('n', &NestedIntersectSmall);
    idx.registerMappingFunction('N', &NestedIntersectLarge);
    idx.registerMappingFunction('a', &MinApart2Keys);
    idx.registerMappingFunction('i', &AnalyzeAllPos);
    idx.registerMappingFunction('m', &MovWinIntersectSmall); // Hay que arreglar el Mov Window
    idx.registerMappingFunction('M', &MovWinIntersectLarge);
    idx.registerMappingFunction('w', &MovWinIntersectMicro);
    idx.registerMappingFunction('W', &MovWinIntersectHuge);
    idx.registerMappingFunction('P', &MapPosMatching);
    idx.registerMappingFunction('1', &HalfMatchingSmall);
    idx.registerMappingFunction('2', &HalfMatchingLarge);
    idx.registerMappingFunction('3', &HalfMatchingHuge);
    idx.registerMappingFunction('4', &ExtremeMatchingSmall);
    idx.registerMappingFunction('5', &ExtremeMatchingLarge);
    idx.registerMappingFunction('6', &ExtremeMatchingHuge);
    idx.registerMappingFunction('b', &BinnedHistogramMatchingSmall);
    idx.registerMappingFunction('B', &BinnedHistogramMatchingLarge);   
    idx.registerMappingFunction('V', &BinnedHistogramMatchingHuge);   
    idx.registerMappingFunction('v', &BinnedHistogramMatchingAll); 
    idx.registerMappingFunction('x', &blastLikeSmall);
    idx.registerMappingFunction('X', &blastLikeLarge);
    idx.registerMappingFunction('o', &blastLikeHuge);
    idx.registerMappingFunction('O', &blastLikeAll);
    idx.registerMappingFunction('u', &Intersect2R_Small);
    idx.registerMappingFunction('U', &Intersect2R_Large);
    idx.registerMappingFunction('c', &Intersect2R_Huge);
    idx.registerMappingFunction('C', &Intersect2R_All);
    //idx.registerMappingFunction('X', &blastLikeBinary);
    idx.registerMappingFunction('7', &CountGeneReadSmall);
    idx.registerMappingFunction('8', &CountGeneReadLarge);
    idx.registerMappingFunction('9', &CountGeneReadHuge);
    idx.registerMappingFunction('k', &CountGeneReadBlastSmall);
    idx.registerMappingFunction('K', &CountGeneReadBlastLarge);
    idx.registerMappingFunction('0', &CountGeneReadBlastHuge);
    idx.registerMappingFunction('q', &SweepAllPosSmall);
    idx.registerMappingFunction('Q', &SweepAllPosLarge);
    idx.registerMappingFunction('G', &SweepAllPosHuge);
    
    
    
    
    
    if (strcmp(argv[1], "checksequences") == 0) {
        gzFile zFile = gzopen(argv[2], "r");
        ogFastAQGZreader fqRd;
        fqRd.setFileReader(&zFile, 0, 0);
        ogFastAQsequence *faqs = (ogFastAQsequence *) calloc(1, sizeof(ogFastAQsequence));
        uint64_t nr = 0;
        while (fqRd.isEOF() == 0 && nr < 10000000) {
            fqRd.getNextSequenceOriginal(faqs);
            if (faqs->seq.l == 0) break; 
            if (++nr % 100000 == 0) fprintf(stderr, "%llu%c", nr, (nr % 500000 == 0 ? '\n' : ' '));
            fprintf(stderr, "Read=%llu; Lengths: Name=%d [%s], Comment=%d, Sequence=%d, Quality=%d\n", nr, faqs->name.l, faqs->name.s, faqs->comment.l, faqs->seq.l, faqs->qual.l);
        }
        free(faqs);
        gzclose(zFile);
        return 0;
    }
    
    if (strcmp(argv[1], "checkseqnaive") == 0) {
        gzFile zFile = gzopen(argv[2], "r");
        ogFastAQGZreader fqRd;
        fqRd.setFileReader(&zFile, 0, 0);
        ogFastAQsequence *faqs = (ogFastAQsequence *) calloc(1, sizeof(ogFastAQsequence));
        uint64_t nr = 0;
        while (fqRd.isEOF() == 0  && nr < 10000000) {
            fqRd.getNextSequenceNaive(faqs);
            if (faqs->seq.l == 0) break; 
            if (++nr % 100000 == 0) fprintf(stderr, "%llu%c", nr, (nr % 500000 == 0 ? '\n' : ' '));
            fprintf(stderr, "Read=%llu; Lengths: Name=%d [%s], Comment=%d, Sequence=%d, Quality=%d\n", nr, faqs->name.l, faqs->name.s, faqs->comment.l, faqs->seq.l, faqs->qual.l);
        }
        free(faqs);
        gzclose(zFile);
        return 0;
    }
    
    if (strcmp(argv[1], "checkseqbytype") == 0) {
        gzFile zFile = gzopen(argv[2], "r");
        ogFastAQGZreader fqRd;
        fqRd.setFileReader(&zFile, 0, 0);
        ogFastAQsequence *faqs = (ogFastAQsequence *) calloc(1, sizeof(ogFastAQsequence));
        uint64_t nr = 0;
        while (fqRd.isEOF() == 0  && nr < 10000000) {
            fqRd.getNextSequenceByType(faqs);
            if (faqs->seq.l == 0) break; 
            if (++nr % 100000 == 0) fprintf(stderr, "%llu%c", nr, (nr % 500000 == 0 ? '\n' : ' '));
            //if (faqs->name.l > 100 || faqs->comment.l > 100 || faqs->seq.l > 200 || faqs->qual.l > 200) {
            //    fprintf(stderr, "Read=%llu; Lengths: Name=%d [%s], Comment=%d, Sequence=%d, Quality=%d\n", nr, faqs->name.l, faqs->name.s, faqs->comment.l, faqs->seq.l, faqs->qual.l);
            //}
            fprintf(stderr, "Read=%llu; Lengths: Name=%d [%s], Comment=%d, Sequence=%d, Quality=%d\n", nr, faqs->name.l, faqs->name.s, faqs->comment.l, faqs->seq.l, faqs->qual.l);
        }
        free(faqs);
        gzclose(zFile);
        return 0;
    }
    
    if (strcmp(argv[1], "checkkseq") == 0) {
        idx.checkKSEQ(argv[2]);
    }
    
    if (strcmp(argv[1], "getChr20_22") == 0) {
        idx.getChr20_22(argv[2+k]);
    }
    
    if (strcmp(argv[1], "testRep") == 0) {
        if (strcmp(argv[2+k], "-maxreads") == 0) {
            uint64_t mr;
            sscanf(argv[3+k], "%llu", &mr);
            idx.setMaxReads(mr);
            k+=2;
        }
        uint16_t repes = 2;
        if (strcmp(argv[2+k], "-rep") == 0) {
            sscanf(argv[3+k], "%hu", &repes);
            k+=2;
        }
        fprintf(stderr, "==== TESTING File: %s ====\n",argv[2+k]);
        idx.testRepes(argv[2+k], repes);
        fprintf(stderr, "==== TESTING %s Done ====\n", argv[2+k]);

    }

    if (strcmp(argv[1], "cigar") == 0) {
        char pCIGAR[1000];
        char *pC = (char *)&pCIGAR;
        sscanf(argv[2], "%s", pC);
        fprintf(stderr, "Read CIGAR: /%s/\n", pCIGAR);
        ogCigarOperations cigarro;
        fprintf(stderr, "CIGAR operations added: %hu\n", cigarro.pushOperationsFromCIGAR(pCIGAR));
        fprintf(stderr, "CIGAR referred size: %u\n", cigarro.totalSizeConsuming());
        fprintf(stderr, "CIGAR sum sizes: %u\n", cigarro.totalSize());
        cigarro.reduceOperations();
        cigarro.buildCIGAR(pCIGAR, 1000);
        fprintf(stderr, "Reduced CIGAR: /%s/\n",pCIGAR);        
        fprintf(stderr, "CIGAR referred size: %u\n", cigarro.totalSizeConsuming());
        fprintf(stderr, "CIGAR sum sizes: %u\n", cigarro.totalSize());
    }

    if (strcmp(argv[1], "testStateMachine") == 0) {
        fprintf(stderr, "testStateMachine <states file> <seq file> <nLen>\n");
        fprintf(stderr, "States file:%s\n", argv[2]);
        fprintf(stderr, "Seq file:%s\n", argv[3]);
        fprintf(stderr, "len:%s\n", argv[4]);
        uint32_t maxLenNo;
        sscanf(argv[4], "%u", &maxLenNo);
        idx.testStateMachine(argv[3], argv[2], maxLenNo);
    }
    
    if (strcmp(argv[1], "index-GTF") == 0) {
        k = 2;
        char *gtfNames = NULL;
        if (strcmp(argv[k], "-o") == 0) { gtfNames = argv[k+1]; k += 2; }
        idx.indexFromGTF(argv[k], argv[k+1], gtfNames);
    }
    
    if (strcmp(argv[1], "index") == 0) {
        if (totalRam < 20179869184) { // This value was obtain in a 16GB Mac plus ~3GB
            fprintf(stderr, "\n\n===============================================================\n");
            fprintf(stderr, "LOW MEMORY WARNING! INDEXING USE ~5 TIMES GENOME SIZE MEMORY\n");
            fprintf(stderr, "IT IS RECOMMENDED TO USE -m 1 TO AVOID LARGE INDEXING TIMES\n");
            fprintf(stderr, "FOR >2GB GENOMES\n");
            fprintf(stderr, "\n===============================================================\n");
        }
        *samout1 = 0;
        for (k=2; k < argc; k += 1) {
            if (strcmp(argv[k], "-k") == 0) {
                sscanf(argv[k+1], "%hd", &keySize);
                k++;
            }
            if (strcmp(argv[k], "-e") == 0) {
                idx.selectKeyEncoding(argv[k+1]);
                k++;
            }
            if (strcmp(argv[k], "-g") == 0) {
                idx.selectGuider(argv[k+1]);
                k++;
            }
            if (strcmp(argv[k], "-m") == 0) {
                idx.setMemoryMode(argv[k+1][0]-'0');
                k++;
            }
            if (strcmp(argv[k], "-gtf") == 0) {
                idx.isExonMode = 1;
                strncpy(idx.countInfo.GTFfile, argv[k+1], MAX_FILENAME);
                idx.countInfo.GTFfile[strlen(argv[k+1])] = 0;
                k++;
            }
            if (strcmp(argv[k], "-ogx") == 0) {
                strncpy(idx.countInfo.OGXfile, argv[k+1], MAX_FILENAME);
                k++;
            }
            if (strcmp(argv[k], "-o") == 0) {
                strncpy(samout1, argv[k+1], 512);
                strncat(samout1, ".ogi", 10);
                k++;
            }
        }
        fprintf(stderr, "==== INDEXING File: %s ====\n",argv[argc-1]);
        if (idx.get_file_size(argv[argc-1]) == -1) {
            fprintf(stderr, "==== ERROR FILE %s DOES NOT EXIST? ====\n", argv[argc-1]);            
        } else {
            idx.generate(argv[argc-1], keySize, samout1, argc, argv);
            fprintf(stderr, "==== INDEXING %s Done ====\n", argv[argc-1]);
        }

    }
    
    if ((strcmp(argv[1], "map") == 0) || (strcmp(argv[1], "count") == 0)) {
        char isMap = (strcmp(argv[1], "map") == 0);
        char schedOptions[6] = "x4s"; //"x4l" "slSL";
        idx.setSchedule(schedOptions);
        
        if (argc < 4) {  
            fprintf(stderr, "Usage: %s map [-k <num>] [-maxreadlen <num>] [-maxreads <num>] [-startread <num>] [-s <sched>] [-t threads] [-q maxQueue] [-UC <0/1>] [-p <0/1>] -i <genome-encoding-no-extension> [-d <n>] [[-o <output.sam>|stdout] [-1 <read-1.gz>] [-2 <read-1.gz> <read-2.gz>]]...\n", argv[0]);  
            return 1;  
        }

        idx.pMapParams->cmdLine = cmdLine;        
        
        char output = 0;
        
        for (k=2; k < argc; k += 1) { // argv[k][0] == '-'
            //fprintf(stderr, "Processing parameter k=%d : %s\n", k, argv[k]); fflush(stderr);
            if (strcmp(argv[k], "-k") == 0) {
                sscanf(argv[k+1], "%hd", &keySize);
                k++;
                
            } else if (strcmp(argv[k], "-maxreadlen") == 0) {
                uint16_t mrl;
                sscanf(argv[k+1], "%hu", &mrl);
                idx.setMaxReadLength(mrl);
                k++;
                
            } else if (strcmp(argv[k], "-maxreads") == 0) {
                uint64_t mr;
                sscanf(argv[k+1], "%llu", &mr);
                idx.setMaxReads(mr);
                k++;
                
            } else if (strcmp(argv[k], "-startread") == 0) {
                uint64_t sr;
                sscanf(argv[k+1], "%llu", &sr);
                idx.setStartingRead(sr);
                k++;
                
            } else if (strcmp(argv[k], "-kseq") == 0) {
                idx.setKSEQmode(argv[k+1][0]-'0');
                k++;
                
            } else if (strcmp(argv[k], "-s") == 0) {
                sscanf(argv[k+1], "%s", schedule);
                idx.setSchedule(schedule);
                k++;
                
            } else if (strcmp(argv[k], "-t") == 0) {
                sscanf(argv[k+1], "%u", &threads);
                if (threads >= 0 && threads < 100000) idx.setThreads(threads);
                k++;
                
            } else if (strcmp(argv[k], "-q") == 0) {
                sscanf(argv[k+1], "%u", &queue);
                if (queue >= 0) idx.setQueue(queue);
                k++;
                
            } else if (strcmp(argv[k], "-UC") == 0) {
                idx.setUpperCase(argv[k+1][0]-'0');
                k++;
                
            } else if (strcmp(argv[k], "-p") == 0) {
                idx.setProductionMode(argv[k+1][0]-'0');
                k++;
                
            } else if (strcmp(argv[k], "-i") == 0) {
                fprintf(stderr, "==== Loading Index File : %s ====\n",argv[k+1]);
                idx.load(argv[k+1], &keySize);
                //fprintf(stderr, "==== LOADED Index File : %s ====\n",argv[k+1]);
                k++;
                
            } else if (strcmp(argv[k], "-o") == 0) {
                if (output == 1) idx.closeOutputSAM();
                output = 0;
                idx.pMapParams->outputSAMcount = 1;
                //idx.openOutputSAM(argv[k+1], NULL, NULL);
                //output = 1;
                strncpy(samout1, argv[k+1], 512);
                strncpy(samout2, argv[k+1], 512);
                k++;

            } else if (strcmp(argv[k], "-kd") == 0) {
                sscanf(argv[k+1], "%u", &idx.pMapParams->MAX_KEY_DISTANCE);
                k++;

            } else if (strcmp(argv[k], "-ecount") == 0) {
                idx.pMapParams->countExons = 1;

            } else if (strcmp(argv[k], "-exons") == 0) {
                // genes is always outputed to "genes-counts.txt"
                idx.pMapParams->countExons = 1;
                strncpy(idx.exonsCountFile, argv[k+1], sizeof(idx.exonsCountFile));
                k++;

            } else if (strcmp(argv[k], "-tcount") == 0) {
                idx.pMapParams->countTranscripts = 1;

            } else if (strcmp(argv[k], "-transcripts") == 0) {
                // genes is always outputed to "transcripts-counts.txt"
                idx.pMapParams->countTranscripts = 1;
                strncpy(idx.transcriptsCountFile, argv[k+1], sizeof(idx.transcriptsCountFile));
                k++;

            } else if (strcmp(argv[k], "-genes") == 0) {
                // genes is always outputed to "genes-counts.txt"
                strncpy(idx.genesCountFile, argv[k+1], sizeof(idx.genesCountFile));
                k++;

            } else if (strcmp(argv[k], "-keycount") == 0) {
                sscanf(argv[k+1], "%hu", &idx.pMapParams->MIN_KEY_COUNTS);
                k++;

            } else if (strcmp(argv[k], "-scorecount") == 0) {
                sscanf(argv[k+1], "%hhu", &idx.pMapParams->MIN_SCORE_COUNT);
                k++;

            } else if (strcmp(argv[k], "-tie1") == 0) {
                sscanf(argv[k+1], "%hhu", &idx.pMapParams->BREAK_TIE_COUNT1);
                k++;
                
            } else if (strcmp(argv[k], "-f") == 0) {
                idx.filePerThread = 1;
                
            } else if (strcmp(argv[k], "-bs") == 0) {
                uint32_t bufSize = 0;
                sscanf(argv[k+1], "%u", &bufSize);
                ogSamWriter::samFileBuffer = bufSize;
                k++;
                
            } else if (strcmp(argv[k], "-d") == 0) {
                uint32_t removeSizes = 0;
                sscanf(argv[k+1], "%u", &removeSizes);
                k++;
                if (idx.pKeys == NULL) {
                    fprintf(stderr, "-d option should be specified after index (-i).\n");
                } else {
                    if (removeSizes > 0) {
                        fprintf(stderr, "))) Removing genome positions for keys having >= %u ... ", removeSizes);
                        removeSizes = idx.pKeys->removeKeySizesLargerThan(removeSizes);
                        fprintf(stderr, "%u keys removed (((\n", removeSizes);
                    }
                }
                
            } else if (strcmp(argv[k], "-unmapped") == 0) {
                idx.unmappedMode = *(argv[k+1]) - '0';
                k++;
                
            } else if (strcmp(argv[k], "-R") == 0) {
                char *pd = idx.readGroupHeader;
                for (char *s = argv[k+1]; *s != 0; s++) {
                    if (*s == '\\' && *(s+1) == 't') {
                        *pd++ = '\t';
                        s++;
                    } else {
                        *pd++ = *s;
                    }
                }
                *pd = 0;
                //strncat(idx.readGroupHeader, argv[k+1], MAX_FILENAME);
                k++;
                
            } else if (strcmp(argv[k], "-1") == 0) {
                if (isMap) {
                    if (output == 0) { idx.openOutputSAM(samout1, argv[k+1], NULL); output = 1; }
                    fprintf(stderr, "==== MAPPING Single: %s ====\n",argv[k+1]);
                    idx.map(argv[k+1], NULL);
                    fprintf(stderr, "\n==== MAPPED : %s ====\n",argv[k+1]);
                    k++;
                } else {
                    if (output == 0 && idx.pMapParams->outputSAMcount) { idx.openOutputSAM(samout1, argv[k+2], NULL); output = 1; }
                    fprintf(stderr, "==== COUNTING Single %s w/ogx %s ====\n",argv[k+2],argv[k+1]);
                    idx.count(argv[k+1], argv[k+2], NULL);
                    fprintf(stderr, "\n==== COUNTED : %s ====\n",argv[k+2]);
                    k += 2;
                }
                
            } else if (strcmp(argv[k], "-2") == 0) {
                if (isMap) {
                    if (output == 0) { idx.openOutputSAM(samout2,argv[k+1],argv[k+2]); output = 1; }
                    fprintf(stderr, "==== MAPPING Paired: %s + %s ====\n",argv[k+1],argv[k+2]);
                    idx.map(argv[k+1], argv[k+2]);
                    fprintf(stderr, "\n==== MAPPED : %s + %s ====\n",argv[k+1],argv[k+2]);
                    k += 3;
                } else {
                    if (output == 0 && idx.pMapParams->outputSAMcount) { idx.openOutputSAM(samout2,argv[k+2],argv[k+3]); output = 1; }
                    fprintf(stderr, "==== COUNTING Paired %s + %s w/ogx %s ====\n",argv[k+2],argv[k+3],argv[k+1]);
                    idx.count(argv[k+1], argv[k+2], argv[k+3]);
                    fprintf(stderr, "\n==== COUNTED : %s + %s ====\n",argv[k+2],argv[k+3]);
                    k += 3;
                }
                
            } else if (strcmp(argv[k], "-prepare") == 0) {
                idx.doPrepare = 0;
                
            } else if (strcmp(argv[k], "-process") == 0) {
                idx.doProcess = 0;
                
//            } else if (strcmp(argv[k], "-b") == 0) {
//                sscanf(argv[k+1], "%u", &idx.maxReadsBuffer);
//                k++;
            } else {
                fprintf(stderr, "OPTION [%s] NOT RECOGNIZED\n.QUIT.", argv[k]);
                return 0;
            }
        }
        if (output == 1) idx.closeOutputSAM();
                
    }

    fprintf(stderr, "**** PROGRAM ENDING ****\n");
    return 0;
}

