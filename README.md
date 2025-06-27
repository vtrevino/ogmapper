# ogmapper
A Fast and Light Genomic Mapper for short reads.

ogMapper generates a tiny index file and memory operations, typically as a fraction of the genome size. 
For example, 1.7 GiB for T2T human genome v2 plus 0.7 GiB for the compressed genome (file ~ 2.4 GiB) and similar amounts of memory in run time, irrespective of the number of threads. This is for 24-bit keys.
The common process starts by creating an index file (*.ogi), then mapping short reads (in pairs or not) to generate .sam files.

ogmapper is written in c++, originally developed in NetBeans 13 on an Intel Mac. I provide binaries for selected operating systems and source files for compilation.

# Binaries
- MacOS/Intel : <a href="ogmapper-mac-intel-clang">ogmapper-mac-intel-clang</a> or <a href="ogmapper-mac-intel-gnu">ogmapper-mac-intel-gnu</a>
- MacOS/M : <a href="ogmapper-mac-m-clang">ogmapper-mac-m-clang</a>
- Linux/x86 : <a href="ogmapper-linux-x86">ogmapper-linux-x86</a>
- Windows :(

# Installation
Download the binary and test it without any arguments. If it works, it should display the options as shown in the section [Running ogMapper] below. If you want to build from source, be my guest and follow the next sections.

# Compilation
If one of the above binaries runs in your system, I recommend using it if they do not cause problems. Still, if you prefer to build it yourself, two components need compilation, the WFA library and ogMapper. Once WFA library has been built, we may proceed to compile ogMapper.

## Compiling WFA
ogMapper uses the WFA2 library to perform alignment operations when needed. I provide the latest version used to build the binaries shown above. Users may opt to download the latest WFA version from https://github.com/smarco/WFA2-lib. You may follow the original WFA2 instructions or follow the steps below as a guide. The goal is to build a library file suitable for ogMapper (libwfacpp.a) that needs to be placed in the /lib folder to be able to compile ogMapper.

### Compiling WFA2 from .zip
- MacOSX: Assuming Xcode command line tools and brew are already installed, download <a href="WFA2-lib-main.zip">WFA2-lib-main.zip</a> and then 

        # Install llvm and libomp
        brew install llvm libomp
  
        # Make WFA2
        unzip WFA2-lib-main.zip
        cd WFA2-lib-main
        make clean
        make build
        # This should generate libwfacpp.a in the lib/ folder
  
        # Copy libwfacpp.a into ogMapper folder under lib/ 
        cp lib/libwfacpp.a [your ogMapper/lib]
  
If you have problems with "VERSION" at compile time, rename it to "__VERSION".
If you have problems with LIBOMP you may need one of the compile flags: "-lomp", "-fopenmp", "-openmp", or "-gomp", depending on the library installed.

- Linux: Assuming c/c++ compilers are installed, download <a href="WFA2-lib-main.zip">WFA2-lib-main.zip</a> and then 
tools/align_benchmark/Makefile

        # 
        unzip WFA2-lib-main.zip
        cd WFA2-lib-main
        make clean
  
        # Edit all Makefiles replacing -lomp flag to -fopenmp (or -openmp, -gomp)
        # Files : Makefile,  wavefront/Makefile, tools/align_benchmark/Makefile, examples/Makefile
        # This "make" should generate libwfacpp.a in the lib/ folder
        make build
        cp lib/libwfacpp.a [your ogMapper/lib]
        # Copy libwfacpp.a into ogMapper folder under lib/ 

If you have problems with "VERSION" at compile time, rename it to "__VERSION".
If you have problems with LIBOMP you may need one of the compile flags: "-lomp", "-fopenmp", "-openmp", or "-gomp", depending on the library installed.


### Compiling WFA2 from https://github.com/smarco/WFA2-lib
Follow the instructions provided in WFA2-lib.

## Compiling ogMapper

### Compiling ogMapper on MacOS
Assuming Xcode command line tools and brew are already installed and WFA2 library is copied in the ogMapper/lib folder:

        cd ogMapper
        make clean
        make Release
        # This should generate a binary in dist/Release

### Compiling ogMapper on Linux
Assuming c/c++ compilers are installed and WFA2 library is copied in the ogMapper/lib folder:

        cd ogMapper
        make clean
        make Release
        # This should generate a binary in dist/Release

# Introduction to Encodings, Guiders, Keys, and Mapping functions
The first step to map or count reads is generating an index of the reference genome. For this task, ogMapper introduced some concepts and also used novel explorative functions. ogMapper does not index all possible DNA subsequences. It first sweeps the DNA sequence until a specific sequence pattern is found; this pattern is known as a "guide". There are two guiders implemented up to now. Once a guide is found, the following DNA sequence of length k (provided by -k argument) is used to build a key, which will be indexed. The key is a binary representation of the DNA sequence, which is provided by an encoding function. There are four encoding functions implemented so far. Once a guide-key is processed, the DNA sweeping continues after the last guide until the whole sequence has been analyzed. In index-time all keys generated from all chromosomes are included in the index file as seed keys. In mapping-time, the keys generated from reads are used to search possible matches in the index. If the guide pattern is not found after L nt, an exception ocurrs indexing following DNA sequences as keys and skipping S nt. This process is better documented in the manuscript. The L parameter can be specified in guiders files as shown below.

For mapping, ogMapper runs a list of mapping functions until one of them maps the reads properly. There are many mapping functions. See the Mapping sections for details.

### Guiders (-g option)
- TupleGuider: The tuple guider is a simple map of all possible n-nt combinations specifying which of them will be indexed and which will not. ogMapper includes a folder of pre-configured Guider files (guiders/ folder). The file format needs first 4 lines for parameters then the needed lines for tuples.

        # Next lines start with "small extension name", maxLenNoPattern, isSymmetric, GuideLen
        XYXt
        24
        1
        3
        # Next lines start tuples and guiders
        #Tuple	isGuider	NT_to_advance	Comment
        AAA	0	1	x
        AAC	0	1	x
        AAG	0	1	x
        AAT	0	1	x
        ACA	1	3	Index
        ACC	0	1	x
        ...
        

To use the TuppleGuider use the option -g TupleGuider:&lt;tuple-file&gt;.

I recommend using tuple guiders since they are easier to understand and maintain.
  
- StateMachineGuider: State machine guiders are used when not simple mappings can be specified using Tuples. For instance, if you want a variable-length guiders, it cannot be done with Tuples.
  
To use the StateMachineGuider use the option -g StateMachineGuider:&lt;state-file&gt;.

### Encodings
The encoding transforms the DNA sequence to a binary key. I recommend BitWise and Plain.

The **BitwiseAT1GC0Encoding** uses 1 bit per nt, transforming A or T/U to 1 and G or C to 0. Any other letter is treated as A. Case-insensitive. So, -k 24 will use 24 nt to generate keys of 24 bits for a total of 16,777,216 different keys. Increasing k would have an important impact on the index size and memory needed.

The **PlainEncoding** uses 2 bits per nt, transforming A to 00, C to 01, G to 10, and T/U to 11. Any other letter is treated as A. So, -k 12 will use 12 nt generating keys of 24 bits. 

The **GappedBitwiseAT1GC0Encoding** is similar to BitwiseAT1GC0Encoding but the nt used for indexing are chosen as [left][gap][right]=k where [left] and [right] are estimated by k/3. So, here -k=36 is equivalent in bits to -k=24 in BitwiseAT1GC0Encoding.

The **SwapBitwiseAT1GC0Encoding** is similar to BitwiseAT1GC0Encoding but it uses non-continuous nt indexing, one nt and skipping one until k is reached.

The **HPCEncoding** (homo-polymer compressed) ignores consecutive repetitions of a nucleotide and then uses 'plain' encoding as the **PlainEncoding** above.


### Mapping Functions
The schedule option (-s) specifies the order of the functions used to attempt mapping each read. In ogMapper, the algorithm for mapping runs a list of mapping functions until one of them maps the reads properly. There are many mapping functions. The list is controlled by the argument -s (schedule). The default is x4s where 'x','4', and 's' are mapping functions. The '4' function is used only when 'x' is unable to map, and so on the next functions. 

The list of functions are:

x – blast_like (schedules x, X, o, O): Check all positions for each key. Process in “x”: for each key in low-to-high size-ordered keys, pick the position of maximum matching. Stop if the identity is higher than 98%. Forward keys and reverse keys are tested. Only use the keys with less than LowKeyCountLimit/2 (LowKeyCountLimit = 53 for CHM13v2 using BitwiseAT1GC0Encoding and XXX guide). All positions with the same best score are added. “X” uses 99% identity and keys sizes < LowKeyCountLimit. “o” uses 100% identity and key sizes < HighKeyCountLimit = 79 for CHM13v2 and BitwiseAT1GC0Encoding on XXX guide. “O” uses 100% identity and all key sizes. The speed order is xXoO while the accuracy is OoXx. In general, highly accurate and fast for most today small reads. Not good for indels around keys.

4 – extreme keys (schedules 4, 5, 6): Check intersections of positions from keys in the left part of the read with those positions of the right part of the read. Process for “4”: for each key within the first half of keys, check intersections with positions +/- 3 nt from keys of the right part of the read. Intersections are added. Only key sizes <  LowKeyCountLimit are considered. Forward keys and reverse keys are tested. “5” uses key sizes < HighKeyCountLimit. “6” uses key sizes < HighKeyCountLimit * 10. Small indels are well handled.

s – two consistent keys (schedules s, S, Z): Intersect positions of two non-overlapping keys. Process for “s”: for unmarked keys, take one key and search for a non-overlapping unmarked key. If they have shared a position (+/- 3 indels), add it to possible matches; otherwise, mark them to avoid re-check. Only key sizes <  LowKeyCountLimit. Forward keys and reverse keys are tested. “S” uses key sizes < HighKeyCountLimit while “Z” uses 10 * HighKeyCountLimit. 

l – left consistent keys (schedules l, L, Y): Intersect positions of two not close keys (> 10 bp) starting from the leftmost keys. Process for “l”: starting from the leftmost key, take a key and get the next not close key, add intersecting positions (+/-3 indels). Forward keys and reverse keys are tested. Only key sizes <  LowKeyCountLimit are tested. “L” uses key sizes < HighKeyCountLimit while “Y” uses < 10 * HighKeyCountLimit.

h – (schedules h, H, j, J): Sorted-size Histogram: An initial 1000-bin histogram is built from key positions. If a position re-occurs in a bin, a second pseudo-histogram is built also on different positions. Forward keys and reverse keys are tested. “h” uses keys sizes <  LowKeyCountLimit, “H” < HighKeyCountLimit, “j” < 10 * HighKeyCountLimit, and “J” uses all keys.

n – (schedules n, N): Nested intersect positions: Intersect positions of a key's “i” and “i+2” considering small indels (+/- 3 nt). Forward keys and reverse keys are tested. “n” uses key sizes <  LowKeyCountLimit while “N” uses HighKeyCountLimit.

a – (a): MinApart2Keys : *pending description*

i – (i): *pending description*

m – (m, M, w, W): *pending description*

P – (P) : *pending description*

1 – (1, 2, 3): *pending description*

b – (b, B, V, v): *pending description*

u – (u, U, c, C): *pending description*

7 – (7, 8, 9): *pending description*

k – (k, K, 0): *pending description*

q – (q, Q, G): *pending description*

If more than 1 position is added to the pool of possibilities or raw alignment is lower than 90%, alignment is performed and the best positions are chosen (unless rules for alternative alignments apply).


# Running ogMapper
The typical operations consist of creating an index file and mapping reads.

### ogMapper options
Running ogMapper without any arguments shows:

      ogMapper version v0.5.0-20-Jan-2025
      Reported Memory=17179869184, argc=1
      Usage: Indexing, Mapping, and Counting sequencing reads in genomes.
      
      INDEXING:
      
          ogmapper index [-k <keysize>] [-g <guider>] [-e <encoding>]
              [-m 0/1] [-gtf <file>] [-ogx <ogx-prefix-file-name>]
              [-o <index-file-no-ext>] <genome-fasta.gz>
      
              Prepare files for mapping (step 1) or counting (step 2).
                 Encodings (-e option):
                     BitwiseAT1GC0Encoding
                     PlainEncoding
                     GappedBitwiseAT1GC0Encoding
                     SwapBitwiseAT1GC0Encoding
                     HPCEncoding
                 Guiders (-g option):
                     StateMachineGuider:<state-file>
                     DefaultGuider
                     TupleGuider:<tuple-file>
                 GTF/Counting (-gtx option used after index-GTF):
                     Enables index for pseudo-counting mode.
                     File name used to generate index-GTF.
                     Must be provided for counting.
      
          ogmapper index-GTF [-o <index-names>] <.gtf[.gz]|.gff[.gz]> <genome.fa.gz>
      
               Prepare files for pseudo-counting (step 1).
               Output:
               - <genome>-GTF.fq.gz file.
               - <gtf>.genes.ogx file.
               - <gtf>.exons.ogx file.
               - <gtf>.transcripts.ogx file.
      
      MAPPING:
      
           ogmapper map [-k <keysize>] [-s <sched>] [-q <nQueue>] [-UC <0/1>]
               [-t <nThreads>] [-p <0/1>] -i <index-file> [-d <n>]
               [-maxreadlen <length>] [-maxreads <n>] [-startread <n>]
               [-o <out.sam>|stdout] [-kseq <0/1>] [-unmapped 0/1/2/3/4] [-R <str>]
               [-1 <reads.gz>] [-2 <read-1.gz> <read-2.gz>]
      
               Perform read mapping (step 2).
               NEEDS an index generated with ogmapper index (step 1).
      
      COUNTING:
      
           ogmapper count [-k <keysize>] [-s <sched>] [-q <nQueue>] [-UC <0/1>]
               [-t <nThreads>] [-p <0/1>] -i <complete-index-file> [-d <n>]
               [-maxreadlen <length>] [-maxreads <n>] [-startread <n>]
               [-genes <outfile>] [-exons <outfile>] [-transcripts <outfile>]
               [-o <out.sam>|stdout] [-kseq <0/1>] [-unmapped 0/1/2/3/4] [-R <str>]
               [-1 <ogx> <reads.gz>] [-2 <ogx> <read-1.gz> <read-2.gz>]
      
               Perform read pseudo-counting (step 3).
               NEEDS a indexes generated with ogmapper index-GTF (step 1).
               Also NEEDS an index generated with ogmapper index with -gtf option (step 2).
      
      OPTIONS:
              -o      Output file. Should appear before -1/-2 parameter.
                      File may end with .gz but is quite slower.
              -k      Key size, in nt depending on the encoding.
              -m      Enables low memory access (-m 1) to save memory when indexing.
                      Default to 0.
              -s      'schedule' functions to call for read mappings. Default 'x4s'.
              -UC     Force uppercase reads.
              -kd     Specifies key distance.
              -q      Queue size for reads. Default 10000. Removal 0. Recommended >= 1000.
              -t      Threads used for processing.
              -f      Generate one output file per thread.
              -p      'Production' mode -p 1 (faster, default). Production mode 0
                      designed to test scheduled mapping functions and times.
              -kseq n Use of the kseq.h library  (n=1, default) or customized (n=0) for
                      reading reads. The customized reader can be faster in some systems.
              -gtf    Specify file for regions 'pseudo'-counting.
              -genes  Output file for gene counting.
              -exons  Output file for exon counting.
              -transcripts Output file for transcript counting.
              -keycount Sets minimum key counts.
              -scorecount Sets minimum score counts.
              -tie1   Sets tie 1 parameter.
              -prepare Do not prepare nor process reads. Used for checking reading times.
              -process Do not process reads. Used for checking reading and preparing times.
              -d <n>   Delete key positions more frequent than <n>. This helps to avoid
                       wasting time searching over large populated keys.
              -unmapped Activate unmapping mode in which unmapped reads are stored in
                       corresponding files adding '.unmap' in RAW/TEXT/PLAIN format.
                       These reads will not be included in the .sam output.
                       This helps to use a second tool for unmapped reads.
                       The mode '-unmapped 1' specify that consider pair unmapped if any read
                       is unmapped, and both reads are not included in sam output.
                       The mode '-unmapped 2' needs both reads unmapped, otherwise both reads
                       will be included in sam output.
                       The mode '-unmapped 3' is similar to 1 but not generating .unmap file.
                       The mode '-unmapped 4' is similar to 2 but not generating .unmap files.
              -maxreads Maximum number of reads to be processed.
              -startread Starting read to be processed.
              -maxreadlen Read length used to estimate initial memory allocations.
              -R       Read group header to add in sam output. Within this, the ID:<id>
                       will be used in every read. Example: ’@RG	ID:foo	SM:bar	PL:ILLUMINA’.
                       This will add the above as header and the tag 'RG:Z:foo' to every read.
                       In the case that ID is numeric, 'RG:i:0' will be added instead.
              -bs <n>  Sets buffer size sam output file.


## Processing DNA reads

### Index Generation for DNA
For indexing DNA the valid options are:

      ogmapper index [-k <keysize>] [-g <guider>] [-e <encoding>]
        [-m 0/1] [-o <index-file-no-ext>] <genome-fasta.gz>

The guiders/ folder contains the state/ and tuple/ folders, which include pre-configured guiders files. Copy the file needed before running the indexing command, for example, the XXX-Tuple-og.txt file used below.

A typical run:

    ogmapper index -k 24 -g TupleGuider:XXX-Tuple-og.txt -e BitwiseAT1GC0Encoding -o chm13v2-XXX-BW chm13v2.0.fa.gz 

The above command will generate the index file "chm13v2-XXX-BW.ogi" of ~2.7 GiB. <a href="https://github.com/marbl/CHM13" target="_blank">T2T-CHM13</a> publish <a href="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz">chmv13v2.0.fa.gz</a>. There are also versions from NCBI. The -g option specifies a TupleGuider. The file XXX-Tuple-og.txt must be available in the same folder than data. This file is included in ogMapper in the guiders/ folder. The -e option specifies the 1-bit-1-nt encoding. See Encodings and Guiders sections above for details.

## Mapping short DNA reads
For mapping, a typical run looks like:

        ogmapper map -s x4s -t 8 -p 1 -i chm13v2-XXX-BW.ogi -o og.sam -2 U0a_R1.fq.gz U0a_R2.fq.gz

Most parameters are self-explanatory in the default output (without arguments). See the section ogMapper options (above). See also the Mapping Functions sections for details in the -s argument.

## Processing RNA reads (for pseudo-counts)
Overall, using ogMapper, one can pseudo-quantify the gene expression from RNA-Seq experiments using the fast mapping. For this, the sequence of exons and transcripts is extracted from DNA to generate fasta files with the sequence of each transcript or exon. Then, these sequences are indexed, and RNA reads are mapped to them. This section explains how to use ogMapper to estimate fast pseudo-counts from RNA-Seq experiments.

#### Step 1: Generate GTF indexes
This step extracts the sequences for genes, exons, and transcripts from a fasta genome and gtf/gff files.

Syntax:

      ogmapper index-GTF [-o <index-names>] <.gtf[.gz]|.gff[.gz]> <genome.fa.gz>

Output:
1) &lt;genome&gt;-GTF.fq.gz file.
2) &lt;genome&gt;-GTF.exons.ogx file.
3) &lt;genome&gt;-GTF.transcripts.ogx file.

Example:

      ogmapper index-GTF -o chm13v2-GTF GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz chm13v2.0.fa.gz

Links:
1) <a href="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz">.gtf.gz</a> file can be downloaded from NCBI.
2) <a href="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz">chmv13v2.0.fa.gz</a> can be downloaded from CHM13 - T2T.


This command generated the files:
1) chm13v2-GTF-GTF.fq.gz (size 52417961)
2) chm13v2-GTF.transcripts.ogx (size 8788712)
3) chm13v2-GTF.exons.ogx (size 13814360)


#### Step 2: Generate exon/transcript-based indexes

Syntax:

          ogmapper index [-k <keysize>] [-g <guider>] [-e <encoding>]
              [-m 0/1] [-gtf <file>] [-ogx <ogx-prefix-file-name>]
              [-o <index-file-no-ext>] <genome-fasta.gz>
      
The -gtf argument specifies the .gtf/.gff file used in Step 1. The -ogx argument specifies the name used in -o argument in Step 1.

Output:
1) .ogi file (index).

Example:

        ogmapper index -k 24 -g TupleGuider:XYX-Tuple-og.txt -e BitwiseAT1GC0Encoding -m 1 
                -gtf GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz -ogx chm13v2-GTF 
                -o chm13v2-GTF chm13v2-GTF-GTF.fq.gz

Output:
1) chm13v2-GTF.ogi (size 281751943)
2) chm13v2-GTF-GTF_gTplXYXt_eBWat1gc0-24.TO_REMOVE.ogX (size 33759692)


#### Step 3: Counting RNA reads

Syntax:

           ogmapper count [-k <keysize>] [-s <sched>] [-q <nQueue>] [-UC <0/1>]
               [-t <nThreads>] [-p <0/1>] -i <complete-index-file> [-d <n>]
               [-maxreadlen <length>] [-maxreads <n>] [-startread <n>]
               [-genes <outfile>] [-exons <outfile>] [-transcripts <outfile>]
               [-o <out.sam>|stdout] [-kseq <0/1>] [-unmapped 0/1/2/3/4] [-R <str>]
               [-1 <ogx> <reads.gz>] [-2 <ogx> <read-1.gz> <read-2.gz>]

