# ogmapper
A Fast and Light Genomic Mapper for short reads.

ogMapper generates a tiny index file and memory operations, typically as a fraction of the genome size. 
For example, 1.7 GiB for T2T human genome v2 plus 0.7 GiB for the compressed genome (file ~ 2.4 GiB) and similar amounts of memory in run time, irrespective of the number of threads. This is for 24-bit keys.
The common process starts by creating an index file (*.ogi), then mapping short reads (in pairs or not) to generate .sam files.

ogmapper is written in c++. I provide binaries for selected operating systems and source files for compilation.

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
  
        # Edit all Makefiles replacing -lomp flag to -fopenmp (Makefile,  wavefront/Makefile, tools/align_benchmark/Makefile, examples/Makefile)
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


### Compiling ogMapper on Linux
Ok

# Introduction to Encodings, Guiders, and Keys
ogMapper introduced some concepts and also used novel explorative functions. ogMapper does not index all possible DNA subsequences. It first sweeps the DNA sequence until a specific sequence pattern is found; this pattern is known as a "guide". Thus, there are two guiders implemented up to now. Once a guide is found, the following DNA sequence of length k (provided by -k argument) is used to build a key, which will be indexed. The key is a binary representation of the DNA sequence, which is provided by an encoding function. There are four encoding functions implemented so far. Once a guide-key is processed, the DNA sweeping continues after the last guide until the whole sequence has been analyzed.

### Guiders (-g option)
- TupleGuider : The tuple guider is a simple map of all possible n-nt combinations specifying which of them will be indexed and which will not. The file format is as follows:

To use the TuppleGuider use the option -g TupleGuider:&lt;tuple-file&gt;.
  
- StateMachineGuider:&lt;state-file&gt;

### Encodings
The encoding transforms the DNA sequence to a binary key. I recommend BitWise and Plain.

The **BitwiseAT1GC0Encoding** uses 1 bit per nt, transforming A or T/U to 1 and G or C to 0. Any other letter is treated as A. Case-insensitive. So, -k 24 will use 24 nt to generate keys of 24 bits for a total of 16,777,216 different keys. Increasing k would have an important impact on the index size and memory needed.

The **PlainEncoding** uses 2 bits per nt, transforming A to 00, C to 01, G to 10, and T/U to 11. Any other letter is treated as A. So, -k 12 will use 12 nt generating keys of 24 bits. 

The **GappedBitwiseAT1GC0Encoding** is similar to BitwiseAT1GC0Encoding but the nt used for indexing are chosen as [left][gap][right]=k where [left] and [right] are estimated by k/3. So, here -k=36 is equivalent in bits to -k=24 in BitwiseAT1GC0Encoding.

The **SwapBitwiseAT1GC0Encoding** is similar to BitwiseAT1GC0Encoding but it uses non-continuous nt indexing, one nt and skipping one until k is reached.

The **HPCEncoding** (homo-polymer compressed) ignores consecutive repetitions of a nucleotide and then uses 'plain' encoding as the **PlainEncoding** above.


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
      Default 0.
                      Default to 0.
              -s      'schedule' functions to call for read mappings.
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


## Index Generation for DNA
For indexing DNA the valid options are:

      ogmapper index [-k <keysize>] [-g <guider>] [-e <encoding>]
        [-m 0/1] [-o <index-file-no-ext>] <genome-fasta.gz>


The typical run

## Mapping short DNA reads
Ok
## Index Generation for RNA
Ok
## Mapping/Counting RNA reads
Ok
## Other Options
Ok
