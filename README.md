# ogmapper
A Fast and Light Genomic Mapper for short reads.

ogMapper generates a tiny index file and memory operations, typically as a fraction of the genome size. 
For example, 1.7 GiB for T2T human genome v2 plus 0.7 GiB for the compressed genome (file ~ 2.4 GiB) and similar amounts of memory in run time, irrespective of the number of threads.
The common process starts by creating an index file (*.ogi), then mapping short reads (in pairs or not) to generate .sam files.

ogmapper is written in c++. I provide binaries for selected operating systems and source files for compilation.

# Binaries
- MacOS/Intel
- MacOS/M
- Linux/x86
- Windows :(

# Installation
- Download

# Compilation
ogMapper uses the WFA2 library to perform alignment operations when needed. I provide the latest version used to compile the above binaries. Users may opt to download the lastest WFA version from https://github.com/smarco/WFA2-lib. You may follow the WFA2 instructions or follow the steps below as a guide. The idea is to build a library file suitable for ogMapper (libwfacpp.a) that is needed in the /lib folder to be able to compile ogMapper.

Once WFA library is ok. We may proceed to compile ogMapper.

## Compiling WFA2 from .zip
Ok
## Compiling WFA2 from https://github.com/smarco/WFA2-lib
Ok
## Compiling ogMapper in MacOS
Ok
## Compiling ogMapper in Linux
Ok

# Running ogMapper

## Index Generation for DNA
Ok
## Mapping short DNA reads
Ok
## Index Generation for RNA
Ok
## Mapping/Counting RNA reads
Ok
## Other Options
Ok
