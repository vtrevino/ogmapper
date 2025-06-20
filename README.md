# ogmapper
A Fast and Light Genomic Mapper for short reads.

ogMapper generates a tiny index file and memory operations, typically as a fraction of the genome size. 
For example, 1.7 GiB for T2T human genome v2 plus 0.7 GiB for the compressed genome (file ~ 2.4 GiB) and similar amounts of memory in run time, irrespective of the number of threads.
The common process starts by creating an index file (*.ogi), then mapping short reads (in pairs or not) to generate .sam files.

# Installation
- download

# Compilation
ogMapper uses the WFA2 library to perform alignment operations when needed.
