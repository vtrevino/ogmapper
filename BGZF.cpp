// Downloaded from https://github.com/ekg/hapviz
// ***************************************************************************
// BGZF.cpp (c) 2009 Derek Barnett, Michael Str�mberg
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 16 August 2010 (DB)
// ---------------------------------------------------------------------------
// BGZF routines were adapted from the bgzf.c code developed at the Broad
// Institute.
// ---------------------------------------------------------------------------
// Provides the basic functionality for reading & writing BGZF files
// ***************************************************************************

//#include <algorithm>
#include "BGZF.h"
using namespace BamTools;
using std::string;
using std::min;

BgzfData::BgzfData(void)
    : UncompressedBlockSize(DEFAULT_BLOCK_SIZE)
    , CompressedBlockSize(MAX_BLOCK_SIZE)
    , BlockLength(0)
    , BlockOffset(0)
    , BlockAddress(0)
    , IsOpen(false)
    , IsWriteOnly(false)
    , IsWriteUncompressed(false)
    , Stream(NULL)
    , UncompressedBlock(NULL)
    , CompressedBlock(NULL)
{
    try {
        CompressedBlock   = new char[CompressedBlockSize];
        UncompressedBlock = new char[UncompressedBlockSize];
    } catch( std::bad_alloc& ba ) {
        fprintf(stderr, "BGZF ERROR: unable to allocate memory for our BGZF object.\n");
        exit(1);
    }
}

// destructor
BgzfData::~BgzfData(void) {
    if( CompressedBlock   ) delete[] CompressedBlock;
    if( UncompressedBlock ) delete[] UncompressedBlock;
}

// closes BGZF file
void BgzfData::Close(void) {

    // skip if file not open, otherwise set flag
    if ( !IsOpen ) return;

    // if writing to file, flush the current BGZF block,
    // then write an empty block (as EOF marker)
    if ( IsWriteOnly ) {
        FlushBlock();
        int blockLength = DeflateBlock();
        fwrite(CompressedBlock, 1, blockLength, Stream);
    }
    
    // flush and close
    fflush(Stream);
    fclose(Stream);
    IsWriteUncompressed = false;
    IsOpen = false;
}

// compresses the current block
int BgzfData::DeflateBlock(void) {

    // initialize the gzip header
    char* buffer = CompressedBlock;
    memset(buffer, 0, 18);
    buffer[0]  = GZIP_ID1;
    buffer[1]  = (char)GZIP_ID2;
    buffer[2]  = CM_DEFLATE;
    buffer[3]  = FLG_FEXTRA;
    buffer[9]  = (char)OS_UNKNOWN;
    buffer[10] = BGZF_XLEN;
    buffer[12] = BGZF_ID1;
    buffer[13] = BGZF_ID2;
    buffer[14] = BGZF_LEN;

    // set compression level
    const int compressionLevel = ( IsWriteUncompressed ? 0 : Z_DEFAULT_COMPRESSION );
    
    // loop to retry for blocks that do not compress enough
    int inputLength = BlockOffset;
    int compressedLength = 0;
    unsigned int bufferSize = CompressedBlockSize;

    while ( true ) {
        
        // initialize zstream values
        z_stream zs;
        zs.zalloc    = NULL;
        zs.zfree     = NULL;
        zs.next_in   = (Bytef*)UncompressedBlock;
        zs.avail_in  = inputLength;
        zs.next_out  = (Bytef*)&buffer[BLOCK_HEADER_LENGTH];
        zs.avail_out = bufferSize - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;

        // initialize the zlib compression algorithm
        if ( deflateInit2(&zs, compressionLevel, Z_DEFLATED, GZIP_WINDOW_BITS, Z_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY) != Z_OK ) {
            fprintf(stderr, "BGZF ERROR: zlib deflate initialization failed.\n");
            exit(1);
        }

        // compress the data
        int status = deflate(&zs, Z_FINISH);
        if ( status != Z_STREAM_END ) {

            deflateEnd(&zs);

            // reduce the input length and try again
            if ( status == Z_OK ) {
                inputLength -= 1024;
                if( inputLength < 0 ) {
                    fprintf(stderr, "BGZF ERROR: input reduction failed.\n");
                    exit(1);
                }
                continue;
            }

            fprintf(stderr, "BGZF ERROR: zlib::deflateEnd() failed.\n");
            exit(1);
        }

        // finalize the compression routine
        if ( deflateEnd(&zs) != Z_OK ) {
            fprintf(stderr, "BGZF ERROR: zlib::deflateEnd() failed.\n");
            exit(1);
        }

        compressedLength = zs.total_out;
        compressedLength += BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;
        if ( compressedLength > MAX_BLOCK_SIZE ) {
            fprintf(stderr, "BGZF ERROR: deflate overflow.\n");
            exit(1);
        }

        break;
    }

    // store the compressed length
    BgzfData::PackUnsignedShort(&buffer[16], (unsigned short)(compressedLength - 1));

    // store the CRC32 checksum
    unsigned int crc = crc32(0, NULL, 0);
    crc = crc32(crc, (Bytef*)UncompressedBlock, inputLength);
    BgzfData::PackUnsignedInt(&buffer[compressedLength - 8], crc);
    BgzfData::PackUnsignedInt(&buffer[compressedLength - 4], inputLength);

    // ensure that we have less than a block of data left
    int remaining = BlockOffset - inputLength;
    if ( remaining > 0 ) {
        if ( remaining > inputLength ) {
            fprintf(stderr, "BGZF ERROR: after deflate, remainder too large.\n");
            exit(1);
        }
        memcpy(UncompressedBlock, UncompressedBlock + inputLength, remaining);
    }

    BlockOffset = remaining;
    return compressedLength;
}

// flushes the data in the BGZF block
void BgzfData::FlushBlock(void) {

    // flush all of the remaining blocks
    while ( BlockOffset > 0 ) {

        // compress the data block
        int blockLength = DeflateBlock();

        // flush the data to our output stream
        int numBytesWritten = fwrite(CompressedBlock, 1, blockLength, Stream);

        if ( numBytesWritten != blockLength ) {
          fprintf(stderr, "BGZF ERROR: expected to write %u bytes during flushing, but wrote %u bytes.\n", blockLength, numBytesWritten);
          exit(1);
        }
              
        BlockAddress += blockLength;
    }
}

// de-compresses the current block
int BgzfData::InflateBlock(const int& blockLength) {

    // Inflate the block in m_BGZF.CompressedBlock into m_BGZF.UncompressedBlock
    z_stream zs;
    zs.zalloc    = NULL;
    zs.zfree     = NULL;
    zs.next_in   = (Bytef*)CompressedBlock + 18;
    zs.avail_in  = blockLength - 16;
    zs.next_out  = (Bytef*)UncompressedBlock;
    zs.avail_out = UncompressedBlockSize;

    int status = inflateInit2(&zs, GZIP_WINDOW_BITS);
    if ( status != Z_OK ) {
        fprintf(stderr, "BGZF ERROR: could not decompress block - zlib::inflateInit() failed\n");
        return -1;
    }

    status = inflate(&zs, Z_FINISH);
    if ( status != Z_STREAM_END ) {
        inflateEnd(&zs);
        fprintf(stderr, "BGZF ERROR: could not decompress block - zlib::inflate() failed\n");
        return -1;
    }

    status = inflateEnd(&zs);
    if ( status != Z_OK ) {
        fprintf(stderr, "BGZF ERROR: could not decompress block - zlib::inflateEnd() failed\n");
        return -1;
    }

    return zs.total_out;
}

// opens the BGZF file for reading (mode is either "rb" for reading, or "wb" for writing)
bool BgzfData::Open(const string& filename, const char* mode, bool isWriteUncompressed ) {

    // determine open mode
    if ( strcmp(mode, "rb") == 0 )
        IsWriteOnly = false;
    else if ( strcmp(mode, "wb") == 0) 
        IsWriteOnly = true;
    else {
        fprintf(stderr, "BGZF ERROR: unknown file mode: %s\n", mode);
        return false; 
    }

    // ----------------------------------------------------------------
    // open Stream to read to/write from file, stdin, or stdout
    // stdin/stdout option contributed by Aaron Quinlan (2010-Jan-03)
    
    // read/write BGZF data to/from a file
    if ( (filename != "stdin") && (filename != "stdout") )
        Stream = fopen(filename.c_str(), mode);
    
    // read BGZF data from stdin
    else if ( (filename == "stdin") && (strcmp(mode, "rb") == 0 ) )
        Stream = freopen(NULL, mode, stdin);
    
    // write BGZF data to stdout
    else if ( (filename == "stdout") && (strcmp(mode, "wb") == 0) )
        Stream = freopen(NULL, mode, stdout);

    if ( !Stream ) {
        fprintf(stderr, "BGZF ERROR: unable to open file %s\n", filename.c_str() );
        return false;
    }
    
    // set flags, return success
    IsOpen = true;
    IsWriteUncompressed = isWriteUncompressed;
    return true;
}

// reads BGZF data into a byte buffer
int BgzfData::Read(char* data, const unsigned int dataLength) {

   if ( !IsOpen || IsWriteOnly || dataLength == 0 ) return 0;

   char* output = data;
   unsigned int numBytesRead = 0;
   while ( numBytesRead < dataLength ) {

       int bytesAvailable = BlockLength - BlockOffset;
       if ( bytesAvailable <= 0 ) {
           if ( !ReadBlock() ) return -1; 
           bytesAvailable = BlockLength - BlockOffset;
           if ( bytesAvailable <= 0 ) break;
       }

       char* buffer   = UncompressedBlock;
       int copyLength = min( (int)(dataLength-numBytesRead), bytesAvailable );
       memcpy(output, buffer + BlockOffset, copyLength);

       BlockOffset  += copyLength;
       output       += copyLength;
       numBytesRead += copyLength;
   }

   if ( BlockOffset == BlockLength ) {
       BlockAddress = ftell64(Stream);
       BlockOffset  = 0;
       BlockLength  = 0;
   }

   return numBytesRead;
}

// reads a BGZF block
bool BgzfData::ReadBlock(void) {

    char    header[BLOCK_HEADER_LENGTH];
    int64_t blockAddress = ftell64(Stream);
    
    int count = fread(header, 1, sizeof(header), Stream);
    if ( count == 0 ) {
        BlockLength = 0;
        return true;
    }

    if ( count != sizeof(header) ) {
        fprintf(stderr, "BGZF ERROR: read block failed - could not read block header\n");
        return false;
    }

    if ( !BgzfData::CheckBlockHeader(header) ) {
        fprintf(stderr, "BGZF ERROR: read block failed - invalid block header\n");
        return false;
    }

    int blockLength = BgzfData::UnpackUnsignedShort(&header[16]) + 1;
    char* compressedBlock = CompressedBlock;
    memcpy(compressedBlock, header, BLOCK_HEADER_LENGTH);
    int remaining = blockLength - BLOCK_HEADER_LENGTH;

    count = fread(&compressedBlock[BLOCK_HEADER_LENGTH], 1, remaining, Stream);
    if ( count != remaining ) {
        fprintf(stderr, "BGZF ERROR: read block failed - could not read data from block\n");
        return false;
    }

    count = InflateBlock(blockLength);
    if ( count < 0 ) { 
      fprintf(stderr, "BGZF ERROR: read block failed - could not decompress block data\n");
      return false;
    }

    if ( BlockLength != 0 )
        BlockOffset = 0;

    BlockAddress = blockAddress;
    BlockLength  = count;
    return true;
}

// seek to position in BGZF file
bool BgzfData::Seek(int64_t position) {

    if ( !IsOpen ) return false;
  
    int     blockOffset  = (position & 0xFFFF);
    int64_t blockAddress = (position >> 16) & 0xFFFFFFFFFFFFLL;

    if ( fseek64(Stream, blockAddress, SEEK_SET) != 0 ) {
        fprintf(stderr, "BGZF ERROR: unable to seek in file\n");
        return false;
    }

    BlockLength  = 0;
    BlockAddress = blockAddress;
    BlockOffset  = blockOffset;
    return true;
}

// get file position in BGZF file
int64_t BgzfData::Tell(void) {
    if ( !IsOpen ) 
        return false;
    else 
        return ( (BlockAddress << 16) | (BlockOffset & 0xFFFF) );
}

// writes the supplied data into the BGZF buffer
unsigned int BgzfData::Write(const char* data, const unsigned int dataLen) {

    if ( !IsOpen || !IsWriteOnly ) return false;
  
    // initialize
    unsigned int numBytesWritten = 0;
    const char* input = data;
    unsigned int blockLength = UncompressedBlockSize;

    // copy the data to the buffer
    while ( numBytesWritten < dataLen ) {
      
        unsigned int copyLength = min(blockLength - BlockOffset, dataLen - numBytesWritten);
        char* buffer = UncompressedBlock;
        memcpy(buffer + BlockOffset, input, copyLength);

        BlockOffset     += copyLength;
        input           += copyLength;
        numBytesWritten += copyLength;

        if ( BlockOffset == blockLength )
            FlushBlock();
    }

    return numBytesWritten;
}
