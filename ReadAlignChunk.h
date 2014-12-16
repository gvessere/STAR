#ifndef READ_ALIGN_CHUNK_DEF
#define READ_ALIGN_CHUNK_DEF

#include "IncludeDefine.h"
#include "Parameters.h"
#include "ReadAlign.h"
#include "OutSJ.h"

class ReadAlignChunk {//chunk of reads and alignments
public:
    Parameters* P;
    ReadAlign* RA;
    
    char** chunkIn; //space for the chunk of input reads
    char*  chunkOutSAM;//space for the chunk of output SAM
    OutSJ *chunkOutSJ, *chunkOutSJ1;
    
    istringstream** readInStream;
    ostringstream*  chunkOutSAMstream;
    ofstream chunkOutSAMfile;
    string chunkOutSAMfileName;
    
    bool noReadsLeft;
    uint iChunkIn; //current chunk # as read from .fastq
    uint iChunkOutSAM; //current chunk # writtedn to Aligned.out.sam
    int iThread; //current thread
    
    ReadAlignChunk(Parameters* Pin, Genome &genomeIn, int iChunk);
    void processChunks();
    void mapChunk();
    void chunkFstreamOpen(string filePrefix, int iChunk, fstream &fstreamOut);
    void chunkFstreamCat (fstream &chunkOut, ofstream &allOut, bool mutexFlag, pthread_mutex_t &mutexVal);
    void chunkFilesCat(ostream *allOut, string filePrefix, uint &iC);
};
#endif
