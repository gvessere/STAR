#include "ReadAlignChunk.h"
#include "GlobalVariables.h"
#include "ThreadControl.h"
#include "ErrorWarning.h"

void ReadAlignChunk::mapChunk() {//map one chunk. Input reads stream has to be setup in RA->readInStream[ii]
    RA->statsRA.resetN();       
    
    for (uint ii=0;ii<P->readNmates;ii++) {//clear eof and rewind the input streams
        RA->readInStream[ii]->clear();
        RA->readInStream[ii]->seekg(0,ios::beg);
    };
    
    if ( P->outSAMorder == "PairedKeepInputOrder" && P->runThreadN>1 ) {//open chunk file
        ostringstream name1("");
        name1 << P->outFileTmp + "/Aligned.tmp.sam.chunk"<<iChunkIn;
        chunkOutSAMfileName = name1.str();
        chunkOutSAMfile.open(chunkOutSAMfileName.c_str());
    };    
    
    int readStatus=0;
    uint chunkOutSAMtotal;
    while (readStatus==0) {//main cycle over all reads

        readStatus=RA->oneRead(); //map one read

        if (readStatus==0) RA->iRead++;

        //write SAM aligns to chunk buffer 
        chunkOutSAMtotal=(uint) RA->outSAMstream->tellp();
        if ( chunkOutSAMtotal > P->chunkOutSAMsizeBytes ) {//this should not happen!
            ostringstream errOut;
            errOut <<"EXITING because of fatal error: buffer size for SAM output is too small\n";
            errOut <<"Solution: increase input parameter --limitOutSAMoneReadBytes\n";
            exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);                      
        } else if ( chunkOutSAMtotal + P->limitOutSAMoneReadBytes > P->chunkOutSAMsizeBytes || (readStatus==-1 && noReadsLeft) ) {//write buffer to disk because it's almost full, or all reads are mapped
            if ( P->outSAMorder == "PairedKeepInputOrder" && P->runThreadN>1 ) {//output chunks into separate files
                chunkOutSAMfile.write(chunkOutSAM,chunkOutSAMtotal);
                chunkOutSAMfile.clear(); //in case 0 bytes were written which could set fail bit
            } else {//standard way, directly into Aligned.out.sam file
                if (P->runThreadN>1) pthread_mutex_lock(&g_threadChunks.mutexOutSAM);    
                P->inOut->outSAM->write(chunkOutSAM,chunkOutSAMtotal);
                P->inOut->outSAM->clear();//in case 0 bytes were written which could set fail bit
                if (P->runThreadN>1) pthread_mutex_unlock(&g_threadChunks.mutexOutSAM);
            };
            RA->outSAMstream->seekp(0,ios::beg); //rewind the chunk storage
        }; 

        //collapse SJ buffer if needed
        if ( chunkOutSJ->N > P->limitOutSJcollapsed ) {//this means the number of collapsed junctions is larger than the chunks size
            ostringstream errOut;
            errOut <<"EXITING because of fatal error: buffer size for SJ output is too small\n";
            errOut <<"Solution: increase input parameter --limitOutSJoneRead\n";
            exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);                      
        } else if ( chunkOutSJ->N + P->limitOutSJoneRead > P->limitOutSJcollapsed || (readStatus==-1 && noReadsLeft) ) {//write buffer to disk because it's almost full, or all reads are mapped
            chunkOutSJ->collapseSJ();
            if ( chunkOutSJ->N + 2*P->limitOutSJoneRead > P->limitOutSJcollapsed ) {
                ostringstream errOut;
                errOut <<"EXITING because of fatal error: buffer size for SJ output is too small\n";
                errOut <<"Solution: increase input parameter --limitOutSJcollapsed\n";
                exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);                      
            };
        };            

        //collapse SJ1 buffer if needed
        if ( chunkOutSJ1->N > P->limitOutSJcollapsed ) {//this means the number of collapsed junctions is larger than the chunks size
            ostringstream errOut;
            errOut <<"EXITING because of fatal error: buffer size for SJ output is too small\n";
            errOut <<"Solution: increase input parameter --limitOutSJoneRead\n";
            exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);                      
        } else if ( chunkOutSJ1->N + P->limitOutSJoneRead > P->limitOutSJcollapsed || (readStatus==-1 && noReadsLeft) ) {//write buffer to disk because it's almost full, or all reads are mapped
            chunkOutSJ1->collapseSJ();
            if ( chunkOutSJ1->N + 2*P->limitOutSJoneRead > P->limitOutSJcollapsed ) {
                ostringstream errOut;
                errOut <<"EXITING because of fatal error: buffer size for SJ output is too small\n";
                errOut <<"Solution: increase input parameter --limitOutSJcollapsed\n";
                exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);                      
            };
        };            

    }; //reads cycle

    if ( P->outSAMorder == "PairedKeepInputOrder" && P->runThreadN>1 ) {//write the remaining part of the buffer, close and rename chunk files
        chunkOutSAMfile.write(chunkOutSAM,chunkOutSAMtotal);
        chunkOutSAMfile.clear(); //in case 0 bytes were written which could set fail bit
        chunkOutSAMfile.close();
        RA->outSAMstream->seekp(0,ios::beg); //rewind the chunk storage
        ostringstream name2("");
        name2 << P->outFileTmp + "/Aligned.out.sam.chunk"<<iChunkIn;                
        rename(chunkOutSAMfileName.c_str(),name2.str().c_str());//marks files as completedly written
    };    
    
    //add stats, write progress if needed
    if (P->runThreadN>1) pthread_mutex_lock(&g_threadChunks.mutexStats);
    g_statsAll.addStats(RA->statsRA);
    g_statsAll.progressReport(P->inOut->logProgress);
    if (P->runThreadN>1) pthread_mutex_unlock(&g_threadChunks.mutexStats); 
};            
