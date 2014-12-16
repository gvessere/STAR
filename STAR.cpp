#include "IncludeDefine.h"
#include "Parameters.h"
#include "SequenceFuns.h"
#include "Genome.h"
#include "ReadAlignChunk.h"
#include "ReadAlign.h"
#include "Stats.h"
#include "genomeGenerate.h"
#include "outputSJ.h"
#include "ThreadControl.h"
#include "GlobalVariables.cpp"
#include "TimeFunctions.h"
#include <sys/types.h>
#include <sys/stat.h>
#include "ErrorWarning.h"
#include "sysRemoveDir.h"

int main(int argInN, char* argIn[]) {
   
    time(&g_statsAll.timeStart);
   
    Parameters *P = new Parameters; //all parameters
       
    P->inputParameters(argInN, argIn);
    
    *(P->inOut->logStdOut) << timeMonthDayTime(g_statsAll.timeStart) << " ..... Started STAR run\n" <<flush;        

    //generate genome
    if (P->runMode=="genomeGenerate") {
        genomeGenerate(P);
        (void) sysRemoveDir (P->outFileTmp);        
        P->inOut->logMain << "DONE: Genome generation, EXITING\n" << flush;
        exit(0);
    } else if (P->runMode!="alignReads") {
        P->inOut->logMain << "EXITING because of INPUT ERROR: unknown value of input parameter runMode=" <<P->runMode<<endl<<flush;
        exit(1);
    };
    
    Genome mainGenome(P);
    mainGenome.genomeLoad();
    //calculate genome-related parameters
    P->winBinN = P->nGenome/(1LLU << P->winBinNbits)+1;

/////////////////////////////////////////////////////////////////////////////////////////////////START
    
    if (P->outSAMmode != "None") {//open SAM file and write header
        *P->inOut->outSAM << "@PG\tID:STAR\tPN:STAR\tVN:" << SVN_VERSION_COMPILED <<"\tCL:" << P->commandLineFull <<  "\tcl:" << P->commandLine <<endl;
        for (uint ii=0;ii<P->nChrReal;ii++) {
            *P->inOut->outSAM << "@SQ\tSN:"<< P->chrName.at(ii) <<"\tLN:"<<P->chrLength[ii]<<endl;
        };
    };
    
    if (P->chimSegmentMin>0) {
        P->inOut->outChimJunction.open((P->outFileNamePrefix + "Chimeric.out.junction").c_str());
        P->inOut->outChimSAM.open((P->outFileNamePrefix + "Chimeric.out.sam").c_str());
        P->inOut->outChimSAM << "@PG\tID:STAR\tPN:STAR\tVN:" << SVN_VERSION_COMPILED <<"\tCL:" << P->commandLineFull <<  "\tcl:" << P->commandLine <<endl;
        
        for (uint ii=0;ii<P->nChrReal;ii++) {
            P->inOut->outChimSAM << "@SQ\tSN:"<< P->chrName.at(ii) <<"\tLN:"<<P->chrLength[ii]<<endl;
        };        
        pthread_mutex_init(&g_threadChunks.mutexOutChimSAM, NULL);   
        pthread_mutex_init(&g_threadChunks.mutexOutChimJunction, NULL);
    };
         
    // P->inOut->logMain << "mlock value="<<mlockall(MCL_CURRENT|MCL_FUTURE) <<"\n"<<flush;
    
   
    ReadAlignChunk *RAchunk[P->runThreadN];
    for (int ii=0;ii<P->runThreadN;ii++) {
        RAchunk[ii]=new ReadAlignChunk(P, mainGenome, ii);
        RAchunk[ii]->RA->iRead=0;
        RAchunk[ii]->iThread=ii;
    };
    
    if (P->runThreadN>1) {
        g_threadChunks.threadArray=new pthread_t[P->runThreadN];
        pthread_mutex_init(&g_threadChunks.mutexInRead, NULL);
        pthread_mutex_init(&g_threadChunks.mutexOutSAM, NULL);
        pthread_mutex_init(&g_threadChunks.mutexOutUnmappedFastx, NULL);
        pthread_mutex_init(&g_threadChunks.mutexOutFilterBySJout, NULL);
        pthread_mutex_init(&g_threadChunks.mutexStats, NULL);
    };
    
    
    ///////////////////////////////////////////////////////////////////
    g_statsAll.progressReportHeader(P->inOut->logProgress);    
    time(&g_statsAll.timeStartMap);
    *P->inOut->logStdOut << timeMonthDayTime(g_statsAll.timeStartMap) << " ..... Started mapping\n" <<flush;
    
    g_statsAll.timeLastReport=g_statsAll.timeStartMap;
    
    for (int ithread=1;ithread<P->runThreadN;ithread++) {//spawn threads
        pthread_create(&g_threadChunks.threadArray[ithread], NULL, &g_threadChunks.threadRAprocessChunks, (void *) RAchunk[ithread]);
    };
    
    RAchunk[0]->processChunks(); //start main thread
    
    for (int ithread=1;ithread<P->runThreadN;ithread++) {//wait for all threads to complete
        int threadJoinStatus = pthread_join(g_threadChunks.threadArray[ithread], NULL);
        if (threadJoinStatus) {//something went wrong with one of threads
                ostringstream errOut;
                errOut << "EXITING because of FATAL ERROR: phtread error while joining thread # " << ithread <<", error code: "<<threadJoinStatus ;
                exitWithError(errOut.str(),std::cerr, P->inOut->logMain, 1, *P);
        };
    };    
    
    if (P->outFilterBySJoutStage==1) {//completed stage 1, go to stage 2
        outputSJ(RAchunk,P);//collapse novel junctions
        
        P->outFilterBySJoutStage=2;
        
        for (int ithread=1;ithread<P->runThreadN;ithread++) {//spawn threads
            pthread_create(&g_threadChunks.threadArray[ithread], NULL, &g_threadChunks.threadRAprocessChunks, (void *) RAchunk[ithread]);
        };

        RAchunk[0]->processChunks(); //start main thread

        for (int ithread=1;ithread<P->runThreadN;ithread++) {//wait for all threads to complete
            int threadJoinStatus = pthread_join(g_threadChunks.threadArray[ithread], NULL);
            if (threadJoinStatus) {//something went wrong with one of threads
                ostringstream errOut;
                errOut << "EXITING because of FATAL ERROR: phtread error while joining thread # " << ithread <<", error code: "<<threadJoinStatus ;
                exitWithError(errOut.str(),std::cerr, P->inOut->logMain, 1, *P);
            };
        };  
    };
    
    if (P->runThreadN>1 && P->outSAMorder=="PairedKeepInputOrder") {//concatenate Aligned.* files
        RAchunk[0]->chunkFilesCat(P->inOut->outSAM, P->outFileTmp + "/Aligned.out.sam.chunk", g_threadChunks.chunkOutN);
    };    
    
    //aggregate output (junctions, signal, etc)
    //collapse splice junctions from different threads/chunks, and output them
    outputSJ(RAchunk,P);
    
    g_statsAll.progressReport(P->inOut->logProgress);
    P->inOut->logProgress  << "ALL DONE!\n"<<flush;
    P->inOut->logFinal.open((P->outFileNamePrefix + "Log.final.out").c_str());
    g_statsAll.reportFinal(P->inOut->logFinal,P);
    *P->inOut->logStdOut << timeMonthDayTime(g_statsAll.timeFinish) << " ..... Finished successfully\n" <<flush;
    
    P->inOut->logMain  << "ALL DONE!\n"<<flush;
   sysRemoveDir (P->outFileTmp);
    
    delete P->inOut; //to close files
    
    return 0;    
};
