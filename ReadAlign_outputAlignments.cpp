#include "ReadAlign.h"
#include "GlobalVariables.h"
#include "ErrorWarning.h"

void ReadAlign::outputAlignments() {
    

    bool mateMapped[2]={false,false};
    
    if ( nW==0 ) {//no good windows
        statsRA.unmappedOther++;
        unmapType=0;
    } else if ( (trBest->maxScore < P->outFilterScoreMin) | (trBest->maxScore < (intScore) (P->outFilterScoreMinOverLread*(Lread-1))) \
              | (trBest->nMatch < P->outFilterMatchNmin)  | (trBest->nMatch < (uint) (P->outFilterMatchNminOverLread*(Lread-1))) ) {//too short
        statsRA.unmappedShort++;
        unmapType=1;
    } else if ( (trBest->nMM > P->outFilterMismatchNmax) | (double(trBest->nMM)/double(trBest->rLength)>P->outFilterMismatchNoverLmax) ) {//too many mismatches
        statsRA.unmappedMismatch++;
        unmapType=2;
    } else if (nTr > P->outFilterMultimapNmax){//too multi
        statsRA.unmappedMulti++;
        unmapType=3;
    } else {//output transcripts

        bool outFilterPassed(true);
        if (P->outFilterBySJoutStage==1) {//no filtering by SJout
            for (uint iTr=0;iTr<nTr;iTr++) {//check transcript for unannotated junctions
                for (uint iex=0;iex<trMult[iTr]->nExons-1;iex++) {//check all junctions
                    if (trMult[iTr]->canonSJ[iex]>=0 && trMult[iTr]->sjAnnot[iex]==0) {
                        outFilterPassed=false;
                        break;
                    };
                };
                if (!outFilterPassed) break;
            };
            if (!outFilterPassed) {//this read is held for further filtering BySJout, record fastq
                unmapType=-3; //the read is not conisddred unmapped
                statsRA.readN--;
                statsRA.readBases -= readLength[0]+readLength[1];
                
//                 if (P->runThreadN>1) pthread_mutex_lock(&g_threadChunks.mutexOutFilterBySJout);
                for (uint im=0;im<P->readNmates;im++) {
                   chunkOutFilterBySJoutFiles[im] << readNameMates[im] << "/" <<im+1;
                   chunkOutFilterBySJoutFiles[im] <<"\n";
                   chunkOutFilterBySJoutFiles[im] << Read0[im] <<"\n";
                    if (readFileType==2) {//fastq
                        chunkOutFilterBySJoutFiles[im] << "+\n";
                        chunkOutFilterBySJoutFiles[im] << Qual0[im] <<"\n";
                    };
                };
//                 if (P->runThreadN>1) pthread_mutex_unlock(&g_threadChunks.mutexOutFilterBySJout);  
            };
        };

        if (P->outSJfilterReads=="All" || nTr==1) {
            uint sjReadStartN=chunkOutSJ1->N;        
            for (uint iTr=0;iTr<nTr;iTr++) {//write all transcripts
                outputTranscriptSJ (*(trMult[iTr]), nTr, chunkOutSJ1, sjReadStartN);            
            };
        };        

        if (outFilterPassed) {
            if (nTr>1) {//multimappers
                statsRA.mappedReadsM++;
                unmapType=-1;
            } else if (nTr==1) {//unique mappers
                statsRA.mappedReadsU++;
                statsRA.transcriptStats(*(trMult[0]),Lread);
                unmapType=-2;
            } else {//cannot be
                ostringstream errOut;
                errOut  << "EXITING because of a BUG: nTr=0 in outputAlignments.cpp";
                exitWithError(errOut.str(), std::cerr, P->inOut->logMain, EXIT_CODE_BUG, *P);                    
            };            
            
            for (uint iTr=0;iTr<nTr;iTr++) {//write all transcripts
                outputTranscriptSAM(*(trMult[iTr]), nTr, iTr, (uint) -1, (uint) -1, 0, outSAMstream);
            };
            if (P->outSJfilterReads=="All" || nTr==1) {
                uint sjReadStartN=chunkOutSJ->N;        
                for (uint iTr=0;iTr<nTr;iTr++) {//write all transcripts
                    outputTranscriptSJ (*(trMult[iTr]), nTr, chunkOutSJ, sjReadStartN);            
                };
            };
            mateMapped[trBest->exons[0][EX_iFrag]]=true;
            mateMapped[trBest->exons[trBest->nExons-1][EX_iFrag]]=true;        
            if (P->readNmates>1 && !(mateMapped[0] && mateMapped[1]) ) unmapType=4;            
        };



    };

    if (unmapType>=0 && P->outSAMunmapped=="Within") {//unmapped read, at least one mate
        for (uint im=0; im<P->readNmates; im++) {   
            if (!mateMapped[im]) {
                int samFLAG=0x4;
                if (P->readNmates==2) {//paired read
                    samFLAG+=0x1 + (im==0 ? 0x40 : 0x80);
                    if (mateMapped[1-im]) {//mate mapped
                        if (trBest->Str!=1-im) samFLAG+=0x20;//mate strand reverted
                    } else {//mate unmapped
                        samFLAG+=0x8;
                    };
                };
                
                *outSAMstream << readName+1 <<"\t"<< samFLAG \
                        <<"\t"<< '*' <<"\t"<< '0' <<"\t"<< '0' <<"\t"<< '*';
                if (mateMapped[1-im]) {//mate is mapped
                    *outSAMstream <<"\t"<< P->chrName[trBest->Chr] <<"\t"<< trBest->exons[0][EX_G] + 1 - P->chrStart[trBest->Chr];
                } else {
                    *outSAMstream <<"\t"<< '*' <<"\t"<< '0';
                };
                *outSAMstream <<"\t"<< '0' <<"\t"<< Read0[im] <<"\t"<< (readFileType==2 ? Qual0[im]:"*") \
                        <<"\tNH:i:0" <<"\tHI:i:0" <<"\tAS:i:"<<trBest->maxScore <<"\tnM:i:"<<trBest->nMM<<"\tuT:A:" <<unmapType <<"\n";
            };
        };
    };
    if (unmapType>=0 && P->outReadsUnmapped=="Fastx" ){//output to fasta/q files
           for (uint im=0;im<P->readNmates;im++) {
               chunkOutUnmappedReadsStream[im] << readNameMates[im] << "/" <<im+1;
               if (P->readNmates>1) chunkOutUnmappedReadsStream[im] <<"\t"<< int(mateMapped[0]) <<  int(mateMapped[1]);
               chunkOutUnmappedReadsStream[im] <<"\n";
               chunkOutUnmappedReadsStream[im] << Read0[im] <<"\n";
                if (readFileType==2) {//fastq
                    chunkOutUnmappedReadsStream[im] << "+\n";
                    chunkOutUnmappedReadsStream[im] << Qual0[im] <<"\n";
                };
           };
    }; 
};



