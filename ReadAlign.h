#ifndef READALIGN_DEF
#define READALIGN_DEF

#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "Genome.h"
#include "Stats.h"
#include "OutSJ.h"
#include <time.h>

class ReadAlign : public Genome {
    public:
        Parameters* P; //pointer to the parameters, will be initialized on construction
          
        //mapping statistics
        Stats statsRA;
        
        //mapping time
        time_t timeStart, timeFinish;
        
        
        //input,output
        istream* readInStream[MAX_N_MATES];
        ostream* outSAMstream;
        OutSJ *chunkOutSJ, *chunkOutSJ1;
        fstream chunkOutChimSAM, chunkOutChimJunction, chunkOutUnmappedReadsStream[MAX_N_MATES], chunkOutFilterBySJoutFiles[MAX_N_MATES];
        
        ostringstream samStreamCIGAR, samStreamSJmotif, samStreamSJintron,samStreamSJannot;
        
        intScore maxScoreMate[MAX_N_MATES];
        intScore *scoreSeedToSeed, *scoreSeedBest;
        uint *scoreSeedBestInd, *seedChain, *scoreSeedBestMM;
        
//         StatsAll *statsRA;
        
        //transcript
        Transcript* trArray; //linear array of transcripts to store all of them from all windows
        Transcript** trArrayPointer; //linear array of transcripts to store all of them from all windows            
        
        //read
        uint iRead, iMate;
        bool revertStrand; //what to do with the strand, according to strandType and iMate
        uint Lread, readLength[MAX_N_MATES], readLengthOriginal[MAX_N_MATES], readLengthPair, readLengthPairOriginal;
        uint clip3pNtotal[MAX_N_MATES], clip5pNtotal[MAX_N_MATES], clip3pAdapterN[MAX_N_MATES]; //total number of trimmed bases from 5p,3p
        int readFileType; //file type: 1=fasta; 2=fastq
        
        char dummyChar[4096];
        char** Read0;
        char** Qual0;
        char** readNameMates;
        char* readName;
        char** Read1;
        char** Qual1; //modified QSs for scoring
        
        //split            
        uint** splitR;
        uint Nsplit;
        
//         uint fragLength[MAX_N_FRAG], fragStart[MAX_N_FRAG]; //fragment Lengths and Starts in read space
        
        //binned alignments
        uintWinBin **winBin; //binned genome: window ID (number) per bin
        
        //alignments
        uiPC *PC; //pieces coordinates
        uiWC *WC; //windows coordinates        
        uiWA **WA; //aligments per window
        
        int unmapType; //marker for why a read is unmapped
        
        uint mapMarker; //alignment marker (typically, if there is something wrong)
        uint nA, nP, nW, nWall, nUM[2]; //number of all alignments,  pieces, windows, U/M, 
        uint *nWA, *nWAP, *WALrec, *WlastAnchor; //number of alignments per window, per window per piece, min recordable length per window
        bool *WAincl; //alginment inclusion mask
        
        uint *swWinCov, *swWinGleft, *swWinGright, *swWinRleft, *swWinRright; //read coverage per window
        char *swT;
        
        uint storedLmin, uniqLmax, uniqLmaxInd, multLmax, multLmaxN, multNmin, multNminL, multNmax, multNmaxL;
        uint nTr, nTrMate; // number of transcripts called
        intScore maxScore, nextWinScore;//maximum alignment score, next best score
        
        uint chimN, chimRepeat, chimStr, chimMotif;
        Transcript trChim[MAX_N_CHIMERAS];
        
        Transcript trA, trA1, *trBest, *trNext, *trInit; //transcript, best tr, next best tr, initialized tr
        Transcript ***trAll; //all transcripts for all windows
        uint *nWinTr; //number of recorded transcripts per window
        
        Transcript *alignC, *extendC, *polyAtailC; //alignment rules/conditions
        
        intScore trMultScores[MAX_N_MULTMAP];//scores for the multiple mappers
        Transcript* trMult[MAX_N_MULTMAP];//multimapping transcripts
               
        ReadAlign (Parameters* Pin, Genome &genomeIn);//allocate arrays
        void resetN();//resets the counters to 0
        void multMapSelect();
        int mapOneRead();
        uint maxMappableLength2strands(uint pieceStart, uint pieceLength, uint iDir, uint iSA1, uint iSA2, uint& maxL, uint iFrag);
        void storeAligns (uint iDir, uint Shift, uint Nrep, uint L, uint indStartEnd[2], uint iFrag);
        bool outputTranscript(Transcript *trOut, uint nTrOut, ofstream *outBED);
        void outputTranscriptSAM(Transcript const &trOut, uint nTrOut, uint iTrOut, uint mateChr, uint mateStart, char mateStrand, ostream *outStream);
        void outputTranscriptSJ(Transcript const &trOut, uint nTrOut, OutSJ *outStream, uint sjReadStartN );
        string outputTranscriptCIGARp(Transcript const &trOut);
        void outTxtMain(ofstream*,Transcript&);
        int createExtendWindowsWithAlign(uint a1, uint aStr); //extends and windows with one alignment
        void assignAlignToWindow(uint a1, uint aLength, uint aStr, uint aNrep, uint aFrag, uint aRstart,bool aAnchor, uint sjA); //assigns one alignment to a window
        void stitchPieces(char **R, char **Q, char *G, PackedArray& SA, uint Lread);
        void outputAlignments();
        void stitchWindowSeeds (uint iW, uint iWrec, char* R, char* Q, char* G);//stitches all seeds in one window: iW
        
        
        int oneRead();
        
};

#endif


