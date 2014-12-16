#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "ReadAlign.h"
#include "SequenceFuns.h"
#include "stitchWindowAligns.h"
#include "sjSplitAlign.cpp"
#include "PackedArray.h"
#include "alignSmithWaterman.h"
#include "GlobalVariables.h"
#include <time.h>

void ReadAlign::stitchPieces(char **R, char **Q, char *G, PackedArray& SA, uint Lread) {
    
       //zero-out winBin
    memset(winBin[0],255,sizeof(winBin[0][0])*P->winBinN);
    memset(winBin[1],255,sizeof(winBin[0][0])*P->winBinN);

    
//     for (uint iWin=0;iWin<nWall;iWin++) {//zero out winBin
//         if (WC[iWin][WC_gStart]<=WC[iWin][WC_gEnd]) {//otherwise the window is dead
//             memset(&(winBin[WC[iWin][WC_Str]][WC[iWin][WC_gStart]]),255,sizeof(winBin[0][0])*(WC[iWin][WC_gEnd]-WC[iWin][WC_gStart]+1));
//         };
// //         for (uint ii=C[iWin][WC_gStart]; ii<WC[iWin][WC_gEnd]; ii++) {
// //             winBin[WC[WC_Str]
// //         };
//     };
    
//     //debug
//     for (uint ii=0;ii<P->winBinN;ii++){
//         if (winBin[0][ii]!=uintWinBinMax || winBin[1][ii]!=uintWinBinMax) {
//             cerr<< "BUG in stitchPieces: ii="<<ii<<"   "<< winBin[0][ii] <<"   "<<winBin[1][ii] <<"   iRead="<<iRead<<"   nW="<<nW<<endl;
//             for (uint iWin=0;iWin<nW;iWin++) {
//                 cerr <<WC[iWin][WC_gStart]<<"   " <<WC[iWin][WC_gEnd] <<"   "<<WC[iWin][WC_Str] <<endl;
//             };
//             exit(1);
//         };
//     };
    
    
    nW=0; //number of windows
    for (uint iP=0; iP<nP; iP++) {//scan through all anchor pieces, create alignment windows

//          if (PC[iP][PC_Nrep]<=P->winAnchorMultimapNmax || PC[iP][PC_Length]>=readLength[PC[iP][PC_iFrag]] ) {//proceed if piece is an anchor, i.e. maps few times or is long enough
       if (PC[iP][PC_Nrep]<=P->winAnchorMultimapNmax ) {//proceed if piece is an anchor, i.e. maps few times
            
            uint aDir   = PC[iP][PC_Dir];     
            uint aLength= PC[iP][PC_Length];            

            for (uint iSA=PC[iP][PC_SAstart]; iSA<=PC[iP][PC_SAend]; iSA++) {//scan through all alignments of this piece

                uint a1 = SA[iSA];
                uint aStr = a1 >> P->GstrandBit;           
                a1 &= P->GstrandMask; //remove strand bit

                //convert to positive strand
                if (aDir==1 && aStr==0) {
                    aStr=1;
                } else if (aDir==0 && aStr==1) {
                    a1 = P->nGenome - (aLength+a1);
                } else if (aDir==1 && aStr==1) {
                    aStr=0;
                    a1 = P->nGenome - (aLength+a1);         
                };

                //final strand            
                if (revertStrand) { //modified strand according to user input CHECK!!!!
                    aStr=1-aStr;
                };   

                if (a1>=P->sjGstart) {//this is sj align
                    uint a1D, aLengthD, a1A, aLengthA, sj1;              
                    if (sjAlignSplit(a1, aLength, P, a1D, aLengthD, a1A, aLengthA, sj1)) {//align crosses the junction

                        int addStatus=createExtendWindowsWithAlign(a1D, aStr);//add donor piece
                        if (addStatus==EXIT_createExtendWindowsWithAlign_TOO_MANY_WINDOWS) {//too many windows
                            break;
                        };         
                        addStatus=createExtendWindowsWithAlign(a1A, aStr);//add acceptor piece
                        if (addStatus==EXIT_createExtendWindowsWithAlign_TOO_MANY_WINDOWS) {//too many windows
                            break;
                        };                       
                    };
                } else {//this is a normal genomic read
                    int addStatus=createExtendWindowsWithAlign(a1, aStr);
                    if (addStatus==EXIT_createExtendWindowsWithAlign_TOO_MANY_WINDOWS) {//too many windows
                        break;
                    };                
                };           
            }; //for (uint iSA=PC[iP][PC_SAstart]; iSA<=PC[iP][PC_SAend]; iSA++) //scan through all alignments of this piece
        };//if (PC[iP][PC_Nrep]<=P->winAnchorMultimapNmax) //proceed if anchor
    };//for (uint iP=0; iP<nP; iP++) //scan through all anchor pieces, create alignment windows
    
    
    for (uint iWin=0;iWin<nW;iWin++) {//extend windows with flanks
        if (WC[iWin][WC_gStart]<=WC[iWin][WC_gEnd]) {//otherwise the window is dead
                       
            uint wb=WC[iWin][WC_gStart];
            for (uint ii=0; ii<P->winFlankNbins && wb>0 && P->chrBin[(wb-1) >> P->winBinChrNbits]==WC[iWin][WC_Chr];ii++) {
                wb--;
                winBin[ WC[iWin][WC_Str] ][ wb ]=(uintWinBin) iWin;
            };
            WC[iWin][WC_gStart] = wb;
            
            wb=WC[iWin][WC_gEnd];
            for (uint ii=0; ii<P->winFlankNbins && wb+1<P->winBinN && P->chrBin[(wb+1) >> P->winBinChrNbits]==WC[iWin][WC_Chr];ii++) {
                wb++;
                winBin[ WC[iWin][WC_Str] ][ wb ]=(uintWinBin) iWin;
            };
            WC[iWin][WC_gEnd] = wb;
            
          
        };
        nWA[iWin]=0; //initialize nWA
        WALrec[iWin]=0; //initialize rec-length        
        WlastAnchor[iWin]=-1;
    };
    
    nWall=nW;
    
    for (uint iP=0; iP<nP; iP++) {//scan through all pieces/aligns, add them to alignment windows, create alignment coordinates
        uint aNrep=PC[iP][PC_Nrep];
        uint aFrag=PC[iP][PC_iFrag];  
        uint aLength=PC[iP][PC_Length];      
        uint aDir=PC[iP][PC_Dir];     
        
        bool aAnchor=(aNrep<=P->winAnchorMultimapNmax); //this align is an anchor or not            

        for (uint ii=0;ii<nW;ii++) {//initialize nWAP
            nWAP[ii]=0;
        };
        
        for (uint iSA=PC[iP][PC_SAstart]; iSA<=PC[iP][PC_SAend]; iSA++) {//scan through all alignments

            uint a1 = SA[iSA];
            uint aStr = a1 >> P->GstrandBit;           
            a1 &= P->GstrandMask; //remove strand bit
            uint aRstart=PC[iP][PC_rStart];

            //convert to positive strand
            if (aDir==1 && aStr==0) {
                aStr=1;
                aRstart = Lread - (aLength+aRstart);                
            } else if (aDir==0 && aStr==1) {
                aRstart = Lread - (aLength+aRstart);                                
                a1 = P->nGenome - (aLength+a1);
            } else if (aDir==1 && aStr==1) {
                aStr=0;
                a1 = P->nGenome - (aLength+a1);         
            };
     
            //final strand            
            if (revertStrand) { //modified strand according to user input CHECK!!!!
                aStr=1-aStr;
            };             

            
            if (a1>=P->sjGstart) {//this is sj read
                uint a1D, aLengthD, a1A, aLengthA, isj1;              
                if (sjAlignSplit(a1, aLength, P, a1D, aLengthD, a1A, aLengthA, isj1)) {//align crosses the junction

                        assignAlignToWindow(a1D, aLengthD, aStr, aNrep, aFrag, aRstart, aAnchor, isj1);
                        assignAlignToWindow(a1A, aLengthA, aStr, aNrep, aFrag, aRstart+aLengthD, aAnchor, isj1);
                        
                    } else {//align does not cross the junction
                        continue; //do not check this align, continue to the next one
                    };
                    
                } else {//this is a normal genomic read
                    assignAlignToWindow(a1, aLength, aStr, aNrep, aFrag, aRstart, aAnchor, -1);
                };
        };
        
//         for (uint ii=0;ii<nW;ii++) {//check of some pieces created too many aligns in some windows, and remove those from WA (ie shift nWA indices
//             if (nWAP[ii]>P->seedNoneLociPerWindow) nWA[ii] -= nWAP[ii];
//         };
    };
    
    //TODO remove windows that have too many alignments
    //aligns are still sorted by original read coordinates, change direction for negative strand
    // DOES NOT HELP!!!
//     for ( uint iW=0;iW<nW;iW++ ) {
//         if (WA[iW][0][WA_rStart]>WA[iW][nWA[iW]-1][WA_rStart]) {//swap
//             for (uint iA=0;iA<nWA[iW]/2;iA++) {
//                 for (uint ii=0;ii<WA_SIZE;ii++) {
//                     uint dummy=WA[iW][iA][ii];
//                     WA[iW][iA][ii]=WA[iW][nWA[iW]-1-iA][ii];
//                     WA[iW][nWA[iW]-1-iA][ii]=dummy;
//                 };
//             };
//         };
//     };
    
#define PacBio    
#if defined PacBio
    if (P->swMode==1) {//stitching is done with Smith-Waterman against the windows
    uint swWinCovMax=0;
        for (uint iWin=0;iWin<nW;iWin++) {//check each window
            swWinCov[iWin]=0;
            if (nWA[iWin]>0) {
                //select good windows by coverage
                uint rLast=0;
                swWinGleft[iWin]=P->chrStart[P->nChrReal]; swWinGright[iWin]=0; swWinRleft[iWin]=0; swWinRright[iWin]=0;
                
                for (uint ia=0; ia<nWA[iWin]; ia++) {//calculate coverage from all aligns
                    uint L1=WA[iWin][ia][WA_Length];
                    uint r1=WA[iWin][ia][WA_rStart];
                    
                    //record ends
                    swWinRright[iWin]=max(swWinRright[iWin],r1+L1-1);
                    swWinGleft[iWin]=min(swWinGleft[iWin],WA[iWin][ia][WA_gStart]);
                    swWinGright[iWin]=max(swWinGright[iWin],WA[iWin][ia][WA_gStart]+L1-1);;

                    
                    if (r1+L1>rLast+1) {
                        if (r1>rLast) {
                            swWinCov[iWin] += L1;
                        } else {
                            swWinCov[iWin] += r1+L1-(rLast+1);
                        };
                        rLast=r1+L1-1;
                    };                    
                };//for (uint ia=0; ia<nWA[iWin]; ia++)
                
                swWinCovMax=max(swWinCovMax,swWinCov[iWin]);
            };//if (nWA[iWin]>0)
        };//for (uint iWin=0;iWin<nW;iWin++)
        
        //debug: read correct loci
        uint trStart,trEnd,trStr,trChr;
        char oneChar;
        istringstream stringStream1;
        stringStream1.str(readName);
        stringStream1 >> oneChar >> trChr >> oneChar >>trStart >> oneChar >> trEnd >> oneChar >>trStr;
        trStart += P->chrStart[trChr];
        trEnd   += P->chrStart[trChr];
        
        uint trNtotal=0, iW1=0;
        trBest = trNext = trInit; //initialize next/best
        for (uint iWin=0;iWin<nW;iWin++) {//check each window
                if (swWinCov[iWin]*100/swWinCovMax >= P->swWinCoverageMinP) {//S-W on all good windows, record the transcripts
                    //full S-W against the window
                    trA=*trInit; //that one is initialized
                    trA.Chr = WC[iWin][WC_Chr];
                    trA.Str = WC[iWin][WC_Str];
                    trA.roStr = revertStrand ? 1-trA.Str : trA.Str; //original strand of the read
                    trA.maxScore=0;

                    uint winLeft =swWinGleft[iWin] -5000;
                    uint winRight=swWinGright[iWin]+5000;
                    
                    //debug: process only correct windows
                    if (!( winLeft<trStart && winRight>trEnd && trA.Str==trStr) ) continue;
                    
                    intSWscore swScore=alignSmithWaterman(R[trA.roStr==0 ? 0:2],Lread,G+winLeft,winRight-winLeft,\
                            (intSWscore) 200, (intSWscore) 200, (intSWscore) 200, (intSWscore) 1, swT, P->swHsize, trA);
                    
                    trA.maxScore = (uint) swScore;
                    
                    trA.cStart  = trA.exons[0][EX_G] + winLeft - P->chrStart[trA.Chr];                     
                    trA.gLength = trA.exons[trA.nExons-1][EX_G]+1;
                    for (uint ii=0;ii<trA.nExons;ii++) {
                        trA.exons[ii][EX_G]+=winLeft;
                    };
                    
//                         uint gg=trA.exons[ii][EX_G]-(trA.exons[ii-1][EX_G]+trA.exons[ii-1][EX_L]);                        
//                         uint rg=trA.exons[ii][EX_R]-trA.exons[ii-1][EX_R]-trA.exons[ii-1][EX_L];
//                         if (gg>P->alignIntronMin) {
//                             trA.canonSJ[ii-1]=0; //sj
//                         } else if (gg>0) {
//                             trA.canonSJ[ii-1]=-1;//deletion
//                         };
//                         if (rg>0) trA.canonSJ[ii-1]=-2;                    
                    
                    trA.rLength=1;
                    trA.nMatch=1;

                    trAll[iW1]=trArrayPointer+trNtotal;
                    *(trAll[iW1][0])=trA;
                    
                    if (trAll[iW1][0]->maxScore > trBest->maxScore || (trAll[iW1][0]->maxScore == trBest->maxScore && trAll[iW1][0]->gLength < trBest->gLength ) ) {
                        trNext=trBest;
                        trBest=trAll[iW1][0];
                    };                    
                    
                    
                    nWinTr[iW1]=1;
                    trNtotal++;
                    iW1++;
                    //output all windows
                    P->inOut->logMain << iRead <<"\t"<< swWinCov[iWin]*100/Lread <<"\t"<< WC[iWin][WC_Str]<<"\t"<< WC[iWin][WC_gStart] \
                            <<"\t"<< WC[iWin][WC_gEnd] <<"\t"<< WA[iWin][0][WA_rStart] <<"\t"<< swWinRright[iWin] \
                            <<"\t"<<swWinGleft[iWin] <<"\t"<< swWinGright[iWin]<<"\t"<< swWinGright[iWin]-swWinGleft[iWin]<<"\t"<<Lread<<"\t"<<trA.maxScore<<endl;                
//                     outputTranscript(&trA, nW, &P->inOut->outBED);                    
                };            
        }; 
        nW=iW1;//number of windows with recorded transcripts
        return;
    };//if (P->swMode==1)
#endif //#if defined PacBio

#ifdef COMPILE_FOR_LONG_READS
uint swWinCovMax=0;
for (uint iW=0;iW<nW;iW++) {//check each window
    swWinCov[iW]=0;
    if (nWA[iW]>0) {
        //select good windows by coverage
        uint rLast=0;

        for (uint ia=0; ia<nWA[iW]; ia++) {//calculate coverage from all aligns
            uint L1=WA[iW][ia][WA_Length];
            uint r1=WA[iW][ia][WA_rStart];

            if (r1+L1>rLast+1) {
                if (r1>rLast) {
                    swWinCov[iW] += L1;
                } else {
                    swWinCov[iW] += r1+L1-(rLast+1);
                };
                rLast=r1+L1-1;
            };                    
        };//for (uint ia=0; ia<nWA[iW]; ia++)

        if (swWinCov[iW]>swWinCovMax) swWinCovMax=swWinCov[iW];
    };//if (nWA[iW]>0)
};//for (uint iW=0;iW<nW;iW++)
for (uint iW=0;iW<nW;iW++) {
    if (swWinCov[iW]<swWinCovMax*10/10) {//remove windows that are not good enough
        nWA[iW]=0;
    } else {//merge pieces that are adjacent in R- and G-spaces
        uint ia1=0;
        for (uint ia=1; ia<nWA[iW]; ia++) {
            if ( WA[iW][ia][WA_rStart] == (WA[iW][ia1][WA_rStart]+WA[iW][ia1][WA_Length]) \
              && WA[iW][ia][WA_gStart] == (WA[iW][ia1][WA_gStart]+WA[iW][ia1][WA_Length]) \
              && WA[iW][ia][WA_iFrag]  ==  WA[iW][ia1][WA_iFrag]     ) {//merge
                
                WA[iW][ia1][WA_Length] += WA[iW][ia][WA_Length];
                WA[iW][ia1][WA_Anchor]=max(WA[iW][ia1][WA_Anchor],WA[iW][ia][WA_Anchor]);
                //NOTE: I am not updating sjA and Nrep fields - this could cause trouble in some cases
                
            } else {//do not merge
                ia1++;
                if (ia1!=ia) {//move from ia to ia1
                    for (uint ii=0; ii<WA_SIZE; ii++) {
                        WA[iW][ia1][ii]=WA[iW][ia][ii];
                    };
                };
            };
        };
        nWA[iW]=ia1+1;
    };
};

//mapping time initialize
std::time(&timeStart);
#endif
    
    
    //generate transcript for each window, choose the best
    trInit->nWAmax=0;
    trBest = trNext = trInit; //initialize next/best
    uint iW1=0;//index of non-empty windows
    uint trNtotal=0; //total number of recorded transcripts
   
    for (uint iW=0; iW<nW; iW++) {//transcripts for all windows

        if (nWA[iW]==0) continue; //the window does not contain any aligns because it was merged with other windows        

//         {//debug
//             for (uint ii=0;ii<nWA[iW];ii++) {
//                         cout << iRead+1 <<"\t"<< WA[iW][ii][WA_Length] <<"\t"<< WA[iW][ii][WA_rStart] <<"\t"<< WA[iW][ii][WA_gStart] << "\n";
//             };
//             continue;
//             cout << nWA[iW]<<" "<<swWinCov[iW]*100/Lread <<"   "<<flush;
//     //         continue;
//         
//         };
        
        if (WlastAnchor[iW]<nWA[iW]) {
            WA[ iW ][ WlastAnchor[iW] ][ WA_Anchor]=2; //mark the last anchor
        };
        
        for (uint ii=0;ii<nWA[iW];ii++) WAincl[ii]=false; //initialize mask
        
        trInit->nWAmax=max(nWA[iW],trInit->nWAmax);        
        trA=*trInit; //that one is initialized
        trA.Chr = WC[iW][WC_Chr];
        trA.Str = WC[iW][WC_Str];
        trA.roStr = revertStrand ? 1-trA.Str : trA.Str; //original strand of the read
        trA.maxScore=0;
        
        trAll[iW1]=trArrayPointer+trNtotal;
        if (trNtotal+P->alignTranscriptsPerWindowNmax > P->alignTranscriptsPerReadNmax) {
            P->inOut->logMain << "WARNING: not enough space allocated for transcript. Did not process all windows for read "<< readName+1 <<endl;
            P->inOut->logMain <<"   SOLUTION: increase alignTranscriptsPerReadNmax and re-run\n" << flush;
            break;
        };
        *(trAll[iW1][0])=trA;
        nWinTr[iW1]=0; //initialize number of transcripts per window
        
        
    #ifdef COMPILE_FOR_LONG_READS
        stitchWindowSeeds(iW, iW1, R[trA.roStr==0 ? 0:2], Q[trA.roStr], G);
    #else
        stitchWindowAligns(0, nWA[iW], 0, WAincl, 0, 0, trA, Lread, WA[iW], R[trA.roStr==0 ? 0:2], Q[trA.roStr], G, sigG, P, trAll[iW1], nWinTr+iW1, this);
    #endif
        trAll[iW1][0]->nextTrScore= nWinTr[iW1]==1 ? 0 : trAll[iW1][1]->maxScore;        
        
        if (trAll[iW1][0]->maxScore > trBest->maxScore || (trAll[iW1][0]->maxScore == trBest->maxScore && trAll[iW1][0]->gLength < trBest->gLength ) ) {
            trNext=trBest;
            trBest=trAll[iW1][0];
        };

        trNtotal += nWinTr[iW1];        
        iW1++;
    };
    
    nW=iW1;//only count windows that had alignments
    
//     {//debug
//         std::time(&timeFinish);
//         double timeDiff=difftime(timeFinish,timeStart);
//         cout << "     "<< timeDiff << "     "<<trBest->maxScore*100/Lread<<"   "<<iRead<<endl;;
//     };
    
    if (trBest->maxScore==0) {//no window was aligned (could happen if for all windows too many reads are multiples)
        mapMarker = MARKER_NO_GOOD_WINDOW;
        nW=0;
        nTr=0;
        return;
    };
            
    nextWinScore=trNext->maxScore;
    
    //output chains for out-of-STAR chimeric detection
    #ifdef OUTPUT_localChains
    {
        P->inOut->outLocalChains << readName <<"\t"<< Read0[0] <<"\t"<< Read0[1] << "\n";
        for (uint iw=0; iw<nW; iw++) {
            for (uint itr=0;itr<nWinTr[iw];itr++) {
                P->inOut->outLocalChains << trAll[iw][itr]->maxScore<<"\t"<< trAll[iw][itr]->Chr<<"\t"<<trAll[iw][itr]->Str<<"\t"<<trAll[iw][itr]->nExons;
                for (uint ib=0;ib<trAll[iw][itr]->nExons;ib++) {                    
                    P->inOut->outLocalChains <<"\t"<< trAll[iw][itr]->exons[ib][EX_G]-P->chrStart[trAll[iw][itr]->Chr] \
                                             <<"\t"<< trAll[iw][itr]->exons[ib][EX_R] <<"\t"<< trAll[iw][itr]->exons[ib][EX_L];
                };
                P->inOut->outLocalChains <<"\n";
            };
        };
    };
    #endif
    //////////////////// chimeras
    //stich windows => chimeras
    //stich only the best window with one of the lower score ones for now - do not stich 2 lower score windows
    //stitch only one window on each end of the read
    
    if (P->chimSegmentMin>0 && nW>1 && trBest->rLength >= P->chimSegmentMin \
            && ( trBest->exons[trBest->nExons-1][EX_R] + trBest->exons[trBest->nExons-1][EX_L] + P->chimSegmentMin <= Lread \
              || trBest->exons[0][EX_R] >= P->chimSegmentMin ) \
             && trBest->nextTrScore+P->outFilterMultimapScoreRange < trBest->maxScore \
             && trBest->intronMotifs[0]==0 && (trBest->intronMotifs[1]==0 || trBest->intronMotifs[2]==0) ) { 
            //there is unmapped space at the start/end, and the main window is not a multimapping window, and non non-canonical junctions, and consistend junction motif
        int chimScoreBest=0,chimScoreNext=0;
        trChim[0]=*trBest;

        uint roStart1=trBest->Str==0 ? trBest->exons[0][EX_R] : Lread - trBest->exons[trBest->nExons-1][EX_R] - trBest->exons[trBest->nExons-1][EX_L];
        uint roEnd1=trBest->Str==0 ? trBest->exons[trBest->nExons-1][EX_R] + trBest->exons[trBest->nExons-1][EX_L] - 1 : Lread - trBest->exons[0][EX_R] - 1;
        if (roStart1>readLength[0]) roStart1--;
        if (roEnd1>readLength[0]) roEnd1--;
        
        uint chimStrBest=0;
        if (trBest->intronMotifs[1]==0 && trBest->intronMotifs[2]==0) {//strand is undefined
            chimStr=0;
        } else if ( (trBest->Str==0) == (trBest->intronMotifs[1]>0)) {//strand the same as RNA
            chimStr=1;
        } else {//strand opposite to RNA
            chimStr=2;
        };
        
        for (uint iW=0; iW<nW; iW++) {//check all other windows for chimeras
            for (uint iWt=0; iWt<nWinTr[iW]; iWt++){//cycl over transcripts in the window    
                if (trBest!=trAll[iW][0] && iWt>0) break; //for all windows except that of the best transcript - hceck only iWt=0 (best trnascripts)
                if (trBest==trAll[iW][0] && iWt==0) continue;
//                 {//same window
//                     if (iWt==0) continue; //do not check the best transcript itself
//                     if (trBest->exons[0][EX_R]<=trAll[iW][iWt]->exons[0][EX_R]) {
//                         //start of the last Best exon is before end of the first Chim exon
//                         if (trBest->exons[trBest->nExons-1][EX_G]<trAll[iW][iWt]->exons[0][EX_G]+trAll[iW][iWt]->exons[0][EX_L]) continue;
//                     } else {
//                         if (trAll[iW][iWt]->exons[trAll[iW][iWt]->nExons-1][EX_G]<trBest->exons[0][EX_G]+trBest->exons[0][EX_L]) continue;                        
//                     };
//                 };
                
                if (trAll[iW][iWt]->intronMotifs[0]>0) continue; //do not stitch a window to itself, or to a window with non-canonical junctions
                uint chimStr1;
                if (trAll[iW][iWt]->intronMotifs[1]==0 && trAll[iW][iWt]->intronMotifs[2]==0) {//strand is undefined
                    chimStr1=0;
                } else if ( (trAll[iW][iWt]->Str==0) == (trAll[iW][iWt]->intronMotifs[1]>0)) {//strand the same as RNA
                    chimStr1=1;
                } else {//strand opposite to RNA
                    chimStr1=2;
                };            

                if (chimStr!=0 && chimStr1!=0 && chimStr!=chimStr1) continue; //chimeric segments have to have consitent strands

                uint roStart2=trAll[iW][iWt]->Str==0 ? trAll[iW][iWt]->exons[0][EX_R] : Lread - trAll[iW][iWt]->exons[trAll[iW][iWt]->nExons-1][EX_R] - trAll[iW][iWt]->exons[trAll[iW][iWt]->nExons-1][EX_L];
                uint roEnd2=trAll[iW][iWt]->Str==0 ? trAll[iW][iWt]->exons[trAll[iW][iWt]->nExons-1][EX_R] + trAll[iW][iWt]->exons[trAll[iW][iWt]->nExons-1][EX_L] - 1 : Lread - trAll[iW][iWt]->exons[0][EX_R] - 1;
                if (roStart2>readLength[0]) roStart2--;
                if (roEnd2>readLength[0]) roEnd2--;          

                uint chimOverlap = roStart2>roStart1 ?  (roStart2>roEnd1 ? 0 : roEnd1-roStart2+1) : (roEnd2<roStart1 ? 0 : roEnd2-roStart1+1);
                bool diffMates=(roEnd1 < readLength[0] && roStart2 >= readLength[0]) || (roEnd2 < readLength[0] && roStart1 >= readLength[0]);

                //segment lengths && (different mates || small gap between segments)
                if (roEnd1 > P->chimSegmentMin + roStart1 + chimOverlap && roEnd2> P->chimSegmentMin + roStart2 + chimOverlap  \
                    && ( diffMates || ( (roEnd1 + P->maxChimReadGap + 1) >= roStart2 && (roEnd2 + P->maxChimReadGap + 1) >= roStart1 ) ) ) {
                                           //maxChimReadGap=0 in Parameters.cpp

                    int chimScore=trBest->maxScore + trAll[iW][iWt]->maxScore - (int)chimOverlap; //subtract overlap to avoid double counting

                    if (chimScore > chimScoreBest && chimScore >= P->chimScoreMin && chimScore+P->chimScoreDropMax >= (int) (readLength[0]+readLength[1]) ) {
                        trChim[1]=*trAll[iW][iWt];                                      
                        chimScoreNext=chimScoreBest;
                        chimScoreBest=chimScore;
                        trChim[1].roStart = trChim[1].roStr ==0 ? trChim[1].rStart : Lread - trChim[1].rStart - trChim[1].rLength;
                        trChim[1].cStart  = trChim[1].gStart - P->chrStart[trChim[1].Chr];      
                        chimStrBest=chimStr1;
                    } else if (chimScore>chimScoreNext) {//replace the nextscore if it's not the best one and is higher than the previous one
                        chimScoreNext=chimScore;              
                    };
                };
            };//cycle over window transcripts
        };//cyecl over windows
        
        if (chimStr==0) chimStr=chimStrBest;
        
        chimN=0;
        if (chimScoreNext + P->chimScoreSeparation < chimScoreBest) {//report only if chimera is unique

            if (trChim[0].roStart > trChim[1].roStart) swap (trChim[0],trChim[1]);
                        
            uint e0 = trChim[0].Str==1 ? 0 : trChim[0].nExons-1;
            uint e1 = trChim[1].Str==0 ? 0 : trChim[1].nExons-1;
                
            uint chimRepeat0=0,chimRepeat1=0,chimJ0=0,chimJ1=0;
            int chimMotif=0;
            chimN=2;            
            if ( trChim[0].exons[e0][EX_iFrag] != trChim[1].exons[e1][EX_iFrag] ) {//mates bracket the chimeric junction
                chimN=2;
                chimRepeat=0;
                chimMotif=-1;
                if (trChim[0].Str==1) {                    
                    chimJ0=trChim[0].exons[e0][EX_G]-1;
                } else {
                    chimJ0=trChim[0].exons[e0][EX_G]+trChim[0].exons[e0][EX_L];                            
                };            
                if (trChim[1].Str==1) {                    
                    chimJ1=trChim[1].exons[e1][EX_G]-1;
                } else {
                    chimJ1=trChim[1].exons[e1][EX_G]+trChim[1].exons[e1][EX_L];                            
                };                    
            } else {//check and shift chimeric junction if necessary
                if (trChim[0].exons[e0][EX_L]>=P->chimJunctionOverhangMin && trChim[1].exons[e1][EX_L]>=P->chimJunctionOverhangMin ) {//large enough overhanh required
                    uint roStart0 = trChim[0].Str==0 ? trChim[0].exons[e0][EX_R] : Lread - trChim[0].exons[e0][EX_R] - trChim[0].exons[e0][EX_L];
                    uint roStart1 = trChim[1].Str==0 ? trChim[1].exons[e1][EX_R] : Lread - trChim[1].exons[e1][EX_R] - trChim[1].exons[e1][EX_L];
                    
                    uint jR, jRbest=0;
                    int jScore=0,jMotif=0,jScoreBest=-999999,jScoreJ=0;
                    for (jR=0; jR<roStart1+trChim[1].exons[e1][EX_L]-roStart0; jR++) {//scan through the exons to find a canonical junction, and check for mismatches
                        
                        if (jR==readLength[0]) jR++; //skip the inter-mate base
                        
                        char bR=Read1[0][roStart0+jR];
                       
                        char b0,b1;
                        if (trChim[0].Str==0) {
                            b0=G[trChim[0].exons[e0][EX_G]+jR];
                        } else {
                            b0=G[trChim[0].exons[e0][EX_G]+trChim[0].exons[e0][EX_L]-1-jR];
                            if (b0<4) b0=3-b0;
                        };
                        
                        if (trChim[1].Str==0) {
                            b1=G[trChim[1].exons[e1][EX_G]-roStart1+roStart0+jR];
                        } else {
                            b1=G[trChim[1].exons[e1][EX_G]+trChim[1].exons[e1][EX_L]-1+roStart1-roStart0-jR];
                            if (b1<4) b1=3-b1;
                        };
                        
                        if (b0>3 || b1>3 || bR>3) {//chimera is not called if there are Ns in the genome or in the read
                            chimN=0;
                            break;
                        };
                        
                        char b01,b02,b11,b12;
                        if (trChim[0].Str==0) {
                            b01=G[trChim[0].exons[e0][EX_G]+jR+1];
                            b02=G[trChim[0].exons[e0][EX_G]+jR+2];                                
                        } else {
                            b01=G[trChim[0].exons[e0][EX_G]+trChim[0].exons[e0][EX_L]-1-jR-1];
                            if (b01<4) b01=3-b01;
                            b02=G[trChim[0].exons[e0][EX_G]+trChim[0].exons[e0][EX_L]-1-jR-2];
                            if (b02<4) b02=3-b02;                                
                        };      
                        if (trChim[1].Str==0) {
                            b11=G[trChim[1].exons[e1][EX_G]-roStart1+roStart0+jR-1];
                            b12=G[trChim[1].exons[e1][EX_G]-roStart1+roStart0+jR];                                
                        } else {
                            b11=G[trChim[1].exons[e1][EX_G]+trChim[1].exons[e1][EX_L]-1+roStart1-roStart0-jR+1];
                            if (b11<4) b11=3-b11;
                            b12=G[trChim[1].exons[e1][EX_G]+trChim[1].exons[e1][EX_L]-1+roStart1-roStart0-jR];
                            if (b12<4) b12=3-b12;                                
                        };        

                        jMotif=0;                        
                        if (b01==2 && b02==3 && b11==0 && b12==2) {//GTAG
                            if (chimStr!=2) {
                                jMotif=1;
                            };
                        } else if(b01==1 && b02==3 && b11==0 && b12==1) {//CTAC
                            if (chimStr!=1) {
                                jMotif=2;
                            };            
                        };  
                                                    
                        if (bR==b0 && bR!=b1) {
                            jScore++;
                        } else if (bR!=b0 && bR==b1) {
                            jScore--;
                        };
                        
                        jScoreJ =jMotif==0 ? jScore +  P->chimScoreJunctionNonGTAG : jScore ;
                        
                        if ( jScoreJ > jScoreBest || (jScoreJ == jScoreBest && jMotif>0) ) {
                            chimMotif=jMotif;
                            jRbest=jR;
                            jScoreBest=jScoreJ;
                        };
                    };//jR cycle
                    if (chimN>0) {//else the chimera was rejected because of mismatches
                        
                        //shift junction in trChim
                        if (trChim[0].Str==1) {
                            trChim[0].exons[e0][EX_R] +=trChim[0].exons[e0][EX_L]-jRbest-1;
                            trChim[0].exons[e0][EX_G] +=trChim[0].exons[e0][EX_L]-jRbest-1;
                            trChim[0].exons[e0][EX_L]=jRbest+1;                        
                            chimJ0=trChim[0].exons[e0][EX_G]-1;
                        } else {
                            trChim[0].exons[e0][EX_L]=jRbest+1;
                            chimJ0=trChim[0].exons[e0][EX_G]+trChim[0].exons[e0][EX_L];                            
                        };                            
                         
                        if (trChim[1].Str==0) {
                            trChim[1].exons[e1][EX_R] +=roStart0+jRbest+1-roStart1;
                            trChim[1].exons[e1][EX_G] +=roStart0+jRbest+1-roStart1;
                            trChim[1].exons[e1][EX_L]=roStart1+trChim[1].exons[e1][EX_L]-roStart0-jRbest-1;  
                            chimJ1=trChim[1].exons[e1][EX_G]-1;
                        } else {
                            trChim[1].exons[e1][EX_L]=roStart1+trChim[1].exons[e1][EX_L]-roStart0-jRbest-1;  
                            chimJ1=trChim[1].exons[e1][EX_G]+trChim[1].exons[e1][EX_L];  
                        };
                        //find repeats
                        char b0,b1;
                        uint jR;
                        for (jR=0;jR<100;jR++) {//forward check
                            if (trChim[0].Str==0) {
                                b0=G[chimJ0+jR];
                            } else {
                                b0=G[chimJ0-jR];
                                if (b0<4) b0=3-b0;
                            };

                            if (trChim[1].Str==0) {
                                b1=G[chimJ1+1+jR];
                            } else {
                                b1=G[chimJ1-1-jR];
                                if (b1<4) b1=3-b1;
                            };
                            if (b0!=b1) break;
                        };
                        chimRepeat1=jR;
                        for (jR=0;jR<100;jR++) {//reverse check
                            if (trChim[0].Str==0) {
                                b0=G[chimJ0-1-jR];
                            } else {
                                b0=G[chimJ0+1+jR];
                                if (b0<4) b0=3-b0;
                            };

                            if (trChim[1].Str==0) {
                                b1=G[chimJ1-jR];
                            } else {
                                b1=G[chimJ1+jR];
                                if (b1<4) b1=3-b1;
                            };
                            if (b0!=b1) break;
                        };                        
                        chimRepeat0=jR;
                    };//chimN>0
                };//large enough overhang
            };//chimeric junction is within a mate
            
            //debug
//             cout << readName <<"\t"<< (trChim[0].Str==0 ? chimJ1-chimJ0 : chimJ0-chimJ1) << "\t"<< (chimMotif>=0 ? P->alignIntronMax :  P->alignMatesGapMax)<<"\n";
//             cout <<  chimRepeat0 <<"\t"<<trChim[0].exons[e0][EX_L]<<"\n";
            //chimeric alignments output
            if ( chimN==2 && trChim[0].exons[e0][EX_L]>=P->chimJunctionOverhangMin+chimRepeat0 \
                    && trChim[1].exons[e1][EX_L]>=P->chimJunctionOverhangMin+chimRepeat1 \
                    && ( trChim[0].Str!=trChim[1].Str ||  trChim[0].Chr!=trChim[1].Chr \
                    || (trChim[0].Str==0 ? chimJ1-chimJ0+1LLU : chimJ0-chimJ1+1LLU) > (chimMotif>=0 ? P->alignIntronMax :  P->alignMatesGapMax) ) )
            {//unique chimeras only && minOverhang1 
             //&& minOverhang2
             //&& (diff str || diff chr || 
             //|| gap > (alignIntronMax,alignMatesGapMax) ) negative gap = very large # because of uint
                                        
                if (trChim[0].exons[0][EX_iFrag]!=trChim[0].exons[trChim[0].nExons-1][EX_iFrag]) {
                    trChim[0].primaryFlag=true;//paired portion is primary
                    trChim[1].primaryFlag=false;
                } else if (trChim[1].exons[0][EX_iFrag]!=trChim[1].exons[trChim[1].nExons-1][EX_iFrag]) {
                    trChim[1].primaryFlag=true;//paired portion is primary
                    trChim[0].primaryFlag=false;                    
                } else if (trChim[0].exons[0][EX_iFrag]!=trChim[1].exons[0][EX_iFrag]) {
                    trChim[1].primaryFlag=true;//paired portion is primary
                    trChim[0].primaryFlag=true;
                } else  {//two chimeric segments are on the same mate, another mate not mapped
                    trChim[0].primaryFlag=true;//paired portion is primary
                    trChim[1].primaryFlag=false;
                };

                for (uint iTr=0;iTr<chimN;iTr++) {//write all chimeric pieces
                    if (P->readNmates==2) {
                        outputTranscriptSAM(trChim[iTr], chimN, iTr, trChim[1-iTr].Chr, trChim[1-iTr].exons[0][EX_G], (int) (trChim[1-iTr].Str!=trChim[1-iTr].exons[0][EX_iFrag]), &chunkOutChimSAM);
                    } else {
                        outputTranscriptSAM(trChim[iTr], chimN, iTr, -1, -1, -1, &chunkOutChimSAM);
                    };                        
                };        
                //junction + SAMp
                chunkOutChimJunction << P->chrName[trChim[0].Chr] <<"\t"<< chimJ0 - P->chrStart[trChim[0].Chr]+1 <<"\t"<< (trChim[0].Str==0 ? "+":"-") \
                        <<"\t"<< P->chrName[trChim[1].Chr] <<"\t"<< chimJ1 - P->chrStart[trChim[1].Chr]+1 <<"\t"<< (trChim[1].Str==0 ? "+":"-") \
                        <<"\t"<< chimMotif <<"\t"<< chimRepeat0  <<"\t"<< chimRepeat1 <<"\t"<< readName+1 \
                        <<"\t"<< trChim[0].exons[0][EX_G] - P->chrStart[trChim[0].Chr]+1 <<"\t"<< outputTranscriptCIGARp(trChim[0]) \
                        <<"\t"<< trChim[1].exons[0][EX_G] - P->chrStart[trChim[1].Chr]+1 <<"\t"<<  outputTranscriptCIGARp(trChim[1]) <<"\n"; //<<"\t"<< trChim[0].exons[0][EX_iFrag]+1 --- no need for that, since trChim[0] is always on the first mate
            };
        };//chimeric score
    };//chimeric search
};//end of function

