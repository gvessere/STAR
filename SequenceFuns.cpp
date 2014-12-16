#include "SequenceFuns.h"

void complementSeqNumbers(char* ReadsIn, char* ReadsOut, uint Lread) {//complement the numeric sequences
    for (uint jj=0;jj<Lread;jj++) {
        switch (int(ReadsIn[jj])){
            case (3): ReadsOut[jj]=char(0);break;
            case (2): ReadsOut[jj]=char(1);break;              
            case (1): ReadsOut[jj]=char(2);break;                          
            case (0): ReadsOut[jj]=char(3);break;
            default:  ReadsOut[jj]=ReadsIn[jj];
        };
    };
};

void revComplementNucleotides(char* ReadsIn, char* ReadsOut, uint Lread) {//complement the numeric sequences
    for (uint jj=0;jj<Lread;jj++) {
        switch (int(ReadsIn[Lread-1-jj])){
            case (65): ReadsOut[jj]=char(84);break;
            case (67): ReadsOut[jj]=char(71);break;              
            case (71): ReadsOut[jj]=char(67);break;                          
            case (84): ReadsOut[jj]=char(65);break;
            
            case (97):  ReadsOut[jj]=char(116);break;
            case (99):  ReadsOut[jj]=char(103);break;
            case (103): ReadsOut[jj]=char(99);break;
            case (116): ReadsOut[jj]=char(97);break;
            
            default:   ReadsOut[jj]=ReadsIn[Lread-1-jj];
        };
    };
};

void convertNucleotidesToNumbers(const char* R0, char* R1, uint Lread) {//transform sequence  from ACGT into 0-1-2-3 code    
    for (uint jj=0;jj<Lread;jj++) {
                    switch (int(R0[jj])){
                        case (65): case(97):  R1[jj]=char(0);break;//A
                        case (67): case(99):  R1[jj]=char(1);break;//C           
                        case (71): case(103): R1[jj]=char(2);break;//G                       
                        case (84): case(116): R1[jj]=char(3);break;//T                                
//                         case (78): R1[jj]=char(9);break;//N
                        default:   R1[jj]=char(9);//anything else
                    };
                };
};

uint chrFind(uint Start, uint i2, uint* chrStart) {// find chromosome from global locus
    uint i1=0, i3;
    while (i1+1<i2) {
        i3=(i1+i2)/2;
        if ( chrStart[i3] > Start ) {
            i2=i3;
        } else {
            i1=i3;
        };
    };
    return i1;
};

uint localSearch(const char *x, uint nx, const char *y, uint ny, double pMM){
    //find the best alignment of two short sequences x and y
    //pMM is the maximum percentage of mismatches
    uint nMatch=0, nMM=0, nMatchBest=0, nMMbest=0, ixBest=nx;            
    for (uint ix=0;ix<nx;ix++) {
        nMatch=0; nMM=0;
        for (uint iy=0;iy<min(ny,nx-ix);iy++) {
            if (x[ix+iy]>3) continue;
            if (x[ix+iy]==y[iy]) {
                nMatch++;
            } else {
                nMM++;
            };
        };
        
        if ( ( nMatch>nMatchBest || (nMatch==nMatchBest && nMM<nMMbest) ) && double(nMM)/double(nMatch)<=pMM) {
            ixBest=ix;
            nMatchBest=nMatch;
            nMMbest=nMM;
        };
    };
    return ixBest;
};

uint qualitySplit(char* r, char* q, uint L, char Qsplit, uint maxNsplit, uint  minLsplit, uint** splitR) {
    //splits the read r[L] by quality scores q[L], outputs in splitR - split coordinate/length - per base
    //returns number of good split regions
    uint iR=0,iS=0,iR1,LgoodMin=0, iFrag=0;
    while ( (iR<L) & (iS<maxNsplit) ) { //main cycle
        //find next good base
        while ( ( (q[iR]<Qsplit) || (r[iR]>3) ) && (iR<L) ) {
            if (r[iR]==MARK_FRAG_SPACER_BASE) iFrag++; //count read fragments
            iR++;
        };
        
        if (iR==L) break; //exit when reached end of read
        
        iR1=iR;
        
        //find the next bad base
        while ( ( (q[iR]>=Qsplit) && (r[iR]<=3) ) && (iR<L) ) {
            iR++;
        };        
               
        if ( (iR-iR1)>LgoodMin ) LgoodMin=iR-iR1;
        if ( (iR-iR1)<minLsplit ) continue; //too short for a good region
        
        splitR[0][iS]=iR1;      //good region start
        splitR[1][iS]=iR-iR1;   //good region length
        splitR[2][iS]=iFrag;    //good region fragment
        iS++;
    };
    
    if (iS==0) splitR[1][0]=LgoodMin; //output min good piece length
    
    return iS;
};

