#include "IncludeDefine.h"
#include "Parameters.h"
#include "SuffixArraysFuns.h"
#include "PackedArray.h"
#include <math.h>
#include "TimeFunctions.h"
#include "ErrorWarning.h"
#include "loadGTF.h"
#include "SjdbClass.h"

#include "serviceFuns.cpp"
#include "streamFuns.h"

char* globalG;
uint globalL;


inline int funCompareSuffixes ( const void *a, const void *b){
	uint jj=0LLU;
        
	uint *ga=(uint*)((globalG-7LLU)+(*((uint*)a)));
	uint *gb=(uint*)((globalG-7LLU)+(*((uint*)b)));
    uint va=0,vb=0;

	while (va==vb && jj<globalL) {
		va=*(ga-jj);
		vb=*(gb-jj);
		jj++;
	};

	if (va>vb) {
		return 1;
	} else if (va==vb) {
		return 0;
	} else {
		return -1;
	};
};

inline bool funCompareSuffixesBool ( const void *a, const void *b) 
{
	uint jj=0LLU;
        
	uint *ga=(uint*)((globalG-7LLU)+(*((uint*)a)));
	uint *gb=(uint*)((globalG-7LLU)+(*((uint*)b)));
    uint va=0,vb=0;

	while (va==vb && jj<globalL) {
		va=*(ga-jj);
		vb=*(gb-jj);
		jj++;
	};

	if (va<vb) {
        return true;
	} else {
		return false;
	};
};


inline uint funG2strLocus (uint SAstr, uint const N, char const GstrandBit, uint const GstrandMask) {
    bool strandG = (SAstr>>GstrandBit) == 0;
    SAstr &= GstrandMask;
    if ( !strandG ) SAstr += N;
    return SAstr;
};

uint genomeScanFastaFiles (Parameters *P, char* G, bool flagRun) {//scans fasta files. flagRun=false: check and find full size, flaRun=true: collect all the data
    uint N=0; //total number of bases in the genome, including chr "spacers"
    ifstream fileIn;
    for (uint ii=0;ii<P->genomeFastaFiles.size();ii++) {//all the input files
        fileIn.open(P->genomeFastaFiles.at(ii).c_str());
        if ( !fileIn.good() ) {//
            ostringstream errOut;
            errOut << "EXITING because of INPUT ERROR: could not open genomeFastaFile: " <<P->genomeFastaFiles.at(ii) <<endl;
            exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);            
        };
        while(!fileIn.eof()) {//read each file until eof
            string lineIn (4096,'.');
            getline(fileIn,lineIn);
            if (lineIn[0]=='>') {//new chromosome
                if (!flagRun) {
                    istringstream lineInStream (lineIn);
                    lineInStream.ignore(1,' ');
                    string chrName1;
                    lineInStream >> chrName1;
                    P->chrName.push_back(chrName1);
                };
                
                if (!flagRun && P->chrStart.size()>0) P->chrLength.push_back(N-P->chrStart.at(P->chrStart.size()-1)); //true length of the chr  
                
                if (N>0) {//pad the chromosomes to bins boudnaries
                    N = ( (N+1)/P->genomeChrBinNbases+1 )*P->genomeChrBinNbases;
                };

                if (!flagRun) {
                    P->chrStart.push_back(N);    
                    P->inOut->logMain << P->genomeFastaFiles.at(ii)<<" : chr # " << P->chrStart.size()-1 << "  \""<<P->chrName.at(P->chrStart.size()-1)<<"\" chrStart: "<<N<<"\n"<<flush;
                };
            } else {//char lines
                if (flagRun) lineIn.copy(G+N,lineIn.size(),0);
                N += lineIn.size();
            };
        };
        fileIn.close();        
    };
    
   
    if (!flagRun) P->chrLength.push_back(N-P->chrStart.at(P->chrStart.size()-1)); //true length of the chr  

    N = ( (N+1)/P->genomeChrBinNbases+1)*P->genomeChrBinNbases;
        
    if (!flagRun) { 
        P->nChrReal=P->chrStart.size();
        P->chrStart.push_back(N); //last chromosome end
        for (uint ii=0;ii<P->nChrReal;ii++) {
            P->chrNameIndex[P->chrName[ii]]=ii;
        };
    };
    
    return N;
};        

void radixPass(PackedArray &SA, char* G, uint N, uint K, uint shiftG, uint* groupStart, uint* SAtemp, uint* c, uint* cc) 
{ // count occurrences 
    for (uint ii = 0;  ii < K;  ii++) cc[ii] = 0;         // reset counters
    
    for (uint ii = 0;  ii < N;  ii++) cc[int(G[SA[ii]+shiftG])]++;    // count occurences; -s for reverse, +s for forward!
    
    c[0]=0;
    for (uint ii = 0;  ii < K;  ii++) c[ii+1] = c[ii]+cc[ii]; // exclusiive prefix sums
    
    for (uint ii = 0;  ii < K+1;  ii++) groupStart[ii] = c[ii]; // copy c for output
    
    for (uint ii = 0;  ii < N;  ii++) SAtemp[c[int(G[SA[ii]+shiftG])]++] = SA[ii];      // sort; -s for reverse, +s for forward!
    
    for (uint ii = 0;  ii < N;  ii++) SA.writePacked(ii,SAtemp[ii]);
}


void genomeGenerate(Parameters *P) {
    
    //check parameters
    if (P->sjdbOverhang==0 && (P->sjdbFileChrStartEnd!="-" || P->sjdbGTFfile!="-")) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL INPUT PARAMETER ERROR: for generating genome with annotations (--sjdbFileChrStartEnd or --sjdbGTFfile options)\n";
        errOut << "you need to specify non-zero --sjdbOverhang\n";
        errOut << "SOLUTION: re-run genome generation specifying non-zero --sjdbOverhang, which ideally should be equal to OneMateLength-1, or could be chosen generically as ~100\n";        
        exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);
    };

    
    //time
    time_t rawTime;
    string timeString;
    
    time(&rawTime);
    P->inOut->logMain     << timeMonthDayTime(rawTime) <<" ... Starting to generate Genome files\n" <<flush;
    *P->inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... Starting to generate Genome files\n" <<flush;
    
    //define some parameters from input parameters
    P->genomeChrBinNbases=1LLU << P->genomeChrBinNbits;
    
    
    //write genome parameters file
    ofstream genomePar((P->genomeDir+("/genomeParameters.txt")).c_str());
    if (genomePar.fail()) {//
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR: could not create output file "<< P->genomeDir+("/genomeParameters.txt") << endl;
        errOut << "Solution: check that the genomeDir directory exists and you have write permission for it: genomeDir="<< P->genomeDir << endl <<flush;        
        exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);
    };    
    
    genomePar << "versionGenome\t" << P->versionSTAR <<endl;
    genomePar << "genomeFastaFiles\t";
    for (uint ii=0;ii<P->genomeFastaFiles.size();ii++) genomePar << P->genomeFastaFiles.at(ii) << " ";
    genomePar << endl;
    genomePar << "genomeSAindexNbases\t" << P->genomeSAindexNbases << endl;
    genomePar << "genomeChrBinNbits\t" << P->genomeChrBinNbits << endl;
    genomePar << "genomeSAsparseD\t" << P->genomeSAsparseD <<endl;
    genomePar << "sjdbOverhang\t" << P->sjdbOverhang <<endl;
    genomePar << "sjdbFileChrStartEnd\t" << P->sjdbFileChrStartEnd <<endl;
    
    genomePar.close();    
    
    //add the sjdb sequences to the genome
    SjdbClass sjdbLoci;
    
   
    if (P->sjdbOverhang>0 && P->sjdbFileChrStartEnd!="-") {       
        ifstream sjdbStreamIn ( P->sjdbFileChrStartEnd.c_str() );   
        if (sjdbStreamIn.fail()) {
            ostringstream errOut;
            errOut << "FATAL error, could not open file sjdbFileChrStartEnd=" << P->sjdbFileChrStartEnd <<"\n";
            exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);
        };

        while (sjdbStreamIn.good()) {
            string oneLine,str1;
            uint u1,u2;
            char c1;
            getline(sjdbStreamIn,oneLine);
            istringstream oneLineStream (oneLine);
            oneLineStream >> str1 >> u1 >> u2 >> c1;
            if (str1!="") {
                sjdbLoci.chr.push_back(str1);
                sjdbLoci.start.push_back(u1);
                sjdbLoci.end.push_back(u2);
                sjdbLoci.str.push_back(c1);
            };
        };         
        
        P->inOut->logMain << "Loaded database junctions from file: " << P->sjdbFileChrStartEnd <<": "<<sjdbLoci.chr.size()<<" junctions\n\n";
        
    }; //if (P->sjdbFileChrStartEnd!="-")

    char *G=NULL, *G1=NULL;        
    uint NbasesChrReal=genomeScanFastaFiles(P,G,false);//first scan the fasta file to fins all the sizes  
    P->chrBinFill();
            
    loadGTF(sjdbLoci, P);    
    
    uint L=10000;//maximum length of genome suffix    
    uint nG1alloc=(NbasesChrReal + sjdbLoci.chr.size()*P->sjdbLength+L)*2;
    G1=new char[nG1alloc];
    G=G1+L;
    
    char Kchr=5; // full alphabet size, indexable alphabet size, chromosome end mark
    memset(G1,Kchr,nG1alloc);//initialize to K-1 all bytes
 
    genomeScanFastaFiles(P,G,true);    //load the genome sequence
     
    //convert the genome to 0,1,2,3,4
    for (uint jj=0;jj<NbasesChrReal;jj++) {
        switch (int(G[jj])){
            case(65): case(97):  G[jj]=char(0);break;//A
            case(67): case(99):  G[jj]=char(1);break;//C           
            case(71): case(103): G[jj]=char(2);break;//G                       
            case(84): case(116): G[jj]=char(3);break;//T                                
            case(78): case(110): G[jj]=char(4);break;//N
            case(48):            G[jj]=Kchr;break;//chromosomal breaks within the sequences
            default:              //anything else
                if (G[jj]!=Kchr) {
//                     P->inOut->logMain << "Unexpected character: char="<< G[jj] << "   int="<<int(G[jj])<<"   at " << jj << " , replacing with N\n";
                     G[jj]=char(4);                                 
                };
        };
    };    

    if (sjdbLoci.chr.size()>0) {//prepare sjdb
        uint *sjdbS=new uint [sjdbLoci.chr.size()];
        uint *sjdbE=new uint [sjdbLoci.chr.size()];
        
        uint8 *sjdbMotif=new uint8 [sjdbLoci.chr.size()];
        uint8 *sjdbShiftLeft=new uint8 [sjdbLoci.chr.size()];
        uint8 *sjdbShiftRight=new uint8 [sjdbLoci.chr.size()];        
        
        
        string chrOld="";
        uint iChr=0;
        for (uint ii=0;ii<sjdbLoci.chr.size();ii++) {
            if (chrOld!=sjdbLoci.chr.at(ii)) {//find numeric value of the chr
                for (iChr=0;iChr<P->nChrReal;iChr++) {
                    if (sjdbLoci.chr.at(ii)==P->chrName[iChr]) break;
                };
                if (iChr>=P->nChrReal) {
                    ostringstream errOut;                    
                    errOut << "EXITING because of FATAL error, the sjdb chromosome " << sjdbLoci.chr.at(ii) << " is not found among the genomic chromosomes\n";
                    errOut << "SOLUTION: fix your file sjdbFileChrStartEnd=" << P->sjdbFileChrStartEnd <<" at line #" <<ii+1<<"\n";
                    exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);
                };
                chrOld=sjdbLoci.chr.at(ii);
            };
            
            sjdbS[ii] = sjdbLoci.start.at(ii) + P->chrStart[iChr] - 1;//sj names contain 1-based intron loci
            sjdbE[ii] = sjdbLoci.end.at(ii)   + P->chrStart[iChr] - 1;

            //motifs
            if ( G[sjdbS[ii]]==2 && G[sjdbS[ii]+1]==3 && G[sjdbE[ii]-1]==0 && G[sjdbE[ii]]==2 ) {//GTAG
                sjdbMotif[ii]=1;
            } else if ( G[sjdbS[ii]]==1 && G[sjdbS[ii]+1]==3 && G[sjdbE[ii]-1]==0 && G[sjdbE[ii]]==1 ) {//CTAC
                sjdbMotif[ii]=2;
            } else if ( G[sjdbS[ii]]==2 && G[sjdbS[ii]+1]==1 && G[sjdbE[ii]-1]==0 && G[sjdbE[ii]]==2 ) {//GCAG
                sjdbMotif[ii]=3;
            } else if ( G[sjdbS[ii]]==1 && G[sjdbS[ii]+1]==3 && G[sjdbE[ii]-1]==2 && G[sjdbE[ii]]==1 ) {//CTGC
                sjdbMotif[ii]=4;
            } else if ( G[sjdbS[ii]]==0 && G[sjdbS[ii]+1]==3 && G[sjdbE[ii]-1]==0 && G[sjdbE[ii]]==1 ) {//ATAC
                sjdbMotif[ii]=5;
            } else if ( G[sjdbS[ii]]==2 && G[sjdbS[ii]+1]==3 && G[sjdbE[ii]-1]==0 && G[sjdbE[ii]]==3 ) {//GTAT
                sjdbMotif[ii]=6;             
            } else {
                sjdbMotif[ii]=0;
            };
            //repeat length: go back and forth around jR to find repeat length
            uint jjL=0,jjR=0;
            while ( jjL <= sjdbS[ii]-1 && G[sjdbS[ii]-1-jjL]==G[sjdbE[ii]-jjL] && G[sjdbS[ii]-1-jjL]<4 && jjL<255) {//go back
                jjL++;
            };
            sjdbShiftLeft[ii]=jjL;
            
            while ( sjdbS[ii]+jjR < NbasesChrReal && G[sjdbS[ii]+jjR]==G[sjdbE[ii]+1+jjR] && G[sjdbS[ii]+jjR]<4 && jjR<255) {//go forward
                jjR++;
            };
            sjdbShiftRight[ii]=jjR;
            
            
            if (jjR==255 || jjL==255) {
                P->inOut->logMain << "WARNING: long repeat for junction # " << ii+1 <<" : " \
                        << sjdbLoci.chr.at(ii) <<" "<<sjdbS[ii] - P->chrStart[iChr] + 1 <<" "<< sjdbE[ii] - P->chrStart[iChr] + 1 \
                        << "; left shift = "<< (int) sjdbShiftLeft[ii] <<"; right shift = "<< (int) sjdbShiftRight[ii] <<"\n";
            };            
            
            sjdbS[ii]-=sjdbShiftLeft[ii];
            sjdbE[ii]-=sjdbShiftLeft[ii];
        };
        
        //sort sjdb
        uint *sjdbSort=new uint [sjdbLoci.chr.size()*3];
        for (uint ii=0;ii<sjdbLoci.chr.size();ii++) {   
            sjdbSort[ii*3]=sjdbS[ii]+(sjdbLoci.str.at(ii)=='-' ? NbasesChrReal : 0); //separate sorting of +/- strand
            sjdbSort[ii*3+1]=sjdbE[ii]+(sjdbLoci.str.at(ii)=='-' ? NbasesChrReal : 0);
            sjdbSort[ii*3+2]=ii;
        };
        
        qsort((void *) sjdbSort, sjdbLoci.chr.size(), sizeof(uint)*3, funCompareUint2);
        
        uint *I=new uint [sjdbLoci.chr.size()];
        uint nsj=0;
        for (uint ii=0;ii<sjdbLoci.chr.size();ii++) {
            uint isj=sjdbSort[ii*3+2];//index of the next sorted junction            
            if (nsj==0 || sjdbS[isj]!=sjdbS[I[nsj-1]] || sjdbE[isj]!=sjdbE[I[nsj-1]]) {//add new junction
                I[nsj++]=isj;
            } else if ( (sjdbMotif[isj]>0 && sjdbMotif[I[nsj-1]]==0) \
                      ||( ((sjdbMotif[isj]>0) == (sjdbMotif[I[nsj-1]]>0)) && sjdbShiftLeft[isj]<sjdbShiftLeft[I[nsj-1]])) {//replace the old junctions
                //canonical or left-most junctions junction win
                I[nsj-1]=isj;
            };
        };
        
        //sort again, after returning canonical junctions back to original loci:
        for (uint ii=0;ii<nsj;ii++) {   
            sjdbSort[ii*3]  =sjdbS[I[ii]] + (sjdbMotif[I[ii]]==0 ? 0 : sjdbShiftLeft[I[ii]]);
            sjdbSort[ii*3+1]=sjdbE[I[ii]] + (sjdbMotif[I[ii]]==0 ? 0 : sjdbShiftLeft[I[ii]]);
            sjdbSort[ii*3+2]=I[ii];
        };
        
        qsort((void *) sjdbSort, nsj, sizeof(uint)*3, funCompareUint2);
        
        P->sjdbStart=new uint [nsj];
        P->sjdbEnd=new uint [nsj];
        P->sjdbMotif=new uint8 [nsj];
        P->sjdbShiftLeft=new uint8 [nsj];
        P->sjdbShiftRight=new uint8 [nsj];    
        P->sjdbStrand=new uint8 [nsj];  
        
        uint nsj1=0;
        for (uint ii=0;ii<nsj;ii++) {
            bool sjReplace=false;
            uint isj=sjdbSort[ii*3+2];
            if ( nsj1>0 && P->sjdbStart[nsj1-1]==sjdbSort[ii*3] && P->sjdbEnd[nsj1-1]==sjdbSort[ii*3+1] ) {//same loci on opposite strands
                if (P->sjdbMotif[nsj1-1]>0 || (P->sjdbMotif[nsj1-1]==0 && sjdbMotif[isj]==0)) {//old sj is canonical, or both are non-canonical (on opposite strand)
                    P->sjdbStrand[nsj1-1]=0;
                    continue;
                } else {//replace the junction
                    nsj1--;
                    sjReplace=true;
                };
            };
            P->sjdbStart[nsj1]=sjdbSort[ii*3];
            P->sjdbEnd[nsj1]=sjdbSort[ii*3+1];
            P->sjdbMotif[nsj1]=sjdbMotif[isj];
            P->sjdbShiftLeft[nsj1]=sjdbShiftLeft[isj];                    
            P->sjdbShiftRight[nsj1]=sjdbShiftRight[isj];
            if (sjdbLoci.str.at(isj)=='+') {
                P->sjdbStrand[nsj1]=1;
            } else if (sjdbLoci.str.at(isj)=='-') {
                P->sjdbStrand[nsj1]=2;
            } else {
                if (P->sjdbMotif[nsj1]==0) {//strand un-defined
                    P->sjdbStrand[nsj1]=0;
                } else {
                    P->sjdbStrand[nsj1]=2-P->sjdbMotif[nsj1]%2;
                };
            };
            if (sjReplace) P->sjdbStrand[nsj1]=0;
            nsj1++;
        };            
        P->sjdbN=nsj1;       
        P->sjDstart = new uint [P->sjdbN];
        P->sjAstart = new uint [P->sjdbN];

        ofstream sjdbInfo((P->genomeDir+"/sjdbInfo.txt").c_str());
        //first line is some general useful information
        sjdbInfo << P->sjdbN <<"\t"<< P->sjdbOverhang <<"\n";
        uint sjGstart=P->chrStart[P->nChrReal];
        for (uint ii=0;ii<P->sjdbN;ii++) {            
            //add sjdb sequence to genome   
            P->sjDstart[ii]   = P->sjdbStart[ii]  - P->sjdbOverhang; 
            P->sjAstart[ii]   = P->sjdbEnd[ii] + 1;     
            if (P->sjdbMotif[ii]==0) {//shinon-canonical junctions back to their true coordinates
                P->sjDstart[ii] += P->sjdbShiftLeft[ii];
                P->sjAstart[ii] += P->sjdbShiftLeft[ii];
            };            
            memcpy(G+sjGstart,G+P->sjDstart[ii],P->sjdbOverhang);//sjdbStart contains 1-based intron loci
            memcpy(G+sjGstart+P->sjdbOverhang,G+P->sjAstart[ii],P->sjdbOverhang);//sjdbStart contains 1-based intron loci
            sjGstart += P->sjdbLength;     
            sjdbInfo << P->sjdbStart[ii] <<"\t"<< P->sjdbEnd[ii] <<"\t"<<(int) P->sjdbMotif[ii] <<"\t"<<(int) P->sjdbShiftLeft[ii] <<"\t"<<(int) P->sjdbShiftRight[ii]<<"\t"<<(int) P->sjdbStrand[ii] <<"\n";
        };
        sjdbInfo.close();
        time ( &rawTime );
        P->inOut->logMain     << timeMonthDayTime(rawTime) <<" ... finished processing splice junctions database ...\n" <<flush;   
        *P->inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... finished processing splice junctions database ...\n" <<flush;
        
    };
    
    uint N = NbasesChrReal + P->sjdbN*P->sjdbLength;
    P->nGenome=N;
    uint N2 = N*2;     

    ofstream chrN((P->genomeDir+("/chrName.txt")).c_str());
    ofstream chrS((P->genomeDir+("/chrStart.txt")).c_str());
    ofstream chrL((P->genomeDir+("/chrLength.txt")).c_str());
    ofstream chrNL((P->genomeDir+("/chrNameLength.txt")).c_str());
    
    for (uint ii=0;ii<P->nChrReal;ii++) {//output names, starts, lengths               
        chrN<<P->chrName[ii]<<endl;
        chrS<<P->chrStart[ii]<<endl;
        chrL<<P->chrLength.at(ii)<<endl;
        chrNL<<P->chrName[ii]<<"\t"<<P->chrLength.at(ii)<<endl;        
    };
    chrS<<P->chrStart[P->nChrReal]<<endl;//size of the genome
    chrN.close();chrL.close();chrS.close(); chrNL.close();   
    
    if (P->limitGenomeGenerateRAM < (nG1alloc+nG1alloc/3)) {//allocate nG1alloc/3 for SA generation
        ostringstream errOut;                            
        errOut <<"EXITING because of FATAL PARAMETER ERROR: limitGenomeGenerateRAM="<< (P->limitGenomeGenerateRAM) <<"is too small for your genome\n";
        errOut <<"SOLUTION: please specify limitGenomeGenerateRAM not less than"<< nG1alloc+nG1alloc/3 <<" and make that much RAM available \n";
        exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);
    };    
    
    //write genome to disk
    ofstream genomeOut((P->genomeDir+("/Genome")).c_str());    
    if (genomeOut.fail()) {//
        ostringstream errOut;                    
        errOut << "FATAL ERROR: could not create output file=Genome, EXITING\n"<<"\n";
        exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);
    };
    P->inOut->logMain << "Writing genome to disk...";
    fstreamWriteBig(genomeOut,G,N);
    genomeOut.close();    
    P->inOut->logMain << " done.\n" <<flush;
      
    
    //preparing to generate SA
    
    for (uint ii=0;ii<N;ii++) {//- strand
        G[N2-1-ii]=G[ii]<4 ? 3-G[ii] : G[ii];
    };      
    
    P->nSA=0;
    for (uint ii=0;ii<N2;ii+=P->genomeSAsparseD) {
        if (G[ii]<4) {
            P->nSA++;
        };
    };     
    
    P->GstrandBit = (uint) floor(log(N)/log(2))+1; 
    if (P->GstrandBit<32) P->GstrandBit=32; //TODO: use simple access function for SA
    
    P->GstrandMask = ~(1LLU<<P->GstrandBit);
    P->nSAbyte=P->nSA*(P->GstrandBit+1)/8+1;
    PackedArray SA1;    
    SA1.defineBits(P->GstrandBit+1,P->nSA);
        
    P->inOut->logMain  << "Number of SA indices: "<< P->nSA << endl<<flush;    
    P->inOut->logMain  << "SA size in bytes: "<< P->nSAbyte << endl<<flush;
    

    //sort SA
        
    time ( &rawTime );
    P->inOut->logMain     << timeMonthDayTime(rawTime) <<" ... starting to sort  Suffix Array. This may take a long time...\n" <<flush;   
    *P->inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... starting to sort  Suffix Array. This may take a long time...\n" <<flush;

    for (uint ii=0;ii<N;ii++) {//re-fill the array backwards for sorting
        swap(G[N2-1-ii],G[ii]);
    };          
    globalG=G;
    globalL=L/sizeof(uint);        

    {//not enough RAM, split into chunks          
        //count the number of indices with 4nt prefix
        uint indPrefN=1LLU << 16;
        uint* indPrefCount = new uint [indPrefN];
        memset(indPrefCount,0,indPrefN*sizeof(indPrefCount[0]));
        P->nSA=0;
        for (uint ii=0;ii<N2;ii+=P->genomeSAsparseD) {
            if (G[ii]<4) {
                uint p1=(G[ii]<<12) + (G[ii-1]<<8) + (G[ii-2]<<4) + G[ii-3];
                indPrefCount[p1]++;
                P->nSA++;
            };
        };

        uint saChunkSize=(P->limitGenomeGenerateRAM-nG1alloc)/8/P->runThreadN; //number of SA indexes per chunk
        saChunkSize=saChunkSize*6/10; //allow extra space for qsort            
        //uint saChunkN=((P->nSA/saChunkSize+1)/P->runThreadN+1)*P->runThreadN;//ensure saChunkN is divisible by P->runThreadN
        //saChunkSize=P->nSA/saChunkN+100000;//final chunk size
        if (P->runThreadN>1) saChunkSize=min(saChunkSize,P->nSA/(P->runThreadN-1));

        uint saChunkN=P->nSA/saChunkSize;//estimate
        uint* indPrefStart = new uint [saChunkN*2]; //start and stop, *2 just in case
        uint* indPrefChunkCount = new uint [saChunkN*2];
        indPrefStart[0]=0;
        saChunkN=0;//start counting chunks
        uint chunkSize1=indPrefCount[0];
        for (uint ii=1; ii<indPrefN; ii++) {
            chunkSize1 += indPrefCount[ii];
            if (chunkSize1 > saChunkSize) {
                saChunkN++;
                indPrefStart[saChunkN]=ii;
                indPrefChunkCount[saChunkN-1]=chunkSize1-indPrefCount[ii];                    
                chunkSize1=indPrefCount[ii];
            };
        };
        saChunkN++;
        indPrefStart[saChunkN]=indPrefN+1;
        indPrefChunkCount[saChunkN-1]=chunkSize1;

        P->inOut->logMain  << "Number of chunks: " << saChunkN <<";   chunks size limit: " << saChunkSize*8 <<" bytes\n" <<flush;

        time ( &rawTime );
        P->inOut->logMain     << timeMonthDayTime(rawTime) <<" ... sorting Suffix Array chunks and saving them to disk...\n" <<flush;   
        *P->inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... sorting Suffix Array chunks and saving them to disk...\n" <<flush;

        #pragma omp parallel for num_threads(P->runThreadN) ordered schedule(dynamic,1)
        for (int iChunk=0; (uint)iChunk < saChunkN; iChunk++) {//start the chunk cycle: sort each chunk with qsort and write to a file
            uint* saChunk=new uint [indPrefChunkCount[iChunk]];//allocate local array for each chunk
            for (uint ii=0,jj=0;ii<N2;ii+=P->genomeSAsparseD) {//fill the chunk with SA indices
                if (G[ii]<4) {
                    uint p1=(G[ii]<<12) + (G[ii-1]<<8) + (G[ii-2]<<4) + G[ii-3];
                    if (p1>=indPrefStart[iChunk] && p1<indPrefStart[iChunk+1]) {
                        saChunk[jj]=ii;
                        jj++;
                    };
                    //TODO: if (jj==indPrefChunkCount[iChunk]) break;
                };
            };


            //sort the chunk
            qsort(saChunk,indPrefChunkCount[iChunk],sizeof(saChunk[0]),funCompareSuffixes);
            for (uint ii=0;ii<indPrefChunkCount[iChunk];ii++) {    
                saChunk[ii]=N2-1-saChunk[ii];
            };  
            //wrtie files
            ostringstream saChunkFileNameStream("");
            saChunkFileNameStream<< P->genomeDir << "/SA_" << iChunk;
            ofstream saChunkFile(saChunkFileNameStream.str().c_str());
            fstreamWriteBig(saChunkFile, (char*) saChunk, sizeof(saChunk[0])*indPrefChunkCount[iChunk]);
            saChunkFile.close();
            delete [] saChunk;
            saChunk=NULL;
        };

        time ( &rawTime );
        P->inOut->logMain     << timeMonthDayTime(rawTime) <<" ... loading chunks from disk, packing SA...\n" <<flush;   
        *P->inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... loading chunks from disk, packing SA...\n" <<flush;    

        //read chunks and pack into full SA1
        SA1.charArray=new char[P->nSAbyte];
        uint N2bit= 1LLU << P->GstrandBit;          
        uint packedInd=0;

        #define SA_CHUNK_BLOCK_SIZE 10000000
        uint* saIn=new uint[SA_CHUNK_BLOCK_SIZE]; //TODO make adjustable
        
        #ifdef genenomeGenerate_SA_textOutput
                ofstream SAtxtStream ((P->genomeDir + "/SAtxt").c_str());
        #endif

        for (uint iChunk=0;iChunk<saChunkN;iChunk++) {//load files one by one and convert to packed
            ostringstream saChunkFileNameStream("");
            saChunkFileNameStream<< P->genomeDir << "/SA_" << iChunk;
            ifstream saChunkFile(saChunkFileNameStream.str().c_str());
            while (! saChunkFile.eof()) {//read blocks from each file
                uint chunkBytesN=fstreamReadBig(saChunkFile,(char*) saIn,SA_CHUNK_BLOCK_SIZE*sizeof(saIn[0]));
                for (uint ii=0;ii<chunkBytesN/sizeof(saIn[0]);ii++) {
                    SA1.writePacked( packedInd+ii, (saIn[ii]<N) ? saIn[ii] : ( (saIn[ii]-N) | N2bit ) );
                    
                    #ifdef genenomeGenerate_SA_textOutput
                        SAtxtStream << saIn[ii] << "\n";
                    #endif
                };
                packedInd += chunkBytesN/sizeof(saIn[0]);
            };
            saChunkFile.close();
            remove(saChunkFileNameStream.str().c_str());//remove the chunk file
        };

        #ifdef genenomeGenerate_SA_textOutput
                SAtxtStream.close();
        #endif        
        delete [] saIn;

        if (packedInd != P->nSA ) {//
            ostringstream errOut;                            
            errOut << "EXITING because of FATAL problem while generating the suffix array\n";
            errOut << "The number of indices read from chunks = "<<packedInd<<" is not equal to expected nSA="<<P->nSA<<endl;
            errOut << "SOLUTION: try to re-run suffix array generation, if it still does not work, report this problem to the author\n"<<flush;
            exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);
        };
                    
        time ( &rawTime );
        P->inOut->logMain     << timeMonthDayTime(rawTime) <<" ... writing Suffix Array to disk ...\n" <<flush;   
        *P->inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... writing Suffix Array to disk ...\n" <<flush;   
        
        ofstream SAout;
        SAout.open((P->genomeDir+("/SA")).c_str());
        fstreamWriteBig(SAout,(char*) SA1.charArray, (streamsize) P->nSAbyte);
        SAout.close();
        
        //DONE with suffix array generation
        
        for (uint ii=0;ii<N;ii++) {//return to normal order for future use
            swap(G[N2-1-ii],G[ii]);
        };         
        
    };    

    time ( &rawTime );
    timeString=asctime(localtime ( &rawTime ));
    timeString.erase(timeString.end()-1,timeString.end());
    P->inOut->logMain     << timeMonthDayTime(rawTime) <<" ... Finished generating suffix array\n" <<flush;  
    *P->inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ... Finished generating suffix array\n" <<flush;          
    
    //check SA size on disk, must agree with P->nSAbyte
    ifstream genomeIn;
    genomeIn.open((P->genomeDir+("/SA")).c_str());
    genomeIn.seekg (0, ios::end);
    uint nSAbyte1=genomeIn.tellg();
    genomeIn.close();
    if (nSAbyte1>P->nSAbyte) {
        ostringstream errOut;
        errOut << "WARNING: in genomeGenerate: the size of the SA file on disk, "<<nSAbyte1<<", is bigger than expected "<<P->nSAbyte<<"\n";
        errOut << "         Will try to cut the file to a correct size and check the SA\n";
        errOut << "         Please report this error to dobin@cshl.edu\n";
        
        ostringstream sysCom;
        sysCom << "cd " <<P->genomeDir<<"; mv SA SA.old; head -c " << P->nSAbyte << " SA.old > SA; rm -f SA.old";
        errOut << "         Executing system command: " <<sysCom.str() <<"\n";
        system(sysCom.str().c_str());
        
        *P->inOut->logStdOut <<errOut.str();
        P->inOut->logMain <<errOut.str();

        genomeIn.open((P->genomeDir+("/SA")).c_str());
        genomeIn.seekg (0, ios::end);
        uint nSAbyte2=genomeIn.tellg();
        genomeIn.close();    
        if (nSAbyte2!=P->nSAbyte) {
                ostringstream errOut;                                            
                errOut << "EXITING: FATAL ERROR in genomeGenerate: could not write correctly sized SA to disk\n";
                errOut << "SOLUTION: Please report this error to dobin@cshl.edu\n";
                exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);            
        };
        
        genomeIn.open((P->genomeDir+("/SA")).c_str());
        fstreamReadBig(genomeIn,SA1.charArray,P->nSAbyte);
        genomeIn.close();
    } else if (nSAbyte1<P->nSAbyte) {
        ostringstream errOut;                                            
        errOut << "EXITING: FATAL ERROR in genomeGenerate: the size of the SA file on disk, "<<nSAbyte1<<", is smaller than expected "<<P->nSAbyte<<"\n";
        errOut << "SOLUTION: Please try to generate the genome files agains in an empty directory. If the error persists, please report it to dobin@cshl.edu\n";
        exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);            
    };
    
    if (nSAbyte1>P->nSAbyte) {//check SA
        time(&rawTime);            
        timeString=asctime(localtime ( &rawTime ));
        timeString.erase(timeString.end()-1,timeString.end());
        P->inOut->logMain    << timeString <<" ... starting to check Suffix Array...\n" <<flush;   
        *P->inOut->logStdOut << timeString <<" ... starting to check Suffix Array...\n" <<flush;   
          
        uint* g1=new uint;
        uint* g2=new uint;        
        for (uint isa=0;isa<P->nSA-1;isa++) {//check SA      
            if (isa%100000000==0) P->inOut->logMain  << isa*100/P->nSA << " " << flush;


            *g1=funG2strLocus(SA1[isa  ],N,P->GstrandBit,P->GstrandMask);
            *g2=funG2strLocus(SA1[isa+1],N,P->GstrandBit,P->GstrandMask);        
            
            uint jj=0;
            while (G[*g1+jj]==G[*g2+jj] && jj<L) ++jj;

            if ( jj<L && G[*g1+jj]>G[*g2+jj] ) {
                ostringstream errOut;                                            
                errOut << "EXITING: FATAL ERROR in genomeGenerate: Suffix Array is not properly sorted\n";
                errOut << "SOLUTION: re-run genomeGenerate from scratch, with an empty genomeDir directory\n";
                exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);
            }
        };
        P->inOut->logMain << "  done: all suffixes ordered correctly\n"<<flush;
        *P->inOut->logStdOut << " done: all suffixes ordered correctly\n"<<flush;
    };
////////////////////////////////////////
//          SA index
//
    time(&rawTime);    
    timeString=asctime(localtime ( &rawTime ));
    timeString.erase(timeString.end()-1,timeString.end());
    P->inOut->logMain    << timeMonthDayTime(rawTime) <<" ... starting to generate Suffix Array index...\n" <<flush;   
    *P->inOut->logStdOut << timeMonthDayTime(rawTime) <<" ... starting to generate Suffix Array index...\n" <<flush; 
        
    P->genomeSAindexStart = new uint [P->genomeSAindexNbases+1];
    P->genomeSAindexStart[0]=0;
    for (uint ii=1;ii<=P->genomeSAindexNbases;ii++) {//L-mer indices starts
        P->genomeSAindexStart[ii] = P->genomeSAindexStart[ii-1] + ( 1LLU<<(2*ii) );
    };
    P->nSAi = P->genomeSAindexStart[P->genomeSAindexNbases];
    
    uint* SAi=new uint[P->nSAi];
//     for (uint isa=0; isa<P->nSAi; isa++) {//initialize
//         SAi[isa]=P->nSA; //if the suffix is not found in the genome, it's location will be marked with this value
//     };
    
    uint* ind0=new uint[P->genomeSAindexNbases];
    uint* indSAlast=new uint[P->genomeSAindexNbases];

    for (uint ii=0; ii<P->genomeSAindexNbases; ii++) {
        ind0[ii]=-1;//this is needed in case "AAA...AAA",i.e. indPref=0 is not present in the genome for some lengths
        indSAlast[ii]=P->nSA;//that's probably not needed
    };

    P->SAiMarkNbit=P->GstrandBit+1;
    P->SAiMarkAbsentBit=P->GstrandBit+2;
    
    P->SAiMarkNmaskC=1LLU << P->SAiMarkNbit;
    P->SAiMarkNmask=~P->SAiMarkNmaskC;
    P->SAiMarkAbsentMaskC=1LLU << P->SAiMarkAbsentBit;
    P->SAiMarkAbsentMask=~P->SAiMarkAbsentMaskC;
       
    
    for (uint isa=0; isa<P->nSA; isa++) {//for all suffixes
        if (isa%100000000==0) P->inOut->logMain  << isa*100/P->nSA << "% " << flush;         
        
        uint SAstr=SA1[isa];
        bool dirG = (SAstr>>P->GstrandBit) == 0; //forward or reverse strand of the genome
        SAstr &= P->GstrandMask;
        if (!dirG) SAstr=P->nGenome-1-SAstr;

        uint indPref=0;
        for (uint iL=0; iL < P->genomeSAindexNbases; iL++) {//calculate index

            indPref <<= 2;
            
            uint g1= (uint) G[dirG ? SAstr+iL : SAstr-iL]; //reverese if (-) strand

            if (g1>3) {//if N, this suffix does not belong in SAi
                for (uint iL1=iL; iL1 < P->genomeSAindexNbases; iL1++) {
                    SAi[P->genomeSAindexStart[iL1]+ind0[iL1]] |= P->SAiMarkNmaskC;
                };
                break;
            };

            if (!dirG) g1=3-g1; //complement if (-) strand

            indPref += (uint) g1;
            
            if ( indPref > ind0[iL] || isa==0 ) {//new && good index, record it
                SAi[P->genomeSAindexStart[iL]+indPref]=isa;
                for (uint ii=ind0[iL]+1; ii<indPref; ii++) {//index is not present, record to the last present suffix
                    SAi[P->genomeSAindexStart[iL]+ii] = isa | P->SAiMarkAbsentMaskC; 
                };
                ind0[iL]=indPref;
//                 indSAlast[iL]=isa;
//             } else if (indPref==ind0[iL]) {
//                 indSAlast[iL]=isa;//last SA index with the same prefix 
            } else if ( indPref < ind0[iL] ) {
                ostringstream errOut;
                errOut << "BUG: next index is smaller than previous, EXITING\n" <<flush;
                exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);
            };
        };
    };//for (uint isa=0; isa<P->nSA; isa++)
    P->inOut->logMain << " done\n"<<flush;
   
    //pack SAi
    PackedArray SAip;
    SAip.defineBits(P->GstrandBit+3,P->nSAi);//SAi uses an extra bit compared to SA because it needs to store values > nSA
    SAip.pointArray((char*) SAi);
    for (uint ii=0;ii<SAip.length;ii++) {
        SAip.writePacked(ii,SAi[ii]);
    };
    
    time(&rawTime);    
    P->inOut->logMain    << timeMonthDayTime(rawTime) <<" ... writing SAindex to disk\n" <<flush;   
    *P->inOut->logStdOut << timeMonthDayTime(rawTime) <<" ... writing SAindex to disk\n" <<flush;   
    
    
    //write SAi to disk
    genomeOut.open((P->genomeDir+("/SAindex")).c_str());    
    if (genomeOut.fail()) {//
        ostringstream errOut;
        errOut << "FATAL ERROR: could not create output file=SAindex, EXITING\n"<<flush;
        exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);
    };

    fstreamWriteBig(genomeOut, (char*) &P->genomeSAindexNbases, sizeof(P->genomeSAindexNbases));
    fstreamWriteBig(genomeOut, (char*) P->genomeSAindexStart, sizeof(P->genomeSAindexStart[0])*(P->genomeSAindexNbases+1));        
    fstreamWriteBig(genomeOut,  SAip.charArray, SAip.lengthByte);
    genomeOut.close();    
    
    time(&rawTime);
    timeString=asctime(localtime ( &rawTime ));
    timeString.erase(timeString.end()-1,timeString.end());
    
    
    time(&rawTime);        
    P->inOut->logMain    << timeMonthDayTime(rawTime) << " ..... Finished successfully\n" <<flush;    
    *P->inOut->logStdOut << timeMonthDayTime(rawTime) << " ..... Finished successfully\n" <<flush;
    
    
};
