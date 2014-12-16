#include "IncludeDefine.h"
#include "Parameters.h"
#include "ErrorWarning.h"
#include "serviceFuns.cpp"
#include "SjdbClass.h"

#include <map>


#define GTF_exonLoci_size 3
#define GTF_exonTrID(ii) ((ii)*GTF_exonLoci_size)
#define GTF_exonStart(ii) ((ii)*GTF_exonLoci_size+1)
#define GTF_exonEnd(ii) ((ii)*GTF_exonLoci_size+2)

uint loadGTF(SjdbClass &sjdbLoci, Parameters *P) {//load gtf file, add junctions to P->sjdb
    //returns number of added junctions
    if (P->sjdbOverhang>0 && P->sjdbGTFfile!="-") {       
        ifstream sjdbStreamIn ( P->sjdbGTFfile.c_str() );   
        if (sjdbStreamIn.fail()) {
            ostringstream errOut;
            errOut << "FATAL error, could not open file sjdbGTFfile=" << P->sjdbGTFfile <<"\n";
            exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_INPUT_FILES, *P);
        };    
        
        std::map <string,uint> transcriptIDnumber;

        uint exonN=0;
        while (sjdbStreamIn.good()) {//count the number of exons
            string chr1,ddd2,featureType;            
            sjdbStreamIn >> chr1 >> ddd2 >> featureType;
            if (chr1.substr(0,1)!="#" && featureType==P->sjdbGTFfeatureExon) {
                exonN++;
            };
            sjdbStreamIn.ignore(1000000000,'\n'); //ignore the rest of the line
        };
        uint* exonLoci=new uint [exonN*GTF_exonLoci_size];
        char* transcriptStrand = new char [exonN];
        vector <string> transcriptID;

        exonN=0;//re-calculate
        sjdbStreamIn.clear();
        sjdbStreamIn.seekg(0,ios::beg);
        while (sjdbStreamIn.good()) {

            string oneLine,chr1,ddd2,featureType;
            getline(sjdbStreamIn,oneLine);
            istringstream oneLineStream (oneLine);
            
            oneLineStream >> chr1 >> ddd2 >> featureType;
            if (chr1.substr(0,1)!="#" && featureType==P->sjdbGTFfeatureExon) {//exonic line, process
                uint ex1,ex2;
                char str1;
                oneLineStream >> ex1 >> ex2 >> ddd2 >> str1 >> ddd2; //read all fields except the last

                string oneLine1;
		getline(oneLineStream, oneLine1);//get the last field
		replace(oneLine1.begin(),oneLine1.end(),';',' ');//to separate attributes
		replace(oneLine1.begin(),oneLine1.end(),'=',' ');//for GFF3 processing
		oneLineStream.str(oneLine1);
                oneLineStream.clear();

                string trID(""), attr1("");
                while (oneLineStream.good()) {
                    oneLineStream >> attr1;
                    if (attr1==P->sjdbGTFtagExonParentTranscript) {
                        oneLineStream >> trID;
                        trID.erase(remove(trID.begin(),trID.end(),'"'),trID.end());
                        trID.erase(remove(trID.begin(),trID.end(),';'),trID.end());
//                         cout <<trID<<endl;
                    };
                };
                if (trID=="") {//no transcript ID
                    P->inOut->logMain << "WARNING: while processing sjdbGTFfile=" << P->sjdbGTFfile <<": no transcript_id for exon feature for line:\n";
                    P->inOut->logMain << oneLine <<"\n"<<flush;
                } else {
                    transcriptIDnumber.insert(std::pair <string,uint> (trID,(uint) transcriptIDnumber.size()));//insert new element if necessary with a new numeric value
                    if (transcriptID.size() < transcriptIDnumber.size()) transcriptID.push_back(trID);
                    if (str1=='+') {
                       transcriptStrand[transcriptIDnumber.size()-1]=1;
                    } else if (str1=='-') {
                       transcriptStrand[transcriptIDnumber.size()-1]=2;
                    } else {
                       transcriptStrand[transcriptIDnumber.size()-1]=0;
                    };
                };
                
                if (P->sjdbGTFchrPrefix!="-") chr1=P->sjdbGTFchrPrefix + chr1;
                if (P->chrNameIndex.count(chr1)==0) {//chr not in Genome
                    P->inOut->logMain << "WARNING: while processing sjdbGTFfile=" << P->sjdbGTFfile <<": chromosome '"<<chr1<<"' not found in Genome fasta files for line:\n";
                    P->inOut->logMain << oneLine <<"\n"<<flush;          
                } else {//record the exon
                    exonLoci[GTF_exonTrID(exonN)]=transcriptIDnumber[trID];
                    exonLoci[GTF_exonStart(exonN)]=ex1+P->chrStart[P->chrNameIndex[chr1]]-1;
                    exonLoci[GTF_exonEnd(exonN)]=ex2+P->chrStart[P->chrNameIndex[chr1]]-1;
                    exonN++;
                };
            };//if (chr1.substr(0,1)!="#" && featureType=="exon")
        };//
        
        //sort exonLoci by transcript ID and exon coordinates
        qsort((void*) exonLoci, exonN, sizeof(uint)*GTF_exonLoci_size, funCompareUint2);
        
        //make junctions
        uint* sjLoci = new uint [exonN*3];
        uint trIDn=exonLoci[0];
        uint sjN=0;
        for (uint exI=1; exI<exonN; exI++) {
            if (trIDn==exonLoci[GTF_exonTrID(exI)]) {
                uint chr1=P->chrBin[exonLoci[GTF_exonStart(exI)] >> P->genomeChrBinNbits];
                if ( exonLoci[GTF_exonStart(exI)]<=exonLoci[GTF_exonEnd(exI-1)] ) {
                    P->inOut->logMain << "WARNING: while processing sjdbGTFfile=" << P->sjdbGTFfile <<": overlapping exons:\n";
                    P->inOut->logMain << P->chrName[chr1] <<"\t"<< exonLoci[GTF_exonStart(exI-1)]+1-P->chrStart[chr1] << "\t"<< exonLoci[GTF_exonEnd(exI-1)]+1-P->chrStart[chr1]  <<"\n";
                    P->inOut->logMain << P->chrName[chr1] <<"\t"<< exonLoci[GTF_exonStart(exI)]+1-P->chrStart[chr1] << "\t"<< exonLoci[GTF_exonEnd(exI)]+1-P->chrStart[chr1]  <<"\n";                    
                } else {
                    sjLoci[sjN*3]=exonLoci[GTF_exonEnd(exI-1)]+1;
                    sjLoci[sjN*3+1]=exonLoci[GTF_exonStart(exI)]-1;
                    sjLoci[sjN*3+2]=(uint) transcriptStrand[trIDn];
                    sjN++;
                };
            } else {
                trIDn=exonLoci[GTF_exonTrID(exI)];
            };
        };
        
        qsort((void*) sjLoci, sjN, sizeof(uint)*3, funCompareUint2);
        
        char strandChar[3]={'.','+','-'};                
        uint sjdbN1=sjdbLoci.chr.size();
        for (uint ii=0;ii<sjN;ii++) {
            if ( ii==0 || (sjLoci[ii*3]!=sjLoci[(ii-1)*3]) || (sjLoci[ii*3+1]!=sjLoci[(ii-1)*3+1]) || (sjLoci[ii*3+2]!=sjLoci[(ii-1)*3+2]) ) {
                uint chr1=P->chrBin[sjLoci[ii*3] >> P->genomeChrBinNbits];
                sjdbLoci.chr.push_back(P->chrName[chr1]);
                sjdbLoci.start.push_back(sjLoci[ii*3]+1-P->chrStart[chr1]);
                sjdbLoci.end.push_back(sjLoci[ii*3+1]+1-P->chrStart[chr1]);
                sjdbLoci.str.push_back(strandChar[sjLoci[ii*3+2]]);
            };
        };
        
        ofstream sjdbList ((P->genomeDir+"/sjdbList.out.tab").c_str());
        for (uint ii=sjdbN1;ii<sjdbLoci.chr.size(); ii++) {
            sjdbList << sjdbLoci.chr.at(ii)<<"\t"<< sjdbLoci.start.at(ii) << "\t"<< sjdbLoci.end.at(ii)  <<"\t"<< sjdbLoci.str.at(ii)<<"\n";
        };
        sjdbList.close();
        
        P->inOut->logMain << "Processing sjdbGTFfile=" << P->sjdbGTFfile <<", found:\n";
        P->inOut->logMain << "\t\t"  << transcriptIDnumber.size() <<" transcripts\n" << "\t\t"  << exonN << " exons (non-collapsed)\n" << "\t\t"  << sjdbLoci.chr.size()-sjdbN1 << " collapsed junctions\n";
        
        return sjdbLoci.chr.size()-sjdbN1;
    } else {
        return 0;
    };
};
