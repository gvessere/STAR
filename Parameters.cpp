#include "Parameters.h"
#include "ErrorWarning.h"
#include "SequenceFuns.h"
#include "OutSJ.h"
#include "sysRemoveDir.h"
#include "stringSubstituteAll.h"

//for mkfifo
#include <sys/stat.h>

#define PAR_NAME_PRINT_WIDTH 30

Parameters::Parameters() {//initalize parameters info
    
    inOut = new InOutStreams;
    
    //versions
    parArray.push_back(new ParameterInfoScalar <uint> (-1, -1, "versionSTAR", &versionSTAR));
    parArray.push_back(new ParameterInfoVector <uint> (-1, -1, "versionGenome", &versionGenome));
    
    //parameters
    parArray.push_back(new ParameterInfoVector <string> (-1, 2, "parametersFiles", &parametersFiles));
    
    //run
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "runMode", &runMode));
    parArray.push_back(new ParameterInfoScalar <int> (-1, -1, "runThreadN", &runThreadN));        

    //genome
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "genomeDir", &genomeDir));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "genomeLoad", &genomeLoad));        
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "genomeFastaFiles", &genomeFastaFiles));
    parArray.push_back(new ParameterInfoScalar <uint> (-1, -1, "genomeSAindexNbases", &genomeSAindexNbases));        
    parArray.push_back(new ParameterInfoScalar <uint> (-1, -1, "genomeChrBinNbits", &genomeChrBinNbits));        
    parArray.push_back(new ParameterInfoScalar <uint> (-1, -1, "genomeSAsparseD", &genomeSAsparseD));        

    //read
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "readFilesIn", &readFilesIn));
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "readFilesCommand", &readFilesCommand));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "readMatesLengthsIn", &readMatesLengthsIn));
    
    //limits
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "limitGenomeGenerateRAM", &limitGenomeGenerateRAM));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "limitIObufferSize", &limitIObufferSize));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "limitOutSAMoneReadBytes", &limitOutSAMoneReadBytes));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "limitOutSJcollapsed", &limitOutSJcollapsed));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "limitOutSJoneRead", &limitOutSJoneRead));

    //output
    parArray.push_back(new ParameterInfoScalar <string>     (-1, 2, "outFileNamePrefix", &outFileNamePrefix)); 
    parArray.push_back(new ParameterInfoScalar <string>     (-1, 2, "outStd", &outStd));        
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "outReadsUnmapped", &outReadsUnmapped));
    parArray.push_back(new ParameterInfoScalar <int>        (-1, -1, "outQSconversionAdd", &outQSconversionAdd));
    
    //outSAM
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "outSAMmode", &outSAMmode));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "outSAMstrandField", &outSAMstrandField));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "outSAMattributes", &outSAMattributes));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "outSAMunmapped", &outSAMunmapped));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "outSAMorder", &outSAMorder));
    parArray.push_back(new ParameterInfoScalar <string>     (-1, -1, "outSAMprimaryFlag", &outSAMprimaryFlag));
    

   //output SJ filtering
    parArray.push_back(new ParameterInfoScalar <string>  (-1, -1, "outSJfilterReads", &outSJfilterReads));
    parArray.push_back(new ParameterInfoVector <int32>   (-1, -1, "outSJfilterCountUniqueMin", &outSJfilterCountUniqueMin));
    parArray.push_back(new ParameterInfoVector <int32>   (-1, -1, "outSJfilterCountTotalMin", &outSJfilterCountTotalMin));    
    parArray.push_back(new ParameterInfoVector <int32>   (-1, -1, "outSJfilterOverhangMin", &outSJfilterOverhangMin));
    parArray.push_back(new ParameterInfoVector <int32>   (-1, -1, "outSJfilterDistToOtherSJmin", &outSJfilterDistToOtherSJmin));
    parArray.push_back(new ParameterInfoVector <int32>   (-1, -1, "outSJfilterIntronMaxVsReadN", &outSJfilterIntronMaxVsReadN));

    //output filtering
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "outFilterType", &outFilterType) );

    parArray.push_back(new ParameterInfoScalar <uint>     (-1, -1, "outFilterMultimapNmax", &outFilterMultimapNmax));
    parArray.push_back(new ParameterInfoScalar <intScore> (-1, -1, "outFilterMultimapScoreRange", &outFilterMultimapScoreRange));
    
    parArray.push_back(new ParameterInfoScalar <intScore> (-1, -1, "outFilterScoreMin", &outFilterScoreMin));
    parArray.push_back(new ParameterInfoScalar <double>   (-1, -1, "outFilterScoreMinOverLread", &outFilterScoreMinOverLread));
    
    parArray.push_back(new ParameterInfoScalar <uint>     (-1, -1, "outFilterMatchNmin", &outFilterMatchNmin));    
    parArray.push_back(new ParameterInfoScalar <double>   (-1, -1, "outFilterMatchNminOverLread", &outFilterMatchNminOverLread));    
    
    parArray.push_back(new ParameterInfoScalar <uint>     (-1, -1, "outFilterMismatchNmax", &outFilterMismatchNmax));
    parArray.push_back(new ParameterInfoScalar <double>   (-1, -1, "outFilterMismatchNoverLmax", &outFilterMismatchNoverLmax));
    parArray.push_back(new ParameterInfoScalar <string>   (-1, -1, "outFilterIntronMotifs", &outFilterIntronMotifs));
  
    //clipping
    parArray.push_back(new ParameterInfoVector <uint>   (-1, -1, "clip5pNbases", &clip5pNbases));
    parArray.push_back(new ParameterInfoVector <uint>   (-1, -1, "clip3pNbases", &clip3pNbases));
    parArray.push_back(new ParameterInfoVector <uint>   (-1, -1, "clip3pAfterAdapterNbases", &clip3pAfterAdapterNbases));
    parArray.push_back(new ParameterInfoVector <string> (-1, -1, "clip3pAdapterSeq", &clip3pAdapterSeq));
    parArray.push_back(new ParameterInfoVector <double> (-1, -1, "clip3pAdapterMMp", &clip3pAdapterMMp));

    //binning, anchors, windows
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "winBinNbits", &winBinNbits));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "winAnchorDistNbins", &winAnchorDistNbins));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "winFlankNbins", &winFlankNbins));
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "winAnchorMultimapNmax", &winAnchorMultimapNmax));
    
    //scoring
    parArray.push_back(new ParameterInfoScalar <intScore>   (-1, -1, "scoreGap", &scoreGap));
    parArray.push_back(new ParameterInfoScalar <intScore>   (-1, -1, "scoreGapNoncan", &scoreGapNoncan));
    parArray.push_back(new ParameterInfoScalar <intScore>   (-1, -1, "scoreGapGCAG", &scoreGapGCAG));
    parArray.push_back(new ParameterInfoScalar <intScore>   (-1, -1, "scoreGapATAC", &scoreGapATAC));
    parArray.push_back(new ParameterInfoScalar <intScore>   (-1, -1, "scoreStitchSJshift", &scoreStitchSJshift));
    parArray.push_back(new ParameterInfoScalar <double>     (-1, -1, "scoreGenomicLengthLog2scale", &scoreGenomicLengthLog2scale));

    parArray.push_back(new ParameterInfoScalar <intScore>   (-1, -1, "scoreDelBase", &scoreDelBase));
    parArray.push_back(new ParameterInfoScalar <intScore>   (-1, -1, "scoreDelOpen", &scoreDelOpen));
    parArray.push_back(new ParameterInfoScalar <intScore>   (-1, -1, "scoreInsOpen", &scoreInsOpen));
    parArray.push_back(new ParameterInfoScalar <intScore>   (-1, -1, "scoreInsBase", &scoreInsBase));   
    
    //alignment    
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "seedSearchLmax", &seedSearchLmax));
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "seedSearchStartLmax", &seedSearchStartLmax));    
    parArray.push_back(new ParameterInfoScalar <double>     (-1, -1, "seedSearchStartLmaxOverLread", &seedSearchStartLmaxOverLread));
    
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "seedPerReadNmax", &seedPerReadNmax));  
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "seedPerWindowNmax", &seedPerWindowNmax));  
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "seedNoneLociPerWindow", &seedNoneLociPerWindow));  
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "seedMultimapNmax", &seedMultimapNmax)); 
    
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "alignIntronMin", &alignIntronMin));        
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "alignIntronMax", &alignIntronMax));        
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "alignMatesGapMax", &alignMatesGapMax));        

    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "alignTranscriptsPerReadNmax", &alignTranscriptsPerReadNmax));    
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "alignSJoverhangMin", &alignSJoverhangMin));    
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "alignSJDBoverhangMin", &alignSJDBoverhangMin));    

    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "alignSplicedMateMapLmin", &alignSplicedMateMapLmin));       
    parArray.push_back(new ParameterInfoScalar <double>     (-1, -1, "alignSplicedMateMapLminOverLmate", &alignSplicedMateMapLminOverLmate));       
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "alignWindowsPerReadNmax", &alignWindowsPerReadNmax));  
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "alignTranscriptsPerWindowNmax", &alignTranscriptsPerWindowNmax));  

    //chimeric
    parArray.push_back(new ParameterInfoScalar <uint>      (-1, -1, "chimSegmentMin", &chimSegmentMin));    
    parArray.push_back(new ParameterInfoScalar <int>       (-1, -1, "chimScoreMin", &chimScoreMin));
    parArray.push_back(new ParameterInfoScalar <int>       (-1, -1, "chimScoreDropMax", &chimScoreDropMax));    
    parArray.push_back(new ParameterInfoScalar <int>       (-1, -1, "chimScoreSeparation", &chimScoreSeparation));    
    parArray.push_back(new ParameterInfoScalar <int>       (-1, -1, "chimScoreJunctionNonGTAG", &chimScoreJunctionNonGTAG));    
    parArray.push_back(new ParameterInfoScalar <uint>       (-1, -1, "chimJunctionOverhangMin", &chimJunctionOverhangMin));    
    
    //sjdb
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "sjdbFileChrStartEnd", &sjdbFileChrStartEnd));
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "sjdbGTFfile", &sjdbGTFfile));    
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "sjdbGTFchrPrefix", &sjdbGTFchrPrefix)); 
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "sjdbGTFfeatureExon", &sjdbGTFfeatureExon)); 
    parArray.push_back(new ParameterInfoScalar <string> (-1, -1, "sjdbGTFtagExonParentTranscript", &sjdbGTFtagExonParentTranscript)); 
    parArray.push_back(new ParameterInfoScalar <uint>   (-1, -1, "sjdbOverhang", &sjdbOverhang));
    parArray.push_back(new ParameterInfoScalar <int>    (-1, -1, "sjdbScore", &sjdbScore));
    
    //SW parameters
    parArray.push_back(new ParameterInfoScalar <uint> (-1, -1, "swMode", &swMode));
    parArray.push_back(new ParameterInfoScalar <uint> (-1, -1, "swWinCoverageMinP", &swWinCoverageMinP));
    
    parameterInputName.push_back("Default");
    parameterInputName.push_back("Command-Line-Initial");
    parameterInputName.push_back("Command-Line");    
    parameterInputName.push_back("genomeParameters.txt");
    
};


void Parameters::inputParameters (int argInN, char* argIn[]) {//input parameters: default, from files, from command line
    
///////// Default parameters
    
    #include "parametersDefault.xxd"
    string parString( (const char*) parametersDefault,parametersDefault_len);
    stringstream parStream (parString);

    scanAllLines(parStream, 0, -1);
    for (uint ii=0; ii<parArray.size(); ii++) {
        if (parArray[ii]->inputLevel<0) {
            ostringstream errOut;
            errOut <<"BUG: DEFAULT parameter value not defined: "<<parArray[ii]->nameString;
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
    };    
    
///////// Initial parameters from Command Line
    
    commandLine=""; 
    string commandLineFile="";    

    if (argInN>1) {//scan parameters from command line
        commandLine += string(argIn[0]);
        for (int iarg=1; iarg<argInN; iarg++) {
            string oneArg=string(argIn[iarg]);
            if (oneArg.at(0)!='-' || oneArg.at(1)!='-') {
                if (oneArg.find_first_of(" \t")!=std::string::npos) {//there is white space in the argument, put "" around
                    oneArg ='\"'  + oneArg +'\"';
                };
                commandLineFile +=' ' + oneArg;
            } else {
                commandLineFile +='\n' + oneArg.substr(2);
            };            
            commandLine += ' ' + oneArg;
        };
        istringstream parStreamCommandLine(commandLineFile);   
        scanAllLines(parStreamCommandLine, 1, 2); //read only initial Command Line parameters 
    };
    
    
    inOut->logMain.open((outFileNamePrefix + "Log.out").c_str());
    if (inOut->logMain.fail()) {
        ostringstream errOut;
        errOut <<"EXITING because of FATAL ERROR: could not create output file: "<<outFileNamePrefix + "Log.out"<<"\n";
        errOut <<"Check if the path " << outFileNamePrefix << " exists and you have permissions to write there\n";
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };
    
    inOut->logMain << "STAR svn revision compiled=" << SVN_VERSION_COMPILED << "\n"<<flush;
    inOut->logMain << "STAR compilation time,server,dir=" << COMPILATION_TIME_PLACE << "\n"<<flush;   
    
    
    
    //define what goes to cout
    if (outStd=="Log") {
        inOut->logStdOut=& std::cout;
    } else if (outStd=="SAM") {
        inOut->logStdOutFile.open((outFileNamePrefix + "Log.std.out").c_str());
        inOut->logStdOut= & inOut->logStdOutFile;
    } else {
        ostringstream errOut;
        errOut <<"EXITING because of FATAL PARAMETER error: outStd="<<outStd <<" is not a valid value of the parameter\n";
        errOut <<"SOLUTION: provide a valid value fot outStd: Log or SAM";
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };  
    
    inOut->logMain << "##### DEFAULT parameters:\n" <<flush;
    for (uint ii=0; ii<parArray.size(); ii++) {
        if (parArray[ii]->inputLevel==0) {
            inOut->logMain << setw(PAR_NAME_PRINT_WIDTH) << parArray[ii]->nameString <<"    "<< *(parArray[ii]) << endl;
        };
    };     
    
    inOut->logMain <<"##### Command Line:\n"<<commandLine <<endl ;

    inOut->logMain << "##### Initial USER parameters from Command Line:\n";                
    for (uint ii=0; ii<parArray.size(); ii++) {
        if (parArray[ii]->inputLevel==1) {
            inOut->logMain << setw(PAR_NAME_PRINT_WIDTH) << parArray[ii]->nameString <<"    "<< *(parArray[ii]) << endl;
        };
    };            
    
///////// Parameters files
    
    if (parametersFiles.at(0) != "-") {//read parameters from a user-defined file
        for (uint ii=0; ii<parametersFiles.size(); ii++) {
            parameterInputName.push_back(parametersFiles.at(ii));
            ifstream parFile(parametersFiles.at(ii).c_str());
            if (parFile.good()) {
                inOut->logMain << "##### USER parameters from user-defined parameters file " <<parametersFiles.at(ii)<< ":\n" <<flush;        
                scanAllLines(parFile, parameterInputName.size()-1, -1);
                parFile.close();        
            } else {
                ostringstream errOut;
                errOut <<"EXITING because of fatal input ERROR: could not open user-defined parameters file " <<parametersFiles.at(ii)<< "\n" <<flush;
                exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
            };
        };
    };
    
///////// Command Line Final
    
    if (argInN>1) {//scan all parameters from command line and override previuos values
        inOut->logMain << "###### All USER parameters from Command Line:\n" <<flush;   
        istringstream parStreamCommandLine(commandLineFile);           
        scanAllLines(parStreamCommandLine, 2, -1);
    };       
    
    inOut->logMain << "##### Finished reading parameters from all sources\n\n" << flush;

    inOut->logMain << "##### Final user re-defined parameters-----------------:\n" << flush;

    ostringstream clFull;
    clFull << argIn[0];
    for (uint ii=0; ii<parArray.size(); ii++) {
        if (parArray[ii]->inputLevel>0) {
            inOut->logMain << setw(PAR_NAME_PRINT_WIDTH) << parArray[ii]->nameString <<"    "<< *(parArray[ii]) << endl;
            if (parArray[ii]->nameString != "parametersFiles" ) {
                clFull << "   --" << parArray[ii]->nameString << " " << *(parArray[ii]);
            };
        };
    };    
    commandLineFull=clFull.str();
    inOut->logMain << "\n-------------------------------\n##### Final effective command line:\n" <<  clFull.str() << "\n";
    
//     parOut.close();
    inOut->logMain << "\n##### Final parameters after user input--------------------------------:\n" << flush;    
//     parOut.open("Parameters.all.out");
    for (uint ii=0; ii<parArray.size(); ii++) {
        inOut->logMain << setw(PAR_NAME_PRINT_WIDTH) << parArray[ii]->nameString <<"    "<< *(parArray[ii]) << endl;
    }; 
//     parOut.close();
   
    inOut->logMain << "----------------------------------------\n\n" << flush;
    
    
    ///////////////////////////////////////// Old variables
    outputBED[0]=0; outputBED[1]=0; //controls the format of starBED output. Will be replaced with HDF5 output
    //chimeric
    maxChimReadGap=0;
    //annot scoring
    annotScoreScale=0;
    annotSignalFile="-";
    //sj-db
    sjdbLength=sjdbOverhang*2+1;


    //splitting
    Qsplit=0;
    maxNsplit=10;
    minLsplit=12;
    minLmap=5;
  
    
    
////////////////////////////////////////////////////// Calculate and check parameters
    outFileTmp=outFileNamePrefix +"_tmp/";
    sysRemoveDir (outFileTmp);
    mkdir (outFileTmp.c_str(),S_IRWXU);
    
    if (runMode=="alignReads") inOut->logProgress.open((outFileNamePrefix + "Log.progress.out").c_str());

    if (outStd=="Log") {
        if (runMode=="alignReads" && outSAMmode != "None") {//open SAM file and write header
            inOut->outSAMfile.open((outFileNamePrefix + "Aligned.out.sam").c_str());
            inOut->outSAM = & inOut->outSAMfile;
        };
    } else if (outStd=="SAM") {
        inOut->outSAM = & std::cout;
    };
    
    //versions
    for (uint ii=0;ii<2;ii++) {
        if (parArray[ii]->inputLevel>0) {
            ostringstream errOut;
            errOut <<"EXITING because of fatal input ERROR: the version parameter "<< parArray[ii]->nameString << " cannot be re-defined by the user\n";
            errOut <<"SOLUTION: please remove this parameter from the command line or input files and re-start STAR\n"<<flush;
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
    };
    
    //run
    if (runThreadN<=0) {
        ostringstream errOut;
        errOut <<"EXITING: fatal input ERROR: runThreadN must be >0, user-defined runThreadN="<<runThreadN<<"\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };
    
    //
    if (outFilterType=="BySJout" && outSAMorder=="PairedKeepInputOrder") {
        ostringstream errOut;
        errOut <<"EXITING: fatal input ERROR: --outFilterType=BySJout is not presently compatible with --outSAMorder=PairedKeepInputOrder\n";
        errOut <<"SOLUTION: re-run STAR without setting one of those parameters. Send a feature request to the Authors\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);        
    };
    //SJ filtering
    for (int ii=0;ii<4;ii++) {
        if (outSJfilterOverhangMin.at(ii)<0) outSJfilterOverhangMin.at(ii)=numeric_limits<int32>::max();
        if (outSJfilterCountUniqueMin.at(ii)<0) outSJfilterCountUniqueMin.at(ii)=numeric_limits<int32>::max();
    };
    
    if (limitGenomeGenerateRAM==0) {//must be >0
        inOut->logMain <<"EXITING because of FATAL PARAMETER ERROR: limitGenomeGenerateRAM=0\n";
        inOut->logMain <<"SOLUTION: please specify a >0 value for limitGenomeGenerateRAM\n"<<flush;
        exit(1);
    } else if (limitGenomeGenerateRAM>1000000000000) {//
        inOut->logMain <<"WARNING: specified limitGenomeGenerateRAM="<<limitGenomeGenerateRAM<<" bytes appears to be too large, if you do not have enough memory the code will crash!\n"<<flush;
    };
            
    readNmates=readFilesIn.size(); //for now the number of mates is defined by the number of input files
    if (runMode=="alignReads" && genomeLoad!="Remove" && genomeLoad!="LoadAndExit") {//open reads files to check if they are present
        
        string readFilesCommandString("");
        if (readFilesCommand.at(0)=="-") {
            if (readFilesIn.at(0).find(',')<readFilesIn.at(0).size()) readFilesCommandString="cat";//concatenate multiple files
        } else {
            for (uint ii=0; ii<readFilesCommand.size(); ii++) readFilesCommandString+=readFilesCommand.at(ii)+"   ";
        };
        
        if (readFilesCommandString=="") {//read from file
            for (uint ii=0;ii<readNmates;ii++) {//open readIn files
                inOut->readIn[ii].open(readFilesIn.at(ii).c_str()); //try to open the Sequences file right away, exit if failed
                if (inOut->readIn[ii].fail()) {
                    ostringstream errOut;
                    errOut <<"EXITING because of fatal input ERROR: could not open readInFile=" << readFilesIn.at(ii) <<"\n";
                    exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);                
                };  
            };    
        } else {//create fifo files, execute pre-processing command
            
            for (uint ii=0;ii<readNmates;ii++) {//open readIn files
                ostringstream sysCom;
                sysCom << outFileTmp <<"tmp.fifo.read"<<ii+1;
                readFilesInTmp.push_back(sysCom.str());
                remove(readFilesInTmp.at(ii).c_str());                
                mkfifo(readFilesInTmp.at(ii).c_str(), S_IRUSR | S_IWUSR );
                string readFilesInString(readFilesIn.at(ii));
                stringSubstituteAll(readFilesInString,",","\" \"");
                string systemString=readFilesCommandString + "   \"" + readFilesInString + "\" > "  + readFilesInTmp.at(ii) + " & ";
                system(systemString.c_str());                
                inOut->readIn[ii].open(readFilesInTmp.at(ii).c_str());
            };
        };
        
        //check sizes of the mate files, if not the same, assume mates are not the same length
        if (readNmates==1) {
            readMatesEqualLengths=true;
        } else if (readNmates > 2){
            inOut->logMain <<"EXITING: because of fatal input ERROR: number of read mates files > 2\n" <<flush;
            exit(1);
        } else if (readMatesLengthsIn=="Equal") {
            readMatesEqualLengths=true;
        } else if (readMatesLengthsIn=="NotEqual") {
            readMatesEqualLengths=false;
        } else {
            inOut->logMain << "EXITING because of FATAL INPUT ERROR: the value of the parameter readMatesLengthsIn=" << readMatesLengthsIn <<" is not among the allowed values: Equal or NotEqual\n"<<flush;
            inOut->logMain << "Solution: restart STAR with correct parameter values\n"<<flush;
            exit(1);
        };

        if ( runMode=="alignReads" && outReadsUnmapped=="Fastx" ) {//open unmapped reads file
            for (uint ii=0;ii<readNmates;ii++) {
                ostringstream ff;
                ff << outFileNamePrefix << "Unmapped.out.mate" << ii+1;
                inOut->outUnmappedReadsStream[ii].open(ff.str().c_str());
            };
        };        
        
        
        if (outFilterType=="Normal") {
            outFilterBySJoutStage=0;
        } else if (outFilterType=="BySJout") {
            outFilterBySJoutStage=1;
        } else {
            ostringstream errOut;
            errOut <<"EXITING because of FATAL input ERROR: unknown value of parameter outFilterType: " << outFilterType <<"\n";
            errOut <<"SOLUTION: specify one of the allowed values: Normal | BySJout\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);                
        };
    };
    
    //in/out buffers
    #define BUFFER_InSizeFraction 0.5
    chunkInSizeBytesArray=(uint) int((limitIObufferSize-limitOutSJcollapsed*Junction::dataSize)*BUFFER_InSizeFraction)/2;
    chunkOutSAMsizeBytes= (uint) int((1.0/BUFFER_InSizeFraction-1.0)*chunkInSizeBytesArray*2.0);
    chunkInSizeBytes=chunkInSizeBytesArray-2*(DEF_readSeqLengthMax+1)-2*DEF_readNameLengthMax;//to prevent overflow
    
    //basic trimming
    if (clip5pNbases.size()==1 && readNmates==2) clip5pNbases.push_back(clip5pNbases[0]);    
    if (clip3pNbases.size()==1 && readNmates==2) clip3pNbases.push_back(clip3pNbases[0]);    
    if (clip3pAfterAdapterNbases.size()==1 && readNmates==2) clip3pAfterAdapterNbases.push_back(clip3pAfterAdapterNbases[0]);    

    //adapter clipping
    if (clip3pAdapterSeq.size()==1 && readNmates==2) clip3pAdapterSeq.push_back(clip3pAdapterSeq[0]);    
    if (clip3pAdapterMMp.size()==1 && readNmates==2) clip3pAdapterMMp.push_back(clip3pAdapterMMp[0]); 
    for (uint ii=0;ii<readNmates;ii++) {
        if (clip3pAdapterSeq.at(ii).at(0)=='-') {// no clipping
            clip3pAdapterSeq.at(ii).assign(""); 
        } else {//clipping
            clip3pAdapterSeqNum[ii]=new char [clip3pAdapterSeq.at(ii).length()];
            convertNucleotidesToNumbers(clip3pAdapterSeq.at(ii).data(),clip3pAdapterSeqNum[ii],clip3pAdapterSeq.at(ii).length());
            //inOut->fastaOutSeqs.open("Seqs.out.fasta");
        };
    };
    
//     #ifdef COMPILE_NO_SHM
//         if (genomeLoad!="NoSharedMemory") {
//             ostringstream errOut;
//             errOut <<"EXITING because of FATAL INPUT ERROR: The code was compiled with NO SHARED MEMORY support, but genomeLoad="<<genomeLoad<<"\n";
//             errOut <<"SOLUTION: run STAR with    --genomeLoad NoSharedMemory    option\n";
//             exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
//         };
//     #endif
    
    //open compilation-dependent streams
    #ifdef OUTPUT_localChains
            inOut->outLocalChains.open((outFileNamePrefix + "LocalChains.out.tab").c_str());
    #endif
    
    inOut->logMain << "Finished loading and checking parameters\n" <<flush;
};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Parameters::scanAllLines (istream &streamIn, int inputLevel,  int inputLevelRequested) {//scan
//     istringstream stringInStream (stringIn);
    string lineIn;
    while (getline(streamIn,lineIn)) {
        scanOneLine(lineIn, inputLevel, inputLevelRequested);
    };
};

int Parameters::scanOneLine (string &lineIn, int inputLevel, int inputLevelRequested) {//scan one line and load the parameters, 
                                                             //0 if comment, 1 if OK
    if (lineIn=="") return 0; //empty line

    istringstream lineInStream (lineIn);
    
    if (inputLevel==0 && ( lineIn.substr(0,1)==" " || lineIn.substr(0,1)=="\t" ) ) return 0;//for Default input spaces also mark comments, for nice formatting    
    
    string parIn("");
    lineInStream >> parIn;
    if (parIn=="" || parIn.substr(0,2)=="//" || parIn.substr(0,1)=="#") return 0; //this is a comment

    uint iPar;
    for (iPar=0; iPar<parArray.size(); iPar++) {
        if (parIn==parArray[iPar]->nameString) {//
            if (inputLevelRequested < 0 || inputLevelRequested == parArray[iPar]->inputLevelAllowed) {
                break;//will read this parameter values                
            } else {
                return 1; //do not read inputs not requested at this level                
            };
        };
    };
    
    string parV("");
    lineInStream >> parV;
    if (parV=="") {//parameter value cannot be empty
        ostringstream errOut;
        errOut << "EXITING: FATAL INPUT ERROR: empty value for paramter \""<< parIn << "\" in input \"" << parameterInputName.at(inputLevel) <<"\"\n";
        errOut << "SOLUTION: use non-empty value for this parameter\n"<<flush;
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };
    
    lineInStream.str(lineIn); lineInStream.clear(); lineInStream >> parIn; //get the correct state of stream, past reading parIn
    
    if (iPar==parArray.size()) {//string is not identified
        ostringstream errOut;
        errOut << "EXITING: FATAL INPUT ERROR: unrecoginzed parameter name \""<< parIn << "\" in input \"" << parameterInputName.at(inputLevel) <<"\"\n";
        errOut << "SOLUTION: use correct parameter name (check the manual)\n"<<flush;
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    } else {//found the corresponding parameter
        if (inputLevel==0 && parArray[iPar]->inputLevel>0) {//this is one of the initial parameters, it was read from Command Line and should not be re-defined
            getline(lineInStream,parV);
            inOut->logMain << setiosflags(ios::left) << setw(PAR_NAME_PRINT_WIDTH) << parArray[iPar]->nameString <<parV<<" ... is RE-DEFINED on Command Line as: " << *(parArray[iPar]) <<"\n";
        } else if (parArray[iPar]->inputLevelAllowed>0 && parArray[iPar]->inputLevelAllowed < inputLevel) {//this is initial parameter and cannot be redefined
            ostringstream errOut;
            errOut << "EXITING: FATAL INPUT ERROR: parameter \""<< parIn << "\" cannot be defined at the input level \"" << parameterInputName.at(inputLevel) << "\"\n";
            errOut << "SOLUTION: define parameter \""<< parIn << "\" in \"" << parameterInputName.at(parArray[iPar]->inputLevelAllowed) <<"\"\n" <<flush;
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        } else if (parArray[iPar]->inputLevel==inputLevel) {//this parameter was already defined at this input level
            ostringstream errOut;
            errOut << "EXITING: FATAL INPUT ERROR: duplicate parameter \""<< parIn << "\" in input \"" << parameterInputName.at(inputLevel) << "\"\n";
            errOut << "SOLUTION: keep only one definition of input parameters in each input source\n"<<flush;
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        } else {//read values
            parArray[iPar]->inputValues(lineInStream);         
            parArray[iPar]->inputLevel=inputLevel;
            if ( inOut->logMain.good() ) {
                inOut->logMain << setiosflags(ios::left) << setw(PAR_NAME_PRINT_WIDTH) << parArray[iPar]->nameString << *(parArray[iPar]);
                if ( parArray[iPar]->inputLevel > 0 ) inOut->logMain <<"     ~RE-DEFINED";
                inOut->logMain << endl;    
            };
        };
    };
    return 0;
};

//////////////////////////////////////////////////////////////////////////////////////////
void Parameters::chrInfoLoad() {//find chrStart,Length,nChr from Genome G
    
    //load chr names
    ifstream chrStreamIn ( (genomeDir+"/chrName.txt").c_str() );   
    if (chrStreamIn.fail()) {
        ostringstream errOut;                            
        errOut << "EXITING because of FATAL error, could not open file " << (genomeDir+"/chrName.txt") <<"\n";
        errOut << "SOLUTION: re-generate genome files with STAR --runMode genomeGenerate\n";
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_INPUT_FILES, *this);
    };
    
    while (chrStreamIn.good()) {
        string chrIn;
        char chrInChar[1000];
        chrStreamIn.getline(chrInChar,1000);
        chrIn=chrInChar;
        if (chrIn=="") break;
        chrName.push_back(chrIn);
    };        
    chrStreamIn.close();
    nChrReal=chrName.size();

    inOut->logMain << "Number of real (reference) chromosmes= " << nChrReal <<"\n"<<flush;
    chrStart.resize(nChrReal+1);
    chrLength.resize(nChrReal);
  
    //load chr lengths
    chrStreamIn.open( (genomeDir+"/chrLength.txt").c_str() );   
    if (chrStreamIn.fail()) {
        ostringstream errOut;                            
        errOut << "EXITING because of FATAL error, could not open file " << (genomeDir+"/chrLength.txt") <<"\n";
        errOut << "SOLUTION: re-generate genome files with STAR --runMode genomeGenerate\n";
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_INPUT_FILES, *this);  
    };
    
    for  (uint ii=0;ii<nChrReal;ii++) {
        chrStreamIn >> chrLength[ii];
    };    
    chrStreamIn.close();    
    
    //load chr starts
    chrStreamIn.open( (genomeDir+"/chrStart.txt").c_str() );   
    if (chrStreamIn.fail()) {
        ostringstream errOut;                            
        errOut << "EXITING because of FATAL error, could not open file " << (genomeDir+"/chrStart.txt") <<"\n";
        errOut << "SOLUTION: re-generate genome files with STAR --runMode genomeGenerate\n";        
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_INPUT_FILES, *this);
    };   
    
    for  (uint ii=0;ii<=nChrReal;ii++) {
        chrStreamIn >> chrStart[ii];
    };    
    chrStreamIn.close();    
  
    //log
    for (uint ii=0; ii<nChrReal;ii++) {
        inOut->logMain << ii+1 <<"\t"<< chrName[ii] <<"\t"<<chrLength[ii]<<"\t"<<chrStart[ii]<<"\n"<<flush;
    };
};

//////////////////////////////////////////////////////////
void Parameters::chrBinFill() {
    genomeChrBinNbases=1LLU<<genomeChrBinNbits;
    chrBinN = chrStart[nChrReal]/genomeChrBinNbases+1;    
    chrBin = new uint [chrBinN];
    for (uint ii=0, ichr=1; ii<chrBinN; ++ii) {
        if (ii*genomeChrBinNbases>=chrStart[ichr]) ichr++;
        chrBin[ii]=ichr-1;
    };
};

