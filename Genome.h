#ifndef GENOME_DEF
#define GENOME_DEF

#include "IncludeDefine.h"
#include "Parameters.h"
#include "PackedArray.h"
class Genome {
    public:
        char *G, *sigG;
        PackedArray SA;
        PackedArray SAi;
        void genomeLoad(Parameters*);
    
    private:
        char *G1; //pointer -200 of G
};
#endif
