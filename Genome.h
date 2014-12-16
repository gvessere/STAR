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
        void genomeLoad();

        Genome (Parameters* Pin);
	Genome() {};    
    private:
	Parameters* P;
        key_t shmKey;  
        char *shmStart;
        int shmID;

        char *G1; //pointer -200 of G

        bool GetSharedObjectByKey(key_t shmKey, int * shmID);
        int CreateSharedObject(key_t shmKey, uint64 shmSize);
        int SharedObjectsUseCount(int shmID);
        void * MapSharedObjectToMemory(int shmID);
        const char * GetPosixObjectKey(key_t shmKey);
        struct stat GetSharedObjectInfo(int shmID);
        void RemoveSharedObject(int shmID, void * * ptr, key_t shmKey);
};
#endif
