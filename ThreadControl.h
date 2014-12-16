#ifndef THREAD_CONTROL_DEF
#define THREAD_CONTROL_DEF

#include "ReadAlignChunk.h"
#include <pthread.h>

#define MAX_chunkOutSAMposition 100000

class ThreadControl {
public:
    pthread_t *threadArray;
    pthread_mutex_t mutexInRead, mutexOutSAM, mutexOutChimSAM, mutexOutChimJunction, mutexOutUnmappedFastx, mutexOutFilterBySJout, mutexStats;
    
    uint chunkInN,chunkOutN;
    
    ThreadControl();
    
    static void* threadRAprocessChunks(void *RAchunk) {
        ( (ReadAlignChunk*) RAchunk )->processChunks();
        return NULL;
    };
};

#endif

