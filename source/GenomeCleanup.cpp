#include "GenomeCleanup.h"

Genome * _genomeReference;


void GenomeCleanUp()
{
    if (_genomeReference != NULL)
        _genomeReference->~Genome();

    _genomeReference=NULL;
};

void (*GetGenomeCleanupFunction(Genome & genome)) ()
{
	if (_genomeReference != NULL)
		_genomeReference = &genome;

	return &GenomeCleanUp;
};