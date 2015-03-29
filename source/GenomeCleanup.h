#include "Genome.h"

#ifndef GENOME_H
#define GENOME_H

extern Genome * _genomeReference;
extern void GenomeCleanUp();
extern void (*GetGenomeCleanupFunction(Genome & genome)) ();

#endif

