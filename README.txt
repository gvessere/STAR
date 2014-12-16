STAR 2.3.0
Spliced Transcripts Alignment to a Reference
© Alexander Dobin, 2009-2013

AUTHOR/SUPPORT
Alex Dobin, dobin@cshl.edu
http://code.google.com/p/rna-star/
https://groups.google.com/d/forum/rna-star


CHANGES 2.3.0 vs 2.2.0

1. STAR can use annotations in the form of GTF and GFF3 files with --sjdbGTFfile /path/to/annotation/file
   GFF3 requires --sjdbGTFtagExonParentTranscript Parent

2. Input files can be processed (e.g. uncompressed) with --readFilesCommand <UncompressCommand>. For example, 
   --readFilesIn Read1.gz Read2.gz  --readFilesCommand zcat 
   will take gzipped files and uncompress them on the fly.

3. Multiple input files can be specified separated by comma:
   --readFilesIn A.read1,B.read1,C.read1 A.read2,B.read2,C.read2

4. Output alignments (Aligned.out.sam) can be filterted by the same criteria (outSJfilter* options) as the junctions in the SJ.out.tab file using 
   --outFilterType BySJout
   This is useful to remove from Aligned.out.sam spurious alignments with junctions supported by too few reads, junctions with large gaps or non-canonical junctions.

5. --outSJfilterIntronMaxVsReadN Gmax1 Gmax2 Gmax3 ... 
   specifies max allowed gap for junctions supported by 1,2,3... reads (for SJ.out.tab) file.
   This is useful to remove large gap junctions supported by too few reads.

6. Multi-mapping read counts per junction are output into column 8 of SJ.out.tab

7. Unmapped reads can be output into separate .fastq (.fasta) files with   
   --outReadsUnmapped Fastx.

8. Quality scores in the SAM output can be transformed with --outQSconversionAdd. For example, 
   --outQSconversionAdd -31 
   will transform Sanger+64 into Sanger+33 representation.

9. The default shared memory option is changed to
   --genomeLoad NoSharedMemory
   This is a safer and more compatible mode. 
   If you want to revert to previous version default, please use --genomeLoad LoadAndKeep

10. SAM output is compatible with Picard tools.

11. STAR can be run on Mac OSX. 
   To compile STAR for Mac OSX, please run 'make STARforMac' inside the source directory.



CHANGES 2.2.0 vs 2.1.4

1. Alignment algoirthm which utilizes annotated splice junctions database (sjdb) was improved 
to increase sensitivity to annotated junctions.
Genomes with sjdb have to be re-generated.

2. The sjdb input file (--sjdbFileChrStartEnd <file>) may include the 4th column that indicates strand (+ or -) of the annotated junctions.
If you have access to splice junction strand infromation it is recommended to include it in the sjdb file.
If no strand information is provided, the strand of annotated junctions will be decided from their motifs.

3. New parameter --alignSJDBoverhangMin (=3 by default) was introduced to define a minimum splice junction overhang for annotated (ajdb) junctions.

4. Alignments with inconsistent intron motif strands are now always filtered out.

5. --outFilterIntronMotifs has the following options:
	None
	RemoveNoncanonical (same as KeepCanonical in the old version)
	RemoveNoncanonicalUnannotated 
The option --outFilterIntronMotifs RemoveNoncanonicalUnannotated is highly recommended if you need to use the STAR alignments with Cufflinks.
Unlike the old version, this option is not set automatically when --outSAMstrandField intronMotif is used.

6. Default --sjdbScore is 2 (instead of 1 in the old version). It allows for better recovery of annotated (sjdb) junctions.



PLANS
I am actively working on the following features:
1. Converting to BAM and sorting the output alignments.
2. Counting reads per transcript/gene.
3. Mapping to personal genomes.

INSTALLATION
To compile, unzip into a directory, and run 
$ make
Standard gnu c++ installation is required, with all the paths to libraries and include files supplied.
For Mac OS X run
$ make STARforMac


HARDWARE/SOFTWARE REQUIREMENTS 
x86-64 compatible processors
64 bit Linux or Mac OS X 
27GB of RAM for human genome 


LIMITATIONS:
----- This release was tested with the default parameters for human and mouse genomes.
Please contact the author for a list of recommended parameters for much larger or much smaller genomes.

