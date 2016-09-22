# FRC Benchmarking automation
---
**by Sergey Koren**

Usage:

You need to have two files in the run directory in addition to the asms, a **mapping.fofn** and **genomesize.** The genomesize tells you Y-axis size for FRC calculations. The mapping.fofn gives you full path to sequence reads and sequence type. It looks like:

>/data/projects/phillippy/seq/goat/m130131_224754_42132_c100461212550000001523059505101380_s1_p0.fastq **pacbio**
>/data/projects/phillippy/seq/goat/illumina/Goat400_NoIndex_L008_R1_016.fastq.gz **illuminaR1**
>/data/projects/phillippy/seq/goat/BGI/SRR488816 **illumina**

Pacbio is self-explanatory, illuminaR1 means files are named *_R1_* and *_R2_* so it will replace R1 in the name with R2. illumina means named _1 and _2 so it will replace _1 with _2.

1. You run index.sh for each fasta file in your folder so if you have 5 files you'd run *index.sh 1, index.sh 2, .. index.sh 5*

2. Then, you have to map all data: *map.sh \<fasta_file\> \<line number from 1 num in mapping.fofn\>*

3. Then: *runAnalysis.sh \<fasta file\>*

4. Finally: *summarizeStats.sh*

You can submit all the above to SGEgrid so index.sh with array job count == num fasta files and map.sh with array count == num seq files.