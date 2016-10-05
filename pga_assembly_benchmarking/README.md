# PGA Assembly Benchmarking

This folder contains a script written by Phase Genomics that compares the output of Proximity Guided Assembly (PGA) scaffolding to the RH map. 

Here is a listing of the files in this directory:
* compareGoats.py: compares contig order of the PGA scaffolds against RH map assignments to determine overall scaffold statistics
* papadum-v5bng-pilon-split2kbgap.condense.rhorder.out: an example file showing RH map assignment for contig inputs to PGA scaffolding.
* lachesis_ordering (directory)
	* Example files showing contig orders determined by [Lachesis](https://github.com/shendurelab/LACHESIS) software.
	* Clusters are ordered by descending size