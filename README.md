# Long-Reads_Analyzer
ONT Analyzer README

ONT Analyzer is a script which takes as input a fastq file containing all the reads coming from the ONT sequencing and performs all the steps of the analysis until the variant calling.

In order to make it run you must have installed conda and all the tools required, which are: 

- PycoQC
- Guppy_barcoder or Porechop
- Filtlong
- Fastqc
- Minimap2 or Flye
- Samtools
- Picard
- CuteSV or Sniffles

Inside the script there are the references to the paths of files, you have to modify them to pass the files correctly. The files needed are :

- Fastq file
- Reference Genome
- Summary file of the run

Inside the Snakemake pipeline folder is present a different version of the script that is recommended because of the enhanced efficiency. It makes usage of snakemake tool to run the same tools, so you still have the same requirements. It comes with a config file where you have to modify the paths and also in the snakefile you have to personalize them.
