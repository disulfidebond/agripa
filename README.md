# agripa
AGRIPA: Annotation and Genomic Reference fIle PArsing tool

AGRIPA is a python based tools for merging and formatting annotation files.  It can also be used to create a homologous list of genes in a GTF format, but this functionality is still under development and should be considered unstable.

#####
Requirements

AGRIPA is a command-line tools, that requires python 2.7 to be on the $PATH, with the following python 2.7 modules:
- > sys, subprocess, argparse, re, md5, time
AGRIPA also requires blastn version 2.2.31+ to be on the $PATH, even if homology functionality will not be used.

#####
Installation

Install by copying the python script agripa.py to the $PATH or work folder

#####
Usage

format -i FILENAME.GTF -o FILENAME.GFF
The format option will reformat the GTF file as a GFF version 3 file.  The first -i argument must be the GTF file, with no additional arguments.  GFF to GTF conversion will be added in a future update.  The -o argument with a filename FILENAME must be provided as the unique output filename for the gff file.

homology -i FILENAME.GTF -i FILENAME.fasta -o FILENAME.GFF
The homology option will conduct recursive homology for the desired organism with Homo sapiens.  The first -i argument must be the organism annotation gene list file, and the second -i argument must be the organism *multi-fasta* rna fasta file (important: this is different, and not the same as a nucleotide fasta file, and has the format:
  >geneName Description (name)
  AATGAGACCTA...
This file is *required*, but can be downloaded from ensemble or ncbi, and is usually called an rna fasta file.
The -o option is required, and is the output unique file name for the homologous list of genes.

merge -i ORGANISM1.GTF -i ORGANISM1.FASTA -i ORGANISM2.GTF -i ORGANISM2.FASTA -o FILENAME.GTF
The merge option merges gtf annotation files from the annotation files of two organisms.  The first -i option must be the organism annotation gene list file as a GTF format, and the second -i argument must be the organism fasta file.  The third -i argument must be the next organism annotation file list as a GTF format that will be merged to the first organism, and the fourth -i argument must be the next organism fasta file.  The -o option is the resulting unique output file name.
Note: only two files can be merged at a time with this release.
