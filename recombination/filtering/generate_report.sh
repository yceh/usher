#!/bin/bash -ex
# Get the raw sequences for all the descendents we want

# Pass in current updated fasta file with all raw sequences
all_sequences_fasta=$1
descendants="allDescendants.txt"
out_fasta="fastas/extractedSeqs.fa"
cores=`grep -c ^processor /proc/cpuinfo`

# Extract raw sequences for all needed descendants 
xzcat "${all_sequences_fasta}" | ./build/faSomeRecords stdin $descendants $out_fasta

python3 analyzerecomb.py -a

# Align raw sequences using mafft
cd fastas/OrderedRecombs
ls . |  parallel -j $cores "mafft --auto {} > ../AlignedRecombs/{} "

# Generate report
cp empty_report.txt report.txt
ls fastas/AlignedRecombs | parallel -j $cores python3 checkmutant.py {.} -r 

cp report.txt data/ 
awk '$19 == "False"' data/report.txt | awk '$14 == "False"' | awk '$11 == "False"' > data/final_report.txt
