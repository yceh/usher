#!/bin/bash
#
# Get the raw sequences for all the descendents

bucket_id="$1"
date="$2"
reference="$3"
startDir=$PWD
cores=`grep -c ^processor /proc/cpuinfo`

# Check correct number of args passed
if [ "$#" -ne 3 ]; then
    echo "ERROR: Incorrect number of arguments passed."
fi

mkdir -p filtering/fastas
cp $reference filtering/fastas/reference.fa

chmod +x filtering/get_raw_sequences.sh
./filtering/get_raw_sequences.sh -b $bucket_id -d $date -i filtering/data -o filtering/data
cp filtering/data/allDescendants.fa filtering/fastas/extractedSeqs.fa

mkdir -p filtering/fastas/OrderedRecombs
mkdir -p filtering/fastas/AlignedRecombs
python3 filtering/analyzerecomb.py -a

# Align raw sequences using mafft
cd filtering/fastas/OrderedRecombs
ls . |  parallel -j $cores "mafft --auto {} > ../AlignedRecombs/{} "
cd $startDir

# Generate report
cp filtering/empty_report.txt filtering/data/report.txt
ls filtering/fastas/AlignedRecombs | parallel -j $cores python3 filtering/checkmutant.py {.} -r 

awk '$19 == "False"' filtering/data/report.txt | awk '$14 == "False"' | awk '$11 == "False"' > filtering/data/final_report.txt
echo "DONE"
