#!/bin/bash -ex

MAT=$1

#  Inputs: (outputs from RIPPLES)
#	  descendants.txt
#	  recombination.txt
python3 recombination/filtering/final_scripts/combineAndGetPVals.py

build/ripplesUtils $MAT
echo "getAllNodes Completed.  Retrieved all relevant nodes."

# Generates allRelevantNodees.vcf
build/matUtils extract -i $MAT -s recombination/filtering/data/allRelevantNodeNames.txt -v recombination/filtering/data/allRelevantNodes.vcf -T 10

python3 recombination/filtering/final_scripts/getABABA.py

python3 recombination/filtering/final_scripts/makeMNK.py   

python3 recombination/filtering/final_scripts/getDescendants.py

python3 recombination/filtering/final_scripts/makeSampleInfo.py

# Get raw sequences for all descendant nodes, align them to reference and generate final_report.txt (more details in report_meaning.txt)
./generate_report.sh
echo "Successfully generated final_report.txt"

# Run 3seq program on mnk_no_duplicates.txt values
( cd recombination/filtering && ./3seq_build/3seq -c 3seq_build/my3seqTable700 )
echo "mnk.log output from 3seq program written to recombination/filtering/data"

python3 recombination/filtering/final_scripts/finish_MNK.py

python3 recombination/filtering/final_scripts/checkClusters.py  
awk '$21 <= .20 {print}' recombination/filtering/data/combinedCatOnlyBestWithPValsFinalReportWithInfSitesNoClusters.txt > recombination/filtering/data/combinedCatOnlyBestWithPValsFinalReportWithInfSitesNoClusters3SeqP02.txt  

#TODO: Generate this file from scratch
# Make sure "optimized-large-radius-pruneCatExcludeB30.usher.no-long-branches.leaves.txt.gz" is present in recombination/filtering for these next two scripts to run 
python3 recombination/filtering/final_scripts/doNewTieBreakers.py 

mkdir -p recombination/filtering/results
python3 recombination/filtering/final_scripts/removeRedundant.py   

echo "Pipeline finished. List of recombinants detected in recombination/filtering/results"
