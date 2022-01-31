#!/bin/bash
# Temporary script to help with running simulations through ripples filtration pipeline

MAT=$1

#  Inputs: 
#	  catDescendants.txt
#	  catRecombination.txt
python3 recombination/filtering/combineAndGetPVals.py

#  This second script is in place of getAllNodes.py and the awk command after
#  Inputs:
#   combinedCatOnlyBestWithPVals.txt
#		MAT
#  Outputs:
#   nodeToParent.txt
#	  allRelevantNodeNames.txt

build/ripplesUtils $MAT
echo "getAllNodes Completed.  Retrieved all relevant nodes."

# Generate .vcf file still using matUtils
build/matUtils extract -i $MAT -s recombination/filtering/data/allRelevantNodeNames.txt -v recombination/filtering/data/allRelevantNodes.vcf -T 10

python3 recombination/filtering/getABABA.py

python3 recombination/filtering/makeMNK.py   # mnk_no_dups.txt gets automatically written into 3seq_build directory 

#recombination/filtering/3seq_build/3seq -c recombination/filtering/3seq_build/my3seqTable700 # mnk.log (3seq output) written into recombination/filtering 
( cd recombination/filtering && ./3seq_build/3seq -c 3seq_build/my3seqTable700 )

echo "All data (including mnk.log output from 3seq) will be written to recombination/filtering/data"


###################################################################################

#TODO: These scripts rely on final_report.txt -> Need to modify these for simulations to not include final_report.txt
#python3 recombination/filtering/finish_MNK.py

#python3 recombination/filtering checkClusters.py  
#awk '$21 <= .20 {print}' recombination/filtering/data/combinedCatOnlyBestWithPValsFinalReportWithInfSitesNoClusters.txt > recombination/filtering/data/combinedCatOnlyBestWithPValsFinalReportWithInfSitesNoClusters3SeqP02.txt  

# Make sure "optimized-large-radius-pruneCatExcludeB30.usher.no-long-branches.leaves.txt.gz" is present in recombination/filtering for these next two scripts to run 
#python3 recombination/filtering/doNewTieBreakers.py 
#python3 recombination/filtering/removeRedundant.py   

#mv recombination/filtering/data/combinedCatOnlyBestWithPValsFinalReportWithInfSitesNoClustersNewTiebreak3seqP02RussPval005RemoveCircular.txt recombination/filtering/results

echo "Pipeline finished. List of recombinants detected in recombination/filtering/results"



