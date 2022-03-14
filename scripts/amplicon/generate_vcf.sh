#!/bin/bash -ex
# Script to generate VCF files for given FASTQ amplicon reads

path_to_amplicon_reads=$1
aligned_amplicons="amplicons.fa"
output_vcf="amplicons.vcf"

# Align amplicon reads to reference to generate SAM file   (single-end alignment)
./build/minimap2 -ax sr test/amplicons/wuhCor1.fa $path_to_amplicon_reads > test/amplicons/aligned_amplicons.sam

# Convert aligned SAM -> BAM file
#./build/samtools-1.9/samtools view -S -b aligned_amplicons.sam > aln.bam
python3 samtomsa.py

# Call variants from aligned amplicon multifasta
./build/faToVcf -maskSites=problematic_sites_sarsCov2.vcf $aligned_amplicons $output_vcf

echo "Finished generating amplicon vcf."
