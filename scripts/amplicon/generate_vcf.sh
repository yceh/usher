#!/bin/bash -ex
# Script to generate VCF files for given FASTQ amplicon reads

amplicon_reads=$1

covid_reference="wuhCor1.fa"
sam="amplicons.sam"
aligned_amplicons="amplicons.fa"
output_vcf="amplicons.vcf"

# Align amplicon reads to reference to generate SAM file   (single-end alignment)
./build/minimap2 -ax sr $covid_reference $amplicon_reads > $sam

# Convert sam -> msa
python3 samtomsa.py

# Call variants from aligned amplicon multifasta
./build/faToVcf -maskSites=problematic_sites_sarsCov2.vcf $aligned_amplicons $output_vcf

echo "Finished generating amplicon vcf."
