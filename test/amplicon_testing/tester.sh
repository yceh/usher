./vcftofasta ../../scripts/amplicon/wuhCor1.fa $1 > spiked_in.fa
cut -f 2 $1 |tail -n +3 > positions
cut -f 2,5 $1 |tail -n +3 >../../scripts/amplicon/position_muts 
head -n 2 spiked_in.fa |tail -n 1 | ./kmerizer positions |sort -u |awk '!/^$/{print "@" NR "\n" $0 "\n+"; for (c=0;c<length($0);c++) printf "F"; printf "\n"}' >../../scripts/amplicon/fake_reads
pushd ../../scripts/amplicon/
./generate_vcf.sh ./fake_reads
../../build/amplicon -i $2 -v amplicons.vcf -f amplicons.fa -d . | grep DEBUG| cut -f 3-4| sort -u -k1,1n >out
diff out position_muts
