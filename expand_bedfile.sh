#This script is used to expand the STR regions to account for reads that need to be masked to prevent STR detection with insert size calculation.
#estimate insertsize with something like the following command, that uses properly paired reads with high mapQ
#samtools view -q 30 -f66 $1 chr1|head -n100000|awk '{ if ($9 > 0) { N+=1; S+=$9; S2+=$9*$9 }} END { M=S/N; print "n="N" mean="M" stdev="sqrt ((S2-M*M*N)/(N-1))}' > $1.ins

#Pass the bed file to be expanded in variable $2..
./expand_bed.pl $1.ins $2 |sort -k1,1 -k2,2n|bedtools merge -i stdin >$1.expanded.bed
