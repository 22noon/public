#Hide selected regions from a BAM file..
#Bam files and BED files should be sorted ..
#First parameter is a BAM file
#second one bed regions
#Third is the output file..

samtools view -b $1 chr1|bedtools intersect -a stdin -b $2 -v >$3.not.bam
samtools view -b $1 chr1|bedtools intersect -a stdin -b $2 |samtools view - >$3.yes.sam
./hide_marker.pl $3.yes.sam mask.bed >$3 2>$3.anon.err

#Convert and merge $3 with $3.not.bam
