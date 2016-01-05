#Split bam file into STR and non STR regions..
#First parameter bam file
#second one is the expanded bed region of STR file (run expand_bedfile.sh to get this) 
#Third is the prefix of the splits..
#fourth is the STR regions to strip..

bedtools intersect -a $1 -b $2 -v >$3.not.bam
bedtools intersect -a $1 -b $2 |samtools view - >$3.yes.sam

#Remove STR information..
./unstr.pl $3.yes.sam $4 >$3.unstr.sam 2>$3.unstr.err
#Convert anonymized sam into bam..
samtools view -H $3.not.bam |grep "@SQ"|cut -f2-|sed 's/SN://'|sed 's/LN://' >$3.header.txt
samtools view -t $3.header.txt -bS $3.unstr.sam >$3.unstr.bam
#$3.unstr.bam and $3.not.bam can be merged together using for example, SAMTools now..
