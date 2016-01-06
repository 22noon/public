This contains several prototype programs used to test the algorithms in the submission.
1) Program to hide the Quasi-Identifiers
hide_marker.sh <BAM file> <BED file> <Output file>
This program takes a BAM file as the input and a set of regions in the BED format and removes sequencing information from reads crossing the regions specified in the BED file. The output is a SAM file containig masked reads and a BAM file containing reads that do not cross the regions of interest. These two files can be merged for further downstream analysis.

2) Program to hide STR information
unstr.sh <BAM file> <Expanded bed file> <Prefix> <STR regions>
This program requires STR regions to be hidden in the BED format and a BAM file. The <Expanded bed file> is created by running the expand_bedfile.sh command on the STR regions, and requires the calculation of the insert size of the BAM file. The insert size can be estimated using a large number of highquality paired reads and example commands are given inside the expand_bedfile.sh file. 
After running the program two files are created, Prefix.unstr.bam and Prefix.not.bam. The first contains the anonymized version of STR regions, and the second those reads falling outside it. These two can be joined together for further downstream analysis. 
