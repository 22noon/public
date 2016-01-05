#!/usr/bin/perl
#expand bed file ranges to mask insert sizes..
use warnings;
use strict;

open INS,$ARGV[0] or die;
open BED,$ARGV[1] or die;

if(my $L=<INS>)
{
	chomp $L;
	my ($n,$Count,$avg,$Avg,$std,$STD)=split(/[= ]/,$L);
	$Avg=int($Avg);$STD=int($STD);
	my $XSTD=$ARGV[2]*$STD;
	while(my $Line=<BED>)
	{
		if(substr($Line,0,1) eq "#")
		{
			print $Line;
		}
		else
		{
			my ($Chr,$St,$Ed,$Rest)=split(/\t/,$Line,4);
			$St-=$XSTD;$St=1 if($St<=0);
			$Ed+=$XSTD;
			print join("\t",$Chr,$St,$Ed,$Rest);
		}
	}
}
else 
{
	die "Empty file..";
}

