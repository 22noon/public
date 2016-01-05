#!/usr/bin/perl
#0-based range masked..

use strict;
use warnings;
use Switch;

open FILE,$ARGV[0] or die;
open BED,$ARGV[1] or die;

my $DEBUG=1;
my $Clip_St=0;
my $Clip_Ed=0;
my $Chromosome;
my %STR_Sites=();

while (my $L=<BED>)
{
	chomp $L;
	($Chromosome,$Clip_St,$Clip_Ed)=split(/\t/,$L);
	print STDERR join("\t",$Chromosome,$Clip_St,$Clip_Ed) if($DEBUG);print STDERR "\n" if($DEBUG); 
	push (@{$STR_Sites{$Chromosome}},join("\t",$Clip_St,$Clip_Ed));
}

$Chromosome="";
my @Current_STR_List=();
my $Skip_Processing=0;
my $ID=0;
my $Extra_Tag="";
my $DB_Tag="";
my $Read_Count=0;
my $STR_Switched=0;

Loop:while (my $L=<FILE>)
{
	chomp $L;
	my $RSkip=0;
	my ($Des,$Flag,$Chr,$Loc,$MapQ,$Cig,$Chr2,$Loc2,$Ins,$Seq,$Qual,$Rest)=split(/\t/,$L,12);

	my $RG="";
	if($Rest=~ m/RG:Z:([0-9a-zA-Z;_.]*)/g)
	{
		$RG=$1;
	}

	$Read_Count++;
	print STDERR "READ_COUNT $Read_Count\n" if(($Read_Count%10000)==0);

	if(!$Chromosome)#Init
	{
		if(exists $STR_Sites{$Chr})
		{
			@Current_STR_List=@{$STR_Sites{$Chr}};
			$Chromosome=$Chr;
		}
		else
		{
			@Current_STR_List=();
		}
		if(@Current_STR_List)
		{
			($Clip_St,$Clip_Ed)=split(/\t/,shift @Current_STR_List);
			print STDERR "Processing $Clip_St\t$Clip_Ed\n";
			$Skip_Processing=0;
			$STR_Switched=1;
		}
		else
		{
			$Skip_Processing=1;#No more STR in the chromosome..
		}
	}
	elsif($Chromosome ne $Chr)#New chromosome..
	{
		if(exists $STR_Sites{$Chr})
		{
			@Current_STR_List=@{$STR_Sites{$Chr}}
		}
		else
		{
			@Current_STR_List=();
		}
		if(@Current_STR_List)
		{
			($Clip_St,$Clip_Ed)=split(/\t/,shift @Current_STR_List);
			print STDERR "Processing $Clip_St\t$Clip_Ed\n";
			$Skip_Processing=0;
			$STR_Switched=1;
		}
		else
		{
			$Skip_Processing=1;
		}
	}
	Bed_Loop:while(($Loc >$Clip_Ed) && (!$Skip_Processing))#SAM has gone past STR..
	{
		print STDERR "L Processed $Clip_St\t$Clip_Ed\t$Loc\n";
		if(@Current_STR_List)
		{
			($Clip_St,$Clip_Ed)=split(/\t/,shift @Current_STR_List);
			print STDERR "Processing $Clip_St\t$Clip_Ed\n";
			$Skip_Processing=0;
			$STR_Switched=1;
		}
		else
		{
			$Skip_Processing=1;
			last Bed_Loop;
		}
	}

	if($Skip_Processing)
	{
		print "$L\n"; 
		next Loop;
	}

	my $Loc_Skip=$Loc;
	my $Cigar_Pos=0;
	my $Revealed_Read="";
	print STDERR "$Des\n" if ($DEBUG);
Loop:	while($Cig=~/([0-9]*)([MDIS])/g)
	{
		switch($2)
		{
			case "M"
			{
				my $Old_Loc=$Loc_Skip;
				$Loc_Skip+=$1;
				if($Loc_Skip>=$Clip_St && $Old_Loc<$Clip_Ed)
				{
					if($Old_Loc<$Clip_St)#Match starts before bdry..
					{
						my $Over=$Loc_Skip-$Clip_St;
						my $Gap=$1-$Over;die"$Loc_Skip $Clip_St\n" if($Gap<0);
						print STDERR "M $Over $Loc_Skip Gap $Gap O $Old_Loc $Revealed_Read\n" if($DEBUG);
						$Revealed_Read.=substr($Seq,$RSkip,$Gap);
						print STDERR "M $Over $Loc_Skip $Old_Loc $Revealed_Read\n" if($DEBUG);

						my $Left_Over=0;
						if($Loc_Skip>$Clip_Ed)
						{
							$Left_Over=$Over-($Clip_Ed-$Clip_St);
							$Over=$Clip_Ed-$Clip_St;
						}
						print STDERR "M $Left_Over $Loc_Skip Over $Over $Revealed_Read\n" if($DEBUG);
						$Revealed_Read.="N"x$Over;
						print STDERR "M $Over $Loc_Skip $Old_Loc $Revealed_Read\n" if($DEBUG);
						if($Left_Over)
						{
							$Revealed_Read.=substr($Seq,$RSkip+$Gap+$Over,$Left_Over);
						}
						print STDERR "M $Over $Loc_Skip $Old_Loc $Revealed_Read\n" if($DEBUG);

					}
					else
					{
						if($Loc_Skip<=$Clip_Ed)
						{
							print STDERR ">M $Loc_Skip $Old_Loc $Revealed_Read\n" if($DEBUG);
							$Revealed_Read.="N"x($Loc_Skip-$Old_Loc);
							print STDERR ">M $Loc_Skip $Old_Loc $Revealed_Read\n" if($DEBUG);
						}	
						else
						{
							my $Over=$Clip_Ed-$Old_Loc;
							print STDERR ">>M $Over $Loc_Skip $Old_Loc $Revealed_Read\n" if($DEBUG);
							$Revealed_Read.="N"x$Over;
							print STDERR ">>M $Over $Loc_Skip $Old_Loc $Revealed_Read\n" if($DEBUG);
							die "$1 $Over" if ($1<$Over);
							$Revealed_Read.=substr($Seq,$RSkip+$Over,$1-$Over);
							print STDERR ">>M $Over $Loc_Skip $Old_Loc $Revealed_Read\n" if($DEBUG);
						}
					}

				}
				else
				{
					print STDERR "*M $RSkip $Revealed_Read\n" if($DEBUG);
					$Revealed_Read.=substr($Seq,$RSkip,$1);
					print STDERR "*M $Revealed_Read\n" if($DEBUG);
				}
				$RSkip+=$1;
			}
			case "S"
			{
				$Revealed_Read.=substr($Seq,$RSkip,$1);
				$RSkip+=$1;
			}
			case "I"
			{
				if($Loc_Skip>=$Clip_St && $Loc_Skip<$Clip_Ed)#Ins is in the region..
				{
					print STDERR ">I $Revealed_Read $RSkip $1\n" if($DEBUG);
					$Revealed_Read.="N"x$1;
					print STDERR ">I $Revealed_Read $RSkip $1\n" if($DEBUG);
				}
				else
				{
					print STDERR "I $Revealed_Read $RSkip $1\n" if($DEBUG);
					$Revealed_Read.=substr($Seq,$RSkip,$1);
					print STDERR "I $Revealed_Read $RSkip $1\n" if($DEBUG);
				}
				$RSkip+=$1;
			}
			case "D"
			{
				$Loc_Skip+=$1;
			}
		}

#		print STDERR "LOC_SKIP $Loc_Skip\n";
		$Cigar_Pos+= length($1)+length($2);

		if(($Loc_Skip>=$Clip_St) && ($2 ne "S"))#Passed the STR?
		{
#			print ">>$L"; 
		}
	}
	print join("\t",$Des,$Flag,$Chr,$Loc,$MapQ,$Cig,$Chr2,$Loc2,$Ins,$Revealed_Read,$Qual,$Rest);print "\n";

}


sub Find_Skip
{
	my $Cig=shift @_;
	my $Skip=0;
	while($Cig=~/([0-9]*)([MDSI])/g)
	{
		if(($2 eq "M")||($2 eq 'D'))
		{
			$Skip+=$1;
		}
	}
	return $Skip;
}
