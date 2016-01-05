#!/usr/bin/perl

use strict;
use warnings;
open FILE,$ARGV[0] or die;
open BED,$ARGV[1] or die;

my $DEBUG=1;
my $Clip_St=0;
my $Clip_Ed=0;
my $Clip_Ed_Real=0;
my $Motif="";
my $Chromosome;
my $Next_Gap;
my $Motif_Length=0;
my %STR_Sites=();
my $Pref_Chop= -1;
my $Motif_Units_Chopped=-1;

while (my $L=<BED>)
{
	chomp $L;
	($Chromosome,$Clip_St,$Clip_Ed,$Motif,$Next_Gap)=split(/\t/,$L);
	print STDERR join("\t",$Chromosome,$Clip_St,$Clip_Ed,$Motif,$Next_Gap);print STDERR "\n"; 
	$Motif_Length=length $Motif;
	push (@{$STR_Sites{$Chromosome}},join("\t",$Clip_St,$Clip_Ed,$Motif_Length,$Next_Gap));
}

$Chromosome="";
my @Current_STR_List=();
my $Skip_Processing=0;
my $ID=0;
my $Copies=0;
my $Max_Chop=0;
my $Max_Trim=0;
my $Chop=0;
my $STR_Switched=0;
my $Read_Count=0;
my $Extra_Tag="";
my $DB_Tag="";

Loop:while (my $L=<FILE>)
{
	my $RSkip=0;
	my ($Des,$Flag,$Chr,$Loc,$MapQ,$Cig,$Chr2,$Loc2,$Ins,$Seq,$Qual,$Rest)=split(/\t/,$L,12);

	my $RG="";
	if($Rest=~ m/RG:Z:([0-9a-zA-Z;_.]*)/g)
	{
		$RG=$1;
	}

	$Read_Count++;
	print STDERR "READ_COUNT $Read_Count\n" if(($Read_Count%10000)==0);

	if(!$Chromosome)
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
			($Clip_St,$Clip_Ed,$Motif_Length,$Next_Gap)=split(/\t/,shift @Current_STR_List);
			print STDERR "Processing $Clip_St\t$Clip_Ed\n";
			$Skip_Processing=0;
			$STR_Switched=1;
		}
		else
		{
			$Skip_Processing=1;#No more STR in the chromosome..
		}
	}
	elsif($Chromosome ne $Chr)
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
			($Clip_St,$Clip_Ed,$Motif_Length,$Next_Gap)=split(/\t/,shift @Current_STR_List);
			print STDERR "Processing $Clip_St\t$Clip_Ed\n";
			$Skip_Processing=0;
			$STR_Switched=1;
		}
		else
		{
			$Skip_Processing=1;
		}
	}
#print STDERR "LOC $Loc   CLIP_ED $Clip_Ed\n";
	Bed_Loop:while(($Loc >$Clip_Ed) && (!$Skip_Processing))#SAM has gone past STR..
	{
		print STDERR "L Processed $Clip_St\t$Clip_Ed\t$Loc\n";
		if(@Current_STR_List)
		{
			($Clip_St,$Clip_Ed,$Motif_Length,$Next_Gap)=split(/\t/,shift @Current_STR_List);
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
	if($STR_Switched)
	{
		$Clip_Ed_Real=$Clip_Ed;
		$Copies=($Clip_Ed-$Clip_St/$Motif_Length);
		$Max_Chop=4;$Max_Chop=$Copies if($Max_Chop>$Copies);#Amount of copies to chop from read..
		$Motif_Units_Chopped=1+int(rand($Max_Chop-1));
		$Chop=$Motif_Units_Chopped*$Motif_Length;#Amount removed from STR
		$Max_Trim=($Max_Chop*$Motif_Length)-$Chop;#Amount removed from the rest of the read..
		$Pref_Chop=int(rand($Max_Trim));#amount of read to be trimmed..

		$Clip_Ed-=$Chop;#Make fake boundary of clip ed to be start..
		$STR_Switched=0;
	}
	if($Skip_Processing)
	{
		print join("\t",&Generate_ID($Des),&Mask_Flag($Flag),$Chr,$Loc,$MapQ,$Cig,"*",0,0,$Seq,$Qual,"RG:Z:$RG");print "\n";
#print "$L\n";
		next Loop;
	}

	print STDERR "SKIP PROCESSING $Skip_Processing\n";
	my $Loc_Skip=$Loc;my $TLoc=$Loc;my $HLoc=$Loc;
	print STDERR "M* TLoc $TLoc \n";
	my $HCig="";
	my $Deleted_Bases=0;my $SC_Bases=0;my $Cig_Pos=0;
	my $Pref_Clip=0;
	my $Pref_Adjust=0;#How much to adjust head location to account for indels..
	my $Suff_Clip=0;
	my $Cigar_Parsed_Length=0;
	my $Suff_Chop=0;
	my $STR_Crossed=0;#Intersecting read or no
	my $Left_Border_Skipped=1;
	$Left_Border_Skipped=0 if ($Loc<$Clip_St);
	print STDERR "REST $Rest $Des $Cig\n";
	die "Error in var : $Pref_Chop"  if($Pref_Chop<0);
Loop:	while($Cig=~/([0-9]*)([MDIS])/g)
	{
		$Cigar_Parsed_Length+= length($1)+length($2);
		if($2 eq "M")
		{
			$RSkip+=$1;
			$Loc_Skip+=$1;
			$TLoc+=$1;
			print STDERR "M& TLoc $TLoc $1\n";
		}
		elsif($2 eq "S")
		{
			if($Cig_Pos==0)
			{
				$Pref_Clip=$1;
				$Pref_Chop=$Pref_Clip if($Pref_Clip>$Pref_Chop);
				$Pref_Adjust-=$1;
			}
			else
			{
				$Suff_Clip=$1;
			}
			$RSkip+=$1;
		}
		elsif($2 eq "I")
		{
			if($RSkip<$Pref_Chop)#Still start location is not fixed..
			{
				if($1<=($Pref_Chop-$RSkip))
				{
					$Pref_Adjust-=$1;
				}
				else
				{
					$Pref_Adjust-=($Pref_Chop-$RSkip);
				}

			}
			$RSkip+=$1;
		}
		elsif($2 eq "D")
		{
			$Loc_Skip+=$1;
			$TLoc+=$1;
		}
		print STDERR "LOC_SKIP $Loc_Skip\n";

		if(($Loc_Skip>=$Clip_St) && ($2 ne "S") && (!$Left_Border_Skipped))#Passed the STR?
		{
			$STR_Crossed=1;	

			my $Cig_Adjust=0;#amount of match in Cigar ..
			my $Gap=$Loc_Skip-$Clip_St;$Gap=0 if($2 eq 'D');
			my $Leftover_Bases=0;
			if(($2 eq "M")||($2 eq "I"))
			{
				$Cig_Adjust=$1-$Gap; 
#die "$Cig_Adjust $Gap $1 $Cig" if ($Cig_Adjust<0);
				next Loop if ($Cig_Adjust<0);
			}
			my $Head=substr($Seq,0,$RSkip-$Gap);my $HeadQ=substr($Qual,0,$RSkip-$Gap);#Chop before STR/After STR..
			my $Tail=substr($Seq,$RSkip-$Gap);my $TailQ=substr($Qual,$RSkip-$Gap);
			my ($Des1,$Flag1,$Chr1,$Loc1,$MapQ1,$Chr21,$Loc21,$Ins1,$Rest1)=($Des,$Flag,$Chr,$Loc_Skip+$Chop,$MapQ,$Chr2,$Loc2,$Ins,$Rest);
			
			$Tail=substr($Tail,$Chop);$TailQ=substr($TailQ,$Chop);
#Mask readlength..
			if($Pref_Chop>length($Head))
			{
				$Head="";$Pref_Chop=0;
			}
			if($Pref_Chop<$Max_Trim)
			{
				$Suff_Chop=int ($Max_Trim-$Pref_Chop);
			}
			$DB_Tag="CT:Z:${Chop}_(${Pref_Chop}_${Suff_Chop})_${Motif_Length}_${Motif_Units_Chopped}_${Max_Trim}";
#Make cigars..
			$Head=substr($Head,$Pref_Chop);$HLoc+=($Pref_Chop+$Pref_Adjust);
			$HeadQ=substr($HeadQ,$Pref_Chop);
			if($Head)
			{
				my $Head_Cig=substr($Cig,0,$Cigar_Parsed_Length);
				my $Last_Match=0;my $Cig_Left=$Cigar_Parsed_Length;
				while($Head_Cig=~/([0-9]*)([MDIS])/g)
				{
					$Cig_Left-=((length $1)+(length $2));
					if(($2 eq 'M')||($2 eq 'I')||($2 eq 'S'))
					{
						my $Count=$1;
						print STDERR "HH $Count $2\n";
						if(!$Cig_Left)
						{
							$Count=$Cig_Adjust;
							print STDERR "HH -- $Count\n";
						}
						if($Count<$Pref_Chop)
						{
							$Pref_Chop-=$Count;
							print STDERR "HH ++ $Count\n";
						}
						else
						{
							print STDERR "HH .. PChop $Pref_Chop $HCig\n";
							$Pref_Chop=$Count-$Pref_Chop;
							$HCig.=$Pref_Chop.$2 if ($Pref_Chop);
							print STDERR "HH .. PChop $Pref_Chop $HCig\n";
							$Pref_Chop=0;
						}
					}
					else
					{
						$HCig.=($1.$2);
						print STDERR "HH ** $HCig\n";
					}

				}
			}

			my $TCig="";
			$TLoc-=$Gap;print STDERR "$TLoc $Gap\n";
			($TCig,my $SC,my $Skip,my $Last_Pos)=&Parse_TCig(substr($Cig,$Cigar_Parsed_Length),$Gap,$Chop,$Suff_Chop);#SC=soft clip, Skip=amount to advance tail..
			print STDERR "CHOP $Chop SKIP $Skip LAST_POS $Last_Pos\n";
			$Last_Pos+=$TLoc;
			my $Extend=$Last_Pos-$Clip_Ed;
			print STDERR "LAST_POS $Last_Pos EXTEND $Extend NEXT GAP $Next_Gap\n";
			if($Extend>$Next_Gap)
			{
				$Head="";$Tail="";
			}
			$Tail="" if(!$TCig);
			$Skip=$Chop+$Skip;   #Account for bases in reference in the chopped region..
			print STDERR "T& Chop $Chop Skip $Skip TLoc $TLoc\n";
			$TLoc+=$Skip; 
			print STDERR "SufCh $Suff_Chop\n";
			if($SC>=$Suff_Chop)
			{
				$Suff_Chop=$SC; 
			}

			print STDERR "SufCh $Suff_Chop\n";
			if($Tail)
			{
				$Tail=substr($Tail,0,length($Tail)-$Suff_Chop);
				$TailQ=substr($TailQ,0,length($TailQ)-$Suff_Chop);
			}
			print STDERR "$Tail\tC:$Cig\tT:$TCig\tSC:$SC TL:$TLoc\n";
			print STDERR "H $Head T $Tail\n$Seq $Cig LOC $Loc\n";
			print STDERR "$Seq CH $Chop HCIG:$HCig HLOC $Loc1 SC:$Suff_Chop PC:$Pref_Chop\n";

#		print $L;
			$Extra_Tag="CT:Z:${Chop}_(${Pref_Chop}_${Suff_Chop})_${Motif_Length}_${Max_Trim}";
			if ($Head) 
			{
				print join("\t",&Generate_ID($Des),&Mask_Flag($Flag),$Chr,$HLoc,$MapQ,$HCig,"*",0,0,$Head,$HeadQ,"RG:Z:$RG",$Extra_Tag,$DB_Tag);print "\n";
			}
			if ($Tail) 
			{
				print join("\t",&Generate_ID($Des),&Mask_Flag($Flag),$Chr,$TLoc,$MapQ,$TCig,"*",0,0,$Tail,$TailQ,"RG:Z:$RG",$Extra_Tag,$DB_Tag);print "\n";
			}
			if((!$Head)&&(!$Tail))
			{
				print join("\t",&Generate_ID($Des),4,$Chr,$Loc,$MapQ,$Cig,"*",0,0,"*","*","RG:Z:$RG",$Extra_Tag,$DB_Tag);print "\n";
			}
			last Loop;
		}
		elsif(($Loc_Skip>=$Clip_Ed) && ($2 ne "S"))#Passed the STR?
		{
			print STDERR "Clip-CHop $Clip_Ed - $Chop\n";
			$STR_Crossed=1;	

			my $Cig_Adjust=0;#amount of match in Cigar ..
			my $Gap=$Loc_Skip-$Clip_Ed;$Gap=0 if($2 eq 'D');
			my $Leftover_Bases=0;
			if(($2 eq "M")||($2 eq "I"))
			{
				$Cig_Adjust=$1-$Gap; 
#die "$Cig_Adjust $Gap $1 $Cig" if ($Cig_Adjust<0);
				next Loop if ($Cig_Adjust<0);
			}
			my $Head=substr($Seq,0,$RSkip-$Gap);my $HeadQ=substr($Qual,0,$RSkip-$Gap);#Chop before STR/After STR..
			my $Tail=substr($Seq,$RSkip-$Gap);my $TailQ=substr($Qual,$RSkip-$Gap);
			my ($Des1,$Flag1,$Chr1,$Loc1,$MapQ1,$Chr21,$Loc21,$Ins1,$Rest1)=($Des,$Flag,$Chr,$Loc_Skip+$Chop,$MapQ,$Chr2,$Loc2,$Ins,$Rest);

			$Tail=substr($Tail,$Chop);$TailQ=substr($TailQ,$Chop);
#Mask readlength..
			if($Pref_Chop>length($Head))
			{
				$Head="";$Pref_Chop=0;
			}
			if($Pref_Chop<$Max_Trim)
			{
				$Suff_Chop=int ($Max_Trim-$Pref_Chop);
			}
			$DB_Tag="CT:Z:${Chop}_${Pref_Chop}_${Suff_Chop}_${Motif_Length}_${Max_Trim}";
#Make cigars..
			$Head=substr($Head,$Pref_Chop);$HLoc+=($Pref_Chop+$Pref_Adjust);
			$HeadQ=substr($HeadQ,$Pref_Chop);
			if($Head)
			{
				my $Head_Cig=substr($Cig,0,$Cigar_Parsed_Length);
				my $Last_Match=0;my $Cig_Left=$Cigar_Parsed_Length;
				print STDERR "HEAD CIG* $Head_Cig\n";
				while($Head_Cig=~/([0-9]*)([MDIS])/g)
				{
					$Cig_Left-=((length $1)+(length $2));
					if(($2 eq 'M')||($2 eq 'I')||($2 eq 'S'))
					{
						my $Count=$1;
						print STDERR "HH* $Count $2\n";
						if(!$Cig_Left)
						{
							$Count=$Cig_Adjust;
							print STDERR "HH* -- $Count\n";
						}
						if($Count<$Pref_Chop)
						{
							$Pref_Chop-=$Count;
							print STDERR "HH* ++ $Count\n";
						}
						else
						{
							print STDERR "HH* .. PChop $Pref_Chop $HCig\n";
							$Pref_Chop=$Count-$Pref_Chop;
							$HCig.=$Pref_Chop.$2 if ($Pref_Chop);
							print STDERR "HH* .. PChop $Pref_Chop $HCig\n";
							$Pref_Chop=0;
						}
					}
					else
					{
						$HCig.=($1.$2);
						print STDERR "HH* ** $HCig\n";
					}

				}
			}

			my $TCig="";
			$TLoc-=$Gap;print STDERR "$TLoc $Gap\n";
			($TCig,my $SC,my $Skip,my $Last_Pos)=&Parse_TCig(substr($Cig,$Cigar_Parsed_Length),$Gap,$Chop,$Suff_Chop);#SC=soft clip, Skip=amount to advance tail..
			print STDERR "*CHOP $Chop SKIP $Skip LAST_POS $Last_Pos\n";
			$Last_Pos+=$TLoc;
			my $Extend=$Last_Pos-$Clip_Ed;
			print STDERR "*LAST_POS $Last_Pos EXTEND $Extend NEXT GAP $Next_Gap\n";
			if($Extend>$Next_Gap)
			{
				$Head="";$Tail="";
			}
			$Tail="" if(!$TCig);
			$Skip=$Chop+$Skip;   #Account for bases in reference in the chopped region..
			print STDERR "*T& Chop $Chop Skip $Skip TLoc $TLoc\n";
			$TLoc+=$Skip; 
			print STDERR "*SufCh $Suff_Chop\n";
			if($SC>=$Suff_Chop)
			{
				$Suff_Chop=$SC; 
			}

			print STDERR "*SufCh $Suff_Chop\n";
			if($Tail)
			{
				$Tail=substr($Tail,0,length($Tail)-$Suff_Chop);
				$TailQ=substr($TailQ,0,length($TailQ)-$Suff_Chop);
			}
			print STDERR "*$Tail\tC:$Cig\tT:$TCig\tSC:$SC TL:$TLoc\n";
			print STDERR "*H $Head T $Tail\n$Seq $Cig LOC $Loc\n";
			print STDERR "*$Seq CH $Chop HCIG:$HCig HLOC $Loc1 SC:$Suff_Chop PC:$Pref_Chop\n";

#	print $L;
			$Extra_Tag="CT:Z:${Chop}_(${Pref_Chop}_${Suff_Chop})_${Motif_Length}_${Max_Trim}";
			if ($Head) 
			{
				print join("\t",&Generate_ID($Des),&Mask_Flag($Flag),$Chr,$HLoc,$MapQ,$HCig,"*",0,0,$Head,$HeadQ,"RG:Z:$RG",$Extra_Tag,$DB_Tag);print "\n";
			}
			if ($Tail) 
			{
				print join("\t",&Generate_ID($Des),&Mask_Flag($Flag),$Chr,$TLoc,$MapQ,$TCig,"*",0,0,$Tail,$TailQ,"RG:Z:$RG",$Extra_Tag,$DB_Tag);print "\n";
			}
			if((!$Head)&&(!$Tail))
			{
				print join("\t",&Generate_ID($Des),4,$Chr,$Loc,$MapQ,$Cig,"*",0,0,"*","*","RG:Z:$RG",$Extra_Tag,$DB_Tag);print "\n";
			}
			last Loop;
		}
	}
	if(!$STR_Crossed)
	{
		print join("\t",&Generate_ID($Des),&Mask_Flag($Flag),$Chr,$Loc,$MapQ,$Cig,"*",0,0,$Seq,$Qual,"RG:Z:$RG");print "\n";
	}

	print STDERR "\n";


}

sub Parse_TCig
{
	my $Cig=shift @_;
	my $Left_Over=shift @_;
	my $Chop=shift @_;
	my $Suff_Chop=shift @_;
	my $Skip=0;
	my $Ref_Length_Travelled=$Left_Over;#Amount skipped along the reference from here..
	print STDERR join ("\t",$Cig,$Left_Over,$Chop,$Suff_Chop);print STDERR "\n";

	my $Cigar_Length=length $Cig;
	if($Cigar_Length==0)
	{
		print STDERR "****** $Left_Over $Chop $Suff_Chop\n";
		if($Left_Over>$Chop+$Suff_Chop)
		{
			$Left_Over-=($Chop+$Suff_Chop);
			return ($Left_Over."M",0,$Skip,$Ref_Length_Travelled);
		}
		else
		{
			return ("",0,0,$Ref_Length_Travelled);
		}


	}
	my $Final_Cigar="";
	my $Last_Match=0;
	my $No_Parse=1;
	my $Cig_Parsed=0;


	while($Cig=~/([0-9]*)([MDSI])/g)
	{
		$Cigar_Length-=(length($1)-length($2));
		$Cig_Parsed=1 if(($2 eq 'S') || ($Cigar_Length==0));
		print STDERR "$1$2\n";
		$No_Parse=0;
		if($2 eq "M")
		{
			$Ref_Length_Travelled+=$1;
			if($Chop)
			{
				if($Left_Over>$Chop)
				{
					print STDERR "+=+= $Left_Over $Chop\n";
					$Last_Match=$1+$Left_Over-$Chop;
					$Chop=0;
				}
				else
				{
					if($Left_Over)#If leftover chop from it..
					{
						print STDERR "+++++++ $Left_Over $Chop\n";
						$Chop=$Chop-$Left_Over;
					}
					else#Otherwise try to chop from match..
					{
						if($1>$Chop)
						{
							print STDERR "----- $Left_Over $Chop\n";
							$Last_Match=$1-$Chop;
							$Chop=0;
						}
						else
						{
							print STDERR "LO $Last_Match CH $Chop C $1\n";
							$Chop-=$1;
							$Last_Match=0;
						}
					}
				}
			}
			else
			{
				print STDERR "@@@@ $Last_Match $1\n";
				$Last_Match=$1;
			}
			$Left_Over=0;
		}
		else
		{
			if ($Left_Over)
			{
				if($Chop)
				{
					if($Chop>=$Left_Over)
					{
						print STDERR "XX CHOP $Chop SKIP $Skip LEFT $Left_Over\n";
						$Chop-=$Left_Over;
						$Left_Over=0;
					}
					else
					{
						print STDERR "XX& CHOP $Chop SKIP $Skip LEFT $Left_Over\n";
						$Left_Over-=$Chop;
						$Chop=0;
					}
				}
			}
			if($2 eq 'I')
			{
				if($Chop)
				{
					if(!$Left_Over)
					{
						if($Chop>$1)#How many location skips to add..
						{
							$Skip-=$1;
							print STDERR "SI $Chop $Skip $1\n";
						}
						else
						{
							$Skip-=$Chop;
							print STDERR "SI& $Chop $Skip $1\n";
						}
					}
				}
			}
			if($2 eq 'D')
			{
				$Ref_Length_Travelled+=$1;
				if($Chop)
				{
					$Skip+=$1;
					print STDERR "SD& $Chop $Skip $1\n";
				}
			}
			if($2 eq 'S')
			{
				if($Chop>0)
				{
					return ("",0,$Skip,$Ref_Length_Travelled);#chopping, but only SC left..
				}
				print STDERR ">>>> $1$2 SC $Suff_Chop\t$Last_Match\n";
				if($1>=$Suff_Chop)
				{
					$Suff_Chop=0; 
				}
				else
				{
					$Suff_Chop-=$1; 
				}
				print STDERR ">>>> $1$2 SC $Suff_Chop\n";die if ($Suff_Chop<0);
			}

			if($Left_Over)
			{
				if($Chop)
				{
					if($Left_Over>$Chop)
					{
						print STDERR "^^^ $Left_Over $Chop\n";
						$Last_Match=$Left_Over-$Chop;
						$Chop=0;
					}
					else
					{
						if($Left_Over)
						{
							print STDERR "^!^ $Left_Over $Chop\n";
							$Last_Match=$Left_Over;
						}
						$Chop=$Chop-$Left_Over;
					}
				}
				else
				{
					print STDERR ":::::$Left_Over\n";
					$Last_Match=$Left_Over;
				}
				$Left_Over=0;
			}

			if($Last_Match<$Suff_Chop)
			{
				if($Cigar_Length)
				{
					my $Remaining_Length=substr($Cig,length($Cig)-$Cigar_Length);
					$Ref_Length_Travelled+=&Find_Skip($Remaining_Length);
					print STDERR "CIG_LEN $Remaining_Length $1$2 \n";

				}
				return ("",0,$Skip,$Ref_Length_Travelled);
			}
			$Last_Match-=$Suff_Chop if($Last_Match && $Cig_Parsed);
			$Final_Cigar.=$Last_Match."M" if($Last_Match);
			if($2 eq "S")
			{
				return ($Final_Cigar,$1,$Skip,$Ref_Length_Travelled);
			}
			$Final_Cigar.="$1$2" if (($2 ne "D") || ($Cigar_Length));
			$Last_Match=0;
		}

	}


	$Last_Match-=$Suff_Chop if($Last_Match);
#die if($Chop);
#die "LM $Last_Match CH $Chop" if($Last_Match<=0);
	$Final_Cigar="" if($Chop);
	$Final_Cigar="" if($Last_Match<=0);
	$Final_Cigar.=$Last_Match."M" if($Last_Match>0);
	return ($Final_Cigar,0,$Skip,$Ref_Length_Travelled);
}

sub Generate_ID
{
	my $Des=shift @_;
	if($DEBUG)
	{
		return $Des;
	}
	$ID++;
	return "$ID";

}

sub Mask_Flag
{
	my $Flag=shift @_;
	if($Flag & 16)
	{
		return 16;
	}
	else
	{
		return 0;
	}
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
