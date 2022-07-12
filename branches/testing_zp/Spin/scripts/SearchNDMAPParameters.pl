#!/usr/bin/perl -w

use strict 'vars';

my $ProgramCommand = "./src/Programs/Spin/NDMAPSpinChain";
my $ScalingFactor = 2.8549;
my $NbrValue = 3;
my $NbrSpin = 8;
my $Momentum = $NbrSpin / 2;

if (!defined($ARGV[5]))
  {
    die ("usage: SearchNDMAPParameters min_d max_d sub_d min_e max_e sub_e\n")
  }
my $MinD = $ARGV[0];
my $MaxD = $ARGV[1];
my $IncD = ($MaxD - $MinD) / $ARGV[2];
if ($IncD == 0)
  {
    $IncD = 1.0;
  }
my $MinE = $ARGV[3];
my $MaxE = $ARGV[4];
my $IncE = ($MaxE - $MinE) / $ARGV[5];
if ($IncE == 0)
  {
    $IncE = 1.0;
  }
my $CurrentD = $MinD;
my $CurrentE;
my $Ground;
my $Excited1;
my $Excited2;
my $Excited3;

while ($CurrentD <= $MaxD)
  {
    $CurrentE = $MinE;
    while ($CurrentE <= $MaxE)
      {
	my $OutputLog = `$ProgramCommand -d $CurrentD -e $CurrentE -z 0 -a 0 --fixed-momentum --momentum 0 $NbrSpin -n 1`;
	print "$ProgramCommand -d $CurrentD -e $CurrentE -z 0 -a 0 --fixed-momentum --momentum 0 $NbrSpin -n 1\n";
	print $OutputLog;
	my @Lines = split (/\n/, $OutputLog);
	my $Flag = 0;
	my $TmpLine;
	foreach $TmpLine (@Lines)
	  {
	    if ($Flag == 1)
	      {
		chomp ($TmpLine);
		my @Values = split (/ /, $TmpLine);
		$Ground = $Values[0];
		last;
	      }
	    if ($TmpLine =~ /Nbr of iterations \=/)
	      {
		$Flag = 1;
	      }
	  }
	$OutputLog = `$ProgramCommand -d $CurrentD -e $CurrentE -z 0 -a 0 --fixed-momentum --momentum $Momentum $NbrSpin -n 3`;
	print "$ProgramCommand -d $CurrentD -e $CurrentE -z 0 -a 0 --fixed-momentum --momentum $Momentum $NbrSpin -n 3\n";
	print $OutputLog;
	@Lines = split (/\n/, $OutputLog);
	$Flag = 0;
	foreach $TmpLine (@Lines)
	  {
	    chomp ($TmpLine);
	    if ($Flag == 1)
	      {
		my @Values = split (/ /, $TmpLine);
		$Excited1 = $Values[0];
		$Excited2 = $Values[1];
		$Excited3 = $Values[2];
		last;
	      }
	    if ($TmpLine =~ /Nbr of iterations \=/)
	      {
		$Flag = 1;
	      }
	  }
	print ($Ground." ".$Excited1." ".$Excited2." ".$Excited3."\n");
	my $GapScalingFactor = 0.42 / ($Excited1 - $Ground);
	open (OUTFILE, ">>ndmapgap.dat");
	print OUTFILE ($CurrentD." ".$CurrentE." 0.42 ".(($Excited2 - $Ground) * $GapScalingFactor)." ".(($Excited3 - $Ground) * $GapScalingFactor)."\n");
	print ("0.42 ".(($Excited2 - $Ground) * $GapScalingFactor)." ".(($Excited3 - $Ground) * $GapScalingFactor)."\n");
	close (OUTFILE);
	print "---------------------------------------------------------------------------\n";
	$CurrentE += $IncE;
      }
    $CurrentD += $IncD;    
  }

