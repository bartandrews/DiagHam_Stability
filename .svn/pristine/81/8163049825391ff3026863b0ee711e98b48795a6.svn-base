#!/usr/bin/perl -w

use strict 'vars';

if (!(defined($ARGV[4])))
  {
    die "usage: FindGapGeneric StartN StartS NInc SInc Caption [PrintFlag]";
  }
my $PrintFlag = 0;
my $NbrBosons = $ARGV[0];
my $S = $ARGV[1];
my $NbrBosonsInc = $ARGV[2];
my $SInc = $ARGV[3];
my $Caption = $ARGV[4];
if (defined($ARGV[5]))
  {
    $PrintFlag = 1;
  }
my %MinArray;
my $TmpFile;
while ($NbrBosons <= 40)
  {
    $TmpFile = "n_".$NbrBosons."/bosons_delta_n_".$NbrBosons."_2s_".$S."_l.dat";
    if (-e $TmpFile)
      {
	print ($TmpFile."\n");
	$MinArray{$NbrBosons} = &FindGap($TmpFile);
      }
    $NbrBosons += $NbrBosonsInc;
    $S += $SInc;
  }
&CreatePostScript(\%MinArray, $Caption, $PrintFlag);

# find gap in a file
#
# $_[0] = file name
# return value = ground state energy

sub FindGap
  {
    my $FileName = $_[0];
    my $Min;
    my $Min2;
    my $Flag = 0;
    open (INFILE, $FileName);
    my $TmpLine;
    foreach $TmpLine (<INFILE>)
      {
	chomp ($TmpLine);
	my @TmpArray = split (/ /, $TmpLine);
	if ($Flag == 0)
	  {
	    $Min = $TmpArray[1];
	    $Flag = 1;
	  }
	else
	  {
	    if ($Flag == 2)
	      {
		if ($TmpArray[1] < $Min)
		  {
		    $Min2 = $Min;
		    $Min = $TmpArray[1];
		  }
		else
		  {
		    if ($TmpArray[1] < $Min2)
		      {
			$Min2 = $TmpArray[1];
		      }
		  }
	      }
	    else
	      {
		if ($TmpArray[1] < $Min)
		  {
		    $Min2 = $Min;
		    $Min = $TmpArray[1];
		  }
		else
		  {
		    $Min2 = $TmpArray[1];
		  }
		$Flag = 2;
	      }
	  }
      }
    close (INFILE);
    return ($Min2 - $Min);
  }

# create postscript graph from data file
#
# $_[0] = hash table containing datas
# $_[1] = print flag (1 if true)
# $_[2] = number of bosons

sub CreatePostScript
  {
    my $Datas = $_[0];
    my $Caption = $_[1];
    my $PrintFlag = $_[2];
    my $N;
    my $E;
    my $FileName = "bosons_delta_gap_".$Caption.".dat";
    open (OUTFILE, ">$FileName");
    my $MinN = 200;
    my $MaxN = 0;
    my $MinGap = 400;
    my $MaxGap = 0;
    while (($N, $E) = each (%$Datas))
      {
	if ($MinN > $N)
	  {
	    $MinN = $N;
	  }
	if ($MaxN < $N)
	  {
	    $MaxN = $N;
	  }
	if ($MinGap > $E)
	  {
	    $MinGap = $E;
	  }
	if ($MaxGap < $E)
	  {
	    $MaxGap = $E;
	  }
	$N = 1.0 / $N;
	print ($N." ".$E."\n");
	print OUTFILE ($N." ".$E."\n");
      }
    close (OUTFILE);
    $MinGap = 0;
    my $Delta = ($MaxGap - $MinGap) / 20.0;
    $MaxGap += $Delta;
    $MinGap -= $Delta;
    $MinGap = 0;
    $MinN--;
    $MaxN++;
    my $Tmp = 1.0 / $MinN;
    $MinN = 1.0 / $MaxN;
    $MaxN = $Tmp;
    my $TmpFileName = "tmp".time().".p";
    my $OutputFile = "bosons_delta_gap_".$Caption.".ps";
    my @TmpArray = split (/_/,  $OutputFile);
    my $Title = "gap ".$Caption;
    open (OUTFILE, ">$TmpFileName");
    print OUTFILE ("set xrange [".$MinN.":".$MaxN."]
set yrange [".$MinGap.":".$MaxGap."]
set xlabel \"1/N\"
set ylabel \"E\"
set size 1.0, 0.6
set terminal postscript portrait enhanced \"Helvetica\" 14
set output \"".$OutputFile."\"
plot \"".$FileName."\" using 1:2 title \"".$Title."\"
");
    close (OUTFILE);
    `gnuplot $TmpFileName`;
    if ($PrintFlag == 1)
      {
	`lpr $OutputFile`;
      }
    `rm -f $TmpFileName`;
  }



