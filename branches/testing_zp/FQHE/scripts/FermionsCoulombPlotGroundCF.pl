#!/usr/bin/perl -w

use strict 'vars';

if (!(defined($ARGV[4])))
  {
    die "usage: FindCoulombPlotGroundCF StartN StartS NInc SInc Caption [1_N_Factor] [Shift] [PrintFlag]";
  }
my $PrintFlag = 0;
my $NbrFermions = $ARGV[0];
my $S = $ARGV[1];
my $NbrFermionsInc = $ARGV[2];
my $SInc = $ARGV[3];
my $Caption = $ARGV[4];
my $FactorInvN = 0.0;
my $FactorShift = 0.0;
if (defined($ARGV[5]))
  {
    $FactorInvN = $ARGV[5];
  }
if (defined($ARGV[6]))
  {
    $FactorShift = $ARGV[6];
  }
if (defined($ARGV[7]))
  {
    $PrintFlag = 1;
  }
my %MinArray;
my $TmpFile;
while ($NbrFermions <= 40)
  {
    $TmpFile = "n_".$NbrFermions."/fermions_coulomb_n_".$NbrFermions."_2s_".$S."_lz.dat";
    if (-e $TmpFile)
      {
	print ($TmpFile."\n");
	my $Scaling = sqrt(($S * $NbrFermionsInc) / ($NbrFermions * $SInc));
	$MinArray{$NbrFermions} = ((&FindMinimum($TmpFile)) * $Scaling) / $NbrFermions;
	$MinArray{$NbrFermions} -= ($FactorInvN / $NbrFermions) + $FactorShift;
      }
    $NbrFermions += $NbrFermionsInc;
    $S += $SInc;
  }
&CreatePostScript(\%MinArray, $Caption, $PrintFlag);

# find a minimum in a file
#
# $_[0] = file name
# return value = ground state energy

sub FindMinimum
  {
    my $FileName = $_[0];
    my $Min;
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
	    if ($TmpArray[1] < $Min)
	      {
		$Min = $TmpArray[1];		
	      }
	  }
      }
    close (INFILE);
    return $Min;
  }

# create postscript graph from data file
#
# $_[0] = hash table containing datas
# $_[1] = print flag (1 if true)
# $_[2] = number of fermions

sub CreatePostScript
  {
    my $Datas = $_[0];
    my $Caption = $_[1];
    my $PrintFlag = $_[2];
    my $N;
    my $E;
    my $FileName = "fermions_coulomb_ground_.dat";
    my $FileName2 = "fermions_coulomb_ground_".$Caption."_FilledShells.dat";

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
    open (OUTFILE, ">$FileName2");
    while (($N, $E) = each (%$Datas))
      {
	my $TmpN = 1;
	my $TmpI = 0;
	while ($TmpN < $N)
	  {
	    $TmpI++;
	    $TmpN += ((2 * $TmpI) + 1);
	  }
	if ($TmpN == $N)
	  {
	    $N = 1.0 / $N;
	    print OUTFILE ($N." ".$E."\n");
	  }
      }
    close (OUTFILE);

    my $Delta = ($MaxGap - $MinGap) / 20.0;
    $MaxGap += $Delta;
    $MinGap -= $Delta;
#    $MinGap = 0.6;
#    $MaxGap = 0.9;
    $MinN--;
    $MaxN++;
    my $Tmp = 1.0 / $MinN;
    $MinN = 1.0 / $MaxN;
    $MaxN = $Tmp;
    $MinN = 0.0;
    my $TmpFileName = "tmp".time().".p";
    my $OutputFile = "fermions_coulomb_ground_".$Caption.".ps";
    my @TmpArray = split (/_/,  $OutputFile);
    my $Title = "gap ".$Caption;
    open (OUTFILE, ">$TmpFileName");
    print OUTFILE ("set xrange [".$MinN.":".$MaxN."]
set yrange [".$MinGap.":".$MaxGap."]
set xlabel \"1/N\"
set ylabel \"E/N\"
set size 1.0, 0.6
set nokey
set terminal postscript portrait enhanced \"Helvetica\" 14
set output \"".$OutputFile."\"
f(x)= a*x+b
fit f(x) \"".$FileName2."\" using 1:2 via a,b
plot \"".$FileName."\" using 1:2 title \"".$Title."\", f(x) with lines 1
");
    close (OUTFILE);
    `gnuplot $TmpFileName`;
    if ($PrintFlag == 1)
      {
	`lpr $OutputFile`;
      }
    `rm -f $TmpFileName`;
  }



