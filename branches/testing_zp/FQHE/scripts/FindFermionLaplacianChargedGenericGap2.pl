#!/usr/bin/perl -w

use strict 'vars';

if (!(defined($ARGV[1])))
  {
    die "usage: FindGapGeneric ground.dat Caption [PrintFlag]";
  }
my $PrintFlag = 0;
my $GapFile = $ARGV[0];
my $Caption = $ARGV[1];
if (defined($ARGV[2]))
  {
    $PrintFlag = 1;
  }
unless (open (INFILE2, $GapFile))
  {
    die ("can't open $GapFile\n")
  }

#my %MinArray;
my %MinArrayElectron;
my %MinArrayHole;
my $TmpFile;
my $TmpLine; 
$TmpLine = <INFILE2>;
chomp ($TmpLine);
my $PFactor = 1;
my $QFactor = 1;
my $Shift = 1;
($PFactor, $QFactor, $Shift) = split (/ /, $TmpLine);
while ($TmpLine = <INFILE2>)
  {
    chomp ($TmpLine);
    if (($TmpLine ne "") && (!($TmpLine =~ /^\#/)))
      {
	print "--------------------------------------------------\n";
	my @Values = split (/ /, $TmpLine);
	my $NbrFermions = $Values[0];
	my $S = $Values[1];
	$TmpFile = "n_".$NbrFermions."/fermions_laplaciandelta_n_".$NbrFermions."_2s_".$S."_lz.dat";
	if (!(-e $TmpFile))
	  {
	    die ("file n_".$NbrFermions."/fermions_laplaciandelta_n_".$NbrFermions."_2s_".$S."_lz.dat does not exist\n")
	  }
	my $Quasi1 = &FindGround($TmpFile);
	$S = ($QFactor * $NbrFermions / $PFactor) - $Shift;
	my $Scaling = ($S * $PFactor) / ($NbrFermions * $QFactor);
	$Quasi1 *= $Scaling * $Scaling;	
	my $Ground = (($Values[4] / ($NbrFermions * $NbrFermions)) + ($Values[5] / $NbrFermions) + $Values[6]);	
	print $Quasi1." ".$Ground."\n";
	$Quasi1 -=  $Ground;
	$MinArrayElectron{$Values[0]} = $Quasi1;
	$NbrFermions = $Values[2];
	$S = $Values[3];
	$TmpFile = "n_".$NbrFermions."/fermions_laplaciandelta_n_".$NbrFermions."_2s_".$S."_lz.dat";
	if (!(-e $TmpFile))
	  {
	    die ("file n_".$NbrFermions."/fermions_laplaciandelta_n_".$NbrFermions."_2s_".$S."_lz.dat does not exist\n")
	  }
	my $Quasi2 = &FindGround($TmpFile);
	$S = ($QFactor * $NbrFermions / $PFactor) - $Shift;
	$Scaling = ($S * $PFactor) / ($NbrFermions * $QFactor);
	$Quasi2 *= $Scaling * $Scaling;		
	$Ground = (($Values[4] / ($NbrFermions * $NbrFermions)) + ($Values[5] / $NbrFermions) + $Values[6]);	
	print $Quasi2." ".$Ground."\n";
	$Quasi2 -= $Ground;	
	$MinArrayHole{$Values[2]} = $Quasi2;
      }
  }
close (INFILE2);

#&CreatePostScript(\%MinArray, $Caption, $PrintFlag);
&CreatePostScript(\%MinArrayElectron, \%MinArrayHole, $Caption, $PrintFlag);

# find ground state energy in a file
#
# $_[0] = file name
# $_[1] = L value where to find the excited state used to obtain the gap value
# return value = ground state energy

sub FindGround
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
    my $DataElectrons = $_[0];
    my $DataHoles = $_[1];
    my $Caption = $_[2];
    my $PrintFlag = $_[3];
    my $NHole;
    my $EHole;
    my $NElectron;
    my $EElectron;
    my $FileName = "fermions_laplaciandelta_chargedgap_".$Caption.".dat";
    open (OUTFILE, ">$FileName");
    my $MinN = 200;
    my $MaxN = 0;
    my $MinGap = 400;
    my $MaxGap = 0;
    while ((($NHole, $EHole) = each (%$DataHoles)) && (($NElectron, $EElectron) = each (%$DataElectrons)))
      {
	if ($MinN > $NHole)
	  {
	    $MinN = $NHole;
	  }
	if ($MaxN < $NHole)
	  {
	    $MaxN = $NHole;
	  }
	if ($MinGap > $EHole)
	  {
	    $MinGap = $EHole;
	  }
	if ($MaxGap < $EHole)
	  {
	    $MaxGap = $EHole;
	  }
	$NHole = 1.0 / $NHole;
	if ($MinN > $NElectron)
	  {
	    $MinN = $NElectron;
	  }
	if ($MaxN < $NElectron)
	  {
	    $MaxN = $NElectron;
	  }
	if ($MinGap > $EElectron)
	  {
	    $MinGap = $EElectron;
	  }
	if ($MaxGap < $EElectron)
	  {
	    $MaxGap = $EElectron;
	  }
	$NElectron = 1.0 / $NElectron;
	print ($NHole." ".$EHole." ".$NElectron." ".$EElectron."\n");
	print OUTFILE ($NHole." ".$EHole." ".$NElectron." ".$EElectron."\n");
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
    $MinN = 0.0;
    my $TmpFileName = "tmp".time().".p";
    my $OutputFile = "fermions_laplaciandelta_chargedgap_".$Caption.".ps";
    my @TmpArray = split (/_/,  $OutputFile);
    my $Title = "gap ".$Caption;
    open (OUTFILE, ">$TmpFileName");
    print OUTFILE ("set xrange [".$MinN.":".$MaxN."]
set yrange [".$MinGap.":".$MaxGap."]
set xlabel \"1/N\"
set ylabel \"E\"
set size 1.0, 0.6
set nokey
set terminal postscript portrait enhanced \"Helvetica\" 14
set output \"".$OutputFile."\"
f(x)= a*x+b
g(x)= m*x+n
fit f(x) \"".$FileName."\" using 1:2 via a,b
fit g(x) \"".$FileName."\" using 3:4 via m,n,p
plot \"".$FileName."\" using 1:2 title \"".$Title."\", \"".$FileName."\" using 3:4 title \"".$Title."\" with points 2, f(x) with lines 1, g(x) with lines 2
");
    close (OUTFILE);
    `gnuplot $TmpFileName`;
    if ($PrintFlag == 1)
      {
	`lpr $OutputFile`;
      }
    `rm -f $TmpFileName`;
  }



