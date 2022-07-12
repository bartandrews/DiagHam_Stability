#!/usr/bin/perl -w

use strict 'vars';

if (!(defined($ARGV[1])))
  {
    die "usage: FindFermionLaplacianChargedGenericGap.pl ground.dat Caption [PrintFlag]";
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

my %MinArrayElectron;
my %MinArrayHole;
my $TmpFile;
my $TmpFileGround1;
my $TmpFileGround2;
my $TmpLine; 
my $PFactor = 1;
my $QFactor = 1;
my $Shift = 1;
my $Flag = 0;
while (($Flag == 0) && ($TmpLine = <INFILE2>))
  {
    chomp ($TmpLine);
    $TmpLine =~ s/^\s*//;
    $TmpLine =~ s/\s*$//;
    if (($TmpLine ne "") && (!($TmpLine =~ /^\#/)))
      {
	($PFactor, $QFactor, $Shift) = split (/ /, $TmpLine);
	$Flag = 1;
      }
  }
$Flag = 0;
while ($TmpLine = <INFILE2>)
  {
    chomp ($TmpLine);
    $TmpLine =~ s/^\s*//;
    $TmpLine =~ s/\s*$//;
    if (($TmpLine ne "") && (!($TmpLine =~ /^\#/)))
      {
	print "--------------------------------------------------\n";
	my @Values = split (/ /, $TmpLine);
	my $NbrFermions = $Values[0];
	my $S = $Values[1];
	my $Ground;
	$TmpFile = "n_".$NbrFermions."/fermions_laplaciandelta_n_".$NbrFermions."_2s_".$S."_lz.dat";
	if ((-e $TmpFile) && (&FindApproximativeGround($NbrFermions, $PFactor, $QFactor, $Shift, \$Ground) == 0))
	  {
	    my $Quasi1 = ((&FindGround($TmpFile) + &ShiftEnergy($S, $NbrFermions)));
	    $S = ($QFactor * $NbrFermions / $PFactor) - $Shift;
	    my $Scaling = ($S * $PFactor) / ($NbrFermions * $QFactor);
	    $Scaling *= $Scaling;
	    $Quasi1 *= $Scaling;
	    $Quasi1 -= $Ground;
	    print ($Quasi1." ".$Ground."\n");
	    $MinArrayElectron{$Values[0]} = $Quasi1;
	  }
	$NbrFermions = $Values[2];
	$S = $Values[3];
	$TmpFile = "n_".$NbrFermions."/fermions_laplaciandelta_n_".$NbrFermions."_2s_".$S."_lz.dat";
	if ((-e $TmpFile) && (&FindApproximativeGround($NbrFermions, $PFactor, $QFactor, $Shift, \$Ground) == 0))
	  {
	    my $Quasi2 = ((&FindGround($TmpFile) + &ShiftEnergy($S, $NbrFermions)));
	    $S = ($QFactor * $NbrFermions / $PFactor) - $Shift;
	    my $Scaling = ($S * $PFactor) / ($NbrFermions * $QFactor);
	    $Scaling *= $Scaling;	
	    $Quasi2 *= $Scaling;
	    $Quasi2 -= $Ground;	
	    print ($Quasi2." ".$Ground."\n");
	    if (($Flag == 1) || ($Values[2] == $Values[0]))
	      {
		$Flag = 1;
		$MinArrayHole{$Values[2]} = $MinArrayElectron{$Values[2]} + $Quasi2;
	      }
	    else
	      {
		$MinArrayHole{$Values[2]} = $Quasi2;
	      }
	  }
      }
  }
close (INFILE2);

if ($Flag == 1)
  {
    &CreateSumPostScript(\%MinArrayHole, $Caption, $PrintFlag);
  }
else
  {
    &CreatePostScript(\%MinArrayElectron, \%MinArrayHole, $Caption, $PrintFlag);
  }

# find ground state energy in a file
#
# $_[0] = file name
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

# find ground state energy in a file using linear fit if needed
#
# $_[0] = number of fermions of the corresponding state
# $_[1] = fraction numerator
# $_[2] = fraction denominator
# $_[3] = shift to the filling factor definition
# $_[4] = reference on the ground state energy
# return value = 0 if the value has be evaluated 

sub FindApproximativeGround
  {
    my $NbrFermions = $_[0];
    my $PFactor = $_[1];
    my $QFactor = $_[2];
    my $Shift = $_[3];
    my $Ground = $_[4];
    my $SMin = (($QFactor * $NbrFermions) / $PFactor) - $Shift;
    my $SMax = $SMin;
    my $NbrFermionsMin = $NbrFermions;
    my $NbrFermionsMax = $NbrFermions;
    if ($SMin != int($SMin))
      {
	$SMax = int($SMin) + 1;
	$NbrFermionsMax = (($SMax + $Shift) * $PFactor) / $QFactor;
	while ($NbrFermionsMax != int($NbrFermionsMax))
	  {
	    $SMax++;
	    $NbrFermionsMax = (($SMax + $Shift) * $PFactor) / $QFactor;
	  }
	$SMin = int($SMin);
	$NbrFermionsMin = (($SMin + $Shift) * $PFactor) / $QFactor;
	while ($NbrFermionsMin != int($NbrFermionsMin))
	  {
	    $SMin--;
	    $NbrFermionsMin = (($SMin + $Shift) * $PFactor) / $QFactor;
	  }
      }
    else
      {
	if (-e "n_".$NbrFermionsMin."/fermions_laplaciandelta_n_".$NbrFermionsMin."_2s_".$SMin."_lz.dat")
	  {
	    $$Ground = &FindGround("n_".$NbrFermionsMin."/fermions_laplaciandelta_n_".$NbrFermionsMin."_2s_".$SMin."_lz.dat") 
	      * sqrt(($SMin * $PFactor) / ($NbrFermionsMin * $QFactor)) / $NbrFermionsMin;
	    return 0;
	  }
	else
	  {
	    return -1;
	  }
      }
    $TmpFileGround1 = "n_".$NbrFermionsMin."/fermions_laplaciandelta_n_".$NbrFermionsMin."_2s_".$SMin."_lz.dat";
    $TmpFileGround2 = "n_".$NbrFermionsMax."/fermions_laplaciandelta_n_".$NbrFermionsMax."_2s_".$SMax."_lz.dat";
    if (-e $TmpFileGround2)
      {
	if ((!(-e $TmpFileGround1)))
	  {
	    $TmpFileGround1 = $TmpFileGround2;
	    $NbrFermionsMin = $NbrFermionsMax;
	    $SMin = $SMax;
	    $NbrFermionsMax += $PFactor;
	    $SMax += $QFactor;
	    $TmpFileGround2 = "n_".$NbrFermionsMax."/fermions_laplaciandelta_n_".$NbrFermionsMax."_2s_".$SMax."_lz.dat";
	    if (!(-e $TmpFileGround2))
	      {
		return -1;
	      }
	  }
      }
    else
      {
	if (!(-e $TmpFileGround1))
	  {
	    return -1;
	  }
	$TmpFileGround2 = $TmpFileGround1;
	$NbrFermionsMax = $NbrFermionsMin;
	$SMax = $SMin;
	$NbrFermionsMin -= $PFactor;
	$SMin -= $QFactor;
	$TmpFileGround1 = "n_".$NbrFermionsMin."/fermions_laplaciandelta_n_".$NbrFermionsMin."_2s_".$SMin."_lz.dat";
	if (!(-e $TmpFileGround1))
	  {
	    return -1;
	  }
      }
    my $Scaling = ($SMin * $PFactor) / ($NbrFermionsMin * $QFactor);	
    my $Ground1 = (((&FindGround($TmpFileGround1) + &ShiftEnergy($SMin, $NbrFermionsMin))
		   * $Scaling * $Scaling / $NbrFermionsMin);	
    $Scaling = ($SMax * $PFactor) / ($NbrFermionsMax * $QFactor);
    my $Ground2 = (((&FindGround($TmpFileGround2) + &ShiftEnergy($SMax, $NbrFermionsMax))
		   * $Scaling * $Scaling / $NbrFermionsMax);	
    $$Ground = (((($Ground1 - $Ground2) / ((1.0 / $NbrFermionsMin) - (1.0 / $NbrFermionsMax))) * 
		 ((1.0 / $NbrFermions) - (1.0 / $NbrFermionsMax)) + $Ground2) * $NbrFermions);
    return 0;
  }

# evaluate the energy shift
#
# $_[0] = 2S value
# $_[1] = number of value
# return value = energy shift

sub ShiftEnergy
  {
    my $S = $_[0];
    my $NbrFermions = $_[1];
    return 0.0;
  }
# create postscript graph from data file
#
# $_[0] = hash table containing electron datas
# $_[1] = hash table containing hole datas
# $_[2] = print flag (1 if true)
# $_[3] = number of fermions

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
    $MinGap = 0.0;
    $MinN--;
    $MaxN++;
    my $Tmp = 1.0 / $MinN;
    $MinN = 1.0 / $MaxN;
    $MaxN = $Tmp;
    $MinN = 0;
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
fit g(x) \"".$FileName."\" using 3:4 via m,n
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



# create postscript graph from data file summin quasihole and quasielectron energy
#
# $_[0] = hash table containing datas
# $_[1] = print flag (1 if true)
# $_[2] = number of fermions

sub CreateSumPostScript
  {
    my $Datas = $_[0];
    my $Caption = $_[1];
    my $PrintFlag = $_[2];
    my $N;
    my $E;
    my $FileName = "fermions_laplaciandelta_chargedgap_".$Caption.".dat";
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
g(x)= m*x*x+n*x+p
fit f(x) \"".$FileName."\" using 1:2 via a,b
fit g(x) \"".$FileName."\" using 1:2 via m,n,p
plot \"".$FileName."\" using 1:2 title \"".$Title."\", f(x) with lines 1, g(x) with lines 2
");
    close (OUTFILE);
    `gnuplot $TmpFileName`;
    if ($PrintFlag == 1)
      {
	`lpr $OutputFile`;
      }
    `rm -f $TmpFileName`;
  }



