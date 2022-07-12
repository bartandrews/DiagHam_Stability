#!/usr/bin/perl -w

use strict 'vars';

if (!(defined($ARGV[4])))
  {
    die "usage: FindGapGeneric StartN StartS NInc SInc Caption [PrintFlag]";
  }
my $PrintFlag = 0;
my $NbrFermions = $ARGV[0];
my $S = $ARGV[1];
my $NbrFermionsInc = $ARGV[2];
my $SInc = $ARGV[3];
my $Caption = $ARGV[4];
if (defined($ARGV[5]))
  {
    $PrintFlag = 1;
  }
my %MinArray;
my $TmpFile;
my $TmpFileElectron;
my $TmpFileHole;
my $ChargeQ = 1.0;
my $ChargeP = 3.0;
while ($NbrFermions <= 40)
  {
    $TmpFile = "n_".$NbrFermions."/fermions_coulomb_n_".$NbrFermions."_2s_".$S."_lz.dat";
    $TmpFileElectron = "n_".$NbrFermions."/fermions_coulomb_n_".$NbrFermions."_2s_".($S - 1)."_lz.dat";
    $TmpFileHole = "n_".$NbrFermions."/fermions_coulomb_n_".$NbrFermions."_2s_".($S + 1)."_lz.dat";
    if ((-e $TmpFile) && (-e $TmpFileElectron) && (-e $TmpFileHole))
      {
	my $Scaling = sqrt((($S - 1) * $NbrFermionsInc) / ($NbrFermions * $SInc));
	$MinArray{$NbrFermions} = &FindGround($TmpFileElectron) * $Scaling;
	$MinArray{$NbrFermions} += (($ChargeQ * $ChargeQ) / ($ChargeP * $ChargeP * sqrt (2.0 * ($S - 1)))) * $Scaling;
	$Scaling = sqrt((($S + 1) * $NbrFermionsInc) / ($NbrFermions * $SInc));
	$MinArray{$NbrFermions} += &FindGround($TmpFileHole) * $Scaling;
	$MinArray{$NbrFermions} += (($ChargeQ * $ChargeQ) / ($ChargeP * $ChargeP * sqrt (2.0 * ($S + 1)))) * $Scaling;
	$Scaling = sqrt(($S * $NbrFermionsInc) / ($NbrFermions * $SInc));
	$MinArray{$NbrFermions} -= 2.0 * &FindGround($TmpFile) * $Scaling;
      }
    $NbrFermions += $NbrFermionsInc;
    $S += $SInc;
  }
&CreatePostScript(\%MinArray, $Caption, $PrintFlag);

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
    my $Datas = $_[0];
    my $Caption = $_[1];
    my $PrintFlag = $_[2];
    my $N;
    my $E;
    my $FileName = "fermions_coulomb_gap_".$Caption.".dat";
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
    $MinGap = 0.045;
    $MinN--;
    $MaxN++;
    my $Tmp = 1.0 / $MinN;
    $MinN = 1.0 / $MaxN;
    $MaxN = $Tmp;
    $MinN = 0.0;
    my $TmpFileName = "tmp".time().".p";
    my $OutputFile = "fermions_coulomb_chargedgap_".$Caption.".ps";
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



