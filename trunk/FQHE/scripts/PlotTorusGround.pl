#!/usr/bin/perl -w

use strict 'vars';

if (!(defined($ARGV[0])))
  {
    die "usage: PlotTorusGround nbr_quantum print_flag\n";
  }

my $NbrQuantum = $ARGV[0];
my $PrintFlag = 0;
if (defined($ARGV[1]))
  {
    $PrintFlag = 1;
  }
my $TmpFile;
my $MinNbrParticle = 3;
while ((!(-e "n_$MinNbrParticle")) && ($MinNbrParticle < 60))
  {
    ++$MinNbrParticle;
  }
if ($MinNbrParticle == 60)
  {
    die ("can't find any directory containing datas\n")
  }
my $Flag = 0;
while ((-e "n_$MinNbrParticle") && ($MinNbrParticle < 60) && ($Flag == 0))
  {
    chdir ("n_$MinNbrParticle");
    foreach $TmpFile (<*>)
      {
	if ($TmpFile =~ /[^\_]*\_torus\_[^\_]*\_n\_$MinNbrParticle\_2s\_$NbrQuantum\_.*\.dat/)
	  {
	    $Flag = 1;
	  }
      }
    if ($Flag == 0)
      {
	++$MinNbrParticle
      }
    chdir ("..");
  }

if ($MinNbrParticle == 60)
  {
    die ("can't find any data atatched to 2s=$NbrQuantum\n");
  }

$Flag = 0;
my $MaxNbrParticle = $MinNbrParticle;
while ((-e "n_$MaxNbrParticle") && ($MaxNbrParticle < 60) && ($Flag == 0))
  {
    chdir ("n_$MaxNbrParticle");
    $Flag = 1;
    foreach $TmpFile (<*>)
      {
	if ($TmpFile =~ /[^\_]*\_torus\_[^\_]*\_n\_$MaxNbrParticle\_2s\_$NbrQuantum\_.*\.dat/)
	  {
	    $Flag = 0;
	  }
      }
    if ($Flag == 0)
      {
	++$MaxNbrParticle
      }
    chdir ("..");
  }
$MaxNbrParticle--;
print ($MinNbrParticle." ".$MaxNbrParticle."\n");
if (($MaxNbrParticle - $MinNbrParticle) < 1)
  {
    die ("not enough datas\n");
  }

my $NbrGroundState = 0;
my $NbrParticle = $MinNbrParticle;
my @AbsoluteGroundStates;
while ($NbrParticle <= $MaxNbrParticle)
  {
    print "processing n=".$NbrParticle." p=".$NbrQuantum."...\n";
    chdir ("n_$NbrParticle");
    my $FileName;
    foreach $TmpFile (<*>)
      {
	if ($TmpFile =~ /[^\_]*\_torus\_[^\_]*\_n\_$NbrParticle\_2s\_$NbrQuantum\_.*\.dat/)
	  {
	    $FileName = $TmpFile;
	  }
      }
    print ($FileName."\n");
    my $Min;
    my $Max;
    &FindMinMax ($FileName, 1, \$Min, \$Max);
    print ($Min." ".$Max."\n");
    push (@AbsoluteGroundStates, $Min);
    chdir ("..");
    $NbrParticle++;
  }

$NbrParticle = $MinNbrParticle + 1;
my %GroundStates;
while ($NbrParticle < $MaxNbrParticle)
  {
    %GroundStates = (%GroundStates, ($NbrParticle/$NbrQuantum),  
		     (($AbsoluteGroundStates[$NbrParticle - $MinNbrParticle + 1] / ($NbrParticle + 1)) + 
		      ($AbsoluteGroundStates[$NbrParticle - $MinNbrParticle - 1] / ($NbrParticle - 1)) -
		      2.0 * ($AbsoluteGroundStates[$NbrParticle - $MinNbrParticle] / $NbrParticle)) * $NbrParticle);
    $NbrParticle++;
  }

&CreatePostScript(\%GroundStates, 0, $NbrQuantum);

# find minimum and maximum values in a file
#
# $_[0] = file name
# $_[1] = column where to search
# $_[2] = reference on min value
# $_[3] = reference on max value
# $_[4] = additonnal constraint on another column
# $_[5] = max value for the constraint
# $_[6] = min value for the constraint

sub FindMinMax
  {
    my $FileName = $_[0];
    my $Column = $_[1];
    my $Min = $_[2];
    my $Max = $_[3];
    if (defined($_[4]))
      {
	my $ColumnConstraint = $_[4];
	my $MinConstraint = $_[5];
	my $MaxConstraint = $_[6];
	my $Flag = 0;
	open (INFILE, $FileName);
	my $TmpLine;
	my @TmpArray;
	foreach $TmpLine (<INFILE>)
	  {
	    chomp ($TmpLine);
	    @TmpArray = split (/ /, $TmpLine);
	    if (($TmpArray[$ColumnConstraint] <= $MaxConstraint) && ($TmpArray[$ColumnConstraint] >= $MinConstraint))
	      {
		if ($Flag == 0)
		  {
		    $$Min = $TmpArray[$Column];
		    $$Max = $$Min;
		    $Flag = 1;
		  }
		else
		  {
		    if ($TmpArray[$Column] < $$Min)
		      {
			$$Min = $TmpArray[$Column];
		      }
		    if ($TmpArray[$Column] > $$Max)
		      {
			$$Max = $TmpArray[$Column];
		      }
		  }
	      }
	  }
	close (INFILE);
      }
    else
      {
	open (INFILE, $FileName);
	my $TmpLine;
	my @TmpArray;
	$TmpLine = <INFILE>;
	chomp ($TmpLine);
	@TmpArray = split (/ /, $TmpLine);
	$$Min = $TmpArray[$Column];
	$$Max = $$Min;
	foreach $TmpLine (<INFILE>)
	  {
	    chomp ($TmpLine);
	    @TmpArray = split (/ /, $TmpLine);
	    if ($TmpArray[$Column] < $$Min)
	      {
		$$Min = $TmpArray[$Column];
	      }
	    if ($TmpArray[$Column] > $$Max)
	      {
		$$Max = $TmpArray[$Column];
	      }
	  }
	close (INFILE);
      }
}

# create postscript graph from data file
#
# $_[0] = hash table containing datas
# $_[1] = print flag (1 if true)
# $_[2] = number of flux quantum

sub CreatePostScript
  {
    my $Datas = $_[0];
    my $PrintFlag = $_[1];
    my $N = $_[2];
    my $S;
    my $E;
    my @FillingFactors = sort(keys(%$Datas));
    my $TmpDataFileName = "tmp".time().".dat";
    open (OUTFILE, ">$TmpDataFileName");
    my $EMax = 0;
    my $EMin = 0;
    my $SMax = 0;
    my $SMin = 0;
    my $Flag = 0;
    foreach $S (@FillingFactors)
#    while (($S, $E) = each (%$Datas))
      {
	$E = $$Datas{$S};
	print OUTFILE ($S." ".$E."\n");
	if ($Flag == 0)
	  {
	    $EMax = $E;
	    $EMin = $E;
	    $SMax = $S;
	    $SMin = $S;
	    $Flag = 1;
	  }
	else
	  {
	    if ($EMax < $E)
	      {
		$EMax = $E;
	      }
	    else
	      {
		if ($EMin > $E)
		  {
		    $EMin = $E;
		  }
	      }
	    if ($SMax < $S)
	      {
		$SMax = $S;
	      }
	    else
	      {
		if ($SMin > $S)
		  {
		    $SMin = $S;
		  }
	      }
	  }
      }
    print ($EMax." ".$EMin."\n") ;
    print ($SMax." ".$SMin."\n") ;
    close (OUTFILE);
    my $Delta = ($EMax - $EMin) / 20.0;
    $EMax += $Delta;
    $EMin -= $Delta;
    $SMax++;
    $SMin--;
    my $TmpFileName = "tmp".time().".p";
    my $Title = "Ground state energy for Nv = ".$N;
    my $OutputFile = "bosons_torus_ground_state_n_".$N.".ps";
    open (OUTFILE, ">$TmpFileName");
    print OUTFILE ("set xrange [".$SMin.":".$SMax."]
set yrange [".$EMin.":".$EMax."]
set xlabel \"filling factor\"
set ylabel \"E(L)\"
set size 1.0, 0.6
set terminal postscript portrait enhanced \"Helvetica\" 14
set output \"".$OutputFile."\"
plot \"".$TmpDataFileName."\" using 1:2 title \"".$Title."\" with lines 1
");
    close (OUTFILE);

    `gnuplot $TmpFileName`;
    if (($PrintFlag == 1) && (-e $OutputFile))
      {
	`lpr $OutputFile`;
      }
    `rm -f $TmpFileName`;
    `rm -f $TmpDataFileName`;
  }

