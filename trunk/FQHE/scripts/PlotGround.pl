#!/usr/bin/perl -w

use strict 'vars';

if (!(defined($ARGV[0])))
  {
    die "usage: BosonsDeltaGraph nbr_bosons print_flag\n";
  }
my $NbrBosons = $ARGV[0];
my $PrintFlag = 0;
if (defined($ARGV[1]))
  {
    $PrintFlag = 1;
  }
my %MinArray;
my $TmpFile;
foreach $TmpFile (<*>)
  {
    if (($TmpFile =~ /bosons\_delta\_n\_$NbrBosons.*\_l\./ ) ||  ($TmpFile =~ /bosons\_coulomb\_n\_$NbrBosons.*\_l\./))
      {
	my $TmpFile2 = $TmpFile;
	$TmpFile2 =~ s/.*2s\_(\d*).*/$1/;
        $MinArray{$TmpFile2} = &FindMin($TmpFile);
	print ("$TmpFile2 ".$MinArray{$TmpFile2}."\n");
      }
  }
&CreatePostScript(\%MinArray, 0, $NbrBosons);

# find minimum in a file
#
# $_[0] = file name
# return value = ground state energy

sub FindMin
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
# $_[2] = number of bosons

sub CreatePostScript
  {
    my $Datas = $_[0];
    my $PrintFlag = $_[1];
    my $N = $_[2];
    my $S;
    my $E;
    my $TmpDataFileName = "tmp".time().".dat";
    open (OUTFILE, ">$TmpDataFileName");
    my $EMax = 0;
    my $EMin = 0;
    my $SMax = 0;
    my $SMin = 0;
    my $Flag = 0;
    while (($S, $E) = each (%$Datas))
      {
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
    my $Title = "Ground state energy for N = ".$N;
    my $OutputFile = "bosons_delta_ground_state_n_".$N.".ps";
    open (OUTFILE, ">$TmpFileName");
    print OUTFILE ("set xrange [".$SMin.":".$SMax."]
set yrange [".$EMin.":".$EMax."]
set xlabel \"Angular Momentum l\"
set ylabel \"E(L)\"
set size 1.0, 0.6
set terminal postscript portrait enhanced \"Helvetica\" 14
set output \"".$OutputFile."\"
plot \"".$TmpDataFileName."\" using 1:2 title \"".$Title."\"
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



