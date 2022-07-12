#!/usr/bin/perl -w

use strict 'vars';

if (!(defined($ARGV[0])))
  {
    die "usage: PlotEigenvalueDensity spectrum [xfig,eps,ps] [print] [maxE]\n";
  }
my $Spectrum = $ARGV[0];
my $ColumnIndex = 0;
my $XFigFlag = 0;
if( (defined($ARGV[1])) && ($ARGV[1] eq "xfig"))
  {
    $XFigFlag = 1;
  }
if( (defined($ARGV[1])) && ($ARGV[1] eq "eps"))
  {
    $XFigFlag = 2;
  }
my $PrintFlag = 0;
my $MaxE;
if (defined($ARGV[2]))
  {
    if ($ARGV[2] eq "print")
      {
	$PrintFlag = 1;
	if (defined($ARGV[3]))
	    {
	      $MaxE = $ARGV[3];
	    }
      }
    else
      {
	$MaxE = $ARGV[2];
      }
  }

my $TmpFile = $Spectrum;
$TmpFile =~ s/\.dat/\_de\.dat/;
ConvertToHistogram($Spectrum, $TmpFile, $ColumnIndex);
&CreatePostScript($TmpFile, $XFigFlag, $PrintFlag, 14, $MaxE);


# convert a spectrum file into a histogram file
#
# $_[0] = input file name
# $_[1] = output file name
# $_[2] = column corresponding to datas to sort

sub ConvertToHistogram
  {
    my $InputFileName = $_[0];
    my $OutputFileName = $_[1];
    my $Column = $_[2];
    my $Min;
    my $Max;
    FindMinMax($InputFileName, $Column, \$Min, \$Max);
    my $TmpLine;
    unless (open (INFILE, $InputFileName))
      {
	die ("can't open ".$InputFileName."\n");
      }
    my $NbrBoxes = 200;
    my $Step  = ($Max - $Min) / $NbrBoxes;    
    my @Histogram;
    my $Pos = 0;
    while ($Pos < $NbrBoxes)
      {
	$Histogram[$Pos] = 0;
	$Pos++;
      }
    while (defined($TmpLine = <INFILE>))
      {
	chomp ($TmpLine);
	if ($TmpLine ne "")
	  {
	    my @TmpArray = split (/ /, $TmpLine);
	    $Pos = ($TmpArray[$Column] - $Min) / $Step;
	    $Histogram[$Pos]++;
	  }
      }
    close (INFILE);
    $Pos = 0;
    unless (open (OUTFILE, ">".$OutputFileName))
      {
	die ("can't open ".$OutputFileName."\n");
      }
    my $ScalingFactor = 10.0 / $NbrBoxes;
    while ($Pos < $NbrBoxes)
      {
#	print OUTFILE ($Min + ($Pos * $Step))." ".$Histogram[$Pos]."\n";
	print OUTFILE ($Pos * $ScalingFactor)." ".$Histogram[$Pos]."\n";
	$Pos++;
      }  
    close (OUTFILE);
}

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
      }
    else
      {
	open (INFILE, $FileName);
	my $TmpLine;
	my @TmpArray;
	$TmpLine = <INFILE>;
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
# $_[0] = file name
# $_[1] = xfig output flag: 1 if true, 0 if none (use ps instead), 2 if xfig+eps
# $_[2] = print flag (1 if true)
# $_[3] = maximum value of L to display
# $_[4] = maximum energy to display

sub CreatePostScript
  {
    my $FileName = $_[0];
    my $XFigFlag = $_[1];
    my $PrintFlag = $_[2];
    my $MaxL = $_[3];
    my $Max;
    my $Min;
#    &FindMinMax($FileName, 1, \$Min, \$Max, 0, 0, $MaxL);
#    if (defined ($_[4]))
#      {
#	$Max = $_[4];
#      }
#    my $Delta = ($Max - $Min) / 20.0;
#    $Max += $Delta;
#    $Min -= $Delta;
    my $TmpFileName = "tmp".time().".p";
    my $OutputFile = $FileName;
    open (OUTFILE, ">$TmpFileName");
#    print OUTFILE ("set xrange [-1:".($MaxL + 1)."]
#set yrange [".$Min.":".$Max."]\n");
    if ($XFigFlag >= 1)
      {
	$OutputFile =~ s/\.dat/\.fig/;
	print OUTFILE ("set xlabel \"E\" font \"default,14\"
set ylabel \"D(E)\" font \"default,14\"
set size 1.0, 1.5
set terminal fig
set key bottom right
set output \"".$OutputFile."\"
plot \"".$FileName."\" using 1:2 notitle with histeps\n");
      }
    else
      {
	$OutputFile =~ s/\.dat/\.ps/;
print OUTFILE ("set xlabel \"E\"
set ylabel \"D(E)\"
set size 1.0, 0.6
set terminal postscript portrait enhanced \"Helvetica\" 14
set output \"".$OutputFile."\"
plot \"".$FileName."\" using 1:2 notitle with histeps\n");
     }
    close (OUTFILE);
    open (INFILE, $TmpFileName);
    my $TmpLine;
    foreach $TmpLine (<INFILE>)
      {
	print $TmpLine;
      }
    close (INFILE);

    `gnuplot $TmpFileName`;
    if ($PrintFlag == 1)
      {
	`lpr $OutputFile`;
      }
    `rm -f $TmpFileName`;
  }



