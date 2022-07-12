#!/usr/bin/perl -w

use strict 'vars';

if (!(defined($ARGV[0])))
  {
    die "usage: BosonsDeltaGraph nbr_bosons [xfig,eps,ps] [print] [maxE]\n";
  }
my $NbrBosons = $ARGV[0];
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
my @ListFiles;
my $TmpFile;
foreach $TmpFile (<*>)
  {
    if (($TmpFile =~ /bosons\_delta\_n\_$NbrBosons.*\_l\.dat/ ) || ($TmpFile =~ /bosons\_coulomb\_n\_$NbrBosons.*\_l\.dat/) || 
	($TmpFile =~ /bosons\_dipolar\_n\_$NbrBosons.*\_l\.dat/) || ($TmpFile =~ /bosons\_v2\_n\_$NbrBosons.*\_l\.dat/) || ($TmpFile =~ /bosons\_.*\_n\_$NbrBosons.*\_l\.dat/))
      {
	push (@ListFiles, $TmpFile);
      }
  }

foreach $TmpFile (@ListFiles)
  {
    &CreatePostScript($TmpFile, $XFigFlag, $PrintFlag, 14, $MaxE);
    print ("\n\n");
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
    &FindMinMax($FileName, 1, \$Min, \$Max, 0, 0, $MaxL);
    if (defined ($_[4]))
      {
	$Max = $_[4];
      }
    my $Delta = ($Max - $Min) / 20.0;
    $Max += $Delta;
    $Min -= $Delta;
    my $TmpFileName = "tmp".time().".p";
    my $OutputFile = $FileName;
    $FileName =~ /n\_(\d+)\_2s\_(\d*)\_/;
    my $Title = "N = ".$1."  2S = ".$2;
    open (OUTFILE, ">$TmpFileName");
    print OUTFILE ("set xrange [-1:".($MaxL + 1)."]
set yrange [".$Min.":".$Max."]\n");
    if ($XFigFlag >= 1)
      {
	$OutputFile =~ s/\_l\.dat/\.fig/;
	print OUTFILE ("set xlabel \"L\" font \"default,14\"
set ylabel \"energy[g]\" font \"default,14\"
set size 1.0, 1.5
set terminal fig
set key bottom right
set output \"".$OutputFile."\"
plot \"".$FileName."\" using 1:2 title \"".$Title."\" with points pt 2\n");
      }
    else
      {
	$OutputFile =~ s/\_l\.dat/\.ps/;
print OUTFILE ("set xlabel \"L\"
set ylabel \"energy[g]\"
set size 1.0, 0.6
set terminal postscript portrait enhanced \"Helvetica\" 14
set output \"".$OutputFile."\"
plot \"".$FileName."\" using 1:2 title \"".$Title."\"\n");
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
    if (($PrintFlag == 1) && ($XFigFlag == 0))
      {
	`lpr $OutputFile`;
      }
    `rm -f $TmpFileName`;
    if ($XFigFlag >= 1)
      {
	open (INFILE, $OutputFile);
	my $XFigFile = "";
	while ((defined ($TmpLine = <INFILE>)) && (!($TmpLine =~ /N \= \d*/)))
	  {
	    $XFigFile .= $TmpLine;
	  }
	$XFigFile .= $TmpLine;	
	while (defined ($TmpLine = <INFILE>))
	  {
	    chomp $TmpLine;
	    $TmpLine =~ s/     / /;
	    my @Tmp = split (/ /, $TmpLine);
	    $Tmp[3]++;
	    $XFigFile .= $Tmp[0]." ".$Tmp[1]." ".$Tmp[2]." ".$Tmp[3]." ".$Tmp[4]." ".$Tmp[5]." ".$Tmp[6]." ".
	      $Tmp[7]." ".$Tmp[8]."     ".$Tmp[9]." ".$Tmp[10]." ".$Tmp[11]." ".$Tmp[12]." ".$Tmp[13]." ".$Tmp[14]." ".$Tmp[15]."\n";
	    $TmpLine = <INFILE>;
	    chomp $TmpLine;
	    $TmpLine =~ s/^\s*//;
	    if ($TmpLine ne "")
	      {
		@Tmp = split (/ /, $TmpLine);
		my $Shift = int (($Tmp[4] - $Tmp[0]) / 2);
		$Tmp[0] -= $Shift;
		$Tmp[4] += $Shift;
		$XFigFile .= "	 ".$Tmp[0]." ".$Tmp[1]." ".$Tmp[0]." ".$Tmp[1]." ".$Tmp[4]." ".$Tmp[1]."\n";
	      }
	    $TmpLine = <INFILE>;
	    $TmpLine = <INFILE>;
	  }
	open (OUTFILE , ">$OutputFile");
	print OUTFILE $XFigFile;
	close (OUTFILE);
	if ($XFigFlag == 2)
	  {
	    my $EpsOutputFile = $OutputFile;
	    $EpsOutputFile =~ s/\.fig/\.eps/;
	    `fig2dev -L eps $OutputFile $EpsOutputFile`;
	  }
      }    
  }



