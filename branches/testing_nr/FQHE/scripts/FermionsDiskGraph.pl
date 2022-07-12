#!/usr/bin/perl -w

use strict 'vars';

if (!(defined($ARGV[0])))
  {
    die "usage: FermionsDiskGraph nbr_fermions [xfig,eps,ps] [print_flag]\n";
  }
my $NbrFermions = $ARGV[0];
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
if (defined($ARGV[2]))
  {
    $PrintFlag = 1;
  }
my @ListFiles;
my $TmpFile;
foreach $TmpFile (<*>)
  {
    if ($TmpFile =~ /fermions\_disk\_.*\_n\_$NbrFermions.*\.dat/ )
      {
	push (@ListFiles, $TmpFile);
	print ($TmpFile."\n");
      }
  }

foreach $TmpFile (@ListFiles)
  {
    &CreatePostScript($TmpFile, $XFigFlag, $PrintFlag);
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
#	$TmpLine = <INFILE>;
#	$TmpLine = <INFILE>;
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
#	$TmpLine = <INFILE>;
#	$TmpLine = <INFILE>;
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

sub CreatePostScript
  {
    my $FileName = $_[0];
    my $XFigFlag = $_[1];
    my $PrintFlag = $_[2];
    my $Max;
    my $Min;
    &FindMinMax($FileName, 1, \$Min, \$Max);
#    $Max = 10.0;
    my $Delta = ($Max - $Min) / 20.0;
    $Max += $Delta;
    $Min -= $Delta;
    my $TmpFileName = "tmp".time().".p";
    my $OutputFile = $FileName;
    my $MaxL = $OutputFile;
    $MaxL =~ s/^.*\_lz\_(\d*).*$/$1/;
    my $N = $OutputFile;;
    $N =~ s/^.*\_n\_(\d*).*$/$1/;
    my $Title = "N = ".$N;
    my $MinL = ($N * ($N - 1)) / 2 - 1;
    ++$MaxL;
#    $MinL = 17;
#    $MaxL = 33;
#    $Min = -0.15;
#    $Max = 1.5;
    open (OUTFILE, ">$TmpFileName");
    print OUTFILE ("set xrange [".$MinL.":".$MaxL."]
set yrange [".$Min.":".$Max."]
set xlabel \"Total angular momentum Lz\"
set ylabel \"Energy\"\n");
    if ($XFigFlag >= 1)
      {
	$OutputFile =~ s/\.dat/\.fig/;
	print OUTFILE ("set size 1.0, 1.5
set terminal fig
set key bottom left
set output \"".$OutputFile."\"
plot \"".$FileName."\" using 1:2 title \"".$Title."\" with points pt 2\n");
      }
    else
      {
	$OutputFile =~ s/\.dat/\.ps/;
print OUTFILE ("set size 1.0, 0.6
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
	    $XFigFile .= $TmpLine;
	    $TmpLine = <INFILE>;
	    $XFigFile .= $TmpLine;
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



