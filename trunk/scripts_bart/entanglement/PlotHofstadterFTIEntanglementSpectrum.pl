#!/usr/bin/perl -w
#
# script for plotting output of FTIEntanglementSpectrumParticlePartition using gnuplot
#
use strict 'vars';


# if (!defined($ARGV[0]))
#   {
#     print("usage PlotFTIEntanglementSpectrum.pl -n NA [Directory|single file]\n");
#     exit(1);
#   }

my $Machine="g";   # choose label to use one of the preset paths
my $BuildDirectory;
my $MPIBuildDirectory;
my $MPIRun;
my $MPIRunPrefix="";
my $FTIPrograms="FTI/src/Programs/FTI/";
my $FCIPrograms="FTI/src/Programs/FCI/";
if ($Machine =~ m/s1/)
  {
    $BuildDirectory = "/scratch/gm360/DiagHam/buildLA/";
    $MPIBuildDirectory = "/scratch/gm360/DiagHam/buildMPI/";
    $MPIRun ="mpirun";
  }
elsif ($Machine =~ m/g/)
  {
    $BuildDirectory = "/Users/gunnar/DiagHam/build1/";
    $MPIBuildDirectory = "/Users/gunnar/DiagHam/buildMPI/";
    $MPIRun ="mpirun";
  } 
elsif ($Machine =~ m/t/)
  {
    $BuildDirectory = "~/DiagHam2/build/";
    $MPIBuildDirectory = "~/DiagHam2/build/";
    $MPIRun ="mpirun";
  }
else
  {
    die ("Machine name not recognized\n");
  }

my $NA=-1;
my $OverrideMax=-1;
while( (defined($ARGV[0])&&$ARGV[0] =~ /^-/ ))
  {
    if ( $ARGV[0] =~ /-n/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $NA = $ARGV[0];
	    $NA =~ s/-n//;
	  }
	else
	  {
	    shift(@ARGV);
	    $NA = $ARGV[0];
	  }
      }
    if ( $ARGV[0] =~ /-m/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $OverrideMax = $ARGV[0];
	    $OverrideMax =~ s/-m//;
	  }
	else
	  {
	    shift(@ARGV);
	    $OverrideMax = $ARGV[0];
	  }
      }
    shift(@ARGV);
  }


if (!defined($ARGV[0]))
  {
    print("usage PlotFTIEntanglementSpectrum.pl [-n N_A] [-m OverrideMax] Directory|'single file'\n");
    exit(1);
  }

my @ListFiles;
my $ReadMore = 1;
if ( defined($ARGV[0]) )
  {
    if ( -d$ARGV[0] )
      {
	chdir($ARGV[0]);
      }
    else # single file:
      {
	@ListFiles = $ARGV[0];
	$ReadMore=0;
      }
  }
my $TmpFile;
if ($ReadMore == 1) # array still empty?
  {
    foreach $TmpFile (<*>)
      {
	if ($TmpFile =~ /.*\.full\.parent/)
	  {
	    push (@ListFiles, $TmpFile);
	  }	
      }
  }


foreach $TmpFile (@ListFiles)
  {
    my $Max;
    my $Min;
    &CreatePostScript($TmpFile, $NA, $OverrideMax);
    print ("\n\n");
#    &FindMinMax($TmpFile, 1, \$Min, \$Max, 0, 0, 15);
#    print ($TmpFile." ".$Min." ".$Max."\n");
  }


# find minimum and maximum values in a file
#
# $_[0] = file name
# $_[1] = column where to search - starting to number from zero
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
	    if ($TmpLine =~ m/^\s*#+/)
	      {
		# print("ignoring comment line ".$TmpLine."\n");
	      }
	    else
	      {
		@TmpArray = split (/ /, $TmpLine);
		if (($TmpArray[$ColumnConstraint] <= $MaxConstraint) && ($TmpArray[$ColumnConstraint] >= $MinConstraint)) {
		  if ($Flag == 0) {
		    $$Min = $TmpArray[$Column];
		    $$Max = $$Min;
		    $Flag = 1;
		  } else {
		    if ($TmpArray[$Column] < $$Min) {
		      $$Min = $TmpArray[$Column];
		    }
		    if ($TmpArray[$Column] > $$Max) {
		      $$Max = $TmpArray[$Column];
		    }
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
# $_[0] = density matrix file name
# $_[1] = N_A value for entanglement subspace to consider

sub CreatePostScript
  {
    my $FileName = $_[0];
    my $NA = $_[1];
    my $OverrideMax = $_[2];
    my $PrintFlag = 0;
    my $Max;
    my $Min;

    $FileName =~ /X\_(\d+)\_Y\_(\d+)\_q\_(\d+)\_n\_(\d+)\_x\_(\d+)\_y\_(\d*)\_/;
    my $FileX = $1;
    my $FileY = $2;
    my $Fileq = $3;
    my $Filen = $4;
    my $Filex = $5;
    my $Filey = $6;

    if ($NA<=0)
      {
	use integer;
	$NA = $Filen / 2;
	print ("Defaulting to NA=n/2=".$NA."\n");
      }

    my $PESFileName = $FileName;
    $PESFileName =~ s/\.full\.parent/.na_${NA}.parentspec/;

    print ("Searching file $PESFileName\n");
    if (! -e $PESFileName)
      {
	system($BuildDirectory.$FTIPrograms."FTIEntanglementSpectrum --show-counting --show-minmaxkya -n $NA --particle-entanglement $FileName");
	if (! -e $PESFileName)
	  {
	    die ("Failed creating entanglement spectrum file ".$PESFileName."\n");
	  }
      }
    else
      {
	print ("Using existing file $PESFileName\n");
      }
    &FindMinMax($PESFileName, 5, \$Min, \$Max, 0, 0, 14);

    if ($OverrideMax>0.0)
      {
	$Max = $OverrideMax;
      }
    my $Delta = ($Max - $Min) / 20.0;
    $Max += $Delta;
    $Min -= $Delta;
    my $PlotFileName =  $FileName; # previously random: "tmp".time().".p";
    $PlotFileName =~ s/\.full\.parent/.na_${NA}.parentspec.gp/;
    #my @TmpArray = split (/_/,  $OutputFile);
  
    my $Title = "N = ".$Filen." N_A = ${NA} $Filex x $Filey MUC with ($Fileq; $FileX x $FileY)";

    my $MaxX = $Filex * $Filey;

    my $OutputFile = $FileName;
    $OutputFile =~ s/\.full\.parent/.na_${NA}.parentspec.ps/;

    open (OUTFILE, ">$PlotFileName");
    print OUTFILE ("set xrange [-1:".($MaxX+1)."]
set yrange [".$Min.":".$Max."]
set xlabel \"k_x*L_y + k_y\"
set ylabel \"{/Symbol x}(k)\"
set size 1.0, 0.6
set terminal postscript portrait enhanced \"Helvetica\" 14
set output \"".$OutputFile."\"
plot \"".$PESFileName."\" using 4:6 title \"".$Title."\" with points lc rgbcolor \"red\" lt 6 pt 1 ps 2
");
    close (OUTFILE);
    # print plotfile on screen
    open (INFILE, $PlotFileName);
    my $TmpLine;
    foreach $TmpLine (<INFILE>)
      {
	print $TmpLine;
      }
    close (INFILE);

    `gnuplot $PlotFileName`;
    if ($PrintFlag == 1)
      {
	`lpr $OutputFile`;
      }
    #`rm -f $PlotFileName`;
  }


