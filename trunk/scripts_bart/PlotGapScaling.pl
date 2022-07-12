#!/usr/bin/perl -w

use strict 'vars';

my $Degeneracy=3; # Laughlin nu = 1/3
my $UseFermions = 0;

while( (defined($ARGV[0])&&$ARGV[0] =~ /^-/ ))
  {
    if ( $ARGV[0] =~ /-d/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $Degeneracy = $ARGV[0];
	    $Degeneracy =~ s/-d//;
	  }
	else
	  {
	    shift(@ARGV);
	    $Degeneracy = $ARGV[0];
	  }
      }
    if ( $ARGV[0] =~ /-f/ )
      {
	print ("Using scaling for fermions (q^2N)\n");
	$UseFermions=1;
      }
    shift(@ARGV);
  }

if (!defined($ARGV[0]))
  {
    print("usage PlotGapScaling.pl gap files'\n");
    print("input: gap files from FindLatticeGap.pl\n");
    exit(1);
  }

my $File;

foreach $File (@ARGV)
  {
    &AnalyseScaling($File);
  }

# subroutine PlotValue
# $_[0] gap value
# $_[1] number of sublattices
sub PlotValue
  {
    my $Gap = $_[0];
    my $Q = $_[1];
    if ($UseFermions==1)
      {
	return $Gap*$Q*$Q;
      }
    else
      {
	return $Gap*$Q;
      }
  }

sub AnalyseScaling
  {
    my $FileName = $_[0];
    my $MinGS=0.0;
    my $MaxGS=1.0;
    my $Exc=2.0;
    open (INFILE, $FileName);
    my $TmpLine;
    my $MaxVal=0.0;
    my $OutputFileName="$FileName.q";
    open (OUTFILE, ">$OutputFileName") or die("Could not open output $OutputFileName");

    print OUTFILE ("# 1/Q 1/N Gap Spread Ratio MinGS MaxGS Exc N FileName\n");

    foreach $TmpLine (<INFILE>)
      {
	chomp ($TmpLine);
	if ($TmpLine =~ m/^\s*#+/)
	  {
	    # print("ignoring comment line ".$TmpLine."\n");
	  }
	else
	  {
	    my @TmpArray = split (/ /, $TmpLine);
	    if ($TmpArray[1] =~ m/N\/A/)
	      {
		print "# skipping line $TmpLine\n";
	      }
	    else
	      {
		my $DataSourceFile=$TmpArray[8];
		my $Gap=$TmpArray[1];
		$DataSourceFile =~ /\_q\_(\d+)\_n\_(\d+)\_/;
		my $QInv = 1.0/$1;
		my $N = $2;
		if ($MaxVal < PlotValue($Gap,$1) )
		  { $MaxVal = PlotValue($Gap,$1) }
		print OUTFILE "$QInv $TmpLine\n";
	      }
	  }
      }
    close (INFILE);
    close (OUTFILE);

    open (GPFILE, ">$FileName.gp");
    my $OutputFile = "$FileName.ps";
    print GPFILE "set size 1.0, 0.6
set terminal postscript portrait enhanced \"Helvetica\" 14 color
set output \"".$OutputFile."\"
set xlabel \"1/q\"\n";
if ($UseFermions==1)
  {
    print GPFILE "set ylabel \"Delta q^2\"\n";
  }
else
  {
    print GPFILE "set ylabel \"Delta q\"\n";
  }
    print GPFILE "set xrange [0:*]
set yrange [0:".1.25*$MaxVal."]\n";
    print GPFILE "plot \"$OutputFileName\" u 1:(\$3/\$1**".(1+$UseFermions),") title \"".$FileName."\" with points lc rgbcolor \"red\" lt 6 pt 1 ps 2\n";
    close (GPFILE);
    `gnuplot $FileName.gp`;
    print "Created plot $OutputFile\n";
    
  }



# read in input file to array
# $_[0] = input file
sub get_data 
    { 
      my $FileName = $_[0];
      my @links; 
      open(my $Input, "<", $FileName) or die "cannot open < $FileName: $!"; 
      while(<$Input>)
	{ 
	  chomp;
	  if ($_ =~ m/^\s*#+/)
	    {
	      # print("ignoring comment line ".$TmpLine."\n");
	    }
	  else
	    {	      
	      my $tuple; 
	      @$tuple = split(/\s+/, $_); 
	      push(@links, $tuple); 
	    }
	} 
      return \@links; 
    }


# read in input file and sort according to the given field
# $_[0] = input file
# $_[1] = field index
sub sort_and_store 
  { 
    my $refToLinks = get_data($_[0]);
    my $field = $_[1];
    my @sorted_links = sort { $a->[$field] <=> $b->[$field] } @$refToLinks; 
    return \@sorted_links; 
  }


# read in input file and sort according to the given field
# $_[0] = input file
sub SortAndPlot
  {    
    my $SpectrumName = $_[0];
    $SpectrumName =~ /\_n\_(\d+)\_x\_(\d+)\_y\_(\d*)\_/;
    my $FileN = $1;
    my $FileLx = $2;
    my $FileLy = $3;
    #print ("Geometry: $FileLx x $FileLy, Degeneracy = $Degeneracy\n");
    my $refToSpectrum = sort_and_store($SpectrumName, 2);
    
    # record first line as ground state:
    my $line = $refToSpectrum->[0];
   
    my $MinGS=$line->[2];
    my $MaxGS;
    my $Exc;
    
    my $Index=0; 
    my $N=0;
    while ($N < $Degeneracy)
      {
	$line = $refToSpectrum->[$Index];
	#print ("$line->[0] $line->[1] $line->[2]\n");
	my $kx = $line->[0];
	my $ky = $line->[1];
	my $Value = $line->[2];
	if (($kx<=$FileLx/2) &&  ($ky<=$FileLy/2)) ## only count states, not their degeneracy partners
	  {
	    $MaxGS=$line->[2];
	    my $deg=0; # calculate degeneracy of the given level
	    if ( ($kx==0) || (( $FileLx%2 == 0) && ($kx== $FileLx/2)))
	      {
		$deg=1;
	      }
	    else
	      {
		$deg=2; # count extra state for symmetry factor       
	      }
	    if ( ($ky==0) || (( $FileLy%2 == 0) && ($ky== $FileLy/2)))
	      {
		$deg *=1;
	      }
	    else
	      {
		$deg *=2; # count extra state for symmetry factor
	      }
	    $N+=$deg;
	  }
	++$Index;
      }
    if ($N != $Degeneracy)
      {
	print ("$SpectrumName N/A\n");
	print ("#Error: ground-state degeneracy exceeded! Go check file $SpectrumName\n");
	$Exc=$MaxGS;
      }
    else
      {
	$line = $refToSpectrum->[$Index];
	$Exc=$line->[2];
	my $Gap = $Exc-$MaxGS;
	my $Spread = $MaxGS-$MinGS;
	my $Ratio = $Gap / ($Spread + 1e-10);
	print ((1.0/$FileN)." $Gap $Spread $Ratio $MinGS $MaxGS $Exc $FileN $SpectrumName\n");
      }    
  }


sub AnalyseGapAllInOne
  {
    my $FileName = $_[0];

    open (my $inFH, q{<}, ${FileName})
      or die qq{open: < $FileName: $!\n};
    open my $outFH, q{>}, q{spw964070.out}
      or die qq{open: > spw964070.out: $!\n};
 
    print $outFH
      map  { $_->[ 0 ] }
      sort { $a->[ 2 ] <=> $b->[ 2 ] }
      map  { [ $_, split ] }
      <$inFH>;
 
    close $inFH 
      or die qq{close: < spw964070.in: $!\n};
    close $outFH
      or die qq{close: > spw964070.out: $!\n};
  }
