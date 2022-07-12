#!/usr/bin/perl -w

use strict 'vars';

my $Degeneracy=3; # Laughlin nu = 1/3
my $UseSymmetry = 0;

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
    if ( $ARGV[0] =~ /-s/ )
      {
	print ("Assuming symmetry k_x <-> -k_x, k_y <-> -k_y for k!=0 sectors\n");
	$UseSymmetry=1;
      }
    shift(@ARGV);
  }

if (!defined($ARGV[0]))
  {
    print("usage FindGap.pl [-d Deg(=3)] [-s] spectrum files'\n");
    print("-s: assume symmetry k <-> -k in counting degeneracy\n");
    exit(1);
  }

my $File;

print ("# 1/N Gap Spread Ratio MinGS MaxGS Exc N FileName\n");

foreach $File (@ARGV)
  {
    &AnalyseGap2($File);
  }

sub AnalyseGap
  {
    my $FileName = $_[0];
    my $MinGS=0.0;
    my $MaxGS=1.0;
    my $Exc=2.0;
    open (INFILE, $FileName);
    my $TmpLine;
    my $LineIdx=0;

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
	    if ($LineIdx==0)
	      {
		$MinGS = $TmpArray[2];
	      }
	    if ($LineIdx==$Degeneracy-1)
	      {
		$MaxGS = $TmpArray[2];
	      }
	    if ($LineIdx==$Degeneracy)
	      {
		$Exc = $TmpArray[2];
		last;
	      }
	    ++$LineIdx;
	  }
      }
    close (INFILE);
    if ($LineIdx<$Degeneracy)
      {
	print ("$FileName N/A\n");
      }
    else
      {
	my $Gap = $Exc-$MaxGS;
	my $Spread = $MaxGS-$MinGS;
	my $Ratio = $Gap / ($Spread + 1e-10);
	print ("$FileName $Gap $Spread $Ratio $MinGS $MaxGS $Exc\n");
      }
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
sub AnalyseGap2
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
