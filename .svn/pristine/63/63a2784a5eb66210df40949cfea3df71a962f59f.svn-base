#!/usr/bin/perl -w

use strict 'vars';

if (!(defined($ARGV[1])))
  {
    die "usage: MergePartialSpectrum [--help --require-overlap] file1 file2 [file 3] ...\n";
  }



my @Files;
my $TmpFile;
my %EnableOptions;
foreach $TmpFile (@ARGV)
  {
    if (!($TmpFile =~ /^\-\-\w/))
      {
	if (!(-e $TmpFile))
	  {
	    die ("can't find file ".$TmpFile."\n");
	  }
	else
	  {
	    push (@Files, $TmpFile);
	  }
      }
    else
      {
	$TmpFile =~ s/^\-\-//;
	$EnableOptions{$TmpFile} = "";
      }
  }
if (defined($EnableOptions{"help"}))
  {
    die "usage: MergePartialSpectrum [--help --require-overlap] file1 file2 [file 3] ...\n";
  }
my $RequireSpectrumOverlap = 0;
if (defined($EnableOptions{"require-overlap"}))
  {
    $RequireSpectrumOverlap = 1;
  }
&OrderFiles (\@Files);

my $TmpSpectrum = "";
my $CurrentLzValue = 0;
my $Index = 0;
my $PreviousFile = $Files[0];

foreach $TmpFile (@Files)
{
  print $TmpFile."\n";
  if ($Index == 0)
    {
      open (INFILE, $TmpFile);
      my $TmpLine;
      foreach $TmpLine (<INFILE>)
	{
	  $TmpSpectrum .= $TmpLine;
	  chomp ($TmpLine);
	  $TmpLine =~ s/\ .*//;
	  if ($TmpLine > $CurrentLzValue)
	    {
	      $CurrentLzValue = $TmpLine;
	    }
	}
      close (INFILE);
    }
  else
    {
      my $MinLz = &TestOverlap ($PreviousFile, $TmpFile);
      if (($MinLz < 0) && ($RequireSpectrumOverlap == 1))
	{
	  die ("no overlap between files ".$PreviousFile." and ".$TmpFile."\n");
	}
      else
	{
	  open (INFILE, $TmpFile);
	  my $TmpLine;
	  foreach $TmpLine (<INFILE>)
	    {
	      my $TmpValue = $TmpLine;
	      chomp ($TmpValue);
	      $TmpValue =~ s/\ .*//;
	      if ($TmpValue >= $MinLz)
		{
		  $TmpSpectrum .= $TmpLine;
		}
	    }
	  close (INFILE);
	}
    }

  $PreviousFile = $TmpFile;
  $Index++;
}

print $TmpSpectrum;

# order partial spectrum files
#
# $_[0] = reference on partial spectrum file name array

sub OrderFiles
  {
    my $Files = $_[0];
    my %MinimumLzValues;
    my $TmpFile;
    foreach $TmpFile (@$Files)
      {
	open (INFILE, $TmpFile);
	my $TmpLine = <INFILE>;
	chomp ($TmpLine);
	$TmpLine =~ s/\ .*//;
	print $TmpLine." ".$TmpFile."\n";
	%MinimumLzValues = (%MinimumLzValues, $TmpFile, $TmpLine);
	close (INFILE);	
      }
    @$Files = sort { $MinimumLzValues{$a} <=> $MinimumLzValues{$b} } keys %MinimumLzValues
  }

# test if an overlap between two files occurs
#
# $_[0] = first file name 
# $_[1] = second file name (must have greater lz values than the lowest lz value of the first file
# return value = first lz value of the second file that is not in the first file (negative is the files overlap but with different spectrum)

sub TestOverlap
  {
    my $File1 = $_[0];
    my $File2 = $_[1];
    open (INFILE2, $File2);
    my $TmpLine = <INFILE2>;
    chomp ($TmpLine);
    $TmpLine =~ s/\ .*//;
    my $MinLz = $TmpLine;
    close (INFILE2);
    open (INFILE1, $File1);
    my $MaxLz = 0;
    foreach $TmpLine (<INFILE1>)
      {
	chomp ($TmpLine);
	$TmpLine =~ s/\ .*//;
	if ($TmpLine > $MaxLz)
	  {
	    $MaxLz = $TmpLine;
	  }
      }
    close (INFILE1);
    if (($MaxLz - $MinLz) < (-1))
      {
	return -1;
      }
    while (($MinLz <= $MaxLz) && ($MinLz >= 0))
      {
	my @Spectrum1 =  &FindSpectrum($File1, $MinLz);
	my @Spectrum2 =  &FindSpectrum($File2, $MinLz);
	my $Index  = 0;
	my $TmpValue;
	foreach $TmpValue (@Spectrum1)
	  {
	    if (abs(($TmpValue - $Spectrum2[$Index]) / $TmpValue) > 1e-13)
	      {
		print ("possible error at ".$MinLz.": ".$TmpValue." ".$Spectrum2[$Index]." ".abs(($TmpValue - $Spectrum2[$Index]) / $TmpValue)."\n");
	      }
	    $Index++;
	  }
	$MinLz++;
      }
    return $MinLz;
  }

# find the part of the spectrum with a given lz value
#
# $_[0] =  name of the file containing the spectrum 
# $_[1] =  lz value
# return value = array containing the spectrum

sub FindSpectrum
  {
    my $FileName = $_[0];
    my $Lz = $_[1];
    open (INFILE, $FileName);
    my $TmpLine;
    my @TmpSpectrum;
    foreach $TmpLine (<INFILE>)
      {
	chomp ($TmpLine);
	my $Value = $TmpLine;
	$Value =~ s/^\d*\ //;
	$TmpLine =~ s/\ .*//;
	if ($TmpLine == $Lz)
	  {
	    push (@TmpSpectrum, $Value);
	  }
      }
    close (INFILE);
    return @TmpSpectrum;
  }

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
    my $N = $_[1];
    my $S;
    my $E;
    my $OutputFile = "bosons_delta_ground_state_n_".$N.".log";
    open (OUTFILE, ">$OutputFile");
    while (($S, $E) = each (%$Datas))
      {
	my $ShiftedE = $E + (0.5 * $N * $N / (sqrt(0.5 * $S)));
	print OUTFILE ($S." ".$E." ".$ShiftedE."\n");
      }
    close (OUTFILE);
  }



