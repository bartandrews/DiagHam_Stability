#!/usr/bin/perl -w

use strict 'vars';

my $InputFile = $ARGV[0];
$InputFile =~ /\_n\_(\d+)\_/;
my $NbrParticles = $1;
my $OutputName = "laughlin";
my $SpectrumType = "L";
my $DimensionFile = "bosons_n_8_2s_14.dim";
$OutputName .= "_".$NbrParticles;

unless (open(INFILE, $InputFile))
  {
    die "can't open ".$InputFile."\n";
  }
my %Spectrum;
my $TmpLine;
while (defined($TmpLine = <INFILE>))
  {
    chomp ($TmpLine);
    if ($TmpLine ne "")
      {
	 my @TmpArray = split (/ /, $TmpLine);
	 if (defined($Spectrum{$TmpArray[0]}))
	   {
	     my $Tmp = $Spectrum{$TmpArray[0]};
	     push (@$Tmp, $TmpArray[1]);	     
	   }
	 else
	   {
	     my @Tmp;
	     push (@Tmp, $TmpArray[1]);
	     $Spectrum{$TmpArray[0]} = \@Tmp;
	   }
      }
  }
close (INFILE);

my @Dimensions;
unless (open(INFILE, $DimensionFile))
  {
    die "can't open ".$InputFile."\n";
  }
if ($SpectrumType eq "Lz")
  {
    while (defined($TmpLine = <INFILE>))
      {
	chomp($TmpLine);
	if ($TmpLine =~ /^Lz \= /)
	  {
	    $TmpLine =~ s/^Lz \= //;
	    @Dimensions = split(/ /, $TmpLine);
	  }
      }
  }
else
  {
    while (defined($TmpLine = <INFILE>))
      {
	chomp($TmpLine);
	if ($TmpLine =~ /^L \= /)
	  {
	    $TmpLine =~ s/^L \= //;
	    @Dimensions = split(/ /, $TmpLine);
	  }
      }
  }
close (INFILE);

mkdir ($OutputName);
chdir ($OutputName);
my $Momentum;
my $PartialSpectrum;
while (($Momentum, $PartialSpectrum) = each(%Spectrum))
  {
    my $NbrEigenvalues = $#$PartialSpectrum + 1;
    if ((!(defined($Dimensions[$Momentum]))) || ($Dimensions[$Momentum] != $NbrEigenvalues))
      {
	die ("dimension mismatch for ".$SpectrumType."=".$Momentum." (get ".$NbrEigenvalues." need ".$Dimensions[$Momentum].")");
      }
    my $TmpOutputFile = "l_".$Momentum."_n_".$NbrEigenvalues;
    unless (open(OUTFILE, ">".$TmpOutputFile))
      {
	die "can't create ".$TmpOutputFile."\n";
      }
    my $Eigenvalue;
    foreach $Eigenvalue (sort {$a <=> $b} (@$PartialSpectrum))
      {
	print OUTFILE $Eigenvalue."\n";
      }
    close OUTFILE;
  }
chdir ("..");
my $Command = "tar -czvf ".$OutputName.".tar.gz ".$OutputName."/*";
system ($Command);
`rm -rf $OutputName`;
