#!/usr/bin/perl -w

use strict 'vars';

if (!(defined($ARGV[0])))
  {
    die "usage: CleanGroundSpectrum spectrum\n";
  }
my $FileName = $ARGV[0];
open (INFILE, "$FileName");
my $TmpLine;
my $Flag = 0;
my $GroundLine = "";
my $WholeLz = "";
my $CurrentLz = 0;
my $NbrValues = 0;
my $FullSpectrum = "";
while (($Flag == 0) && (defined ($TmpLine = <INFILE>)))
  {    
    chomp ($TmpLine);
    my @TmpArray = split (/ /, $TmpLine);
    if ($CurrentLz == $TmpArray[0])
      {
	if ($GroundLine eq "")
	  {
	    $GroundLine = $TmpArray[0]." ".$TmpArray[1]."\n";
	    $WholeLz = $GroundLine;
	    $NbrValues++;
	  }
	else
	  {
	    $WholeLz .= $TmpArray[0]." ".$TmpArray[1]."\n";
	    $NbrValues++;
	  }
      }
    else
      {
	if ($NbrValues > 2)
	  {
	    $FullSpectrum .= $WholeLz;
	    $FullSpectrum .= $TmpArray[0]." ".$TmpArray[1]."\n";
	    $Flag = 1;
	  }
	else
	  {
	    $FullSpectrum .= $GroundLine;
	    $GroundLine = $TmpArray[0]." ".$TmpArray[1]."\n";
	    $WholeLz = $GroundLine;
	    $NbrValues = 1;
	    $CurrentLz = $TmpArray[0];
	  }
      }
  }
while (defined ($TmpLine = <INFILE>))
  {
    $FullSpectrum .= $TmpLine;
  }
close (INFILE);
print $FullSpectrum;
