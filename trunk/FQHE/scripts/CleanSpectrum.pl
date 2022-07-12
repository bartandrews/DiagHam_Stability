#!/usr/bin/perl -w

use strict 'vars';

if (!(defined($ARGV[0])))
  {
    die "usage: CleanSpectrum spectrum\n";
  }
my $FileName = $ARGV[0];
open (INFILE, "$FileName");
my $TmpLine;
foreach $TmpLine (<INFILE>)
  {
    chomp ($TmpLine);
    my @TmpArray = split (/ /, $TmpLine);
    if (defined($TmpArray[1]))
      {
	print $TmpArray[1]."\n";
      }
  }
close (INFILE);
