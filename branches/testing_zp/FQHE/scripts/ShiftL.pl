#!/usr/bin/perl -w

use strict 'vars';

if (!(defined($ARGV[0])))
  {
    die "usage: BosonsDeltaGraph bosons_delta_n_x_2s_y_l.dat\n";
  }
my $FileName = $ARGV[0];
my $NbrBosons = $FileName;
$NbrBosons =~ s/.*n\_(\d*).*/$1/;
my $S = $FileName;
$S =~ s/.*2s\_(\d*).*/$1/;

my $Flag = 0;
open (INFILE, $FileName);
my $TmpLine;
my $TmpFile = "";
foreach $TmpLine (<INFILE>)
  {
    chomp ($TmpLine);
    my @TmpArray = split (/ /, $TmpLine);
    $TmpArray[1] += (0.5 * $NbrBosons * $NbrBosons) / (0.5 * $S);
    $TmpFile .= $TmpArray[0]." ".$TmpArray[1]."\n";
  }
close (INFILE);

open (OUTFILE, ">$FileName");
print OUTFILE $TmpFile;
close (OUTFILE);
