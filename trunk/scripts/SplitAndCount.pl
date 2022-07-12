#!/usr/bin/perl -w

use strict 'vars';

if (!(defined($ARGV[0])))
  {
    die "usage: SplitAndCount file\n";
  }

my $FileName = $ARGV[0];

open (INFILE, $FileName);
my $TmpLine;
while (defined($TmpLine = <INFILE>))
{
    chomp ($TmpLine);
    my @TmpArray = split (/ /, $TmpLine);
    my $TmpValue;
    my $Count = 0;
    foreach $TmpValue (@TmpArray)
    {
	print $Count." ".$TmpValue."\n";
	$Count++;
    }
}
close (INFILE);
