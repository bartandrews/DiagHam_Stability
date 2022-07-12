#!/usr/bin/perl -w

use strict 'vars';

my $ScalingFactor = 2.8549;
my $NbrValue = 3;
unless (open(INFILE, $ARGV[0]))
  {
    die ("can't open ".$ARGV[0]."\n");
  }

my $TmpLine;
while (defined($TmpLine = <INFILE>))
  {
    chomp ($TmpLine);
    if ($TmpLine =~ /\-d \d*\.\d* \-e \d*\.\d* \-z \d*\.\d* \-a \d*\.\d*/)
      {
	print $TmpLine."\n";
      }
    else
      {
	if ($TmpLine =~ /Nbr of iterations \=/)
	  {
	    $TmpLine = <INFILE>;
	    chomp ($TmpLine);
	    my @Values = split (/ /, $TmpLine);
	    my $Index = 1;
	    while ($Index < $NbrValue)
	      {
		print (($Values[$Index] - $Values[$Index - 1]) * $ScalingFactor);
		print " ";
		++$Index;
	      }
	    print "\n";
	  }
      }
  }
close (INFILE);

