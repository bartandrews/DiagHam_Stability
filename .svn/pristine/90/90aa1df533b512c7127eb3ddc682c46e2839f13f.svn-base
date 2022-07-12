#!/usr/bin/perl -w

use strict 'vars';

if (!(defined($ARGV[0])))
  {
    die "usage: FigCrossToDot file.fig radius_mag red green blue\n";
  }

my $FileName = $ARGV[0];
my $Magnification = 1.0;
if (defined($ARGV[1]))
  {
    $Magnification = $ARGV[1];
  }
my $Red = 1.0;
if (defined($ARGV[2]))
  {
    $Red = $ARGV[2];
  }
my $Green = 1.0;
if (defined($ARGV[3]))
  {
    $Green = $ARGV[3];
  }
my $Blue = 1.0;
if (defined($ARGV[4]))
  {
    $Blue = $ARGV[4];
  }

unless (open (INFILE, $FileName))
  {
    die "can't open $FileName\n";
  }
my $XFigFile = "";
my $TmpLine;
while (defined ($TmpLine = <INFILE>))
  {
    if ($TmpLine =~ /^2 1 0 1 -1 -1 10 0 -1/)
      {
	my $TmpLine2 = <INFILE>;
	my $TmpLine3 = $TmpLine2;
	chomp ($TmpLine3);
	$TmpLine3 =~ s/^\s*//;
	my @TmpValues = split (/ /, $TmpLine3);	
	if ((!defined($TmpValues[4])) && ($TmpValues[0] != $TmpValues[2]) && ($TmpValues[1] != $TmpValues[3]))
	  {
	    $TmpLine2 = <INFILE>;
	    $TmpLine2 = <INFILE>;
	    $TmpValues[0] += $TmpValues[2];
	    $TmpValues[0] *= 0.5;
	    $TmpValues[1] += $TmpValues[3];
	    $TmpValues[1] *= 0.5;
	    $XFigFile .= "1 4 0 1 4 4 50 -1 20 1 0.0000 ".$TmpValues[0]." ".$TmpValues[1]." 53 53 2775 -150 2950 -75\n"

	  }
	else
	  {
	    $XFigFile .= $TmpLine;
	    $XFigFile .= $TmpLine2;
	  }
      }
    else
      {
	if ($TmpLine =~ /^4 2 -1 0 -1 0 10/)
	  {
	    $TmpLine =~ s/^4 2 -1 0 -1 0 10/4 2 -1 0 -1 0 12/;
	    $XFigFile .= $TmpLine;
	  }
	else
	  {
	    $XFigFile .= $TmpLine;
	  }
      }
  }
$XFigFile .= $TmpLine;	
while (defined ($TmpLine = <INFILE>))
  {
  }
close (INFILE);
print $XFigFile
