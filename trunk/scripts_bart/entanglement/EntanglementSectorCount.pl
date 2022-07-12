#!/usr/bin/perl -w

use strict 'vars';

use Getopt::Long;

my $SpectrumFile="";
my $RowIndex = 3;
my $RowData = 5;
my $Threshold;
my $LatexFlag = 0;

my $Result = GetOptions ("spectrum=s" => \$SpectrumFile, "threshold=s" => \$Threshold, "index:i" => \$RowIndex, "data:i" => \$RowData);

#if (($SpectrumFile eq "") || (!(-e $SpectrumFile)) || (!($Eigenvalue =~ /^[\+\-]?\d*\.?\d*e?\d+$/)) || (!()))
if ($SpectrumFile eq "")
  {
    die ("usage: EntanglementSectorCount.pl --spectrum file_name --threshold max_xi [--index ROW_INDEX --data ROW_DATA]\n");
  }

my %Degeneracy;

unless (open (INFILE, $SpectrumFile))
  {
    die ("can't open ".$SpectrumFile."\n");
  }
my $TmpLine;
while (defined($TmpLine = <INFILE>))
  {
    chomp ($TmpLine);
    if (! ($TmpLine =~ /^#/ ))
      {
	my @TmpArray = split (/ /, $TmpLine);
	if ($TmpArray[$RowData] < $Threshold)
	  {
	    if (defined($Degeneracy{$TmpArray[$RowIndex]}))
	      {
		$Degeneracy{$TmpArray[$RowIndex]}++;
	      }
	    else
	      {
		$Degeneracy{$TmpArray[$RowIndex]} = 1;
	      }
	  }
      }
  }
close (INFILE);

print ("State count up the threshold of $Threshold\nIndex\tCount\n");

foreach $TmpLine (sort {$a <=> $b} (keys(%Degeneracy)))
  {
    print $TmpLine." ".$Degeneracy{$TmpLine}."\n";
  }
