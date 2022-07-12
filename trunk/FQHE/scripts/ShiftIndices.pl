#!/usr/bin/perl -w

if (!(defined($ARGV[1])))
  {
    die "usage: ShiftIndices.pl NbrColumns SourceFile";
  }

my $IndexShift=-1;
my $NbrColumns=$ARGV[0];
open (INFILE, $ARGV[1]);
my $TmpLine;
my $OutFileName = $ARGV[1]."DH";
open (OUTFILE, ">$OutFileName");

foreach $TmpLine (<INFILE>)
  {
    chomp ($TmpLine);
    my @TmpArray = split (/\s/, $TmpLine);
    for (my $c=0; $c<$NbrColumns; ++$c)
      {
	print OUTFILE (($TmpArray[$c]+$IndexShift)."\t");
      }
    for (my $c=$NbrColumns; $c<=$#TmpArray-1; ++$c)
      {
	print OUTFILE ($TmpArray[$c]."\t");
      }
    print OUTFILE ($TmpArray[$#TmpArray]."\n");
  }
close(OUTFILE);
