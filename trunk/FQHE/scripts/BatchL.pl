#!/usr/bin/perl -w

use strict 'vars';

my $NbrFermions = $ARGV[0];
my $SValue = $ARGV[1];
my $LProgram = "/home/regnault/development/Physics/DiagHam/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereLValue";


my $TmpFile;
my %LValues;
my $OutputFileName;
foreach  $TmpFile (<*>)
  {
    if ((-f $TmpFile) && ($TmpFile=~ /\.vec$/))
      {
	print "processing ".$TmpFile."\n";
	my $Command = $LProgram." ".$TmpFile;
	my $LValue = `$Command`;
	my $LzValue = $TmpFile;
	$LzValue =~ s/^.*\_lz\_(-?\d+)\..*$/$1/;
	my @TmpArray = split (/\n/, $LValue);
	$TmpArray[1] =~ s/^.*\<L\> \= //;
	my $Index  = $TmpFile;
	$Index  =~ s/^.*\_lz\_-?\d+\.(\d+)\..*$/$1/;	    
	if (defined($LValues{$LzValue}))
	  { 
	    my $TmpArray2 = $LValues{$LzValue};
	    $$TmpArray2{$Index} = $TmpFile." ".$TmpArray[1];
	  }
	else
	  {
	    my %TmpArray2;
	    $TmpArray2{$Index} = $TmpFile." ".$TmpArray[1];
	    $LValues{$LzValue} = \%TmpArray2;
	  }
	$OutputFileName = $TmpFile;
      }
  }

$OutputFileName =~ s/\_lz\_.*/\.lvalues\.dat/;
unless (open (OUTFILE, ">".$OutputFileName))
  {
    die ("can't open ".$OutputFileName."\n");
  }
foreach $TmpFile(sort {$a <=> $b} (keys(%LValues)))
  {
    my $LzValue = $LValues{$TmpFile};
    my $TmpLz;
    foreach $TmpLz (sort {$a <=> $b} (keys(%$LzValue)))
      {
	print OUTFILE $$LzValue{$TmpLz}."\n";
      }
  }
close (OUTFILE);
