#!/usr/bin/perl -w

use strict 'vars';

my $TmpFile;
foreach $TmpFile (<*>)
  {
    if (($TmpFile =~ /^bosons\_disk\_hardcore\_nbody\_2\_n\_(\d+)\_2s\_(\d+)\_lz\_(\d+)\.0\.vec$/) && ($1 eq $2) && ($1 eq $3))
      {
	print $TmpFile."\n";
	my $NbrParticles = $1;
	unless (open(OUTFILE, ">tmp.dat"))
	  {
	    die ("can't create tmp.dat\n");
	  }
	print OUTFILE ("coefficients = 1.0\ndescriptions=0 ".$NbrParticles);
	my $Pos = 2;
	while ($Pos <= $NbrParticles)
	  {
	    print OUTFILE (" 0");
	    $Pos++;
	  }
	print OUTFILE ("\n");
	close(OUTFILE);
	print "/home/regnault/development/Physics/DiagHam/build/src/Programs/QHE/QHEOnDisk/QHEForgeEigenstate --bosons -p ".$NbrParticles." -l ".$NbrParticles." -i tmp.dat -o tmp.vec\n";
	system("/home/regnault/development/Physics/DiagHam/build/src/Programs/QHE/QHEOnDisk/QHEForgeEigenstate --bosons -p ".$NbrParticles." -l ".$NbrParticles." -i tmp.dat -o tmp.vec");
	system("/home/regnault/development/Physics/DiagHam/build/src/Programs/QHE/QHEBosonsDeltaOverlap --use-exact tmp.vec --exact-state ".$TmpFile);
#	unlink ("tmp.dat");
#	unlink ("tmp.vec");
      }
  }

