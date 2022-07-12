#!/usr/bin/perl -w

use strict 'vars';

use Getopt::Long;


my $PathToProgram = "";
my $PathToDatas = ".";
my $FractionDescriptionFile = "";
my $NbrNBody = 3;

my $Result = GetOptions ("prog:s" => \$PathToProgram, "data:s" => \$PathToDatas, "fraction=s" => \$FractionDescriptionFile,  "nbrnbody:i" => \$NbrNBody); 


if ($PathToProgram eq "")
  {
    $PathToProgram = "/home/regnault/development/DMRG/DiagHam/src/Programs/QHE/QHEOnSphere";
  }

my @Fractions;
unless (open (INFILE, $FractionDescriptionFile))
  {
    die ("can't open ".$FractionDescriptionFile."\n");
  }
my $TmpLine;
foreach $TmpLine (<INFILE>)
  {
    chomp ($TmpLine);
    push (@Fractions, $TmpLine);
  }
close (INFILE);

unless (chdir ($PathToDatas))
  {
    die ("directory ".$PathToDatas." does not exist or is not accessible\n");
  }

foreach $TmpLine (@Fractions)
  {
    my @TmpArray = split (/\|/, $TmpLine);
    my $NbrParticles = $TmpArray[3];
    while ($NbrParticles < 4)
      {
	$NbrParticles += $TmpArray[3];
      }
    my $FillingFactor = $TmpArray[0] / $TmpArray[1];
    my $SValue = ($NbrParticles / $TmpArray[0]) * $TmpArray[1] - $TmpArray[2];
    while ($NbrParticles < 200)
      {
	$SValue = ($NbrParticles / $TmpArray[0]) * $TmpArray[1] - $TmpArray[2];
	if (-e "n_".$NbrParticles."/2s_".$SValue)
	  {
	    chdir ("n_".$NbrParticles."/2s_".$SValue);
	    my $TmpFileName;
	    foreach $TmpFileName (<*.vec>)
	      {
		if ($TmpFileName =~ /\_n_$NbrParticles\_2s\_$SValue\_lz\_0\.0\.vec$/)
		  {
		    my $Command = $PathToProgram."/QHEBosonsNBodyCorrelation -n ".$NbrNBody." -s ".$TmpFileName;		    
		    my $Correlation = `$Command`;
		    chomp ($Correlation);
		    print $NbrParticles." ".$FillingFactor." ".$Correlation."\n";
		  }
	      }
	    chdir ("../..");
	  }
	$NbrParticles += $TmpArray[3];	
      }
  }
