#!/usr/bin/perl -w

use strict 'vars';

use Getopt::Long;


my $PathToLzToLProgram = "/home/regnault/development/Physics/DiagHam/build/FQHE/src/Programs/FQHEOnSphere/LzToL";
my $CheckValidity = 0;
my $ValidtyError = 1e-12;
my $DegeneracyError = 1e-12;
my $FilePattern;

my $Result = GetOptions ("progdiag:s" => \$PathToLzToLProgram, 
			 "file-pattern=s" => \$FilePattern,
			 "check-validity" => \$CheckValidity,
			 "validity-error:s" => \$ValidtyError,
			 "degeneracy-error:s" => \$DegeneracyError);



$FilePattern =~ /\_n\_(\d+)\_/;
my $NbrParticles = $1;
if (!(defined($NbrParticles)))
  {
    die ("can't find number of particles in file pattern ".$FilePattern."\n");
  }
$FilePattern =~ /^(.+)\_sz\_\d+\_iz\_\d+\_pz\_\d+\_lz\.dat$/;
my $Prefix = $1;

my $TmpFile;
my @Sprectra;
foreach $TmpFile (<*>)
  {
    if ($TmpFile =~ /^$Prefix\_sz\_(\d+)\_iz\_(\d+)\_pz\_(\d+)\_lz\.dat$/)
      {
	if (!(-f $Prefix."_sz_".$1."_iz_".$2."_pz_".$3."_l.dat"))
	  {
	    my $Command = $PathToLzToLProgram." -f ".$TmpFile." > ".$Prefix."_sz_".$1."_iz_".$2."_pz_".$3."_l.dat";
	    `$Command`;
	    print "converting Lz->L ".$TmpFile."\n";
	  }
	push (@Sprectra, $Prefix."_sz_".$1."_iz_".$2."_pz_".$3."_l.dat");
      }
  }

my $MinFile = shift(@Sprectra);
my $MinEnergy;
my $MinL;
&FindMinEnergy ($MinFile, \$MinEnergy, \$MinL);
my $Degeneracy = &FindGroundStateDegeneracy($MinFile, $MinEnergy, $MinL, $DegeneracyError);
print "local minimum energy in file ".$MinFile." at L=".$MinL." (".$MinEnergy.") with degeneracy ".$Degeneracy."\n";

foreach $TmpFile (@Sprectra)
  {
    my $CurrentMinEnergy ;
    my $CurrentMinL;    
    &FindMinEnergy($TmpFile, \$CurrentMinEnergy, \$CurrentMinL);
    $Degeneracy = &FindGroundStateDegeneracy($TmpFile, $CurrentMinEnergy, $CurrentMinL, $DegeneracyError);
    print "local minimum energy in file ".$TmpFile." at L=".$CurrentMinL." (".$CurrentMinEnergy.") with degeneracy ".$Degeneracy."\n";
    if ($CurrentMinEnergy < $MinEnergy)
      {
	$MinEnergy = $CurrentMinEnergy;
	$MinL = $CurrentMinL;
	$MinFile = $TmpFile;
      }
  }

print "Minimum energy found in file ".$MinFile." at L=".$MinL." (".$MinEnergy.")\n";



sub FindMinEnergy ()
  {
    my $FileName = $_[0];
    my $MinEnergy = $_[1];
    my $CurrentMinL = $_[2];
    unless (open(INFILE ,$FileName))
      {
	die ("can't open ".$FileName."\n");	
      }
    my $TmpLine = <INFILE>;
    chomp($TmpLine);    
    my @TmpArray = split (/ /, $TmpLine);
    $$MinEnergy = $TmpArray[1];
    $$CurrentMinL = $TmpArray[0];
    while (defined($TmpLine = <INFILE>))
      {
	chomp($TmpLine);    
	@TmpArray = split (/ /, $TmpLine);
	if ($$MinEnergy > $TmpArray[1])
	  {
	    $$MinEnergy = $TmpArray[1];
	    $$CurrentMinL = $TmpArray[0];
	  }
      }
    close (INFILE);
  }

sub FindGroundStateDegeneracy ()
  {
    my $FileName = $_[0];
    my $GroundStateEnergy = $_[1];
    my $GroundStateL = $_[2];
    my $DegeneracyError = $_[3];
    my $Degeneracy = 0;
    unless (open(INFILE ,$FileName))
      {
	die ("can't open ".$FileName."\n");	
      }
    my $TmpLine;
    while (defined($TmpLine = <INFILE>))
      {
	chomp($TmpLine);    
	my @TmpArray = split (/ /, $TmpLine);
	if ($TmpArray[0] == $GroundStateL)
	  {
	    if (((abs($TmpArray[1]) < $DegeneracyError) && (abs($GroundStateEnergy) < $DegeneracyError)) ||
		(abs($TmpArray[1] - $GroundStateEnergy) < ($DegeneracyError * abs($GroundStateEnergy))))
	      {
		$Degeneracy++;
	      }
	  }
      }
    close (INFILE);
    return $Degeneracy;
  }
