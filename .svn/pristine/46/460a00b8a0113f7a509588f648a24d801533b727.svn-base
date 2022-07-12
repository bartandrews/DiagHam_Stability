#!/usr/bin/perl -w
# Batch for calculations of the combined V0-V1 interaction for fermions with spin
use strict 'vars';
use Math::Trig;

# global settings
my $DiagonalizationProgram="/home/gunnar/DiagHam/build1/FQHE/src/Programs/FQHEOnSphere/FQHESphereFermionsWithSpin";
my $memory = 0;
my $NumStates = 10;

if (!(defined($ARGV[1])))
  {
    print ("Run FQHESphereFermionsWithSpin with pseudopotentials V1=1\n");
    die "usage: UnitPseudoPot.pl nbr_fermions nbr_flux [GS] [Exc]\n";
  }

my $NbrFermions = $ARGV[0];
my $NbrFlux = $ARGV[1];

my $OnlyGS="";
if (defined($ARGV[2]))
  {
   $OnlyGS=$ARGV[2];
  }
my $DoExc="";
if (defined($ARGV[3]))
  {
   $DoExc=$ARGV[3];
  }

my $PseudopotentialFile = "pseudopotentials_UnitV1_2s_".$NbrFlux.".dat";
unless (open (OUTFILE, ">".$PseudopotentialFile))
  {
    die ("can't create pseudopotential file\n");
  }
print OUTFILE ("# Pseudopotentials for unit V1 potential for 2s=".$NbrFlux."\n");
print OUTFILE ("Pseudopotentials =0 1");
for ( my $i=2; $i<=$NbrFlux; $i++)
  {
    print OUTFILE (" 0");
  }
close (OUTFILE);
# run calculation
my $Command = $DiagonalizationProgram." -p ".$NbrFermions." -l ".$NbrFlux." -s 0 --show-itertime --memory ".$memory." --interaction-file ".$PseudopotentialFile." --interaction-name UnitV1";
my $Command2=$Command;
#$Command = $Command." --szsymmetrized-basis";
if (($OnlyGS =~ "GS") && !($DoExc =~ "Exc"))
  {
    $Command = $Command." --nbr-lz 1 --eigenstate -n 1 --force-reorthogonalize";
  }
else
  {
    if ($DoExc =~ "Exc")
      {
	$Command = $Command." --nbr-lz 2 --eigenstate -n 1 --force-reorthogonalize";
      }
    else
      {
	$Command = $Command."  -n".$NumStates;
      }
  }
print ("#To run, type: \n".$Command."\n");

