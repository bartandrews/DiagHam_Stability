#!/usr/bin/perl -w
# Batch for calculations of the combined V0-V1 interaction for fermions with spin
use strict 'vars';
use Math::Trig;

# global settings
my $DiagonalizationProgram="/home/gunnar/DiagHam/build1/FQHE/src/Programs/FQHEOnSphere/QHEFermionsSphereWithSpin";
my $memory = 0;
my $NumStates = 10;

if (!(defined($ARGV[4])))
  {
    print ("Run QHEFermionsSphereWithSpin for with given pseudopotentials V0,V1\n");
    die "usage: MapClockModel.pl nbr_fermions min_flux max_flux V0 V1  [GS] [Exc]\n";
  }

my $NbrFermions = $ARGV[0];
my $MinNbrFlux = $ARGV[1];
my $MaxNbrFlux = $ARGV[2];
my $V0 = $ARGV[3];
my $V1 = $ARGV[4];

if ( ($NbrFermions % 2) == 1 )
  {
    die ("Even particle number required");
  }

if  ( $MinNbrFlux < $NbrFermions/2 - 1)
  {
    $MinNbrFlux = $NbrFermions/2 - 1;
  }

my $OnlyGS="";
if (defined($ARGV[5]))
  {
   $OnlyGS=$ARGV[5];
  }
my $DoExc="";
if (defined($ARGV[6]))
  {
   $DoExc=$ARGV[6];
  }
my $Flux;
for ( $Flux=$MinNbrFlux; $Flux<=$MaxNbrFlux; $Flux++)
  {
    # run calculation
    my $Command = $DiagonalizationProgram." -p ".$NbrFermions." -l ".$Flux." -s 0 --show-itertime --memory ".$memory." -v".$V0." -w".$V1;
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
    print ("#running now: \n".$Command."\n");
    system($Command);
  }
