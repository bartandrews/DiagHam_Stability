#!/usr/bin/perl -w
#
# script for composing a bilayer pseudopotential file
#
use strict 'vars';
use File::stat;

my $Program="GenerateEffectiveInteraction";
my $MCProgram="FQHESphereMonteCarlo";

my $SimulateMC=0;
my $NbrParticles=0;
my $NbrIterations=10000000;
my $Exponent=2;
my $BackgroundStr="";
my $RefineFlag="";

while( (defined($ARGV[0])) && ( $ARGV[0] =~ /^-/ ))
  {
    if ( $ARGV[0] =~ /^--/ )
      {
	if ( $ARGV[0] =~ /--refine/ )
	  {
	    $RefineFlag="-r";
	  }
      }
    if ( $ARGV[0] =~ /-r/ )
      {
	$SimulateMC=1;
      }
    if ( $ARGV[0] =~ /-b/ )
      {
	$BackgroundStr="-b";
      }
    if ( $ARGV[0] =~ /-p/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $NbrParticles = $ARGV[0];
	    $NbrParticles =~ s/-p//;
	  }
	else
	  {
	    shift(@ARGV);
	    $NbrParticles = $ARGV[0];
	  }
      }
    if ( $ARGV[0] =~ /-i/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $NbrIterations = $ARGV[0];
	    $NbrIterations =~ s/-i//;
	  }
	else
	  {
	    shift(@ARGV);
	    $NbrIterations = $ARGV[0];
	  }
      }
    if ( $ARGV[0] =~ /-e/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $Exponent = $ARGV[0];
	    $Exponent =~ s/-e//;
	  }
	else
	  {
	    shift(@ARGV);
	    $Exponent = $ARGV[0];
	  }
      }
    shift(@ARGV);
  }

if (!defined($ARGV[1]))
  {
    print("usage MakeBilayerInteraction.pl [-r -p nbrParticle [-e exponent(=2)] [-i iterations]] flux pseudopotentials [datafile]\n");
    print("options\n  -r: run Monte-Carlo simulation using the generated potentials\n");
    print("  -p number of particles (required for simulation)\n");
    print("  -i number of MC iterations (default: 10000000)\n");
    print("  -e exponent of the Jastrow-factor in the Pfaffian state \n");
    print("  -b print background energy only, then exit\n");
    print("  if datafile is given, a line with the state energy will be added, otherwise, write to MCSummary.dat\n");
    exit(1);
  }

if (($SimulateMC>0)&&($NbrParticles==0))
  {
    print("To run a simulation, please indicate the number of particles with -p\n");
    exit(1);
  }

my $Flux = $ARGV[0];
my $PseudoPotFile = $ARGV[1];
my $Datafile = "MCSummarySL.dat";
if (length($BackgroundStr)>0)
  {
    $Datafile = "MCSummarySL-b.dat";
  }

if (defined($ARGV[2]))
  {
    $Datafile=$ARGV[2];
  }

my $PseudoPotBase = $PseudoPotFile;
$PseudoPotBase =~ s/.dat//;

print("Using BaseName = $PseudoPotBase\n");

system("$Program -s $Flux -i $PseudoPotFile -o $PseudoPotBase.int $RefineFlag");

my $Options="--interaction-params $PseudoPotBase.int";
print("Options: $Options\n");

if ($SimulateMC>0)
  {
    my $TmpFileName = "/tmp/tmp".time().".p";
    system("$MCProgram -p $NbrParticles -l $Flux -i $NbrIterations --test-wavefunction pairedcf --pair-coeff 0 --nbr-flux $Exponent --sampler laughlin --laughlin-exponent $Exponent $Options -o $TmpFileName $BackgroundStr");
    open (INFILE, $TmpFileName);
    open (OUTFILE, ">>$Datafile");
    my $TmpLine;
    foreach $TmpLine (<INFILE>)
      {
	if (defined($TmpLine) && (! ( $TmpLine =~ m/#/ )))
	  {
	    chomp($TmpLine);
	    print ($PseudoPotBase." ".$Exponent." ".$TmpLine."\n");
	    print OUTFILE ($PseudoPotBase."\t".$Exponent."\t".$TmpLine."\n");
	  }
      }
    close(OUTFILE);
    close(INFILE);
    system("rm $TmpFileName");
  }
