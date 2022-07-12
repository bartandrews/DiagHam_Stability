#!/usr/bin/perl -w

use strict 'vars';

my $Count = 0;;

my $NbrFermions = 10;
my $SValue = 18;
my $LandauLevel = 1;
my $ScaleFactor = 1.0 / $SValue;# / sqrt($SValue);
my $NbrNBody = 3;

#my $InteractionName = "coulomb_".$LandauLevel;
my $InteractionName = "hardcore_nbody_3";

my $PathToProgram = "/home/regnault/development/Physics/DiagHam/build/src/Programs/QHE/QHEOnSphere";
#my $DiagonalizationProgram = "QHEFermionsTwoBodyGeneric";
my $DiagonalizationProgram = "QHEFermionsNBodyHardCore";
my $DiagonalizationOptions = "--use-lapack --full-diag 700";
my $CorrelationProgram = "QHEFermionsCorrelation";
my $PseudopotentialProgram = "CoulombPseudopotentials";

#my @Values = (-0.001, -0.1, -0.3, -0.6);
my @Values = (0, 0.05, 0.1, 0.15);
#my @Values = (0, -0.01, -0.02, -0.04, -0.08, -0.16);
#my @Values = (0, 0.01, 0.02, 0.03,0.04, 0.05);
#my @Values = (0.15, 0.05, 0.035, 0.0175, 0, -0.0175, -0.035, -0.05, -0.15);

my @GnuplotLines;
my $TmpValues;
my $Index = 0;
foreach $TmpValues (@Values)
  {
    if ($NbrNBody == 2)
      {
	system($PathToProgram."/".$PseudopotentialProgram." -l ".$LandauLevel." -s ".($SValue - (2 * $LandauLevel))." --add-impurities --north-potential ".($TmpValues * $ScaleFactor)." --south-potential ".($TmpValues * $ScaleFactor));
	system($PathToProgram."/".$DiagonalizationProgram." -p ".$NbrFermions." -l ".$SValue." -n 2 --eigenstate --interaction-name ".$InteractionName." --interaction-file pseudopotential_coulomb_l_".$LandauLevel."_2s_".($SValue - (2 * $LandauLevel)).".dat ".$DiagonalizationOptions);
      }
    else
      {
	system($PathToProgram."/".$DiagonalizationProgram." -p ".$NbrFermions." -l ".$SValue." -n 1 --eigenstate --nbr-nbody ".$NbrNBody."  --add-impurities --impurity-potential ".($TmpValues * $ScaleFactor)." --landau-level ".$LandauLevel." ".$DiagonalizationOptions);
	print $PathToProgram."/".$DiagonalizationProgram." -p ".$NbrFermions." -l ".$SValue." -n 1 --eigenstate --nbr-nbody ".$NbrNBody."  --add-impurities --impurity-potential ".($TmpValues * $ScaleFactor)." --landau-level ".$LandauLevel." ".$DiagonalizationOptions."\n";
      }
    my $LzValue = &FindLzGround("fermions_".$InteractionName."_n_".$NbrFermions."_2s_".$SValue."_lz.dat");
    system($PathToProgram."/".$CorrelationProgram." -z ".(2 * $LzValue)." --density -p ".$NbrFermions." -l ".$SValue." -i ".$InteractionName." -r -s fermions_".$InteractionName."_n_".$NbrFermions."_2s_".$SValue."_lz_".(2 * $LzValue).".0.vec --landau-level ".$LandauLevel);
    rename("fermions_".$InteractionName."_n_".$NbrFermions."_2s_".$SValue.".rho.dat", "fermions_".$InteractionName."_n_".$NbrFermions."_2s_".$SValue.".rho.".$TmpValues.".dat");
    push (@GnuplotLines, "\"fermions_".$InteractionName."_n_".$NbrFermions."_2s_".$SValue.".rho.".$TmpValues.".dat\" using 1:2 title \"V=".$TmpValues." (L_z=".$LzValue.")\" with lines lt ".$Index);
    $Index++;
  }

open (OUTFILE, ">tmp.p");
#print OUTFILE "set yrange [0.5:0.9]
print OUTFILE "set xlabel \"{/Symbol q}\"
set ylabel \"{/Symbol r}\"
set size 1, 0.5
set terminal postscript landscape color solid enhanced \"Helvetica\" 14
set output \"fermions_".$InteractionName."_n_".$NbrFermions."_2s_".$SValue."_vp.rho.ps\"
plot ".join(", ", @GnuplotLines)."\n";
close (OUTFILE);

system ("gnuplot tmp.p");


# Find the lz value at which the ground state appears
#
# $_[0] = name of the file that contains the spectrum
# return value = Lz value (minus 1/2 if real Lz value is an half integer)

sub FindLzGround ()
  {
    my $FileName = $_[0];
    my $Min;
    my $MinLz;
    my $Flag = 0;
    unless (open (INFILE, $FileName))
      {
	die ("can't open ".$FileName."\n");
      }
    my $TmpLine;
    foreach $TmpLine (<INFILE>)
      {
	chomp ($TmpLine);
	my @TmpArray = split (/ /, $TmpLine);
	if ($Flag == 0)
	  {
	    $Min = $TmpArray[1];
	    $MinLz = $TmpArray[0];
	    $Flag = 1;
	  }
	else
	  {
	    if ($TmpArray[1] < $Min)
	      {
		$Min = $TmpArray[1];
		$MinLz = $TmpArray[0];
	      }
	  }
      }
    close (INFILE);
    return $MinLz;
  }
