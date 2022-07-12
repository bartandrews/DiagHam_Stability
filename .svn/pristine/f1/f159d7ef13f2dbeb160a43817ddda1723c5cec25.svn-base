#!/usr/bin/perl -w

use strict 'vars';

my $NbrPoints = 51;
my $NbrParticles= 14;
my $SValue = 25;

my $Index = 0;
my $Step = 1.0 / ($NbrPoints - 1);
my $Value = 0.0;


my $ProgramPath = " /home/regnault/development/Physics/DiagHam/build64mpi/FQHE/src/Programs/FQHEOnSphere";
my $ReferenceState = "/scratch/regnault/nu_5_2_gap/n_14/fermions_mixed_coulomb_1_dv1_0.05_nbody_3_1_n_14_2s_25_lz_0.0.vec";
my $ReferenceState2 = "/scratch/regnault/nu_5_2_gap/n_14/fermions_mixed_coulomb_1_dv1_0.05_nbody_3_1_n_14_2s_25_lz_2.0.vec";

#my $ThreeBodyWeight = 0.306582478478;
#my $ThreeBodyWeight = 0.644340439407;
#my $ThreeBodyWeight = 0.4821806525841742527;
#my $ThreeBodyWeight = 0.4397786440410844119;
my $ThreeBodyWeight = 0.665260458701285005;
my $InteractionName = "mixed_coulomb_1_dv1_0.05_nbody_3";
#my $InteractionName = "mixed_coulomb_1_nbody_3";
#my $InteractionName = "mixed_v1_nbody_3";
#my $PseudopotentialFileName = "pseudopotential_v1_2s_15.dat";
my $PseudopotentialFileName = "pseudopotential_coulomb_l_1_dv1_0.05_2s_".($SValue - 2).".dat";
#my $PseudopotentialFileName = "pseudopotential_coulomb_l_1_2s_".($SValue - 2).".dat";


my $GapFileName = "fermions_".$InteractionName."_n_".$NbrParticles."_2s_".$SValue.".gap.dat";


my $DiagonalizationProgram = "mpirun -C ".$ProgramPath."/QHEFermionsNBodyHardCore --mpi";
my $OverlapProgram = $ProgramPath."/QHEBosonsDeltaOverlap";
my $LProgram = $ProgramPath."/FQHESphereLValue";

my @Pseudopotentials;

unless (open (INFILE, $PseudopotentialFileName))
  {
    die ("can't open ".$PseudopotentialFileName."\n");
  }
my $TmpLine;
while (defined($TmpLine = <INFILE>))
  {
    if ($TmpLine =~ /^\s*Pseudopotentials\s*\=/)
      {
	chomp($TmpLine);
	$TmpLine =~ s/^\s*Pseudopotentials\s*\=\s*//;
	$TmpLine =~ s/\s*$//;
	@Pseudopotentials = split (/\s+/, $TmpLine);
      }
  }
close (INFILE);

unless (open (OUTFILE, ">".$GapFileName))
  {
    die ("can't open ".$GapFileName."\n");
  }
close (OUTFILE);

my $Previous0 = "";
my $Previous2 = "";

while ($Index < $NbrPoints)
  {
    unless (open (OUTFILE, ">/home/regnault/tmp/pseudopotential_mixed_".$Value."_2s_".$SValue.".dat"))
      {
	die ("can't create file /home/regnault/tmp/pseudopotential_mixed_".$Value."_2s_".$SValue.".dat\n");
      }
    print OUTFILE "# pseudopotentials on the sphere for coulomb interaction 
# in the Landau level N=1 for 2S=".$SValue." flux quanta
#
# Pseudopotentials = V_0 V_1 ...

Pseudopotentials =";
    my $TmpPseudopotential;
    foreach $TmpPseudopotential (@Pseudopotentials)
      {
	print OUTFILE " ".((1.0 - $Value) * $TmpPseudopotential);
      }
    print OUTFILE "\n\nNbrNBody=3

Weights=0 0 0 ".($ThreeBodyWeight * $Value)."\n";
    close (OUTFILE);

    my $OutputName = "fermions_".$InteractionName."_".$Value."_n_".$NbrParticles."_2s_".$SValue."_lz";
    
    if ((!(-e $OutputName."_0.0.vec")) || (!(-e $OutputName.".0.dat")))
      {
	my $Command = $DiagonalizationProgram." -p ".$NbrParticles." -l ".$SValue."  --eigenstate --memory 1500 --force-reorthogonalize --nbr-lz 1 -n 1 --nbody-file /home/regnault/tmp/pseudopotential_mixed_".$Value."_2s_".$SValue.".dat";
	if ($Previous0 ne "")
	  {
	    $Command .= " --initial-vector ".$Previous0;
	  }
	system ($Command);
	
	rename ("fermions_hardcore_nbody_3_n_".$NbrParticles."_2s_".$SValue."_lz_0.0.vec", $OutputName."_0.0.vec");    
	rename ("fermions_hardcore_nbody_3_n_".$NbrParticles."_2s_".$SValue."_lz.dat", $OutputName.".0.dat");
	$Previous0 = "/home/regnault/tmp/tmp0.vec";
	$Command = "cp ".$OutputName."_0.0.vec ".$Previous0;
	system ($Command);
      }
	
    if ((!(-e $OutputName."_2.0.vec")) || (!(-e $OutputName.".2.dat")))
      {
	my $Command = $DiagonalizationProgram." -p ".$NbrParticles." -l ".$SValue."  --eigenstate --memory 1500 --force-reorthogonalize --nbr-lz 1 --initial-lz 2 -n 1 --nbody-file /home/regnault/tmp/pseudopotential_mixed_".$Value."_2s_".$SValue.".dat";
	if ($Previous2 ne "")
	  {
	    $Command .= " --initial-vector ".$Previous2;
	  }
	
	system ($Command);
	rename ("fermions_hardcore_nbody_3_n_".$NbrParticles."_2s_".$SValue."_lz_2.0.vec", $OutputName."_2.0.vec");    
	rename ("fermions_hardcore_nbody_3_n_".$NbrParticles."_2s_".$SValue."_lz.dat", $OutputName.".2.dat");
	$Previous2 = "/home/regnault/tmp/tmp2.vec";
	$Command = "cp ".$OutputName."_2.0.vec ".$Previous2;
	system ($Command);
      }
    
    unless (open (INFILE, $OutputName.".0.dat"))
      {
	die ("can't open ".$OutputName.".0.dat\n");
      }
    $TmpLine = <INFILE>;
    chomp ($TmpLine);
    my @TmpArray = split (/ /, $TmpLine);
    my $Gap = -$TmpArray[1];
    close (INFILE);
    unless (open (INFILE, $OutputName.".2.dat"))
      {
	die ("can't open ".$OutputName.".2.dat\n");
      }
    $TmpLine = <INFILE>;
    chomp ($TmpLine);
    @TmpArray = split (/ /, $TmpLine);
    $Gap += $TmpArray[1];
    close (INFILE);

    unless (open (OUTFILE, ">>".$GapFileName))
      {
	die ("can't open ".$GapFileName."\n");
      }
    print OUTFILE $Value." ".$Gap;
    
    if ($ReferenceState ne "")
      {
	my $Command = $OverlapProgram." --use-exact ".$OutputName."_0.0.vec --exact-state ".$ReferenceState;
	my $Overlap = `$Command`;
	chomp ($Overlap);
	$Overlap =~ s/^\s*overlap\s*=\s*//; 
	print OUTFILE " ".($Overlap * $Overlap);
      }	
    if ($ReferenceState2 ne "")
      {
	my $Command = $OverlapProgram." --use-exact ".$OutputName."_2.0.vec --exact-state ".$ReferenceState2;
	my $Overlap = `$Command`;
	chomp ($Overlap);
	$Overlap =~ s/^\s*overlap\s*=\s*//;
	print OUTFILE " ".($Overlap * $Overlap);
	$Command = $LProgram." ".$OutputName."_2.0.vec";
	$Overlap = `$Command`;
	my @TmpArray = split (/\n/, $Overlap);
	my $TmpLine;
	foreach $TmpLine (@TmpArray)
	  {
	    if ($TmpLine =~ /^\s*\<L\>\s*\=/)
	      {
		$TmpLine =~ s/^\s*\<L\>\s*\=\s*//;
		print OUTFILE " ".$TmpLine;
	      }
	  }
      }
    
    
    print OUTFILE "\n";
    close (OUTFILE);
    
    $Value += $Step;
    $Index++;
  }
