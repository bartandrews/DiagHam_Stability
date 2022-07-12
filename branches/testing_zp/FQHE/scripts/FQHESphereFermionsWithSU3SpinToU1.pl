#!/usr/bin/perl -w

use strict 'vars';
use POSIX;

my $Error = 1e-10;
my $SpectrumFileName = $ARGV[0];
my $Lz = 0;
my $L2DiagonalizationProgram = "";
my $SUNConvertProgram = "";
my $ProgramPath = "/home/regnault/development/Physics/DiagHam/build/";
my $SUNToU1Program = $ProgramPath."FQHE/src/Programs/FQHEOnSphere/FQHESphereSUKToU1";
my $SU3Flag = 0;
my $SU2Flag = 0;
my $SU4Flag = 0;
my $InteractionName = 0;
my $NbrParticles = 0;
my $NbrFluxQuanta = 0;
my $InputSymmetrizedL2VectorPrefix = "" ;
my $InputUnsymmetrizedL2VectorPrefix = "";
if ($SpectrumFileName =~ /\_su3\_/)
  {
    $SU3Flag = 1;
  }
elsif ($SpectrumFileName =~ /\_su2\_/)
  {
    $SU2Flag = 1;
  }
elsif($SpectrumFileName =~ /\_su4\_/)
  {
    $SU4Flag = 1;
  }
  else
  {
    die ("unsuppored SU(N) symmetry\n");
  }

if ($SU3Flag == 1)
  {
    $L2DiagonalizationProgram = $ProgramPath."FQHE/src/Programs/FQHEOnSphere/FQHESphereWithSU3SpinL2Diagonalize";
    $SUNConvertProgram =$ProgramPath."FQHE/src/Programs/FQHEOnSphere/FQHESphereWithSU3SpinConvertSymmetrizedState";
    $SpectrumFileName =~ /^fermions\_sphere\_su3\_tzpz3sym_(.*)\_n\_(\d+)\_2s\_(\d+)\_tz\_(\d+)\_y\_(\d+)\_lz\.dat$/;
    $InteractionName = $1;
    $NbrParticles = $2;
    $NbrFluxQuanta = $3;
    my $TotalTz = $4;
    my $TotalY = $5;
    $L2DiagonalizationProgram .= " -p ".$NbrParticles." -l ".$NbrFluxQuanta." -z ".$Lz." -y ".$TotalY." -t ".$TotalTz." --z3symmetrized-basis --tzsymmetrized-basis --use-hilbert basis.tmp --eigenstate --interaction-name tzpz3sym_".$InteractionName."_l2";
    $InputSymmetrizedL2VectorPrefix = "fermions_sphere_su3_tzpz3sym_".$InteractionName."_l2_n_".$NbrParticles."_2s_".$NbrFluxQuanta."_tz_".$TotalTz."_y_".$TotalY."_lz";
    $InputUnsymmetrizedL2VectorPrefix = "fermions_sphere_su3_".$InteractionName."_l2_n_".$NbrParticles."_2s_".$NbrFluxQuanta."_tz_".$TotalTz."_y_".$TotalY."_lz";
  }
elsif ($SU2Flag == 1)
  {
    $L2DiagonalizationProgram = $ProgramPath."FQHE/src/Programs/FQHEOnSphere/FQHESphereWithSpinL2Diagonalize";
    $SUNConvertProgram =$ProgramPath."FQHE/src/Programs/FQHEOnSphere/FQHESphereWithSpinConvertSymmetrizedState";
    $SpectrumFileName =~ /^fermions\_sphere\_su2\_szpsym_(.*)\_n\_(\d+)\_2s\_(\d+)\_sz\_(\d+)\_lz\.dat$/;
    $InteractionName = $1;
    $NbrParticles = $2;
    $NbrFluxQuanta = $3;
    my $TotalSz = $4;
    $L2DiagonalizationProgram .= " -p ".$NbrParticles." -l ".$NbrFluxQuanta." -z ".$Lz." -s ".$TotalSz." --szsymmetrized-basis --use-hilbert basis.tmp --eigenstate --interaction-name szpsym_".$InteractionName."_l2";
    $InputSymmetrizedL2VectorPrefix = "fermions_sphere_su2_szpsym_".$InteractionName."_l2_n_".$NbrParticles."_2s_".$NbrFluxQuanta."_sz_".$TotalSz."_lz";
    $InputUnsymmetrizedL2VectorPrefix = "fermions_sphere_su2_".$InteractionName."_l2_n_".$NbrParticles."_2s_".$NbrFluxQuanta."_sz_".$TotalSz."_lz";
  }
else
  {
    $L2DiagonalizationProgram = $ProgramPath."FQHE/src/Programs/FQHEOnSphere/FQHESphereWithSU4SpinL2Diagonalize";
    $SUNConvertProgram =$ProgramPath."FQHE/src/Programs/FQHEOnSphere/FQHESphereWithSU4SpinConvertSymmetrizedState";
  }
my $LinearlyIndependentProgram = $ProgramPath."src/Programs/ExtractLinearlyIndependentVectors";

my $NbrZeroEnergyStates = &ExtractNbrZeroEnergyStates($SpectrumFileName, $Lz, $Error);
$SpectrumFileName =~ /^fermions\_sphere\_su3\_tzmz3sym_(.*)\_n\_(\d+)\_2s\_(\d+)\_tz\_(\d+)\_y\_(\d+)\_lz\.dat$/;
my $InputVectorPrefix = $SpectrumFileName;
$InputVectorPrefix =~ s/\.dat$//;

unless (open (OUTFILE, ">basis.tmp"))
  {
    die ("can't open basis.tmp\n");
  }
print OUTFILE "Basis =";
my $Count = 0;
while ($Count < $NbrZeroEnergyStates)
  {
    print OUTFILE " ".$InputVectorPrefix."_".$Lz.".".$Count.".vec";
    $Count++;
  }
print OUTFILE "\n";
close(OUTFILE);

my $Command = $L2DiagonalizationProgram;
system($Command);
unlink ("basis.tmp");

$Count = $NbrZeroEnergyStates - 1;
while ($Count >= 0)
  {
    my $TmpVectorName = $InputSymmetrizedL2VectorPrefix."_".$Lz.".".$Count.".vec";
    my $TmpVectorName2 = $InputUnsymmetrizedL2VectorPrefix."_".$Lz.".".$Count.".vec";
    $Command = $SUNConvertProgram." ".$TmpVectorName." -o ".$TmpVectorName2;
    system($Command);
    $Command = $SUNToU1Program." -s ".$TmpVectorName2." --interaction-name ".$InteractionName;
    system($Command);
    if ($Count != 0)
      {
	rename("fermions_sphere_".$InteractionName."_n_".$NbrParticles."_2s_".$NbrFluxQuanta."_lz_".$Lz.".0.vec",
	       "fermions_sphere_".$InteractionName."_n_".$NbrParticles."_2s_".$NbrFluxQuanta."_lz_".$Lz.".".$Count.".vec");
      }
    unlink ($TmpVectorName);
    $Count--;
  }

my @SUNLDegeneracy;
my @U1LDegeneracy;
&ExtractLDegeneracy($InputSymmetrizedL2VectorPrefix.".dat", \@SUNLDegeneracy);
my $MaxL = $#SUNLDegeneracy;
my $TmpL = 0;
my $Shift = 0;
while ($TmpL <= $MaxL)
  {
    if ($SUNLDegeneracy[$TmpL] > 1)
      {
	unless (open (OUTFILE, ">basis.tmp"))
	  {
	    die ("can't open basis.tmp\n");
	  }
	print OUTFILE "Basis =";
	$Count = 0;
	while ($Count < $SUNLDegeneracy[$TmpL])
	  {
	    print OUTFILE " fermions_sphere_".$InteractionName."_n_".$NbrParticles."_2s_".$NbrFluxQuanta."_lz_".$Lz.".".($Shift + $Count).".vec";
	    $Count++;
	  }
	print OUTFILE "\n";
	close(OUTFILE);
	$Command = $LinearlyIndependentProgram." -b basis.tmp ";
	my $Tmp = `$Command`;
	$Tmp =~ /(\d+) linearly independent vectors/gm;
	$U1LDegeneracy[$TmpL] = $1;
	$Count = 0;
	print $Tmp;
	while ($Count < $SUNLDegeneracy[$TmpL])
	  {
	    unlink ("fermions_sphere_".$InteractionName."_n_".$NbrParticles."_2s_".$NbrFluxQuanta."_lz_".$Lz.".".($Shift + $Count).".vec");
	    $Count++;
	  }
	$Count = 0;
	while ($Count < $U1LDegeneracy[$TmpL])
	  {
	    rename("vector_".$Count.".vec", "fermions_sphere_".$InteractionName."_n_".$NbrParticles."_2s_".$NbrFluxQuanta."_l_".$TmpL.".".$Count.".vec");
	    $Count++;
	  }
	unlink ("basis.tmp");
	$Shift += $SUNLDegeneracy[$TmpL];
      }
    else
      {
	if ($SUNLDegeneracy[$TmpL] == 1)
	  {
	    $U1LDegeneracy[$TmpL] = 1;
	    rename("fermions_sphere_".$InteractionName."_n_".$NbrParticles."_2s_".$NbrFluxQuanta."_lz_".$Lz.".".$Shift.".vec", 
		   "fermions_sphere_".$InteractionName."_n_".$NbrParticles."_2s_".$NbrFluxQuanta."_l_".$TmpL.".0.vec");
	    $Shift++;
	  }
	else
	  {
	    $U1LDegeneracy[$TmpL] = 0;
	  }
      }
    $TmpL++;
  }
print join(" ", @SUNLDegeneracy)."\n";
print join(" ", @U1LDegeneracy)."\n";

# return the number of zero energy states in a given Lz sector
#
# $_[0] = spectrum file name
# $_[1] = twice the Lz value
# $_[2] = error on the definition of a zero energy state
# return value = number of zero energy states

sub ExtractNbrZeroEnergyStates ()
  {
    my $Spectrum = $_[0];
    my $Lz = $_[1] / 2;    
    my $Error = $_[2];
    my $Count = 0;
    unless (open (INFILE, $Spectrum))
      {
	die ("can't open ".$Spectrum."\n");
      }
    my $TmpLine;
    while (defined($TmpLine = <INFILE>))
      {
	chomp($TmpLine);
	$TmpLine =~ s/\#.*$//;
	$TmpLine =~ s/^\s+//;
	if ($TmpLine ne "")
	  {
	    my @TmpArray = split (/\s+/, $TmpLine);
	    if (($TmpArray[0] == $Lz) && (abs($TmpArray[1]) < $Error))
	      {
		$Count++;
	      }
	    if ($TmpArray[0] > $Lz)
	      {
		last;
	      }
	  }
      }
    close (INFILE);
    return $Count;
  }

# extract the L degeneracy decomposition 
#
# $_[0] = L2 spectrum file name
# $_[1] = reference on the array where the L degeneracy decomposition has to be stored

sub ExtractLDegeneracy()
  {
    my $Spectrum = $_[0];
    my $Degeneracy = $_[1];
    unless (open (INFILE, $Spectrum))
      {
	die ("can't open ".$Spectrum."\n");
      }
    my $TmpLine;
    my %TmpDegeneracy;
    my $MaxL = 0;
    while (defined($TmpLine = <INFILE>))
      {
	chomp($TmpLine);
	$TmpLine =~ s/\#.*$//;
	$TmpLine =~ s/^\s+//;
	if ($TmpLine ne "")
	  {
	    my @TmpArray = split (/\s+/, $TmpLine);
	    my $LValue = 0;
	    if ($TmpArray[1] > 1e-9)
	      {
		$LValue = POSIX::ceil(0.5 * (sqrt((4.0 * $TmpArray[1]) + 1.0) - 1.0) - 0.1);
		if ($MaxL < $LValue)
		  {
		    $MaxL = $LValue;
		  }
	      }
	    if (defined($TmpDegeneracy{$LValue}))
	      {
		$TmpDegeneracy{$LValue}++;
	      }
	    else
	      {
		$TmpDegeneracy{$LValue} = 1;
	      }
	  }
      }
    close (INFILE);
    while ($MaxL >= 0)
      {
	if (defined($TmpDegeneracy{$MaxL}))
	  {
	    $$Degeneracy[$MaxL] = $TmpDegeneracy{$MaxL};
	  }
	else
	  {
	     $$Degeneracy[$MaxL] = 0;
	  }
	$MaxL--;
      }
  }
