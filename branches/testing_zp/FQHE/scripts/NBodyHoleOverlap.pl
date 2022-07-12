#!/usr/bin/perl -w

use strict 'vars';

use Getopt::Long;


my $Error = 1.0e-12;
my $KValue = 0;
my $NbrHoles = 0;
my $RowFlag = 0;
my $SumFlag = 0;
my $LatexFlag = 0;
my $DebugFlag = 0;
my $FermionFlag = 0;
my $Precision = 4;
my $MaxNbrParticles = 100;
my $DataPath = "/home/regnault/public_html/qhe/QHEOnSphere";
my $DiagHamPath = "/home/regnault/development/DMRG/DiagHam";
my $TrueIntercationName = "delta";
my $Result = GetOptions ("kvalue=s" => \$KValue, "nbrholes=s" => \$NbrHoles, "error:s" => \$Error, 
			 "datapath:s" => \$DataPath, "diaghampath:s" => \$DiagHamPath, "latex" => \$LatexFlag, "precision" => \$Precision,
			 "maxparticles:s" => \$MaxNbrParticles, "debug" => \$DebugFlag, "fermion" => \$FermionFlag, "interaction:s" => \$TrueIntercationName); 

if ((!($KValue =~ /^\d+$/)) || ($KValue == 0) || (!($NbrHoles =~ /^\d+$/)) || ($NbrHoles == 0) || (!($Precision =~ /^\d+$/)) || ($Precision <= 2) ||
    (!($MaxNbrParticles =~ /^\d+$/)) || ($MaxNbrParticles <= 4))
  {
    die ("usage: NBodyHoleOverlap.pl --kvalue k --nbrholes h [--datapath --diaghampath --error 1e-12 --latex --precision 4 --maxparticles 100 --debug --fermion --interaction delta]\n");
  }

my $QHENBodyQuasiHoleOverlap = $DiagHamPath."/build/src/Programs/QHE/QHEOnSphere/QHENBodyQuasiHoleOverlap";
my $SphereSpectrumDegeneracy = $DiagHamPath."/scripts/QHE/SphereSpectrumDegeneracy.pl";
if (!(-x $QHENBodyQuasiHoleOverlap))
  {
    die ("can't find or execute ".$QHENBodyQuasiHoleOverlap."\n");
  }
if (!(-x $SphereSpectrumDegeneracy))
  {
    die ("can't find or execute ".$SphereSpectrumDegeneracy."\n");
  }

my $NbrParticles = 2 * $KValue;
my $SValue = 2 + $NbrHoles;
my $Statistics = "bosons";
if ($FermionFlag == 1)
  {
     $Statistics = "fermions";
     $NbrParticles = 2 * $KValue;
     $SValue = (2 * $KValue) + 1 + $NbrHoles;
  }
$DataPath .= "/".$Statistics;
my $FinalOutput = "";
my $MaxLz = 0;

while ($NbrParticles < $MaxNbrParticles)
  {
    my $HardcorePrefix = $DataPath."/hardcore_".($KValue + 1)."/n_".$NbrParticles."/2s_".$SValue."/".$Statistics."_hardcore_nbody_".($KValue + 1)."_n_".$NbrParticles."_2s_".$SValue."_lz";
    my $SpectrumName = $HardcorePrefix.".dat";
    print $SpectrumName."\n";
    if (-e $SpectrumName)
      {
	my $SpectrumCommand = $SphereSpectrumDegeneracy." -spectrum ".$SpectrumName." --row --error ".$Error;
	my $SpectrumOutput = `$SpectrumCommand`;
	chomp ($SpectrumOutput);
	my $DeltaPrefix = $DataPath."/".$TrueIntercationName."/n_".$NbrParticles."/2s_".$SValue."/".$Statistics."_".$TrueIntercationName."_n_".$NbrParticles."_2s_".$SValue."_lz";
	my $TemporaryOverlapDefinition = "NbrParticles = ".$NbrParticles."
Statistics = ".$Statistics."
LzMax = ".$SValue."
Degeneracy = ".$SpectrumOutput."
InputVectors = ".$DeltaPrefix."_
OutputVectors = ".$HardcorePrefix."_
Spectrum = ".$DeltaPrefix.".dat\n";
	my $TemporaryFileName = &GetTemporaryFileName().".def";
	unless (open (OUTFILE, ">".$TemporaryFileName))
	  {
	    die "can't create file ".$TemporaryFileName."\n";
	  }
	print OUTFILE $TemporaryOverlapDefinition;
	close (OUTFILE);
	if ($DebugFlag == 1)
	  {
	    print "===============================================================\nNew input file:\n".$TemporaryOverlapDefinition."\n\n";
	  }
	my $OverlapCommand = $QHENBodyQuasiHoleOverlap." --global-overlap --input-file ".$TemporaryFileName." --output-precision ".$Precision;
	if ($LatexFlag == 1)
	  {
	    $OverlapCommand .= " --latex-output";
	  }
	
	my $OverlapOutput = `$OverlapCommand`;	
	if ($? != 0)
	  {
	    my $ErrorCode = $?;
#	    `rm -f $TemporaryFileName`;
	    die ("command ".$OverlapCommand." failed with error code ".$ErrorCode." while processing N=".$NbrParticles." 2S=".$SValue."\n");
	  }
	if ($DebugFlag == 1)
	  {
	    print $OverlapOutput."\n===============================================================\n\n";
	  }
	if ($LatexFlag == 1)
	  {
	    chomp ($OverlapOutput);
	    my @TmpArray = split (/\n/, $OverlapOutput);
	    my $LatexValue = pop (@TmpArray);
	    my $TmpLine = pop (@TmpArray);
	    while ((defined($TmpLine)) && ($TmpLine ne "latex output:"))
	      {
		if (!($TmpLine =~ /^\s*$/))
		  {
		    $LatexValue = $TmpLine;
		  }
		$TmpLine = pop (@TmpArray);
	      }
	    $FinalOutput .= "\$".$NbrParticles;
	    if ((($NbrParticles % 2) == 1) && (($SValue % 2) == 1))
	      {
		$FinalOutput .= "^*";
	      }
	    $FinalOutput .= "\$"." & ".$LatexValue."\n";
	    my @TmpArray2 = split (/ /, $SpectrumOutput);
	    if (($#TmpArray2 - 1) > $MaxLz)
	      {
		$MaxLz = $#TmpArray2 - 1;
	      }
	  }	
	else
	  {
	    $OverlapOutput =~ s/^.*\-{2,}\n\n//mg;
	    $FinalOutput .= "N = ".$NbrParticles."\n".$OverlapOutput;
	  }
	`rm -f $TemporaryFileName`;
      }
    $NbrParticles += $KValue;
    if ($FermionFlag == 1)
      {
	$SValue += $KValue + 2;
      }
    else
      {
	$SValue += 2;
      }
  }

if ($LatexFlag == 1)
  {
    print "\\begin{table*}
\\begin{ruledtabular}
\\begin{tabular}{cc";
    my $Index = 0;
    while ($Index <= $MaxLz)
      {
	print "c";
	$Index++;
      }
    print "}
\$N\$ & \$\\cal{O}^{".$KValue.",".$NbrHoles."}\$ & \$\\cal{O}_0^{".$KValue.",".$NbrHoles."}\$ ";
    $Index = 1;
    while ($Index <= $MaxLz)
      {
	print " & \$\\cal{O}_".$Index."^{".$KValue.",".$NbrHoles."}\$";
	$Index++;
      }
    print " \\\\\\hline\n";
  }
print $FinalOutput;
if ($LatexFlag == 1)
  {
    my $Fraction = "";
    if (($KValue % 2) == 0)
      {
	$Fraction = ($KValue / 2);
      }
    else
      {
	$Fraction = $KValue."/ 2";
      }
    my $RealNbrHoles = $KValue * $NbrHoles;
    print "\\end{tabular}
\\end{ruledtabular}
\\caption{overlap of the lowest energy exact wave functions and the \$\\nu=".$Fraction."\$ ".$RealNbrHoles." quasi-hole excitation wave functions. \$L\$ values for rows with a star, have to be understood as \$L-1/2\$}
\\end{table*}\n";
  }

# get temporary file name
#
# return value = file name

sub GetTemporaryFileName
  {
    my $FileName = "tmp".time();
    my $Pattern1 = "^".$FileName.".*";
    my $Pattern2 = "^".$FileName."(\\d*).*";
    my $TmpFile;
    my $Maximum = -1;
    foreach $TmpFile (<tmp*>)
      {
	if ((-f $TmpFile) && ($TmpFile =~ /\b$Pattern1\b/))
	  {
	    $TmpFile =~ /\b$Pattern2\b/;
	    if (($1 ne "") && ($1 > $Maximum))
	      {
		$Maximum = $1;
	      }
	  }
      }
    if ($Maximum != -1)
      {
	$Maximum++;
	$FileName .= $Maximum;
      }
    else
      {
	$FileName .= "0";
      }
    return $FileName;
  }
