#!/usr/bin/perl -w

use strict 'vars';

use Getopt::Long;


my $PathToDiagonalizationProgram = "/home/regnault/development/DMRG/DiagHam/build/src/Programs/QHE/QHEOnSphere/QHEBosonsTwoBodyGeneric";
my $PathToOverlapProgram = "/home/regnault/development/DMRG/DiagHam/build/src/Programs/QHE/QHEOnSphere/QHENBodyQuasiHoleOverlap";
my $PathToDegeneracyProgram = "/home/regnault/development/DMRG/DiagHam/build/src/Programs/QHE/QHEOnSphere/ParafermionQuasiholeDimension";
my $Step = 0;
my $MinValue = 0.0;
my $NbrValues = 0;
my $NbrBosons = 0;
my $LzMax = 0;
my $ReferenceVector = "";
my $PlotFlag = 0;
my $DiagonalizationProgramOptions = "";
my $GridFile = "";
my $KeepFlag = 0;
my $ShowBestOverlaps = 1;
my $KValue = 0;
my $MaxDegeneracy = 0;
my $OverlapOnlyFlag = 0;

my $Result = GetOptions ("progdiag:s" => \$PathToDiagonalizationProgram, "progoverlap:s" => \$PathToOverlapProgram, 
			 "diagoption:s" => \$DiagonalizationProgramOptions, "progdeg:s" => \$PathToDegeneracyProgram,
			 "ref=s" => \$ReferenceVector, "nbrbosons:i" => \$NbrBosons, "lzmax:i" => \$LzMax,
			 "step:f" => \$Step, "min:f" => \$MinValue, "nbrvalues:i" => \$NbrValues, "plot" => \$PlotFlag,
			 "kvalue:f" => \$KValue, "keep" => \$KeepFlag, "grid:s" => \$GridFile, "maxdeg:i" => \$MaxDegeneracy,
			 "overlaponly" => \$OverlapOnlyFlag);



if ($GridFile eq "")
  {
    if ($NbrValues <= 0)
      {
	die ("invalid value for nbrvalues option\n");
      }
    if ($Step <= 0)
      {
	die ("invalid value for step option\n");
      }
  }

if (($ReferenceVector eq ""))
  {
    die ("invalid value for ref option, or file ".$ReferenceVector." does not exist\n");
  }
if ($NbrBosons == 0)
  {
    $ReferenceVector =~ /\_n\_(\d+)/;
    $NbrBosons = $1;
    if (!defined($NbrBosons))
      {
	die ("can't find number of bosons from file name ".$ReferenceVector.", use the nbrbosons option\n");
      }
  }
if ($LzMax == 0)
  {
    $ReferenceVector =~ /\_2s\_(\d+)/;
    $LzMax = $1;
    if (!defined($LzMax))
      {
	die ("can't find lzmax from file name ".$ReferenceVector.", use the lzmax option\n");
      }
  }
if ($KValue == 0)
  {
    $ReferenceVector =~ /\_nbody\_(\d+)/;
    $KValue = $1;
    if (!defined($KValue))
      {
	die ("can't find k value  from file name ".$ReferenceVector.", use the kvalue option\n");
      }
    $KValue--;
  }


my $DataFileName = "bosons_v2v0_n_".$NbrBosons."_2s_".$LzMax;
unless (open (OUTFILE, ">".$DataFileName.".overlap"))
  {
    die ("can't open ".$DataFileName.".overlap\n");
  }
close(OUTFILE);

my $Command = $PathToDegeneracyProgram." --lz-values -p ".$NbrBosons." -k ".$KValue." -q ".int(($LzMax + 2) * $KValue - (2 * $NbrBosons));
my @TmpArrayDegeneracy = split (/\n/, `$Command`);
my $Degeneracy = shift (@TmpArrayDegeneracy);
my @TmpArrayDegeneracy2 = split (/ /, $Degeneracy);
my $NbrLz = $#TmpArrayDegeneracy2;
if ($MaxDegeneracy < $TmpArrayDegeneracy2[0])
  {
    $MaxDegeneracy = $TmpArrayDegeneracy2[0];
  }
#$MaxDegeneracy = 2;

if ($GridFile eq "")
  {
    while ($NbrValues > 0)
      {
	unless (open (OUTFILE, ">tmppseudopotential.dat"))
	  {
	    die ("can't open tmppseudopotential.dat\n");
	  }
	print OUTFILE "Pseudopotentials = 1 0 ".$MinValue;
	my $TmpLz = 3;
	while ($TmpLz <= $LzMax)
	  {
	    print OUTFILE " 0";
	    $TmpLz++;
	  }
	print OUTFILE "\n";
	close(OUTFILE);
	
	my $DiagOutputFileName = $MinValue;
	$DiagOutputFileName =~ s/\./\_/;
	$DiagOutputFileName = "v2v0".$DiagOutputFileName;
	my $Overlap = &EvaluateOverlap($PathToDiagonalizationProgram, $PathToOverlapProgram, $NbrBosons, $LzMax, $ReferenceVector, $DiagOutputFileName, $KeepFlag,
				      $Degeneracy, $MaxDegeneracy, $NbrLz, $OverlapOnlyFlag);
	unless (open (OUTFILE, ">>".$DataFileName.".overlap"))
	  {
	    die ("can't open ".$DataFileName.".overlap\n");
	  }
	print OUTFILE $MinValue." ".$Overlap."\n";
	close(OUTFILE);

	unlink("tmppseudopotential.dat");
	$NbrValues--;
	$MinValue += $Step;
      }
    
    if ($PlotFlag == 1)
      {
	$MinValue += $Step;
	unless (open (OUTFILE, ">".$DataFileName.".p"))
	  {
	    die ("can't create file ".$DataFileName.".p\n");
	  }
	print OUTFILE ("set xrange [0:".$MinValue."]
set yrange [0:1.1]
set xlabel \"V2/V0\"
set ylabel \"overlap\"
set size 1, 0.9
set terminal postscript landscape enhanced \"Helvetica\" 14
set output \"".$DataFileName.".ps\"
plot \"".$DataFileName.".overlap\" using 1:2 title \"N=".$NbrBosons." 2S=".$LzMax."\", \"".$DataFileName.".overlap\" using 1:2 notitle with lines
");  
	my $Command = "gnuplot ".$DataFileName.".p";
	`$Command`;
	unlink ($DataFileName.".p");
      }
  }
else
  {
    my @GridVertices;
    &ParseGridDefinition ($GridFile, \@GridVertices, $LzMax);
    if ($#GridVertices == 0)
      {
	die ($GridFile." is an invalid grid file\n");
      }
    my $VertexCount = 0;
    my $VertexPosition;
    foreach $VertexPosition (@GridVertices)
      {
	print "processing ".$VertexCount."/".$#GridVertices."\r";
	unless (open (OUTFILE, ">tmppseudopotential.dat"))
	  {
	    die ("can't open tmppseudopotential.dat\n");
	  }
	print OUTFILE "Pseudopotentials = ".$VertexPosition."\n";
	close(OUTFILE);

	my $DiagOutputFileName = "v2v0".$VertexCount;
	my $Overlap = &EvalauteOverlap($PathToDiagonalizationProgram, $PathToOverlapProgram, $NbrBosons, $LzMax, $ReferenceVector, $DiagOutputFileName, $KeepFlag,
				      $Degeneracy, $MaxDegeneracy, $NbrLz, $OverlapOnlyFlag);

	unless (open (OUTFILE, ">>".$DataFileName.".overlap"))
	  {
	    die ("can't open ".$DataFileName.".overlap\n");
	  }
	print OUTFILE $VertexPosition."|".$Overlap."\n";
	close(OUTFILE);

	unlink("tmppseudopotential.dat");
	$VertexCount++;
      }
    if ($ShowBestOverlaps == 1)
      {
	my $Index = 0;
	my %BestOverlaps;
	my @BestOverlapPseudoPotentials;
	unless (open (INFILE, $DataFileName.".overlap"))
	  {
	    die ("can't open ".$DataFileName.".overlap\n");
	  }
	my $TmpLine;
	while (defined($TmpLine = <INFILE>))
	  {
	    chomp ($TmpLine);
	    if ($TmpLine ne "")
	      {
		my @TmpArray = split (/\|/, $TmpLine);
		push (@BestOverlapPseudoPotentials, $TmpArray[0]);
		$BestOverlaps{$Index} = $TmpArray[1];
		$Index++;
	      }
	  }
	close(INFILE);
	my @SortedIndices= sort {$BestOverlaps{$b} <=> $BestOverlaps{$a}}(keys(%BestOverlaps));
	$Index = 0;
	while (($Index <= $#SortedIndices) && ($Index < 10))
	  {
	    print $Index." : ".$BestOverlaps{$SortedIndices[$Index]}." (".$BestOverlapPseudoPotentials[$SortedIndices[$Index]].")\n";
	    $Index++;
	  }
      }
  }


# evaluate overlap betwwen maximally symmetric ground state evalauted for a given pseudopotential configuration, and a given eigenstate
#
# $_[0] = diagonalization program (with full path if any)
# $_[1] = overlap program (with full path if any)
# $_[2] = number of bosons
# $_[3] = maximum Lz value
# $_[4] = file that contains the reference vector (with path)
# $_[5] = interaction name
# $_[6] = flag to indicate if partial datas have to be kept (0 if false)
# $_[7] = quasihole ground state L-degeneracy
# $_[8] = maximum degeneracy in a given Lz sector
# $_[9] = number of Lz value to evaluate
# $_[10] = 1 if only overlaps have to be evaluated (eigenvectores and eigenvalues are retrieved from current directory instead of being evaluated)
# return value = square overlap

sub EvaluateOverlap()
  {
    my $PathToDiagonalizationProgram = $_[0];
    my $PathToOverlapProgram = $_[1];
    my $NbrBosons = $_[2];
    my $LzMax = $_[3];
    my $ReferenceVector = $_[4];
    my $DiagOutputFileName = $_[5];
    my $KeepFlag = $_[6];
    my $Degeneracy = $_[7];
    my $MaxDegeneracy = $_[8];
    my $NbrLz = $_[9];
    my $OverlapOnlyFlag = $_[10];

    my $ParityValue = 0;
    my $MaxLz = 2 * ($NbrLz - 1);
    if ((($NbrBosons % 2) == 1) && (($LzMax % 2) == 1))
      {
	$ParityValue = 1;
	$MaxLz += 1;
      }

    if ($OverlapOnlyFlag != 1)
      {
	my $StartLz = 0;
	if ($ParityValue == 1)
	  {
	    $StartLz++;
	  }
	my $ReducedNbrLz = $NbrLz;
	while (($ReducedNbrLz >= 0) && (-e $DiagOutputFileName2."_".$StartLz.".".($MaxDegeneracy - 1).".vec"))
	  {
	    $StartLz += 2;
	    $ReducedNbrLz--;
	  }
	if ($ReducedNbrLz >= 0)
	  {
	    if ($StartLz > 0)
	      {
		rename($DiagOutputFileName2.".dat", $DiagOutputFileName2.".dat.0")
	      }
	    my $Command = $PathToDiagonalizationProgram." -p ".$NbrBosons." -l ".$LzMax." --initial-lz ".$StartLz." --nbr-lz ".($ReducedNbrLz + 1)." -n ".$MaxDegeneracy." --force-reorthogonalize --eigenstate --interaction-name ".$DiagOutputFileName." --interaction-file tmppseudopotential.dat  --full-diag 3000 -S ".$DiagonalizationProgramOptions;
	    system($Command);
	    if ($StartLz > 0)
	      {
		$Command = "cat ".$DiagOutputFileName2.".dat >> ".$DiagOutputFileName2.".dat.0";
		`$Command`;
		unlink ($DiagOutputFileName2.".dat");
		rename($DiagOutputFileName2.".dat.0", $DiagOutputFileName2.".dat")
	      }
	  }
      }

    my $DiagOutputFileName2 = "bosons_".$DiagOutputFileName."_n_".$NbrBosons."_2s_".$LzMax."_lz";

    unless (open(OUTFILE, ">tmpoverlap.dat"))
      {	
	die ("can't create  tmpoverlap.dat\n");
      }
    print OUTFILE ("NbrParticles=".$NbrBosons."
LzMax=".$LzMax."
Degeneracy=".$Degeneracy."
Spectrum=".$DiagOutputFileName2.".dat
QuasiholeStates=".$ReferenceVector."
ExactStates=".$DiagOutputFileName2."_\n");
    close (OUTFILE);

    my $Overlap = "-1";
    if (-e $DiagOutputFileName2."_".$ParityValue.".0.vec")
      {
	$Command = $PathToOverlapProgram." --input-file tmpoverlap.dat --global-overlap --latex-output";
	$Overlap = `$Command`;
	print $Overlap;
	my @TmpArray = split (/\n/, $Overlap);
	my $TmpLine = shift(@TmpArray);
	while (($TmpLine ne "latex output:") && ($TmpLine ne "no possible overlap calculation"))
	  {
	    $TmpLine = shift(@TmpArray);
	  }
	if ($TmpLine eq "latex output:")
	  {
	    $Overlap = shift(@TmpArray);
	    $Overlap =~ s/\&  \&/ 0 /mg; 
	    $Overlap =~ s/[^\d\. e\-]//mg; 
	    $Overlap =~ s/ +/ /mg;
	  }
	else
	  {
	    $Overlap = $TmpLine;
	  }
      }

    if (($KeepFlag == 0) && ($OverlapOnlyFlag != 1))
      {
	my $Lz = $ParityValue;
	while ($Lz <= $MaxLz)
	  {
	    my $Index = 0;
	    while ($Index <= $MaxDegeneracy)
	      {
		if (-e $DiagOutputFileName2."_".$Lz.".".$Index.".vec")
		  {
		    unlink($DiagOutputFileName2."_".$Lz.".".$Index.".vec")
		  }
		$Index++;
	      }
	    $Lz += 2;
	  }
	if (-e $DiagOutputFileName2.".dat")
	  {
	    unlink ($DiagOutputFileName2.".dat");
	  }
      }
    return $Overlap;
  }


# parse a grid definition file
#
# $_[0] = name file that contains the grid definition
# $_[1] = reference on the array where vertex positions will be stored (each position is a string that gives pseudo-potential definition)
# $_[2] = maximum Lz value

sub ParseGridDefinition()
  {
    my $GridFile = $_[0];
    my $GridVertices = $_[1];
    my $LzMax = $_[2];

    unless (open (INFILE, $GridFile))
      {
	die ("can't open ".$GridFile."\n");
      }
    my %Pseudopotentials;
    my $TmpLine;
    while (defined($TmpLine = <INFILE>))
      {
	chomp ($TmpLine);
	$TmpLine =~ s/^\s+//;
	$TmpLine =~ s/\s+$//;
	$TmpLine =~ s/^\#.*//;
	if ($TmpLine ne "")
	  {
	    my @TmpArray = split (/\s+/, $TmpLine);
	    if (($#TmpArray != 3) || (!($TmpArray[0] =~ /^\d+$/)) || (!($TmpArray[1] =~ /^\-?\d*\.?\d+$/)) || (!($TmpArray[2] =~ /^\-?\d*\.?\d+$/))
		|| (!($TmpArray[3] =~ /^\d+$/)))
	      {
		close (INFILE);
		return;
	      }
	    my $Index = shift (@TmpArray);
	    $Pseudopotentials{$Index} = \@TmpArray;
	  }
      }
    close (INFILE);

    my $Index = 0;
    my @PseudopotentialIndices = sort {$a <=> $b} (keys (%Pseudopotentials));

    my @Strife;
    $Strife[0] = 1;
    $Index = 0;
    while ($Index <= $#PseudopotentialIndices)
      {
	my $TmpArray = $Pseudopotentials{$PseudopotentialIndices[$Index]};
	$Strife[$Index + 1] = $Strife[$Index] * $$TmpArray[2];
	$Index++;
      }

    $Index = 1;
    my $Total = $Strife[$#Strife];
    my $TmpString = "1";
    while ($Index < $PseudopotentialIndices[0])
      {
	$TmpString .= " 0";
	$Index++;
      }
    $Index = 0;
    while ($Index < $Total)
      {
	push (@$GridVertices, $TmpString);
	$Index++;
      }

    $Index = 0;
    while ($Index <= $#PseudopotentialIndices)
      {
	my $TmpArray2 = $Pseudopotentials{$PseudopotentialIndices[$Index]};
	my $Value = $$TmpArray2[0];
	my $Step = $$TmpArray2[1];
	my $NbrValues = $$TmpArray2[2];
	my @TmpStrings;
	my $Max = $LzMax + 1;
	if ($Index < $#PseudopotentialIndices)
	  {
	    $Max = $PseudopotentialIndices[$Index + 1];
	  }
	my $TmpString = "";
	my $Index2 = $PseudopotentialIndices[$Index] + 1;
	while ($Index2 < $Max)
	  {
	    $TmpString .= " 0";
	    $Index2++;
	  }
	$Index2 = 0;
	while ($Index2 < $NbrValues)
	  {
	    push (@TmpStrings, " ".$Value.$TmpString);
	    $Value += $Step;
	    $Index2++;
	  }
	
	$Index2 = 0;
	my $Index3 = 0;
	my $Index4 = 0;
	my $TmpStrife = $Strife[$Index];

	while ($Index2 < $Total)
	  {
	    $Index3 = 0;
	    while ($Index3 < $NbrValues)
	      {
		$TmpString = $TmpStrings[$Index3];
		$Index4 = 0;
		while ($Index4 < $TmpStrife)
		  {
		    $$GridVertices[$Index2] .= $TmpString;
		    $Index2++;
		    $Index4++;
		  }
		$Index3++;
	      }
	  }
	$Index++;
      }
  }
