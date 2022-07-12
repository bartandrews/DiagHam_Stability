#!/usr/bin/perl -w

use strict 'vars';

use Getopt::Long;


my $PathToDiagonalizationProgram = "/home/regnault/development/Physics/DiagHam/build/FQHE/src/Programs/FQHEOnSphere/QHEFermionsTwoBodyGeneric";
my $NbrPoints = 0;
my $NbrParticles = 0;
my $LzMax = 0;
my $DiagonalizationProgramOptions = "--show-itertime -S --processors 4 --memory 6500";
my $InitialPseudopotentialFilename;
my $FinalPseudopotentialFilename;

my $Result = GetOptions ("progdiag:s" => \$PathToDiagonalizationProgram, 
			 "initial-pseudopotentials=s" => \$InitialPseudopotentialFilename,
			 "final-pseudopotentials=s" => \$FinalPseudopotentialFilename,
			 "nbrparticles=i" => \$NbrParticles, "lzmax=i" => \$LzMax,
			 "nbrpoints=i" => \$NbrPoints);

if (!defined($InitialPseudopotentialFilename))
  {
    die ("no initial pseudopotentials, use the initial-pseudopotentials option\n");
  }
if (!defined($FinalPseudopotentialFilename))
  {
    die ("no final pseudopotentials, use the final-pseudopotentials option\n");
  }

if ($NbrParticles == 0)
  {
    die ("can't find number of particles, use the nbrparticles option\n");
  }
if ($LzMax == 0)
  {
    die ("can't find lzmax, use the lzmax option\n");
  }


my @InitialPseudopotentials;
&ParsePseudopotentials($InitialPseudopotentialFilename, \@InitialPseudopotentials);
my @FinalPseudopotentials;
&ParsePseudopotentials($FinalPseudopotentialFilename, \@FinalPseudopotentials);


my $Lambda = 0.0;
my $Step = 1.0 / ($NbrPoints - 1);

my $OutputFile = "gap_n_".$NbrParticles."_2s_".$LzMax.".dat";
unless (open (OUTFILE, ">".$OutputFile))
  {
    die ("can't open ".$OutputFile."\n");
  }
close (OUTFILE);
my $DataFile = "fermions";
if ($PathToDiagonalizationProgram =~ /boson[^\/]*$/i)
  {
    $DataFile = "bosons";
  }
$DataFile .= "_mixed_n_".$NbrParticles."_2s_".$LzMax."_lz.dat";
while ($NbrPoints > 0)
  {
    unless (open (OUTFILE, ">tmppseudopotential.dat"))
      {
	die ("can't open tmppseudopotential.dat\n");
      }
    my $Index = 0;
    my @TmpPseudopotentials;
    while ($Index <= $LzMax)
      {
	$TmpPseudopotentials[$Index] = ((1.0 - $Lambda) * $InitialPseudopotentials[$Index]) + ($Lambda * $FinalPseudopotentials[$Index]);
	$Index++;
      }
    print OUTFILE "Pseudopotentials = ".join(" ", @TmpPseudopotentials)."\n";
    close(OUTFILE);
    system ($PathToDiagonalizationProgram." -p ".$NbrParticles." -l ".$LzMax." -n 1 --nbr-lz 2 --interaction-name mixed --interaction-file tmppseudopotential.dat ".$DiagonalizationProgramOptions);
    my $Gap = 0.0;
    unless (open (INFILE, $DataFile))
      {
	die ("can't open ".$DataFile."\n");
      }
    my $TmpLine = <INFILE>;
    chomp ($TmpLine);
    $TmpLine =~ s/^\d+ //;
    $Gap = -$TmpLine;
    $TmpLine = <INFILE>;
    chomp ($TmpLine);
    $TmpLine =~ s/^\d+ //;
    $Gap += $TmpLine;
    if ($Gap < 0.0)
      {
	$Gap = 0.0
      }
    close (INFILE);
    rename ($DataFile, $DataFile.".".$Lambda);
    unless (open (OUTFILE, ">>".$OutputFile))
      {
	die ("can't open ".$OutputFile."\n");
      }
    print OUTFILE ($Lambda." ".$Gap."\n");
    close (OUTFILE);    
    unlink("tmppseudopotential.dat");
    $Lambda += $Step;
    $NbrPoints--;
  }


#    if ($PlotFlag == 1)
#      {
#	$MinValue += $Step;
#	unless (open (OUTFILE, ">".$DataFileName.".p"))
#	  {
#	    die ("can't create file ".$DataFileName.".p\n");
#	  }
#	print OUTFILE ("set xrange [0:".$MinValue."]
#set yrange [0:1.1]
#set xlabel \"V2/V0\"
#set ylabel \"overlap\"
#set size 1, 0.9
#set terminal postscript landscape enhanced \"Helvetica\" 14
#set output \"".$DataFileName.".ps\"
#plot \"".$DataFileName.".overlap\" using 1:2 title \"N=".$NbrFermions." 2S=".$LzMax."\", \"".$DataFileName.".overlap\" using 1:2 notitle with lines
#");  
#	my $Command = "gnuplot ".$DataFileName.".p";
#	`$Command`;
#	unlink ($DataFileName.".p");
#      }
#  }

# get pseudo-potentials from description file
#
# $_[0] = name of the file that contains the pseudo-potentials
# $_[1] = reference on the array where pseudo-potentials have to be stored

sub ParsePseudopotentials
  {
    my $FileName = $_[0];
    my $Pseudopotentials = $_[1];
    unless (open(INFILE, $FileName))
      {
	die ("can't open pseudo-potential description ".$FileName."\n");
      }
    my $TmpLine;
    while (defined($TmpLine = <INFILE>))
      {
	chomp($TmpLine);
	if ($TmpLine =~ /^\s*Pseudopotentials\s*\=/)
	  {
	    $TmpLine =~ s/^\s*Pseudopotentials\s*\=\s*//;
	    @$Pseudopotentials = split (/ /, $TmpLine);
	  }
      }
    close (INFILE);
  }
