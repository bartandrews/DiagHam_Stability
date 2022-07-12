#!/usr/bin/perl -w

use strict 'vars';

use Getopt::Long;


my $PathToLzToLProgram = "/home/regnault/development/Physics/DiagHam/build/FQHE/src/Programs/FQHEOnSphere/LzToL";
my $CheckValidity = 0;
my $ValidtyError = 1e-12;
my $FilePattern;
my $FullFlag = 0;

my $Result = GetOptions ("progdiag:s" => \$PathToLzToLProgram, 
			 "file-pattern=s" => \$FilePattern,
			 "check-validity" => \$CheckValidity,
			 "validity-error:s" => \$ValidtyError,
			 "full" => \$FullFlag);



$FilePattern =~ /\_n\_(\d+)\_2s\_(\d+)/;
my $NbrParticles = $1;
my $LzMax = $2;

if (!(defined($NbrParticles)))
  {
    die ("can't find number of particles in file pattern ".$FilePattern."\n");
  }
$FilePattern =~ /^(.+)\_sz\_\d+\_iz\_\d+\_pz\_\d+\_lz\.dat$/;
my $Prefix = $1;

my $TmpFile;
my %Sprectra;
my $GroundStateEnergy = 1000000.0;
foreach $TmpFile (<*>)
  {
    if ($TmpFile =~ /^$Prefix\_sz\_(\d+)\_iz\_(\d+)\_pz\_(\d+)\_l\.dat$/)
      {
	if (!defined($Sprectra{"sz_".$1."_iz_".$2}))
	  {
	    my @TmpArray = ($Prefix."_sz_".$1."_iz_".$2."_pz_".$3."_l.dat");
	    $Sprectra{"sz_".$1."_iz_".$2} = \@TmpArray;
	  }
	else
	  {
	    my $TmpArray = $Sprectra{"sz_".$1."_iz_".$2};
	    push (@$TmpArray, $Prefix."_sz_".$1."_iz_".$2."_pz_".$3."_l.dat");
	  }
	unless (open(INFILE, $TmpFile))
	  {
	    die ("can't open ".$TmpFile."\n");
	  }
	my $TmpLine;
	while (defined($TmpLine = <INFILE>))
	  {
	    chomp ($TmpLine);
	    my @TmpArray = split (/ /, $TmpLine);
	    if ($TmpArray[1] < $GroundStateEnergy)
	      {
		$GroundStateEnergy = $TmpArray[1];		
	      }
	  }
	close (INFILE);
      }
  }

my $TmpArray;
my $SzIz;
my @PointSytles = (70, 11, 31, 63, 2);
while (($SzIz, $TmpArray) = each (%Sprectra))
  {
    &FindGroundState($TmpArray, $GroundStateEnergy, $FullFlag);
    my @SortedFiles = sort {$a =~ /\_pz\_(\d+)\_/; my $Pz1 = $1; $b =~ /\_pz\_(\d+)\_/; my $Pz2 = $1; return ($Pz1 <=> $Pz2);} (@$TmpArray);
    my @GnuplotPlot;
    my $PointStyleIndex = 0;
    foreach $TmpFile (@SortedFiles)
      {
	$TmpFile =~ /\_pz\_(\d+)\_/;
	my $Pz = $1;
	if ($FullFlag == 0)
	  {
	    $TmpFile =~ s/l\.dat$/l\.ground.dat/;	
	  }
	else
	  {
	    $TmpFile =~ s/l\.dat$/l\.full.dat/;	
	  }
	push (@GnuplotPlot,  "\"".$TmpFile."\" using 1:2 title \"2P_z=".$Pz."\" with points ".$PointSytles[$PointStyleIndex]);
	$PointStyleIndex++;
      }    
    my $OutputPSFile = shift (@SortedFiles);
    $OutputPSFile =~ s/\_pz\_(\d+)\_/\_/;
    if ($FullFlag == 0)
      {
	$OutputPSFile =~ s/\_l\.ground\.dat$/.ground.eps/;
      }
    else
      {
	$OutputPSFile =~ s/\_l\.full\.dat$/.full.eps/;
      }
	
    $SzIz =~ /^sz\_(\d+)\_iz\_(\d+)$/;
    my $Sz = $1;
    my $Iz = $2;
    unless (open (OUTFILE, ">tmp.p"))
      {
	die ("can't create tmp.p file\n");
     }
    print OUTFILE "set xrange [0:10]
set yrange [0:]
set xlabel \"total angular momentum\"
set ylabel \"energy (arb. unit.)\"
set size 0.5, 1.0
set key top left height 1 box
set title \"N=".$NbrParticles." 2S=".$LzMax.", 2S_z=".$Sz." 2I_z=".$Iz."\"
set terminal postscript landscape enhanced eps \"Helvetica\" 18
set output \"".$OutputPSFile."\"
plot ".join (", ", @GnuplotPlot)."\n";
    close (OUTFILE);
    `gnuplot tmp.p`;
  }



sub FindGroundState()
  {
    my $Files = $_[0];
    my $GroundStateEnergy = $_[1];
    my $FullFlag = $_[2];
    my @SortedFiles = sort {$a =~ /\_pz\_(\d+)\_/; my $Pz1 = $1; $b =~ /\_pz\_(\d+)\_/; my $Pz2 = $1; return ($Pz1 <=> $Pz2);} (@$Files);
    my $TmpFile;
    my $MaxL = 0;
    my @TmpSpectra;
    my @TmpSpectrumMaxL;
    foreach $TmpFile (@SortedFiles)
      {
	print $TmpFile."\n";
	my $TmpSpectrum = "";
	unless (open(INFILE, $TmpFile))
	  {
	    die ("can't open ".$TmpFile."\n");
	  }
	my $TmpLine;
	my $CurrentMaxL = 10000;
	while (defined($TmpLine = <INFILE>))
	  {
	    chomp ($TmpLine);
	    my @TmpArray = split (/ /, $TmpLine);
	    if ($FullFlag == 0)
	      {
		if ($TmpArray[0] < $CurrentMaxL)
		  {
		    $CurrentMaxL = $TmpArray[0];		
		    $TmpSpectrum .= $TmpArray[0]." ".($TmpArray[1] - $GroundStateEnergy)."\n";
		  }
	      }
	    else
	      {
		$TmpSpectrum .= $TmpArray[0]." ".($TmpArray[1] - $GroundStateEnergy)."\n";
	      }
	  }
	close (INFILE);
	if ($FullFlag == 0)
	  {
	    $TmpFile =~ s/l\.dat$/l\.ground.dat/;
	  }
	else
	  {
	    $TmpFile =~ s/l\.dat$/l\.full.dat/;
	  }
	unless (open(OUTFILE, ">".$TmpFile))
	  {
	    die ("can't open ".$TmpFile."\n");
	  }
	print OUTFILE $TmpSpectrum;
	close (OUTFILE);
      }
  }
