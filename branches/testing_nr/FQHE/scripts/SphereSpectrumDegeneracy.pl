#!/usr/bin/perl -w

use strict 'vars';

use Getopt::Long;


my $SpectrumFile = "";
my $Eigenvalue = 0.0;
my $Error = 1.0e-12;
my $LFlag = 0;
my $RowFlag = 0;
my $SumFlag = 0;
my $LatexFlag = 0;
my $NbrParticles = -1;
my $NbrFlux = -1;
my $Result = GetOptions ("spectrum=s" => \$SpectrumFile, "eigenvalue:s" => \$Eigenvalue, "error:s" => \$Error, "lsort" => \$LFlag,
			"row" => \$RowFlag, "sum" => \$SumFlag, "latex:s" => \$LatexFlag, "particles:s" => \$NbrParticles,
			"flux" => \$NbrFlux); 

#if (($SpectrumFile eq "") || (!(-e $SpectrumFile)) || (!($Eigenvalue =~ /^[\+\-]?\d*\.?\d*e?\d+$/)) || (!()))
if ($SpectrumFile eq "")
  {
    die ("usage: SphereSpectrumDegenracy.pl --spectrum file_name [--eigenvalue 0.0 --error 1.0e-12 --lsort --sum --latex 0 --particles -1 --flux -1]\n");
  }

if (abs($Eigenvalue) < $Error)
  {
    $Eigenvalue = 0.0;
  }
my %Degeneracy;

unless (open (INFILE, $SpectrumFile))
  {
    die ("can't open ".$SpectrumFile."\n");
  }
my $TmpLine;
while (defined($TmpLine = <INFILE>))
  {
    chomp ($TmpLine);
    my @TmpArray = split (/ /, $TmpLine);
    if ((abs ($TmpArray[1] - $Eigenvalue) < abs($Error * $Eigenvalue)) || (($Eigenvalue == 0.0) && (abs ($TmpArray[1]) < $Error)))
      {
	if (defined($Degeneracy{$TmpArray[0]}))
	  {
	    $Degeneracy{$TmpArray[0]}++;
	  }
	else
	  {
	    $Degeneracy{$TmpArray[0]} = 1;
	  }
      }
  }
close (INFILE);

my $Sum = 0;
if ($SumFlag != 0)
  {
    if ($NbrParticles == -1)
      {
	$SpectrumFile =~ /\_n\_(\d+)\_/;
	$NbrParticles = $1;
      }
    if ($NbrFlux == -1)
      {
	$SpectrumFile =~ /\_2s\_(\d+)\_/;
	$NbrFlux = $1;
      }
    if ((($NbrParticles * $NbrFlux) % 2) == 0)
      {
	foreach $TmpLine (keys(%Degeneracy))
	  {
	    if ($TmpLine != 0)
	      {
		$Sum += 2 * $Degeneracy{$TmpLine};
	      }
	    else
	      {
		$Sum += $Degeneracy{$TmpLine};
	      }
	  }
      }
    else
      {
	foreach $TmpLine (keys(%Degeneracy))
	  {
	    $Sum += 2 * $Degeneracy{$TmpLine};
	  }
      }
  }

if ($LFlag != 0)
  {
    my $CurrentNbr = 0;
    foreach $TmpLine (sort {$b <=> $a} (keys(%Degeneracy)))
      {	
	$Degeneracy{$TmpLine} -= $CurrentNbr;
	$CurrentNbr += $Degeneracy{$TmpLine};
      }
    
  }

if (($LatexFlag == 0) && ($RowFlag == 0))
  {
    foreach $TmpLine (sort {$a <=> $b} (keys(%Degeneracy)))
      {
	print $TmpLine." ".$Degeneracy{$TmpLine}."\n";
      }
  }
else
  {
    if ($LatexFlag == 0)
      {
	foreach $TmpLine (sort {$a <=> $b} (keys(%Degeneracy)))
	  {
	    print $Degeneracy{$TmpLine}." ";
	  }
	print "\n";
      }
    else
      {
	my $Tmp = 0;
	foreach $TmpLine (sort {$a <=> $b} (keys(%Degeneracy)))
	  {
	    print "\$".$Degeneracy{$TmpLine}."\$";
	    if ($Tmp < $LatexFlag)
	      {
		print " & "; 
	      }
	    $Tmp++;
	  }
	while ($Tmp < $LatexFlag)
	  {
	    print " & "; 
	    $Tmp++;
	  }
	print "\\\\\n";
      }
  }

if ($SumFlag != 0)
  {
    print "Sum = ".$Sum."\n";
  }
