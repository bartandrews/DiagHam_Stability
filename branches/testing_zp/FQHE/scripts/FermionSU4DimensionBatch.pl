#!/usr/bin/perl -w

use strict 'vars';

use Getopt::Long;


my $PathToDimensionProgram = "/home/regnault/development/Physics/DiagHam/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereGetDimension";
my $MinN = 0;
my $MaxN = 0;
my $MaximumDimension = 0;
my $MaxNbrFlux = 50;

my $Result = GetOptions ("progdim:s" => \$PathToDimensionProgram, 
			 "minnbrfermions=i" => \$MinN, "maxnbrfermions=i" => \$MaxN,
			 "maxdim=i" => \$MaximumDimension, "maxnbrflux:i");

while ($MinN <= $MaxN)
  {
    my $NbrFluxQuanta = int($MinN / 4) - 1;
    if (($NbrFluxQuanta < 0) || ($MinN != ($NbrFluxQuanta * 4)))
      {
	$NbrFluxQuanta++;
      }
    my $TmpMaxDim = 0;    
    while (($TmpMaxDim < $MaximumDimension) && ($NbrFluxQuanta <= $MaxNbrFlux))
      {
        my $FileName = "fermions_sphere_su4_n_".$MinN."_2s_".$NbrFluxQuanta.".dim";
	if (!(-e $FileName))
	  {
	    my $Command = $PathToDimensionProgram." -n ".$MinN." -s " .$NbrFluxQuanta." --su4-spin --save-disk";
	    print ("processing N=".$MinN." 2S=" .$NbrFluxQuanta."\n");
	    system($Command);
	  }
	$TmpMaxDim = &GetMaximumDimension($FileName);
	print "  maximum dimension = ".$TmpMaxDim."\n";
	$NbrFluxQuanta++;
      }
    $MinN++;
  }


# retrieve maximum dimension in a given Lz subspace
#
# $_[0] = file name
# return value = maximum dimension

sub GetMaximumDimension()
  {
    my $FileName = $_[0];
    my $MaxDim = 0;
    unless (open (INFILE, $FileName))
      {
	die ("error, can't open ".$FileName."\n");
      }
    my $TmpLine;
    while (defined($TmpLine = <INFILE>))
      {
	chomp($TmpLine);
	$TmpLine =~ s/\#.*$//;
	$TmpLine =~ s/\s+$//;
	$TmpLine =~ s/^\s+//;
	if ($TmpLine ne "")
	  {
	    my @TmpArray = split (/ /, $TmpLine);
	    if ($MaxDim < $TmpArray[4])
	      {
		$MaxDim = $TmpArray[4];
	      }
	  }
      }
    close (INFILE);
    return $MaxDim;
  }

