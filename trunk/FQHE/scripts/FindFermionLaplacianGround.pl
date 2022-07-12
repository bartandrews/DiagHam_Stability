#!/usr/bin/perl -w

use strict 'vars';

my $NbrFermions = 3;
while ($NbrFermions < 40)
  {
#    my $SValue = 3 * ($NbrFermions - 1);
      my $SValue = $NbrFermions + 1;
      while ($SValue < 40)
{
    my $TmpFile = "n_".$NbrFermions."/fermions_laplaciandelta_n_".$NbrFermions."_2s_".$SValue."_lz.dat";
    if (-e  $TmpFile)
      {
	my $Energy = &FindMin($TmpFile);
	if ($Energy < 0)
	{
#       print ($NbrFermions." ".$SValue." ".$Energy."\n");
	print ("./QHEFermionsLaplacianDelta -S -p ".$NbrFermions." -l ".$SValue." -n 40\n");
    }
      }
$SValue++;    
}
    $NbrFermions++;
  }

# find minimum in a file
#
# $_[0] = file name
# return value = ground state energy

sub FindMin
  {
    my $FileName = $_[0];
    my $Min;
    my $Flag = 0;
    open (INFILE, $FileName);
    my $TmpLine;
    foreach $TmpLine (<INFILE>)
      {
	chomp ($TmpLine);
	my @TmpArray = split (/ /, $TmpLine);
	if ($Flag == 0)
	  {
	    $Min = $TmpArray[1];
	    $Flag = 1;
	  }
	else
	  {
	    if ($TmpArray[1] < $Min)
	    {
	      $Min = $TmpArray[1];
	    }
	  }
      }
    close (INFILE);
    return $Min;
  }
