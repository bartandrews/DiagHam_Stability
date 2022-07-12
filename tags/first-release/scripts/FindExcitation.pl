#!/usr/bin/perl -w

use strict 'vars';


my $NbrBosons = 4;
while ($NbrBosons <= 40)
  {
    if (-e "n_".$NbrBosons)
      {
	chdir ("n_".$NbrBosons);
	my $SValue = 2;
	while ($SValue <= 20)
	  {
	    if (-e "bosons_delta_n_".$NbrBosons."_2s_".$SValue."_lz.dat")
	      {
		my $E = &FindMin("bosons_delta_n_".$NbrBosons."_2s_".$SValue."_lz.dat");
		my $ShiftedE = $E + (0.5 * $NbrBosons * $NbrBosons / ((0.5 * $SValue)));
		if (abs($ShiftedE) < 1e-12)
		  {
		    $ShiftedE = 0;
		  }
		print $ShiftedE." ";
	      }
	    else
	      {
		print "nan ";
	      }
	    ++$SValue;
	  }
	chdir ("..");
      }
    print "\n";
    ++$NbrBosons
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

# create postscript graph from data file
#
# $_[0] = hash table containing datas
# $_[1] = print flag (1 if true)
# $_[2] = number of bosons

sub CreatePostScript
  {
    my $Datas = $_[0];
    my $N = $_[1];
    my $S;
    my $E;
    my $OutputFile = "bosons_delta_ground_state_n_".$N.".log";
    open (OUTFILE, ">$OutputFile");
    while (($S, $E) = each (%$Datas))
      {
	my $ShiftedE = $E + (0.5 * $N * $N / (sqrt(0.5 * $S)));
	print OUTFILE ($S." ".$E." ".$ShiftedE."\n");
      }
    close (OUTFILE);
  }



