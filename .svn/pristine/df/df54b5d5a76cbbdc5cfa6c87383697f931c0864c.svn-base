#!/usr/bin/perl -w

use strict 'vars';

if (!(defined($ARGV[0])))
  {
    die "usage: BosonsDeltaGraph nbr_bosons print_flag\n";
  }
my $NbrBosons = $ARGV[0];
my $PrintFlag = 0;
if (defined($ARGV[1]))
  {
    $PrintFlag = 1;
  }
my %MinArray;
my $TmpFile;
foreach $TmpFile (<*>)
  {
    if (($TmpFile =~ /bosons\_delta\_n\_$NbrBosons.*\_l\./ ) ||  ($TmpFile =~ /bosons\_coulomb\_n\_$NbrBosons.*\_l\./))
      {
	my $TmpFile2 = $TmpFile;
	$TmpFile2 =~ s/.*2s\_(\d*).*/$1/;
        $MinArray{$TmpFile2} = &FindMin($TmpFile);
	print ("$TmpFile2 ".$MinArray{$TmpFile2}."\n");
      }
  }
&CreatePostScript(\%MinArray, $NbrBosons);

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



