#!/usr/bin/perl -w

use strict 'vars';

if (!(defined($ARGV[1])))
  {
    die "usage: FindGapBilayer NbrFermions NbrFlux [Caption] [PrintFlag]";
  }
my $NbrFermions = $ARGV[0];
my $S = $ARGV[1];
my $Directory = "n_".$NbrFermions;
my $Caption = "";
if (defined($ARGV[2]))
  {
    my $Caption = $ARGV[2];
  }
my $PrintFlag = 0;
if (defined($ARGV[3]))
  {
    $PrintFlag = 1;
  }
my %MinArray;
my %LArray;
my $TmpFile;
my $L;
my $LayerSeparation;
# search for the case of entire spectra:
my $List=`ls $Directory/fermions_sphere_su2_coulomb_d_*_2s_${S}_*lz.dat $Directory/2S_$S/fermions_sphere_su2_coulomb_d_*_2s_${S}_*lz.dat`;
my @TmpFiles = split (/\n/, $List);
foreach $TmpFile (@TmpFiles)
  {
    print ("next file: ".$TmpFile." ");
    $TmpFile =~ /.*_d_([\d\.]*)_n.*/;
    $LayerSeparation = $1;
    print "d= ".$LayerSeparation."\n";
    $MinArray{$LayerSeparation} = &FindGap($TmpFile,\$L);
    print ("L= ".$L."\n");
    $LArray{$LayerSeparation} = $L;
  }
# search for the case of separately calculated lz=0 and lz=2:
$List=`ls $Directory/fermions_sphere_su2_coulomb_d_*_2s_${S}_*lz_0.dat $Directory/2S_$S/fermions_sphere_su2_coulomb_d_*_2s_${S}_*lz_0.dat`;
@TmpFiles = split (/\n/, $List);
my $TmpFile2;
foreach $TmpFile (@TmpFiles)
  {
    print ("next file: ".$TmpFile."\n");
    $TmpFile =~ /.*_d_([\d\.]*)_n.*/;
    $TmpFile2 = $TmpFile;
    $TmpFile2 =~ s/lz_0/lz_2/;    
    if (-e $TmpFile2)
      {
        $TmpFile =~ /.*_d_([\d\.]*)_n.*/;
        $LayerSeparation = $1;
        $MinArray{$LayerSeparation} = &CalculateGap($TmpFile,$TmpFile2,\$L);
        print ("L= ".$L."\n");
        $LArray{$LayerSeparation} = $L;
      }
  }
&CreatePostScript(\%MinArray, $Caption, $PrintFlag);
&PrintLs(\%LArray);

# find gap in a file
#
# $_[0] = file name
# $_[1] = value of groundstate angular momentum
# return value = ground state energy

sub FindGap
  {
    my $FileName = $_[0];
    my $L = $_[1];
    my $Min;
    my $Min2;
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
	    if ($Flag == 2)
	      {
		if ($TmpArray[1] < $Min)
		  {
		    $Min2 = $Min;
		    $Min = $TmpArray[1];
		  }
		else
		  {
		    if ($TmpArray[1] < $Min2)
		      {
			$Min2 = $TmpArray[1];
		      }
		  }
	      }
	    else
	      {
		if ($TmpArray[1] < $Min)
		  {
		    $Min2 = $Min;
		    $Min = $TmpArray[1];
		  }
		else
		  {
		    $Min2 = $TmpArray[1];
		  }
		$Flag = 2;
	      }
	  }
      }
    close (INFILE);
    my $gap=($Min2 - $Min);
    open (INFILE, $FileName);
    my $MaxLz=0;
    foreach $TmpLine (<INFILE>)
      {
	chomp ($TmpLine);
	my @TmpArray = split (/ /, $TmpLine);
	if ( abs($Min-$TmpArray[1]) < 1e-11)
	  {
	    if ( $TmpArray[0] > $MaxLz)
	      {
		$MaxLz = $TmpArray[0];
	      }
	  }
      }
    close(INFILE);
    if ( $MaxLz > 0 )
      {
	print ("     /vvv\\     \nGS at L=".$MaxLz." with gap ".$gap." in ".$FileName."\n     \\^^^/     \n");
	$gap=0.0;
      }
    $$L=$MaxLz;
    return $gap;
  }

# find gap from lz_0 and lz_2 files
#
# $_[0] = file name lz_0
# $_[1] = file name lz_2
# $_[2] = return value of groundstate angular momentum
# return value = ground state energy

sub CalculateGap
  {
    my $FileName = $_[0];
    my $FileName2 = $_[1];
    my $L = $_[2];
    my $E0;
    my $E2;
    my $Flag = 0;
    open (INFILE, $FileName);
    my $TmpLine;
    $TmpLine = <INFILE>;
    chomp ($TmpLine);
    my @TmpArray = split (/ /, $TmpLine);
    $E0 = $TmpArray[1];
    close(INFILE);
    open (INFILE, $FileName2);
    $TmpLine = <INFILE>;
    chomp ($TmpLine);
    @TmpArray = split (/ /, $TmpLine);
    $E2 = $TmpArray[1];
    close(INFILE);
    my $gap =  ($E2 - $E0);
    if (abs($gap) < 1e-11)
      {
	$gap=0.0;
	$$L=99;
      }
    else
      {
	$$L=0;
      }
    return $gap;
  }

# create postscript graph from data file
#
# $_[0] = hash table containing datas
# $_[1] = print flag (1 if true)
# $_[2] = number of bosons

sub CreatePostScript
  {
    my $Datas = $_[0];
    my $Caption = $_[1];
    my $PrintFlag = $_[2];
    my $D;
    my $E;
    my $basename = "fermions_sphere_su2_coulomb_n_".$NbrFermions."_2s_".$S;
    if ( ! $Caption eq "" )
      {
	$basename = $basename ."_".$Caption;
      }
    my $FileName = $basename.".dat";
    open (OUTFILE, ">$FileName");
    my $MinD = 200;
    my $MaxD = 0;
    my $MinGap = 4000;
    my $MaxGap = 0;
    while (($D, $E) = each (%$Datas))
      {
	if ($MinD > $D)
	  {
	    $MinD = $D;
	  }
	if ($MaxD < $D)
	  {
	    $MaxD = $D;
	  }
	if ($MinGap > $E)
	  {
	    $MinGap = $E;
	  }
	if ($MaxGap < $E)
	  {
	    $MaxGap = $E;
	  }
	print ($D." ".$E."\n");
	print OUTFILE ($D." ".$E."\n");
      }
    close (OUTFILE);
    $MinGap = 0;
    my $Delta = ($MaxGap - $MinGap) / 20.0;
    $MaxGap += $Delta;
    $MinGap -= $Delta;
    $MinGap = 0;
    $MinD-=0.25;
    $MaxD+=0.25;
    my $TmpFileName = "tmp".time().".p";
    my $OutputFile = $basename.".ps";
    my @TmpArray = split (/_/,  $OutputFile);
    my $Title = "gap: N=".$NbrFermions." S=".$S." ".$Caption;
    open (OUTFILE, ">$TmpFileName");
    print OUTFILE ("set xrange [".$MinD.":".$MaxD."]
set yrange [".$MinGap.":".$MaxGap."]
set xlabel \"d [l_0]\"
set ylabel \"E\"
set size 1.0, 0.6
set terminal postscript portrait enhanced \"Helvetica\" 14
set output \"".$OutputFile."\"
plot \"".$FileName."\" using 1:2 title \"".$Title."\"
");
    close (OUTFILE);
    `gnuplot $TmpFileName`;
    if ($PrintFlag == 1)
      {
	`lpr $OutputFile`;
      }
    `rm -f $TmpFileName`;
  }


# print table with GS angular momentum
#
# $_[0] = hash table containing datas

sub PrintLs
  {
    my $Datas = $_[0];
    my $D;
    my $L0;
    my $basename = "fermions_sphere_su2_coulomb_n_".$NbrFermions."_2s_".$S;
    if ( ! $Caption eq "" )
      {
	$basename = $basename ."_".$Caption;
      }
    my $FileName = $basename.".L0";
    open (OUTFILE, ">$FileName");
    while (($D, $L0) = each (%$Datas))
      {
	print ($D." ".$L0."\n");
	print OUTFILE ($D." ".$L0."\n");
      }
    close (OUTFILE);
    }
