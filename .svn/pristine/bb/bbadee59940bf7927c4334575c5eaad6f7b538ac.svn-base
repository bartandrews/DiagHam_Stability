#!/usr/bin/perl -w
#
# script for analysis of eigenvector files designed to run on TCM group PC's
#
use strict 'vars';

# set to zero to suppress some output
my $debug = 0;

# do not include files with tag _rh_
my $noRH = 1;

if (!defined($ARGV[2]))
  {
    print("usage:   LatticePlotGapInSeries.pl directory offset slope [lowest_U]\n");
    print("=========================================================\n");
    print("script searches for all protocol files *.gs in 'directory' and subdir's\n");
    print("extracts the gap for states on line n_phi = offset + slope\n");
    print("please use fractions p/q for both slope and offset\n");
    print("lowest_U: do not include points where interaction is less than U, default: hardcore!\n");
    exit(1);
  }

my ($DirName,$offStr,$slopeStr) = @ARGV;

my ($P,$Q)=split(/\//,$offStr);
my ($R,$S)=split(/\//,$slopeStr);

my $Mod = &FindGCD($Q,$S);



my $tmp=$Q;
$P *= $S/$Mod;
$Q *= $S/$Mod;

$R *= $tmp/$Mod;
$S *= $tmp/$Mod;

print("analyzing n_phi=$P/$Q+($R/$S)n\n");

if ($S!=$Q)
  {
    die ("Error calculating smallest common denominator");
  }

my $lowestU="hc";
my $requireHardCore=1;
if (defined($ARGV[3]))
  {
    $requireHardCore=0;
    $lowestU = $ARGV[3];
  }

my $TmpFileName = "tmp".time().".p";
	
system("find $DirName -type f -name \"*.gs\" > ".$TmpFileName);
open (INFILE, $TmpFileName);

$DirName =~ s/\/$//;
$DirName =~ s/\//-/;

my $OutputName= "LatticeGap_".$DirName."_".$lowestU.".gap";

open (INFILE, $TmpFileName);
open (OUTFILE, ">$OutputName");
print OUTFILE ("# n_phi\tGap\tGap2\tN\tNPhi\tNs\tinteract\tsource\n");
my $FileName;
my $CountEntries=0;
my %CurveLabels = ();
my %ReverseCurveLabels = ();
my $CountLabels=0;
my %Curves = ();
foreach $FileName (<INFILE>)
  {
    chomp ($FileName);
    my $HardCore;
    my $N;
    my $x;
    my $y;
    my $u;
    my $q;
    if ($FileName =~ m/hardcore/)
      {
	$FileName =~ /\_n\_(\d+)\_x\_(\d*)\_y\_(\d*)\_/;
	$N = $1;
	$x = $2;
	$y = $3;
	$u = 0;
	$HardCore=1;
	$q = -1;
      }
    else
      {
	$FileName =~ /\_n\_(\d+)\_x\_(\d*)\_y\_(\d*)\_u\_(-*\d*[\.]*\d*)\_/;
	$N = $1;
	$x = $2;
	$y = $3;
	$u = $4;
	$HardCore=0;
	$q = -1;
      }
    if ( ( ( $HardCore == 1) || ( ($requireHardCore==0) && ($u >= $lowestU))) &&
	 (! ( ( $noRH==1) && ($FileName =~ m/\_rh\_/))))
      {
	my $labelU;
	my $textU;
	if ( $HardCore == 1)
	  {
	    print "including n = ".$N."  x = ".$x."  y = ".$y."  (hardcore bosons) from $FileName\n";
	    $textU = "hard";
	    $labelU = "hardcore";
	  }
	else
	  {
	    print "including n = ".$N."  x = ".$x."  y = ".$y."  u = ".$u." from $FileName\n";
	    $textU = "$u";
	    $labelU = "u=$u";
	  }
	++$CountEntries;
	my $Ns=$x*$y;

	if (($P*$Ns + $R*$N) % $Q == 0)
	  {
	    my $TargetNPhi = ($P*$Ns + $R*$N)/$Q;

	    my $TmpLine = `grep ^${TargetNPhi} ${FileName}`;
	    chomp ($TmpLine);
	    if ($TmpLine ne "")
	      {
		print("Target $TargetNPhi found: $TmpLine\n");

		my @Entries = split(/\t/,$TmpLine);
		my $SpecName=$FileName;
		$SpecName =~ s/.gs/\_${TargetNPhi}.dat/;
		if (! -e $SpecName)
		  {
		    $SpecName =~ s/\_${TargetNPhi}\.dat/.dat/;
		  }
		my @TmpArray = `grep ^${TargetNPhi} ${SpecName}`;
		my $Gap2=0.0;
		if (defined($TmpArray[2]))
		  {
		    my @Array2=split(/ /,$TmpArray[2]);
		    my @Array0=split(/ /,$TmpArray[0]);
		    $Gap2 = $Array2[1]-$Array0[1];
		  }
		my $OutLine = sprintf  ("%.6f\t%.10f\t%.10f\t%d\t%d\t%d\t%s\t%s\n",
					$TargetNPhi/$Ns, $Entries[3], $Gap2,
					$N, $TargetNPhi, $Ns, $textU, $SpecName);
		print OUTFILE ($OutLine."\n");
		my $CurveLabel = "N=$N, $labelU";
		if (!exists($CurveLabels{$CurveLabel}))
		  {
		    $CurveLabels{$CurveLabel}=$CountLabels;
		    $ReverseCurveLabels{$CountLabels} = "$N\t$textU";
		    ++$CountLabels;
		  }
		push( @{$Curves{$CurveLabels{$CurveLabel}}},$OutLine);
		
	      }
	    else
	      {
		print "Missing flux value: $TargetNPhi in ${FileName}";
	      }
	  }
	else
	  {
	    print("Not commensurate\n");
	  }
	
      }
    else
      {
	if ($debug)
	  {
	    print "excluding n = ".$N."  x = ".$x."  y = ".$y."  u = ".$u."\n";
	  }
      }


  }

system ("rm ".$TmpFileName);
close(OUTFILE);
close(INFILE);

if (! -e $OutputName."_files" )
  {
    system ("mkdir ".$OutputName."_files");
  }

$TmpFileName = "tmp".time().".p";
	
open (OUTFILE, ">$TmpFileName");

my @OutFileN;
my $CurveLabel;
foreach $CurveLabel (keys %CurveLabels)
  {
    my $numCurve=$CurveLabels{$CurveLabel};
    print ("Label: ".$CurveLabel." has value ".$CurveLabels{$CurveLabel}." and entries\n");
    my @sorted = sort @{$Curves{$CurveLabels{$CurveLabel}}};
    my @params = split("\t",$ReverseCurveLabels{$CurveLabels{$CurveLabel}});
    $OutFileN[$CurveLabels{$CurveLabel}] = $OutputName."_files/points_n_".$params[0]."_u_".$params[1].".dat";
    open (OUTFILE2, ">".$OutFileN[$CurveLabels{$CurveLabel}]);
    foreach my $Line (@sorted)
      {	
	print OUTFILE2 ($Line);
      }
    close (OUTFILE2);
    print OUTFILE ("READ BLOCK \"".$OutFileN[$CurveLabels{$CurveLabel}]."\"\n");
    print OUTFILE ("BLOCK xy \"1:2\"\n");
    print OUTFILE ("S_ legend \"${CurveLabel}\"\n");
    print OUTFILE ("S_ comment \"".$OutFileN[$CurveLabels{$CurveLabel}]." gap\"\n");
    print OUTFILE ("S_ linestyle 0\n");
    print OUTFILE ("S_ symbol color ".(($numCurve+1)%15+1)."\n");
    print OUTFILE ("S_ symbol ".(($numCurve+1)%11+1)."\n");
    print OUTFILE ("S_ symbol fill pattern 1\n");
    print OUTFILE ("S_ symbol fill color ".(($numCurve+1)%15+1)."\n");
    print OUTFILE ("BLOCK xy \"1:3\"\n");
    print OUTFILE ("S_ comment \"".$OutFileN[$CurveLabels{$CurveLabel}]." gap 2\"\n");
    print OUTFILE ("S_ linestyle 0\n");
    print OUTFILE ("S_ symbol color ".(($numCurve+1)%15+1)."\n");
    print OUTFILE ("S_ symbol ".(($numCurve+1)%11+1)."\n");
    print OUTFILE ("S_ symbol fill pattern 0\n");
    print OUTFILE ("S_ symbol fill color ".(($numCurve+1)%15+1)."\n");
  }

print OUTFILE ("xaxis label \"n\\s\\xf\"\n");
print OUTFILE ("yaxis label \"\\xD\"\n");
close (OUTFILE);

system("xmgrace -batch ${TmpFileName} -nosafe -noask &");
system("sleep 2");

system ("rm $TmpFileName");



# find greatest common divider
#
# m = first integer
# n = second integer
# return value = GCD

sub FindGCD
  {
    my ($m, $n) = @_;
    if ($m < $n)
      {
	return RecursiveFindGCD ($m, $n);
      }
    else
      {
	return RecursiveFindGCD ($n, $m);
      }
    return $n;
  }


# find greatest common divider (recurisive part of the method)
#
# m = first integer
# n = second integer (must be greater than m)
# return value = GCD

sub RecursiveFindGCD
  {
    my ( $m, $n) = @_;
    if ($m == 0)
      {
	return $n;
      }
    else
      {
	return RecursiveFindGCD (($n % $m), $m);
      }
}
