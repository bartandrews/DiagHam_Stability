#!/usr/bin/perl -w
#
# script for optimizing overlaps for vectors associated with given protocol files
#
use strict 'vars';
use File::stat;

# hardwire which state to look at
my $attach=+1;  # direction of flux attachement +/- 1
my $attachStr="";
#my $Program_32="/rscratch/gm360/bin/FQHELatticeDensityMatrix";
my $Program_32="FQHELatticeCompositeFermions";
my $Program_64="FQHELatticeCompositeFermions_64";

if (!defined($ARGV[0]))
  {
    print("usage EvaluateBosonsOnLattice.pl [-m] (protocol files)\n");
    exit(1);
  }

my $attachCmd="-f 1";

while( $ARGV[0] =~ /^-/ )
  {
    if ( $ARGV[0] =~ /-m/ )
      {
	$attach=-1;
	$attachStr="-1";
	$attachCmd="-f -1";
	print ("Calculating for negative flux attachment\n");
      }
    shift(@ARGV);
  }

my $Program;
my $Memory=0;
my $Have64Bits=0;
my $tmp = `status`;
if ( $tmp =~ /x86_64/ )
  {
    $Program = $Program_64;
  }
else
  {
    $Program = $Program_32;	
  }

my $TmpFile;
foreach $TmpFile (@ARGV)
  {
    if ($TmpFile =~ /bosons\_lattice.*.dat/)
      {
	print ("Analyzing vectors for protocol ".$TmpFile."\n");
	&AnalyzeVectors($TmpFile);
      }
  }


# test for available eigenvectors, and analyze them.
#
# $_[0] = spectrum file name

sub AnalyzeVectors
  {
    my $FileName = $_[0];
    my $HardCore;
    my $N;
    my $x;
    my $y;
    my $u;
    my $q;
    if ($FileName =~ m/hardcore/)
      {

	$FileName =~ /n\_(\d+)\_x\_(\d*)\_y\_(\d*)\_/;
	$N = $1;
	$x = $2;
	$y = $3;
	$u = 0;
	$HardCore=1;
	$q = -1;
      }
    else
      {
	$FileName =~ /n\_(\d+)\_x\_(\d*)\_y\_(\d*)\_u\_(-*\d*[\.]*\d*)\_/;
	$N = $1;
	$x = $2;
	$y = $3;
	$u = $4;
	$HardCore=0;
	$q = -1;
      }
    my $TotalSolenoid="";
    my $SolenoidX=0.0;
    my $SolenoidY=0.0;
    if ($FileName =~ m/\_s\_/)
      {	
	$FileName =~ /\_s\_(-*\d*[\.]*\d*e*-*\d*)\_(-*\d*[\.]*\d*e*-*\d*)/;
	$SolenoidX = $1;
	$SolenoidY = $2;
	$TotalSolenoid="--solenoid-flux $SolenoidX,$SolenoidY";
      }
    my $BaseName = $FileName;
    if ($FileName =~ /bosons\_lattice.*\_q\_(\d*).dat/)
      {
	$q = $1;
	$BaseName =~ s/q\_$q.dat/q/;	
      }
    else
      {
	$BaseName =~ s/q.dat/q/;	
      }
    my $Interaction;
    if ( $HardCore == 1)
      {
	print "n = ".$N."  x = ".$x."  y = ".$y."  (hardcore bosons)\n";
	$Interaction ="-c";
      }
    else
      {
	print "n = ".$N."  x = ".$x."  y = ".$y."  u = ".$u."\n";
	$Interaction ="";
      }
    my $MinQ;
    my $MaxQ;
    if ( $q==-1 )
      {
	$MinQ=0;
	$MaxQ=$x*$y;	
      }
    else
      {
	$MinQ=$q;
	$MaxQ=$q;
      }
    my $TmpFileName = "tmp".time().".p";
    for ($q=$MinQ; $q<=$MaxQ;++$q)
      {
	print ("Analyzing q=".$q."\n");
	system ("grep ^$q $FileName | sed -e 's/$q //'> $TmpFileName");
	
	open (INFILE, $TmpFileName);
	my $TmpLine;
	my $CountEV=0;
	my $VectorFile;
	my @EigenValues;
	my @EigenVectors;
	my $ProtocolName = $BaseName."\_$q\.Ovl".$attachStr;	
	open (OUTFILE, ">$ProtocolName");
	print OUTFILE ("# ID\tE\tO\tSolenoid_x\tSolenoid_y\n");
	foreach $TmpLine (<INFILE>)
	  {
	    $VectorFile = $BaseName."\_$q\.$CountEV\.vec";
	    if ( -e $VectorFile )
	      {
		chomp($TmpLine);
		my @Tmp = split (/ /, $TmpLine);
		push (@EigenValues, $Tmp[0]);
		push (@EigenVectors, $VectorFile);
	      }
	    ++$CountEV;
	  }
	close(INFILE);
	print ("found ".($#EigenVectors+1)." vector files:\n");
	if ($#EigenVectors > -1)
	  {
	    for ($CountEV=0;$CountEV<$#EigenVectors+1;++$CountEV)
	      {
		my $ProtocolName2 = $BaseName."\_$q.$CountEV.ovl$attachStr";		
		if (! -e $ProtocolName2 )
		  {
		    system("$Program -p $N -q $q -x $x -y $y $Interaction $TotalSolenoid --optimize ".$EigenVectors[$CountEV]." > $ProtocolName2");
		    print("running: $Program -p $N -q $q -x $x -y $y $Interaction $attachCmd $TotalSolenoid --optimize ".$EigenVectors[$CountEV]." > $ProtocolName2\n");
		  }
		else
		  {
		    print("found existing: N=$N, N_phi=$q x=$x y=$y $Interaction for EV ".$EigenVectors[$CountEV]. "($ProtocolName2)\n");
		  }
		my $ResultStr = `grep ^Final $ProtocolName2`;
		$ResultStr =~ /Final overlap at theta=\( (-*\d*[\.]*\d*e*-*\d*),(-*\d*[\.]*\d*e*-*\d*) \), O=(-*\d*[\.]*\d*e*-*\d*)/;
		my $Sx = $1;
		my $Sy = $2;
		my $Ovl = $3;
		print OUTFILE ("$CountEV\t$EigenValues[$CountEV]\t$Ovl\t$Sx\t$Sy\n");
	      }
	    close(OUTFILE);
	  }
	else
	  {
	    # no vector files present: clean up and delete otherwise empty protocol file!
	    close(OUTFILE);
	    system("rm $ProtocolName");
	  }
      }
    system ("rm $TmpFileName");
  }


