#!/usr/bin/perl -w
#
# script for optimizing overlaps for vectors associated with given protocol files
#
use strict 'vars';
use File::stat;
use Math::Complex;

# flag allowing to switch on some extra output
my $Verbose=0;

# hardwire which state to look at
my $Program="FQHELatticeBosonsGeneric";
my $MatrixProgram="MatrixElement";
my $Program64Suffix="_64";

my $DiceLattice="NbrSites = 6\n"
  ."Dimension = 2\n"
  ."LatticeVector0 = 1.73205080756888,0\n"
  ."LatticeVector1 = 0,3\n"
  ."SubLatticeVector0 = 0,0\n"
  ."SubLatticeVector1 = 0,1\n"
  ."SubLatticeVector2 = 0,2\n"
  ."SubLatticeVector3 = 0.866025403784439,0.5\n"
  ."SubLatticeVector4 = 0.866025403784439,1.5\n"
  ."SubLatticeVector5 = 0.866025403784439,2.5\n"
  ."NeighborsInCell = 0,1 | 1,2 | 1,3 | 1,4 | 2,5 | 4,5\n"
  ."NeighborCells = 0,1 | 1,0 | 0,-1 | -1,0 | 1,1 | -1,-1\n"
  ."NeighborsAcrossBoundary0_1 = 5,0 | 5,3\n"
  ."NeighborsAcrossBoundary0_-1 = 0,5 | 3,5\n"
  ."NeighborsAcrossBoundary1_1 = 5,0\n"
  ."NeighborsAcrossBoundary-1_-1 = 0,5\n"
  ."NeighborsAcrossBoundary1_0 = 3,1 | 4,1 | 5,2\n"
  ."NeighborsAcrossBoundary-1_0 = 1,3 | 1,4 | 2,5\n"
  ."UseGauge = yes\n"
  ."GaugeAyx = 1.15470053837925\n";

my $DiceLatticeReal="NbrSites = 6\n"
  ."Dimension = 2\n"
  ."LatticeVector0 = 1.73205080756888,0\n"
  ."LatticeVector1 = 0,3\n"
  ."SubLatticeVector0 = 0,0\n"
  ."SubLatticeVector1 = 0,1\n"
  ."SubLatticeVector2 = 0,2\n"
  ."SubLatticeVector3 = 0.866025403784439,0.5\n"
  ."SubLatticeVector4 = 0.866025403784439,1.5\n"
  ."SubLatticeVector5 = 0.866025403784439,2.5\n"
  ."NeighborsInCell = 0,1 | 1,2 | 1,3,1 | 1,4 | 2,5 | 4,5,1\n"
  ."NeighborCells = 0,1 | 1,0 | 0,-1 | -1,0 | 1,1 | -1,-1\n"
  ."NeighborsAcrossBoundary0_1 = 5,0 | 5,3\n"
  ."NeighborsAcrossBoundary0_-1 = 0,5 | 3,5\n"
  ."NeighborsAcrossBoundary1_1 = 5,0,1\n"
  ."NeighborsAcrossBoundary-1_-1 = 0,5,1\n"
  ."NeighborsAcrossBoundary1_0 = 3,1 | 4,1 | 5,2\n"
  ."NeighborsAcrossBoundary-1_0 = 1,3 | 1,4 | 2,5\n"
  ."UseGauge = no\n"
  ."NbrFlux = 3\n";

my $EffectiveTriangularLattice="Descriptor = eff_triang\n"
  ."NbrSites = 2\n"
  ."Dimension = 2\n"
  ."LatticeVector0 = 1.73205080756888,0\n"
  ."LatticeVector1 = 0,3\n"
  ."SubLatticeVector0 = 0,0\n"
  ."SubLatticeVector1 = 0.866025403784439,1.5\n"
  ."NeighborsInCell = 0,1\n"
  ."NeighborCells = 0,1 | 1,0 | 0,-1 | -1,0 | 1,1 | -1,-1\n"
  ."NeighborsAcrossBoundary0_1 = 1,0\n"
  ."NeighborsAcrossBoundary0_-1 = 0,1\n"
  ."NeighborsAcrossBoundary1_1 = 1,0\n"
  ."NeighborsAcrossBoundary-1_-1 = 0,1\n"
  ."NeighborsAcrossBoundary1_0 =  0,0 | 1,0 | 1,1\n"
  ."NeighborsAcrossBoundary-1_0 = 0,0 | 0,1 | 1,1\n"
  ."UseGauge = yes\n"
  ."GaugeAyx = 1.15470053837925\n";

my $Descriptor = "dice_doubled";
my $Trapping = 0.001;
my $RatioU = 1.0;

my $Directory="";
my @UnitCells;
$UnitCells[0]=3;
$UnitCells[1]=2;
my $DisplayHelp=0;
my $AllElements=0;
my $KeepFiles=0;
my $UseExtrapolateAnalytical=1;
my $UseRealRepresentation=0;
my $RealVectorString="";
my $NonZeroThreshold = 1e-6;
my $SolenoidCmd="";
my $SolenoidStr="";
my $HaveSolenoid=0;
my $SolenoidX=0;
my $SolenoidY=0;
my $VectorsOnly=0;

while( (defined($ARGV[0])&&$ARGV[0] =~ /^-/ ))
  {
    if ( $ARGV[0] =~ /-a/ )
      {
	$AllElements=1;
      }
    if ( $ARGV[0] =~ /-d/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $Directory = $ARGV[0];
	    $Directory =~ s/-d//;
	  }
	else
	  {
	    shift(@ARGV);
	    $Directory = $ARGV[0];
	  }
      }
    if ( $ARGV[0] =~ /-h/ )
      {
	$DisplayHelp=1;
      }
    if ( $ARGV[0] =~ /-k/ )
      {
	$KeepFiles=1;
      }
    if ( $ARGV[0] =~ /-v/ )
      {
	$VectorsOnly=1;
      }

    if ( $ARGV[0] =~ /-u/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $RatioU = $ARGV[0];
	    $RatioU =~ s/-u//;
	  }
	else
	  {
	    shift(@ARGV);
	    $RatioU = $ARGV[0];
	  }
      }
    if ( $ARGV[0] =~ /-z/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $NonZeroThreshold = $ARGV[0];
	    $NonZeroThreshold =~ s/-z//;
	  }
	else
	  {
	    shift(@ARGV);
	    $NonZeroThreshold = $ARGV[0];
	  }
      }
    if ( $ARGV[0] =~ /-r/ )
      {
	print("Using raw matrix elements with no rounding\n");
	$UseExtrapolateAnalytical=0;
      }
    if ( $ARGV[0] =~ /-p/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $Trapping = $ARGV[0];
	    $Trapping =~ s/-p//;
	  }
	else
	  {
	    shift(@ARGV);
	    $Trapping = $ARGV[0];
	  }
	print ("set trapping to ".$Trapping."\n");
      }
    if ( $ARGV[0] =~ /-C/ )
      {
	my $TmpStr;
	if (length($ARGV[0])>2)
	  {
	    $TmpStr = $ARGV[0];
	    $TmpStr =~ s/-C//;
	  }
	else
	  {
	    shift(@ARGV);
	    $TmpStr = $ARGV[0];
	  }
	@UnitCells = split(/,/,$TmpStr);	
	my $TmpLength = $#UnitCells+1;
	if ( $TmpLength != 2)
	  {
	    die ("need a lattice dimension of two: -C Lx,Ly\n");
	  }
	print ("Evaluating matrix elements for $UnitCells[0]x$UnitCells[1] unit cells\n");
      }
    if ( $ARGV[0] =~ /-R/ )
      {
	print("Using real representation of dice unit cell\n");
	$UseRealRepresentation=1;
	$RealVectorString="_re"
      }
    if ( $ARGV[0] =~ /-s/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $SolenoidCmd = $ARGV[0];
	    $SolenoidCmd =~ s/-s//;
	  }
	else
	  {
	    shift(@ARGV);
	    $SolenoidCmd = $ARGV[0];
	  }
	$HaveSolenoid=1;
	my @TmpArray = split(/,/,$SolenoidCmd);	
	my $TmpLength = $#TmpArray+1;
	if ( $TmpLength >= 1)
	  {
	    $SolenoidX=$TmpArray[0];
	  }
	if ( $TmpLength >= 2)
	  {
	    $SolenoidY=$TmpArray[1];
	  }
	if (($SolenoidX!=0.0)||($SolenoidY!=0.0))
	  {
	    $SolenoidCmd="-s ".$SolenoidCmd;
	    $SolenoidStr=sprintf ("_s_%g_%g", $SolenoidX, $SolenoidY);
	  }
	else
	  {
	    $SolenoidCmd="";
	    $SolenoidStr="";
	  }
	print ("Using solenoid fluxes $SolenoidStr $SolenoidCmd\n");
      }
    shift(@ARGV);
  }

if ($DisplayHelp)
  {
    print("usage DiceLatticeModel.pl -C Lx,Ly [-p V_p] [-d dir]\n");
    print("option -C: indicate length in x- and y-directions\n");
    print("       -p: strength of pinning potential (default: 0.001)\n");
    print("       -d: name of directory to generate files in (default: ./dice_Lx_Ly)\n");
    print("       -k: keep intermediate files on disk\n");
    print("       -u: ratio of U6 to U3 on 6-fold and 3-fold connected sites of dice lattice\n");
    print("       -a: include all matrix elements, ignoring obvious symmetries\n");
    print("       -r: write raw amplitudes without using rounding to expected WF amplitudes\n");
    print("       -R: use real gauge for dice lattice\n");
    print("       -s: use solenoid fluxes (s_x,s_y)\n");
    print("       -v: calculate vectors only\n");
    print("       -z: threshold above which elements are considered as non-zero\n");
    exit(1);
  }


# set default directory name, if not given
if ($Directory eq "")
  {
    $Directory = "./dice_".($UnitCells[0])."_".($UnitCells[1]);
  }

my $NbrCells = $UnitCells[0]*$UnitCells[1];
my $NbrSites = 2*$NbrCells;

# append lattice geometry to lattice definition
if ($UseRealRepresentation==1)
  {
    $DiceLattice=$DiceLatticeReal;
  }
$DiceLattice.="PeriodicRepeat = ".($UnitCells[0]).",".($UnitCells[1])."\n";
$EffectiveTriangularLattice.="PeriodicRepeat = ".($UnitCells[0]).",".($UnitCells[1])."\n";
$EffectiveTriangularLattice.="NbrFlux = ".($NbrSites/2)."\n";
my $Have64Bits=0;
my $tmp = "";
$tmp = `status`;
if ( $tmp =~ /x86_64/ )
  {
    $Program .= $Program64Suffix;
    $MatrixProgram .= $Program64Suffix;
  }

# create directory if not existent
if ( ! -e $Directory )
  {
    system ("mkdir -p $Directory");
  }

chdir($Directory);

#calculate localized single-particle states:
for (my $i=0; $i<$NbrSites; ++$i)
  {
    # generate lattice definition
    my $LatticeFile = "DiceDoubledPhases".$RealVectorString."_S".$i."_V_-".abs($Trapping)."_on_$UnitCells[0]x$UnitCells[1].dat";
    open (DEFINITION, ">$LatticeFile");
    my $CurrentDescriptor=$Descriptor.$RealVectorString."_S".$i."_V_-".abs($Trapping);
    print DEFINITION ("Descriptor = ".$CurrentDescriptor."\n");
    print DEFINITION ("LocalPotentials = ".GetSiteIndex($i).",-".abs($Trapping)."\n");
    print DEFINITION ($DiceLattice);
    close(DEFINITION);
    # run single particle calculation
    my $Command = "$Program --cmdlog-off --use-lapack -p 1 -L $LatticeFile -q ".(3*$NbrCells)." -c --eigenstate -n 1 $SolenoidCmd";
    if ($Verbose==1)
      {
	print ("running: $Command\n");
      }
    system($Command);
  }

#write definition of effective lattice
my $LatticeFile = "DiceDoubledEffective_$UnitCells[0]x$UnitCells[1].dat";
open (DEFINITION, ">$LatticeFile");
print DEFINITION ($EffectiveTriangularLattice);
close(DEFINITION);

#generate description of interaction, here: delta interaction
my $Interaction="HilbertSpaceDimension = ".(6*$NbrCells)."\n"
  ."NbrMatrixElements = 2\n"
  ."# Matrix element 0 -> U3, Matrix element 1 -> U6\n"
  ."MatrixElements = 1 1\n"
  ."Sparse = false\n"
  ."HaveOffDiagonal = false\n"
  ."DiagonalEntries = ";
for (my $i=0; $i<$NbrCells; ++$i)
  {
    $Interaction.="0 1 0 0 0 1 ";
  }
$Interaction.="\n";
my $InteractionFile = "DeltaInteraction_$UnitCells[0]x$UnitCells[1].dat";
open (DEFINITION, ">$InteractionFile");
print DEFINITION ($Interaction);
close (DEFINITION);


# test if we are done after calculating vectors
if ( $VectorsOnly == 1)
  {
    if ($KeepFiles==0)
      {
	system ("rm DiceDoubledPhases_S*_V_-".abs($Trapping)."_on_$UnitCells[0]x$UnitCells[1].dat");
	system ("rm bosons_lattice_dice_doubled_S*_V_-".abs($Trapping)."_".$UnitCells[0]."x".$UnitCells[1]."_n_1_hardcore".$SolenoidStr."_q_"
		.(3*$NbrCells).".dat");
      }
    print ("Done calculating vectors\n");
    exit(0);
  }

my $MatrixFile = "MatrixElements".$RealVectorString."_Delta_u_".$RatioU."_".$UnitCells[0]."x".$UnitCells[1].$SolenoidStr.".dat";
open (MATRIX, ">$MatrixFile");

# calculate all matrix elements
for (my $Index1=0; $Index1<$NbrSites; ++$Index1)
  {
    my $UpperBound2=$NbrSites;
    if ($AllElements==0)
      {
	$UpperBound2=$Index1+1;
      }
    for (my $Index2=0; $Index2<$UpperBound2; ++$Index2)
      {
	for (my $Index3=0; $Index3<$NbrSites; ++$Index3)
	  {
	    my $UpperBound4=$NbrSites;
	    if ($AllElements==0)
	      {
		$UpperBound4=$Index3+1;
	      }
	    for (my $Index4=0; $Index4<$UpperBound4; ++$Index4)
	      {
		# test whether matrix element can be non-zero:
		my $X1;
		my $Y1;
		my $X2;
		my $Y2;
		my $X3;
		my $Y3;
		my $X4;
		my $Y4;
		$X1 = ($Index1/2)%$UnitCells[0];
		$Y1 = ($Index1/2)/$UnitCells[0];
		$X2 = ($Index2/2)%$UnitCells[0];
		$Y2 = ($Index2/2)/$UnitCells[0];
		$X3 = ($Index3/2)%$UnitCells[0];
		$Y3 = ($Index3/2)/$UnitCells[0];
		$X4 = ($Index4/2)%$UnitCells[0];
		$Y4 = ($Index4/2)/$UnitCells[0];
		if ( ((abs($X1-$X2)<2)||(abs($X1-$X2)>$UnitCells[0]-2))
		     && ((abs($X1-$X3)<2)||(abs($X1-$X3)>$UnitCells[0]-2))
		     && ((abs($X1-$X4)<2)||(abs($X1-$X4)>$UnitCells[0]-2))
		     && ((abs($X2-$X3)<2)||(abs($X2-$X3)>$UnitCells[0]-2))
		     && ((abs($X2-$X4)<2)||(abs($X2-$X4)>$UnitCells[0]-2))
		     && ((abs($X3-$X4)<2)||(abs($X3-$X4)>$UnitCells[0]-2))
		     && ((abs($Y1-$Y2)<2)||(abs($Y1-$Y2)>$UnitCells[1]-2))
		     && ((abs($Y1-$Y3)<2)||(abs($Y1-$Y3)>$UnitCells[1]-2))
		     && ((abs($Y1-$Y4)<2)||(abs($Y1-$Y4)>$UnitCells[1]-2))
		     && ((abs($Y2-$Y3)<2)||(abs($Y2-$Y3)>$UnitCells[1]-2))
		     && ((abs($Y2-$Y4)<2)||(abs($Y2-$Y4)>$UnitCells[1]-2))
		     && ((abs($Y3-$Y4)<2)||(abs($Y3-$Y4)>$UnitCells[1]-2)) )
		  {
		    # call MatrixElement evaluator: 4,3: creation operators, 2,1 annihilation operators
		    my $Command = "$MatrixProgram --cmdlog-off -c --gauge --quiet --interaction $InteractionFile ".GetLocalWavefunction($Index4)." "
		      .GetLocalWavefunction($Index3)." ".GetLocalWavefunction($Index2)." ".GetLocalWavefunction($Index1);
		    my $Output = `$Command`;
		    my @OutputLines = split(/\n/,$Output);
		    chomp (@OutputLines);
		    my $TmpLength = $#OutputLines+1;

		    if ($Verbose==1)
		      {
			print ("Output for element $Index1 $Index2 $Index3 $Index4:\ngenerated by command: $Command");
			for (my $i=0; $i<$TmpLength; ++$i)
			  {
			    print ($OutputLines[$i]."\n");
			  }
		      }
		
		    #get matrix element for U3:
		    my @Complex3 = split(/ /,$OutputLines[0]);
		    if ($Verbose==1)
		      {
			print ("$OutputLines[0] read as $Complex3[0]+I*$Complex3[1]\n");
		      }
		    my $RoundedMultiples3="";
		    if ($UseExtrapolateAnalytical==1)
		      {
			ExtrapolateToAnalyticalValue($Complex3[0],$Complex3[1],$RoundedMultiples3);
		      }
		    #get matrix element for U6:
		    my @Complex6 = split(/ /,$OutputLines[1]);
		    my $RoundedMultiples6="";
		    if ($UseExtrapolateAnalytical==1)
		      {
			ExtrapolateToAnalyticalValue($Complex6[0],$Complex6[1],$RoundedMultiples6);
		      }
		    my $Real = $Complex3[0]+$RatioU*$Complex6[0];
		    if (abs($Real)<1e-8)
		      {
			$Real=0.0;
		      }
		    my $Imag = $Complex3[1]+$RatioU*$Complex6[1];
		    if (abs($Imag)<1e-8)
		      {
			$Imag=0.0;
		      }
		    if (((abs($Real)>0.0)||(abs($Imag)>0.0))&&($Real*$Real+$Imag*$Imag>1e-12))
		      {
			printf("$Index1 $Index2 $Index3 $Index4 ".$Real." ".$Imag." "
			       .$RoundedMultiples3." ".$RoundedMultiples6."\n");
			printf MATRIX ("$Index1 $Index2 $Index3 $Index4 ".$Real." ".$Imag."\n"); #." ".$RoundedMultiples3." ".$RoundedMultiples6."\n");
		      }
		  }
	      }
	  }
      }
  }
close (MATRIX);

if ($KeepFiles==0)
  {
    system ("rm DiceDoubledPhases_S*_V_-".abs($Trapping)."_on_$UnitCells[0]x$UnitCells[1].dat");
    system ("rm bosons_lattice_dice_doubled".$RealVectorString."_S*_V_-".abs($Trapping)."_".$UnitCells[0]."x".$UnitCells[1]."_n_1_hardcore".$SolenoidStr."_q_"
	    .(3*$NbrCells)."*");
  }


# get coordinate of trapping site
#
# $_[0] = index of trapping site

sub GetSiteIndex
  {
    my $SiteIndex = $_[0];
    if ( $SiteIndex % 2 == 0)
      {
	return (6*($SiteIndex/2)+1);
      }
    else
      {
	return (6*(($SiteIndex-1)/2)+5)
      }
  }

# get filename of wavefunction at given trapping site
#
# $_[0] = index of trapping site

sub GetLocalWavefunction
  {
    my $SiteIndex = $_[0];
    return "bosons_lattice_dice_doubled".$RealVectorString."_S".$SiteIndex."_V_-".abs($Trapping)."_".$UnitCells[0]."x".$UnitCells[1]."_n_1_hardcore".$SolenoidStr."_q_"
      .(3*$NbrCells).".0.vec";
  }

# try to guess a known exact value from an approximate numerical matrix element
#
# $_[0] = real part of matrix element
# $_[1] = imaginary part of matrix element
# $_[2] = multiplicity that was guessed
sub ExtrapolateToAnalyticalValue
  {
    my $Re=$_[0];
    my $Im=$_[1];
    my @KnownPrefactors;
    # global prefactor
    my $GlobalPrefactor = 1.0/144;
    # insert known additional prefactors here:
    my $NbrKnown = push (@KnownPrefactors,1.0);
    $NbrKnown = push (@KnownPrefactors,sqrt(2.0));
    if (($HaveSolenoid==0)&&(abs($Im)>1e-5)&&(abs(abs($Re/$Im)-1.0)<1e-5))
      {
	if ($Re*$Im>0.0)
	  {
	    $Re=$Im;
	  }
	else
	  {
	    $Re=-$Im;
	  }
      }
    my $AbsVal = sqrt($Re*$Re+$Im*$Im);
    if ($AbsVal>$NonZeroThreshold)
      {
	my $RoundedMultiples=0;
	my $UsedPrefactor=-1;
	for (my $NbrValue=0; ($NbrValue<$NbrKnown) && ($UsedPrefactor<0); ++$NbrValue)
	  {
	    $RoundedMultiples = sprintf("%.0f", $AbsVal/($GlobalPrefactor*$KnownPrefactors[$NbrValue]));
	    my $Factor = $RoundedMultiples*($GlobalPrefactor*$KnownPrefactors[$NbrValue])/$AbsVal;
	    if (abs($Factor-1.0)<0.03)
	      {
		$Re*=$Factor;
		$Im*=$Factor;
		$UsedPrefactor = $NbrValue;
	      }
	  }
	if ($UsedPrefactor<0)
	  {
	    print("Did not find appropriate factor for matrix element ".$Re."+I*".$Im."\n");
	    exit;
	  }
	$_[0]=$Re;
	$_[1]=$Im;
	$_[2]=$KnownPrefactors[$UsedPrefactor]*$RoundedMultiples;
      }
    else
      {
	$_[0]=0.0;
	$_[1]=0.0;
	$_[2]=0;
      }
  }
