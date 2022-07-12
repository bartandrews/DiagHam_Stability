#!/usr/bin/perl -w
#
# script for optimizing overlaps for vectors associated with given protocol files
#
use strict 'vars';
use File::stat;
use Math::Complex;


# hardwire which state to look at
my $Program_32="FQHELatticeBosons";
my $Program_64="FQHELatticeBosons_64";
my $OverlapExe="GenericOverlap";

my $CalculateVectors=0;
my $NbrGrid=20;
my $GridString="";
my $GridString2="";
my $ReferenceString="0,0,1,1";
my $Degeneracy=1;
my @Multiplet;
$Multiplet[0]=0;
my $Memory=1000;
my $Options="";
my $NbrCalculate=1;

while( (defined($ARGV[0])&&$ARGV[0] =~ /^-/ ))
  {
    if ( $ARGV[0] =~ /-d/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $Degeneracy = $ARGV[0];
	    $Degeneracy =~ s/-d//;
	  }
	else
	  {
	    shift(@ARGV);
	    $Degeneracy = $ARGV[0];
	  }
	for (my $i=1; $i<$Degeneracy; ++$i)
	  {
	    $Multiplet[$i]=$i;
	  }
      }
    if ( $ARGV[0] =~ /-g/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $GridString = $ARGV[0];
	    $GridString =~ s/-g//;
	  }
	else
	  {
	    shift(@ARGV);
	    $GridString = $ARGV[0];
	  }
      }
    if ( $ARGV[0] =~ /-h/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $GridString2 = $ARGV[0];
	    $GridString2 =~ s/-h//;
	  }
	else
	  {
	    shift(@ARGV);
	    $GridString2 = $ARGV[0];
	  }
      }
    if ( $ARGV[0] =~ /-m/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $Memory = $ARGV[0];
	    $Memory =~ s/-m//;
	  }
	else
	  {
	    shift(@ARGV);
	    $Memory = $ARGV[0];
	  }
      }
    if ( $ARGV[0] =~ /-n/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $NbrGrid = $ARGV[0];
	    $NbrGrid =~ s/-n//;
	  }
	else
	  {
	    shift(@ARGV);
	    $NbrGrid = $ARGV[0];
	  }
      }
    if ( $ARGV[0] =~ /-o/ )
      {
        if (length($ARGV[0])>2)
          {
            $Options = $ARGV[0];
            $Options =~ s/-o//;
          }
        else
          {
            shift(@ARGV);
            $Options = $ARGV[0];
          }
      }

    if ( $ARGV[0] =~ /-q/ )
      {
	my $TmpStr;
	if (length($ARGV[0])>2)
	  {
	    $TmpStr = $ARGV[0];
	    $TmpStr =~ s/-q//;
	  }
	else
	  {
	    shift(@ARGV);
	    $TmpStr = $ARGV[0];
	  }
	@Multiplet = split(/,/,$TmpStr);	
	sort(@Multiplet);
	$Degeneracy = $#Multiplet+1;
	print ("Analysing multiplet  [$Multiplet[0]");
	for (my $i=1; $i<$Degeneracy; ++$i)
	  {
	    print (", $Multiplet[$i]");
	  }
	print ("]\n");
	print ("Multiplet degeneracy: ".$Degeneracy."\n");
      }
    if ( $ARGV[0] =~ /-r/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $ReferenceString = $ARGV[0];
	    $ReferenceString =~ s/-r//;
	  }
	else
	  {
	    shift(@ARGV);
	    $ReferenceString = $ARGV[0];
	  }
      }
    if ( $ARGV[0] =~ /-s/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $NbrCalculate = $ARGV[0];
	    $NbrCalculate =~ s/-s//;
	  }
	else
	  {
	    shift(@ARGV);
	    $NbrCalculate = $ARGV[0];
	  }
      }
    if ( $ARGV[0] =~ /-c/ )
      {
	$CalculateVectors=1;
	print("Will calculate missing vectors!\n");
      }
    shift(@ARGV);
  }

if (!defined($ARGV[0]))
  {
    print("usage ChernNumberLattice.pl [-g GRIDPOINTS] [-n NbrPoints] [-c] [-r s1x,s1y,s2x,s2y] basename_*SX_SY*\n");
    print("option -g: list of discrete points used for both x- and y- directions\n");
    print("       -h: list of discrete points used for both y- directions (if different from x)\n");
    print("       -d: degeneracy of groundstate multiplet (or indicate states with -q)\n");
    print("       -n: number of gridpoints to be used\n");
    print("       -r: reference points A,B as theta1A,theta2A,theta1B,theta2B (default 0,0,1,1)\n");
    print("       -c: optionally calculate missing vector files\n");
    print("       -m: memory for precalculations when calculating vectors\n");
    print("       -q: quantum numbers of states to be considered part of multiplet (-q q1,q2,q3,...)\n");
    print("       -s: number of states to be calculated at each point\n");
    exit(1);
  }


my @GridPointsX;
my @GridPointsY;
my $NbrGridY;

if (length($GridString)>0)
  {
    @GridPointsX = split (/,/, $GridString);
    $NbrGrid = $#GridPointsX + 1;
  }
else
  {
    my $Sep=2.0/$NbrGrid;
    for (my $i=0; $i<$NbrGrid; ++$i)
      {
	$GridPointsX[$i]=-1.0+($i+1)*$Sep;
      }
  }
if ( length($GridString2)>0)
  {
    @GridPointsY = split (/,/, $GridString2);
    $NbrGridY = $#GridPointsY + 1;
  }
else
  {
    @GridPointsY=@GridPointsX;
    $NbrGridY=$NbrGrid;
  }

if ($NbrCalculate<$Multiplet[$#Multiplet]+1)
  {
    $NbrCalculate=$Multiplet[$#Multiplet]+1;
  }

if ($NbrCalculate<$Degeneracy)
  {
    $NbrCalculate=$Degeneracy;
  }

my $Program;
my $Have64Bits=0;
my $tmp = "";
$tmp = `status`;
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
    if ($TmpFile =~ m/bosons\_lattice/)
      {
	print ("Analyzing states with base name ".$TmpFile."\n");
	&AnalyzeChern($TmpFile);
      }
  }


# calculate gauge and prepare output for plotting
#
# $_[0] = base file name

sub AnalyzeChern
  {
    my $BaseName = $_[0];
    my $HardCore;
    my $N;
    my $x;
    my $y;
    my $u;
    my $q;
    if ($BaseName =~ m/hardcore/)
      {

	$BaseName =~ /n\_(\d+)\_x\_(\d*)\_y\_(\d*)\_.*\_q\_(\d*)/;
	$N = $1;
	$x = $2;
	$y = $3;
	$u = 0;
	$HardCore=1;
	$q = $4;
      }
    else
      {
	$BaseName =~ /n\_(\d+)\_x\_(\d*)\_y\_(\d*)\_u\_(-*\d*[\.]*\d*)\_.*\_q\_(\d*)/;
	$N = $1;
	$x = $2;
	$y = $3;
	$u = $4;
	$HardCore=0;
	$q = $5;
      }
    my $TotalSolenoid="";
    my $Interaction;
    if ( $HardCore == 1)
      {
	print "n = ".$N."  x = ".$x."  y = ".$y."  q = ".$q."  (hardcore bosons)\n";
	$Interaction ="-c";
      }
    else
      {
	print "n = ".$N."  x = ".$x."  y = ".$y."  q = ".$q."  u = ".$u."\n";
	$Interaction ="";
      }
    
    my @ReferenceVals = split (/,/, $ReferenceString);

    my $RefS1x = $ReferenceVals[0];
    my $RefS1y = $ReferenceVals[1];
    my $RefS2x = $ReferenceVals[2];
    my $RefS2y = $ReferenceVals[3];

    my $CommandLine = "$Program -p $N -x $x -y $y $Interaction -q $q -n $NbrCalculate -m $Memory $Options";

    TestVectors ($BaseName, $RefS1x, $RefS1y, $Degeneracy, $CalculateVectors, $CommandLine, \@Multiplet);
    TestVectors ($BaseName, $RefS2x, $RefS2y, $Degeneracy, $CalculateVectors, $CommandLine, \@Multiplet);

    my $LogFileName = GetVectorName($BaseName, 0, 0, 0);
    $LogFileName =~ s/0\.vec/cn/;

    SavePrevious($LogFileName);
    open (LOGFILE, ">$LogFileName");
    
    my $PlotFileName = GetVectorName($BaseName, 0, 0, 0);
    $PlotFileName =~ s/0\.vec/cnp1/;
    open (PLOTFILE1, ">$PlotFileName");
    
    $PlotFileName =~ s/cnp1/cnp2/;
    open (PLOTFILE2, ">$PlotFileName");
    
    print LOGFILE ("# Evaluation of Chern-Number for state at N = $N, x = $x, y = $y, q = $q  ($Interaction)\n");
    print LOGFILE ("# Theta_x\tTheta_y\tOmega_x\tOmega_y\n");

    my $GridX;
    my $GridY;
    for ($GridX=0; $GridX<$NbrGrid;++$GridX)
      {
	my $SolenoidX = $GridPointsX[$GridX];
	printf ("%8g\t",$SolenoidX);
	for ($GridY=0; $GridY<$NbrGridY;++$GridY)
	  {
	    my $SolenoidY = $GridPointsY[$GridY];
	    #print ("Analysing point $SolenoidX,$SolenoidY \n");
	    print ("x");

	    # make sure we have all necessary data
	    TestVectors ($BaseName, $SolenoidX, $SolenoidY, $Degeneracy, $CalculateVectors, $CommandLine, \@Multiplet);
	    my $Alpha;
	    my $Beta;
	    my $OvlCommand="GenericOverlap --quiet -s -c ";
	    my @Matrix1;
	    my @Matrix2;
	    my @Overlaps;
	    for ($Alpha=0; $Alpha<$Degeneracy;++$Alpha)
	      {
		my @Column1;
		my @Column2;
		my $VectorA1=GetVectorName($BaseName, $RefS1x, $RefS1y, $Multiplet[$Alpha]);
		my $VectorA2=GetVectorName($BaseName, $RefS2x, $RefS2y, $Multiplet[$Alpha]);
		for ($Beta=0; $Beta<$Degeneracy;++$Beta)
		  {
		    my $VectorB=GetVectorName($BaseName, $SolenoidX, $SolenoidY, $Multiplet[$Beta]);
		    my $TmpCmd = $OvlCommand." ".$VectorA1." ".$VectorB;
		    my $OvlString= `$TmpCmd`;
		    chomp($OvlString);
		    #print("run $TmpCmd\n$OvlString\n");
		    @Overlaps = split(/ /,$OvlString);
		    my $z = cplx($Overlaps[0],$Overlaps[1]);
		    push(@Column1,$z);
		    $TmpCmd = $OvlCommand." ".$VectorA2." ".$VectorB;
		    $OvlString= `$TmpCmd`;
		    chomp($OvlString);
		    #print("run $TmpCmd\n$OvlString\n");
		    @Overlaps = split(/ /,$OvlString);
		    $z = cplx($Overlaps[0],$Overlaps[1]);
		    push(@Column2,$z);
		  }
		push(@Matrix1,\@Column1);
		push(@Matrix2,\@Column2);
	      }
	    my $Det1 = det (\@Matrix1);
	    my $Det2 = det (\@Matrix2);
	    my $Arg1 = arg($Det1);
	    my $Arg2 = arg($Det2);
	    my $Sqr1 = abs($Det1*$Det1);
	    my $Sqr2 = abs($Det2*$Det2);
	    my $TotalX = sin($Arg1-$Arg2);
	    my $TotalY = cos($Arg1-$Arg2);
	    print LOGFILE ("$SolenoidX\t$SolenoidY\t$TotalX\t$TotalY\t$Sqr1\t$Sqr2\n");
	    print PLOTFILE1 ("$SolenoidX\t$SolenoidY\t$Sqr1\n");
	    print PLOTFILE2 ("$SolenoidX\t$SolenoidY\t$Sqr2\n");
	  }
	print PLOTFILE1 ("\n");
	print PLOTFILE2 ("\n");
	print ("\n");
      }
    close(LOGFILE);
  } # end of AnalyzeChern



# test if all required vectors at a given point are present
# if not, calculate
sub TestVectors {
  my $BaseName = $_[0];
  my $SolenoidX = $_[1];
  my $SolenoidY = $_[2];
  my $Degeneracy = $_[3];
  my $Calculate = $_[4];
  my $Command = $_[5];
  my $Multiplet = $_[6];

  while ($SolenoidX<=-1.0)
    {
      $SolenoidX+=2.0;
    }
  while ($SolenoidX>2.0)
    {
      $SolenoidX-=2.0;
    }
  while ($SolenoidY<=-1.0)
    {
      $SolenoidY+=2.0;
    }
  while ($SolenoidY>2.0)
    {
      $SolenoidY-=2.0;
    }
  
  my $HavePoint=1;
  for (my $i=0; $i<$Degeneracy; ++$i)
    {
      my $VectorName = GetVectorName($BaseName, $SolenoidX, $SolenoidY, $Multiplet->[$i]);
      if ( ! -e $VectorName )
	{
	  print ("State $VectorName not found!\n");
	  $HavePoint=0;
	}
#      else
#	{
#	  print ("State $VectorName found!\n");
#	}
    }
  if ( $HavePoint == 0)
    {
      if ( $Calculate == 1 )
	{
	  print ("Missing vectors at $SolenoidX,$SolenoidY ... recalculating\n");
	  my $Instruction = sprintf("%s --eigenstate --show-itertime --solenoid-flux %g,%g", $Command, $SolenoidX,$SolenoidY);
	  print ("Command executed: ".$Instruction."\n");

	  system($Instruction);
	}
      else
	{
	  die("Missing vectors at $SolenoidX,$SolenoidY ...\nAborting - please supply the missing files manually!\n");
	}
    }
}

# generate the name of a vector-file
#
sub GetVectorName
  {
    my $BaseName=$_[0];
    my $SolenoidX = $_[1];
    my $SolenoidY = $_[2];
    my $ID = $_[3];

    while ($SolenoidX<=-1.0)
      {
	$SolenoidX+=2.0;
      }
    while ($SolenoidX>2.0)
      {
	$SolenoidX-=2.0;
      }
    while ($SolenoidY<=-1.0)
      {
	$SolenoidY+=2.0;
      }
    while ($SolenoidY>2.0)
      {
	$SolenoidY-=2.0;
      }

    my $VectorName = $BaseName.".".$ID.".vec";
    my $XStr = sprintf("%g",$SolenoidX);
    my $YStr = sprintf("%g",$SolenoidY);
    $VectorName =~ s/XX/$XStr/;
    $VectorName =~ s/YY/$YStr/;
    if (($SolenoidX==0)&&($SolenoidY==0))
      {
	$VectorName =~ s/\_s\_0\_0//;
      }
    return $VectorName;
  }

#  see if a file with the requested name already exists. If so, rename it to increase a counter
#
sub SavePrevious {
  my $MaxSave = 5;   # parameter that can be adjusted to manage number of files kept on record
  my $FileName = $_[0];
  my $EffectiveFileName = $FileName;
  my $Level = -1;
  if ( defined($_[1]) )
    {
      $Level = $_[1];
      $EffectiveFileName = $FileName.".".$Level;
    }
  print("SavePrevious ($FileName,$Level) -> effective name : $EffectiveFileName \n");
  if ( -e $EffectiveFileName )
    {
      if ( $Level < $MaxSave )
	{
	  SavePrevious($FileName,$Level+1);
	}
      print("move $EffectiveFileName $FileName".".".($Level+1)."\n");
      system ("mv $EffectiveFileName $FileName".".".($Level+1));
    }
}

# calculate the determinant of a matrix given as a double list
# call this function by reference: det (\@matrix)
sub det {
    my $matrix = shift;
    my $size   = $#{ $matrix } + 1;

    foreach (@$matrix) {
#        print ("Line ".$_.", size ".($#{ $_ } + 1)."\n");
#	foreach (@$_)
#	  {
#	    print("Entry".$_."\n");
#	  }
        die "det(Matrix) requires n x n matrix!" if @$_ != $size;
        }

    return $matrix->[0][0] if $size == 1;
    return $matrix->[0][0] * $matrix->[1][1] - $matrix->[1][0] * $matrix->[0][1]
      if $size == 2;
    return _det_helper( $matrix, $size );
}

sub _det_helper {
    my $matrix = shift;
    my $size   = shift;

    return $matrix->[0][0] * $matrix->[1][1] * $matrix->[2][2] + $matrix->[1][0]
      * $matrix->[2][1] * $matrix->[0][2] + $matrix->[2][0] * $matrix->[0][1] *
      $matrix->[1][2] - $matrix->[0][2] * $matrix->[1][1] * $matrix->[2][0] -
      $matrix->[1][2] * $matrix->[2][1] * $matrix->[0][0] - $matrix->[2][2] *
      $matrix->[0][1] * $matrix->[1][0]
      if $size == 3;

    my $det;
    foreach ( 0 .. $size - 1 ) {
        if ( $_ % 2 ) {
            $det -=
              $matrix->[0][$_] *
              _det_helper( _matrix_slice( $matrix, 0, $_ ), $size - 1 );
        }
        else {
            $det +=
              $matrix->[0][$_] *
              _det_helper( _matrix_slice( $matrix, 0, $_ ), $size - 1 );
        }
    }
    return $det;
}

sub _matrix_slice {
    my $matrix = shift;
    my $x      = shift;
    my $y      = shift;

    return [ map { [ @{$_}[ 0 .. $y - 1, $y + 1 ... $#$_ ] ] }
          @{$matrix}[ 0 .. $x - 1, $x + 1 .. $#$matrix ] ];
}
