#!/usr/bin/perl -w

use strict 'vars';

use Getopt::Long;

my $WorkDir="./"; # working directory for storage of all output files
my $Path="1";
my $Machine="s2";   # choose label to use one of the preset paths
my $BuildDirectory;
my $MPIBuildDirectory;
my $MPIRun;
my $MPIRunPrefix="";
my $FTIPrograms="FTI/src/Programs/FTI/";
my $FCIPrograms="FTI/src/Programs/FCI/";
if ($Machine =~ m/s2/)
  {
    $BuildDirectory = "/scratch/ba340/DiagHam/build/";
    $MPIBuildDirectory = "/scratch/ba340/DiagHam/buildMPI/";
    $MPIRun ="mpirun";
  }
elsif ($Machine =~ m/g/)
  {
    $BuildDirectory = "/scratch/ba340/DiagHam/build/";
    $MPIBuildDirectory = "/scratch/ba340/DiagHam/buildMPI/";
    $MPIRun ="mpirun";
  }
elsif ($Machine =~ m/t/)
  {
    $BuildDirectory = "~/DiagHam/build/";
    $MPIBuildDirectory = "~/DiagHam/build/";
    $MPIRun ="mpirun";
  }
else
  {
    die ("Machine name not recognized\n");
  }

my $Processors=8;
my $UseDisk=0;
my $ExecOptions="";
my $ExecPrefix="";
my $MPIOptions="";

while( (defined($ARGV[0])&&$ARGV[0] =~ /^-/ ))
  {
    if ( $ARGV[0] =~ /--mpi/ )
      {
	my $NbrMPI;
	if (length($ARGV[0])>5)
	  {
	    $NbrMPI = $ARGV[0];
	    $NbrMPI =~ s/--mpi//;
	  }
	else
	  {
	    shift(@ARGV);
	    $NbrMPI = $ARGV[0];
	  }
	$ExecPrefix=$MPIRun." -np $NbrMPI ";
	$BuildDirectory=$MPIBuildDirectory;
	$MPIOptions = "--mpi-smp cluster.dat"
      }
    if ( $ARGV[0] =~ /--disk/ )
      {
	$UseDisk=1;
      }
    if ( $ARGV[0] =~ /-o/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $ExecOptions = $ARGV[0];
	    $ExecOptions =~ s/-P//;
	  }
	else
	  {
	    shift(@ARGV);
	    $ExecOptions = $ARGV[0];
	  }
      }
    if ( $ARGV[0] =~ /-P/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $Processors = $ARGV[0];
	    $Processors =~ s/-P//;
	  }
	else
	  {
	    shift(@ARGV);
	    $Processors = $ARGV[0];
	  }
      }
    shift(@ARGV);
  }
if (!defined $ARGV[0])
  {
    printf("No input file given.\n");
    printf("Input file syntax: 
GetEntanglementSpectrum [-P NbrProcessors] [-o ExecOptions] [--mpi mpiNbrThread] [--disk] TaskFile
where Taskfile has lines
spectrum.dat d kx-0 ky-0 ... kx-(d-1) ky-(d-1)
for subspaces with degeneracy d
NbrProcessors = number of proc to run on
ExecOptions = options to be used for calculation of eigenstates
mpiNbrThread = number of MPI instances\n"); 
    exit();
  }

# remaining files are data input
my $TaskFile = $ARGV[0];


my $ExecutableEntanglementParticlePartition =  $BuildDirectory.$FTIPrograms."FTIEntanglementEntropyParticlePartition";
my $ExecutableEntanglementSpectrum =  $BuildDirectory.$FTIPrograms."FTIEntanglementSpectrum";

# hand-coded options for entanglement spectrum
my $Options =  "--use-lapack --show-time  -S --processors $Processors";

unless (open (TASKFILE ,"$TaskFile"))
  {
    die ("can't open ".$TaskFile."\n");	
  }


my $TmpTask;
my @Params;
foreach $TmpTask (<TASKFILE>)
  {
    if (!($TmpTask =~ m/^#/))
      {
	chomp($TmpTask);
	@Params = split(" ", $TmpTask);

	if (scalar (@Params)<2)
	  {
	    print ("Skipping task $TmpTask\n");
	  }
	else
	  {
	    my $TmpSpectrum = $Params[0];
	    shift(@Params);
	    my $Degeneracy = $Params[0];
	    shift(@Params);

	    # define ground state sectors for simplicity
	    my @GroundStateSectorsX;
	    my @GroundStateSectorsY;
	
	    my $kx; 
	    my $ky;
	
	    my $GroundStates = $TmpSpectrum;
	    $GroundStates =~ s/\.dat/\.gs/;
	
	    open (GSFILE, ">$GroundStates");
	    my $count=0;
	    for (my $state=0; $state<$Degeneracy; ++$state) {
	      if (!defined($Params[0])) {
		die "incorrect degeneracy or number of k-points given for spectrum $TmpSpectrum\n";
	      }
	      $kx= $Params[0];
	      shift(@Params);

	      if (!defined($Params[0])) {
		die "incorrect degeneracy or number of k-points given for spectrum $TmpSpectrum\n";
	      }
	      $ky= $Params[0];
	      shift(@Params);

	      push(@GroundStateSectorsX,$kx);
	      push(@GroundStateSectorsY,$ky);
	    
	      if (($state>0)&&($kx==$GroundStateSectorsX[$state-1])&&($ky==$GroundStateSectorsY[$state-1])) {
		++$count;
	      } else {
		$count=0;
	      }

	      my $SectorFile =  $TmpSpectrum;
	      $SectorFile =~ s/\.dat/\_kx\_$kx\_ky\_$ky\.$count\.vec/;
	      print GSFILE ("$SectorFile\n");

	      if (! -e $SectorFile) {
		print("#Attention, file for ground state sector not available: ".$SectorFile."\n");

		# generate command
		&GenerateHofstadterCommand($TmpSpectrum, $kx, $ky, $count);
	      }
	    }
	
	    my $DensityMatrix = $TmpSpectrum;
	    $DensityMatrix =~ s/\.dat/\.full.parent/;
	    print ("$ExecPrefix$ExecutableEntanglementParticlePartition $MPIOptions $Options --degenerated-groundstate $GroundStates --density-matrix $DensityMatrix\n");
	    
	    close (GSFILE);
	  }
      }
  }



sub FindMinEnergy ()
  {
    my $FileName = $_[0];
    my $MinEnergy = $_[1];
    my $CurrentMinL = $_[2];
    unless (open(INFILE ,$FileName))
      {
	die ("can't open ".$FileName."\n");	
      }
    my $TmpLine = <INFILE>;
    chomp($TmpLine);    
    my @TmpArray = split (/ /, $TmpLine);
    $$MinEnergy = $TmpArray[1];
    $$CurrentMinL = $TmpArray[0];
    while (defined($TmpLine = <INFILE>))
      {
	chomp($TmpLine);    
	@TmpArray = split (/ /, $TmpLine);
	if ($$MinEnergy > $TmpArray[1])
	  {
	    $$MinEnergy = $TmpArray[1];
	    $$CurrentMinL = $TmpArray[0];
	  }
      }
    close (INFILE);
  }


sub FindGroundStateDegeneracy ()
  {
    my $FileName = $_[0];
    my $GroundStateEnergy = $_[1];
    my $GroundStateL = $_[2];
    my $DegeneracyError = $_[3];
    my $Degeneracy = 0;
    unless (open(INFILE ,$FileName))
      {
	die ("can't open ".$FileName."\n");	
      }
    my $TmpLine;
    while (defined($TmpLine = <INFILE>))
      {
	chomp($TmpLine);    
	my @TmpArray = split (/ /, $TmpLine);
	if ($TmpArray[0] == $GroundStateL)
	  {
	    if (((abs($TmpArray[1]) < $DegeneracyError) && (abs($GroundStateEnergy) < $DegeneracyError)) ||
		(abs($TmpArray[1] - $GroundStateEnergy) < ($DegeneracyError * abs($GroundStateEnergy))))
	      {
		$Degeneracy++;
	      }
	  }
      }
    close (INFILE);
    return $Degeneracy;
  }

# generate the command for calling FCIHofstadterModel
# $_[0] = spectrum name $TmpSpectrum
# $_[1] = $kx
# $_[2] = $ky
# $_[3] = vector index within momentum subsector

sub GenerateHofstadterCommand
  {
    my $HofstadterExecutable =  $BuildDirectory.$FCIPrograms."FCIHofstadterModel";

    my $SpectrumName = $_[0];
    my $KxSector = $_[1];
    my $KySector = $_[2];
    my $Index = $_[3];

    my $Statistics = "";

    if ($SpectrumName =~ m/boson/)
      {
	$Statistics = "--boson";
      }
    
    # geometry parameters
    $SpectrumName =~ /X\_(\d+)\_Y\_(\d+)\_q\_(\d+)\_n\_(\d+)\_x\_(\d+)\_y\_(\d*)\_/;
    my $FileX = $1;
    my $FileY = $2;
    my $Fileq = $3;
    my $Filen = $4;
    my $Filex = $5;
    my $Filey = $6;

    # interaction / projection
    my $Interaction;
    if ($SpectrumName =~ m/\_u\_\d+/)
      {
	$SpectrumName =~ /\_u\_(-*\d*[\.]*\d*)\_/;
	$Interaction = "-u $1";
      }
    else
      {
	$Interaction = "--flat-band";
      }
    

    # gauge fluxes
    $SpectrumName =~ /\_gx\_(-*\d*[\.]*\d*)\_gy\_(-*\d*[\.]*\d*)/;
    my $GaugeX="", my $GaugeY="";
    if (abs($1)>0.0)
      {
	$GaugeX = "--gauge-x $1";
      }
    if (abs($2)>0.0)
      {
	$GaugeY = "--gauge-y $1";
      }
    if ($Index==0)
      {
	if ($UseDisk>0)
	  {
	    print ("$ExecPrefix$HofstadterExecutable $MPIOptions $Statistics -X $FileX -Y $FileY -q $Fileq -p $Filen -x $Filex -y $Filey $Interaction $GaugeX $GaugeY $ExecOptions --only-kx $KxSector --only-ky $KySector --eigenstate -n ".($Index+1)." --fast-disk --show-itertime\n");
	  }
	else
	  {
	    print ("$ExecPrefix$HofstadterExecutable $MPIOptions $Statistics -X $FileX -Y $FileY -q $Fileq -p $Filen -x $Filex -y $Filey $Interaction $GaugeX $GaugeY $ExecOptions --only-kx $KxSector --only-ky $KySector --eigenstate -n ".($Index+1)." --force-reorthogonalize --show-itertime\n");
	  }
      }
    else
      {
	print ("$ExecPrefix$HofstadterExecutable $MPIOptions $Statistics -X $FileX -Y $FileY -q $Fileq -p $Filen -x $Filex -y $Filey $Interaction $GaugeX $GaugeY $ExecOptions --only-kx $KxSector --only-ky $KySector --eigenstate -n ".($Index+1)." --force-reorthogonalize --show-itertime\n");
      }
  }
