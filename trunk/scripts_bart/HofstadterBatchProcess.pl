#!/usr/bin/perl -w
#
# HofstadterBatchProcess.pl: run jobs for Hofstadter spectrum.
#
use strict 'vars';


# if (!defined($ARGV[0]))
#   {
#     print("usage HofstadterBatchProcess.pl -o "Options" BatchParameters.dat\n");

#     exit(1);
#   }

my $Machine="pc42";   # choose label to use one of the preset paths
my $BuildDirectory;
my $MPIBuildDirectory;
my $MPIRun;
my $MPIRunPrefix="";
my $FTIPrograms="FTI/src/Programs/FTI/";
my $FCIPrograms="FTI/src/Programs/FCI/";
if ($Machine =~ m/pc42/)
  {
    $BuildDirectory = "/scratch/ba340/DiagHam/build/";
    $MPIBuildDirectory = "/scratch/ba340/DiagHam/buildMPI/";
    $MPIRun ="mpirun";
  }
else
  {
    die ("Machine name not recognized\n");
  }

my $Processors=4;
my $UseDisk=0;
my $ExecOptions="";
my $ExecPrefix="";
my $MPIOptions="";
my $NbrState=1;
my $NbrMPI=0;
my $GetEigenstates=0;
my $UseBlockLanczos=0;
my $RunByN=1;

while( (defined($ARGV[0])&&$ARGV[0] =~ /^-/ ))
  {
    if ( $ARGV[0] =~ /--mpi/ )
      {
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
    if ( $ARGV[0] =~ /--run-by/ )
      {
	if (length($ARGV[0])>8)
	  {
	    $RunByN = $ARGV[0];
	    $RunByN =~ s/--run-by//;
	  }
	else
	  {
	    shift(@ARGV);
	    $RunByN = $ARGV[0];
	  }
      }
    if ( $ARGV[0] =~ /--block/ )
      {
	$UseBlockLanczos=1;
      }
    if ( $ARGV[0] =~ /--disk/ )
      {
	$UseDisk=1;
      }
    if ( $ARGV[0] =~ /--state/ )
      {
	$GetEigenstates=1;
      }
    if ( $ARGV[0] =~ /-n/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $NbrState = $ARGV[0];
	    $NbrState =~ s/-P//;
	  }
	else
	  {
	    shift(@ARGV);
	    $NbrState = $ARGV[0];
	  }
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


if (!defined($ARGV[0]))
  {
    print("usage HofstadterBatchProcess.pl [-P NbrProcessors] [-o ExecOptions] [--mpi mpiNbrThread] [--block] [--disk] [--state] [-n nbrState] -o \"Options\" BatchParameters.dat [...]
NbrProcessors = number of proc to run on
ExecOptions = options to be used for calculation of eigenstates
mpiNbrThread = number of MPI instances
block = use block-lanczos
disk = use disk storage
state = calculate eigenstates
nbrState = number of states to calculate per sector
BatchParameters can be generated by Mathematica notebook SquareHofstadterCells.nb\n");
    exit(1);
  }

# set executable
my $HofstadterExecutable =  $BuildDirectory.$FCIPrograms."FCIHofstadterModel";

my @ListFiles;
my $ReadMore = 1;
if ( defined($ARGV[0]) )
  {
    if ( -d$ARGV[0] )
      {
	chdir($ARGV[0]);
      }
    else # single file:
      {
	@ListFiles = $ARGV[0];
	$ReadMore=0;
      }
  }

my $Job=0;

my $TmpFile;
foreach $TmpFile (@ListFiles)
  {
     unless (open (TASKFILE ,"$TmpFile"))
       {
         die ("can't open ".$TmpFile."\n");	
       }
     my $TmpTask;
     my @Params;
     foreach $TmpTask (<TASKFILE>)
     {
        if (!($TmpTask =~ m/^#/))
        {
	  chomp($TmpTask);
	  &runHofstadterCode($TmpTask);
        }    
      }
    close(TASKFILE);
  }

# subroutine runHofstadterCode
# $_[0] Task parameters, as a string
sub runHofstadterCode
  {    
    my $TaskParameters = $_[0];
    my @Params = split("\t", $TaskParameters);
    if (scalar (@Params)<11)
      {
	print ("Skipping task $TaskParameters\n");
	return;
      }    
    # format: N\tN_s\tNx\tNy\tp\tq\tsgn\tX\tY\tx\ty
    
    my $N=$Params[0];
    my $p=$Params[4];
    my $q=$Params[5];
    my $sgn=$Params[6];
    my $X=$Params[7];
    my $Y=$Params[8];
    my $x=$Params[9];
    my $y=$Params[10];

    my $Command = "$ExecPrefix$HofstadterExecutable $MPIOptions -p $N -x $x -y $y -X $X -Y $Y -q $p --band-min 0 --band-max 0 --boson --u-potential 10 --show-itertime -n $NbrState";
    if ($NbrMPI==0)
      { $Command .= " -S --processors $Processors"; }
    if ($GetEigenstates>0)
      { $Command .= " --eigenstate"; }

    if ($NbrState==1)
      {
	if (($UseDisk>0) && ( ($GetEigenstates>0) || ($ExecOptions =~ m/--eigenstate/)))
	  {
	    $Command .= " --fast-disk";
	  }
      }
    else
      {
	if ($UseBlockLanczos)
	  {
	    if ($UseDisk>0)
	      {
		print("#Not using fast-disk mode for Block-Lanczos - Algorithm currently broken\n");
	      }
	    $Command .= " --block-lanczos --block-size $NbrState --use-lapack --force-reorthogonalize";
	  }
	else
	  {
	    $Command .= " --force-reorthogonalize";
	  }
      }
    { 
      use integer;
    if ($Job % $RunByN != $RunByN-1)
      {
	$Command .= " &";
      }
      print("$Command\n");
    }
    ++$Job;
  }