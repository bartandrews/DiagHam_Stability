#!/usr/bin/perl -w
#
# script for optimizing overlaps for vectors associated with given protocol files
#
use strict 'vars';
use File::stat;
use Math::Complex;


# hardwire which state to look at
my $Program="FQHESphereImportVector";
my $Extension="";

my $Format="";
my $OutputFormat="";
my $Reverse="--reverse";
my $InputType="--raw-state";
my $TargetDirectory="./";
my $DisplayHelp=0;

while( (defined($ARGV[0])&&$ARGV[0] =~ /^-/ ))
  {
    if ( $ARGV[0] =~ /-d/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $TargetDirectory = $ARGV[0];
	    $TargetDirectory =~ s/-d//;
	  }
	else
	  {
	    shift(@ARGV);
	    $TargetDirectory = $ARGV[0];
	  }
	print ("As of yet, the output directory cannot be changed from ./");
      }
    if ( $ARGV[0] =~ /-f/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $Format = $ARGV[0];
	    $Format =~ s/-f//;
	  }
	else
	  {
	    shift(@ARGV);
	    $Format = $ARGV[0];
	  }
      }
    if ( $ARGV[0] =~ /-o/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $OutputFormat = $ARGV[0];
	    $OutputFormat =~ s/-o//;
	  }
	else
	  {
	    shift(@ARGV);
	    $OutputFormat = $ARGV[0];
	  }
      }
    if ( $ARGV[0] =~ /-n/ )
      {
	$Reverse = "";
      }
    if ( $ARGV[0] =~ /-a/ )
      {
	$InputType="";
      }
    if ( $ARGV[0] =~ /-h/ )
      {
	$DisplayHelp=1;
      }
    shift(@ARGV);
  }

if (!defined($ARGV[0]))
  {
    print("usage ImportAssistant.pl -f format [-a -d target_directory -n -h] INPUT-VECTORS\n");
    print("option -d: option to indicate target directory\n");
    print("       -a: input file in ASCII format\n");
    print("       -f: format of input filename, containing directives NN and SS for particle and flux numbers\n");
    print("       -h: display help\n");
    print("       -n: no reversal of order in vector\n");
    exit(1);
  }

if ($Format eq "")
  {
    print("Error: please indicate a format for the input vectors\n");
    exit(1);
  }

my $tmp = "";
my $tmp2 = `which status`;

if (length($tmp2)>0)
  {
    $tmp = `status`;
    if ( $tmp =~ /x86_64/ )
      {
	$Extension = "_64";
      }
  }

my $TmpFile;
foreach $TmpFile (@ARGV)
  {
    print ("Importing vector ".$TmpFile."\n");
    &ImportVector($TmpFile, $Format, $OutputFormat);
  }

# calculate gauge and prepare output for plotting
#
# $_[0] = base file name

sub ImportVector
  {
    my $ImportFile = $_[0];
    my $Format = $_[1];
    my $OutputFormat = $_[2];

    $Format =~ /(.*)NN(.*)SS(.*)/;
    my $FormatBeg = $1;
    my $FormatMid = $2;
    my $FormatEnd = $3;

    my $Pos = rindex($ImportFile,"/");
    my $InputFileName = $ImportFile;
    if ($Pos>0)
      {
	$InputFileName = substr($ImportFile,$Pos+1);
      }

    if ($InputFileName =~ m/$FormatBeg\d+$FormatMid\d+$FormatEnd/ )
      {
	$InputFileName =~ /$FormatBeg(\d+)$FormatMid(\d+)$FormatEnd/;
	my $NbrParticles = $1;
	my $NbrFlux = $2;
	my $OutputString="";
	if (($OutputFormat =~ m/NN/ )&&($OutputFormat =~ m/SS/ ))
	  {
	    $OutputString = "-o $OutputFormat";
	    $OutputString =~ s#NN#$NbrParticles#;
	    $OutputString =~ s#SS#$NbrFlux#;
	  }
	my $Command = "$Program$Extension $Reverse $InputType $ImportFile -p $NbrParticles -l $NbrFlux $OutputString";
	print ("Import command: ".$Command."\n");
	system($Command);
      }
    else
      {
	print ("Could not match pattern for ".$ImportFile."\n");
      }
  }
