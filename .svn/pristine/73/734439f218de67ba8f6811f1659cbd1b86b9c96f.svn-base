#!/usr/bin/perl -w

use strict 'vars';

use Getopt::Long;

my $MainDirectory = ".";
my $PngFlag = 0;
my $PdfFlag = 0;
my $PsFlag = 0;
my $Temperature = -1.0;
my $BinaryDirectory = "/home/regnault/development/Physics/DiagHam/build/src/Tools/QuantumDot/Analysis";

my $Result = GetOptions ("directory:s" => \$MainDirectory, "png" => \$PngFlag, "pdf" => \$PdfFlag, "ps" => \$PsFlag, "bindir:s" => \$BinaryDirectory,
			 "temperature:s" => \$Temperature);

if (!(-d $MainDirectory))
  {
    die ("directorty ".$MainDirectory." does not exist\n");
  }

my $InitialDirectory = $MainDirectory."/initial";
if (!(-d $InitialDirectory))
  {
    die ("directorty ".$InitialDirectory." does not exist\n");
  }
my $FinalDirectory = $MainDirectory."/final";
if (!(-d $FinalDirectory))
  {
    die ("directorty ".$FinalDirectory." does not exist\n");
  }
my %InitialDirectories;
&GetBFieldDirectories ($InitialDirectory, \%InitialDirectories);
my %FinalDirectories;
&GetBFieldDirectories ($FinalDirectory, \%FinalDirectories);

my $TemperatureCaption = "";
my $TemperatureSuffix = "";
if ($Temperature > 0)
  {
    $TemperatureCaption = " T=".$Temperature."K ";
    $TemperatureSuffix = "_t".$Temperature;
  }
my $BField;
my $TmpDirectory1;
my $TmpDirectory2;
my $BroadeningOuptut = "";
my @TmpBFileds = sort {$a <=> $b} (keys(%InitialDirectories));
foreach $BField (@TmpBFileds)
  {    
    if (defined($FinalDirectories{$BField}))
      {
	$TmpDirectory1 = $InitialDirectories{$BField};
	print "parsing B=".$BField."T :\n";
	my $InitialDOS = "initialdos".$BField.".dat";
        my $FinalDOS = "finaldos".$BField.".dat";
	my $Absorption = "absorption".$BField.$TemperatureSuffix.".dat";
	my $Input = sprintf("%.6f", $BField);
	my $Command = $BinaryDirectory."/SumDOS2 --gamma 0.2 --min 0 --max 250 --step 0.1 --directory ".$TmpDirectory1."/run_ --input eigenvalues".$Input.". --output ".$InitialDOS;
	system($Command);
	$Command = $BinaryDirectory."/SumDOS2 --gamma 0.2 --min 0 --max 250 --step 0.1 --directory ".$FinalDirectories{$BField}."/run_ --input eigenvalues".$Input.". --output ".$FinalDOS;
        system($Command);
	$Command = $BinaryDirectory."/AbsorptionInQuantumWellBField --gamma 0.2 --min 0 --max 250 --step 0.1 --initial-dir ".$TmpDirectory1."/run_ --final-dir ".$FinalDirectories{$BField}."/run_ --initial-input eigenvalues".$Input.". -z 125.6 --output ".$Absorption. " --temperature ".$Temperature;
	system($Command);
        $Command = $BinaryDirectory."/EvaluateBroadening --half-height 0.36787944 ".$Absorption;
	$BroadeningOuptut .= "B = ".$BField." T \n";
	$BroadeningOuptut .= `$Command`;
	if ($PsFlag == 1)
          {
            &PlotDOS($InitialDOS, "initial B=".$BField."T", "DOS", 0, 250, "ps");
            &PlotDOS($FinalDOS, "final B=".$BField."T", "DOS", 0, 250, "ps");
	    &PlotDOS($Absorption, "B=".$BField."T".$TemperatureCaption, "absorption (a.u.)", 100, 200, "ps");
          }
	if ($PngFlag == 1)
	  {
	    &PlotDOS($InitialDOS, "initial B=".$BField."T", "DOS", 0, 250, "png");
            &PlotDOS($FinalDOS, "final B=".$BField."T", "DOS", 0, 250, "png");
	    &PlotDOS($Absorption, "B=".$BField."T".$TemperatureCaption, "absorption (a.u.)", 100, 200, "png");
	  }
	if ($PdfFlag == 1)
	  {
            &PlotDOS($InitialDOS, "initial B=".$BField."T", "DOS", 0, 250, "pdf");
            &PlotDOS($FinalDOS, "final B=".$BField."T", "DOS", 0, 250, "pdf");
	    &PlotDOS($Absorption, "B=".$BField."T".$TemperatureCaption, "absorption (a.u.)", 100, 200, "pdf");
          }
      }
  }
unless (open(OUTFILE, ">broadening".$TemperatureSuffix.".log"))
  {
    die ("can't create file broadening".$TemperatureSuffix.".log\n");
  }
print OUTFILE $BroadeningOuptut;
close (OUTFILE);


# get all directories for different b field (directory names have to be of the form bfield_xxxx)
#
# $_[0] = directory where to look at
# $_[1] = reference on theassociative where found directories will be stored (key is b field value and value is with directory name relative path)

sub GetBFieldDirectories()
  {
    my $Directory = $_[0];
    my $Directories = $_[1];
    my $TmpFile;
    my $CurrentDirectory = `pwd`;
    chomp ($CurrentDirectory);
    chdir ($Directory);
    foreach $TmpFile (<*>)
      {
	if ((-d $TmpFile) && ($TmpFile =~ /^bfield\_(\d+\.?\d*)$/))
	  {
	    $$Directories{$1} = $Directory."/".$TmpFile;
	  }
      }
    chdir($CurrentDirectory);
  }


# plot grepah
#
# $_[0] = name of the file containing dos (must have .dat extension)
# $_[1] = caption
# $_[2] = y axis label
# $_[3] = min x value
# $_[4] = max x value 
# $_[5] = output format (png, pdf, fig or ps)

sub PlotDOS()
  {
    my $FileName = $_[0];
    my $Caption = $_[1];
    my $YAxis = $_[2];
    my $MinX = $_[3];
    my $MaxX = $_[4];
    my $Format = $_[5];
    my $TmpOutput = $FileName;
    $TmpOutput =~ s/\.dat$/\.p/;
    my $Output = $FileName;
    $Output =~ s/\.dat$/\.$Format/;
    unless (open (OUTFILE, ">".$TmpOutput))
      {
	dir ("can't create ".$TmpOutput."\n");	  
      }
    print OUTFILE ("set xrange [".$MinX.":".$MaxX."]
set xlabel \"E(meV)\"
set ylabel \"".$YAxis."\"\n");
    if ($Format eq "ps")
      {
	print OUTFILE ("set size 1, 0.9
set terminal postscript landscape enhanced \"Helvetica\" 14\n");
      }
    elsif ($Format eq "png")
      {
	print OUTFILE ("set terminal png medium size 640,480\n");
      }
    elsif ($Format eq "pdf")
      {
        print OUTFILE ("set size 1, 0.9
set terminal pdf enhanced fname \"Helvetica\" fsize 14\n");
      }
  print OUTFILE ("set output \"".$Output."\"
plot \"".$FileName."\" using 1:2 title \"".$Caption."\" with lines 1
");
    close (OUTFILE);
    `gnuplot $TmpOutput`;
    unlink ($TmpOutput);
  }
