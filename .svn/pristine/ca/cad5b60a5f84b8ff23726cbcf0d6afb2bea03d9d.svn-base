#!/usr/bin/perl -w
#
# script for analysis of eigenvector files designed to run on TCM group PC's
#
use strict 'vars';


if (!defined($ARGV[0]))
  {
    print("usage SummarizeBosonsOnLattice.pl [Directory|single file]\n");
    exit(1);
  }


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
my $TmpFile;
if ($ReadMore == 1) # array still empty?
  {
    foreach $TmpFile (<*>)
      {
	if ($TmpFile =~ /bosons\_lattice.*q_\d*\.eval/)
	  {
	    push (@ListFiles, $TmpFile);
	  }	
      }
  }

foreach $TmpFile (@ListFiles)
  {
    print ("Analyzing series of ".$TmpFile."\n");
    &AnalyzeProtocols($TmpFile);
  }


# summarize ground-state properties at different values of the flux
#
# $_[0] = spectrum file name

sub AnalyzeProtocols
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
    my $BaseName = $FileName;
    if ($FileName =~ /bosons\_lattice.*\_q\_(\d*).eval/)
      {
	$q = $1;
	$BaseName =~ s/q\_$q.eval/q/;	
      }
    else
      {
	$BaseName =~ s/q.eval/q/;	
      }
    if ( $HardCore == 1)
      {
	print "n = ".$N."  x = ".$x."  y = ".$y."  (hardcore bosons)\n";
      }
    else
      {
	print "n = ".$N."  x = ".$x."  y = ".$y."  u = ".$u."\n";
      }
    my $Ns=$x*$y;
    my $TmpLine;
    my $ProtocolName = $BaseName."\.gs";	
    open (OUTFILE, ">$ProtocolName");
    my $CountEntries=0;
    print OUTFILE ("# q\tN\tN_s\tGap\trho_0\trho_0/rho_1\trho_ave\tEVCount\t|K|\tDeg\n");
    print ("q=");
    for ($q=0; $q<=$Ns;++$q)
      {
	$FileName = $BaseName."_${q}.eval";
	print ("${q}..");
	if ( -e $FileName )
	  {
	    $TmpLine = `grep ^0 $FileName`;
	    if ($TmpLine gt "")
	      {
		++$CountEntries;
		chomp($TmpLine);
		my @Data = split(/\t/,$TmpLine);
		chomp(@Data);
		my $Energy0 = $Data[1];
		my $Rho0 = $Data[2];
		my $Rho01 = $Data[3];
		my $RhoBar = $Data[4];
		my $EVCount = $Data[5];
		my $AbsK = $Data[6];
		my $Degeneracy = $Data[7];
		my $Gap = -1;
		$TmpLine = `grep ^1 $FileName`;
		if ($TmpLine gt "")
		  {
		    chomp($TmpLine);
		    @Data = split(/\t/,$TmpLine);
		    chomp(@Data);
		    my $Energy1 = $Data[1];
		    $Gap=$Energy1-$Energy0;
		  }
		print OUTFILE ("${q}\t${N}\t${Ns}\t${Gap}\t${Rho0}\t${Rho01}\t${RhoBar}\t${EVCount}\t${AbsK}\t${Degeneracy}\n");
	      }
	  }
      }
    print ("done!\n");
    close(OUTFILE);
    if ($CountEntries == 0)
      {
	print ("No data found for present series\n");
	system ("rm $ProtocolName");
      }
  }


