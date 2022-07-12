#!/usr/bin/perl -w
#
# script for analysis of eigenvector files designed to run on TCM group PC's
#
use strict 'vars';
use File::stat;

#my $Program_32="/rscratch/gm360/bin/FQHELatticeDensityMatrix";
my $Program_32="FQHELatticeDensityMatrix";
my $Program_64="/rscratch/gm360/bin/FQHELatticeDensityMatrix_64";

if (!defined($ARGV[0]))
  {
    print("usage EvaluateBosonsOnLattice.pl Directory|single file\n");
    exit(1);
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
	if ($TmpFile =~ /bosons\_lattice.*.dat/)
	  {
	    push (@ListFiles, $TmpFile);
	  }
      }
  }

foreach $TmpFile (@ListFiles)
  {
    print ("Analyzing ".$TmpFile."\n");
    &AnalyzeVectors($TmpFile);
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
    if ( $HardCore == 1)
      {
	print "n = ".$N."  x = ".$x."  y = ".$y."  (hardcore bosons)\n";
      }
    else
      {
	print "n = ".$N."  x = ".$x."  y = ".$y."  u = ".$u."\n";
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
	my $ProtocolName = $BaseName."\_$q\.eval";	
	open (OUTFILE, ">$ProtocolName");
	print OUTFILE ("# ID\tEnergy\trho_0\trho_0/rho_1\trho_ave\tEVCount\t|K|\tDeg\tKx\tKy\n");
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
	    $CountEV=0;
	    my $BeginMultiplet=0;
	    my $EndMultiplet=0;
	    my $SizeMultiplet;
	    my $Rank=0;
	    while ($EndMultiplet <= $#EigenVectors)
	      {
		print ("Multiplet ".$Rank."\n");
		while( ($EndMultiplet < $#EigenVectors) &&
		       ( abs($EigenValues[$EndMultiplet+1]-$EigenValues[$EndMultiplet]) < 1e-10 ) )
		  {
		    ++$EndMultiplet;
		  }
		$SizeMultiplet = $EndMultiplet - $BeginMultiplet + 1;
		my $Vectors="";		
		for (my $i=$BeginMultiplet; $i<=$EndMultiplet; ++$i)
		  {
		    $Vectors=$Vectors." ".$EigenVectors[$i];
		  }
		my $ProtocolName2 = $BaseName."\_$q\_r\_$Rank\.dm";
		my $WantToCompute = 1;
		if ( -e $ProtocolName2 )
		  {
		    my $stat1 = stat($ProtocolName2);
		    my $mtime2=0;
		    for (my $i=$BeginMultiplet; $i<=$EndMultiplet; ++$i)
		      {
			my $stat2 = stat($EigenVectors[$i]);
			if ( $stat2->mtime > $mtime2 )
			  {
			    $mtime2 =$stat2->mtime;
			  }
		      }
		    if ($mtime2 > $stat1->mtime)
		      {
			$WantToCompute = 1;
		      }
		    else
		      {
			$WantToCompute = 0;
		      }
		  }
		
		if ( $WantToCompute == 1)
		  {
		    system("$Program $Vectors > $ProtocolName2");
		  }
		else
		  {
		    print ($ProtocolName2." exists, skipping...\n");
		  }
		system("grep ^EV $ProtocolName2 > $TmpFileName");

		open (INFILE2, $TmpFileName);
		my $Threshold = $SizeMultiplet*1.0/($x*$y);
		my @DensityEVs;		
		foreach $TmpLine (<INFILE2>)
		  {
		    chomp($TmpLine);
		    my @Tmp = split (/ /, $TmpLine);
		    push(@DensityEVs,$Tmp[2]);
		  }
		my $DensityRatio=1e6;
		if ($#DensityEVs > 0)
		  {
		    $DensityRatio = $DensityEVs[0]/$DensityEVs[1];
		  }
		my $DensityRho0 = $DensityEVs[0];
		my $AverageDensity=0.0;
		my $LargeEVCount=0;
		while (($LargeEVCount<=$#DensityEVs) && ($DensityEVs[$LargeEVCount] > $Threshold))
		  {
		    $AverageDensity+=$DensityEVs[$LargeEVCount];
		    ++$LargeEVCount;
		  }
		if ( $LargeEVCount > 0.0 )
		  {
		    $AverageDensity /= $LargeEVCount;
		  }
		close(INFILE2);

		system("grep -A255 ^#i $ProtocolName2 | grep -A254 ^0 > $TmpFileName");
		open (INFILE2, $TmpFileName);
		my $Kx=0.0;
		my $Ky=0.0;
		my $NewAbsK;
		my $AbsK=-1.0;
		my $QuestionMark="";
		foreach $TmpLine (<INFILE2>)
		  {
		    chomp($TmpLine);
		    my @Tmp = split (/\t/, $TmpLine);
		    $Kx = $Tmp[1];
		    $Ky = $Tmp[2];
		    $NewAbsK = sqrt($Kx*$Kx + $Ky*$Ky);
		    if (($AbsK>0.0) && ($NewAbsK!=$AbsK))
		      {
			$QuestionMark="\t?";
			print("Unresolved degeneracies or other problems in multiplet $Rank\n");
		      }
		    else
		      {
			if ($AbsK<0.0)
			  {
			    $AbsK = $NewAbsK;
			  }
		      }
		  }
		close(INFILE2);

		if ( $WantToCompute == 1)
		  {
		    system("echo '# Corresponding Energy levels' >> $ProtocolName2");
		    for (my $i=$BeginMultiplet; $i<=$EndMultiplet; ++$i)
		      {
			system("echo 'E${i} = ".$EigenValues[$i]."' >> $ProtocolName2");
		      }
		  }		
		print OUTFILE ($Rank."\t".$EigenValues[$BeginMultiplet]."\t".$DensityRho0."\t".$DensityRatio."\t".$AverageDensity."\t".$LargeEVCount."\t".$AbsK."\t".$SizeMultiplet."\t".$Kx."\t".$Ky.$QuestionMark."\n");
		++$Rank;
		$BeginMultiplet=$EndMultiplet+1;
		$EndMultiplet=$BeginMultiplet;		
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


