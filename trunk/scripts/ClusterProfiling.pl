#!/usr/bin/perl -w

use strict 'vars';


my $InputFile = $ARGV[0];
my $ClusterDescriptionFile = $ARGV[1];
my $OptimizerOperation = "VectorHamiltonianMultiply core";
if (defined($ARGV[2]))
{
   $OptimizerOperation = $ARGV[2]." core";    
}
my $NbrMPINodes;
my $TotalNbrNodes;
my %ClusterDescription;
my %OperationDetails;

unless (open(INFILE, $InputFile))
  {
    die ("can't open ".$InputFile."\n");
  }

my $TmpLine;

while ((!defined($TotalNbrNodes)) && (defined($TmpLine = <INFILE>)))
  {
    chomp ($TmpLine);
    $TmpLine =~ s/^\s+//;
    $TmpLine =~ s/\s+$//;
    if ($TmpLine =~ /number of nodes\s*\=\s*(\d+)/)
      {
	$TotalNbrNodes = $1;
      }
  }

while ((defined($TmpLine = <INFILE>)) && ($TmpLine !~ /cluster description/))
  {
  }

while ((defined($TmpLine = <INFILE>)) && ($TmpLine !~ /\-\-\-\-+/))
  {
    chomp ($TmpLine);
    $TmpLine =~ s/^\s+//;
    $TmpLine =~ s/\s+$//;
    if (($TmpLine ne "") && ($TmpLine =~ /^node\s+\d+/))
      {
	$TmpLine =~ /^(node\s+\d+)\s*\:\s+hostname\=([^\s]+)\s*cpu\=(\d+)\s+perf\.\s+index\=(\d+\.?\d*)\s+mem\=(\-?\d+)Mb/;
	my %TmpDescription;
	$TmpDescription{"hostname"} = $2;
	$TmpDescription{"nbrcpu"} = $3;
	$TmpDescription{"perf"} = $4;
	$TmpDescription{"mem"} = $5;
	$ClusterDescription{$1} = \%TmpDescription;
      }
  }

while ((defined($TmpLine = <INFILE>)) && ($TmpLine !~ /\-\-\-\-+/))
  {
  }

while (defined($TmpLine = <INFILE>))
  {
    chomp ($TmpLine);
    $TmpLine =~ s/^\s+//;
    $TmpLine =~ s/\s+$//;
    if ($TmpLine ne "")
      {
	$TmpLine =~ /^(node\s+\d+)\s*\:\s*(.*)\s+operation\s+.*\s+(\d+\.\d+)\s+seconds/;	
	my $TmpNodeID = $1;
	my $TmpOperation = $2;
	my $TmpTime = $3;
	if (defined($OperationDetails{$TmpOperation}))
	  {
	    my $TmpNodeInfo = $OperationDetails{$TmpOperation};
	    if (defined($$TmpNodeInfo{$TmpNodeID}))
	      {		
		my $TmpNodeInfo2 = $$TmpNodeInfo{$TmpNodeID};
		$$TmpNodeInfo2{"nbrcall"}++;
		$$TmpNodeInfo2{"totaltime"} += $TmpTime;
		$$TmpNodeInfo2{"totalsqrtime"} += $TmpTime * $TmpTime;
		if ($TmpTime < $$TmpNodeInfo2{"mintime"})
		  {
		    $$TmpNodeInfo2{"mintime"} = $TmpTime;
		  }
		if ($TmpTime > $$TmpNodeInfo2{"maxtime"})
		  {
		    $$TmpNodeInfo2{"maxtime"} = $TmpTime;
		  }
	      }
	    else
	      {
		my %TmpNodeInfo2;
		$TmpNodeInfo2{"nbrcall"} = 1;
		$TmpNodeInfo2{"totaltime"} = $TmpTime;
		$TmpNodeInfo2{"totalsqrtime"} = $TmpTime * $TmpTime;
		$TmpNodeInfo2{"maxtime"} = $TmpTime;
		$TmpNodeInfo2{"mintime"} = $TmpTime;
		$$TmpNodeInfo{$TmpNodeID} = \%TmpNodeInfo2;		
	      }
	  }
	else
	  {
	    my %TmpNodeInfo;
	    my %TmpNodeInfo2;
	    $TmpNodeInfo2{"nbrcall"} = 1;
	    $TmpNodeInfo2{"totaltime"} = $TmpTime;
	    $TmpNodeInfo2{"totalsqrtime"} = $TmpTime * $TmpTime;
	    $TmpNodeInfo2{"mintime"} = $TmpTime;
	    $TmpNodeInfo2{"maxtime"} = $TmpTime;
	    $TmpNodeInfo{$TmpNodeID} = \%TmpNodeInfo2;	       
	    $OperationDetails{$TmpOperation} = \%TmpNodeInfo;
	  }
      }
  }

close(INFILE);

my $OperationName;
my $NodeDependentOperationInfo;
my %OptimizationFactors;

while (($OperationName, $NodeDependentOperationInfo) = each(%OperationDetails))
  {
    print $OperationName." operation :\n";
    my $TmpNodeName;
    my $TmpNodeInfo;
    my $NbrCalls = 0;
    my $TotalNbrCalls = 0;
    my $MaxTotalTime = 0;
    my $TotalCPUTime = 0;
    my $MeanTime = 0;
    while (($TmpNodeName, $TmpNodeInfo) = each(%$NodeDependentOperationInfo))
      {
	if ($$TmpNodeInfo{"nbrcall"} > $NbrCalls)
	  {
	    $NbrCalls = $$TmpNodeInfo{"nbrcall"};
	  }
	if ($$TmpNodeInfo{"totaltime"} > $MaxTotalTime)
	  {
	    $MaxTotalTime = $$TmpNodeInfo{"totaltime"};
	  }
	$TotalCPUTime += $$TmpNodeInfo{"totaltime"};
	$TotalNbrCalls += $$TmpNodeInfo{"nbrcall"};
	$MeanTime += $$TmpNodeInfo{"totaltime"};
      }
    $MeanTime /= $TotalNbrCalls;
    print "  global total time = ".$MaxTotalTime."s\n";
    print "  total CPU time = ".$TotalCPUTime."s\n";
    print "  mean time = ".$MeanTime."s\n";
    print "  number of calls = ".$NbrCalls."\n";
    print "  total number of calls = ".$TotalNbrCalls."\n";
    print "  per node informations :\n";
    foreach $TmpNodeName (sort (keys(%$NodeDependentOperationInfo)))
      {
	$TmpNodeInfo = $$NodeDependentOperationInfo{$TmpNodeName};
	print "    ".$TmpNodeName." : \n";
	print "      total time = ".$$TmpNodeInfo{"totaltime"}."s\n";
	print "      number of calls = ".$$TmpNodeInfo{"nbrcall"}."s\n";
	print "      mean time per call = ".($$TmpNodeInfo{"totaltime"} / $$TmpNodeInfo{"nbrcall"})."+/-".(sqrt(($$TmpNodeInfo{"totalsqrtime"} * $$TmpNodeInfo{"nbrcall"}) - ($$TmpNodeInfo{"totaltime"} * $$TmpNodeInfo{"totaltime"})) / $$TmpNodeInfo{"nbrcall"})."s\n";
	print "      min time per call = ".$$TmpNodeInfo{"mintime"}."s\n";	
	print "      max time per call = ".$$TmpNodeInfo{"maxtime"}."s\n";	
      }
    if ($OperationName eq $OptimizerOperation)
      {
	foreach $TmpNodeName (sort (keys(%$NodeDependentOperationInfo)))
	  {
	    $TmpNodeInfo = $$NodeDependentOperationInfo{$TmpNodeName};
	    $OptimizationFactors{$TmpNodeName} = $MeanTime / ($$TmpNodeInfo{"totaltime"} / $$TmpNodeInfo{"nbrcall"});
	  }	
      }
    print "------------\n";
  }

my $TotalPerFactor = 0;
my $TmpNodeName;
foreach $TmpNodeName (sort (keys(%OptimizationFactors)))
  {
    my $NodeDescription = $ClusterDescription{$TmpNodeName};
    $$NodeDescription{"optperf"} = $$NodeDescription{"perf"} / $$NodeDescription{"nbrcpu"} * $OptimizationFactors{$TmpNodeName};
    $TotalPerFactor += $$NodeDescription{"optperf"};
  }
foreach $TmpNodeName (sort (keys(%OptimizationFactors)))
  {
    my $NodeDescription = $ClusterDescription{$TmpNodeName};
    $$NodeDescription{"optperf"} /= $TotalPerFactor;
  }

print "------------
optimized performance index
------------\n";
foreach $TmpNodeName (sort (keys(%ClusterDescription)))
  {
    my $NodeDescription = $ClusterDescription{$TmpNodeName};
    print "  ".$TmpNodeName." :  perf. index=".$$NodeDescription{"optperf"}."\n";
  }


if (defined($ClusterDescriptionFile))
{
    unless(open (INFILE, $ClusterDescriptionFile))
    {
	die "can't open ".$ClusterDescriptionFile."\n";
    }
    while (defined($TmpLine = <INFILE>))
    {
	my $TmpLine2 = $TmpLine;
	chomp ($TmpLine);
	$TmpLine =~ s/^\s+//;
	$TmpLine =~ s/\s+$//;
	$TmpLine =~ s/^\#.*//;
	if ($TmpLine eq "")
	{
	    print $TmpLine2;	    
	}
	else
	{
	    my @TmpArray = split (/\s+/, $TmpLine);
	    foreach $TmpNodeName (sort (keys(%ClusterDescription)))
	    {	
		my $NodeDescription = $ClusterDescription{$TmpNodeName};
		if ($$NodeDescription{"hostname"} eq $TmpArray[0])
		{
		    $TmpArray[2] = $$NodeDescription{"optperf"};
		}
	    }
	    print join(" ", @TmpArray)."\n";
	}
    }
    close (INFILE);
}
