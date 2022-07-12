#!/usr/bin/perl -w

use strict 'vars';

if (!(defined($ARGV[1])))
  {
    die "usage: BosonsBench min_n max_n\n";
  }

my $NbrBosons = $ARGV[0];
my $MaxNbrBosons = $ARGV[1];
my $NbrIter = 2;

my %Benchmarks;
my %RealTimeBenchmarks;
my %SMPBenchmarks;
my %SMPRealTimeBenchmarks;
while ($NbrBosons <= $MaxNbrBosons)
  {
    my $S = 2 * ($NbrBosons - 1);
    $Benchmarks{$NbrBosons} = 0;
    my $TmpOutputFile = "tmp".time().".out";
    my $Iter = 0;
    while ($Iter < $NbrIter)
      {
	`time -f \"\%U \%e\" -o time_smp.out /home/regnault/development/DMRG/DiagHam/src/Programs/QHEBosonsCoulombSMP $NbrBosons $S > $TmpOutputFile`;
	`time -f \"\%U \%e\" -o time.out /home/regnault/development/DMRG/DiagHam/src/Programs/QHEBosonsCoulombMono $NbrBosons $S > $TmpOutputFile`;
	my $LzFileName = "bosons_coulomb_n_".$NbrBosons."_2s_".$S."_lz.dat";
	`rm -f $TmpOutputFile`;
	`rm -f $LzFileName`;
	open (INFILE , "time.out");
	my $TmpLine = <INFILE>;
	chomp ($TmpLine);
	my @TimeValues = split (/ /, $TmpLine);
	print ($NbrBosons." : u = ".$TimeValues[0]."s t = ".$TimeValues[1]."s\n");
	$Benchmarks{$NbrBosons} += $TimeValues[0];
	$RealTimeBenchmarks{$NbrBosons} += $TimeValues[1];
	close (INFILE);
	`rm -f time.out`;

	open (INFILE , "time_smp.out");
	my $TmpLine = <INFILE>;
	chomp ($TmpLine);
	my @TimeValues = split (/ /, $TmpLine);
	print ($NbrBosons." : u = ".$TimeValues[0]."s t = ".$TimeValues[1]."s\n");
	$SMPBenchmarks{$NbrBosons} += $TimeValues[0];
	$SMPRealTimeBenchmarks{$NbrBosons} += $TimeValues[1];
	close (INFILE);
	`rm -f time_smp.out`;
	$Iter++;
      }
    $Benchmarks{$NbrBosons} /= $NbrIter;
    $RealTimeBenchmarks{$NbrBosons} /= $NbrIter;
    $SMPBenchmarks{$NbrBosons} /= $NbrIter;
    $SMPRealTimeBenchmarks{$NbrBosons} /= $NbrIter;
    print ($NbrBosons." : mono mean values : u = ".$Benchmarks{$NbrBosons}."s t = ".$RealTimeBenchmarks{$NbrBosons}."s
     SMP mean values : u = ".$SMPBenchmarks{$NbrBosons}."s t = ".$SMPRealTimeBenchmarks{$NbrBosons}."s
---------------------------------------------------------------\n");
    $NbrBosons += 1;
  }
