#!/usr/bin/perl -w

use strict 'vars';

if (!(defined($ARGV[2])))
  {
    die "usage: BosonsBatch nbr_bosons min_2s max_2s\n";
  }

my $NbrBosons = $ARGV[0];
my $MaxS = $ARGV[2];

my $S = $ARGV[1];

while ($S <= $MaxS)
  {
    print ("proceeding 2S=".$S."\n");
    `/home/regnault/development/DMRG/DiagHam/src/Programs/QHEBosonsDelta $NbrBosons $S`;
    my $LzFileName = "bosons_delta_n_".$NbrBosons."_2s_".$S."_lz.dat";
    my $LFileName = "bosons_delta_n_".$NbrBosons."_2s_".$S."_l.dat";
    `/home/regnault/development/DMRG/SpectrumTools/LzToL $LzFileName > $LFileName`;
    $S += 1;
  }
`/home/regnault/development/DMRG/DiagHam/BosonsDeltaGraphFromL.pl $NbrBosons`;
`/home/regnault/development/DMRG/DiagHam/PlotGround.pl $NbrBosons`;
`/home/regnault/development/DMRG/DiagHam/FindGround.pl $NbrBosons`;
