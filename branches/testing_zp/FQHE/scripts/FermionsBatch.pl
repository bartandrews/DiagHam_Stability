#!/usr/bin/perl -w

use strict 'vars';

if (!(defined($ARGV[2])))
  {
    die "usage: FermionsBatch nbr_bosons min_2s max_2s\n";
  }

my $NbrFermions = $ARGV[0];
my $MaxS = $ARGV[2];
my $S = $ARGV[1];

if (!(-e "n_$NbrFermions"))
  {
    mkdir("n_$NbrFermions");
  }
chdir("n_$NbrFermions");

while ($S <= $MaxS)
  {
    print ("proceeding N=".$NbrFermions." 2S=".$S."\n");
    `/home/regnault/development/DMRG/DiagHam/src/Programs/QHE/QHEFermionsV3 -S --processors 2 -p $NbrFermions -l $S`;
    my $LzFileName = "fermions_v3_n_".$NbrFermions."_2s_".$S."_lz.dat";
    my $LFileName = "fermions_v3_n_".$NbrFermions."_2s_".$S."_l.dat";
    `/home/regnault/development/DMRG/SpectrumTools/LzToL $LzFileName > $LFileName`;
    $S += 1;
  }
`/home/regnault/development/DMRG/DiagHam/scripts/FermionsGraphFromL.pl $NbrFermions`;
chdir ("..");
