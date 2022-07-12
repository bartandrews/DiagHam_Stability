#!/usr/bin/perl -w

use strict 'vars';

my $Shift = 0.001;
while ($Shift <= )
  {
    print ("proceeding 2S=".$S."\n");
    `/home/regnault/development/DMRG/DiagHam/src/Programs/QHEBosonsDelta $NbrBosons $S`;
    my $LzFileName = "bosons_delta_n_".$NbrBosons."_2s_".$S."_lz.dat";
    my $LFileName = "bosons_delta_n_".$NbrBosons."_2s_".$S."_l.dat";
    `/home/regnault/development/DMRG/SpectrumTools/LzToL $LzFileName > $LFileName`;
    $Shift += 0.001;
  }

