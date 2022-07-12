#!/usr/bin/perl -w

use strict 'vars';

if (!(defined($ARGV[2])))
  {
    die "usage: BosonsTorusBatch nbr_bosons min_2s max_2s\n";
  }

my $NbrBosons = $ARGV[0];
my $MaxS = $ARGV[2];
my $S = $ARGV[1];
my $NbrEigenvalue = "";
if (defined($ARGV[3]))
  { 
    $NbrEigenvalue = "-n ".$ARGV[3]." ";
  }
if (!(-e "n_$NbrBosons"))
  {
    mkdir("n_$NbrBosons");
  }
chdir("n_$NbrBosons");

while ($S <= $MaxS)
  {
    my $Ratio = 2.0 * sqrt(3) / $S;
    print ("proceeding N=".$NbrBosons." 2S=".$S."\n");
    print ("/home/regnault/development/DMRG/DiagHam/src/Programs/QHEBosonsTorus -S --processors 2 -p $NbrBosons -l $S -r $Ratio $NbrEigenvalue\n");
    `/home/regnault/development/DMRG/DiagHam/src/Programs/QHEBosonsTorus -S --processors 2 -p $NbrBosons -l $S -r $Ratio $NbrEigenvalue`;
    $S += 1;
  }
`/home/regnault/development/DMRG/DiagHam/scripts/TorusGraph.pl $NbrBosons`;
chdir ("..");
