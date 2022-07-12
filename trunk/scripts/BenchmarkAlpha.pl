#!/usr/bin/perl -w

use strict 'vars';

if (!(defined($ARGV[0])))
  {
    die "usage: Benchmark log_file\n";
  }
my $LogFile = $ARGV[0];
unless (open (INFILE, $LogFile))
  {
    die ("can't open ".$ARGV[0]."\n");
  }

my $Line1;
my $Line2;
my $Line3;
my $Line4;
my $Line5;
my %BosonSizeMono;
my @BosonNbrPerSizeMono;
my @BosonTimeMono;
my %FermionSizeMono;
my @FermionNbrPerSizeMono;
my @FermionTimeMono;
my %BosonSizeBipro;
my @BosonNbrPerSizeBipro;
my @BosonTimeBipro;
my %FermionSizeBipro;
my @FermionNbrPerSizeBipro;
my @FermionTimeBipro;
while (defined($Line1 = <INFILE>) && defined($Line2 = <INFILE>) && defined($Line3 = <INFILE>) && defined($Line4 = <INFILE>) && defined($Line5 = <INFILE>))
  {
    chomp ($Line1);
    chomp ($Line2);
    chomp ($Line3);
    chomp ($Line4);
    chomp ($Line5);
    my @TmpArray = split (/ /, $Line1);
    if ($TmpArray[0] eq "bosons")
      {
	if ($TmpArray[4] eq "mono")
	  {
	    &ParseResults ($TmpArray[3], $Line3." ".$Line5, \%BosonSizeMono, \@BosonNbrPerSizeMono, \@BosonTimeMono);
	  }
	else
	  {
	    &ParseResults ($TmpArray[3], $Line3." ".$Line5, \%BosonSizeBipro, \@BosonNbrPerSizeBipro, \@BosonTimeBipro);
	  }
      }
    else
      {
	if ($TmpArray[0] eq "fermions")
	  {
	    if ($TmpArray[4] eq "mono")
	      {
		&ParseResults ($TmpArray[3], $Line3." ".$Line5, \%FermionSizeMono, \@FermionNbrPerSizeMono, \@FermionTimeMono);
	      }
	    else
	      {
		&ParseResults ($TmpArray[3], $Line3." ".$Line5, \%FermionSizeBipro, \@FermionNbrPerSizeBipro, \@FermionTimeBipro);
	      }
	  }
      }
  }

my $TmpName = $LogFile;
$TmpName =~ s/\.log//;
print ("bosons short mono\n");
print &PrintResults(\%BosonSizeMono, \@BosonNbrPerSizeMono, \@BosonTimeMono);
open (OUTFILE, ">".$TmpName."_bosons_mono.dat");
print OUTFILE &PrintResults(\%BosonSizeMono, \@BosonNbrPerSizeMono, \@BosonTimeMono);
close(OUTFILE);
print ("bosons short bipro\n");
print &PrintResults(\%BosonSizeBipro, \@BosonNbrPerSizeBipro, \@BosonTimeBipro);
open (OUTFILE, ">".$TmpName."_bosons_bipro.dat");
print OUTFILE &PrintResults(\%BosonSizeBipro, \@BosonNbrPerSizeBipro, \@BosonTimeBipro);
close(OUTFILE);
print ("fermions short mono\n");
print &PrintResults(\%FermionSizeMono, \@FermionNbrPerSizeMono, \@FermionTimeMono);
open (OUTFILE, ">".$TmpName."_fermions_mono.dat");
print OUTFILE &PrintResults(\%FermionSizeMono, \@FermionNbrPerSizeMono, \@FermionTimeMono);
close(OUTFILE);
print ("fermions short bipro\n");
print &PrintResults(\%FermionSizeBipro, \@FermionNbrPerSizeBipro, \@FermionTimeBipro);
open (OUTFILE, ">".$TmpName."_fermions_bipro.dat");
print OUTFILE &PrintResults(\%FermionSizeBipro, \@FermionNbrPerSizeBipro, \@FermionTimeBipro);
close(OUTFILE);

close(INFILE);

# parse time result line
#
# $_[0] = id corresponding to the run
# $_[1] = time result line (chomped line)
# $_[2] = reference on hash table containing the id associated to the run type (data attach to each key gives position of the datas corresponding id in others arrays)
# $_[3] = reference on array containing the number of runs per id
# $_[4] = reference on array containing total time of run for each id

sub ParseResults
  {
    my $Id = $_[0];
    my $Line = $_[1];
    my $IdHash = $_[2];
    my $NbrPerId = $_[3];
    my $TimePerId = $_[4];
    my $Position = $$IdHash{$Id};
#    print $Line."\n";
    my @TmpArray = split (/  */, $Line);
    $TmpArray[1] =~ s/s//;
    $TmpArray[3] =~ s/s//;
#    print ($TmpArray[1]." ".$TmpArray[3]."\n");
    my @TmpArray2 = split (/m/, $TmpArray[1]);
    my @TmpArray3 = split (/m/, $TmpArray[3]);
    my $RunTime = (($TmpArray2[0] - $TmpArray3[0]) * 60) + $TmpArray2[1] - $TmpArray3[1];
    if (defined ($Position))
      {
	$$NbrPerId[$Position]++;
	$$TimePerId[$Position] += $RunTime;
      }
    else
      {
	$Position = $#{@$TimePerId} + 1;
	%$IdHash = (%$IdHash, $Id, $Position);
	$$NbrPerId[$Position]++;
	$$TimePerId[$Position] += $RunTime;	
      }
  }

# print bechmark results
#
# $_[1] = reference on hash table containing the id associated to the run type (data attach to each key gives position of the datas corresponding id in others arrays)
# $_[2] = reference on array containing the number of runs per id
# $_[3] = reference on array containing total time of run for each id

sub PrintResults
  {
    my $TmpString = "";
    my $IdHash = $_[0];
    my $NbrPerId = $_[1];
    my $TimePerId = $_[2];
    my @SortedArray = sort {$a <=> $b}keys %$IdHash;
    my $Id;
    foreach $Id (@SortedArray)
      {
	my $Position = $$IdHash{$Id};
	$TmpString .= $Id." ".$$NbrPerId[$Position]." ".$$TimePerId[$Position]." ".($$TimePerId[$Position] / $$NbrPerId[$Position])."\n";
      }
    return $TmpString;
  }
