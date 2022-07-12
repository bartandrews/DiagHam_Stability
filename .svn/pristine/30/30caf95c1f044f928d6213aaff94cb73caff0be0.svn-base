#!/usr/bin/perl -w

if (!(defined($ARGV[0])))
  {
    die "usage: Prepare3DPlot.pl SourceFile";
  }

my $TmpFileName= "tmp".time();

system ("sort -g $ARGV[0] > $TmpFileName");
system ("mv $TmpFileName $ARGV[0]");

open (INFILE, $ARGV[0]);
my $TmpLine;
my $OldValue=-999.94262457514;
my $OutFileName = $ARGV[0]."-3D.dat";
open (OUTFILE, ">$OutFileName");

my @TmpBlock;

foreach $TmpLine (<INFILE>)
  {
    chomp ($TmpLine);
    my @TmpArray = split (/\s/, $TmpLine);
    if ($TmpArray[0] != $OldValue)
      {
	ProcessBlock(\@TmpBlock);
	@TmpBlock=();
	push (@TmpBlock,\@TmpArray);
	$OldValue=$TmpArray[0];
      }
    else
      {
	push (@TmpBlock,\@TmpArray);
      }
  }
ProcessBlock(\@TmpBlock);
print ("to plot, run splot \"$OutFileName\" u 1:2:3 \n");

close(OUTFILE);


# sort block on second row index and print
#
# $_[0] = reference to block that shall be sorted
sub ProcessBlock
  {
    my $Ref = $_[0];
    my @Block = @{$Ref};
    my @SortedBlock;
    if ($#Block >= 0)
      {
	my $ArrayRef;
	my $Item;
#	print "Block0=@{$Block[0]}";
#	print "Block01=${$Block[0]}[1]";
	@Block = sort { $a->[1] <=> $b->[1] } @Block;
#	print "Block0=@{$Block[0]}\n";

	foreach $ArrayRef (@Block)
	  {
	    foreach $Item (@{$ArrayRef})
	      {
		print OUTFILE ($Item."\t");
	      }
	    print OUTFILE ("\n");
	  }
	print OUTFILE ("\n");
      }
  }
