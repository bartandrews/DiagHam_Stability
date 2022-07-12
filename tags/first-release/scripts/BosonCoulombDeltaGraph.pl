#!/usr/bin/perl -w

use strict 'vars';

if (!(defined($ARGV[0])))
  {
    die "usage: BosonsDeltaGraph nbr_bosons print_flag\n";
  }
my $NbrBosons = $ARGV[0];
my $PrintFlag = 0;
if (defined($ARGV[1]))
  {
    $PrintFlag = 1;
  }
my @ListLzFiles;
my $TmpFile;
foreach $TmpFile (<*>)
  {
    if (($TmpFile =~ /bosons\_delta\_.*\_n\_$NbrBosons.*\_lz\./ ) || ($TmpFile =~ /bosons\_coulomb\_.*\_n\_$NbrBosons.*\_lz\./ ))
      {
        push (@ListLzFiles, $TmpFile);
      }
  }
my @ListFiles;
foreach $TmpFile (@ListLzFiles)
  {
    my $TmpFile2 = $TmpFile;
    $TmpFile2 =~ s/\_lz/\_l/;
    print "$TmpFile\n";
    `/home/regnault/development/DMRG/SpectrumTools/LzToL $TmpFile > $TmpFile2`;
    push (@ListFiles, $TmpFile2);
  }

foreach $TmpFile (@ListFiles)
  {
    my $Max;
    my $Min;
    &CreatePostScript($TmpFile, $PrintFlag);
    print ("\n\n");
#    &FindMinMax($TmpFile, 1, \$Min, \$Max, 0, 0, 15);
#    print ($TmpFile." ".$Min." ".$Max."\n");
  }


# find minimum and maximum values in a file
#
# $_[0] = file name
# $_[1] = column where to search
# $_[2] = reference on min value
# $_[3] = reference on max value
# $_[4] = additonnal constraint on another column
# $_[5] = max value for the constraint
# $_[6] = min value for the constraint

sub FindMinMax
  {
    my $FileName = $_[0];
    my $Column = $_[1];
    my $Min = $_[2];
    my $Max = $_[3];
    if (defined($_[4]))
      {
        my $ColumnConstraint = $_[4];
        my $MinConstraint = $_[5];
        my $MaxConstraint = $_[6];
        my $Flag = 0;
        open (INFILE, $FileName);
        my $TmpLine;
        my @TmpArray;
#       $TmpLine = <INFILE>;
#       $TmpLine = <INFILE>;
        foreach $TmpLine (<INFILE>)
          {
            chomp ($TmpLine);
            @TmpArray = split (/ /, $TmpLine);
            if (($TmpArray[$ColumnConstraint] <= $MaxConstraint) && ($TmpArray[$ColumnConstraint] >= $MinConstraint))
              {
                if ($Flag == 0)
                  {
                    $$Min = $TmpArray[$Column];
                    $$Max = $$Min;
                    $Flag = 1;
                  }
                else
                  {
                    if ($TmpArray[$Column] < $$Min)
                      {
                        $$Min = $TmpArray[$Column];
                      }
                    if ($TmpArray[$Column] > $$Max)
                      {
                        $$Max = $TmpArray[$Column];
                      }
                  }
              }
          }
      }
    else
      {
        open (INFILE, $FileName);
        my $TmpLine;
        my @TmpArray;
#       $TmpLine = <INFILE>;
#       $TmpLine = <INFILE>;
        $TmpLine = <INFILE>;
        @TmpArray = split (/ /, $TmpLine);
        $$Min = $TmpArray[$Column];
        $$Max = $$Min;
        foreach $TmpLine (<INFILE>)
          {
            chomp ($TmpLine);
            @TmpArray = split (/ /, $TmpLine);
            if ($TmpArray[$Column] < $$Min)
              {
                $$Min = $TmpArray[$Column];
              }
            if ($TmpArray[$Column] > $$Max)
              {
                $$Max = $TmpArray[$Column];
              }
          }
        close (INFILE);
      }
}

# create postscript graph from data file
#
# $_[0] = file name
# $_[1] = print flag (1 if true)

sub CreatePostScript
  {
    my $FileName = $_[0];
    my $PrintFlag = $_[1];
    my $Max;
    my $Min;
    &FindMinMax($FileName, 1, \$Min, \$Max, 0, 0, 14);
    my $Delta = ($Max - $Min) / 20.0;
    $Max += $Delta;
    $Min -= $Delta;
    my $TmpFileName = "tmp".time().".p";
    my $OutputFile = $FileName;
    my @TmpArray = split (/_/,  $OutputFile);
    my $Title = "N = ".$TmpArray[3]."  2S = ".$TmpArray[5];
    $OutputFile =~ s/\_l\.dat/\.ps/;
    open (OUTFILE, ">$TmpFileName");
    print OUTFILE ("set xrange [-1:15]
set yrange [".$Min.":".$Max."]
set xlabel \"Angular Momentum l\"
set ylabel \"E(L)\"
set size 1.0, 0.6
set terminal postscript portrait enhanced \"Helvetica\" 14
set output \"".$OutputFile."\"
plot \"".$FileName."\" using 1:2 title \"".$Title."\"
");
    close (OUTFILE);
    open (INFILE, $TmpFileName);
    my $TmpLine;
    foreach $TmpLine (<INFILE>)
      {
        print $TmpLine;
      }
    close (INFILE);

    `gnuplot $TmpFileName`;
    if ($PrintFlag == 1)
      {
        `lpr $OutputFile`;
      }
    `rm -f $TmpFileName`;
  }

