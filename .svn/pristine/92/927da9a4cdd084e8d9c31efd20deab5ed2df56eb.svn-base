#!/usr/bin/perl -w

my $Column = $ARGV[0];
my @ReplacementString;
my $DisplayHelp=0;
my $Divisor=1.0;
my $Multiplier=1.0;


while( (defined($ARGV[0])&&$ARGV[0] =~ /^-/ ))
  {
    if ( $ARGV[0] =~ /-c/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $Column = $ARGV[0];
	    $Column =~ s/-c//;
	  }
	else
	  {
	    shift(@ARGV);
	    $Column = $ARGV[0];
	  }
      }
    if ( $ARGV[0] =~ /-d/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $Divisor = $ARGV[0];
	    $Divisor =~ s/-d//;
	  }
	else
	  {
	    shift(@ARGV);
	    $Divisor = $ARGV[0];
	  }
      }
    if ( $ARGV[0] =~ /-h/ )
      {
	$DisplayHelp=1;
      }
    if ( $ARGV[0] =~ /-i/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $TmpString = $ARGV[0];
	    $TmpString =~ s/-i//;
	  }
	else
	  {
	    shift(@ARGV);
	    $TmpString = $ARGV[0];
	  }
	push(@ReplacementString,$TmpString);
      }
    if ( $ARGV[0] =~ /-m/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $Multiplier = $ARGV[0];
	    $Multiplier =~ s/-m//;
	  }
	else
	  {
	    shift(@ARGV);
	    $Multiplier = $ARGV[0];
	  }
      }
    shift(@ARGV);
  }

if ((!defined($ARGV[0]))||($DisplayHelp==1))
  {
    die "usage: FindMinMax.pl -c Column -i Replacements [-i ...] [-d Divisor] [-m Multiplier] -h SourceFiles";
  }


my @Replacements=();
my @NbrReplacements=();
&ParseExpressions(\@ReplacementString,\@Replacements,\@NbrReplacements);




foreach $FileName (@ARGV)
  {
    open (INFILE, $FileName);
    my $TmpLine;
    my $MinValue=1e99;
    my $MaxValue=-1e99;
    my $Average;
    my @TmpBlock;
    my $Sum=0.0;
    my $NbrLine=0;
    foreach $TmpLine (<INFILE>)
      {
	chomp ($TmpLine);
	if (!( $TmpLine =~ m/^\#/ ))
	  {
	    my @TmpArray = split (/\s/, $TmpLine);
	    if ($#TmpArray < $Column )
	      {
		die ("File ".$FileName." has less than the requested ".($Column+1)." columns");
	      }
	    if ( $TmpArray[$Column] < $MinValue )
	      {
		$MinValue=$TmpArray[$Column];
	      }
	    if ( $TmpArray[$Column] > $MaxValue )
	      {
		$MaxValue=$TmpArray[$Column];
	      }
	    $Sum+=$TmpArray[$Column];
	    ++$NbrLine;
	  }
      }
    $Average=$Sum/$NbrLine;

    my @RepValues;
    my $NbrRep = $#Replacements + 1;

    for (my $i=0; $i<$NbrRep; ++$i)
      {
	my $Value = $FileName;
	my $Format = $Replacements[$i];
	if ($Value =~ m/$Format/)
	  {
	    $Value =~ /$Format/;
	    for (my $j=1; $j<=$NbrReplacements[$i]; ++$j)
	      {
		push(@RepValues, $$j);
	      }
	  }
	else
	  {
	    for (my $j=1; $j<=$NbrReplacements[$i]; ++$j)
	      {
		push(@RepValues, 0);
	      }
	  }
	#print ("Replacement $i : ".$Format." value: ".$Value."\n");
      }
    my $item;
    foreach $item (@RepValues)
      {
	print ("$item\t");
      }

    print ($Average*$Multiplier/$Divisor."\t".($MaxValue-$Average)*$Multiplier/$Divisor."\t".($Average-$MinValue)*$Multiplier/$Divisor."\t".$MaxValue*$Multiplier/$Divisor."\t".$MinValue*$Multiplier/$Divisor."\t".$FileName."\n");
  }

# analyse filename and add select parts of filename as columns
#
# $_[0] = reference to an array of string value of expressions in printf conventions (input)
# $_[1] = reference to an array of string value of expressions in perl conventions (output)
# $_[2] = reference to an array where the number of replacements per string shall be returned

sub ParseExpressions
  {
    my $ReplacementString = $_[0];
    my $Replacements = $_[1];
    my $NbrReplacements = $_[2];
    my $NbrFields = $#$ReplacementString + 1;

    for (my $i=0; $i<$NbrFields; ++$i)
      {
	my $Format = @$ReplacementString[$i];
	my $Count = 0;
	while ( $Format =~ m/%/ )
	  {
	    my $Recognized=0;
	    if ( $Format =~ m/%f/ )
	      {
		# print ("Recognized floating point replacement ".@$ReplacementString[$i]." ");
		$Format =~ s/%f/\(-*\\d*\[\\.\]*\\d*e*-*\\d*\)/;
		# print ("perl format: ".$Format."\n");
		$Recognized=1;
		++$Count;
	      }
	    if (($Recognized==0) && ($Format =~ m/%d/ ))
	      {
		# print ("Recognized integer replacement ".@$ReplacementString[$i]." ");
		$Format =~ s/%d/\(\\d*\)/;
		# print ("perl format: ".$Format."\n");
		$Recognized=1;
		++$Count;
	      }
	    if (($Recognized==0) && ($Format =~ m/%s/ ))
	      {
		# print ("Recognized string replacement ".@$ReplacementString[$i]." ");
		$Format =~ s/%s/\(\.*\)/;
		# print ("perl format: ".$Format."\n");
		$Recognized=1;
		++$Count;
	      }
	    if ($Recognized==0)
	      {
		print("Attention: format ".@$ReplacementString[$i]." not recognized\n");
		exit(1);
	      }
	  }
	push(@$Replacements,$Format);
	push(@$NbrReplacements,$Count);
      }
  }
