#!/usr/bin/perl -w
#
# script which extracts parts of a filename to add as additional columns to the data
#
use strict 'vars';
use File::stat;


my @ReplacementString;
my $DisplayHelp=0;
my $TmpString;
my $Separator="\t";
my $ExcludeComments=0;
my $DefaultValue="--";

while( (defined($ARGV[0])&&$ARGV[0] =~ /^-/ ))
  {
    if ( $ARGV[0] =~ /-f/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $TmpString = $ARGV[0];
	    $TmpString =~ s/-f//;
	  }
	else
	  {
	    shift(@ARGV);
	    $TmpString = $ARGV[0];
	  }
	push(@ReplacementString,$TmpString);
      }
    if ( $ARGV[0] =~ /-d/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $Separator = $ARGV[0];
	    $Separator =~ s/-d//;
	  }
	else
	  {
	    shift(@ARGV);
	    $Separator = $ARGV[0];
	  }
      }
    if ( $ARGV[0] =~ /-e/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $DefaultValue = $ARGV[0];
	    $DefaultValue =~ s/-e//;
	  }
	else
	  {
	    shift(@ARGV);
	    $DefaultValue = $ARGV[0];
	  }
      }
    if ( $ARGV[0] =~ /-h/ )
      {
	$DisplayHelp=1;
      }
    shift(@ARGV);
  }

if ((!defined($ARGV[0]))||($DisplayHelp==1))
  {
    print("usage FileNameToColumn.pl [-f format] Filename\n");
    print("option -f: format descriptor of piece to extract from filename\n");
    print("       -d: delimiter to use in output instead of TAB\n");
    print("       -e: value to use for empty replacements (default: --)\n");
    print("       -h: display this help\n");
    exit(1);
  }

my @Replacements;
my @NbrReplacements;

&ParseExpressions(\@ReplacementString,\@Replacements,\@NbrReplacements);


my $TmpFile =$ARGV[0];

&ProcessFileName($TmpFile, \@Replacements, \@NbrReplacements);

exit(0);



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



# analyse filename and add select parts of filename as columns
#
# $_[0] = input file name
# $_[1] = replacement for final columns
# $_[2] = nbr of replacement for final columns

sub ProcessFileName
  {
    my $FileName = $_[0];
    my $Replacements = $_[1];
    my $NbrReplacements = $_[2];
    my $NbrItems = $#$Replacements + 1;

    my @InitialValues;
    my @Values;

    for (my $i=0; $i<$NbrItems; ++$i)
      {
	my $Value = $FileName;
	my $Format = @$Replacements[$i];
	if ($Value =~ m/$Format/)
	  {
	    $Value =~ /$Format/;
	    for (my $j=1; $j<=@$NbrReplacements[$i]; ++$j)
	      {
		push(@Values, $$j);
	      }
	  }
	else
	  {
	    push(@Values, $DefaultValue);
	  }
	# print (" Replacement $i : ".$Format." value: ".$Value."\n");
      }

    my $MyItem;
    my @UserInput=<STDIN>;
    my $InputItem;
    foreach $InputItem (@UserInput)
      {
	chomp ($InputItem);
	foreach $MyItem (@Values)
	  {
	    print ("$MyItem$Separator");
	  }
	print ("$InputItem\n");
      }
  }
