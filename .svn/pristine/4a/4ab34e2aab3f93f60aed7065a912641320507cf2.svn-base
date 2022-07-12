#!/usr/bin/perl -w
#
# script which extracts parts of a filename to add as additional columns to the data
#
use strict 'vars';
use File::stat;


my @InitialReplacementString;
my @FinalReplacementString;
my @LineNumbers = ();
my $DisplayHelp=0;
my $TmpString;
my $Separator="\t";
my $ExcludeComments=0;
my $RawInput=0;
my $DefaultValue="--";

while( (defined($ARGV[0])&&$ARGV[0] =~ /^-/ ))
  {
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
	push(@InitialReplacementString,$TmpString);
      }
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
	push(@FinalReplacementString,$TmpString);
      }
    if ( $ARGV[0] =~ /-l/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $TmpString = $ARGV[0];
	    $TmpString =~ s/-l//;
	  }
	else
	  {
	    shift(@ARGV);
	    $TmpString = $ARGV[0];
	  }
	@LineNumbers = split(/,/,$TmpString);
	@LineNumbers = sort(@LineNumbers);
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
    if ( $ARGV[0] =~ /-c/ )
      {
	$ExcludeComments=1;
      }
    if ( $ARGV[0] =~ /-h/ )
      {
	$DisplayHelp=1;
      }
    if ( $ARGV[0] =~ /-r/ )
      {
	$RawInput=1;
      }
    shift(@ARGV);
  }

if ((!defined($ARGV[0]))||($DisplayHelp==1))
  {
    print("usage FileNameToColumn.pl [-i InitialColumn] [-f FinalColumn] Filenames\n");
    print("option -i: format descriptor of columns to add before existing data\n");
    print("       -f: format descriptor of columns to add after existing data\n");
    print("       -d: delimiter to use in output instead of TAB\n");
    print("       -e: value to use for empty replacements (default: --)\n");
    print("       -c: exclude comment lines entirely\n");
    print("       -r: print existing lines in raw format without analysing delimiters\n");
    print("       -l: only print a given list (separated by commas) of (uncommented) line numbers (starting at 0) rather than the entire file\n");
    print("       -h: display this help\n");
    exit(1);
  }

my @InitialReplacements;
my @FinalReplacements;
my @NbrInitialReplacements;
my @NbrFinalReplacements;

&ParseExpressions(\@InitialReplacementString,\@InitialReplacements,\@NbrInitialReplacements);
&ParseExpressions(\@FinalReplacementString,\@FinalReplacements,\@NbrFinalReplacements);


my $TmpFile;
foreach $TmpFile (@ARGV)
  {
    if ( -e $TmpFile)
      {
	&ProcessFile($TmpFile, \@InitialReplacements, \@NbrInitialReplacements, \@FinalReplacements, \@NbrFinalReplacements);
      }
    else
      {
	print ("Did not find file ".$TmpFile."\n");
      }
  }

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
# $_[1] = replacement for initial columns
# $_[2] = nbr of replacement for initial columns
# $_[3] = replacement for final columns
# $_[4] = nbr of replacement for final columns

sub ProcessFile
  {
    my $FileName = $_[0];
    my $InitialReplacements = $_[1];
    my $NbrInitialReplacements = $_[2];
    my $FinalReplacements = $_[3];
    my $NbrFinalReplacements = $_[4];
    my $NbrInitial = $#$InitialReplacements + 1;
    my $NbrFinal = $#$FinalReplacements + 1;

    my @InitialValues;
    my @FinalValues;

    for (my $i=0; $i<$NbrInitial; ++$i)
      {
	my $Value = $FileName;
	my $Format = @$InitialReplacements[$i];
	if ($Value =~ m/$Format/)
	  {
	    $Value =~ /$Format/;
	    for (my $j=1; $j<=@$NbrInitialReplacements[$i]; ++$j)
	      {
		push(@InitialValues, $$j);
	      }
	  }
	else
	  {
	    for (my $j=1; $j<=@$NbrInitialReplacements[$i]; ++$j)
	      {
		push(@InitialValues, $DefaultValue);
	      }
	  }
	# print ("Initial Replacement $i : ".$Format." value: ".$Value."\n");
      }

    for (my $i=0; $i<$NbrFinal; ++$i)
      {
	my $Value = $FileName;
	my $Format = @$FinalReplacements[$i];
	if ($Value =~ m/$Format/)
	  {
	    $Value =~ /$Format/;
	    for (my $j=1; $j <= @$NbrFinalReplacements[$i]; ++$j)
	      {
		push(@FinalValues, $$j);
	      }
	  }
	else
	  {
	    push(@FinalValues, $DefaultValue);
	  }
	# print ("Final Replacement $i : ".$Format." value: ".$Value."\n");
      }

    open (INPUTFILE, "$FileName");

    my $TmpLine;
    my $NbrLine=0;
    my $NextLine=0;
    my $Continue=1;
    my $PrintLine;
    my $item;
    while ((defined($TmpLine = <INPUTFILE> )) && ($Continue==1))
      {
	chomp ($TmpLine);
	if ( $TmpLine =~ m/\s*#/ )
	  {
	    if ( $ExcludeComments == 0 )
	      {
		print ("#".$FileName.": ".$TmpLine."\n");
	      }
	  }
	else
	  {
	    if (defined $LineNumbers[$NextLine])
	      {
		if ($NbrLine==$LineNumbers[$NextLine])
		  {
		    $PrintLine=1;
		    ++$NextLine;
		    if (!defined $LineNumbers[$NextLine])
		      {
			$Continue=0;
		      }
		  }
		else
		  {
		    $PrintLine=0;
		  }
	      }
	    else
	      {
		$PrintLine=1;
	      }
	    if ($PrintLine==1)
	      {
		foreach $item (@InitialValues)
		  {
		    print ("$item$Separator");
		  }
		if ($RawInput == 1)
		  {
		    print ($TmpLine);
		  }
		else
		  {
		    my @LineItems = split(/\s/,$TmpLine);
		    my $ItemNbr = 0;
		    foreach $item (@LineItems)
		      {
			if ($ItemNbr==0)
			  {
			    print ($item);
			  }
			else
			  {
			    print ("$Separator$item");
			  }
			++$ItemNbr;
		      }
		  }
		foreach $item (@FinalValues)
		  {
		    print ("$Separator$item");
		  }
		print ("\n");
	      }
	    ++$NbrLine;
	  }
      }
  }
