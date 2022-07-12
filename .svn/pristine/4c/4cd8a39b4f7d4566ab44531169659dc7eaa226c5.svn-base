#!/usr/bin/perl -w

use strict 'vars';

my $TmpFileName = GetTemporaryFileName();

my $CVSOutput = `cvs log > $TmpFileName`;

unless (open(INFILE, $TmpFileName))
  {
    die ("error occured while retrieving cvs log\n");
  }

my $TmpLine;
while ((defined ($TmpLine = <INFILE>)) && (index($TmpLine, "RCS file") != 0))
  {
  }

my %ChangeLog;
my $Flag = 0;
while ($Flag == 0)
  {
    while ((defined ($TmpLine = <INFILE>)) && (index($TmpLine, "description:") != 0))
      {	
      }
    if (index($TmpLine, "description:") == 0)
      {
	while ((defined ($TmpLine = <INFILE>)) && (!($TmpLine =~ /^\=+/)))
	  {	
	    while ((defined ($TmpLine = <INFILE>)) && (index($TmpLine, "date: ") != 0))
	      {	
	      }
	    if (index($TmpLine, "date: ") == 0)
	      {
		$TmpLine =~ /^date\: (\d*\/\d*\/\d*) (\d*)\:(\d*)\:(\d*)/;
		my $TmpDate = $1."/".$2."/".$3."/".$4;
		my $Comment = "";
		while ((defined ($TmpLine = <INFILE>)) && (!($TmpLine =~ /^[\-\=]{2,}/)))
		  {
		    chomp ($TmpLine);
		    $TmpLine =~ s/^\s*//;
		    $TmpLine =~ s/\s*$//;
		    if ($TmpLine ne "")
		      {
			$Comment .= $TmpLine."\n";
		      }
		  }
		if ($TmpLine =~ /^[\-\=]+/)
		  {	    
		    if (!(defined($ChangeLog{$TmpDate})))
		      {
			$ChangeLog{$TmpDate} = $Comment;
		      }
		  }
		else
		  {
		    $Flag = 1;
		  }
	      }
	    else
	      {
		$Flag = 1;
	      }
	  }
	if ((!(defined($TmpLine))) || (!($TmpLine =~ /^[\-\=]{2,}/)))
	  {
	    $Flag = 1;	    
	  }
      }
    else
      {
	$Flag = 1;	
      }
  }

close (INFILE);

my $TmpDate;
my %ChangeLog2;

foreach $TmpDate (sort {$b cmp $a} (keys(%ChangeLog)))
  {
    if ((!($ChangeLog{$TmpDate} =~ /empty log message/)) && ($ChangeLog{$TmpDate} ne ""))
      {
	my $TmpDate2 = $TmpDate;
	$TmpDate2 =~ s/^(\d*\/\d*\/\d*).*/$1/;
	$ChangeLog{$TmpDate} =~ s/^\n*//;
	my @TmpArray = split (/\n/, $ChangeLog{$TmpDate});
	my $Author = shift(@TmpArray);
	if (!($Author =~ /\:\s*$/))
	  {
	    @TmpArray = ($Author, @TmpArray);
	    $Author = "Anonymous Coward:"
	  }
	while ((defined($TmpArray[0])) && ($TmpArray[0] =~ /^\s*\n/))
	  {
	    shift(@TmpArray);
	  }
	if (defined($TmpArray[0]))
	  {
	    $TmpArray[0] =~ s/^([^\+\-\*])/\* $1/;
	    $ChangeLog{$TmpDate} = $TmpDate2."  ".$Author;
	    my $TmpLine;
	    foreach $TmpLine (@TmpArray)
	      {		
		$TmpLine =~ s/^[\+\-\*]\s*/\* /;
		if ($TmpLine =~ /^\*/)
		  {
		    $ChangeLog{$TmpDate} .= "\n\n".$TmpLine;
		  }
		else
		  {
		    $ChangeLog{$TmpDate} .= $TmpLine;
		  }
	      }
	    $ChangeLog{$TmpDate} .= "\n\n\n";
	    if (defined($ChangeLog2{$TmpDate2}))
	      {
		if (index ($ChangeLog2{$TmpDate2}, $ChangeLog{$TmpDate}) < 0)
		  {
		    $ChangeLog2{$TmpDate2} .= $ChangeLog{$TmpDate};	    
		  }
	      }
	    else
	      {
		$ChangeLog2{$TmpDate2} = $ChangeLog{$TmpDate};
	      }
	  }
      }
  }

foreach $TmpDate (sort {$b cmp $a} (keys(%ChangeLog2)))
  {
    print $ChangeLog2{$TmpDate};
  }

`rm -f $TmpFileName`;


# get temporary file name
#
# return value = file name

sub GetTemporaryFileName
  {
    my $FileName = "tmp".time();
    my $Pattern1 = "^".$FileName.".*";
    my $Pattern2 = "^".$FileName."(\\d*).*";
    my $TmpFile;
    my $Maximum = -1;
    foreach $TmpFile (<tmp*>)
      {
	if ((-f $TmpFile) && ($TmpFile =~ /\b$Pattern1\b/))
	  {
	    $TmpFile =~ /\b$Pattern2\b/;
	    if (($1 ne "") && ($1 > $Maximum))
	      {
		$Maximum = $1;
	      }
	  }
      }
    if ($Maximum != -1)
      {
	$Maximum++;
	$FileName .= $Maximum;
      }
    else
      {
	$FileName .= "0";
      }
    return $FileName.".log";
  }
