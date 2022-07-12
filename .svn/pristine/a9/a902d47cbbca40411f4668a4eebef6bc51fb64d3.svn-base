#!/usr/bin/perl -w

use strict 'vars';

my $TmpFile;
foreach $TmpFile (<*>)
  {
    if ($TmpFile =~ /\_l\.dat/ )
      {
	print ("proceed ".$TmpFile."...\n");
	open (INFILE, $TmpFile);
	my $TmpData = "";
	my $TmpLine = <INFILE>;
	if ($TmpLine =~ /boson/)
	  {
	    $TmpLine = <INFILE>;
	  }
	else
	  {
	    $TmpData .= $TmpLine;
	  }
	while (defined ($TmpLine = <INFILE>))
	  {
	    $TmpData .= $TmpLine;
	  }
	close (INFILE);
	open (OUTFILE, ">$TmpFile");
	print OUTFILE $TmpData;
	close (OUTFILE);
      }
  }

