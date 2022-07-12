#!/usr/bin/perl -w

use strict 'vars';
my $Pattern = &CleanPattern($ARGV[0]);
&SearchPattern($ARGV[0]);


# search pattern in .cc files in surrent and sub directories
#
# $_[0] = pattern to find

sub SearchPattern
  {
    my $TmpFile;
    foreach $TmpFile (<*>)
      {
	if (-d $TmpFile)
	  {
	    chdir ($TmpFile);
	    &SearchPattern();
	    chdir("..");
	  }
	if ((-f $TmpFile) && ($TmpFile =~ /\.cc$/))
	  {
	    open (INFILE, $TmpFile);
	    my $TmpLine;
	    my $Line = 1;
	    while (defined ($TmpLine = <INFILE>))
	      {
		chomp ($TmpLine);
		if ($TmpLine =~ /$Pattern/)
		  {
		    print "in ".$TmpFile. " line ".$Line.": ".$TmpLine."\n";
		  }
		$Line++;
	      }
	    close (INFILE);
	    
	  }
      }
  }



# clean pattern to allow use with regular expression
#
# $_[0] = pattern to clean
# return value = cleaned pattern

sub CleanPattern()
  {
    my $Pattern = $_[0];
    $Pattern =~ s/\./\\\./;
    $Pattern =~ s/\-/\\\-/;
    $Pattern =~ s/\>/\\\>/;
    return $Pattern;
  }
