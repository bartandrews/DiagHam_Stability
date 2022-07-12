#!/usr/bin/perl -w

use strict 'vars';

unless (open (INFILE, "configure"))
  {
    if (-e "configure")
      {
	die "can't open configure\n"
      }
    else
      {
	die "can't find configure, please run AIXPatch in the top directory\n"
      }
  }

my $TmpLine;
my $TmpFile = "";
foreach $TmpLine (<INFILE>)
  {
    if ($TmpLine =~ /CXXFLAGS\=\"\-O2 \-Wall \"/)
      {
	$TmpLine =~ s/CXXFLAGS\=\"\-O2 \-Wall \"/CXXFLAGS\=\"\-O2 \-bbigtoc\"/;
      }
    $TmpLine =~ s/mpiCC/\"mpCC\_r \-cpp\"/g;
    $TmpFile .= $TmpLine;
  }
close(INFILE);
unless (open (OUTFILE, ">configure"))
  {
    die "can't overwrite configure\n"
  }
print OUTFILE $TmpFile;
close (OUTFILE);

&PatchMakefileIn();

# patch the Makefile.in in the current directory and all subdirectories
#

sub PatchMakefileIn()
  {
    my $TmpFileName;
    foreach $TmpFileName (<*>)
      {
	if (-d $TmpFileName)
	  {
	    chdir ($TmpFileName);
	    &PatchMakefileIn();
	    chdir "..";
	  }
      }    
     if (-e "Makefile.in")
      {
	unless (open (INFILE, "Makefile.in"))
	  {
	    my $CurrentDirectory = `pwd`;
	    chomp ($CurrentDirectory);
	    die "can't open Makefile.in in $CurrentDirectory\n"
	  }
	my $TmpLine;
	my $TmpFile = "";
	foreach $TmpLine (<INFILE>)
	  {
	    if ($TmpLine =~ /\$\(COMPILE\) \-Wp\,\-MD\,\.deps\/\$\(\*F\)\.pp \-c \$\</)
	      {
		$TmpLine =~ s/\$\(COMPILE\) \-Wp\,\-MD\,\.deps\/\$\(\*F\)\.pp \-c \$\</\$\(COMPILE\) \-M \-c \$\</;
		$TmpFile .= $TmpLine;
		$TmpFile .= "	mv \$(*F).u .deps/\$(*F).pp\n";
	      }
	    else
	      {
		if ($TmpLine =~ /\$\(CXXCOMPILE\) \-Wp\,\-MD\,\.deps\/\$\(\*F\)\.pp \-c \$\</)
		  {
		    $TmpLine =~ s/\$\(CXXCOMPILE\) \-Wp\,\-MD\,\.deps\/\$\(\*F\)\.pp \-c \$\</\$\(CXXCOMPILE\) \-M \-c \$\</;
		    $TmpFile .= $TmpLine;
		    $TmpFile .= "	mv \$(*F).u .deps/\$(*F).pp\n";
		  }
		else
		  {
		    $TmpFile .= $TmpLine;
		  }
	      }
	  }
	close (INFILE);
	unless (open (OUTFILE, ">Makefile.in"))
	  {
	    my $CurrentDirectory = `pwd`;
	    chomp ($CurrentDirectory);
	    die "can't overwrite Makefile.in in $CurrentDirectory\n"
	  }
	print OUTFILE $TmpFile;
	close (OUTFILE);
      }
  }

