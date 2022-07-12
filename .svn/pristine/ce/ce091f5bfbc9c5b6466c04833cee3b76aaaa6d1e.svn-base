#!/usr/bin/perl -w

open(INFILE, "src/config.h");
my $TmpLine;
my $TmpFile = "";
foreach $TmpLine (<INFILE>)
  {
    if ($TmpLine =~ /^\s*\#define MACHINE\_PRECISION 1e\-14/)
      {
	$TmpFile .= "#define MACHINE_PRECISION 1e-18\n";
      }
    else
      {
	$TmpFile .= $TmpLine;
      }
  }
close (INFILE);
open(OUTFILE, ">src/config.h");
print OUTFILE $TmpFile;
close(OUTFILE);

&ApplySwitchToDirectoryTree();

# switch each occurrence of double to long double in a file
#
# $_[0] = name of the file

sub SwitchToLongDouble
  {
    open(INFILE, "$_[0]");
    my $TmpLine;
    my $TmpFile = "";
    foreach $TmpLine (<INFILE>)
      {
	if (!($TmpLine =~ /^\s*\/\//))
	  {
	    $TmpLine =~ s/double/long double/mg;
	    $TmpLine =~ s/long long double/long double/mg;
	    $TmpLine =~ s/pow\s*\(/powl\(/mg;
	    $TmpLine =~ s/cos\s*\(/cosl\(/mg;
	    $TmpLine =~ s/sin\s*\(/sinl\(/mg;
	    $TmpLine =~ s/tan\s*\(/tanl\(/mg;
	    $TmpLine =~ s/acos\s*\(/acosl\(/mg;
	    $TmpLine =~ s/asin\s*\(/asinl\(/mg;
	    $TmpLine =~ s/atan\s*\(/atanl\(/mg;
	    $TmpLine =~ s/atan2\s*\(/atan2l\(/mg;
	    $TmpLine =~ s/cos\s*\(/cosl\(/mg;
	    $TmpLine =~ s/exp\s*\(/expl\(/mg;
	    $TmpLine =~ s/exp10\s*\(/exp10l\(/mg;
	    $TmpLine =~ s/exp2\s*\(/exp2l\(/mg;
	    $TmpLine =~ s/log\s*\(/logl\(/mg;
	    $TmpLine =~ s/sqrt\s*\(/sqrtl\(/mg;
	  }
	$TmpFile .= $TmpLine;
      }
    close(INFILE);
    open(OUTFILE, ">$_[0]");
    print OUTFILE $TmpFile;
    close(OUTFILE);
  }

# switch each occurrence of double to long double for each source files of a directory tree
#

sub ApplySwitchToDirectoryTree
  {
    my $TmpFile;
    foreach $TmpFile (<*>)
      {
	if ((-d $TmpFile) && ($TmpFile ne "Output"))
	  {
	    chdir ($TmpFile);
	    &ApplySwitchToDirectoryTree();
	    chdir ("..");
	  }
	else
	  {
	    if (($TmpFile =~ /\.cc$/) || ($TmpFile =~ /\.h$/))
	      {
		&SwitchToLongDouble($TmpFile);
	      }
	  }
      }
  }
