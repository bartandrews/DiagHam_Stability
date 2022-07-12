#!/usr/bin/perl -w

use strict 'vars';

# parse arguments 

my @MultipleArguments = ("-x");
my %Arguments;
my $Option = "";
foreach (@ARGV)
  {
    if (/^\-h/)
      {
	&ShowUsage();
      }
    if (($Option eq "") && (/^\-./))
      {
	if (/^\-.$/)
	  {
	    $Option = $_;
	  }
	else
	  {
	    my $Value = $_;
	    $Value =~ s/^\-.//;
	    s/$Value//;
	    %Arguments = (%Arguments, $_, $Value);
	  }
	next;
      }
    if (($Option ne "") && !(/^\-./))
      {
	if (&GetPosition(\@MultipleArguments, $Option) == -1)
	  {
	    %Arguments = (%Arguments, $Option, $_);	
	    $Option = "";
	    next;
	  }
	else
	  {
	    s/\s//;
	    my $TmpValue = $Arguments{$Option};
	    if (!defined($TmpValue))
	      {
		%Arguments = (%Arguments, $Option, $_);	
	      }
	    else
	      {
		$TmpValue .= " ".$_;
		%Arguments = (%Arguments, $Option, $TmpValue);	
	      }
	    next;
	  }
      }
    if (($Option ne "") && (/^\-./))
      {
	if (&GetPosition(\@MultipleArguments, $Option) == -1)
	  {
	    %Arguments = (%Arguments, $Option, "");
	  }
	if (/^\-.$/)
	  {
	    $Option = $_;
	  }
	else
	  {
	    my $Value = $_;
	    $Value =~ s/^\-.//;
	    s/$Value//;
	    %Arguments = (%Arguments, $_, $Value);
	    $Option = "";
	  }
	next;
      }
    die "error : wrong arguments, type genmake -h to show usage\n";
  }

# initialize variables
my $CurrentDir = `pwd`;
chomp ($CurrentDir);
my $LibraryOnlyFlag = 0;
my $Source = "";
my $Executable= "";
my $BaseDirectory = "";
my @ExcludedPackages;

# proceed arguments
my $Arg;
my $Value;
while (($Arg, $Value) = each (%Arguments))
  {
    if ($Arg eq "-l")
      {
	$LibraryOnlyFlag = 1;
	next;
      }
    if ($Arg eq "-e")
      {
	$Executable = $Value;
	next;
      }
    if ($Arg eq "-d")
      {
	$BaseDirectory = $Value;
	next;
      }
    if ($Arg eq "-s")
      {
	$Source = $Value;
	next;
      }
    if ($Arg eq "-x")
      {
	@ExcludedPackages = split (/ /, $Value);
	next;
      }
    die "unknown option $Arg, type genmake -h to show usage\n";
  }

#look if default values are needed
if ($BaseDirectory eq "")
  {
    $BaseDirectory = $CurrentDir;
  }
else
  {
    unless (chdir ($BaseDirectory))
      {
	die "error : can't find $BaseDirectory directory\n";
      }
    $BaseDirectory = `pwd`;
    chomp ($BaseDirectory);
  }
if ($Executable eq "")
  {
    $Executable = $BaseDirectory;
    $Executable =~ s/.*\///;
  }
if ($Source eq "")
  {
    my @Sources = &ListSourceFile();
    if ($#Sources < 0)
      {
	die "error : no source file found\n";
      }
    $Source = $Sources[0];
  }
else
  {
    if (!(-e $Source))
      {
	die "error : file $Source does not exist\n";
      }
  }

#generate makefiles
if ($LibraryOnlyFlag == 1)
  {
    &GenerateMakefile($BaseDirectory, \@ExcludedPackages);
  }
else
  {
    &GenerateBaseMakefile($BaseDirectory, $Source, $Executable, \@ExcludedPackages);
  }

# Show Usage
#

sub ShowUsage()
  {
    print 
"usage : genmake [OPTIONS]

    -e Executable : use Executable as name for executable file (default value = current directory name)
    -s Source : use Source as executable source file (default value = first .cc file in the current directory)
    -d Directory : use Directory as base directory
    -l : indicate base package is a library (do not create executable)
    -x Package1 ... : exclude enumerated packages from list of packages to built
    -h : display help

";
    exit;
  }



# return a list of all directories
#
# return value = list of all sudirectories in current directory

sub ListDir
  {
    foreach (<*>)
      {
	if (-d)
	  {
	    @_=(@_,$_);
	  }
      }    
    return @_;
  }



# return a list of all source files
#
# return value = list of all source files in current directory

sub ListSourceFile
  {
    foreach (<*>)
      {
	if ((-f) && (/\.cc$/))
	  {
	    @_=(@_,$_);
	  }
      }    
    return @_;
  }



# return a list of all header files
#
# return value = list of all headers in current directory

sub ListHeaderFile
  {
    foreach (<*>)
      {
	if ((-f) && (/\.h$/))
	  {
	    @_=(@_,$_);
	  }
      }    
    return @_;
  }



# return a list of all dependence files
#
# return value = list of all dependence files in current directory

sub ListDependenceFile
  {
    foreach (<*>)
      {
	if ((-f) && (/\.d$/))
	  {
	    @_=(@_,$_);
	  }
      }    
    return @_;
  }



# generate base makefile (which produces the executable file)
#
# $_[0] = path to the directory contenaining the config file
# $_[1] = name of the source file needed to create executable file
# $_[2] = executable file name
# $_[3] = reference to the list of excluded packages

sub GenerateBaseMakefile
  {
    # prepare data
    my $CurrentDir = `pwd`;
    chomp ($CurrentDir);
    my $RootDir = $_[0];
    my $Source = $_[1];
    my $ExecutableName = $_[2];
    my $ExcludedPackages = $_[3];

    # erase previous file
    if (-e "Makefile.am")
      {
	`rm -f Makefile.am`;
      }

    my $TmpMakefileAm = "SUBDIRS=";

    # look for subdirectories and generate makefile in them
    my @SubDir = ListDir();
    my $FlagSubDir = 0;
    my $CurrentSubDir;
    my $FirstFlag = 0;
    my @UsedPackages;
    foreach $CurrentSubDir (@SubDir)
      {
 	chdir $CurrentSubDir;
	my $TmpSubDir = &GenerateMakefile($RootDir, $ExcludedPackages, \@UsedPackages);
	if (-e "Makefile.am")
	  {
	    $FlagSubDir = 1;
	  }
	chdir "..";
	if ($TmpSubDir ne "")
	  {
	    if ($FirstFlag == 0)
	      {
		$FirstFlag = 1;
	      }
	    else
	      {
		$TmpMakefileAm .= " ";
	      }
	    $TmpMakefileAm .= $TmpSubDir;
	  }
      }
    
    # fill Makefile.am
    $TmpMakefileAm .= "\n\nbin_PROGRAMS=".$_[2]."\n\n".$_[2]."_SOURCES=";
    my $TmpSrcFile;
    $FirstFlag = 0;
    foreach $TmpSrcFile (&ListSourceFile())
      {
	if ($FirstFlag == 0)
	  {
	    $FirstFlag = 1;
	  }
	else
	  {
	    $TmpMakefileAm .= " \\\n";
	  }
	$TmpMakefileAm .= $TmpSrcFile;
      }
    $TmpMakefileAm .= "\n\n";
    
    $TmpMakefileAm .= "LDADD=";
    foreach $CurrentSubDir (@UsedPackages)
      {
	$TmpMakefileAm .= "-L\@srcdir\@/".$CurrentSubDir." ";
      }
    foreach $CurrentSubDir (@UsedPackages)
      {
	$CurrentSubDir =~ s/^.*\///g;
	$TmpMakefileAm .= "-l".$CurrentSubDir." ";
      }
    $TmpMakefileAm .= "\n";

    unless (open(OUTFILE, ">Makefile.am"))
      {
	die "error : can't create Makefile.am\n";
      }
    print OUTFILE $TmpMakefileAm;
    close (OUTFILE);	
  }



# generate makefile for a given package
#
# $_[0] = path to the directory contenaining the config file
# $_[1] = reference to the list of excluded packages
# $_[2] = reference to the list of used packages
# return value = string associated to the package (empty if package is excluded)

sub GenerateMakefile
  {
# prepare data
    my $CurrentDir = `pwd`;
    chomp ($CurrentDir);
    my $RootDir = $_[0];
    my $ExcludedPackages = $_[1];
    my $UsedPackages = $_[2];
# erase previous files
    if (-e "Makefile.am")
      {
	`rm -f Makefile.am`;
      }

# look for subdirectories and generate makefile in them
    my $TmpMakefileAm = "SUBDIRS=";
    my @SubDir = &ListDir();
    my $FlagSubDir = 0;
    my $CurrentSubDir;
    my $FirstFlag = 0;
    foreach $CurrentSubDir (@SubDir)
      {
	chdir $CurrentSubDir;
	my $TmpSubDir = &GenerateMakefile($RootDir, $ExcludedPackages, $UsedPackages);
	if (-e "Makefile.am")
	  {
	    $FlagSubDir = 1;
	  }
	chdir "..";
	if ($TmpSubDir ne "")
	  {
	    if ($FirstFlag == 0)
	      {
		$FirstFlag = 1;
	      }
	    else
	      {
		$TmpMakefileAm .= " ";
	      }
	    $TmpMakefileAm .= $TmpSubDir;
	  }
     }

# find package name and test if it has to be built
    my $PackageName = $CurrentDir;
    $PackageName =~ s/.*\///g;
    my $ExcludedPackageName;
    foreach $ExcludedPackageName (@{$ExcludedPackages})
      {
	if ($PackageName eq $ExcludedPackageName)
	  {
	    return "";
	  }
      }

# test if a Makefile.am fas to be created
   my @SourceFiles = &ListSourceFile(); 
   if (($#SourceFiles < 0) && ($FlagSubDir == 0))
      {
	return "";
      }

# begin to create Makefile.am
    $TmpMakefileAm .= "\n\nnoinst_LIBRARIES=lib".$PackageName.".a\n
lib".$PackageName."_a_SOURCES=";
    $FirstFlag = 0;
    my $TmpSrcFile;
    foreach $TmpSrcFile (&ListSourceFile(