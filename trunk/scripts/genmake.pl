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

sub GenerateBaseSrcMakefile
  {
    # prepare data
    my $CurrentDir = `pwd`;
    chomp ($CurrentDir);
    my $RootDir = $_[0];
    my $Source = $_[1];
    my $ExecutableName = $_[2];
    my $ExcludedPackages = $_[3];

    # erase previous files
    if (-e "makefile")
      {
	`rm -f makefile`;
      }
    if (-e "Dependence")
      {
	`rm -f Dependence`;
      }
    my $TmpDepFile;
    foreach $TmpDepFile (&ListDependenceFile())
      {
	 `rm -f $TmpDepFile`;
      }

    # look for subdirectories and generate makefile in them
    my @SubDir = ListDir();
    my $FlagSubDir = 0;
    my $CurrentSubDir;
    foreach $CurrentSubDir (@SubDir)
      {
 	chdir $CurrentSubDir;
	GenerateMakefile($RootDir, $ExcludedPackages);
	if (-e "makefile")
	  {
	    $FlagSubDir = 1;
	  }
	chdir "..";
     }

    # find config path from current directory
    my $ConfigDir = $CurrentDir;
    $ConfigDir =~ s/$RootDir//;
    if ($ConfigDir =~ /\//)
      {
	my @ParentDir = split (/\//, $ConfigDir);
	my $NbrBack = $#ParentDir;
	$ConfigDir = "..";
	while ($NbrBack-- > 1)
	  {
	    $ConfigDir .= "/..";
	  }
      }
    else
      {
	$ConfigDir = ".";
      }

    # test if a makefile fas to be created
    my @SourceFiles = &ListSourceFile(); 
    if (($#SourceFiles < 0) && ($FlagSubDir == 0))
      {
	return;
      }

    # begin to create makefile
    my $Makefile = "";
    $Makefile .= "include $ConfigDir/config
EXE = \$(SrcDir)/$ExecutableName
SRC = $Source
OBJ = \$(SRC:.cc=.o)
include ./Dependence\n";	

    my $DepFile = $Source;
    &GenerateDependences($DepFile, $RootDir);
    $DepFile =~ s/\.cc$/\.d/;	
    $Makefile .= "include ./$DepFile\n";	
    
    &GenerateFullDependenceFile($RootDir);
    my %Packages;
    &GetPackages($BaseDirectory, \%Packages);
    my $LibraryLinkFlag = &GetLibraryFlag(\%Packages);

    $Makefile .= "\n\n\$(EXE): \$(OBJ) \$(DEP)";

    $Makefile .= "\n";
    if ($#SubDir >= 0)
      {
	my $Dir;
	foreach $Dir (@SubDir)
	  {
	    chdir ($Dir);
	    if (-e "makefile")
	      {
		$Makefile .= "\tcd $Dir && \$(MAKE)\n";
	      }
	    chdir ("..");
	  }
      }
    $Makefile .= "\n	\$(C++) \$(CCFLAGS) -o \$(EXE) \$(OBJ) -L\$(LibDir) \$(EXTLIBPATH) $LibraryLinkFlag \$(EXTLIB)\n\n";

    my $DepName = $Source;
    $DepName =~ s/\.cc$//;
    my $ObjName = $Source;
    $ObjName =~ s/\.cc$/\.o/;
    $Makefile .= "$ObjName : $Source \$(Dep$DepName)
\t\$(C++) \$(CCFLAGS) -I\$(SrcDir) -c $Source

";
    
    $Makefile .= "clean:\n";
    if ($#SubDir != 0)
      {
	my $Dir;
	foreach $Dir (@SubDir)
	  {
	    chdir ($Dir);
	    if (-e "makefile")
	      {
		$Makefile .="\tcd $Dir && \$(MAKE) clean\n";
	      }
	    chdir ("..");
	  }
      }
    $Makefile .= "	-rm \$(OBJ)	
	-rm \$(EXE)\n";
    
    # save makefile
    unless (open(OUTFILE, ">makefile"))
      {
	die "error : can't create makefile \n";
      }
    print OUTFILE $Makefile;
    close (OUTFILE);	
  }



# generate makefile for a given package
#
# $_[0] = path to the directory contenaining the config file
# $_[1] = reference to the list of excluded packages

sub GenerateMakefile
  {
# prepare data
    my $CurrentDir = `pwd`;
    chomp ($CurrentDir);
    my $RootDir = $_[0];
    my $ExcludedPackages = $_[1];

# erase previous files
    if (-e "makefile")
      {
	`rm -f makefile`;
      }
    if (-e "Dependence")
      {
	`rm -f Dependence`;
      }
    my $TmpDepFile;
    foreach $TmpDepFile (&ListDependenceFile())
      {
	 `rm -f $TmpDepFile`;
      }

# look for subdirectories and generate makefile in them
    my @SubDir = ListDir();
    my $FlagSubDir = 0;
    my $CurrentSubDir;
    foreach $CurrentSubDir (@SubDir)
      {
 	unless (chdir $CurrentSubDir)
	  {
	    die ("can't open $CurrentDir/$CurrentSubDir\n");
	  }
	GenerateMakefile($RootDir, $ExcludedPackages);
	if (-e "makefile")
	  {
	    $FlagSubDir = 1;
	  }
	chdir "..";
     }

# find config path from current directory
    my $ConfigDir = $CurrentDir;
    $ConfigDir =~ s/$RootDir//;
    if ($ConfigDir =~ /\//)
      {
	my @ParentDir = split (/\//, $ConfigDir);
	my $NbrBack = $#ParentDir;
	$ConfigDir = "..";
	while ($NbrBack-- > 1)
	  {
	    $ConfigDir .= "/..";
	  }
      }
    else
      {
	$ConfigDir = ".";
      }

# find package name and test if it has to be built
    my $PackageName = $CurrentDir;
    $PackageName =~ s/.*\///;
    my $ExcludedPackageName;
    foreach $ExcludedPackageName (@{$ExcludedPackages})
      {
	if ($PackageName eq $ExcludedPackageName)
	  {
	    return;
	  }
      }

# test if a makefile fas to be created
   my @SourceFiles = &ListSourceFile(); 
   if (($#SourceFiles < 0) && ($FlagSubDir == 0))
      {
	return;
      }

# begin to create makefile
    my $Makefile = "";
    $Makefile .= "include $ConfigDir/config
LIB = \$(LibDir)/lib$PackageName.a
SRC = ";
    foreach (@SourceFiles)
      {
	$Makefile .= $_." ";
      }
    $Makefile .="
OBJ = \$(SRC:.cc=.o)
";
    $Makefile .= "include ./Dependence\n\n";	

    foreach (@SourceFiles)
      {
	my $DepFile = $_;
	&GenerateDependences($DepFile, $RootDir);
	$DepFile =~ s/\.cc$/\.d/;	
	$Makefile .= "include ./$DepFile\n";	
      }
    
    &GenerateFullDependenceFile($RootDir);

    $Makefile .= "\n\n\$(LIB): \$(OBJ) \$(DEP)";

    $Makefile .= "\n";
    if ($#SubDir >= 0)
      {
	my $Dir;
	foreach $Dir (@SubDir)
	  {
	    chdir ($Dir);
	    if (-e "makefile")
	      {
		$Makefile .= "\tcd $Dir && \$(MAKE)\n";
	      }
	    chdir ("..");
	  }
      }
    $Makefile .= "\n	\$(AR) \$(ARFLAGS) \$(LIB) \$(OBJ)\n\n";

    foreach (@SourceFiles)
      {
	my $DepName = $_;
	$DepName =~ s/\.cc$//;
	my $ObjName = $_;
	$ObjName =~ s/\.cc$/\.o/;
	$Makefile .= "$ObjName : $_ \$(Dep$DepName)
\t\$(C++) \$(CCFLAGS) -I\$(SrcDir) \$(INCPATH) -c $_

";
      }

$Makefile .= "clean:\n";
  if ($#SubDir != 0)
    {
      my $Dir;
      foreach $Dir (@SubDir)
	{
	  chdir ($Dir);
	  if (-e "makefile")
	    {
	      $Makefile .="\tcd $Dir && \$(MAKE) clean\n";
	    }
	  chdir ("..");
	}
    }
  $Makefile .= "	-rm \$(OBJ)	
	-rm \$(LIB)\n";

# save makefile
    unless (open(OUTFILE, ">makefile"))
      {
	die "error : can't open makefile \n";
      }
    print OUTFILE $Makefile;
    close (OUTFILE);	
  }



# generate header dependence file of a given source file
#
# $_[0] = Source File Name
# $_[1] = source root directory

sub GenerateDependences
  {
    my $SourceFileName = $_[0];
    my $SourceDir = $_[1];
    my @ExplicitDepList = &GenerateExplicitSourceDepList($SourceFileName);
    my $Header;
    foreach $Header (@ExplicitDepList)
      {
#	print (" $Header \n");
	my @TmpDepList = &GenerateExplicitHeaderDepList($Header, $SourceDir);
	&MergeDepList(\@ExplicitDepList, \@TmpDepList);
      }
    my @DepList = @ExplicitDepList;
    my $DepFile = "Dep$SourceFileName";
    $DepFile =~ s/\.cc$//;
    $DepFile .= " =";
    my $Flag = 0;
    foreach (@DepList)
      {
	if ($Flag == 1)
	  {
	    $DepFile .= " \\\n";
	  }
	$DepFile .= "\t\$(SrcDir)/$_";
	$Flag = 1;
      }
    $DepFile .= "\n";
    $SourceFileName =~ s/\.cc$/\.d/;
    unless (open(OUTFILE, ">$SourceFileName"))
      {
	die "error : can't open $SourceFileName \n";
      }
    print OUTFILE $DepFile;
    close (OUTFILE);	    
  }



# generate explicit (from class source) dependence list
#
# $_[0] = Source File Name
# return value = List of all explicit header dependence of the source file

sub GenerateExplicitSourceDepList
  {
    my $SourceFileName = $_[0];
    unless (open(INFILE, $SourceFileName))
      {
	die "error : can't open  $SourceFileName\n";
      }
    my @ListDep;
    my $CurrentLine;
    while (defined($CurrentLine = <INFILE>))
      {
	if (($CurrentLine =~ /\#include *\"/) && !($CurrentLine =~ /\/\//))
	  {
	    $CurrentLine = substr($CurrentLine, index ($CurrentLine, "\"") + 1);				
	    $CurrentLine = substr($CurrentLine, 0, index ($CurrentLine, "\""));				
	    if ($CurrentLine ne "")
	      {
		push (@ListDep, $CurrentLine);
	      }
	  }
     }
    close (INFILE);
    return @ListDep;
  }



# generate explicit (from class header) dependence list
#
# $_[0] = Header File Name (with path from SrcDir)
# $_[1] = Source Root directory
# return value = List of all explicit header dependence of the header file

sub GenerateExplicitHeaderDepList
  {
    my $HeaderFileName = $_[0];
    $HeaderFileName =~ s/^.*\///;   
    my $Directory = $_[0];
    if ($Directory =~ /\//)
      {
	$Directory =~ s/\/[^\/]*$//;
	$Directory = $_[1]."/".$Directory; 	
      }
    else
      {
	$Directory = $_[1];
      }
    my $CurrentDir = `pwd`;
    chomp ($CurrentDir);   
    chdir $Directory;
    unless (open(INFILE, $HeaderFileName))
      {
	die "error : can't open  $Directory/$HeaderFileName\n";
      }
    my @ListDep;
    my $CurrentLine;
    while (defined($CurrentLine = <INFILE>))
      {
	if (($CurrentLine =~ /\#include *\"/) && !($CurrentLine =~ /\/\//))
	  {
	    $CurrentLine = substr($CurrentLine, index ($CurrentLine, "\"") + 1);				
	    $CurrentLine = substr($CurrentLine, 0, index ($CurrentLine, "\""));				
	    if ($CurrentLine ne "")
	      {
		push (@ListDep, $CurrentLine);
	      }
	  }
     }
    close (INFILE);
    chdir $CurrentDir;
#    my $Header;
#    print "$HeaderFileName :\n";
#    foreach $Header (@ListDep)
#      {
#	print "$Header\n";
#      }
#    print "\n\n";
    return @ListDep;
  }



# generate full dependence file for the current package
#
# $_[0] = source root directory

sub GenerateFullDependenceFile()
  {
    my $RootDir = $_[0];
    my @DepList;
    my @SubDir = ListDir();
    my $TmpDependenceFile = "DEP = ";
# find dependences in subdirectory 
    foreach (@SubDir)
      {
 	chdir $_;
	if (-e "Dependence")
	  {
	    unless (open (INFILE, "Dependence"))
	      {
		die ("error : can't open $_/Dependence\n");
	      }
	    my $CurrentLine;
	    while (defined($CurrentLine = <INFILE>))
	      {
		chomp ($CurrentLine);
		$CurrentLine =~ s/\\//;
		$CurrentLine =~ s/.*\$/\$/;		
		my $TmpDep;
		my $Flag = 0;
		foreach $TmpDep (@DepList)
		  {
		    if ($TmpDep eq $CurrentLine)
		      {
			$Flag = 1;
			last;
		      }
		  }
		if ($Flag == 0)
		  {
		    push(@DepList, $CurrentLine);
		  }
	      }
	    close (INFILE)
	  }
	chdir "..";
     }

# add header dependences of the current package 
    my $SrcFile;
    foreach $SrcFile (&ListDependenceFile())
      {
	unless (open (INFILE, $SrcFile))
	  {
	    my $Dir = `pwd`;
	    die ("error : can't open $Dir/$SrcFile\n");
	  }
	my $CurrentLine;
	while (defined($CurrentLine = <INFILE>))
	  {
	    chomp ($CurrentLine);
	    $CurrentLine =~ s/.*\$/\$/;
	    $CurrentLine =~ s/\\//;
	    my $TmpDep;
	    my $Flag = 0;
	    foreach $TmpDep (@DepList)
	      {
		if ($TmpDep eq $CurrentLine)
		  {
		    $Flag = 1;
		    last;
		  }
	      }
	    if ($Flag == 0)
	      {
		push(@DepList, $CurrentLine);
	      }
	  }
      }

#prepare file
    my $RelDir = `pwd`;
    chomp($RelDir);
    $RelDir =~ s/^$RootDir//;
    foreach $SrcFile (&ListSourceFile())
      {
	$TmpDependenceFile .= "\t\$(SrcDir)$RelDir/$SrcFile \\\n"
      }
    foreach $SrcFile (@DepList)
      {
	$TmpDependenceFile .= "\t".$SrcFile."\\\n"
      }
    chomp ($TmpDependenceFile);
    $TmpDependenceFile =~ s/\\$/\n/;

#write file
    my $Dir2 = `pwd`;
    print "$Dir2\n";
    unless (open(OUTFILE, ">Dependence"))
      {
	my $Dir = `pwd`;
	die "error : can't open $Dir/Dependence \n";
      }
    print OUTFILE $TmpDependenceFile;
    close (OUTFILE);	    
 
  }



# merge list $_[1] into $_[0] with no duplication
#
# $_[0] = reference to destination list
# $_[1] = reference to the list to merge

sub MergeDepList
  {    
    my $CurrentElement;
    foreach $CurrentElement(@{$_[1]})
      {
	my $flag = 0;
	foreach (@{$_[0]})
	  {
	    if ($_ eq $CurrentElement)
	      {
		$flag = 1;
		last;
	      }
	  }
	if ($flag == 0)
	  {
	    push (@{$_[0]}, $CurrentElement);
	  }
      }
  }

# find all packages and their dependences
#
# $_[0] = base directory
# $_[1] = reference on the hqsh table to use 

sub GetPackages()
  {
    my $BaseDir = $_[0];
    my $Packages = $_[1];
    my $SubDir;
    foreach $SubDir (&ListDir())
      {
	chdir $SubDir;
	&GetPackages($BaseDir, \%{$Packages});
	chdir ("..");
      }
    my $CurrentDir = `pwd`;
    chomp ($CurrentDir);
    if ($CurrentDir eq $BaseDir)
      {
	return;
      }
    my $PackageName = $CurrentDir;
    $PackageName =~ s/^.*\///;
    print "$PackageName\n";
    my @PackageDependences = &FindPackageDependences($BaseDir);
    if ($#PackageDependences >= 0)
      {
	%{$Packages} = (%{$Packages}, $PackageName, \@PackageDependences);
      }
  }



# find list of all package dependences for a given package
#
# $_[0] = base directory
# return value = list of package dependences

sub FindPackageDependences()
  {
    my $BaseDir = $_[0];
    my @ListDep;
    my $DepFile;
    my $CurrentDir = `pwd`;
    chomp ($CurrentDir);
    my $PackageName = $CurrentDir;
    $PackageName =~ s/^.*\///;
    foreach $DepFile (&ListDependenceFile())
      {
	unless (open(INFILE, $DepFile))
	  {	    
	    die ("error : can't open $CurrentDir/$DepFile\n");
	  }
	my $CurrentLine;
	while ($CurrentLine = <INFILE>)
	  {
	    chomp ($CurrentLine);
	    $CurrentLine =~ s/.*\$\(SrcDir\)//;
	    $CurrentLine =~ s/ *\\$//;
	    $CurrentLine =~ s/\/[^\/]*$//;
	    $CurrentLine =~ s/^\///;
	    if ($CurrentLine ne "")
	      {
		my $Flag = 0;
		my $TmpPackage;
		foreach $TmpPackage (@ListDep)
		  {
		    if ($TmpPackage eq $CurrentLine)
		      {
			$Flag = 1;
			last;
		      }
		  }
		if ($Flag == 0)
		  {		
		    push (@ListDep, $CurrentLine);
		  }
	      }
	  }
	close (INFILE);
      }

# check if package is really a librairy
    my @RealDepList;
    my $TmpPackage;
    foreach $TmpPackage (@ListDep)
      {
	if (-e "$BaseDir/$TmpPackage/makefile")
	  {
	    $TmpPackage =~ s/^.*\///;
	    push (@RealDepList, $TmpPackage);
	  }
      }
    return @RealDepList;
  }



# return library flag for linker
#
# $_[0] = reference to te hash table of packages and their dependences
# return value = library flag for linker

sub GetLibraryFlag
  {
    my %Packages = %{$_[0]};
    my @Libraries;
    
    # find all package with no external dependence
    my $PackageName;
    my $PackageDependences;
    while (($PackageName, $PackageDependences) = each (%Packages))
      {
	if ($#{@{$PackageDependences}} == 0)
	  {
	    push (@Libraries, $PackageName);
	    delete $Packages{$PackageName};
	  }
      }

    # find maximal number of dependences
    my $MaxNbrDep = 0;
    while (($PackageName, $PackageDependences) = each (%Packages))
      {
	if ($#{@{$PackageDependences}} > $MaxNbrDep)
	  {
	    $MaxNbrDep = $#{@{$PackageDependences}};
	  }
      }

    # find packages with no loop dependences
    my $CurrentNbrDep =1;
    my $MainFlag = 0;
    while ($CurrentNbrDep <= $MaxNbrDep)
      {
	while (($MainFlag == 0) && (($PackageName, $PackageDependences) = each (%Packages)))
	  {
	    if ($#{@{$PackageDependences}} ==  $CurrentNbrDep)
	      {
		my $CurrentDep;
		my $Flag = 0;
		foreach $CurrentDep (@{$PackageDependences})
		  {
		    if (($CurrentDep ne $PackageName) && (GetPosition(\@Libraries, $CurrentDep) == -1))
		      {
			$Flag = 1;
			last;
		      }
		  }
		if ($Flag == 0)
		  {
		    $MainFlag = 1;
		    push (@Libraries, $PackageName);
		    delete $Packages{$PackageName};		    
		  }
	      }
	  }
	if ($MainFlag == 1)
	  {
	    $MainFlag = 0;
	  }
	else
	  {
	    $CurrentNbrDep++;
	  }
      }
    
    # sort left packages with respect to the number of non-sorted dependences 
    my @TmpList = keys(%Packages);
    while ($#TmpList >= 0)
      {
	my $MaxOccurrence = -1;
	my @MaxOccurrenceList;
	foreach $PackageName (@TmpList)
	  {
	    my $TmpOccurrence = &GetPackageOccurrence(\%Packages, $PackageName);
	    print "$TmpOccurrence \n";
	    if ($TmpOccurrence > $MaxOccurrence)
	      {
		@MaxOccurrenceList = ($PackageName);
		$MaxOccurrence = $TmpOccurrence;
	      }
	    else
	      {
		if ($TmpOccurrence == $MaxOccurrence)
		  {
		    push (@MaxOccurrenceList, $PackageName);
		  }
	      }
	  }
	print "$MaxOccurrence @MaxOccurrenceList\n";
	my $NbrDepTable = &GetNbrDependences(\@Libraries, \%Packages);
	my $MinDep = $#{@{$Packages{$MaxOccurrenceList[0]}}} - ${%{$NbrDepTable}}{$MaxOccurrenceList[0]};
	my $MinPackageName = $MaxOccurrenceList[0];
	foreach $PackageName (@MaxOccurrenceList)
	  {
	    my $TmpDep = $#{@{$Packages{$PackageName}}} - ${%{$NbrDepTable}}{$PackageName};
	    if ($TmpDep < $MinDep)
	      {
		$MinDep = $TmpDep;
		$MinPackageName = $PackageName;
	      }
	  }
	push (@Libraries, $MinPackageName);
	delete $Packages{$MinPackageName};		    	
	@TmpList = keys(%Packages);
      }

    # write flag
    my $LinkerFlag = "";
    my $TmpLibrary;
    foreach $TmpLibrary (@Libraries)
      {
#	$LinkerFlag = "\$(LibDir)/$TmpLibrary.a ".$LinkerFlag;
	$LinkerFlag = "-l$TmpLibrary ".$LinkerFlag;
      }
    $LinkerFlag =~ s/ $//;
    return $LinkerFlag;
  }



# create a hash table with number of remaining dependences for each package
#
# $_[0] = reference to the list of already use packages
# $_[1] = reference to the hash table of packages and their dependences
# return value = reference to hash table with number of remaining dependences for each package

sub GetNbrDependences()
  {
    my %NbrDep;
    my $PackageName;
    my $PackageDependences;
    while (($PackageName, $PackageDependences) = each (%{$_[1]}))
      {
	my $TmpNbr = 0;
	my $TmpPackageName;
	foreach $TmpPackageName (@{$PackageDependences})
	  {
	    if (GetPosition($_[0], $TmpPackageName) != -1)
	      {
		$TmpNbr++;
	      }
	  }
	%NbrDep = (%NbrDep, $PackageName, $TmpNbr);
      }
    return \%NbrDep;
  }


# get occurence of package in a hash table of package dependences
#
# $_[0] = hash table reference
# $_[1] = package to check
# return value = package occurence

sub GetPackageOccurrence()
  {
    my $Occurrence = 0;
    my $PackageName;
    my $PackageDependences;
    while (($PackageName, $PackageDependences) = each (%{$_[0]}))
      {
	if (GetPosition(\@{$PackageDependences}, $_[1]) != -1)
	    {
	      $Occurrence++;
	      next;
	    }
      }
    return $Occurrence;
  }



# return position of an element in a list of strings
#
# $_[0] = list reference
# $_[1] = string to check
# return value = element position in the list (-1 if it doesn't belong to)

sub GetPosition()
  {
    my $TmpString;
    my $Pos = 0;    
    foreach $TmpString (@{$_[0]})
      {
	if ($TmpString eq $_[1])
	  {
	    last;
	  }
	else
	  {
	    $Pos++;
	  }
      }
    if ($Pos == @{$_[0]})
      {
	return -1;
      }
    else
      {
	return $Pos;
      }
  }
